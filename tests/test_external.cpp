// atomforge external-potential tests: the native Lennard-Jones potential, the ExternalPotential bridge against it (through
// quantum/worker.py lj), FIRE relaxation through the bridge, and the bridge's failure paths.
//
//   test_external <group>    groups: lj external external_relax external_errors
//   test_external external_search  the binary search pipeline through the bridge vs the native LJ (reduced)
//   test_external external_cli <path-to-atomforge>   `atomforge relax --potential external:...` end to end
//
// Groups that need python3 exit with status 77 (CTest SKIP_RETURN_CODE) when it is not on PATH.
#include "atomforge/atomforge.hpp"
#include "atomforge/external.hpp"
#include "atomforge/lj.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <clocale>
#include <csignal>
#include <filesystem>
#include <unistd.h>
#include <sys/wait.h>
#include <filesystem>
#include <fstream>
#include <sstream>

#ifndef ATOMFORGE_SOURCE_DIR
#define ATOMFORGE_SOURCE_DIR "."
#endif

using namespace atomforge;
static int g_fail = 0, g_pass = 0;
static bool g_skip = false;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static const std::string kWorker = std::string(ATOMFORGE_SOURCE_DIR) + "/quantum/worker.py";
static LJParams make_lj(bool per_species) {
    LJParams p; p.epsilon = 0.0104; p.sigma = 3.40; p.rc = 8.5; p.mass = 39.948;
    if (per_species) { p.epsilon_by_species = {0.0104, 0.0141}; p.sigma_by_species = {3.40, 3.65}; }   // Ar-like, Kr-like
    return p;
}
static const LJParams kLJ = make_lj(false);
static const LJParams kLJ2 = make_lj(true);   // species-dependent: catches species mix-ups in the bridge

static std::string lj_command(const std::string& extra = "", const LJParams& p = kLJ) {
    char b[320];
    if (p.epsilon_by_species.empty())
        std::snprintf(b, sizeof b, " lj --epsilon %.17g --sigma %.17g --rc %.17g --mass %.17g", p.epsilon, p.sigma, p.rc, p.mass);
    else
        std::snprintf(b, sizeof b, " lj --epsilon %.17g,%.17g --sigma %.17g,%.17g --rc %.17g --mass %.17g", p.epsilon_by_species[0],
                      p.epsilon_by_species[1], p.sigma_by_species[0], p.sigma_by_species[1], p.rc, p.mass);
    return "python3 '" + kWorker + "'" + b + (extra.empty() ? "" : " " + extra);
}
// Exact-PID liveness (never pgrep -f: that matches unrelated command lines). A zombie counts as dead.
static bool pid_alive(long pid) {
    if (pid <= 0) return false;
    std::FILE* f = std::fopen(("/proc/" + std::to_string(pid) + "/stat").c_str(), "r");
    if (!f) return ::kill(pid_t(pid), 0) == 0 && !std::filesystem::exists("/proc/self");   // no procfs: kill(0) fallback
    char buf[512] = {0};
    const std::size_t n = std::fread(buf, 1, sizeof buf - 1, f);
    std::fclose(f);
    const char* rp = std::strrchr(buf, ')');
    return n > 0 && rp && rp[1] == ' ' && rp[2] != 'Z' && rp[2] != 'X';
}
static bool wait_dead(long pid, double secs) {
    for (int k = 0; k < int(secs * 20); ++k) { if (!pid_alive(pid)) return true; ::usleep(50000); }
    return !pid_alive(pid);
}
static long read_pid(const std::string& path) {
    for (int k = 0; k < 100; ++k) {   // up to 5 s for the process to write it
        std::FILE* f = std::fopen(path.c_str(), "r");
        long pid = 0;
        if (f) { if (std::fscanf(f, "%ld", &pid) != 1) pid = 0; std::fclose(f); }
        if (pid > 0) return pid;
        ::usleep(50000);
    }
    return 0;
}
static bool have_python() {
    if (std::system("command -v python3 >/dev/null 2>&1") == 0) return true;
    std::puts("  SKIP: python3 not found on PATH");
    g_skip = true;
    return false;
}
static std::string log_path(const char* tag) {
    const char* t = std::getenv("TMPDIR");
    return std::string(t && *t ? t : "/tmp") + "/atomforge-test-" + tag + ".log";
}

static Structure proto(const std::string& name) {
    for (const auto& p : prototypes()) if (p.name.rfind(name, 0) == 0) return p.build(0, 1);
    throw std::runtime_error("no prototype " + name);
}
static Structure rattle(Structure s, double amp, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(-amp, amp);
    for (auto& r : s.pos) r += Vec3{U(rng), U(rng), U(rng)};
    return s;
}
// A general triclinic cell: L1_1 (rhombohedral) sheared further, rattled.
static Structure triclinic(unsigned seed) {
    Structure s = proto("L1_1").supercell(2, 1, 1);
    s.deform(Mat3::diag(5.4, 5.4, 5.4));
    Mat3 F = Mat3::identity(); F(0, 1) = 0.11; F(1, 2) = -0.07; F(2, 0) = 0.05; F(0, 0) = 1.03;
    s.deform(F);
    return rattle(s, 0.12, seed);
}
static Structure fcc_rattled(unsigned seed) {
    Structure s = proto("A1").supercell(2, 2, 2);
    s.deform(Mat3::diag(5.30, 5.30, 5.30));
    return rattle(s, 0.15, seed);
}

// ---------------------------------------------------------------- native LJ derivatives
static void test_lj() {
    std::puts("[lj] native Lennard-Jones: forces and virial vs finite differences, minimum position");
    LennardJones pot({"Ar", "Kr"}, kLJ);
    for (const Structure& s : {fcc_rattled(1), triclinic(2)}) {
        const auto ev = evaluate(s, pot);
        const double h = 1e-5;
        double worst = 0, fmax = 0;
        Vec3 sum{};
        for (const auto& f : ev.forces) { sum += f; fmax = std::max(fmax, norm(f)); }
        for (std::size_t i = 0; i < s.size(); i += 3) for (int k = 0; k < 3; ++k) {
            Structure p = s, m = s; p.pos[i][k] += h; m.pos[i][k] -= h;
            const double fd = -(evaluate(p, pot).energy - evaluate(m, pot).energy) / (2 * h);
            worst = std::max(worst, std::abs(fd - ev.forces[i][k]));
        }
        double wworst = 0, wmax = 0;
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
            Mat3 ep = Mat3::identity(), em = Mat3::identity();
            ep(a, b) += h / 2; ep(b, a) += h / 2; em(a, b) -= h / 2; em(b, a) -= h / 2;
            if (a == b) { ep(a, a) = 1 + h; em(a, a) = 1 - h; }
            Structure p = s, m = s; p.deform(ep); m.deform(em);
            const double fd = (evaluate(p, pot).energy - evaluate(m, pot).energy) / (2 * h);
            wworst = std::max(wworst, std::abs(fd - ev.virial(a, b)));
            wmax = std::max(wmax, std::abs(ev.virial(a, b)));
        }
        std::printf("  N=%-3zu max|F|=%.3f eV/A |F-FD|=%.2e |sumF|=%.2e  max|W|=%.3f eV |W-FD|=%.2e\n", s.size(), fmax, worst, norm(sum), wmax, wworst);
        CHECK(fmax > 0.01, "structure not displaced enough (max|F| %.3e)", fmax);
        CHECK(worst < 1e-7, "LJ force vs FD error %.3e", worst);
        CHECK(norm(sum) < 1e-12, "LJ net force %.3e", norm(sum));
        CHECK(wworst < 1e-7 * std::max(1.0, wmax), "LJ virial vs finite strain error %.3e", wworst);
    }
    // Dimer in a big box: minimum at 2^(1/6) sigma, depth eps (less the cutoff shift).
    Structure d; d.cell = Mat3::diag(30, 30, 30);
    const double rmin = std::pow(2.0, 1.0 / 6) * kLJ.sigma;
    d.pos = {{0, 0, 0}, {rmin, 0, 0}}; d.species = {0, 0};
    const auto ev = evaluate(d, pot);
    const double sr6 = std::pow(kLJ.sigma / kLJ.rc, 6), shift = 4 * kLJ.epsilon * (sr6 * sr6 - sr6);
    CHECK(std::abs(ev.energy - (-kLJ.epsilon - shift)) < 1e-14, "dimer energy %.15g", ev.energy);
    CHECK(norm(ev.forces[0]) < 1e-14, "dimer force at r_min %.3e", norm(ev.forces[0]));
}

// ---------------------------------------------------------------- bridge vs native
static void test_external() {
    std::puts("[external] ExternalPotential(worker.py lj) vs native LennardJones: energy, forces, virial");
    if (!have_python()) return;
    LennardJones nat({"Ar", "Kr"}, kLJ);
    ExternalOptions o; o.log_path = log_path("external"); o.hello_timeout_s = 30; o.eval_timeout_s = 60;
    ExternalPotential ext(lj_command(), {"Ar", "Kr"}, o);
    CHECK(ext.name().find("external:lj(") == 0, "name '%s'", ext.name().c_str());
    CHECK(std::abs(ext.cutoff() - kLJ.rc) < 1e-12, "cutoff %.6g", ext.cutoff());
    CHECK(std::abs(ext.mass(1) - kLJ.mass) < 1e-12, "mass %.6g", ext.mass(1));
    struct Case { const char* name; Structure s; };
    Structure shifted = fcc_rattled(5);
    shifted.pos[3] += shifted.cell.row(0) * 2.0 - shifted.cell.row(2) * 1.0;   // unwrapped atom
    const Case cases[] = {{"fcc 2x2x2 rattled (32 atoms)", fcc_rattled(3)},
                          {"triclinic L1_1 2x1x1 rattled (4 atoms)", triclinic(4)},
                          {"fcc rattled, one atom 2 cells outside", shifted},
                          {"fcc 2x2x2 rattled, 2 species", [] { Structure s = fcc_rattled(6); for (std::size_t i = 0; i < s.size(); i += 2) s.species[i] = 1; return s; }()}};
    for (const auto& c : cases) {
        NeighborList nl;
        const auto a = nat.evaluate(c.s, nl);
        const auto b = evaluate(c.s, ext);
        double df = 0, dw = 0, fscale = 0, wscale = 0;
        CHECK(b.forces.size() == c.s.size(), "%s: %zu forces", c.name, b.forces.size());
        for (std::size_t i = 0; i < c.s.size() && i < b.forces.size(); ++i) { df = std::max(df, norm(a.forces[i] - b.forces[i])); fscale = std::max(fscale, norm(a.forces[i])); }
        for (int k = 0; k < 9; ++k) { dw = std::max(dw, std::abs(a.virial.m[k] - b.virial.m[k])); wscale = std::max(wscale, std::abs(a.virial.m[k])); }
        const double de = std::abs(a.energy - b.energy);
        std::printf("  %-42s E=%.6f eV |dE|=%.1e  max|F|=%.3f |dF|=%.1e  max|W|=%.3f |dW|=%.1e\n", c.name, a.energy, de, fscale, df, wscale, dw);
        CHECK(fscale > 0.01 && wscale > 0.01, "%s: trivial test structure", c.name);
        CHECK(de < 1e-9, "%s: energy differs by %.3e eV", c.name, de);
        CHECK(df < 1e-9, "%s: forces differ by %.3e eV/A", c.name, df);
        CHECK(dw < 1e-9, "%s: virial differs by %.3e eV", c.name, dw);
    }
    CHECK(ext.evaluations() == 4 && ext.alive(), "worker should have served 4 requests and still be alive (%ld)", ext.evaluations());

    // Species-dependent parameters on both sides (Lorentz-Berthelot), so a species mix-up in the bridge would show; once
    // with "virial" replies and once with "stress" replies (W/V, the form the pyscf and mace backends send).
    LennardJones nat2({"Ar", "Kr"}, kLJ2);
    auto two = [](Structure x) { for (std::size_t i = 0; i < x.size(); i += 2) x.species[i] = 1; return x; };
    const Structure mixed[] = {two(fcc_rattled(6)), two(triclinic(9))};
    for (const char* mode : {"virial", "stress"}) {
        ExternalOptions o2 = o; o2.log_path = log_path(mode);
        ExternalPotential e2(lj_command(std::string("--reply ") + mode, kLJ2), {"Ar", "Kr"}, o2);
        for (const auto& st : mixed) {
            const auto a = evaluate(st, nat2), b = evaluate(st, e2);
            double df = 0, dw = 0;
            for (std::size_t i = 0; i < st.size(); ++i) df = std::max(df, norm(a.forces[i] - b.forces[i]));
            for (int k = 0; k < 9; ++k) dw = std::max(dw, std::abs(a.virial.m[k] - b.virial.m[k]));
            std::printf("  two species, per-pair LJ, %-6s replies, N=%-3zu E=%.6f |dE|=%.1e |dF|=%.1e |dW|=%.1e\n", mode, st.size(), a.energy,
                        std::abs(a.energy - b.energy), df, dw);
            CHECK(std::abs(a.energy - b.energy) < 1e-9 && df < 1e-9 && dw < 1e-9, "%s replies: dE %.3e dF %.3e dW %.3e", mode, a.energy - b.energy, df, dw);
        }
        // the species really matter: swapping them changes the energy
        Structure sw = mixed[0]; for (auto& sp : sw.species) sp = 1 - sp;
        CHECK(std::abs(evaluate(sw, nat2).energy - evaluate(mixed[0], nat2).energy) > 1e-3, "per-species LJ does not depend on species");
    }

    // A worker whose backend prints to stdout (Python print and raw writes to fd 1, at start-up and per request, as
    // mace_mp does) still speaks clean protocol: worker.py moves the channel off fd 1.
    {
        ExternalOptions on = o; on.log_path = log_path("noisy"); std::remove(on.log_path.c_str());
        ExternalPotential noisy(lj_command("--noisy"), {"Ar", "Kr"}, on);
        const Structure st = fcc_rattled(3);
        const double d = std::abs(evaluate(st, noisy).energy - evaluate(st, nat).energy);
        std::ifstream lf(on.log_path); std::stringstream ls; ls << lf.rdbuf();
        CHECK(d < 1e-9, "noisy worker: energy differs by %.3e", d);
        CHECK(ls.str().find("noisy: raw fd 1 write during evaluate") != std::string::npos, "noisy output did not reach the log");
    }
}

// ---------------------------------------------------------------- FIRE through the bridge
static void test_external_relax() {
    std::puts("[external_relax] FIRE (positions + full cell) through the bridge matches FIRE on the native potential");
    if (!have_python()) return;
    LennardJones nat({"Ar"}, kLJ);
    ExternalOptions o; o.log_path = log_path("external_relax"); o.eval_timeout_s = 60;
    ExternalPotential ext(lj_command(), {"Ar"}, o);
    Structure s = proto("A1");
    s.deform(Mat3::diag(5.6, 5.6, 5.6));
    Mat3 F = Mat3::identity(); F(0, 1) = 0.04; F(2, 2) = 0.97;
    s.deform(F);
    s = rattle(s, 0.08, 7);
    for (auto& sp : s.species) sp = 0;
    const RelaxOptions ro{CellMode::Full, 1e-4, 2000, 0.0, 0.15};
    const auto t0 = std::chrono::steady_clock::now();
    const auto rn = relax(s, nat, ro);
    const auto re = relax(s, ext, ro);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double a = std::cbrt(re.s.volume() * 4 / double(re.s.size()));
    std::printf("  native: %d steps E=%.8f eV/atom   external: %d steps E=%.8f eV/atom  V=%.4f A^3/atom (a_cubic=%.4f A)  %ld worker calls, %.2f s\n",
                rn.steps, rn.energy / s.size(), re.steps, re.energy / s.size(), re.s.volume() / s.size(), a, ext.evaluations(), secs);
    CHECK(re.converged, "relaxation through the bridge did not converge (%d steps)", re.steps);
    CHECK(rn.steps == re.steps, "native and bridged FIRE took different paths (%d vs %d steps)", rn.steps, re.steps);
    CHECK(std::abs(rn.energy - re.energy) < 1e-9, "relaxed energies differ by %.3e eV", rn.energy - re.energy);
    CHECK(std::abs(rn.s.volume() - re.s.volume()) < 1e-6, "relaxed volumes differ by %.3e A^3", rn.s.volume() - re.s.volume());
    // relaxed residuals, independently
    const auto ev = evaluate(re.s, ext);
    double fm = 0; for (const auto& f : ev.forces) fm = std::max(fm, norm(f));
    const double pGPa = (ev.virial(0, 0) + ev.virial(1, 1) + ev.virial(2, 2)) / 3 / re.s.volume() * units::eV_per_A3_to_GPa;
    std::printf("  residual max|F| = %.2e eV/A, mean stress = %.2e GPa\n", fm, pGPa);
    CHECK(fm < 1e-3, "residual force %.3e", fm);
    CHECK(std::abs(pGPa) < 1e-2, "residual stress %.3e GPa", pGPa);
    CHECK(re.energy / s.size() < -0.06, "relaxed LJ fcc energy per atom %.4f eV is not bound", re.energy / s.size());
}

// ---------------------------------------------------------------- failure paths
template <class Fn>
static void expect_error(const char* what, Fn&& fn, const char* needle, double max_s = 20) {
    const auto t0 = std::chrono::steady_clock::now();
    bool threw = false; std::string msg;
    try { fn(); } catch (const ExternalPotentialError& e) { threw = true; msg = e.what(); }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string first = msg.substr(0, msg.find('\n'));
    std::printf("  %-34s %.2f s  -> %s\n", what, secs, threw ? first.c_str() : "(no exception)");
    CHECK(threw, "%s: no ExternalPotentialError", what);
    CHECK(msg.find(needle) != std::string::npos, "%s: message lacks '%s': %s", what, needle, msg.c_str());
    CHECK(secs < max_s, "%s: took %.1f s", what, secs);
}
static void test_external_errors() {
    std::puts("[external_errors] missing worker, crash, timeout, NaN, wrong shape, bad JSON, wrong id: clear exceptions, no hang");
    // JSON parser corner cases (no python3 needed)
    CHECK(std::isnan(json::parse("{\"e\": NaN}").get("e")->num), "NaN token");
    CHECK(json::parse("[1, -2.5e-3, \"a\\\"b\", true, null]").arr.size() == 5, "array parse");
    auto rejects = [](const std::string& t) { try { (void)json::parse(t); } catch (const ExternalPotentialError&) { return true; } return false; };
    CHECK(rejects("{\"e\": 1,}"), "trailing comma must be rejected");
    CHECK(rejects(std::string("[1\0 5]", 6)), "an embedded NUL must be rejected");
    CHECK(rejects("[1.5.2]") && rejects("[01x]"), "malformed numbers must be rejected");
    CHECK(json::num(0.1) == "0.1" && json::parse("[" + json::num(1.0 / 3) + "]").arr[0].num == 1.0 / 3, "round trip of doubles");
    // Locale independence: with a decimal-comma locale, numbers are still written and read with a point.
    if (std::setlocale(LC_ALL, "de_DE.UTF-8") || std::setlocale(LC_ALL, "de_DE") || std::setlocale(LC_ALL, "fr_FR.UTF-8")) {
        CHECK(json::num(8.5) == "8.5" && json::parse("[8.5]").arr[0].num == 8.5, "number formatting depends on the locale");
        std::setlocale(LC_ALL, "C");
    } else std::puts("  (no decimal-comma locale installed: locale check limited to the C locale)");
    if (!have_python()) return;
    const Structure s = fcc_rattled(8);
    ExternalOptions o; o.hello_timeout_s = 20; o.eval_timeout_s = 20;
    auto with = [&](const char* tag) { ExternalOptions x = o; x.log_path = log_path(tag); std::remove(x.log_path.c_str()); return x; };

    expect_error("worker missing", [&] { ExternalPotential p("/nonexistent/atomforge-worker --x", {"Ar"}, with("missing")); }, "exit status 127");
    expect_error("worker missing (stderr in message)", [&] { ExternalPotential p("/nonexistent/atomforge-worker", {"Ar"}, with("missing2")); }, "worker stderr");
    expect_error("worker exits before handshake", [&] { ExternalPotential p("sh -c 'echo boom >&2; exit 4'", {"Ar"}, with("early")); }, "exit status 4");
    // The worker exits but a process it started keeps the socket open: noticed by the exit, not by end of file.
    expect_error("exit with a child holding the socket", [&] { ExternalPotential p("sh -c 'sleep 47.123 & exit 4'", {"Ar"}, with("holder")); }, "exit status 4", 5);
    expect_error("unusable log path", [&] { ExternalOptions x = o; x.log_path = "/nonexistent-atomforge-dir/w.log"; ExternalPotential p(lj_command(), {"Ar"}, x); }, "cannot open worker log");
    expect_error("silent worker (handshake timeout)", [&] { ExternalOptions x = with("silent"); x.hello_timeout_s = 1;
                                                             ExternalPotential p("sleep 30", {"Ar"}, x); }, "timed out", 5);
    expect_error("wrong protocol version", [&] { ExternalPotential p("sh -c 'read l; echo \"{\\\"protocol\\\": 99}\"; sleep 5'", {"Ar"}, with("proto")); }, "protocol 99");

    {
        ExternalPotential p(lj_command("--fault crash --fault-after 1"), {"Ar"}, with("crash"));
        evaluate(s, p);   // first request is fine
        expect_error("worker crashes mid-run", [&] { evaluate(s, p); }, "exit status 3");
        CHECK(!p.alive(), "crashed worker still marked alive");
        expect_error("call after a crash", [&] { evaluate(s, p); }, "not running");
    }
    {
        ExternalPotential p(lj_command("--fault spawn-exit"), {"Ar"}, with("spawnexit"));
        expect_error("exit mid-request, child keeps socket", [&] { evaluate(s, p); }, "exit status 5", 5);
    }
    {
        ExternalOptions x = with("hang"); x.eval_timeout_s = 1.5;
        ExternalPotential p(lj_command("--fault hang"), {"Ar"}, x);
        expect_error("worker hangs (eval timeout 1.5 s)", [&] { evaluate(s, p); }, "timed out", 6);
    }
    {
        ExternalOptions x = with("inf"); x.eval_timeout_s = HUGE_VAL; x.hello_timeout_s = 1e300;   // no deadline, no overflow
        ExternalPotential p(lj_command(), {"Ar"}, x);
        CHECK(evaluate(s, p).forces.size() == s.size(), "infinite timeouts");
        Structure bad = s; bad.pos[2].y = std::nan("");
        expect_error("NaN position (not sent)", [&] { evaluate(bad, p); }, "not finite");
        CHECK(p.alive(), "a rejected input must not kill the worker");
        Structure badc = s; badc.cell(1, 1) = HUGE_VAL;
        expect_error("infinite cell (not sent)", [&] { evaluate(badc, p); }, "not finite");
        CHECK(evaluate(s, p).forces.size() == s.size(), "worker still serves after rejected inputs");
    }
    {
        ExternalPotential p(lj_command("--fault errnull"), {"Ar"}, with("errnull"));
        bool ok = true;
        try { evaluate(s, p); } catch (const std::exception& e) { ok = false; std::printf("  %s\n", e.what()); }
        CHECK(ok, "a reply with \"error\": null must count as success");
    }
    {
        ExternalPotential p(lj_command("--fault nan"), {"Ar"}, with("nan"));
        expect_error("NaN energy", [&] { evaluate(s, p); }, "not finite");
    }
    {
        ExternalPotential p(lj_command("--fault shape"), {"Ar"}, with("shape"));
        expect_error("too few force vectors", [&] { evaluate(s, p); }, "must have 32 entries, got 31");
    }
    {
        ExternalPotential p(lj_command("--fault garbage"), {"Ar"}, with("garbage"));
        expect_error("non-JSON reply", [&] { evaluate(s, p); }, "malformed JSON");
    }
    {
        ExternalPotential p(lj_command("--fault badid"), {"Ar"}, with("badid"));
        expect_error("reply with the wrong id", [&] { evaluate(s, p); }, "does not match");
    }
    {
        // A worker-side error is reported, and the worker keeps serving.
        ExternalPotential p(lj_command("--fault error"), {"Ar"}, with("werr"));
        expect_error("worker-side error", [&] { evaluate(s, p); }, "injected worker-side error");
        bool ok = true;
        try { evaluate(s, p); } catch (const std::exception& e) { ok = false; std::printf("  %s\n", e.what()); }
        CHECK(ok && p.alive(), "worker should keep serving after a reported error");
    }
    {
        // stderr goes to the log file
        ExternalOptions x = with("log");
        { ExternalPotential p(lj_command(), {"Ar"}, x); evaluate(s, p); }
        std::FILE* f = std::fopen(x.log_path.c_str(), "r");
        std::string txt; char buf[256];
        while (f && std::fgets(buf, sizeof buf, f)) txt += buf;
        if (f) std::fclose(f);
        CHECK(txt.find("[worker] started lj(") != std::string::npos, "worker stderr not in log %s: '%s'", x.log_path.c_str(), txt.c_str());
    }
    // A failed worker's children die with it (the worker runs in its own process group). Checked by exact PID.
    {
        const std::string pf = log_path("childpid") + ".pid";
        std::remove(pf.c_str());
        expect_error("worker with a background child", [&] {
            ExternalPotential p("sh -c 'sleep 7.3125 & echo $! > " + pf + "; read l; echo garbage; sleep 9'", {"Ar"}, with("child")); }, "malformed JSON");
        const long child = read_pid(pf);
        CHECK(child > 0 && wait_dead(child, 3), "the worker's background child (pid %ld) survived the worker", child);
        std::remove(pf.c_str());
    }
    // A worker killed without a word: the message quotes this instance's stderr only, not an earlier run's.
    {
        ExternalOptions x = with("sections");
        { std::FILE* f = std::fopen(x.log_path.c_str(), "w"); if (f) { std::fputs("OLD-RUN-MARKER\n", f); std::fclose(f); } }
        std::string msg;
        try { ExternalPotential p("sh -c 'echo NEW-RUN-MARKER >&2; exit 3'", {"Ar"}, x); } catch (const ExternalPotentialError& e) { msg = e.what(); }
        CHECK(msg.find("NEW-RUN-MARKER") != std::string::npos && msg.find("OLD-RUN-MARKER") == std::string::npos, "log tail mixes runs: %s", msg.c_str());
        ExternalOptions d1 = o, d2 = o;
        ExternalPotential p1(lj_command(), {"Ar"}, d1), p2(lj_command(), {"Ar"}, d2);
        CHECK(p1.log_path() != p2.log_path(), "two potentials share the default log %s", p1.log_path().c_str());
    }
    // Commands with shell operators run without the "exec" prefix (a leading builtin works).
    {
        ExternalPotential p("cd '" + std::string(ATOMFORGE_SOURCE_DIR) + "/quantum' && python3 worker.py lj", {"Ar"}, with("cd"));
        CHECK(evaluate(s, p).forces.size() == s.size(), "cd ... && worker");
    }

}

// ---------------------------------------------------------------- the search pipeline through the bridge
static void test_external_search() {
    std::puts("[external_search] binary_search through ExternalPotential(lj) reproduces the native LJ search; report names the backend");
    if (!have_python()) return;
    // Reduced: two prototypes, no 48-site solid-solution / MC cells (about 50 ms per call in pure Python), short MD,
    // Gamma-only phonons. The native run below covers the full reduced pipeline including MC.
    SearchOptions o;
    o.prototypes = {"L1_2", "B2"};
    o.ss_parents = {}; o.mc_parents = {};
    o.md_steps = 50; o.volume_per_atom = {38, 38}; o.verbose = false;
    o.phonons = true; o.phonon_opt.n[0] = o.phonon_opt.n[1] = o.phonon_opt.n[2] = 1;
    LennardJones nat({"Ar", "Kr"}, kLJ);
    const auto Rn = binary_search("Ar", "Kr", nat, o);
    ExternalOptions xo; xo.log_path = log_path("external_search"); xo.eval_timeout_s = 120;
    ExternalPotential ext(lj_command(), {"Ar", "Kr"}, xo);
    const auto t0 = std::chrono::steady_clock::now();
    const auto Re = binary_search("Ar", "Kr", ext, o);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(Rn.candidates.size() == Re.candidates.size() && !Re.candidates.empty(), "candidate counts %zu vs %zu", Rn.candidates.size(), Re.candidates.size());
    double dE = 0, dV = 0; int nscreen = 0, same_verdict = 0, nph = 0;
    for (std::size_t k = 0; k < std::min(Rn.candidates.size(), Re.candidates.size()); ++k) {
        const auto &a = Rn.candidates[k], &b = Re.candidates[k];
        dE = std::max(dE, std::abs(a.e_per_atom - b.e_per_atom));
        dV = std::max(dV, std::abs(a.vol_per_atom - b.vol_per_atom));
        if (b.screened) { ++nscreen; same_verdict += (a.verdict == b.verdict); nph += b.phonons.computed; }
    }
    std::printf("  %zu candidates, %d screened, %ld worker calls, %.1f s: max |dE| %.1e eV/atom, max |dV| %.1e A^3/atom, potential '%s'\n",
                Re.candidates.size(), nscreen, ext.evaluations(), secs, dE, dV, Re.potential.c_str());
    CHECK(dE < 1e-9 && dV < 1e-9, "bridge and native searches differ (dE %.3e, dV %.3e)", dE, dV);
    CHECK(nscreen > 0 && same_verdict == nscreen && nph == nscreen, "screen: %d screened, %d same verdicts, %d with phonons", nscreen, same_verdict, nph);
    CHECK(Re.potential.rfind("external:lj(", 0) == 0, "SearchResult::potential '%s'", Re.potential.c_str());
    const std::string path = log_path("external_search_report") + ".json";
    write_report_json(path, Re);
    std::ifstream f(path); std::stringstream ss; ss << f.rdbuf();
    const std::string js = ss.str();
    CHECK(js.find("\"potential\": \"external:lj(") != std::string::npos, "report.json does not name the external backend");
    CHECK(js.find("\"phonons\": {\"supercell\": [1, 1, 1]") != std::string::npos && js.find("\"dynamically_stable\": ") != std::string::npos, "report.json has no screen.phonons block");
    bool parsed = true;
    try { (void)json::parse(js); } catch (const std::exception& e) { parsed = false; std::printf("  %s\n", e.what()); }
    CHECK(parsed, "report.json with phonons is not valid JSON");
    std::remove(path.c_str());

    // the full reduced pipeline (solid solutions, MC ordering) on the native LJ, with a non-Zhou element pair
    SearchOptions full; full.verbose = false; full.ss_seeds = 1; full.mc_sweeps = 2; full.md_steps = 50;
    full.mc_parents = {"fcc"}; full.mc_fractions = {0.5}; full.volume_per_atom = {38, 38};
    const auto Rf = binary_search("Ar", "Kr", nat, full);
    int nmc = 0; for (const auto& c : Rf.candidates) nmc += c.source == "mc-ordering";
    CHECK(nmc == 1 && Rf.candidates.size() > Re.candidates.size(), "full native LJ search: %zu candidates, %d MC", Rf.candidates.size(), nmc);
    bool threw = false;
    try { SearchOptions bad = o; bad.volume_per_atom.clear(); (void)binary_search("Ar", "Kr", nat, bad); } catch (const std::runtime_error& e) { threw = std::string(e.what()).find("volume_per_atom") != std::string::npos; }
    CHECK(threw, "an element outside the Zhou04 table without volume_per_atom must be rejected");
}

// ---------------------------------------------------------------- the CLI path
static std::string g_cli;
static int run(const std::string& cmd, std::string& out) {
    out.clear();
    FILE* fp = popen((cmd + " 2>&1").c_str(), "r");
    if (!fp) return 127;
    char buf[512];
    while (std::fgets(buf, sizeof buf, fp)) out += buf;
    const int st = pclose(fp);
    return WIFEXITED(st) ? WEXITSTATUS(st) : 128;
}
static double field(const std::string& out, const char* key) {
    const auto k = out.find(key);
    return k == std::string::npos ? std::nan("") : std::atof(out.c_str() + k + std::strlen(key));
}
static void test_external_cli() {
    std::puts("[external_cli] atomforge relax --potential external:\"python3 quantum/worker.py lj\" matches --potential lj");
    if (g_cli.empty()) { std::puts("  SKIP: no atomforge path given"); g_skip = true; return; }
    std::string a, b, c;
    const std::string cli = "/usr/bin/timeout 60s '" + g_cli + "'";
    const std::string base = cli + " relax A1 Ar --volume 40 --no-elastic";
    // Exit codes without python3: 1 for usage / argument errors, 2 for failures while running.
    int sc = run(base + " --potential bogus", c);
    CHECK(sc == 2 && c.find("unknown --potential") != std::string::npos, "bogus spec: status %d", sc);
    sc = run(cli + " relax A1 Ar --potential lj", c);
    CHECK(sc == 2 && c.find("--volume") != std::string::npos, "element without a reference volume: status %d", sc);
    sc = run(cli + " relax L1_2 Cu Xx", c);
    CHECK(sc == 2 && c.find("no parameters for Xx") != std::string::npos, "unsupported element: status %d '%s'", sc, c.c_str());
    sc = run(cli + " relax L1_2 Cu Au extra", c);
    CHECK(sc == 1 && c.find("unexpected argument extra") != std::string::npos, "extra argument: status %d", sc);
    sc = run(cli + " relax A1 Cu --fmax nope", c);
    CHECK(sc == 1 && c.find("--fmax expects") != std::string::npos, "bad option value: status %d", sc);
    sc = run(base + " --potential external:true --worker-log /nonexistent-atomforge-dir/w.log", c);
    CHECK(sc == 2 && c.find("cannot open worker log") != std::string::npos, "unusable --worker-log: status %d '%s'", sc, c.c_str());
    if (!have_python()) return;

    const int sa = run(base + " --potential lj", a);
    const int sb = run(base + " --potential \"external:" + lj_command() + "\" --worker-log '" + log_path("cli") + "'", b);
    std::printf("%s%s", a.c_str(), b.c_str());
    CHECK(sa == 0 && sb == 0, "exit status %d / %d", sa, sb);
    CHECK(b.find("potential: external:lj(") != std::string::npos && b.find("worker evaluations") != std::string::npos, "external run output");
    const double ea = field(a, "E = "), eb = field(b, "E = "), va = field(a, "V = "), vb = field(b, "V = ");
    CHECK(std::abs(ea - eb) < 1e-12 && std::abs(va - vb) < 1e-12 && ea < -0.07, "E %.6f vs %.6f, V %.4f vs %.4f", ea, eb, va, vb);
    sc = run(base + " --potential external:/nonexistent/worker", c);
    CHECK(sc == 2 && c.find("error:") != std::string::npos && c.find("exit status 127") != std::string::npos, "missing worker: status %d '%s'", sc, c.c_str());
    sc = run(base + " --worker-timeout inf --potential \"external:" + lj_command() + "\"", c);
    CHECK(sc == 0, "--worker-timeout inf: status %d", sc);

    // atomforge killed (Ctrl-C = SIGINT, or SIGKILL) while the worker is busy: the worker must not outlive it. Once with
    // worker.py (which also has its own watchdog), once with a plain shell worker that relies on atomforge's supervisor.
    struct Kill { const char* sig; bool shell; };
    for (const Kill k : {Kill{"INT", false}, Kill{"KILL", false}, Kill{"INT", true}}) {
        const char* sig = k.sig;
        const std::string pf = log_path("orphan") + "." + sig + (k.shell ? ".sh" : "") + ".pid";
        std::remove(pf.c_str());
        const std::string worker = k.shell ? "sh -c 'echo $$ > " + pf + "; exec sleep 30'" : lj_command("--fault hang --pid-file " + pf);
        sc = run(std::string("/usr/bin/timeout -s ") + sig + " 2 '" + g_cli + "' relax A1 Ar --volume 40 --potential \"external:" + worker + "\"", c);
        const long wpid = read_pid(pf);
        const bool gone = wpid > 0 && wait_dead(wpid, 5);
        std::printf("  SIG%s to atomforge mid-%s (%s worker): status %d, worker pid %ld %s\n", sig, k.shell ? "handshake" : "evaluation",
                    k.shell ? "shell" : "worker.py", sc, wpid, gone ? "exited" : "STILL RUNNING");
        CHECK(wpid > 0 && gone, "worker %ld outlived atomforge after SIG%s", wpid, sig);
        if (wpid > 0 && !gone) ::kill(pid_t(wpid), SIGKILL);   // exact pid, created by this test
        std::remove(pf.c_str());
    }

    // search end to end through the worker (reduced: one prototype, no solid solutions or MC, Gamma-only phonons)
    const std::string dir = log_path("cli_search") + ".d";
    sc = run("/usr/bin/timeout 120s '" + g_cli + "' search Ar Kr --potential \"external:" + lj_command() + "\" --volume 38,38"
             " --prototypes L1_2 --no-ss --no-mc --md-steps 20 --phonon-supercell 1 --out '" + dir + "'", c);
    std::ifstream rf(dir + "/report.json");
    std::stringstream rs; rs << rf.rdbuf();
    const std::string rep = rs.str();
    std::printf("  search: status %d, report %zu bytes\n", sc, rep.size());
    CHECK(sc == 0 && c.find("worker evaluations") != std::string::npos && c.find("dynamically stable") != std::string::npos, "search via the worker: status %d\n%s", sc, c.c_str());
    CHECK(rep.find("\"potential\": \"external:lj(") != std::string::npos, "report.json does not name the external backend");
    CHECK(rep.find("\"phonons\": {") != std::string::npos && rep.find("L1_2 Ar3Kr") != std::string::npos, "report.json lacks the screened L1_2 candidate with phonons");
    bool ok = true;
    try { (void)json::parse(rep); } catch (const std::exception&) { ok = false; }
    CHECK(ok, "report.json is not valid JSON");
    std::filesystem::remove_all(dir);
    sc = run("'" + g_cli + "' search Ar Kr --potential lj", c);
    CHECK(sc == 2 && c.find("volume_per_atom") != std::string::npos, "search outside the Zhou04 table without --volume: status %d", sc);
    sc = run("'" + g_cli + "' search Cu Au --volume 12", c);
    CHECK(sc == 1 && c.find("--volume expects VA,VB") != std::string::npos, "bad --volume: status %d", sc);
}

int main(int argc, char** argv) {
    struct Group { const char* name; void (*fn)(); };
    const Group groups[] = {{"lj", test_lj}, {"external", test_external}, {"external_relax", test_external_relax}, {"external_errors", test_external_errors}, {"external_cli", test_external_cli}, {"external_search", test_external_search}};
    const char* want = argc > 1 ? argv[1] : nullptr;
    if (argc > 2) g_cli = argv[2];
    bool ran = false;
    for (const auto& g : groups) {
        if (want && std::strcmp(want, g.name)) continue;
        const auto t0 = std::chrono::steady_clock::now();
        try { g.fn(); } catch (const std::exception& e) { ++g_fail; std::printf("  FAIL %s: uncaught exception: %s\n", g.name, e.what()); }
        std::printf("  (%s: %.2f s)\n", g.name, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ran = true;
    }
    if (!ran) { std::fprintf(stderr, "unknown test group '%s'\n", want); return 2; }
    std::printf("\n%d checks passed, %d failed%s\n", g_pass, g_fail, g_skip ? " (some groups skipped: python3 not found)" : "");
    if (g_fail) return 1;
    return (g_skip && g_pass == 0) ? 77 : 0;
}
