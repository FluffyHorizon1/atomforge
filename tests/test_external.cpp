// atomforge external-potential tests: the native Lennard-Jones potential, the ExternalPotential bridge against it (through
// quantum/worker.py lj), FIRE relaxation through the bridge, and the bridge's failure paths.
//
//   test_external <group>    groups: lj external external_relax external_errors
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
#include <sys/wait.h>

#ifndef ATOMFORGE_SOURCE_DIR
#define ATOMFORGE_SOURCE_DIR "."
#endif

using namespace atomforge;
static int g_fail = 0, g_pass = 0;
static bool g_skip = false;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

static const std::string kWorker = std::string(ATOMFORGE_SOURCE_DIR) + "/quantum/worker.py";
static const LJParams kLJ{0.0104, 3.40, 8.5, 39.948};

static std::string lj_command(const std::string& extra = "") {
    char b[160];
    std::snprintf(b, sizeof b, " lj --epsilon %.17g --sigma %.17g --rc %.17g --mass %.17g", kLJ.epsilon, kLJ.sigma, kLJ.rc, kLJ.mass);
    return "python3 '" + kWorker + "'" + b + (extra.empty() ? "" : " " + extra);
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
    if (!have_python()) return;
    const Structure s = fcc_rattled(8);
    ExternalOptions o; o.hello_timeout_s = 20; o.eval_timeout_s = 20;
    auto with = [&](const char* tag) { ExternalOptions x = o; x.log_path = log_path(tag); std::remove(x.log_path.c_str()); return x; };

    expect_error("worker missing", [&] { ExternalPotential p("/nonexistent/atomforge-worker --x", {"Ar"}, with("missing")); }, "exit status 127");
    expect_error("worker missing (stderr in message)", [&] { ExternalPotential p("/nonexistent/atomforge-worker", {"Ar"}, with("missing2")); }, "worker stderr");
    expect_error("worker exits before handshake", [&] { ExternalPotential p("sh -c 'echo boom >&2; exit 4'", {"Ar"}, with("early")); }, "exit status 4");
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
        ExternalOptions x = with("hang"); x.eval_timeout_s = 1.5;
        ExternalPotential p(lj_command("--fault hang"), {"Ar"}, x);
        expect_error("worker hangs (eval timeout 1.5 s)", [&] { evaluate(s, p); }, "timed out", 6);
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
    // JSON parser corner cases
    CHECK(std::isnan(json::parse("{\"e\": NaN}").get("e")->num), "NaN token");
    CHECK(json::parse("[1, -2.5e-3, \"a\\\"b\", true, null]").arr.size() == 5, "array parse");
    bool threw = false;
    try { json::parse("{\"e\": 1,}"); } catch (const ExternalPotentialError&) { threw = true; }
    CHECK(threw, "trailing comma must be rejected");
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
    if (!have_python()) return;
    std::string a, b, c;
    const std::string base = "/usr/bin/timeout 60s '" + g_cli + "' relax A1 Ar --volume 40 --no-elastic";
    const int sa = run(base + " --potential lj", a);
    const int sb = run(base + " --potential \"external:" + lj_command() + "\" --worker-log '" + log_path("cli") + "'", b);
    std::printf("%s%s", a.c_str(), b.c_str());
    CHECK(sa == 0 && sb == 0, "exit status %d / %d", sa, sb);
    CHECK(b.find("potential: external:lj(") != std::string::npos && b.find("worker evaluations") != std::string::npos, "external run output");
    const double ea = field(a, "E = "), eb = field(b, "E = "), va = field(a, "V = "), vb = field(b, "V = ");
    CHECK(std::abs(ea - eb) < 1e-12 && std::abs(va - vb) < 1e-12 && ea < -0.07, "E %.6f vs %.6f, V %.4f vs %.4f", ea, eb, va, vb);
    int sc = run(base + " --potential external:/nonexistent/worker", c);
    CHECK(sc == 2 && c.find("error:") != std::string::npos && c.find("exit status 127") != std::string::npos, "missing worker: status %d '%s'", sc, c.c_str());
    sc = run(base + " --potential bogus", c);
    CHECK(sc == 2 && c.find("unknown --potential") != std::string::npos, "bogus spec: status %d", sc);
    sc = run("/usr/bin/timeout 60s '" + g_cli + "' relax A1 Ar --potential lj", c);
    CHECK(sc == 2 && c.find("--volume") != std::string::npos, "element without a reference volume: status %d", sc);
}

int main(int argc, char** argv) {
    struct Group { const char* name; void (*fn)(); };
    const Group groups[] = {{"lj", test_lj}, {"external", test_external}, {"external_relax", test_external_relax}, {"external_errors", test_external_errors}, {"external_cli", test_external_cli}};
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
