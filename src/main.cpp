// atomforge CLI
//   atomforge elements
//   atomforge validate
//   atomforge search <A> <B> [--out DIR] [--T 600] [--md-steps 3000] [--mc-sweeps 400] [--quick] [--materials]
//     --materials also writes <out>/materials-contract-<A><B>.json (forge.materials.v0 sidecar; see include/atomforge/materials.hpp)
//   atomforge relax <prototype> <A> [B] [--potential eam|lj[:eps,sigma,rc]|external:CMD] ...   e.g. atomforge relax L1_2 Cu Au
#include "atomforge/atomforge.hpp"
#include "atomforge/external.hpp"
#include "atomforge/lj.hpp"
#include "atomforge/phonons.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <memory>

using namespace atomforge;
namespace fs = std::filesystem;

static void usage() {
    std::puts("atomforge — atomistic alloy structure search\n"
              "  atomforge elements\n"
              "  atomforge validate\n"
              "  atomforge relax <prototype> <A> [B] [--potential eam|lj[:eps,sigma,rc]|external:CMD] [--volume V]\n"
              "                 [--cell full|hydrostatic|fixed] [--fmax F] [--elastic|--no-elastic] [--phonons] [--phonon-supercell N]\n"
              "                 [--worker-timeout S] [--worker-hello-timeout S] [--worker-log PATH]\n"
              "                                          prototypes: A1 A2 A3 L1_2 L1_0 L1_1 B2 D0_3 C15 B1\n"
              "  atomforge search <A> <B> [--out DIR] [--T K] [--md-steps N] [--mc-sweeps N] [--quick] [--materials]\n"
              "                 [--potential eam|lj[:eps,sigma,rc]|external:CMD] [--volume VA,VB] [--prototypes P1,P2,..]\n"
              "                 [--no-ss] [--no-mc] [--phonons] [--phonon-supercell N] [--worker-timeout S] [--worker-log PATH]");
}

static int cmd_elements() {
    std::printf("Zhou-Johnson-Wadley 2004 EAM elements:\n");
    for (const auto& p : zhou04::table) std::printf("  %-3s Z=%-3d mass=%8.3f  fitted ground state: %s\n", p.name, p.Z, p.mass, p.fitted_ground_state);
    return 0;
}

// Reference data and the relaxation loop live in validate.hpp so the unit tests can inject bad references.
// Exit non-zero if any ground state is wrong (as on ci/actions-ctest) or any lattice constant is more than
// validate_lattice_tol_pct from experiment. A passing run's stdout is unchanged.
static int cmd_validate() {
    const auto S = validate_elements(reference_elements());
    return S.passed() ? 0 : 1;
}

// --potential SPEC:  eam (default) | lj[:eps,sigma,rc] | external:<command>   (see docs/EXTERNAL_PROTOCOL.md)
static std::unique_ptr<Potential> make_potential(const std::string& spec, const std::vector<std::string>& els, const ExternalOptions& xo) {
    if (spec == "eam") return std::make_unique<ZhouEAM>(els);
    if (spec == "lj" || spec.rfind("lj:", 0) == 0) {
        LJParams p;
        if (spec.size() > 3) {
            double v[3];
            char tail = 0;
            if (std::sscanf(spec.c_str() + 3, "%lf,%lf,%lf%c", &v[0], &v[1], &v[2], &tail) != 3)
                throw std::runtime_error("--potential lj:EPS,SIGMA,RC expects three numbers (eV, A, A), got '" + spec + "'");
            p.epsilon = v[0]; p.sigma = v[1]; p.rc = v[2];
        }
        return std::make_unique<LennardJones>(els, p);
    }
    if (spec.rfind("external:", 0) == 0) {
        const std::string cmd = spec.substr(9);
        if (cmd.find_first_not_of(" \t") == std::string::npos) throw std::runtime_error("--potential external: needs a worker command");
        return std::make_unique<ExternalPotential>(cmd, els, xo);
    }
    throw std::runtime_error("unknown --potential '" + spec + "' (expected eam, lj[:eps,sigma,rc] or external:<command>)");
}

static int cmd_relax(int argc, char** argv) {
    std::vector<std::string> pos;
    std::string spec = "eam";
    ExternalOptions xo;
    double volume = 0, fmax = 1e-3;
    bool phonons = false;
    PhononOptions po;
    int elastic = -1;   // -1: default (on for eam/lj, off for external, where it costs at least 24 more evaluations)
    CellMode cell = CellMode::Full;
    // Exit status, as for `search`: 1 for usage and argument errors (unknown option, bad value, unknown prototype, extra
    // arguments), 2 for failures while running (unsupported element, worker failure); see main().
    try {
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto val = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(a + " requires a value");
            return argv[++i];
        };
        auto num = [&](double lo) {
            const std::string v = val();
            char* end = nullptr;
            const double x = std::strtod(v.c_str(), &end);
            if (end == v.c_str() || *end != '\0' || !(x > lo)) throw std::runtime_error(a + " expects a number > " + std::to_string(lo) + ", got '" + v + "'");
            return x;
        };
        if (a == "--potential") spec = val();
        else if (a == "--volume") volume = num(0);
        else if (a == "--fmax") fmax = num(0);
        else if (a == "--worker-timeout") xo.eval_timeout_s = num(0);            // inf (or >= 1e9 s): no deadline
        else if (a == "--worker-hello-timeout") xo.hello_timeout_s = num(0);
        else if (a == "--worker-log") xo.log_path = val();
        else if (a == "--elastic") elastic = 1;
        else if (a == "--no-elastic") elastic = 0;
        else if (a == "--phonons") phonons = true;
        else if (a == "--phonon-supercell") {
            const double n = num(0);
            if (n != std::floor(n) || n > 16) throw std::runtime_error("--phonon-supercell expects an integer 1..16");
            po.n[0] = po.n[1] = po.n[2] = int(n);
            phonons = true;
        }
        else if (a == "--cell") {
            const std::string c = val();
            if (c == "full") cell = CellMode::Full;
            else if (c == "hydrostatic") cell = CellMode::Hydrostatic;
            else if (c == "fixed") cell = CellMode::Fixed;
            else throw std::runtime_error("--cell expects full, hydrostatic or fixed, got '" + c + "'");
        } else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "error: unknown option %s\n", a.c_str());
            usage();
            return 1;
        } else pos.push_back(a);
    }
    } catch (const std::runtime_error& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    if (pos.size() > 3) { std::fprintf(stderr, "error: unexpected argument %s\n", pos[3].c_str()); usage(); return 1; }
    if (pos.size() < 2) { usage(); return 1; }
    const std::string proto = pos[0], A = pos[1], B = pos.size() > 2 ? pos[2] : pos[1];
    const std::vector<std::string> els = A == B ? std::vector<std::string>{A} : std::vector<std::string>{A, B};
    const int ib = A == B ? 0 : 1;
    const auto protos = prototypes();
    const Prototype* P = nullptr;
    for (const auto& p : protos) if (p.name.rfind(proto, 0) == 0) { P = &p; break; }
    if (!P) { std::fprintf(stderr, "error: unknown prototype %s\n", proto.c_str()); return 1; }
    // The potential first: for eam an unsupported element is reported by ZhouEAM itself ("no parameters for Xx").
    const auto pot = make_potential(spec, els, xo);
    if (volume == 0 && (!ZhouEAM::lookup(A) || !ZhouEAM::lookup(B)))
        throw std::runtime_error("no Zhou04 reference volume for " + (ZhouEAM::lookup(A) ? B : A) + ": pass --volume (A^3/atom)");
    Structure s = P->build(0, ib);
    detail::scale_to_volume(s, volume > 0 ? volume : detail::vegard_volume(A, B, s.fraction(ib)));
    const auto t0 = std::chrono::steady_clock::now();
    auto r = relax(s, *pot, RelaxOptions{cell, fmax, 3000, 0.0, 0.15});
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (spec != "eam") std::printf("potential: %s\n", pot->name().c_str());
    std::printf("%s %s%s: E = %.5f eV/atom, V = %.4f A^3/atom, %d steps %s\n", P->name.c_str(), A.c_str(), A == B ? "" : B.c_str(),
                r.energy / s.size(), r.s.volume() / s.size(), r.steps, r.converged ? "(converged)" : "(NOT converged)");
    for (int i = 0; i < 3; ++i) std::printf("  a%d = %9.5f %9.5f %9.5f\n", i + 1, r.s.cell(i, 0), r.s.cell(i, 1), r.s.cell(i, 2));
    if (const auto* ext = dynamic_cast<const ExternalPotential*>(pot.get()))
        std::printf("  %ld worker evaluations in %.1f s (worker log: %s)\n", ext->evaluations(), secs, ext->log_path().c_str());
    if (elastic == 1 || (elastic == -1 && !dynamic_cast<const ExternalPotential*>(pot.get()))) {
        auto el = elastic_tensor(r.s, *pot);
        std::printf("  C11=%.1f C12=%.1f C44=%.1f GPa | B=%.1f G=%.1f GPa | min eig %.1f -> %s\n", el.C[0][0], el.C[0][1], el.C[3][3],
                    el.bulk_voigt, el.shear_voigt, el.eigen[0], el.mechanically_stable ? "Born-stable" : "Born-UNSTABLE");
    }
    if (phonons) {
        const auto t1 = std::chrono::steady_clock::now();
        const auto FC = force_constants(r.s, *pot, po);
        const auto P = phonon_check(FC, po.imag_tol_THz);
        std::printf("  phonons (%dx%dx%d supercell, %ld evaluations, %.1f s): %d commensurate q, min %.3f max %.3f THz, "
                    "%d modes below -%.2f THz -> %s\n", P.n[0], P.n[1], P.n[2], P.evaluations,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count(), P.nq, P.min_THz, P.max_THz,
                    P.n_imag, P.imag_tol_THz, P.stable ? "dynamically stable" : "dynamically UNSTABLE");
    }
    return 0;
}

static int cmd_search(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const std::string A = argv[2], B = argv[3];
    SearchOptions opt;
    std::string out = "results/" + A + B;
    bool materials = false;
    std::string spec = "eam";
    ExternalOptions xo;
    int i = 4;
    auto need = [&](const char* optname) -> const char* {
        if (i + 1 >= argc) {
            std::fprintf(stderr, "error: %s requires a value\n", optname);
            return nullptr;
        }
        return argv[++i];
    };
    auto parse_int = [&](const char* optname, int lo, int* dest) -> bool {
        const char* s = need(optname); if (!s) return false;
        char* end = nullptr;
        const long v = std::strtol(s, &end, 10);
        if (end == s || *end != '\0' || v < lo || v > 1000000000L) {
            std::fprintf(stderr, "error: %s expects an integer >= %d, got '%s'\n", optname, lo, s);
            return false;
        }
        *dest = int(v);
        return true;
    };
    auto parse_double = [&](const char* optname, double lo, double* dest) -> bool {
        const char* s = need(optname); if (!s) return false;
        char* end = nullptr;
        const double v = std::strtod(s, &end);
        if (end == s || *end != '\0' || !(v > lo)) {
            std::fprintf(stderr, "error: %s expects a number > %g, got '%s'\n", optname, lo, s);
            return false;
        }
        *dest = v;
        return true;
    };
    for (; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--out")) {
            const char* s = need("--out"); if (!s) return 1;
            out = s;
        } else if (!std::strcmp(argv[i], "--T")) {
            if (!parse_double("--T", 0.0, &opt.md_T)) return 1;
        } else if (!std::strcmp(argv[i], "--md-steps")) {
            if (!parse_int("--md-steps", 1, &opt.md_steps)) return 1;
        } else if (!std::strcmp(argv[i], "--mc-sweeps")) {
            if (!parse_int("--mc-sweeps", 1, &opt.mc_sweeps)) return 1;
        } else if (!std::strcmp(argv[i], "--materials")) {
            materials = true;
        } else if (!std::strcmp(argv[i], "--potential")) {
            const char* v = need("--potential"); if (!v) return 1;
            spec = v;
        } else if (!std::strcmp(argv[i], "--volume")) {   // VA,VB in A^3/atom (needed outside the Zhou04 table)
            const char* v = need("--volume"); if (!v) return 1;
            double va = 0, vb = 0; char tail = 0;
            if (std::sscanf(v, "%lf,%lf%c", &va, &vb, &tail) != 2 || !(va > 0) || !(vb > 0)) {
                std::fprintf(stderr, "error: --volume expects VA,VB (two positive volumes per atom), got '%s'\n", v);
                return 1;
            }
            opt.volume_per_atom = {va, vb};
        } else if (!std::strcmp(argv[i], "--prototypes")) {   // comma-separated name prefixes, e.g. L1_2,B2
            const char* v = need("--prototypes"); if (!v) return 1;
            opt.prototypes.clear();
            std::string cur;
            for (const char* c = v;; ++c) {
                if (*c == ',' || *c == '\0') { if (!cur.empty()) opt.prototypes.push_back(cur); cur.clear(); if (!*c) break; }
                else cur += *c;
            }
            if (opt.prototypes.empty()) { std::fprintf(stderr, "error: --prototypes needs at least one name\n"); return 1; }
        } else if (!std::strcmp(argv[i], "--no-mc")) {
            opt.mc_parents.clear();
        } else if (!std::strcmp(argv[i], "--no-ss")) {
            opt.ss_parents.clear();
        } else if (!std::strcmp(argv[i], "--phonons")) {
            opt.phonons = true;
        } else if (!std::strcmp(argv[i], "--phonon-supercell")) {
            int n = 0;
            if (!parse_int("--phonon-supercell", 1, &n)) return 1;
            opt.phonons = true;
            opt.phonon_opt.n[0] = opt.phonon_opt.n[1] = opt.phonon_opt.n[2] = n;
        } else if (!std::strcmp(argv[i], "--worker-timeout")) {
            if (!parse_double("--worker-timeout", 0.0, &xo.eval_timeout_s)) return 1;
        } else if (!std::strcmp(argv[i], "--worker-log")) {
            const char* v = need("--worker-log"); if (!v) return 1;
            xo.log_path = v;
        } else if (!std::strcmp(argv[i], "--quick")) {
            opt.mc_sweeps = 100; opt.md_steps = 1000; opt.ss_seeds = 2;
        } else if (argv[i][0] == '-') {
            std::fprintf(stderr, "error: unknown option %s\n", argv[i]);
            usage();
            return 1;
        } else {
            std::fprintf(stderr, "error: unexpected argument %s\n", argv[i]);
            usage();
            return 1;
        }
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto pot = make_potential(spec, {A, B}, xo);
    auto R = binary_search(A, B, *pot, opt);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    fs::create_directories(out + "/structures");
    write_report_json(out + "/report.json", R);
    if (materials) write_materials_contract_json(out + "/materials-contract-" + A + B + ".json", R);
    std::vector<std::string> els{A, B};
    for (const auto& c : R.candidates) {
        const std::string base = out + "/structures/" + std::to_string(c.id);
        write_poscar(base + ".POSCAR", c.relaxed, els);
        write_extxyz(base + ".xyz", c.relaxed, els, c.e_per_atom * c.relaxed.size());
    }

    auto sorted = R.candidates;
    std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b) { return a.x < b.x || (a.x == b.x && a.ef < b.ef); });
    std::printf("\n%s-%s  (%s)   refs: %s (%.4f eV), %s (%.4f eV)\n", A.c_str(), B.c_str(), R.potential.c_str(), R.refA.c_str(), R.muA, R.refB.c_str(), R.muB);
    std::printf("%4s  %-34s %-15s %6s %10s %10s\n", "id", "structure", "source", "x_B", "Ef meV/at", "Ehull");
    for (const auto& c : sorted) {
        if (c.source == "element") continue;
        std::printf("%4d  %-34s %-15s %6.3f %10.1f %10.1f %s\n", c.id, c.label.c_str(), c.source.c_str(), c.x, 1000 * c.ef, 1000 * c.e_hull,
                    c.e_hull < 1e-4 ? "<- hull" : "");
    }
    std::printf("\nStability screen:\n");
    for (const auto& c : sorted) {
        if (!c.screened) continue;
        std::printf("  %-34s Ef=%7.1f meV  B=%6.1f GPa  minCeig=%6.1f GPa  MD %4.0fK: dEq=%6.1f meV shift=%.2f A  => %s\n", c.label.c_str(),
                    1000 * c.ef, c.elastic.bulk_voigt, c.elastic.eigen[0], c.thermal.T, 1000 * c.thermal.dE_quench_per_atom,
                    c.thermal.max_site_shift, c.verdict.c_str());
        if (c.phonons.computed)
            std::printf("  %-34s phonons %dx%dx%d: %d q, min %.3f max %.3f THz, %d imaginary -> %s\n", "", c.phonons.n[0], c.phonons.n[1],
                        c.phonons.n[2], c.phonons.nq, c.phonons.min_THz, c.phonons.max_THz, c.phonons.n_imag,
                        c.phonons.stable ? "dynamically stable" : "dynamically UNSTABLE");
    }
    if (const auto* ext = dynamic_cast<const ExternalPotential*>(pot.get()))
        std::printf("\n%ld worker evaluations (worker log: %s)\n", ext->evaluations(), ext->log_path().c_str());
    std::printf("\nwrote %s/report.json and %zu structures (POSCAR + extxyz) in %.1f s\n", out.c_str(), R.candidates.size(), secs);
    if (materials) std::printf("wrote %s/materials-contract-%s%s.json (forge.materials.v0; Q_dielectric, loss_tangent, junction_gap_ueV are placeholders)\n", out.c_str(), A.c_str(), B.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    const std::string cmd = argv[1];
    try {
        if (cmd == "elements") return cmd_elements();
        if (cmd == "validate") return cmd_validate();
        if (cmd == "relax") return cmd_relax(argc, argv);
        if (cmd == "search") return cmd_search(argc, argv);
    } catch (const std::exception& e) { std::fprintf(stderr, "error: %s\n", e.what()); return 2; }
    usage();
    return 1;
}
