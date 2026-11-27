// atomforge CLI
//   atomforge elements
//   atomforge validate
//   atomforge search <A> <B> [--out DIR] [--T 600] [--md-steps 3000] [--mc-sweeps 400] [--quick]
//   atomforge relax <prototype> <A> [B]        e.g. atomforge relax L1_2 Cu Au
#include "atomforge/atomforge.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>

using namespace atomforge;
namespace fs = std::filesystem;

static void usage() {
    std::puts("atomforge — atomistic alloy structure search\n"
              "  atomforge elements\n"
              "  atomforge validate\n"
              "  atomforge relax <prototype> <A> [B]     prototypes: A1 A2 A3 L1_2 L1_0 L1_1 B2 D0_3 C15 B1\n"
              "  atomforge search <A> <B> [--out DIR] [--T K] [--md-steps N] [--mc-sweeps N] [--quick]");
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

static int cmd_relax(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const std::string proto = argv[2], A = argv[3], B = argc > 4 ? argv[4] : argv[3];
    ZhouEAM pot(A == B ? std::vector<std::string>{A} : std::vector<std::string>{A, B});
    const int ib = A == B ? 0 : 1;
    for (const auto& p : prototypes()) {
        if (p.name.rfind(proto, 0) != 0) continue;
        Structure s = p.build(0, ib);
        detail::scale_to_volume(s, detail::vegard_volume(A, B, s.fraction(ib)));
        auto r = relax(s, pot);
        std::printf("%s %s%s: E = %.5f eV/atom, V = %.4f A^3/atom, %d steps %s\n", p.name.c_str(), A.c_str(), A == B ? "" : B.c_str(),
                    r.energy / s.size(), r.s.volume() / s.size(), r.steps, r.converged ? "(converged)" : "(NOT converged)");
        for (int i = 0; i < 3; ++i) std::printf("  a%d = %9.5f %9.5f %9.5f\n", i + 1, r.s.cell(i, 0), r.s.cell(i, 1), r.s.cell(i, 2));
        auto el = elastic_tensor(r.s, pot);
        std::printf("  C11=%.1f C12=%.1f C44=%.1f GPa | B=%.1f G=%.1f GPa | min eig %.1f -> %s\n", el.C[0][0], el.C[0][1], el.C[3][3],
                    el.bulk_voigt, el.shear_voigt, el.eigen[0], el.mechanically_stable ? "Born-stable" : "Born-UNSTABLE");
        return 0;
    }
    std::fprintf(stderr, "unknown prototype %s\n", proto.c_str());
    return 1;
}

static int cmd_search(int argc, char** argv) {
    if (argc < 4) { usage(); return 1; }
    const std::string A = argv[2], B = argv[3];
    SearchOptions opt;
    std::string out = "results/" + A + B;
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
    auto R = binary_search(A, B, opt);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    fs::create_directories(out + "/structures");
    write_report_json(out + "/report.json", R);
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
    }
    std::printf("\nwrote %s/report.json and %zu structures (POSCAR + extxyz) in %.1f s\n", out.c_str(), R.candidates.size(), secs);
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
