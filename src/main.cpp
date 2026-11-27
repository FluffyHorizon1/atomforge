// atomforge CLI
//   atomforge elements
//   atomforge validate
//   atomforge search <A> <B> [--out DIR] [--T 600] [--md-steps 3000] [--mc-sweeps 400] [--quick] [--materials]
//     --materials also writes <out>/materials-contract-<A><B>.json (forge.materials.v0 sidecar; see include/atomforge/materials.hpp)
//   atomforge relax <prototype> <A> [B]        e.g. atomforge relax L1_2 Cu Au
#include "atomforge/atomforge.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>

using namespace atomforge;
namespace fs = std::filesystem;

static void usage() {
    std::puts("atomforge — atomistic alloy structure search\n"
              "  atomforge elements\n"
              "  atomforge validate\n"
              "  atomforge relax <prototype> <A> [B]     prototypes: A1 A2 A3 L1_2 L1_0 L1_1 B2 D0_3 C15 B1\n"
              "  atomforge search <A> <B> [--out DIR] [--T K] [--md-steps N] [--mc-sweeps N] [--quick] [--materials]");
}

static int cmd_elements() {
    std::printf("Zhou-Johnson-Wadley 2004 EAM elements:\n");
    for (const auto& p : zhou04::table) std::printf("  %-3s Z=%-3d mass=%8.3f  fitted ground state: %s\n", p.name, p.Z, p.mass, p.fitted_ground_state);
    return 0;
}

// Reference data: experimental lattice constants (A) and cohesive energies (eV/atom).
static int cmd_validate() {
    struct Ref { const char* el; const char* gs; double a; double ec; };
    const Ref refs[] = {{"Cu", "fcc", 3.615, 3.49}, {"Ag", "fcc", 4.086, 2.95}, {"Au", "fcc", 4.078, 3.81}, {"Ni", "fcc", 3.524, 4.44},
                        {"Pd", "fcc", 3.891, 3.89}, {"Pt", "fcc", 3.924, 5.84}, {"Al", "fcc", 4.050, 3.39}, {"Pb", "fcc", 4.950, 2.03},
                        {"Fe", "bcc", 2.867, 4.28}, {"Mo", "bcc", 3.147, 6.82}, {"Ta", "bcc", 3.301, 8.10}, {"W", "bcc", 3.165, 8.90},
                        {"Cr", "bcc", 2.885, 4.10}, {"Mg", "hcp", 3.209, 1.51}, {"Co", "hcp", 2.507, 4.39}, {"Ti", "hcp", 2.951, 4.85},
                        {"Zr", "hcp", 3.232, 6.25}};
    std::printf("%-3s %-4s | %-9s %-9s %-6s | %-9s %-9s | %s\n", "El", "GS", "a_calc", "a_expt", "err%", "Ec_calc", "Ec_expt", "predicted GS");
    int ok = 0;
    for (const auto& r : refs) {
        ZhouEAM pot({r.el});
        double best = 1e9; std::string bestp; double a_gs = 0;
        for (const auto& p : prototypes()) {
            if (p.nB) continue;
            Structure s = p.build(0, 0);
            detail::scale_to_volume(s, 0.74 * std::pow(ZhouEAM::lookup(r.el)->re, 3));
            auto rr = relax(s, pot);
            const double e = rr.energy / s.size();
            if (p.parent == r.gs) a_gs = (p.parent == "hcp") ? norm(rr.s.cell.row(0)) : std::cbrt(rr.s.volume());
            if (e < best - 1e-6) { best = e; bestp = p.parent; }
        }
        const bool gs_ok = bestp == r.gs;
        ok += gs_ok;
        std::printf("%-3s %-4s | %-9.4f %-9.4f %+6.2f | %-9.4f %-9.4f | %s %s\n", r.el, r.gs, a_gs, r.a, 100 * (a_gs - r.a) / r.a,
                    -best, r.ec, bestp.c_str(), gs_ok ? "OK" : "MISMATCH");
    }
    const int nref = int(sizeof(refs) / sizeof(refs[0]));
    std::printf("ground state correct for %d/%d elements\n", ok, nref);
    return ok == nref ? 0 : 1;
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
    bool materials = false;
    for (int i = 4; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
        else if (!std::strcmp(argv[i], "--T") && i + 1 < argc) opt.md_T = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--md-steps") && i + 1 < argc) opt.md_steps = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--mc-sweeps") && i + 1 < argc) opt.mc_sweeps = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--materials")) materials = true;
        else if (!std::strcmp(argv[i], "--quick")) { opt.mc_sweeps = 100; opt.md_steps = 1000; opt.ss_seeds = 2; }
    }
    const auto t0 = std::chrono::steady_clock::now();
    auto R = binary_search(A, B, opt);
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
    }
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
