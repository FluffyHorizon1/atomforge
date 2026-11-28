// atomforge phonon tests: finite-displacement force constants, dynamical matrix and dynamical stability.
//
//   test_phonons phonons      fcc Cu (Zhou04 EAM) stable, acoustic sum rule, X-point frequencies vs experiment;
//                             bcc Cu dynamically unstable (imaginary modes), consistent with the Born test
//   test_phonons phonons_lj   native Lennard-Jones fcc vs the analytic lattice-sum dispersion; the same phonons through
//                             ExternalPotential(worker.py lj) match the native ones (skipped part without python3)
#include "atomforge/atomforge.hpp"
#include "atomforge/external.hpp"
#include "atomforge/lj.hpp"
#include "atomforge/phonons.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifndef ATOMFORGE_SOURCE_DIR
#define ATOMFORGE_SOURCE_DIR "."
#endif

using namespace atomforge;
static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// Experiment: Cu at 80 K, inelastic neutron scattering, Svensson, Brockhouse and Rowe, Phys. Rev. 155, 619 (1967):
// X point, longitudinal 7.25 THz, transverse 5.10 THz.
static constexpr double kCuXL = 7.25, kCuXT = 5.10;

static void test_phonons() {
    std::puts("[phonons] fcc Cu stable with zero acoustic modes at Gamma; bcc Cu dynamically unstable");
    ZhouEAM cu({"Cu"});
    // fcc: relax the primitive cell, phonons on the automatic supercell
    const auto rf = relax(primitive_fcc(3.61), cu, RelaxOptions{CellMode::Hydrostatic, 1e-6, 3000});
    const double a = std::cbrt(4 * rf.s.volume());
    const auto FC = force_constants(rf.s, cu);
    const auto P = phonon_check(FC);
    const double k = 2 * units::pi / a;
    const auto fX = frequencies_THz(FC, {k, 0, 0});
    const auto fL = frequencies_THz(FC, {k / 2, k / 2, k / 2});
    std::printf("  fcc Cu a=%.4f A, supercell %dx%dx%d, %d q: min %.3f max %.3f THz, |acoustic(Gamma)| %.1e THz, %d imaginary\n",
                a, P.n[0], P.n[1], P.n[2], P.nq, P.min_THz, P.max_THz, P.gamma_acoustic_THz, P.n_imag);
    std::printf("  X: T %.3f (exp %.2f)  L %.3f (exp %.2f) THz    L point: %.3f %.3f %.3f THz\n", fX[0], kCuXT, fX[2], kCuXL, fL[0], fL[1], fL[2]);
    CHECK(P.stable && P.n_imag == 0, "fcc Cu must be dynamically stable (%d imaginary modes, min %.3f THz)", P.n_imag, P.min_THz);
    CHECK(P.gamma_acoustic_THz < 1e-3, "acoustic modes at Gamma %.3e THz, expected 0 (sum rule)", P.gamma_acoustic_THz);
    CHECK(std::abs(fX[0] - fX[1]) < 1e-6, "X-point transverse modes not degenerate (%.6f %.6f)", fX[0], fX[1]);
    CHECK(std::abs(fX[0] / kCuXT - 1) < 0.05, "X transverse %.3f THz vs experiment %.2f", fX[0], kCuXT);
    // The Zhou04 EAM is stiffer than experiment for longitudinal modes: L(X) = 8.34 THz, 15 % above 7.25 THz (Svensson et al.
    // 1967), and the highest mode 8.39 THz (+16 %). It is the potential, not the supercell (6x6x6 and 8x8x8 agree, checked
    // below), hence 20 % here against 5 % for the transverse mode.
    CHECK(std::abs(fX[2] / kCuXL - 1) < 0.20, "X longitudinal %.3f THz vs experiment %.2f", fX[2], kCuXL);
    CHECK(std::abs(P.max_THz / kCuXL - 1) < 0.20, "max frequency %.3f THz vs experiment %.2f", P.max_THz, kCuXL);
    CHECK(P.assessed && P.asr_violation < 1e-9, "assessed %d, sum-rule violation before correction %.2e eV/A^2", int(P.assessed), P.asr_violation);
    // Convergence in the supercell. At commensurate points (X, L) 6x6x6 and 8x8x8 agree by construction, so compare an
    // interpolated, non-commensurate q instead.
    PhononOptions o8; o8.n[0] = o8.n[1] = o8.n[2] = 8;
    const auto FC8 = force_constants(rf.s, cu, o8);
    const Vec3 qi{0.31 * k, 0.12 * k, 0.05 * k};
    const auto f6 = frequencies_THz(FC, qi), f8 = frequencies_THz(FC8, qi);
    double dconv = 0; for (int m = 0; m < 3; ++m) dconv = std::max(dconv, std::abs(f6[m] - f8[m]));
    std::printf("  non-commensurate q = (0.31, 0.12, 0.05) 2pi/a: 6x6x6 %.4f %.4f %.4f, 8x8x8 %.4f %.4f %.4f THz (max diff %.1e)\n",
                f6[0], f6[1], f6[2], f8[0], f8[1], f8[2], dconv);
    CHECK(dconv < 0.02, "interpolated frequencies not converged in the supercell: %.3e THz", dconv);
    // A Gamma-only (1x1x1) check cannot see anything (the sum rule zeroes the acoustic modes): it must not claim stability.
    PhononOptions o1; o1.n[0] = o1.n[1] = o1.n[2] = 1;

    // bcc: imaginary transverse modes along Gamma-N, as C' = (C11 - C12)/2 < 0 in the Born test
    const auto rb = relax(primitive_bcc(2.87), cu, RelaxOptions{CellMode::Hydrostatic, 1e-6, 3000});
    const double ab = std::cbrt(2 * rb.s.volume());
    const auto FB = force_constants(rb.s, cu);
    const auto PB = phonon_check(FB);
    const double kb = 2 * units::pi / ab;
    const auto fN = frequencies_THz(FB, {kb / 2, kb / 2, 0});
    std::printf("  bcc Cu a=%.4f A: %d imaginary modes on %d q, min %.3f THz; N point %.3f %.3f %.3f THz\n", ab, PB.n_imag, PB.nq, PB.min_THz, fN[0], fN[1], fN[2]);
    CHECK(!PB.stable && PB.n_imag > 0, "bcc Cu must be dynamically unstable");
    CHECK(fN[0] < -0.5, "lowest N-point mode %.3f THz should be imaginary", fN[0]);
    CHECK(PB.gamma_acoustic_THz < 1e-3, "bcc acoustic modes at Gamma %.3e THz", PB.gamma_acoustic_THz);
    const auto PB1 = phonon_check(rb.s, cu, o1);
    std::printf("  bcc Cu, 1x1x1: %s\n", phonon_status(PB1).c_str());
    CHECK(!PB1.assessed && !PB1.stable && phonon_status(PB1).rfind("not assessed", 0) == 0, "a Gamma-only check must report 'not assessed'");
}

// Multi-atom cells. The conventional 4-atom fcc cell holds the primitive cell's bands folded: at q, its 12 frequencies are
// the primitive ones at q + G for G = 0 and the three conventional reciprocal vectors that are not primitive ones. Then an
// unequal-mass compound (L1_2 Cu3Au) against the same cell doubled along z, which folds q and q + (0, 0, pi/a).
static std::vector<double> sorted_union(std::vector<std::vector<double>> parts) {
    std::vector<double> all;
    for (auto& p : parts) all.insert(all.end(), p.begin(), p.end());
    std::sort(all.begin(), all.end());
    return all;
}
static double max_diff(std::vector<double> a, std::vector<double> b) {
    std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
    if (a.size() != b.size()) return 1e9;
    double d = 0; for (std::size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}
static void test_phonons_basis() {
    std::puts("[phonons_basis] 4-atom cubic fcc vs 1-atom primitive (q + G); L1_2 Cu3Au (unequal masses) vs its doubled cell");
    ZhouEAM cu({"Cu"});
    const double a = 3.6150;
    Structure conv;
    for (const auto& p : prototypes()) if (p.name.rfind("A1", 0) == 0) { conv = p.build(0, 0); break; }
    conv.deform(Mat3::diag(a, a, a));
    PhononOptions oc; oc.n[0] = oc.n[1] = oc.n[2] = 4;            // 16 x 16 x 16 A cube
    PhononOptions op; op.n[0] = op.n[1] = op.n[2] = 8;            // commensurate with every (i, j, k) / 4 of 2 pi / a
    const auto FCc = force_constants(conv, cu, oc);
    const auto FCp = force_constants(primitive_fcc(a), cu, op);
    const double k = 2 * units::pi / a;
    const Vec3 G[4] = {{0, 0, 0}, {k, 0, 0}, {0, k, 0}, {0, 0, k}};
    double worst = 0;
    int nq = 0;
    for (const Vec3 q : {Vec3{0, 0, 0}, Vec3{k / 4, 0, 0}, Vec3{k / 2, k / 4, 0}, Vec3{k / 4, k / 4, k / 4}, Vec3{3 * k / 4, k / 2, k / 4}}) {
        std::vector<std::vector<double>> parts;
        for (const auto& g : G) parts.push_back(frequencies_THz(FCp, q + g));
        worst = std::max(worst, max_diff(frequencies_THz(FCc, q), sorted_union(parts)));
        ++nq;
    }
    std::printf("  fcc Cu: conventional 4-atom cell vs primitive q+G at %d q-points: worst |diff| %.1e THz\n", nq, worst);
    CHECK(worst < 1e-6, "conventional and primitive cells disagree by %.3e THz", worst);

    ZhouEAM cuau({"Cu", "Au"});
    Structure l12;
    for (const auto& p : prototypes()) if (p.name.rfind("L1_2", 0) == 0) { l12 = p.build(0, 1); break; }
    l12.deform(Mat3::diag(3.75, 3.75, 3.75));
    Structure dbl = l12.supercell(1, 1, 2);
    std::swap(dbl.pos[0], dbl.pos[5]); std::swap(dbl.species[0], dbl.species[5]);   // a different basis order on purpose
    PhononOptions o4; o4.n[0] = o4.n[1] = o4.n[2] = 4;
    PhononOptions o2; o2.n[0] = 4; o2.n[1] = 4; o2.n[2] = 2;   // the same 15 x 15 x 15 A supercell lattice
    const auto F1 = force_constants(l12, cuau, o4);
    const auto F2 = force_constants(dbl, cuau, o2);
    const double kz = 2 * units::pi / 3.75;
    double w2 = 0, gam = 0;
    for (const Vec3 q : {Vec3{0, 0, 0}, Vec3{kz / 4, 0, kz / 4}, Vec3{kz / 2, kz / 4, 0}, Vec3{kz / 4, kz / 2, kz / 2}}) {
        w2 = std::max(w2, max_diff(frequencies_THz(F2, q), sorted_union({frequencies_THz(F1, q), frequencies_THz(F1, q + Vec3{0, 0, kz / 2})})));
    }
    const auto f0 = frequencies_THz(F1, {0, 0, 0});
    for (int m = 0; m < 3; ++m) gam = std::max(gam, std::abs(f0[m]));
    std::printf("  L1_2 Cu3Au (m = 63.5, 197.0): 4-atom vs 8-atom cell: worst |diff| %.1e THz; acoustic at Gamma %.1e THz; optic %.3f..%.3f THz\n",
                w2, gam, f0[3], f0.back());
    CHECK(w2 < 1e-6, "Cu3Au: 4-atom and doubled 8-atom cells disagree by %.3e THz", w2);
    CHECK(gam < 1e-4 && f0[3] > 1.0, "Cu3Au Gamma: acoustic %.3e, first optic %.3f THz", gam, f0[3]);
}

// A potential with an "egg-box" error, as real-space DFT grids have: an extra force -kappa (g/2pi) sin(2 pi x / g) on each
// coordinate, so translating the crystal is no longer free. The sum rule then has to correct sum_j Phi = kappa cos(...),
// which the check reports (kappa for an atom at the origin).
struct EggBox final : Potential {
    const Potential& base; double kappa, g;
    EggBox(const Potential& b, double k, double gg) : base(b), kappa(k), g(gg) {}
    EvalResult evaluate(const Structure& s, NeighborList& nl) const override {
        EvalResult r = base.evaluate(s, nl);
        for (std::size_t i = 0; i < s.size(); ++i) for (int a = 0; a < 3; ++a)
            r.forces[i][a] -= kappa * g / (2 * units::pi) * std::sin(2 * units::pi * s.pos[i][a] / g);
        return r;
    }
    double cutoff() const override { return base.cutoff(); }
    const std::vector<std::string>& elements() const override { return base.elements(); }
    double mass(int sp) const override { return base.mass(sp); }
    std::string name() const override { return "eggbox"; }
};
static void test_phonons_asr() {
    std::puts("[phonons_asr] the sum-rule violation before correction is reported (egg-box potential)");
    ZhouEAM cu({"Cu"});
    const auto clean = phonon_check(primitive_fcc(3.615), cu);
    EggBox egg(cu, 0.05, 0.37);
    const auto dirty = phonon_check(primitive_fcc(3.615), egg);
    std::printf("  plain EAM: %.1e eV/A^2; with a 0.05 eV/A^2 egg box: %.6f eV/A^2\n", clean.asr_violation, dirty.asr_violation);
    CHECK(clean.asr_violation < 1e-9, "plain EAM violation %.3e", clean.asr_violation);
    // the central difference with h = 0.01 A sees kappa sin(x)/x, x = 2 pi h / g
    const double x = 2 * units::pi * 0.01 / 0.37, expect = 0.05 * std::sin(x) / x;
    CHECK(std::abs(dirty.asr_violation - expect) < 1e-6, "egg-box violation reported as %.6f, expected %.6f", dirty.asr_violation, expect);
}

// The stability screen's verdicts with phonons, on synthetic candidates (every branch), and one real EAM search.
static void test_phonons_search() {
    std::puts("[phonons_search] screen verdicts with phonons; a reduced Cu-Au search screens L1_2 with assessed phonons");
    Candidate c;
    c.e_hull = 0; c.ef = -0.09; c.elastic.mechanically_stable = true; c.thermal.retained = true;
    c.phonons.computed = true; c.phonons.assessed = true; c.phonons.stable = true;
    CHECK(screen_verdict(c, true) == "STABLE (on hull, mechanically + thermally + dynamically stable)", "%s", screen_verdict(c, true).c_str());
    CHECK(screen_verdict(c, false) == "STABLE (on hull, mechanically + thermally stable)", "%s", screen_verdict(c, false).c_str());
    Candidate na = c; na.phonons.assessed = false; na.phonons.stable = false;
    CHECK(screen_verdict(na, true).find("phonons not assessed") != std::string::npos, "%s", screen_verdict(na, true).c_str());
    Candidate im = c; im.phonons.stable = false; im.phonons.n_imag = 3;
    CHECK(screen_verdict(im, true) == "UNSTABLE (imaginary phonons)", "%s", screen_verdict(im, true).c_str());
    Candidate imna = im; imna.phonons.assessed = false;   // imaginary modes count even on a small grid
    CHECK(screen_verdict(imna, true) == "UNSTABLE (imaginary phonons)", "%s", screen_verdict(imna, true).c_str());
    Candidate meta = c; meta.e_hull = 0.01;
    CHECK(screen_verdict(meta, true) == "METASTABLE (near hull; dynamically stable)", "%s", screen_verdict(meta, true).c_str());
    Candidate born = im; born.elastic.mechanically_stable = false;
    CHECK(screen_verdict(born, true) == "UNSTABLE (fails Born criteria)", "%s", screen_verdict(born, true).c_str());

    SearchOptions o;
    o.prototypes = {"L1_2"}; o.ss_parents = {}; o.mc_parents = {}; o.md_steps = 200; o.verbose = false; o.phonons = true;
    const auto R = binary_search("Cu", "Au", o);
    int n = 0, stable = 0;
    for (const auto& k : R.candidates) {
        if (!k.screened) continue;
        ++n;
        std::printf("  %-14s %s | phonons %dx%dx%d %s, max %.2f THz\n", k.label.c_str(), k.verdict.c_str(), k.phonons.n[0], k.phonons.n[1],
                    k.phonons.n[2], phonon_status(k.phonons).c_str(), k.phonons.max_THz);
        stable += k.verdict == "STABLE (on hull, mechanically + thermally + dynamically stable)" && k.phonons.assessed;
    }
    CHECK(n >= 1 && stable >= 1, "no screened Cu-Au L1_2 candidate reached STABLE with assessed phonons (%d screened)", n);
}

// Analytic dynamical matrix of a monatomic Bravais lattice with a pair potential: D(q) = 1/m sum_R K(R) (1 - cos q.R),
// K_ab = V'' R_a R_b / R^2 + V'/R (delta_ab - R_a R_b / R^2).
static std::array<double, 3> lj_analytic(const LJParams& p, double a, const Vec3& q) {
    std::array<std::array<double, 3>, 3> D{};
    const Structure c = primitive_fcc(a);
    const int R = int(std::ceil(p.rc / (a / std::sqrt(3.0)))) + 1;
    for (int i = -R; i <= R; ++i) for (int j = -R; j <= R; ++j) for (int k = -R; k <= R; ++k) {
        if (!i && !j && !k) continue;
        const Vec3 r = c.cell.row(0) * i + c.cell.row(1) * j + c.cell.row(2) * k;
        const double d = norm(r);
        if (d >= p.rc) continue;
        const double s6 = std::pow(p.sigma / d, 6), s12 = s6 * s6;
        const double v1 = -24 * p.epsilon * (2 * s12 - s6) / d, v2 = 24 * p.epsilon * (26 * s12 - 7 * s6) / (d * d);
        const double cq = 1 - std::cos(dot(q, r));
        for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y)
            D[x][y] += cq * ((v2 - v1 / d) * r[x] * r[y] / (d * d) + (x == y ? v1 / d : 0)) / p.mass;
    }
    std::array<double, 3> f{};
    const auto ev = sym_eigenvalues<3>(D);
    for (int k = 0; k < 3; ++k) f[k] = (ev[k] >= 0 ? 1 : -1) * std::sqrt(std::abs(ev[k])) * units::THz_per_sqrt_eV_A2_amu;
    return f;
}

static void test_phonons_lj() {
    std::puts("[phonons_lj] Lennard-Jones fcc: finite displacements vs the analytic dispersion; native vs external worker");
    LJParams p; p.epsilon = 0.0104; p.sigma = 3.40; p.rc = 8.5; p.mass = 39.948;
    LennardJones lj({"Ar"}, p);
    const auto r = relax(primitive_fcc(5.3), lj, RelaxOptions{CellMode::Hydrostatic, 1e-7, 3000});
    const double a = std::cbrt(4 * r.s.volume());
    const auto FC = force_constants(r.s, lj);
    double worst = 0, fmax = 0;
    int nq = 0;
    for (const auto& q : commensurate_qpoints(FC)) {
        const auto f = frequencies_THz(FC, q);
        const auto g = lj_analytic(p, a, q);
        for (int k = 0; k < 3; ++k) { worst = std::max(worst, std::abs(f[k] - g[k])); fmax = std::max(fmax, g[k]); }
        ++nq;
    }
    const double k = 2 * units::pi / a;
    const auto fX = frequencies_THz(FC, {k, 0, 0});
    const auto gX = lj_analytic(p, a, {k, 0, 0});
    std::printf("  a=%.4f A, %dx%dx%d supercell, %d commensurate q: max %.4f THz, worst |FD - analytic| %.2e THz; X: %.4f %.4f vs %.4f %.4f\n",
                a, FC.n[0], FC.n[1], FC.n[2], nq, fmax, worst, fX[0], fX[2], gX[0], gX[2]);
    CHECK(fmax > 1.0, "LJ dispersion implausibly soft (%.3f THz)", fmax);
    CHECK(worst < 1e-3, "finite-displacement phonons differ from the analytic LJ dispersion by %.3e THz", worst);
    // off the commensurate grid (interpolated): still close because the 6x6x6 supercell covers the 8.5 A cutoff
    const Vec3 qg{0.31 * k, 0.12 * k, 0.05 * k};
    const auto fg = frequencies_THz(FC, qg);
    const auto gg = lj_analytic(p, a, qg);
    double wg = 0; for (int i = 0; i < 3; ++i) wg = std::max(wg, std::abs(fg[i] - gg[i]));
    CHECK(wg < 5e-3, "interpolated q: |FD - analytic| %.3e THz", wg);

    if (std::system("command -v python3 >/dev/null 2>&1") != 0) { std::puts("  SKIP external part: python3 not found"); return; }
    char cmd[512];
    std::snprintf(cmd, sizeof cmd, "python3 '%s/quantum/worker.py' lj --epsilon %.17g --sigma %.17g --rc %.17g --mass %.17g", ATOMFORGE_SOURCE_DIR, p.epsilon, p.sigma, p.rc, p.mass);
    ExternalOptions xo; xo.eval_timeout_s = 120;
    ExternalPotential ext(cmd, {"Ar"}, xo);
    PhononOptions o; o.n[0] = o.n[1] = o.n[2] = 4;
    const auto FN = force_constants(r.s, lj, o);
    const auto FE = force_constants(r.s, ext, o);
    double d = 0;
    for (std::size_t i = 0; i < FN.phi.size(); ++i) d = std::max(d, std::abs(FN.phi[i] - FE.phi[i]));
    std::printf("  through the bridge (4x4x4, %ld worker calls): max |Phi_native - Phi_external| = %.2e eV/A^2\n", FE.evaluations, d);
    CHECK(d < 1e-8, "external force constants differ by %.3e", d);
}

int main(int argc, char** argv) {
    struct Group { const char* name; void (*fn)(); };
    const Group groups[] = {{"phonons", test_phonons}, {"phonons_lj", test_phonons_lj}, {"phonons_basis", test_phonons_basis},
                            {"phonons_asr", test_phonons_asr}, {"phonons_search", test_phonons_search}};
    const char* want = argc > 1 ? argv[1] : nullptr;
    bool ran = false;
    for (const auto& g : groups) {
        if (want && std::strcmp(want, g.name)) continue;
        const auto t0 = std::chrono::steady_clock::now();
        try { g.fn(); } catch (const std::exception& e) { ++g_fail; std::printf("  FAIL %s: uncaught exception: %s\n", g.name, e.what()); }
        std::printf("  (%s: %.2f s)\n", g.name, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ran = true;
    }
    if (!ran) { std::fprintf(stderr, "unknown test group '%s'\n", want); return 2; }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
