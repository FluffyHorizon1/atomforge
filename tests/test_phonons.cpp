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
    // A larger supercell changes nothing at commensurate points (force constants are converged).
    PhononOptions o8; o8.n[0] = o8.n[1] = o8.n[2] = 8;
    const auto P8 = phonon_check(rf.s, cu, o8);
    CHECK(std::abs(P8.max_THz - P.max_THz) < 1e-3, "max frequency changes with supercell: %.4f vs %.4f", P8.max_THz, P.max_THz);

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
    const Group groups[] = {{"phonons", test_phonons}, {"phonons_lj", test_phonons_lj}};
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
