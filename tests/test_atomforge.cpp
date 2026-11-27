// atomforge test suite — analytic derivatives against finite differences, symmetry invariances, relaxation, hull, Born
// stability and the validate path. Same custom CHECK style as qpuforge.
//
//   test_atomforge            run every group
//   test_atomforge <group>    run one group: forces virial invariance relax hull born validate_unit cli io unwrapped md_invariance
//   The CLI `atomforge validate` path is the CTest test named `validate` from ci/actions-ctest, not a group here.
#include "atomforge/atomforge.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <filesystem>
#include <unistd.h>
#include <sys/wait.h>
#include <stdexcept>
namespace fs = std::filesystem;

using namespace atomforge;
static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// ---------------------------------------------------------------- helpers
static Structure proto(const std::string& name, int a, int b) {
    for (const auto& p : prototypes()) if (p.name.rfind(name, 0) == 0) return p.build(a, b);
    throw std::runtime_error("no prototype " + name);
}
// Cu/Au structure at a sensible lattice parameter, with seeded random displacements (breaks all symmetry).
static Structure rattled(const std::string& name, int n1, int n2, int n3, double a0, double amp, unsigned seed) {
    Structure s = proto(name, 0, 1).supercell(n1, n2, n3);
    s.deform(Mat3::diag(a0, a0, a0));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(-amp, amp);
    for (auto& r : s.pos) r += Vec3{U(rng), U(rng), U(rng)};
    return s;
}
static double energy(const Structure& s, const Potential& pot) { return evaluate(s, pot).energy; }
static double max_force(const std::vector<Vec3>& f) { double m = 0; for (const auto& v : f) m = std::max(m, norm(v)); return m; }
static Mat3 rotation(const Vec3& axis, double th) {
    const Vec3 u = axis * (1.0 / norm(axis));
    const double c = std::cos(th), s = std::sin(th), t = 1 - c;
    Mat3 R;
    R(0, 0) = t * u.x * u.x + c;       R(0, 1) = t * u.x * u.y - s * u.z; R(0, 2) = t * u.x * u.z + s * u.y;
    R(1, 0) = t * u.x * u.y + s * u.z; R(1, 1) = t * u.y * u.y + c;       R(1, 2) = t * u.y * u.z - s * u.x;
    R(2, 0) = t * u.x * u.z - s * u.y; R(2, 1) = t * u.y * u.z + s * u.x; R(2, 2) = t * u.z * u.z + c;
    return R;
}

// ---------------------------------------------------------------- tests
static void test_forces() {
    std::puts("[forces] analytic forces vs central finite differences of the energy");
    ZhouEAM pot({"Cu", "Au"});
    struct Case { const char* name; Structure s; };
    const Case cases[] = {{"fcc Cu 2x2x1 (16 atoms, rattled)", rattled("A1", 2, 2, 1, 3.61, 0.08, 1)},
                          {"L1_2 CuAu3 2x1x1 (8 atoms, rattled)", rattled("L1_2", 2, 1, 1, 3.75, 0.08, 2)},
                          {"L1_1 CuAu triclinic 2x2x2 (16 atoms, rattled)", rattled("L1_1", 2, 2, 2, 3.75, 0.06, 3)}};
    const double h = 1e-4;
    for (const auto& c : cases) {
        const auto ev = evaluate(c.s, pot);
        Vec3 fsum{}; for (const auto& f : ev.forces) fsum += f;
        double worst = 0, fscale = max_force(ev.forces);
        const int N = int(c.s.size());
        for (int i = 0; i < N; i += (N > 8 ? 3 : 1)) for (int k = 0; k < 3; ++k) {
            Structure p = c.s, m = c.s;
            p.pos[i][k] += h; m.pos[i][k] -= h;
            const double fd = -(energy(p, pot) - energy(m, pot)) / (2 * h);
            worst = std::max(worst, std::abs(fd - ev.forces[i][k]));
        }
        std::printf("  %-48s max|F|=%.3f eV/A  worst |F_analytic - F_FD| = %.2e  |sum F| = %.2e\n", c.name, fscale, worst, norm(fsum));
        CHECK(fscale > 0.05, "%s: structure is not displaced enough to test forces (max|F| %.3e)", c.name, fscale);
        CHECK(worst < 1e-6, "%s: force vs FD error %.3e", c.name, worst);
        CHECK(norm(fsum) < 1e-10, "%s: net force %.3e (Newton's third law)", c.name, norm(fsum));
    }
}

static void test_virial() {
    std::puts("[virial] analytic virial W = dE/d(eps) vs central finite strain");
    ZhouEAM pot({"Cu", "Au"});
    struct Case { const char* name; Structure s; };
    const Case cases[] = {{"fcc Cu 2x1x1 (8 atoms, rattled)", rattled("A1", 2, 1, 1, 3.55, 0.06, 4)},
                          {"L1_1 CuAu triclinic 2x2x2 (16 atoms, rattled)", rattled("L1_1", 2, 2, 2, 3.70, 0.05, 5)}};
    const double h = 1e-5;
    for (const auto& c : cases) {
        const auto ev = evaluate(c.s, pot);
        double worst = 0, wmax = 0, asym = 0;
        for (int a = 0; a < 3; ++a) for (int b = a; b < 3; ++b) {
            auto strained = [&](double e) {
                Mat3 eps; eps(a, b) += e; if (a != b) eps(b, a) += e;
                Structure s = c.s; s.deform(Mat3::identity() + eps); return s;
            };
            // symmetric strain eps_ab = eps_ba = e changes E by (W_ab + W_ba) e = 2 W_ab e (off-diagonal)
            const double fd = (energy(strained(h), pot) - energy(strained(-h), pot)) / (2 * h) / (a == b ? 1.0 : 2.0);
            worst = std::max(worst, std::abs(fd - ev.virial(a, b)));
            wmax = std::max(wmax, std::abs(ev.virial(a, b)));
            asym = std::max(asym, std::abs(ev.virial(a, b) - ev.virial(b, a)));
        }
        std::printf("  %-48s max|W|=%.3f eV  worst |W_analytic - W_FD| = %.2e  asymmetry %.1e\n", c.name, wmax, worst, asym);
        CHECK(wmax > 0.5, "%s: virial is trivially small (%.3e)", c.name, wmax);
        CHECK(worst < 1e-5 * std::max(1.0, wmax), "%s: virial vs finite strain error %.3e", c.name, worst);
        CHECK(asym < 1e-10, "%s: virial not symmetric (%.3e)", c.name, asym);
    }
    // Stress sign convention: compressed fcc Cu is under compression, i.e. negative-pressure tensor (stress = W/V < 0 means
    // the lattice wants to expand) -> trace of W must be negative when compressed, positive when stretched.
    ZhouEAM cu({"Cu"});
    Structure s = proto("A1", 0, 0);
    s.deform(Mat3::diag(3.45, 3.45, 3.45));
    const double w_comp = evaluate(s, cu).virial(0, 0);
    s.deform(Mat3::diag(3.80 / 3.45, 3.80 / 3.45, 3.80 / 3.45));
    const double w_tens = evaluate(s, cu).virial(0, 0);
    std::printf("  fcc Cu W_xx: %.3f eV at a=3.45 (compressed), %.3f eV at a=3.80 (stretched)\n", w_comp, w_tens);
    CHECK(w_comp < 0 && w_tens > 0, "virial sign convention (tensile positive)");
}

static void test_invariance() {
    std::puts("[invariance] energy under translation, rotation and periodic-image shifts");
    ZhouEAM pot({"Cu", "Au"});
    const Structure s = rattled("L1_1", 2, 2, 2, 3.75, 0.07, 6);
    const auto ref = evaluate(s, pot);
    const double tol = 1e-9 * s.size();

    Structure t = s;
    for (auto& r : t.pos) r += Vec3{0.37, -1.21, 2.05};
    const auto et = evaluate(t, pot);
    double df = 0; for (std::size_t i = 0; i < s.size(); ++i) df = std::max(df, norm(et.forces[i] - ref.forces[i]));
    std::printf("  translation: dE = %.2e eV, max dF = %.2e\n", et.energy - ref.energy, df);
    CHECK(std::abs(et.energy - ref.energy) < tol, "translation dE %.3e", et.energy - ref.energy);
    CHECK(df < 1e-9, "translation changes forces by %.3e", df);

    // rigid rotation of cell + atoms: E invariant, F rotates, virial transforms as R W R^T
    const Mat3 R = rotation({1, 2, -0.5}, 0.83);
    Structure q = s; q.deform(R);
    const auto eq = evaluate(q, pot);
    double dfr = 0; for (std::size_t i = 0; i < s.size(); ++i) dfr = std::max(dfr, norm(eq.forces[i] - R * ref.forces[i]));
    const Mat3 Wr = R * ref.virial * transpose(R);
    double dW = 0; for (int k = 0; k < 9; ++k) dW = std::max(dW, std::abs(eq.virial.m[k] - Wr.m[k]));
    std::printf("  rotation:    dE = %.2e eV, max |F - R F0| = %.2e, max |W - R W0 R^T| = %.2e\n", eq.energy - ref.energy, dfr, dW);
    CHECK(std::abs(eq.energy - ref.energy) < tol, "rotation dE %.3e", eq.energy - ref.energy);
    CHECK(dfr < 1e-9, "forces do not rotate covariantly: %.3e", dfr);
    CHECK(dW < 1e-9, "virial does not rotate covariantly: %.3e", dW);

    // move atoms by whole lattice vectors (PBC images of the same crystal), then wrap back into the cell.
    // (Evaluating the shifted, unwrapped structure directly is covered by test_unwrapped.)
    Structure im = s;
    for (std::size_t i = 0; i < im.size(); ++i) {
        const int a = int(i % 5) - 2, b = int((i / 5) % 3) - 1, c = int(i % 2) ? 3 : -2;
        im.pos[i] += s.cell.row(0) * a + s.cell.row(1) * b + s.cell.row(2) * c;
    }
    im.wrap();
    const auto ei = evaluate(im, pot);
    double dfi = 0; for (std::size_t i = 0; i < s.size(); ++i) dfi = std::max(dfi, norm(ei.forces[i] - ref.forces[i]));
    std::printf("  image shift + wrap(): dE = %.2e eV, max dF = %.2e\n", ei.energy - ref.energy, dfi);
    CHECK(std::abs(ei.energy - ref.energy) < tol, "image-shift dE %.3e", ei.energy - ref.energy);
    CHECK(dfi < 1e-9, "image shift changes forces by %.3e", dfi);

    // supercell consistency: E(2x1x1 of a primitive-ish cell) = 2 E(cell), and the same cell expressed with a different lattice basis
    Structure p = proto("A1", 0, 0); p.deform(Mat3::diag(3.60, 3.60, 3.60));
    ZhouEAM cu({"Cu"});
    const double e1 = energy(p, cu) / p.size(), e2 = energy(p.supercell(2, 1, 1), cu) / (2 * p.size()), e3 = energy(p.supercell(2, 2, 2), cu) / (8 * p.size());
    std::printf("  supercells: E/atom = %.8f (1x1x1)  %.8f (2x1x1)  %.8f (2x2x2)\n", e1, e2, e3);
    CHECK(std::abs(e1 - e2) < 1e-9 && std::abs(e1 - e3) < 1e-9, "energy per atom depends on supercell size");
}

// Regression test for a fixed bug: NeighborList::build used to scan only images up to ceil((rc+skin)/width) around the home
// cell, so it silently dropped pairs when atoms were not in (or very near) the home cell: offsets of 2+ cells in a 7 A cell,
// or of just 0/1 cells between atoms in a small triclinic cell. The energy then depended on how far an atom had drifted.
// The image range now grows with the spread of the atoms' fractional coordinates.
static void test_unwrapped() {
    std::puts("[unwrapped] energy and forces of atoms displaced by whole lattice vectors (no wrap)");
    ZhouEAM cu({"Cu"});
    Structure s = proto("A1", 0, 0).supercell(2, 2, 2);
    s.deform(Mat3::diag(3.6148, 3.6148, 3.6148));
    const auto ref = evaluate(s, cu);
    const double e0 = ref.energy;
    for (int k : {-3, -2, 2, 3, 7}) {
        Structure t = s; t.pos[0] += s.cell.row(0) * k;
        const auto ev = evaluate(t, cu);
        const double d = ev.energy - e0;
        double df = 0; for (std::size_t i = 0; i < s.size(); ++i) df = std::max(df, norm(ev.forces[i] - ref.forces[i]));
        std::printf("  atom 0 shifted by %2d lattice vectors: dE = %.3e eV, max dF = %.2e eV/A\n", k, d, df);
        CHECK(std::abs(d) < 1e-9, "energy changed by %.3e eV after shifting one atom by %d lattice vectors", d, k);
        CHECK(df < 1e-9, "forces changed by %.3e eV/A after shifting one atom by %d lattice vectors", df, k);
    }
    // also in a small triclinic cell, where even a single-cell offset between atoms is not enough to be safe
    const Structure q = rattled("L1_1", 2, 2, 2, 3.75, 0.07, 6);
    ZhouEAM pot({"Cu", "Au"});
    Structure im = q;
    for (std::size_t i = 0; i < im.size(); ++i) im.pos[i] += q.cell.row(0) * int(i % 2) - q.cell.row(1) * int((i / 2) % 2) + q.cell.row(2) * int((i / 4) % 2);
    const double d = energy(im, pot) - energy(q, pot);
    std::printf("  triclinic CuAu, atoms offset by 0/1 lattice vectors: dE = %.3e eV\n", d);
    CHECK(std::abs(d) < 1e-9, "energy changed by %.3e eV after offsetting atoms by 0/1 lattice vectors", d);
    // every atom shifted by a different, larger offset, in both a cubic and the triclinic cell
    Structure big = q;
    for (std::size_t i = 0; i < big.size(); ++i) big.pos[i] += q.cell.row(0) * (int(i % 5) - 2) + q.cell.row(1) * (int(i % 3) - 1) + q.cell.row(2) * ((i % 2) ? 4 : -3);
    const double db = energy(big, pot) - energy(q, pot);
    std::printf("  triclinic CuAu, atoms offset by up to 4 lattice vectors: dE = %.3e eV\n", db);
    CHECK(std::abs(db) < 1e-9, "energy changed by %.3e eV after offsetting atoms by up to 4 lattice vectors", db);
    // wrapped configurations must be unaffected: wrap() of the shifted structures gives the same energy
    Structure bw = big; bw.wrap();
    CHECK(std::abs(energy(bw, pot) - energy(q, pot)) < 1e-9, "wrapped copy of the shifted structure has a different energy");
}

// MD energy must not depend on which periodic image each atom started in. Shifts of +-1 cell stay in place; shifts of 3-6 cells
// make run_md wrap atoms back (it wraps atoms that are more than one cell outside the home cell).
static void test_md_invariance() {
    std::puts("[md_invariance] run_md results are invariant under per-atom lattice-vector shifts of the start structure");
    ZhouEAM cu({"Cu"});
    Structure s = proto("A1", 0, 0).supercell(2, 2, 2);
    s.deform(Mat3::diag(3.6148, 3.6148, 3.6148));
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(-0.05, 0.05);
    for (auto& r : s.pos) r += Vec3{U(rng), U(rng), U(rng)};
    MDOptions o; o.T = 600; o.steps = 300; o.seed = 3;
    const auto ref = run_md(s, cu, o);
    struct Case { const char* name; int lo, span; };
    for (const Case& c : {Case{"shifts -1..1 cells (kept unwrapped)", -1, 3}, Case{"shifts 3..6 cells (wrapped by run_md)", 3, 4}}) {
        Structure t = s;
        for (std::size_t i = 0; i < t.size(); ++i)
            t.pos[i] += s.cell.row(0) * (c.lo + int(i % c.span)) + s.cell.row(1) * (c.lo + int((i / 2) % c.span)) - s.cell.row(2) * (c.lo + int((i / 3) % c.span));
        const auto md = run_md(t, cu, o);
        // final positions agree modulo a lattice vector
        const Mat3 inv_t = transpose(inverse(ref.final.cell));
        double dpos = 0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            Vec3 f = inv_t * (md.final.pos[i] - ref.final.pos[i]);
            for (int k = 0; k < 3; ++k) f[k] -= std::round(f[k]);
            dpos = std::max(dpos, norm(transpose(ref.final.cell) * f));
        }
        std::printf("  %-40s <Epot> %.8f vs %.8f eV/atom, <T> %.3f vs %.3f K, <P> %.4f vs %.4f GPa, rms_disp %.6f vs %.6f, max dpos %.1e A\n",
                    c.name, md.mean_epot_per_atom, ref.mean_epot_per_atom, md.mean_T, ref.mean_T, md.mean_P_GPa, ref.mean_P_GPa, md.rms_disp, ref.rms_disp, dpos);
        CHECK(std::abs(md.mean_epot_per_atom - ref.mean_epot_per_atom) < 1e-8, "%s: mean Epot differs by %.3e eV/atom", c.name, md.mean_epot_per_atom - ref.mean_epot_per_atom);
        CHECK(std::abs(md.mean_T - ref.mean_T) < 1e-6, "%s: mean T differs by %.3e K", c.name, md.mean_T - ref.mean_T);
        CHECK(std::abs(md.mean_P_GPa - ref.mean_P_GPa) < 1e-6, "%s: mean P differs by %.3e GPa", c.name, md.mean_P_GPa - ref.mean_P_GPa);
        CHECK(std::abs(md.rms_disp - ref.rms_disp) < 1e-6, "%s: rms displacement differs by %.3e A", c.name, md.rms_disp - ref.rms_disp);
        CHECK(dpos < 1e-6, "%s: final positions differ by %.3e A (mod lattice)", c.name, dpos);
    }
    // sanity: the control must actually move atoms, otherwise the comparison proves nothing
    CHECK(ref.rms_disp > 0.01, "MD barely moved the atoms (rms displacement %.3e A)", ref.rms_disp);
}

static void test_relax() {
    std::puts("[relax] FIRE convergence on small cells");
    ZhouEAM cu({"Cu"});
    // 1. perfect fcc Cu, started 3% off in lattice constant: full cell relaxation must find the equilibrium a0
    {
        Structure s = proto("A1", 0, 0);
        s.deform(Mat3::diag(3.72, 3.72, 3.72));
        const auto t0 = std::chrono::steady_clock::now();
        const auto r = relax(s, cu);
        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const double a = std::cbrt(r.s.volume());
        double smax = 0; for (int k = 0; k < 9; ++k) smax = std::max(smax, std::abs(r.stress_GPa.m[k]));
        std::printf("  fcc Cu from a=3.72: converged=%d in %d steps (%.2f s), a=%.5f A, E=%.5f eV/atom, max|stress|=%.2e GPa\n", r.converged,
                    r.steps, secs, a, r.energy / s.size(), smax);
        CHECK(r.converged, "fcc relax did not converge in %d steps", r.steps);
        CHECK(std::abs(a - 3.6148) < 2e-3, "a0 = %.5f, expected 3.6148 (Zhou04 Cu)", a);
        CHECK(smax < 0.05, "residual stress %.3e GPa", smax);
        CHECK(std::abs(r.energy / s.size() + 3.54) < 0.01, "Ec = %.5f eV/atom, expected -3.54", r.energy / s.size());
        // 2. lower energy than the starting point, and independent re-evaluation confirms the force criterion
        CHECK(r.energy < energy(s, cu) - 1e-3, "relaxation did not lower the energy");
    }
    // 2. rattled L1_2 Cu3Au: atoms-only (fixed cell) then full relax; check the final forces independently
    {
        ZhouEAM pot({"Cu", "Au"});
        Structure s = proto("L1_2", 0, 1);
        s.deform(Mat3::diag(3.78, 3.78, 3.78));
        std::mt19937 rng(11);
        std::uniform_real_distribution<double> U(-0.1, 0.1);
        for (auto& r : s.pos) r += Vec3{U(rng), U(rng), U(rng)};
        const double e0 = energy(s, pot);
        RelaxOptions fixed; fixed.cell = CellMode::Fixed; fixed.fmax = 1e-4;
        const auto rf = relax(s, pot, fixed);
        const double fres = max_force(evaluate(rf.s, pot).forces);
        std::printf("  L1_2 CuAu3 rattled, fixed cell: E %.5f -> %.5f eV, %d steps, residual max|F| = %.2e eV/A\n", e0, rf.energy, rf.steps, fres);
        CHECK(rf.converged, "fixed-cell relax did not converge");
        CHECK(rf.energy < e0 - 1e-3, "no energy decrease");
        CHECK(fres < 1.5e-4, "independent force check: max|F| = %.3e after relax to 1e-4", fres);
        // atoms must be back on the perfect L1_2 sites: nearest-neighbour distance a/sqrt(2) for all pairs of the fcc sublattice
        const double a = std::cbrt(rf.s.volume());
        double dmin = 1e9;
        for (std::size_t i = 0; i < rf.s.size(); ++i) for (std::size_t j = i + 1; j < rf.s.size(); ++j) {
            Vec3 d = rf.s.to_frac(rf.s.pos[j] - rf.s.pos[i]);
            for (int k = 0; k < 3; ++k) d[k] -= std::round(d[k]);
            dmin = std::min(dmin, norm(rf.s.to_cart(d)));
        }
        CHECK(std::abs(dmin - a / std::sqrt(2.0)) < 1e-3, "nearest-neighbour distance %.5f, expected %.5f", dmin, a / std::sqrt(2.0));
        const auto rc = relax(rf.s, pot);
        std::printf("  L1_2 CuAu3 full cell:  a=%.5f A (expt/paper-quoted 3.750), E=%.5f eV/atom, %d steps, converged=%d\n", std::cbrt(rc.s.volume()),
                    rc.energy / rc.s.size(), rc.steps, rc.converged);
        CHECK(rc.converged, "full relax did not converge");
        CHECK(std::abs(std::cbrt(rc.s.volume()) - 3.750) < 0.02, "Cu3Au-type a = %.4f", std::cbrt(rc.s.volume()));
    }
    // 3. external pressure: at +5 GPa the cell must shrink, enthalpy must be reported with the pV term
    {
        Structure s = proto("A1", 0, 0); s.deform(Mat3::diag(3.62, 3.62, 3.62));
        const auto r0 = relax(s, cu);
        RelaxOptions po; po.pressure_GPa = 5.0;
        const auto rp = relax(s, cu, po);
        const double pV = 5.0 / units::eV_per_A3_to_GPa * rp.s.volume();
        std::printf("  fcc Cu at 5 GPa: V %.4f -> %.4f A^3, H-E = %.4f eV (pV = %.4f)\n", r0.s.volume(), rp.s.volume(), rp.enthalpy - rp.energy, pV);
        CHECK(r0.converged && rp.converged, "pressure relax did not converge");
        CHECK(rp.s.volume() < r0.s.volume(), "cell did not shrink under pressure");
        CHECK(std::abs(rp.enthalpy - rp.energy - pV) < 1e-6, "enthalpy != E + pV");
    }
}

static void test_hull() {
    std::puts("[hull] lower convex hull and a reduced Cu-Au search");
    // synthetic points with a known answer
    std::vector<HullPoint> pts = {{0, 0, 0}, {1, 0, 1}, {0.25, -0.10, 2}, {0.5, -0.05, 3} /* above hull */, {0.75, -0.08, 4}, {0.5, -0.12, 5}, {0.25, -0.02, 6}};
    const auto h = lower_hull(pts);
    std::printf("  synthetic hull ids:");
    for (const auto& p : h) std::printf(" %d", p.id);
    std::printf("\n");
    // (0.75, -0.08) lies below the chord (0.5,-0.12)-(1,0) (= -0.06), so it is a vertex; ids 3 and 6 are above the hull
    CHECK(h.size() == 5 && h[0].id == 0 && h[1].id == 2 && h[2].id == 5 && h[3].id == 4 && h[4].id == 1, "synthetic hull vertices wrong (size %zu)", h.size());
    bool convex = true;
    for (std::size_t k = 1; k + 1 < h.size(); ++k)
        convex &= (h[k].x - h[k - 1].x) * (h[k + 1].ef - h[k - 1].ef) - (h[k].ef - h[k - 1].ef) * (h[k + 1].x - h[k - 1].x) > 0;
    CHECK(convex, "hull is not strictly convex");
    CHECK(std::abs(hull_energy_at(h, 0.25) + 0.10) < 1e-12 && std::abs(hull_energy_at(h, 0.0)) < 1e-12 && std::abs(hull_energy_at(h, 1.0)) < 1e-12,
          "hull energy at vertices");
    CHECK(std::abs(hull_energy_at(h, 0.375) + 0.11) < 1e-12, "hull energy interpolation: %.6f", hull_energy_at(h, 0.375));
    for (const auto& p : pts) CHECK(p.ef - hull_energy_at(h, p.x) > -1e-12, "point %d below the hull", p.id);

    // reduced real search (fcc only, one composition, short MC/MD) — pipeline sanity, not a physics benchmark
    SearchOptions o;
    o.verbose = false; o.mc_parents = {"fcc"}; o.mc_fractions = {0.25, 0.5}; o.mc_sweeps = 20; o.ss_seeds = 1; o.md_steps = 200;
    const auto t0 = std::chrono::steady_clock::now();
    const auto R = binary_search("Cu", "Au", o);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double min_ehull = 1e9, min_ef = 1e9, el_min[2] = {1e9, 1e9}, el_lowest = 1e9; int nscreened = 0; bool l12 = false;
    for (const auto& c : R.candidates) {
        min_ehull = std::min(min_ehull, c.e_hull); min_ef = std::min(min_ef, c.ef);
        if (c.source == "element") { el_min[c.x > 0.5] = std::min(el_min[c.x > 0.5], c.ef); el_lowest = std::min(el_lowest, c.ef); }
        nscreened += c.screened;
        if (c.source == "prototype" && c.label.rfind("L1_2 Cu3Au", 0) == 0) { l12 = true; CHECK(c.ef < -0.03, "Cu3Au L1_2 Ef = %.4f eV/atom, expected clearly negative", c.ef); }
    }
    std::printf("  Cu-Au search: %zu candidates, %zu hull vertices, %d screened, min Ehull %.2e, min Ef %.4f eV/atom (%.1f s)\n", R.candidates.size(),
                R.hull.size(), nscreened, min_ehull, min_ef, secs);
    CHECK(l12, "L1_2 Cu3Au not among candidates");
    CHECK(min_ehull > -1e-9, "a candidate lies below the hull (%.3e)", min_ehull);
    CHECK(R.hull.size() >= 3, "hull has no intermediate compound");
    CHECK(R.hull.front().x == 0 && R.hull.back().x == 1 && std::abs(R.hull.front().ef) < 1e-12 && std::abs(R.hull.back().ef) < 1e-12, "hull end points");
    // each element has three trial lattices; the lowest defines the reference (Ef = 0), the others are >= 0
    CHECK(std::abs(el_min[0]) < 1e-9 && std::abs(el_min[1]) < 1e-9 && el_lowest > -1e-9, "elemental references: lowest Ef = %.3e / %.3e", el_min[0], el_min[1]);
    CHECK(nscreened > 0, "no candidate was stability-screened");
}

static void test_born() {
    std::puts("[born] elastic tensor and Born stability");
    ZhouEAM cu({"Cu"});
    Structure s = proto("A1", 0, 0);
    s.deform(Mat3::diag(3.62, 3.62, 3.62));
    const auto r = relax(s, cu, RelaxOptions{CellMode::Full, 1e-5, 6000});
    CHECK(r.converged, "tight fcc relax did not converge");
    const auto el = elastic_tensor(r.s, cu);
    const auto& C = el.C;
    std::printf("  fcc Cu: C11=%.1f C12=%.1f C44=%.1f GPa  B=%.1f G=%.1f  min eig %.1f  stable=%d  (expt 168 / 121 / 75)\n", C[0][0], C[0][1], C[3][3],
                el.bulk_voigt, el.shear_voigt, el.eigen[0], el.mechanically_stable);
    CHECK(el.mechanically_stable, "fcc Cu must be Born-stable");
    CHECK(C[0][0] - C[0][1] > 0 && C[0][0] + 2 * C[0][1] > 0 && C[3][3] > 0, "cubic Born criteria C11>|C12|, C11+2C12>0, C44>0");
    CHECK(C[0][0] > 150 && C[0][0] < 210, "C11 = %.1f GPa", C[0][0]);
    CHECK(C[0][1] > 90 && C[0][1] < 150, "C12 = %.1f GPa", C[0][1]);
    CHECK(C[3][3] > 60 && C[3][3] < 100, "C44 = %.1f GPa", C[3][3]);
    // cubic symmetry of the tensor
    CHECK(std::abs(C[0][0] - C[1][1]) < 0.5 && std::abs(C[1][1] - C[2][2]) < 0.5 && std::abs(C[0][1] - C[1][2]) < 0.5 && std::abs(C[3][3] - C[4][4]) < 0.5 &&
          std::abs(C[4][4] - C[5][5]) < 0.5, "tensor lacks cubic symmetry");
    CHECK(std::abs(C[0][3]) < 0.5 && std::abs(C[0][4]) < 0.5 && std::abs(C[3][4]) < 0.5, "spurious cubic off-block terms");
    const double Bk = (C[0][0] + 2 * C[0][1]) / 3;
    CHECK(std::abs(el.bulk_voigt - Bk) < 0.5, "Voigt bulk modulus %.2f vs (C11+2C12)/3 = %.2f", el.bulk_voigt, Bk);
    // bulk modulus independently from the energy-volume curve E(V): B = V d2E/dV2
    {
        const double V0 = r.s.volume(), dv = 0.01 * V0;
        auto Ev = [&](double V) { Structure t = r.s; const double f = std::cbrt(V / V0); t.deform(Mat3::diag(f, f, f)); return energy(t, cu); };
        const double B_ev = V0 * (Ev(V0 + dv) - 2 * Ev(V0) + Ev(V0 - dv)) / (dv * dv) * units::eV_per_A3_to_GPa;
        std::printf("  bulk modulus: %.2f GPa (elastic tensor) vs %.2f GPa (E-V curve)\n", Bk, B_ev);
        CHECK(std::abs(B_ev - Bk) / Bk < 0.02, "B from E(V) %.2f vs from C_ij %.2f", B_ev, Bk);
    }
    // negative control: bcc Cu (a known mechanically unstable phase for fcc metals) must fail the Born test
    Structure b = proto("A2", 0, 0);
    b.deform(Mat3::diag(2.87, 2.87, 2.87));
    const auto rb = relax(b, cu, RelaxOptions{CellMode::Full, 1e-5, 6000});
    const auto eb = elastic_tensor(rb.s, cu);
    std::printf("  bcc Cu: C11=%.1f C12=%.1f C44=%.1f GPa  C'=(C11-C12)/2=%.1f  min eig %.1f  stable=%d\n", eb.C[0][0], eb.C[0][1], eb.C[3][3],
                (eb.C[0][0] - eb.C[0][1]) / 2, eb.eigen[0], eb.mechanically_stable);
    CHECK(!eb.mechanically_stable && eb.C[0][0] < eb.C[0][1], "bcc Cu should be Born-unstable (C11 < C12)");
}

static void test_validate_unit() {
    std::puts("[validate_unit] element validation logic, pass and fail paths");
    // pass path on a subset (the full 17-element table is exercised by the CLI test, 'atomforge validate')
    std::vector<ElementRef> good;
    for (const char* e : {"Cu", "Fe", "Mg"}) for (const auto& r : reference_elements()) if (!std::strcmp(r.el, e)) good.push_back(r);
    const auto S = validate_elements(good, nullptr);
    std::printf("  subset Cu Fe Mg: gs %d/%d, lattice %d/%d, max |err| %.2f%%, passed=%d\n", S.gs_ok, S.n, S.lattice_ok, S.n, S.max_abs_err_pct, S.passed());
    CHECK(S.n == 3 && S.passed(), "good references must pass");
    CHECK(S.max_abs_err_pct < 0.75, "lattice error %.3f%% larger than the 0.75%% claimed in the README", S.max_abs_err_pct);
    // wrong ground state -> fails
    auto bad_gs = good; bad_gs[0].gs = "bcc";
    const auto G = validate_elements(bad_gs, nullptr);
    CHECK(!G.passed() && G.gs_ok == 2, "wrong ground state must fail (gs_ok=%d)", G.gs_ok);
    // wrong lattice constant (+3%) -> fails, ground state still fine
    auto bad_a = good; bad_a[1].a *= 1.03;
    const auto L = validate_elements(bad_a, nullptr);
    CHECK(!L.passed() && L.gs_ok == 3 && L.lattice_ok == 2, "wrong lattice constant must fail (gs_ok=%d lattice_ok=%d)", L.gs_ok, L.lattice_ok);
    CHECK(L.max_abs_err_pct > 2.5, "reported lattice error %.2f%%", L.max_abs_err_pct);
    // just inside the tolerance still passes
    auto edge = good; edge[0].a /= 1 + 0.5 * validate_lattice_tol_pct / 100;
    CHECK(validate_elements(edge, nullptr).passed(), "error inside tolerance must pass");
    // reference table sanity
    CHECK(reference_elements().size() == 17, "reference table has %zu elements, expected 17", reference_elements().size());
    for (const auto& r : reference_elements()) CHECK(ZhouEAM::lookup(r.el) != nullptr, "%s has no Zhou04 parameters", r.el);
}


// ---------------------------------------------------------------- cli: unknown / invalid search arguments must fail fast (a01)
static std::string atomforge_exe() {
    char buf[4096];
    const ssize_t n = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (n > 0) {
        buf[n] = 0;
        const fs::path p = fs::path(buf).parent_path() / "atomforge";
        if (fs::exists(p)) return p.string();
    }
    if (fs::exists("./atomforge")) return "./atomforge";
    return "atomforge";
}
static int run_cli(const std::string& args, std::string& out, int timeout_s = 3) {
    const std::string cmd = "/usr/bin/timeout " + std::to_string(timeout_s) + "s " + atomforge_exe() + " " + args + " 2>&1";
    out.clear();
    FILE* fp = popen(cmd.c_str(), "r");
    if (!fp) return 127;
    char buf[512];
    while (fgets(buf, sizeof(buf), fp)) out += buf;
    const int st = pclose(fp);
    if (st == -1) return 127;
    if (WIFEXITED(st)) return WEXITSTATUS(st);
    return 128;
}
static void test_cli() {
    std::puts("[cli] unknown, missing, non-numeric and out-of-range search arguments exit non-zero without starting a search");
    std::string out;
    int st = run_cli("", out);
    CHECK(st == 1 && out.find("search") != std::string::npos, "no args: status %d", st);
    st = run_cli("search Cu Au --quik", out);
    CHECK(st == 1 && out.find("unknown option") != std::string::npos, "--quik typo: status %d out='%s'", st, out.c_str());
    CHECK(out.find("[1/") == std::string::npos, "--quik must not start a search");
    st = run_cli("search Cu Au --T abc", out);
    CHECK(st == 1 && out.find("error:") != std::string::npos, "--T abc: status %d out='%s'", st, out.c_str());
    st = run_cli("search Cu Au --md-steps", out);
    CHECK(st == 1 && out.find("requires a value") != std::string::npos, "--md-steps missing: status %d out='%s'", st, out.c_str());
    st = run_cli("search Cu Au --T", out);
    CHECK(st == 1 && out.find("requires a value") != std::string::npos, "--T missing: status %d", st);
    st = run_cli("search Cu Au --out", out);
    CHECK(st == 1 && out.find("requires a value") != std::string::npos, "--out missing: status %d", st);
    st = run_cli("search Cu Au --md-steps -5", out);
    CHECK(st == 1 && out.find("error:") != std::string::npos, "--md-steps -5: status %d out='%s'", st, out.c_str());
    st = run_cli("search Cu Au --mc-sweeps 0", out);
    CHECK(st == 1 && out.find("error:") != std::string::npos, "--mc-sweeps 0: status %d out='%s'", st, out.c_str());
    st = run_cli("search Cu Au --T 0", out);
    CHECK(st == 1, "--T 0 (not > 0): status %d", st);
}

// ---------------------------------------------------------------- io: writers must report stream errors (a02)
static void test_io() {
    std::puts("[io] write_poscar / write_extxyz / write_report_json throw when the path cannot be written");
    const Structure s = proto("A1", 0, 0);
    const std::vector<std::string> els{"Cu"};
    const char* bad = "/nonexistent-dir-atomforge-a02/x";
    auto expect_throw = [&](auto&& fn, const char* what) {
        bool threw = false; std::string msg;
        try { fn(); } catch (const std::runtime_error& e) { threw = true; msg = e.what(); }
        CHECK(threw && msg.find("cannot write") != std::string::npos, "%s: threw=%d msg='%s'", what, threw, msg.c_str());
    };
    expect_throw([&]{ write_poscar(std::string(bad) + ".POSCAR", s, els); }, "write_poscar");
    expect_throw([&]{ write_extxyz(std::string(bad) + ".xyz", s, els); }, "write_extxyz");
    SearchResult R; R.A = "Cu"; R.B = "Au";
    expect_throw([&]{ write_report_json(std::string(bad) + ".json", R); }, "write_report_json");
    // a writable path still succeeds and leaves a file
    const fs::path ok = fs::temp_directory_path() / "atomforge-a02-ok.POSCAR";
    write_poscar(ok.string(), s, els);
    CHECK(fs::file_size(ok) > 0, "successful POSCAR write");
    fs::remove(ok);
}

int main(int argc, char** argv) {
    struct Group { const char* name; void (*fn)(); };
    const Group groups[] = {{"forces", test_forces}, {"virial", test_virial}, {"invariance", test_invariance}, {"relax", test_relax},
                            {"hull", test_hull}, {"born", test_born}, {"validate_unit", test_validate_unit}, {"cli", test_cli}, {"io", test_io}, {"unwrapped", test_unwrapped}, {"md_invariance", test_md_invariance}};
    const char* want = argc > 1 ? argv[1] : nullptr;
    bool ran = false;
    for (const auto& g : groups) {
        if (want && std::strcmp(want, g.name)) continue;
        const auto t0 = std::chrono::steady_clock::now();
        g.fn();
        std::printf("  (%s: %.2f s)\n", g.name, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        ran = true;
    }
    if (!ran) { std::fprintf(stderr, "unknown test group '%s'\n", want); return 2; }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
