// atomforge — stability analysis: elastic tensor (Born), MD heat-and-quench, binary convex hull
#pragma once
#include "relax.hpp"
#include <random>
#include <vector>

namespace atomforge {

// ---------------------------------------------------------------- elastic constants
struct ElasticResult {
    std::array<std::array<double, 6>, 6> C{};   // GPa, Voigt
    std::array<double, 6> eigen{};              // GPa
    double bulk_voigt{0}, shear_voigt{0};       // GPa
    bool mechanically_stable{false};            // Born: C positive definite
};

// Relaxed-ion elastic tensor by central finite differences of stress (structure must be relaxed).
inline ElasticResult elastic_tensor(const Structure& relaxed, const Potential& pot, double delta = 0.005) {
    ElasticResult R;
    static const int va[6] = {0, 1, 2, 1, 0, 0}, vb[6] = {0, 1, 2, 2, 2, 1};
    RelaxOptions ro; ro.cell = CellMode::Fixed; ro.fmax = 2e-4;
    auto stress_voigt = [&](const Mat3& eps) {
        Structure s = relaxed;
        s.deform(Mat3::identity() + eps);
        auto r = relax(s, pot, ro);
        auto ev = evaluate(r.s, pot);
        std::array<double, 6> sv{};
        for (int i = 0; i < 6; ++i) sv[i] = ev.virial(va[i], vb[i]) / r.s.volume() * units::eV_per_A3_to_GPa;
        return sv;
    };
    for (int j = 0; j < 6; ++j) {
        Mat3 e;
        if (j < 3) e(j, j) = delta; else { e(va[j], vb[j]) = delta / 2; e(vb[j], va[j]) = delta / 2; }
        auto sp = stress_voigt(e), sm = stress_voigt(e * -1.0);
        for (int i = 0; i < 6; ++i) R.C[i][j] = (sp[i] - sm[i]) / (2 * delta);
    }
    for (int i = 0; i < 6; ++i) for (int j = i + 1; j < 6; ++j) { const double m = 0.5 * (R.C[i][j] + R.C[j][i]); R.C[i][j] = R.C[j][i] = m; }
    R.eigen = sym_eigenvalues<6>(R.C);
    R.mechanically_stable = R.eigen[0] > 0.5;   // GPa tolerance against FD noise
    const auto& C = R.C;
    R.bulk_voigt = (C[0][0] + C[1][1] + C[2][2] + 2 * (C[0][1] + C[0][2] + C[1][2])) / 9.0;
    R.shear_voigt = ((C[0][0] + C[1][1] + C[2][2]) - (C[0][1] + C[0][2] + C[1][2]) + 3 * (C[3][3] + C[4][4] + C[5][5])) / 15.0;
    return R;
}

// ---------------------------------------------------------------- molecular dynamics
struct MDOptions {
    double T = 300;             // K
    int steps = 4000;
    double dt_fs = 2.0;
    double gamma_per_ps = 1.0;  // Langevin friction
    bool barostat = true;       // Berendsen isotropic, P=0
    double tau_p_fs = 500;
    double compressibility_per_GPa = 0.01;
    unsigned seed = 7;
};
struct MDResult {
    Structure final;
    double mean_T{0}, mean_P_GPa{0}, mean_epot_per_atom{0};
    double rms_disp{0};         // A, time-averaged RMS displacement from (drifting) reference
};

// BAOAB Langevin integrator + optional Berendsen barostat.
inline MDResult run_md(Structure s, const Potential& pot, const MDOptions& o) {
    const int N = int(s.size());
    const double dt = o.dt_fs / units::time_fs;
    const double gamma = o.gamma_per_ps / 1000.0 * units::time_fs;
    const double kT = units::kB * o.T;
    std::mt19937 rng(o.seed);
    std::normal_distribution<double> nd(0, 1);
    std::vector<double> m(N);
    std::vector<Vec3> v(N);
    for (int i = 0; i < N; ++i) {
        m[i] = pot.mass(s.species[i]);
        const double sd = std::sqrt(kT / m[i]);
        v[i] = {sd * nd(rng), sd * nd(rng), sd * nd(rng)};
    }
    auto remove_com = [&]() {
        Vec3 P{}; double M = 0;
        for (int i = 0; i < N; ++i) { P += v[i] * m[i]; M += m[i]; }
        const Vec3 vc = P * (1.0 / M);
        for (int i = 0; i < N; ++i) v[i] -= vc;
    };
    remove_com();
    NeighborList nl;
    auto ev = pot.evaluate(s, nl);
    const double c1 = std::exp(-gamma * dt), c2 = std::sqrt(1 - c1 * c1);
    std::vector<Vec3> ref = s.pos;          // unwrapped tracking
    MDResult R;
    double sumT = 0, sumP = 0, sumE = 0, sumD = 0; int nsamp = 0;
    for (int step = 0; step < o.steps; ++step) {
        for (int i = 0; i < N; ++i) v[i] += ev.forces[i] * (0.5 * dt / m[i]);
        for (int i = 0; i < N; ++i) s.pos[i] += v[i] * (0.5 * dt);
        for (int i = 0; i < N; ++i) {
            const double sd = std::sqrt(kT / m[i]) * c2;
            v[i] = v[i] * c1 + Vec3{sd * nd(rng), sd * nd(rng), sd * nd(rng)};
        }
        remove_com();
        for (int i = 0; i < N; ++i) s.pos[i] += v[i] * (0.5 * dt);
        // Keep coordinates bounded so the neighbour-list image range (and cost) does not grow as atoms diffuse. Only atoms that
        // have left the home cell by more than one full cell are wrapped, so trajectories that stay near the cell are untouched
        // bit for bit. `ref` moves with the atom, so the displacement s.pos - ref is unchanged.
        {
            const Mat3 inv_t = transpose(inverse(s.cell)), cell_t = transpose(s.cell);
            for (int i = 0; i < N; ++i) {
                Vec3 f = inv_t * s.pos[i];
                if (f.x >= -1 && f.x < 2 && f.y >= -1 && f.y < 2 && f.z >= -1 && f.z < 2) continue;
                for (int k = 0; k < 3; ++k) f[k] -= std::floor(f[k]);
                const Vec3 wrapped = cell_t * f;
                ref[i] += wrapped - s.pos[i];
                s.pos[i] = wrapped;
            }
        }
        ev = pot.evaluate(s, nl);
        for (int i = 0; i < N; ++i) v[i] += ev.forces[i] * (0.5 * dt / m[i]);

        double K = 0; for (int i = 0; i < N; ++i) K += 0.5 * m[i] * dot(v[i], v[i]);
        const double V = s.volume();
        const double P = (2 * K - (ev.virial(0, 0) + ev.virial(1, 1) + ev.virial(2, 2))) / (3 * V) * units::eV_per_A3_to_GPa;
        if (o.barostat) {
            double mu = 1.0 - o.compressibility_per_GPa * (o.dt_fs / o.tau_p_fs) * (0.0 - P) / 3.0;
            mu = std::clamp(mu, 0.995, 1.005);
            s.deform(Mat3::diag(mu, mu, mu));
            for (auto& r : ref) r = r * mu;
        }
        if (step >= o.steps / 2) {
            sumT += 2 * K / (3 * N * units::kB); sumP += P; sumE += ev.energy / N;
            Vec3 com{}; for (int i = 0; i < N; ++i) com += s.pos[i] - ref[i];
            com = com * (1.0 / N);
            double d2 = 0; for (int i = 0; i < N; ++i) { const Vec3 d = s.pos[i] - ref[i] - com; d2 += dot(d, d); }
            sumD += std::sqrt(d2 / N); ++nsamp;
        }
    }
    R.mean_T = sumT / nsamp; R.mean_P_GPa = sumP / nsamp; R.mean_epot_per_atom = sumE / nsamp; R.rms_disp = sumD / nsamp;
    s.wrap();
    R.final = s;
    return R;
}

struct ThermalCheck {
    double T{0};
    double rms_disp{0};
    double dE_quench_per_atom{0};   // E(quenched after MD) - E(start); << 0 => transformed to lower-energy state
    double max_site_shift{0};       // A, after quench, vs original sites (fractional mapping)
    bool retained{false};
};

// Heat to T, then quench and compare to the starting relaxed structure.
inline ThermalCheck thermal_check(const Structure& relaxed, const Potential& pot, double T, int steps = 3000, unsigned seed = 7) {
    ThermalCheck C; C.T = T;
    // use a supercell with >= ~100 atoms so thermal fluctuations are meaningful
    int n = 1; while (int(relaxed.size()) * n * n * n < 96) ++n;
    Structure s0 = relaxed.supercell(n, n, n);
    const double e0 = evaluate(s0, pot).energy;
    MDOptions o; o.T = T; o.steps = steps; o.seed = seed;
    auto md = run_md(s0, pot, o);
    auto q = relax(md.final, pot, RelaxOptions{CellMode::Full, 2e-3, 4000});
    C.rms_disp = md.rms_disp;
    C.dE_quench_per_atom = (q.energy - e0) / double(s0.size());
    std::vector<Vec3> df(s0.size());
    Vec3 mean{};
    for (std::size_t i = 0; i < s0.size(); ++i) {
        Vec3 d = q.s.to_frac(q.s.pos[i]) - s0.to_frac(s0.pos[i]);
        for (int k = 0; k < 3; ++k) d[k] -= std::round(d[k]);
        df[i] = d; mean += d;
    }
    mean = mean * (1.0 / double(s0.size()));
    double mx = 0;
    for (auto& d : df) {
        Vec3 e = d - mean;
        for (int k = 0; k < 3; ++k) e[k] -= std::round(e[k]);
        mx = std::max(mx, norm(q.s.to_cart(e)));
    }
    C.max_site_shift = mx;
    C.retained = (C.dE_quench_per_atom > -0.005) && (mx < 0.5);
    return C;
}

// ---------------------------------------------------------------- binary convex hull
struct HullPoint { double x; double ef; int id; };

// Lower convex hull of (x, Ef); returns ids of vertices sorted by x. Endpoints (x=0,1, Ef=0) assumed present.
inline std::vector<HullPoint> lower_hull(std::vector<HullPoint> pts) {
    std::sort(pts.begin(), pts.end(), [](auto& a, auto& b) { return a.x < b.x || (a.x == b.x && a.ef < b.ef); });
    std::vector<HullPoint> h;
    for (const auto& p : pts) {
        if (!h.empty() && std::abs(h.back().x - p.x) < 1e-9) continue;   // keep the lowest at each x
        while (h.size() >= 2) {
            const auto& a = h[h.size() - 2]; const auto& b = h.back();
            const double cr = (b.x - a.x) * (p.ef - a.ef) - (b.ef - a.ef) * (p.x - a.x);
            if (cr <= 1e-12) h.pop_back(); else break;
        }
        h.push_back(p);
    }
    return h;
}
inline double hull_energy_at(const std::vector<HullPoint>& h, double x) {
    for (std::size_t k = 0; k + 1 < h.size(); ++k)
        if (x >= h[k].x - 1e-12 && x <= h[k + 1].x + 1e-12) {
            const double t = (h[k + 1].x - h[k].x) < 1e-12 ? 0 : (x - h[k].x) / (h[k + 1].x - h[k].x);
            return h[k].ef + t * (h[k + 1].ef - h[k].ef);
        }
    return 0.0;
}

} // namespace atomforge
