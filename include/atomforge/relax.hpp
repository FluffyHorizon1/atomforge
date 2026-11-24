// atomforge — FIRE relaxation of atoms and (optionally) the cell, at external pressure
//
// Cell DOF follow the "unit cell filter" idea: positions r = D r0, generalized cell coordinate
// X = c*D with c = N, generalized cell force G = -(W + pV I) D^-T / c (W = virial dE/d eps).
#pragma once
#include "potential.hpp"
#include <vector>

namespace atomforge {

enum class CellMode { Fixed, Hydrostatic, Full };

struct RelaxOptions {
    CellMode cell = CellMode::Full;
    double fmax = 1e-3;         // eV/A on atoms; eV per atom on cell DOF
    int max_steps = 3000;
    double pressure_GPa = 0.0;
    double max_move = 0.15;     // A per step
};

struct RelaxResult {
    Structure s;
    double energy{0};           // eV total
    double enthalpy{0};         // eV total (E + pV)
    int steps{0};
    bool converged{false};
    Mat3 stress_GPa;            // W/V converted; tensile positive
};

inline RelaxResult relax(Structure s, const Potential& pot, const RelaxOptions& opt = {}) {
    const int N = int(s.size());
    const double cf = double(N);
    const double p = opt.pressure_GPa / units::eV_per_A3_to_GPa;   // eV/A^3
    const Mat3 H0 = s.cell;
    std::vector<Vec3> r0 = s.pos;
    Mat3 D = Mat3::identity();
    NeighborList nl;

    const int ncell = (opt.cell == CellMode::Fixed) ? 0 : 9;
    const int ndof = 3 * N + ncell;
    std::vector<double> v(ndof, 0.0), g(ndof, 0.0);

    auto apply = [&]() {
        s.cell = transpose(D * transpose(H0));
        for (int i = 0; i < N; ++i) s.pos[i] = D * r0[i];
    };
    EvalResult ev;
    double H = 0;
    auto compute = [&]() {
        apply();
        ev = pot.evaluate(s, nl);
        const double V = s.volume();
        H = ev.energy + p * V;
        const Mat3 DT = transpose(D);
        for (int i = 0; i < N; ++i) {
            const Vec3 f = DT * ev.forces[i];
            g[3 * i] = f.x; g[3 * i + 1] = f.y; g[3 * i + 2] = f.z;
        }
        if (ncell) {
            Mat3 Wp = ev.virial + Mat3::identity() * (p * V);
            Mat3 G = (Wp * transpose(inverse(D))) * (-1.0 / cf);
            if (opt.cell == CellMode::Hydrostatic) {
                const double tr = (G(0, 0) + G(1, 1) + G(2, 2)) / 3.0;
                G = Mat3::identity() * tr;
            }
            for (int k = 0; k < 9; ++k) g[3 * N + k] = G.m[k];
        }
    };
    auto fmax_now = [&]() {
        double m = 0;
        for (int i = 0; i < N; ++i) m = std::max(m, std::sqrt(g[3*i]*g[3*i] + g[3*i+1]*g[3*i+1] + g[3*i+2]*g[3*i+2]));
        for (int k = 0; k < ncell; ++k) m = std::max(m, std::abs(g[3 * N + k]));
        return m;
    };

    // FIRE (Bitzek et al., PRL 97, 170201 (2006))
    double dt = 0.1, alpha = 0.1;
    const double dtmax = 1.0, finc = 1.1, fdec = 0.5, alpha0 = 0.1, falpha = 0.99;
    int npos = 0;
    RelaxResult res;
    compute();
    int step = 0;
    for (; step < opt.max_steps; ++step) {
        if (fmax_now() < opt.fmax) { res.converged = true; break; }
        double P = 0, vn = 0, fn = 0;
        for (int k = 0; k < ndof; ++k) { P += g[k] * v[k]; vn += v[k] * v[k]; fn += g[k] * g[k]; }
        vn = std::sqrt(vn); fn = std::sqrt(fn);
        if (P > 0) {
            for (int k = 0; k < ndof; ++k) v[k] = (1 - alpha) * v[k] + alpha * g[k] / (fn + 1e-300) * vn;
            if (++npos > 5) { dt = std::min(dt * finc, dtmax); alpha *= falpha; }
        } else {
            std::fill(v.begin(), v.end(), 0.0);
            dt *= fdec; alpha = alpha0; npos = 0;
        }
        for (int k = 0; k < ndof; ++k) v[k] += dt * g[k];
        // step, capped
        double dmax = 0;
        for (int k = 0; k < ndof; ++k) dmax = std::max(dmax, std::abs(dt * v[k]));
        const double scale = dmax > opt.max_move ? opt.max_move / dmax : 1.0;
        for (int i = 0; i < N; ++i) {
            r0[i].x += scale * dt * v[3 * i]; r0[i].y += scale * dt * v[3 * i + 1]; r0[i].z += scale * dt * v[3 * i + 2];
        }
        if (ncell) {
            for (int k = 0; k < 9; ++k) D.m[k] += scale * dt * v[3 * N + k] / cf;
        }
        compute();
    }
    res.steps = step;
    res.energy = ev.energy;
    res.enthalpy = H;
    res.stress_GPa = ev.virial * (units::eV_per_A3_to_GPa / s.volume());
    s.wrap();
    res.s = s;
    return res;
}

// Single-point evaluation helper.
inline EvalResult evaluate(const Structure& s, const Potential& pot) { NeighborList nl; return pot.evaluate(s, nl); }

} // namespace atomforge
