// atomforge — truncated-and-shifted Lennard-Jones pair potential
//
//   V(r) = 4 eps [ (sigma/r)^12 - (sigma/r)^6 ] - V_LJ(rc)   for r < rc, 0 otherwise
//
// One (eps, sigma) for every species pair. The shift makes the energy continuous at rc; the force has a small jump there,
// which is harmless for relaxation and for the cross-checks it exists for. Its main job is to cross-check the external
// worker bridge (quantum/worker.py lj implements the same formula), and to give a cheap non-EAM model for tests.
#pragma once
#include "potential.hpp"
#include <stdexcept>
#include <string>
#include <vector>

namespace atomforge {

struct LJParams {
    double epsilon = 0.0104;   // eV (argon-like)
    double sigma = 3.40;       // A
    double rc = 8.5;           // A
    double mass = 39.948;      // amu, used for every species
};

class LennardJones final : public Potential {
public:
    LennardJones(std::vector<std::string> els, const LJParams& p = {}) : names_(std::move(els)), p_(p) {
        if (names_.empty()) throw std::runtime_error("LennardJones: no elements");
        if (!(p_.epsilon > 0) || !(p_.sigma > 0) || !(p_.rc > p_.sigma)) throw std::runtime_error("LennardJones: need eps > 0, sigma > 0, rc > sigma");
        shift_ = raw(p_.rc);
    }
    double cutoff() const override { return p_.rc; }
    const std::vector<std::string>& elements() const override { return names_; }
    double mass(int) const override { return p_.mass; }
    std::string name() const override { return "LJ(eps=" + std::to_string(p_.epsilon) + ",sigma=" + std::to_string(p_.sigma) + ",rc=" + std::to_string(p_.rc) + ")"; }
    const LJParams& params() const { return p_; }

    EvalResult evaluate(const Structure& s, NeighborList& nl) const override {
        if (nl.needs_rebuild(s) || std::abs(nl.cutoff - p_.rc) > 1e-12) nl.build(s, p_.rc);
        const double rc2 = p_.rc * p_.rc;
        EvalResult R;
        R.forces.assign(s.size(), Vec3{});
        for (const auto& pr : nl.pairs) {
            const Vec3 d = s.pos[pr.j] + pr.shift - s.pos[pr.i];
            const double r2 = dot(d, d);
            if (r2 >= rc2) continue;
            const double sr2 = p_.sigma * p_.sigma / r2, sr6 = sr2 * sr2 * sr2, sr12 = sr6 * sr6;
            R.energy += 4 * p_.epsilon * (sr12 - sr6) - shift_;
            const double dEdr_over_r = -24 * p_.epsilon * (2 * sr12 - sr6) / r2;   // (dV/dr) / r
            const Vec3 fvec = d * dEdr_over_r;                                     // force on i is +fvec, on j is -fvec
            R.forces[pr.i] += fvec;
            R.forces[pr.j] -= fvec;
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) R.virial(a, b) += dEdr_over_r * d[a] * d[b];
        }
        return R;
    }

private:
    double raw(double r) const { const double sr6 = std::pow(p_.sigma / r, 6); return 4 * p_.epsilon * (sr6 * sr6 - sr6); }
    std::vector<std::string> names_;
    LJParams p_;
    double shift_{0};
};

} // namespace atomforge
