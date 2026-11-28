// atomforge — truncated-and-shifted Lennard-Jones pair potential
//
//   V(r) = 4 eps [ (sigma/r)^12 - (sigma/r)^6 ] - V_LJ(rc)   for r < rc, 0 otherwise
//
// One (eps, sigma) for every species pair, or per-species values mixed by the Lorentz-Berthelot rule
// eps_ab = sqrt(eps_a eps_b), sigma_ab = (sigma_a + sigma_b) / 2, each pair with its own cutoff shift. The shift makes the
// energy continuous at rc; the force has a small jump there, which is harmless for relaxation and for the cross-checks it exists for. Its main job is to cross-check the external
// worker bridge (quantum/worker.py lj implements the same formula), and to give a cheap non-EAM model for tests.
#pragma once
#include "potential.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace atomforge {

struct LJParams {
    double epsilon = 0.0104;   // eV (argon-like)
    double sigma = 3.40;       // A
    double rc = 8.5;           // A
    double mass = 39.948;      // amu, used for every species
    std::vector<double> epsilon_by_species, sigma_by_species;   // optional, one entry per element (override the above)
};

class LennardJones final : public Potential {
public:
    LennardJones(std::vector<std::string> els, const LJParams& p = {}) : names_(std::move(els)), p_(p) {
        if (names_.empty()) throw std::runtime_error("LennardJones: no elements");
        const std::size_t n = names_.size();
        auto per = [&](const std::vector<double>& v, double one, const char* what) {
            if (v.empty()) return std::vector<double>(n, one);
            if (v.size() != n) throw std::runtime_error(std::string("LennardJones: ") + what + " needs one value per element");
            return v;
        };
        const auto e1 = per(p_.epsilon_by_species, p_.epsilon, "epsilon_by_species");
        const auto s1 = per(p_.sigma_by_species, p_.sigma, "sigma_by_species");
        eps_.assign(n * n, 0); sig_.assign(n * n, 0); shift_.assign(n * n, 0);
        for (std::size_t a = 0; a < n; ++a) for (std::size_t b = 0; b < n; ++b) {
            const double e = std::sqrt(e1[a] * e1[b]), sg = 0.5 * (s1[a] + s1[b]);
            if (!(e > 0) || !(sg > 0) || !(p_.rc > sg)) throw std::runtime_error("LennardJones: need eps > 0, sigma > 0, rc > sigma");
            eps_[a * n + b] = e; sig_[a * n + b] = sg;
            const double sr6 = std::pow(sg / p_.rc, 6);
            shift_[a * n + b] = 4 * e * (sr6 * sr6 - sr6);
        }
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
            const std::size_t k = std::size_t(s.species[pr.i]) * names_.size() + std::size_t(s.species[pr.j]);
            const double eps = eps_[k], sg = sig_[k];
            const double sr2 = sg * sg / r2, sr6 = sr2 * sr2 * sr2, sr12 = sr6 * sr6;
            R.energy += 4 * eps * (sr12 - sr6) - shift_[k];
            const double dEdr_over_r = -24 * eps * (2 * sr12 - sr6) / r2;   // (dV/dr) / r
            const Vec3 fvec = d * dEdr_over_r;                                     // force on i is +fvec, on j is -fvec
            R.forces[pr.i] += fvec;
            R.forces[pr.j] -= fvec;
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) R.virial(a, b) += dEdr_over_r * d[a] * d[b];
        }
        return R;
    }

private:
    std::vector<std::string> names_;
    LJParams p_;
    std::vector<double> eps_, sig_, shift_;   // per species pair (row-major)
};

} // namespace atomforge
