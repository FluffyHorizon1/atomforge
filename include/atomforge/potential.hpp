// atomforge — potential interface + Zhou/Johnson/Wadley (2004) analytic alloy EAM
#pragma once
#include "structure.hpp"
#include "zhou04_params.hpp"
#include <memory>
#include <string>
#include <vector>

namespace atomforge {

struct EvalResult {
    double energy{0};               // eV (total)
    std::vector<Vec3> forces;       // eV/A
    Mat3 virial;                    // W_ab = dE/d(eps_ab) for symmetric strain eps (eV). stress = W / V
};

// Any energy model. Swap this for an ML potential (MACE/CHGNet) or a DFT bridge later.
class Potential {
public:
    virtual ~Potential() = default;
    virtual EvalResult evaluate(const Structure& s, NeighborList& nl) const = 0;
    virtual double cutoff() const = 0;
    virtual const std::vector<std::string>& elements() const = 0;
    virtual double mass(int species) const = 0;
    virtual std::string name() const = 0;
    int index_of(const std::string& el) const {
        const auto& e = elements();
        for (std::size_t i = 0; i < e.size(); ++i) if (e[i] == el) return int(i);
        throw std::runtime_error("element not in potential: " + el);
    }
};

// Analytic EAM, X.W. Zhou, R.A. Johnson, H.N.G. Wadley, Phys. Rev. B 69, 144113 (2004).
// Cross pair terms use the paper's mixing rule:
//   phi_ab(r) = 1/2 [ f_b(r)/f_a(r) phi_aa(r) + f_a(r)/f_b(r) phi_bb(r) ]
// NOTE: cross terms are NOT fitted to alloy data -> trends are qualitative.
class ZhouEAM final : public Potential {
public:
    explicit ZhouEAM(const std::vector<std::string>& els) : names_(els) {
        for (const auto& el : els) {
            bool found = false;
            for (const auto& p : zhou04::table) if (el == p.name) { P_.push_back(p); found = true; break; }
            if (!found) throw std::runtime_error("ZhouEAM: no parameters for " + el + " (supported: " + supported() + ")");
        }
        double blat_max = 0;
        for (const auto& p : P_) blat_max = std::max(blat_max, std::sqrt(2.0) * p.re);
        rc_ = std::sqrt(10.0) / 2.0 * blat_max;   // same cutoff convention as LAMMPS setfl generator
    }

    static std::string supported() {
        std::string s; for (const auto& p : zhou04::table) { if (!s.empty()) s += ' '; s += p.name; } return s;
    }
    static const zhou04::Params* lookup(const std::string& el) {
        for (const auto& p : zhou04::table) if (el == p.name) return &p;
        return nullptr;
    }

    double cutoff() const override { return rc_; }
    const std::vector<std::string>& elements() const override { return names_; }
    double mass(int s) const override { return P_[s].mass; }
    std::string name() const override { return "Zhou04-EAM"; }

    // Generic radial form g(r) = c exp(-k (r/re - 1)) / (1 + (r/re - lam)^20), with derivative.
    static void radial(double r, double re, double c, double k, double lam, double& g, double& dg) {
        const double x = r / re;
        const double e = c * std::exp(-k * (x - 1.0));
        const double t = x - lam;
        const double t19 = std::pow(t, 19), t20 = t19 * t;
        const double den = 1.0 + t20;
        g = e / den;
        dg = g * (-k / re - 20.0 * t19 / (re * den));
    }
    void density(int s, double r, double& f, double& df) const {
        const auto& p = P_[s]; radial(r, p.re, p.fe, p.beta1, p.ramda1, f, df);
    }
    void pair_self(int s, double r, double& phi, double& dphi) const {
        const auto& p = P_[s];
        double g1, d1, g2, d2;
        radial(r, p.re, p.A, p.alpha, p.cai, g1, d1);
        radial(r, p.re, p.B, p.beta, p.ramda, g2, d2);
        phi = g1 - g2; dphi = d1 - d2;
    }
    void pair(int a, int b, double r, double& phi, double& dphi) const {
        if (a == b) { pair_self(a, r, phi, dphi); return; }
        double pa, dpa, pb, dpb, fa, dfa, fb, dfb;
        pair_self(a, r, pa, dpa); pair_self(b, r, pb, dpb);
        density(a, r, fa, dfa); density(b, r, fb, dfb);
        const double u = fb / fa, du = (dfb * fa - fb * dfa) / (fa * fa);
        const double iu = 1.0 / u, diu = -du / (u * u);
        phi = 0.5 * (u * pa + iu * pb);
        dphi = 0.5 * (du * pa + u * dpa + diu * pb + iu * dpb);
    }
    void embed(int s, double rho, double& F, double& dF) const {
        const auto& p = P_[s];
        if (rho <= 0) { F = 0; dF = 0; return; }
        const double rhoin = p.rhol * p.rhoe, rhoout = p.rhoh * p.rhoe;
        if (rho < rhoin) {
            const double d = rho / rhoin - 1.0;
            F = p.Fi0 + p.Fi1 * d + p.Fi2 * d * d + p.Fi3 * d * d * d;
            dF = (p.Fi1 + 2 * p.Fi2 * d + 3 * p.Fi3 * d * d) / rhoin;
        } else if (rho < rhoout) {
            const double d = rho / p.rhoe - 1.0;
            const double F3 = rho < p.rhoe ? p.Fm3 : p.Fm4;
            F = p.Fm0 + p.Fm1 * d + p.Fm2 * d * d + F3 * d * d * d;
            dF = (p.Fm1 + 2 * p.Fm2 * d + 3 * F3 * d * d) / p.rhoe;
        } else {
            const double x = rho / p.rhos, lx = std::log(x), xn = std::pow(x, p.fnn);
            F = p.Fn * (1.0 - p.fnn * lx) * xn;
            dF = -p.Fn * p.fnn * p.fnn * lx * xn / x / p.rhos;
        }
    }

    EvalResult evaluate(const Structure& s, NeighborList& nl) const override {
        if (nl.needs_rebuild(s) || std::abs(nl.cutoff - rc_) > 1e-12) nl.build(s, rc_);
        const int N = int(s.size());
        const double rc2 = rc_ * rc_;
        EvalResult R;
        R.forces.assign(N, Vec3{});
        std::vector<double> rho(N, 0.0), dF(N, 0.0);

        // pass 1: densities
        for (const auto& pr : nl.pairs) {
            const Vec3 d = s.pos[pr.j] + pr.shift - s.pos[pr.i];
            const double r2 = dot(d, d);
            if (r2 >= rc2) continue;
            const double r = std::sqrt(r2);
            double f, df;
            density(s.species[pr.j], r, f, df); rho[pr.i] += f;
            density(s.species[pr.i], r, f, df); rho[pr.j] += f;
        }
        for (int i = 0; i < N; ++i) { double F; embed(s.species[i], rho[i], F, dF[i]); R.energy += F; }

        // pass 2: pair energy, forces, virial
        for (const auto& pr : nl.pairs) {
            const Vec3 d = s.pos[pr.j] + pr.shift - s.pos[pr.i];
            const double r2 = dot(d, d);
            if (r2 >= rc2) continue;
            const double r = std::sqrt(r2);
            const int si = s.species[pr.i], sj = s.species[pr.j];
            double phi, dphi, fi, dfi, fj, dfj;
            pair(si, sj, r, phi, dphi);
            density(si, r, fi, dfi);
            density(sj, r, fj, dfj);
            R.energy += phi;
            const double dEdr = dphi + dF[pr.i] * dfj + dF[pr.j] * dfi;
            const Vec3 fvec = d * (dEdr / r);           // force on j is -fvec, on i is +fvec
            R.forces[pr.i] += fvec;
            R.forces[pr.j] -= fvec;
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) R.virial(a, b) += dEdr * d[a] * d[b] / r;
        }
        return R;
    }

private:
    std::vector<std::string> names_;
    std::vector<zhou04::Params> P_;
    double rc_{0};
};

} // namespace atomforge
