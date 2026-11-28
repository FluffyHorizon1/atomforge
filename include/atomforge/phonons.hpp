// atomforge — harmonic phonons by finite displacements
//
//   1. Build an n1 x n2 x n3 supercell of the (relaxed) cell. Displace each atom of the original cell by +-h along x, y, z
//      and take central differences of the forces on every supercell atom: Phi(i a; j b) = -dF_jb / du_ia  (eV/A^2).
//   2. Acoustic sum rule: the self term is reset so that sum_j Phi(i a; j b) = 0 exactly (rigid translations cost nothing).
//   3. Dynamical matrix D_{ia,kb}(q) = sum_{j -> k} Phi(i a; j b) w_ij exp(i q.(r_j - r_i)) / sqrt(m_i m_k), with r_j - r_i
//      taken over the shortest supercell images of the pair (weights w_ij = 1/number of equally short images), then made
//      Hermitian. Eigenvalues are omega^2; frequencies are reported in THz, negative for imaginary modes.
//
// At q-points commensurate with the supercell (q = sum_i m_i b_i / n_i) the result is exact for the given force constants,
// whatever their range. Elsewhere (a band path) it is an interpolation that improves with the supercell size.
// Dynamical stability: no frequency below -imag_tol_THz on the commensurate grid. The check counts as assessed only if the
// supercell has at least two cells along every lattice vector and faces at least one potential cutoff apart: a Gamma-only
// (1x1x1) check would always pass, because the sum rule zeroes the acoustic modes there. Imaginary modes found on any
// grid do mean "unstable".
//
// Every displacement is one Potential::evaluate on the supercell, called sequentially: 6 N evaluations for N atoms in the
// cell. That is cheap for EAM and LJ; through an ExternalPotential (DFT) it is the dominant cost.
#pragma once
#include "potential.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace atomforge {

namespace units {
inline constexpr double pi = 3.14159265358979323846;
// sqrt(eV / (A^2 amu)) in rad/s, divided by 2 pi and 1e12: THz per sqrt(eV/A^2/amu) = 15.633302...
inline const double THz_per_sqrt_eV_A2_amu = std::sqrt(1.602176634e-19 / 1.66053906660e-27) * 1e10 / (2 * units::pi) / 1e12;
} // namespace units

struct PhononOptions {
    int n[3] = {0, 0, 0};        // supercell repeats; 0 = automatic (each face width >= min_width)
    double min_width = 0;        // A; 0 -> 2 * potential cutoff, or 10 A if the potential reports no cutoff
    int max_repeat = 8;
    double h = 0.01;             // displacement, A
    double imag_tol_THz = 0.1;   // modes above -imag_tol count as real (finite-difference noise near Gamma)
};

struct ForceConstants {
    Structure cell, super;
    int n[3] = {1, 1, 1};
    std::vector<double> mass;            // per atom of `cell`
    std::vector<double> phi;             // (3 Np) x (3 Ns), row-major
    long evaluations = 0;
    double cutoff = 0;                   // the potential's cutoff (0 if it reports none)
    double asr_violation = 0;            // max |sum_j Phi(i a; j b)| before the sum rule is imposed, eV/A^2 (egg-box errors)
    // pair geometry: for each (i in cell, j in super) the shortest image vectors, with weight 1/count
    std::vector<std::vector<Vec3>> images;
    std::size_t np() const { return cell.size(); }
    std::size_t ns() const { return super.size(); }
};

inline ForceConstants force_constants(const Structure& cell, const Potential& pot, const PhononOptions& opt = {}) {
    ForceConstants FC;
    FC.cell = cell;
    const Vec3 w = face_widths(cell.cell);
    const double minw = opt.min_width > 0 ? opt.min_width : (pot.cutoff() > 0 ? 2 * pot.cutoff() : 10.0);
    for (int k = 0; k < 3; ++k) {
        FC.n[k] = opt.n[k] > 0 ? opt.n[k] : std::max(1, int(std::ceil(minw / w[k] - 1e-9)));
        FC.n[k] = std::min(FC.n[k], std::max(opt.max_repeat, opt.n[k]));
    }
    FC.super = cell.supercell(FC.n[0], FC.n[1], FC.n[2]);   // atoms 0..Np-1 are the original cell's, at the same positions
    const std::size_t Np = cell.size(), Ns = FC.super.size();
    FC.mass.resize(Np);
    for (std::size_t i = 0; i < Np; ++i) FC.mass[i] = pot.mass(cell.species[i]);
    FC.phi.assign(3 * Np * 3 * Ns, 0.0);
    NeighborList nl;
    for (std::size_t i = 0; i < Np; ++i) for (int a = 0; a < 3; ++a) {
        Structure p = FC.super, m = FC.super;
        p.pos[i][a] += opt.h; m.pos[i][a] -= opt.h;
        const auto fp = pot.evaluate(p, nl).forces;
        const auto fm = pot.evaluate(m, nl).forces;
        FC.evaluations += 2;
        if (fp.size() != Ns || fm.size() != Ns) throw std::runtime_error("force_constants: potential returned the wrong number of forces");
        for (std::size_t j = 0; j < Ns; ++j) for (int b = 0; b < 3; ++b)
            FC.phi[(3 * i + a) * 3 * Ns + 3 * j + b] = -(fp[j][b] - fm[j][b]) / (2 * opt.h);
    }
    // acoustic sum rule on the self term; the size of the correction is kept as a diagnostic
    FC.cutoff = pot.cutoff();
    for (std::size_t i = 0; i < Np; ++i) for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
        double s = 0;
        for (std::size_t j = 0; j < Ns; ++j) s += FC.phi[(3 * i + a) * 3 * Ns + 3 * j + b];
        FC.asr_violation = std::max(FC.asr_violation, std::abs(s));
        FC.phi[(3 * i + a) * 3 * Ns + 3 * i + b] -= s;
    }
    // shortest supercell images of every (i, j) pair
    FC.images.assign(Np * Ns, {});
    const Mat3& L = FC.super.cell;
    for (std::size_t i = 0; i < Np; ++i) for (std::size_t j = 0; j < Ns; ++j) {
        const Vec3 d0 = FC.super.pos[j] - FC.super.pos[i];
        std::vector<Vec3> cand;
        double best = 1e300;
        for (int x = -2; x <= 2; ++x) for (int y = -2; y <= 2; ++y) for (int z = -2; z <= 2; ++z) {
            const Vec3 d = d0 + L.row(0) * x + L.row(1) * y + L.row(2) * z;
            const double r = norm(d);
            if (r < best - 1e-6) { best = r; cand.assign(1, d); }
            else if (std::abs(r - best) <= 1e-6) cand.push_back(d);
        }
        FC.images[i * Ns + j] = cand;
    }
    return FC;
}

namespace detail {
// Eigenvalues of a real symmetric matrix (row-major n x n) by cyclic Jacobi.
inline std::vector<double> sym_eig(std::vector<double> a, int n) {
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0, diag = 0;
        for (int p = 0; p < n; ++p) { diag += a[p * n + p] * a[p * n + p]; for (int q = p + 1; q < n; ++q) off += a[p * n + q] * a[p * n + q]; }
        if (off <= 1e-30 * std::max(diag, 1e-300)) break;
        for (int p = 0; p < n; ++p) for (int q = p + 1; q < n; ++q) {
            const double apq = a[p * n + q];
            if (std::abs(apq) < 1e-300) continue;
            const double theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
            const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
            const double c = 1 / std::sqrt(t * t + 1), s = t * c;
            for (int k = 0; k < n; ++k) {
                const double akp = a[k * n + p], akq = a[k * n + q];
                a[k * n + p] = c * akp - s * akq; a[k * n + q] = s * akp + c * akq;
            }
            for (int k = 0; k < n; ++k) {
                const double apk = a[p * n + k], aqk = a[q * n + k];
                a[p * n + k] = c * apk - s * aqk; a[q * n + k] = s * apk + c * aqk;
            }
        }
    }
    std::vector<double> ev(n);
    for (int i = 0; i < n; ++i) ev[i] = a[i * n + i];
    std::sort(ev.begin(), ev.end());
    return ev;
}
} // namespace detail

// Hermitian dynamical matrix at q (Cartesian, 1/A, including the 2 pi), as a 3Np x 3Np complex matrix (row-major).
inline std::vector<std::complex<double>> dynamical_matrix(const ForceConstants& FC, const Vec3& q) {
    const std::size_t Np = FC.np(), Ns = FC.ns(), n = 3 * Np;
    std::vector<std::complex<double>> D(n * n, 0.0);
    for (std::size_t i = 0; i < Np; ++i) for (std::size_t j = 0; j < Ns; ++j) {
        const std::size_t k = j % Np;
        const auto& im = FC.images[i * Ns + j];
        std::complex<double> ph = 0;
        for (const auto& d : im) ph += std::polar(1.0, dot(q, d));
        ph /= double(im.size());
        const double s = 1.0 / std::sqrt(FC.mass[i] * FC.mass[k]);
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
            D[(3 * i + a) * n + 3 * k + b] += FC.phi[(3 * i + a) * 3 * Ns + 3 * j + b] * s * ph;
    }
    for (std::size_t r = 0; r < n; ++r) for (std::size_t c = r; c < n; ++c) {
        const auto h = 0.5 * (D[r * n + c] + std::conj(D[c * n + r]));
        D[r * n + c] = h; D[c * n + r] = std::conj(h);
    }
    return D;
}

// Frequencies at q in THz, ascending; imaginary modes are returned as negative numbers.
inline std::vector<double> frequencies_THz(const ForceConstants& FC, const Vec3& q) {
    const auto D = dynamical_matrix(FC, q);
    const int n = int(3 * FC.np()), m = 2 * n;
    // Hermitian H = A + iB  ->  real symmetric [[A, -B], [B, A]]; every eigenvalue appears twice.
    std::vector<double> R(std::size_t(m) * m);
    for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c) {
        const double A = D[r * n + c].real(), B = D[r * n + c].imag();
        R[r * m + c] = A; R[(r + n) * m + c + n] = A;
        R[r * m + c + n] = -B; R[(r + n) * m + c] = B;
    }
    const auto ev = detail::sym_eig(std::move(R), m);
    std::vector<double> f(n);
    for (int k = 0; k < n; ++k) {
        const double l = 0.5 * (ev[2 * k] + ev[2 * k + 1]);
        f[k] = (l >= 0 ? 1.0 : -1.0) * std::sqrt(std::abs(l)) * units::THz_per_sqrt_eV_A2_amu;
    }
    return f;
}

// Reciprocal lattice vectors (rows), including 2 pi.
inline Mat3 reciprocal(const Mat3& cell) { return transpose(inverse(cell)) * (2 * units::pi); }

// The q-points at which the supercell force constants are exact: q = sum_i (m_i / n_i) b_i, m_i = 0..n_i-1.
inline std::vector<Vec3> commensurate_qpoints(const ForceConstants& FC) {
    const Mat3 B = reciprocal(FC.cell.cell);
    std::vector<Vec3> q;
    for (int a = 0; a < FC.n[0]; ++a) for (int b = 0; b < FC.n[1]; ++b) for (int c = 0; c < FC.n[2]; ++c)
        q.push_back(B.row(0) * (double(a) / FC.n[0]) + B.row(1) * (double(b) / FC.n[1]) + B.row(2) * (double(c) / FC.n[2]));
    return q;
}

struct PhononCheck {
    bool computed = false;
    int n[3] = {0, 0, 0};
    int nq = 0;                  // commensurate q-points sampled
    double min_THz = 0, max_THz = 0;
    double gamma_acoustic_THz = 0;   // largest |frequency| of the three lowest modes at Gamma (after the sum rule: ~0)
    double asr_violation = 0;        // eV/A^2, before the sum rule (see ForceConstants::asr_violation)
    double min_width = 0;            // smallest face width of the supercell, A
    int n_imag = 0;              // modes below -imag_tol over all sampled q
    double imag_tol_THz = 0;
    bool assessed = false;       // >= 2 cells along every lattice vector and faces >= one cutoff apart
    bool stable = false;         // assessed and no imaginary mode
    long evaluations = 0;
};

// "dynamically stable", "dynamically UNSTABLE" (imaginary modes found, whatever the grid), or "not assessed (...)".
inline std::string phonon_status(const PhononCheck& P) {
    if (!P.computed) return "not computed";
    if (P.n_imag > 0) return "dynamically UNSTABLE";
    if (!P.assessed) return "not assessed (" + std::to_string(P.n[0]) + "x" + std::to_string(P.n[1]) + "x" + std::to_string(P.n[2]) +
                            " supercell too small; use >= 2 per direction and faces >= one cutoff apart)";
    return "dynamically stable";
}

inline PhononCheck phonon_check(const ForceConstants& FC, double imag_tol_THz = 0.1) {
    PhononCheck P;
    P.computed = true;
    for (int k = 0; k < 3; ++k) P.n[k] = FC.n[k];
    P.evaluations = FC.evaluations;
    P.imag_tol_THz = imag_tol_THz;
    P.min_THz = 1e300; P.max_THz = -1e300;
    for (const auto& q : commensurate_qpoints(FC)) {
        const auto f = frequencies_THz(FC, q);
        ++P.nq;
        if (norm(q) < 1e-12) for (int k = 0; k < 3 && k < int(f.size()); ++k) P.gamma_acoustic_THz = std::max(P.gamma_acoustic_THz, std::abs(f[k]));
        for (double x : f) { P.min_THz = std::min(P.min_THz, x); P.max_THz = std::max(P.max_THz, x); P.n_imag += (x < -imag_tol_THz); }
    }
    P.asr_violation = FC.asr_violation;
    const Vec3 w = face_widths(FC.super.cell);
    P.min_width = std::min({w.x, w.y, w.z});
    P.assessed = FC.n[0] >= 2 && FC.n[1] >= 2 && FC.n[2] >= 2 && (FC.cutoff <= 0 || P.min_width >= FC.cutoff);
    P.stable = P.assessed && P.n_imag == 0;
    return P;
}
inline PhononCheck phonon_check(const Structure& s, const Potential& pot, const PhononOptions& opt = {}) {
    return phonon_check(force_constants(s, pot, opt), opt.imag_tol_THz);
}

// Primitive cells of the cubic Bravais lattices (lattice parameter a of the conventional cube), one atom of species sp.
inline Structure primitive_fcc(double a, int sp = 0) {
    Structure s; s.label = "fcc-primitive";
    s.cell = Mat3::rows({0, a / 2, a / 2}, {a / 2, 0, a / 2}, {a / 2, a / 2, 0});
    s.pos = {{0, 0, 0}}; s.species = {sp};
    return s;
}
inline Structure primitive_bcc(double a, int sp = 0) {
    Structure s; s.label = "bcc-primitive";
    s.cell = Mat3::rows({-a / 2, a / 2, a / 2}, {a / 2, -a / 2, a / 2}, {a / 2, a / 2, -a / 2});
    s.pos = {{0, 0, 0}}; s.species = {sp};
    return s;
}

} // namespace atomforge
