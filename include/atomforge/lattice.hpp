// atomforge — crystal prototype builders (header-only, C++17)
//
// Every builder makes a cell in units of a=1, then scales it to the requested volume/atom.
// Species roles: 'A' -> species index a, 'B' -> species index b.
#pragma once
#include "structure.hpp"
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace atomforge {

struct Prototype {
    std::string name;      // Strukturbericht / common name
    std::string formula;   // stoichiometry pattern, e.g. "A3B"
    int nA, nB;            // atoms of each role in the cell
    std::string parent;    // parent lattice: fcc / bcc / hcp / other
    std::function<Structure(int a, int b)> build;   // unit-a cell
};

namespace detail {
inline Structure cubic(double a = 1.0) { Structure s; s.cell = Mat3::diag(a, a, a); return s; }
inline void scale_to_volume(Structure& s, double vol_per_atom) {
    const double f = std::cbrt(vol_per_atom * s.size() / s.volume());
    s.deform(Mat3::diag(f, f, f));
}
} // namespace detail

inline std::vector<Prototype> prototypes() {
    using detail::cubic;
    std::vector<Prototype> P;
    P.push_back({"A1 (fcc)", "A", 4, 0, "fcc", [](int a, int) {
        Structure s = cubic(); s.label = "A1-fcc";
        for (auto f : {Vec3{0, 0, 0}, Vec3{0, .5, .5}, Vec3{.5, 0, .5}, Vec3{.5, .5, 0}}) s.add_frac(f, a);
        return s; }});
    P.push_back({"A2 (bcc)", "A", 2, 0, "bcc", [](int a, int) {
        Structure s = cubic(); s.label = "A2-bcc";
        s.add_frac({0, 0, 0}, a); s.add_frac({.5, .5, .5}, a); return s; }});
    P.push_back({"A3 (hcp)", "A", 2, 0, "hcp", [](int a, int) {
        Structure s; s.label = "A3-hcp";
        const double c = std::sqrt(8.0 / 3.0);
        s.cell = Mat3::rows({1, 0, 0}, {-0.5, std::sqrt(3.0) / 2, 0}, {0, 0, c});
        s.add_frac({0, 0, 0}, a); s.add_frac({1.0 / 3, 2.0 / 3, 0.5}, a); return s; }});
    P.push_back({"L1_2 (Cu3Au)", "A3B", 3, 1, "fcc", [](int a, int b) {
        Structure s = cubic(); s.label = "L1_2";
        s.add_frac({0, 0, 0}, b);
        for (auto f : {Vec3{0, .5, .5}, Vec3{.5, 0, .5}, Vec3{.5, .5, 0}}) s.add_frac(f, a);
        return s; }});
    P.push_back({"L1_0 (CuAu)", "AB", 2, 2, "fcc", [](int a, int b) {
        Structure s = cubic(); s.label = "L1_0";
        s.add_frac({0, 0, 0}, b); s.add_frac({.5, .5, 0}, b);
        s.add_frac({.5, 0, .5}, a); s.add_frac({0, .5, .5}, a); return s; }});
    P.push_back({"L1_1 (CuPt)", "AB", 1, 1, "fcc", [](int a, int b) {
        // alternating (111) planes; cell from fcc primitive vectors p1-p2, p2-p3, p1+p2
        Structure s; s.label = "L1_1";
        const Vec3 p1{0, .5, .5}, p2{.5, 0, .5}, p3{.5, .5, 0};
        s.cell = Mat3::rows(p1 - p2, p2 - p3, p1 + p2);
        s.pos = {Vec3{0, 0, 0}, p1}; s.species = {a, b}; s.wrap(); return s; }});
    P.push_back({"B2 (CsCl)", "AB", 1, 1, "bcc", [](int a, int b) {
        Structure s = cubic(); s.label = "B2";
        s.add_frac({0, 0, 0}, a); s.add_frac({.5, .5, .5}, b); return s; }});
    P.push_back({"D0_3 (Fe3Al)", "A3B", 12, 4, "bcc", [](int a, int b) {
        Structure s = cubic(); s.label = "D0_3";
        for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) {
            s.add_frac({i * .5, j * .5, k * .5}, ((i + j + k) % 2 == 0) ? b : a);
            s.add_frac({i * .5 + .25, j * .5 + .25, k * .5 + .25}, a);
        }
        return s; }});
    P.push_back({"C15 (MgCu2 Laves)", "AB2", 8, 16, "other", [](int a, int b) {
        // A on 8a (diamond), B on 16d; Fd-3m origin choice 1. Here A = larger atom role.
        Structure s = cubic(); s.label = "C15";
        const Vec3 T[4] = {{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}};
        const Vec3 s8a[2] = {{0, 0, 0}, {.75, .25, .75}};
        const Vec3 s16d[4] = {{.625, .625, .625}, {.375, .875, .125}, {.875, .125, .375}, {.125, .375, .875}};
        for (auto t : T) { for (auto f : s8a) s.add_frac(f + t, a); for (auto f : s16d) s.add_frac(f + t, b); }
        s.wrap(); return s; }});
    P.push_back({"B1 (NaCl)", "AB", 4, 4, "other", [](int a, int b) {
        Structure s = cubic(); s.label = "B1";
        for (auto f : {Vec3{0, 0, 0}, Vec3{0, .5, .5}, Vec3{.5, 0, .5}, Vec3{.5, .5, 0}}) {
            s.add_frac(f, a); s.add_frac(f + Vec3{.5, 0, 0}, b);
        }
        s.wrap(); return s; }});
    return P;
}

// Build a parent lattice (fcc/bcc/hcp) supercell with all sites = species a.
inline Structure parent_supercell(const std::string& parent, int n, int sp) {
    for (const auto& p : prototypes())
        if (p.parent == parent && p.nB == 0) return p.build(sp, sp).supercell(n, n, n);
    throw std::runtime_error("unknown parent lattice " + parent);
}

// Random substitutional alloy: exactly round(x*N) sites become species b.
inline void randomize_occupancy(Structure& s, int a, int b, double xb, std::mt19937& rng) {
    const int N = int(s.size()), nb = int(std::lround(xb * N));
    std::vector<int> idx(N);
    for (int i = 0; i < N; ++i) idx[i] = i;
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int i = 0; i < N; ++i) s.species[idx[i]] = (i < nb) ? b : a;
}

} // namespace atomforge
