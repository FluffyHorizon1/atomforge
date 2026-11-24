// atomforge — binary alloy structure search pipeline
//
//   1. elemental references (lowest of fcc/bcc/hcp for this potential)
//   2. known ordered prototypes in both A/B role assignments (L1_2, L1_0, L1_1, B2, D0_3, C15, B1)
//   3. random solid solutions on fcc/bcc/hcp (averaged over seeds)
//   4. simulated-annealing Monte Carlo of chemical ordering on fcc/bcc/hcp supercells
//      (finds orderings that are not in the prototype list)
//   5. formation energies -> convex hull -> energy above hull
//   6. stability screen of on/near-hull candidates: Born elastic + MD heat/quench
#pragma once
#include "analysis.hpp"
#include "lattice.hpp"
#include <functional>
#include <numeric>
#include <iostream>
#include <set>
#include <sstream>

namespace atomforge {

struct Candidate {
    int id{0};
    std::string source;         // "prototype" | "solid-solution" | "mc-ordering" | "element"
    std::string label;          // e.g. "L1_2 Cu3Au"
    std::string parent;
    Structure relaxed;
    double x{0};                // fraction of element B
    double e_per_atom{0};       // eV/atom
    double ef{0};               // formation energy eV/atom
    double e_hull{0};           // eV/atom above hull
    double vol_per_atom{0};
    bool converged{false};
    // stability screen (only for screened candidates)
    bool screened{false};
    ElasticResult elastic;
    ThermalCheck thermal;
    std::string verdict;
};

struct SearchOptions {
    std::vector<double> mc_fractions{0.125, 0.25, 1.0 / 3.0, 0.5, 2.0 / 3.0, 0.75, 0.875};
    std::vector<std::string> mc_parents{"fcc", "bcc", "hcp"};
    int mc_sweeps = 400;           // swap attempts per atom
    int ss_seeds = 3;
    double screen_window = 0.030;  // eV/atom above hull eligible for stability screening
    double md_T = 600;             // K
    int md_steps = 3000;
    bool verbose = true;
    unsigned seed = 2026;
};

struct SearchResult {
    std::string A, B, potential;
    double muA{0}, muB{0};
    std::string refA, refB;
    std::vector<Candidate> candidates;
    std::vector<HullPoint> hull;
};

namespace detail {
inline double vegard_volume(const std::string& A, const std::string& B, double xb) {
    const auto* pa = ZhouEAM::lookup(A); const auto* pb = ZhouEAM::lookup(B);
    const double va = 0.74 * std::pow(pa->re, 3), vb = 0.74 * std::pow(pb->re, 3);
    return (1 - xb) * va + xb * vb;
}
inline std::string formula(const std::string& A, const std::string& B, const Structure& s, int ia, int ib) {
    auto c = s.composition();
    int na = c.count(ia) ? c[ia] : 0, nb = c.count(ib) ? c[ib] : 0;
    int g = std::gcd(na, nb); if (g == 0) g = 1;
    std::ostringstream o;
    if (na) { o << A; if (na / g > 1) o << na / g; }
    if (nb) { o << B; if (nb / g > 1) o << nb / g; }
    return o.str();
}
// 48-site supercells on every parent lattice, so x = 1/8, 1/4, 1/3, 1/2, 2/3, 3/4, 7/8 are all exact.
inline Structure mc_supercell(const std::string& parent, int sp) {
    if (parent == "fcc") return parent_supercell("fcc", 1, sp).supercell(2, 2, 3);   // 4 x 12
    if (parent == "bcc") return parent_supercell("bcc", 1, sp).supercell(2, 3, 4);   // 2 x 24
    return parent_supercell("hcp", 1, sp).supercell(4, 3, 2);                        // 2 x 24
}
} // namespace detail

// Simulated annealing over site occupancies on a fixed (Vegard-scaled, isotropically relaxed) lattice.
inline Structure mc_ordering(Structure s, const Potential& pot, int ia, int ib, double xb, int sweeps, std::mt19937& rng) {
    randomize_occupancy(s, ia, ib, xb, rng);
    // isotropic volume relaxation at random occupancy gives a sensible lattice for the alloy
    {   // keep ideal sites; only adopt the relaxed volume
        const double v1 = relax(s, pot, RelaxOptions{CellMode::Hydrostatic, 5e-3, 500}).s.volume();
        const double f = std::cbrt(v1 / s.volume());
        s.deform(Mat3::diag(f, f, f));
    }
    const int N = int(s.size());
    NeighborList nl;
    double E = pot.evaluate(s, nl).energy;
    std::vector<int> As, Bs;
    for (int i = 0; i < N; ++i) (s.species[i] == ib ? Bs : As).push_back(i);
    if (As.empty() || Bs.empty()) return s;
    std::uniform_int_distribution<int> pa(0, int(As.size()) - 1), pb(0, int(Bs.size()) - 1);
    std::uniform_real_distribution<double> U(0, 1);
    const long total = long(sweeps) * N;
    const double T0 = 0.15, T1 = 0.002;   // eV per atom-swap scale (annealing schedule)
    Structure best = s; double Ebest = E;
    for (long step = 0; step < total; ++step) {
        const double kT = T0 * std::pow(T1 / T0, double(step) / total);
        const int a = pa(rng), b = pb(rng);
        std::swap(s.species[As[a]], s.species[Bs[b]]);
        const double En = pot.evaluate(s, nl).energy;
        if (En <= E || U(rng) < std::exp(-(En - E) / kT)) {
            E = En; std::swap(As[a], Bs[b]);
            if (E < Ebest) { Ebest = E; best = s; }
        } else {
            std::swap(s.species[As[a]], s.species[Bs[b]]);
        }
    }
    return best;
}

inline SearchResult binary_search(const std::string& A, const std::string& B, const SearchOptions& opt = {}) {
    ZhouEAM pot({A, B});
    const int ia = 0, ib = 1;
    SearchResult R; R.A = A; R.B = B; R.potential = pot.name();
    std::mt19937 rng(opt.seed);
    int next_id = 0;
    auto log = [&](const std::string& m) { if (opt.verbose) std::cerr << m << std::flush; };

    auto add = [&](const std::string& source, const std::string& label, const std::string& parent, Structure s) {
        auto r = relax(std::move(s), pot);
        Candidate c;
        c.id = next_id++; c.source = source; c.parent = parent;
        c.relaxed = r.s; c.converged = r.converged;
        c.x = r.s.fraction(ib);
        c.e_per_atom = r.energy / double(r.s.size());
        c.vol_per_atom = r.s.volume() / double(r.s.size());
        c.label = label + " " + detail::formula(A, B, r.s, ia, ib);
        c.relaxed.label = c.label;
        R.candidates.push_back(std::move(c));
        return R.candidates.size() - 1;
    };

    // 1. elements
    log("[1/6] elemental references\n");
    for (int which : {ia, ib}) {
        const std::string& el = (which == ia) ? A : B;
        double best = 1e9; std::string bestname;
        for (const auto& p : prototypes()) {
            if (p.nB) continue;
            Structure s = p.build(which, which);
            detail::scale_to_volume(s, detail::vegard_volume(A, B, which == ib ? 1.0 : 0.0));
            auto k = add("element", p.name, p.parent, s);
            if (R.candidates[k].e_per_atom < best) { best = R.candidates[k].e_per_atom; bestname = p.name; }
        }
        (which == ia ? R.muA : R.muB) = best;
        (which == ia ? R.refA : R.refB) = el + " " + bestname;
    }

    // 2. prototypes
    log("[2/6] ordered prototypes\n");
    for (const auto& p : prototypes()) {
        if (p.nB == 0) continue;
        for (int swapped = 0; swapped < 2; ++swapped) {
            const int a = swapped ? ib : ia, b = swapped ? ia : ib;
            if (p.nA == p.nB && swapped) continue;   // AB-symmetric prototype
            Structure s = p.build(a, b);
            detail::scale_to_volume(s, detail::vegard_volume(A, B, s.fraction(ib)));
            add("prototype", p.name.substr(0, p.name.find(' ')), p.parent, s);
        }
    }

    // 3. random solid solutions
    log("[3/6] random solid solutions\n");
    for (const std::string parent : {"fcc", "bcc", "hcp"})
        for (double x : {0.25, 0.5, 0.75}) {
            double sumE = 0; std::size_t best_k = 0; double bestE = 1e9;
            for (int sd = 0; sd < opt.ss_seeds; ++sd) {
                Structure s = detail::mc_supercell(parent, ia);
                randomize_occupancy(s, ia, ib, x, rng);
                detail::scale_to_volume(s, detail::vegard_volume(A, B, x));
                auto k = add("solid-solution", "SS-" + parent, parent, s);
                sumE += R.candidates[k].e_per_atom;
                if (R.candidates[k].e_per_atom < bestE) { bestE = R.candidates[k].e_per_atom; best_k = k; }
            }
            // keep a single averaged entry (the representative is the lowest seed; energy = mean)
            Candidate avg = R.candidates[best_k];
            avg.e_per_atom = sumE / opt.ss_seeds;
            avg.label = "SS-" + parent + " (avg of " + std::to_string(opt.ss_seeds) + ") " + detail::formula(A, B, avg.relaxed, ia, ib);
            R.candidates.resize(R.candidates.size() - opt.ss_seeds);
            avg.id = next_id++ ;
            R.candidates.push_back(avg);
        }

    // 4. MC ordering
    log("[4/6] Monte Carlo chemical-ordering search");
    for (const auto& parent : opt.mc_parents)
        for (double x : opt.mc_fractions) {
            Structure s = detail::mc_supercell(parent, ia);
            const double nb = x * s.size();
            if (std::abs(nb - std::round(nb)) > 1e-6) continue;
            detail::scale_to_volume(s, detail::vegard_volume(A, B, x));
            Structure best = mc_ordering(s, pot, ia, ib, x, opt.mc_sweeps, rng);
            add("mc-ordering", "MC-" + parent, parent, best);
            log(".");
        }
    log("\n");

    // 5. formation energies + hull
    log("[5/6] convex hull\n");
    std::vector<HullPoint> pts;
    for (auto& c : R.candidates) {
        c.ef = c.e_per_atom - ((1 - c.x) * R.muA + c.x * R.muB);
        pts.push_back({c.x, c.ef, c.id});
    }
    R.hull = lower_hull(pts);
    for (auto& c : R.candidates) c.e_hull = c.ef - hull_energy_at(R.hull, c.x);

    // 6. stability screen: lowest structure per composition within window, skipping elements
    log("[6/6] stability screen (elastic + MD at " + std::to_string(int(opt.md_T)) + " K)\n");
    std::map<long, int> best_at_x;
    for (std::size_t k = 0; k < R.candidates.size(); ++k) {
        const auto& c = R.candidates[k];
        if (c.source == "element" || c.e_hull > opt.screen_window) continue;
        const long key = std::lround(c.x * 1e6);
        if (!best_at_x.count(key) || c.e_hull < R.candidates[best_at_x[key]].e_hull) best_at_x[key] = int(k);
    }
    for (auto& [key, k] : best_at_x) {
        auto& c = R.candidates[k];
        log("   screening " + c.label + "\n");
        c.screened = true;
        c.elastic = elastic_tensor(c.relaxed, pot);
        c.thermal = thermal_check(c.relaxed, pot, opt.md_T, opt.md_steps);
        const bool on_hull = c.e_hull < 1e-4;
        const bool bound = c.ef < -1e-3;
        if (on_hull && bound && c.elastic.mechanically_stable && c.thermal.retained) c.verdict = "STABLE (on hull, mechanically + thermally stable)";
        else if (bound && c.elastic.mechanically_stable && c.thermal.retained) c.verdict = "METASTABLE (near hull)";
        else if (!c.elastic.mechanically_stable) c.verdict = "UNSTABLE (fails Born criteria)";
        else if (!c.thermal.retained) c.verdict = "UNSTABLE (transforms under MD)";
        else c.verdict = "NOT BOUND (Ef >= 0, phase-separates)";
    }
    return R;
}

} // namespace atomforge
