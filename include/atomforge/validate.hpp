// atomforge — validation of the EAM element set against experiment (ground state + lattice constant)
#pragma once
#include "lattice.hpp"
#include "relax.hpp"
#include "search.hpp"
#include <cstdio>
#include <string>
#include <vector>

namespace atomforge {

// Reference data: experimental ground state, lattice constant (A) and cohesive energy (eV/atom).
struct ElementRef { const char* el; const char* gs; double a; double ec; };

// A lattice constant further than this from experiment (in %) counts as a validation failure,
// in addition to a wrong predicted ground state.
inline constexpr double validate_lattice_tol_pct = 1.0;

inline const std::vector<ElementRef>& reference_elements() {
    static const std::vector<ElementRef> refs = {
        {"Cu", "fcc", 3.615, 3.49}, {"Ag", "fcc", 4.086, 2.95}, {"Au", "fcc", 4.078, 3.81}, {"Ni", "fcc", 3.524, 4.44},
        {"Pd", "fcc", 3.891, 3.89}, {"Pt", "fcc", 3.924, 5.84}, {"Al", "fcc", 4.050, 3.39}, {"Pb", "fcc", 4.950, 2.03},
        {"Fe", "bcc", 2.867, 4.28}, {"Mo", "bcc", 3.147, 6.82}, {"Ta", "bcc", 3.301, 8.10}, {"W", "bcc", 3.165, 8.90},
        {"Cr", "bcc", 2.885, 4.10}, {"Mg", "hcp", 3.209, 1.51}, {"Co", "hcp", 2.507, 4.39}, {"Ti", "hcp", 2.951, 4.85},
        {"Zr", "hcp", 3.232, 6.25}};
    return refs;
}

struct ValidateSummary {
    int n{0};
    int gs_ok{0};
    int lattice_ok{0};
    double max_abs_err_pct{0};
    bool passed() const { return gs_ok == n && lattice_ok == n; }
};

// Relax every elemental prototype of each reference element, pick the lowest-energy parent lattice and compare with
// experiment. Prints the table to `out` (nullptr = silent). A passing run's printed table matches the historical
// `atomforge validate` output, including "ground state correct for N/N elements".
inline ValidateSummary validate_elements(const std::vector<ElementRef>& refs, std::FILE* out = stdout) {
    ValidateSummary S;
    S.n = int(refs.size());
    if (out) std::fprintf(out, "%-3s %-4s | %-9s %-9s %-6s | %-9s %-9s | %s\n", "El", "GS", "a_calc", "a_expt", "err%", "Ec_calc", "Ec_expt", "predicted GS");
    for (const auto& r : refs) {
        ZhouEAM pot({r.el});
        double best = 1e9; std::string bestp; double a_gs = 0;
        for (const auto& p : prototypes()) {
            if (p.nB) continue;
            Structure s = p.build(0, 0);
            detail::scale_to_volume(s, 0.74 * std::pow(ZhouEAM::lookup(r.el)->re, 3));
            auto rr = relax(s, pot);
            const double e = rr.energy / s.size();
            if (p.parent == r.gs) a_gs = (p.parent == "hcp") ? norm(rr.s.cell.row(0)) : std::cbrt(rr.s.volume());
            if (e < best - 1e-6) { best = e; bestp = p.parent; }
        }
        const bool gs_ok = bestp == r.gs;
        const double err_pct = 100 * (a_gs - r.a) / r.a;
        const bool lat_ok = a_gs > 0 && std::abs(err_pct) <= validate_lattice_tol_pct;
        S.gs_ok += gs_ok;
        S.lattice_ok += lat_ok;
        S.max_abs_err_pct = std::max(S.max_abs_err_pct, std::abs(err_pct));
        if (out) std::fprintf(out, "%-3s %-4s | %-9.4f %-9.4f %+6.2f | %-9.4f %-9.4f | %s %s%s\n", r.el, r.gs, a_gs, r.a, err_pct,
                              -best, r.ec, bestp.c_str(), gs_ok ? "OK" : "MISMATCH", lat_ok ? "" : " LATTICE-ERR");
    }
    if (out) std::fprintf(out, "ground state correct for %d/%d elements\n", S.gs_ok, S.n);
    return S;
}

} // namespace atomforge
