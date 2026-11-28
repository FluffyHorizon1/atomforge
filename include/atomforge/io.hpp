// atomforge — structure + report output
#pragma once
#include "search.hpp"
#include <fstream>
#include <stdexcept>
#include <iomanip>

namespace atomforge {

inline void check_write(std::ofstream& f, const std::string& path) {
    f.flush();
    if (!f) throw std::runtime_error("cannot write " + path);
    f.close();
    if (!f) throw std::runtime_error("cannot write " + path);
}


// VASP POSCAR (feed straight into DFT for verification: VASP / Quantum ESPRESSO via ASE / GPAW)
inline void write_poscar(const std::string& path, const Structure& s, const std::vector<std::string>& els) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << s.label << "\n1.0\n" << std::fixed << std::setprecision(10);
    for (int i = 0; i < 3; ++i) f << "  " << s.cell(i, 0) << " " << s.cell(i, 1) << " " << s.cell(i, 2) << "\n";
    std::vector<int> present;
    for (std::size_t e = 0; e < els.size(); ++e) for (int sp : s.species) if (sp == int(e)) { present.push_back(int(e)); break; }
    for (int e : present) f << " " << els[e];
    f << "\n";
    for (int e : present) { int n = 0; for (int sp : s.species) n += (sp == e); f << " " << n; }
    f << "\nDirect\n";
    for (int e : present)
        for (std::size_t i = 0; i < s.size(); ++i)
            if (s.species[i] == e) { Vec3 fr = s.to_frac(s.pos[i]); f << "  " << fr.x << " " << fr.y << " " << fr.z << "\n"; }
    check_write(f, path);
}

// Extended XYZ (opens in OVITO / ASE / VESTA)
inline void write_extxyz(const std::string& path, const Structure& s, const std::vector<std::string>& els, double energy = 0) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << s.size() << "\n" << std::fixed << std::setprecision(8) << "Lattice=\"";
    for (int i = 0; i < 9; ++i) f << s.cell.m[i] << (i < 8 ? " " : "");
    f << "\" Properties=species:S:1:pos:R:3 energy=" << energy << " pbc=\"T T T\" label=\"" << s.label << "\"\n";
    for (std::size_t i = 0; i < s.size(); ++i)
        f << els[s.species[i]] << " " << s.pos[i].x << " " << s.pos[i].y << " " << s.pos[i].z << "\n";
    check_write(f, path);
}

inline std::string json_escape(const std::string& s) {
    std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; } return o;
}

inline void write_report_json(const std::string& path, const SearchResult& R) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << std::setprecision(6) << "{\n  \"system\": [\"" << R.A << "\", \"" << R.B << "\"],\n"
      << "  \"potential\": \"" << R.potential << "\",\n"
      << "  \"reference\": {\"" << R.A << "\": {\"phase\": \"" << json_escape(R.refA) << "\", \"mu\": " << R.muA << "}, \""
      << R.B << "\": {\"phase\": \"" << json_escape(R.refB) << "\", \"mu\": " << R.muB << "}},\n  \"hull\": [";
    for (std::size_t k = 0; k < R.hull.size(); ++k) f << (k ? ", " : "") << "{\"x\": " << R.hull[k].x << ", \"ef\": " << R.hull[k].ef << ", \"id\": " << R.hull[k].id << "}";
    f << "],\n  \"candidates\": [\n";
    for (std::size_t k = 0; k < R.candidates.size(); ++k) {
        const auto& c = R.candidates[k];
        f << "    {\"id\": " << c.id << ", \"label\": \"" << json_escape(c.label) << "\", \"source\": \"" << c.source
          << "\", \"parent\": \"" << c.parent << "\", \"natoms\": " << c.relaxed.size() << ", \"x\": " << c.x
          << ", \"e_per_atom\": " << c.e_per_atom << ", \"ef\": " << c.ef << ", \"e_hull\": " << c.e_hull
          << ", \"vol_per_atom\": " << c.vol_per_atom << ", \"converged\": " << (c.converged ? "true" : "false");
        if (c.screened) {
            f << ", \"screen\": {\"verdict\": \"" << json_escape(c.verdict) << "\", \"bulk_GPa\": " << c.elastic.bulk_voigt
              << ", \"shear_GPa\": " << c.elastic.shear_voigt << ", \"min_C_eig_GPa\": " << c.elastic.eigen[0]
              << ", \"C11\": " << c.elastic.C[0][0] << ", \"C12\": " << c.elastic.C[0][1] << ", \"C44\": " << c.elastic.C[3][3]
              << ", \"md_T\": " << c.thermal.T << ", \"md_rms_disp\": " << c.thermal.rms_disp
              << ", \"dE_quench\": " << c.thermal.dE_quench_per_atom << ", \"max_site_shift\": " << c.thermal.max_site_shift;
            if (c.phonons.computed) {
                const auto& P = c.phonons;
                f << ", \"phonons\": {\"supercell\": [" << P.n[0] << ", " << P.n[1] << ", " << P.n[2] << "], \"n_q\": " << P.nq
                  << ", \"min_THz\": " << P.min_THz << ", \"max_THz\": " << P.max_THz << ", \"gamma_acoustic_THz\": " << P.gamma_acoustic_THz
                  << ", \"n_imaginary\": " << P.n_imag << ", \"imag_tol_THz\": " << P.imag_tol_THz
                  << ", \"asr_violation_eV_A2\": " << P.asr_violation << ", \"min_width_A\": " << P.min_width
                  << ", \"assessed\": " << (P.assessed ? "true" : "false")
                  << ", \"dynamically_stable\": " << (P.stable ? "true" : "false") << ", \"evaluations\": " << P.evaluations << "}";
            }
            f << "}";
        }
        f << "}" << (k + 1 < R.candidates.size() ? "," : "") << "\n";
    }
    f << "  ]\n}\n";
    check_write(f, path);
}

} // namespace atomforge
