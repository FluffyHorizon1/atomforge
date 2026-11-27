// atomforge — forge.materials.v0 sidecar writer
//
// Implements the sidecar proposed in qpuforge docs/DATA_CONTRACT.md (branch docs/materials-data-contract,
// commit f64238f, section 3.2), plus one additive optional block (`atomforge_computed`) that the contract
// allows on v0 ("Additive optional fields stay on v0").
//
// What is REAL (computed by atomforge's EAM energy model from a SearchResult):
//   produced_by, system, structures_ref, and atomforge_computed (elemental references, hull phases:
//   energies, formation energies, cell lengths/angles, Born elastic constants).
// What is a PLACEHOLDER (atomforge cannot compute it today): quantities.Q_dielectric,
//   quantities.loss_tangent, quantities.junction_gap_ueV. These carry value=null, status="placeholder"
//   and are listed in the top-level `placeholders` array. No physics number is invented for them;
//   the qpuforge default is quoted only under `qpuforge_default_for_reference`.
#pragma once
#include "search.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef ATOMFORGE_GIT_COMMIT
#define ATOMFORGE_GIT_COMMIT "unknown"
#endif

namespace atomforge {

struct MaterialsOptions {
    std::string commit = ATOMFORGE_GIT_COMMIT;   // atomforge git commit (set by CMake; "unknown" if built outside git)
    std::string timestamp;                       // ISO-8601 with offset; empty = now (UTC, "+00:00")
    std::string report_path = "report.json";     // path of the search report, relative to the sidecar
};

namespace materials_detail {

inline std::string now_iso8601_utc() {
    const std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{};
#ifdef _WIN32
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[40];
    std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%S+00:00", &tm);
    return buf;
}

inline std::string num(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream o; o.imbue(std::locale::classic());
    o.precision(10);   // 10 significant digits: round-trips the values atomforge actually computes to well below FD noise
    o << v;
    return o.str();
}
inline std::string str(const std::string& s) {
    std::string o = "\"";
    for (unsigned char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); o += b; }
                else o += char(c);   // UTF-8 passes through unchanged
        }
    }
    return o + "\"";
}

struct CellGeom { double a, b, c, alpha, beta, gamma; };
inline CellGeom cell_geometry(const Mat3& cell) {
    const Vec3 v0 = cell.row(0), v1 = cell.row(1), v2 = cell.row(2);
    const double la = norm(v0), lb = norm(v1), lc = norm(v2);
    const double rad = 180.0 / 3.14159265358979323846;
    auto ang = [&](const Vec3& p, const Vec3& q, double lp, double lq) {
        double c = dot(p, q) / (lp * lq);
        c = c > 1 ? 1 : (c < -1 ? -1 : c);
        return std::acos(c) * rad;
    };
    return {la, lb, lc, ang(v1, v2, lb, lc), ang(v0, v2, la, lc), ang(v0, v1, la, lb)};
}

inline std::string cell_json(const Structure& s, const char* ind) {
    const CellGeom g = cell_geometry(s.cell);
    std::ostringstream o;
    o << "{\n"
      << ind << "  \"natoms\": " << s.size() << ",\n"
      << ind << "  \"cell_lengths_A\": [" << num(g.a) << ", " << num(g.b) << ", " << num(g.c) << "],\n"
      << ind << "  \"cell_angles_deg\": [" << num(g.alpha) << ", " << num(g.beta) << ", " << num(g.gamma) << "],\n"
      << ind << "  \"note\": \"Edges and angles of the stored (relaxed) cell. For prototypes and Monte Carlo supercells this is the simulation cell, not a primitive cell.\"\n"
      << ind << "}";
    return o.str();
}

inline const Candidate* find_by_id(const SearchResult& R, int id) {
    for (const auto& c : R.candidates) if (c.id == id) return &c;
    return nullptr;
}

// Elemental reference structure: the "element" candidate at the pure composition whose energy defines mu.
inline const Candidate* find_reference(const SearchResult& R, double x, double mu) {
    for (const auto& c : R.candidates)
        if (c.source == "element" && std::abs(c.x - x) < 1e-9 && c.e_per_atom == mu) return &c;
    return nullptr;
}

inline std::string placeholder_quantity(const std::string& unit, const std::string& maps_to, const std::string& extra,
                                        const std::string& default_ref, const std::string& reason) {
    std::ostringstream o;
    o << "{\n"
      << "      \"value\": null,\n"
      << "      \"unit\": " << str(unit) << ",\n"
      << "      \"maps_to\": " << maps_to << ",\n"
      << "      \"status\": \"placeholder\",\n"
      << "      \"provenance\": \"placeholder\",\n"
      << "      \"reason\": " << str(reason) << ",\n"
      << extra
      << "      \"qpuforge_default_for_reference\": " << default_ref << "\n"
      << "    }";
    return o.str();
}

} // namespace materials_detail

// Serialise the forge.materials.v0 sidecar for a finished search. Deterministic for a fixed
// (SearchResult, MaterialsOptions) pair; with an empty `timestamp` the only varying field is the timestamp.
inline std::string materials_contract_json(const SearchResult& R, const MaterialsOptions& opt = {}) {
    using namespace materials_detail;
    std::ostringstream o;
    o << "{\n  \"schema\": \"forge.materials.v0\",\n";
    o << "  \"produced_by\": {\n"
      << "    \"tool\": \"atomforge\",\n"
      << "    \"commit\": " << str(opt.commit) << ",\n"
      << "    \"potential_or_method\": " << str(R.potential) << ",\n"
      << "    \"timestamp\": " << str(opt.timestamp.empty() ? now_iso8601_utc() : opt.timestamp) << "\n  },\n";
    o << "  \"system\": [" << str(R.A) << ", " << str(R.B) << "],\n";
    o << "  \"structures_ref\": {\n"
      << "    \"report\": " << str(opt.report_path) << ",\n"
      << "    \"format_notes\": \"POSCAR in Å, extxyz energy in eV (total), report energies in eV/atom\"\n  },\n";

    // --- contract quantities: all three are placeholders today -----------------------------------
    const std::string eam_reason = "atomforge computes structures and energies with an analytic EAM potential; it has no dielectric-loss model, "
                                   "and EAM cannot represent O, N or Si. Do not use this value.";
    o << "  \"quantities\": {\n";
    o << "    \"Q_dielectric\": "
      << placeholder_quantity("dimensionless", "\"ChipDesign.materials.Q_dielectric\"",
                              "      \"conditions\": {\"T_K\": null, \"f_GHz\": null},\n", "3000000", eam_reason)
      << ",\n";
    o << "    \"loss_tangent\": "
      << placeholder_quantity("dimensionless", "null",
                              "      \"relation\": \"If the uniform-loss model in device.hpp holds, Q_dielectric = 1 / loss_tangent.\",\n",
                              "null", eam_reason)
      << ",\n";
    o << "    \"junction_gap_ueV\": "
      << placeholder_quantity("ueV", "\"ChipDesign.materials.junction_gap_ueV\"", "", "180",
                              "The superconducting gap is not a ground-state total-energy quantity and atomforge does not compute it. "
                              "Reserved in forge.materials.v0.")
      << "\n  },\n";

    o << "  \"placeholders\": [\"quantities.Q_dielectric\", \"quantities.loss_tangent\", \"quantities.junction_gap_ueV\"],\n";
    o << "  \"contains_placeholders\": true,\n";

    // --- additive optional block: real atomforge results ----------------------------------------
    o << "  \"atomforge_computed\": {\n"
      << "    \"status\": \"computed\",\n"
      << "    \"method\": " << str(R.potential) << ",\n"
      << "    \"fidelity\": \"classical-eam\",\n"
      << "    \"units\": {\"energy\": \"eV/atom\", \"volume\": \"A^3/atom\", \"length\": \"A\", \"angle\": \"deg\", \"elastic\": \"GPa\"},\n"
      << "    \"caveat\": \"Classical EAM results, not DFT. Not a QPU-stack dielectric calculation and not valid as DFT input without verification (see results/NiAl for an EAM failure case).\",\n";

    o << "    \"references\": {\n";
    const std::string elA = R.A, elB = R.B;
    for (int w = 0; w < 2; ++w) {
        const std::string& el = w ? elB : elA;
        const double mu = w ? R.muB : R.muA;
        const Candidate* ref = find_reference(R, w ? 1.0 : 0.0, mu);
        o << "      " << str(el) << ": {\"phase\": " << str(w ? R.refB : R.refA) << ", \"mu_eV_per_atom\": " << num(mu);
        if (ref) o << ", \"vol_per_atom_A3\": " << num(ref->vol_per_atom) << ", \"structure\": " << cell_json(ref->relaxed, "      ");
        o << "}" << (w ? "\n" : ",\n");
    }
    o << "    },\n";

    o << "    \"hull_phases\": [";
    for (std::size_t k = 0; k < R.hull.size(); ++k) {
        const Candidate* c = find_by_id(R, R.hull[k].id);
        o << (k ? ",\n" : "\n") << "      {\"id\": " << R.hull[k].id << ", \"x\": " << num(R.hull[k].x) << ", \"ef_eV_per_atom\": " << num(R.hull[k].ef);
        if (c) {
            o << ", \"label\": " << str(c->label) << ", \"source\": " << str(c->source)
              << ", \"e_per_atom_eV\": " << num(c->e_per_atom) << ", \"e_hull_eV_per_atom\": " << num(c->e_hull)
              << ", \"vol_per_atom_A3\": " << num(c->vol_per_atom) << ", \"relax_converged\": " << (c->converged ? "true" : "false")
              << ", \"structure\": " << cell_json(c->relaxed, "      ");
            if (c->screened) {
                o << ", \"elastic_GPa\": {\"C11\": " << num(c->elastic.C[0][0]) << ", \"C12\": " << num(c->elastic.C[0][1])
                  << ", \"C44\": " << num(c->elastic.C[3][3]) << ", \"bulk_voigt\": " << num(c->elastic.bulk_voigt)
                  << ", \"shear_voigt\": " << num(c->elastic.shear_voigt) << ", \"min_C_eig\": " << num(c->elastic.eigen[0])
                  << ", \"born_stable\": " << (c->elastic.mechanically_stable ? "true" : "false") << "}"
                  << ", \"verdict\": " << str(c->verdict);
            }
        }
        o << "}";
    }
    o << (R.hull.empty() ? "" : "\n    ") << "]\n  },\n";

    o << "  \"not_in_contract\": {\n"
      << "    \"Rn_spread_pct\": \"Process spread used by fab_monte_carlo(); not an atomistic output.\",\n"
      << "    \"environment\": \"T_eff_mK, flux_noise_uPhi0, charge_noise_e, Tphi_other_us stay in chip JSON.\"\n  },\n";
    o << "  \"open_stack\": [],\n";
    o << "  \"notes\": [\n"
      << "    \"Written by atomforge's forge.materials.v0 writer. Quantities listed in `placeholders` are NOT computed (value is null); copy nothing from them into a chip JSON.\",\n"
      << "    \"`atomforge_computed` is an additive optional block (the contract allows additive optional fields on v0); qpuforge does not read it.\"\n"
      << "  ]\n}\n";
    return o.str();
}

inline void write_materials_contract_json(const std::string& path, const SearchResult& R, const MaterialsOptions& opt = {}) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path + " for writing");
    f << materials_contract_json(R, opt);
    if (!f) throw std::runtime_error("write failed: " + path);
}

} // namespace atomforge
