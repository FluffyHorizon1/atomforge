// Tests for the forge.materials.v0 sidecar writer (include/atomforge/materials.hpp).
//   test_materials_contract <group>      groups: schema placeholders determinism roundtrip real_search
// Exit code 0 = pass. Each group prints PASS/FAIL lines; any failed check makes the group fail.
#include "atomforge/atomforge.hpp"
#include "json_mini.hpp"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>

using namespace atomforge;
using mini_json::Value;

static int g_fail = 0, g_checks = 0;
#define CHECK(cond) do { ++g_checks; if (!(cond)) { ++g_fail; std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)
#define CHECK_THROWS(expr) do { ++g_checks; bool t_ = false; try { (void)(expr); } catch (...) { t_ = true; } \
    if (!t_) { ++g_fail; std::printf("FAIL %s:%d: expected throw: %s\n", __FILE__, __LINE__, #expr); } } while (0)

static const std::string kTs = "2026-11-26T12:00:00+00:00";
static MaterialsOptions fixed_opts() {
    MaterialsOptions o; o.commit = "0123456789abcdef0123456789abcdef01234567"; o.timestamp = kTs; o.report_path = "results/CuAu/report.json";
    return o;
}

// Hand-built, fast, fully deterministic SearchResult (no relaxation). Values are arbitrary test data.
static SearchResult synthetic() {
    SearchResult R; R.A = "Cu"; R.B = "Au"; R.potential = "Zhou04-EAM";
    R.muA = -3.5; R.muB = -3.9; R.refA = "Cu A1 (fcc)"; R.refB = "Au A1 (fcc)";
    auto mk = [](int id, const char* src, const char* label, double x, double e, double ef, double eh, double a, bool scr) {
        Candidate c; c.id = id; c.source = src; c.label = label; c.parent = "fcc"; c.x = x; c.e_per_atom = e; c.ef = ef; c.e_hull = eh;
        c.converged = true; c.relaxed.cell = Mat3::diag(a, a, a); c.relaxed.label = label;
        for (int i = 0; i < 4; ++i) { c.relaxed.pos.push_back({0, 0, 0}); c.relaxed.species.push_back(0); }
        c.vol_per_atom = a * a * a / 4.0;
        if (scr) {
            c.screened = true; c.verdict = "STABLE (on hull, mechanically + thermally stable)";
            c.elastic.C[0][0] = 170.5; c.elastic.C[0][1] = 122.25; c.elastic.C[3][3] = 76.125;
            c.elastic.bulk_voigt = 138.3; c.elastic.shear_voigt = 51.5; c.elastic.eigen[0] = 48.25; c.elastic.mechanically_stable = true;
        }
        return c;
    };
    R.candidates.push_back(mk(0, "element", "A1 (fcc) Cu", 0.0, -3.5, 0.0, 0.0, 3.615, false));
    R.candidates.push_back(mk(1, "element", "A1 (fcc) Au", 1.0, -3.9, 0.0, 0.0, 4.078, false));
    R.candidates.push_back(mk(2, "prototype", "L1_0 CuAu", 0.5, -3.83, -0.13, 0.0, 3.9, true));
    R.hull = {{0.0, 0.0, 0}, {0.5, -0.13, 2}, {1.0, 0.0, 1}};
    return R;
}

static Value parsed(const SearchResult& R, const MaterialsOptions& o = fixed_opts()) { return mini_json::parse(materials_contract_json(R, o)); }

static void test_schema() {
    const Value v = parsed(synthetic());
    CHECK(v.type == Value::Object);
    // required top-level fields from docs/DATA_CONTRACT.md section 3.2
    for (const char* k : {"schema", "produced_by", "system", "structures_ref", "quantities", "not_in_contract", "open_stack", "notes"}) CHECK(v.has(k));
    CHECK(v.at("schema").type == Value::String && v.at("schema").s == "forge.materials.v0");
    // produced_by
    const Value& pb = v.at("produced_by");
    CHECK(pb.at("tool").s == "atomforge");
    CHECK(pb.at("commit").s == "0123456789abcdef0123456789abcdef01234567");
    CHECK(pb.at("potential_or_method").s == "Zhou04-EAM");
    CHECK(pb.at("timestamp").s == kTs);
    {   // ISO-8601 with explicit offset, e.g. 2026-11-26T12:00:00+00:00
        const std::string& t = pb.at("timestamp").s;
        CHECK(t.size() == 25 && t[4] == '-' && t[10] == 'T' && t[19] == '+' && t[22] == ':');
    }
    // system, structures_ref
    CHECK(v.at("system").type == Value::Array && v.at("system").a.size() == 2);
    CHECK(v.at("system").a[0].s == "Cu" && v.at("system").a[1].s == "Au");
    CHECK(v.path("structures_ref.report").s == "results/CuAu/report.json");
    CHECK(v.path("structures_ref.format_notes").s == "POSCAR in \xC3\x85, extxyz energy in eV (total), report energies in eV/atom");
    // quantities: the three contract entries, each with value/unit/maps_to
    const Value& q = v.at("quantities");
    for (const char* k : {"Q_dielectric", "loss_tangent", "junction_gap_ueV"}) {
        CHECK(q.has(k));
        CHECK(q.at(k).has("value")); CHECK(q.at(k).has("unit")); CHECK(q.at(k).has("maps_to"));
    }
    CHECK(q.path("Q_dielectric.unit").s == "dimensionless");
    CHECK(q.path("Q_dielectric.maps_to").s == "ChipDesign.materials.Q_dielectric");
    CHECK(q.path("Q_dielectric.conditions").has("T_K") && q.path("Q_dielectric.conditions").has("f_GHz"));
    CHECK(q.path("loss_tangent.unit").s == "dimensionless");
    CHECK(q.path("loss_tangent.maps_to").is_null());
    CHECK(q.path("loss_tangent.relation").s.find("Q_dielectric = 1 / loss_tangent") != std::string::npos);
    CHECK(q.path("junction_gap_ueV.unit").s == "ueV");
    CHECK(q.path("junction_gap_ueV.maps_to").s == "ChipDesign.materials.junction_gap_ueV");
    CHECK(v.path("not_in_contract.Rn_spread_pct").type == Value::String);
    CHECK(v.path("not_in_contract.environment").type == Value::String);
    CHECK(v.at("open_stack").type == Value::Array);
    CHECK(v.at("notes").type == Value::Array);
    // additive block: real atomforge quantities
    const Value& c = v.at("atomforge_computed");
    CHECK(c.at("status").s == "computed");
    CHECK(c.path("units.elastic").s == "GPa");
    CHECK(c.path("references.Cu.phase").s == "Cu A1 (fcc)");
    CHECK(std::abs(c.path("references.Cu.mu_eV_per_atom").n - (-3.5)) < 1e-12);
    CHECK(std::abs(c.path("references.Au.structure.cell_lengths_A").a[0].n - 4.078) < 1e-12);
    CHECK(std::abs(c.path("references.Cu.structure.cell_angles_deg").a[2].n - 90.0) < 1e-9);
    CHECK(c.at("hull_phases").a.size() == 3);
    const Value& ord = c.at("hull_phases").a[1];
    CHECK(ord.at("label").s == "L1_0 CuAu");
    CHECK(std::abs(ord.at("ef_eV_per_atom").n - (-0.13)) < 1e-12);
    CHECK(std::abs(ord.path("elastic_GPa.C11").n - 170.5) < 1e-12);
    CHECK(std::abs(ord.path("elastic_GPa.C12").n - 122.25) < 1e-12);
    CHECK(std::abs(ord.path("elastic_GPa.C44").n - 76.125) < 1e-12);
    CHECK(ord.path("elastic_GPa.born_stable").b);
    CHECK(!c.at("hull_phases").a[0].has("elastic_GPa"));   // unscreened phases carry no elastic block
}

static void test_placeholders() {
    const Value v = parsed(synthetic());
    CHECK(v.at("contains_placeholders").b);
    std::set<std::string> listed;
    for (const auto& p : v.at("placeholders").a) listed.insert(p.s);
    const std::set<std::string> want{"quantities.Q_dielectric", "quantities.loss_tangent", "quantities.junction_gap_ueV"};
    CHECK(listed == want);
    for (const auto& p : listed) {
        const Value& q = v.path(p);
        // not a number that could be mistaken for physics, and explicitly flagged
        CHECK(q.at("value").is_null());
        CHECK(q.at("status").s == "placeholder");
        CHECK(q.at("provenance").s == "placeholder");
        CHECK(!q.at("reason").s.empty());
        CHECK(q.has("qpuforge_default_for_reference"));
    }
    // every quantity in the contract block is a listed placeholder (no unflagged quantity can sneak in)
    for (const auto& kv : v.at("quantities").o) CHECK(listed.count("quantities." + kv.first) == 1);
    // the qpuforge defaults appear only under the clearly-named reference key, never as `value`
    CHECK(v.path("quantities.Q_dielectric.qpuforge_default_for_reference").n == 3000000);
    CHECK(v.path("quantities.junction_gap_ueV.qpuforge_default_for_reference").n == 180);
    CHECK(v.path("quantities.loss_tangent.qpuforge_default_for_reference").is_null());
    // the raw text must not contain a bare "value": <number> for any placeholder
    const std::string txt = materials_contract_json(synthetic(), fixed_opts());
    CHECK(txt.find("\"value\": 0") == std::string::npos);
    CHECK(txt.find("\"value\": 3") == std::string::npos);
    CHECK(txt.find("\"value\": 1") == std::string::npos);
    // computed quantities are not flagged as placeholders
    CHECK(v.path("atomforge_computed.status").s == "computed");
    for (const auto& p : listed) CHECK(p.rfind("atomforge_computed", 0) != 0);
}

static void test_determinism() {
    const SearchResult R = synthetic();
    const std::string a = materials_contract_json(R, fixed_opts()), b = materials_contract_json(R, fixed_opts());
    CHECK(a == b);
    CHECK(!a.empty() && a.back() == '\n');
    // file writer produces byte-identical output to the string builder
    const auto path = std::filesystem::temp_directory_path() / "atomforge_materials_determinism.json";
    write_materials_contract_json(path.string(), R, fixed_opts());
    std::ifstream f(path, std::ios::binary); std::stringstream ss; ss << f.rdbuf();
    CHECK(ss.str() == a);
    std::filesystem::remove(path);
    // different timestamp / commit -> only those fields change
    MaterialsOptions o2 = fixed_opts(); o2.timestamp = "2026-11-27T00:00:00+00:00";
    CHECK(materials_contract_json(R, o2) != a);
    // default timestamp is generated in ISO-8601 UTC with offset
    MaterialsOptions o3 = fixed_opts(); o3.timestamp.clear();
    const Value v3 = parsed(R, o3);
    const std::string t = v3.path("produced_by.timestamp").s;
    CHECK(t.size() == 25 && t.substr(19) == "+00:00");
    // unwritable path throws rather than silently succeeding
    CHECK_THROWS(write_materials_contract_json("/nonexistent-dir-atomforge/x.json", R, fixed_opts()));
}

static void test_roundtrip() {
    // numbers survive write -> parse within 1e-9 relative; strings survive escaping; non-finite -> null
    SearchResult R = synthetic();
    R.candidates[2].label = "L1_0 \"quoted\" \\ back\tslash CuAu";
    R.candidates[2].e_per_atom = -3.8312345678912;
    R.candidates[2].ef = -0.1298765432109;
    R.candidates[2].elastic.C[0][0] = 170.123456789012;
    R.candidates[2].vol_per_atom = std::nan("");   // must serialise as null, not "nan"
    R.hull[1].ef = -0.1298765432109;
    const std::string txt = materials_contract_json(R, fixed_opts());
    for (const char* bad : {": nan", ": -nan", ": inf", ": -inf", ", nan", ", inf"}) CHECK(txt.find(bad) == std::string::npos);
    const Value v = mini_json::parse(txt);
    const Value& ph = v.path("atomforge_computed").at("hull_phases").a[1];
    CHECK(ph.at("label").s == R.candidates[2].label);
    CHECK(std::abs(ph.at("e_per_atom_eV").n - R.candidates[2].e_per_atom) < 1e-9 * 3.9);
    CHECK(std::abs(ph.at("ef_eV_per_atom").n - R.hull[1].ef) < 1e-9);
    CHECK(std::abs(ph.path("elastic_GPa.C11").n - 170.123456789012) < 1e-7);
    CHECK(ph.at("vol_per_atom_A3").is_null());
    // re-serialising the parsed document's key set is stable: parse(write(R)) has the same key order as parse(write(R))
    auto keys = [](const Value& x) { std::vector<std::string> k; for (auto& kv : x.o) k.push_back(kv.first); return k; };
    CHECK(keys(v) == keys(parsed(R)));
    const std::vector<std::string> want{"schema", "produced_by", "system", "structures_ref", "quantities", "placeholders",
                                        "contains_placeholders", "atomforge_computed", "not_in_contract", "open_stack", "notes"};
    CHECK(keys(v) == want);
    // empty hull still yields valid JSON
    SearchResult E; E.A = "Cu"; E.B = "Au"; E.potential = "Zhou04-EAM";
    const Value ve = parsed(E);
    CHECK(ve.path("atomforge_computed.hull_phases").a.empty());
    CHECK(ve.at("schema").s == "forge.materials.v0");
    // the parser itself rejects malformed input (guards the guard)
    CHECK_THROWS(mini_json::parse("{\"a\": }"));
    CHECK_THROWS(mini_json::parse(txt.substr(0, txt.size() / 2)));
}

// Real (tiny) atomforge run: no Monte Carlo, 1 solid-solution seed, short MD. Exercises relax + Born + hull end to end.
static void test_real_search() {
    SearchOptions so; so.verbose = false; so.mc_parents.clear(); so.ss_seeds = 1; so.md_steps = 200; so.md_T = 300;
    const SearchResult R = binary_search("Cu", "Au", so);
    const MaterialsOptions mo = fixed_opts();
    const std::string a = materials_contract_json(R, mo);
    const Value v = mini_json::parse(a);
    CHECK(v.at("schema").s == "forge.materials.v0");
    CHECK(v.path("produced_by.potential_or_method").s == R.potential);
    const Value& comp = v.at("atomforge_computed");
    // elemental references reproduce the search's own chemical potentials
    CHECK(std::abs(comp.path("references.Cu.mu_eV_per_atom").n - R.muA) < 1e-8);
    CHECK(std::abs(comp.path("references.Au.mu_eV_per_atom").n - R.muB) < 1e-8);
    // Cu fcc lattice constant sanity: cell edge of the 4-atom conventional cell (expt 3.615 A)
    const double aCu = comp.path("references.Cu.structure.cell_lengths_A").a[0].n;
    std::printf("Cu reference: %s, cell edge %.4f A\n", comp.path("references.Cu.phase").s.c_str(), aCu);
    CHECK(aCu > 3.4 && aCu < 3.8);
    CHECK(comp.path("references.Cu.structure.natoms").n == 4);
    // hull phases mirror the SearchResult; screened ones carry elastic constants equal to the computed tensor
    CHECK(comp.at("hull_phases").a.size() == R.hull.size());
    int screened = 0;
    for (std::size_t k = 0; k < R.hull.size(); ++k) {
        const Value& hp = comp.at("hull_phases").a[k];
        CHECK(int(hp.at("id").n) == R.hull[k].id);
        CHECK(std::abs(hp.at("ef_eV_per_atom").n - R.hull[k].ef) < 1e-8);
        const Candidate* c = materials_detail::find_by_id(R, R.hull[k].id);
        CHECK(c != nullptr);
        if (c && c->screened) {
            ++screened;
            CHECK(std::abs(hp.path("elastic_GPa.C11").n - c->elastic.C[0][0]) < 1e-6 * std::max(1.0, std::abs(c->elastic.C[0][0])));
            CHECK(std::abs(hp.path("elastic_GPa.C44").n - c->elastic.C[3][3]) < 1e-6 * std::max(1.0, std::abs(c->elastic.C[3][3])));
            CHECK(hp.at("verdict").s == c->verdict);
        } else CHECK(!hp.has("elastic_GPa"));
    }
    std::printf("real search: %zu candidates, %zu hull phases, %d with elastic constants\n", R.candidates.size(), R.hull.size(), screened);
    CHECK(screened > 0);
    // placeholders still flagged on a real run
    CHECK(v.at("quantities").path("Q_dielectric.value").is_null());
    CHECK(v.at("quantities").path("Q_dielectric.status").s == "placeholder");
    // determinism: the same search re-run gives byte-identical sidecar text
    const SearchResult R2 = binary_search("Cu", "Au", so);
    CHECK(materials_contract_json(R2, mo) == a);
}

int main(int argc, char** argv) {
    const std::vector<std::pair<std::string, std::function<void()>>> groups{
        {"schema", test_schema}, {"placeholders", test_placeholders}, {"determinism", test_determinism},
        {"roundtrip", test_roundtrip}, {"real_search", test_real_search}};
    if (argc < 2) { std::fprintf(stderr, "usage: test_materials_contract <group|all>\n"); return 2; }
    const std::string want = argv[1];
    bool ran = false;
    for (const auto& [name, fn] : groups) {
        if (want != "all" && want != name) continue;
        ran = true;
        const int before = g_fail;
        try { fn(); } catch (const std::exception& e) { ++g_fail; std::printf("FAIL %s: exception: %s\n", name.c_str(), e.what()); }
        std::printf("%s %s\n", g_fail == before ? "PASS" : "FAIL", name.c_str());
    }
    if (!ran) { std::fprintf(stderr, "unknown group %s\n", want.c_str()); return 2; }
    std::printf("%d checks, %d failed\n", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
