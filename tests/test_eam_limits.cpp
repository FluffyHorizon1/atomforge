// EAM limits: tests that pin what the Zhou-Johnson-Wadley EAM (the only built-in potential for search) cannot do.
//
//   test_eam_limits eam_nial    Ni-Al is a KNOWN FAILURE (README "Validation"): the ordered prototypes come out at zero or positive
//                               formation energy, where experiment gives about -0.6 eV/atom for B2 NiAl. This group asserts the
//                               current known-bad numbers, so it passes while the failure is unchanged and FAILS (with a message
//                               saying what to update) if the failure changes in either direction: a fixed potential, or a
//                               regression that moves the numbers. It is not a statement that the numbers are right.
//   test_eam_limits eam_scope   Which elements the EAM can treat at all: 17 metals. O, N, Si and Nb (the QPU stack's oxides, nitride,
//                               substrate and a second superconductor) have no parameters; asking for them is an error, not a guess.
//
// The pins are the values recorded in results/NiAl/report.json (potential Zhou04-EAM); eam_nial also checks that file against a fresh
// computation, so the committed record cannot drift from the code unnoticed.
#include "atomforge/atomforge.hpp"
#include "json_mini.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifndef ATOMFORGE_SOURCE_DIR
#define ATOMFORGE_SOURCE_DIR "."
#endif

using namespace atomforge;
static int g_fail = 0, g_pass = 0;
#define CHECK(cond, ...) do { if (cond) { ++g_pass; } else { ++g_fail; std::printf("  FAIL %s:%d  ", __FILE__, __LINE__); std::printf(__VA_ARGS__); std::printf("\n"); } } while (0)

// What to do if a pin below fails. Printed with every failure of the group.
static const char* kUpdate =
    "the EAM Ni-Al known failure has CHANGED. If that is intended (new potential, fixed mixing), update together: the pins in "
    "tests/test_eam_limits.cpp, README.md 'Ni-Al: fails', docs/EAM_LIMITATIONS.md, and results/NiAl/ (regenerate with the PROVENANCE command). "
    "If it is not intended, a change to the potential, the relaxation or the prototypes moved a number that was meant to stay put.";

// Experiment, README "Validation": B2 NiAl forms with about -0.6 eV/atom.
static constexpr double kExperimentB2 = -0.6;
// Pins, eV/atom. kEfB2 and the elemental energies are the values in results/NiAl/report.json (Zhou04-EAM, default search options).
// kEhullB2 is relative to the hull of the ordered prototypes only (report.json's 0.0746 includes the Monte Carlo point at NiAl3,
// which this fast test does not run). kEfMinProto is L1_2 NiAl3, the most strongly bound prototype: -9 meV/atom, "near zero".
static constexpr double kEfB2 = 0.0676393, kEhullB2 = 0.0738018, kMuNi = -4.45, kMuAl = -3.58, kEfMinProto = -0.00924;
// Ordering failure (equiatomic ordered phases): the EAM puts B2 (the observed NiAl structure) 55 meV ABOVE L1_1 and 51 meV above L1_0.
static constexpr double kEfL10 = 0.0162, kEfL11 = 0.0121;
static constexpr double kBoundFloor = -0.05;   // no Ni-Al prototype or hull vertex may bind by more than 50 meV/atom (experiment: ~600)
static constexpr double kTolPin = 2e-3;       // eV/atom: tight enough that any change to the potential or relaxation shows
static constexpr double kTolRecord = 1e-6;    // fresh computation vs the committed report.json: the same code, so equal

static SearchResult nial_prototypes() {
    SearchOptions o;
    o.verbose = false;
    o.ss_parents.clear();         // no random solid solutions
    o.mc_parents.clear();         // no Monte Carlo ordering: the ordered prototypes decide the verdict and are deterministic
    o.screen_window = -1;         // no stability screen: its MD is not what is pinned here
    return binary_search("Ni", "Al", o);
}

static void test_eam_nial() {
    std::puts("[eam_nial] Ni-Al with the Zhou04 EAM: the documented known failure, pinned");
    const SearchResult R = nial_prototypes();
    const Candidate* b2 = nullptr; const Candidate* l10 = nullptr; const Candidate* l11 = nullptr;
    double min_proto_ef = 1e9; std::string min_proto;
    int n_proto = 0;
    for (const auto& c : R.candidates) {
        if (c.source != "prototype") continue;
        ++n_proto;
        if (c.label.rfind("B2 ", 0) == 0) b2 = &c;
        if (c.label.rfind("L1_0 NiAl", 0) == 0) l10 = &c;
        if (c.label.rfind("L1_1 NiAl", 0) == 0) l11 = &c;
        if (c.ef < min_proto_ef) { min_proto_ef = c.ef; min_proto = c.label; }
    }
    std::printf("  references: %s (%.4f eV), %s (%.4f eV); %d ordered prototypes, lowest Ef %+.4f eV/atom (%s)\n",
                R.refA.c_str(), R.muA, R.refB.c_str(), R.muB, n_proto, min_proto_ef, min_proto.c_str());
    CHECK(R.refA == "Ni A1 (fcc)" && R.refB == "Al A1 (fcc)", "references are %s / %s, expected Ni A1 (fcc) / Al A1 (fcc); %s", R.refA.c_str(), R.refB.c_str(), kUpdate);
    CHECK(std::abs(R.muA - kMuNi) < kTolPin && std::abs(R.muB - kMuAl) < kTolPin, "elemental energies %.4f / %.4f eV, pinned %.4f / %.4f; %s", R.muA, R.muB, kMuNi, kMuAl, kUpdate);
    CHECK(b2 != nullptr, "no B2 NiAl prototype among the candidates");
    if (!b2) return;
    std::printf("  B2 NiAl: Ef %+.5f eV/atom, Ehull %.5f eV/atom, %.3f A^3/atom (experiment: about %.1f eV/atom)\n", b2->ef, b2->e_hull, b2->vol_per_atom, kExperimentB2);
    CHECK(b2->converged, "B2 NiAl relaxation did not converge");
    CHECK(std::abs(b2->x - 0.5) < 1e-12 && b2->relaxed.size() == 2, "B2 NiAl should be the 2-atom 50:50 cell");

    // 1. The known-bad numbers, pinned.
    CHECK(std::abs(b2->ef - kEfB2) < kTolPin, "B2 NiAl Ef = %+.5f eV/atom, pinned %+.5f +- %.0e; %s", b2->ef, kEfB2, kTolPin, kUpdate);
    CHECK(std::abs(b2->e_hull - kEhullB2) < kTolPin, "B2 NiAl Ehull = %.5f eV/atom, pinned %.5f +- %.0e; %s", b2->e_hull, kEhullB2, kTolPin, kUpdate);

    // 2. The failure itself, stated as inequalities so that the message says which direction it moved.
    CHECK(b2->ef > 0, "B2 NiAl Ef = %+.5f: no longer positive. The EAM now binds Ni-Al; %s", b2->ef, kUpdate);
    CHECK(b2->ef - kExperimentB2 > 0.5, "B2 NiAl Ef = %+.5f eV/atom is within 0.5 eV of experiment (%.1f); the failure shrank; %s", b2->ef, kExperimentB2, kUpdate);
    // Ordering failure: of the equiatomic ordered phases the EAM ranks B2 above L1_1 and L1_0 (the wrong way round for NiAl).
    CHECK(l10 && l11, "L1_0 / L1_1 NiAl prototypes missing");
    if (l10 && l11) {
        std::printf("  equiatomic Ef: L1_1 %+.4f, L1_0 %+.4f, B2 %+.4f eV/atom\n", l11->ef, l10->ef, b2->ef);
        CHECK(std::abs(l10->ef - kEfL10) < kTolPin && std::abs(l11->ef - kEfL11) < kTolPin, "L1_0 / L1_1 NiAl Ef = %+.5f / %+.5f, pinned %+.4f / %+.4f; %s", l10->ef, l11->ef, kEfL10, kEfL11, kUpdate);
        CHECK(b2->ef > l10->ef + 0.03 && b2->ef > l11->ef + 0.03, "B2 NiAl is no longer 30 meV or more above L1_0 and L1_1 (the ordering failure changed); %s", kUpdate);
    }
    // The most strongly bound prototype is still only -9 meV/atom ("near zero"): nothing comes close to a real intermetallic.
    CHECK(std::abs(min_proto_ef - kEfMinProto) < kTolPin && min_proto == "L1_2 NiAl3", "lowest prototype Ef = %+.5f (%s), pinned %+.5f (L1_2 NiAl3); %s", min_proto_ef, min_proto.c_str(), kEfMinProto, kUpdate);
    CHECK(min_proto_ef > kBoundFloor, "an ordered Ni-Al prototype now binds by %.3f eV/atom (%s), more than %.2f; %s", -min_proto_ef, min_proto.c_str(), -kBoundFloor, kUpdate);
    // Prototype-only hull: the two elements plus that one weakly bound NiAl3 point; no strongly bound compound.
    double hull_min = 0; for (const auto& h : R.hull) hull_min = std::min(hull_min, h.ef);
    CHECK(R.hull.size() == 3 && hull_min > kBoundFloor, "prototype hull has %zu vertices (pinned 3: Ni, L1_2 NiAl3, Al), deepest Ef %+.5f; %s", R.hull.size(), hull_min, kUpdate);

    // 3. The committed record (results/NiAl/report.json) says the same as the code.
    const std::string path = std::string(ATOMFORGE_SOURCE_DIR) + "/results/NiAl/report.json";
    std::ifstream f(path);
    CHECK(bool(f), "cannot open %s", path.c_str());
    if (!f) return;
    std::stringstream ss; ss << f.rdbuf();
    const mini_json::Value J = mini_json::parse(ss.str());
    CHECK(J.at("potential").s == "Zhou04-EAM", "results/NiAl/report.json potential is %s", J.at("potential").s.c_str());
    bool found = false;
    for (const auto& c : J.at("candidates").a) {
        if (c.at("label").s != "B2 NiAl") continue;
        found = true;
        CHECK(std::abs(c.at("ef").n - b2->ef) < kTolRecord, "results/NiAl/report.json B2 NiAl ef %.7f vs computed %.7f", c.at("ef").n, b2->ef);
    }
    CHECK(found, "results/NiAl/report.json has no B2 NiAl candidate");
}

static void test_eam_scope() {
    std::puts("[eam_scope] elements the Zhou04 EAM can and cannot treat (QPU materials stack)");
    CHECK(zhou04::table.size() == 17, "the Zhou04 table has %zu elements, documented 17 (docs/EAM_LIMITATIONS.md lists them)", zhou04::table.size());
    // Present: metals only. Al, Ta and Ti are elements of the QPU stack that the EAM can describe as bulk metals.
    for (const char* el : {"Al", "Ta", "Ti"}) CHECK(ZhouEAM::lookup(el) != nullptr, "%s should be in the Zhou04 table", el);
    // Absent: the oxygen, nitrogen, silicon and niobium the stack needs (AlOx barriers, Nb/Ta oxides, TiN, Si substrate).
    for (const char* el : {"O", "N", "Si", "Nb"}) {
        CHECK(ZhouEAM::lookup(el) == nullptr, "%s is now in the Zhou04 table: update docs/EAM_LIMITATIONS.md and this test", el);
        bool threw = false; std::string msg;
        try { ZhouEAM p({"Al", el}); } catch (const std::runtime_error& e) { threw = true; msg = e.what(); }
        CHECK(threw && msg.find(std::string("no parameters for ") + el) != std::string::npos, "ZhouEAM({Al, %s}) should refuse with 'no parameters for %s' (got: %s)", el, el, threw ? msg.c_str() : "no error");
        threw = false; msg.clear();
        try { binary_search("Al", el, SearchOptions{}); } catch (const std::runtime_error& e) { threw = true; msg = e.what(); }
        CHECK(threw, "binary_search(Al, %s) should fail, not produce numbers", el);
    }
    // The compounds the stack is made of cannot be set up: each needs at least one absent element.
    const char* compounds[][3] = {{"AlOx barrier", "Al", "O"}, {"Nb2O5 / Ta2O5 surface oxide", "Ta", "O"}, {"TiN", "Ti", "N"}, {"Si substrate", "Si", "Si"}};
    for (const auto& c : compounds) {
        bool treatable = ZhouEAM::lookup(c[1]) && ZhouEAM::lookup(c[2]);
        std::printf("  %-30s %s + %s: %s\n", c[0], c[1], c[2], treatable ? "TREATABLE" : "not treatable (no EAM parameters)");
        CHECK(!treatable, "%s became treatable; the README 'Toward the quantum level' paragraph and docs/EAM_LIMITATIONS.md need updating", c[0]);
    }
}

int main(int argc, char** argv) {
    const char* g = argc > 1 ? argv[1] : "all";
    const bool all = !std::strcmp(g, "all");
    bool ran = false;
    if (all || !std::strcmp(g, "eam_nial")) { test_eam_nial(); ran = true; }
    if (all || !std::strcmp(g, "eam_scope")) { test_eam_scope(); ran = true; }
    if (!ran) { std::fprintf(stderr, "usage: test_eam_limits [eam_nial|eam_scope|all]\n"); return 2; }
    std::printf("\n%d checks passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
