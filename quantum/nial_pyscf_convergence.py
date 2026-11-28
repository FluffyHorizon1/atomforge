#!/usr/bin/env python3
"""Convergence study: formation energy of B2 NiAl with the pyscf backend of quantum/worker.py.

Periodic Kohn-Sham DFT, PBE, GTH pseudopotentials (Ni q18, Al q3), Gaussian basis, Fermi smearing. Three total energies
per setting, all at the experimental lattice constants:

    fcc Ni   a = 3.524 A, 4-atom cubic cell, spin-polarised (KUKS)
    fcc Al   a = 4.050 A, 4-atom cubic cell, spin-restricted
    B2 NiAl  a = 2.887 A, 2-atom cubic cell, spin-restricted (one setting is repeated spin-polarised)

    E_f = E(NiAl)/2 - E(Ni)/8 - E(Al)/8      eV per atom

The fcc elements are computed in the cubic cell, not the 1-atom primitive one, because PySCF's multigrid integrator is
far slower in a non-orthogonal cell (fcc Ni, primitive cell, 2x2x2 k-points, 45^3 grid: 485 s for one SCF). The same
k-mesh "k" is therefore a denser sampling for the larger fcc cells than for B2 NiAl; each energy has to converge on
its own.

The SCF is the one the bridge runs (worker.PySCF.scf with --multigrid); only the forces and the stress, which a
fixed-geometry energy does not need, are left out. quantum/worker.py refuses a spin-restricted run whose electron count
over the k-mesh is odd (Al has 3 valence electrons, NiAl 21), so the meshes here are even.

    python3 quantum/nial_pyscf_convergence.py --cache results/NiAl_pyscf_convergence.json [--threads 4] [--set NAME ...]

Every finished SCF is stored in the cache file at once, so an interrupted run continues where it stopped. The table is
printed from the cache; --table prints it without computing anything.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import worker  # noqa: E402

A_NI, A_AL, A_B2, A_L12 = 3.524, 4.050, 2.887, 3.570


def structures(scale=1.0):
    """The cells at the experimental lattice constants times `scale`."""
    def cubic4(a, els):   # conventional cubic cell of the fcc lattice: corner, then the three face centres
        return ([[a, 0, 0], [0, a, 0], [0, 0, a]], [[0, 0, 0], [0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]], els)
    b = A_B2 * scale
    return {
        "Ni": cubic4(A_NI * scale, ["Ni"] * 4),
        "Al": cubic4(A_AL * scale, ["Al"] * 4),
        "NiAl": ([[b, 0, 0], [0, b, 0], [0, 0, b]], [[0, 0, 0], [b / 2, b / 2, b / 2]], ["Ni", "Al"]),
        "Ni3Al": cubic4(A_L12 * scale, ["Al", "Ni", "Ni", "Ni"]),
    }


# formation energy per atom = E(compound cell) * c - E(Ni cell) * n - E(Al cell) * m  (cells as in structures())
FORMATION = {"NiAl": (1 / 2, 1 / 8, 1 / 8), "Ni3Al": (1 / 4, 3 / 16, 1 / 16)}
POLARISED = {"Ni": True, "Al": False, "NiAl": False, "Ni3Al": True}


# name -> list of (basis, k, sigma in Hartree, ke_cutoff in Hartree or None). The sets are sized for a 4-core machine.
SZV, DZVP, TZVP, TZV2P = "gth-szv-molopt-sr", "gth-dzvp-molopt-sr", "TZVP-MOLOPT-SR-GTH", "TZV2P-MOLOPT-SR-GTH"
CUT = 400.0     # Hartree; the cutoff table shows the formation energy converged to 0.1 meV/atom there
BASE = (DZVP, 4, 0.005, CUT)
SETS = {
    "1-cutoff": [(SZV, 2, 0.005, c) for c in (200.0, 400.0, 800.0, None)],
    "2-smearing": [(SZV, 4, s, CUT) for s in (0.01, 0.005, 0.0025)],
    "3-kpoints-szv": [(SZV, k, 0.005, CUT) for k in (2, 4, 6, 8, 10)],
    "4-kpoints-dzvp": [(DZVP, k, 0.005, CUT) for k in (2, 4, 6, 8, 10)],
    "5-kpoints-tzvp": [(TZVP, k, 0.005, CUT) for k in (4, 6, 8)],
    "6-basis": [(b, 6, 0.005, CUT) for b in (SZV, DZVP, TZVP, TZV2P)],
    "7-smearing-dzvp": [(DZVP, 8, s, CUT) for s in (0.005, 0.0025)],
}
# L1_2 Ni3Al (spin-polarised), same references
SETS_NI3AL = {
    "9-ni3al": [(DZVP, k, 0.005, CUT) for k in (6, 8, 10)] + [(TZVP, 6, 0.005, CUT)],
}
# energy-volume scans for the relaxed lattice constants: lattice constant = experiment * scale
EOS = (DZVP, 6, 0.005, CUT)
SCALES = (0.98, 0.99, 1.0, 1.01, 1.02)


def key(system, basis, k, sigma, cut, polarised, scale=1.0):
    kk = "%s|%s|k%d|s%g|c%s|%s" % (system, basis, k, sigma, "auto" if cut is None else "%g" % cut, "u" if polarised else "r")
    return kk if scale == 1.0 else kk + "|x%.3f" % scale


def run(system, basis, k, sigma, cut, polarised, threads, scale=1.0):
    args = ["pyscf", "--basis", basis, "--kpts", "%d,%d,%d" % (k, k, k), "--smearing", "%g" % sigma, "--multigrid",
            "--conv-tol", "1e-8", "--max-cycle", "150", "--threads", str(threads)]
    if cut is not None:
        args += ["--ke-cutoff", "%g" % cut]
    if polarised:
        args += ["--spin-polarised"]
    model = worker.PySCF(worker.build_parser().parse_args(args))
    model.kpoint_symmetry = True       # energies only: the SCF runs on the irreducible k-points
    cell, pos, sym = structures(scale)[system]
    t0 = time.time()
    _, info = model.scf(cell, pos, sym)
    info["wall"] = time.time() - t0
    info["name"] = model.name
    info["a"] = cell[0][0]
    return info


def need(system, basis, k, sigma, cut, polarised, cache, path, threads, compute, scale=1.0):
    kk = key(system, basis, k, sigma, cut, polarised, scale)
    if kk not in cache and compute:
        print("# running %s" % kk, file=sys.stderr, flush=True)
        try:
            cache[kk] = run(system, basis, k, sigma, cut, polarised, threads, scale)
        except Exception as e:  # an SCF that does not converge is a result too
            cache[kk] = {"error": str(e)}
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(cache, f, indent=1, sort_keys=True)
        os.replace(tmp, path)
    return cache.get(kk)


def table(name, points, cache, path, threads, compute, compound="NiAl"):
    ev = worker.HARTREE_EV
    wc, wn, wa = FORMATION[compound]
    print("\n## %s  (%s)" % (name, compound))
    print("%-19s %2s %7s %6s %4s %9s | %10s %10s | %7s %7s | %8s %8s %8s" % (
        "basis", "k", "sigma", "cut", "nao", "mesh", "Ef(F)", "Ef(s->0)", "m(Ni)", "m(cmp)", "t(Ni)", "t(Al)", "t(cmp)"))
    print("%-19s %2s %7s %6s %4s %9s | %10s %10s | %7s %7s | %8s %8s %8s" % (
        "", "", "Ha", "Ha", "cmp", "cmp", "eV/atom", "eV/atom", "muB/at", "muB/at", "s", "s", "s"))
    for basis, k, sigma, cut in points:
        r = {s: need(s, basis, k, sigma, cut, POLARISED[s], cache, path, threads, compute) for s in ("Ni", "Al", compound)}
        cuts = "auto" if cut is None else "%g" % cut
        if any(v is None for v in r.values()):
            print("%-19s %2d %7g %6s   (not computed)" % (basis, k, sigma, cuts))
            continue
        bad = [s for s, v in r.items() if "error" in v]
        if bad:
            print("%-19s %2d %7g %6s   failed: %s" % (basis, k, sigma, cuts, "; ".join("%s: %s" % (s, r[s]["error"]) for s in bad)))
            continue
        ef = (r[compound]["f"] * wc - r["Ni"]["f"] * wn - r["Al"]["f"] * wa) * ev
        ez = (r[compound]["e_zero"] * wc - r["Ni"]["e_zero"] * wn - r["Al"]["e_zero"] * wa) * ev
        nat = len(structures()[compound][2])
        print("%-19s %2d %7g %6s %4d %9s | %10.4f %10.4f | %7.3f %7.3f | %8.0f %8.0f %8.0f" % (
            basis, k, sigma, cuts, r[compound]["nao"], "%d^3" % r[compound]["mesh"][0], ef, ez, r["Ni"]["moment"] / 4,
            r[compound]["moment"] / nat, r["Ni"]["wall"], r["Al"]["wall"], r[compound]["wall"]))


def eos_fit(vols, energies):
    """Third-order Birch-Murnaghan fit (a cubic polynomial in V^(-2/3)). Returns V0, E0, B0 (eV/A^3) and whether V0
    lies inside the scanned range."""
    import numpy as np
    x = np.asarray(vols, float) ** (-2.0 / 3)
    c = np.polyfit(x, np.asarray(energies, float), 3)
    best = None
    for root in np.roots([3 * c[0], 2 * c[1], c[2]]):
        if abs(root.imag) > 1e-12 or root.real <= 0:
            continue
        x0 = root.real
        curv = 6 * c[0] * x0 + 2 * c[1]            # d2E/dx2; a minimum in V needs it positive
        if curv > 0 and (best is None or abs(x0 - x.mean()) < abs(best - x.mean())):
            best = x0
    if best is None:
        return None
    v0 = best ** -1.5
    e0 = float(np.polyval(c, best))
    b0 = v0 * (6 * c[0] * best + 2 * c[1]) * (2.0 / 3 * v0 ** (-5.0 / 3)) ** 2
    return v0, e0, b0, min(vols) <= v0 <= max(vols)


def eos(cache, path, threads, compute, compounds=("NiAl",)):
    ev = worker.HARTREE_EV
    basis, k, sigma, cut = EOS
    print("\n## 8-eos: energy-volume scans (%s, k=%d, sigma=%g Ha, cut %g Ha), lattice constant = experiment x %s"
          % (basis, k, sigma, cut, ", ".join("%g" % x for x in SCALES)))
    print("# Third-order Birch-Murnaghan fit of E(s->0) per cell. The cubic cells have no internal degrees of freedom,")
    print("# so the minimum is the relaxed structure.")
    print("%-6s %10s %10s %9s %13s" % ("", "a0 (A)", "a_exp (A)", "B0 (GPa)", "E0-E(a_exp)"))
    print("%-6s %10s %10s %9s %13s" % ("", "", "", "", "meV/atom"))
    fit = {}
    for system in ("Ni", "Al") + tuple(compounds):
        pts = [need(system, basis, k, sigma, cut, POLARISED[system], cache, path, threads, compute, x) for x in SCALES]
        if any(p is None or "error" in p for p in pts):
            print("%-6s (not computed)" % system)
            continue
        cell0 = structures()[system]
        nat = len(cell0[2])
        a_exp = cell0[0][0][0]
        vols = [(a_exp * x) ** 3 for x in SCALES]
        res = eos_fit(vols, [p["e_zero"] * ev for p in pts])
        if res is None:
            print("%-6s no minimum" % system)
            continue
        v0, e0, b0, inside = res
        fixed = pts[SCALES.index(1.0)]["e_zero"] * ev
        fit[system] = e0
        print("%-6s %10.4f %10.3f %9.1f %13.2f%s" % (system, v0 ** (1.0 / 3), a_exp, b0 * 160.21766, (e0 - fixed) / nat * 1000,
                                                     "" if inside else "   (outside the scanned range)"))
    for compound in compounds:
        if all(s in fit for s in ("Ni", "Al", compound)):
            wc, wn, wa = FORMATION[compound]
            r = {s: cache[key(s, basis, k, sigma, cut, POLARISED[s])]["e_zero"] * ev for s in ("Ni", "Al", compound)}
            ef_rel = fit[compound] * wc - fit["Ni"] * wn - fit["Al"] * wa
            ef_fix = r[compound] * wc - r["Ni"] * wn - r["Al"] * wa
            print("# %s: Ef(s->0) relaxed %.4f eV/atom, at the experimental lattice constants %.4f eV/atom; relaxation changes it by %+.4f eV/atom"
                  % (compound, ef_rel, ef_fix, ef_rel - ef_fix))


def main():
    allsets = sorted(SETS) + sorted(SETS_NI3AL) + ["8-eos", "8-eos-ni3al"]
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default="results/NiAl_pyscf_convergence.json")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--set", action="append", choices=allsets, help="which table(s) to compute (default: all)")
    ap.add_argument("--table", action="store_true", help="print the tables from the cache, compute nothing")
    o = ap.parse_args()
    cache = {}
    if os.path.exists(o.cache):
        with open(o.cache) as f:
            cache = json.load(f)
    import numpy
    import pyscf
    import platform
    print("# Formation energies of B2 NiAl (and L1_2 Ni3Al) with periodic DFT through quantum/worker.py pyscf:")
    print("# PBE, GTH-PBE pseudopotentials (Ni q18, Al q3), Gaussian basis sets, Fermi smearing, multigrid SCF, symmetry-reduced k-mesh.")
    print("# pyscf %s, numpy %s, Python %s, %d threads." % (pyscf.__version__, numpy.__version__, platform.python_version(), o.threads))
    print("# Cells: fcc Ni a = %.3f (4-atom cubic, spin-polarised), fcc Al a = %.3f (4-atom cubic), B2 NiAl a = %.3f (2 atoms),"
          % (A_NI, A_AL, A_B2))
    print("# L1_2 Ni3Al a = %.3f (4 atoms, spin-polarised); experimental lattice constants unless scanned."
          % A_L12)
    print("# Ef(F): from the Mermin free energies. Ef(s->0): from E0 = (E + F)/2, the usual zero-smearing estimate.")
    print("# m: magnetic moment per atom. t: wall time of one SCF, energy only.")
    for name in sorted(SETS):
        if o.set is None or name in o.set:
            table(name, SETS[name], cache, o.cache, o.threads, not o.table)
    # NiAl spin-polarised at the base setting: does a moment survive, does the energy change?
    basis, k, sigma, cut = BASE
    u = need("NiAl", basis, k, sigma, cut, True, cache, o.cache, o.threads, not o.table and (o.set is None or "4-kpoints-dzvp" in o.set))
    r = cache.get(key("NiAl", basis, k, sigma, cut, False))
    if u and r and "error" not in u and "error" not in r:
        print("\n## B2 NiAl from a spin-polarised start (%s, k=%d, sigma=%g): moment %.4f muB per cell, F(unrestricted) - F(restricted) = %.2e eV per cell"
              % (basis, k, sigma, u["moment"], (u["f"] - r["f"]) * worker.HARTREE_EV))
    if o.set is None or "8-eos" in o.set:
        eos(cache, o.cache, o.threads, not o.table)
    for name in sorted(SETS_NI3AL):
        if o.set is None or name in o.set:
            table(name, SETS_NI3AL[name], cache, o.cache, o.threads, not o.table, compound="Ni3Al")
    if o.set is None or "8-eos-ni3al" in o.set:
        eos(cache, o.cache, o.threads, not o.table, compounds=("Ni3Al",))


if __name__ == "__main__":
    main()
