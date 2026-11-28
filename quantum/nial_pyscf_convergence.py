#!/usr/bin/env python3
"""Convergence study: formation energy of B2 NiAl at fixed geometry with the pyscf backend of quantum/worker.py.

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

A_NI, A_AL, A_B2 = 3.524, 4.050, 2.887


def structures():
    def fcc(a, el):   # conventional cubic cell, 4 atoms
        return ([[a, 0, 0], [0, a, 0], [0, 0, a]], [[0, 0, 0], [0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]], [el] * 4)
    b = A_B2
    return {
        "Ni": fcc(A_NI, "Ni"),
        "Al": fcc(A_AL, "Al"),
        "NiAl": ([[b, 0, 0], [0, b, 0], [0, 0, b]], [[0, 0, 0], [b / 2, b / 2, b / 2]], ["Ni", "Al"]),
    }


# name -> list of (basis, k, sigma in Hartree, ke_cutoff in Hartree or None). The sets are sized for a 4-core machine.
SZV, DZVP = "gth-szv-molopt-sr", "gth-dzvp-molopt-sr"
BASE = (SZV, 4, 0.005, None)
SETS = {
    "1-kpoints-szv": [(SZV, k, 0.005, None) for k in (2, 4, 6)],
    "2-kpoints-dzvp": [(DZVP, k, 0.005, None) for k in (2, 4)],
    "3-basis": [(b, k, 0.005, None) for k in (2, 4) for b in (SZV, DZVP)],
    "4-smearing": [(SZV, 4, s, None) for s in (0.01, 0.005)],
    "5-cutoff": [(SZV, 2, 0.005, c) for c in (200.0, 400.0, 800.0, None)],
}


def key(system, basis, k, sigma, cut, polarised):
    return "%s|%s|k%d|s%g|c%s|%s" % (system, basis, k, sigma, "auto" if cut is None else "%g" % cut, "u" if polarised else "r")


def run(system, basis, k, sigma, cut, polarised, threads):
    args = ["pyscf", "--basis", basis, "--kpts", "%d,%d,%d" % (k, k, k), "--smearing", "%g" % sigma, "--multigrid",
            "--conv-tol", "1e-8", "--max-cycle", "150", "--threads", str(threads)]
    if cut is not None:
        args += ["--ke-cutoff", "%g" % cut]
    if polarised:
        args += ["--spin-polarised"]
    model = worker.PySCF(worker.build_parser().parse_args(args))
    model.kpoint_symmetry = True       # energies only: the SCF runs on the irreducible k-points
    cell, pos, sym = structures()[system]
    t0 = time.time()
    _, info = model.scf(cell, pos, sym)
    info["wall"] = time.time() - t0
    info["name"] = model.name
    return info


def need(system, basis, k, sigma, cut, polarised, cache, path, threads, compute):
    kk = key(system, basis, k, sigma, cut, polarised)
    if kk not in cache and compute:
        print("# running %s" % kk, file=sys.stderr, flush=True)
        try:
            cache[kk] = run(system, basis, k, sigma, cut, polarised, threads)
        except Exception as e:  # an SCF that does not converge is a result too
            cache[kk] = {"error": str(e)}
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(cache, f, indent=1, sort_keys=True)
        os.replace(tmp, path)
    return cache.get(kk)


def table(name, points, cache, path, threads, compute):
    ev = worker.HARTREE_EV
    print("\n## %s" % name)
    print("%-20s %2s %7s %6s %4s %9s | %10s %10s | %7s | %8s %8s %8s" % (
        "basis", "k", "sigma", "cut", "nao", "mesh", "Ef(F)", "Ef(s->0)", "m(Ni)", "t(Ni)", "t(Al)", "t(NiAl)"))
    print("%-20s %2s %7s %6s %4s %9s | %10s %10s | %7s | %8s %8s %8s" % (
        "", "", "Ha", "Ha", "NiAl", "NiAl", "eV/atom", "eV/atom", "muB", "s", "s", "s"))
    for basis, k, sigma, cut in points:
        r = {s: need(s, basis, k, sigma, cut, s == "Ni", cache, path, threads, compute) for s in ("Ni", "Al", "NiAl")}
        cuts = "auto" if cut is None else "%g" % cut
        if any(v is None for v in r.values()):
            print("%-20s %2d %7g %6s   (not computed)" % (basis, k, sigma, cuts))
            continue
        bad = [s for s, v in r.items() if "error" in v]
        if bad:
            print("%-20s %2d %7g %6s   failed: %s" % (basis, k, sigma, cuts, "; ".join("%s: %s" % (s, r[s]["error"]) for s in bad)))
            continue
        ef = (r["NiAl"]["f"] / 2 - r["Ni"]["f"] / 8 - r["Al"]["f"] / 8) * ev
        ez = (r["NiAl"]["e_zero"] / 2 - r["Ni"]["e_zero"] / 8 - r["Al"]["e_zero"] / 8) * ev
        print("%-20s %2d %7g %6s %4d %9s | %10.4f %10.4f | %7.3f | %8.0f %8.0f %8.0f" % (
            basis, k, sigma, cuts, r["NiAl"]["nao"], "%d^3" % r["NiAl"]["mesh"][0], ef, ez, r["Ni"]["moment"] / 4,
            r["Ni"]["wall"], r["Al"]["wall"], r["NiAl"]["wall"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default="results/NiAl_pyscf_convergence.json")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--set", action="append", choices=sorted(SETS), help="which table(s) to compute (default: all)")
    ap.add_argument("--table", action="store_true", help="print the tables from the cache, compute nothing")
    o = ap.parse_args()
    cache = {}
    if os.path.exists(o.cache):
        with open(o.cache) as f:
            cache = json.load(f)
    import pyscf
    print("# B2 NiAl formation energy at fixed geometry: PBE, GTH-PBE (Ni q18, Al q3), Fermi smearing, pyscf %s" % pyscf.__version__)
    print("# a(Ni) = %.3f, a(Al) = %.3f, a(NiAl) = %.3f A. Ef(F): from the free energies; Ef(s->0): extrapolated to zero smearing." % (A_NI, A_AL, A_B2))
    print("# t: wall time of one SCF, %d threads. m(Ni): moment of fcc Ni per atom." % o.threads)
    for name in (o.set or sorted(SETS)):
        table(name, SETS[name], cache, o.cache, o.threads, not o.table)
    # NiAl spin-polarised at the base setting: does a moment survive, does the energy change?
    basis, k, sigma, cut = BASE
    u = need("NiAl", basis, k, sigma, cut, True, cache, o.cache, o.threads, not o.table and not o.set)
    r = cache.get(key("NiAl", basis, k, sigma, cut, False))
    if u and r and "error" not in u and "error" not in r:
        print("\n## B2 NiAl spin-polarised start (%s, k=%d, sigma=%g): moment %.4f muB per cell, F(unrestricted) - F(restricted) = %.2e eV"
              % (basis, k, sigma, u["moment"], (u["f"] - r["f"]) * worker.HARTREE_EV))


if __name__ == "__main__":
    main()
