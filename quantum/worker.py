#!/usr/bin/env python3
"""atomforge external-potential worker (protocol 1).

Reads one JSON request per line on stdin and writes one JSON reply per line on stdout. Diagnostics go to stderr, which
ExternalPotential appends to its log file. The protocol is documented in docs/EXTERNAL_PROTOCOL.md.

    python3 quantum/worker.py lj    [--epsilon EV] [--sigma A] [--rc A] [--mass AMU]
    python3 quantum/worker.py pyscf [--basis gth-szv] [--pseudo gth-pbe] [--xc pbe] [--kpts 2,2,2] [--smearing HA]
    python3 quantum/worker.py mace  [--model medium] [--device cpu]
    python3 quantum/worker.py --check BACKEND     exit 0 if BACKEND can run here, 77 if its packages are missing

Units on the wire: A, eV, eV/A. "virial" is W_ab = dE/d(eps_ab) in eV (tensile positive); "stress" is W/V in eV/A^3.

The lj backend is pure Python (no numpy) so the C++ tests can run it anywhere python3 exists. Its --fault options exist
only to exercise the error paths of the bridge.
"""
import argparse
import json
import math
import sys

PROTOCOL = 1
SKIP = 77  # CTest SKIP_RETURN_CODE

HARTREE_EV = 27.211386245988
BOHR_A = 0.529177210903


def log(*a):
    print("[worker]", *a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- small 3x3 helpers (pure Python)
def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def inv3(m):
    d = det3(m)
    if abs(d) < 1e-14:
        raise ValueError("singular cell")
    return [[(m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d, (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d,
             (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d],
            [(m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d, (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d,
             (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d],
            [(m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d, (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d,
             (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d]]


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


# ---------------------------------------------------------------- backends
class LJ:
    """Truncated-and-shifted Lennard-Jones; the same formula as atomforge::LennardJones (include/atomforge/lj.hpp)."""

    def __init__(self, o):
        self.eps, self.sig, self.rc, self.m = o.epsilon, o.sigma, o.rc, o.mass
        if not (self.eps > 0 and self.sig > 0 and self.rc > self.sig):
            raise ValueError("lj: need epsilon > 0, sigma > 0, rc > sigma")
        sr6 = math.pow(self.sig / self.rc, 6)
        self.shift = 4 * self.eps * (sr6 * sr6 - sr6)
        self.name = "lj(eps=%g,sigma=%g,rc=%g)" % (self.eps, self.sig, self.rc)
        self.cutoff = self.rc

    def masses(self, elements):
        return [self.m] * len(elements)

    def evaluate(self, cell, pos, species, symbols):
        n = len(pos)
        rc2 = self.rc * self.rc
        a, b, c = cell
        vol = abs(dot(a, cross(b, c)))
        widths = [vol / math.sqrt(dot(x, x)) for x in (cross(b, c), cross(c, a), cross(a, b))]
        # fractional coordinates: f = r . inv(cell)  (cell rows are lattice vectors)
        inv = inv3(cell)
        frac = [[sum(r[k] * inv[k][j] for k in range(3)) for j in range(3)] for r in pos]
        nimg = []
        for k in range(3):
            fk = [f[k] for f in frac] or [0.0]
            spread = max(fk) - min(fk)
            x = self.rc / widths[k]
            nimg.append(max(int(math.ceil(x)), int(math.floor(x + spread))))
        e = 0.0
        F = [[0.0, 0.0, 0.0] for _ in range(n)]
        W = [[0.0] * 3 for _ in range(3)]
        for ia in range(-nimg[0], nimg[0] + 1):
            for ib in range(-nimg[1], nimg[1] + 1):
                for ic in range(-nimg[2], nimg[2] + 1):
                    sh = [ia * a[k] + ib * b[k] + ic * c[k] for k in range(3)]
                    zero = ia == 0 and ib == 0 and ic == 0
                    positive = ia > 0 or (ia == 0 and (ib > 0 or (ib == 0 and ic > 0)))
                    for i in range(n):
                        pi = pos[i]
                        for j in range(i + 1 if zero else i, n):
                            if i == j and not positive:
                                continue
                            pj = pos[j]
                            d0 = pj[0] + sh[0] - pi[0]
                            d1 = pj[1] + sh[1] - pi[1]
                            d2 = pj[2] + sh[2] - pi[2]
                            r2 = d0 * d0 + d1 * d1 + d2 * d2
                            if r2 >= rc2:
                                continue
                            sr2 = self.sig * self.sig / r2
                            sr6 = sr2 * sr2 * sr2
                            sr12 = sr6 * sr6
                            e += 4 * self.eps * (sr12 - sr6) - self.shift
                            g = -24 * self.eps * (2 * sr12 - sr6) / r2  # (dV/dr)/r
                            d = (d0, d1, d2)
                            for k in range(3):
                                F[i][k] += g * d[k]
                                F[j][k] -= g * d[k]
                                for l in range(3):
                                    W[k][l] += g * d[k] * d[l]
        return {"energy": e, "forces": F, "virial": W}


class PySCF:
    """Periodic Kohn-Sham DFT (PySCF), as in quantum/pyscf_dft_probe.py: GTH pseudopotentials, k-point mesh, Fermi smearing."""

    def __init__(self, o):
        from pyscf import lib  # noqa: F401  (ImportError -> clean refusal)
        from pyscf.pbc import gto, dft  # noqa: F401
        from pyscf.pbc.grad import krks_stress  # noqa: F401
        if o.threads:
            lib.num_threads(o.threads)
        self.o = o
        self.kmesh = [int(x) for x in o.kpts.split(",")]
        if len(self.kmesh) != 3 or min(self.kmesh) < 1:
            raise ValueError("pyscf: --kpts must be three positive integers, e.g. 2,2,2")
        self.name = "pyscf(%s/%s/%s,k=%s,sigma=%gHa)" % (o.xc.upper(), o.basis, o.pseudo, o.kpts, o.smearing)
        self.cutoff = 0.0
        self.dm = None
        self.last_key = None

    def masses(self, elements):
        from pyscf.data import elements as el
        return [float(el.MASSES[el.charge(s)]) for s in elements]

    def evaluate(self, cell_a, pos, species, symbols):
        import numpy as np
        from pyscf.pbc import gto, dft
        from pyscf.pbc.scf.addons import smearing_
        from pyscf.pbc.grad import krks_stress
        o = self.o
        cell = gto.Cell()
        cell.atom = [(s, tuple(p)) for s, p in zip(symbols, pos)]
        cell.a = cell_a
        cell.unit = "A"
        cell.basis = o.basis
        cell.pseudo = o.pseudo
        cell.verbose = 0
        cell.build()
        kpts = cell.make_kpts(self.kmesh)
        mf = smearing_(dft.KRKS(cell, kpts).set(xc=o.xc, conv_tol=o.conv_tol), sigma=o.smearing, method="fermi")
        key = tuple(symbols)
        dm0 = self.dm if (self.dm is not None and key == self.last_key) else None
        e = mf.kernel(dm0=dm0)
        if not mf.converged:
            raise RuntimeError("SCF did not converge")
        self.dm, self.last_key = mf.make_rdm1(), key
        g = mf.nuc_grad_method()
        grad = np.asarray(g.kernel())                      # dE/dR, Ha/Bohr
        sigma = np.asarray(krks_stress.kernel(g))          # (1/V) dE/d(eps), Ha/Bohr^3
        forces = (-grad * HARTREE_EV / BOHR_A).tolist()
        stress = (sigma * HARTREE_EV / BOHR_A ** 3)
        stress = (0.5 * (stress + stress.T)).tolist()
        log("E = %.8f Ha  max|g| = %.3e  cycles ok" % (e, float(np.abs(grad).max())))
        return {"energy": float(e) * HARTREE_EV, "forces": forces, "stress": stress}


class MACE:
    """MACE foundation model through ASE (optional)."""

    def __init__(self, o):
        from mace.calculators import mace_mp
        self.calc = mace_mp(model=o.model, device=o.device, default_dtype="float64")
        self.name = "mace(%s)" % o.model
        self.cutoff = 0.0

    def masses(self, elements):
        from ase.data import atomic_masses, atomic_numbers
        return [float(atomic_masses[atomic_numbers[s]]) for s in elements]

    def evaluate(self, cell, pos, species, symbols):
        from ase import Atoms
        at = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True)
        at.calc = self.calc
        e = float(at.get_potential_energy())
        f = at.get_forces().tolist()
        s = at.get_stress(voigt=False).tolist()  # ASE: (1/V) dE/d(eps), eV/A^3 -- same convention
        return {"energy": e, "forces": f, "stress": s}


BACKENDS = {"lj": LJ, "pyscf": PySCF, "mace": MACE}
REQUIRES = {"lj": [], "pyscf": ["numpy", "pyscf"], "mace": ["numpy", "ase", "mace"]}


def check(backend):
    import importlib.util
    missing = [m for m in REQUIRES[backend] if importlib.util.find_spec(m) is None]
    if missing:
        print("%s backend unavailable: missing %s" % (backend, ", ".join(missing)))
        return SKIP
    print("%s backend available" % backend)
    return 0


def send(obj):
    sys.stdout.write(json.dumps(obj) + "\n")  # json writes NaN for non-finite floats; the bridge rejects it
    sys.stdout.flush()


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("backend", nargs="?", choices=sorted(BACKENDS))
    ap.add_argument("--check", metavar="BACKEND", choices=sorted(BACKENDS))
    ap.add_argument("--epsilon", type=float, default=0.0104)
    ap.add_argument("--sigma", type=float, default=3.40)
    ap.add_argument("--rc", type=float, default=8.5)
    ap.add_argument("--mass", type=float, default=39.948)
    ap.add_argument("--basis", default="gth-szv")
    ap.add_argument("--pseudo", default="gth-pbe")
    ap.add_argument("--xc", default="pbe")
    ap.add_argument("--kpts", default="2,2,2")
    ap.add_argument("--smearing", type=float, default=0.01, help="Fermi smearing width in Hartree")
    ap.add_argument("--conv-tol", type=float, default=1e-9)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--model", default="medium")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--fault", choices=["crash", "hang", "nan", "shape", "garbage", "badid", "error"],
                    help="testing aid: misbehave on evaluate request number --fault-after+1")
    ap.add_argument("--fault-after", type=int, default=0)
    o = ap.parse_args(argv)
    if o.check:
        return check(o.check)
    if not o.backend:
        ap.error("a backend is required")

    try:
        model = BACKENDS[o.backend](o)
    except ImportError as e:
        log("%s backend unavailable: %s" % (o.backend, e))
        line = sys.stdin.readline()
        if line:
            send({"protocol": PROTOCOL, "error": "%s backend unavailable: %s" % (o.backend, e)})
        return SKIP
    log("started %s" % model.name)

    nevals = 0
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            req = json.loads(line)
        except ValueError as e:
            log("bad request: %s" % e)
            return 2
        op = req.get("op")
        if op == "exit":
            break
        if op == "hello":
            if req.get("protocol") != PROTOCOL:
                send({"protocol": PROTOCOL, "error": "unsupported protocol %r" % req.get("protocol")})
                return 2
            try:
                masses = model.masses(req.get("elements", []))
            except Exception as e:  # unknown element symbols etc.
                send({"protocol": PROTOCOL, "error": "masses: %s" % e})
                return 2
            send({"protocol": PROTOCOL, "name": model.name, "cutoff": model.cutoff, "masses": masses})
            continue
        if op != "evaluate":
            send({"id": req.get("id"), "error": "unknown op %r" % op})
            continue
        rid = req.get("id")
        n = len(req["positions"])
        fault = o.fault if nevals == o.fault_after else None
        nevals += 1
        if fault == "crash":
            log("fault injection: exiting with status 3")
            return 3
        if fault == "hang":
            log("fault injection: hanging")
            import time
            while True:
                time.sleep(1)
        if fault == "garbage":
            sys.stdout.write("this is not json\n")
            sys.stdout.flush()
            continue
        try:
            out = model.evaluate(req["cell"], req["positions"], req["species"], req["symbols"])
        except Exception as e:
            log("evaluate failed: %s" % e)
            send({"id": rid, "error": str(e)})
            continue
        if fault == "error":
            send({"id": rid, "error": "injected worker-side error"})
            continue
        if fault == "nan":
            out["energy"] = float("nan")
        if fault == "shape":
            out["forces"] = out["forces"][:max(0, n - 1)]
        if fault == "badid":
            rid = (rid or 0) + 1000
        out["id"] = rid
        send(out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
