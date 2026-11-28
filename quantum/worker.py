#!/usr/bin/env python3
"""atomforge external-potential worker (protocol 1).

Reads one JSON request per line on stdin and writes one JSON reply per line on stdout. Diagnostics go to stderr, which
ExternalPotential appends to its log file. The protocol is documented in docs/EXTERNAL_PROTOCOL.md.

    python3 quantum/worker.py lj    [--epsilon EV] [--sigma A] [--rc A] [--mass AMU]
    python3 quantum/worker.py pyscf [--basis gth-szv] [--pseudo gth-pbe] [--xc pbe] [--kpts 2,2,2] [--smearing HA]
                                    [--ke-cutoff HA] [--spin-polarised] [--multigrid] [--threads N]
    python3 quantum/worker.py mace  [--model medium|FILE] [--device cpu] [--threads N]
    python3 quantum/worker.py --check BACKEND     exit 0 if BACKEND can run here, 77 if its packages are missing

Units on the wire: A, eV, eV/A. "virial" is W_ab = dE/d(eps_ab) in eV (tensile positive); "stress" is W/V in eV/A^3.

The lj backend is pure Python (no numpy) so the C++ tests can run it anywhere python3 exists. Its --fault options exist
only to exercise the error paths of the bridge.
"""
import argparse
import json
import math
import os
import sys
import threading
import time

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
        # --epsilon / --sigma: one value for every species, or one per element (in handshake order), mixed by the
        # Lorentz-Berthelot rule, eps_ab = sqrt(eps_a eps_b), sigma_ab = (sigma_a + sigma_b) / 2, as in lj.hpp.
        self.eps1 = [float(x) for x in str(o.epsilon).split(",")]
        self.sig1 = [float(x) for x in str(o.sigma).split(",")]
        self.rc, self.m = o.rc, o.mass
        if min(self.eps1) <= 0 or min(self.sig1) <= 0 or not self.rc > max(self.sig1):
            raise ValueError("lj: need epsilon > 0, sigma > 0, rc > sigma")
        self.stress_reply = o.reply == "stress"
        self.cutoff = self.rc
        self.name = "lj(eps=%s,sigma=%s,rc=%g)" % (o.epsilon, o.sigma, self.rc)
        if len(self.eps1) == 1 and len(self.sig1) == 1:
            self.name = "lj(eps=%g,sigma=%g,rc=%g)" % (self.eps1[0], self.sig1[0], self.rc)
        self.setup(max(len(self.eps1), len(self.sig1)))   # redone for the real element count at the handshake

    def setup(self, nel):
        def per(v):
            if len(v) == 1:
                return v * nel
            if len(v) != nel:
                raise ValueError("lj: %d per-element values for %d elements" % (len(v), nel))
            return v
        e1, s1 = per(self.eps1), per(self.sig1)
        self.E = [[math.sqrt(e1[a] * e1[b]) for b in range(nel)] for a in range(nel)]
        self.S = [[0.5 * (s1[a] + s1[b]) for b in range(nel)] for a in range(nel)]
        self.SH = [[0.0] * nel for _ in range(nel)]
        for a in range(nel):
            for b in range(nel):
                sr6 = math.pow(self.S[a][b] / self.rc, 6)
                self.SH[a][b] = 4 * self.E[a][b] * (sr6 * sr6 - sr6)

    def masses(self, elements):
        self.setup(len(elements))
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
                            sa, sb = species[i], species[j]
                            sig, eps = self.S[sa][sb], self.E[sa][sb]
                            sr2 = sig * sig / r2
                            sr6 = sr2 * sr2 * sr2
                            sr12 = sr6 * sr6
                            e += 4 * eps * (sr12 - sr6) - self.SH[sa][sb]
                            g = -24 * eps * (2 * sr12 - sr6) / r2  # (dV/dr)/r
                            d = (d0, d1, d2)
                            for k in range(3):
                                F[i][k] += g * d[k]
                                F[j][k] -= g * d[k]
                                for l in range(3):
                                    W[k][l] += g * d[k] * d[l]
        if self.stress_reply:  # the form the pyscf and mace backends use: W / V in eV/A^3
            return {"energy": e, "forces": F, "stress": [[w / vol for w in row] for row in W]}
        return {"energy": e, "forces": F, "virial": W}


class PySCF:
    """Periodic Kohn-Sham DFT (PySCF), as in quantum/pyscf_dft_probe.py: GTH pseudopotentials, k-point mesh, Fermi smearing.

    Spin-restricted (KRKS) by default. --spin-polarised runs KUKS from an initial guess with more spin-up than spin-down
    density; with smearing there is one Fermi level for both spins, so the moment is free to settle or to vanish."""

    def __init__(self, o):
        import pyscf
        from pyscf import lib  # noqa: F401  (ImportError -> clean refusal)
        from pyscf.pbc import gto, dft  # noqa: F401
        from pyscf.pbc.grad import krks_stress, kuks_stress  # noqa: F401
        if o.threads:
            lib.num_threads(o.threads)
        self.o = o
        self.kmesh = [int(x) for x in o.kpts.split(",")]
        if len(self.kmesh) != 3 or min(self.kmesh) < 1:
            raise ValueError("pyscf: --kpts must be three positive integers, e.g. 2,2,2")
        if o.ke_cutoff is not None and not o.ke_cutoff > 0:
            raise ValueError("pyscf: --ke-cutoff must be positive (Hartree)")
        self.name = "pyscf(%s/%s/%s,k=%s,sigma=%gHa" % (o.xc.upper(), o.basis, o.pseudo, o.kpts, o.smearing)
        if o.ke_cutoff is not None:
            self.name += ",cut=%gHa" % o.ke_cutoff
        if o.spin_polarised:
            self.name += ",spin-polarised"
        if o.multigrid:
            self.name += ",multigrid"
        self.name += ")"
        self.cutoff = 0.0
        self.dm = None
        self.last_key = None
        log("pyscf %s, %d thread(s)" % (pyscf.__version__, lib.num_threads()))

    def masses(self, elements):
        from pyscf.data import elements as el
        return [float(el.MASSES[el.charge(s)]) for s in elements]

    def scf(self, cell_a, pos, symbols):
        """Run the SCF and return (mf, info); info has the energies in Hartree, the cell moment and the timings."""
        import numpy as np
        from pyscf.pbc import gto, dft
        from pyscf.pbc.scf.addons import smearing_
        o = self.o
        t0 = time.time()
        cell = gto.Cell()
        cell.atom = [(s, tuple(p)) for s, p in zip(symbols, pos)]
        cell.a = cell_a
        cell.unit = "A"
        cell.basis = o.basis
        cell.pseudo = o.pseudo
        cell.verbose = 0
        if o.ke_cutoff is not None:
            cell.ke_cutoff = o.ke_cutoff
        cell.build()
        kpts = cell.make_kpts(self.kmesh)
        if o.spin_polarised:
            # PySCF counts the unpaired electrons of KUKS over all k-points. Only the parity has to be consistent here:
            # the smearing sets the occupations, with one Fermi level for both spins.
            cell.spin = (cell.nelectron * len(kpts)) % 2
        elif (cell.nelectron * len(kpts)) % 2:
            # Restricted orbitals hold electron pairs. PySCF's smeared KRKS then rounds the count up to the next pair
            # without an error: fcc Al (3 electrons) on 3x3x3 k-points is computed with 82 electrons on 27 k-points,
            # 3.037 per cell, a charged cell with an energy about 0.3 eV too high.
            raise RuntimeError("%d electrons per cell on %d k-points is an odd number: spin-restricted DFT would run a "
                               "charged cell. Use --spin-polarised, or a k-mesh with an even number of points"
                               % (cell.nelectron, len(kpts)))
        ks = dft.KUKS(cell, kpts) if o.spin_polarised else dft.KRKS(cell, kpts)
        ks.set(xc=o.xc, conv_tol=o.conv_tol, max_cycle=o.max_cycle)
        self.numint = ks._numint
        if o.multigrid:
            # The multigrid integrator builds the density and the Coulomb and exchange-correlation matrices on a
            # hierarchy of grids. For Ni (semi-core GTH potential, 85^3 grid points in B2 NiAl) the SCF is about 17
            # times faster and the energy agrees with the single-grid one to 1e-9 Ha. It has no k-point gradients, so
            # evaluate() puts the standard integrator back for the forces and the stress.
            from pyscf.pbc.dft import multigrid
            ks._numint = multigrid.MultiGridNumInt(cell)
        mf = smearing_(ks, sigma=o.smearing, method="fermi")
        key = (tuple(symbols), bool(o.spin_polarised), tuple(int(x) for x in cell.mesh))
        dm0 = self.dm if (self.dm is not None and key == self.last_key) else None
        if dm0 is None and o.spin_polarised:
            dm0 = np.array(mf.get_init_guess())            # (2, nk, nao, nao), equal spins
            dm0[0] *= 1.0 + o.init_polarisation
            dm0[1] *= 1.0 - o.init_polarisation
        ncyc = [0]
        mf.callback = lambda env: ncyc.__setitem__(0, ncyc[0] + 1)
        e = mf.kernel(dm0=dm0)
        if not mf.converged:
            raise RuntimeError("SCF did not converge in %d cycles" % o.max_cycle)
        dm = mf.make_rdm1()
        self.dm, self.last_key = dm, key
        moment = 0.0
        if o.spin_polarised:
            s1e = mf.get_ovlp()
            moment = float(sum(np.einsum("ij,ji->", dm[0][k] - dm[1][k], s1e[k]).real for k in range(len(kpts))) / len(kpts))
        # With smearing the SCF minimises the Mermin free energy F = E - sigma*S, and the analytic forces and stress are
        # derivatives of F, not of E. F is the energy the bridge gets, so that energy, forces and stress are consistent.
        f = getattr(mf, "e_free", None)
        f = float(e) if f is None else float(f)
        ez = getattr(mf, "e_zero", None)
        ez = float(e) if ez is None else float(ez)
        info = {"e": float(e), "f": f, "e_zero": ez, "moment": moment, "cycles": ncyc[0], "nao": int(cell.nao),
                "mesh": [int(x) for x in cell.mesh], "nk": len(kpts), "t_scf": time.time() - t0}
        return mf, info

    def evaluate(self, cell_a, pos, species, symbols):
        import numpy as np
        from pyscf.pbc.grad import krks_stress, kuks_stress
        mf, info = self.scf(cell_a, pos, symbols)
        t0 = time.time()
        mf._numint = self.numint
        g = mf.nuc_grad_method()
        grad = np.asarray(g.kernel())                      # dF/dR, Ha/Bohr
        t1 = time.time()
        stress_mod = kuks_stress if self.o.spin_polarised else krks_stress
        sigma = np.asarray(stress_mod.kernel(g))           # (1/V) dF/d(eps), Ha/Bohr^3
        t2 = time.time()
        forces = (-grad * HARTREE_EV / BOHR_A).tolist()
        stress = (sigma * HARTREE_EV / BOHR_A ** 3)
        stress = (0.5 * (stress + stress.T)).tolist()
        log("E = %.8f Ha  F = %.8f Ha  max|g| = %.3e  moment = %.3f muB  %d cycles  scf %.1f s  grad %.1f s  stress %.1f s"
            % (info["e"], info["f"], float(np.abs(grad).max()), info["moment"], info["cycles"], info["t_scf"], t1 - t0, t2 - t1))
        out = {"energy": info["f"] * HARTREE_EV, "energy_sigma0": info["e_zero"] * HARTREE_EV, "forces": forces, "stress": stress}
        if self.o.spin_polarised:
            out["moment"] = info["moment"]
        return out


class MACE:
    """MACE foundation model through ASE (optional). A machine-learned interatomic potential fitted to DFT data; not DFT."""

    def __init__(self, o):
        import ase
        import mace
        import torch
        from mace.calculators import mace_mp
        if o.threads:
            torch.set_num_threads(o.threads)
        # --model is a name known to mace_mp ("medium", "medium-mpa-0", ...), which downloads the file on first use, or
        # the path of a model file. The file's SHA-256 goes to the log, so that a recorded run identifies its model.
        self.calc = mace_mp(model=o.model, device=o.device, default_dtype="float64")
        label = o.model
        if os.path.isfile(o.model):
            import hashlib
            h = hashlib.sha256()
            with open(o.model, "rb") as f:
                for block in iter(lambda: f.read(1 << 20), b""):
                    h.update(block)
            label = os.path.basename(o.model)
            log("mace model file %s  sha256 %s  %d bytes" % (label, h.hexdigest(), os.path.getsize(o.model)))
        log("mace-torch %s, torch %s, ase %s, %d thread(s), device %s"
            % (mace.__version__, torch.__version__, ase.__version__, torch.get_num_threads(), o.device))
        self.name = "mace(%s)" % label
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
REQUIRES = {"lj": [], "pyscf": ["numpy", "pyscf", "pyscf.pbc.grad.krks_stress", "pyscf.pbc.grad.kuks_stress"], "mace": ["numpy", "ase", "mace.calculators"]}


def check(backend):
    import importlib
    missing = []
    for m in REQUIRES[backend]:
        try:
            importlib.import_module(m)
        except Exception:  # ImportError, or a broken install
            missing.append(m)
    if missing:
        print("%s backend unavailable: missing %s" % (backend, ", ".join(missing)))
        return SKIP
    print("%s backend available" % backend)
    return 0


# The protocol channel. main() moves it off fd 1 before any backend is imported, so that library output (mace_mp prints
# to stdout, C extensions may write to fd 1) goes to stderr, i.e. the log, and cannot corrupt the replies.
PROTO = None


def claim_stdout():
    global PROTO
    fd = os.dup(1)
    os.dup2(2, 1)
    sys.stdout = sys.stderr
    PROTO = os.fdopen(fd, "w", buffering=1)


def send(obj):
    PROTO.write(json.dumps(obj) + "\n")  # json writes NaN for non-finite floats; the bridge rejects it
    PROTO.flush()


def watchdog():
    """Exit when the parent (atomforge) goes away, even in the middle of a long evaluation. Covers systems without
    PR_SET_PDEATHSIG and commands where a shell sits between atomforge and this process. It is a Python thread, so a
    call that holds the GIL delays it until the call returns; atomforge's supervisor kills the process group regardless."""
    parent = os.getppid()

    def run():
        while True:
            time.sleep(0.5)
            if os.getppid() != parent:
                os._exit(70)
    threading.Thread(target=run, daemon=True).start()


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("backend", nargs="?", choices=sorted(BACKENDS))
    ap.add_argument("--check", metavar="BACKEND", choices=sorted(BACKENDS))
    ap.add_argument("--epsilon", default="0.0104", help="eV; one value, or one per element separated by commas")
    ap.add_argument("--sigma", default="3.40", help="A; one value, or one per element separated by commas")
    ap.add_argument("--reply", choices=["virial", "stress"], default="virial", help="lj: reply with the virial or with W/V")
    ap.add_argument("--rc", type=float, default=8.5)
    ap.add_argument("--mass", type=float, default=39.948)
    ap.add_argument("--basis", default="gth-szv")
    ap.add_argument("--pseudo", default="gth-pbe")
    ap.add_argument("--xc", default="pbe")
    ap.add_argument("--kpts", default="2,2,2")
    ap.add_argument("--smearing", type=float, default=0.01, help="Fermi smearing width in Hartree")
    ap.add_argument("--conv-tol", type=float, default=1e-9)
    ap.add_argument("--max-cycle", type=int, default=100, help="pyscf: SCF cycle limit")
    ap.add_argument("--ke-cutoff", type=float, default=None, help="pyscf: plane-wave cutoff of the density grid in Hartree "
                    "(default: PySCF's estimate from the basis)")
    ap.add_argument("--multigrid", action="store_true", help="pyscf: run the SCF with PySCF's multigrid integrator "
                    "(same energy, much faster for hard pseudopotentials); forces and stress use the standard one")
    ap.add_argument("--spin-polarised", action="store_true", help="pyscf: unrestricted (KUKS) with a polarised initial guess")
    ap.add_argument("--init-polarisation", type=float, default=0.2,
                    help="pyscf --spin-polarised: the initial guess has (1+p)/2 of the density spin-up")
    ap.add_argument("--threads", type=int, default=0, help="pyscf, mace: number of threads (0: the library default)")
    ap.add_argument("--model", default="medium")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--fault", choices=["crash", "hang", "nan", "shape", "garbage", "badid", "error", "spawn-exit", "errnull"],
                    help="testing aid: misbehave on evaluate request number --fault-after+1")
    ap.add_argument("--fault-after", type=int, default=0)
    ap.add_argument("--delay", type=float, default=0.0, help="testing aid: sleep this many seconds before each reply")
    ap.add_argument("--noisy", action="store_true", help="testing aid: print to stdout (Python and fd 1) at start and per request")
    ap.add_argument("--pid-file", help="testing aid: write this process's pid to the file")
    return ap


def main(argv):
    ap = build_parser()
    o = ap.parse_args(argv)
    if o.check:
        return check(o.check)
    if not o.backend:
        ap.error("a backend is required")
    claim_stdout()
    watchdog()
    if o.pid_file:
        with open(o.pid_file, "w") as f:
            f.write("%d\n" % os.getpid())
    if o.noisy:
        print("noisy: Python print to stdout during construction")
        os.write(1, b"noisy: raw write to fd 1 during construction\n")

    try:
        model = BACKENDS[o.backend](o)
    except ImportError as e:
        log("%s backend unavailable: %s" % (o.backend, e))
        line = sys.stdin.readline()
        if line:
            send({"protocol": PROTOCOL, "error": "%s backend unavailable: %s" % (o.backend, e)})
        return SKIP
    except Exception as e:  # bad options
        log("%s backend: %s" % (o.backend, e))
        if sys.stdin.readline():
            send({"protocol": PROTOCOL, "error": "%s backend: %s" % (o.backend, e)})
        return 2
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
        if o.delay > 0:
            time.sleep(o.delay)
        if fault == "crash":
            log("fault injection: exiting with status 3")
            return 3
        if fault == "hang":
            log("fault injection: hanging")
            while True:
                time.sleep(1)
        if o.noisy:
            print("noisy: print during evaluate %d" % nevals)
            os.write(1, b"noisy: raw fd 1 write during evaluate\n")
        if fault == "spawn-exit":  # leave a child holding the socket, then die without replying
            import subprocess
            subprocess.Popen(["sleep", "30"], stdin=sys.stdin.fileno(), stdout=PROTO.fileno(), stderr=subprocess.DEVNULL)
            log("fault injection: started a child and exiting with status 5")
            os._exit(5)
        if fault == "garbage":
            PROTO.write("this is not json\n")
            PROTO.flush()
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
        if fault == "errnull":
            out["error"] = None
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
