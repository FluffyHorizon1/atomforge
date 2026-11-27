# External worker protocol, version 1

`atomforge::ExternalPotential` (`include/atomforge/external.hpp`) gets energies, forces and virials from a separate
process. The reference worker is `quantum/worker.py`. Anything that follows this page can be a worker: a Python script around
PySCF, GPAW, Quantum ESPRESSO or a machine-learned potential, or a compiled program.

## Process model

- atomforge starts the worker once, with `/bin/sh -c "exec <command>"`, when the `ExternalPotential` is constructed. The worker
  stays alive for every evaluation of a relaxation or search and is told to exit when the potential is destroyed.
- The worker's stdin and stdout are one end of a Unix socket pair. Each message is one JSON object on one line, terminated by
  `\n`. The worker must flush stdout after every reply.
- The worker's stderr is appended to a log file (`ExternalOptions::log_path`, `--worker-log` on the CLI; default
  `$TMPDIR/atomforge-worker-<pid>.log`). Error messages from atomforge quote the last 1.5 kB of it.
- Requests are strictly sequential: atomforge sends one request and waits for its reply before sending the next.

## Messages

### Handshake

```json
-> {"op":"hello","protocol":1,"elements":["Cu","Au"]}
<- {"protocol":1,"name":"pyscf(PBE/gth-szv/gth-pbe,k=2,2,2,sigma=0.01Ha)","cutoff":0,"masses":[63.546,196.97]}
```

`protocol` is required and must equal 1. `name` (string), `cutoff` (Å, informational; 0 if none) and `masses` (amu, one per
element, used by molecular dynamics) are optional. Without `masses`, atomforge uses the Zhou04 table for elements it knows,
or `ExternalOptions::masses`. A worker that cannot serve the elements replies `{"protocol":1,"error":"..."}`.

### Evaluate

```json
-> {"op":"evaluate","id":7,"natoms":2,
    "cell":[[3.6,0,0],[0,3.6,0],[0,0,3.6]],
    "positions":[[0,0,0],[1.8,1.8,0.1]],
    "species":[0,1],"symbols":["Cu","Au"]}
<- {"id":7,"energy":-7.12,"forces":[[0,0,0.01],[0,0,-0.01]],"virial":[[..],[..],[..]]}
```

- `cell`: rows are the lattice vectors, in Å. `positions`: Cartesian, Å, not necessarily inside the home cell.
- `species` indexes the handshake's `elements`; `symbols` repeats the same information as element symbols.
- `id` increases by one per request, starting at 1. The reply must carry the same `id`.
- `energy`: total energy of the cell, eV.
- `forces`: one `[fx,fy,fz]` per atom, in the order of `positions`, eV/Å, `F = -dE/dR`.
- `virial`: 3×3, eV, `W_ab = dE/d(eps_ab)` for a symmetric strain `eps` applied to cell and positions (`r -> (I + eps) r`).
  This is the convention of `ZhouEAM`: W is positive along a direction in which stretching raises the energy, so a cell
  under compression (wanting to expand) has a negative trace. The stress is `W / V` (tensile positive).
- Instead of `virial` a worker may send `stress`, 3×3 in eV/Å³, equal to `W / V` (same sign; this is the ASE and PySCF
  convention `(1/V) dE/d(eps)`). atomforge multiplies by the cell volume it sent.

A worker that cannot evaluate a structure (for example an SCF that does not converge) replies

```json
<- {"id":7,"error":"SCF did not converge"}
```

atomforge then throws `ExternalPotentialError` with that message, and the worker stays in use for later requests.

### Exit

```json
-> {"op":"exit"}
```

atomforge also closes the worker's stdin. A worker should exit when it reads `exit` or end of file. A worker that has not
exited after 2 s is killed.

Numbers are written with 17 significant digits (round-trip exact for IEEE doubles).

## Failure handling

`evaluate()` throws `atomforge::ExternalPotentialError` (a `std::runtime_error`) and never blocks indefinitely when:

| condition | message contains |
|-----------|------------------|
| the command cannot be started or exits early (e.g. not found: shell status 127) | `exit status N` |
| no handshake within `hello_timeout_s` (default 60 s) | `timed out` |
| a protocol version other than 1 | `worker speaks protocol N` |
| the worker exits or closes stdout during a request | `worker closed its output`, plus the exit status or signal |
| no reply within `eval_timeout_s` (default 3600 s; `--worker-timeout`) | `timed out` |
| a reply that is not a JSON object | `malformed JSON` |
| a reply with another `id` | `does not match` |
| a missing field, the wrong number of force vectors, a non-3×3 virial | the field name and the expected shape |
| `NaN` or `±Infinity` anywhere in energy, forces or virial | `not finite` |
| an `error` reply | `worker reported an error` |

In every case except an `error` reply, the worker is killed (SIGKILL) and later calls throw `external worker is not
running`; atomforge does not restart workers. Messages include the tail of the stderr log.

## The reference worker

```
python3 quantum/worker.py lj    [--epsilon 0.0104] [--sigma 3.40] [--rc 8.5] [--mass 39.948]
python3 quantum/worker.py pyscf [--basis gth-szv] [--pseudo gth-pbe] [--xc pbe] [--kpts 2,2,2] [--smearing 0.01] [--threads N]
python3 quantum/worker.py mace  [--model medium] [--device cpu]
python3 quantum/worker.py --check pyscf     # exit 0 if usable here, 77 if its packages are missing
```

- `lj` is the truncated-and-shifted Lennard-Jones potential of `include/atomforge/lj.hpp`, in pure Python with periodic
  images, so it runs wherever `python3` does. It exists to test the bridge; `--fault crash|hang|nan|shape|garbage|badid|error`
  makes it misbehave on purpose for the error-path tests.
- `pyscf` runs periodic Kohn–Sham DFT with the settings of `quantum/pyscf_dft_probe.py` (GTH pseudopotentials, k-point
  mesh, Fermi smearing in Hartree), returns analytic forces and the analytic stress, and reuses the previous density matrix as
  the initial guess when the atoms are the same species in the same order.
- `mace` uses the MACE foundation model through ASE, in double precision.

A missing optional package makes the worker refuse the handshake with `<backend> backend unavailable: ...` and exit 77.
