# External worker protocol, version 1

`atomforge::ExternalPotential` (`include/atomforge/external.hpp`) gets energies, forces and virials from a separate
process. The reference worker is `quantum/worker.py`. Anything that follows this page can be a worker: a Python script around
PySCF, GPAW, Quantum ESPRESSO or a machine-learned potential, or a compiled program.

## Process model

- Platforms: POSIX systems with `/bin/sh` (Linux and macOS; tested on Linux). It uses fork/exec, a Unix socket pair and
  pipes. `close_range` is used where the kernel has it and a loop over the descriptor limit elsewhere; `SO_NOSIGPIPE`
  replaces `MSG_NOSIGNAL` where the latter is missing. Numbers go through floating-point `std::to_chars`/`from_chars` when
  the standard library has them, and through classic-locale streams otherwise.
- atomforge starts the worker once, when the `ExternalPotential` is constructed, as `/bin/sh -c "exec <command>"`. The
  `exec` is left out when the command contains shell operators (`; & | < > ( ) $` or a backtick) or starts with a builtin
  such as `cd`, so `cd quantum && python3 worker.py lj` works. A command that starts with variable assignments
  (`OMP_NUM_THREADS=1 python3 ...`, the usual way to set DFT thread counts) is run as `exec env ...`. The worker stays
  alive for every evaluation of a relaxation or search and is told to exit when the potential is destroyed.
- atomforge forks once and immediately execs a supervisor: `/bin/sh` running a short script. The supervisor shares no
  memory, descriptors or signal handlers with the host. Before the exec the child resets every signal disposition to the
  default, clears the signal mask and closes every descriptor except the five it passes on: /dev/null, the log, the socket
  and a "lifeline" pipe whose only write end the host keeps. The supervisor leads a process group of its own, which
  contains the worker and everything the worker starts.
  - When the worker ends, the supervisor ends the same way: with the same exit status, or by re-raising the same signal
    with core dumps turned off. atomforge therefore reports the worker's own status (`exit status N` or
    `killed by signal N`). A worker that exits normally with a status between 129 and 159 is indistinguishable from one
    killed by a signal.
  - When the host goes away for any reason (normal exit, Ctrl-C, `timeout`, SIGHUP, SIGKILL, a crash), the lifeline reads
    end of file and the supervisor's watcher kills the whole group at once.
  - The group is not the terminal's, so Ctrl-C reaches the worker only this way. As a background job of a non-interactive
    shell, the worker also starts with SIGINT and SIGQUIT ignored.
  - `quantum/worker.py` additionally exits by itself when its parent process disappears. Its watchdog is a Python thread,
    so a worker stuck inside a call that holds the GIL notices only when the call returns (about 3.5 s in the reviewer's
    probe). The supervisor's group kill does not depend on it.
- The worker's stdin and stdout are one end of a Unix socket pair. Each message is one JSON object on one line, terminated by
  `\n`. The worker must flush after every reply and must not write anything else to that channel. `worker.py` keeps a
  private copy of fd 1 for replies and points fd 1 and `sys.stdout` at stderr before importing a backend, so library output
  (`mace_mp` prints to stdout, for example) ends up in the log.
- The worker's stderr is appended to a log file (`ExternalOptions::log_path`, `--worker-log` on the CLI). The default is
  `$TMPDIR/atomforge-worker-<pid>-<n>.log`, a new file for each potential in the process. Each start writes a header line,
  and error messages quote at most the last 1.5 kB written after it, never an earlier run's output. A log path that cannot
  be opened is reported as such.
- Requests are strictly sequential: atomforge sends one request and waits for its reply before sending the next.
- Timeouts of 1e9 s or more, or infinite ones, mean no deadline.

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
- `energy`: total energy of the cell, eV. It must be the quantity that `forces` and `virial` are derivatives of: with
  electronic smearing that is the Mermin free energy `E - sigma*S`, not `E`.
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

Replies may carry extra fields (the `pyscf` backend adds `energy_sigma0`, the energy extrapolated to zero smearing);
atomforge ignores them. `"error": null` counts as no error. atomforge writes numbers in the shortest form that reads back
exactly (`std::to_chars`) and reads them with `std::from_chars`, so the C locale (a decimal comma, say) does not matter.

## Failure handling

`evaluate()` throws `atomforge::ExternalPotentialError` (a `std::runtime_error`) and never blocks indefinitely when:

| condition | message contains |
|-----------|------------------|
| the log file cannot be opened | `cannot open worker log` |
| the command cannot be started or exits early (e.g. not found: shell status 127) | `exit status N` |
| no handshake within `hello_timeout_s` (default 60 s) | `timed out` |
| a protocol version other than 1 | `worker speaks protocol N` |
| the worker exits or closes stdout during a request | `worker closed its output` or `worker exited`, plus the exit status or signal; noticed within about 1 s even when a process the worker started still holds the socket |
| no reply within `eval_timeout_s` (default 3600 s; `--worker-timeout`) | `timed out` |
| a reply that is not a JSON object | `malformed JSON` |
| a reply with another `id` | `does not match` |
| a missing field, the wrong number of force vectors, a non-3×3 virial | the field name and the expected shape |
| `NaN` or `±Infinity` anywhere in energy, forces or virial | `not finite` |
| an `error` reply | `worker reported an error` |
| a non-finite cell or position in the structure (not sent; the worker keeps running) | `not finite` |

In every case except an `error` reply and a rejected input, the worker is killed (SIGKILL to its process group: the worker is started in a
group of its own, so processes it started die with it) and later calls throw `external worker is not
running`; atomforge does not restart workers. Messages include the tail of the stderr log.

## The reference worker

```
python3 quantum/worker.py lj    [--epsilon 0.0104[,EPS_B,...]] [--sigma 3.40[,SIG_B,...]] [--rc 8.5] [--mass 39.948] [--reply virial|stress]
python3 quantum/worker.py pyscf [--basis gth-szv] [--pseudo gth-pbe] [--xc pbe] [--kpts 2,2,2] [--smearing 0.01] [--threads N]
python3 quantum/worker.py mace  [--model medium] [--device cpu]
python3 quantum/worker.py --check pyscf     # exit 0 if usable here, 77 if its packages are missing
```

- `lj` is the truncated-and-shifted Lennard-Jones potential of `include/atomforge/lj.hpp`, in pure Python with periodic
  images, so it runs wherever `python3` does. Per-element `--epsilon`/`--sigma` lists (in handshake order) are mixed by the
  Lorentz–Berthelot rule, as in `LJParams::epsilon_by_species`. `--reply stress` sends `stress` = W/V instead of `virial`.
  The backend exists to test the bridge. Several options make it misbehave on purpose for the error-path tests:
  `--fault crash|hang|nan|shape|garbage|badid|error|errnull|spawn-exit`, `--noisy` (prints to stdout) and `--pid-file`.
- `pyscf` runs periodic Kohn–Sham DFT with the settings of `quantum/pyscf_dft_probe.py` (GTH pseudopotentials, k-point
  mesh, Fermi smearing in Hartree), returns the free energy with its analytic forces and stress, and reuses the previous density matrix as
  the initial guess when the atoms are the same species in the same order.
  It is spin-restricted (KRKS): magnetic elements such as Fe, Co and Ni are computed non-magnetic.
- `mace` uses the MACE foundation model through ASE, in double precision. The model is downloaded on first use, inside
  the handshake, so the first start can need more than the default 60 s handshake timeout (`--worker-hello-timeout S`).

`--check pyscf` imports `pyscf.pbc.grad.krks_stress` too, and `--check mace` imports `mace.calculators`. A missing optional
package makes the worker refuse the handshake with `<backend> backend unavailable: ...` and exit 77.
