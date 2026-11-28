# atomforge — atomistic structure search for alloys (materials layer)

Header-only C++17. Mix elements, generate candidate crystal structures, relax them, and decide which are
stable: thermodynamically (convex hull), mechanically (Born criteria) and thermally (MD heat and quench).

```
elements ──► prototypes (L1₂, L1₀, L1₁, B2, D0₃, C15, B1) ──┐
         ──► random solid solutions (fcc/bcc/hcp)          ├─► FIRE relax (atoms + cell) ──► formation energy ──► convex hull
         ──► Monte Carlo chemical-ordering search          ┘                                                           │
                                                          stability screen: elastic tensor (Born) + MD heat/quench ◄──┘
```

## Build & run

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/atomforge validate          # 17 elements vs experiment; exits non-zero on any failure
./build/atomforge relax L1_2 Cu Au  # one structure + elastic constants
./build/atomforge search Cu Au      # full binary search -> results/CuAu/ (report.json, POSCAR, extxyz)
```

Options: `-DATOMFORGE_NATIVE=ON` adds `-march=native` (off by default, so binaries are portable);
`-DATOMFORGE_BUILD_TESTS=OFF` skips the `tests/` suite (the CLI `validate` CTest entry is always registered). Tests run with `ctest --test-dir build --output-on-failure`
(see [Testing](#testing)).

Every candidate is written as **POSCAR** (VASP / Quantum ESPRESSO via ASE) and **extended XYZ** (OVITO, ASE).
Hull structures can go straight to DFT for verification.

## Materials sidecar (`forge.materials.v0`)

```bash
./build/atomforge search Cu Au --materials   # also writes results/CuAu/materials-contract-CuAu.json
```

Writes the interchange sidecar proposed in qpuforge
[`docs/DATA_CONTRACT.md`](https://github.com/FluffyHorizon1/qpuforge/blob/f64238f/docs/DATA_CONTRACT.md)
(branch `docs/materials-data-contract`, commit `f64238f`); the writer is `include/atomforge/materials.hpp`.

- **Real (computed by the EAM model):** `produced_by`, `system`, `structures_ref`, and an additive optional block
  `atomforge_computed` with elemental references and hull phases (energies, formation energies, cell edges and angles, Born
  elastic constants). These are classical-EAM numbers, not DFT.
- **Placeholders (not computed, `value: null`, `status: "placeholder"`):** `Q_dielectric`, `loss_tangent`,
  `junction_gap_ueV`. They are listed in the top-level `placeholders` array. Do not copy them into a chip JSON.

The contract is still a proposal, so the format may change.

## Energy model

The Zhou–Johnson–Wadley analytic EAM (Phys. Rev. B 69, 144113, 2004) covers 17 metals: Cu Ag Au Ni Pd Pt Al Pb Fe
Mo Ta W Cr Mg Co Ti Zr. Alloy cross terms come from the paper's mixing rule. Parameters were generated from the
LAMMPS `tools/eam_database` table (GPL-2.0), see `tools/gen_params.py`. Forces and virial are analytic and are checked
against central finite differences of the energy and of the strain in `ctest` (worst errors about 8e-9 eV/Å and 4e-8 eV
on rattled Cu/Au cells, including a triclinic one). The `Potential` interface is abstract, so an ML potential or a
DFT backend can replace the EAM without touching the search, relaxation, hull or stability code.

## Validation

`atomforge validate` relaxes every elemental prototype for each of the 17 elements and compares with experiment. It
exits non-zero if any predicted ground state is wrong or any lattice constant is more than 1% from experiment.

* All 17 elements: lattice constants within 0.75% of experiment (largest deviation: Zr, −0.74%), and the correct ground
  state (fcc/bcc/hcp) for all 17. This is checked by the CTest `validate` test (the CLI itself). Cohesive energies are printed for
  comparison but are not pass/fail criteria.
* **Cu–Au:** the search independently rediscovers the experimentally known ordered phases Cu₃Au (L1₂, a = 3.750 Å,
  matching experiment), CuAu (L1₀) and CuAu₃ (L1₂). All are on the hull and pass the mechanical and thermal checks.
  It also predicts extra hcp-based phases, which are likely artefacts of the potential. (The test suite only checks the
  Cu₃Au lattice constant and, on a reduced search, that L1₂ Cu₃Au is bound and the hull is consistent. The full
  Cu–Au result above is a recorded run in `results/`, not a test.)
* **Cu–Ag:** all formation energies positive, so immiscible, as observed (recorded run in `results/`, not a test).
* **Ni–Al: fails.** The generic mixing rule gives formation energies near zero or positive (B2 NiAl: +0.07 eV/atom),
  where experiment gives about −0.6 eV/atom.
  Strongly bonding intermetallics need quantum-level (DFT) energies. That is the next rung.

## Testing

`ctest` runs `tests/test_atomforge.cpp` and `tests/test_external.cpp` (custom `CHECK` macros, one CTest entry per group, labels shown in brackets):

| group | what it checks |
|-------|----------------|
| `forces` [unit, derivatives] | analytic forces vs central finite differences on rattled fcc, L1₂ and triclinic L1₁ cells; net force = 0 |
| `virial` [unit, derivatives] | analytic virial vs finite strain (all six components, triclinic cell too); symmetry and sign convention |
| `invariance` [unit, symmetry] | energy/forces/virial under translation, rotation, lattice-vector shifts + wrap, and supercell size |
| `relax` [physics] | FIRE convergence on fcc Cu and rattled Cu₃Au, independent residual-force check, relaxation at 5 GPa |
| `born` [physics] | fcc Cu elastic constants and Born stability (bulk modulus cross-checked against E(V)); bcc Cu must be unstable |
| `hull` [physics, slow] | convex hull on synthetic points with a known answer, plus a reduced Cu–Au search |
| `validate` | `atomforge validate` (the CLI); exits 0 and prints `ground state correct for 17/17 elements`. Registered by the CI branch. |
| `validate_unit` [unit] | pass and fail paths of the element validation logic (wrong ground state, wrong lattice constant) |
| `unwrapped` [physics, symmetry] | energy and forces unchanged when atoms are displaced by whole lattice vectors (2-7 cells, cubic and triclinic), no wrapping |
| `lj` [unit, derivatives] | native Lennard-Jones forces and virial vs finite differences; dimer minimum at 2^(1/6) σ |
| `external` [unit, external] | `ExternalPotential` with the `lj` worker vs the native LJ: energy, forces, virial to 1e-9 (rattled, triclinic, unwrapped, two species) |
| `external_relax` [physics, external] | FIRE with full cell through the bridge: same steps, energy and volume as the native run; residual force and stress |
| `external_errors` [unit, external] | missing worker, early exit, handshake timeout, wrong protocol, crash, hang, NaN, wrong shape, bad JSON, wrong id, worker-side error, stderr log |
| `external_cli` [cli, external] | `atomforge relax --potential external:...` end to end against `--potential lj`, and its error exits |
| `md_invariance` [physics, symmetry] | `run_md` averages and final positions unchanged when the start structure's atoms are shifted by lattice vectors (including shifts that make `run_md` wrap) |
| `structure_misc` [unit] | composition, fraction and supercell bookkeeping; errors for unknown elements, prototypes and parent lattices |
| `io` [unit] | writers throw `cannot write` on unwritable paths; POSCAR, extended XYZ and `report.json` writers checked against the in-memory structure, `report.json` parsed as strict JSON (escaping, empty and real search results) |
| `cli` [unit, cli, slow] | the `atomforge` executable: unknown, missing, non-numeric and out-of-range `search` arguments fail fast; usage and error exit codes, `elements`, `relax` (stable and Born-unstable cases), and a reduced `search` whose `report.json` and structure files are cross-checked. Runs as `test_atomforge cli <path-to-atomforge>`; the executable checks are skipped without the path |

It does not cover the Monte Carlo ordering search beyond the reduced hull and CLI runs, or MD thermodynamics beyond the
shift-invariance check.
The `external*` groups need `python3` and are reported as skipped (exit 77) without it; they do not run the `pyscf` or
`mace` backends.

## Toward the quantum level

`quantum/pyscf_dft_probe.py` runs periodic Kohn–Sham DFT (PBE, GTH pseudopotentials, k-points) with PySCF.
It returns energy, analytic forces and the stress tensor, which is everything the FIRE relaxer needs.
On a 2-core cloud box, fcc Al with 27 k-points took about 60 s per evaluation (`results/pyscf_probe_log.txt`).

### External potentials (the DFT bridge)

`ExternalPotential` (`include/atomforge/external.hpp`) is a `Potential` that streams structures to a persistent worker
process and reads back the energy, forces and virial. The worker is started once and kept for the whole relaxation or search;
messages are line-delimited JSON, versioned (protocol 1) and specified in [docs/EXTERNAL_PROTOCOL.md](docs/EXTERNAL_PROTOCOL.md).
The worker's stderr goes to a log file. A worker that cannot start, exits, times out, or answers with malformed JSON, the
wrong shape or a NaN raises `ExternalPotentialError` with the end of that log; it never hangs.

`quantum/worker.py` is the reference worker, with three backends:

| backend | what it is | needs |
|---------|------------|-------|
| `lj` | the truncated-and-shifted Lennard-Jones of `include/atomforge/lj.hpp`, pure Python with periodic images | python3 |
| `pyscf` | periodic Kohn–Sham DFT with the probe's settings; analytic forces and stress | numpy, pyscf |
| `mace` | the MACE foundation model through ASE | numpy, ase, mace-torch |

From the command line:

```bash
./build/atomforge relax A1 Ar --volume 40 --potential lj            # native C++ Lennard-Jones
./build/atomforge relax A1 Ar --volume 40 \
    --potential external:"python3 quantum/worker.py lj"           # the same model through the bridge
./build/atomforge relax A1 Al --cell hydrostatic --fmax 0.01 \
    --potential external:"python3 quantum/worker.py pyscf --kpts 2,2,2 --threads 2" --worker-log al.log
```

Both LJ runs above give E = -0.07762 eV/atom and V = 36.5108 Å³/atom in 33 FIRE steps (34 worker calls). For external
potentials the elastic tensor is off by default (`--elastic` turns it on), because it costs at least 24 more evaluations (12 strained cells, each with an internal relaxation; 24 for fcc Ar with the `lj` worker).
`--worker-timeout S` bounds each evaluation (default 3600 s; `inf` for none). `--worker-hello-timeout S` bounds the
start-up and handshake (default 60 s; MACE downloads its model then). The `pyscf` backend is spin-restricted, so Fe, Co and
Ni come out non-magnetic. `relax` exits with 0 on success, 1 on usage and argument errors (unknown option or prototype, bad
value, extra arguments), and 2 on failures while running (unsupported element, a worker that fails), as `search` does.

In the tests, the bridge with the `lj` worker reproduces the native LJ energy, forces and virial to better than 1e-9 (the
differences were exactly zero in our runs) on rattled fcc, triclinic and unwrapped cells, and FIRE through the bridge
follows the native path step for step (79 steps, full cell).

The `pyscf` backend was checked against finite differences of its own energy (`quantum/check_pyscf_fd.py`, fcc Al,
gth-szv, 2×2×2 k-points, log in `results/pyscf_bridge_fd_log.txt`): the trace of the virial was -0.6362 eV analytic and
-0.6367 eV by finite strain, and a force on a displaced atom -0.11778 eV/Å analytic and -0.11777 eV/Å by finite
displacement. With Fermi smearing the worker returns the free energy E - σS, because that is what PySCF's forces and
stress are derivatives of; returning E instead put the virial off by a factor of 1.9. One evaluation took 26 s for the
1-atom cell and 131 s for the 4-atom cell on 2 threads. The CI tests do not install PySCF or MACE.

For QPU materials (Al/AlOx junction barriers, Nb and Ta surface oxides, TiN, the Si substrate) EAM is not enough.
Oxides, nitrides and covalent solids need DFT or reactive potentials.

## Limitations

The Monte Carlo search runs on 48-site cells. Dynamical stability is checked by MD heat-and-quench, not phonons.
There is no magnetism and no charge transfer (EAM). Alloy energies are qualitative.

The neighbour list accepts atoms anywhere in space: the periodic-image range grows with how far atoms lie outside the home
cell (earlier versions silently dropped pairs for atoms two or more cells away, or one cell away in a small triclinic cell).
The cost of building the list grows with that distance, so long simulations should keep positions wrapped; `run_md` wraps atoms
that drift more than one cell outside the home cell.

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This is proprietary software of Secured Cyber Solutions Ltd, licensed under the
[Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**:
no use, copying, modification or distribution is permitted without written permission from SCS.
Third-party components are listed in [NOTICE](NOTICE).

Licensing enquiries: https://securedcybersolutions.co.uk
