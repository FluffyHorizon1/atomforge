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
(72 s for 32 entries run one at a time on the 4-core cloud box used for this README, 27 s for `ctest -LE slow`; measured before the two sub-second `eam_*` entries below were added, 34 in all; see [Testing](#testing)).

Every candidate is written as **POSCAR** (VASP / Quantum ESPRESSO via ASE) and **extended XYZ** (OVITO, ASE).
Hull structures can go straight to DFT for verification.

## Materials sidecar (`forge.materials.v0`)

```bash
./build/atomforge search Cu Au --materials   # also writes results/CuAu/materials-contract-CuAu.json
```

Writes the interchange sidecar defined in qpuforge
[`docs/DATA_CONTRACT.md`](https://github.com/FluffyHorizon1/qpuforge/blob/main/docs/DATA_CONTRACT.md) (on qpuforge
`main`; first proposed on branch `docs/materials-data-contract`, commit `f64238f`); the writer is
`include/atomforge/materials.hpp`.

- **Real (computed by the EAM model):** `produced_by`, `system`, `structures_ref`, and an additive optional block
  `atomforge_computed` with elemental references and hull phases (energies, formation energies, cell edges and angles, Born
  elastic constants). These are classical-EAM numbers, not DFT.
- **Placeholders (not computed, `value: null`, `status: "placeholder"`):** `Q_dielectric`, `loss_tangent`,
  `junction_gap_ueV`. They are listed in the top-level `placeholders` array. Do not copy them into a chip JSON.

The contract is still a proposal, so the format may change.

**Cross-repo check (qpuforge reader).** `.github/workflows/qpuforge-crosscheck.yml` builds qpuforge's sidecar reader at a
pinned ref and runs `tools/qpuforge_crosscheck.sh`: `atomforge search Cu Au --quick --materials` writes the sidecar,
`qpuforge analyze designs/alpha1-v11/alpha1_v11.json --materials <it>` must exit 0 and skip all three quantities as
placeholders, and, when that qpuforge has `materials apply`, the saved chip copy must record the sidecar's SHA-256. It is
the one check that fails when either repository changes the format alone (decision D3, qpuforge
`docs/design/q16-materials-ingestion.md`). The job is **non-blocking** (`continue-on-error`); `ci.yml` stays the merge
gate. The qpuforge ref is the `qpuforge_ref` input of a manual run (Actions → "qpuforge reader cross-check" → Run
workflow), else the repository variable `QPUFORGE_REF`, else `main`; pin a full commit SHA to reproduce a run. A private
qpuforge needs a read-only token in the secret `QPUFORGE_READ_TOKEN`. Locally:

```bash
bash tools/qpuforge_crosscheck.sh build/atomforge ../qpuforge/build/qpuforge ../qpuforge
```

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
* **Ni–Al: fails with the EAM.** The generic mixing rule gives formation energies near zero or positive (B2 NiAl: +0.07
  eV/atom), where experiment gives about −0.6 eV/atom.
  Strongly bonding intermetallics need quantum-level energies. Through the external-potential bridge the sign and the size
  come out right: a MACE foundation model (a machine-learned potential fitted to PBE data, not DFT itself) gives B2 NiAl
  −0.66 eV/atom and L1₂ Ni₃Al −0.43 eV/atom, both on the hull and stable in every screen. See
  [Ni–Al through the bridge](#nial-through-the-bridge) for the runs, and for what the direct DFT calculation did and did not reach.
  The failure is pinned, not hidden: `eam_nial` asserts the current known-bad numbers (and fails if they change), and
  [docs/EAM_LIMITATIONS.md](docs/EAM_LIMITATIONS.md) lists it with the other EAM limitations and the QPU-stack materials the
  EAM cannot treat at all (no O, N, Si or Nb parameters).

## Testing

`ctest` runs `tests/test_atomforge.cpp`, `tests/test_external.cpp`, `tests/test_phonons.cpp` and `tests/test_materials_contract.cpp` (custom `CHECK` macros, one CTest entry per group, labels shown in brackets):

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
| `external_cli` [cli, external] | `atomforge relax --potential external:...` end to end against `--potential lj`, a reduced `search` through the worker (report names the backend, has `screen.phonons`), and the error exits |
| `external_json` [unit, external] | the protocol's JSON reader and number formatting (no python3 needed) |
| `external_locale` [unit, external] | JSON numbers and a whole handshake under a de_DE decimal-comma locale, generated with `localedef` into a temporary `LOCPATH`; skipped if that is not possible |
| `external_cli_args` [cli, external] | `relax` argument and failure exit codes: 1 for usage and argument errors (including a bad `--potential`), 2 for failures while running (no python3 needed) |
| `external_search` [physics, external, slow] | `binary_search` through the `lj` worker reproduces the native LJ search exactly; `report.json` names the backend and carries `screen.phonons`; full reduced pipeline (MC included) on the native LJ |
| `external_pyscf` [physics, external, pyscf, slow] | periodic DFT through the bridge on the 1-atom cell of fcc Al (PBE, gth-szv, 2×2×2 k-points): regression energy, zero force, isotropic virial, virial against a finite strain (1%), a hydrostatic FIRE relaxation driven by the DFT stress, the spin-polarised path against the restricted one, and the refusal of an odd electron count. About 4 minutes on 2 threads |
| `external_mace` [physics, external, mace, slow] | a MACE model through the bridge: forces and virial against finite differences on rattled, sheared B2 NiAl; lattice constants of Ni, Al and NiAl; B2 NiAl formation energy between −0.85 and −0.5 eV/atom |
| `phonons` [physics, phonons] | fcc Cu (EAM) stable, acoustic modes 0 at Γ, X-point frequencies vs experiment, supercell convergence; bcc Cu imaginary at N |
| `phonons_basis` [unit, phonons] | 4-atom cubic fcc vs 1-atom primitive at q + G; L1₂ Cu₃Au (unequal masses) vs its doubled cell with a reordered basis |
| `phonons_asr` [unit, phonons] | the sum-rule violation before correction: 1e-13 for the EAM, the finite-difference value for a synthetic egg box |
| `phonons_search` [physics, phonons, slow] | every screen verdict with phonons on synthetic candidates; a reduced Cu–Au search reaches `STABLE (... + dynamically stable)` |
| `phonons_lj` [unit, phonons] | Lennard-Jones fcc phonons vs the analytic dispersion (commensurate and interpolated q); force constants through the bridge equal the native ones |
| `md_invariance` [physics, symmetry] | `run_md` averages and final positions unchanged when the start structure's atoms are shifted by lattice vectors (including shifts that make `run_md` wrap) |
| `eam_nial` [physics, eam] | the Ni–Al known failure, pinned: B2 NiAl Ef = +0.0676 eV/atom (experiment about −0.6), B2 above L1₁/L1₀, no prototype bound by more than 50 meV, and `results/NiAl/report.json` equal to a fresh computation. Fails, saying what to update, if the failure changes in either direction |
| `eam_scope` [unit, eam] | the Zhou04 table is the 17 metals; O, N, Si and Nb have no parameters and `ZhouEAM` / `binary_search` refuse them, so AlOx, Ta/Nb oxides, TiN and Si cannot be set up |
| `structure_misc` [unit] | composition, fraction and supercell bookkeeping; errors for unknown elements, prototypes and parent lattices |
| `io` [unit] | writers throw `cannot write` on unwritable paths; POSCAR, extended XYZ and `report.json` writers checked against the in-memory structure, `report.json` parsed as strict JSON (escaping, empty and real search results) |
| `cli` [unit, cli, slow] | the `atomforge` executable: unknown, missing, non-numeric and out-of-range `search` arguments fail fast; usage and error exit codes, `elements`, `relax` (stable and Born-unstable cases), and a reduced `search` whose `report.json` and structure files are cross-checked. Runs as `test_atomforge cli <path-to-atomforge>`; the executable checks are skipped without the path |

It does not cover the Monte Carlo ordering search beyond the reduced hull and CLI runs, or MD thermodynamics beyond the
shift-invariance check.
The `external*` groups need `python3` and are reported as skipped (exit 77) without it. `external_pyscf` is also skipped
without the `pyscf` module, and `external_mace` without `mace-torch` or without a model file named by
`ATOMFORGE_MACE_MODEL`. CI runs these two in jobs of their own (`pyscf backend`, `mace backend`), which install the
packages (`pip install pyscf==2.15.1`; torch 2.15.0 CPU and `mace-torch==0.3.17`), fetch `mace-mpa-0-medium.model` from
the ACEsuit GitHub release, verify its SHA-256, cache it, and fail if the test is skipped.

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
| `pyscf` | periodic Kohn–Sham DFT with the probe's settings; analytic forces and stress; `--spin-polarised`, `--multigrid`, `--ke-cutoff` | numpy, pyscf |
| `mace` | a MACE foundation model through ASE (`--model` name or file); a machine-learned potential, not DFT | numpy, ase, mace-torch |

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
start-up and handshake (default 60 s; MACE downloads its model then; `relax` only). The `pyscf` backend is spin-restricted
unless `--spin-polarised` is given, so Fe, Co and Ni come out non-magnetic by default. `relax` exits with 0 on success, 1 on usage and argument errors (unknown option or prototype, bad
value, extra arguments), and 2 on failures while running (unsupported element, a worker that fails), as `search` does.

In the tests, the bridge with the `lj` worker reproduces the native LJ energy, forces and virial to better than 1e-9 (the
differences were exactly zero in our runs) on rattled fcc, triclinic and unwrapped cells, and FIRE through the bridge
follows the native path step for step (79 steps, full cell).

The `pyscf` backend was checked against finite differences of its own energy (`quantum/check_pyscf_fd.py`, fcc Al,
gth-szv, 2×2×2 k-points, log in `results/pyscf_bridge_fd_log.txt`): the trace of the virial was -0.6362 eV analytic and
-0.6367 eV by finite strain, and a force on a displaced atom -0.11778 eV/Å analytic and -0.11777 eV/Å by finite
displacement. With Fermi smearing the worker returns the free energy E - σS, because that is what PySCF's forces and
stress are derivatives of; returning E instead put the virial off by a factor of 1.9. One evaluation took 26 s for the
1-atom cell and 131 s for the 4-atom cell on 2 threads. CI runs both optional backends in jobs of their own (see
[Testing](#testing)).

### Ni–Al through the bridge

<a id="nial-through-the-bridge"></a>
Two different things were run through `ExternalPotential` for Ni–Al, and they must not be confused. **MACE** is a
machine-learned interatomic potential trained on PBE calculations; it is cheap enough for relaxations, phonons and a
search, and its numbers are predictions of a fitted model. **PySCF** is the DFT calculation itself (PBE, GTH
pseudopotentials); it is exact for its functional only when converged in basis set and k-points.

**MACE** (`quantum/run_nial_mace.sh`; records `results/NiAl_mace_mpa0_relax.txt` and `results/NiAl_mace_mp0_relax.txt`;
mace-torch 0.3.17, torch 2.15.0, ASE 3.30.0, CPU, 4 threads, 2026-11-28). Each line is
`atomforge relax <prototype> ... --fmax 0.001 --elastic --phonons --potential external:"python3 quantum/worker.py mace --model FILE --threads 4"`.
The two models are `mace-mpa-0-medium.model` (SHA-256 `75428afe…4fb638`) and MACE-MP-0 medium,
`2023-12-03-mace-128-L1_epoch-199.model` (SHA-256 `01bfe221…d35845`), both from the ACEsuit `mace-mp` GitHub releases; the
full hashes are in the records.

| | mace-mpa-0 | MACE-MP-0 | PBE, published | experiment |
|---|---|---|---|---|
| fcc Ni, a (Å) | 3.521 | 3.510 | | 3.524 |
| fcc Al, a (Å) | 4.039 | 4.061 | | 4.050 |
| B2 NiAl, a (Å) | 2.885 | 2.883 | 2.898 | 2.887 |
| L1₂ Ni₃Al, a (Å) | 3.561 | 3.556 | 3.582 | 3.57 |
| B2 NiAl, formation energy (eV/atom) | −0.659 | −0.692 | −0.677 | −0.685 (−0.60 to −0.70) |
| L1₂ Ni₃Al, formation energy (eV/atom) | −0.426 | −0.476 | −0.441 | −0.428 (−0.39 to −0.43) |
| reference energies Ni, Al (eV/atom) | −5.787, −3.740 | −5.732, −3.728 | | |

The published PBE values are from Lechermann, Fähnle and Sanchez, Intermetallics 13, 1096 (2005) (arXiv:cond-mat/0410285,
tables 1 and 2: mixed-basis pseudopotential calculations, Ni₃Al ferromagnetic). The experimental formation energies are
enthalpies of formation measured by aluminium-solution calorimetry at 298–300 K: B2 AlNi −66.1 ± 1.3 and L1₂ AlNi₃
−41.3 ± 1.3 kJ/mol of atoms (direct-reaction calorimetry: −65.7 ± 1.1 and −39.4 ± 0.6). They are taken from tables 7.1.2,
7.1.3, 7.3.2 and 7.3.3 of K. Rzyman, Z. Moser and J.-C. Gachon, "Calorimetric studies of the intermetallic compounds from
Al-Ti and Al-Ni systems", monograph chapter, pp. 93–146, read as the PDF in the RCIN repository of the Polish Academy of
Sciences (https://rcin.org.pl/ippt/Content/67610/PDF/WA727_25873_55602_Rzyman-Calorimetric.pdf; the enclosing volume
could not be confirmed from the PDF; its newest references are from 2004). The AlNi₃ solution values, −41.3 to −42.3
kJ/mol between 300 and 1123 K, were first published by Rzyman et al., J. Phase Equilibria 17, 173 (1996). The earlier
measurements the chapter tabulates (tables 7.1.4 and 7.3.4) span −58.3 to −67.4 kJ/mol of atoms for AlNi and −37.6 to
−40.6 for AlNi₃, the ranges in brackets above; that spread across methods, not the ± of one method, is the experimental
uncertainty. Room-temperature enthalpies are compared here with 0 K energies without zero-point motion.
The experimental lattice constants of Ni and Al are those of `atomforge validate`. The two models differ by 0.03 eV/atom
for NiAl and 0.05 eV/atom for Ni₃Al, and both lie within 0.04 eV/atom of the published PBE values and of the
calorimetric values. That spread, about 0.05 eV/atom, is the honest error bar on the MACE numbers; the Zhou04 EAM is off
by 0.7 eV/atom.

All four phases are Born-stable and have no imaginary phonon on the commensurate grid with both models (mace-mpa-0, in GPa:
Ni C11 = 290, C12 = 178, C44 = 136; Al 115, 65, 29; NiAl 214, 149, 117; Ni₃Al 199, 136, 124; phonons on 3×3×3 supercells
of the 4-atom cells and 4×4×4 of the B2 cell, highest mode 8.0, 7.9, 10.4 and 9.3 THz). The four relaxations with elastic
constants and phonons take about 2 minutes per model.

The search (`results/NiAl_mace_run.txt`, `results/NiAl_mace/`):

```bash
./build/atomforge search Ni Al --out results/NiAl_mace --no-ss --no-mc --md-steps 300 --phonons \
    --potential external:"python3 quantum/worker.py mace --model mace-mpa-0-medium.model --threads 4"
```

It made 1540 worker evaluations in 24 minutes. The hull has L1₂ Ni₃Al (−426 meV/atom) and NiAl (−659 meV/atom), both
`STABLE (on hull, mechanically + thermally + dynamically stable)`. Its limits: it is reduced to the ordered prototypes,
because the 27 solid-solution relaxations on 48-site cells had not finished after 70 minutes on this machine and the Monte
Carlo stage would need about 4e5 evaluations; the MD screen is 300 steps; the L1₀ prototype relaxes into B2 (c/a → the
bcc value, same energy), so the hull point is listed as `L1_0 NiAl`; and the prototype list has none of the Al-rich Ni–Al
structures (D0₁₁ NiAl₃, D5₁₃ Ni₂Al₃), so nothing can be concluded for x(Al) > 0.5.

### The whole search through any potential

`atomforge search` takes the same `--potential` option, and every stage then uses that potential: the elemental
references, prototypes, solid solutions, Monte Carlo ordering, relaxation, hull, Born and MD screen, and the optional phonons.
All evaluations are made one at a time on a single `Potential`, with no threading over candidates. One worker process
therefore serves the whole search, and `ExternalPotential` refuses a concurrent call from a second thread rather than
interleaving requests. `report.json` names the backend in `"potential"` (for example `"external:lj(eps=0.0104,sigma=3.4,rc=8.5)"`).

```bash
./build/atomforge search Ar Kr --potential external:"python3 quantum/worker.py lj" --volume 38,38 \
    --prototypes L1_2,B2 --no-ss --no-mc --md-steps 50 --phonons --phonon-supercell 1
```

That run made 785 worker calls in 14.8 s, and the native `--potential lj` search gives identical energies and volumes
(test `external_search`). `--volume VA,VB` sets the starting volumes per atom for elements outside the Zhou04 table.
`--prototypes`, `--no-ss` and `--no-mc` cut the candidate list. That matters at DFT cost, because a full search makes
several thousand evaluations (2347 for an Ar–Kr search with all prototypes, one fcc solid solution and one fcc MC run, LJ cut off at 5 Å, 18 candidates).

For QPU materials (Al/AlOx junction barriers, Nb and Ta surface oxides, TiN, the Si substrate) EAM is not enough. (Which of these the EAM can treat at all: [docs/EAM_LIMITATIONS.md](docs/EAM_LIMITATIONS.md).)
Oxides, nitrides and covalent solids need DFT or reactive potentials.

## Phonons

`include/atomforge/phonons.hpp` computes harmonic phonons by finite displacements. Each atom of the cell is displaced
by ±0.01 Å in a supercell whose faces are at least twice the cutoff apart, which gives the force constants by central
differences. The acoustic sum rule is imposed on the self term, and the dynamical matrix is built with minimum-image
weights. Frequencies at the q-points commensurate with the supercell are exact for those force constants, so the
stability check samples that grid, not a band path. A structure is dynamically stable when no mode is below −0.1 THz on
that grid. The check counts only when the supercell has at least two cells along every lattice vector and faces at least
one cutoff apart. Otherwise it is reported as `not assessed`: a Γ-only (1×1×1) check always passes, because the sum rule
zeroes the acoustic modes there, so bcc Cu at 1×1×1 is `not assessed`, not stable. Imaginary modes on any grid still mean
`dynamically UNSTABLE`. The size of the sum-rule correction, max |Σⱼ Φ|, is reported as `asr_violation`. It is 1e-13
eV/Å² for the EAM, and it shows egg-box errors of a DFT grid: a synthetic 0.05 eV/Å² egg box is reported as 0.0498, the
finite-difference value. `atomforge relax ... --phonons` prints the check. `search --phonons` adds it to the stability
screen (`STABLE (on hull, mechanically + thermally + dynamically stable)`, `UNSTABLE (imaginary phonons)`, or `phonons
not assessed` in the verdict) and to `report.json` as `screen.phonons`.

| | this code (Zhou04 EAM) | experiment, Cu at 80 K (Svensson, Brockhouse and Rowe, Phys. Rev. 155, 619 (1967)) |
|---|---|---|
| fcc Cu, X transverse | 5.11 THz | 5.10 THz |
| fcc Cu, X longitudinal | 8.34 THz | 7.25 THz |
| fcc Cu, highest mode (L point) | 8.39 THz | |
| fcc Cu, acoustic modes at Γ | 1e-7 THz | 0 |
| bcc Cu | imaginary, −1.12 THz at N | (bcc Cu is not stable) |

The longitudinal modes come out 15% stiff with this potential (16% for the highest mode); the transverse X mode is within
0.2%. That is the Zhou04 potential's stiffness, not a supercell effect: at a q that is commensurate with neither supercell,
(0.31, 0.12, 0.05) 2π/a, the 6×6×6 and 8×8×8 supercells agree to 4.5e-4 THz. (At X and L they agree exactly by
construction.) The `phonons` test therefore allows 20% for the longitudinal X mode and the highest frequency, and 5% for
the transverse X mode. Multi-atom cells are checked by folding: the 4-atom cubic fcc cell reproduces the primitive cell's
frequencies at q + G to 1e-7 THz, and L1₂ Cu₃Au (masses 63.5 and 197.0) matches its doubled 8-atom cell to 3e-10 THz.
For the native Lennard-Jones fcc lattice, the finite-displacement
frequencies agree with the analytic lattice-sum dispersion to 2e-4 THz, and the force constants through the `lj` worker
are identical to the native ones. The cost is 6 evaluations per atom of the cell, on the supercell, which is small for EAM
but dominant for DFT.

## Limitations

The Monte Carlo search runs on 48-site cells. Dynamical stability is checked by MD heat-and-quench, and by phonons when
`--phonons` is given; the phonon check samples only the q-points commensurate with its supercell, so a small supercell
(`--phonon-supercell 1` is Γ only, and is reported as `not assessed`) can miss an instability.
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
