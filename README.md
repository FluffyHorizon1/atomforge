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
(about 4 s, see [Testing](#testing)).

Every candidate is written as **POSCAR** (VASP / Quantum ESPRESSO via ASE) and **extended XYZ** (OVITO, ASE).
Hull structures can go straight to DFT for verification.

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

`ctest` runs `tests/test_atomforge.cpp` (custom `CHECK` macros, one CTest entry per group, labels shown in brackets):

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
| `unwrapped` [known-bug, xfail] | **expected to fail**: see Limitations |

The suite takes about 4 s. It does not cover the Monte Carlo ordering search beyond the reduced hull run, the MD
integrator, file output (`io.hpp`) or the `search`/`relax` command-line paths.

## Toward the quantum level

`quantum/pyscf_dft_probe.py` runs periodic Kohn–Sham DFT (PBE, GTH pseudopotentials, k-points) with PySCF.
It returns energy, analytic forces and the stress tensor, which is everything the FIRE relaxer needs.
On a 2-core cloud box, fcc Al with 27 k-points took about 60 s per evaluation (`results/pyscf_probe_log.txt`).
The planned bridge is an `ExternalPotential` that streams structures to a persistent Python worker
(PySCF / GPAW / Quantum ESPRESSO / MACE), so the same search runs at DFT accuracy on your own hardware.

For QPU materials (Al/AlOx junction barriers, Nb and Ta surface oxides, TiN, the Si substrate) EAM is not enough.
Oxides, nitrides and covalent solids need DFT or reactive potentials.

## Limitations

The Monte Carlo search runs on 48-site cells. Dynamical stability is checked by MD heat-and-quench, not phonons.
There is no magnetism and no charge transfer (EAM). Alloy energies are qualitative.

Known bug: the neighbour list assumes atoms are in (or very close to) the home cell. Atoms displaced by two or more
cell vectors, or by one vector in a small triclinic cell, silently lose neighbours and the energy is wrong. `relax` and
the search keep atoms close to the cell, but `run_md` does not wrap positions, so long or hot MD runs on small cells
can be affected. Tracked by the expected-fail `unwrapped` test.

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This is proprietary software of Secured Cyber Solutions Ltd, licensed under the
[Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**:
no use, copying, modification or distribution is permitted without written permission from SCS.
Third-party components are listed in [NOTICE](NOTICE).

Licensing enquiries: https://securedcybersolutions.co.uk
