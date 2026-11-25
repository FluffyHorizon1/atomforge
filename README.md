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
./build/atomforge validate          # 17 elements vs experiment
./build/atomforge relax L1_2 Cu Au  # one structure + elastic constants
./build/atomforge search Cu Au      # full binary search -> results/CuAu/ (report.json, POSCAR, extxyz)
```

Every candidate is written as **POSCAR** (VASP / Quantum ESPRESSO via ASE) and **extended XYZ** (OVITO, ASE).
Hull structures can go straight to DFT for verification.

## Energy model

The Zhou–Johnson–Wadley analytic EAM (Phys. Rev. B 69, 144113, 2004) covers 17 metals: Cu Ag Au Ni Pd Pt Al Pb Fe
Mo Ta W Cr Mg Co Ti Zr. Alloy cross terms come from the paper's mixing rule. Parameters were generated from the
LAMMPS `tools/eam_database` table (GPL-2.0), see `tools/gen_params.py`. Forces and virial are analytic and checked
against finite differences (errors 6e-9 and 4e-8). The `Potential` interface is abstract, so an ML potential or a
DFT backend can replace the EAM without touching the search, relaxation, hull or stability code.

## Validation

* All 17 elements: lattice constants within 0.75% of experiment, and the correct ground state (fcc/bcc/hcp) for all 17.
* **Cu–Au:** the search independently rediscovers the experimentally known ordered phases Cu₃Au (L1₂, a = 3.750 Å,
  matching experiment), CuAu (L1₀) and CuAu₃ (L1₂). All are on the hull and pass the mechanical and thermal checks.
  It also predicts extra hcp-based phases, which are likely artefacts of the potential.
* **Cu–Ag:** all formation energies positive, so immiscible, as observed.
* **Ni–Al: fails.** The generic mixing rule gives formation energies near zero or positive (B2 NiAl: +0.07 eV/atom),
  where experiment gives about −0.6 eV/atom.
  Strongly bonding intermetallics need quantum-level (DFT) energies. That is the next rung.

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

## Licence

Copyright (c) 2026 Secured Cyber Solutions Ltd. All rights reserved.

This is proprietary software of Secured Cyber Solutions Ltd, licensed under the
[Secured Cyber Solutions Proprietary Software Licence v1.0](LICENSE). It is **not open source**:
no use, copying, modification or distribution is permitted without written permission from SCS.
Third-party components are listed in [NOTICE](NOTICE).

Licensing enquiries: https://securedcybersolutions.co.uk
