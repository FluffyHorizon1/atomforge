# What the built-in EAM cannot do

atomforge's built-in energy model is the Zhou–Johnson–Wadley (2004) embedded-atom method (`include/atomforge/potential.hpp`,
parameters in `include/atomforge/zhou04_params.hpp`, generated from the LAMMPS EAM database). It is fast and it gets
elemental metals and weakly-bonded alloys qualitatively right. This page states where it does not, with the numbers, so
that the DFT / `ExternalPotential` work has an explicit baseline to beat rather than a footnote to find. Everything here is
software screening; nothing on it is a fabrication claim.

The two groups in `tests/test_eam_limits.cpp` (`eam_nial`, `eam_scope`; `ctest -L eam`) assert the statements below, so
they stay true: they pass while the limitation is unchanged and fail, with a message listing what to update, if it moves.

## 1. Ni–Al is a known failure

`README.md` ("Validation") says the search fails for Ni–Al. The committed record is `results/NiAl/` (`report.json`, run log
`results/NiAl_run.txt`, potential `Zhou04-EAM`). The numbers, in eV/atom, relative to the lowest elemental reference
(Ni fcc −4.4500 eV, Al fcc −3.5800 eV):

| structure | x(Al) | formation energy Ef | above hull | note |
|---|---|---|---|---|
| B2 NiAl | 0.50 | **+0.0676** | +0.0746 (full search) | experiment: about −0.6 (README) |
| L1₁ NiAl | 0.50 | +0.0121 | | lowest equiatomic ordered prototype |
| L1₀ NiAl | 0.50 | +0.0162 | | |
| L1₂ Ni₃Al | 0.25 | +0.0334 | | |
| L1₂ NiAl₃ | 0.75 | −0.0092 | | the most strongly bound ordered prototype |
| MC-bcc NiAl₃ (Monte Carlo ordering) | 0.75 | −0.0104 | 0 (on hull) | the deepest point in the whole search |

What this shows:

- **Wrong sign and size.** B2 NiAl, the structure NiAl is observed in, is predicted *unbound* by +0.07 eV/atom against about
  −0.6 experimentally: an error of roughly 0.7 eV/atom. The deepest hull point in the whole search is −0.010 eV/atom.
- **Wrong ordering.** Among the equiatomic ordered phases the EAM puts B2 about 55 meV *above* L1₁ and 51 meV above L1₀.
- **Spurious "stable" verdicts.** The stability screen is relative to the model's own hull. In `results/NiAl_run.txt` the
  only two candidates it calls `STABLE (on hull, mechanically + thermally …)` are bcc-ordered Al-rich cells (NiAl₃ at
  −10.4 meV and NiAl₇ at −9.8 meV/atom). Those are artefacts of an almost flat energy landscape, not predictions; do not
  read a Ni–Al `STABLE` from this potential as evidence a compound exists.
- The README gives the reason: strongly bonding intermetallics need quantum-level energies (DFT). The EAM's alloy energies
  come from combining single-element fits with a generic mixing rule, which has no term for the charge transfer and
  directional d-band bonding that bind Ni–Al.

`eam_nial` pins B2 Ef (+0.0676 ± 0.002), the lowest prototype (−0.0092, L1₂ NiAl₃), the equiatomic ordering and the
prototype-only hull, states the failure as inequalities (Ef > 0, more than 0.5 eV from experiment, no prototype binding by
more than 50 meV), and checks `results/NiAl/report.json` against a fresh computation. It runs in well under a second
(prototypes only: no Monte Carlo, no screening), so it is in the default `ctest`.

**When a DFT or ML backend makes Ni–Al work**, this test is the one to invert: B2 NiAl Ef should land near the experimental
value, and the pin becomes an acceptance test for that backend (planned milestone: DFT / ExternalPotential). Until then a
failure of `eam_nial` means the EAM failure *changed*, and the message says which README paragraph, which page (this one) and
which result files to update together.

That has now happened for a machine-learned backend. Through `ExternalPotential`, a MACE foundation model gives B2 NiAl
−0.66 eV/atom and puts NiAl and L1₂ Ni₃Al on the hull (README, "Ni–Al through the bridge"; records under `results/NiAl_mace*`),
and the test `external_mace` asserts the B2 formation energy for that backend. `eam_nial` stays as it is: it pins the
built-in EAM, which has not changed.

## 2. Other limitations of the model and the search

From `README.md` ("Limitations", "Validation") and the committed runs:

- **No magnetism and no charge transfer.** Fe, Co and Ni are treated as non-magnetic metals; ionic and mixed ionic–covalent
  bonding is out of reach. Alloy energies are qualitative.
- **Spurious phases.** The Cu–Au search predicts extra hcp-based phases that are likely artefacts of the potential
  (README, "Cu–Au"). Treat the on-hull phases that experiment does not know as suspect until checked with DFT.
- **Stiff longitudinal phonons.** For fcc Cu the EAM's X-point longitudinal mode is about 15% above experiment (8.34 vs 7.25
  THz); the `phonons` test allows 20% for it and 5% for the transverse mode (`README.md`, Phonons).
- **Small cells.** The Monte Carlo ordering search runs on 48-site cells, so it can only see orderings commensurate with them.
- **Stability screen.** Dynamical stability is the MD heat-and-quench check, and the finite-displacement phonon check only when
  `--phonons` is given and only at q-points commensurate with the phonon supercell (`--phonon-supercell 1` is Γ only and
  is reported as `not assessed`).

These are requirements on a DFT (or ML) backend, not caveats to carry silently: charge transfer for oxides and
intermetallics, and a route to phonon stability, are what the `ExternalPotential` bridge (`docs/EXTERNAL_PROTOCOL.md`) is for.

## 3. Which QPU-stack materials the EAM cannot treat

The Zhou04 table covers **17 metals only**: Cu, Ag, Au, Ni, Pd, Pt, Al, Pb, Fe, Mo, Ta, W, Mg, Co, Ti, Zr, Cr. There is no
oxygen, nitrogen, silicon or niobium. `ZhouEAM` refuses them (`ZhouEAM: no parameters for O (supported: …)`), and the
search cannot be started for a pair containing one; it never falls back to a guess. For the materials stack named in
qpuforge's `ROADMAP.md` Phase 0 and in the README ("Toward the quantum level"):

| QPU-stack material | needs | EAM can treat it? |
|---|---|---|
| Al metal | Al | **yes**, as an elemental fcc metal only (all 17 elements pass `validate`) |
| Al/AlOx junction barrier (amorphous or crystalline oxide) | Al, O | **no**: no O |
| Ta metal | Ta | **yes**, as an elemental bcc metal only |
| Ta surface oxide | Ta, O | **no**: no O |
| Nb metal | Nb | **no**: Nb is not in the table |
| Nb surface oxide | Nb, O | **no**: no Nb, no O |
| TiN | Ti, N | **no**: no N (Ti alone is in the table) |
| Si substrate | Si | **no**: no Si (covalent; needs DFT or a reactive potential) |

Even where an element is in the table, "treatable" means *bulk elemental metal structure and metal–metal alloy energies at
the qualitative level of section 2*: not oxides, nitrides, interfaces, surfaces, dielectric loss, or superconducting
properties. Oxides, nitrides and covalent solids need DFT or reactive potentials (README, "Toward the quantum level"), which
is the planned "Materials loop" work after `ExternalPotential`. `eam_scope` pins the 17-element table and that O, N, Si and
Nb are absent, and that each compound row above cannot be set up, so a table change cannot silently make this page wrong.

Note on Nb: the backlog text for this item lists O, N and Si as missing; Nb is missing from the table too (checked against
`zhou04_params.hpp`), which matters because Nb is one of the stack's superconductors.
