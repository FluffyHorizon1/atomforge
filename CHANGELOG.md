# Changelog

## [Unreleased]

### Added
- Ni–Al through the bridge. Recorded MACE runs (`results/NiAl_mace_*`, `quantum/run_nial_mace.sh`): fcc Ni, fcc Al, B2
  NiAl and L1₂ Ni₃Al relaxed with elastic constants and phonons for two models, and a reduced `search Ni Al` with its
  hull. A PySCF convergence study of the B2 NiAl formation energy at fixed geometry
  (`quantum/nial_pyscf_convergence.py`, `results/NiAl_pyscf_convergence.*`).
- `quantum/worker.py pyscf`: `--spin-polarised` (KUKS, reply field `moment`), `--multigrid`, `--ke-cutoff`,
  `--max-cycle`. `quantum/worker.py mace`: `--model` accepts a file, `--threads`; versions and the model's SHA-256 go to
  the log. The protocol is unchanged.
- CTest groups `external_pyscf` and `external_mace`, and CI jobs `pyscf backend` and `mace backend`.

- CI: `.github/workflows/qpuforge-crosscheck.yml` and `tools/qpuforge_crosscheck.sh` (qpuforge decision D3,
  2026-11-28). The job builds qpuforge's sidecar reader at a pinned ref (`qpuforge_ref` input, `QPUFORGE_REF` variable,
  default `main`), runs `search Cu Au --quick --materials` and requires qpuforge to read the sidecar: exit 0, all three
  quantities skipped as placeholders, nothing applied, and the SHA-256 in a `materials apply --out` copy when that
  qpuforge has the command. Non-blocking (`continue-on-error`); `ci.yml` is unchanged.
- README: the contract link points at qpuforge `main` instead of commit `f64238f`; the cross-check is documented.
- `ExternalPotential` (`include/atomforge/external.hpp`): a `Potential` backed by a persistent worker process, speaking
  line-delimited JSON (protocol 1, [docs/EXTERNAL_PROTOCOL.md](docs/EXTERNAL_PROTOCOL.md)). Energy (eV), forces (eV/Å) and
  the virial or stress in the ZhouEAM sign convention. Crashes, timeouts, malformed replies, wrong shapes and non-finite
  numbers raise `ExternalPotentialError` with the tail of the worker's stderr log.
- `quantum/worker.py`: reference worker with `lj` (pure Python), `pyscf` (periodic Kohn–Sham DFT) and `mace` backends;
  `--check BACKEND` exits 77 when a backend's packages are missing.
- `LennardJones` (`include/atomforge/lj.hpp`): native truncated-and-shifted Lennard-Jones, used to cross-check the bridge.
- `atomforge relax`: `--potential eam|lj[:eps,sigma,rc]|external:CMD`, `--volume`, `--cell`, `--fmax`,
  `--elastic`/`--no-elastic`, `--worker-timeout`, `--worker-log`. The default EAM output is unchanged.
- `atomforge search --potential eam|lj[:...]|external:CMD` and `binary_search(A, B, pot, opt)`: the whole search through
  any `Potential`. Also `--volume VA,VB`, `--prototypes`, `--no-ss`, `--no-mc`, and `SearchOptions::volume_per_atom`,
  `prototypes`, `ss_parents`, `ss_fractions`, `phonons`. `ExternalPotential` refuses concurrent `evaluate()` calls.
- Finite-displacement phonons (`include/atomforge/phonons.hpp`): force constants, acoustic sum rule, dynamical matrix,
  dynamical stability on the commensurate q grid, reported as `not assessed` when the supercell is too small (Γ only),
  and the sum-rule violation before correction. `relax --phonons`, `search --phonons` (verdicts `STABLE (... +
  dynamically stable)`, `UNSTABLE (imaginary phonons)`, `phonons not assessed`), and an additive `screen.phonons` block
  in `report.json` (`assessed`, `asr_violation_eV_A2`, `min_width_A`).
- CTest groups `external_search`, `phonons`, `phonons_lj`, `phonons_basis`, `phonons_asr` and `phonons_search`.
- CTest groups `lj`, `external`, `external_relax`, `external_errors` and `external_cli` (`tests/test_external.cpp`).

### Fixed
- `quantum/worker.py pyscf` refuses a spin-restricted run whose electron count over the k-mesh is odd. PySCF's smeared
  KRKS computed a charged cell there without an error (fcc Al on 3×3×3 k-points: 3.037 electrons per cell, energy
  0.33 eV too high); the k = 3 line of `results/pyscf_probe_log.txt` is such a run and is annotated.
