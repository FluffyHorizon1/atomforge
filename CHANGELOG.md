# Changelog

## [Unreleased]

### Added
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
