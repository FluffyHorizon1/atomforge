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
- CTest groups `lj`, `external`, `external_relax`, `external_errors` and `external_cli` (`tests/test_external.cpp`).
