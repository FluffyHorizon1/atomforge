#!/bin/sh
# Cost of one full evaluation (energy, forces and stress) of B2 NiAl through the bridge with the pyscf backend.
#   sh quantum/run_nial_pyscf_bridge_cost.sh OUT.txt [THREADS]
# Each setting runs `atomforge relax B2 Ni Al --cell fixed` at the experimental lattice constant, 2.887 A
# (12.0315 A^3/atom). The atoms sit on symmetry positions, so the forces vanish and FIRE stops after the first
# evaluation. The worker log gives the wall time of the SCF, the forces and the stress separately. PySCF computes
# forces and stress on the full k-mesh (no symmetry reduction), so these runs measure what a relaxation through the
# bridge costs per step. Needs ./build/atomforge and a python3 with pyscf.
set -eu
OUT=$1
THREADS=${2:-4}
LOG=$(mktemp)
{
  echo "# B2 NiAl, one evaluation through the bridge (pyscf backend). Date: $(date -u +%Y-%m-%dT%H:%MZ). atomforge commit $(git rev-parse --short HEAD)."
  echo "# Host: $(nproc) cores ($(grep -m1 'model name' /proc/cpuinfo | sed 's/.*: //')), $(uname -sr). $(python3 --version)."
  for setting in "gth-szv-molopt-sr 2" "gth-szv-molopt-sr 4" "gth-dzvp-molopt-sr 2"; do
    # shellcheck disable=SC2086
    set -- $setting
    POT="external:python3 quantum/worker.py pyscf --basis $1 --kpts $2,$2,$2 --smearing 0.005 --ke-cutoff 400 --multigrid --threads $THREADS"
    echo
    echo "\$ ./build/atomforge relax B2 Ni Al --volume 12.0315 --cell fixed --no-elastic --worker-timeout inf --potential \"$POT\""
    : > "$LOG"
    ./build/atomforge relax B2 Ni Al --volume 12.0315 --cell fixed --no-elastic --worker-timeout inf --worker-log "$LOG" --potential "$POT"
    grep '^\[worker\] \(pyscf \|E = \)' "$LOG"
  done
} > "$OUT" 2>&1
rm -f "$LOG"
