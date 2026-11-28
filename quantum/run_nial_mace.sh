#!/bin/sh
# Recorded Ni-Al runs with a MACE foundation model through the external-potential bridge.
#   sh quantum/run_nial_mace.sh MODEL_FILE OUT.txt [THREADS]
# MODEL_FILE is a MACE model file, for example mace-mpa-0-medium.model from
# https://github.com/ACEsuit/mace-mp/releases/download/mace_mpa_0/mace-mpa-0-medium.model
# Needs ./build/atomforge and a python3 with mace-torch. MACE is a machine-learned potential, not DFT.
set -eu
MODEL=$1
OUT=$2
THREADS=${3:-4}
LOG=$(mktemp)
POT="external:python3 quantum/worker.py mace --model $MODEL --threads $THREADS"
{
  echo "# Ni-Al with MACE through the bridge. Date: $(date -u +%Y-%m-%dT%H:%MZ). atomforge commit $(git rev-parse --short HEAD)."
  echo "# The L1_2 line below is Ni3Al (three Ni, one Al); the CLI prints the element pair."
  echo "# Host: $(nproc) cores, $(uname -sr). $(python3 --version)."
  for spec in "A1 Ni" "A1 Al" "B2 Ni Al" "L1_2 Ni Al"; do
    echo
    echo "\$ ./build/atomforge relax $spec --fmax 0.001 --elastic --phonons --worker-hello-timeout 300 --potential \"$POT\""
    # shellcheck disable=SC2086
    ./build/atomforge relax $spec --fmax 0.001 --elastic --phonons --worker-hello-timeout 300 --worker-log "$LOG" --potential "$POT"
  done
  echo
  echo "# worker log (versions and model hash):"
  grep '^\[worker\]' "$LOG" | sort -u
} > "$OUT" 2>&1
# Formation energies per atom from the four relaxations above.
awk 'function e(   i){for(i=1;i<NF;i++) if($i=="E" && $(i+1)=="=") return $(i+2); return 0}
     /^A1 \(fcc\) Ni:/{ni=e()} /^A1 \(fcc\) Al:/{al=e()} /^B2 /{b2=e()} /^L1_2 /{l12=e()}
     END{printf "\n# formation energies from the lines above (eV/atom): B2 NiAl %.4f, L1_2 Ni3Al %.4f\n", b2-(ni+al)/2, l12-(3*ni+al)/4}' "$OUT" >> "$OUT"
rm -f "$LOG"
