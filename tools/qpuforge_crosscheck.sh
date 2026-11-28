#!/usr/bin/env bash
# Cross-repo check of the forge.materials.v0 contract (qpuforge docs/DATA_CONTRACT.md, decision D3 of qpuforge
# docs/design/q16-materials-ingestion.md): atomforge writes a sidecar, qpuforge's own reader reads it.
#
#   tools/qpuforge_crosscheck.sh <atomforge-binary> <qpuforge-binary> <qpuforge-source-dir> [work-dir]
#
# 1. `atomforge search Cu Au --quick --materials` writes materials-contract-CuAu.json (three null placeholders today).
# 2. `qpuforge analyze <Alpha-1 v1.1> --materials <that file>` must exit 0, report every quantity as a skipped
#    placeholder and apply nothing (atomforge computes none of the three quantities yet).
# 3. When the qpuforge build has `materials apply` (branch feat/materials-apply-out and later), the sidecar is also
#    saved into a chip copy, whose materials_source must carry the sidecar's SHA-256.
# Exit 0 when every check passes, 1 otherwise. Needs bash, grep and sha256sum.
set -u
AF=${1:?atomforge binary}; QF=${2:?qpuforge binary}; QSRC=${3:?qpuforge source dir}; WORK=${4:-qpuforge-crosscheck}
CHIP="$QSRC/designs/alpha1-v11/alpha1_v11.json"
mkdir -p "$WORK"
fails=0
check() { if eval "$2"; then echo "PASS  $1"; else echo "FAIL  $1"; fails=$((fails + 1)); fi; }

"$AF" search Cu Au --quick --materials --out "$WORK/af" > "$WORK/atomforge.log" 2>&1
check "atomforge search Cu Au --quick --materials exits 0" "[ $? -eq 0 ]"
SIDECAR="$WORK/af/materials-contract-CuAu.json"
check "atomforge wrote $SIDECAR" "[ -s '$SIDECAR' ]"
check "the sidecar says forge.materials.v0" "grep -q '\"schema\": \"forge.materials.v0\"' '$SIDECAR'"

"$QF" analyze "$CHIP" --materials "$SIDECAR" > "$WORK/analyze.log" 2>&1
check "qpuforge analyze --materials <atomforge sidecar> exits 0" "[ $? -eq 0 ]"
for q in Q_dielectric loss_tangent junction_gap_ueV; do
  check "qpuforge skips $q as a placeholder" "grep -Eq '^  skipped +$q +placeholder' '$WORK/analyze.log'"
done
check "qpuforge applies nothing from a placeholder sidecar" "! grep -Eq '^  applied ' '$WORK/analyze.log'"
check "qpuforge reports no error" "! grep -q 'error:' '$WORK/analyze.log'"

if "$QF" materials --help > /dev/null 2>&1; then
  "$QF" materials apply "$CHIP" "$SIDECAR" --out "$WORK/chip-applied.json" > "$WORK/apply.log" 2>&1
  check "qpuforge materials apply --out exits 0" "[ $? -eq 0 ]"
  sha=$(sha256sum "$SIDECAR" | cut -d' ' -f1)
  check "materials_source records the sidecar SHA-256 ($sha)" "grep -q '\"sidecar_sha256\": \"$sha\"' '$WORK/chip-applied.json'"
else
  echo "SKIP  qpuforge materials apply (not in this qpuforge build)"
fi

echo "--- qpuforge analyze output (materials lines)"
grep -E '^(materials sidecar|  (applied|skipped) )' "$WORK/analyze.log" || cat "$WORK/analyze.log"
if [ "$fails" -ne 0 ]; then echo "$fails check(s) failed"; exit 1; fi
echo "all cross-repo checks passed"
