#!/usr/bin/env bash
# Phase-1 driver: the ideal-fluid collapse / migration runs (octant, Z4c), sequentially.
D="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/code/athenak-bdnk"
EXE="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/ref-code/athenak/build-win/src/athena.exe"
export PATH="/c/msys64/ucrt64/bin:$PATH" OMP_NUM_THREADS=8
TAGS=("$@"); [ ${#TAGS[@]} -eq 0 ] && TAGS=(p1a_collapse_eps0003_K64 p1b_collapse_font_K64 p1b_migration_font_K64)
for tag in "${TAGS[@]}"; do
  R="$D/runs/$tag"; mkdir -p "$R"; cd "$R" || exit 1
  echo "=== $tag start $(date '+%F %T')" | tee -a "$D/runs/run_p1.log"
  T0=$(date +%s); "$EXE" -i "$D/inputs/$tag.athinput" > run.log 2>&1; EC=$?
  LAST=$(grep -c "^elapsed" run.log)
  echo "=== $tag exit=$EC wall=$(( $(date +%s) - T0 ))s cycles~$((LAST*100)) $(grep zone-cycles run.log)" | tee -a "$D/runs/run_p1.log"
done
echo "P1 DONE" | tee -a "$D/runs/run_p1.log"
