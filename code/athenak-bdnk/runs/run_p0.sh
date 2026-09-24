#!/usr/bin/env bash
# Phase-0 driver: the three 64^3 baseline runs, sequentially, 8 OpenMP threads (Windows/MinGW build).
# Usage (Git Bash): bash code/athenak-bdnk/runs/run_p0.sh [tags...]
D="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/code/athenak-bdnk"
EXE="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/ref-code/athenak/build-win/src/athena.exe"
export PATH="/c/msys64/ucrt64/bin:$PATH" OMP_NUM_THREADS=8
TAGS=("$@"); [ ${#TAGS[@]} -eq 0 ] && TAGS=(p0a_tov_cowling_K64_static p0a_tov_cowling_K64_kick p0b_tov_z4c_K64_static)
for tag in "${TAGS[@]}"; do
  R="$D/runs/$tag"; mkdir -p "$R"; cd "$R" || exit 1
  echo "=== $tag start $(date '+%F %T')" | tee -a "$D/runs/run_p0.log"
  T0=$(date +%s); "$EXE" -i "$D/inputs/$tag.athinput" > run.log 2>&1; EC=$?
  echo "=== $tag exit=$EC wall=$(( $(date +%s) - T0 )) s  $(grep zone-cycles run.log)" | tee -a "$D/runs/run_p0.log"
done
echo "P0 DONE" | tee -a "$D/runs/run_p0.log"
