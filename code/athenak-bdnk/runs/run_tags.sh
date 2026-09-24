#!/usr/bin/env bash
# Generic Phase-1 runner: run the named AthenaK tags in sequence, logging to runs/run_tags.log.
# Honours OMP_NUM_THREADS (default 4) so a run can share the box with another one.
#   bash run_tags.sh p1b_collapse_font_K64_L1
D="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/code/athenak-bdnk"
EXE="C:/Users/nicff/OneDrive/Desktop/Research/BDNK/ref-code/athenak/build-win/src/athena.exe"
export PATH="/c/msys64/ucrt64/bin:$PATH" OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
TAGS=("$@"); [ ${#TAGS[@]} -eq 0 ] && { echo "usage: run_tags.sh <tag> [tag...]"; exit 2; }
for tag in "${TAGS[@]}"; do
  R="$D/runs/$tag"; mkdir -p "$R"; cd "$R" || exit 1
  echo "=== $tag start $(date '+%F %T') threads=$OMP_NUM_THREADS" | tee -a "$D/runs/run_tags.log"
  T0=$(date +%s); "$EXE" -i "$D/inputs/$tag.athinput" > run.log 2>&1; EC=$?
  LAST=$(grep -c "^elapsed" run.log)
  echo "=== $tag exit=$EC wall=$(( $(date +%s) - T0 ))s cycles~$((LAST*100)) $(grep zone-cycles run.log)" | tee -a "$D/runs/run_tags.log"
done
echo "RUN_TAGS DONE" | tee -a "$D/runs/run_tags.log"
