#!/usr/bin/env bash
# ralph_stop_guard.sh — progress-aware Stop-hook guard (BDNK consumer).
#
# Consults phys-agentic-loop/_common/loop_gate.py. Blocks the stop (forces the
# Ralph loop to continue) ONLY while the gate verdict is `continue` — i.e. the
# loop is making verified progress within its iteration / wall-clock budget. On
# any `halt:*` verdict it stays silent and lets the agent stop; loop_gate writes
# .claude/HUMAN_REVIEW_REQUIRED.md for halts needing human attention. Falls back
# to the legacy `iteration < max_iterations` rule if loop_gate is unavailable.

set -u
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || printf '%s' "$(cd "$(dirname "$0")/.." && pwd)")"
STATE="$REPO_ROOT/.claude/ralph-loop.local.md"
[ -f "$STATE" ] || exit 0

emit_block() {
    local reason="$1" encoded
    if command -v python3 >/dev/null 2>&1; then
        encoded="$(printf '%s' "$reason" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))')"
    else
        encoded="\"${reason//\"/\'}\""
    fi
    printf '{"decision":"block","reason":%s}\n' "$encoded"
}

legacy_guard() {
    local active iter max
    active=$(awk '/^---$/{n++; next} n==1 && /^active:/ {print $2}' "$STATE" | head -1)
    iter=$(awk '/^---$/{n++; next} n==1 && /^iteration:/ {print $2}' "$STATE" | head -1)
    max=$(awk '/^---$/{n++; next} n==1 && /^max_iterations:/ {print $2}' "$STATE" | head -1)
    [ "$active" = "true" ] || exit 0
    [ -n "$iter" ] && [ -n "$max" ] || exit 0
    case "$iter$max" in *[!0-9]*) exit 0 ;; esac
    if [ "$iter" -lt "$max" ]; then
        emit_block "Ralph loop active — iteration ${iter}/${max} (legacy guard; loop_gate.py unavailable). Advance the counter in .claude/ralph-loop.local.md, ship the next verified reproduction, then continue."
    fi
    exit 0
}

GATE=""
for cand in \
    "$REPO_ROOT/phys-agentic-loop/_common/loop_gate.py" \
    "$REPO_ROOT/_common/loop_gate.py" \
    "${PHYS_AGENTIC_LOOP:-}/_common/loop_gate.py"; do
    if [ -n "$cand" ] && [ -f "$cand" ]; then GATE="$cand"; break; fi
done

if [ -z "$GATE" ] || ! command -v python3 >/dev/null 2>&1; then legacy_guard; fi

out="$(python3 "$GATE" decide --repo-root "$REPO_ROOT" --state-file "$STATE" --write-gate 2>/dev/null)"
decision="$(printf '%s' "$out" | python3 -c 'import json,sys
try: print(json.load(sys.stdin)["decision"])
except Exception: pass' 2>/dev/null)"
[ -n "$decision" ] || legacy_guard

case "$decision" in
    continue)
        reason="$(printf '%s' "$out" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("reason",""))' 2>/dev/null)"
        emit_block "Ralph loop progress gate: ${reason}. Advance the iteration counter in .claude/ralph-loop.local.md, reproduce the next DAG node (append knowledge/error rows + a VLM-checked figure, commit per substage), then continue. The gate halts automatically if progress stalls."
        ;;
    halt:*)
        : ;;  # graceful stop; loop_gate wrote .claude/HUMAN_REVIEW_REQUIRED.md
esac
exit 0
