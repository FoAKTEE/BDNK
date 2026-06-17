#!/usr/bin/env bash
# inject_infra.sh — BDNK consumer bootstrap (adapted from phys-agentic-loop).
#
# Emits a <session-start-briefing> with the phys-agentic-loop methodology (read
# from the phys-agentic-loop/ subdir) plus this mission's live state: the Ralph
# loop file, research-state note, the project DAG, and the exact-target ledger.
# Wired as the SessionStart hook in .claude/settings.json; also a manual bootstrap.

set -u
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || printf '%s' "$(cd "$(dirname "$0")/.." && pwd)")"
INFRA="$REPO_ROOT/phys-agentic-loop"

emit_file() {
    local label="$1" path="$2"
    printf '\n----- BEGIN %s (%s) -----\n' "$label" "$path"
    if [ -f "$path" ]; then cat "$path"; else printf '(missing: %s)\n' "$path"; fi
    printf '\n----- END %s -----\n' "$label"
}

printf '<session-start-briefing enforcement="MANDATORY">\n'
printf '\n=== PHYS-AGENTIC-LOOP INFRA (source: %s) ===\n' "$INFRA"
emit_file "INDEX.md"                                   "$INFRA/INDEX.md"
emit_file "_common/contracts/markers.md"               "$INFRA/_common/contracts/markers.md"
emit_file "_common/contracts/progress_principles.md"   "$INFRA/_common/contracts/progress_principles.md"
emit_file "notes/multi_timescale_tracking_template.md" "$INFRA/notes/multi_timescale_tracking_template.md"

printf '\n=== BDNK-HMNS MISSION STATE ===\n'
emit_file "RALPH LOOP"      "$REPO_ROOT/.claude/ralph-loop.local.md"
emit_file "RESEARCH STATE"  "$REPO_ROOT/progress/research_state.md"
emit_file "CURRENT ITER"    "$REPO_ROOT/progress/loop_notes/current_iter.md"
emit_file "PROJECT DAG"     "$REPO_ROOT/progress/dag/project_dag.md"

printf '\n=== LIVE DAG PROGRESS (knowledge + error ledgers, doubly linked) ===\n'
if command -v python3 >/dev/null 2>&1; then
    python3 "$INFRA/_common/visualization/dag_mermaid.py" progress --papers bdnk-hmns --repo-root "$REPO_ROOT" 2>/dev/null \
      | python3 -c 'import sys,json
try:
    d=json.load(sys.stdin)
    order=sorted(d,key=lambda x:(x["status"]!="blocking",x["status"]!="preliminary",x["node_id"]))
    for r in order:
        print("  %-11s %-22s k%s t%s(pass %s)" % (r["status"],r["node_id"],r["n_knowledge"],r["n_trials"],r["pass"]))
except Exception: pass'
fi
printf '\nRULES: reproduce EXACTLY (no fabrication); each iter ship ONE verified DAG\n'
printf 'node (append knowledge/error rows, generate+VLM a figure); commit per substage\n'
printf '(no Claude co-author); all commits -> FoAKTEE/BDNK only; <20 CPU, no GPU; Julia.\n'
printf '</session-start-briefing>\n'
