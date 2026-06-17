---
active: true
iteration: 28
session_id: "bdnk-hmns-repro"
max_iterations: 1000
no_progress_limit: 8
stuck_counter_limit: 3
max_wall_seconds: 0
paper: "bdnk-hmns"
completion_promise: "BDNK_HMNS_ALL_REFPAPER_FIGURES_REPRODUCED"
started_at: "2026-06-17T08:00:00Z"
---

BDNK-HMNS reproduction loop. Project code `bdnk-hmns`; Julia package
`code/BDNKStar/`; branch `main`; commits → FoAKTEE/BDNK only (no Claude
co-author). Canonical plan: `progress/research_state.md` (research-state note),
DAG: `progress/dag/project_dag.md`, exact targets: `progress/reproduction/LEDGER.md`.
Constraints: < 20 CPU threads, no GPU, Julia only.

Each iteration runs a tight loop:

1. **PLAN** — Read `progress/research_state.md` + `progress/dag/project_dag.md`.
   Pick the SINGLE highest-priority `[BLOCKING]`/`[FUTURE]`/`[PRELIMINARY]` DAG
   node. State the exact reproduction target (number + source from
   `progress/reproduction/LEDGER.md`), evidence type, and the verifier command.

2. **EDIT** — Implement that one node only into `code/BDNKStar/` as a reusable
   module (one concern per file). Ground every formula in the reference code/.tex
   (no fabrication). Reproduce EXACTLY against the ledger target.

3. **VERIFY** — `cd code/BDNKStar && JULIA_NUM_THREADS=4 julia --project=. test/runtests.jl`;
   generate the validation figure and VLM-inspect it (pre- against the paper
   figure when available, post dual-compare). Numerical-relativity comparison in
   log scale. A failed gate blocks advancing.

4. **COMMIT (per substage)** — Conventional-commit grammar, NO Claude co-author,
   verification output in the body. Push to `origin main` when the credential is
   available (currently queued locally on auth).

5. **UPDATE** — Append knowledge/error ledger rows under the DAG node (doubly
   linked); re-render `project_dag.md`; rewrite `current_iter.md`; extend
   `research_state.md`; rewrite `nodal_note.md` every 10 iters.

6. **ESCALATION** — If a node has not advanced for 3 iters or 30 min, run
   `phys-agentic-loop/pipelines/6-escalation/spec.md` (acquire missing
   sources / run the reference code NSO.jl for an exact cross-check). Sub-agents
   only of the allowed type (dynamical workflow claude 4.8 ultracode / claude 4.8 max).

7. **PROGRESS-GATED TERMINATION.** The Stop guard consults
   `phys-agentic-loop/_common/loop_gate.py decide`. Keep working while it returns
   `continue`; it halts (writing `.claude/HUMAN_REVIEW_REQUIRED.md`) on
   no_progress / stuck_counter / time_budget / max_iterations. Do not defeat the
   gate. On halt: read the gate, fix the root cause, `loop_gate.py reset`, resume.

8. **COMPLETION.** Output the promise ONLY when ALL hold:
   - Every reproducible figure of all reference papers reproduced + VLM dual-checked
     (Caballero-Yunes 2506.09149; Bussières 2604.13208; Shum 2509.15303;
     Redondo-Yuste 2411.16841; Pandya 2201.12317/2209.09265; BDN 2009.11388;
     Kovtun 1907.08191; Chabanov-Rezzolla 2311.13027).
   - Every benchmark in `progress/reproduction/LEDGER.md` matched within tolerance.
   - No `[OPEN]`/`[BLOCKING]`/`[PRELIMINARY]` markers remain in `research_state.md`.
   - Every commit pushed to FoAKTEE/BDNK.

Completion promise (output ONLY when unequivocally TRUE):
<promise>BDNK_HMNS_ALL_REFPAPER_FIGURES_REPRODUCED</promise>
