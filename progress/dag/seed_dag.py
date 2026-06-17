#!/usr/bin/env python3
"""Seed the BDNK-HMNS project knowledge DAG + error ledger from the source
synthesis (progress/understanding/synthesis.json), doubly linking knowledge and
error rows under each DAG node. Idempotent (append-batch dedups on
status+summary). Run from the repo root:

    python3 progress/dag/seed_dag.py

Renders the Mermaid DAG to progress/dag/project_dag.md and the live HTML board.
"""
from __future__ import annotations
import json, sys, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]          # repo root (…/BDNK)
PAL = ROOT / "phys-agentic-loop"
sys.path.insert(0, str(PAL))
from _common.ledgers import knowledge_database as kdb   # noqa: E402
from _common.ledgers import error_database as edb        # noqa: E402
from _common.visualization import dag_mermaid as dm      # noqa: E402

PAPER = "bdnk-hmns"
SYNTH = json.loads((ROOT / "progress/understanding/synthesis.json").read_text())

# STEP 0 completion state (substage 1, commit 835faae): the EOS, ideal con2prim
# and causality monitor are SOLID; the BDNK gradient-frozen recovery is
# PRELIMINARY (barotropic done; general-EOS linear solve + conformal explicit
# pending). s1a.tov_background is the immediate next buildable node.
SOLID   = {"step0.eos", "step0.con2prim_ideal", "step0.causality"}
PRELIM  = {"step0.bdnk_recovery"}
BLOCK   = {"s1a.tov_background"}
EVID_835 = ("code/BDNKStar test/runtests.jl 189/189 pass; prim<->cons round-trip "
            "<=1e-10 (machine precision); commit 835faae")
DOMAIN_OK = {"symbolic", "numerical", "proof"}


def status_of(n):
    nid = n["node_id"]
    if nid in SOLID:  return "solid"
    if nid in PRELIM: return "preliminary"
    if nid in BLOCK:  return "blocking"
    return "future"


def build_rows():
    rows = []
    for n in SYNTH["dag_nodes"]:
        st = status_of(n)
        dom = n["domain"] if n["domain"] in DOMAIN_OK else "numerical"
        row = {
            "paper": PAPER,
            "node_id": n["node_id"],
            "task_id": n.get("stage", "?"),
            "domain": dom,
            "status": st,
            "summary": n["summary"],
            "predecessors": n.get("predecessors", []),
            "equation_labels": n.get("equation_labels", []) or [],
            "notes": "GATE: " + n.get("gate", "(see plan)"),
        }
        if st == "solid":
            row["evidence"] = EVID_835
            row["metric_at_landing"] = {
                "name": "prim<->cons round-trip max rel err",
                "value": 1.2e-14, "threshold": 1e-10, "pass": True}
            row["code_block_refs"] = ["code/BDNKStar/src"]
        rows.append(row)
    return rows


def seed_error_trial():
    """The STEP 0 round-trip validation trial, anchored under
    step0.con2prim_ideal (the doubly-linked error<->DAG edge)."""
    trial = {
        "paper": PAPER, "node_id": "step0.con2prim_ideal",
        "task_id": "STEP0", "stage": "validation", "domain": "numerical",
        "change_type": "structural", "iteration": 1, "wall_clock_seconds": 5.0,
        "pass_fail": "pass",
        "change_summary": "STEP 0 prim<->cons round-trip across polytrope/ideal/tabulated + BDNK frozen",
        "metric": {"name": "round_trip_max_rel_err", "value": 1.2e-14,
                   "threshold": 1e-10, "pass": True},
    }
    edb.append_row(dict(trial), repo_root=ROOT)


def main():
    res = kdb.append_batch(build_rows(), repo_root=ROOT)
    print("knowledge append-batch:", res)
    seed_error_trial()
    print("error trial appended under step0.con2prim_ideal")

    # render the project DAG as Mermaid (single paper view)
    md = dm.render_single(PAPER, repo_root=ROOT)
    out = ROOT / "progress/dag/project_dag.md"
    out.write_text(md, encoding="utf-8")
    print("rendered", out)
    # per-node progress readout
    prog = dm.node_progress([PAPER], repo_root=ROOT)
    n_solid = sum(1 for p in prog if p["status"] == "solid")
    print(f"DAG: {len(prog)} nodes; {n_solid} solid")


if __name__ == "__main__":
    main()
