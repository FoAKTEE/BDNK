# Source library — downloaded + reformulated

Every reference paper's arXiv **source** (LaTeX) is downloaded to
`ref-paper/sources/arXiv-<id>/src/` (gitignored: copyright / local-only) and
**reformulated** into a structured note here (tracked). The notes feed the
project knowledge DAG (`progress/dag/project_dag.md`), whose `equation_labels`
back-reference these papers; the whole-project synthesis is `SYNTHESIS.md`.

| arXiv / code | reformulated note | role | DAG nodes |
|---|---|---|---|
| 2209.09265 Pandya–Most–Pretorius (recovery, ideal-gas microphysics) | `arxiv-2209.09265.md` | STEP 0 recovery + causality | step0.eos, step0.con2prim_ideal, step0.bdnk_recovery |
| 2201.12317 Pandya–Most–Pretorius (conformal numerical method) | `arxiv-2201.12317.md` | 1C method, conformal frame | step0.bdnk_recovery, s2.gr_coupling |
| 2009.11388 Bemfica–Disconzi–Noronha (foundational BDNK / PRX) | `arxiv-2009.11388.md` | causality/stability inequalities | step0.causality |
| 1907.08191 Kovtun (first-order stable hydro) | `arxiv-1907.08191.md` | general-frame constitutive relations | step0.eos, step0.causality |
| 2506.09149 Caballero–Yunes (radial perturbations) | `arxiv-2506.09149.md` | STAGE 1A target | s1a.radial_eig, s1a.heat_criterion |
| 2411.16841 Redondo-Yuste (dissipative-star perturbations) | `arxiv-2411.16841.md` | STAGE 1B axial wave eqs | s1b.axial_wave_eqs |
| 2604.13208 Bussières et al. (axial oscillations of viscous NS) | `arxiv-2604.13208.md` | STAGE 1B target (benchmarks) | s1b.axial_wave_eqs, s1b.qnm_freqdomain |
| 2509.15303 Shum et al. (nonlinear BDNK NS, Cowling) | `arxiv-2509.15303.md` | STAGE 1C target | s1c.hrsc_core, s1c.imex, s1c.qnm_extract |
| 2311.13027 Chabanov–Rezzolla (IS bulk viscosity) | `arxiv-2311.13027.md` | STAGE 2 IS contrast | s2.is_contrast |
| Fantini–Rubio 2025 (constraint evolution) | `fantini-rubio-2025.md` | STAGE 2 constraint damping | s2.gr_coupling |
| Keeble–Pretorius 2025 (two-sphere viscous hydro) | `keeble-pretorius-2025.md` | STAGE 3 angular sector | s3.cowling_3p1 |
| **code** 1D_conformal_bdnk (Pandya, C) | `c-1D_conformal_bdnk.md` | recovery + WENO/KT/Heun reference | step0.bdnk_recovery, s1c.hrsc_core |
| **code** BDNK_frame_constraints (Mathematica) | `nb-BDNK_frame_constraints.md` | frame inequalities | step0.causality |
| **code** NeutronStarOscillations.jl (Keeble–Redondo-Yuste, Julia) | `jl-NeutronStarOscillations.md` | TOV + radial/axial linear, in-language | s1a.*, s1b.* |

## Reformulation passes done

1. **Per-source extraction** — equations, algorithms, benchmark numbers, code
   pointers (14 notes above; workflow `bdnk-understand`).
2. **Cross-paper notation cross-walk** — Kovtun ↔ BDN ↔ PMP ↔ Shum unified into
   one canonical general-frame stress tensor (`SYNTHESIS.md` §bdnk_formalism).
3. **Project DAG decomposition** — 16 build nodes STEP0→STAGE4, doubly linked to
   the knowledge + error ledgers (`progress/dag/project_dag.md`).

## Pending reformulation (next iterations)

- Per-paper *equation* DAGs (each paper → its own `logic.md` node graph) and the
  cross-paper duplicate-collapse pass (`dag_mermaid duplicates`).
