# Nodal note — window iters 001–003 (full-rewrite every 10 iters; keeps 10)

## 10-iter window
- **error-DB pass/fail** (paper bdnk-hmns): 6 trials, 6 pass / 0 fail in window.
- **node coverage delta**: 0 → 16 DAG nodes created; 5 promoted to SOLID,
  1 to PRELIMINARY. Reproductions landed: STEP-0 round-trip gate, conformal
  recovery + RH shock, TOV (Bussières + Shum), radial Cowling spectrum.
- **simplification cycles consumed**: 0 (greenfield build; no refactor trigger yet).
- **strategic redirects**: scope expanded by user to "reproduce EXACTLY ALL
  results + ALL figures of ALL ref-papers" → reproduction program R1–R5 created.

## Logic-DAG snapshot  (mirrors progress/dag/project_dag.md)
- ● SOLID  step0.eos · step0.con2prim_ideal · step0.causality ·
  step0.bdnk_recovery · s1a.tov_background
- ◐ PRELIMINARY  s1a.radial_eig (converged+stable; NSO exact cross-check [OPEN])
- ✗/□ remaining  s1a.heat_criterion · s1b.axial_wave_eqs · s1b.qnm_freqdomain ·
  s1c.hrsc_core · s1c.imex · s1c.qnm_extract · s2.gr_coupling · s2.is_contrast ·
  s3.cowling_3p1 · s4.production
- external deps: NSO.jl run (Arpack/DifferentialEquations) for exact radial/axial
  cross-validation [OPEN].

## Accepted-results snapshot
| claim | evidence type | verifier | status |
|---|---|---|---|
| STEP-0 prim↔cons round-trip ≤1e-10 (max 4.3e-15) | Convergence | test/runtests.jl 189/189; figures/step0_validation.png | [SOLID] |
| conformal BDNK recovery = solver.c (1.8e-16) | CrossCheck | test_conformal.jl | [SOLID] |
| RH steady shock εR=4.40741, vR=0.41667 | Analytic | test_conformal.jl | [SOLID] |
| TOV Bussières EOS1 M=1.266/R=8.861 (t 1.27/8.86) | NumericalSim+VLM | test_tov.jl; tov_reproduction.png | [SOLID] |
| TOV Shum M_T=1.40016 (t 1.4) | NumericalSim | test_tov.jl | [SOLID] |
| radial Cowling f0=5.503 kHz converged+stable | Convergence+VLM | test_radial.jl; radial_modes.png | [PRELIMINARY] |

## Simplification cycle
- none this window (initial build). Code already modular (one concern/file:
  eos/ tov/ transport/ recovery/ conformal/ perturbations/). Reuse maintained:
  EOS + TOV consumed by RadialModes; recovery shared.

## Failure-mode drift
- two transient build failures, both fixed same-iter (no amended rows needed):
  (1) ShumPolytrope export missing from top module → added re-export;
  (2) radial c_s'∝1/c_s³ NaN at surface (ε=0) → atmosphere floor + barotrope p(ε).
- no new failure_mode enum extensions required.
