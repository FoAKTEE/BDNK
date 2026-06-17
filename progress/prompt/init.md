
## **CENTRAL TASK**
* it is now demonstrated that should be stored in a directed acyclic graph in knowledge database, use Mermaid
* check prompt message
* deep search about other ref-code and ref-paper on weband run decomposition and escalation frequently!
* use fewer than 20 CPUs and no GPUs
* all code should be written in julia only


**ref-paper**
* prompt message + ref-paper/paper/*

**ref-code**
https://github.com/aapandy2/1D_conformal_bdnk
https://github.com/aapandy2/BDNK_frame_constraints

**management**
* audit, modify, utilize phys-agentic-loop for management and monitor
* git repo - substage progress doc + commit - https://github.com/FoAKTEE/BDNK
    - remove claude coauthor
    - all commits should go to https://github.com/FoAKTEE/BDNK only!


**modification to knowledge and error database structure**
* both databases should be directly connected to the new DAG
    - labeling: knowledge/error under which node 
      - under the same node: knowledge/error list sequentially
    - knowledge/error database and DAG should be doubly linked for debug

* prompt message:
This is an executable specification for an autonomous research-engineering agent, or for a multi-node AI pipeline. Run it end-to-end with a single agent, or run each STEP/STAGE block as an independent node — every block is self-contained (TASK → METHOD → DELIVERABLES → GATE → ON FAILURE). Do not start a block until its upstream GATE has passed. Fill the {{...}} parameters before launching.
Parameters to set (orchestrator/human supplies)

{{EOS_TABLE}} — finite-temperature tabulated EOS, p(ρ, T, Yₑ), with source and citation.
{{TRANSPORT}} — parametrized (dimensionless η, ζ, κ, τ) orphysical (named dense-matter prescriptions).
{{TOL}} — acceptance tolerance for each benchmark match (e.g. 2% on mode frequencies, 5% on decay rates); may differ per gate.
{{FRAME_SET}} — ≥2 distinct BDNK hydrodynamic-frame parametrizations.
{{GR_BACKEND}} — the existing IS/MIS infrastructure to fork (Valencia + BSSN/CCZ4, or Misner–Sharp radial gauge).
{{COMPUTE}} — resolution ceiling and wall-clock/GPU budget per run.
Role and objective
You are building a causal, first-order viscous (BDNK) neutron-star code to quantify how bulk viscosity, shear viscosity, and heat conduction affect a hot, out-of-equilibrium post-merger hypermassive remnant — specifically its radial/collapse stability and its non-radial oscillation spectrum. Work in two forks off one shared core, validating every step against published linear and nonlinear results before extending.
Global rules (enforce in every block)

Reproduce before extend. Never pass a GATE you have not numerically met within {{TOL}}. On failure, iterate; after two failed iterations, escalate and stop — do not proceed on an unmet gate.
No fabrication. Validate against the Benchmarks below by computing, not asserting. You may fetch and read the cited papers (arXiv IDs given). If a benchmark value you need cannot be extracted, HALT and request it — never invent a published number or report a match you did not verify.
Reproducibility. For every run, persist the full parameter set, code commit SHA, RNG seeds, grid, and EOS/coefficient choices. Emit self-convergence plots and constraint-monitor time series as artifacts.
Self-consistency. The primitive:left_right_arrow:conservative round-trip must close to ≤ 1e-10 in smooth regions. Constraint monitors are first-class outputs, not afterthoughts.
Frame discipline. A claimed mode or instability is physical only if it persists across {{FRAME_SET}} and converges under refinement at small Knudsen number Kn ~ τω. Report frame dependence explicitly.
Novelty discipline. Do not assert any "first ___" without a literature check; flag every such claim for human confirmation.
Per-block reporting. Each block writes REPORT.md: what ran, GATE result (pass/fail + the numbers), plots, and open issues.


STEP 0 — Shared EOS + primitive-recovery module (the trunk; build first)
TASK. One abstracted module reused by every later stage. The two "not done" items in the plan — realistic-EOS recovery (Stage 2) and 3+1D polytropic recovery (Stage 3) — are this one module with different inputs; build it once, cleanly. METHOD.

EOS driver: polytropic, Γ-law ideal gas, and tabulated {{EOS_TABLE}}. Expose thermodynamically consistent derivatives (cₛ², ∂p/∂ε|_ρ, ∂p/∂ρ|_ε); enforce monotonicity; implement a density floor / atmosphere.
BDNK recovery: the conserved densities carry the dissipative gradient corrections, so do not use the algebraic ideal-hydro inversion. Evaluate the first-derivative dissipative terms from grid data, freeze them, and solve the ideal-like inversion shifted by those known source terms with a 1D/2D Newton iteration; add an outer iteration if gradient self-consistency is required (cf. Pandya–Most–Pretorius).
Causality/stability: implement the BDNK coefficient inequalities; monitor them pointwise; expose a violation flag. DELIVERABLES. Module + unit tests. GATE. Round-trip ≤ 1e-10 across polytropic → ideal → tabulated; tabulated recovery convergent; causality monitor operative. ON FAILURE. Localize (the stellar surface/atmosphere is the usual culprit); tighten Newton globalization; if tabulated derivatives are inconsistent, escalate the {{EOS_TABLE}} choice.
STAGE 1A — Radial linear benchmark (Caballero–Yunes)
DEPENDS ON: STEP 0 (EOS only). TASK. A frequency-domain / shooting eigenvalue solver for radial perturbations of a TOV star under the Eckart, BDNK, and Müller–Israel–Stewart closures, including bulk, shear, and heat conduction. GATE. Recover, within {{TOL}}: stability to bulk and shear; the conditional heat-conduction instability; and the heat-conduction stability criterion.
STAGE 1B — Axial linear benchmark (Redondo-Yuste / Bussières)
TASK. A linear axial-sector solver; reduce the problem to the coupled wave equations of Redondo-Yuste; compute the QNM spectrum and identify the viscosity-driven mode families that have no perfect-fluid counterpart. GATE. Reproduce the Bussières axial spectrum and the new families within {{TOL}}.
STAGE 1C — Nonlinear core (Shum reproduction)
DEPENDS ON: STEP 0. TASK. Nonlinear, spherically symmetric BDNK hydrodynamics under the Cowling approximation, simplified EOS. This is the engine both forks inherit. METHOD. High-resolution shock capturing (HLL-type or central / Kurganov–Tadmor with high-order reconstruction) on the advective operator; IMEX Runge–Kutta for the stiff relaxation-time terms; positivity limiting; hardened atmosphere. Extract QNM frequency content and the fundamental-mode decay rate. GATE. Stable evolution within the restricted parameter window; QNM frequencies and decay rate match Shum within {{TOL}}; self-convergence at the design order.
STAGE 2 — Fork A: 1+1D dynamical GR + realistic EOS + collapse
DEPENDS ON: 1C. TASK. Drop Cowling; couple BDNK hydro to a dynamical spacetime; add baryon number; run viscous spherical collapse. METHOD.

GR coupling: fork {{GR_BACKEND}} and substitute the BDNK stress tensor for the IS one; carry constraint damping (cf. Fantini–Rubio). Apparent-horizon finding is a 1D root find for the marginally trapped surface.
Add baryon-number conservation + dissipative heat/number flux; realistic EOS via STEP 0.
Numerics: HRSC + IMEX + positivity + hardened atmosphere.
Frames: evolve in {{FRAME_SET}}; physical observables must converge across frames at small Kn.
Seed collapse: pressure depletion, inward velocity, or a configuration just beyond M_max. SCIENCE. Q1 — does full GR coupling change radial stability? (focus the heat-conduction channel; check against the CY criterion). Q2 — how does viscosity affect the collapse? GATE. Linear radial spectrum matches the CY criterion across {{FRAME_SET}}; convergent collapse with apparent-horizon formation; constraints preserved under refinement.
STAGE 3 — Fork B: 3+1D Cowling, non-radial modes
DEPENDS ON: 1C. TASK. Lift the core to a 3+1D angular grid (Cowling retained) to resolve non-radial modes. METHOD. Extend the STEP 0 recovery with the angular-derivative gradient terms (this is the "3+1D polytropic recovery"); the two-sphere formulation (Keeble–Pretorius) is a useful precursor. Validate the linear axial sector against Bussières: seed small ℓ ≥ 2 axial perturbations, extract the QNMs, recover the new viscous families. Then enable the polar sector and mode coupling. SCIENCE. (i) new modes from non-equilibrium effects; (ii) modified stability; (iii) excitation of the non-radial modes. GATE. Linear axial spectrum matches the Bussières new families within {{TOL}}; self-convergence.
STAGE 4 — Production
TASK. Run the science questions to convergence. Include IS comparison runs (same code, swapped closure) to separate physical effects from BDNK frame artifacts. Produce final plots/tables and a synthesis report. GATE. Every science claim convergence-tested, frame-checked, and IS-cross-checked; all novelty claims human-confirmed.

Decision points (escalate to human)

{{TRANSPORT}} parametrized vs physical, and which prescriptions — this determines whether Stage 2 is a methods paper or the astrophysics result.
{{EOS_TABLE}} source.
per-gate {{TOL}}.
any "first ___" novelty claim.
{{COMPUTE}} resolution/budget ceilings.
Benchmarks (validate against; never fabricate)

Caballero & Yunes, Phys. Rev. D 112, 063050 (2025), arXiv:2506.09149 — radial; stable to bulk/shear, conditional heat-conduction instability + stability criterion.
Bussières, Redondo-Yuste, Ortega Gómez, Cardoso, arXiv:2604.13208 (2026) — axial; new viscosity-driven mode families with no perfect-fluid counterpart.
Shum, Abalos, Bea, Bezares, Figueras, Palenzuela, Phys. Rev. D 113, 084029 (2026), arXiv:2509.15303 — first nonlinear spherically symmetric BDNK under Cowling; QNMs + fundamental decay rate; restricted stable window.
Redondo-Yuste, Class. Quantum Grav. 42, 075012 (2025), arXiv:2411.16841 — axial sector → two coupled wave equations, one a novel viscous mode.
Method references: Pandya–Most–Pretorius, PRD 106, 123036 (2022); Keeble–Pretorius, PRD 112, 124034 (2025); Fantini–Rubio, PRD 112, 063038 (2025); Chabanov–Rezzolla, arXiv:2311.13027 (IS comparison + bulk-viscous collapse precedent).