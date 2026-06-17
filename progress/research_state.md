# BDNK-HMNS — research state (append/modify long-memory note)

**Mission.** Build a causal first-order (BDNK) viscous neutron-star code in Julia
to quantify how bulk/shear viscosity + heat conduction affect a hot
out-of-equilibrium HMNS remnant (radial/collapse stability + non-radial
spectrum). Reproduce before extend; validate every gate numerically.

**Repo.** github.com/FoAKTEE/BDNK (commits here only; no Claude co-author).
Branch `main`. Compute: < 20 CPU threads, no GPU. Language: Julia only (physics);
phys-agentic-loop (Python) is the management/ledger tool.

## Source library (downloaded + reformulated)
- 9 arXiv sources (`.tex` in `ref-paper/sources/`, gitignored) + 3 ref-code repos
  (`ref-code/`, gitignored): 1D_conformal_bdnk (C), BDNK_frame_constraints (NB),
  NeutronStarOscillations.jl (Julia, in-language radial/axial reference).
- Reformulated: 14 structured notes `progress/understanding/*.md`, cross-walk +
  benchmark table + build spec `SYNTHESIS.md`, index `INDEX.md`.

## Working context (canonical conventions)
- Geometrized G=c=1, lengths km (Units.jl). Mostly-plus metric.
- General-frame BDNK: T^{ab}=(ε+A)u^a u^b+(P+Π)Δ^{ab}+Q^a u^b+Q^b u^a−2η σ^{ab};
  frame coeffs {τ_ε,τ_P,τ_Q,η,ζ,σ,β_ε,β_n} are FIRST-ORDER (no IS evolution eqn).
- Ideal-gas micro: P=(Γ−1)ρϵ, c_s²=ΓP/(ε+P)=Γ(Γ−1)ϵ/(1+Γϵ). T=p/ρ.
- Shum 1C EOS: cold Γ=2, p(ε)=[1+2κε−√(1+4κε)]/(2κ), c_s²=1−1/√(1+4κε).

## DAG status  (knowledge-database/paper_bdnk-hmns, 16 nodes; Mermaid: progress/dag/project_dag.md)
- ● SOLID: step0.eos, step0.con2prim_ideal, step0.causality
- ◐ PRELIMINARY: step0.bdnk_recovery (barotropic gradient-frozen done; general-EOS
  linear-solve + conformal explicit compute_xiD/uxD port pending)
- ✗ BLOCKING (next): s1a.tov_background
- □ FUTURE: s1a.radial_eig, s1a.heat_criterion, s1b.axial_wave_eqs,
  s1b.qnm_freqdomain, s1c.hrsc_core, s1c.imex, s1c.qnm_extract, s2.gr_coupling,
  s2.is_contrast, s3.cowling_3p1, s4.production

## Accepted results (gates met)
1. **STEP 0 round-trip gate** [SOLID]. prim↔cons closes to ≤1e-10 across
   polytrope / Shum-Γ2 / ideal-gas / tabulated / BDNK-frozen. Measured max
   4.3e-15. Causality monitor solves the char-speed biquadratic exactly.
   Evidence: `code/BDNKStar` 189/189 tests; figure `figures/step0_validation.png`.
   Commits 835faae, (substage-4).

## Open questions / human decision points (defaults adopted, flagged for confirm)
- {{EOS_TABLE}}: default Shum analytic Γ=2 κ=100 for 1C trunk; real tabulated EOS
  is a Stage-2 decision.
- {{TRANSPORT}}: default dimensionless parametrized hats for STEP0/1A/1C; physical
  η_c[g cm⁻¹ s⁻¹] for Stage-4 production.
- {{FRAME_SET}}: Shum (ŝ,â,q̂)=(1,1,0.999) for 1C; conformal luminal (25/4,25/7)η₀
  for STEP-0; ≥2 frames at Stage 2.
- {{TOL}}: Δ_visc PF-fallback OFF (TOL<0) for benchmarks.
- {{GR_BACKEND}}: Cowling for 1A/1B/1C/3; dynamical 1+1D GR for Stage 2.

## Benchmark targets (from SYNTHESIS.md benchmark_table — to reproduce)
- 1A: TOV n=1 κ=100; heat-conduction stability criterion c_s²−c_n²≥0.
- 1B: Bussières EOS1 M=1.27M☉ R=8.86km; ℓ=2 w-mode (10.4884 kHz, 29.587 µs).
- 1C: Shum M_T=1.4M☉ (ρ0c=0.00128); QNM F=2.69 kHz, H1≈4.55, H2≈6.36 kHz;
  1/τ≈0.00157 M☉⁻¹; stable window τ_ε=(4/3)η̂+ζ̂ ≲ 0.1.

## Next 3
1. step0.bdnk_recovery → solid: port conformal compute_xiD/uxD (solver.c), verify
   CrossCheck round-trip; scaffold general-EOS linear solve.
2. s1a.tov_background: TOV RK4 + implicit-CN cross-check (NSO ref); gate
   e^ν→1−2M/R exterior, cross-solver 1e-6.
3. s1a.radial_eig + heat_criterion: radial Cowling eigensolver; CY c_s²−c_n²≥0.
