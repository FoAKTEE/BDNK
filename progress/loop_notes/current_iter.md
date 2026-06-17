# Iter 003 — radial Cowling eigensolver (STAGE 1A)

**Paper anchor.** Caballero–Yunes 2506.09149 radial sector / NSO.jl
PerfectFluidCowling A0/A1/A2 operator; DAG node s1a.radial_eig.

**What shipped this iter.**
- commit 6c2a86b feat(radial): radial Cowling eigensolver.
- src/perturbations/RadialModes.jl (matrix eigensolver, ξ(0)=0 + Δp=0 BCs);
  EOS d2pde2; test_radial.jl; viz/radial_modes.jl.
- error-DB: 1 pass trial under s1a.radial_eig; knowledge: s1a.radial_eig→preliminary.
- figure radial_modes.png (VLM-ok).

**Next-3 roadmap.**
1. s1a.heat_criterion — Caballero–Yunes c_s²−c_n²≥0 (needs 2-param EOS for c_n²).
2. Cross-validate radial f0 vs NSO.jl run (install Arpack/DiffEq) → s1a.radial_eig solid.
3. s1b axial wave eqs + QNM (Redondo-Yuste/Bussières); s1c nonlinear Cowling core.

**Simplification flag.** not required (greenfield; modular).

**Verifier output.**
`julia --project=. test/runtests.jl` → all pass (EOS/recovery/causality/tov/
conformal/radial). radial spectrum f[kHz]=[5.5032,9.0961,12.5707,16.0148],
converged 2nd-order (|Δf0|∝N⁻²), all ω²>0 (stable).
