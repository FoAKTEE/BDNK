# Iter 004 — conformal evolution engine + CY heat criterion

**Paper anchors.** Pandya 2201.12317 / solver.c (flat-space BDNK evolution);
Caballero–Yunes 2506.09149 (heat-conduction criterion).

**What shipped this iter (local commits 77bef47, f77f026).**
- src/conformal/ConformalEvolution.jl — full WENO5+KT+Heun engine (solver.c port):
  constant state exact (Δε=0), Gaussian→2 pulses, steady shock RH-consistent.
- EOS cn2 + heat_conduction_stable; reproduce CY Δ=c_s²−c_n²=−(Γ-1)/(1+Γϵ)<0.
- figures conformal_evolution.png, heat_criterion.png (both VLM-ok).
- test suite 271→ (EOS/recovery/causality/tov/conformal/conf-evolution/radial/heat).

**Push status.** BLOCKED — VS Code git-askpass token expired (IDE/MCP dropped).
2 commits queued locally; push on credential restore (user: "keep building locally").

**Next-3 roadmap.**
1. R4 axial QNM (Bussières 2604.13208 w-modes 10.4884 kHz/29.587 µs; Redondo-Yuste
   coupled wave eqs 14a/14b; shooting + Leaver continued fraction).
2. R5 nonlinear stellar Cowling (Shum): areal→isotropic, BDNK evolve on TOV
   background w/ general EOS, FFT QNM (F=2.69/H1=4.55/H2=6.36) + decay 0.00157.
3. NSO.jl run to lock exact radial/axial cross-checks.

**Verifier output.** `julia test/runtests.jl` → 271/271 pass; constant Δε=0.0;
CY Δ exact to 1e-12.
