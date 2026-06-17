# arXiv:2311.13027 — Chabanov & Rezzolla, "Numerical modelling of bulk viscosity in neutron stars"

## One line
Full GR implementation of **Müller-Israel-Stewart (MIS)** second-order bulk-viscous hydrodynamics in the FIL/ECHO code, with a robust conservative-to-primitive recovery that treats the bulk pressure Π as an evolved variable. This is the IS comparison baseline (not BDNK) and contains a bulk-viscous migration precedent.

## Relevance to BDNK project
- **Stage 2 (recovery):** Their cons2prim scheme is the cleanest published recipe for recovering Π from a conserved variable D_Π and is directly adaptable. BDNK has no Π as a fundamental field (it is algebraic in gradients), but the *limiting / causality enforcement* and 1D root-find structure are reusable.
- **Stage 4 (collapse/migration tests):** Provides IS baseline numbers for migration and linear damping that BDNK runs can be cross-checked against. Note: their migration does NOT collapse to BH; viscosity stabilizes to a less compact state.

## Formalism (MIS / Maxwell-Cattaneo, NOT BDNK)

Stress-energy (Landau frame, bulk only):
- T^{μν} = e u^μ u^ν + (p + Π) h^{μν}   (Eq. 3)
- ∇_μ J^μ = 0 (baryon), ∇_μ T^{μν} = 0   (Eqs. 4, 5)

Relaxation equation (Maxwell-Cattaneo):
- τ_Π u^μ ∇_μ Π = Π_NS − Π    (Eq. 6)
- Π_NS := −ζ Θ    (Eq. 7)  [Navier-Stokes / first-order limit]
- Θ = ∇_μ u^μ = ϑ + Λ − K W   (Eq. 16, 3+1 split)

Flux-conservative system:
- ∂_t U + ∂_i F^i(U) = S(U, ∂_μ U)   (Eq. 8)

Conservative variables (Eq. 9), U := √γ [D, S_j, τ, D_Π]^T:
- D = ρ W
- S_j = (e + p + Π) W² v_j
- τ = (e + p + Π) W² − (p + Π) − ρ W
- D_Π = ρ W Π

Effective (viscous-modified) quantities (Eq. 23):
- a' = (p + Π) / [ρ(1+ε)]
- h' = h + Π/ρ
- z = W v

## Equation of state
- Hybrid / Γ-law thermal split: **p = κ ρ^Γ + (Γ_th − 1) ρ ε_th**  (Eq. 58)
- Binary merger: cold part = β-equilibrated slice of **TNTYST** EOS; **Γ_th = 1.7**.
- TOV / migration tests use polytrope (κ, Γ).

## Transport coefficients (spatially varying)
- ζ(ρ): constant ζ_l for ρ ≤ ρ_l, **cubic-polynomial** transition for ρ_l < ρ < ρ_h, constant ζ_h for ρ ≥ ρ_h  (Eq. 48).
- τ_Π: **linear** interpolation across transition zone (Eq. 56). For stability τ_h ≈ 1.1 Δt_min, τ_l ≈ 1.1 Δt_max.
- ζ, τ_Π held *constant* in time per run.
- Physical motivation: effective bulk viscosity from violations of weak (β/Urca) chemical equilibrium.

## Causality / sound speed (KEY for BDNK cross-check)
Viscous (effective) sound speed (Eq. 30):
- c_s'² = (ζ / τ_Π)(1 / (ρ h')) + (∂p/∂e)_ρ + (1/h')(∂p/∂ρ)_e

Causality enforcement (Eq. 31): if c_s'² > c_max², reset τ_Π to
- τ_Π = [ζ/(ρ h')] · [ c_max² − (∂p/∂e)_ρ − (1/h')(∂p/∂ρ)_e ]^{−1}

This makes the *gradient* term in c_s'² (the ζ/τ_Π piece) the analog of the BDNK causality constraint: subluminal propagation requires the dissipative contribution to the characteristic speed be bounded. c_max is set just below 1.

## Conservative-to-primitive recovery (Stage-2 relevant)
Treat Π as an extra primitive recovered from D_Π. Uses a Kastaun-style 1D root-find on z = W v with the *primed* enthalpy h':
- f(z) = z − r / h'(z)    (root in z)
- Bracket: z₋ = (k/2)/√(1 − k²/4),  z₊ = k/√(1 − k²)
Limiting applied **inside** the root-find each iteration:
1. If Π < 0 and Π < α p (α ≈ −0.9, with −1 ≤ α ≤ 0): set Π = α p.
2. If Π > 0 and a' > 1: set Π = e − p.
3. If c_s'² > c_max²: adjust τ_Π via Eq. 31.
Caveat: Π is only truly "constant" if the limiter never fires.

## Numerical method
- Code: **FIL** (4th-order finite-difference, Cartesian) on top of IllinoisGRMHD/ETK.
- Hydro: **ECHO** scheme, 4th-order shock-capturing.
- Time integration: **Heun's** 2nd-order, 2-stage. Temporal source derivatives (Θ etc.) via 1st-order backward differencing vs previous full timestep (Eqs. 21-22).
- Spacetime: **Z4c** constraint-damping.

## Benchmark numbers (cross-check targets)

### Linear density-oscillation damping test (Sec III.1)
- TOV: M = 1.4 M_⊙, R = 14.2 km, ρ_c ≈ 7.91×10¹⁴ g/cm³.
- ~2% amplitude fundamental radial mode.
- ζ_h ∈ {0, 9.42×10²⁵, 1.98×10²⁶, 8.20×10²⁶} g s⁻¹ cm⁻¹.
- Resolutions Δx ∈ {207, 281, 369, 487} m.
- Convergence fit: ζ(Δx) = ζ_a + ζ_s (Δx/λ)^p, λ = 2π c_s ω⁻¹.
- Recovered asymptotic ζ_a vs input: ≲6% (low), ≲20% (medium), ≲35% (high).

### Migration test (Sec III.2) — bulk-viscous collapse precedent
- Unstable TOV: ρ_c ≈ 4.94×10¹⁵ g/cm³, M = 1.447 M_⊙.
- ζ ∈ {0, 8×10²⁷, 4×10²⁸, 2×10²⁹, 10³⁰} g s⁻¹ cm⁻¹.
- Δx = 0.25 M_⊙ ≈ 370 m (finest level).
- Outcome: **does NOT collapse to BH**; viscosity damps radial oscillations → heat → less compact asymptotic state, lower central ρ. Shock fronts weaken; highest-ζ second wave stalls at ~60 km.

### Binary merger (Sec III.3)
- ζ_h ∈ {ζ_0, ζ_0/2, ζ_0/5, 0}, ζ_0 = 10³⁰ g cm⁻¹ s⁻¹.
- Transition densities: ρ_h ≈ 4.52×10¹⁴, ρ_l ≈ 1.13×10¹² g/cm³.
- τ_h ≈ 2.7×10⁻⁴ ms.
- Δx ≈ 0.17 M_⊙ ≈ 260 m (6th level).
- Post-merger: Π ≈ 20% of equilibrium p; |R⁻¹| (inverse Reynolds) ~1% peak, → ~10⁻⁵ by t−t_mer ≈ 5 ms.
- High viscosity suppresses dynamically ejected mass by factor ~5.

## Open questions
- No closed-form ζ(ρ,T) functional form given (ζ held constant per run); the microphysical β-equilibrium ζ must come from elsewhere.
- Exact polytropic κ, Γ for the TOV/migration tests not extracted.
- Precise definition of `r` and `k` in the root-find brackets not fully quoted (standard Kastaun: k = |S|/(τ+D+p+Π)-type), should confirm from Eqs. 24-29 if needed for exact reuse.
