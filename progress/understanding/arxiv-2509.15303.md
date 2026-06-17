# arXiv:2509.15303 — Shum, Abalos, Bea, Bezares, Figueras, Palenzuela
"Neutron star evolution with the Bemfica-Disconzi-Noronha-Kovtun viscous hydrodynamics framework"
PRD 113, 084029 (2026).

**Role in our project: STAGE 1C reproduction target (the nonlinear trunk / acceptance test).**
First nonlinear, spherically-symmetric BDNK neutron-star simulation under the Cowling approximation. Reproducing their (i) stable parameter window, (ii) QNM frequencies, and (iii) fundamental-mode decay rate is the Stage-1 gate for the BDNK Cowling core.

> NOTE: numbers below extracted via ar5iv full text. Cross-checked between two independent fetches. A few coefficient sign/label details should be re-verified against the PDF before final implementation (flagged in open questions).

---

## 1. Equation of State (EOS)

A "simplified" EOS built by combining a polytrope `p = κ ρ₀^Γ` with an ideal-gas relation `p = (Γ-1) ε₀ ρ₀`, eliminating ρ₀ to give a closed `p(ε)`.

- **Γ = 2**, **κ = 100** (note: this κ is the polytropic constant, NOT a heat-conduction coefficient).
- Closed form (Eq. 53):
  ```
  p(ε) = [ 1 + 2εκ − √(1 + 4εκ) ] / (2κ)
  ```
- Sound speed `c_s² = p'(ε)` follows by differentiation.

**Star configuration:**
- central rest-mass density ρ₀,c = 0.00128 M_⊙⁻²
- total gravitational mass M_T = 1.4 M_⊙
- Initial data: TOV / hydrostatic equilibrium (dissipation vanishes in equilibrium).

---

## 2. BDNK stress tensor, frame, transport coefficients

General-frame BDNK stress-energy tensor (Eq. 6):
```
T^{μν} = (ε + A) u^μ u^ν + (p + Π) Δ^{μν} + Q^μ u^ν + u^μ Q^ν − 2 η σ^{μν}
```
with projector Δ^{μν} = g^{μν} + u^μ u^ν.

Dissipative corrections (general / Bemfica-Kovtun frame):
```
A   = τ_ε [ u^μ ∇_μ ε + (ε+p) ∇_μ u^μ ]
Π   = −ζ ∇_μ u^μ + τ_p [ u^μ ∇_μ ε + (ε+p) ∇_μ u^μ ]
Q^μ = τ_Q [ (ε+p) u^ν ∇_ν u^μ + p'(ε) Δ^{μν} ∇_ν ε ]
```
Shear tensor (Eq. 7):
```
σ^{μν} = (1/2)[ Δ^{μα}Δ^{νβ}(∇_α u_β + ∇_β u_α) − (2/3) Δ^{μν} Δ^{αβ} ∇_α u_β ]
```
Heat-flux coupling closure (Eq. 8): β_ε = τ_Q p'(ε).

Transport coefficients: η (shear visc.), ζ (bulk visc.), and relaxation-time-like
coefficients τ_ε, τ_p, τ_Q (these are the BDNK "frame" coefficients; they
multiply first derivatives, NOT second-order relaxation terms as in Israel-Stewart).

### Dimensionless parametrization (Eq. 65)
With ρ ≡ ε + p, length scale L = 1, and hatted dimensionless controls:
```
η   ≡ q̂ L c_s² ρ η̂
ζ   ≡ q̂ L c_s² ρ ζ̂
V̂   ≡ (4/3) η̂ + ζ̂
β_ε ≡ c_s² â V̂ L
τ_p ≡ ŝ c_s² L V̂
τ_Q ≡ â L V̂
τ_ε ≡ V̂ L
```
All hatted quantities > 0; ζ̂ ≥ 0. The frame is fixed by (ŝ, â, q̂); viscosity
magnitude by (η̂, ζ̂). Note τ_ε = V̂ L = (4/3)η̂ + ζ̂ when L=1 — this combination
is the empirical stability control parameter.

---

## 3. Causality / stability inequalities

Characteristic speeds (Eqs. 67-68):
```
c_0 = c_s √( q̂ η̂ / (â V̂) )

c_± = c_s √{ [ â(1+ŝ) + q̂ ± √( q̂² + â²(4q̂ + (ŝ−1)²) + 2 â q̂ (1+ŝ) ) ] / (2 â) }
```
Well-posedness / strong hyperbolicity (Eq. 69), also the linear-stability condition:
```
0 < q̂ < ŝ
```
Causality (Eqs. 70-71):
```
q̂ < [ (1 − c_s²)/c_s² ] · [ (1 − ŝ c_s²) / (c_s² + â⁻¹) ]
ŝ < 1 / c_s²
```

---

## 4. Restricted STABLE parameter window

Empirical bound found in the simulations (≈ Eq. 73, Sec. III.1):
```
τ_ε = (4/3) η̂ + ζ̂  ≲ 0.1
```
i.e. stable evolutions (to the explored simulation times) require the combined
viscous/relaxation control below ~0.1.

**Four parameter cases run (all with frame ŝ = 1, â = 1, q̂ = 0.999):**

| Case        | τ_ε   | η̂      | ζ̂     |
|-------------|-------|--------|-------|
| smallSB-F2  | 0.023 | 0.01   | 0.01  |
| medS-F2     | 0.023 | 0.01725| 0     |
| highB-F9    | 0.092 | 0.0015 | 0.09  |
| medSB-F9    | 0.092 | 0.03525| 0.045 |

(F2 / F9 labels track the two τ_ε families ~0.023 and ~0.092.)
For smallSB-F2, η ≈ 0.00999 c_s² (ε+p) in geometric (M_⊙) units.

---

## 5. QNM frequencies and fundamental-mode decay rate

**Table 1 — mode frequencies (kHz):**

| Mode | PF (perfect fluid) | smallSB-F2 | highB-F9 |
|------|--------------------|-----------|----------|
| F  (fundamental) | 2.69 | 2.69 | 2.67 |
| H1 (1st overtone)| 4.55 | 4.60 | 4.60 |
| H2 (2nd overtone)| 6.36 | 6.36 | 6.30 |

**Table 2 — fundamental-mode decay rate (geometric units, M_⊙⁻¹), Δr = 0.002 M_⊙:**

| Case        | 1/τ_l (M_⊙⁻¹) | 1/τ_nl (M_⊙⁻¹) | ω_nl (M_⊙⁻¹) |
|-------------|---------------|----------------|--------------|
| smallSB-F2  | 0.00157       | 0.00157        | 0.0834       |
| medS-F2     | 0.00150       | 0.00150        | 0.0834       |
| highB-F9    | 0.00215       | 0.00215        | 0.0834       |
| medSB-F9    | 0.00182       | 0.00182        | 0.0834       |

Fit model (damped sinusoid) for central-energy-density oscillation:
```
ε̃_c(t) = A exp(−t/τ) cos(ω t + φ_0) + C
```
Subscripts l / nl = linear / nonlinear extracted damping; ω_nl ≈ 0.0834 M_⊙⁻¹
corresponds to f ≈ 2.71 kHz (consistent with the F mode). Higher viscosity (highB)
gives the largest decay rate, as expected physically.

Unit conversion: 1 M_⊙ (time) ≈ 4.925e-6 s, so M_⊙⁻¹ ≈ 203 kHz angular-rate scale.
1/τ_l = 0.00157 M_⊙⁻¹ ⇒ damping time τ ≈ 637 M_⊙ ≈ 3.14 ms.

---

## 6. Numerical method

- Conservation-law form (Eq. 19): ∂_t q + ∂_i F^i(q) = S(q).
- Spatial: third-order finite-volume (FDOC, finite-difference Osher-Chakrabarthy);
  "equivalent to 4th-order finite difference with 3rd-order dissipation."
- Time integrator: third-order Strong-Stability-Preserving Runge-Kutta (SSP-RK3).
- CFL: Δt/Δr = 0.25.
- Numerical dissipation set by max characteristic speed (Eq. 72), with a floor
  of 0.1 c to stabilize the stellar surface / atmosphere.
- Grid: staggered (avoids r=0), Δr ∈ [0.001, 0.0032] M_⊙; r_max = 20 M_⊙;
  outflow / outer boundary conditions.
- Atmosphere: cells with p < κ ρ₀,atm^Γ (ρ₀,atm = 1e-12 M_⊙⁻²) reset to
  ρ₀ = 1e-13 M_⊙⁻², v=0, and dissipative vars ε̂=0, v̄̂=0.
- Evolution variables (spherical sym): { γ̃E, γ̃S_r, ε, ∂_r ε, ṽ^r, ∂_r ṽ^r },
  with regularized radial velocity ṽ^r = v^r / r.
- Final times: t_f = 8000 M_⊙ (4500 M_⊙ at the highest resolution).

---

## 7. Cowling approximation setup

Fluid evolved on a FIXED background spacetime (perfect-fluid TOV equilibrium metric).
- Solve hydrostatic equations (Eqs. 56-58) for α(r), a(r), p(r) in areal-polar
  (Schwarzschild-like) coordinates.
- Transform to isotropic coordinates (Eq. 61):
  ```
  ds² = −α²(r) dt² + ψ⁴(r) ( dr² + r² dΩ² )
  ```
- Initial fluid: v^i = 0, ε̂ = 0, v̄̂ = 0. Perturbation/oscillation seeded by the
  truncation error / initial-data mismatch; QNMs read off the resulting ringdown.

---

## 8. Convergence

- "Qualitative convergence" demonstrated across Δr ∈ [0.001, 0.0032] M_⊙ (Fig. 2);
  deviations in ε(r) decrease with resolution.
- Quantitative convergence study in Appendix B.
- Formal order: 3rd-order spatial FV (4th-order FD-equivalent), 3rd-order SSP-RK3
  in time.

---

## 9. Conservative → primitive recovery (KEY for our Stage 0/1C)

BDNK key feature: the conserved densities carry the first-derivative (dissipative)
gradient corrections, but in the chosen formulation those corrections are LINEAR in
the auxiliary time-derivative variables, so recovery is a LINEAR inversion (no
nonlinear Newton root-find here, unlike Pandya-Most-Pretorius's general approach).

Variables:
```
primitives  p₀ = (ε, v^i),   p₁ = (ε̂, v̄̂^i)
   with  ε̂  = −n^μ ∇_μ ε      (negative normal derivative of ε)
         v̄̂^i = γ^i_α n^ν ∇_ν v^α
conserved   q = (γ E, γ S_i)
```
Relations (Eqs. 28-30, schematic):
```
E    = W² ε − p(1 − W²) + [viscous corrections linear in ε̂, v̄̂]
S^i  = − v^i W² (ε + p) + [viscous corrections]
S^{ij} = p γ^{ij} + W² (ε+p) v_i v_j + [viscous corrections]
```
Algorithm:
1. Given p₀ at current step and q, the viscous terms (A, Π, Q^μ) are linear in
   (ε̂, v̄̂), so solving for the state is a linear (≤4×4) system solved analytically.
2. After advancing to the next timestep and obtaining p₀ there, reconstruct p₁ =
   (ε̂, v̄̂) from linear combinations of the E, S^i equations.
3. No nonlinear root-finding required.

Spherical-symmetry evolution of the auxiliary/derivative variables (Eqs. 43-48):
```
∂_t ε        = −α ε̂
∂_t ṽ^r      = α(−v̄̂^r / r + K^r_r ṽ^r)
∂_t(∂_r ε)   = −∂_r(α ε̂)
∂_t(∂_r ṽ^r) = ∂_r[ α(−v̄̂^r / r + K^r_r ṽ^r) ]
```

---

## Open questions / to verify against the PDF
- Exact signs/index placement in A, Π, Q^μ and the precise definition of n^μ vs u^μ
  in ε̂, v̄̂ (transcribed from ar5iv summary; confirm against Eqs. 6-8, 28-30).
- Whether β_ε appears as an independent coefficient or is fully fixed by τ_Q p'(ε).
- Exact form of the hydrostatic equations 56-58 and the areal→isotropic transform.
- Precise Table-2 separation between "linear" (1/τ_l) and "nonlinear" (1/τ_nl)
  damping fits and the resolution dependence (Δr=0.002 quoted).
- Detailed convergence-order number from Appendix B (only "qualitative" stated in body).
