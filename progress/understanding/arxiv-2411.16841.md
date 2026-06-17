# arXiv:2411.16841 — Redondo-Yuste, "Perturbations of relativistic dissipative stars" (CQG 42, 075012, 2025)

**Stage focus: 1B — axial sector reduction to two coupled wave equations + novel viscous mode.**

## One-line
Linearized BDNK perturbation theory of spherically symmetric self-gravitating viscous fluids; the axial (odd-parity) sector reduces to **two coupled wave equations** (one a damped GW mode, one a NOVEL rotational viscous mode with no perfect-fluid analogue); polar sector reduces to 5 coupled wave equations + 1 constraint. Causal structure recovers BDNK causality bounds. **No QNM numbers, no specific stellar model — pure formalism, ready for numerical study.**

## What is and is NOT in this paper
- **IS:** Full covariant linearized perturbation equations for BDNK self-gravitating fluid; frame-basis (comoving) projection; axial → 2 coupled wave eqs; polar → 5 coupled wave eqs + constraint; characteristic/causality discussion; companion Wolfram package.
- **IS NOT:** No numerical QNM frequencies. No specific EoS chosen (only barotropic p=p(e) assumed). No specific compactness / star model integrated. No explicit numeric values for transport coefficients. The paper itself states a first-principles QNM/damping-time calculation "is still lacking."

## BDNK formalism (general frame)
Stress-energy tensor (eq. 16):
```
T_μν = ε u_μ u_ν + 𝒫 ⊥_μν + 2 u_(μ Q_ν) + 𝒯_μν,    ⊥_μν = g_μν + u_μ u_ν
```
with u^μ Q_μ = 0 (heat current transverse), 𝒯_μν transverse-traceless shear.

Constitutive relations (eqs. 17–19), barotropic p=p(e):
```
ε  = e  + τ_e [ u^μ ∇_μ e + ρ ∇_μ u^μ ]
𝒫  = p  - ζ ∇_μ u^μ + τ_p [ u^μ ∇_μ e + ρ ∇_μ u^μ ]
Q_μ = τ_𝒬 [ ρ u^ν ∇_ν u_μ + c_s² ⊥_μν ∇^ν e ]
𝒯_μν = -2 η σ_μν
```
where ρ = e + p (enthalpy density), c_s² = ∂p/∂e (sound speed squared), σ_μν = shear tensor.

Five transport coefficients: ζ (bulk visc), η (shear visc), τ_e, τ_p (energy/pressure relaxation times), τ_𝒬 (heat-current relaxation time).

Useful combinations (eq. 18):
```
V   = ζ + (4/3) η
τ_± = τ_e ± τ_p
V_± = V ± ρ τ_∓
```

**Frame choice:** comoving — the fluid 4-velocity IS one of the orthonormal frame vectors. Two-vector frame basis in the (t,r) plane:
- l_A = u_A (timelike, l^A l_A = -1)
- n_A (spacelike, n^A n_A = +1), with n_A = ε_AB u^B
Derivative shorthand: ḟ = l^A ∇_A f (along u), f' = n^A ∇_A f (radial).
Background scalar projections: v_A = ∇_A r, U = l^A v_A, W = n^A v_A, μ = ∇_A l^A, ν = ∇_A n^A.

Background stress components in frame basis (eqs. ~20):
```
t_g = p - e - μ V - (μ ρ + ė)(τ_e - τ_p)
t_p = p + e - μ V + (μ ρ + ė)(τ_e + τ_p)
t_q = -2 τ_𝒬 (ν ρ + c_s² e')
t_S = p - μ(ζ - (2/3)η) + τ_p(μ ρ + ė)
```
(V = ζ + 4η/3.)

## Axial / odd sector (Section 5) — STAGE 1B CORE
- **Gauge:** Regge-Wheeler.
- **Master fields:** two odd-parity METRIC perturbation scalars **k_n** and **k_l**. (Both are gravitational; the fluid enters only via the odd-parity stress-energy source terms ϑ_n, ϑ_l, which carry the fluid velocity perturbations.)
- λ² = ℓ(ℓ+1); GW modes ℓ≥2, fluid modes ℓ≥1.
- One master variable rescaling noted: ψ = r⁻¹ e^{Φ/2} k_n (Table 1).

The two coupled wave equations (eq. 14), as extracted:

**𝒪_l equation:**
```
k'_ṅ - k''_l + 2(μ - U) k'_n + (2W - ν) k̇_n - ν k'_l
 + (4π t_q + μ' - 2 U W) k_n
 + (4π(t_g - t_p + 4 t_S) + 2(W² + μ U - ν W) - v² - ν' + (λ²-1)/r²) k_l
 = 16π ϑ_l
```

**𝒪_n equation:**
```
k̈_n - k'̇_l + μ k̇_n - 2U k'_l + (2W - ν) k̇_l + (2 U W - 4π t_q - ν') k_l
 + (4π(2 t_g + t_p + 2 t_S) + 2(μ U - ν W - W²) + ν' + ν² - μ² + λ²/r) k_n
 = 16π ϑ_n
```
(Transcribed from ar5iv; verify exact coefficient/derivative placement against the published eq. 14 before coding — some index/derivative orderings are at the edge of OCR reliability.)

Source terms ϑ_n, ϑ_l are built from the odd-parity fluid perturbations (velocity/density variations α, β, γ, ω). 16π factor = gravitational coupling (G=c=1).

### The novel viscous mode (key physics, Stage 1B)
Verbatim from the paper:
> "One of them describes the propagation of GWs inside the star, and includes explicit dissipation effects due to the shear viscosity. The second equation corresponds to a viscous mode, with no analogue in the perfect fluid case."
> "This mode corresponds to a rotational motion in the star, which oscillates with a frequency proportional to the shear viscosity."

Interpretation:
- Mode 1: gravitational-wave propagation inside the star, **damped by shear viscosity η** (explicit dissipative terms in the wave operator).
- Mode 2: **NEW rotational (shear) fluid mode**, absent for a perfect fluid. In a perfect fluid the axial fluid sector is non-dynamical (only a stationary differential-rotation / w-mode-type structure); shear viscosity η promotes the axial rotational motion to a genuinely PROPAGATING/oscillating dynamical mode whose oscillation frequency scales with η. Analogy drawn by the author: "akin to the appearance of a second sound mode for superfluid stars."
- No closed-form dispersion relation ω(k) is printed in the paper; the existence + η-scaling is the stated result.

## Polar / even sector (for context)
Reduces to **five coupled wave equations + one constraint** (Section on polar sector). Carries GW modes plus a viscous mode analogous to the axial one.

## Causality / stability
- The paper does NOT reprint explicit inequalities. It states the characteristic/principal-part analysis of the linearized curved-background system "recovers immediately the BDNK causality constraints" of Bemfica–Disconzi–Noronha–Kovtun (Bemfica et al. 2022, PRX-type general-frame conditions).
- Net claim: theory is (i) causal given transport-coefficient constraints, (ii) linearly stable, (iii) only first-order gradients in T_μν.
- **ACTION for our build:** the explicit causality inequalities must be pulled from Bemfica et al. 2022 (general-frame conditions on τ_e, τ_p, τ_𝒬, η, ζ, c_s²), NOT from this paper. This paper only confirms they re-emerge in the stellar setting.

## Numerical method / tooling
- No PDE solver / discretization in the paper itself — the equations are presented "ready to be studied numerically for particular stellar models."
- Companion: a **Wolfram Language package + scripts** (Redondo-Yuste 2024) released to derive the perturbation equations symbolically. Useful as a cross-check for our coefficient definitions, but not a hydro evolution code.

## Relevance to our Julia BDNK NS project
- **Direct use:** the general-frame BDNK constitutive relations (eqs. 16–19) and coefficient combinations V, τ_±, V_± are exactly the constitutive forms we need to encode.
- **Stage 1B specific:** this is the reference establishing the axial 2-wave-equation structure and the existence of the η-driven rotational viscous mode — a target physical phenomenon for validation, but with no numbers to reproduce here.
- **Gaps to fill elsewhere:** explicit BDNK causality inequalities (→ Bemfica 2022), any QNM benchmark numbers (→ other refs), EoS choice (→ project decision).

## Open questions
- Exact eq. (14) coefficients need verification against the published PDF (OCR-level uncertainty in derivative orderings like k'_ṅ vs k̇'_n).
- No printed dispersion relation for the viscous mode — only the qualitative "frequency ∝ η" statement.
- Frame-basis scalars (U, W, μ, ν) defined relative to background metric functions Φ, Λ — need the explicit background TOV map to those before coding.
