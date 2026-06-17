# arXiv:2201.12317 — Pandya, Most, Pretorius: Conservative Finite-Volume BDNK

**Title:** "Conservative finite volume scheme for first-order viscous relativistic hydrodynamics" (conformal BDNK, μ=0).
This is the *method paper* behind the C reference code. Target stage: 1C-method.
Source obtained from arXiv LaTeX e-print (`paper.tex`, 1872 lines). Equations below copied exactly.

> Note: the ar5iv WebFetch summarizer mislabeled some equation numbers. The numbers below come from the actual `.tex` labels and are authoritative. The "(eqs 19-20)" in the task brief correspond to the **EOS specialization block** (`paper.tex` after `\label{eq:J_grad_correction}`) and the **frame-coefficient block** `eq:frame_coeffs_def` / `eq:frame_coeffs`.

---

## 1. BDNK stress tensor & constitutive relations (general frame)

Conserved currents (`eq:Tab_cov`, `eq:J_cov`):
```
T^{ab}_1 = (ε + 𝒜) u^a u^b + (P + Π) Δ^{ab} + 𝒬^a u^b + 𝒬^b u^a − 2 η σ^{ab}
J^a_1   = 𝒩 u^a + 𝒥^a
```
Δ^{ab} = g^{ab} + u^a u^b (mostly-plus metric).

Gradient corrections (`eq:A_cov`, then the align block):
```
𝒜   = τ_ε [ u^c ∇_c ε + (ε+P) ∇_c u^c ]
Π   = −ζ ∇_c u^c + τ_P [ u^c ∇_c ε + (ε+P) ∇_c u^c ]
𝒬^a = τ_Q (ε+P) u^c ∇_c u^a + β_ε Δ^{ac} ∇_c ε + β_n Δ^{ac} ∇_c n
```
Shear (`eq:shear_cov`):
```
σ^{ab} = ½ ( Δ^{ac}Δ^{bd}∇_c u_d + Δ^{ac}Δ^{bd}∇_d u_c − (2/3) Δ^{ab}Δ^{cd}∇_c u_d )
```
Particle-current gradient corrections are dropped (`eq:J_grad_correction`):
```
𝒩 = n ,  𝒥^a = 0
```
So J^a is identical to the ideal-fluid current.

Thermal transport coefficients:
```
β_ε = τ_Q (∂P/∂ε)_n + (σ T (ε+P)/n) (∂(μ/T)/∂ε)_n
β_n = τ_Q (∂P/∂n)_ε + (σ T (ε+P)/n) (∂(μ/T)/∂n)_ε
```
Free coefficients in general: η (shear), ζ (bulk), σ (thermal cond.), and three relaxation times τ_ε, τ_Q, τ_P.

---

## 2. Conformal EOS specialization (μ=0, g_{ab}T^{ab}=0) — "eqs 19"

```
P(ε,n) = ε/3 ,   Π = 𝒜/3 ,   ζ = 0 ,
β_ε = τ_Q/3 ,    β_n = 0 ,    τ_P = τ_ε/3
```
Leaves only η, τ_ε, τ_Q. T^{ab} becomes independent of n, so T^{ab} and J^a **decouple**.
Sound speed c_s = sqrt(∂P/∂ε) = 1/√3.

---

## 3. Frame coefficients — "eqs 20" (`eq:frame_coeffs_def`, `eq:frame_coeffs`) — USED IN C CODE

```
η      ≡ η_0 ε^{3/4}
τ_ε = (3/4ε) χ ≡ (3/4ε) χ_0 ε^{3/4}
τ_Q = (3/4ε) λ ≡ (3/4ε) λ_0 ε^{3/4}
```
Hydrodynamic frame choice (sets all characteristic speeds = 1, i.e. exactly luminal):
```
(χ_0, λ_0) = ( 25/4 η_0 , 25/7 η_0 )
```
η_0 is the single free viscosity knob. η_0 → 0 is the inviscid limit (recovers rel. Euler exactly).
This frame satisfies Bemfica et al. (2018) conditions for existence/uniqueness, causality, linear
stability about equilibrium, and Freistühler-Temple smooth-shock existence.

---

## 4. Conservative form & primitive variables

Conservation form: ∂_t q + ∂_x f^x + ∂_y f^y = ψ. For T^{ab}:
q = (T^{tt}, T^{tx}, T^{ty})^T, fluxes are the corresponding T^{xj}, T^{yj}. For J: q=J^t, f=(J^x,J^y), ψ=0.

**Perfect-fluid primitives** p_0 = (ε, u^x, u^y)^T (`eq:PF_p0`).
Analytic PF recovery (`eq:PF_pvr`) — valid because conformal + Cartesian + flat:
```
ε   = −T^{tt} + sqrt( 6 (T^{tt})^2 + 3[ (T^{tt})^2 − (T^{tx})^2 − (T^{ty})^2 ] )
|v| = sqrt((T^{tx})^2 + (T^{ty})^2) / (T^{tt} + 3ε) ,   u^t = 1/sqrt(1−|v|^2)
u^x = 3 u^t T^{tx} / (3 T^{tt} + ε) ,   u^y = 3 u^t T^{ty} / (3 T^{tt} + ε)
```

**BDNK primitives** (first-order reduction): evolve ξ ≡ ln(ε) for stability, with
p_1 = (ξ̇, u̇^x, u̇^y)^T (`eq:BDNK_p1`), plus trivial evolution ∂_t ξ = ξ̇ etc.
(Evolving ξ=ln ε instead of ε keeps state physical over (−∞,∞); avoids unphysical ε,v.)

---

## 5. BDNK primitive recovery (LINEAR — key advantage)

Because T^{ab}_1 is linear in the gradient (time-derivative) terms (`eq:BDNK_linearity_in_p`):
```
q_1 = q_0(p_0) + η_0 [ A(p_0)·p_1 + b(p_0, ∂_i p_0) ]
```
=> exact analytic inverse (`eq:naive_BDNK_pvr`):
```
p_1 = A^{-1} · [ (1/η_0)(q_1 − q_0) − b ]
```
A^{-1} exists for all physical states in the chosen frame.

**Adaptive recovery for small η_0** (truncation error τ acts like numerical viscosity, blows up as 1/η_0):
Define shifted vars (`eq:q_tilde_defn`), p_1^{PF} = perfect-fluid time-derivs c(p_0,∂_i p_0) (`eq:RE_eqns_nonconservative`):
```
q̃_1 ≡ q_1 − q_1|_{p_1→p_1^{PF}} = η_0 A·(p_1 − p_1^{PF})
p_1 = (1/η_0) A^{-1}·q̃_1 + p_1^{PF}        (eq:shifted_BDNK_pvr)
```
Algorithm per cell, with "viscous tolerance" Δ_η (empirically tuned, ↓ as resolution ↑):
1. estimate numerical viscosity → Δ_η
2. compute q̃_1 via eq:q_tilde_defn
3. if q̃_1 ≥ Δ_η: use shifted BDNK solution; update p_0 via trivial evol eqs.
   if q̃_1 < Δ_η: set q̃_1=0 (cell "inviscid"); update p_0 via perfect-fluid recovery eq:PF_pvr (trivial evol NOT used).

(Non-conformal generalization: A·C replaces η_0 A; p_1 = C^{-1}A^{-1}q̃_1 + p_1^{PF}.)

---

## 6. Reconstruction

**5th-order WENO** for p_0 at interfaces (`sec:reconstruction`, Appendix `sec:WENO`).
Right-side ENO stencils (`eq:WENO_right_stencils`):
```
v0 = −1/6 p̄_{i−2} + 5/6 p̄_{i−1} + 1/3 p̄_i
v1 =  1/3 p̄_{i−1} + 5/6 p̄_i    − 1/6 p̄_{i+1}
v2 = 11/6 p̄_i     − 7/6 p̄_{i+1} + 1/3 p̄_{i+2}
```
Left-side ENO stencils (`eq:WENO_left_stencils`):
```
u0 =  1/3 p̄_i    + 5/6 p̄_{i+1} − 1/6 p̄_{i+2}
u1 = −1/6 p̄_{i−1}+ 5/6 p̄_i     + 1/3 p̄_{i+1}
u2 =  1/3 p̄_{i−2}− 7/6 p̄_{i−1} +11/6 p̄_i
```
Smoothness indicators:
```
β0 = 1/4 (3p̄_i − 4p̄_{i+1} + p̄_{i+2})^2 + 13/12 (p̄_i − 2p̄_{i+1} + p̄_{i+2})^2
β1 = 1/4 (p̄_{i+1} − p̄_{i−1})^2          + 13/12 (p̄_{i−1} − 2p̄_i + p̄_{i+1})^2
β2 = 1/4 (p̄_{i−2} − 4p̄_{i−1} + 3p̄_i)^2  + 13/12 (p̄_{i−2} − 2p̄_{i−1} + p̄_i)^2
```
Nonlinear weights (`eq:WENO_weights`):
```
w_k = α_k / Σ_l α_l ,   α_k = d_k / (ε_W + β^k)^2 ,   d_k = (3/10, 3/5, 1/10)   [WENO]
```
Final: p^+_{i+1/2} = Σ w_k v^k ; p^-_{i+1/2} = Σ w_k u^k.
ε_W is the free WENO sensitivity parameter (ε_W→∞ ⇒ w→d ⇒ plain 5th order).

**CWENO derivative computation** for the ∂_i p_0 gradient terms appearing in BDNK fluxes (`eq:CWENO_derivative`):
```
p̄'_i = w0 (p̄_{i−2} − 4p̄_{i−1} + 3p̄_i)/(2h)
      + w1 (p̄_{i+1} − p̄_{i−1})/(2h)
      + w2 (−3p̄_i + 4p̄_{i+1} − p̄_{i+2})/(2h)
```
weighted backward/centered/forward 2nd-order stencils. Same β^k, but linear weights
`d_k = (1/6, 2/3, 1/6)` → 4th-order accurate derivative.
After PVR: compute ∂_i p_0 over whole grid via CWENO, store, then WENO-reconstruct these
derivative terms to interfaces just like p_0, then feed to flux.
(Caveat: curl constraint ∂_x∂_y ξ − ∂_y∂_x ξ = 0 only held to O(h^2) by CWENO; exact for fixed stencils.)

---

## 7. Numerical flux & time integration

**Kurganov-Tadmor central-upwind flux** (`eq:KT_flux`), example at (x_{i+1/2}, y_j):
```
F_{i+1/2,j} = ½ ( f(p^-) + f(p^+) − a [ q(p^+) − q(p^-) ] )
```
with **a = 1 exactly** (max characteristic speed = c in chosen frame). At a=1 this equals HLL =
local Lax-Friedrichs. Symmetric ⇒ discrete conservation. Same flux (a=1) used for J^a evolution.

**Time integration: Heun's method** (2nd-order TVD-RK2), q̇ = H(q):
```
q̂^{n+1} = q^n + Δt H(q^n)
q^{n+1} = q^n + (Δt/2) [ H(q^n) + H(q̂^{n+1}) ]
```
Time step Δt = λ h / a, a=1, Courant λ ∈ (0,1) (Table `table:courant`):
- 1D Gaussian: max λ 0.5, used 0.1
- 2D viscous rotor: 0.5 / 0.1
- 1D shock tube: 0.5 / 0.1
- 2D oblique shockwave: 0.1 / 0.1
- **1D steady-state shockwave: max 0.5, used 0.1**
- 2D Kelvin-Helmholtz: 0.5 / 0.5

**No artificial dissipation / Kreiss-Oliger is used** in this work (deliberate choice; the FV/WENO
scheme is non-oscillatory). KO is only mentioned as what Pandya 2021 (the semi-FD scheme) used.

**Boundaries:** outermost 3 cells = ghost cells. Ghost (`eq:ghost_cells`): X_{k,j}:=X_{3,j} for k∈[0,2];
X_{k,j}:=X_{N−4,j} for k∈[N−3,N−1]. Periodic variant identifies 3 boundary cells with 3 interior cells.
p_1 in ghost cells: either copy or set 0 (no practical difference for these tests).

---

## 8. Steady-state shock test (the requested benchmark)

Planar smooth shockwave in its rest frame (`sec` 1D steady-state shockwave). Domain x∈[−L,L].

Asymptotic states satisfy ideal-fluid Rankine-Hugoniot in rest frame (`eq:rankine_hugoniot`):
```
ε_R = ε_L (9 v_L^2 − 1) / (3(1 − v_L^2))
v_R = 1 / (3 v_L)
```
Initial data — erf bridge (`eq:steady_state_ID`):
```
ε(0,x)   = (ε_R − ε_L)/2 [ erf(x/w) + 1 ] + ε_L
v^x(0,x) = (v_R − v_L)/2 [ erf(x/w) + 1 ] + v_L ,   w = 10
```
**Chosen left state: ε_L = 1, v_L = 0.8** ⇒ via RH: v_R = 1/(3·0.8) = 0.41667,
ε_R = 1·(9·0.64 − 1)/(3·(1−0.64)) = (5.76−1)/1.08 = 4.4074.
Viscosity **η_0 = 0.2**. erf-blob settles to true BDNK static profile after a small blob leaves the domain.
Reference "true" profile = BDNK solver (eq:shifted_BDNK_pvr) applied everywhere (stable at this η_0, resolution).
Test purpose: validate the adaptive PVR — gray "non-equilibrium" region grows as Δ_η shrinks; error
vs the all-BDNK solution falls with smaller Δ_η. (Fig `steady_state_fig.pdf`.)

Related 1D shock tube (`eq:1D_shock_tube_ID`): L=200, ε=1 (x≤0)/0.1 (x>0), u^x=0, η_0=0.2,
resolutions N_x = 2^9,2^10,2^11. Convergence runs use N_x = 2^7,2^8,2^9.

---

## 9. Causality / stability conditions

Frame (χ_0,λ_0)=(25/4 η_0, 25/7 η_0) chosen to (a) make all characteristic speeds exactly 1, and
(b) satisfy the Bemfica-Disconzi-Noronha 2018/2020 inequalities guaranteeing strong hyperbolicity,
causality, thermodynamic stability (2nd law), existence/uniqueness, and linear stability about
equilibrium; also satisfies Freistühler-Temple smooth-shock existence. The paper cites these
conditions (refs Bemfica_2018, Bemfica:2020zjp) but does NOT reprint the explicit inequality set —
look there (and in arXiv:2005.11632 / 1907.12695) for the literal inequalities on χ_0, λ_0, η_0.

## Open questions for the build
- Explicit component forms of A, b, c, q_0 (relegated to Appendix `sec:coord_eqns`) not transcribed here.
- Explicit numeric value of Δ_η used per test not tabulated (tuned empirically).
- Literal causality inequalities not reprinted in this paper.
