# arXiv:2209.09265 — Pandya, Most, Pretorius (PRD 106, 123036, 2022)
**"Causal, stable first-order viscous relativistic hydrodynamics with ideal gas microphysics"**

THE key paper for STEP 0 of the Julia BDNK neutron-star project. Provides:
the general-EOS BDNK constitutive relations, the explicit ideal-gas frame
coefficient choices, the gradient-frozen "shifted" primitive recovery, and
the causality/stability inequalities. Source TeX downloaded to
`/tmp/arxiv_src/paper.tex` (1463 lines). Companion PDE-numerics paper is
Pandya:2022pif = arXiv:2201.12317, whose simplified conformal C solver is in
`/data/haiyangw/claude/BDNK/ref-code/1D_conformal_bdnk/`.

Conventions: mostly-plus signature (-+++). Notation (Table I, eq. lines 293-307):
- energy density `epsilon`, specific internal energy `e`, pressure `P`
- `rho ≡ epsilon + P` (enthalpy density), rest-mass density `m n`
- adiabatic index `Gamma` (lit. `gamma`), baryon density `n`, 4-velocity `u^a`

---

## 1. Conservation laws & decomposition (eqs. 1-2, 7-10)

```
∇_a T^{ab} = 0            (eq:Tab_cons_law)
∇_a J^a   = 0            (eq:Ja_cons_law)
```

Decomposition w.r.t. flow velocity u^a (eq. 7-8):
```
T^{ab} = E u^a u^b + P_scr Δ^{ab} + Q^a u^b + Q^b u^a + T_scr^{ab}     (eq:Tab)
J^a    = N u^a + J_scr^a                                               (eq:Ja)
Δ^{ab} ≡ g^{ab} + u^a u^b
```
Projections (eq. 9):
```
E = u_a u_b T^{ab},   P_scr = (1/3) Δ_{ab} T^{ab},   Q^a = -Δ^{ab} u^c T_{bc}
T_scr^{ab} = T^{<ab>},  N = -u_a J^a,  J_scr^a = Δ^{ab} J_b
X^{<ab>} ≡ (1/2)(Δ^{ac}Δ^{bd} + Δ^{ad}Δ^{bc} - (2/3)Δ^{ab}Δ^{cd}) X_{cd}   (traceless transverse)
```
Ideal-fluid limit: E0=epsilon, P0=P, N0=n, Q0=T0=J0=0.

---

## 2. BDNK constitutive relations — GENERAL EOS (eqs. 11-19)

These are the build-critical constitutive relations (script quantities defined
in terms of primitives epsilon, n, u^a and their gradients):

```
E      = epsilon + τ_ε [ u^c ∇_c epsilon + ρ ∇_c u^c ]                    (eq:script_E, 11)
P_scr  = P - ζ ∇_c u^c + τ_P [ u^c ∇_c epsilon + ρ ∇_c u^c ]              (eq:script_P, 12)
Q^a    = τ_Q ρ u^c ∇_c u^a + β_ε Δ^{ac} ∇_c epsilon + β_n Δ^{ac} ∇_c n    (eq:Q_a, 13)
T_scr^{ab} = -2 η σ^{ab} ≡ -2 η ∇^{<a} u^{b>}                             (eq:script_T_ab, 14)
N      = n                                                                (eq:script_N, 15)
J_scr^a = 0                                                               (eq:script_J_a, 16)
ρ ≡ epsilon + P                                                           (eq:rho, 17)
```

The β coefficients (heat-flux microphysics) (eqs. 18-19):
```
β_ε = τ_Q p'_ε + (σ/ρ) κ_ε                                               (eq:beta_eps, 18)
β_n = τ_Q p'_n + (σ/n) κ_n                                               (eq:beta_n, 19)
```
where (eqs. 20a-20e, lines 380-385):
```
p'_ε ≡ (∂P/∂epsilon)_n
p'_n ≡ (∂P/∂n)_epsilon
κ_ε  ≡ (ρ²T/n) (∂(μ/T)/∂epsilon)_n
κ_n  ≡ ρ T (∂(μ/T)/∂n)_epsilon
κ_s  ≡ κ_ε + κ_n
```

**8 transport coefficients** that define the "hydrodynamic frame":
`τ_ε, τ_P, τ_Q` (relaxation times); `η, ζ, σ` (shear/bulk visc., thermal cond.);
`β_ε, β_n` (derived from above). 5 EOM evolve epsilon, n, 3 components of u^a;
4th component from u_c u^c = -1.

NOTE the gradient operators used throughout:
- material derivative `u^c ∇_c X` (= "dot" in rest frame)
- expansion `θ = ∇_c u^c`
- shear `σ^{ab} = ∇^{<a} u^{b>}` (traceless transverse symmetric)

---

## 3. Ideal-gas ("gamma-law") microphysics (eqs. 25-37)

EOS (eq. 25):
```
P(epsilon, n) = [Γ-1] m n e(epsilon,n) = n T(epsilon,n),     Γ ∈ (1,2)
epsilon = m n (1 + e)                                         (eq:e_defn, 26)
```
Thermo (Euler/first law) gives entropy & chemical potential (eqs. 30-31):
```
s(epsilon,n) = m n ( (1/((Γ-1)m)) ln[ e/n^{Γ-1} ] + const )
μ(epsilon,n) = m + m e ( Γ - ln[ e/n^{Γ-1} ] + const )
```
Derived microphysics derivatives (eqs. 31-34, lines 423-426):
```
p'_ε = Γ - 1
p'_n = -(Γ-1) m
κ_ε  = -(Γ-1) epsilon ρ² / (n² P)
κ_n  = (ρ/(n² P)) [ (Γ-1) epsilon² + P² ]
```
so (eqs. 35-36):
```
β_ε = (Γ-1) τ_Q - (Γ-1) σ epsilon ρ / (n² P)
β_n = -(Γ-1) m τ_Q + (σ ρ/(n³ P)) [ (Γ-1) epsilon² + P² ]
```
Sound speed (eq. 37) and auxiliaries (38-40):
```
c_s² ≡ (∂P/∂epsilon)_sbar = (Γ-1) + (n/ρ)(-(Γ-1)m) = Γ P / ρ
κ_s = -(Γ-1) m ρ / n
ω   ≡ κ_s/κ_ε = m n P /(epsilon ρ)                          (eq:omega, 39)
α   ≡ p'_ε / c_s² = (Γ-1)/c_s²                              (eq:alpha, 40)
```

---

## 4. Hydrodynamic frame choice — explicit ideal-gas coefficients (eqs. 41-44)

The "new class of frames" used in the paper (eq. 41, lines 464-469):
```
η = ρ c_s² L η̂ ,     ζ = ρ c_s² L ζ̂ ,     σ = (V̂ L ρ c_s² / (-κ_ε)) σ̂
τ_ε = τ_Q = L V̂ τ̂ ,     τ_P = 2(Γ-1) L V̂
```
with combined viscosity & inverse Reynolds number (eqs. 42-43):
```
V  ≡ (4η/3) + ζ
V̂ ≡ V/(ρ c_s² L) ≡ Re^{-1}
```
- Hatted quantities `η̂, ζ̂, σ̂, τ̂` are dimensionless free parameters.
- `L > 0` dimensionful length scale; paper sets **L = 1** throughout.
- Free: `η̂>0, ζ̂≥0, σ̂≥0`; relaxation times share single param `τ̂` which fixes
  characteristic speeds (notably max speed c_+).
- Equivalent general form (eq. 70, App. A): τ_P = 2 α c_s² L V̂  (= 2(Γ-1)LV̂ since α c_s²=Γ-1).

**Frame constraints (simplified, eq. 44, lines 505-507):**
```
σ̂ ≤ 1/3                                          (linear stability)
τ̂ ≥ [ (Γ-1)(2-c_s²) + c_s² ] / (1 - c_s²)        (causality)
```
Note c_s²→1 forces τ̂→∞ to keep speeds subluminal (frame-ansatz limitation;
neutron-star EOS with c_s²≈1 may need a different frame).

---

## 5. Causality, hyperbolicity, stability inequalities (Appendix A, eqs. 60-71)

Assumed parameter ranges for 2nd-law consistency (eq. 60, lines 1312-1318):
```
ρ, n, τ_ε, τ_P, τ_Q, η > 0 ;   m, ζ, σ ≥ 0
0 < c_s² < 1 ;   0 < ω < 3 - 2√2 ≈ 0.2 ;   1 ≤ α
```

**Causality (BDNK eq. 20 of Bemfica:2020zjp), eqs. CAUS A-D (lines 1322-1335):**
```
ρ τ_Q > η                              (CAUS A)
B² ≥ 4 A C ≥ 0                          (CAUS B)
2A > -B ≥ 0                            (CAUS C)
A > -B - C                            (CAUS D)
```
with shorthands (eqs. 61-65):
```
A = ρ τ_ε τ_Q
B = -τ_ε (ρ c_s² τ_Q + V + σ κ_s) - ρ τ_P τ_Q
C = τ_P (ρ c_s² τ_Q + σ κ_s) - β_ε V
D ≡ ρ c_s² (τ_ε + τ_Q) + V + σ κ_ε
E ≡ σ (p'_ε κ_s - c_s² κ_ε)
```

**Linear stability (BDNK eq. 48), eqs. STAB A1-E (lines 1340-1357):**
```
(τ_ε+τ_Q)|B| ≥ τ_ε τ_Q D                                         (STAB A1)
τ_ε τ_Q D ≥ ρ c_s² τ_ε τ_Q (τ_ε+τ_Q)                              (STAB A2)
(τ_ε+τ_Q)|B|D + ρ τ_ε τ_Q (τ_ε+τ_Q) E > τ_ε τ_Q D² + ρ(τ_ε+τ_Q)² C  (STAB B)
c_s² D - E ≥ ρ c_s⁴ (τ_ε+τ_Q)                                      (STAB C)
(τ_ε+τ_Q)[ |B|(c_s²D-2E) + 2 c_s² ρ τ_ε τ_Q E + C D ] >
        2 c_s² ρ (τ_ε+τ_Q)² C + τ_ε τ_Q D (c_s² D - E)             (STAB D)
|B|D[C(τ_ε+τ_Q)+E τ_ε τ_Q] + 2ρ τ_ε τ_Q (τ_ε+τ_Q) C E >
   ρ C²(τ_ε+τ_Q)² + τ_ε τ_Q (C D² + ρ τ_ε τ_Q E²) + B² E (τ_ε+τ_Q)  (STAB E)
```
Rescaled (B̂,Ĉ,D̂,Ê) versions in eqs. 66-67 reduce these to (lines 1383-1390):
```
|B̂|≥D̂ ;  D̂≥1 ;  |B̂|D̂+Ê-D̂²-Ĉ>0 ;  D̂-Ê≥1 ;  (+ two more)
```
For the frame ansatz (eq. 68: τ_ε=τ_Q=LV̂τ̂, τ_P=2α c_s² LV̂) plus σ̂≤1/3,
ALL stability inequalities are satisfied; reduces to `1-(2-αω)σ̂ ≥ 0`,
implied by `σ̂ ≤ 1/2`, strengthened to `σ̂ ≤ 1/3` for the full set.

Causality reduces (eq. 71, lines 1410-1414) to:
```
τ̂ > c_s² η/((4η/3)+ζ)
2τ̂ > c_s² (2α - ω σ̂ + τ̂ + 1)
c_s⁴(-2αω σ̂ + σ̂ + α τ̂) + τ̂² ≥ c_s² τ̂ (2α - ω σ̂ + τ̂ + 1)
```
all implied by the single eq. 44 bound `τ̂ ≥ [(Γ-1)(2-c_s²)+c_s²]/(1-c_s²)`.

**Characteristic speeds for this frame (eqs. 72-73, lines 1424-1430):**
```
c_±² = (c_s²/(2τ̂)) ( 2α - ω σ̂ + τ̂ + 1
        ± [ ω σ̂(4α+ω σ̂) + (2α+1)² - 2(ω+2)σ̂ + τ̂² + τ̂(2-2ω σ̂) ]^{1/2} )
c_1² = c_s² η/(V τ̂)
```
Saturating eq. 44 with σ̂>0 gives |c_+|<1; for strict |c_+|<1 replace ≥ with >.

Physical positivity note: τ_Q > η/(m n) (i.e. ρτ_Q>η is CAUS A) keeps P>0.

---

## 6. Primitive recovery — gradient-frozen "shifted" inversion

This paper's PDE numerics defer to Pandya:2022pif (2201.12317). The concrete
algorithm is in the conformal C solver
`/data/haiyangw/claude/BDNK/ref-code/1D_conformal_bdnk/solver.c`.
Key insight on the BDNK con2prim:

Because E, P_scr, Q^a contain TIME derivatives of primitives (via u^c∇_c),
the conserved T^{tt}, T^{tx} depend on (epsilon, u, ∂_t epsilon, ∂_t u,
∂_x epsilon, ∂_x u). The scheme does NOT do a 2D root-find for (epsilon,u) at
each cell. Instead it EVOLVES THE PRIMITIVES directly:

1. Evolve conserved T^{tt}, T^{tx} one (Heun/RK2 TVD) step from flux divergence.
2. Evolve primitives via their time derivatives:
   `xi[n+1] = xi[n] + dt*xiD[n]`, `ux[n+1] = ux[n] + dt*uxD[n]`
   (xi = log epsilon, ux = u^x).
3. Recompute spatial gradient arrays `xiP=∂_x xi, uxP=∂_x ux` (now FROZEN).
4. Recompute the primitive time derivatives by inverting the BDNK structure
   with the spatial gradients held fixed: `compute_xiD(...)`, `compute_uxD(...)`.

The inversion in `compute_xiD`/`compute_uxD`:
- Build the PERFECT-FLUID stress at current (xi,ux) WITH the frozen spatial
  gradients folded in (TttPF, TtxPF in code).
- The "shift" = mismatch between PF prediction and actual conserved value:
  `dtt = (TttPF - T00)/eta0`, `dtx = (TtxPF - T01)/eta0`
  (this is the gradient-frozen *shifted* part: viscous correction measured
  as deviation from PF).
- Solve the resulting linear system analytically for ∂_t xi, ∂_t ux as closed-
  form rational functions of (xi, ux, xiP, uxP, dtt, dtx, frame coeffs ch=chi0/eta0,
  l=lambda0/eta0). (solver.c lines 447-497 and 287-302 give compute_A, compute_Qx,
  compute_m2sxx for the constitutive pieces.)
- Fallback: if viscous shift below TOL (eta0*|dtt|<TOL), use closed-form PF
  inversion: `eps = -T00 + sqrt(4 T00² - 3 T01²)`,
  `ux = 3 T01 / sqrt((3 T00 + eps)² - (3 T01)²)`. TOL<0 ⇒ always BDNK solve.

So "primitive recovery" here = recover the TIME DERIVATIVES of the primitives
from the conserved variables with SPATIAL gradients frozen, then integrate in
time (method of lines). For a general-EOS Julia build this generalizes to a
small Newton/linear solve per cell on (∂_t epsilon, ∂_t n, ∂_t u^i) given
(T^{tt}, T^{ti}, N^t) and frozen spatial gradients. For neutron stars the
PF closed-form root may not exist analytically (general EOS), so the PF fallback
becomes a 1D/2D Newton root-find.

Code pointers (solver.c):
- `compute_xiD` (line 404), `compute_uxD` (line 455): shifted gradient-frozen inversion.
- `compute_A`/`compute_Qx`/`compute_m2sxx` (lines 287-307): constitutive E-corr, Q^x, shear.
- `T_tt`/`T_tx`/`T_xx` (lines 309-378): full BDNK conserved components.
- `T_tt0`/`T_tx0`/`T_xx0` (lines 329-402): perfect-fluid components.
- `Heun_solve_system` (line 562): predictor/corrector MoL loop, gradient recompute.
- `flux_x` (line 500): KT/CWENO central flux, KT speed a=1 (speed of light).
- WENO reconstruction `WENO_reconst_qRx/qLx` (lines 63-127); deriv `wDx`/`Dx` (139,173).

---

## 7. Numerical method (Appendix B, lines 1433-1457)

- **ODE tests** (Bjorken, shockwave ODE): explicit **RK4**, N = 2^9 … 2^13.
  Convergence factor Q_N = ||R_{N/2}||/||R_N|| → 16 for 4th order (Table III).
- **PDE solver** (from Pandya:2022pif): conservative finite-volume, **method of
  lines**, **TVD RK2 = Heun's method** in time. **WENO/CWENO** spatial recon
  (≤4th order smooth, 2nd order overall). CFL `λ = Δt/Δx = 0.1` (0.01 for stiff
  superluminal). Kurganov-Tadmor central flux with max speed a=1.
- Boundary conditions: ghost cells (outflow) or periodic; ghost-cell interaction
  degrades convergence to between 1st and 2nd order at late times.

---

## 8. Test problems & benchmark numbers

- **Bjorken flow EOM** (eq. 50, Milne coords g=diag(-1,1,1,τ²)):
  ```
  τ_ε ε̈ = -(1/τ)(τ + 2τ_ε + τ_P) ε̇ - (1/τ²)[ ρ(τ+τ_P) - V ]
  n(τ) = n_0/τ ,  u^a=(1,0,0,0)
  ```
  Inviscid: ε(τ) = m n_0 τ^{-1}[1 + e_0 τ^{-(Γ-1)}].
  τ_ε→∞: ε = c_1 τ^{-1} + c_2.
  Runs τ=1→20, n_0=0.1, Γ=4/3, m=1, τ̂∈{0.5,1,2}.
  Characteristic speeds: τ̂=0.5 → c_+≈1.3 (superluminal always);
  τ̂=1 → c_+≈1.05 early, ≈0.9 late; τ̂=2 → c_+≈0.7 (always subluminal).
  Convergence Q_N → ~16 (Table III: 16.9, 16.3, 16.1 at N=2^11).

- **Planar shockwave** (rest frame, var. in x; eq. 56 + shockwave_nprime/epsP/velP):
  Asymptotic left state {ε_L, v_L, n_L} = {1, 0.8, 0.1} at x→-∞.
  Tuned-frame condition δ ≡ β_ε ρ + β_n n - ρ c_s² τ_Q - σ κ_s = 0.
  v_L ≥ c_+ ⇒ acausal instability. Shockwave ODE convergence Q_N≈15.9 (N=2^13).

- **Heat conduction** (1+1)D telegrapher-type tests, σ̂∈{0, 0.15, 1/3}, τ̂∈{1.5,15,75},
  V̂≈2/15; convergence Q_N→4 (2nd order) until boundary interaction (~t=150).

- Frame benchmark values used: L=1, Γ=4/3, m∈{0.1,1}, V̂ (≈1/10 Bjorken,
  2/15 & 4/3 shockwaves), σ̂∈{0,1/3}, τ̂∈{0.5,1,1.5,2,3,15,75}.

---

## 9. Eckart/MIS connection & "purely frame" dynamics (Sec. III.A, lines 564-772)

- Eckart = BDNK with τ_ε=τ_P=0, τ_Q=-κT/ρ, κ ≡ σρ²/(n²T) (eq. 45-46).
- Out-of-equilibrium temperature is frame-dependent (eq. 51):
  `T = (Γ-1)/n (T^{tt} - m n) - τ_ε Ṫ`.
- Relaxation form (eq. 49): ε̇ = (1/τ_ε)(T^{tt} - ε)  ⇒  ε(t)=T^{tt}+(ε_0-T^{tt})e^{-t/τ}.
  T^{tt}=u_a u_b T^{ab}=ε+δε; the regulator terms (eqs. 52-53) are the
  perfect-fluid Euler eqs (O(∇²) on-shell), giving BDNK same physical content
  as Eckart but causal.

---

## Open questions / gaps for the Julia build
- The general-EOS (non-conformal) explicit primitive time-derivative inversion is
  NOT written out in closed form in this paper; only the conformal closed form
  exists in the C code. For a tabulated/neutron-star EOS, STEP 0 must construct
  the per-cell linear solve for (∂_t ε, ∂_t n, ∂_t u^i) symbolically or via Newton.
- Full Appendix B numerical details (exact WENO weights, ghost-cell prescription)
  live in Pandya:2022pif (2201.12317) — see the conformal C solver for the
  reference implementation (WENO_reconst_*, flux_x, Heun_solve_system).
- The frame ansatz fails as c_s²→1; neutron-star matter with stiff EOS
  (c_s²≈1) likely needs an alternative frame — flagged in eq. 44 discussion.
