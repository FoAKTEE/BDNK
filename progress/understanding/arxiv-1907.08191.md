# arXiv:1907.08191 — Kovtun, "First-order relativistic hydrodynamics is stable"

Source: https://arxiv.org/abs/1907.08191 (JHEP 10 (2019) 034). Priority: medium. Stage: foundational.
This paper fixes the **general-frame notation** that BDNK (Bemfica–Disconzi–Noronha) work uses. It is the linearized-stability companion; it derives the general-frame constitutive relations and the equilibrium stability/causality inequalities.

## Conventions
- **Metric signature: mostly-plus**, `g = diag(-1,+1,+1,+1)`, `u^2 = u_mu u^mu = -1`.
- **Projector**: `Delta^{ab} = g^{ab} + u^a u^b` (orthogonal to u^a). Eq. (2.3).
- **Dotted derivative**: `Ẋ = u^lambda ∂_lambda X` (material/comoving derivative). E.g. `Ṫ = u^lambda ∂_lambda T`, `µ̇ = u^lambda ∂_lambda µ`.
- d = number of spatial dims; in 3+1, `1/d` trace uses d=3 for the spatial trace in P.
- Units: ω, k measured in units of `w0/γs` (sound) or `w0/η` (shear), with `w0 ≡ ε0+p0`.

## Decomposition (Eq. 2.2, 2.3)
```
T^{µν} = E u^µ u^ν + P Δ^{µν} + (Q^µ u^ν + Q^ν u^µ) + T^{µν}     (2.2a)
J^µ   = N u^µ + J^µ                                              (2.2b)

E ≡ u_µ u_ν T^{µν},   P ≡ (1/d) Δ_{µν} T^{µν},   Q_µ ≡ -Δ_{µα} u_β T^{αβ}     (2.3a)
T_{µν} ≡ (1/2)( Δ_{µα}Δ_{νβ} + Δ_{να}Δ_{µβ} - (2/d) Δ_{µν}Δ_{αβ} ) T^{αβ}     (2.3b)
N ≡ -u_µ J^µ,   J_µ ≡ Δ_{µα} J^α                                              (2.3c)
```
T^{µν} (script-T) is the transverse, symmetric, traceless part (the shear stress).

## General-frame first-order constitutive relations (Eq. 2.4) — KEY
The first-order data: 3 scalars `Ṫ, ∂_λ u^λ, µ̇`; 3 transverse vectors `Δ^{ρσ}∂_σ T, u̇^ρ, Δ^{ρσ}∂_σ µ`; 1 TT tensor `σ^{µν}`.
Shear tensor: `σ^{µν} = Δ^{µρ}Δ^{νσ}(∂_ρ u_σ + ∂_σ u_ρ - (2/3) g_{ρσ} ∂_λ u^λ)`.

```
E   = ε  + ε1 Ṫ/T + ε2 ∂_λu^λ + ε3 u^λ∂_λ(µ/T) + O(∂²)     (2.4a)
P   = p  + π1 Ṫ/T + π2 ∂_λu^λ + π3 u^λ∂_λ(µ/T) + O(∂²)     (2.4b)
Q^µ = θ1 u̇^µ + (θ2/T) Δ^{µλ}∂_λT + θ3 Δ^{µλ}∂_λ(µ/T) + O(∂²)  (2.4c)
T^{µν} = -η σ^{µν} + O(∂²)                                  (2.4d)
N   = n  + ν1 Ṫ/T + ν2 ∂_λu^λ + ν3 u^λ∂_λ(µ/T) + O(∂²)     (2.4e)
J^µ = γ1 u̇^µ + (γ2/T) Δ^{µλ}∂_λT + γ3 Δ^{µλ}∂_λ(µ/T) + O(∂²)  (2.4f)
```
- **16 a priori coefficients** (charged): ε1,2,3; π1,2,3; θ1,2,3; ν1,2,3; γ1,2,3; η.
- **7 for uncharged**: ε1, ε2, π1, π2, θ1, θ2, η.
- The `1/T` factors are notational convenience only.
- ε1, π1, θ1, ν1, γ1 act as **relaxation times** for energy density, pressure, momentum density, charge density, charge current respectively (Discussion sec).

## Field redefinitions / frame changes (Eq. 2.8, 2.9)
General first-order field redefinition:
```
δT = a1 Ṫ/T + a2 ∂·u + a3 u^λ∂_λ(µ/T)        (2.8a)
δu^µ = b1 u̇^µ + (b2/T)Δ^{µν}∂_νT + b3 Δ^{µλ}∂_λ(µ/T)  (2.8b)
δµ = c1 Ṫ/T + c2 ∂·u + c3 u^λ∂_λ(µ/T)        (2.8c)
```
Coefficient shifts (comma = partial derivative, ε,T ≡ (∂ε/∂T)_µ):
```
εi → εi - ε,T a_i - ε,µ c_i      (2.9a)
πi → πi - p,T a_i - p,µ c_i      (2.9b)
νi → νi - n,T a_i - n,µ c_i      (2.9c)
θi → θi - (ε+p) b_i              (2.9d)
γi → γi - n b_i                  (2.9e)
η  → η                           (2.9f)
```
**Frame-invariant combinations** (Eq. 2.10) — the genuine 1-derivative transport coefficients:
```
f_i ≡ π_i - (∂p/∂ε)_n ε_i - (∂p/∂n)_ε ν_i
ℓ_i ≡ γ_i - (n/(ε+p)) θ_i
η   (invariant)
```
There are ultimately only **three genuine 1-derivative transport coefficients**: η, ζ, σ.

## Thermodynamic consistency (extensivity) — Eq. 2.11, 2.12
At O(∂): consistency forces `θ1 = θ2`, `γ1 = γ2`, hence `ℓ1 = ℓ2`. So general-frame uncharged hydro has **6** coeffs: ε1,2, π1,2, **θ ≡ θ1 = θ2**, η.
Zero-derivative: `ε = -p + T(∂p/∂T) + µ(∂p/∂µ)`, `n = ∂p/∂µ`.

## Landau-Lifshitz frame (Eq. 2.13–2.15)
Conditions: `E = ε, N = n, Q^µ = 0`. Reduces to f1,2,3, ℓ1,2,3, η. On-shell:
```
T^{µν} = ε u^µu^ν + [p - ζ ∂·u] Δ^{µν} - η σ^{µν}              (2.14a)
J^α   = n u^α - σT Δ^{αν}∂_ν(µ/T) + χT Δ^{αν}∂_ν T            (2.14b)
```
χT = 0 by consistency (ℓ1=ℓ2). Bulk viscosity / conductivity (Eq. 2.15):
```
ζ = -f2 + [ ((ε+p)(∂n/∂µ) - n(∂ε/∂µ)) f1 + (n(∂ε/∂T)_{µ/T} - (ε+p)(∂n/∂T)_{µ/T}) f3 ]
          / [ T( (∂ε/∂T)(∂n/∂µ) - (∂ε/∂µ)(∂n/∂T) ) ]          (2.15a)
σ = (n/(ε+p)) ℓ1 - (1/T) ℓ3                                    (2.15b)
```
**Uncharged bulk viscosity (footnote 8)** — directly useful for BDNK:
```
ζ = vs² (π1 - vs² ε1) - π2 + vs² ε2 ,   vs² ≡ ∂p/∂ε
```
(equivalently π2 is traded for ζ via this relation).

## Eckart frame (Eq. 2.16–2.19)
Conditions `E=ε, N=n, J^µ=0` (charged only). Heat conductivity `κ ≡ (ε+p)² σ / (n² T)`.

## Stability & causality definitions (Sec. 3)
- Stability: `Im ω(k) ⩽ 0` for all k.   (3.1)
- Causality: `lim_{k→∞} Re ω(k)/k < 1`.  (3.2)
- Dispersion written `H(w,z)=0` with `w ≡ -q·ū`, `z ≡ q² + (q·ū)²`.  (3.3)
- In a Lorentz-covariant theory, causality violation ⇒ instability. Stability must be checked at v0 ≠ 0; but stability+causality at v0=0 imply they persist at v0≠0 for gaps and large-k.

## STABILITY CONDITIONS (uncharged, general frame) — KEY DELIVERABLE

### Shear channel (Sec 4.1)
Eigenmodes are zeros of F_shear (Eq. 4.3). Gapless mode (4.4), gapped mode (4.5):
```
ω = iw0 √(1-v0²) / (η v0² - θ) + ...   (4.5),   w0 ≡ ε0+p0
```
Stability requires:
```
θ > η > 0          (4.6)
```
Landau frame sets θ=0 ⇒ unstable. Large-k shear speed solves Eq. (4.7); real and <1 for θ>η.

### Sound channel (Sec 4.2)
Eigenmodes are zeros of F_sound (Eq. 4.8, quartic in ω). Definitions:
```
vs² ≡ ∂p0/∂ε0      (speed of sound squared)
γs ≡ (4/3) η + ζ   (sound damping combination)
w0 ≡ ε0 + p0
```
Fluid at rest, small k: sound waves `ω = ±vs|k| - (i/2) γs/(ε0+p0) k²` (4.9);
gapped modes `ω = -i (ε0+p0)/(vs² ε1)`, `ω = -i (ε0+p0)/θ` (4.10).
Small-k stability:
```
γs > 0 ,  ε1 > 0 ,  θ > 0          (4.11)
```
At arbitrary k, set ω = iΔ → quartic `aΔ⁴ + bΔ³ + cΔ² + dΔ + e = 0` (4.12) with a>0,b>0,d>0.
**Routh-Hurwitz** (Re Δ < 0) ⇔ `e > 0` AND `ae < (d/b)(ac - db... )` i.e. `a e < d b (a c - d b)/b²` form. First condition:
```
ε2 + π1 > γs/vs² + vs² ε1          (4.13)
```
Second condition (with dimensionless ε̄1≡vs²ε1/γs, ε̄2≡ε2/γs, θ̄≡θ/γs, π̄1≡π1/γs):
```
(ε̄1²/vs²) + vs²(ε̄1-ε̄2)(ε̄1+θ̄)(ε̄1-π̄1)
  + (ε̄1+θ̄)[ 2ε̄1² - ε̄1(ε̄2+π̄1) + (θ̄+ε̄2)(θ̄+π̄1) ] > 0     (4.14)
```
Large-k sound speed solves quadratic (4.15) for c²_sound; require real, 0 < c²<1 (causality).
Causality of a quadratic ax²+bx+c (a>0): `b² > 4ac, 0<c<a, -c-a<b<0` (footnote 14).

### Moving fluid (v0 ≠ 0) — extra necessary condition (Eq. 4.21)
```
vs² ε1 + θ > γs/(1 - vs²)          (4.21)
```
Landau-Lifshitz uncharged (ε1=θ=0) ⇒ unstable at v0≠0.

### Summary of NECESSARY stability conditions (Discussion, uncharged)
```
θ > η > 0
vs² ε1 + θ > γs/(1 - vs²)
ε2 + π1 > γs/vs² + vs² ε1
```
Minimal regulator = 3 dimensionless params: θ/γs, ε1/γs, π1/γs.

## Conformal uncharged fluid (Sec 4, Eq. 4.22–4.23)
Conformal symmetry in 3+1 ⇒ `ε1 = 3π1, ε2 = 3π2, π1 = 3π2`, `vs = 1/√3`, `ζ = 0`.
Three independent coeffs: π1, θ, η. Stability:
```
1 - 3η/θ - η/π1 > 0 ,   π1 > 4η          (4.23)
```
Sufficient: `θ > 4η, π1 > 4η`. (matches Ref. [19].)

## Notes for BDNK build
- This is the **notation source**: BDNK's A (≡ E-ε scalar correction), P (pressure correction), Q^µ (heat-flux vector) map directly onto Kovtun's E, P, Q^µ with coefficients ε_i, π_i, θ_i.
- For an uncharged neutron-star fluid: implement Eqs (2.4a–d) with 6 coeffs {ε1,ε2,π1,π2 (or ζ), θ(=θ1=θ2), η}; use footnote-8 relation to swap π2↔ζ.
- Enforce stability via (4.6), (4.11), (4.13), (4.14), (4.21). For a quick causal/stable parameter choice use conformal-style sufficient `θ > 4η, π1 > 4η` analog (but real NS EOS is non-conformal — must check (4.13),(4.14),(4.21) numerically with actual vs²).

## Open questions
- The exact algebraic form of the Routh-Hurwitz second condition `ae < (d/b)(ac-db)` as printed (4.12 text) is slightly garbled in PDF text extraction; the consolidated invariant form is (4.14) which is unambiguous.
- BDNK papers (Bemfica-Disconzi-Noronha 2017/2019) give a specific 2-parameter subfamily (χ, λ, etc.); cross-map to be done in the BDNK-specific source extraction.
