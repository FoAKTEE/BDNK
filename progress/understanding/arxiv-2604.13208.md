# arXiv:2604.13208 — Axial Oscillations of Viscous Neutron Stars

**Authors:** S. Bussières, J. Redondo-Yuste, J. J. Ortega Gómez, V. Cardoso (submitted 2026-04-14, rev 2026-06-08)
**URL:** https://arxiv.org/abs/2604.13208 ; HTML: https://arxiv.org/html/2604.13208
**Stage relevance:** 1B — axial QNM spectrum + new viscosity-driven mode families ("η-modes").

## One-line
BDNK first-order causal viscous hydro applied to *axial* (odd-parity, shear-dominated) perturbations of TOV neutron stars; computes w-mode shifts vs viscosity AND uncovers a new η-mode family (kHz, ms damping) with no perfect-fluid counterpart, showing avoided crossings with w-modes.

## 1. BDNK stress tensor & constitutive relations
Perfect fluid (Eq.1): `T_ab^(0) = ρ u_a u_b + p Δ_ab`, `Δ_ab = g_ab + u_a u_b`.

First-order correction (Eq.2): `T_ab^(1) = A u_a u_b + Π Δ_ab + 2 u_(a Q_b) - 2η σ_ab`.

Shear (Eq.3): `σ_ab = ½ Δ_ac Δ_bd (∇^c u^d + ∇^d u^c - (2/3)Δ^cd ∇_e u^e)`.

Constitutive relations (Eq.6):
- `A = τ[ u^a ∇_a ρ + (ρ+p) ∇_a u^a ]`
- `Π = c_s² θ [ u^a ∇_a ρ + (ρ+p) ∇_a u^a ] - ζ ∇_a u^a`
- `Q_a = τ[ (ρ+p) u^b ∇_b u_a + c_s² Δ_ab ∇^b ρ ]`

Bulk viscosity ζ=0 for the axial sector (only η enters axial dynamics).

`c_s² = dp/dρ`.

## 2. Transport-coefficient parametrizations (Eq.7), L_0 scale
Two dimensionless frames, coefficients scaled by reference length L_0 (numeric value of L_0 NOT given explicitly; pure scaling param):

Parametrization A:
- `η = η̂ (ρ+p) L_0 c_s²` ; `θ = L_0 η̂` ; `τ = τ̂ L_0 η̂`

Parametrization B:
- `η = η̂ p L_0` ; `θ = L_0 p/ρ` ; `τ = τ̂ L_0 p/ρ`

η̂, τ̂ dimensionless. Four frame labels used in tables: A1, A2, B1, B2.

## 3. Causality / stability inequalities (Eq.8)
Param A (8a): `η̂ ≥ 0 , τ̂ > 0 , 0 ≤ c_s² ≤ τ̂/(2+τ̂)`.
Param B (8b): `0 ≤ η̂ ≤ 3/4 , τ̂ > max(η̂, 2/(1-c_s²))`.
Entropy (Eq.9): `∇_a S^a ≥ 2η σ_ab σ^ab / T + O(∂³)` ⇒ η ≥ 0.

## 4. Background (TOV)
`ds² = -e^ν dt² + e^λ dr² + r² dΩ²`.
- `m' = 4π r² ρ`
- `ν' = (2m + 8π r³ p)/(r(r-2m))`
- `p' = -(ρ+p)(m + 4π r³ p)/(r(r-2m))`

## 5. EOS (Eq.13)
Polytrope `p = κ ρ^(1+1/n)`, both at ρ_c = 3×10^15 g/cm³:
- EOS1: κ=100 km², n=1 → M=1.27 M☉, R=8.86 km
- EOS2: κ=700 km^2.5, n=0.8 → M=1.54 M☉, R=8.78 km
- Also constant-density (ρ=const) stars in Sec IV.3 (ultracompact, stable light ring).

## 6. Axial perturbation equations (coupled QNM system, Eqs.17-18)
Variables: ψ (Regge-Wheeler metric/fluid axial var) and Z (new viscous/shear var). `f² = e^(ν-λ)`.

```
f[ (f ψ')' ] + (ω² - V) ψ = -16π e^(ν/2) iω η ψ + C_1 Z
f[ (f Z')' ] + (c_η^{-2} ω² - U) Z = C_2 Z' + C_3 Z + C_4 ψ' + C_5 ψ
```

Regge-Wheeler potential (Eq.19): `V = e^ν[ ℓ(ℓ+1)/r² - 6m/r³ + 4π(ρ-p) ]`.
Viscous "second-sound" speed (Eq.20): `c_η² = η / [ τ (p+ρ) ]`.

Potential U and couplings C_1..C_5 (Eq.21):
- `U = e^ν[ ℓ(ℓ+1)/r² - 2m/r³ + 8π(2p+ρ) ]`
- `C_1 = (8π e^{ν-λ/2}/r²)[ 2r η' + (e^λ(1+8π r² p) - 1) η ]`
- `C_2 = (f²/(2r))[ e^λ(1+8π r² p) - 1 - 2r η'/η ]`
- `C_3 = -iω(p+ρ) e^{ν/2}(1/η + 6πτ) + 2 f² η'/(r η)`
- `C_4 = r f[ iω + (p+ρ)/η · (e^{ν/2} - iω τ) ]`
- `C_5 = f[ (p+ρ) e^{ν/2}/η - (iω/2)(-7 + e^λ(1+8π r² p)) + (iω/η)(r η' - (p+ρ)τ) ]`

Surface regularity (Eq.24): `B_1 Z(R) + B_2 Z'(R) + B_3 ψ(R) + B_4 ψ'(R) = 0`; B_i depend on frame (Eqs.25 param A, 26 param B).

## 7. Numerical method
- Interior: shoot from r_min≈0 with two regular seeds:
  `ψ^(1)=r^{ℓ+1}+..., Z^(1)=0` and `ψ^(2)=0, Z^(2)=r^{ℓ+1}+...`; impose Eq.24 to fix linear-combination coefficient K (Eq.30).
- Exterior (vacuum, ψ only): Leaver (1985) continued-fraction method. `v=1-a/r`, `ψ = χ(r)φ(v)`, `χ = (r-2M)^{2iωM} e^{iωr}`. Four-term recurrence (Eq.36) → three-term via Gaussian elimination (Eqs.38-39) → continued fraction (Eq.40) gives ψ'(a)/ψ(a).
- QNM condition: match interior/exterior log-derivatives, vanishing Wronskian `Δ(ω)=ψ_in ψ_up' - ψ_up ψ_in' = 0`; complex-ω root search with viscosity continuation tracking.

## 8. Numbers
### w-mode (ℓ=2, EOS1, ρ_c=3e15) — Table 2, format (f[kHz], τ[μs])
| η_c [g cm⁻¹ s⁻¹] | A1 | A2 | B1 | B2 |
|---|---|---|---|---|
| 3×10²⁹ | (10.4884, 29.5870) | (10.4884, 29.5870) | (10.4868, 29.5894) | (10.4868, 29.5891) |
| 1×10³⁰ | (10.4571, 29.6917) | (10.4571, 29.6898) | (10.4523, 29.6938) | (10.4522, 29.6964) |
| 1×10³¹ | (10.0898, 30.8857) | (10.0932, 30.8905) | (10.0608, 30.7400) | (10.1271, 30.8477) |

Perfect-fluid reference w-mode is the η_c→0 limit (≈10.49 kHz, τ≈29.6 μs). Viscosity lowers f and raises τ.

Compactness scaling (Eq.41): `Δf/f₀ ~ C_f - 1.8(M/R)`, `Δτ/τ₀ ~ C_τ - 5.0(M/R)` — less compact stars feel viscosity more.

### η-modes (new family, Fig.2)
- kHz frequencies, ms-scale damping times (much longer-lived than w-modes' μs damping).
- `Im ω → 0` as η_c → 0 — NO perfect-fluid counterpart (undamped/absent in inviscid limit).
- Tracked over η_c ∈ [3×10²⁹, 1×10³¹] g cm⁻¹ s⁻¹; authors could NOT follow them to very small viscosity.
- Frame-sensitive (unlike w-modes which are frame-robust).
- **Mode avoidance:** η-mode and w-mode branches approach but repel (avoided crossing) at η_c ≳ 10³⁰ g cm⁻¹ s⁻¹; destabilizes w-mode freqs in that region.
- NO standalone numeric table for η-modes (only Fig.2 complex-plane trajectories).

### Ultracompact / constant-density (Fig.3, η_c=10³¹)
Viscosity strongly damps long-lived trapped modes of the stable light ring: `|Im ω_ℓ| ≲ 10⁻²` independent of ℓ/compactness.

## Open questions / gaps
- L_0 numeric value not stated in text.
- No tabulated η-mode (f, τ) values — only Fig.2; reproduction needs digitizing or re-derivation.
- Explicit B_1..B_4 forms (Eqs.25-26) not transcribed here (frame-dependent; need PDF).
- Units/normalization of ω in figures vs Table 2 (kHz) need cross-check.
