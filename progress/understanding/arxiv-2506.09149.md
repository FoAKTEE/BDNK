# arXiv:2506.09149 — Caballero & Yunes, "Neutron star radial perturbations for causal, viscous, relativistic fluids"

- **Journal**: PRD 112, 063050 (2025). 25 pages.
- **Authors**: Daniel A. Caballero, Nicolás Yunes.
- **Source read**: ar5iv mirror (https://ar5iv.labs.arxiv.org/abs/2506.09149) + arXiv abstract.
- **Stage relevance**: 1A (radial perturbation / stability target). This is an **analytic, perturbative** paper (small transport-coefficient expansion), NOT a numerical mode-frequency catalog. No specific EOS, central density, or transport-coefficient numbers are tabulated.

## 1. What the paper is

Compares three causal/viscous closures for **radial** perturbations of neutron stars:
Eckart, BDNK (Bemfica-Disconzi-Noronha-Kovtun), Müller-Israel-Stewart (MIS, in Maxwell-Cattaneo form).
Includes bulk viscosity (ζ), shear viscosity (η), and **heat conduction (κ)** — the heat-conduction channel is the novel/dangerous one.

### Main physics results
- For small viscosity, **all three models are always radially stable to bulk and shear viscosity** (independent of the magnitude of ζ, η). The viscous operator 𝔽 is positive-definite.
- All three **can be unstable to heat conduction** unless a thermodynamic condition is satisfied. They derive a **necessary stability criterion common to all three fluids**.
- BDNK: heat-conduction perturbations show oscillatory propagation at a "second sound" speed. MIS: only decaying behavior on fast timescales.

## 2. Background (TOV) and EOS

Line element (their convention):
```
ds² = -e^{2Φ(r)} dt² + e^{2Λ(r)} dr² + r² dΩ²
```
(ar5iv rendered it as e^{-2Φ}; treat Φ sign as a convention detail — Φ is the metric potential, p' = -(ε+p)Φ' giving hydrostatic equilibrium.)

TOV / equilibrium (Eq. 32):
```
dp/dr  = -(ε+p) dΦ/dr
dΦ/dr  = (r/2)[ -e^{2Λ}(8πp + 1/r²) + 1/r² ]      (sign per hydrostatic eq.)
dΛ/dr  = (r/2)[  e^{2Λ}(8πε + 1/r²) - 1/r² ]
```

EOS written generically: p = p(n,s), ε = ε(n,s). Background is **isentropic** (s = s₀ const), so effectively p=p(ε).
First law (Eq. 5):  dε = ((ε+p)/n) dn + n T ds.

Thermodynamic derivatives (Eq. 26, 27):
```
c_s² ≡ ∂p/∂ε|_s     (adiabatic sound speed²)
c_n² ≡ ∂p/∂ε|_n     (sound speed² at fixed baryon number)
γ p  = (ε+p) c_s²     (relativistic adiabatic index)
```
Fugacity (Eq. 14): φ ≡ (ε+p)/(nT) − s  (constant in equilibrium background).

## 3. Constitutive relations (decomposition: ℰ energy, 𝒫 bulk-like, 𝒬^μ heat flux, π^μν shear, 𝒩, 𝒥^μ)

### Eckart (Eq. 15)
```
𝒫_Eck   = -ζ Θ
𝒬^μ_Eck = (κ T (ε+p)/n) Π^{μλ} ∇_λ φ
π^{μν}_Eck = -2 η σ^{μν}
ℰ_Eck = 𝒩_Eck = 𝒥^μ_Eck = 0
```

### BDNK (Eq. 13) — general (Eckart-like) frame, first-order gradients with relaxation-style τ terms
```
ℰ_BDNK   = τ_ℰ [ u^μ ∇_μ ε + (ε+p) Θ ]
𝒫_BDNK   = -ζ Θ + τ_𝒫 [ u^μ ∇_μ ε + (ε+p) Θ ]
𝒬^μ_BDNK = (κ T (ε+p)/n) Π^{μλ} ∇_λ φ + τ_𝒬 [ (ε+p) a^μ + Π^{μλ} ∇_λ p ]
π^{μν}_BDNK = -2 η σ^{μν}
𝒩_BDNK = 𝒥^μ_BDNK = 0
```
(τ_ℰ, τ_𝒫, τ_𝒬 are BDNK frame/transport coefficients; a^μ = u·∇u^μ acceleration, Θ = ∇·u expansion, σ^{μν} shear tensor, Π^{μν} spatial projector.)

### MIS / Maxwell-Cattaneo (Eq. 17) — relaxation equations
```
τ_0 u^ν ∇_ν 𝒫_MC + 𝒫_MC = -ζ Θ
τ_1 Π^μ_ν u^λ ∇_λ 𝒬^ν_MC + 𝒬^μ_MC = (κ T (ε+p)/n) Π^{μλ} ∇_λ φ
τ_2 Π^μ_α Π^ν_β u^λ ∇_λ π^{αβ}_MC + π^{μν}_MC = -2 η σ^{μν}
ℰ_MC = 𝒩_MC = 𝒥^μ_MC = 0
```

## 4. Perturbation equations

### Zeroth-order (ideal) Chandrasekhar radial pulsation, Sturm-Liouville form
Renormalized variable: **Ξ ≡ r² e^Φ ξ**, where ξ is the Lagrangian radial displacement.

```
W(r) ∂_t² Ξ = ∂_r( P(r) ∂_r Ξ ) - Q(r) Ξ        (zeroth-order, ideal)
```
Coefficients (Eq. 40a–40c):
```
W(r) = e^{3Λ-Φ} (ε+p) / r²
P(r) = e^{Λ-3Φ} γ p / r²                  ( = e^{Λ-3Φ}(ε+p)c_s²/r² )
Q(r) = -e^{Λ-3Φ} (ε+p)/r² [ Φ'' - (2/r)Φ' - (Φ')² ]
        + 4π e^{3Λ-3Φ} (ε+p)² / r · ( Φ' - 1/r )
```
Eigenproblem: with Ξ ~ e^{iωt}, L₀ Ξ = ω² Ξ defining real ω²_j, eigenfunctions φ_j^{(0)}.

### Viscous corrections
Per-model perturbation equations: Eq. 46 (Eckart), Eq. 52 (BDNK), Eq. 54 (MIS).
Schematically the full operator adds viscous (𝔽) and heat-conduction (𝔾) pieces to the ideal operator L₀.

## 5. Perturbative stability analysis (small α: (ζ,η,κ)→(αζ,αη,ακ), α≪1)

First-order complex frequency shift (Eq. 62):
```
ω_j^{(1)} = (i/2) ⟨φ_j^{(0)}| 𝔽 |φ_j^{(0)}⟩ / ⟨φ_j^{(0)}|φ_j^{(0)}⟩
          + (i/2) (1/ω_j^{(0)2}) ⟨φ_j^{(0)}| 𝔾 |φ_j^{(0)}⟩ / ⟨φ_j^{(0)}|φ_j^{(0)}⟩
```
Stability ⇔ Im(i ω_j^{(1)}) < 0 (mode decays). 𝔽 = viscous (bulk+shear) part; 𝔾 = heat-conduction part.

### Viscous operator 𝔽 (Eq. 64–65) — positive definite ⇒ bulk/shear always stable
```
W(r) 𝔽 ϕ = -∂_r[ (ζ + 4η/3) e^{Λ-2Φ} ∂_r ϕ / r² ]
            + 4 η e^{Λ-2Φ} (Λ' - 2Φ') ϕ / r³
            - 4 (∂_r η) e^{Λ-2Φ} ϕ / r³
            - (κ T (ε+p)/n) ∂_r[ ((ε+p)/(nT)) (c_s² - c_n²) e^{-Φ} ∂_r ϕ / r² ]
```
(The κ term appears here through Y[ϕ] below; the pure ζ,η part is positive-definite, giving **unconditional bulk/shear stability** — Eq. 69.)

### Heat-conduction source and the stability criterion
Heat-conduction contribution carried by
```
Y[ϕ] = -(κ T (ε+p)/n) ∂_r[ ((ε+p)/(nT)) (c_s² - c_n²) e^{-Φ} ∂_r ϕ / r² ]
```
Stability against heat conduction requires the heat-conduction expectation value to be non-destabilizing, which reduces to a condition on the sign of **(c_s² − c_n²)**.

**Heat-conduction stability criterion (necessary, common to Eckart/BDNK/MIS):**
> the fluid is stable to heat conduction when **c_s² − c_n² ≥ 0** (adiabatic sound speed ≥ fixed-baryon-number sound speed); if c_n² > c_s² over the star, heat conduction can drive a radial instability.
The precise statement is an eigenfunction-weighted integral of (c_s²−c_n²) (via 𝔾) being non-positive; locally c_s² ≥ c_n² guarantees it. This is the inequality to implement/check in Stage 1A.

## 6. Numbers / benchmarks

- **No specific EOS** (no polytrope index, no SLy/APR) — analysis is EOS-generic.
- **No central density / pressure values**, no κ,ζ,η magnitudes, **no tabulated mode frequencies or growth rates**. The work is analytic + first-order perturbative; α is a formal small parameter, never assigned a number.
- For reproduction we must supply our own EOS + TOV solve, build W,P,Q, solve the ideal Chandrasekhar eigenproblem (ω²_j, φ_j), then evaluate 𝔽, 𝔾 matrix elements to get ω^{(1)}.

## 7. Open questions for build
- Exact forms of BDNK τ_ℰ, τ_𝒫, τ_𝒬 and MIS τ_0,τ_1,τ_2 in terms of (ε,p,n,T) — need full Eq. 13/17 context + their causality/stability constraints (likely cross-reference Kovtun 2019, Bemfica et al.).
- Exact 𝔾 operator (Eq. 65) full expression and the precise integral form of the stability functional (only the (c_s²−c_n²) dependence is confirmed here).
- Metric signature/Φ-sign convention vs our TOV solver (ar5iv showed e^{-2Φ} in g_tt; reconcile against p'=-(ε+p)Φ').
- BDNK "second sound" speed expression for the oscillatory heat mode (numeric check target if we add a benchmark).
