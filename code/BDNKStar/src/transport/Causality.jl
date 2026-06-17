#=
    Causality — the pointwise BDNK causality / stability monitor (STEP 0
    "causality monitor operative" gate).

    The principal characteristic speeds of the BDNK system solve the biquadratic
        Λ₂ c⁴ - 2 Λ₁ c² + Λ₀ = 0,
    with coefficients (Λ₀, Λ₁, Λ₂) functions of (p, e, cs, η, ζ, τ_ε, τ_P, τ_Q, L).
    Ported verbatim from the reference solver BDNKCharacteristicSpeeds.jl
    (Keeble & Redondo-Yuste, Zenodo 19207244) — a CrossCheck-grade port.

    Causality + (weak) hyperbolicity at a point require:
      * real characteristic speeds:        Λ₁² - Λ₀ Λ₂ ≥ 0
      * non-negative speeds-squared:        c²₊ , c²₋ ≥ 0
      * subluminal:                          c²₊ ≤ 1
    `causality_flag` returns these as a NamedTuple so the evolution can monitor a
    pointwise violation flag (a first-class output, not an afterthought).
=#
module Causality

using ..Transport
using ..EquationOfState

export characteristic_speeds, causality_flag, is_causal

# Reference biquadratic coefficients (BDNKCharacteristicSpeeds.jl, lines 20–24).
Λ0(p, e, cs, η, ζ, τε, τP, τQ, L) =
    (4*L^4*(3ζ + 4η)^4*(-1 + τP)*τQ^2*τε*(p + e)^2*cs^4) / 81
Λ1(p, e, cs, η, ζ, τε, τP, τQ, L) =
    (L^2*(3ζ + 4η)^2*(τε + τQ*(τP + τε))*(p + e)*cs^2) / 9
Λ2(p, e, cs, η, ζ, τε, τP, τQ, L) =
    (2*L^2*(3ζ + 4η)^2*τQ*τε*(p + e)) / 9

"""
    characteristic_speeds(p, e, cs, tc::TransportCoefficients) -> (c2_minus, c2_plus, disc)

Squared characteristic speeds c²∓ = (Λ₁ ∓ √(Λ₁²-Λ₀Λ₂))/Λ₂ and the discriminant
`disc = Λ₁² - Λ₀Λ₂`. `cs` is the (adiabatic) sound speed = √(cs²).
"""
function characteristic_speeds(p::Real, e::Real, cs::Real, tc::TransportCoefficients)
    λ0 = Λ0(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
    λ1 = Λ1(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
    λ2 = Λ2(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
    disc = λ1^2 - λ0*λ2
    if λ2 == 0
        return (NaN, NaN, disc)
    end
    sq = disc ≥ 0 ? sqrt(disc) : NaN
    c2m = (λ1 - sq) / λ2
    c2p = (λ1 + sq) / λ2
    return (c2m, c2p, disc)
end

"""
    causality_flag(p, e, cs2, tc) -> NamedTuple

Pointwise monitor. `real_speeds` (disc ≥ 0), `nonneg` (c² ≥ 0), `subluminal`
(c²₊ ≤ 1+tol), and the overall `causal` AND of the three, plus the raw speeds.
"""
function causality_flag(p::Real, e::Real, cs2::Real, tc::TransportCoefficients;
                        tol::Real=1e-12)
    cs = sqrt(max(cs2, 0.0))
    c2m, c2p, disc = characteristic_speeds(p, e, cs, tc)
    real_speeds = disc ≥ 0
    nonneg      = (c2m ≥ -tol) && (c2p ≥ -tol)
    subluminal  = (c2p ≤ 1 + tol)
    causal = real_speeds && nonneg && subluminal
    return (causal=causal, real_speeds=real_speeds, nonneg=nonneg,
            subluminal=subluminal, c2_minus=c2m, c2_plus=c2p, disc=disc)
end

is_causal(p, e, cs2, tc::TransportCoefficients; tol=1e-12) =
    causality_flag(p, e, cs2, tc; tol=tol).causal

end # module Causality
