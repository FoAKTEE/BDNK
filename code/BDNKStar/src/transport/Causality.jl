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

    ── RECONCILIATION NOTE (vs Kovtun/Shum char-speed paths) ──────────────────
    Loose end resolved in repro/causality_reconcile.jl.  The companion paths are:
      • repro/kovtun_sound.jl   — full sound-channel dispersion ω(k); its LARGE-k
        slopes ω/k → c∞ ARE the characteristic speeds (numeric ≡ analytic, the
        eq.4.18 biquadratic `largek_sound2` in bdnk_frame_independence.jl).
      • Transport.jl::shum_frame_speeds — scale-free hatted frame, c₊²+c₋²=3cs².
    VERDICT: documented convention/system split, NOT a coefficient bug.
    Two independent witnesses (causality_reconcile.jl Tables 1–4 + the scaling
    test) show:
      (1) These (Λ₀,Λ₁,Λ₂) expect the relaxation times τε,τP,τQ as DIMENSIONLESS
          frame numbers, not the dimensionful relaxation lengths.  Evidence: the
          bare `(-1+τP)` in Λ₀ and the expansion 2Λ₁/Λ₂ = cs²(1/τQ + τP/τε + 1)
          mix dimensions — only a pure-number τ makes `(τP-1)` and the `1/τQ`
          term homogeneous.  Feeding dimensionful τ≈0.02 blows up 1/τQ≈43,
          which is exactly why frame-independence saw disc<0 for ALL frames:
          the speeds are NOT invariant under a uniform τ→λτ rescale (the scaling
          test in causality_reconcile.jl), whereas any true k→∞ slope must be.
      (2) Even fed the dimensionless ratios, this biquadratic does NOT equal
          Kovtun's large-k limit (Table 3): the correct Kovtun large-k object is
          (ε1θ)c⁴ - A2 c² + θ(cs²(ε2+π1-cs²ε1)-γs) = 0 with
          A2 = cs⁴ε1²+γsε1+(ε2+π1)(θ-cs²ε1)+ε2π1 — it carries Kovtun's π1,ε2,
          which this Keeble–Redondo-Yuste form (Λ₀∝(3ζ+4η)⁴, no π1/ε2 cross
          terms, q̂-free) lacks.  Different physical char-system, not a typo.
    ⇒ The MATCHED causality check is the Kovtun eq.4.18 route (route (i) in
      bdnk_frame_independence.jl); only the frame RATIOS map to Shum
      (ŝ=τP/(cs²τε)=ŝ, â=τQ/τε).  This module's flag is retained as the literal
      Keeble–Redondo-Yuste port for cross-reference, NOT as the definitive check.
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
