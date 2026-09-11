#=
    IsraelStewart — second-order Müller–Israel–Stewart (MIS) bulk-viscous
    closure with the Chabanov–Rezzolla causal limiter (arXiv:2311.13027).

    This is the IS half of the STAGE-2 IS-vs-BDNK bulk-viscous contrast (DAG
    node s2.is_contrast): the *1D-reproducible* physics of the otherwise
    compute-infeasible 3+1D GRMHD merger paper. Where BDNK (`Causality`)
    enforces causality through frame-constraint inequalities on the relaxation
    times, the MIS scheme keeps the viscous sound speed subluminal by a pointwise
    limiter that rescales τ_Π — that mechanism contrast is the whole point.

    Maxwell–Cattaneo relaxation (eq:is_bulk / eq:nsvalue):
        τ_Π u·∇Π = Π_NS − Π,   Π_NS = −ζ Θ,   Θ = ∇·u,
    Π an INDEPENDENT evolved variable; the extra conserved density is
        DΠ = ρ W Π                                         (eq:bulk_implement1).

    Effective viscous sound speed (eq:limit_sound, h' = h + Π/ρ):
        c_s'^2 = (ζ/τ_Π)(1/(ρ h')) + (∂p/∂e)_ρ + (1/h')(∂p/∂ρ)_e.
    Causal limiter (iii) (eq:causal_limit): if c_s'^2 > c_max^2 reset
        τ_Π = (ζ/(ρ h')) [ c_max^2 − (∂p/∂e)_ρ − (1/h')(∂p/∂ρ)_e ]^{-1},
    which drives c_s'^2 → c_max^2 EXACTLY (0 ≤ c_max^2 < 1, free parameter).

    Inverse Reynolds number R^{-1} = Π/(p+e) (eq:reynolds); its ideal-EOS
    central-density bounds (Sec. limiting), for p = κ ρ^Γ:
        R^{-1}_min = σ p/(e+p),   R^{-1}_max = (e−p)/(e+p).

    All line/equation labels refer to bulk_vis_in_bns_simulations.tex via
    progress/reproduction/2311.13027.md.
=#
module IsraelStewart

export bulk_pressure_NS, enthalpy_prime, cs2_viscous, tauPi_causal,
       apply_causality_fix, reynolds_inv, reynolds_inv_min, reynolds_inv_max,
       DPi_conserved

"""Navier–Stokes bulk pressure Π_NS = −ζ Θ (eq:nsvalue)."""
@inline bulk_pressure_NS(ζ::Real, Θ::Real) = -ζ * Θ

"""h' = h + Π/ρ — specific enthalpy including bulk pressure (line 481)."""
@inline enthalpy_prime(h::Real, Π::Real, ρ::Real) = h + Π / ρ

"""
    cs2_viscous(cs2_eq, ζ, τΠ, ρ, hp)

Effective MIS viscous sound speed squared c_s'^2 (eq:limit_sound). `cs2_eq =
(∂p/∂e)_ρ + (1/h')(∂p/∂ρ)_e` is the equilibrium (Π→0) value — for a cold
barotrope just the adiabatic c_s². `hp = h' = h + Π/ρ`.
"""
@inline cs2_viscous(cs2_eq::Real, ζ::Real, τΠ::Real, ρ::Real, hp::Real) =
    (ζ / τΠ) * (1 / (ρ * hp)) + cs2_eq

"""
    tauPi_causal(cmax2, cs2_eq, ζ, ρ, hp)

Relaxation time that makes c_s'^2 == cmax2 exactly (eq:causal_limit). Requires
`cmax2 > cs2_eq` (the equilibrium part must itself be subluminal).
"""
@inline function tauPi_causal(cmax2::Real, cs2_eq::Real, ζ::Real, ρ::Real, hp::Real)
    0 ≤ cmax2 < 1 || error("IS limiter: c_max²=$cmax2 must satisfy 0 ≤ c_max² < 1 (subluminal target)")
    cmax2 > cs2_eq || error("IS limiter ill-posed: equilibrium c_s²=$cs2_eq ≥ c_max²=$cmax2")
    return (ζ / (ρ * hp)) / (cmax2 - cs2_eq)
end

"""
    apply_causality_fix(cs2_eq, ζ, τΠ, ρ, hp; cmax2=0.9) -> (τΠ_used, cs2_after, fixed)

Limiter (iii): rescale τ_Π only where the naive c_s'^2 would exceed `cmax2`;
otherwise leave τ_Π untouched. Returns the τ_Π actually used, the resulting
c_s'^2, and whether the fix fired.
"""
function apply_causality_fix(cs2_eq, ζ, τΠ, ρ, hp; cmax2=0.9)
    cs2_naive = cs2_viscous(cs2_eq, ζ, τΠ, ρ, hp)
    if cs2_naive > cmax2
        τf = tauPi_causal(cmax2, cs2_eq, ζ, ρ, hp)
        return (τf, cs2_viscous(cs2_eq, ζ, τf, ρ, hp), true)
    else
        return (τΠ, cs2_naive, false)
    end
end

"""Inverse Reynolds number R^{-1} = Π/(p+e) (eq:reynolds, Denicol2012b)."""
@inline reynolds_inv(Π::Real, p::Real, e::Real) = Π / (p + e)

"""Lower bound R^{-1}_min = σ p/(e+p) (Sec. limiting, ideal-EOS estimate)."""
@inline reynolds_inv_min(σ::Real, p::Real, e::Real) = σ * p / (e + p)

"""Upper bound R^{-1}_max = (e−p)/(e+p) (Sec. limiting, ideal-EOS estimate)."""
@inline reynolds_inv_max(p::Real, e::Real) = (e - p) / (e + p)

"""Extra IS conserved variable DΠ = ρ W Π (eq:bulk_implement1, the variable
that distinguishes true second-order IS evolution from first-order BDNK)."""
@inline DPi_conserved(ρ::Real, W::Real, Π::Real) = ρ * W * Π

end # module IsraelStewart
