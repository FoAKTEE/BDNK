#=
    Transport — the BDNK transport-coefficient container and named frames.

    BDNK first-order viscous hydro is parametrized by shear η, bulk ζ, heat
    conductivity κ_Q, and the hydrodynamic-frame relaxation times (τ_ε, τ_P, τ_Q)
    — the "frame" choice. Causality / stability / hyperbolicity hold only when
    these satisfy nonlinear inequalities (Bemfica–Disconzi–Noronha PRX 12 021044;
    checked numerically in BDNK_frame_constraints.nb), monitored in `Causality`.

    `{{TRANSPORT}}` (parametrized dimensionless vs physical dense-matter
    prescriptions) and `{{FRAME_SET}}` (≥2 distinct frames) are human decision
    points; this struct is the parametrized container both choices fill in.
=#
module Transport

export TransportCoefficients, conformal_frame_PMP, shum_frame_speeds, shum_frame_wellposed

"""
    TransportCoefficients(; η, ζ, κQ, τε, τP, τQ, L=1.0)

BDNK transport coefficients and relaxation times. `L` is the characteristic
gradient length scale entering the characteristic-speed analysis (set 1 for the
dimensionless/parametrized regime). All in geometrized units consistent with the
EOS (energy density km^-2).
"""
Base.@kwdef struct TransportCoefficients
    η::Float64        # shear viscosity
    ζ::Float64        # bulk viscosity
    κQ::Float64       # heat conductivity
    τε::Float64       # energy-frame relaxation time
    τP::Float64       # pressure-frame relaxation time
    τQ::Float64       # heat-flux relaxation time
    L::Float64 = 1.0  # gradient length scale
end

"""
    conformal_frame_PMP(e) -> TransportCoefficients

The conformal hydrodynamic frame used in the 1D reference solver
(Pandya, arXiv:2201.12317 eqs 19–20; 1D_conformal_bdnk/parameters.h):
η = e^{1/4}/(3π), λ = 25/7 η, χ = 25/4 η. Conformal ⇒ bulk ζ = 0 and the
coefficients scale as e^{1/4} (T³ in conformal microphysics). Returned in the
(η, ζ, κQ, τ…) container with χ→τ_Q-like and λ→τ_ε,P-like roles for cross-checks.
"""
function conformal_frame_PMP(e::Real)
    η = e^(0.25) / (3π)
    λ = (25/7) * η
    χ = (25/4) * η
    # Conformal frame: no bulk; map λ,χ into the relaxation slots for the monitor
    # (these reproduce the reference's frame; full mapping finalized vs synthesis).
    return TransportCoefficients(η=η, ζ=0.0, κQ=χ, τε=λ, τP=λ, τQ=χ, L=1.0)
end

"""
    shum_frame_speeds(ŝ, â, q̂, η̂, ζ̂, cs) -> (c0, c_plus, c_minus)

Flat-space BDNK characteristic speeds in the Shum et al. (2509.15303 eqs.67–71)
spherical-Cowling frame parametrization (V̂ = (4/3)η̂ + ζ̂):
  c0  = cs √(q̂ η̂ / (â V̂))
  c±  = cs √{ [â(1+ŝ)+q̂ ± √(q̂² + â²(4q̂+(ŝ-1)²) + 2â q̂(1+ŝ))] / (2â) }
Production frame (ŝ,â,q̂)=(1,1,0.999) gives c₊=√3 cs, c₋≈0.0183 cs.
"""
function shum_frame_speeds(ŝ::Real, â::Real, q̂::Real, η̂::Real, ζ̂::Real, cs::Real)
    V̂ = (4/3)*η̂ + ζ̂
    c0 = cs * sqrt(q̂*η̂/(â*V̂))
    disc = q̂^2 + â^2*(4q̂ + (ŝ-1)^2) + 2*â*q̂*(1+ŝ)
    base = â*(1+ŝ) + q̂
    cp = cs * sqrt((base + sqrt(disc))/(2â))
    cm = cs * sqrt((base - sqrt(disc))/(2â))
    return c0, cp, cm
end

"""Shum well-posedness + linear stability (eq.71): 0 < q̂ < ŝ."""
shum_frame_wellposed(ŝ::Real, â::Real, q̂::Real) = (0 < q̂ < ŝ)

end # module Transport
