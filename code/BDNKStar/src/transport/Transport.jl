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

export TransportCoefficients, conformal_frame_PMP

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

end # module Transport
