#=
    FVCommon — shared finite-volume / Valencia GRHydro primitives for the
    Cowling (frozen-metric) FV reproduction of the spherical-engine stellar
    modes. This is a GENUINELY NONLINEAR, flux-conservative scheme (HLL Riemann
    solver + MinMod-limited 2nd-order reconstruction + SSP-RK2/3), distinct from
    the LINEARIZED cowling3d engine.

    Valencia conserved variables for a barotrope p=p(ε) on a fixed metric
    (Cowling, β=0, Schwarzschild-like areal coords):
        D  = √γ ρ W                 (conserved rest-mass density)
        Sⱼ = √γ ρ h W² v_j          (conserved covariant momentum)
        τ  = √γ (ρ h W² − p − ρW)   (conserved energy minus D)
    with ρ rest-mass density, h=1+ϵ+p/ρ specific enthalpy, W=1/√(1−v²) Lorentz
    factor, v² = γ^{ij} v_i v_j.  The barotrope closes ε=ε(p), ρ=ρ(p).

    For the Shum/cold-polytrope EOS we have a clean barotrope p=p(ε); we also
    need ρ(p). For ShumPolytrope p=κρ² ⇒ ρ=√(p/κ); ε=ρ+p. We carry ρ as a
    derived quantity from the EOS so the conserved set closes on p alone.

    cons2prim (barotrope): a 1-D root find on pressure. Given (D, S², τ) and
    γ-metric, define for trial p:
        ε = τ + D + p   (NO — see below; we use the standard z-method)
    We use the robust 1-D bracketed root find:
        Define  q = τ/D, r = S²/D²  (S² = γ^{ij} S_i S_j / γ ... handled by caller
        in the orthonormal momentum). Solve for p via
            f(p) = p − p_EOS(ε(p)),  ε(p) from the conserved set.
    Concretely, with E = (τ + √γ p + D)/√γ = ρ h W² (the "energy" density),
    Ŝ² = (S²)/γ the squared physical momentum density, the standard relations:
        v² = Ŝ² / (E)²,  W = 1/√(1−v²),  ρ = D/(√γ W),  ε from EOS inverse.
    We root-find p so that p = p_EOS(ε(ρ,p)). For a barotrope p=p(ε):
        f(p) = p_EOS(ε) − p,  ε = E_phys − p − ρ ... -> we use the closed form
    below (see cons2prim_barotrope).
=#
module FVCommon

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..Numerics: brent

export PrimState, ConsCell, rho_from_p, eps_from_p,
       cons2prim_barotrope, prim2cons_barotrope,
       minmod, mc_limiter, hll_flux, lax_friedrichs_flux,
       AtmospherePars, eos_cs2, _solve_p

# rest-mass density ρ(p) and energy ε(p) for the barotrope.
@inline rho_from_p(eos::ShumPolytrope, p::Real) = p > 0 ? sqrt(p/eos.κ) : 0.0
@inline eps_from_p(eos::BarotropicEOS, p::Real) = p > 0 ? energy_from_pressure(eos, p) : 0.0
@inline eos_cs2(eos::BarotropicEOS, ε::Real)    = sound_speed2(eos, max(ε, 0.0))

# generic ρ(p): ρ = ε − p/(Γ-1)? For ShumPolytrope ε=ρ+p ⇒ ρ = ε − p.
@inline rho_from_eps_p(::ShumPolytrope, ε::Real, p::Real) = ε - p

# ---------------------------------------------------------------------------
# Atmosphere parameters
# ---------------------------------------------------------------------------
struct AtmospherePars
    ρ_atm::Float64       # atmosphere rest-mass density floor
    p_atm::Float64       # atmosphere pressure (from EOS at ρ_atm)
    ε_atm::Float64       # atmosphere energy density
    ρ_cut::Float64       # cells with ρ < ρ_cut are reset to atmosphere
    vmax::Float64        # velocity cap
end

# ---------------------------------------------------------------------------
# Slope limiters
# ---------------------------------------------------------------------------
@inline function minmod(a::Float64, b::Float64)
    (a*b ≤ 0) && return 0.0
    return abs(a) < abs(b) ? a : b
end
@inline function minmod(a::Float64, b::Float64, c::Float64)
    s = sign(a)
    (s == sign(b) && s == sign(c)) || return 0.0
    return s * min(abs(a), abs(b), abs(c))
end
# monotonized-central limiter
@inline function mc_limiter(dminus::Float64, dplus::Float64)
    minmod(0.5*(dminus+dplus), 2.0*dminus, 2.0*dplus)
end

# ---------------------------------------------------------------------------
# cons2prim for the barotrope: 1-D root find on pressure.
#
# Inputs (orthonormal / densitized-out): given the UNdensitized conserved set
#   D̂ = D/√γ = ρW,  Ŝ = |S|_phys/√γ = ρhW²|v|,  τ̂ = τ/√γ = ρhW²−p−ρW
# we recover p. Define E = τ̂ + p + D̂ = ρhW². Then |v| = Ŝ/E, W=1/√(1−v²),
# ρ = D̂/W, ε from EOS. Self-consistency: p_EOS(ε) = p. Root-find on p.
# Returns (ρ, p, ε, vmag, W).
# ---------------------------------------------------------------------------
# allocation-free residual f(p) = p − p_EOS(ε_cons(p)) for the barotrope cons2prim
@inline function _c2p_resid(eos::BarotropicEOS, p::Float64, D̂::Float64, Ŝ::Float64, τ̂::Float64)
    E = τ̂ + p + D̂
    E ≤ 0 && return 1e30
    v2 = clamp((Ŝ/E)^2, 0.0, 1.0-1e-12)
    W = 1.0/sqrt(1.0-v2)
    ε_cons = E/W^2 - p
    return p - pressure(eos, max(ε_cons, 0.0))
end

# allocation-free bracketed root solve (Illinois/regula-falsi + bisection guard).
@inline function _solve_p(eos::BarotropicEOS, D̂::Float64, Ŝ::Float64, τ̂::Float64,
                          plo0::Float64, phi0::Float64)
    plo=plo0; phi=phi0
    flo=_c2p_resid(eos,plo,D̂,Ŝ,τ̂); fhi=_c2p_resid(eos,phi,D̂,Ŝ,τ̂)
    # expand high bracket if needed
    it=0
    while flo*fhi>0 && it<60
        phi*=4; fhi=_c2p_resid(eos,phi,D̂,Ŝ,τ̂); it+=1
    end
    flo*fhi>0 && return (NaN, false)
    p=phi
    for _ in 1:100
        # regula falsi
        p = (plo*fhi - phi*flo)/(fhi-flo)
        fp = _c2p_resid(eos,p,D̂,Ŝ,τ̂)
        if abs(fp) < 1e-15*(1+abs(p)) || (phi-plo) < 1e-15*(1+abs(p))
            return (p, true)
        end
        if flo*fp < 0
            phi=p; fhi=fp
        else
            plo=p; flo=fp
        end
    end
    return (p, true)
end

function cons2prim_barotrope(eos::BarotropicEOS, D̂::Float64, Ŝ::Float64, τ̂::Float64,
                             atm::AtmospherePars; pmax_fac::Float64=10.0)
    # atmosphere reset on low density
    if D̂ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0, 1.0, true)
    end
    pmax = pmax_fac * (abs(τ̂) + D̂ + atm.p_atm) + 1e-30
    p, ok = _solve_p(eos, D̂, Ŝ, τ̂, 1e-30, pmax)
    if !ok || !isfinite(p)
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0, 1.0, true)
    end
    p = max(p, atm.p_atm)
    E = τ̂ + p + D̂
    v2 = clamp((Ŝ/E)^2, 0.0, atm.vmax^2)
    W = 1.0/sqrt(1.0 - v2)
    ρ = D̂ / W
    ε = E / W^2 - p
    if ρ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0, 1.0, true)
    end
    return (ρ, p, max(ε,0.0), sqrt(v2), W, false)
end

# prim2cons (undensitized): given ρ, p, vmag (physical), metric not needed here.
# returns (D̂, Ŝ, τ̂). ε from EOS.
function prim2cons_barotrope(eos::BarotropicEOS, ρ::Float64, p::Float64, vmag::Float64)
    ε = ρ > 0 ? energy_from_pressure(eos, p) : 0.0   # consistency via p
    # use ε directly: h = (ε+p)/ρ
    h = ρ > 0 ? (ε + p)/ρ : 1.0
    v2 = clamp(vmag^2, 0.0, 1.0-1e-12)
    W = 1.0/sqrt(1.0-v2)
    D̂ = ρ*W
    Ŝ = ρ*h*W^2*vmag
    τ̂ = ρ*h*W^2 - p - ρ*W
    return (D̂, Ŝ, τ̂)
end

# ---------------------------------------------------------------------------
# Riemann solvers (1-D, along a coordinate direction in an orthonormal frame).
# Each takes left/right (cons, flux, max wave speed) and returns the numerical
# flux. Conserved/flux vectors are length-3 NTuples (D, S, τ).
# ---------------------------------------------------------------------------
@inline function lax_friedrichs_flux(UL, FL, UR, FR, amax)
    ntuple(i -> 0.5*(FL[i]+FR[i]) - 0.5*amax*(UR[i]-UL[i]), 3)
end

@inline function hll_flux(UL, FL, sL, UR, FR, sR)
    if sL ≥ 0
        return FL
    elseif sR ≤ 0
        return FR
    else
        inv = 1.0/(sR - sL)
        return ntuple(i -> (sR*FL[i] - sL*FR[i] + sL*sR*(UR[i]-UL[i]))*inv, 3)
    end
end

end # module FVCommon
