#=
    EquationOfState — the shared EOS driver (STEP 0 trunk).

    One abstracted module reused by every later stage (STAGE 1A radial, 1C
    nonlinear, 2 realistic-EOS recovery, 3 polytropic 3+1D recovery). The plan
    (bdnk_hmns_plan.tex) makes this the highest-leverage decision in the project:
    the two "not done" recovery items are this module with different inputs.

    Two EOS families share one interface:

      BarotropicEOS  — p = p(e)            (cold stars; the reference polytrope
                                            p = κ e^{1+1/n})
      GeneralEOS     — p = p(ρ, ϵ)         (finite-T; ideal-gas microphysics,
                                            the BDNK heat-conduction sector)

    where  ρ  = rest-mass density,  ϵ  = specific internal energy,
           e  = ρ(1+ϵ) total energy density,  h = 1 + ϵ + p/ρ specific enthalpy.

    Every EOS exposes thermodynamically consistent derivatives — the sound speed
    cs² (which sets the characteristic speeds and hence the BDNK causality
    inequalities) and the partials ∂p/∂ρ|_ϵ, ∂p/∂ϵ|_ρ used by the recovery and
    flux Jacobians. Monotonicity (cs² > 0) and a density floor / atmosphere are
    enforced. Tabulated EOS are represented by a monotone cubic-Hermite table
    whose recovery is *convergent* under refinement (the real {{EOS_TABLE}} is a
    downstream human decision; the interface is fixed here).
=#
module EquationOfState

using ..Numerics

export AbstractEOS, BarotropicEOS, GeneralEOS,
       PolytropeEnergy, IdealGas, TabulatedBarotrope, tabulate,
       pressure, sound_speed2, energy_from_pressure,
       dpdrho_eps, dpdeps_rho, specific_enthalpy, total_energy_density,
       temperature, is_thermodynamically_valid, apply_floor

abstract type AbstractEOS end
abstract type BarotropicEOS <: AbstractEOS end   # p = p(e)
abstract type GeneralEOS   <: AbstractEOS end     # p = p(ρ, ϵ)

# ---------------------------------------------------------------------------
# Cold energy polytrope  p = κ e^{1+1/n}   (the NeutronStarOscillations.jl EOS)
# ---------------------------------------------------------------------------
struct PolytropeEnergy <: BarotropicEOS
    κ::Float64
    n::Float64
end

@inline pressure(eos::PolytropeEnergy, e::Real) = eos.κ * e^(1 + 1/eos.n)
# cs² = dp/de = (1 + 1/n) κ e^{1/n} = (1 + 1/n) p/e
@inline sound_speed2(eos::PolytropeEnergy, e::Real) = (1 + 1/eos.n) * eos.κ * e^(1/eos.n)
@inline energy_from_pressure(eos::PolytropeEnergy, p::Real) = (p / eos.κ)^(eos.n/(eos.n + 1))

# ---------------------------------------------------------------------------
# Γ-law ideal gas  p = (Γ-1) ρ ϵ   (finite temperature; BDNK ideal-gas micro-
# physics, cf. Pandya–Most–Pretorius 2209.09265). Temperature in units k_B/m=1.
# ---------------------------------------------------------------------------
struct IdealGas <: GeneralEOS
    Γ::Float64
end

@inline pressure(eos::IdealGas, ρ::Real, ϵ::Real)   = (eos.Γ - 1) * ρ * ϵ
@inline dpdrho_eps(eos::IdealGas, ρ::Real, ϵ::Real) = (eos.Γ - 1) * ϵ
@inline dpdeps_rho(eos::IdealGas, ρ::Real, ϵ::Real) = (eos.Γ - 1) * ρ
@inline temperature(eos::IdealGas, ρ::Real, ϵ::Real) = (eos.Γ - 1) * ϵ   # T = p/ρ

# total energy density and specific enthalpy (generic for any GeneralEOS)
@inline total_energy_density(::GeneralEOS, ρ::Real, ϵ::Real) = ρ * (1 + ϵ)
@inline function specific_enthalpy(eos::GeneralEOS, ρ::Real, ϵ::Real)
    return 1 + ϵ + pressure(eos, ρ, ϵ) / ρ
end

# Relativistic sound speed:  cs² = (1/h)(∂p/∂ρ|_ϵ + (p/ρ²) ∂p/∂ϵ|_ρ)
@inline function sound_speed2(eos::GeneralEOS, ρ::Real, ϵ::Real)
    p = pressure(eos, ρ, ϵ)
    h = 1 + ϵ + p/ρ
    return (dpdrho_eps(eos, ρ, ϵ) + (p/ρ^2) * dpdeps_rho(eos, ρ, ϵ)) / h
end

# ---------------------------------------------------------------------------
# Tabulated barotrope — monotone cubic-Hermite table of p(e), with cs² = dp/de
# read off the same interpolant (thermodynamically consistent by construction).
# `tabulate(base, e_lo, e_hi, N)` builds a log-spaced table from any barotrope.
# ---------------------------------------------------------------------------
struct TabulatedBarotrope <: BarotropicEOS
    loge::Vector{Float64}   # log e nodes (ascending)
    p::Vector{Float64}      # p at nodes
    dpde::Vector{Float64}   # dp/de at nodes (= cs²)
end

"""
    tabulate(base::BarotropicEOS, e_lo, e_hi, N) -> TabulatedBarotrope

Sample a base barotrope on `N` log-spaced energy-density nodes in `[e_lo, e_hi]`,
storing p and dp/de at the nodes for cubic-Hermite reconstruction. Refining `N`
drives the tabulated recovery error to zero (STEP 0 "tabulated convergent").
"""
function tabulate(base::BarotropicEOS, e_lo::Real, e_hi::Real, N::Int)
    loge = collect(range(log(e_lo), log(e_hi); length=N))
    e = exp.(loge)
    p = [pressure(base, ei) for ei in e]
    dpde = [sound_speed2(base, ei) for ei in e]
    return TabulatedBarotrope(loge, p, dpde)
end

@inline function _locate(t::TabulatedBarotrope, le::Float64)
    lo, hi = first(t.loge), last(t.loge)
    le = clamp(le, lo, hi)
    # binary search for interval index i with loge[i] ≤ le ≤ loge[i+1]
    a, b = 1, length(t.loge)
    while b - a > 1
        m = (a + b) >>> 1
        (t.loge[m] ≤ le) ? (a = m) : (b = m)
    end
    return a, le
end

function pressure(t::TabulatedBarotrope, e::Real)
    i, le = _locate(t, log(e))
    Δ = t.loge[i+1] - t.loge[i]
    s = (le - t.loge[i]) / Δ
    # Hermite in log e: store derivatives wrt log e -> d(p)/d(loge) = (dp/de)*e
    e_i  = exp(t.loge[i]);  e_i1 = exp(t.loge[i+1])
    m0 = t.dpde[i]   * e_i  * Δ
    m1 = t.dpde[i+1] * e_i1 * Δ
    h00 = (2s^3 - 3s^2 + 1); h10 = (s^3 - 2s^2 + s)
    h01 = (-2s^3 + 3s^2);    h11 = (s^3 - s^2)
    return h00*t.p[i] + h10*m0 + h01*t.p[i+1] + h11*m1
end

function sound_speed2(t::TabulatedBarotrope, e::Real)
    i, le = _locate(t, log(e))
    Δ = t.loge[i+1] - t.loge[i]
    s = (le - t.loge[i]) / Δ
    e_i  = exp(t.loge[i]);  e_i1 = exp(t.loge[i+1])
    m0 = t.dpde[i]   * e_i  * Δ
    m1 = t.dpde[i+1] * e_i1 * Δ
    dh00 = (6s^2 - 6s); dh10 = (3s^2 - 4s + 1)
    dh01 = (-6s^2 + 6s); dh11 = (3s^2 - 2s)
    dp_dloge = (dh00*t.p[i] + dh10*m0 + dh01*t.p[i+1] + dh11*m1) / Δ
    return dp_dloge / e                      # dp/de = (dp/d loge)/e
end

function energy_from_pressure(t::TabulatedBarotrope, p::Real)
    e_lo = exp(first(t.loge)); e_hi = exp(last(t.loge))
    r = brent(e -> pressure(t, e) - p, e_lo, e_hi; xtol=1e-14)
    return r.root
end

# ---------------------------------------------------------------------------
# Validity / monotonicity and atmosphere floor
# ---------------------------------------------------------------------------
"""
    is_thermodynamically_valid(eos, args...) -> Bool

Monotone (cs² > 0) and subluminal-or-marginal (cs² ≤ 1) — the minimal physical
admissibility used to flag bad table regions / stellar-surface atmosphere.
"""
is_thermodynamically_valid(eos::BarotropicEOS, e::Real) =
    (cs2 = sound_speed2(eos, e); cs2 > 0 && cs2 ≤ 1 + 1e-12)
is_thermodynamically_valid(eos::GeneralEOS, ρ::Real, ϵ::Real) =
    (cs2 = sound_speed2(eos, ρ, ϵ); ρ > 0 && cs2 > 0 && cs2 ≤ 1 + 1e-12)

@inline apply_floor(e::Real, floor::Real) = max(e, floor)

end # module EquationOfState
