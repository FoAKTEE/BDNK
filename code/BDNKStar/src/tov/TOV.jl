#=
    TOV — Tolman–Oppenheimer–Volkoff background star (the shared static
    background for STAGE 1A/1B/1C linear + nonlinear analyses).

    Geometrized G=c=1, lengths km (Units). Barotropic EOS ε=ε(p). Equations
    (matching NeutronStarOscillations.jl TOV.jl):

        m'(r) = 4π r² ε
        ν'(r) = 2(m + 4π r³ p) / [r(r-2m)]
        p'(r) = -(ε+p)(m + 4π r³ p) / [r(r-2m)]

    Integrated outward by RK4 from a regular center (Taylor-seeded at r=h):
        m(h)=(4π/3)ε_c h³,  p(h)=p_c-2π(ε_c+p_c)(ε_c/3+p_c)h².
    Surface R at p→0 (linear interpolation of the last step); M=m(R). The metric
    potential ν is shifted so e^{ν(R)} = 1-2M/R (Schwarzschild exterior match).
=#
module TOV

using ..EquationOfState

export TOVStar, solve_tov, mass_solar, IdealGasStar, solve_tov_idealgas

struct TOVStar
    r::Vector{Float64}
    m::Vector{Float64}
    p::Vector{Float64}
    ε::Vector{Float64}
    ν::Vector{Float64}
    M::Float64      # gravitational mass [km]
    R::Float64      # areal radius [km]
end

@inline function _rhs(eos::BarotropicEOS, r, m, p)
    ε = p > 0 ? energy_from_pressure(eos, p) : 0.0
    denom = r * (r - 2m)
    fac = m + 4π * r^3 * p
    dm = 4π * r^2 * ε
    dν = 2 * fac / denom
    dp = -(ε + p) * fac / denom
    return dm, dp, dν, ε
end

"""
    solve_tov(eos, εc; h=1e-3, ptol_rel=1e-10, rmax=100.0) -> TOVStar

Integrate the TOV star with central energy density `εc` (km⁻²) under barotropic
`eos`. `h` is the RK4 radial step (km); termination at p ≤ ptol_rel·p_c.
"""
function solve_tov(eos::BarotropicEOS, εc::Float64; h::Float64=1e-3,
                   ptol_rel::Float64=1e-10, rmax::Float64=100.0)
    pc = pressure(eos, εc)
    ptol = ptol_rel * pc
    # regular-center Taylor seed at r = h
    r0 = h
    m = (4π/3) * εc * r0^3
    p = pc - 2π*(εc + pc)*(εc/3 + pc) * r0^2
    ν = 0.0
    rs = Float64[r0]; ms = Float64[m]; ps = Float64[p]
    εs = Float64[energy_from_pressure(eos, p)]; νs = Float64[ν]
    r = r0
    while p > ptol && r < rmax
        dm1, dp1, dν1, _ = _rhs(eos, r,       m,            p)
        dm2, dp2, dν2, _ = _rhs(eos, r+h/2,   m+h/2*dm1,    max(p+h/2*dp1, 0.0))
        dm3, dp3, dν3, _ = _rhs(eos, r+h/2,   m+h/2*dm2,    max(p+h/2*dp2, 0.0))
        p4 = max(p + h*dp3, 0.0)
        dm4, dp4, dν4, _ = _rhs(eos, r+h,     m+h*dm3,      p4)
        m_new = m + h/6*(dm1 + 2dm2 + 2dm3 + dm4)
        p_new = p + h/6*(dp1 + 2dp2 + 2dp3 + dp4)
        ν_new = ν + h/6*(dν1 + 2dν2 + 2dν3 + dν4)
        r += h
        if p_new ≤ ptol
            # linear interpolation to the surface p = 0 (between p and p_new)
            frac = p / (p - p_new)
            R = (r - h) + frac * h
            M = m + frac * (m_new - m)
            push!(rs, R); push!(ms, M); push!(ps, 0.0)
            push!(εs, 0.0); push!(νs, ν + frac*(ν_new - ν))
            # shift ν to Schwarzschild exterior at the surface
            νsurf = log(1 - 2M/R)
            @. νs += -νs[end] + νsurf
            return TOVStar(rs, ms, ps, εs, νs, M, R)
        end
        m, p, ν = m_new, p_new, ν_new
        push!(rs, r); push!(ms, m); push!(ps, p)
        push!(εs, energy_from_pressure(eos, p)); push!(νs, ν)
    end
    error("TOV integration did not reach the surface within rmax=$rmax km")
end

# ---------------------------------------------------------------------------
# Finite-temperature ideal-gas star
# ---------------------------------------------------------------------------
#=
    A STATIC star's structure is necessarily barotropic, so we integrate the
    hydrostatic structure with the ideal-gas-on-an-adiabat polytrope

        p = K ρ^{Γ_struct},     ε = ρ + p/(Γ_struct - 1)             (structure)

    (`isentropic_idealgas` with index Γ_struct), then RECONSTRUCT the genuine
    finite-temperature ideal-gas profiles along the resulting star using the
    micro-EOS  p = (Γ-1) ρ ϵ  with the *adiabatic* (perturbation) index Γ:

        ρ(r) = (p(r)/K)^{1/Γ_struct}         (invert the structure polytrope)
        ϵ(r) = p(r) / [(Γ-1) ρ(r)]           (ideal-gas internal energy)
        T(r) = p(r)/ρ(r) = (Γ-1) ϵ(r)        (temperature, units k_B/m = 1)
        cs²(r) = sound_speed2(IdealGas(Γ), ρ, ϵ)   (adiabatic ∂p/∂ε|_s)
        cn²(r) = cn2(IdealGas(Γ), ρ, ϵ) = Γ-1      (fixed-baryon ∂p/∂ε|_n)

    CLOSURE CHOICE (documented):
      * Γ_struct == Γ  → ISENTROPIC star: the hydrostatic adiabat coincides with
        the perturbation adiabat. This is the honest, self-consistent ideal-gas
        star. The Caballero–Yunes difference is STILL nonzero on it:
        cs²−cn² = −(Γ−1)/(1+Γϵ) < 0 (the ideal gas always violates the heat-
        conduction criterion), so the heat-conduction sector is genuinely
        exercised — it is *not* marginal as a cold barotrope would be.
      * Γ_struct ≠ Γ  → composition/entropy STRATIFIED star: the structure
        polytrope and the perturbation adiabat differ. T(r) ∝ ρ^{Γ_struct-1}
        is then a non-isentropic temperature profile; this is the closure that
        (in a full treatment) activates buoyancy g-modes. cs²−cn² has the same
        sign but a different radial profile.

    In both cases T(r) decreases outward (dT/dr < 0), so heat conduction has a
    real temperature gradient to act on. This delivers the finite-T BACKGROUND
    and the heat-conduction sector (criterion + causality); the full thermal-
    perturbation (δT / heat-flux) mode equations are a further extension.
=#
struct IdealGasStar
    r::Vector{Float64}      # areal radius [km]
    ρ::Vector{Float64}      # rest-mass density [km^-2]
    p::Vector{Float64}      # pressure [km^-2]
    ε::Vector{Float64}      # total energy density ρ(1+ϵ) [km^-2]
    ϵ::Vector{Float64}      # specific internal energy [dimensionless]
    T::Vector{Float64}      # temperature p/ρ [units k_B/m = 1]
    cs2::Vector{Float64}    # adiabatic sound speed squared
    cn2::Vector{Float64}    # equilibrium (fixed-baryon) sound speed squared
    m::Vector{Float64}      # enclosed mass [km]
    ν::Vector{Float64}      # metric potential
    M::Float64              # gravitational mass [km]
    R::Float64              # areal radius [km]
    Γ::Float64              # adiabatic / perturbation index
    Γ_struct::Float64       # structure polytrope index
    K::Float64              # polytropic constant
end

"""
    solve_tov_idealgas(; Γ=5/3, Γ_struct=Γ, K=100.0, ρc, h=2e-4,
                       N_tab=1200, kwargs...) -> IdealGasStar

Build a finite-temperature ideal-gas TOV star and reconstruct the finite-T
profiles along it. The hydrostatic structure is integrated with the ideal-gas-
on-an-adiabat polytrope p = K ρ^{Γ_struct} (a barotrope, as a static star must
be), then the ideal-gas micro-EOS p = (Γ-1)ρϵ with the *adiabatic* index Γ is
used to reconstruct ρ, ϵ, T = p/ρ, the adiabatic sound speed cs² and the
fixed-baryon sound speed cn² = Γ-1 (the Caballero–Yunes heat-conduction inputs).

Keyword `ρc` is the central rest-mass density [km^-2]; `Γ_struct` defaults to Γ
(isentropic — see `IdealGasStar` docs for the stratified Γ_struct ≠ Γ closure).
Extra keywords (`ptol_rel`, `rmax`) are forwarded to `solve_tov`.

Returns an `IdealGasStar` with the profile vectors and M, R. By construction
T(r) > 0 with dT/dr < 0, and cs², cn² ∈ (0,1) for Γ ∈ (1,2).
"""
function solve_tov_idealgas(; Γ::Real=5/3, Γ_struct::Real=Γ, K::Real=100.0,
                            ρc::Real, h::Float64=2e-4, N_tab::Int=1200,
                            ρ_lo::Real=1e-9, kwargs...)
    Γ = float(Γ); Γ_struct = float(Γ_struct); K = float(K)
    # central total energy density of the structure polytrope
    pc = K * ρc^Γ_struct
    εc = ρc + pc / (Γ_struct - 1)
    # tabulate the structure barotrope p(ε) over a range bracketing the star
    eos = isentropic_idealgas(Γ=Γ_struct, K=K, ρ_lo=ρ_lo, ρ_hi=2*ρc, N=N_tab)
    star = solve_tov(eos, εc; h=h, kwargs...)

    n = length(star.r)
    ρ   = Vector{Float64}(undef, n)
    ϵ   = Vector{Float64}(undef, n)
    T   = Vector{Float64}(undef, n)
    cs2 = Vector{Float64}(undef, n)
    cn2v= Vector{Float64}(undef, n)
    micro = IdealGas(Γ)
    for i in 1:n
        p = star.p[i]
        if p > 0
            ρi  = (p / K)^(1/Γ_struct)            # invert structure polytrope
            ϵi  = p / ((Γ - 1) * ρi)              # ideal-gas internal energy
            ρ[i] = ρi; ϵ[i] = ϵi
            T[i] = p / ρi                          # = (Γ-1) ϵ
            cs2[i]  = sound_speed2(micro, ρi, ϵi)
            cn2v[i] = cn2(micro, ρi, ϵi)
        else
            # surface: take the limit (ϵ, T, cs², cn² are continuous to the edge)
            ρ[i] = 0.0; ϵ[i] = 0.0; T[i] = 0.0
            cs2[i]  = 0.0
            cn2v[i] = Γ - 1                        # cn² = Γ-1 is ρ,ϵ-independent
        end
    end
    return IdealGasStar(star.r, ρ, star.p, star.ε, ϵ, T, cs2, cn2v,
                        star.m, star.ν, star.M, star.R, Γ, Γ_struct, K)
end

# Msun_to_km from Units (duplicated constant to avoid a cross-module dep cycle)
const _MSUN_TO_KM = 1.988416e30 * (6.6743015e-11 / 299_792_458.0^2) * 1e-3

"""
    mass_solar(star) -> M [M⊙]

Convert the geometric mass `star.M` to solar masses, **assuming the KM-geometric
convention** (lengths in km, ε in km⁻²), i.e. `M[M⊙] = M[km]/1.4766`.

!!! warning "The unit convention is set by the CALLER, not by the solver"
    `solve_tov` is scale-agnostic: it integrates in whatever geometric units the
    EOS + εc are expressed in. Two conventions are in use in this package and
    `mass_solar` is only correct for the first:

    * **km-geometric** — εc in km⁻² (e.g. `piecewise_polytrope` via the CGS map
      `ρ[g/cm³]·7.4237e-19`). Lengths are km, `star.M` is km ⇒ `mass_solar` is
      CORRECT. SLy at ρc=1e15 g/cm³ ⇒ M=1.417 M⊙, R=11.68 km.
    * **M⊙-geometric** — the G=c=M⊙=1 convention standard in the NS-oscillation
      literature. Here `star.M` is ALREADY in M⊙ and `star.R` is in M⊙ units
      (multiply by 1.4766 for km); applying `mass_solar` would wrongly divide
      again by 1.4766.

    Because the Γ=2 (n=1) polytrope is SCALE-INVARIANT (M ∝ √κ), the same raw
    `ShumPolytrope(100.0)` solution legitimately reproduces two DIFFERENT published
    benchmark stars depending on the units assigned to κ — see the note on
    `ShumPolytrope` in eos/EquationOfState.jl. Always state which convention a
    reported mass/radius/frequency uses.
"""
mass_solar(star::TOVStar)     = star.M / _MSUN_TO_KM
mass_solar(star::IdealGasStar)= star.M / _MSUN_TO_KM
mass_solar(M_km::Real)        = M_km / _MSUN_TO_KM

end # module TOV
