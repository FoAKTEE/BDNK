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

export TOVStar, solve_tov, mass_solar

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

# Msun_to_km from Units (duplicated constant to avoid a cross-module dep cycle)
const _MSUN_TO_KM = 1.988416e30 * (6.6743015e-11 / 299_792_458.0^2) * 1e-3
mass_solar(star::TOVStar) = star.M / _MSUN_TO_KM
mass_solar(M_km::Real)    = M_km / _MSUN_TO_KM

end # module TOV
