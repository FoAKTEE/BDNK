#=
    RadialModes — radial (l=0) stellar pulsations in the relativistic Cowling
    approximation (frozen metric), the linear backbone of STAGE 1A
    (Caballero–Yunes). Frequency-domain matrix eigensolver.

    Pulsation ODE for the Lagrangian displacement ξ(r) (verbatim operator from
    NeutronStarOscillations.jl FrequencyDomain.jl PerfectFluidCowling):

        A₂(r) ξ'' + A₁(r) ξ' + A₀(r) ξ = ω² ξ,

    with A₀,A₁,A₂ explicit functions of the TOV background (m,p,ε,ν) and the
    sound speed c_s=√(dp/dε) and c_s'(r). Second-order central differences on a
    uniform grid; ξ(0)=0 at the centre; Lagrangian Δp=0 at the surface. The
    eigenvalues ω² (real, positive ⇒ stable) give the mode frequencies.
=#
module RadialModes

using ..EquationOfState
using ..TOV
using LinearAlgebra

export radial_cowling_spectrum, hz_per_invkm

# kHz <-> km^-1 (geometric): f[kHz] = (√ω²/2π)[km^-1] / kHz_to_km
const _sec_to_km = 299792458.0 * 1e-3
const _kHz_to_km = 1e3 / _sec_to_km
"""frequency in kHz from ω² in km⁻² (geometric, εc supplied in km⁻²)."""
hz_per_invkm(ω2::Real) = sqrt(ω2) / (2π) / _kHz_to_km

# --- pulsation-operator coefficients (NSO PerfectFluidCowling A0/A1/A2) -------
@inline _A0(m,p,ε,ν,cs,csp,r) = (2*exp(ν)*(m^2*(1 + 5*cs*(cs - 2*r*csp)) +
    r*m*(-1 - 8π*r^2*p + cs*((-5 + 8π*r^2*ε)*cs + r*(9 + 8π*r^2*ε)*csp)) +
    r^2*(cs^2 + 2*r*(π*r*(p + ε + 8π*r^2*p*ε - 8π*r^2*ε^2*cs^2) -
    (1 + 2π*r^2*ε)*cs*csp)))) / (r^3*(r - 2*m))
@inline _A1(m,p,ε,ν,cs,csp,r) = (exp(ν)*(4π*r^3*p - 2*r*cs*((1 + 2π*r^2*ε)*cs + r*csp) +
    m*(1 + 5*cs^2 + 4*r*cs*csp))) / r^2
@inline _A2(m,p,ε,ν,cs,csp,r) = -((exp(ν)*(r - 2*m)*cs^2)/r)

@inline _lin(xs, ys, x) = begin            # linear interpolation on ascending xs
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j]) / (xs[j+1] - xs[j])
    ys[j] + t*(ys[j+1] - ys[j])
end

"""
    radial_cowling_spectrum(eos, εc; N=1000, h_tov=5e-5, nmodes=5)
        -> (freqs_kHz, ω²_invkm2, R)

Eigenfrequencies of the lowest `nmodes` radial Cowling modes. `εc` in km⁻².
"""
function radial_cowling_spectrum(eos::BarotropicEOS, εc::Float64;
                                 N::Int=1000, h_tov::Float64=5e-5, nmodes::Int=5)
    star = solve_tov(eos, εc; h=h_tov)
    R = star.R
    rt, mt, pt, et, νt = star.r, star.m, star.p, star.ε, star.ν
    hh = R / N
    rg = collect(range(0.0, R; length=N+1))           # r[1]=0 … r[N+1]=R
    εf = εc * 1e-9                                      # atmosphere floor (cs>0 at surface)
    # background + sound speed on the uniform grid (barotrope: p=p(ε) exactly)
    bg(r) = begin
        m = _lin(rt, mt, r); ν = _lin(rt, νt, r)
        ε = max(_lin(rt, et, r), εf)
        p = pressure(eos, ε)
        cs = sqrt(max(sound_speed2(eos, ε), 0.0))
        pprime = -(ε+p)*(m + 4π*r^3*p)/(r*(r - 2m))     # TOV dp/dr
        csp = d2pde2(eos, ε) * pprime / (2*cs^3)        # dcs/dr
        return m, p, ε, ν, cs, csp
    end
    n = N                                               # matrix on r[2..N+1]
    A = zeros(n, n)
    # first row at r[2] (uses ξ(0)=0)
    m,p,ε,ν,cs,csp = bg(rg[2]); a0=_A0(m,p,ε,ν,cs,csp,rg[2]); a1=_A1(m,p,ε,ν,cs,csp,rg[2]); a2=_A2(m,p,ε,ν,cs,csp,rg[2])
    A[1,1] = a0 - 2a2/hh^2
    A[1,2] = (2a2 + a1*hh)/(2hh^2)
    # interior rows
    for i in 2:n-1
        r = rg[i+1]; m,p,ε,ν,cs,csp = bg(r)
        a0=_A0(m,p,ε,ν,cs,csp,r); a1=_A1(m,p,ε,ν,cs,csp,r); a2=_A2(m,p,ε,ν,cs,csp,r)
        A[i,i-1] = (2a2 - a1*hh)/(2hh^2)
        A[i,i]   = a0 - 2a2/hh^2
        A[i,i+1] = (2a2 + a1*hh)/(2hh^2)
    end
    # surface row at r[end]=R, enforcing Lagrangian Δp=0 (NSO An_minus1/Ann)
    r = rg[end]; m,p,ε,ν,cs,csp = bg(r)
    a0=_A0(m,p,ε,ν,cs,csp,r); a1=_A1(m,p,ε,ν,cs,csp,r); a2=_A2(m,p,ε,ν,cs,csp,r)
    A[n,n-1] = (2a2)/hh^2
    A[n,n] = ((2a2*(5hh+2r) + hh^2*(5a1 - 2a0*r))*m -
              r*(2a2*(2hh+r) + hh^2*(2a1 - a0*r) + 4hh*(2a2 + a1*hh)*π*r^2*ε)) /
             (hh^2*r*(r - 2m))
    vals = eigvals(A)
    ω2 = sort([real(v) for v in vals if abs(imag(v)) < 1e-6*abs(real(v)) + 1e-30 && real(v) > 0])
    k = min(nmodes, length(ω2))
    return [hz_per_invkm(w) for w in ω2[1:k]], ω2[1:k], R
end

end # module RadialModes
