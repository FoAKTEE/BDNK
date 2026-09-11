#=
    SphBackground — axisymmetric (r,θ) BOUNDARY-CONFORMING TOV background for the
    spherical BDNK rewrite (STAGE 3, surface treatment).

    The instability of the Cartesian full-frame BDNK engine is the BDNK velocity
    recovery (den = τ_Q w₀ − η → 0) evaluated on the ragged, staircased stellar
    surface of a sphere embedded in a Cartesian grid. The cure (as in the working
    1-D radial BDNK codes and the Pretorius two-sphere paper, arXiv:2508.20998) is
    a boundary-conforming grid: the stellar surface sits on the coordinate line
    r = R, so the Δp=0 surface condition is applied cleanly and the near-surface
    radial structure is resolved rather than staircased.

    Grid: r ∈ (0, R] (cell-centred, last point at the surface R), θ ∈ (0, π)
    (cell-centred, avoiding the pole coordinate singularity). Axisymmetric ⇒ no φ;
    captures the ℓ=2, m=0 f-mode (degenerate in m for a non-rotating star).

    Metric (geometrised, ν,Λ full exponents):
        ds² = −e^{ν(r)}dt² + e^{Λ(r)}dr² + r²(dθ² + sin²θ dφ²),
        e^{Λ}=(1−2m/r)^{-1},  α=e^{ν/2},  √γ = e^{Λ/2} r² sinθ.
    The background is spherically symmetric ⇒ all matter/metric profiles depend on
    r only (stored as length-Nr vectors); the θ grid carries only sinθ/cosθ.
=#
module SphBackground

using ..EquationOfState
using ..TOV

export SphGrid, SphStar, build_sphstar

struct SphGrid
    Nr::Int; Nθ::Int
    r::Vector{Float64};  θ::Vector{Float64}
    dr::Float64;         dθ::Float64
    sinθ::Vector{Float64}; cosθ::Vector{Float64}
end

# radial background profiles (length Nr) + grid; θ-dependence is only sinθ/cosθ
struct SphStar
    grid::SphGrid
    R::Float64; M::Float64
    ρ0::Vector{Float64};  p0::Vector{Float64};  ε0::Vector{Float64};  cs2::Vector{Float64}
    α::Vector{Float64}        # lapse e^{ν/2}
    eΛ::Vector{Float64}       # g_rr = (1−2m/r)^{-1}
    m::Vector{Float64}        # enclosed gravitational mass
    νp::Vector{Float64}       # ν'(r) = 2Φ' (redshift gradient)
    sqγr::Vector{Float64}     # √γ without sinθ:  e^{Λ/2} r²
end

@inline function _lin(xs, ys, x)
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j]) / (xs[j+1] - xs[j])
    ys[j] + t * (ys[j+1] - ys[j])
end

"""
    build_sphstar(eos, εc; Nr=200, Nθ=32, h_tov=2e-4) -> SphStar

Solve the TOV star and sample its frozen background on a cell-centred (r,θ) grid
with r ∈ (0,R] (surface ON the last radial point) and θ ∈ (0,π).
"""
function build_sphstar(eos::BarotropicEOS, εc::Float64;
                       Nr::Int=200, Nθ::Int=32, h_tov::Float64=2e-4)
    star = solve_tov(eos, εc; h=h_tov)
    R, M = star.R, star.M
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε

    dr = R / Nr
    rg = collect(range(dr/2, R - dr/2; length=Nr))     # cell centres, last just inside R
    dθ = π / Nθ
    θg = collect(range(dθ/2, π - dθ/2; length=Nθ))      # cell centres, avoid poles
    grid = SphGrid(Nr, Nθ, rg, θg, dr, dθ, sin.(θg), cos.(θg))

    ρ0=similar(rg); p0=similar(rg); ε0=similar(rg); cs2=similar(rg)
    α=similar(rg); eΛ=similar(rg); mm=similar(rg); νp=similar(rg); sqγr=similar(rg)
    @inbounds for i in 1:Nr
        r = rg[i]
        m = _lin(rt, mt, r); ν = _lin(rt, νt, r)
        ε = max(_lin(rt, et, r), εc*1e-12)
        p = pressure(eos, ε)
        mm[i]=m; ε0[i]=ε; p0[i]=p
        ρ0[i] = max(ε - p, εc*1e-12)               # ρ = ε−p for the Γ=2 (e=ρ+p) family
        cs2[i] = max(sound_speed2(eos, ε), 0.0)
        α[i] = exp(ν/2)
        eΛ[i] = 1.0 / (1.0 - 2m/r)
        νp[i] = 2*(m + 4π*r^3*p) / (r*(r - 2m))    # ν' = 2Φ'
        sqγr[i] = sqrt(eΛ[i]) * r^2                  # √γ / sinθ = e^{Λ/2} r²
    end
    SphStar(grid, R, M, ρ0, p0, ε0, cs2, α, eΛ, mm, νp, sqγr)
end

end # module SphBackground
