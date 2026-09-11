#=
    Background3D — frozen TOV background on a 3D Cartesian grid (the static
    spacetime for the STAGE-3 Phase-2 non-radial Cowling evolution).

    A spherical TOV star (areal/Schwarzschild coordinates, metric
        ds² = −e^{ν}dt² + (1−2m/r)^{-1}dr² + r²dΩ² )
    is embedded in a cubic Cartesian box x,y,z ∈ [−L,L]³. At each cell the 1-D
    TOV profile is interpolated by radius r=√(x²+y²+z²); outside the surface the
    Schwarzschild exterior (m=M, e^{ν}=1−2M/r) and an atmosphere floor are used.

    Lapse  α = e^{ν/2}  (−g_tt = e^{ν} = α²).
    Spatial metric (areal→Cartesian):  γ_ij = δ_ij + (e^{2λ}−1) n_i n_j ,
        n_i = x_i/r,  e^{2λ} = (1−2m/r)^{-1}  ⇒  det γ = e^{2λ},  √γ = e^{λ}.
    Eigenvalues of γ_ij are (e^{2λ}, 1, 1): radial stretch, transverse unchanged.

    This module owns ONLY the static background + grid (verifiable against the
    1-D TOV); the evolution scheme lives alongside it.
=#
module Background3D

using ..EquationOfState
using ..TOV

export CartesianGrid, Star3D, build_star3d, radius_field

struct CartesianGrid
    N::Int                 # cells per dimension
    L::Float64             # half-box size (same geometric units as the TOV star)
    dx::Float64
    x::Vector{Float64}     # cell-centre coordinates (shared across dims)
end

function CartesianGrid(N::Int, L::Float64)
    dx = 2L / N
    x = collect(range(-L + dx/2, L - dx/2; length=N))
    CartesianGrid(N, L, dx, x)
end

# All background fields live on the N³ grid (areal-Cartesian).
struct Star3D
    grid::CartesianGrid
    R::Float64; M::Float64
    r::Array{Float64,3}                 # areal radius
    ρ0::Array{Float64,3}                # rest-mass density (atmosphere floor outside)
    p0::Array{Float64,3}                # pressure
    ε0::Array{Float64,3}                # total energy density
    cs2::Array{Float64,3}              # sound speed squared dp/dε
    α::Array{Float64,3}                 # lapse e^{ν/2}
    e2λ::Array{Float64,3}              # g_rr = (1−2m/r)^{-1}
    sqrtγ::Array{Float64,3}            # √det γ = e^{λ}
    nx::Array{Float64,3}; ny::Array{Float64,3}; nz::Array{Float64,3}  # radial unit vector
    interior::BitArray{3}              # r ≤ R (inside the star)
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
    build_star3d(eos, εc; N=64, Lfac=1.3, h_tov=2e-4, ρ_atm_rel=1e-9) -> Star3D

Solve the TOV star and sample its frozen background on an N³ Cartesian grid of
half-width `Lfac·R`. `εc` central energy density (EOS geometric units).
"""
function build_star3d(eos::BarotropicEOS, εc::Float64;
                      N::Int=64, Lfac::Float64=1.3, h_tov::Float64=2e-4,
                      ρ_atm_rel::Float64=1e-9)
    star = solve_tov(eos, εc; h=h_tov)
    R, M = star.R, star.M
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε
    grid = CartesianGrid(N, Lfac * R)
    xs = grid.x
    ρ_atm = εc * ρ_atm_rel

    r   = Array{Float64}(undef, N, N, N)
    ρ0  = similar(r); p0 = similar(r); ε0 = similar(r); cs2 = similar(r)
    α   = similar(r); e2λ = similar(r); sqγ = similar(r)
    nx  = similar(r); ny = similar(r); nz = similar(r)
    interior = falses(N, N, N)

    @inbounds for k in 1:N, j in 1:N, i in 1:N
        X, Y, Z = xs[i], xs[j], xs[k]
        rr = sqrt(X^2 + Y^2 + Z^2); rr = max(rr, 1e-12)
        r[i,j,k] = rr
        nx[i,j,k] = X/rr; ny[i,j,k] = Y/rr; nz[i,j,k] = Z/rr
        if rr ≤ R
            interior[i,j,k] = true
            m = _lin(rt, mt, rr); ν = _lin(rt, νt, rr)
            ε = max(_lin(rt, et, rr), ρ_atm)
            p = pressure(eos, ε)
            ε0[i,j,k]  = ε
            p0[i,j,k]  = p
            cs2[i,j,k] = max(sound_speed2(eos, ε), 0.0)
            α[i,j,k]   = exp(ν/2)
            e2λ[i,j,k] = 1.0 / (1.0 - 2m/rr)
        else
            # Schwarzschild exterior + atmosphere
            ε = ρ_atm
            interior[i,j,k] = false
            ρ0[i,j,k]  = ρ_atm
            ε0[i,j,k]  = ε
            p0[i,j,k]  = pressure(eos, ε)
            cs2[i,j,k] = max(sound_speed2(eos, ε), 0.0)
            α[i,j,k]   = sqrt(max(1 - 2M/rr, 1e-12))
            e2λ[i,j,k] = 1.0 / max(1 - 2M/rr, 1e-12)
        end
        sqγ[i,j,k] = sqrt(e2λ[i,j,k])           # √γ = e^{λ}
    end
    # rest-mass density from the barotrope (ε = ρ(1+ε_int); for Shum Γ=2: ρ = ε − p)
    @inbounds for I in eachindex(ρ0)
        ρ0[I] = max(ε0[I] - p0[I], ρ_atm)       # exact for the Γ=2 (e=ρ+p) family
    end

    Star3D(grid, R, M, r, ρ0, p0, ε0, cs2, α, e2λ, sqγ, nx, ny, nz, interior)
end

"""radius of every cell (diagnostic helper)."""
radius_field(s::Star3D) = s.r

end # module Background3D
