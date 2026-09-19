#=
    DGSubcell — the DG / finite-difference subcell machinery of Dumbser et al. (2014) as
    implemented in SpECTRE (Deppe et al., PRD 105, 123031 (2022), arXiv:2109.12033; source
    `src/Evolution/DgSubcell/`): a troubled-cell indicator decides, element by element, whether
    the solution is evolved with the high-order DG scheme or with a robust finite-difference
    scheme on a subcell grid carrying the SAME degrees of freedom, and the two representations
    are exchanged by a conservative projection / constrained-least-squares reconstruction pair.

    This module holds the grid-agnostic pieces; the 1D star (DGStarFD) and the 3D Cartesian
    star apply them per element (in 3D dimension by dimension, as a tensor product).

    SUBCELL GRID. A DG element with N Legendre–Gauss–Lobatto points per direction is paired
    with M = 2N−1 equidistant finite-volume cells of width h = 2/M in the reference coordinate
    (SpECTRE `evolution::dg::subcell::fd::mesh`). M > N, so the FD grid is finer than the DG
    one wherever it matters and the DOF count grows only mildly.

    PROJECTION P (DG → FD), M×N. The subcell value is the exact cell average of the DG
    polynomial,  P_ij = (1/h) ∫_{cell i} ℓ_j(ξ) dξ,  computed with Gauss–Legendre quadrature
    (exact: ℓ_j has degree N−1). Conservative by construction: h Σ_i (Pu)_i = ∫u dξ.

    RECONSTRUCTION R (FD → DG), N×M. There are more subcells than DG points, so R solves the
    constrained least-squares problem
        minimise ‖P u − v‖²   subject to   Σ_j w_j u_j = h Σ_i v_i,
    i.e. the reconstructed DG polynomial matches the subcell data as well as possible while
    carrying exactly the same integral. With A = 2PᵀP,
        R = A⁻¹(2Pᵀ) − (A⁻¹w / wᵀA⁻¹w) (wᵀA⁻¹(2Pᵀ) − h 1ᵀ),
    which is SpECTRE's formula (`subcell::fd::reconstruction_matrix`). Because the LGL rule is
    exact for the degree-(N−1) DG polynomial the constraint is automatically satisfied by the
    unconstrained minimiser when v = Pu, so R P = I exactly — projecting and reconstructing a
    DG solution returns it unchanged (verified to 1e-15 in test_dgsubcell.jl).

    TROUBLED-CELL INDICATORS.
      • Persson (`subcell::persson_tci`, Persson & Peraire 2006): with the element's modal
        (orthonormal-Legendre) coefficients ĉ, let Û be the nodal field of the top `nmodes`
        modes alone. The element is troubled when ‖Û‖₂/‖U‖₂ > (N − nmodes)^(−α), the default
        being α = 4 and nmodes = 1. This measures how much power sits in the modes that carry
        the Gibbs oscillations of an under-resolved feature.
      • Relaxed discrete maximum principle (`subcell::rdmp_tci`): the candidate solution must
        satisfy  min_N(u^n) − δ ≤ u* ≤ max_N(u^n) + δ  with the neighbourhood N being the
        element and its neighbours and δ = max(δ₀, ε (max_N − min_N)), δ₀ = 1e-4, ε = 1e-3.
      • Physical admissibility (SpECTRE's `TciOptions`) is checked by the caller, which knows
        the equation of state: density above the atmosphere cutoff, positive pressure proxy.

    FINITE-DIFFERENCE RECONSTRUCTION. `mc_slope` is the monotonised-central slope used to
    reconstruct primitives to the subcell faces; it is SpECTRE's low-order fallback
    (`MonotonisedCentral`), second order and positivity-friendly.
=#
module DGSubcell

using ..DGCommon: LGLBasis
using LinearAlgebra: SymTridiagonal, eigen, dot

export SubcellOps, subcell_ops, persson_ratio, persson_threshold, persson_troubled,
       rdmp_troubled, minmod2, mc_slope

# ----------------------------------------------------------------------------------
# Gauss–Legendre nodes/weights on [−1,1] (Golub–Welsch)
# ----------------------------------------------------------------------------------
function gauss_legendre(n::Int)
    n == 1 && return ([0.0], [2.0])
    β = [k/sqrt(4k^2 - 1.0) for k in 1:n-1]
    F = eigen(SymTridiagonal(zeros(n), β))
    x = F.values; w = [2*F.vectors[1, i]^2 for i in 1:n]
    perm = sortperm(x)
    return x[perm], w[perm]
end

@inline function _lagrange(ξ::Vector{Float64}, j::Int, x::Float64)
    ℓ = 1.0
    @inbounds for m in eachindex(ξ)
        m == j && continue
        ℓ *= (x - ξ[m])/(ξ[j] - ξ[m])
    end
    ℓ
end

"""
    SubcellOps

Projection/reconstruction pair for one DG basis: `N` LGL points ↔ `M = 2N−1` equidistant
subcells of reference width `h`, centres `ξc`. `P` (M×N) projects a DG nodal field to subcell
averages, `R` (N×M) reconstructs a DG nodal field from subcell averages; `R*P = I`.
"""
struct SubcellOps
    N::Int
    M::Int
    h::Float64
    ξc::Vector{Float64}
    P::Matrix{Float64}
    R::Matrix{Float64}
end

"""
    subcell_ops(basis) -> SubcellOps

Build the conservative projection and the constrained-least-squares reconstruction for an
LGL basis (SpECTRE `subcell::fd::projection_matrix` / `reconstruction_matrix`).
"""
function subcell_ops(b::LGLBasis)
    N = b.N; M = 2N - 1; h = 2/M
    ξc = [-1 + (i - 0.5)*h for i in 1:M]
    xg, wg = gauss_legendre(max(N, 2))
    P = zeros(M, N)
    @inbounds for i in 1:M, j in 1:N
        s = 0.0
        for q in eachindex(xg)
            s += wg[q]*_lagrange(b.ξ, j, ξc[i] + 0.5*h*xg[q])
        end
        P[i, j] = 0.5*s                      # (1/h)∫ = (1/h)(h/2)Σ w f
    end
    A = 2*(transpose(P)*P)
    Ai = inv(A)
    w = collect(b.w)
    Aiw = Ai*w
    denom = dot(w, Aiw)
    twoPt = 2*transpose(P)
    R = Ai*twoPt .- (Aiw/denom)*(transpose(w)*Ai*twoPt .- h*ones(1, M))
    return SubcellOps(N, M, h, ξc, P, R)
end

# ----------------------------------------------------------------------------------
# troubled-cell indicators
# ----------------------------------------------------------------------------------
"""
    persson_ratio(u, basis; nmodes=1) -> ‖Û‖₂/‖U‖₂

Fraction of the element's nodal L2 norm carried by its top `nmodes` spectral modes.
"""
function persson_ratio(u::AbstractVector{Float64}, b::LGLBasis; nmodes::Int=1)
    N = b.N
    û = b.invV*collect(u)
    top = zeros(N)
    @inbounds for i in max(1, N - nmodes + 1):N; top[i] = û[i]; end
    ûf = b.V*top
    den = sqrt(sum(abs2, u))
    den ≤ 1e-300 && return 0.0
    return sqrt(sum(abs2, ûf))/den
end

"""Persson threshold (N − nmodes)^(−α); the element is troubled above it."""
@inline persson_threshold(b::LGLBasis; α::Float64=4.0, nmodes::Int=1) =
    Float64(max(b.N - nmodes, 1))^(-α)

persson_troubled(u::AbstractVector{Float64}, b::LGLBasis; α::Float64=4.0, nmodes::Int=1) =
    persson_ratio(u, b; nmodes=nmodes) > persson_threshold(b; α=α, nmodes=nmodes)

"""
    rdmp_troubled(cand_min, cand_max, past_min, past_max; δ0=1e-4, ε=1e-3)

Relaxed discrete maximum principle: `true` when the candidate leaves the neighbourhood range
of the previous step by more than δ = max(δ₀, ε(past_max − past_min)).
"""
@inline function rdmp_troubled(cand_min::Float64, cand_max::Float64, past_min::Float64,
                               past_max::Float64; δ0::Float64=1e-4, ε::Float64=1e-3)
    δ = max(δ0, ε*(past_max - past_min))
    return (cand_min < past_min - δ) || (cand_max > past_max + δ)
end

# ----------------------------------------------------------------------------------
# finite-difference reconstruction
# ----------------------------------------------------------------------------------
@inline minmod2(a::Float64, b::Float64) = (a*b ≤ 0.0) ? 0.0 : (abs(a) < abs(b) ? a : b)
"""Monotonised-central slope (SpECTRE's `MonotonisedCentral`) from three cell averages."""
@inline mc_slope(um::Float64, u0::Float64, up::Float64) =
    minmod2(2*(u0 - um), minmod2(0.5*(up - um), 2*(up - u0)))

end # module DGSubcell
