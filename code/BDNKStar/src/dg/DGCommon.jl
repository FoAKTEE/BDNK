#=
    DGCommon — nodal Discontinuous-Galerkin (RKDG) building blocks shared by the
    1D special-relativistic shock tube (Stage 1), the 1D radial TOV star
    (Stage 2) and the 2D Cartesian non-conforming surface (Stage 3).

    This is the FALLBACK scheme for the verified FV-Cartesian verdict
    'fv-insufficient-need-dg' (see src/fvcartesian/FVCartesian.jl docstring): a
    2nd-order FV + MinMod staircases the non-conforming stellar surface and
    drives an ℓ=2-perturbed secular instability that WORSENS with resolution.
    A nodal DG scheme carries a sub-cell polynomial in each element, so a steep
    surface is resolved INSIDE a cell; a troubled-cell limiter + a
    positivity-preserving (Zhang–Shu) limiter keep ρ,p>0 monotone without the
    1/Δx noise amplification that kills the FV surface.

    Nodal DG on Legendre–Gauss–Lobatto (LGL) nodes:
      * solution in element e is a degree-p polynomial sampled at the (p+1) LGL
        nodes ξ_a ∈ [-1,1];
      * the LGL nodes INCLUDE the endpoints ±1, so the interface (face) values
        needed by the numerical flux are nodal values — no extra interpolation;
      * the diagonal LGL mass matrix (mass-lumping) makes the scheme cheap and is
        exactly the spectral-element collocation form.

    Semidiscrete DG for ∂_t u + ∂_x f(u) = s on element [x_L,x_R], Jacobian
    J=Δx/2:
        M (du/dt) = (1/J) Dᵀ M f  −  [ φ_a f̂ ]_face  + M s
    With mass lumping (M=diag(w)) and the strong/weak form we use the standard
        du_a/dt = (1/J) [ Σ_b D̃_{ab} f_b ]  − boundary correction + s_a
    where the boundary correction injects the numerical (Rusanov) flux at the two
    faces through the inverse-mass-weighted Lagrange endpoint values.
=#
module DGCommon

export lgl_nodes_weights, vandermonde_legendre, diff_matrix,
       legendre_P, rusanov_flux, tvb_minmod, cell_average,
       LGLBasis, build_lgl_basis

# ---------------------------------------------------------------------------
# Legendre polynomial P_n(x) and derivative (Bonnet recursion).
# ---------------------------------------------------------------------------
@inline function legendre_P(n::Int, x::Float64)
    n == 0 && return 1.0
    n == 1 && return x
    p0 = 1.0; p1 = x
    for k in 1:n-1
        p2 = ((2k+1)*x*p1 - k*p0)/(k+1)
        p0 = p1; p1 = p2
    end
    return p1
end

# P_n and dP_n/dx together.
@inline function legendre_P_dP(n::Int, x::Float64)
    n == 0 && return (1.0, 0.0)
    n == 1 && return (x, 1.0)
    p0=1.0; p1=x; dp0=0.0; dp1=1.0
    for k in 1:n-1
        p2 = ((2k+1)*x*p1 - k*p0)/(k+1)
        dp2 = dp0 + (2k+1)*p1
        p0=p1; p1=p2; dp0=dp1; dp1=dp2
    end
    return (p1, dp1)
end

# ---------------------------------------------------------------------------
# Legendre–Gauss–Lobatto nodes & weights on [-1,1] for polynomial degree p
# (N=p+1 nodes). Interior nodes are roots of P_p'(x); endpoints are ±1.
# Newton iteration on (1-x²)P_p'(x) using the relation with P_p.
# ---------------------------------------------------------------------------
function lgl_nodes_weights(p::Int)
    N = p + 1
    x = zeros(N); w = zeros(N)
    x[1] = -1.0; x[N] = 1.0
    # Lobatto weight at endpoints: 2/(p(p+1))
    we = 2.0/(p*(p+1))
    w[1] = we; w[N] = we
    # interior nodes: initial guess = Chebyshev–Gauss–Lobatto, Newton on P_p'(x)=0
    for i in 2:N-1
        xi = -cos(π*(i-1)/p)
        for _ in 1:100
            # we need root of P_p'(x). Use d/dx of P_p via legendre_P_dP, and the
            # second derivative from the ODE (1-x²)P'' - 2xP' + p(p+1)P = 0 ⇒
            # P'' = (2x P' - p(p+1)P)/(1-x²).
            Pp, dPp = legendre_P_dP(p, xi)
            ddPp = (2*xi*dPp - p*(p+1)*Pp)/(1-xi^2)
            δ = dPp/ddPp
            xi -= δ
            abs(δ) < 1e-15 && break
        end
        x[i] = xi
        Pp = legendre_P(p, xi)
        w[i] = 2.0/(p*(p+1)*Pp^2)
    end
    # sort (Newton may land slightly out of order for high p)
    perm = sortperm(x)
    return x[perm], w[perm]
end

# ---------------------------------------------------------------------------
# Vandermonde V_{a,n}=P_n(x_a) (normalized Legendre) and the nodal
# differentiation matrix D = Vr * inv(V) with Vr_{a,n}=P_n'(x_a).
# We use ORTHONORMAL Legendre so V is well conditioned.
# ---------------------------------------------------------------------------
@inline _norm_leg(n) = sqrt((2n+1)/2)   # ∫_{-1}^1 P_n² = 2/(2n+1)

function vandermonde_legendre(x::Vector{Float64}, p::Int)
    N = length(x)
    V = zeros(N, p+1); Vr = zeros(N, p+1)
    for n in 0:p
        c = _norm_leg(n)
        for a in 1:N
            P, dP = legendre_P_dP(n, x[a])
            V[a, n+1] = c*P
            Vr[a, n+1] = c*dP
        end
    end
    return V, Vr
end

# Nodal differentiation matrix on the given nodes: D u gives du/dξ at the nodes.
function diff_matrix(x::Vector{Float64}, p::Int)
    V, Vr = vandermonde_legendre(x, p)
    return Vr / V    # D = Vr V^{-1}
end

# ---------------------------------------------------------------------------
# Packaged LGL basis for an element of polynomial degree p.
# ---------------------------------------------------------------------------
struct LGLBasis
    p::Int
    N::Int                       # = p+1
    ξ::Vector{Float64}           # LGL nodes on [-1,1] (ξ[1]=-1, ξ[N]=+1)
    w::Vector{Float64}           # LGL quadrature weights (Σw=2)
    D::Matrix{Float64}           # nodal differentiation matrix d/dξ
    V::Matrix{Float64}           # orthonormal-Legendre Vandermonde
    invV::Matrix{Float64}        # V^{-1} (modal ↔ nodal)
end

function build_lgl_basis(p::Int)
    ξ, w = lgl_nodes_weights(p)
    V, _ = vandermonde_legendre(ξ, p)
    D = diff_matrix(ξ, p)
    LGLBasis(p, p+1, ξ, w, D, V, inv(V))
end

# cell average of a nodal field u_a using LGL weights: ū = (1/2)Σ w_a u_a
@inline function cell_average(u, w)
    s = 0.0; @inbounds for a in eachindex(u); s += w[a]*u[a]; end
    return 0.5*s
end

# ---------------------------------------------------------------------------
# Local Lax–Friedrichs (Rusanov) numerical flux for a conserved system.
# UL,FL,UR,FR are NTuples; amax is the max signal speed across the interface.
# ---------------------------------------------------------------------------
@inline function rusanov_flux(UL::NTuple{K,Float64}, FL::NTuple{K,Float64},
                              UR::NTuple{K,Float64}, FR::NTuple{K,Float64},
                              amax::Float64) where {K}
    ntuple(i -> 0.5*(FL[i]+FR[i]) - 0.5*amax*(UR[i]-UL[i]), K)
end

# ---------------------------------------------------------------------------
# TVB-modified minmod (Cockburn–Shu). m̄(a1,a2,a3) returns a1 if |a1|≤M·h², else
# the standard minmod. Used by the troubled-cell limiter.
# ---------------------------------------------------------------------------
@inline function tvb_minmod(a1::Float64, a2::Float64, a3::Float64, Mh2::Float64)
    abs(a1) ≤ Mh2 && return a1
    s = sign(a1)
    (s == sign(a2) && s == sign(a3)) || return 0.0
    return s*min(abs(a1), abs(a2), abs(a3))
end

end # module DGCommon
