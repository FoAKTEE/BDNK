#=
    DGSRHydro — STAGE 1 (the MVP, the user's explicit ask "if FV does not work,
    try pseudospectral / discontinuous-Galerkin to capture shocks").

    1D special-relativistic hydrodynamics on flat space with a Γ-law ideal gas
    p=(Γ-1)ρϵ, solved with a NODAL RKDG scheme:
      * Legendre–Gauss–Lobatto nodal basis, degree p (default 3);
      * local Lax–Friedrichs (Rusanov) interface flux;
      * SSP-RK3 time stepping;
      * a TVB-minmod TROUBLED-CELL limiter (component-wise on the characteristic-
        free conserved variables, via cell means + linear slope) flags & limits
        shocked cells;
      * a Zhang–Shu POSITIVITY-PRESERVING limiter squeezes the nodal polynomial
        toward its cell mean so D>0 and p>0 at every node and quadrature point.

    Conserved (flat space, Valencia / SRHD):
        D  = ρW
        S  = ρhW² v
        τ  = ρhW² − p − D
    with W=1/√(1−v²), h=1+ϵ+p/ρ. Fluxes F=(Dv, Sv+p, S−Dv)=(Dv, Sv+p, (τ+p)v).
    cons2prim: 1-D root find on pressure (same z-style residual as the FV code).

    EXACT relativistic Riemann solution is provided by DGExactRiemann for the
    DG-vs-exact validation (no spurious oscillations, positivity, convergence).
=#
module DGSRHydro

using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis, rusanov_flux, tvb_minmod, cell_average

export SRGrid, SRDG, setup_srdg, set_initial!, evolve_srdg!, srdg_primitives,
       blast_initial!, shocktube_initial!, srdg_cell_means

# ---------------------------------------------------------------------------
# Ideal-gas SRHD primitive recovery: given (D,S,τ), Γ, find p by 1-D root find.
#   E = τ + p + D = ρhW² ;  v = S/E ;  W=1/√(1−v²) ;  ρ=D/W ;
#   ϵ = (E/W² − p − ρ)/ρ = (E/W² − p)/ρ − 1 ;  residual = p − (Γ−1)ρϵ.
# ---------------------------------------------------------------------------
@inline function _c2p_resid(p, D, S, τ, Γ)
    E = τ + p + D
    E ≤ 0 && return 1e30
    v2 = clamp((S/E)^2, 0.0, 1.0-1e-13)
    W = 1.0/sqrt(1.0-v2)
    ρ = D/W
    ρ ≤ 0 && return 1e30
    ϵ = (E/W^2 - p)/ρ - 1.0
    return p - (Γ-1.0)*ρ*ϵ
end

@inline function cons2prim_sr(D, S, τ, Γ)
    # bracket p in (plo, phi) and regula-falsi/bisection
    plo = 1e-14
    phi = max(1.0, abs(τ)+abs(S)+D)
    flo = _c2p_resid(plo, D, S, τ, Γ); fhi = _c2p_resid(phi, D, S, τ, Γ)
    it=0
    while flo*fhi > 0 && it < 80
        phi *= 3.0; fhi = _c2p_resid(phi, D, S, τ, Γ); it += 1
    end
    if flo*fhi > 0
        # fallback: cold/atmosphere
        return (max(D,1e-14), 0.0, 1e-13, 1.0)
    end
    p = phi
    for _ in 1:100
        p = (plo*fhi - phi*flo)/(fhi-flo)
        fp = _c2p_resid(p, D, S, τ, Γ)
        if abs(fp) < 1e-14*(1+abs(p)) || (phi-plo) < 1e-15
            break
        end
        if flo*fp < 0; phi=p; fhi=fp else plo=p; flo=fp end
    end
    E = τ + p + D
    v = S/E
    v2 = clamp(v^2, 0.0, 1.0-1e-13)
    W = 1.0/sqrt(1.0-v2)
    ρ = D/W
    return (ρ, v, p, W)
end

@inline function prim2cons_sr(ρ, v, p, Γ)
    ϵ = p/((Γ-1.0)*ρ)
    h = 1.0 + ϵ + p/ρ
    W = 1.0/sqrt(1.0 - v^2)
    D = ρ*W
    S = ρ*h*W^2*v
    τ = ρ*h*W^2 - p - D
    return (D, S, τ)
end

@inline function flux_sr(ρ, v, p, D, S, τ)
    return (D*v, S*v + p, (τ + p)*v)
end

# relativistic acoustic signal speed bound (for Rusanov amax)
@inline function _signal_speed(ρ, v, p, Γ)
    ϵ = p/((Γ-1.0)*ρ)
    h = 1.0 + ϵ + p/ρ
    cs2 = clamp(Γ*p/(ρ*h), 0.0, 1.0)
    cs = sqrt(cs2)
    v2 = v^2
    den = 1.0 - v2*cs2
    disc = sqrt(max(cs2*(1-v2)*(1-v2*cs2 - v2*(1-cs2)), 0.0))
    λp = (v*(1-cs2) + disc)/den
    λm = (v*(1-cs2) - disc)/den
    return max(abs(λp), abs(λm))
end

# ---------------------------------------------------------------------------
# Grid & state: K elements, each carrying (p+1) LGL nodes.
# ---------------------------------------------------------------------------
struct SRGrid
    K::Int            # number of elements
    xL::Float64
    xR::Float64
    Δx::Float64
    J::Float64        # Jacobian = Δx/2
    xc::Vector{Float64}        # element centers
    xnode::Matrix{Float64}     # (N, K) physical node coordinates
end

mutable struct SRDG
    g::SRGrid
    b::LGLBasis
    Γ::Float64
    cfl::Float64
    M_tvb::Float64    # TVB constant (0 ⇒ pure minmod TVD)
    # conserved nodal fields: (N, K)
    D::Matrix{Float64}
    S::Matrix{Float64}
    τ::Matrix{Float64}
    bc::Symbol        # :outflow (copy) for shock tubes
end

function SRGrid(K::Int, xL::Float64, xR::Float64, b::LGLBasis)
    Δx = (xR-xL)/K
    J = Δx/2
    xc = [xL + (k-0.5)*Δx for k in 1:K]
    N = b.N
    xnode = zeros(N, K)
    for k in 1:K, a in 1:N
        xnode[a,k] = xc[k] + J*b.ξ[a]
    end
    SRGrid(K, xL, xR, Δx, J, xc, xnode)
end

"""
    setup_srdg(; K=200, xL=0.0, xR=1.0, p=3, Γ=5/3, cfl=0.15, M_tvb=0.0) -> SRDG

Build a 1D special-relativistic hydro RKDG solver: `K` elements of degree `p`
LGL nodal basis on [xL,xR], ideal-gas Γ-law. `M_tvb` is the TVB troubled-cell
constant (0 = strict TVD minmod).
"""
function setup_srdg(; K::Int=200, xL::Float64=0.0, xR::Float64=1.0, p::Int=3,
                    Γ::Float64=5/3, cfl::Float64=0.15, M_tvb::Float64=0.0,
                    bc::Symbol=:outflow)
    b = build_lgl_basis(p)
    g = SRGrid(K, xL, xR, b)
    N = b.N
    z() = zeros(N, K)
    SRDG(g, b, Γ, cfl, M_tvb, z(), z(), z(), bc)
end

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
"""set_initial!(dg, ρ0, v0, p0) — ρ0,v0,p0 are functions of x."""
function set_initial!(dg::SRDG, ρ0, v0, p0)
    g = dg.g; N = dg.b.N
    @inbounds for k in 1:g.K, a in 1:N
        x = g.xnode[a,k]
        ρ = ρ0(x); v = v0(x); p = p0(x)
        D,S,τ = prim2cons_sr(ρ, v, p, dg.Γ)
        dg.D[a,k]=D; dg.S[a,k]=S; dg.τ[a,k]=τ
    end
    return dg
end

"""Standard relativistic shock tube: (ρL,vL,pL) for x<x0, (ρR,vR,pR) for x>x0."""
function shocktube_initial!(dg::SRDG; x0::Float64=0.5,
        ρL=1.0, vL=0.0, pL=1.0, ρR=0.125, vR=0.0, pR=0.1)
    set_initial!(dg, x->x<x0 ? ρL : ρR, x->x<x0 ? vL : vR, x->x<x0 ? pL : pR)
end

"""Mildly/strongly relativistic blast (high pL): default the classic
   Martí–Müller relativistic test (ρL=1,pL=1000 ... )."""
function blast_initial!(dg::SRDG; x0::Float64=0.5, ρL=1.0, vL=0.0, pL=1000.0,
                        ρR=1.0, vR=0.0, pR=0.01)
    set_initial!(dg, x->x<x0 ? ρL : ρR, x->x<x0 ? vL : vR, x->x<x0 ? pL : pR)
end

# ---------------------------------------------------------------------------
# Primitive recovery over the whole grid (for diagnostics / output)
# ---------------------------------------------------------------------------
function srdg_primitives(dg::SRDG)
    g = dg.g; N = dg.b.N
    ρ = similar(dg.D); v = similar(dg.D); pr = similar(dg.D)
    @inbounds for k in 1:g.K, a in 1:N
        ρa,va,pa,_ = cons2prim_sr(dg.D[a,k], dg.S[a,k], dg.τ[a,k], dg.Γ)
        ρ[a,k]=ρa; v[a,k]=va; pr[a,k]=pa
    end
    return ρ, v, pr
end

# cell means of the conserved (and primitive) fields (for plotting / error)
function srdg_cell_means(dg::SRDG)
    g=dg.g; w=dg.b.w; N=dg.b.N
    ρ,v,pr = srdg_primitives(dg)
    ρm=zeros(g.K); vm=zeros(g.K); pm=zeros(g.K); xm=copy(g.xc)
    @inbounds for k in 1:g.K
        ρm[k]=cell_average(view(ρ,:,k), w)
        vm[k]=cell_average(view(v,:,k), w)
        pm[k]=cell_average(view(pr,:,k), w)
    end
    return xm, ρm, vm, pm
end

# ---------------------------------------------------------------------------
# Positivity-preserving (Zhang–Shu) limiter: rescale nodal polynomial toward its
# cell mean so D and p stay positive at all nodes. Operates per element on the
# conserved variables (the admissible set {D>0, τ+D−√(D²+S²)>0 ⇒ p>0} is convex).
# ---------------------------------------------------------------------------
@inline _pfloor(D,S,τ) = (τ + D) - sqrt(D^2 + S^2)   # >0 ⇔ p>0 region proxy (q)

function _positivity_limit!(dg::SRDG, Dε::Float64, qε::Float64)
    g=dg.g; w=dg.b.w; N=dg.b.N
    @inbounds for k in 1:g.K
        D̄ = cell_average(view(dg.D,:,k), w)
        S̄ = cell_average(view(dg.S,:,k), w)
        τ̄ = cell_average(view(dg.τ,:,k), w)
        D̄ = max(D̄, Dε)
        # --- step 1: enforce D>0 at nodes ---
        Dmin = Inf
        for a in 1:N; Dmin=min(Dmin, dg.D[a,k]); end
        if Dmin < Dε
            θ = (D̄-Dε)/(D̄-Dmin + 1e-300)
            θ = clamp(θ, 0.0, 1.0)
            for a in 1:N; dg.D[a,k] = θ*(dg.D[a,k]-D̄)+D̄; end
        end
        # --- step 2: enforce q=p-proxy>0 at nodes ---
        q̄ = _pfloor(D̄,S̄,τ̄)
        if q̄ ≤ qε
            # mean itself non-positive: reset element to a tiny floor state
            for a in 1:N; dg.D[a,k]=max(D̄,Dε); dg.S[a,k]=S̄; dg.τ[a,k]=τ̄; end
            continue
        end
        θq = 1.0
        for a in 1:N
            q = _pfloor(dg.D[a,k], dg.S[a,k], dg.τ[a,k])
            if q < qε
                t = (q̄-qε)/(q̄-q + 1e-300)
                θq = min(θq, clamp(t,0.0,1.0))
            end
        end
        if θq < 1.0
            for a in 1:N
                dg.D[a,k]=θq*(dg.D[a,k]-D̄)+D̄
                dg.S[a,k]=θq*(dg.S[a,k]-S̄)+S̄
                dg.τ[a,k]=θq*(dg.τ[a,k]-τ̄)+τ̄
            end
        end
    end
end

# ---------------------------------------------------------------------------
# TVB-minmod troubled-cell limiter. For each element and each conserved field,
# compare the element's nodal slope against the minmod of neighbour-mean
# differences; if limited, REPLACE the nodal polynomial by the limited LINEAR
# profile (mean + limited slope·ξ). This is the classic Cockburn–Shu ΛΠ limiter.
# ---------------------------------------------------------------------------
function _tvb_limit_field!(U::Matrix{Float64}, dg::SRDG)
    g=dg.g; w=dg.b.w; ξ=dg.b.ξ; N=dg.b.N
    Mh2 = dg.M_tvb*g.Δx^2
    means = zeros(g.K)
    @inbounds for k in 1:g.K; means[k]=cell_average(view(U,:,k), w); end
    # nodal linear-slope coefficient: project U onto ξ. With LGL weights,
    # slope u1 = (3/2)Σ w_a ξ_a U_a (since ∫ξ²=2/3 over [-1,1]).
    newU = copy(U)
    @inbounds for k in 1:g.K
        ū = means[k]
        ūm = k>1   ? means[k-1] : means[k]
        ūp = k<g.K ? means[k+1] : means[k]
        # internal slope (value at ξ=+1 minus mean) in mean-difference units:
        # represent the cell's own undivided slope as u_right_face - mean.
        # use the LGL linear projection:
        s_int = 0.0
        for a in 1:N; s_int += 1.5*w[a]*ξ[a]*U[a,k]; end   # coeff of ξ
        # candidate differences to neighbours (forward/backward mean diffs)
        a1 = s_int
        a2 = ūp - ū
        a3 = ū - ūm
        sl = tvb_minmod(a1, a2, a3, Mh2)
        if abs(sl - s_int) > 1e-13*(1+abs(s_int))
            # flagged troubled: replace by limited linear profile
            for a in 1:N
                newU[a,k] = ū + sl*ξ[a]
            end
        end
    end
    copyto!(U, newU)
end

# detect troubled cells via the conserved-D field, then limit ALL fields there.
function _limit!(dg::SRDG)
    _tvb_limit_field!(dg.D, dg)
    _tvb_limit_field!(dg.S, dg)
    _tvb_limit_field!(dg.τ, dg)
end

# ---------------------------------------------------------------------------
# DG spatial operator: du/dt = L(u). Strong-form nodal DG with LGL mass lumping.
#   On element: (du/dt)_a = -(1/J) (D̃ f)_a + (1/(J w_a)) [endpoint flux jumps]
# where D̃ = D (differentiation), and the face correction injects (f̂ - f) at the
# two endpoints (a=1 left, a=N right) weighted by 1/w_a.
# ---------------------------------------------------------------------------
function _rhs!(rD,rS,rτ, dg::SRDG)
    g=dg.g; b=dg.b; N=b.N; K=g.K; J=g.J; Γ=dg.Γ
    fill!(rD,0.0); fill!(rS,0.0); fill!(rτ,0.0)
    # precompute nodal primitives & physical fluxes
    fD=zeros(N,K); fS=zeros(N,K); fτ=zeros(N,K)
    ρn=zeros(N,K); vn=zeros(N,K); pn=zeros(N,K)
    @inbounds for k in 1:K, a in 1:N
        ρ,v,p,_ = cons2prim_sr(dg.D[a,k],dg.S[a,k],dg.τ[a,k],Γ)
        ρn[a,k]=ρ; vn[a,k]=v; pn[a,k]=p
        f1,f2,f3 = flux_sr(ρ,v,p, dg.D[a,k],dg.S[a,k],dg.τ[a,k])
        fD[a,k]=f1; fS[a,k]=f2; fτ[a,k]=f3
    end
    # volume term: -(1/J) D̃ f   (matrix-vector per element)
    Dm = b.D
    @inbounds for k in 1:K
        for a in 1:N
            sD=0.0; sS=0.0; sτ=0.0
            for c in 1:N
                d = Dm[a,c]
                sD += d*fD[c,k]; sS += d*fS[c,k]; sτ += d*fτ[c,k]
            end
            rD[a,k] = -sD/J; rS[a,k] = -sS/J; rτ[a,k] = -sτ/J
        end
    end
    # surface term: numerical (Rusanov) flux at each interior + boundary face.
    # face k+1/2 between element k (right node N) and k+1 (left node 1).
    @inbounds for kf in 0:K
        kL = kf; kR = kf+1
        # boundary handling (outflow = copy)
        if kL == 0
            kL2 = 1; aL = 1
        else
            kL2 = kL; aL = N
        end
        if kR == K+1
            kR2 = K; aR = N
        else
            kR2 = kR; aR = 1
        end
        # left/right interface states
        if kL == 0
            # left boundary: ghost = copy of first node (outflow)
            DL=dg.D[1,1]; SL=dg.S[1,1]; τL=dg.τ[1,1]
        else
            DL=dg.D[N,kL]; SL=dg.S[N,kL]; τL=dg.τ[N,kL]
        end
        if kR == K+1
            DR=dg.D[N,K]; SR=dg.S[N,K]; τR=dg.τ[N,K]
        else
            DR=dg.D[1,kR]; SR=dg.S[1,kR]; τR=dg.τ[1,kR]
        end
        ρLf,vLf,pLf,_ = cons2prim_sr(DL,SL,τL,Γ)
        ρRf,vRf,pRf,_ = cons2prim_sr(DR,SR,τR,Γ)
        FL = flux_sr(ρLf,vLf,pLf, DL,SL,τL)
        FR = flux_sr(ρRf,vRf,pRf, DR,SR,τR)
        amax = max(_signal_speed(ρLf,vLf,pLf,Γ), _signal_speed(ρRf,vRf,pRf,Γ))
        F̂ = rusanov_flux((DL,SL,τL),FL,(DR,SR,τR),FR, amax)
        # inject into the two adjacent elements at their endpoints.
        # right node of element kL gets -(F̂ - F_int)/ (J w_N) with outward normal +1
        if kL ≥ 1
            wend = b.w[N]
            rD[N,kL] -= (F̂[1]-FL[1])/(J*wend)
            rS[N,kL] -= (F̂[2]-FL[2])/(J*wend)
            rτ[N,kL] -= (F̂[3]-FL[3])/(J*wend)
        end
        if kR ≤ K
            w1 = b.w[1]
            rD[1,kR] += (F̂[1]-FR[1])/(J*w1)
            rS[1,kR] += (F̂[2]-FR[2])/(J*w1)
            rτ[1,kR] += (F̂[3]-FR[3])/(J*w1)
        end
    end
    return nothing
end

@inline function _maxspeed(dg::SRDG)
    g=dg.g; N=dg.b.N; Γ=dg.Γ; a=0.0
    @inbounds for k in 1:g.K, i in 1:N
        ρ,v,p,_ = cons2prim_sr(dg.D[i,k],dg.S[i,k],dg.τ[i,k],Γ)
        a = max(a, _signal_speed(ρ,v,p,Γ))
    end
    return a
end

"""
    evolve_srdg!(dg; tmax, Dε=1e-13, qε=1e-13) -> nsteps

Advance the SRHD DG solution to time `tmax` with SSP-RK3 + TVB troubled-cell
limiter + Zhang–Shu positivity limiter applied after every stage. CFL is scaled
by 1/(2p+1) (the DG inverse-CFL).
"""
function evolve_srdg!(dg::SRDG; tmax::Float64, Dε::Float64=1e-13, qε::Float64=1e-13,
                      maxsteps::Int=2_000_000)
    g=dg.g; N=dg.b.N; K=g.K; p=dg.b.p
    rD=zeros(N,K); rS=zeros(N,K); rτ=zeros(N,K)
    D0=zeros(N,K); S0=zeros(N,K); τ0=zeros(N,K)
    cfl_dg = dg.cfl/(2p+1)
    t=0.0; nstep=0
    _limit!(dg); _positivity_limit!(dg, Dε, qε)
    while t<tmax && nstep<maxsteps
        amax=max(_maxspeed(dg),1e-6)
        dt = cfl_dg*g.Δx/amax
        dt = min(dt, tmax-t)
        @inbounds copyto!(D0,dg.D); copyto!(S0,dg.S); copyto!(τ0,dg.τ)
        # stage 1
        _rhs!(rD,rS,rτ, dg)
        @inbounds for I in eachindex(dg.D)
            dg.D[I]=D0[I]+dt*rD[I]; dg.S[I]=S0[I]+dt*rS[I]; dg.τ[I]=τ0[I]+dt*rτ[I]
        end
        _limit!(dg); _positivity_limit!(dg, Dε, qε)
        # stage 2
        _rhs!(rD,rS,rτ, dg)
        @inbounds for I in eachindex(dg.D)
            dg.D[I]=0.75*D0[I]+0.25*(dg.D[I]+dt*rD[I])
            dg.S[I]=0.75*S0[I]+0.25*(dg.S[I]+dt*rS[I])
            dg.τ[I]=0.75*τ0[I]+0.25*(dg.τ[I]+dt*rτ[I])
        end
        _limit!(dg); _positivity_limit!(dg, Dε, qε)
        # stage 3
        _rhs!(rD,rS,rτ, dg)
        @inbounds for I in eachindex(dg.D)
            dg.D[I]=(1/3)*D0[I]+(2/3)*(dg.D[I]+dt*rD[I])
            dg.S[I]=(1/3)*S0[I]+(2/3)*(dg.S[I]+dt*rS[I])
            dg.τ[I]=(1/3)*τ0[I]+(2/3)*(dg.τ[I]+dt*rτ[I])
        end
        _limit!(dg); _positivity_limit!(dg, Dε, qε)
        t+=dt; nstep+=1
    end
    return nstep
end

end # module DGSRHydro
