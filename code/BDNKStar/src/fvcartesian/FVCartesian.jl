#=
    FVCartesian — STAGE 2: 3-D Cartesian finite-volume Valencia GRHydro on the
    frozen TOV (Cowling) metric, with an ATMOSPHERE. The non-conforming spherical
    surface on the Cartesian grid is handled by the atmosphere (the standard cure
    for the staircase instability that destabilized the prior linear cowling3d
    engine). Goal: reproduce the ℓ=2 f-mode of the (r,θ) spherical engine from a
    totally different (Cartesian, nonlinear FV) discretization.

    Metric (frozen TOV, isotropic-in-angle areal Schwarzschild-like; Cowling):
        α = e^{ν/2},  γ_ij = δ_ij + (e^{λ}−1) n_i n_j,  n_i = x_i/r,
        √γ = e^{λ/2},  γ^{ij} = δ_ij + (e^{−λ}−1) n_i n_j.
    Physical squared velocity v² = γ_ij v^i v^j; conserved momentum is COVARIANT
    S_i = √γ ρhW² v_i. Conserved set per cell (densitized by √γ=e^{λ/2}):
        D = √γ ρW,  S_i = √γ ρhW² v_i,  τ = √γ(ρhW² − p − ρW).

    Flux along axis d: ∂_t U + ∂_d( α F^d ) = √γ S,  with F^d the standard
    Valencia fluxes built from the contravariant v^d = γ^{dk} v_k. We use a
    directional HLL solver with MinMod reconstruction, dimension-by-dimension,
    and SSP-RK2. Octant symmetry (x,y,z ≥ 0) with reflection BCs cuts the cost 8×.
    Well-balanced: the equilibrium flux+source is subtracted so the static star
    has exactly zero RHS (the same lake-at-rest trick as FVRadial).

    HONEST ASSESSMENT (verdict: fv-insufficient-need-dg). The atmosphere DOES cure
    the staircase for the UNPERTURBED star — the static Cartesian star is exactly
    well-balanced (machine-zero ρc drift). But under an ℓ=2 perturbation the non-
    conforming staircased surface drives a SECULAR drift whose growth ACCELERATES
    with resolution (ρc drift @t=300: ≈1.8% at N=32³ → 7.4% at N=48³), so the f-
    mode signal is swamped before even ~3 mode periods (one f-period ≈108 M⊙). At
    affordable octant resolution (N≲48) the 2nd-order FV cannot hold the perturbed
    surface cleanly enough to reproduce the 1.883 kHz ℓ=2 f-mode. A high-order
    DISCONTINUOUS-GALERKIN scheme (sub-cell resolution + a shock-capturing limiter
    that keeps the surface monotone without 1/Δx noise amplification) is the
    indicated next step — exactly mirroring why the boundary-conforming (r,θ)
    SphBDNK engine succeeded where the Cartesian one failed.
=#
module FVCartesian

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, minmod, rho_from_p, eos_cs2, _solve_p
using ..Numerics: brent
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram

export FVCartGrid, FVCartEngine, FVCartState, setup_fvcart, evolve_fvcart!,
       seed_l2_velocity!, fvcart_quadrupole, fvcart_central_density,
       fvcart_periodogram_freqs

# ---------------------------------------------------------------------------
# Grid: octant Cartesian, cell centers x_i=(i-1/2)Δ for i=1..N (so x≥0), NG ghosts.
# ---------------------------------------------------------------------------
struct FVCartGrid
    N::Int
    NG::Int
    Δ::Float64
    x::Vector{Float64}    # cell centers incl ghosts, length N+2NG; ghost on low (mirror) & high (atm)
    L::Float64
end
function FVCartGrid(N::Int, L::Float64; NG::Int=2)
    Δ = L/N
    ntot = N+2NG
    x = [ (i-NG-0.5)*Δ for i in 1:ntot ]
    FVCartGrid(N, NG, Δ, x, L)
end

# ---------------------------------------------------------------------------
# Background fields on cell centers (3D arrays, array-indexed incl ghosts)
# ---------------------------------------------------------------------------
struct FVCartBG
    eos::BarotropicEOS
    α::Array{Float64,3}
    elam::Array{Float64,3}      # e^{λ} = γ_rr
    sqrtγ::Array{Float64,3}     # √γ = e^{λ/2}
    nx::Array{Float64,3}; ny::Array{Float64,3}; nz::Array{Float64,3}
    r::Array{Float64,3}
    Φp::Array{Float64,3}        # Φ' = ν'/2
    interior::BitArray{3}       # r < R (star)
    R::Float64; M::Float64
end

mutable struct FVCartState
    D::Array{Float64,3}
    Sx::Array{Float64,3}; Sy::Array{Float64,3}; Sz::Array{Float64,3}
    τ::Array{Float64,3}
    # cached primitives
    ρ::Array{Float64,3}; p::Array{Float64,3}; ε::Array{Float64,3}
    vx::Array{Float64,3}; vy::Array{Float64,3}; vz::Array{Float64,3}
end
function FVCartState(ntot::Int)
    z()=zeros(ntot,ntot,ntot)
    FVCartState(z(),z(),z(),z(),z(), z(),z(),z(), z(),z(),z())
end

mutable struct FVCartEngine
    g::FVCartGrid
    bg::FVCartBG
    atm::AtmospherePars
    cfl::Float64
    Seq_D::Array{Float64,3}; Seq_Sx::Array{Float64,3}; Seq_Sy::Array{Float64,3}
    Seq_Sz::Array{Float64,3}; Seq_τ::Array{Float64,3}
end

@inline _lininterp(xs, ys, x) = begin
    n=length(xs)
    x≤xs[1] && return ys[1]
    x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j])
    ys[j]+t*(ys[j+1]-ys[j])
end

@inline _p_of_rho(eos::ShumPolytrope, ρ) = eos.κ*ρ^2

# ---------------------------------------------------------------------------
# cons2prim for the 3D Cartesian conserved set. Given densitized (D,S_i,τ) and
# metric (√γ, e^λ, n_i), recover primitives. The physical squared momentum:
#   Ŝ_i = S_i/√γ ; Ŝ² = γ^{ij}Ŝ_iŜ_j = |Ŝ|² + (e^{-λ}-1)(n·Ŝ)²
# then the same 1-D pressure root find as the radial case (orthonormal).
# Returns (ρ,p,ε, v^x,v^y,v^z) — contravariant velocity.
# ---------------------------------------------------------------------------
function _cons2prim_cell(eos, atm::AtmospherePars, D, Sx, Sy, Sz, τ, sqrtγ, elam, nx, ny, nz)
    if sqrtγ ≤ 0
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0,0.0,0.0)
    end
    D̂ = D/sqrtγ
    Ŝx = Sx/sqrtγ; Ŝy = Sy/sqrtγ; Ŝz = Sz/sqrtγ
    τ̂ = τ/sqrtγ
    if D̂ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0,0.0,0.0)
    end
    nS = nx*Ŝx + ny*Ŝy + nz*Ŝz
    # Ŝ² = γ^{ij}Ŝ_iŜ_j = δ + (e^{-λ}-1)(n·Ŝ)²
    S2 = Ŝx^2+Ŝy^2+Ŝz^2 + (1/elam - 1)*nS^2
    S2 = max(S2, 0.0)
    Ŝmag = sqrt(S2)
    pmax = 10.0*(abs(τ̂)+D̂+atm.p_atm)+1e-30
    p, ok = _solve_p(eos, D̂, Ŝmag, τ̂, 1e-30, pmax)
    if !ok || !isfinite(p)
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0,0.0,0.0)
    end
    p=max(p,atm.p_atm)
    E=τ̂+p+D̂
    v2=clamp((Ŝmag/E)^2,0.0,atm.vmax^2)
    W=1.0/sqrt(1.0-v2)
    ρ=D̂/W
    if ρ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0,0.0,0.0)
    end
    ε=E/W^2-p
    # contravariant velocity: v_i = Ŝ_i/(ρhW²) = Ŝ_i/E ; v^i = γ^{ij}v_j
    invE = E>0 ? 1.0/E : 0.0
    vlx=Ŝx*invE; vly=Ŝy*invE; vlz=Ŝz*invE      # covariant v_i
    nv = nx*vlx+ny*vly+nz*vlz
    q = (1/elam - 1)
    vux = vlx + q*nx*nv; vuy = vly + q*ny*nv; vuz = vlz + q*nz*nv   # contravariant v^i
    return (ρ, p, max(ε,0.0), vux, vuy, vuz)
end

# prim2cons for a cell: given ρ,p and contravariant v^i, build densitized cons.
function _prim2cons_cell(eos, ρ, p, vux, vuy, vuz, sqrtγ, elam, nx, ny, nz)
    ε = ρ>0 ? energy_from_pressure(eos,p) : 0.0
    h = ρ>0 ? (ε+p)/ρ : 1.0
    # covariant v_i = γ_ij v^j ; γ_ij = δ + (e^λ-1)n_in_j
    nv = nx*vux+ny*vuy+nz*vuz
    q = (elam-1)
    vlx = vux + q*nx*nv; vly = vuy + q*ny*nv; vlz = vuz + q*nz*nv
    v2 = vux*vlx+vuy*vly+vuz*vlz
    v2 = clamp(v2,0.0,1.0-1e-12)
    W = 1.0/sqrt(1.0-v2)
    D = sqrtγ*ρ*W
    fac = sqrtγ*ρ*h*W^2
    Sx = fac*vlx; Sy = fac*vly; Sz = fac*vlz
    τ = sqrtγ*(ρ*h*W^2 - p - ρ*W)
    return (D,Sx,Sy,Sz,τ)
end

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
"""
    setup_fvcart(eos, εc; N=48, L_fac=1.3, atm_fac=1e-7, cfl=0.25, h_tov=2e-4)
        -> (engine, state)

Build the octant 3-D Cartesian FV engine on the TOV star. `N` cells per axis over
[0, L_fac·R]. Returns equilibrium-initialized engine+state.
"""
function setup_fvcart(eos::BarotropicEOS, εc::Float64; N::Int=48, L_fac::Float64=1.3,
                      atm_fac::Float64=1e-7, cfl::Float64=0.25, h_tov::Float64=2e-4)
    star = solve_tov(eos, εc; h=h_tov)
    R,M = star.R, star.M
    L = L_fac*R
    g = FVCartGrid(N, L)
    ntot = N+2g.NG
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p
    m_of(r)= r≤R ? _lininterp(rt,mt,r) : M
    ν_of(r)= r≤R ? _lininterp(rt,νt,r) : log(1-2M/r)
    p_of(r)= r≤R ? max(_lininterp(rt,pt,r),0.0) : 0.0

    α=zeros(ntot,ntot,ntot); elam=similar(α); sqrtγ=similar(α)
    nx=similar(α); ny=similar(α); nz=similar(α); rr=similar(α); Φp=similar(α)
    interior=falses(ntot,ntot,ntot)
    @inbounds for k in 1:ntot, j in 1:ntot, i in 1:ntot
        X=g.x[i]; Y=g.x[j]; Z=g.x[k]
        r=sqrt(X^2+Y^2+Z^2); rr[i,j,k]=r
        if r ≤ 1e-12
            α[i,j,k]=exp(0.5*ν_of(0.0)); elam[i,j,k]=1.0; sqrtγ[i,j,k]=1.0
            nx[i,j,k]=ny[i,j,k]=nz[i,j,k]=0.0; Φp[i,j,k]=0.0; continue
        end
        m=m_of(r); fm=max(1-2m/r,1e-12)
        α[i,j,k]=exp(0.5*ν_of(r)); elam[i,j,k]=1/fm; sqrtγ[i,j,k]=1/sqrt(fm)
        nx[i,j,k]=X/r; ny[i,j,k]=Y/r; nz[i,j,k]=Z/r
        denom=r*(r-2m); Φp[i,j,k]= denom>0 ? (m+4π*r^3*p_of(r))/denom : 0.0
        interior[i,j,k] = r < R
    end
    bg=FVCartBG(eos,α,elam,sqrtγ,nx,ny,nz,rr,Φp,interior,R,M)

    ρc=rho_from_p(eos,pressure(eos,εc))
    ρ_atm=atm_fac*ρc; p_atm=_p_of_rho(eos,ρ_atm)
    ε_atm=energy_from_pressure(eos,p_atm)
    atm=AtmospherePars(ρ_atm,p_atm,ε_atm,5*ρ_atm,0.999)

    st=FVCartState(ntot)
    @inbounds for k in 1:ntot, j in 1:ntot, i in 1:ntot
        r=rr[i,j,k]
        ρ = r≤R ? max(rho_from_p(eos,p_of(r)),ρ_atm) : ρ_atm
        p = r≤R ? max(p_of(r),p_atm) : p_atm
        ε = energy_from_pressure(eos,p)
        st.ρ[i,j,k]=ρ; st.p[i,j,k]=p; st.ε[i,j,k]=ε
        st.vx[i,j,k]=0.0; st.vy[i,j,k]=0.0; st.vz[i,j,k]=0.0
        D,Sx,Sy,Sz,τ = _prim2cons_cell(eos,ρ,p,0.0,0.0,0.0,sqrtγ[i,j,k],elam[i,j,k],
                                       nx[i,j,k],ny[i,j,k],nz[i,j,k])
        st.D[i,j,k]=D; st.Sx[i,j,k]=Sx; st.Sy[i,j,k]=Sy; st.Sz[i,j,k]=Sz; st.τ[i,j,k]=τ
    end
    eng=FVCartEngine(g,bg,atm,cfl, zeros(ntot,ntot,ntot),zeros(ntot,ntot,ntot),
                     zeros(ntot,ntot,ntot),zeros(ntot,ntot,ntot),zeros(ntot,ntot,ntot))
    _update_primitives!(st,eng); _fill_ghosts!(st,eng)
    _raw_rhs!(eng.Seq_D,eng.Seq_Sx,eng.Seq_Sy,eng.Seq_Sz,eng.Seq_τ, st, eng)
    return eng, st
end

# ---------------------------------------------------------------------------
# primitive update over the whole grid
# ---------------------------------------------------------------------------
function _update_primitives!(st::FVCartState, eng::FVCartEngine)
    bg=eng.bg; eos=bg.eos; atm=eng.atm
    @inbounds Threads.@threads for k in axes(st.D,3)
        for j in axes(st.D,2), i in axes(st.D,1)
            ρ,p,ε,vx,vy,vz = _cons2prim_cell(eos,atm, st.D[i,j,k],st.Sx[i,j,k],
                st.Sy[i,j,k],st.Sz[i,j,k],st.τ[i,j,k], bg.sqrtγ[i,j,k],bg.elam[i,j,k],
                bg.nx[i,j,k],bg.ny[i,j,k],bg.nz[i,j,k])
            st.ρ[i,j,k]=ρ; st.p[i,j,k]=p; st.ε[i,j,k]=ε
            st.vx[i,j,k]=vx; st.vy[i,j,k]=vy; st.vz[i,j,k]=vz
            # re-densitize on atmosphere reset to keep cons consistent
            if ρ ≤ atm.ρ_cut*1.0000001 && p ≤ atm.p_atm*1.0000001
                D,Sx,Sy,Sz,τ=_prim2cons_cell(eos,ρ,p,0.0,0.0,0.0,bg.sqrtγ[i,j,k],
                    bg.elam[i,j,k],bg.nx[i,j,k],bg.ny[i,j,k],bg.nz[i,j,k])
                st.D[i,j,k]=D; st.Sx[i,j,k]=Sx; st.Sy[i,j,k]=Sy; st.Sz[i,j,k]=Sz; st.τ[i,j,k]=τ
            end
        end
    end
end

# ghost fill: low side mirror (octant symmetry, even scalars, odd normal velocity),
# high side outflow→atmosphere.
function _fill_ghosts!(st::FVCartState, eng::FVCartEngine)
    # operate on primitives, then re-densitize the ghost conserved variables.
    # low side: octant mirror; high side: outflow → atmosphere.
    _ghost_axis!(st,eng,1); _ghost_axis!(st,eng,2); _ghost_axis!(st,eng,3)
end

function _ghost_axis!(st::FVCartState, eng::FVCartEngine, ax::Int)
    g=eng.g; NG=g.NG; N=g.N; bg=eng.bg; eos=bg.eos; atm=eng.atm
    ntot=N+2NG
    rng = 1:ntot
    @inbounds for gc in 1:NG
        lo = NG-gc+1; mir = NG+gc        # low mirror
        hi = NG+N+gc; src = NG+N         # high outflow
        for c in rng, b in rng
            # build index tuples
            Ilo = _idx(ax, lo, b, c); Imir = _idx(ax, mir, b, c)
            Ihi = _idx(ax, hi, b, c)
            # low: octant mirror — scalars even, velocity component ax flips
            st.ρ[Ilo...]=st.ρ[Imir...]; st.p[Ilo...]=st.p[Imir...]; st.ε[Ilo...]=st.ε[Imir...]
            st.vx[Ilo...]=st.vx[Imir...]; st.vy[Ilo...]=st.vy[Imir...]; st.vz[Ilo...]=st.vz[Imir...]
            _flipv!(st, Ilo, ax)
            _redens!(st,eng,Ilo)
            # high: atmosphere
            st.ρ[Ihi...]=atm.ρ_atm; st.p[Ihi...]=atm.p_atm; st.ε[Ihi...]=atm.ε_atm
            st.vx[Ihi...]=0.0; st.vy[Ihi...]=0.0; st.vz[Ihi...]=0.0
            _redens!(st,eng,Ihi)
        end
    end
end

@inline _idx(ax,a,b,c) = ax==1 ? (a,b,c) : ax==2 ? (b,a,c) : (b,c,a)
@inline function _flipv!(st,I,ax)
    if ax==1; st.vx[I...]=-st.vx[I...]
    elseif ax==2; st.vy[I...]=-st.vy[I...]
    else; st.vz[I...]=-st.vz[I...] end
end
@inline function _redens!(st,eng,I)
    bg=eng.bg; eos=bg.eos
    D,Sx,Sy,Sz,τ=_prim2cons_cell(eos,st.ρ[I...],st.p[I...],st.vx[I...],st.vy[I...],st.vz[I...],
        bg.sqrtγ[I...],bg.elam[I...],bg.nx[I...],bg.ny[I...],bg.nz[I...])
    st.D[I...]=D; st.Sx[I...]=Sx; st.Sy[I...]=Sy; st.Sz[I...]=Sz; st.τ[I...]=τ
end

# ---------------------------------------------------------------------------
# Physical (orthonormal-along-axis) fluxes & wave speeds for HLL along axis ax.
# We use the standard Valencia flux in the contravariant velocity component
# vᵃ = v^ax along the coordinate axis (orthonormal value vᵃ_phys = √γ_aa vᵃ).
# For simplicity & robustness we form the directional flux in the orthonormal
# frame using the projected normal speed wᵃ = √(γ_aa) v^ax (physical), since the
# atmosphere/surface dynamics is dominated by the radial direction handled by
# the metric n_i. Here γ_aa = 1 + (e^λ-1)n_a² (diagonal metric component).
# ---------------------------------------------------------------------------
# wave speeds along axis: λ± = α-frame; in orthonormal radial-dominated approx,
# use the relativistic 1D formula with the physical normal speed w and sound cs.
@inline function _axis_speeds(cs2, w)
    cs=sqrt(clamp(cs2,0.0,1.0))
    w2=clamp(w^2,0.0,1.0-1e-12)
    a=1.0/(1.0-w2*cs^2)
    disc=sqrt(max(cs^2*(1-w2)*(1-w2*cs^2-w^2*(1-cs^2)),0.0))
    return a*(w*(1-cs^2)-disc), a*(w*(1-cs^2)+disc)
end

# MinMod reconstruct a 3D primitive array along axis ax to the face between
# cells (left=index a, right=a+1). Returns (qL,qR).
@inline function _recon3(q, ax, i,j,k)
    if ax==1
        sL=minmod(q[i,j,k]-q[i-1,j,k], q[i+1,j,k]-q[i,j,k])
        sR=minmod(q[i+1,j,k]-q[i,j,k], q[i+2,j,k]-q[i+1,j,k])
        return q[i,j,k]+0.5*sL, q[i+1,j,k]-0.5*sR
    elseif ax==2
        sL=minmod(q[i,j,k]-q[i,j-1,k], q[i,j+1,k]-q[i,j,k])
        sR=minmod(q[i,j+1,k]-q[i,j,k], q[i,j+2,k]-q[i,j+1,k])
        return q[i,j,k]+0.5*sL, q[i,j+1,k]-0.5*sR
    else
        sL=minmod(q[i,j,k]-q[i,j,k-1], q[i,j,k+1]-q[i,j,k])
        sR=minmod(q[i,j,k+1]-q[i,j,k], q[i,j,k+2]-q[i,j,k+1])
        return q[i,j,k]+0.5*sL, q[i,j,k+1]-0.5*sR
    end
end

# Build the FULL covariant conserved + flux state at a reconstructed face value.
# Given primitives (ρ,p,vux,vuy,vuz) and metric at the face (sqrtγ_face implicit
# =1 in undensitized orthonormal; we work undensitized and multiply α√γ outside).
# Returns U=(D̂,Ŝx,Ŝy,Ŝz,τ̂) covariant and F=(FD,FSx,FSy,FSz,Fτ) along axis ax,
# plus the physical normal speed w for the wave-speed estimate.
@inline function _face_state(eos, ρ,p,vux,vuy,vuz, elam,nx,ny,nz, ax)
    ε=energy_from_pressure(eos,p)
    h=ρ>0 ? (ε+p)/ρ : 1.0
    nv=nx*vux+ny*vuy+nz*vuz
    q=(elam-1)
    vlx=vux+q*nx*nv; vly=vuy+q*ny*nv; vlz=vuz+q*nz*nv  # covariant v_i
    v2=clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0-1e-12)
    W=1.0/sqrt(1.0-v2)
    D̂=ρ*W
    fac=ρ*h*W^2
    Ŝx=fac*vlx; Ŝy=fac*vly; Ŝz=fac*vlz
    τ̂=fac - p - ρ*W
    # contravariant velocity along axis ax: va = v^ax (the flux advection speed)
    va = ax==1 ? vux : ax==2 ? vuy : vuz
    # fluxes F^ax: F_D = D̂ v^ax ; F_{S_i} = Ŝ_i v^ax + p δ^ax_i ; F_τ=(τ̂+p)v^ax
    FD=D̂*va
    δx = ax==1 ? 1.0 : 0.0; δy=ax==2 ? 1.0 : 0.0; δz=ax==3 ? 1.0 : 0.0
    FSx=Ŝx*va + p*δx; FSy=Ŝy*va + p*δy; FSz=Ŝz*va + p*δz
    Fτ=(τ̂+p)*va
    # physical normal speed (for wave estimate): w = √(γ_aa)·|v^ax| sign; use va
    # times √γ_aa where γ_aa=1+(e^λ-1)n_a²
    na = ax==1 ? nx : ax==2 ? ny : nz
    γaa = 1+(elam-1)*na^2
    w = sqrt(γaa)*va
    return (D̂,Ŝx,Ŝy,Ŝz,τ̂),(FD,FSx,FSy,FSz,Fτ), w, ε
end

# ---------------------------------------------------------------------------
# Raw RHS: directional HLL sweep + analytic geometric source.
# ---------------------------------------------------------------------------
function _raw_rhs!(rD,rSx,rSy,rSz,rτ, st::FVCartState, eng::FVCartEngine)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG; Δ=g.Δ; atm=eng.atm
    fill!(rD,0.0); fill!(rSx,0.0); fill!(rSy,0.0); fill!(rSz,0.0); fill!(rτ,0.0)
    lo=NG+1; hi=NG+N
    # interior faces along each axis. For face between a and a+1, area = α√γ at face
    # (averaged). We add contributions to both cells.
    for ax in 1:3
        @inbounds Threads.@threads for k in lo:hi
            for j in lo:hi
                for i in lo:hi
                    # face on the high side of cell (i,j,k) along ax
                    (ia,ja,ka) = ax==1 ? (i+1,j,k) : ax==2 ? (i,j+1,k) : (i,j,k+1)
                    ia>hi+1 && continue   # only up to last interior face
                    # reconstruct primitives to the face (left=this cell, right=neighbor)
                    ρL,ρR=_recon3(st.ρ,ax,i,j,k); pL,pR=_recon3(st.p,ax,i,j,k)
                    vxL,vxR=_recon3(st.vx,ax,i,j,k); vyL,vyR=_recon3(st.vy,ax,i,j,k); vzL,vzR=_recon3(st.vz,ax,i,j,k)
                    ρL=max(ρL,atm.ρ_atm); ρR=max(ρR,atm.ρ_atm)
                    pL=max(pL,atm.p_atm); pR=max(pR,atm.p_atm)
                    # metric at the face: average of the two cells
                    elamf=0.5*(bg.elam[i,j,k]+bg.elam[ia,ja,ka])
                    nxf=0.5*(bg.nx[i,j,k]+bg.nx[ia,ja,ka]); nyf=0.5*(bg.ny[i,j,k]+bg.ny[ia,ja,ka]); nzf=0.5*(bg.nz[i,j,k]+bg.nz[ia,ja,ka])
                    αf=0.5*(bg.α[i,j,k]+bg.α[ia,ja,ka]); sgf=0.5*(bg.sqrtγ[i,j,k]+bg.sqrtγ[ia,ja,ka])
                    UL,FL,wL,εL=_face_state(eos,ρL,pL,vxL,vyL,vzL,elamf,nxf,nyf,nzf,ax)
                    UR,FR,wR,εR=_face_state(eos,ρR,pR,vxR,vyR,vzR,elamf,nxf,nyf,nzf,ax)
                    csL=eos_cs2(eos,εL); csR=eos_cs2(eos,εR)
                    smL,spL=_axis_speeds(csL,wL); smR,spR=_axis_speeds(csR,wR)
                    sL=min(smL,smR,0.0); sR=max(spL,spR,0.0)
                    A=αf*sgf
                    # HLL flux (5-component)
                    invs = sR>sL ? 1.0/(sR-sL) : 0.0
                    @inline hll(u_l,u_r,f_l,f_r) = sL≥0 ? f_l : sR≤0 ? f_r : (sR*f_l-sL*f_r+sL*sR*(u_r-u_l))*invs
                    fD =A*hll(UL[1],UR[1],FL[1],FR[1])
                    fSx=A*hll(UL[2],UR[2],FL[2],FR[2])
                    fSy=A*hll(UL[3],UR[3],FL[3],FR[3])
                    fSz=A*hll(UL[4],UR[4],FL[4],FR[4])
                    fτ =A*hll(UL[5],UR[5],FL[5],FR[5])
                    # subtract from left cell, add to right cell
                    rD[i,j,k]-=fD/Δ; rSx[i,j,k]-=fSx/Δ; rSy[i,j,k]-=fSy/Δ; rSz[i,j,k]-=fSz/Δ; rτ[i,j,k]-=fτ/Δ
                    rD[ia,ja,ka]+=fD/Δ; rSx[ia,ja,ka]+=fSx/Δ; rSy[ia,ja,ka]+=fSy/Δ; rSz[ia,ja,ka]+=fSz/Δ; rτ[ia,ja,ka]+=fτ/Δ
                end
            end
        end
    end
    # geometric/gravity source: S_{S_i} = √γ[ −(ε+p)W² Φ' α n_i ] (radial gravity);
    # S_τ = −√γ α (ε+p)W² (v·n) Φ'. (Cartesian projection of the spherical source.)
    @inbounds Threads.@threads for k in lo:hi
        for j in lo:hi, i in lo:hi
            α=bg.α[i,j,k]; sg=bg.sqrtγ[i,j,k]; Φp=bg.Φp[i,j,k]
            nx=bg.nx[i,j,k]; ny=bg.ny[i,j,k]; nz=bg.nz[i,j,k]
            ρ=st.ρ[i,j,k]; p=st.p[i,j,k]; ε=st.ε[i,j,k]
            vux=st.vx[i,j,k]; vuy=st.vy[i,j,k]; vuz=st.vz[i,j,k]
            # covariant v_i for v·n and W
            elam=bg.elam[i,j,k]; q=(elam-1); nv=nx*vux+ny*vuy+nz*vuz
            vlx=vux+q*nx*nv; vly=vuy+q*ny*nv; vlz=vuz+q*nz*nv
            v2=clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
            vn = nx*vlx+ny*vly+nz*vlz     # v·n (covariant·normal)
            gforce = -α*sg*(ε+p)*W^2*Φp
            rSx[i,j,k]+=gforce*nx; rSy[i,j,k]+=gforce*ny; rSz[i,j,k]+=gforce*nz
            rτ[i,j,k]+= -α*sg*(ε+p)*W^2*vn*Φp
        end
    end
    return nothing
end

function _rhs!(rD,rSx,rSy,rSz,rτ, st::FVCartState, eng::FVCartEngine)
    _raw_rhs!(rD,rSx,rSy,rSz,rτ, st, eng)
    @inbounds for I in eachindex(rD)
        rD[I]-=eng.Seq_D[I]; rSx[I]-=eng.Seq_Sx[I]; rSy[I]-=eng.Seq_Sy[I]
        rSz[I]-=eng.Seq_Sz[I]; rτ[I]-=eng.Seq_τ[I]
    end
end

# ---------------------------------------------------------------------------
# time stepping (SSP-RK2)
# ---------------------------------------------------------------------------
function _max_speed(st::FVCartState, eng::FVCartEngine)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG
    lo=NG+1; hi=NG+N; a=0.0
    @inbounds for k in lo:hi, j in lo:hi, i in lo:hi
        cs=sqrt(clamp(eos_cs2(eos,st.ε[i,j,k]),0.0,1.0))
        # crude: |v|+cs in physical terms
        elam=bg.elam[i,j,k]; nx=bg.nx[i,j,k]; ny=bg.ny[i,j,k]; nz=bg.nz[i,j,k]
        vux=st.vx[i,j,k]; vuy=st.vy[i,j,k]; vuz=st.vz[i,j,k]
        q=elam-1; nv=nx*vux+ny*vuy+nz*vuz
        vlx=vux+q*nx*nv; vly=vuy+q*ny*nv; vlz=vuz+q*nz*nv
        vmag=sqrt(clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0))
        a=max(a, (vmag+cs)/(1+vmag*cs)*sqrt(elam))   # radial-worst physical speed
    end
    return a
end

"""
    evolve_fvcart!(st, eng; tmax, sample_dt, A_probe...) -> (ts, q_l2, ρc)

Advance the 3-D Cartesian FV star to coordinate time `tmax` (SSP-RK2), recording
the ℓ=2 quadrupole moment and the central density every ~`sample_dt`. Returns
(times, q_ℓ2(t), ρ_c(t)). `drift` available via ρc history.
"""
function evolve_fvcart!(st::FVCartState, eng::FVCartEngine; tmax::Float64,
                        sample_dt::Float64=-1.0, cfl::Float64=-1.0, verbose::Bool=false)
    g=eng.g; bg=eng.bg; N=g.N; NG=g.NG; Δ=g.Δ; ntot=N+2NG
    cfl=cfl>0 ? cfl : eng.cfl
    sample_dt = sample_dt>0 ? sample_dt : 3*Δ
    z()=zeros(ntot,ntot,ntot)
    rD,rSx,rSy,rSz,rτ = z(),z(),z(),z(),z()
    Dn,Sxn,Syn,Szn,τn = z(),z(),z(),z(),z()
    ts=Float64[]; q2=Float64[]; ρch=Float64[]
    ic=NG+1; ρc0=st.ρ[ic,ic,ic]
    t=0.0; last=-1e30; nstep=0
    _update_primitives!(st,eng); _fill_ghosts!(st,eng)
    while t<tmax && nstep<5_000_000
        amax=max(_max_speed(st,eng),1e-3)
        dt=cfl*Δ/amax; dt=min(dt,tmax-t)
        @inbounds for I in eachindex(st.D)
            Dn[I]=st.D[I]; Sxn[I]=st.Sx[I]; Syn[I]=st.Sy[I]; Szn[I]=st.Sz[I]; τn[I]=st.τ[I]
        end
        _rhs!(rD,rSx,rSy,rSz,rτ, st, eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=Dn[I]+dt*rD[I]; st.Sx[I]=Sxn[I]+dt*rSx[I]; st.Sy[I]=Syn[I]+dt*rSy[I]
            st.Sz[I]=Szn[I]+dt*rSz[I]; st.τ[I]=τn[I]+dt*rτ[I]
        end
        _update_primitives!(st,eng); _fill_ghosts!(st,eng)
        _rhs!(rD,rSx,rSy,rSz,rτ, st, eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=0.5*(Dn[I]+st.D[I]+dt*rD[I]); st.Sx[I]=0.5*(Sxn[I]+st.Sx[I]+dt*rSx[I])
            st.Sy[I]=0.5*(Syn[I]+st.Sy[I]+dt*rSy[I]); st.Sz[I]=0.5*(Szn[I]+st.Sz[I]+dt*rSz[I])
            st.τ[I]=0.5*(τn[I]+st.τ[I]+dt*rτ[I])
        end
        _update_primitives!(st,eng); _fill_ghosts!(st,eng)
        t+=dt; nstep+=1
        if t-last≥sample_dt
            push!(ts,t); push!(q2,fvcart_quadrupole(st,eng)); push!(ρch,st.ρ[ic,ic,ic])
            last=t
            verbose && nstep%200==0 && println("  t=",round(t,digits=1)," ρc/ρc0=",round(st.ρ[ic,ic,ic]/ρc0,digits=5))
        end
        if !isfinite(st.ρ[ic,ic,ic]) || st.ρ[ic,ic,ic]>100*ρc0
            return ts,q2,ρch    # blow-up
        end
    end
    return ts,q2,ρch
end

"""seed an ℓ=2,m=0 velocity perturbation: v^r = A (r/R)(2n_z²−n_x²−n_y²)."""
function seed_l2_velocity!(st::FVCartState, eng::FVCartEngine; A::Float64=1e-3)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG; lo=NG+1; hi=NG+N
    @inbounds for k in lo:hi, j in lo:hi, i in lo:hi
        bg.interior[i,j,k] || continue
        r=bg.r[i,j,k]; r≤1e-12 && continue
        nx=bg.nx[i,j,k]; ny=bg.ny[i,j,k]; nz=bg.nz[i,j,k]
        Y20 = 2*nz^2 - nx^2 - ny^2
        vr = A*(r/bg.R)*Y20
        # contravariant v^i = vr * n^i ; n^i (contravariant) = γ^{ij}n_j; since n is
        # radial unit covector, n^i = e^{-λ}? For simplicity set v^i along n_i (covariant
        # direction) and let prim2cons handle the metric — use v^i = vr_phys n_i / √γ_rr.
        # Physical radial velocity vr_phys → contravariant v^r = vr_phys/√(γ_rr)=vr_phys e^{-λ/2}.
        el2=sqrt(bg.elam[i,j,k])
        vux=vr/el2*nx; vuy=vr/el2*ny; vuz=vr/el2*nz
        st.vx[i,j,k]=vux; st.vy[i,j,k]=vuy; st.vz[i,j,k]=vuz
        D,Sx,Sy,Sz,τ=_prim2cons_cell(eos,st.ρ[i,j,k],st.p[i,j,k],vux,vuy,vuz,
            bg.sqrtγ[i,j,k],bg.elam[i,j,k],nx,ny,nz)
        st.D[i,j,k]=D; st.Sx[i,j,k]=Sx; st.Sy[i,j,k]=Sy; st.Sz[i,j,k]=Sz; st.τ[i,j,k]=τ
    end
    _fill_ghosts!(st,eng)
end

"""√γ-weighted ℓ=2,m=0 quadrupole of (ρ−ρ_bg): Σ δρ (2n_z²−n_x²−n_y²) √γ Δ³."""
function fvcart_quadrupole(st::FVCartState, eng::FVCartEngine)
    g=eng.g; bg=eng.bg; N=g.N; NG=g.NG; lo=NG+1; hi=NG+N; acc=0.0
    @inbounds for k in lo:hi, j in lo:hi, i in lo:hi
        bg.interior[i,j,k] || continue
        r=bg.r[i,j,k]; r≤1e-12 && continue
        nx=bg.nx[i,j,k]; ny=bg.ny[i,j,k]; nz=bg.nz[i,j,k]
        acc += st.ρ[i,j,k]*(2nz^2-nx^2-ny^2)*bg.sqrtγ[i,j,k]
    end
    acc*g.Δ^3
end

fvcart_central_density(st::FVCartState, eng::FVCartEngine) =
    (c=eng.g.NG+1; st.ρ[c,c,c])

"""periodogram peak frequencies (kHz) of a quadrupole time series."""
function fvcart_periodogram_freqs(ts::Vector{Float64}, q0::Vector{Float64};
        fmin_kHz::Float64=0.5, fmax_kHz::Float64=10.0, npts::Int=4000,
        Lunit_km::Float64=Msun_to_km, npeaks::Int=4, window::Bool=true)
    νmin=fmin_kHz*Lunit_km*kHz_to_km; νmax=fmax_kHz*Lunit_km*kHz_to_km
    νgrid=range(νmin,νmax;length=npts)
    q=copy(q0)
    if window
        n=length(q); μ=sum(q)/n
        @inbounds for i in 1:n; q[i]=(q[i]-μ)*(0.5-0.5*cos(2π*(i-1)/(n-1))); end
    end
    P=periodogram(ts,q,νgrid)
    peaks=Tuple{Float64,Float64}[]
    for n in 2:length(P)-1
        (P[n]>P[n-1] && P[n]≥P[n+1]) && push!(peaks,(P[n],νgrid[n]))
    end
    sort!(peaks,by=x->-x[1]); k=min(npeaks,length(peaks))
    freqs=[p[2]/(Lunit_km*kHz_to_km) for p in peaks[1:k]]
    return freqs, collect(P), collect(νgrid)
end

end # module FVCartesian
