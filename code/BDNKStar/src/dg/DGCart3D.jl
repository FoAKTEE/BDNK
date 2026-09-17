#=
    DGCart3D — STAGE 3 (3+1D): FULL three-dimensional Cartesian tensor-product
    nodal RKDG GRHydro on the frozen TOV (Cowling) metric, with the
    NON-CONFORMING stellar surface + atmosphere + positivity/TVB limiter.

    This is the genuine 3+1D extension of the validated 2+1D axisymmetric engine
    DGCart2D (which was meridional x–z, m=0 ONLY). Here the solution lives on a
    full 3D (x,y,z) octant grid and carries a FULL 3-momentum (S_x,S_y,S_z), so a
    genuinely NON-AXISYMMETRIC m≠0 perturbation (Y_22, Y_21) can be represented —
    the thing 2D-axisymmetric could not do.

    REUSE: this module reuses the proven 1D LGL kernels (DGCommon), the
    atmosphere/cons2prim root-find machinery (FVCommon._solve_p, AtmospherePars,
    eos_cs2), the TOV background and EOS, exactly as DGCart2D does. Everything is
    lifted dimension-by-dimension: the same densitized Valencia conserved set, the
    same directional Rusanov sweeps (now x THEN y THEN z), the same SSP-RK3, the
    same per-element TVB-troubled / Zhang–Shu positivity limiter, and the same
    well-balanced lake-at-rest (store & subtract the static background DG RHS).

    EQUILIBRIUM PRESERVATION (wb=true, the default for a star). The RHS
    subtraction makes rhs(U_eq) ≡ 0 exactly, but two other pieces of a standard
    RKDG scheme are NOT equilibrium-preserving and, on a star whose surface
    element spans many orders of magnitude in D, each one alone destroys the
    background:
      (i)  the mean-based Zhang–Shu limiter U ← Ū + θ(U−Ū) flattens the
           equilibrium's own intra-element profile whenever θ<1, pushing mass
           from the stellar edge into the atmosphere nodes (measured on the
           unseeded star: thrown into a ±3% radial oscillation with 0.6c surface
           velocities in 50 M⊙; ρ_c −10% in 4.6 f-mode periods once seeded);
      (ii) the Rusanov dissipation a_max(U_R−U_L) acts on the equilibrium jump
           across each face with a STATE-DEPENDENT coefficient, so any change of
           the local wave speed sources a first-order error proportional to the
           (large) equilibrium surface jump.
    Both are fixed by working on the deviation δU = U − U_eq: the limiter
    (_limit_wb!) scales toward a mean-preserving reference that reduces to U_eq
    when δU=0 and shares the mean deviation out in proportion to the background,
    and the flux (_rus5wb) dissipates a_max on δU but only a frozen equilibrium
    sound speed on U_eq. U_eq is then a BITWISE fixed point of limiter+RHS; the
    unseeded star stays static to |ρ_c/ρ_c0−1| ~ 1e-13 and |v| ~ 4e-11 over
    250 M⊙ with the limiter ON (it was −10% eroding before); a 1% velocity seed
    (A=1e-2, K=6) leaves ρ_c at −5.8e-4 by t=600 and −7.3e-4 by t=1000,
    decelerating, with no stellar element ever handed to the fallback limiter
    (dgcart3d_limiter_census). The stored U_eq is the nodal projection of the
    TOV star used for the initial data, so this is exact for that data and
    reduces to the plain scheme (wb=false, used by the shock tube) when no
    background is stored. What the fix does NOT remove is the surface
    systematic of a nodal DG star at 14–23 nodes across R: the ℓ=2 f-mode comes
    out 2–10% below the 1D Cowling value with estimator scatter of the same
    size (VALIDATION.md §7.8), the analogue of the linear engine's masked
    (rigid-wall) surface. For precision use the cut-cell linear engines.

    Conserved (densitized, FULL 3-momentum):
        D=√γ ρW, S_i=√γ ρhW²v_i (i=x,y,z), τ=√γ(ρhW²−p−ρW), √γ=e^{λ/2}.
    Cowling metric γ_ij=δ_ij+(e^λ−1)n_in_j with n=(x,y,z)/r.

    Storage: 6D arrays (N,N,N,Kx,Ky,Kz). Octant symmetry (reflecting x=0,y=0,z=0)
    is used for the m=0 staircase test (cheap); for the m≠0 (Y_22/Y_21) test we
    keep the SAME octant grid but verify the perturbation's allowed reflection
    parity (Y_22 ∝ n_x²−n_y² and Y_21 ∝ n_x n_z are even under x→−x AND y→−y, so
    the octant box captures one lobe consistently).

    Keep resolution MODEST (K=6..10 per dim, p=2..3): 3D at ≤6 threads is
    expensive — short evolutions suffice for stability/trend/oscillation, NOT long
    f-mode frequency extraction.
=#
module DGCart3D

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, rho_from_p, eos_cs2, _solve_p

using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis, tvb_minmod

export DGCart3DEngine, DGCart3DState, setup_dgcart3d, evolve_dgcart3d!,
       seed_dgcart3d_l2!, seed_dgcart3d_Y22!, seed_dgcart3d_Y21!,
       dgcart3d_central_density, dgcart3d_quadrupole, dgcart3d_quadrupole_m2,
       dgcart3d_shocktube_diagonal!, dgcart3d_prim_minmax, dgcart3d_limiter_census

@inline _p_of_rho(eos::ShumPolytrope, ρ) = eos.κ*ρ^2
@inline function _lininterp(xs, ys, x)
    n = length(xs); x ≤ xs[1] && return ys[1]; x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x); t = (x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end

# cell-average of an (N,N,N) nodal block with tensor LGL weights (Σw=2 per dim).
@inline function _cellavg3(blk, w)
    s = 0.0; N = length(w)
    @inbounds for c in 1:N, b in 1:N, a in 1:N
        s += w[a]*w[b]*w[c]*blk[a,b,c]
    end
    return 0.125*s
end

struct DGCart3DEngine
    bs::LGLBasis
    Kx::Int; Ky::Int; Kz::Int
    Δ::Float64; J::Float64
    eos::BarotropicEOS
    atm::AtmospherePars
    R::Float64; M::Float64
    cfl::Float64; M_tvb::Float64
    # per-node background, arrays (N,N,N,Kx,Ky,Kz)
    x::Array{Float64,6}; y::Array{Float64,6}; z::Array{Float64,6}; r::Array{Float64,6}
    α::Array{Float64,6}; elam::Array{Float64,6}; sqrtγ::Array{Float64,6}
    nx::Array{Float64,6}; ny::Array{Float64,6}; nz::Array{Float64,6}; Φp::Array{Float64,6}
    interior::BitArray{6}
    Seq_D::Array{Float64,6}; Seq_Sx::Array{Float64,6}; Seq_Sy::Array{Float64,6}
    Seq_Sz::Array{Float64,6}; Seq_τ::Array{Float64,6}
    # WELL-BALANCED EQUILIBRIUM (wb=true): the static conserved state U_eq and its primitives.
    # Both the positivity limiter and the Rusanov dissipation act on the DEVIATION U − U_eq, so
    # U_eq is an exact fixed point of the whole update, not only of the RHS subtraction.
    wb::Bool
    Deq::Array{Float64,6}; Sxeq::Array{Float64,6}; Syeq::Array{Float64,6}
    Szeq::Array{Float64,6}; τeq::Array{Float64,6}
    ρeq::Array{Float64,6}; peq::Array{Float64,6}
end

mutable struct DGCart3DState
    D::Array{Float64,6}; Sx::Array{Float64,6}; Sy::Array{Float64,6}; Sz::Array{Float64,6}; τ::Array{Float64,6}
    ρ::Array{Float64,6}; p::Array{Float64,6}; ε::Array{Float64,6}
    vx::Array{Float64,6}; vy::Array{Float64,6}; vz::Array{Float64,6}
end

# ---------------------------------------------------------------------------
# cons2prim for the full 3D Cowling metric γ_ij=δ_ij+(e^λ−1)n_in_j.
# γ^{ij}=δ_ij+(1/e^λ−1)n_in_j ; |S|²_phys = γ^{ij}Ŝ_iŜ_j .
# ---------------------------------------------------------------------------
@inline function _c2p(eos, atm, D,Sx,Sy,Sz,τ, sqrtγ,elam,nx,ny,nz)
    if sqrtγ ≤ 0; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0,0.0) end
    D̂=D/sqrtγ; Ŝx=Sx/sqrtγ; Ŝy=Sy/sqrtγ; Ŝz=Sz/sqrtγ; τ̂=τ/sqrtγ
    if D̂ ≤ atm.ρ_cut; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0,0.0) end
    nS=nx*Ŝx+ny*Ŝy+nz*Ŝz
    S2=Ŝx^2+Ŝy^2+Ŝz^2+(1/elam-1)*nS^2; S2=max(S2,0.0); Ŝ=sqrt(S2)
    pmax=10.0*(abs(τ̂)+D̂+atm.p_atm)+1e-30
    p,ok=_solve_p(eos,D̂,Ŝ,τ̂,1e-30,pmax)
    if !ok || !isfinite(p); return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0,0.0) end
    p=max(p,atm.p_atm); E=τ̂+p+D̂
    v2=clamp((Ŝ/E)^2,0.0,atm.vmax^2); W=1.0/sqrt(1.0-v2); ρ=D̂/W
    if ρ ≤ atm.ρ_cut; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0,0.0) end
    ε=E/W^2-p; invE = E>0 ? 1.0/E : 0.0
    vlx=Ŝx*invE; vly=Ŝy*invE; vlz=Ŝz*invE; nv=nx*vlx+ny*vly+nz*vlz; q=(1/elam-1)
    vux=vlx+q*nx*nv; vuy=vly+q*ny*nv; vuz=vlz+q*nz*nv
    return (ρ,p,max(ε,0.0),vux,vuy,vuz)
end

function _p2c(eos, ρ,p,vux,vuy,vuz, sqrtγ,elam,nx,ny,nz)
    ε = ρ>0 ? energy_from_pressure(eos,p) : 0.0; h = ρ>0 ? (ε+p)/ρ : 1.0
    nv=nx*vux+ny*vuy+nz*vuz; q=(elam-1)
    vlx=vux+q*nx*nv; vly=vuy+q*ny*nv; vlz=vuz+q*nz*nv
    v2=clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
    D=sqrtγ*ρ*W; fac=sqrtγ*ρ*h*W^2
    return (D, fac*vlx, fac*vly, fac*vlz, sqrtγ*(ρ*h*W^2-p-ρ*W))
end

# undensitized face state + flux along axis (1=x,2=y,3=z); returns (U,F,w-frame-speed,ε)
@inline function _faceflux(eos, ρ,p,vux,vuy,vuz, elam,nx,ny,nz, ax)
    ε=energy_from_pressure(eos,p); h = ρ>0 ? (ε+p)/ρ : 1.0
    nv=nx*vux+ny*vuy+nz*vuz; q=(elam-1)
    vlx=vux+q*nx*nv; vly=vuy+q*ny*nv; vlz=vuz+q*nz*nv
    v2=clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
    D̂=ρ*W; fac=ρ*h*W^2; Ŝx=fac*vlx; Ŝy=fac*vly; Ŝz=fac*vlz; τ̂=fac-p-ρ*W
    va = ax==1 ? vux : (ax==2 ? vuy : vuz)
    δx = ax==1 ? 1.0 : 0.0; δy = ax==2 ? 1.0 : 0.0; δz = ax==3 ? 1.0 : 0.0
    FD=D̂*va; FSx=Ŝx*va+p*δx; FSy=Ŝy*va+p*δy; FSz=Ŝz*va+p*δz; Fτ=(τ̂+p)*va
    na = ax==1 ? nx : (ax==2 ? ny : nz); γaa=1+(elam-1)*na^2; w=sqrt(γaa)*va
    return (D̂,Ŝx,Ŝy,Ŝz,τ̂),(FD,FSx,FSy,FSz,Fτ), w, ε
end
@inline function _axisspeeds(cs2,w)
    cs=sqrt(clamp(cs2,0.0,1.0)); w2=clamp(w^2,0.0,1.0-1e-12)
    a=1.0/(1.0-w2*cs^2); disc=sqrt(max(cs^2*(1-w2)*(1-w2*cs^2-w^2*(1-cs^2)),0.0))
    return a*(w*(1-cs^2)-disc), a*(w*(1-cs^2)+disc)
end

"""
    setup_dgcart3d(eos, εc; Kx=8, Ky=8, Kz=8, p=2, L_fac=1.3, atm_fac=1e-7,
                   cfl=0.2, M_tvb=50.0, h_tov=2e-4) -> (engine, state)

Build the full-3D octant DG engine on the TOV star (reflecting symmetry on
x=0,y=0,z=0). `Kx,Ky,Kz` elements of degree `p`; effective resolution per axis
≈ K·p. Modest defaults for ≤6 threads.
"""
function setup_dgcart3d(eos::BarotropicEOS, εc::Float64; Kx::Int=8, Ky::Int=8, Kz::Int=8,
        p::Int=2, L_fac::Float64=1.3, atm_fac::Float64=1e-7, cfl::Float64=0.2,
        M_tvb::Float64=50.0, h_tov::Float64=2e-4)
    star=solve_tov(eos,εc;h=h_tov); R,M=star.R,star.M
    L=L_fac*R; Δ=L/Kx; J=Δ/2
    bs=build_lgl_basis(p); N=bs.N
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p
    m_of(r)= r≤R ? _lininterp(rt,mt,r) : M
    ν_of(r)= r≤R ? _lininterp(rt,νt,r) : log(1-2M/r)
    p_of(r)= r≤R ? max(_lininterp(rt,pt,r),0.0) : 0.0

    dims=(N,N,N,Kx,Ky,Kz)
    x=zeros(dims); y=zeros(dims); z=zeros(dims); r=zeros(dims)
    α=zeros(dims); elam=zeros(dims); sqrtγ=zeros(dims)
    nx=zeros(dims); ny=zeros(dims); nz=zeros(dims); Φp=zeros(dims)
    interior=falses(dims)
    @inbounds Threads.@threads for kz in 1:Kz
        for ky in 1:Ky, kx in 1:Kx
            xc=(kx-0.5)*Δ; yc=(ky-0.5)*Δ; zc=(kz-0.5)*Δ
            for c in 1:N, b in 1:N, a in 1:N
                X=xc+J*bs.ξ[a]; Y=yc+J*bs.ξ[b]; Z=zc+J*bs.ξ[c]; rr=sqrt(X^2+Y^2+Z^2)
                x[a,b,c,kx,ky,kz]=X; y[a,b,c,kx,ky,kz]=Y; z[a,b,c,kx,ky,kz]=Z; r[a,b,c,kx,ky,kz]=rr
                if rr ≤ 1e-10
                    α[a,b,c,kx,ky,kz]=exp(0.5*ν_of(0.0)); elam[a,b,c,kx,ky,kz]=1.0; sqrtγ[a,b,c,kx,ky,kz]=1.0
                    continue
                end
                m=m_of(rr); fm=max(1-2m/rr,1e-12)
                α[a,b,c,kx,ky,kz]=exp(0.5*ν_of(rr)); elam[a,b,c,kx,ky,kz]=1/fm; sqrtγ[a,b,c,kx,ky,kz]=1/sqrt(fm)
                nx[a,b,c,kx,ky,kz]=X/rr; ny[a,b,c,kx,ky,kz]=Y/rr; nz[a,b,c,kx,ky,kz]=Z/rr
                den=rr*(rr-2m); Φp[a,b,c,kx,ky,kz]= den>0 ? (m+4π*rr^3*p_of(rr))/den : 0.0
                interior[a,b,c,kx,ky,kz]= rr<R
            end
        end
    end

    ρc=rho_from_p(eos,pressure(eos,εc)); ρ_atm=atm_fac*ρc; p_atm=_p_of_rho(eos,ρ_atm)
    ε_atm=energy_from_pressure(eos,p_atm)
    atm=AtmospherePars(ρ_atm,p_atm,ε_atm,5*ρ_atm,0.999)

    st=DGCart3DState((zeros(dims) for _ in 1:11)...)
    @inbounds Threads.@threads for kz in 1:Kz
        for ky in 1:Ky, kx in 1:Kx, c in 1:N, b in 1:N, a in 1:N
            rr=r[a,b,c,kx,ky,kz]
            ρ = rr≤R ? max(rho_from_p(eos,p_of(rr)),ρ_atm) : ρ_atm
            pp = rr≤R ? max(p_of(rr),p_atm) : p_atm
            st.ρ[a,b,c,kx,ky,kz]=ρ; st.p[a,b,c,kx,ky,kz]=pp; st.ε[a,b,c,kx,ky,kz]=energy_from_pressure(eos,pp)
            D,Sx,Sy,Sz,τ=_p2c(eos,ρ,pp,0.0,0.0,0.0,sqrtγ[a,b,c,kx,ky,kz],elam[a,b,c,kx,ky,kz],
                nx[a,b,c,kx,ky,kz],ny[a,b,c,kx,ky,kz],nz[a,b,c,kx,ky,kz])
            st.D[a,b,c,kx,ky,kz]=D; st.Sx[a,b,c,kx,ky,kz]=Sx; st.Sy[a,b,c,kx,ky,kz]=Sy
            st.Sz[a,b,c,kx,ky,kz]=Sz; st.τ[a,b,c,kx,ky,kz]=τ
        end
    end
    eng=DGCart3DEngine(bs,Kx,Ky,Kz,Δ,J,eos,atm,R,M,cfl,M_tvb,x,y,z,r,α,elam,sqrtγ,nx,ny,nz,Φp,interior,
                       zeros(dims),zeros(dims),zeros(dims),zeros(dims),zeros(dims),
                       true, copy(st.D),copy(st.Sx),copy(st.Sy),copy(st.Sz),copy(st.τ), zeros(dims),zeros(dims))
    _update_prims!(st,eng)
    # the equilibrium PRIMITIVES are the post-floor ones the raw RHS actually sees
    copyto!(eng.ρeq, st.ρ); copyto!(eng.peq, st.p)
    # Seq must be built with the deviation-form flux at U=U_eq (deviation ≡ 0 there), so that
    # raw_rhs(U_eq) − Seq vanishes identically for every later a_max(U)
    _raw_rhs!(eng.Seq_D,eng.Seq_Sx,eng.Seq_Sy,eng.Seq_Sz,eng.Seq_τ, st, eng)
    return eng, st
end

function _update_prims!(st::DGCart3DState, eng::DGCart3DEngine)
    eos=eng.eos; atm=eng.atm; N=eng.bs.N
    @inbounds Threads.@threads for kz in 1:eng.Kz
        for ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
            ρ,p,ε,vx,vy,vz=_c2p(eos,atm,st.D[a,b,c,kx,ky,kz],st.Sx[a,b,c,kx,ky,kz],st.Sy[a,b,c,kx,ky,kz],
                st.Sz[a,b,c,kx,ky,kz],st.τ[a,b,c,kx,ky,kz],eng.sqrtγ[a,b,c,kx,ky,kz],eng.elam[a,b,c,kx,ky,kz],
                eng.nx[a,b,c,kx,ky,kz],eng.ny[a,b,c,kx,ky,kz],eng.nz[a,b,c,kx,ky,kz])
            st.ρ[a,b,c,kx,ky,kz]=ρ; st.p[a,b,c,kx,ky,kz]=p; st.ε[a,b,c,kx,ky,kz]=ε
            st.vx[a,b,c,kx,ky,kz]=vx; st.vy[a,b,c,kx,ky,kz]=vy; st.vz[a,b,c,kx,ky,kz]=vz
        end
    end
end

@inline function _rus5(UL,FL,UR,FR,amax,Af)
    ntuple(i->0.5*Af*(FL[i]+FR[i])-0.5*amax*Af*(UR[i]-UL[i]), 5)
end
# WELL-BALANCED Rusanov: the dissipation acts on the deviation from equilibrium with the
# CURRENT speed, and on the equilibrium's own inter-element jump with the FROZEN equilibrium
# speed. At U=U_eq the second term is all there is, and it is exactly what Seq contains — so
# the subtraction cancels for any a_max(U). With plain Rusanov, a_max changes as soon as the
# star moves and the equilibrium jump (the DG interpolation error of the steep surface) leaks
# into the RHS in proportion to a(U)−a(U_eq): a rectified, secular forcing on the background.
@inline function _rus5wb(UL,FL,UR,FR,amax,Af, ULe,URe,ae)
    ntuple(i->0.5*Af*(FL[i]+FR[i]) - 0.5*Af*( amax*((UR[i]-URe[i])-(UL[i]-ULe[i])) + ae*(URe[i]-ULe[i]) ), 5)
end

function _raw_rhs!(rD,rSx,rSy,rSz,rτ, st::DGCart3DState, eng::DGCart3DEngine)
    bs=eng.bs; N=bs.N; Kx=eng.Kx; Ky=eng.Ky; Kz=eng.Kz; J=eng.J; eos=eng.eos; atm=eng.atm
    Dm=bs.D; w=bs.w
    fill!(rD,0.0); fill!(rSx,0.0); fill!(rSy,0.0); fill!(rSz,0.0); fill!(rτ,0.0)

    # VOLUME term (all 3 directions), element-local flux buffers per thread.
    @inbounds Threads.@threads for kz in 1:Kz
        FxD=zeros(N,N,N); FxSx=zeros(N,N,N); FxSy=zeros(N,N,N); FxSz=zeros(N,N,N); Fxτ=zeros(N,N,N)
        FyD=zeros(N,N,N); FySx=zeros(N,N,N); FySy=zeros(N,N,N); FySz=zeros(N,N,N); Fyτ=zeros(N,N,N)
        FzD=zeros(N,N,N); FzSx=zeros(N,N,N); FzSy=zeros(N,N,N); FzSz=zeros(N,N,N); Fzτ=zeros(N,N,N)
        for ky in 1:Ky, kx in 1:Kx
            for c in 1:N, b in 1:N, a in 1:N
                ρ=st.ρ[a,b,c,kx,ky,kz]; p=st.p[a,b,c,kx,ky,kz]
                vx=st.vx[a,b,c,kx,ky,kz]; vy=st.vy[a,b,c,kx,ky,kz]; vz=st.vz[a,b,c,kx,ky,kz]
                el=eng.elam[a,b,c,kx,ky,kz]
                nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
                A=eng.α[a,b,c,kx,ky,kz]*eng.sqrtγ[a,b,c,kx,ky,kz]
                _,Fx,_,_=_faceflux(eos,ρ,p,vx,vy,vz,el,nxv,nyv,nzv,1)
                _,Fy,_,_=_faceflux(eos,ρ,p,vx,vy,vz,el,nxv,nyv,nzv,2)
                _,Fz,_,_=_faceflux(eos,ρ,p,vx,vy,vz,el,nxv,nyv,nzv,3)
                FxD[a,b,c]=A*Fx[1];FxSx[a,b,c]=A*Fx[2];FxSy[a,b,c]=A*Fx[3];FxSz[a,b,c]=A*Fx[4];Fxτ[a,b,c]=A*Fx[5]
                FyD[a,b,c]=A*Fy[1];FySx[a,b,c]=A*Fy[2];FySy[a,b,c]=A*Fy[3];FySz[a,b,c]=A*Fy[4];Fyτ[a,b,c]=A*Fy[5]
                FzD[a,b,c]=A*Fz[1];FzSx[a,b,c]=A*Fz[2];FzSy[a,b,c]=A*Fz[3];FzSz[a,b,c]=A*Fz[4];Fzτ[a,b,c]=A*Fz[5]
            end
            for c in 1:N, b in 1:N, a in 1:N
                sD=0.0;sSx=0.0;sSy=0.0;sSz=0.0;sτ=0.0
                for d in 1:N
                    dx=Dm[a,d]; dy=Dm[b,d]; dz=Dm[c,d]
                    sD +=dx*FxD[d,b,c] +dy*FyD[a,d,c] +dz*FzD[a,b,d]
                    sSx+=dx*FxSx[d,b,c]+dy*FySx[a,d,c]+dz*FzSx[a,b,d]
                    sSy+=dx*FxSy[d,b,c]+dy*FySy[a,d,c]+dz*FzSy[a,b,d]
                    sSz+=dx*FxSz[d,b,c]+dy*FySz[a,d,c]+dz*FzSz[a,b,d]
                    sτ +=dx*Fxτ[d,b,c] +dy*Fyτ[a,d,c] +dz*Fzτ[a,b,d]
                end
                sg=eng.sqrtγ[a,b,c,kx,ky,kz]; iJg=1.0/(J*sg)
                rD[a,b,c,kx,ky,kz]=-sD*iJg; rSx[a,b,c,kx,ky,kz]=-sSx*iJg
                rSy[a,b,c,kx,ky,kz]=-sSy*iJg; rSz[a,b,c,kx,ky,kz]=-sSz*iJg; rτ[a,b,c,kx,ky,kz]=-sτ*iJg
            end
        end
    end

    # SURFACE term — X faces (loop kz outer for thread safety: each thread owns a kz slab)
    @inbounds Threads.@threads for kz in 1:Kz
        for ky in 1:Ky, kf in 0:Kx
            kL=kf; kR=kf+1
            for c in 1:N, b in 1:N
                if kL==0
                    ρR=st.ρ[1,b,c,1,ky,kz];pR=st.p[1,b,c,1,ky,kz];vxR=st.vx[1,b,c,1,ky,kz];vyR=st.vy[1,b,c,1,ky,kz];vzR=st.vz[1,b,c,1,ky,kz]
                    ρL=ρR;pL=pR;vxL=-vxR;vyL=vyR;vzL=vzR
                    ρeR=eng.ρeq[1,b,c,1,ky,kz];peR=eng.peq[1,b,c,1,ky,kz];ρeL=ρeR;peL=peR
                    elf=eng.elam[1,b,c,1,ky,kz];nxf=eng.nx[1,b,c,1,ky,kz];nyf=eng.ny[1,b,c,1,ky,kz];nzf=eng.nz[1,b,c,1,ky,kz]
                    Af=eng.α[1,b,c,1,ky,kz]*eng.sqrtγ[1,b,c,1,ky,kz]
                elseif kR==Kx+1
                    ρL=st.ρ[N,b,c,Kx,ky,kz];pL=st.p[N,b,c,Kx,ky,kz];vxL=st.vx[N,b,c,Kx,ky,kz];vyL=st.vy[N,b,c,Kx,ky,kz];vzL=st.vz[N,b,c,Kx,ky,kz]
                    ρR=atm.ρ_atm;pR=atm.p_atm;vxR=0.0;vyR=0.0;vzR=0.0
                    ρeL=eng.ρeq[N,b,c,Kx,ky,kz];peL=eng.peq[N,b,c,Kx,ky,kz];ρeR=atm.ρ_atm;peR=atm.p_atm
                    elf=eng.elam[N,b,c,Kx,ky,kz];nxf=eng.nx[N,b,c,Kx,ky,kz];nyf=eng.ny[N,b,c,Kx,ky,kz];nzf=eng.nz[N,b,c,Kx,ky,kz]
                    Af=eng.α[N,b,c,Kx,ky,kz]*eng.sqrtγ[N,b,c,Kx,ky,kz]
                else
                    ρL=st.ρ[N,b,c,kL,ky,kz];pL=st.p[N,b,c,kL,ky,kz];vxL=st.vx[N,b,c,kL,ky,kz];vyL=st.vy[N,b,c,kL,ky,kz];vzL=st.vz[N,b,c,kL,ky,kz]
                    ρR=st.ρ[1,b,c,kR,ky,kz];pR=st.p[1,b,c,kR,ky,kz];vxR=st.vx[1,b,c,kR,ky,kz];vyR=st.vy[1,b,c,kR,ky,kz];vzR=st.vz[1,b,c,kR,ky,kz]
                    ρeL=eng.ρeq[N,b,c,kL,ky,kz];peL=eng.peq[N,b,c,kL,ky,kz];ρeR=eng.ρeq[1,b,c,kR,ky,kz];peR=eng.peq[1,b,c,kR,ky,kz]
                    elf=0.5*(eng.elam[N,b,c,kL,ky,kz]+eng.elam[1,b,c,kR,ky,kz])
                    nxf=0.5*(eng.nx[N,b,c,kL,ky,kz]+eng.nx[1,b,c,kR,ky,kz]); nyf=0.5*(eng.ny[N,b,c,kL,ky,kz]+eng.ny[1,b,c,kR,ky,kz])
                    nzf=0.5*(eng.nz[N,b,c,kL,ky,kz]+eng.nz[1,b,c,kR,ky,kz])
                    Af=0.5*(eng.α[N,b,c,kL,ky,kz]*eng.sqrtγ[N,b,c,kL,ky,kz]+eng.α[1,b,c,kR,ky,kz]*eng.sqrtγ[1,b,c,kR,ky,kz])
                end
                UL,FL,wL,εL=_faceflux(eos,ρL,pL,vxL,vyL,vzL,elf,nxf,nyf,nzf,1)
                UR,FR,wR,εR=_faceflux(eos,ρR,pR,vxR,vyR,vzR,elf,nxf,nyf,nzf,1)
                amax=max(abs.(_axisspeeds(eos_cs2(eos,εL),wL))...,abs.(_axisspeeds(eos_cs2(eos,εR),wR))...)
                if eng.wb
                    ULe,_,_,εeL=_faceflux(eos,ρeL,peL,0.0,0.0,0.0,elf,nxf,nyf,nzf,1)
                    URe,_,_,εeR=_faceflux(eos,ρeR,peR,0.0,0.0,0.0,elf,nxf,nyf,nzf,1)
                    ae=max(sqrt(clamp(eos_cs2(eos,εeL),0.0,1.0)),sqrt(clamp(eos_cs2(eos,εeR),0.0,1.0)))
                    F̂=_rus5wb(UL,FL,UR,FR,amax,Af,ULe,URe,ae)
                else
                    F̂=_rus5(UL,FL,UR,FR,amax,Af)
                end
                FLp=(Af.*FL); FRp=(Af.*FR)
                if kL≥1
                    wend=w[N]; sg=eng.sqrtγ[N,b,c,kL,ky,kz]; f=1.0/(J*wend*sg)
                    rD[N,b,c,kL,ky,kz]-=(F̂[1]-FLp[1])*f; rSx[N,b,c,kL,ky,kz]-=(F̂[2]-FLp[2])*f
                    rSy[N,b,c,kL,ky,kz]-=(F̂[3]-FLp[3])*f; rSz[N,b,c,kL,ky,kz]-=(F̂[4]-FLp[4])*f; rτ[N,b,c,kL,ky,kz]-=(F̂[5]-FLp[5])*f
                end
                if kR≤Kx
                    w1=w[1]; sg=eng.sqrtγ[1,b,c,kR,ky,kz]; f=1.0/(J*w1*sg)
                    rD[1,b,c,kR,ky,kz]+=(F̂[1]-FRp[1])*f; rSx[1,b,c,kR,ky,kz]+=(F̂[2]-FRp[2])*f
                    rSy[1,b,c,kR,ky,kz]+=(F̂[3]-FRp[3])*f; rSz[1,b,c,kR,ky,kz]+=(F̂[4]-FRp[4])*f; rτ[1,b,c,kR,ky,kz]+=(F̂[5]-FRp[5])*f
                end
            end
        end
    end

    # SURFACE term — Y faces (thread over kz; Y-face couples ky/ky+1 within a kz slab)
    @inbounds Threads.@threads for kz in 1:Kz
        for kx in 1:Kx, kf in 0:Ky
            kL=kf; kR=kf+1
            for c in 1:N, a in 1:N
                if kL==0
                    ρR=st.ρ[a,1,c,kx,1,kz];pR=st.p[a,1,c,kx,1,kz];vxR=st.vx[a,1,c,kx,1,kz];vyR=st.vy[a,1,c,kx,1,kz];vzR=st.vz[a,1,c,kx,1,kz]
                    ρL=ρR;pL=pR;vxL=vxR;vyL=-vyR;vzL=vzR
                    ρeR=eng.ρeq[a,1,c,kx,1,kz];peR=eng.peq[a,1,c,kx,1,kz];ρeL=ρeR;peL=peR
                    elf=eng.elam[a,1,c,kx,1,kz];nxf=eng.nx[a,1,c,kx,1,kz];nyf=eng.ny[a,1,c,kx,1,kz];nzf=eng.nz[a,1,c,kx,1,kz]
                    Af=eng.α[a,1,c,kx,1,kz]*eng.sqrtγ[a,1,c,kx,1,kz]
                elseif kR==Ky+1
                    ρL=st.ρ[a,N,c,kx,Ky,kz];pL=st.p[a,N,c,kx,Ky,kz];vxL=st.vx[a,N,c,kx,Ky,kz];vyL=st.vy[a,N,c,kx,Ky,kz];vzL=st.vz[a,N,c,kx,Ky,kz]
                    ρR=atm.ρ_atm;pR=atm.p_atm;vxR=0.0;vyR=0.0;vzR=0.0
                    ρeL=eng.ρeq[a,N,c,kx,Ky,kz];peL=eng.peq[a,N,c,kx,Ky,kz];ρeR=atm.ρ_atm;peR=atm.p_atm
                    elf=eng.elam[a,N,c,kx,Ky,kz];nxf=eng.nx[a,N,c,kx,Ky,kz];nyf=eng.ny[a,N,c,kx,Ky,kz];nzf=eng.nz[a,N,c,kx,Ky,kz]
                    Af=eng.α[a,N,c,kx,Ky,kz]*eng.sqrtγ[a,N,c,kx,Ky,kz]
                else
                    ρL=st.ρ[a,N,c,kx,kL,kz];pL=st.p[a,N,c,kx,kL,kz];vxL=st.vx[a,N,c,kx,kL,kz];vyL=st.vy[a,N,c,kx,kL,kz];vzL=st.vz[a,N,c,kx,kL,kz]
                    ρR=st.ρ[a,1,c,kx,kR,kz];pR=st.p[a,1,c,kx,kR,kz];vxR=st.vx[a,1,c,kx,kR,kz];vyR=st.vy[a,1,c,kx,kR,kz];vzR=st.vz[a,1,c,kx,kR,kz]
                    ρeL=eng.ρeq[a,N,c,kx,kL,kz];peL=eng.peq[a,N,c,kx,kL,kz];ρeR=eng.ρeq[a,1,c,kx,kR,kz];peR=eng.peq[a,1,c,kx,kR,kz]
                    elf=0.5*(eng.elam[a,N,c,kx,kL,kz]+eng.elam[a,1,c,kx,kR,kz])
                    nxf=0.5*(eng.nx[a,N,c,kx,kL,kz]+eng.nx[a,1,c,kx,kR,kz]); nyf=0.5*(eng.ny[a,N,c,kx,kL,kz]+eng.ny[a,1,c,kx,kR,kz])
                    nzf=0.5*(eng.nz[a,N,c,kx,kL,kz]+eng.nz[a,1,c,kx,kR,kz])
                    Af=0.5*(eng.α[a,N,c,kx,kL,kz]*eng.sqrtγ[a,N,c,kx,kL,kz]+eng.α[a,1,c,kx,kR,kz]*eng.sqrtγ[a,1,c,kx,kR,kz])
                end
                UL,FL,wL,εL=_faceflux(eos,ρL,pL,vxL,vyL,vzL,elf,nxf,nyf,nzf,2)
                UR,FR,wR,εR=_faceflux(eos,ρR,pR,vxR,vyR,vzR,elf,nxf,nyf,nzf,2)
                amax=max(abs.(_axisspeeds(eos_cs2(eos,εL),wL))...,abs.(_axisspeeds(eos_cs2(eos,εR),wR))...)
                if eng.wb
                    ULe,_,_,εeL=_faceflux(eos,ρeL,peL,0.0,0.0,0.0,elf,nxf,nyf,nzf,2)
                    URe,_,_,εeR=_faceflux(eos,ρeR,peR,0.0,0.0,0.0,elf,nxf,nyf,nzf,2)
                    ae=max(sqrt(clamp(eos_cs2(eos,εeL),0.0,1.0)),sqrt(clamp(eos_cs2(eos,εeR),0.0,1.0)))
                    F̂=_rus5wb(UL,FL,UR,FR,amax,Af,ULe,URe,ae)
                else
                    F̂=_rus5(UL,FL,UR,FR,amax,Af)
                end
                FLp=(Af.*FL); FRp=(Af.*FR)
                if kL≥1
                    wend=w[N]; sg=eng.sqrtγ[a,N,c,kx,kL,kz]; f=1.0/(J*wend*sg)
                    rD[a,N,c,kx,kL,kz]-=(F̂[1]-FLp[1])*f; rSx[a,N,c,kx,kL,kz]-=(F̂[2]-FLp[2])*f
                    rSy[a,N,c,kx,kL,kz]-=(F̂[3]-FLp[3])*f; rSz[a,N,c,kx,kL,kz]-=(F̂[4]-FLp[4])*f; rτ[a,N,c,kx,kL,kz]-=(F̂[5]-FLp[5])*f
                end
                if kR≤Ky
                    w1=w[1]; sg=eng.sqrtγ[a,1,c,kx,kR,kz]; f=1.0/(J*w1*sg)
                    rD[a,1,c,kx,kR,kz]+=(F̂[1]-FRp[1])*f; rSx[a,1,c,kx,kR,kz]+=(F̂[2]-FRp[2])*f
                    rSy[a,1,c,kx,kR,kz]+=(F̂[3]-FRp[3])*f; rSz[a,1,c,kx,kR,kz]+=(F̂[4]-FRp[4])*f; rτ[a,1,c,kx,kR,kz]+=(F̂[5]-FRp[5])*f
                end
            end
        end
    end

    # SURFACE term — Z faces (thread over kx; Z-face couples kz/kz+1 within a kx slab)
    @inbounds Threads.@threads for kx in 1:Kx
        for ky in 1:Ky, kf in 0:Kz
            kL=kf; kR=kf+1
            for b in 1:N, a in 1:N
                if kL==0
                    ρR=st.ρ[a,b,1,kx,ky,1];pR=st.p[a,b,1,kx,ky,1];vxR=st.vx[a,b,1,kx,ky,1];vyR=st.vy[a,b,1,kx,ky,1];vzR=st.vz[a,b,1,kx,ky,1]
                    ρL=ρR;pL=pR;vxL=vxR;vyL=vyR;vzL=-vzR
                    ρeR=eng.ρeq[a,b,1,kx,ky,1];peR=eng.peq[a,b,1,kx,ky,1];ρeL=ρeR;peL=peR
                    elf=eng.elam[a,b,1,kx,ky,1];nxf=eng.nx[a,b,1,kx,ky,1];nyf=eng.ny[a,b,1,kx,ky,1];nzf=eng.nz[a,b,1,kx,ky,1]
                    Af=eng.α[a,b,1,kx,ky,1]*eng.sqrtγ[a,b,1,kx,ky,1]
                elseif kR==Kz+1
                    ρL=st.ρ[a,b,N,kx,ky,Kz];pL=st.p[a,b,N,kx,ky,Kz];vxL=st.vx[a,b,N,kx,ky,Kz];vyL=st.vy[a,b,N,kx,ky,Kz];vzL=st.vz[a,b,N,kx,ky,Kz]
                    ρR=atm.ρ_atm;pR=atm.p_atm;vxR=0.0;vyR=0.0;vzR=0.0
                    ρeL=eng.ρeq[a,b,N,kx,ky,Kz];peL=eng.peq[a,b,N,kx,ky,Kz];ρeR=atm.ρ_atm;peR=atm.p_atm
                    elf=eng.elam[a,b,N,kx,ky,Kz];nxf=eng.nx[a,b,N,kx,ky,Kz];nyf=eng.ny[a,b,N,kx,ky,Kz];nzf=eng.nz[a,b,N,kx,ky,Kz]
                    Af=eng.α[a,b,N,kx,ky,Kz]*eng.sqrtγ[a,b,N,kx,ky,Kz]
                else
                    ρL=st.ρ[a,b,N,kx,ky,kL];pL=st.p[a,b,N,kx,ky,kL];vxL=st.vx[a,b,N,kx,ky,kL];vyL=st.vy[a,b,N,kx,ky,kL];vzL=st.vz[a,b,N,kx,ky,kL]
                    ρR=st.ρ[a,b,1,kx,ky,kR];pR=st.p[a,b,1,kx,ky,kR];vxR=st.vx[a,b,1,kx,ky,kR];vyR=st.vy[a,b,1,kx,ky,kR];vzR=st.vz[a,b,1,kx,ky,kR]
                    ρeL=eng.ρeq[a,b,N,kx,ky,kL];peL=eng.peq[a,b,N,kx,ky,kL];ρeR=eng.ρeq[a,b,1,kx,ky,kR];peR=eng.peq[a,b,1,kx,ky,kR]
                    elf=0.5*(eng.elam[a,b,N,kx,ky,kL]+eng.elam[a,b,1,kx,ky,kR])
                    nxf=0.5*(eng.nx[a,b,N,kx,ky,kL]+eng.nx[a,b,1,kx,ky,kR]); nyf=0.5*(eng.ny[a,b,N,kx,ky,kL]+eng.ny[a,b,1,kx,ky,kR])
                    nzf=0.5*(eng.nz[a,b,N,kx,ky,kL]+eng.nz[a,b,1,kx,ky,kR])
                    Af=0.5*(eng.α[a,b,N,kx,ky,kL]*eng.sqrtγ[a,b,N,kx,ky,kL]+eng.α[a,b,1,kx,ky,kR]*eng.sqrtγ[a,b,1,kx,ky,kR])
                end
                UL,FL,wL,εL=_faceflux(eos,ρL,pL,vxL,vyL,vzL,elf,nxf,nyf,nzf,3)
                UR,FR,wR,εR=_faceflux(eos,ρR,pR,vxR,vyR,vzR,elf,nxf,nyf,nzf,3)
                amax=max(abs.(_axisspeeds(eos_cs2(eos,εL),wL))...,abs.(_axisspeeds(eos_cs2(eos,εR),wR))...)
                if eng.wb
                    ULe,_,_,εeL=_faceflux(eos,ρeL,peL,0.0,0.0,0.0,elf,nxf,nyf,nzf,3)
                    URe,_,_,εeR=_faceflux(eos,ρeR,peR,0.0,0.0,0.0,elf,nxf,nyf,nzf,3)
                    ae=max(sqrt(clamp(eos_cs2(eos,εeL),0.0,1.0)),sqrt(clamp(eos_cs2(eos,εeR),0.0,1.0)))
                    F̂=_rus5wb(UL,FL,UR,FR,amax,Af,ULe,URe,ae)
                else
                    F̂=_rus5(UL,FL,UR,FR,amax,Af)
                end
                FLp=(Af.*FL); FRp=(Af.*FR)
                if kL≥1
                    wend=w[N]; sg=eng.sqrtγ[a,b,N,kx,ky,kL]; f=1.0/(J*wend*sg)
                    rD[a,b,N,kx,ky,kL]-=(F̂[1]-FLp[1])*f; rSx[a,b,N,kx,ky,kL]-=(F̂[2]-FLp[2])*f
                    rSy[a,b,N,kx,ky,kL]-=(F̂[3]-FLp[3])*f; rSz[a,b,N,kx,ky,kL]-=(F̂[4]-FLp[4])*f; rτ[a,b,N,kx,ky,kL]-=(F̂[5]-FLp[5])*f
                end
                if kR≤Kz
                    w1=w[1]; sg=eng.sqrtγ[a,b,1,kx,ky,kR]; f=1.0/(J*w1*sg)
                    rD[a,b,1,kx,ky,kR]+=(F̂[1]-FRp[1])*f; rSx[a,b,1,kx,ky,kR]+=(F̂[2]-FRp[2])*f
                    rSy[a,b,1,kx,ky,kR]+=(F̂[3]-FRp[3])*f; rSz[a,b,1,kx,ky,kR]+=(F̂[4]-FRp[4])*f; rτ[a,b,1,kx,ky,kR]+=(F̂[5]-FRp[5])*f
                end
            end
        end
    end

    # geometric/gravity source (Cartesian projection of the spherical source)
    @inbounds Threads.@threads for kz in 1:Kz
        for ky in 1:Ky, kx in 1:Kx, c in 1:N, b in 1:N, a in 1:N
            α=eng.α[a,b,c,kx,ky,kz]; Φp=eng.Φp[a,b,c,kx,ky,kz]
            nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
            ρ=st.ρ[a,b,c,kx,ky,kz]; p=st.p[a,b,c,kx,ky,kz]; ε=st.ε[a,b,c,kx,ky,kz]
            vux=st.vx[a,b,c,kx,ky,kz]; vuy=st.vy[a,b,c,kx,ky,kz]; vuz=st.vz[a,b,c,kx,ky,kz]
            el=eng.elam[a,b,c,kx,ky,kz]; q=(el-1)
            nv=nxv*vux+nyv*vuy+nzv*vuz; vlx=vux+q*nxv*nv; vly=vuy+q*nyv*nv; vlz=vuz+q*nzv*nv
            v2=clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
            vn=nxv*vlx+nyv*vly+nzv*vlz
            gforce=-α*(ε+p)*W^2*Φp
            rSx[a,b,c,kx,ky,kz]+=gforce*nxv; rSy[a,b,c,kx,ky,kz]+=gforce*nyv; rSz[a,b,c,kx,ky,kz]+=gforce*nzv
            rτ[a,b,c,kx,ky,kz]+= -α*(ε+p)*W^2*vn*Φp
        end
    end
    return nothing
end

function _rhs!(rD,rSx,rSy,rSz,rτ, st::DGCart3DState, eng::DGCart3DEngine)
    _raw_rhs!(rD,rSx,rSy,rSz,rτ, st, eng)
    @inbounds for I in eachindex(rD)
        rD[I]-=eng.Seq_D[I]; rSx[I]-=eng.Seq_Sx[I]; rSy[I]-=eng.Seq_Sy[I]
        rSz[I]-=eng.Seq_Sz[I]; rτ[I]-=eng.Seq_τ[I]
    end
end

# Zhang–Shu positivity proxy: q = (τ+D) − √(D²+|S|²) > 0 guards p>0.
@inline _q3d(D,Sx,Sy,Sz,τ,elam,nx,ny,nz) = begin
    nS=nx*Sx+ny*Sy+nz*Sz; S2=Sx^2+Sy^2+Sz^2+(1/elam-1)*nS^2
    (τ+D) - sqrt(D^2 + max(S2,0.0))
end

function _limit!(st::DGCart3DState, eng::DGCart3DEngine)
    eng.wb ? _limit_wb!(st, eng) : _limit_plain!(st, eng)
end

# ---- EQUILIBRIUM-PRESERVING Zhang–Shu (wb=true) ----------------------------------------
# Standard Zhang–Shu scales every field toward its CELL MEAN: U ← Ū + θ(U−Ū). On a static star
# the equilibrium's own intra-element variation (D spans orders of magnitude across a surface
# element) is then flattened whenever θ<1, mass is pushed from the stellar edge into the
# atmosphere nodes, and the star evaporates from the surface inward (measured: the UNSEEDED star
# thrown into a ±3% radial oscillation with 0.6c surface velocities within 50 M⊙; ρ_c −10% in
# 4.6 f-mode periods once seeded). Here the scaling acts about a REFERENCE R that (i) equals
# U_eq when the deviation δU=U−U_eq vanishes, so U_eq is a fixed point, and (ii) has the same
# cell mean as U, so the update stays conservative:
#     R_D = D_eq + δD̄·w_D,   R_S = S_eq + δS̄·w_S,   R_τ = τ_eq + ΔK + (δτ̄ − ΔK̄)·w_τ,   U ← R + θ(U−R),
# with mean-one weights w_D = D_eq/D̄_eq, w_τ = τ_eq/τ̄_eq and w_S ∝ (D_eq − D_floor): the mean
# deviation is shared out in proportion to the background, so the atmosphere nodes — where
# D_eq sits on the floor and τ_eq is 1e-14 of the interior — receive essentially none of it.
# ΔK = K(R_D,R_S) − K(D_eq,S_eq), with K(D,S) = √(D²+|S|²) − D the cold kinetic energy implied by
# a momentum S at density D, is the energy the reference momentum MUST carry: without it (the
# first version of this limiter) R_τ ≈ τ_eq at the outermost stellar nodes while R_S ≈ D v̄, and
# since τ_eq ∝ ρ² vanishes at the surface faster than ½Dv̄², the reference had negative pressure
# there for a 1% velocity seed and every surface element fell back to the flattening limiter
# (measured: the three surface elements of the K=6 octant grid on every stage of every step,
# a residual ρ_c drift of −0.7% in 5 periods). With ΔK the pressure of the reference is
# q(R) = τ_eq (1 + (δτ̄ − ΔK̄)/τ̄_eq): feasible unless the element has lost its entire equilibrium
# thermal energy, which a physical state cannot do.
# (A UNIFORM shift R = U_eq + δŪ was tried even earlier: it swamps the atmosphere nodes, R
# becomes infeasible in every surface element on the negative side of the perturbation, and
# the code fell back to the flattening limiter — ρ_c −3.6% at t=600 instead of −10%.)
# θ is chosen so D ≥ D_fl = ½√γρ_atm (strictly BELOW the equilibrium atmosphere, so round-off
# cannot trip it) and the pressure proxy q ≥ q_ε at every node. q is concave in U, so for a node
# with q(R) > q_ε > q(U) the chord root underestimates the true root and q(R+θ(U−R)) ≥ q_ε is
# guaranteed. Only if R itself is infeasible (an element that lost more than half its mass or
# all its thermal energy — in practice ejecta arriving in pure-atmosphere elements with a
# negative τ error) does the plain mean-based limiter run for that element.
struct _WBBufs
    δD::Array{Float64,3}; δSx::Array{Float64,3}; δSy::Array{Float64,3}; δSz::Array{Float64,3}; δτ::Array{Float64,3}
    RD::Array{Float64,3}; RSx::Array{Float64,3}; RSy::Array{Float64,3}; RSz::Array{Float64,3}; Rτ::Array{Float64,3}
    wS::Array{Float64,3}; ΔK::Array{Float64,3}
end
_WBBufs(N::Int) = _WBBufs([zeros(N,N,N) for _ in 1:12]...)

# One element. Returns 0 (state untouched), 1 (scaled about the equilibrium-preserving
# reference), 2 (reference infeasible → plain mean-based limiter). apply=false only classifies.
function _limit_wb_elem!(st::DGCart3DState, eng::DGCart3DEngine, kx::Int, ky::Int, kz::Int, B::_WBBufs; apply::Bool=true)
    bs=eng.bs; w=bs.w; N=bs.N; ρa=eng.atm.ρ_atm
    δD=B.δD; δSx=B.δSx; δSy=B.δSy; δSz=B.δSz; δτ=B.δτ; RD=B.RD; RSx=B.RSx; RSy=B.RSy; RSz=B.RSz; Rτ=B.Rτ; wS=B.wS; ΔK=B.ΔK
    @inbounds begin
        for c in 1:N,b in 1:N,a in 1:N
            δD[a,b,c] =st.D[a,b,c,kx,ky,kz] -eng.Deq[a,b,c,kx,ky,kz]
            δSx[a,b,c]=st.Sx[a,b,c,kx,ky,kz]-eng.Sxeq[a,b,c,kx,ky,kz]
            δSy[a,b,c]=st.Sy[a,b,c,kx,ky,kz]-eng.Syeq[a,b,c,kx,ky,kz]
            δSz[a,b,c]=st.Sz[a,b,c,kx,ky,kz]-eng.Szeq[a,b,c,kx,ky,kz]
            δτ[a,b,c] =st.τ[a,b,c,kx,ky,kz] -eng.τeq[a,b,c,kx,ky,kz]
            wS[a,b,c] =eng.Deq[a,b,c,kx,ky,kz]-eng.sqrtγ[a,b,c,kx,ky,kz]*ρa      # ≥0; 0 in the atmosphere
        end
        mD=_cellavg3(δD,w); mSx=_cellavg3(δSx,w); mSy=_cellavg3(δSy,w); mSz=_cellavg3(δSz,w); mτ=_cellavg3(δτ,w)
        Deqm=_cellavg3(view(eng.Deq,:,:,:,kx,ky,kz),w); τeqm=_cellavg3(view(eng.τeq,:,:,:,kx,ky,kz),w)
        wSm=_cellavg3(wS,w)
        for c in 1:N,b in 1:N,a in 1:N
            wD=eng.Deq[a,b,c,kx,ky,kz]/Deqm
            wsS = wSm>0 ? wS[a,b,c]/wSm : wD                 # pure-atmosphere element: share like D
            RD[a,b,c] =eng.Deq[a,b,c,kx,ky,kz]+mD*wD
            RSx[a,b,c]=eng.Sxeq[a,b,c,kx,ky,kz]+mSx*wsS
            RSy[a,b,c]=eng.Syeq[a,b,c,kx,ky,kz]+mSy*wsS
            RSz[a,b,c]=eng.Szeq[a,b,c,kx,ky,kz]+mSz*wsS
            el=eng.elam[a,b,c,kx,ky,kz]; nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
            # kinetic energy of the reference momentum, relative to the equilibrium's: K = −q(D,S,τ=0)
            ΔK[a,b,c]=_q3d(eng.Deq[a,b,c,kx,ky,kz],eng.Sxeq[a,b,c,kx,ky,kz],eng.Syeq[a,b,c,kx,ky,kz],eng.Szeq[a,b,c,kx,ky,kz],0.0,el,nxv,nyv,nzv) -
                      _q3d(RD[a,b,c],RSx[a,b,c],RSy[a,b,c],RSz[a,b,c],0.0,el,nxv,nyv,nzv)
        end
        ΔKm=_cellavg3(ΔK,w)
        for c in 1:N,b in 1:N,a in 1:N
            wτ=eng.τeq[a,b,c,kx,ky,kz]/τeqm
            Rτ[a,b,c]=eng.τeq[a,b,c,kx,ky,kz]+ΔK[a,b,c]+(mτ-ΔKm)*wτ
        end
        # ---- step 1: D ≥ D_fl = ½√γρ_atm at every node
        θ=1.0
        for c in 1:N,b in 1:N,a in 1:N
            Dfl=0.5*eng.sqrtγ[a,b,c,kx,ky,kz]*ρa
            Dn=st.D[a,b,c,kx,ky,kz]
            if Dn<Dfl
                RD[a,b,c]>Dfl || return 2
                θ=min(θ, clamp((RD[a,b,c]-Dfl)/(RD[a,b,c]-Dn+1e-300),0.0,1.0))
            end
        end
        if θ<1.0 && apply
            for c in 1:N,b in 1:N,a in 1:N
                st.D[a,b,c,kx,ky,kz]=RD[a,b,c]+θ*(st.D[a,b,c,kx,ky,kz]-RD[a,b,c])
            end
        end
        # ---- step 2: pressure proxy q ≥ q_ε at every node, about the same reference
        θq=1.0
        for c in 1:N,b in 1:N,a in 1:N
            el=eng.elam[a,b,c,kx,ky,kz]; nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
            qR=_q3d(RD[a,b,c],RSx[a,b,c],RSy[a,b,c],RSz[a,b,c],Rτ[a,b,c], el,nxv,nyv,nzv)
            qR>0.0 || return 2
            Dn = (θ<1.0 && !apply) ? RD[a,b,c]+θ*(st.D[a,b,c,kx,ky,kz]-RD[a,b,c]) : st.D[a,b,c,kx,ky,kz]
            qn=_q3d(Dn,st.Sx[a,b,c,kx,ky,kz],st.Sy[a,b,c,kx,ky,kz],st.Sz[a,b,c,kx,ky,kz],st.τ[a,b,c,kx,ky,kz], el,nxv,nyv,nzv)
            qε=1e-12*qR
            if qn<qε
                θq=min(θq, clamp((qR-qε)/(qR-qn+1e-300),0.0,1.0))
            end
        end
        if θq<1.0 && apply
            for c in 1:N,b in 1:N,a in 1:N
                st.D[a,b,c,kx,ky,kz] =RD[a,b,c] +θq*(st.D[a,b,c,kx,ky,kz] -RD[a,b,c])
                st.Sx[a,b,c,kx,ky,kz]=RSx[a,b,c]+θq*(st.Sx[a,b,c,kx,ky,kz]-RSx[a,b,c])
                st.Sy[a,b,c,kx,ky,kz]=RSy[a,b,c]+θq*(st.Sy[a,b,c,kx,ky,kz]-RSy[a,b,c])
                st.Sz[a,b,c,kx,ky,kz]=RSz[a,b,c]+θq*(st.Sz[a,b,c,kx,ky,kz]-RSz[a,b,c])
                st.τ[a,b,c,kx,ky,kz] =Rτ[a,b,c] +θq*(st.τ[a,b,c,kx,ky,kz] -Rτ[a,b,c])
            end
        end
    end
    return (θ<1.0 || θq<1.0) ? 1 : 0
end

function _limit_wb!(st::DGCart3DState, eng::DGCart3DEngine)
    N=eng.bs.N; Kx=eng.Kx; Ky=eng.Ky; Kz=eng.Kz
    @inbounds Threads.@threads for kz in 1:Kz
        B=_WBBufs(N)
        for ky in 1:Ky, kx in 1:Kx
            code=_limit_wb_elem!(st,eng,kx,ky,kz,B)
            code==2 && _limit_plain_elem!(st,eng,kx,ky,kz)
        end
    end
end

"""
    dgcart3d_limiter_census(st, eng) -> (free, scaled, fallback, star_fallback, r_fallback)

Classify every element by what the equilibrium-preserving limiter WOULD do to the current
state, without modifying it: `free` (θ=1, untouched), `scaled` (scaled about the
equilibrium-preserving reference), `fallback` (reference infeasible → plain mean-based
limiter). `star_fallback` counts fallback elements that contain stellar (interior) nodes —
the ones whose flattening would erode the star — and `r_fallback` lists the fallback
elements' rms radius in units of R. For wb=false every element is reported as `fallback`.
"""
function dgcart3d_limiter_census(st::DGCart3DState, eng::DGCart3DEngine)
    N=eng.bs.N; free=0; scaled=0; fallback=0; star_fallback=0; rf=Float64[]
    eng.wb || return (free=0, scaled=0, fallback=eng.Kx*eng.Ky*eng.Kz, star_fallback=count(any, (view(eng.interior,:,:,:,kx,ky,kz) for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx)), r_fallback=rf)
    B=_WBBufs(N)
    for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx
        code=_limit_wb_elem!(st,eng,kx,ky,kz,B; apply=false)
        if code==0; free+=1
        elseif code==1; scaled+=1
        else
            fallback+=1
            any(view(eng.interior,:,:,:,kx,ky,kz)) && (star_fallback+=1)
            push!(rf, sqrt(sum(abs2, view(eng.r,:,:,:,kx,ky,kz))/N^3)/eng.R)
        end
    end
    return (free=free, scaled=scaled, fallback=fallback, star_fallback=star_fallback, r_fallback=rf)
end

# ---- the ORIGINAL mean-based Zhang–Shu, one element (fallback, and the wb=false path) ----
function _limit_plain_elem!(st::DGCart3DState, eng::DGCart3DEngine, kx::Int, ky::Int, kz::Int)
    bs=eng.bs; w=bs.w; N=bs.N
    @inbounds begin
        D̄=_cellavg3(view(st.D,:,:,:,kx,ky,kz),w)
        S̄x=_cellavg3(view(st.Sx,:,:,:,kx,ky,kz),w); S̄y=_cellavg3(view(st.Sy,:,:,:,kx,ky,kz),w)
        S̄z=_cellavg3(view(st.Sz,:,:,:,kx,ky,kz),w); τ̄=_cellavg3(view(st.τ,:,:,:,kx,ky,kz),w)
        Dε=eng.atm.ρ_atm*_cellavg3(view(eng.sqrtγ,:,:,:,kx,ky,kz),w)
        D̄=max(D̄,Dε)
        Dmin=Inf; for c in 1:N,b in 1:N,a in 1:N; Dmin=min(Dmin,st.D[a,b,c,kx,ky,kz]); end
        if Dmin<Dε
            θ=clamp((D̄-Dε)/(D̄-Dmin+1e-300),0.0,1.0)
            for c in 1:N,b in 1:N,a in 1:N; st.D[a,b,c,kx,ky,kz]=θ*(st.D[a,b,c,kx,ky,kz]-D̄)+D̄; end
        end
        elf=eng.elam[1,1,1,kx,ky,kz]; nxf=eng.nx[1,1,1,kx,ky,kz]; nyf=eng.ny[1,1,1,kx,ky,kz]; nzf=eng.nz[1,1,1,kx,ky,kz]
        q̄=_q3d(D̄,S̄x,S̄y,S̄z,τ̄,elf,nxf,nyf,nzf)
        if q̄ ≤ 0.0
            for c in 1:N,b in 1:N,a in 1:N
                st.D[a,b,c,kx,ky,kz]=max(D̄,Dε); st.Sx[a,b,c,kx,ky,kz]=S̄x; st.Sy[a,b,c,kx,ky,kz]=S̄y
                st.Sz[a,b,c,kx,ky,kz]=S̄z; st.τ[a,b,c,kx,ky,kz]=τ̄
            end
            return
        end
        qε=1e-12*q̄; θq=1.0
        for c in 1:N,b in 1:N,a in 1:N
            qn=_q3d(st.D[a,b,c,kx,ky,kz],st.Sx[a,b,c,kx,ky,kz],st.Sy[a,b,c,kx,ky,kz],st.Sz[a,b,c,kx,ky,kz],
                    st.τ[a,b,c,kx,ky,kz],eng.elam[a,b,c,kx,ky,kz],eng.nx[a,b,c,kx,ky,kz],
                    eng.ny[a,b,c,kx,ky,kz],eng.nz[a,b,c,kx,ky,kz])
            if qn<qε
                θq=min(θq, clamp((q̄-qε)/(q̄-qn+1e-300),0.0,1.0))
            end
        end
        if θq<1.0
            for c in 1:N,b in 1:N,a in 1:N
                st.D[a,b,c,kx,ky,kz]=θq*(st.D[a,b,c,kx,ky,kz]-D̄)+D̄
                st.Sx[a,b,c,kx,ky,kz]=θq*(st.Sx[a,b,c,kx,ky,kz]-S̄x)+S̄x
                st.Sy[a,b,c,kx,ky,kz]=θq*(st.Sy[a,b,c,kx,ky,kz]-S̄y)+S̄y
                st.Sz[a,b,c,kx,ky,kz]=θq*(st.Sz[a,b,c,kx,ky,kz]-S̄z)+S̄z
                st.τ[a,b,c,kx,ky,kz]=θq*(st.τ[a,b,c,kx,ky,kz]-τ̄)+τ̄
            end
        end
    end
end

function _limit_plain!(st::DGCart3DState, eng::DGCart3DEngine)
    @inbounds Threads.@threads for kz in 1:eng.Kz
        for ky in 1:eng.Ky, kx in 1:eng.Kx
            _limit_plain_elem!(st,eng,kx,ky,kz)
        end
    end
end

@inline function _maxspeed(st::DGCart3DState, eng::DGCart3DEngine)
    eos=eng.eos; N=eng.bs.N; a=0.0
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, aa in 1:N
        cs=sqrt(clamp(eos_cs2(eos,st.ε[aa,b,c,kx,ky,kz]),0.0,1.0))
        el=eng.elam[aa,b,c,kx,ky,kz]; nxv=eng.nx[aa,b,c,kx,ky,kz]; nyv=eng.ny[aa,b,c,kx,ky,kz]; nzv=eng.nz[aa,b,c,kx,ky,kz]
        vux=st.vx[aa,b,c,kx,ky,kz]; vuy=st.vy[aa,b,c,kx,ky,kz]; vuz=st.vz[aa,b,c,kx,ky,kz]
        q=el-1; nv=nxv*vux+nyv*vuy+nzv*vuz
        vlx=vux+q*nxv*nv; vly=vuy+q*nyv*nv; vlz=vuz+q*nzv*nv
        vmag=sqrt(clamp(vux*vlx+vuy*vly+vuz*vlz,0.0,1.0))
        a=max(a,(vmag+cs)/(1+vmag*cs)*sqrt(el))
    end
    return a
end

"""
    evolve_dgcart3d!(st, eng; tmax, sample_dt) -> (ts, q2, ρc)

SSP-RK3 evolution of the full-3D DG star; records the ℓ=2 (m=0) quadrupole and
the central density.
"""
function evolve_dgcart3d!(st::DGCart3DState, eng::DGCart3DEngine; tmax::Float64,
        sample_dt::Float64=4.0, cfl::Float64=-1.0, verbose::Bool=false,
        record_m2::Bool=false, limiter::Bool=true)
    bs=eng.bs; p=bs.p
    cfl=cfl>0 ? cfl : eng.cfl; cfl_dg=cfl/(2p+1)
    dims=size(st.D)
    rD=zeros(dims);rSx=zeros(dims);rSy=zeros(dims);rSz=zeros(dims);rτ=zeros(dims)
    D0=zeros(dims);Sx0=zeros(dims);Sy0=zeros(dims);Sz0=zeros(dims);τ0=zeros(dims)
    ts=Float64[];q2=Float64[];ρch=Float64[]; qm2=Float64[]
    ρc0=dgcart3d_central_density(st,eng)
    t=0.0; last=-1e30; ns=0
    _update_prims!(st,eng)
    while t<tmax && ns<5_000_000
        amax=max(_maxspeed(st,eng),1e-3); dt=cfl_dg*eng.Δ/amax; dt=min(dt,tmax-t)
        copyto!(D0,st.D);copyto!(Sx0,st.Sx);copyto!(Sy0,st.Sy);copyto!(Sz0,st.Sz);copyto!(τ0,st.τ)
        _rhs!(rD,rSx,rSy,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=D0[I]+dt*rD[I];st.Sx[I]=Sx0[I]+dt*rSx[I];st.Sy[I]=Sy0[I]+dt*rSy[I]
            st.Sz[I]=Sz0[I]+dt*rSz[I];st.τ[I]=τ0[I]+dt*rτ[I]
        end
        limiter && _limit!(st,eng);_update_prims!(st,eng)
        _rhs!(rD,rSx,rSy,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=0.75*D0[I]+0.25*(st.D[I]+dt*rD[I]);st.Sx[I]=0.75*Sx0[I]+0.25*(st.Sx[I]+dt*rSx[I])
            st.Sy[I]=0.75*Sy0[I]+0.25*(st.Sy[I]+dt*rSy[I]);st.Sz[I]=0.75*Sz0[I]+0.25*(st.Sz[I]+dt*rSz[I])
            st.τ[I]=0.75*τ0[I]+0.25*(st.τ[I]+dt*rτ[I])
        end
        limiter && _limit!(st,eng);_update_prims!(st,eng)
        _rhs!(rD,rSx,rSy,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=(1/3)*D0[I]+(2/3)*(st.D[I]+dt*rD[I]);st.Sx[I]=(1/3)*Sx0[I]+(2/3)*(st.Sx[I]+dt*rSx[I])
            st.Sy[I]=(1/3)*Sy0[I]+(2/3)*(st.Sy[I]+dt*rSy[I]);st.Sz[I]=(1/3)*Sz0[I]+(2/3)*(st.Sz[I]+dt*rSz[I])
            st.τ[I]=(1/3)*τ0[I]+(2/3)*(st.τ[I]+dt*rτ[I])
        end
        limiter && _limit!(st,eng);_update_prims!(st,eng)
        t+=dt; ns+=1
        if t-last≥sample_dt
            push!(ts,t);push!(q2,dgcart3d_quadrupole(st,eng));push!(ρch,dgcart3d_central_density(st,eng))
            record_m2 && push!(qm2,dgcart3d_quadrupole_m2(st,eng))
            last=t
            verbose && println("  t=",round(t,digits=1)," ρc/ρc0=",round(ρch[end]/ρc0,digits=6))
        end
        if !isempty(ρch) && (!isfinite(ρch[end]) || ρch[end]>100*ρc0)
            return record_m2 ? (ts,q2,ρch,qm2) : (ts,q2,ρch)
        end
    end
    return record_m2 ? (ts,q2,ρch,qm2) : (ts,q2,ρch)
end

dgcart3d_central_density(st::DGCart3DState, eng::DGCart3DEngine) =
    _cellavg3(view(st.ρ,:,:,:,1,1,1), eng.bs.w)

"""min/max of (ρ,p) over all nodes — for positivity / overshoot checks."""
function dgcart3d_prim_minmax(st::DGCart3DState)
    return (minimum(st.ρ), maximum(st.ρ), minimum(st.p), maximum(st.p),
            minimum(st.vx), maximum(st.vx))
end

"""seed ℓ=2,m=0 (Y20) radial velocity: v^r=A(r/R)(3cos²θ−1), here ∝ (2n_z²−n_x²−n_y²)."""
function seed_dgcart3d_l2!(st::DGCart3DState, eng::DGCart3DEngine; A::Float64=1e-4)
    eos=eng.eos; N=eng.bs.N
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        eng.interior[a,b,c,kx,ky,kz] || continue
        rr=eng.r[a,b,c,kx,ky,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
        Y20=2*nzv^2-nxv^2-nyv^2
        vr=A*(rr/eng.R)*Y20
        el2=sqrt(eng.elam[a,b,c,kx,ky,kz])
        vux=vr/el2*nxv; vuy=vr/el2*nyv; vuz=vr/el2*nzv
        st.vx[a,b,c,kx,ky,kz]=vux; st.vy[a,b,c,kx,ky,kz]=vuy; st.vz[a,b,c,kx,ky,kz]=vuz
        D,Sx,Sy,Sz,τ=_p2c(eos,st.ρ[a,b,c,kx,ky,kz],st.p[a,b,c,kx,ky,kz],vux,vuy,vuz,
            eng.sqrtγ[a,b,c,kx,ky,kz],eng.elam[a,b,c,kx,ky,kz],nxv,nyv,nzv)
        st.D[a,b,c,kx,ky,kz]=D; st.Sx[a,b,c,kx,ky,kz]=Sx; st.Sy[a,b,c,kx,ky,kz]=Sy
        st.Sz[a,b,c,kx,ky,kz]=Sz; st.τ[a,b,c,kx,ky,kz]=τ
    end
end

"""
seed a genuinely NON-AXISYMMETRIC ℓ=2,m=2 (Y22 real part) perturbation:
Y22 ∝ sin²θ cos2φ ∝ (n_x²−n_y²). Radial velocity v^r=A(r/R)(n_x²−n_y²). Even
under x→−x and y→−y, so the octant box captures one lobe consistently. This is
the capability the 2D axisymmetric engine CANNOT represent (m=0 only).
"""
function seed_dgcart3d_Y22!(st::DGCart3DState, eng::DGCart3DEngine; A::Float64=1e-4)
    eos=eng.eos; N=eng.bs.N
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        eng.interior[a,b,c,kx,ky,kz] || continue
        rr=eng.r[a,b,c,kx,ky,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
        Y22=nxv^2-nyv^2
        vr=A*(rr/eng.R)*Y22
        el2=sqrt(eng.elam[a,b,c,kx,ky,kz])
        vux=vr/el2*nxv; vuy=vr/el2*nyv; vuz=vr/el2*nzv
        st.vx[a,b,c,kx,ky,kz]=vux; st.vy[a,b,c,kx,ky,kz]=vuy; st.vz[a,b,c,kx,ky,kz]=vuz
        D,Sx,Sy,Sz,τ=_p2c(eos,st.ρ[a,b,c,kx,ky,kz],st.p[a,b,c,kx,ky,kz],vux,vuy,vuz,
            eng.sqrtγ[a,b,c,kx,ky,kz],eng.elam[a,b,c,kx,ky,kz],nxv,nyv,nzv)
        st.D[a,b,c,kx,ky,kz]=D; st.Sx[a,b,c,kx,ky,kz]=Sx; st.Sy[a,b,c,kx,ky,kz]=Sy
        st.Sz[a,b,c,kx,ky,kz]=Sz; st.τ[a,b,c,kx,ky,kz]=τ
    end
end

"""
seed ℓ=2,m=1 (Y21 real part): Y21 ∝ sinθcosθ cosφ ∝ n_x n_z. Even under x→−x AND
z→−z so consistent in the octant. Another non-axisymmetric mode.
"""
function seed_dgcart3d_Y21!(st::DGCart3DState, eng::DGCart3DEngine; A::Float64=1e-4)
    eos=eng.eos; N=eng.bs.N
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        eng.interior[a,b,c,kx,ky,kz] || continue
        rr=eng.r[a,b,c,kx,ky,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
        Y21=nxv*nzv
        vr=A*(rr/eng.R)*Y21
        el2=sqrt(eng.elam[a,b,c,kx,ky,kz])
        vux=vr/el2*nxv; vuy=vr/el2*nyv; vuz=vr/el2*nzv
        st.vx[a,b,c,kx,ky,kz]=vux; st.vy[a,b,c,kx,ky,kz]=vuy; st.vz[a,b,c,kx,ky,kz]=vuz
        D,Sx,Sy,Sz,τ=_p2c(eos,st.ρ[a,b,c,kx,ky,kz],st.p[a,b,c,kx,ky,kz],vux,vuy,vuz,
            eng.sqrtγ[a,b,c,kx,ky,kz],eng.elam[a,b,c,kx,ky,kz],nxv,nyv,nzv)
        st.D[a,b,c,kx,ky,kz]=D; st.Sx[a,b,c,kx,ky,kz]=Sx; st.Sy[a,b,c,kx,ky,kz]=Sy
        st.Sz[a,b,c,kx,ky,kz]=Sz; st.τ[a,b,c,kx,ky,kz]=τ
    end
end

"""√γ-weighted ℓ=2,m=0 quadrupole moment ∫ρ(2n_z²−n_x²−n_y²)√γ dV over interior."""
function dgcart3d_quadrupole(st::DGCart3DState, eng::DGCart3DEngine)
    N=eng.bs.N; w=eng.bs.w; acc=0.0
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        eng.interior[a,b,c,kx,ky,kz] || continue
        rr=eng.r[a,b,c,kx,ky,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]; nzv=eng.nz[a,b,c,kx,ky,kz]
        acc+=st.ρ[a,b,c,kx,ky,kz]*(2nzv^2-nxv^2-nyv^2)*eng.sqrtγ[a,b,c,kx,ky,kz]*w[a]*w[b]*w[c]
    end
    acc*eng.J^3
end

"""√γ-weighted ℓ=2,m=2 quadrupole moment ∫ρ(n_x²−n_y²)√γ dV over interior — the
NON-AXISYMMETRIC moment. Zero for any axisymmetric configuration; nonzero only
when a genuine m=2 deformation is present."""
function dgcart3d_quadrupole_m2(st::DGCart3DState, eng::DGCart3DEngine)
    N=eng.bs.N; w=eng.bs.w; acc=0.0
    @inbounds for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        eng.interior[a,b,c,kx,ky,kz] || continue
        rr=eng.r[a,b,c,kx,ky,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,c,kx,ky,kz]; nyv=eng.ny[a,b,c,kx,ky,kz]
        acc+=st.ρ[a,b,c,kx,ky,kz]*(nxv^2-nyv^2)*eng.sqrtγ[a,b,c,kx,ky,kz]*w[a]*w[b]*w[c]
    end
    acc*eng.J^3
end

"""
    dgcart3d_shocktube_diagonal!(eos, atm, p, K; ...) -> (eng-like NamedTuple, evolve closure)

Standalone FLAT-space 3D DG shock sanity test (Martí–Müller test 1) along a grid
direction, to confirm the 3D directional sweeps + limiters preserve the
1D-validated shock behavior (positivity, no oscillation) in three spatial dims.

Builds a small flat box (γ_ij=δ_ij, no gravity) and initializes a 1D shock tube
along axis `dir` (1=x, 2=y, 3=z, or :diag for the grid diagonal). Returns the
evolved (ρ,p,v) cell-mean profiles and positivity flag. This reuses the same
_faceflux / Rusanov / SSP-RK3 / limiter machinery via a thin flat-metric engine.
"""
function dgcart3d_shocktube_diagonal!(eos::BarotropicEOS; K::Int=24, p::Int=2,
        L::Float64=1.0, tmax::Float64=0.30, cfl::Float64=0.2, dir::Symbol=:x,
        ρL=10.0, pL=13.33, ρR=1.0, pR=1e-3)
    bs=build_lgl_basis(p); N=bs.N
    Δ=L/K; J=Δ/2
    dims=(N,N,N,K,K,K)
    # flat metric
    elam=ones(dims); sqrtγ=ones(dims); α=ones(dims)
    nx=zeros(dims); ny=zeros(dims); nz=zeros(dims); Φp=zeros(dims)
    interior=trues(dims)
    x=zeros(dims);y=zeros(dims);z=zeros(dims);r=zeros(dims)
    ρ_atm=1e-10*ρR; p_atm=max(_p_of_rho(eos, ρ_atm),1e-30)
    ε_atm=energy_from_pressure(eos,p_atm)
    atm=AtmospherePars(ρ_atm,p_atm,ε_atm,5*ρ_atm,0.999)
    eng=DGCart3DEngine(bs,K,K,K,Δ,J,eos,atm,1e30,0.0,cfl,0.0,
        x,y,z,r,α,elam,sqrtγ,nx,ny,nz,Φp,interior,
        zeros(dims),zeros(dims),zeros(dims),zeros(dims),zeros(dims),
        false, zeros(dims),zeros(dims),zeros(dims),zeros(dims),zeros(dims), zeros(dims),zeros(dims))
    st=DGCart3DState((zeros(dims) for _ in 1:11)...)
    # signed coordinate s along the test direction; split at the MIDPOINT of its
    # range so equal halves are L/R. For :diag, s=(x+y+z)/√3 ranges [0, L√3], mid=L√3/2.
    smid = dir==:diag ? L*sqrt(3.0)/2 : L/2
    @inbounds for kz in 1:K, ky in 1:K, kx in 1:K, c in 1:N, b in 1:N, a in 1:N
        X=(kx-0.5)*Δ+J*bs.ξ[a]; Y=(ky-0.5)*Δ+J*bs.ξ[b]; Z=(kz-0.5)*Δ+J*bs.ξ[c]
        eng.x[a,b,c,kx,ky,kz]=X; eng.y[a,b,c,kx,ky,kz]=Y; eng.z[a,b,c,kx,ky,kz]=Z
        s = dir==:x ? X : (dir==:y ? Y : (dir==:z ? Z : (X+Y+Z)/sqrt(3.0)))
        useL = s < smid
        ρ0 = useL ? ρL : ρR; p0 = useL ? pL : pR
        st.ρ[a,b,c,kx,ky,kz]=ρ0; st.p[a,b,c,kx,ky,kz]=p0
        st.ε[a,b,c,kx,ky,kz]=energy_from_pressure(eos,p0)
        D,Sx,Sy,Sz,τ=_p2c(eos,ρ0,p0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0)
        st.D[a,b,c,kx,ky,kz]=D;st.Sx[a,b,c,kx,ky,kz]=Sx;st.Sy[a,b,c,kx,ky,kz]=Sy
        st.Sz[a,b,c,kx,ky,kz]=Sz;st.τ[a,b,c,kx,ky,kz]=τ
    end
    _update_prims!(st,eng)
    # NOTE: Seq_* stay zero (flat, no background subtraction) → pure conservation law.
    evolve_dgcart3d!(st,eng; tmax=tmax, sample_dt=tmax, cfl=cfl)
    mn=dgcart3d_prim_minmax(st)
    return (st=st, eng=eng, ρmin=mn[1], ρmax=mn[2], pmin=mn[3], pmax=mn[4],
            vmin=mn[5], vmax=mn[6], ρL=ρL, ρR=ρR)
end

end # module DGCart3D
