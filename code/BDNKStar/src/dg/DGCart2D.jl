#=
    DGCart2D — STAGE 3: 2D Cartesian (axisymmetric x–z plane, m=0) tensor-product
    nodal RKDG on the frozen TOV (Cowling) metric, with the NON-CONFORMING stellar
    surface + atmosphere + positivity/TVB limiter. This is the actual target test:
    does DG cure the staircase instability that 2nd-order FV could not?

    We work in the (x,z) meridional plane (the ℓ=2,m=0 mode is axisymmetric), so
    the 2D Cartesian DG with reflecting symmetry on x=0,z=0 (quadrant) captures the
    ℓ=2 surface dynamics at a fraction of a full 3D cost — feasible at ≤6 threads.
    The metric is the spherical Cowling one projected on Cartesian axes (same as
    FVCartesian): γ_ij=δ_ij+(e^λ−1)n_in_j with n=(x,z)/r in-plane.

    Conserved (densitized, 2D in-plane momentum):
        D=√γ ρW, S_x=√γ ρhW²v_x, S_z=√γ ρhW²v_z, τ=√γ(ρhW²−p−ρW), √γ=e^{λ/2}.
    Tensor-product LGL nodes per element; dimension-by-dimension Rusanov flux;
    SSP-RK3; per-element TVB troubled-cell + Zhang–Shu positivity limiter. The
    static background DG RHS is stored & subtracted (well-balanced).

    The KEY diagnostic is the RESOLUTION TREND of the ℓ=2-perturbed central-density
    drift: FV got WORSE with refinement; DG should improve or at least not worsen.
=#
module DGCart2D

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, rho_from_p, eos_cs2, _solve_p
using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis, tvb_minmod
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram

export DGCart2DEngine, DGCart2DState, setup_dgcart2d, evolve_dgcart2d!,
       seed_dgcart2d_l2!, dgcart2d_quadrupole, dgcart2d_central_density

@inline _lininterp(xs,ys,x)=begin
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end
@inline _p_of_rho(eos::ShumPolytrope,ρ)=eos.κ*ρ^2

# cell-average of a (N×N) nodal block with tensor LGL weights
@inline function _cellavg2(blk, w)
    s=0.0; N=length(w)
    @inbounds for b in 1:N, a in 1:N; s+=w[a]*w[b]*blk[a,b]; end
    return 0.25*s
end

struct DGCart2DEngine
    bs::LGLBasis
    Kx::Int; Kz::Int
    Δ::Float64; J::Float64
    eos::BarotropicEOS
    atm::AtmospherePars
    R::Float64; M::Float64
    cfl::Float64; M_tvb::Float64
    # per-node background: arrays (N,N,Kx,Kz)
    x::Array{Float64,4}; z::Array{Float64,4}; r::Array{Float64,4}
    α::Array{Float64,4}; elam::Array{Float64,4}; sqrtγ::Array{Float64,4}
    nx::Array{Float64,4}; nz::Array{Float64,4}; Φp::Array{Float64,4}
    interior::BitArray{4}
    Seq_D::Array{Float64,4}; Seq_Sx::Array{Float64,4}; Seq_Sz::Array{Float64,4}; Seq_τ::Array{Float64,4}
end

mutable struct DGCart2DState
    D::Array{Float64,4}; Sx::Array{Float64,4}; Sz::Array{Float64,4}; τ::Array{Float64,4}
    ρ::Array{Float64,4}; p::Array{Float64,4}; ε::Array{Float64,4}
    vx::Array{Float64,4}; vz::Array{Float64,4}
end

# cons2prim for the 2D in-plane metric (n in x–z plane, n_y=0)
function _c2p(eos, atm, D,Sx,Sz,τ, sqrtγ,elam,nx,nz)
    if sqrtγ≤0; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0) end
    D̂=D/sqrtγ; Ŝx=Sx/sqrtγ; Ŝz=Sz/sqrtγ; τ̂=τ/sqrtγ
    if D̂≤atm.ρ_cut; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0) end
    nS=nx*Ŝx+nz*Ŝz
    S2=Ŝx^2+Ŝz^2+(1/elam-1)*nS^2; S2=max(S2,0.0); Ŝ=sqrt(S2)
    pmax=10.0*(abs(τ̂)+D̂+atm.p_atm)+1e-30
    p,ok=_solve_p(eos,D̂,Ŝ,τ̂,1e-30,pmax)
    if !ok||!isfinite(p); return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0) end
    p=max(p,atm.p_atm); E=τ̂+p+D̂
    v2=clamp((Ŝ/E)^2,0.0,atm.vmax^2); W=1.0/sqrt(1.0-v2); ρ=D̂/W
    if ρ≤atm.ρ_cut; return (atm.ρ_atm,atm.p_atm,atm.ε_atm,0.0,0.0) end
    ε=E/W^2-p; invE=E>0 ? 1.0/E : 0.0
    vlx=Ŝx*invE; vlz=Ŝz*invE; nv=nx*vlx+nz*vlz; q=(1/elam-1)
    vux=vlx+q*nx*nv; vuz=vlz+q*nz*nv
    return (ρ,p,max(ε,0.0),vux,vuz)
end

function _p2c(eos, ρ,p,vux,vuz, sqrtγ,elam,nx,nz)
    ε=ρ>0 ? energy_from_pressure(eos,p) : 0.0; h=ρ>0 ? (ε+p)/ρ : 1.0
    nv=nx*vux+nz*vuz; q=(elam-1)
    vlx=vux+q*nx*nv; vlz=vuz+q*nz*nv
    v2=clamp(vux*vlx+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
    D=sqrtγ*ρ*W; fac=sqrtγ*ρ*h*W^2
    return (D, fac*vlx, fac*vlz, sqrtγ*(ρ*h*W^2-p-ρ*W))
end

# undensitized face state + flux along axis (1=x,2=z)
@inline function _faceflux(eos, ρ,p,vux,vuz, elam,nx,nz, ax)
    ε=energy_from_pressure(eos,p); h=ρ>0 ? (ε+p)/ρ : 1.0
    nv=nx*vux+nz*vuz; q=(elam-1)
    vlx=vux+q*nx*nv; vlz=vuz+q*nz*nv
    v2=clamp(vux*vlx+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
    D̂=ρ*W; fac=ρ*h*W^2; Ŝx=fac*vlx; Ŝz=fac*vlz; τ̂=fac-p-ρ*W
    va= ax==1 ? vux : vuz
    δx= ax==1 ? 1.0 : 0.0; δz= ax==2 ? 1.0 : 0.0
    FD=D̂*va; FSx=Ŝx*va+p*δx; FSz=Ŝz*va+p*δz; Fτ=(τ̂+p)*va
    na= ax==1 ? nx : nz; γaa=1+(elam-1)*na^2; w=sqrt(γaa)*va
    return (D̂,Ŝx,Ŝz,τ̂),(FD,FSx,FSz,Fτ), w, ε
end
@inline function _axisspeeds(cs2,w)
    cs=sqrt(clamp(cs2,0.0,1.0)); w2=clamp(w^2,0.0,1.0-1e-12)
    a=1.0/(1.0-w2*cs^2); disc=sqrt(max(cs^2*(1-w2)*(1-w2*cs^2-w^2*(1-cs^2)),0.0))
    return a*(w*(1-cs^2)-disc), a*(w*(1-cs^2)+disc)
end

"""
    setup_dgcart2d(eos, εc; Kx=12, Kz=12, p=2, L_fac=1.3, atm_fac=1e-7,
                   cfl=0.2, M_tvb=50.0, h_tov=2e-4) -> (engine, state)

Build the meridional (x–z) quadrant DG engine on the TOV star. `Kx,Kz` elements
of degree `p`. Effective resolution per axis ≈ Kx·p.
"""
function setup_dgcart2d(eos::BarotropicEOS, εc::Float64; Kx::Int=12, Kz::Int=12,
        p::Int=2, L_fac::Float64=1.3, atm_fac::Float64=1e-7, cfl::Float64=0.2,
        M_tvb::Float64=50.0, h_tov::Float64=2e-4)
    star=solve_tov(eos,εc;h=h_tov); R,M=star.R,star.M
    L=L_fac*R; Δ=L/Kx; J=Δ/2
    bs=build_lgl_basis(p); N=bs.N
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p
    m_of(r)= r≤R ? _lininterp(rt,mt,r) : M
    ν_of(r)= r≤R ? _lininterp(rt,νt,r) : log(1-2M/r)
    p_of(r)= r≤R ? max(_lininterp(rt,pt,r),0.0) : 0.0

    dims=(N,N,Kx,Kz)
    x=zeros(dims); z=zeros(dims); r=zeros(dims); α=zeros(dims); elam=zeros(dims)
    sqrtγ=zeros(dims); nx=zeros(dims); nz=zeros(dims); Φp=zeros(dims)
    interior=falses(dims)
    for kz in 1:Kz, kx in 1:Kx
        xc=(kx-0.5)*Δ; zc=(kz-0.5)*Δ
        for b in 1:N, a in 1:N
            X=xc+J*bs.ξ[a]; Z=zc+J*bs.ξ[b]; rr=sqrt(X^2+Z^2)
            x[a,b,kx,kz]=X; z[a,b,kx,kz]=Z; r[a,b,kx,kz]=rr
            if rr≤1e-10
                α[a,b,kx,kz]=exp(0.5*ν_of(0.0)); elam[a,b,kx,kz]=1.0; sqrtγ[a,b,kx,kz]=1.0
                nx[a,b,kx,kz]=0.0; nz[a,b,kx,kz]=0.0; Φp[a,b,kx,kz]=0.0; continue
            end
            m=m_of(rr); fm=max(1-2m/rr,1e-12)
            α[a,b,kx,kz]=exp(0.5*ν_of(rr)); elam[a,b,kx,kz]=1/fm; sqrtγ[a,b,kx,kz]=1/sqrt(fm)
            nx[a,b,kx,kz]=X/rr; nz[a,b,kx,kz]=Z/rr
            den=rr*(rr-2m); Φp[a,b,kx,kz]= den>0 ? (m+4π*rr^3*p_of(rr))/den : 0.0
            interior[a,b,kx,kz]= rr<R
        end
    end

    ρc=rho_from_p(eos,pressure(eos,εc)); ρ_atm=atm_fac*ρc; p_atm=_p_of_rho(eos,ρ_atm)
    ε_atm=energy_from_pressure(eos,p_atm)
    atm=AtmospherePars(ρ_atm,p_atm,ε_atm,5*ρ_atm,0.999)

    st=DGCart2DState((zeros(dims) for _ in 1:9)...)
    for kz in 1:Kz, kx in 1:Kx, b in 1:N, a in 1:N
        rr=r[a,b,kx,kz]
        ρ= rr≤R ? max(rho_from_p(eos,p_of(rr)),ρ_atm) : ρ_atm
        pp= rr≤R ? max(p_of(rr),p_atm) : p_atm
        st.ρ[a,b,kx,kz]=ρ; st.p[a,b,kx,kz]=pp; st.ε[a,b,kx,kz]=energy_from_pressure(eos,pp)
        D,Sx,Sz,τ=_p2c(eos,ρ,pp,0.0,0.0,sqrtγ[a,b,kx,kz],elam[a,b,kx,kz],nx[a,b,kx,kz],nz[a,b,kx,kz])
        st.D[a,b,kx,kz]=D; st.Sx[a,b,kx,kz]=Sx; st.Sz[a,b,kx,kz]=Sz; st.τ[a,b,kx,kz]=τ
    end
    eng=DGCart2DEngine(bs,Kx,Kz,Δ,J,eos,atm,R,M,cfl,M_tvb,x,z,r,α,elam,sqrtγ,nx,nz,Φp,interior,
                       zeros(dims),zeros(dims),zeros(dims),zeros(dims))
    _update_prims!(st,eng)
    _raw_rhs!(eng.Seq_D,eng.Seq_Sx,eng.Seq_Sz,eng.Seq_τ, st, eng)
    return eng, st
end

function _update_prims!(st::DGCart2DState, eng::DGCart2DEngine)
    eos=eng.eos; atm=eng.atm; N=eng.bs.N
    @inbounds Threads.@threads for kz in 1:eng.Kz
        for kx in 1:eng.Kx, b in 1:N, a in 1:N
            ρ,p,ε,vx,vz=_c2p(eos,atm,st.D[a,b,kx,kz],st.Sx[a,b,kx,kz],st.Sz[a,b,kx,kz],st.τ[a,b,kx,kz],
                eng.sqrtγ[a,b,kx,kz],eng.elam[a,b,kx,kz],eng.nx[a,b,kx,kz],eng.nz[a,b,kx,kz])
            st.ρ[a,b,kx,kz]=ρ; st.p[a,b,kx,kz]=p; st.ε[a,b,kx,kz]=ε; st.vx[a,b,kx,kz]=vx; st.vz[a,b,kx,kz]=vz
        end
    end
end

# Rusanov flux at a face between (UL,FL) and (UR,FR), area Af.
@inline function _rus4(UL,FL,UR,FR,amax,Af)
    ntuple(i->0.5*Af*(FL[i]+FR[i])-0.5*amax*Af*(UR[i]-UL[i]), 4)
end

function _raw_rhs!(rD,rSx,rSz,rτ, st::DGCart2DState, eng::DGCart2DEngine)
    bs=eng.bs; N=bs.N; Kx=eng.Kx; Kz=eng.Kz; J=eng.J; eos=eng.eos; atm=eng.atm
    Dm=bs.D; w=bs.w
    fill!(rD,0.0); fill!(rSx,0.0); fill!(rSz,0.0); fill!(rτ,0.0)
    # VOLUME term (both directions), element-local flux buffers per thread.
    @inbounds Threads.@threads for kz in 1:Kz
        # element-local flux buffers
        FxD=zeros(N,N); FxSx=zeros(N,N); FxSz=zeros(N,N); Fxτ=zeros(N,N)
        FzD=zeros(N,N); FzSx=zeros(N,N); FzSz=zeros(N,N); Fzτ=zeros(N,N)
        for kx in 1:Kx
            for b in 1:N, a in 1:N
                ρ=st.ρ[a,b,kx,kz]; p=st.p[a,b,kx,kz]; vx=st.vx[a,b,kx,kz]; vz=st.vz[a,b,kx,kz]
                el=eng.elam[a,b,kx,kz]; nxv=eng.nx[a,b,kx,kz]; nzv=eng.nz[a,b,kx,kz]
                A=eng.α[a,b,kx,kz]*eng.sqrtγ[a,b,kx,kz]
                _,Fx,_,_=_faceflux(eos,ρ,p,vx,vz,el,nxv,nzv,1)
                _,Fz,_,_=_faceflux(eos,ρ,p,vx,vz,el,nxv,nzv,2)
                FxD[a,b]=A*Fx[1]; FxSx[a,b]=A*Fx[2]; FxSz[a,b]=A*Fx[3]; Fxτ[a,b]=A*Fx[4]
                FzD[a,b]=A*Fz[1]; FzSx[a,b]=A*Fz[2]; FzSz[a,b]=A*Fz[3]; Fzτ[a,b]=A*Fz[4]
            end
            # volume divergence: d/dx along a-index, d/dz along b-index
            for b in 1:N, a in 1:N
                sD=0.0;sSx=0.0;sSz=0.0;sτ=0.0
                for c in 1:N
                    dx=Dm[a,c]; dz=Dm[b,c]
                    sD+=dx*FxD[c,b]+dz*FzD[a,c]
                    sSx+=dx*FxSx[c,b]+dz*FzSx[a,c]
                    sSz+=dx*FxSz[c,b]+dz*FzSz[a,c]
                    sτ+=dx*Fxτ[c,b]+dz*Fzτ[a,c]
                end
                sg=eng.sqrtγ[a,b,kx,kz]
                rD[a,b,kx,kz]=-sD/(J*sg); rSx[a,b,kx,kz]=-sSx/(J*sg)
                rSz[a,b,kx,kz]=-sSz/(J*sg); rτ[a,b,kx,kz]=-sτ/(J*sg)
            end
        end
    end
    # SURFACE term — X faces
    @inbounds Threads.@threads for kz in 1:Kz
        for kf in 0:Kx
            kL=kf; kR=kf+1
            for b in 1:N
                # left/right states at the face
                if kL==0
                    ρR=st.ρ[1,b,1,kz];pR=st.p[1,b,1,kz];vxR=st.vx[1,b,1,kz];vzR=st.vz[1,b,1,kz]
                    ρL=ρR;pL=pR;vxL=-vxR;vzL=vzR
                    elf=eng.elam[1,b,1,kz];nxf=eng.nx[1,b,1,kz];nzf=eng.nz[1,b,1,kz];Af=eng.α[1,b,1,kz]*eng.sqrtγ[1,b,1,kz]
                elseif kR==Kx+1
                    ρL=st.ρ[N,b,Kx,kz];pL=st.p[N,b,Kx,kz];vxL=st.vx[N,b,Kx,kz];vzL=st.vz[N,b,Kx,kz]
                    ρR=atm.ρ_atm;pR=atm.p_atm;vxR=0.0;vzR=0.0
                    elf=eng.elam[N,b,Kx,kz];nxf=eng.nx[N,b,Kx,kz];nzf=eng.nz[N,b,Kx,kz];Af=eng.α[N,b,Kx,kz]*eng.sqrtγ[N,b,Kx,kz]
                else
                    ρL=st.ρ[N,b,kL,kz];pL=st.p[N,b,kL,kz];vxL=st.vx[N,b,kL,kz];vzL=st.vz[N,b,kL,kz]
                    ρR=st.ρ[1,b,kR,kz];pR=st.p[1,b,kR,kz];vxR=st.vx[1,b,kR,kz];vzR=st.vz[1,b,kR,kz]
                    elf=0.5*(eng.elam[N,b,kL,kz]+eng.elam[1,b,kR,kz])
                    nxf=0.5*(eng.nx[N,b,kL,kz]+eng.nx[1,b,kR,kz]); nzf=0.5*(eng.nz[N,b,kL,kz]+eng.nz[1,b,kR,kz])
                    Af=0.5*(eng.α[N,b,kL,kz]*eng.sqrtγ[N,b,kL,kz]+eng.α[1,b,kR,kz]*eng.sqrtγ[1,b,kR,kz])
                end
                UL,FL,wL,εL=_faceflux(eos,ρL,pL,vxL,vzL,elf,nxf,nzf,1)
                UR,FR,wR,εR=_faceflux(eos,ρR,pR,vxR,vzR,elf,nxf,nzf,1)
                amax=max(abs.(_axisspeeds(eos_cs2(eos,εL),wL))...,abs.(_axisspeeds(eos_cs2(eos,εR),wR))...)
                F̂=_rus4(UL,FL,UR,FR,amax,Af)
                FLp=(Af*FL[1],Af*FL[2],Af*FL[3],Af*FL[4]); FRp=(Af*FR[1],Af*FR[2],Af*FR[3],Af*FR[4])
                if kL≥1
                    wend=w[N]; sg=eng.sqrtγ[N,b,kL,kz]
                    rD[N,b,kL,kz]-=(F̂[1]-FLp[1])/(J*wend*sg); rSx[N,b,kL,kz]-=(F̂[2]-FLp[2])/(J*wend*sg)
                    rSz[N,b,kL,kz]-=(F̂[3]-FLp[3])/(J*wend*sg); rτ[N,b,kL,kz]-=(F̂[4]-FLp[4])/(J*wend*sg)
                end
                if kR≤Kx
                    w1=w[1]; sg=eng.sqrtγ[1,b,kR,kz]
                    rD[1,b,kR,kz]+=(F̂[1]-FRp[1])/(J*w1*sg); rSx[1,b,kR,kz]+=(F̂[2]-FRp[2])/(J*w1*sg)
                    rSz[1,b,kR,kz]+=(F̂[3]-FRp[3])/(J*w1*sg); rτ[1,b,kR,kz]+=(F̂[4]-FRp[4])/(J*w1*sg)
                end
            end
        end
    end
    # SURFACE term — Z faces
    @inbounds Threads.@threads for kx in 1:Kx
        for kf in 0:Kz
            kL=kf; kR=kf+1
            for a in 1:N
                if kL==0
                    ρR=st.ρ[a,1,kx,1];pR=st.p[a,1,kx,1];vxR=st.vx[a,1,kx,1];vzR=st.vz[a,1,kx,1]
                    ρL=ρR;pL=pR;vxL=vxR;vzL=-vzR
                    elf=eng.elam[a,1,kx,1];nxf=eng.nx[a,1,kx,1];nzf=eng.nz[a,1,kx,1];Af=eng.α[a,1,kx,1]*eng.sqrtγ[a,1,kx,1]
                elseif kR==Kz+1
                    ρL=st.ρ[a,N,kx,Kz];pL=st.p[a,N,kx,Kz];vxL=st.vx[a,N,kx,Kz];vzL=st.vz[a,N,kx,Kz]
                    ρR=atm.ρ_atm;pR=atm.p_atm;vxR=0.0;vzR=0.0
                    elf=eng.elam[a,N,kx,Kz];nxf=eng.nx[a,N,kx,Kz];nzf=eng.nz[a,N,kx,Kz];Af=eng.α[a,N,kx,Kz]*eng.sqrtγ[a,N,kx,Kz]
                else
                    ρL=st.ρ[a,N,kx,kL];pL=st.p[a,N,kx,kL];vxL=st.vx[a,N,kx,kL];vzL=st.vz[a,N,kx,kL]
                    ρR=st.ρ[a,1,kx,kR];pR=st.p[a,1,kx,kR];vxR=st.vx[a,1,kx,kR];vzR=st.vz[a,1,kx,kR]
                    elf=0.5*(eng.elam[a,N,kx,kL]+eng.elam[a,1,kx,kR])
                    nxf=0.5*(eng.nx[a,N,kx,kL]+eng.nx[a,1,kx,kR]); nzf=0.5*(eng.nz[a,N,kx,kL]+eng.nz[a,1,kx,kR])
                    Af=0.5*(eng.α[a,N,kx,kL]*eng.sqrtγ[a,N,kx,kL]+eng.α[a,1,kx,kR]*eng.sqrtγ[a,1,kx,kR])
                end
                UL,FL,wL,εL=_faceflux(eos,ρL,pL,vxL,vzL,elf,nxf,nzf,2)
                UR,FR,wR,εR=_faceflux(eos,ρR,pR,vxR,vzR,elf,nxf,nzf,2)
                amax=max(abs.(_axisspeeds(eos_cs2(eos,εL),wL))...,abs.(_axisspeeds(eos_cs2(eos,εR),wR))...)
                F̂=_rus4(UL,FL,UR,FR,amax,Af)
                FLp=(Af*FL[1],Af*FL[2],Af*FL[3],Af*FL[4]); FRp=(Af*FR[1],Af*FR[2],Af*FR[3],Af*FR[4])
                if kL≥1
                    wend=w[N]; sg=eng.sqrtγ[a,N,kx,kL]
                    rD[a,N,kx,kL]-=(F̂[1]-FLp[1])/(J*wend*sg); rSx[a,N,kx,kL]-=(F̂[2]-FLp[2])/(J*wend*sg)
                    rSz[a,N,kx,kL]-=(F̂[3]-FLp[3])/(J*wend*sg); rτ[a,N,kx,kL]-=(F̂[4]-FLp[4])/(J*wend*sg)
                end
                if kR≤Kz
                    w1=w[1]; sg=eng.sqrtγ[a,1,kx,kR]
                    rD[a,1,kx,kR]+=(F̂[1]-FRp[1])/(J*w1*sg); rSx[a,1,kx,kR]+=(F̂[2]-FRp[2])/(J*w1*sg)
                    rSz[a,1,kx,kR]+=(F̂[3]-FRp[3])/(J*w1*sg); rτ[a,1,kx,kR]+=(F̂[4]-FRp[4])/(J*w1*sg)
                end
            end
        end
    end
    # geometric/gravity source (Cartesian projection of the spherical source)
    @inbounds Threads.@threads for kz in 1:Kz
        for kx in 1:Kx, b in 1:N, a in 1:N
            α=eng.α[a,b,kx,kz]; sg=eng.sqrtγ[a,b,kx,kz]; Φp=eng.Φp[a,b,kx,kz]
            nxv=eng.nx[a,b,kx,kz]; nzv=eng.nz[a,b,kx,kz]
            ρ=st.ρ[a,b,kx,kz]; p=st.p[a,b,kx,kz]; ε=st.ε[a,b,kx,kz]
            vux=st.vx[a,b,kx,kz]; vuz=st.vz[a,b,kx,kz]; el=eng.elam[a,b,kx,kz]; q=(el-1)
            nv=nxv*vux+nzv*vuz; vlx=vux+q*nxv*nv; vlz=vuz+q*nzv*nv
            v2=clamp(vux*vlx+vuz*vlz,0.0,1.0-1e-12); W=1.0/sqrt(1.0-v2)
            vn=nxv*vlx+nzv*vlz
            gforce=-α*(ε+p)*W^2*Φp     # per unit (undensitized): source/sg cancels sg
            rSx[a,b,kx,kz]+=gforce*nxv; rSz[a,b,kx,kz]+=gforce*nzv
            rτ[a,b,kx,kz]+= -α*(ε+p)*W^2*vn*Φp
        end
    end
    return nothing
end

function _rhs!(rD,rSx,rSz,rτ, st::DGCart2DState, eng::DGCart2DEngine)
    _raw_rhs!(rD,rSx,rSz,rτ, st, eng)
    @inbounds for I in eachindex(rD)
        rD[I]-=eng.Seq_D[I]; rSx[I]-=eng.Seq_Sx[I]; rSz[I]-=eng.Seq_Sz[I]; rτ[I]-=eng.Seq_τ[I]
    end
end

# Zhang–Shu positivity proxy on the densitized cons: q = (τ+D) − √(D²+|S|²) > 0
# guards p>0 (analogue of the 1D SRHD limiter), with the √γ-metric momentum norm.
@inline _q2d(D,Sx,Sz,τ,elam,nx,nz) = begin
    nS=nx*Sx+nz*Sz; S2=Sx^2+Sz^2+(1/elam-1)*nS^2
    (τ+D) - sqrt(D^2 + max(S2,0.0))
end

# limiter: Zhang–Shu positivity (D>0 AND p-proxy>0) squeezing toward the block
# mean — the essential cure for the steep ℓ=2-perturbed surface.
function _limit!(st::DGCart2DState, eng::DGCart2DEngine)
    bs=eng.bs; w=bs.w; N=bs.N; Kx=eng.Kx; Kz=eng.Kz
    @inbounds Threads.@threads for kz in 1:Kz
        for kx in 1:Kx
            D̄=_cellavg2(view(st.D,:,:,kx,kz),w)
            S̄x=_cellavg2(view(st.Sx,:,:,kx,kz),w); S̄z=_cellavg2(view(st.Sz,:,:,kx,kz),w)
            τ̄=_cellavg2(view(st.τ,:,:,kx,kz),w)
            Dε=eng.atm.ρ_atm*_cellavg2(view(eng.sqrtγ,:,:,kx,kz),w)
            # --- step 1: D>0 at nodes ---
            D̄=max(D̄,Dε)
            Dmin=Inf; for b in 1:N,a in 1:N; Dmin=min(Dmin,st.D[a,b,kx,kz]); end
            if Dmin<Dε
                θ=clamp((D̄-Dε)/(D̄-Dmin+1e-300),0.0,1.0)
                for b in 1:N,a in 1:N; st.D[a,b,kx,kz]=θ*(st.D[a,b,kx,kz]-D̄)+D̄; end
            end
            # --- step 2: p-proxy q>0 at nodes ---
            elf=eng.elam[1,1,kx,kz]; nxf=eng.nx[1,1,kx,kz]; nzf=eng.nz[1,1,kx,kz]
            q̄=_q2d(D̄,S̄x,S̄z,τ̄,elf,nxf,nzf)
            if q̄ ≤ 0.0
                for b in 1:N,a in 1:N
                    st.D[a,b,kx,kz]=max(D̄,Dε); st.Sx[a,b,kx,kz]=S̄x; st.Sz[a,b,kx,kz]=S̄z; st.τ[a,b,kx,kz]=τ̄
                end
                continue
            end
            qε=1e-12*q̄; θq=1.0
            for b in 1:N,a in 1:N
                qn=_q2d(st.D[a,b,kx,kz],st.Sx[a,b,kx,kz],st.Sz[a,b,kx,kz],st.τ[a,b,kx,kz],
                        eng.elam[a,b,kx,kz],eng.nx[a,b,kx,kz],eng.nz[a,b,kx,kz])
                if qn<qε
                    θq=min(θq, clamp((q̄-qε)/(q̄-qn+1e-300),0.0,1.0))
                end
            end
            if θq<1.0
                for b in 1:N,a in 1:N
                    st.D[a,b,kx,kz]=θq*(st.D[a,b,kx,kz]-D̄)+D̄
                    st.Sx[a,b,kx,kz]=θq*(st.Sx[a,b,kx,kz]-S̄x)+S̄x
                    st.Sz[a,b,kx,kz]=θq*(st.Sz[a,b,kx,kz]-S̄z)+S̄z
                    st.τ[a,b,kx,kz]=θq*(st.τ[a,b,kx,kz]-τ̄)+τ̄
                end
            end
        end
    end
end

@inline function _maxspeed(st::DGCart2DState, eng::DGCart2DEngine)
    eos=eng.eos; N=eng.bs.N; a=0.0
    @inbounds for kz in 1:eng.Kz, kx in 1:eng.Kx, b in 1:N, aa in 1:N
        cs=sqrt(clamp(eos_cs2(eos,st.ε[aa,b,kx,kz]),0.0,1.0))
        el=eng.elam[aa,b,kx,kz]; nxv=eng.nx[aa,b,kx,kz]; nzv=eng.nz[aa,b,kx,kz]
        vux=st.vx[aa,b,kx,kz]; vuz=st.vz[aa,b,kx,kz]; q=el-1; nv=nxv*vux+nzv*vuz
        vlx=vux+q*nxv*nv; vlz=vuz+q*nzv*nv
        vmag=sqrt(clamp(vux*vlx+vuz*vlz,0.0,1.0))
        a=max(a,(vmag+cs)/(1+vmag*cs)*sqrt(el))
    end
    return a
end

"""
    evolve_dgcart2d!(st, eng; tmax, sample_dt) -> (ts, q2, ρc)

SSP-RK3 evolution of the 2D meridional DG star; records the ℓ=2 quadrupole and
central density."""
function evolve_dgcart2d!(st::DGCart2DState, eng::DGCart2DEngine; tmax::Float64,
        sample_dt::Float64=4.0, cfl::Float64=-1.0, verbose::Bool=false)
    bs=eng.bs; N=bs.N; Kx=eng.Kx; Kz=eng.Kz; p=bs.p
    cfl=cfl>0 ? cfl : eng.cfl; cfl_dg=cfl/(2p+1)
    dims=size(st.D)
    rD=zeros(dims);rSx=zeros(dims);rSz=zeros(dims);rτ=zeros(dims)
    D0=zeros(dims);Sx0=zeros(dims);Sz0=zeros(dims);τ0=zeros(dims)
    ts=Float64[];q2=Float64[];ρch=Float64[]; ρc0=dgcart2d_central_density(st,eng)
    t=0.0; last=-1e30; ns=0
    _update_prims!(st,eng)
    while t<tmax && ns<5_000_000
        amax=max(_maxspeed(st,eng),1e-3); dt=cfl_dg*eng.Δ/amax; dt=min(dt,tmax-t)
        copyto!(D0,st.D);copyto!(Sx0,st.Sx);copyto!(Sz0,st.Sz);copyto!(τ0,st.τ)
        _rhs!(rD,rSx,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D); st.D[I]=D0[I]+dt*rD[I];st.Sx[I]=Sx0[I]+dt*rSx[I];st.Sz[I]=Sz0[I]+dt*rSz[I];st.τ[I]=τ0[I]+dt*rτ[I]; end
        _limit!(st,eng);_update_prims!(st,eng)
        _rhs!(rD,rSx,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=0.75*D0[I]+0.25*(st.D[I]+dt*rD[I]);st.Sx[I]=0.75*Sx0[I]+0.25*(st.Sx[I]+dt*rSx[I])
            st.Sz[I]=0.75*Sz0[I]+0.25*(st.Sz[I]+dt*rSz[I]);st.τ[I]=0.75*τ0[I]+0.25*(st.τ[I]+dt*rτ[I])
        end
        _limit!(st,eng);_update_prims!(st,eng)
        _rhs!(rD,rSx,rSz,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=(1/3)*D0[I]+(2/3)*(st.D[I]+dt*rD[I]);st.Sx[I]=(1/3)*Sx0[I]+(2/3)*(st.Sx[I]+dt*rSx[I])
            st.Sz[I]=(1/3)*Sz0[I]+(2/3)*(st.Sz[I]+dt*rSz[I]);st.τ[I]=(1/3)*τ0[I]+(2/3)*(st.τ[I]+dt*rτ[I])
        end
        _limit!(st,eng);_update_prims!(st,eng)
        t+=dt; ns+=1
        if t-last≥sample_dt
            push!(ts,t);push!(q2,dgcart2d_quadrupole(st,eng));push!(ρch,dgcart2d_central_density(st,eng))
            last=t
            verbose && ns%100==0 && println("  t=",round(t,digits=1)," ρc/ρc0=",round(ρch[end]/ρc0,digits=5))
        end
        if !isfinite(ρch[end]) || ρch[end]>100*ρc0; return ts,q2,ρch end
    end
    return ts,q2,ρch
end

dgcart2d_central_density(st::DGCart2DState, eng::DGCart2DEngine) =
    _cellavg2(view(st.ρ,:,:,1,1), eng.bs.w)

"""seed ℓ=2,m=0 velocity v^r=A(r/R)(3cos²θ−1)/... ; in x–z plane Y20∝(2n_z²−n_x²)."""
function seed_dgcart2d_l2!(st::DGCart2DState, eng::DGCart2DEngine; A::Float64=1e-4)
    eos=eng.eos; N=eng.bs.N
    @inbounds for kz in 1:eng.Kz, kx in 1:eng.Kx, b in 1:N, a in 1:N
        eng.interior[a,b,kx,kz] || continue
        rr=eng.r[a,b,kx,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,kx,kz]; nzv=eng.nz[a,b,kx,kz]
        Y20=2*nzv^2-nxv^2          # axisymmetric ℓ=2,m=0 in x–z plane (n_y=0)
        vr=A*(rr/eng.R)*Y20
        el2=sqrt(eng.elam[a,b,kx,kz])
        vux=vr/el2*nxv; vuz=vr/el2*nzv
        st.vx[a,b,kx,kz]=vux; st.vz[a,b,kx,kz]=vuz
        D,Sx,Sz,τ=_p2c(eos,st.ρ[a,b,kx,kz],st.p[a,b,kx,kz],vux,vuz,
            eng.sqrtγ[a,b,kx,kz],eng.elam[a,b,kx,kz],nxv,nzv)
        st.D[a,b,kx,kz]=D; st.Sx[a,b,kx,kz]=Sx; st.Sz[a,b,kx,kz]=Sz; st.τ[a,b,kx,kz]=τ
    end
end

"""√γ-weighted ℓ=2 quadrupole of ρ over the interior (meridional)."""
function dgcart2d_quadrupole(st::DGCart2DState, eng::DGCart2DEngine)
    N=eng.bs.N; w=eng.bs.w; acc=0.0
    @inbounds for kz in 1:eng.Kz, kx in 1:eng.Kx, b in 1:N, a in 1:N
        eng.interior[a,b,kx,kz] || continue
        rr=eng.r[a,b,kx,kz]; rr≤1e-10 && continue
        nxv=eng.nx[a,b,kx,kz]; nzv=eng.nz[a,b,kx,kz]
        acc+=st.ρ[a,b,kx,kz]*(2nzv^2-nxv^2)*eng.sqrtγ[a,b,kx,kz]*w[a]*w[b]
    end
    acc*eng.J^2
end

end # module DGCart2D
