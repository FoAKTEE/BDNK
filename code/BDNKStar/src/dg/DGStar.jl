#=
    DGStar — STAGE 2: nodal RKDG on the 1D radial TOV (Cowling) star with an
    atmosphere + positivity limiter. The goal mirrors FVRadial but with sub-cell
    polynomial resolution at the steep stellar surface.

    Spherical Valencia GRHydro on the frozen TOV metric (areal Schwarzschild-like,
    Cowling, β=0):
        ds² = −e^{ν}dt² + e^{λ}dr² + r²dΩ²,  α=e^{ν/2}, e^{λ/2}=1/√(1−2m/r).
    We evolve the UNDENSITIZED orthonormal conserved variables (D̂,Ŝ,τ̂) per node
    with the well-balanced flux/source split: the static-background DG RHS is
    stored and subtracted, so the unperturbed star has machine-zero RHS (the same
    lake-at-rest device as FVRadial). The atmosphere floor + Zhang–Shu positivity
    limiter keep ρ,p>0; a TVB troubled-cell limiter monotonizes the surface
    element.

    Reuses the ShumPolytrope barotrope and the FVCommon cons2prim/prim2cons.
=#
module DGStar

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, cons2prim_barotrope, prim2cons_barotrope,
                  rho_from_p, eos_cs2
using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis, tvb_minmod, cell_average
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram

export DGStarEngine, DGStarState, setup_dgstar, evolve_dgstar!,
       seed_dgstar_radial!, dgstar_central_density, dgstar_surface_width,
       dgstar_radial_freq

@inline _lininterp(xs,ys,x) = begin
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end
@inline _p_of_rho(eos::ShumPolytrope, ρ) = eos.κ*ρ^2

struct DGStarEngine
    b::LGLBasis
    K::Int
    Δr::Float64
    J::Float64
    eos::BarotropicEOS
    atm::AtmospherePars
    R::Float64; M::Float64
    cfl::Float64
    M_tvb::Float64
    # node coords (N,K)
    r::Matrix{Float64}
    α::Matrix{Float64}
    elam2::Matrix{Float64}      # e^{λ/2}
    Φp::Matrix{Float64}
    sqrtg::Matrix{Float64}      # e^{λ/2} r²
    # well-balanced equilibrium RHS (subtracted)
    Seq_D::Matrix{Float64}; Seq_S::Matrix{Float64}; Seq_τ::Matrix{Float64}
end

mutable struct DGStarState
    D::Matrix{Float64}; S::Matrix{Float64}; τ::Matrix{Float64}
    ρ::Matrix{Float64}; p::Matrix{Float64}; ε::Matrix{Float64}; v::Matrix{Float64}
end
DGStarState(N,K)=DGStarState((zeros(N,K) for _ in 1:7)...)

@inline function _phys_flux(eos,ρ,p,v)
    ε=energy_from_pressure(eos,p); h=ρ>0 ? (ε+p)/ρ : 1.0
    W=1.0/sqrt(1.0-clamp(v^2,0.0,1.0-1e-12))
    D̂=ρ*W; Ŝ=ρ*h*W^2*v; τ̂=ρ*h*W^2-p-ρ*W
    return (D̂,Ŝ,τ̂),(D̂*v, Ŝ*v+p, (τ̂+p)*v)
end
@inline function _wavespeeds(eos,ρ,p,ε,v)
    cs=sqrt(clamp(eos_cs2(eos,ε),0.0,1.0)); v2=clamp(v^2,0.0,1.0-1e-12)
    a=1.0/(1.0-v2*cs^2); disc=sqrt(max(cs^2*(1-v2)*(1-v2*cs^2-v^2*(1-cs^2)),0.0))
    return a*(v*(1-cs^2)-disc), a*(v*(1-cs^2)+disc)
end

"""
    setup_dgstar(eos, εc; K=64, p=3, rmax_fac=1.3, atm_fac=1e-7, cfl=0.2,
                 M_tvb=50.0, h_tov=2e-4) -> (engine, state)

Build the radial DG star engine on the TOV star of central energy density `εc`.
`K` elements of degree `p` over [0, rmax_fac·R].
"""
function setup_dgstar(eos::BarotropicEOS, εc::Float64; K::Int=64, p::Int=3,
        rmax_fac::Float64=1.3, atm_fac::Float64=1e-7, cfl::Float64=0.2,
        M_tvb::Float64=50.0, h_tov::Float64=2e-4, rin_fac::Float64=0.05)
    star=solve_tov(eos,εc;h=h_tov); R,M=star.R,star.M
    rmax=rmax_fac*R
    # excise a tiny inner core r_in>0 (reflecting boundary there) so NO LGL node
    # sits on the r=0 coordinate singularity — the standard cure for the 1/r
    # noise growth that destabilizes spherical nodal DG at the center.
    rin=rin_fac*R
    Δr=(rmax-rin)/K; J=Δr/2
    b=build_lgl_basis(p); N=b.N
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p
    m_of(r)= r≤R ? _lininterp(rt,mt,r) : M
    ν_of(r)= r≤R ? _lininterp(rt,νt,r) : log(1-2M/r)
    p_of(r)= r≤R ? max(_lininterp(rt,pt,r),0.0) : 0.0

    r=zeros(N,K); α=zeros(N,K); elam2=zeros(N,K); Φp=zeros(N,K); sqrtg=zeros(N,K)
    for k in 1:K
        rc = rin + (k-0.5)*Δr
        for a in 1:N
            rr = rc + J*b.ξ[a]; rr=max(rr,1e-8)
            r[a,k]=rr
            m=m_of(rr); f=max(1-2m/rr,1e-12)
            α[a,k]=exp(0.5*ν_of(rr)); elam2[a,k]=1/sqrt(f)
            den=rr*(rr-2m); Φp[a,k]= den>0 ? (m+4π*rr^3*p_of(rr))/den : 0.0
            sqrtg[a,k]=elam2[a,k]*rr^2
        end
    end

    ρc=rho_from_p(eos,pressure(eos,εc)); ρ_atm=atm_fac*ρc; p_atm=_p_of_rho(eos,ρ_atm)
    ε_atm=energy_from_pressure(eos,p_atm)
    atm=AtmospherePars(ρ_atm,p_atm,ε_atm,5*ρ_atm,0.999)

    st=DGStarState(N,K)
    for k in 1:K, a in 1:N
        rr=r[a,k]
        ρ = rr≤R ? max(rho_from_p(eos,p_of(rr)),ρ_atm) : ρ_atm
        pp= rr≤R ? max(p_of(rr),p_atm) : p_atm
        εe=energy_from_pressure(eos,pp)
        st.ρ[a,k]=ρ; st.p[a,k]=pp; st.ε[a,k]=εe; st.v[a,k]=0.0
        (D̂,Ŝ,τ̂),_=_phys_flux(eos,ρ,pp,0.0)
        sg=sqrtg[a,k]
        st.D[a,k]=sg*D̂; st.S[a,k]=sg*Ŝ; st.τ[a,k]=sg*τ̂   # DENSITIZED
    end
    eng=DGStarEngine(b,K,Δr,J,eos,atm,R,M,cfl,M_tvb,r,α,elam2,Φp,sqrtg,
                     zeros(N,K),zeros(N,K),zeros(N,K))
    _update_prims!(st,eng)
    _raw_rhs!(eng.Seq_D,eng.Seq_S,eng.Seq_τ, st, eng)
    return eng, st
end

function _update_prims!(st::DGStarState, eng::DGStarEngine)
    eos=eng.eos; atm=eng.atm
    @inbounds for k in 1:eng.K, a in 1:eng.b.N
        sg=eng.sqrtg[a,k]
        if sg ≤ 1e-30
            st.ρ[a,k]=atm.ρ_atm; st.p[a,k]=atm.p_atm; st.ε[a,k]=atm.ε_atm; st.v[a,k]=0.0
            continue
        end
        D̂=st.D[a,k]/sg; Ŝ=st.S[a,k]/sg; τ̂=st.τ[a,k]/sg   # UNdensitize
        ρ,pp,ε,vmag,W,isatm = cons2prim_barotrope(eos,D̂,abs(Ŝ),τ̂,atm)
        v= isatm ? 0.0 : sign(Ŝ)*vmag
        st.ρ[a,k]=ρ; st.p[a,k]=pp; st.ε[a,k]=ε; st.v[a,k]=v
        if isatm
            (D̂2,Ŝ2,τ̂2),_=_phys_flux(eos,ρ,pp,0.0)
            st.D[a,k]=sg*D̂2; st.S[a,k]=sg*Ŝ2; st.τ[a,k]=sg*τ̂2
        end
    end
end

# DG RHS (raw): flux divergence with Rusanov interface flux + spherical geometric
# source. Inner face r=0: reflecting (S odd). Outer: atmosphere outflow.
function _raw_rhs!(rD,rS,rτ, st::DGStarState, eng::DGStarEngine)
    b=eng.b; N=b.N; K=eng.K; J=eng.J; eos=eng.eos; atm=eng.atm
    fill!(rD,0.0); fill!(rS,0.0); fill!(rτ,0.0)
    # node fluxes
    fD=zeros(N,K); fS=zeros(N,K); fτ=zeros(N,K)
    @inbounds for k in 1:K, a in 1:N
        ρ=st.ρ[a,k]; pp=st.p[a,k]; v=st.v[a,k]
        _,F=_phys_flux(eos,ρ,pp,v)
        # densitized flux uses area A=α e^{λ/2} r² ; but for well-balanced split
        # we use the area-weighted flux divergence consistent with FVRadial.
        A=eng.α[a,k]*eng.sqrtg[a,k]
        fD[a,k]=A*F[1]; fS[a,k]=A*F[2]; fτ[a,k]=A*F[3]
    end
    Dm=b.D
    # DENSITIZED evolution: ∂_t U = -(1/J) D̃(A F) + sg·source. NO 1/sg division
    # (this is what kills the 1/r noise amplification that broke the undensitized
    # form and the fine-Δr FV surface).
    @inbounds for k in 1:K, a in 1:N
        sD=0.0;sS=0.0;sτ=0.0
        for c in 1:N
            d=Dm[a,c]; sD+=d*fD[c,k]; sS+=d*fS[c,k]; sτ+=d*fτ[c,k]
        end
        rD[a,k]=-sD/J; rS[a,k]=-sS/J; rτ[a,k]=-sτ/J
    end
    # surface (Rusanov) flux at faces
    @inbounds for kf in 0:K
        kL=kf; kR=kf+1
        if kL==0
            # inner r=0 reflecting: mirror first node with S flipped
            ρR=st.ρ[1,1]; pR=st.p[1,1]; vR=st.v[1,1]
            ρL=ρR; pL=pR; vL=-vR
            UL,FL=_phys_flux(eos,ρL,pL,vL); UR,FR=_phys_flux(eos,ρR,pR,vR)
        elseif kR==K+1
            ρL=st.ρ[N,K]; pL=st.p[N,K]; vL=st.v[N,K]
            ρR=atm.ρ_atm; pR=atm.p_atm; vR=0.0
            UL,FL=_phys_flux(eos,ρL,pL,vL); UR,FR=_phys_flux(eos,ρR,pR,vR)
        else
            ρL=st.ρ[N,kL]; pL=st.p[N,kL]; vL=st.v[N,kL]
            ρR=st.ρ[1,kR]; pR=st.p[1,kR]; vR=st.v[1,kR]
            UL,FL=_phys_flux(eos,ρL,pL,vL); UR,FR=_phys_flux(eos,ρR,pR,vR)
        end
        εL=energy_from_pressure(eos,pL); εR=energy_from_pressure(eos,pR)
        λmL,λpL=_wavespeeds(eos,ρL,pL,εL,vL); λmR,λpR=_wavespeeds(eos,ρR,pR,εR,vR)
        amax=max(abs(λmL),abs(λpL),abs(λmR),abs(λpR))
        # area at face: use the element-boundary sqrtg·α (continuous in r)
        if kL≥1; AL=eng.α[N,kL]*eng.sqrtg[N,kL] else AL=eng.α[1,1]*eng.sqrtg[1,1] end
        if kR≤K; AR=eng.α[1,kR]*eng.sqrtg[1,kR] else AR=eng.α[N,K]*eng.sqrtg[N,K] end
        Af=0.5*(AL+AR)
        F̂=ntuple(i->0.5*Af*(FL[i]+FR[i]) - 0.5*amax*Af*(UR[i]-UL[i]), 3)
        FLp=(Af*FL[1],Af*FL[2],Af*FL[3]); FRp=(Af*FR[1],Af*FR[2],Af*FR[3])
        if kL≥1
            wend=b.w[N]
            rD[N,kL]-=(F̂[1]-FLp[1])/(J*wend)
            rS[N,kL]-=(F̂[2]-FLp[2])/(J*wend)
            rτ[N,kL]-=(F̂[3]-FLp[3])/(J*wend)
        end
        if kR≤K
            w1=b.w[1]
            rD[1,kR]+=(F̂[1]-FRp[1])/(J*w1)
            rS[1,kR]+=(F̂[2]-FRp[2])/(J*w1)
            rτ[1,kR]+=(F̂[3]-FRp[3])/(J*w1)
        end
    end
    # geometric source (DENSITIZED, same as FVRadial): rhs for U gets sg·source.
    @inbounds for k in 1:K, a in 1:N
        rr=eng.r[a,k]; α=eng.α[a,k]; el2=eng.elam2[a,k]; Φp=eng.Φp[a,k]
        ρ=st.ρ[a,k]; pp=st.p[a,k]; ε=st.ε[a,k]; v=st.v[a,k]
        W=1.0/sqrt(1.0-clamp(v^2,0.0,1.0-1e-12)); sg=eng.sqrtg[a,k]
        rS[a,k]+= α*sg*(-(ε+pp)*W^2*Φp) + 2.0*α*el2*rr*pp
        rτ[a,k]+= α*sg*(-(ε+pp)*W^2*v*Φp)
    end
    return nothing
end

function _rhs!(rD,rS,rτ, st::DGStarState, eng::DGStarEngine)
    _raw_rhs!(rD,rS,rτ, st, eng)
    @inbounds for I in eachindex(rD)
        rD[I]-=eng.Seq_D[I]; rS[I]-=eng.Seq_S[I]; rτ[I]-=eng.Seq_τ[I]
    end
end

# positivity + TVB limiter on the conserved fields (operate on D̂ mainly)
function _limit!(st::DGStarState, eng::DGStarEngine)
    b=eng.b; w=b.w; ξ=b.ξ; N=b.N; K=eng.K; Mh2=eng.M_tvb*eng.Δr^2
    for U in (st.D, st.S, st.τ)
        means=[cell_average(view(U,:,k),w) for k in 1:K]
        new=copy(U)
        @inbounds for k in 1:K
            ū=means[k]; ūm=k>1 ? means[k-1] : means[k]; ūp=k<K ? means[k+1] : means[k]
            s_int=0.0; for a in 1:N; s_int+=1.5*w[a]*ξ[a]*U[a,k]; end
            sl=tvb_minmod(s_int, ūp-ū, ū-ūm, Mh2)
            if abs(sl-s_int)>1e-13*(1+abs(s_int))
                for a in 1:N; new[a,k]=ū+sl*ξ[a]; end
            end
        end
        copyto!(U,new)
    end
    # positivity: floor densitized D=sg·ρW toward mean. Floor per node uses the
    # local sg·ρ_atm (so the atmosphere is the admissible lower bound).
    @inbounds for k in 1:K
        D̄=cell_average(view(st.D,:,k),w)
        Dε_mean=eng.atm.ρ_atm*cell_average(view(eng.sqrtg,:,k),w)
        D̄=max(D̄, Dε_mean)
        # find the worst node ratio (D - sg·ρ_atm) and squeeze toward mean
        θmin=1.0
        for a in 1:N
            Dε=eng.atm.ρ_atm*eng.sqrtg[a,k]
            if st.D[a,k] < Dε
                θ=clamp((D̄-Dε)/(D̄-st.D[a,k]+1e-300),0.0,1.0)
                θmin=min(θmin,θ)
            end
        end
        if θmin<1.0
            for a in 1:N; st.D[a,k]=θmin*(st.D[a,k]-D̄)+D̄; end
        end
    end
end

@inline function _maxspeed(st::DGStarState, eng::DGStarEngine)
    eos=eng.eos; a=0.0
    @inbounds for k in 1:eng.K, i in 1:eng.b.N
        λm,λp=_wavespeeds(eos,st.ρ[i,k],st.p[i,k],st.ε[i,k],st.v[i,k])
        a=max(a,abs(λm),abs(λp))
    end
    return a
end

"""
    evolve_dgstar!(st, eng; tmax, sample_dt, probe_frac=0.3) -> (ts, probe, ρc, drift)

SSP-RK3 evolution of the radial DG star with limiter. Records a velocity probe
and the central density (mean of first element)."""
function evolve_dgstar!(st::DGStarState, eng::DGStarEngine; tmax::Float64,
        sample_dt::Float64=2.0, probe_frac::Float64=0.3, cfl::Float64=-1.0)
    b=eng.b; N=b.N; K=eng.K; p=b.p
    cfl=cfl>0 ? cfl : eng.cfl
    cfl_dg=cfl/(2p+1)
    rD=zeros(N,K); rS=zeros(N,K); rτ=zeros(N,K)
    D0=zeros(N,K); S0=zeros(N,K); τ0=zeros(N,K)
    kp=clamp(round(Int, probe_frac*K),1,K)
    ρc0=dgstar_central_density(st,eng)
    ts=Float64[]; probe=Float64[]; ρch=Float64[]; drift=0.0; t=0.0; last=-1e30; ns=0
    _update_prims!(st,eng)
    while t<tmax && ns<5_000_000
        amax=max(_maxspeed(st,eng),1e-3); dt=cfl_dg*eng.Δr/amax; dt=min(dt,tmax-t)
        copyto!(D0,st.D); copyto!(S0,st.S); copyto!(τ0,st.τ)
        _rhs!(rD,rS,rτ,st,eng)
        @inbounds for I in eachindex(st.D); st.D[I]=D0[I]+dt*rD[I]; st.S[I]=S0[I]+dt*rS[I]; st.τ[I]=τ0[I]+dt*rτ[I]; end
        _limit!(st,eng); _update_prims!(st,eng)
        _rhs!(rD,rS,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=0.75*D0[I]+0.25*(st.D[I]+dt*rD[I]); st.S[I]=0.75*S0[I]+0.25*(st.S[I]+dt*rS[I]); st.τ[I]=0.75*τ0[I]+0.25*(st.τ[I]+dt*rτ[I])
        end
        _limit!(st,eng); _update_prims!(st,eng)
        _rhs!(rD,rS,rτ,st,eng)
        @inbounds for I in eachindex(st.D)
            st.D[I]=(1/3)*D0[I]+(2/3)*(st.D[I]+dt*rD[I]); st.S[I]=(1/3)*S0[I]+(2/3)*(st.S[I]+dt*rS[I]); st.τ[I]=(1/3)*τ0[I]+(2/3)*(st.τ[I]+dt*rτ[I])
        end
        _limit!(st,eng); _update_prims!(st,eng)
        t+=dt; ns+=1
        if t-last≥sample_dt
            push!(ts,t); push!(probe,cell_average(view(st.v,:,kp),b.w)); push!(ρch,dgstar_central_density(st,eng))
            last=t; drift=max(drift,abs(dgstar_central_density(st,eng)-ρc0)/ρc0)
        end
        if !isfinite(ρch==[] ? ρc0 : ρch[end]) || dgstar_central_density(st,eng)>100*ρc0
            return ts,probe,ρch,Inf
        end
    end
    return ts,probe,ρch,drift
end

dgstar_central_density(st::DGStarState, eng::DGStarEngine) =
    cell_average(view(st.ρ,:,1), eng.b.w)

"""seed homologous radial velocity v=A r/R inside the star."""
function seed_dgstar_radial!(st::DGStarState, eng::DGStarEngine; A::Float64=1e-4)
    eos=eng.eos
    @inbounds for k in 1:eng.K, a in 1:eng.b.N
        rr=eng.r[a,k]
        if rr<eng.R
            v=A*(rr/eng.R); st.v[a,k]=v
            (D̂,Ŝ,τ̂),_=_phys_flux(eos,st.ρ[a,k],st.p[a,k],v)
            sg=eng.sqrtg[a,k]
            st.D[a,k]=sg*D̂; st.S[a,k]=sg*Ŝ; st.τ[a,k]=sg*τ̂
        end
    end
end

"""
    dgstar_surface_width(st, eng) -> (Δr_cells, R)

Number of radial CELLS over which ρ drops from 50% to 5% of central — the
'surface sharpness' diagnostic. A sub-cell-resolved surface gives <~1 cell here.
"""
function dgstar_surface_width(st::DGStarState, eng::DGStarEngine)
    K=eng.K; w=eng.b.w
    ρm=[cell_average(view(st.ρ,:,k),w) for k in 1:K]
    rm=[ (k-0.5)*eng.Δr for k in 1:K ]
    ρc=ρm[1]
    r50=NaN; r05=NaN
    for k in 1:K-1
        if isnan(r50) && ρm[k]≥0.5ρc && ρm[k+1]<0.5ρc; r50=rm[k] end
        if isnan(r05) && ρm[k]≥0.05ρc && ρm[k+1]<0.05ρc; r05=rm[k]; break end
    end
    width = (isnan(r50)||isnan(r05)) ? NaN : (r05-r50)/eng.Δr
    return width, eng.R
end

function dgstar_radial_freq(ts::Vector{Float64}, ρc::Vector{Float64};
        fmin_kHz=2.0, fmax_kHz=4.0, npts=4000, npeaks=1)
    νmin=fmin_kHz*Msun_to_km*kHz_to_km; νmax=fmax_kHz*Msun_to_km*kHz_to_km
    νgrid=range(νmin,νmax;length=npts)
    q=copy(ρc); n=length(q); μ=sum(q)/n
    @inbounds for i in 1:n; q[i]=(q[i]-μ)*(0.5-0.5*cos(2π*(i-1)/(n-1))); end
    P=periodogram(ts,q,νgrid)
    peaks=Tuple{Float64,Float64}[]
    for m in 2:length(P)-1; (P[m]>P[m-1]&&P[m]≥P[m+1])&&push!(peaks,(P[m],νgrid[m])) end
    sort!(peaks,by=x->-x[1]); k=min(npeaks,length(peaks))
    k==0 && return Float64[]
    return [pk[2]/(Msun_to_km*kHz_to_km) for pk in peaks[1:k]]
end

end # module DGStar
