#=
    ConformalEvolution — full flat-space conformal-BDNK evolution (slab symmetry),
    a faithful Julia port of ref-code/1D_conformal_bdnk/solver.c (Pandya;
    arXiv:2201.12317). The evolution engine STAGE 1C inherits:

      * 5th-order WENO face reconstruction (WENO_reconst_qLx/qRx)
      * Kurganov–Tadmor central flux at light speed (a = c = 1)
      * Heun (SSP-RK2) predictor–corrector
      * BDNK primitive solve for ξ̇,u̇ each substage (ConformalBDNK.recover_time_derivs)

    Conserved: T^{tt}, T^{tx}. Primitives: ξ=ln ε, u=u^x (trivially advanced by
    the recovered ξ̇,u̇). Reproduces the steady-shock stationary solution and the
    smooth Gaussian evolution.
=#
module ConformalEvolution

using ..ConformalBDNK
using ..ConformalBDNK: T_tt, T_tx, T_xx

export ConfState, init_gaussian, init_smooth_shock, evolve!, energy_density

const NG = 3   # ghost cells per side

mutable struct ConfState
    fr::ConformalFrame
    x::Vector{Float64}
    dx::Float64
    dt::Float64
    ξ::Vector{Float64};  u::Vector{Float64}
    ξx::Vector{Float64}; ux::Vector{Float64}     # spatial derivatives
    ξt::Vector{Float64}; ut::Vector{Float64}     # time derivatives (recovered)
    Ttt::Vector{Float64}; Ttx::Vector{Float64}
    periodic::Bool
end

energy_density(s::ConfState) = exp.(s.ξ)

# --- WENO5 face reconstruction (solver.c WENO_reconst_qLx / qRx) -------------
# _qLx = left state at face i+1/2 (C reconst_qLx); _qRx = right state (C reconst_qRx)
const _EPSW = 1e-3
@inline function _qLx(q, i)        # solver.c WENO_reconst_qLx
    v0 = (1/3)q[i]   + (5/6)q[i+1] - (1/6)q[i+2]
    v1 = (-1/6)q[i-1] + (5/6)q[i]  + (1/3)q[i+1]
    v2 = (1/3)q[i-2] - (7/6)q[i-1] + (11/6)q[i]
    d0, d1, d2 = 3/10, 3/5, 1/10
    b0 = 0.25*(3q[i]-4q[i+1]+q[i+2])^2 + (13/12)*(q[i]-2q[i+1]+q[i+2])^2
    b1 = 0.25*(q[i+1]-q[i-1])^2        + (13/12)*(q[i-1]-2q[i]+q[i+1])^2
    b2 = 0.25*(q[i-2]-4q[i-1]+3q[i])^2 + (13/12)*(q[i-2]-2q[i-1]+q[i])^2
    w0=d0/(_EPSW+b0)^2; w1=d1/(_EPSW+b1)^2; w2=d2/(_EPSW+b2)^2; s=w0+w1+w2
    return (w0*v0+w1*v1+w2*v2)/s
end
@inline function _qRx(q, i)        # solver.c WENO_reconst_qRx (shifts to i+1)
    j = i+1
    v0 = (11/6)q[j]   - (7/6)q[j+1] + (1/3)q[j+2]
    v1 = (1/3)q[j-1]  + (5/6)q[j]   - (1/6)q[j+1]
    v2 = (-1/6)q[j-2] + (5/6)q[j-1] + (1/3)q[j]
    d0, d1, d2 = 1/10, 3/5, 3/10
    b0 = 0.25*(3q[j]-4q[j+1]+q[j+2])^2 + (13/12)*(q[j]-2q[j+1]+q[j+2])^2
    b1 = 0.25*(q[j+1]-q[j-1])^2        + (13/12)*(q[j-1]-2q[j]+q[j+1])^2
    b2 = 0.25*(q[j-2]-4q[j-1]+3q[j])^2 + (13/12)*(q[j-2]-2q[j-1]+q[j])^2
    w0=d0/(_EPSW+b0)^2; w1=d1/(_EPSW+b1)^2; w2=d2/(_EPSW+b2)^2; s=w0+w1+w2
    return (w0*v0+w1*v1+w2*v2)/s
end

# centered WENO 4th-order first derivative (solver.c wDx); 2nd-order at edges
@inline function _Dx(f, i, dx, N)
    if i > 2 && i < N-1
        b1 = 0.25*(f[i-2]-4f[i-1]+3f[i])^2 + (13/12)*(f[i-2]-2f[i-1]+f[i])^2
        b2 = 0.25*(f[i+1]-f[i-1])^2        + (13/12)*(f[i-1]-2f[i]+f[i+1])^2
        b3 = 0.25*(3f[i]-4f[i+1]+f[i+2])^2 + (13/12)*(f[i]-2f[i+1]+f[i+2])^2
        a1=(1/6)/(_EPSW+b1)^2; a2=(2/3)/(_EPSW+b2)^2; a3=(1/6)/(_EPSW+b3)^2; s=a1+a2+a3
        D1=(f[i-2]-4f[i-1]+3f[i])/(2dx); D2=(f[i+1]-f[i-1])/(2dx); D3=(-3f[i]+4f[i+1]-f[i+2])/(2dx)
        return (a1*D1 + a2*D2 + a3*D3)/s
    end
    return (f[i+1]-f[i-1])/(2dx)
end

# --- ghost cells ------------------------------------------------------------
function _set_ghost!(a::Vector{Float64}, N, periodic)
    if periodic
        a[1]=a[N-5]; a[2]=a[N-4]; a[3]=a[N-3]
        a[N-2]=a[4]; a[N-1]=a[5]; a[N]=a[6]
    else  # outflow
        a[1]=a[4]; a[2]=a[4]; a[3]=a[4]
        a[N-2]=a[N-3]; a[N-1]=a[N-3]; a[N]=a[N-3]
    end
end
function _ghosts!(s::ConfState)
    N=length(s.x)
    for a in (s.ξ,s.u,s.ξx,s.ux,s.ξt,s.ut,s.Ttt,s.Ttx); _set_ghost!(a,N,s.periodic); end
end

# --- KT flux at face i+1/2 (COMP: 0 -> T^{tt} eq, flux T^{tx}; 1 -> T^{tx}, flux T^{xx}) ---
@inline _Tcomp(comp, fr, ξ,u,ξx,ux,ξt,ut) =
    comp==0 ? T_tx(fr,ξ,u,ξx,ux,ξt,ut) : T_xx(fr,ξ,u,ξx,ux,ξt,ut)
@inline _Tjump(comp, fr, ξ,u,ξx,ux,ξt,ut) =
    comp==0 ? T_tt(fr,ξ,u,ξx,ux,ξt,ut) : T_tx(fr,ξ,u,ξx,ux,ξt,ut)
function _flux(s::ConfState, i, comp)
    fr=s.fr
    ξL=_qLx(s.ξ,i);  ξR=_qRx(s.ξ,i)        # L = reconst_qLx, R = reconst_qRx (solver.c)
    uL=_qLx(s.u,i);  uR=_qRx(s.u,i)
    ξxL=_qLx(s.ξx,i);ξxR=_qRx(s.ξx,i)
    uxL=_qLx(s.ux,i);uxR=_qRx(s.ux,i)
    ξtL=_qLx(s.ξt,i);ξtR=_qRx(s.ξt,i)
    utL=_qLx(s.ut,i);utR=_qRx(s.ut,i)
    fL=_Tcomp(comp,fr,ξL,uL,ξxL,uxL,ξtL,utL)
    fR=_Tcomp(comp,fr,ξR,uR,ξxR,uxR,ξtR,utR)
    jump=_Tjump(comp,fr,ξR,uR,ξxR,uxR,ξtR,utR) - _Tjump(comp,fr,ξL,uL,ξxL,uxL,ξtL,utL)
    return 0.5*(fL + fR - 1.0*jump)   # KT central flux, a = c = 1
end
@inline _div(s,i,comp) = (_flux(s,i,comp) - _flux(s,i-1,comp))/s.dx

# recover primitive time-derivatives + spatial derivatives over the interior
function _update_aux!(s::ConfState)
    N=length(s.x)
    _ghosts!(s)
    for i in NG+1:N-NG
        s.ξx[i]=_Dx(s.ξ,i,s.dx,N); s.ux[i]=_Dx(s.u,i,s.dx,N)
    end
    _ghosts!(s)
    for i in NG+1:N-NG
        ξt,ut = recover_time_derivs(s.fr, s.ξ[i], s.u[i], s.ξx[i], s.ux[i], s.Ttt[i], s.Ttx[i])
        s.ξt[i]=ξt; s.ut[i]=ut
    end
    _ghosts!(s)
end

# --- Heun (SSP-RK2) step ----------------------------------------------------
function _step!(s::ConfState)
    N=length(s.x); dt=s.dt
    _update_aux!(s)
    Ttt0=copy(s.Ttt); Ttx0=copy(s.Ttx); ξ0=copy(s.ξ); u0=copy(s.u)
    dTtt1=zeros(N); dTtx1=zeros(N); dξ1=copy(s.ξt); du1=copy(s.ut)
    for i in NG+1:N-NG; dTtt1[i]=-_div(s,i,0); dTtx1[i]=-_div(s,i,1); end
    for i in NG+1:N-NG
        s.Ttt[i]=Ttt0[i]+dt*dTtt1[i]; s.Ttx[i]=Ttx0[i]+dt*dTtx1[i]
        s.ξ[i]=ξ0[i]+dt*dξ1[i];       s.u[i]=u0[i]+dt*du1[i]
    end
    _update_aux!(s)
    dTtt2=zeros(N); dTtx2=zeros(N); dξ2=copy(s.ξt); du2=copy(s.ut)
    for i in NG+1:N-NG; dTtt2[i]=-_div(s,i,0); dTtx2[i]=-_div(s,i,1); end
    for i in NG+1:N-NG
        s.Ttt[i]=Ttt0[i]+0.5dt*(dTtt1[i]+dTtt2[i]); s.Ttx[i]=Ttx0[i]+0.5dt*(dTtx1[i]+dTtx2[i])
        s.ξ[i]=ξ0[i]+0.5dt*(dξ1[i]+dξ2[i]);         s.u[i]=u0[i]+0.5dt*(du1[i]+du2[i])
    end
    _update_aux!(s)
    return s
end

"""evolve!(s, nsteps) advances `nsteps` Heun steps; returns s."""
function evolve!(s::ConfState, nsteps::Int)
    for _ in 1:nsteps; _step!(s); end
    return s
end

# --- initial data -----------------------------------------------------------
function _make_state(fr, x, dx, cfl, periodic)
    N=length(x); dt=cfl*dx
    z=zeros(N)
    ConfState(fr, x, dx, dt, copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z), periodic)
end

"""Gaussian energy clump: ε = A exp(-(x-x0)²/w) + c, u=0."""
function init_gaussian(fr::ConformalFrame; N=257, xmin=-200.0, xmax=200.0,
                       A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)
    x=collect(range(xmin,xmax;length=N)); dx=x[2]-x[1]
    s=_make_state(fr,x,dx,cfl,false)
    for i in 1:N
        ε=A*exp(-(x[i]-x0)^2/w)+c; s.ξ[i]=log(ε); s.u[i]=0.0
        s.Ttt[i]=T_tt(fr,s.ξ[i],0.0,0,0,0,0); s.Ttx[i]=T_tx(fr,s.ξ[i],0.0,0,0,0,0)
    end
    _update_aux!(s); return s
end

"""Steady planar shock (conformal RH) in its rest frame: a stationary solution."""
function init_smooth_shock(fr::ConformalFrame; N=257, xmin=-200.0, xmax=200.0,
                           εL=1.0, vL=0.8, cfl=0.1)
    εR,vR = rankine_hugoniot(εL,vL)
    x=collect(range(xmin,xmax;length=N)); dx=x[2]-x[1]
    s=_make_state(fr,x,dx,cfl,false)
    for i in 1:N
        ε=(εR-εL)/2*(erf_(x[i]/10)+1)+εL
        v=(vL-vR)/2*(1-erf_(x[i]/10))+vR
        u=v/sqrt(1-v^2)
        s.ξ[i]=log(ε); s.u[i]=u
        s.Ttt[i]=T_tt(fr,s.ξ[i],u,0,0,0,0); s.Ttx[i]=T_tx(fr,s.ξ[i],u,0,0,0,0)
    end
    _update_aux!(s); return s
end

# erf without SpecialFunctions (Abramowitz–Stegun 7.1.26, |err|<1.5e-7)
function erf_(x::Real)
    t=1/(1+0.3275911*abs(x))
    y=1-(((((1.061405429t-1.453152027)t)+1.421413741)t-0.284496736)t+0.254829592)t*exp(-x^2)
    return x≥0 ? y : -y
end

end # module ConformalEvolution
