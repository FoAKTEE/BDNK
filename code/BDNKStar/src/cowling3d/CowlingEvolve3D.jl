#=
    CowlingEvolve3D — STAGE-3 Phase-2 engine: 3+1D time-domain evolution of the
    LINEARIZED relativistic fluid perturbations on the frozen TOV background
    (Cowling approximation), on a Cartesian grid. The non-radial spectrum is read
    off by FFT of an ℓ-projected diagnostic and validated against the Phase-1
    eigensolver (NonRadialModes).

    Linearized Cowling system (Gaertig–Kokkotas; static TOV, v₀=0, areal coords;
    α=lapse, γ_ij spatial metric, c_s²=dp/dε, Φ'=½ν'). Evolved variables are the
    energy-density perturbation δε and the COVARIANT momentum δS_i=(ε₀+p₀)δv_i —
    evolving δS_i (not δv^i) makes the factor (ε₀+p₀) CANCEL out of every RHS term,
    removing the 1/(ε₀+p₀) blow-up at the stellar surface (Font/Valencia momentum
    form). With δp=c_s²δε and contravariant δS^i = γ^{ij}δS_j:
        ∂_t δε  = −(1/√γ) ∂_i[ α √γ δS^i ]  −  α Φ' (n_i δS^i)
        ∂_t δS_i = −α[ ∂_i δp  +  Φ' n_i (δp + δε) ]
    (the lapse on the momentum RHS is derived in
     exact_solutions/mathematica/t175_cowling_source.m; the continuity
     equation above is exact as written)
    with γ^{ij} = δ_ij + (e^{-2λ}−1) n_i n_j  (inverse of γ_ij=δ_ij+(e^{2λ}−1)n_in_j).

    Method of lines: central 2nd-order finite differences + Kreiss–Oliger
    dissipation, RK4 in time. Perturbations are confined to the stellar interior
    (atmosphere frozen) — linear ⇒ the exact background never drifts.
=#
module CowlingEvolve3D

using ..Background3D
using ..Background3D: CartesianGrid, Star3D
using ..Units: Msun_to_km, kHz_to_km
using LinearAlgebra: svd, eigvals, pinv

export EvolState, Evo3D, setup_evo3d, seed_l2!, evolve3d!, l2_quadrupole, central_deps,
       periodogram, freq_kHz_cyclic, damping_rate,
       ylm_real, seed_ylm!, ylm_moment, evolve3d_moments!,
       spectral_peak, pencil_modes, envelope_ratio, analyze_qnm

# ---- precomputed background coefficients on the grid ----------------------
struct Evo3D
    s::Star3D
    w0::Array{Float64,3}      # ε₀+p₀
    Φp::Array{Float64,3}      # Φ' = (m+4πr³p)/(r(r−2m))
    p0p::Array{Float64,3}     # p₀' = −w0 Φ'
    qinv::Array{Float64,3}    # 1/e^{2λ} − 1  (inverse-metric radial correction)
    η0::Array{Float64,3}      # BDNK shear viscosity field
    ζ0::Array{Float64,3}      # BDNK bulk viscosity field
    viscous::Bool
    σ_ko::Float64
    cut::Bool                       # cut-cell (aperture-weighted) surface treatment
    κmin::Float64                   # small-cell floor on the volume fraction
    κ::Array{Float64,3}             # volume fraction inside r=R
    axm::Array{Float64,3}; axp::Array{Float64,3}   # face apertures, -x/+x
    aym::Array{Float64,3}; ayp::Array{Float64,3}
    azm::Array{Float64,3}; azp::Array{Float64,3}
end

# --- cut-cell geometry: both kappa and the apertures reduce to 1D quadratures of
# the SAME disc-rectangle area, so no interface plane-reconstruction is needed.
# EXACT TO MACHINE PRECISION: the integrands have kinks where the disc boundary meets a
# rectangle edge or corner (and a sqrt endpoint singularity at the disc rim), so a single
# Gauss rule converges only algebraically — 24 points left ~1e-5 quadrature error, and
# because kappa integrates over x while the faces integrate over y and z, that error was
# NOT symmetric under x<->y and leaked a ~1e-5 cross-moment between E_g partners. Here
# the rim singularity is removed by the substitution y = a sin(theta) (x = R sin(phi) for
# kappa) and every integral is split at its kink points, so each piece is smooth and
# 16-point Gauss is exact to round-off. The geometry then respects the sphere's full
# octahedral symmetry, and the parity cancellations of the Y_lm moments hold exactly.
# Verified against the analytic sphere: volume to 1e-12, and the closure condition
# (face-aperture area vectors balancing the boundary face) reproduces 4 pi R^2 at 2nd order.
function _gauss(ng)
    gx=zeros(ng); gw=zeros(ng)
    for i in 1:ng
        x=cos(pi*(i-0.25)/(ng+0.5))
        for _ in 1:100
            p0=1.0; p1=0.0
            for j in 1:ng; p2=p1; p1=p0; p0=((2j-1)*x*p1-(j-1)*p2)/j; end
            dp=ng*(x*p0-p1)/(x^2-1); dxx=-p0/dp; x+=dxx; abs(dxx)<1e-14 && break
        end
        p0=1.0; p1=0.0
        for j in 1:ng; p2=p1; p1=p0; p0=((2j-1)*x*p1-(j-1)*p2)/j; end
        dp=ng*(x*p0-p1)/(x^2-1); gx[i]=x; gw[i]=2/((1-x^2)*dp^2)
    end
    (gx,gw)
end
const _GX16, _GW16 = _gauss(16)

@inline function _gauss_int(f, a, b, gx, gw)
    b <= a && return 0.0
    c=(a+b)/2; hh=(b-a)/2; acc=0.0
    @inbounds for i in eachindex(gx); acc += gw[i]*f(c+hh*gx[i]); end
    acc*hh
end
# ∫_a^b f split at the breakpoints inside (a,b): exact for piecewise-smooth f
function _piecewise_int(f, a, b, bps, gx, gw)
    pts = sort!(filter(t -> a < t < b, bps)); pushfirst!(pts, a); push!(pts, b)
    acc=0.0
    for i in 1:length(pts)-1; acc += _gauss_int(f, pts[i], pts[i+1], gx, gw); end
    acc
end
# area of the disc of radius a (centred at the origin of the (y,z) plane) ∩ [y0,y1]×[z0,z1]
function _disc_rect(a, y0,y1, z0,z1, gx=_GX16, gw=_GW16)
    a <= 0 && return 0.0
    lo=max(y0,-a); hi=min(y1,a); hi<=lo && return 0.0
    θlo=asin(clamp(lo/a,-1.0,1.0)); θhi=asin(clamp(hi/a,-1.0,1.0))
    F(θ) = begin zh=a*cos(θ); w=min(z1,zh)-max(z0,-zh); w>0 ? w*zh : 0.0 end   # dy = a cosθ dθ
    bps=Float64[]
    for v in (z1, -z0, z0, -z1)                     # kinks where zh(θ) = v
        if 0 < v < a; t=acos(v/a); push!(bps, t); push!(bps, -t); end
    end
    _piecewise_int(F, θlo, θhi, bps, gx, gw)
end
# volume fraction of the sphere of radius R in the cell centred at (cx,cy,cz), half-width h
function _kappa_cell(R, cx,cy,cz, h, gx=_GX16, gw=_GW16)
    xlo=max(cx-h,-R); xhi=min(cx+h,R); xhi<=xlo && return 0.0
    φlo=asin(clamp(xlo/R,-1.0,1.0)); φhi=asin(clamp(xhi/R,-1.0,1.0))
    y0=cy-h; y1=cy+h; z0=cz-h; z1=cz+h
    G(φ) = begin a=R*cos(φ); _disc_rect(a, y0,y1,z0,z1, gx,gw)*a end           # dx = R cosφ dφ
    bps=Float64[]
    for v in (abs(y0),abs(y1),abs(z0),abs(z1), hypot(y0,z0),hypot(y0,z1),hypot(y1,z0),hypot(y1,z1))
        if 0 < v < R; t=acos(v/R); push!(bps, t); push!(bps, -t); end          # disc rim meets an edge/corner
    end
    _piecewise_int(G, φlo, φhi, bps, gx, gw) / (2h)^3
end
function _cutcell_geometry(s::Star3D)
    N=s.grid.N; dx=s.grid.dx; R=s.R; X=s.grid.x
    κ=zeros(N,N,N); axm=zeros(N,N,N); axp=zeros(N,N,N)
    aym=zeros(N,N,N); ayp=zeros(N,N,N); azm=zeros(N,N,N); azp=zeros(N,N,N)
    h=dx/2
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        cx=X[i]; cy=X[j]; cz=X[k]; r=sqrt(cx^2+cy^2+cz^2)
        if r > R + 2dx
            continue                                   # wholly outside
        elseif r < R - 2dx
            κ[i,j,k]=1.0; axm[i,j,k]=axp[i,j,k]=1.0
            aym[i,j,k]=ayp[i,j,k]=azm[i,j,k]=azp[i,j,k]=1.0
            continue
        end
        ax(xf)= _disc_rect(sqrt(max(0.0,R^2-xf^2)), cy-h,cy+h, cz-h,cz+h)/dx^2
        ay(yf)= _disc_rect(sqrt(max(0.0,R^2-yf^2)), cx-h,cx+h, cz-h,cz+h)/dx^2
        az(zf)= _disc_rect(sqrt(max(0.0,R^2-zf^2)), cx-h,cx+h, cy-h,cy+h)/dx^2
        axm[i,j,k]=ax(cx-h); axp[i,j,k]=ax(cx+h)
        aym[i,j,k]=ay(cy-h); ayp[i,j,k]=ay(cy+h)
        azm[i,j,k]=az(cz-h); azp[i,j,k]=az(cz+h)
        κ[i,j,k]=_kappa_cell(R, cx,cy,cz, h)
    end
    # The sphere is invariant under every permutation of the axes and the grid uses the
    # same coordinate vector on all three, so κ[i,j,k] is EXACTLY permutation-symmetric;
    # average over the six permutations to remove the residual ~1e-9 quadrature asymmetry
    # (the apertures already map into each other exactly).
    κs = similar(κ)
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        κs[i,j,k] = (κ[i,j,k]+κ[j,i,k]+κ[i,k,j]+κ[k,j,i]+κ[j,k,i]+κ[k,i,j])/6
    end
    κ = κs
    # Gauss quadrature of a full cell returns 1±few·eps; a volume fraction / aperture
    # outside [0,1] is meaningless and would leak into the κmin clipping, so pin it.
    for A in (κ,axm,axp,aym,ayp,azm,azp), I in eachindex(A)
        A[I] = clamp(A[I], 0.0, 1.0)
    end
    (κ,axm,axp,aym,ayp,azm,azp)
end

mutable struct EvolState
    δε::Array{Float64,3}
    Sx::Array{Float64,3}; Sy::Array{Float64,3}; Sz::Array{Float64,3}   # covariant δS_i
end

# preallocated RHS scratch (avoids per-step allocation in the long evolution)
struct Scratch
    δp::Array{Float64,3}
    Sux::Array{Float64,3}; Suy::Array{Float64,3}; Suz::Array{Float64,3}
    Fx::Array{Float64,3};  Fy::Array{Float64,3};  Fz::Array{Float64,3}
    # viscous: contravariant velocity + symmetric viscous stress π^{ij}
    vx::Array{Float64,3};  vy::Array{Float64,3};  vz::Array{Float64,3}
    πxx::Array{Float64,3}; πyy::Array{Float64,3}; πzz::Array{Float64,3}
    πxy::Array{Float64,3}; πxz::Array{Float64,3}; πyz::Array{Float64,3}
end
Scratch(N::Int) = Scratch((zeros(N,N,N) for _ in 1:16)...)

"""
    setup_evo3d(s; σ_ko=0.02, η̂=0.0, ζ̂=0.0, cut=true, κmin=0.5)

Precompute the evolution coefficients. `η̂`,`ζ̂` are dimensionless BDNK shear/bulk
viscosity knobs; the viscosity fields are η₀=η̂·ε₀, ζ₀=ζ̂·ε₀ (vanishing at the
surface). η̂=ζ̂=0 ⇒ the ideal (Phase-2) engine.

`cut=true` (DEFAULT) is the cut-cell free-surface treatment: the boundary face is
omitted from every aperture-weighted divergence, which imposes Δp=0 / zero mass flux /
zero traction at the true sphere r=R. `cut=false` is the legacy binary interior mask,
whose frozen exterior imposes δS_i=0 on the stair-stepped surface — a RIGID WALL. The
mask gives the ℓ=2 f-mode −6.95% at N=32 AND N=40 (it does not converge away) from a
run that looks perfectly healthy; the cut cells give −0.08% at N=48 (σ_ko=0.01). The
default was flipped to `true` for that reason; pass `cut=false` only for flat-space
tests where there is no surface (e.g. the plane-wave test).
"""
function setup_evo3d(s::Star3D; σ_ko::Float64=0.02, η̂::Float64=0.0, ζ̂::Float64=0.0,
                     cut::Bool=true, κmin::Float64=0.5)
    N = s.grid.N
    w0 = s.ε0 .+ s.p0
    Φp = similar(w0); p0p = similar(w0); qinv = similar(w0)
    @inbounds for I in eachindex(w0)
        r = s.r[I]; e2λ = s.e2λ[I]
        m = r * (1 - 1/e2λ) / 2                       # from e^{2λ}=(1−2m/r)^{-1}
        denom = r * (r - 2m)
        Φp[I]   = denom > 0 ? (m + 4π*r^3*s.p0[I]) / denom : 0.0
        p0p[I]  = -w0[I] * Φp[I]
        qinv[I] = 1/e2λ - 1
    end
    η0 = η̂ .* s.ε0; ζ0 = ζ̂ .* s.ε0
    z3()=zeros(N,N,N)
    κ,axm,axp,aym,ayp,azm,azp = cut ? _cutcell_geometry(s) :
                                (z3(),z3(),z3(),z3(),z3(),z3(),z3())
    Evo3D(s, w0, Φp, p0p, qinv, η0, ζ0, (η̂ > 0 || ζ̂ > 0), σ_ko,
          cut, κmin, κ, axm, axp, aym, ayp, azm, azp)
end

EvolState(N::Int) = EvolState(zeros(N,N,N), zeros(N,N,N), zeros(N,N,N), zeros(N,N,N))

# central first derivative along dim d (2nd order); zero outside [2,N-1]
@inline _dx(A,i,j,k,h) = (A[i+1,j,k]-A[i-1,j,k])/(2h)
@inline _dy(A,i,j,k,h) = (A[i,j+1,k]-A[i,j-1,k])/(2h)
@inline _dz(A,i,j,k,h) = (A[i,j,k+1]-A[i,j,k-1])/(2h)
# 2nd-order Kreiss–Oliger (5-point) dissipation per dim, returns the sum
@inline function _ko(A,i,j,k)
    (A[i-2,j,k]-4A[i-1,j,k]+6A[i,j,k]-4A[i+1,j,k]+A[i+2,j,k]) +
    (A[i,j-2,k]-4A[i,j-1,k]+6A[i,j,k]-4A[i,j+1,k]+A[i,j+2,k]) +
    (A[i,j,k-2]-4A[i,j,k-1]+6A[i,j,k]-4A[i,j,k+1]+A[i,j,k+2])
end

# RHS of the linearized momentum system into (dε,dSx,dSy,dSz)
function _rhs!(dε,dSx,dSy,dSz, st::EvolState, e::Evo3D, scr::Scratch)
    s = e.s; N = s.grid.N; h = s.grid.dx; σ = e.σ_ko
    δp, Sux, Suy, Suz = scr.δp, scr.Sux, scr.Suy, scr.Suz
    Fx, Fy, Fz = scr.Fx, scr.Fy, scr.Fz
    @inbounds for I in eachindex(δp)
        δp[I] = s.cs2[I] * st.δε[I]                       # closure δp = c_s² δε
        nS = s.nx[I]*st.Sx[I] + s.ny[I]*st.Sy[I] + s.nz[I]*st.Sz[I]   # n·δS
        Sux[I] = st.Sx[I] + e.qinv[I]*s.nx[I]*nS          # contravariant δS^i
        Suy[I] = st.Sy[I] + e.qinv[I]*s.ny[I]*nS
        Suz[I] = st.Sz[I] + e.qinv[I]*s.nz[I]*nS
        asg = s.α[I]*s.sqrtγ[I]                            # continuity flux α√γ δS^i
        Fx[I] = asg*Sux[I]; Fy[I] = asg*Suy[I]; Fz[I] = asg*Suz[I]
    end
    # --- BDNK viscous stress (LEADING / Navier–Stokes-limit shear+bulk) -------
    # π^{ij} ≡ −T^{ij}_visc = η(∂_iv^j+∂_jv^i−⅔δ_ijθ)+ζδ_ijθ : the overall sign is
    # absorbed so the momentum force below is +∂_jπ^{ij} = −∂_jT^{ij}_visc (the
    # dissipative RHS). SCOPE — this is the LEADING-order BDNK dissipation: a FLAT
    # coordinate divergence ∂_jπ^{ij} (no √γ/Christoffel weighting) of the shear+
    # bulk stress, WITHOUT the BDNK frame relaxation times (τ_ε,τ_P,τ_Q) or the
    # energy/heat sectors (A, Q^a). It is valid in the long-wavelength regime
    # where BDNK ≈ Navier–Stokes; the BDNK frame's role here is to regularize the
    # UV and render the first-order system causal — verified independently by the
    # Causality biquadratic monitor, NOT enforced by this (parabolic) evolution.
    # The full frame-recovery BDNK is the next refinement.
    if e.viscous
        vx, vy, vz = scr.vx, scr.vy, scr.vz
        # δv^i = δS^i/(ε₀+p₀). The floor MUST be relative: w₀ falls to the
        # atmosphere value ~1e-12·w₀max outside the star, so an absolute 1e-30 floor
        # gives δv ~ 1e12·δS there, and ∂δv in the viscous stress then blows up.
        # (Harmless while only masked interior cells were evolved; fatal once cut
        # cells evolve cells whose centre lies outside r=R.) No fluid ⇒ no velocity.
        wcut = 1e-6 * maximum(e.w0)
        @inbounds for I in eachindex(vx)
            if e.w0[I] < wcut
                vx[I] = 0.0; vy[I] = 0.0; vz[I] = 0.0
            else
                w = e.w0[I]
                vx[I] = Sux[I]/w; vy[I] = Suy[I]/w; vz[I] = Suz[I]/w
            end
        end
        πxx,πyy,πzz,πxy,πxz,πyz = scr.πxx,scr.πyy,scr.πzz,scr.πxy,scr.πxz,scr.πyz
        @inbounds for k in 2:N-1, j in 2:N-1, i in 2:N-1
            η = e.η0[i,j,k]; ζ = e.ζ0[i,j,k]
            if η == 0 && ζ == 0
                πxx[i,j,k]=πyy[i,j,k]=πzz[i,j,k]=πxy[i,j,k]=πxz[i,j,k]=πyz[i,j,k]=0.0
                continue
            end
            dxvx=_dx(vx,i,j,k,h); dyvy=_dy(vy,i,j,k,h); dzvz=_dz(vz,i,j,k,h)
            θ = dxvx+dyvy+dzvz
            dyvx=_dy(vx,i,j,k,h); dzvx=_dz(vx,i,j,k,h)
            dxvy=_dx(vy,i,j,k,h); dzvy=_dz(vy,i,j,k,h)
            dxvz=_dx(vz,i,j,k,h); dyvz=_dy(vz,i,j,k,h)
            b = (ζ - (2/3)*η)*θ
            πxx[i,j,k]=2η*dxvx+b; πyy[i,j,k]=2η*dyvy+b; πzz[i,j,k]=2η*dzvz+b
            πxy[i,j,k]=η*(dyvx+dxvy); πxz[i,j,k]=η*(dzvx+dxvz); πyz[i,j,k]=η*(dzvy+dyvz)
        end
    end
    fill!(dε,0.0); fill!(dSx,0.0); fill!(dSy,0.0); fill!(dSz,0.0)
    @inbounds for k in 3:N-2, j in 3:N-2, i in 3:N-2
        e.cut ? (e.κ[i,j,k] > 0.0 || continue) : (s.interior[i,j,k] || continue)
        α = s.α[i,j,k]; sg = s.sqrtγ[i,j,k]
        nx = s.nx[i,j,k]; ny = s.ny[i,j,k]; nz = s.nz[i,j,k]; Φp = e.Φp[i,j,k]
        # --- continuity:  ∂_t δε = −(1/√γ)∂_iF^i − α Φ' (n_i δS^i)
        # CUT CELL: aperture-weighted FINITE-VOLUME divergence. The boundary face
        # carries ZERO mass flux (free surface: the surface moves WITH the fluid),
        # so it simply does not appear in the sum -- that IS the BC, imposed at the
        # true r=R rather than at a stair-stepped wall. Face values are the 2nd-order
        # average of the two adjacent cell-centred fluxes.
        divF = if e.cut
            κc = max(e.κ[i,j,k], e.κmin)
            ( e.axp[i,j,k]*0.5*(Fx[i,j,k]+Fx[i+1,j,k]) - e.axm[i,j,k]*0.5*(Fx[i-1,j,k]+Fx[i,j,k])
            + e.ayp[i,j,k]*0.5*(Fy[i,j,k]+Fy[i,j+1,k]) - e.aym[i,j,k]*0.5*(Fy[i,j-1,k]+Fy[i,j,k])
            + e.azp[i,j,k]*0.5*(Fz[i,j,k]+Fz[i,j,k+1]) - e.azm[i,j,k]*0.5*(Fz[i,j,k-1]+Fz[i,j,k]) ) / (κc*h)
        else
            (Fx[i+1,j,k]-Fx[i-1,j,k] + Fy[i,j+1,k]-Fy[i,j-1,k] +
             Fz[i,j,k+1]-Fz[i,j,k-1]) / (2h)
        end
        ndotSup = nx*Sux[i,j,k] + ny*Suy[i,j,k] + nz*Suz[i,j,k]
        dε[i,j,k] = -divF/sg - α*Φp*ndotSup - σ*_ko(st.δε,i,j,k)
        # --- Euler:  ∂_t δS_i = −α[ ∂_iδp + Φ' n_i (δp+δε) ]
        #   The LAPSE multiplies the whole RHS. Derived from Div[δT^μ_ν]=0 via
        #     Div[T^μ_ν] = (1/√-g)∂_μ(√-g T^μ_ν) − ½T^{ab}∂_ν g_{ab},
        #   see exact_solutions/mathematica/t175_cowling_source.m: the flux is
        #   −(1/√γ)∂_i(α√γ δp) and the source is α[−δε Φ' + δp ∂_i ln√γ], and the
        #   ∂_i ln√γ pieces CANCEL between them, leaving exactly the line below.
        #   So it is NOT the α√γ-inside-the-divergence form (that leaves an
        #   uncancelled δp ∂ln√γ and overshoots); it is a bare α on the old RHS.
        #   Omitting it leaves the momentum and continuity equations with
        #   different time normalisations and puts a resolution-INDEPENDENT
        #   ~+5% into every mode frequency. The continuity equation above is
        #   already exact and is unchanged.
        # CUT CELL: INT d_j dp dV = OINT dp n_j dA. Coordinate faces contribute
        # aperture*dp; the boundary face contributes dp_b n_j A_b with dp_b = 0 --
        # the free-surface condition (p0'(R)=0 for a polytrope, so Delta p = dp).
        # Faces perpendicular to the other axes have n_j = 0 and drop out.
        local dpx, dpy, dpz
        if e.cut
            κc = max(e.κ[i,j,k], e.κmin)
            dpx = (e.axp[i,j,k]*0.5*(δp[i,j,k]+δp[i+1,j,k]) - e.axm[i,j,k]*0.5*(δp[i-1,j,k]+δp[i,j,k]))/(κc*h)
            dpy = (e.ayp[i,j,k]*0.5*(δp[i,j,k]+δp[i,j+1,k]) - e.aym[i,j,k]*0.5*(δp[i,j-1,k]+δp[i,j,k]))/(κc*h)
            dpz = (e.azp[i,j,k]*0.5*(δp[i,j,k]+δp[i,j,k+1]) - e.azm[i,j,k]*0.5*(δp[i,j,k-1]+δp[i,j,k]))/(κc*h)
        else
            dpx = _dx(δp,i,j,k,h); dpy = _dy(δp,i,j,k,h); dpz = _dz(δp,i,j,k,h)
        end
        src = Φp*(δp[i,j,k] + st.δε[i,j,k])
        dSx[i,j,k] = α*(-dpx - src*nx) - σ*_ko(st.Sx,i,j,k)
        dSy[i,j,k] = α*(-dpy - src*ny) - σ*_ko(st.Sy,i,j,k)
        dSz[i,j,k] = α*(-dpz - src*nz) - σ*_ko(st.Sz,i,j,k)
        # BDNK viscous force +∂_j π^{ij} (dissipative ⇒ damps the modes)
        if e.viscous
            πxx,πyy,πzz = scr.πxx,scr.πyy,scr.πzz; πxy,πxz,πyz = scr.πxy,scr.πxz,scr.πyz
            if e.cut
                # Aperture-weighted viscous divergence, same FV form as the pressure
                # gradient. The boundary face contributes π^{ij}n_j A_b, which is ZERO
                # for a free surface (no viscous traction on it) — so it drops out,
                # exactly as the mass flux and δp terms do.
                κc = max(e.κ[i,j,k], e.κmin); iκh = 1/(κc*h)
                dSx[i,j,k] += ( e.axp[i,j,k]*0.5*(πxx[i,j,k]+πxx[i+1,j,k]) - e.axm[i,j,k]*0.5*(πxx[i-1,j,k]+πxx[i,j,k])
                              + e.ayp[i,j,k]*0.5*(πxy[i,j,k]+πxy[i,j+1,k]) - e.aym[i,j,k]*0.5*(πxy[i,j-1,k]+πxy[i,j,k])
                              + e.azp[i,j,k]*0.5*(πxz[i,j,k]+πxz[i,j,k+1]) - e.azm[i,j,k]*0.5*(πxz[i,j,k-1]+πxz[i,j,k]) )*iκh
                dSy[i,j,k] += ( e.axp[i,j,k]*0.5*(πxy[i,j,k]+πxy[i+1,j,k]) - e.axm[i,j,k]*0.5*(πxy[i-1,j,k]+πxy[i,j,k])
                              + e.ayp[i,j,k]*0.5*(πyy[i,j,k]+πyy[i,j+1,k]) - e.aym[i,j,k]*0.5*(πyy[i,j-1,k]+πyy[i,j,k])
                              + e.azp[i,j,k]*0.5*(πyz[i,j,k]+πyz[i,j,k+1]) - e.azm[i,j,k]*0.5*(πyz[i,j,k-1]+πyz[i,j,k]) )*iκh
                dSz[i,j,k] += ( e.axp[i,j,k]*0.5*(πxz[i,j,k]+πxz[i+1,j,k]) - e.axm[i,j,k]*0.5*(πxz[i-1,j,k]+πxz[i,j,k])
                              + e.ayp[i,j,k]*0.5*(πyz[i,j,k]+πyz[i,j+1,k]) - e.aym[i,j,k]*0.5*(πyz[i,j-1,k]+πyz[i,j,k])
                              + e.azp[i,j,k]*0.5*(πzz[i,j,k]+πzz[i,j,k+1]) - e.azm[i,j,k]*0.5*(πzz[i,j,k-1]+πzz[i,j,k]) )*iκh
            else
            dSx[i,j,k] += (πxx[i+1,j,k]-πxx[i-1,j,k] + πxy[i,j+1,k]-πxy[i,j-1,k] + πxz[i,j,k+1]-πxz[i,j,k-1])/(2h)
            dSy[i,j,k] += (πxy[i+1,j,k]-πxy[i-1,j,k] + πyy[i,j+1,k]-πyy[i,j-1,k] + πyz[i,j,k+1]-πyz[i,j,k-1])/(2h)
            dSz[i,j,k] += (πxz[i+1,j,k]-πxz[i-1,j,k] + πyz[i,j+1,k]-πyz[i,j-1,k] + πzz[i,j,k+1]-πzz[i,j,k-1])/(2h)
            end
        end
    end
    return nothing
end

# one classical RK4 step (allocating scratch reused via closures kept simple)
function _rk4!(st::EvolState, e::Evo3D, dt::Float64, k, tmp, scr::Scratch)
    (k1ε,k1x,k1y,k1z, k2ε,k2x,k2y,k2z, k3ε,k3x,k3y,k3z, k4ε,k4x,k4y,k4z) = k
    (tε,tx,ty,tz) = tmp
    base = st
    _rhs!(k1ε,k1x,k1y,k1z, base, e, scr)
    @. tε = base.δε + dt/2*k1ε; @. tx = base.Sx + dt/2*k1x; @. ty = base.Sy + dt/2*k1y; @. tz = base.Sz + dt/2*k1z
    s2 = EvolState(tε,tx,ty,tz); _rhs!(k2ε,k2x,k2y,k2z, s2, e, scr)
    @. tε = base.δε + dt/2*k2ε; @. tx = base.Sx + dt/2*k2x; @. ty = base.Sy + dt/2*k2y; @. tz = base.Sz + dt/2*k2z
    _rhs!(k3ε,k3x,k3y,k3z, s2, e, scr)
    @. tε = base.δε + dt*k3ε; @. tx = base.Sx + dt*k3x; @. ty = base.Sy + dt*k3y; @. tz = base.Sz + dt*k3z
    _rhs!(k4ε,k4x,k4y,k4z, s2, e, scr)
    @. base.δε += dt/6*(k1ε+2k2ε+2k3ε+k4ε)
    @. base.Sx += dt/6*(k1x+2k2x+2k3x+k4x)
    @. base.Sy += dt/6*(k1y+2k2y+2k3y+k4y)
    @. base.Sz += dt/6*(k1z+2k2z+2k3z+k4z)
    return nothing
end

"""seed an ℓ=2, m=0 density perturbation δε = A ε₀ (r/R) (2z²−x²−y²)/r²."""
function seed_l2!(st::EvolState, e::Evo3D; A::Float64=1e-3)
    s = e.s; xs = s.grid.x; N = s.grid.N
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        s.interior[i,j,k] || continue
        X,Y,Z = xs[i],xs[j],xs[k]; r = s.r[i,j,k]
        st.δε[i,j,k] = A * s.ε0[i,j,k] * (r/s.R) * (2Z^2 - X^2 - Y^2)/r^2
    end
    return nothing
end

"""√γ-weighted ℓ=2,m=0 quadrupole moment of δε (the non-radial diagnostic)."""
function l2_quadrupole(st::EvolState, e::Evo3D)
    s = e.s; xs = s.grid.x; N = s.grid.N; acc = 0.0
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        s.interior[i,j,k] || continue
        X,Y,Z = xs[i],xs[j],xs[k]; r = s.r[i,j,k]
        acc += st.δε[i,j,k] * (2Z^2 - X^2 - Y^2)/r^2 * s.sqrtγ[i,j,k]
    end
    acc * s.grid.dx^3
end

central_deps(st::EvolState, e::Evo3D) = (c = e.s.grid.N÷2 + 1; st.δε[c,c,c])

"""
    evolve3d!(st, e; dt, nsteps, sample=1) -> (ts, q_l2, q_c)

Advance the linear Cowling perturbations `nsteps` RK4 steps of size `dt`,
recording the ℓ=2 quadrupole and central δε every `sample` steps.
"""
function evolve3d!(st::EvolState, e::Evo3D; dt::Float64, nsteps::Int, sample::Int=1)
    N = e.s.grid.N
    k   = ntuple(_ -> zeros(N,N,N), 16)
    tmp = ntuple(_ -> zeros(N,N,N), 4)
    scr = Scratch(N)
    ts = Float64[]; q2 = Float64[]; qc = Float64[]
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt); push!(q2, l2_quadrupole(st,e)); push!(qc, central_deps(st,e))
        end
        n == nsteps && break
        _rk4!(st, e, dt, k, tmp, scr)
    end
    return ts, q2, qc
end

# ---- spectral analysis (stdlib-only direct periodogram) -------------------
"""power spectrum |Σ q(t) e^{−iωt}|² of series `q(t)` at cyclic freqs `νs`."""
function periodogram(ts::Vector{Float64}, q::Vector{Float64}, νs::AbstractVector)
    q̄ = q .- sum(q)/length(q)
    P = similar(collect(float(νs)))
    @inbounds for (m,ν) in enumerate(νs)
        ω = 2π*ν; re = 0.0; im = 0.0
        for n in eachindex(ts)
            s,c = sincos(ω*ts[n]); re += q̄[n]*c; im -= q̄[n]*s
        end
        P[m] = re^2 + im^2
    end
    P
end

"""cyclic frequency ν (geometric, length-unit⁻¹) → kHz (M⊙ units by default)."""
freq_kHz_cyclic(ν::Real; Lunit_km::Real=Msun_to_km) = ν / (Lunit_km * kHz_to_km)

"""
    damping_rate(ts, q) -> 1/τ  (geometric units)

Mode damping rate from a decaying oscillation q(t): log-linear least-squares fit
of the local-maximum envelope. Returns ~0 for an undamped signal.

Returns **NaN** when the record holds fewer than 3 peaks, i.e. when the rate cannot
be MEASURED — the window is too short, or the mode is overdamped and no longer
rings. Deliberately not `0.0`: that is a physically meaningful rate (no damping),
so using it as a failure sentinel silently turns "could not measure" into
"measured zero". That collision made `test_bdnk_fullframe` compare a sentinel
against a real rate and conclude the opposite of the truth.

Callers must guard with `isfinite`. Comparisons against NaN are false, so an
unguarded `d2 > d0` now fails loudly instead of passing by accident.
"""
function damping_rate(ts::Vector{Float64}, q::Vector{Float64})
    a = abs.(q); tp = Float64[]; ap = Float64[]
    for n in 2:length(a)-1
        if a[n] > a[n-1] && a[n] ≥ a[n+1] && a[n] > 0
            push!(tp, ts[n]); push!(ap, a[n])
        end
    end
    length(ap) < 3 && return NaN
    x = tp; y = log.(ap); n = length(x)
    sx = sum(x); sy = sum(y); sxx = sum(abs2, x); sxy = sum(x .* y)
    slope = (n*sxy - sx*sy) / (n*sxx - sx^2)
    return -slope
end

# ═══════════════════════════════════════════════════════════════════════════════════
#  GENERAL (ℓ,m) SEEDING AND MOMENTS — full cube, no octant restriction, so every m
#  (including the T₂g partners m=±1 that an octant grid cannot carry) is admissible.
# ═══════════════════════════════════════════════════════════════════════════════════
"""
    ylm_real(l, m, nx, ny, nz) -> Float64

Real solid-harmonic ANGULAR factor Y_ℓm(n̂) on the unit vector n̂=(nx,ny,nz), ℓ∈{2,3,4},
m∈{−ℓ..ℓ}: m≥0 is the cosine-type harmonic (∝cos mφ), m<0 the sine-type (∝sin|m|φ).
Unnormalised Cartesian polynomials (the normalisation is irrelevant for a frequency
projection). ℓ=2, m=0 is exactly the (2n_z²−n_x²−n_y²) used by `seed_l2!`.
Under the three grid reflections x→−x, y→−y, z→−z the ℓ=2 quintet splits into the
cubic-group irreps E_g = {Y₂₀, Re Y₂₂} (even under all three) and T₂g = {Re Y₂₁,
Im Y₂₁, Im Y₂₂} (odd under two) — the full cube carries both.
"""
function ylm_real(l::Int, m::Int, nx::Float64, ny::Float64, nz::Float64)
    x2 = nx*nx; y2 = ny*ny; z2 = nz*nz; s2 = x2 + y2
    if l == 2
        m ==  0 && return 2z2 - s2
        m ==  1 && return nx*nz
        m == -1 && return ny*nz
        m ==  2 && return x2 - y2
        m == -2 && return 2nx*ny
    elseif l == 3
        m ==  0 && return nz*(2z2 - 3s2)
        m ==  1 && return nx*(4z2 - s2)
        m == -1 && return ny*(4z2 - s2)
        m ==  2 && return nz*(x2 - y2)
        m == -2 && return 2nx*ny*nz
        m ==  3 && return nx*(x2 - 3y2)
        m == -3 && return ny*(3x2 - y2)
    elseif l == 4
        m ==  0 && return 8z2*z2 - 24z2*s2 + 3s2*s2
        m ==  1 && return nx*nz*(4z2 - 3s2)
        m == -1 && return ny*nz*(4z2 - 3s2)
        m ==  2 && return (x2 - y2)*(6z2 - s2)
        m == -2 && return 2nx*ny*(6z2 - s2)
        m ==  3 && return nx*nz*(x2 - 3y2)
        m == -3 && return ny*nz*(3x2 - y2)
        m ==  4 && return x2*x2 - 6x2*y2 + y2*y2
        m == -4 && return 4nx*ny*(x2 - y2)
    end
    throw(ArgumentError("ylm_real: unsupported (l,m)=($l,$m); l∈{2,3,4}, |m|≤l"))
end

# evolved region: cut cells with positive volume fraction, else the interior mask
@inline _active(e::Evo3D, i, j, k) = e.cut ? (e.κ[i,j,k] > 0.0) : e.s.interior[i,j,k]

"""
    seed_ylm!(st, e; l=2, m=0, A=1e-3)

Seed δε = A ε₀ (r/R)^ℓ Y_ℓm(n̂) — the regular solid harmonic r^ℓY_ℓm, a polynomial in
(x,y,z), so the seed is smooth at the centre. NOTE the radial profile is (r/R)^ℓ, not the
(r/R)^1 of the legacy `seed_l2!`; both excite the same modes, this one with more f-mode
and less overtone content. δS_i = 0 (density seed, fluid at rest).
"""
function seed_ylm!(st::EvolState, e::Evo3D; l::Int=2, m::Int=0, A::Float64=1e-3)
    s = e.s; N = s.grid.N
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        _active(e, i, j, k) || continue
        r = s.r[i,j,k]
        st.δε[i,j,k] = A * s.ε0[i,j,k] * (r/s.R)^l * ylm_real(l, m, s.nx[i,j,k], s.ny[i,j,k], s.nz[i,j,k])
    end
    fill!(st.Sx, 0.0); fill!(st.Sy, 0.0); fill!(st.Sz, 0.0)
    return nothing
end

"""√γ-weighted projection ∫ δε Y_ℓm(n̂) √γ dV over the interior — the (ℓ,m) diagnostic.
For (2,0) this is `l2_quadrupole` exactly."""
function ylm_moment(δε::Array{Float64,3}, s::Star3D; l::Int=2, m::Int=0)
    N = s.grid.N; acc = 0.0
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        s.interior[i,j,k] || continue
        acc += δε[i,j,k] * ylm_real(l, m, s.nx[i,j,k], s.ny[i,j,k], s.nz[i,j,k]) * s.sqrtγ[i,j,k]
    end
    acc * s.grid.dx^3
end
ylm_moment(st::EvolState, e::Evo3D; l::Int=2, m::Int=0) = ylm_moment(st.δε, e.s; l=l, m=m)

"""
    evolve3d_moments!(st, e; dt, nsteps, sample=1, moments=[(2,0)]) -> (ts, Q)

As `evolve3d!`, but records an arbitrary list of (ℓ,m) moments; `Q[n, j]` is moment
`moments[j]` at sample `n`. Use several moments in one run to test the m-degeneracy
(or the E_g/T₂g splitting) on ONE background with ONE realisation of every systematic.
"""
function evolve3d_moments!(st::EvolState, e::Evo3D; dt::Float64, nsteps::Int, sample::Int=1,
                           moments::Vector{Tuple{Int,Int}}=[(2,0)])
    N = e.s.grid.N
    k   = ntuple(_ -> zeros(N,N,N), 16)
    tmp = ntuple(_ -> zeros(N,N,N), 4)
    scr = Scratch(N)
    ts = Float64[]; rows = Vector{Vector{Float64}}()
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt)
            push!(rows, [ylm_moment(st, e; l=l, m=m) for (l,m) in moments])
        end
        n == nsteps && break
        _rk4!(st, e, dt, k, tmp, scr)
    end
    Q = Matrix{Float64}(undef, length(ts), length(moments))
    for (n, r) in enumerate(rows); Q[n, :] .= r; end
    return ts, Q
end

# ═══════════════════════════════════════════════════════════════════════════════════
#  FREQUENCY EXTRACTION — two independent estimators that must agree, plus the gate.
#  A record of length T resolves Δν = 1/T and NOTHING finer; everything below Δν is a
#  parametric estimate and is reported as such (see `analyze_qnm`).
# ═══════════════════════════════════════════════════════════════════════════════════
"""
    spectral_peak(ts, q, νs) -> ν_peak

Periodogram peak on the cyclic-frequency grid `νs`, refined by a log-parabolic fit to
the three points around the maximum (a Gaussian-peak interpolation). Geometric units.
"""
function spectral_peak(ts::Vector{Float64}, q::Vector{Float64}, νs::AbstractVector)
    P = periodogram(ts, q, νs); i = argmax(P)
    (1 < i < length(νs)) || return float(νs[i])
    y1, y2, y3 = log(P[i-1] + 1e-300), log(P[i] + 1e-300), log(P[i+1] + 1e-300)
    den = y1 - 2y2 + y3
    d = den < 0 ? (y1 - y3)/(2den) : 0.0
    float(νs[i]) + clamp(d, -1.0, 1.0)*(νs[2] - νs[1])
end

"""
    pencil_modes(ts, q; npoles=14, nmax=600) -> Vector{NamedTuple}

Matrix-pencil (Hankel-SVD) pole estimator. The record is decimated uniformly to ≤`nmax`
samples, a Hankel matrix with pencil parameter L=M÷3 is formed, its SVD truncated to
`npoles`, and the poles z_k are the eigenvalues of pinv(V₁)V₂ built from the leading
right singular vectors. Each pole z=exp(s·Δt) gives (ν=|Im s|/2π, γ=−Re s) — decaying
modes have γ>0. Returns the poles with Im s>0, UNRANKED: residue ranking via a
Vandermonde over hundreds of samples under/overflows for |z|≠1 (a trap already hit in
this project), so the caller selects the pole nearest an independent periodogram peak.
"""
function pencil_modes(ts::Vector{Float64}, q::Vector{Float64}; npoles::Int=14, nmax::Int=600)
    n = length(q); step = max(1, cld(n, nmax))
    y = q[1:step:end]; y = y .- sum(y)/length(y)
    Δt = ts[1+step] - ts[1]
    M = length(y); L = M ÷ 3
    (M - L ≥ 2 && L ≥ 2) || return NamedTuple{(:ν,:γ),Tuple{Float64,Float64}}[]
    Y = Matrix{Float64}(undef, M-L, L+1)
    @inbounds for j in 1:L+1, i in 1:M-L; Y[i,j] = y[i+j-1]; end
    F = svd(Y); p = min(npoles, count(>(1e-12*F.S[1]), F.S))
    V = F.V[:, 1:p]
    z = eigvals(pinv(V[1:end-1, :]) * V[2:end, :])
    out = NamedTuple{(:ν,:γ),Tuple{Float64,Float64}}[]
    for zk in z
        sk = log(Complex(zk)) / Δt
        imag(sk) > 0 || continue
        push!(out, (ν = imag(sk)/(2π), γ = -real(sk)))
    end
    out
end

"""
    envelope_ratio(ts, q; npk=3) -> (ratio, npeaks)

Ratio of the mean of the last `npk` local maxima of |q| to the mean of the first `npk`.
>1 ⇒ the record is GROWING and must not be Fourier-analysed as a mode; this is the gate
that a naive FFT of e.g. an APR4 N=48 run (amplitude up 8e14) would have failed.
"""
function envelope_ratio(ts::Vector{Float64}, q::Vector{Float64}; npk::Int=3)
    a = abs.(q); pk = Float64[]
    for n in 2:length(a)-1
        (a[n] > a[n-1] && a[n] ≥ a[n+1] && a[n] > 0) && push!(pk, a[n])
    end
    length(pk) < 2npk && return (NaN, length(pk))
    (sum(pk[end-npk+1:end]) / sum(pk[1:npk]), length(pk))
end

"""
    analyze_qnm(ts, q; ν_ref, Lunit_km=Msun_to_km, band=(0.3, 3.0), nν=3000) -> NamedTuple

Gated two-estimator QNM read-out of one mode record `q(t)`:
  * `stable`   — finite, and envelope end/start ratio < 1
  * `f_kHz`    — periodogram peak (log-parabolic refined) in a band around `ν_ref`
  * `f_pencil_kHz`, `γ_pencil` — the matrix-pencil pole nearest that peak
  * `γ`        — envelope log-linear damping rate (geometric, amplitude e-folding)
  * `df_kHz`   — the Rayleigh resolution 1/T, which bounds what the record RESOLVES
  * `nperiods` — record length in periods of the measured mode
  * `agree_pct`— |f_pencil − f|/f in percent; the two estimators must agree for a claim
`ν_ref` is a cyclic frequency estimate in geometric units used only to set the search
band. If `stable` is false the frequencies are still returned but should not be quoted.
"""
function analyze_qnm(ts::Vector{Float64}, q::Vector{Float64}; ν_ref::Float64,
                     Lunit_km::Real=Msun_to_km, band::Tuple{Float64,Float64}=(0.3, 3.0), nν::Int=3000)
    T = ts[end] - ts[1]
    finite = all(isfinite, q)
    ratio, npk = finite ? envelope_ratio(ts, q) : (NaN, 0)
    stable = finite && isfinite(ratio) && ratio < 1.0
    νs = range(band[1]*ν_ref, band[2]*ν_ref; length=nν)
    ν = finite ? spectral_peak(ts, q, νs) : NaN
    poles = finite ? pencil_modes(ts, q) : NamedTuple{(:ν,:γ),Tuple{Float64,Float64}}[]
    νp, γp = if isempty(poles) || !isfinite(ν)
        (NaN, NaN)
    else
        j = argmin([abs(pl.ν - ν) for pl in poles]); (poles[j].ν, poles[j].γ)
    end
    γ = finite ? damping_rate(ts, q) : NaN
    (stable = stable, finite = finite, envelope_ratio = ratio, npeaks = npk,
     f_kHz = freq_kHz_cyclic(ν; Lunit_km=Lunit_km),
     f_pencil_kHz = freq_kHz_cyclic(νp; Lunit_km=Lunit_km), γ_pencil = γp, γ = γ,
     df_kHz = freq_kHz_cyclic(1/T; Lunit_km=Lunit_km), nperiods = T*ν,
     agree_pct = isfinite(νp) ? 100*abs(νp - ν)/ν : NaN)
end

end # module CowlingEvolve3D
