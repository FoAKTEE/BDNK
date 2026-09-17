#=
    DGBall3D — 3+1D nodal RKDG GRHydro for the Cowling TOV star on a CUBED-SPHERE grid that
    conforms to the star, after Hébert, Kidder & Teukolsky, PRD 98, 044041 (2018), Secs. II–IV
    and VI.B, Appendices A–B.

    Grid (their Fig. 14 / Table II, rescaled to the star's areal radius R):
      • a filled ball: a ROUNDED central cube (their eqs. A7–A10) surrounded by six-wedge
        cubed-sphere shells (eqs. A1–A6, no equiangular map) whose inner surface conforms to
        the cube (curvature c rises from 0.55 through 0.85 to 1 = spherical);
      • spherical shells inside the star (radial degree p_int), thin QUADRATIC (or linear)
        shells across the surface with the stellar surface EXACTLY on a shell boundary,
        larger shells outside to 3R; each wedge split n_t×n_t tangentially (degree p_t).
    Scheme:
      • strong-form nodal DG on mapped hexahedra (their eq. 12 / Teukolsky 2016): the
        physical flux divergence is formed by the chain rule, ∂_a f^a = (∂x̄^b/∂x^a) D^b f^a,
        with the discrete Jacobian (D applied to the nodal coordinates), so a uniform state
        is preserved EXACTLY; surface terms carry the coordinate area vector of each face
        node; conforming faces are matched geometrically (face centroids, then node
        positions), so wedge rotations need no bookkeeping;
      • Valencia GRHD on the frozen TOV metric in Cartesian areal coordinates,
        γ_ab = δ_ab + (e^λ−1) n_a n_b, α = e^Φ, with the COMPLETE sources
            s_{S_i} = (α/2)√γ T^{ab}∂_iγ_ab − √γ(ρhW²−p)∂_iα,   s_τ = −√γ ρhW² v^a ∂_aα,
        T^{ab}∂_iγ_ab = ρhW²[e^λλ'(v·n)² n_i + (2(e^λ−1)/r)((v·n)v^i − (v·n)²n_i)] + pλ' n_i,
        which reduce to TOV for the static star (DGCart3D omits the metric-derivative term
        and the √γ on the gravity term; its static residual is not TOV);
      • Γ-law ideal-gas EOS for the evolution (the polytrope is initial data and the lower
        entropy bound), atmosphere fixing after Galeazzi et al. 2013 (ρ_atm = 10⁻¹³ρ_c,
        κρ^(Γ−1)/(Γ−1) ≤ ε ≤ 100×), HLL (Davis) or LLF flux, SSP-RK3;
      • limiter :wb (the equilibrium-preserving positivity scaling of DGCart3D /
        VALIDATION.md §7.8, J-weighted means) on the surface shells, :none, or the
        mean-based :mean; optional exponential modal filter on the deviation from
        equilibrium; optional static-residual subtraction (exact fixed point; a fake force
        of the size of the raw residual once the surface moves — VALIDATION.md §7.9).
    Diagnostics: err[D̃], ρ_c (origin), M_b, the real spherical-harmonic moments
    q_ℓm = ∫ D̃ r^ℓ Y_ℓm d³x for (ℓ,m) = (0,0),(2,0),(2,2), max atmosphere speed.

    STATUS (VALIDATION.md §7.10, repro/dgball3d_hkt.jl, 2026-09-17). Coarsest grid
    (nt=2, p=3; 464 elements, 27 776 nodes, 0.03 s/step on 8 threads): exact free
    stream; the projected star is a fixed point to 5e-16 with the subtraction and
    settles to err[D̃] 5e-4 without it (the paper's B1: 6e-4); with the paper's
    momentum filter the seeded runs are stable for 2000 M⊙ and give F = 2.674 kHz
    (−0.45%) and f = 1.845/1.888/1.912 kHz (three estimators; −2..+1.5%, Q≈3).
    Without the filter every seeded run grows an instability with e-fold 250 M⊙.
    REFINEMENT DOES NOT YET PAY: at nt=3 (or nt=2, p_t=5) a spurious quadrupolar
    grid mode grows after the seeded transient (e-fold 190–310 M⊙) and saturates
    at q20/q00 ≈ 7%, with or without the subtraction; the static momentum residual
    is 2–5% of gravity (the geometric error of a cubic interpolant of a 45° patch,
    0.4% at p_t=5). The chain-rule strong form is the suspect; a split-form or
    over-integrated volume term is the next step.
=#
module DGBall3D

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon: rho_from_p
using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis
using ..Units: Msun_to_km, kHz_to_km
using LinearAlgebra: det, inv, norm, cross, dot

export DGBall3DEngine, DGBall3DState, ball_grid, setup_dgball3d, evolve_dgball3d!,
       seed_dgball3d_radial!, seed_dgball3d_l2!, dgball3d_central_density, dgball3d_errD,
       dgball3d_baryon_mass, dgball3d_moment, dgball3d_static_residual, dgball3d_volume,
       dgball3d_rhs_norm

# ----------------------------------------------------------------------------------
# geometry: mappings from the reference cube [−1,1]³ (their Appendix A)
# ----------------------------------------------------------------------------------
# wedge along +x: reference (ξ,η,ζ) → affine radial/tangential box → cubed-sphere → rotation
struct WedgeMap
    xmin::Float64; xmax::Float64; cmin::Float64; cmax::Float64
    ηlo::Float64; ηhi::Float64; ζlo::Float64; ζhi::Float64
    Rot::NTuple{9,Float64}                     # column-major 3×3 rotation taking +x̂ to the wedge axis
end
struct CubeMap
    xmin::Float64; cmin::Float64
    ξlo::Float64; ξhi::Float64; ηlo::Float64; ηhi::Float64; ζlo::Float64; ζhi::Float64
end
@inline _aff(lo,hi,s) = lo + (hi-lo)*(s+1)/2
function mapx(m::WedgeMap, ξ, η, ζ)
    X = _aff(m.xmin, m.xmax, ξ); Ȳ = _aff(m.ηlo, m.ηhi, η); Z̄ = _aff(m.ζlo, m.ζhi, ζ)
    a = 1/sqrt(1 + Ȳ^2 + Z̄^2)
    bmin = m.xmin*(1 + m.cmin*(a-1)); bmax = m.xmax*(1 + m.cmax*(a-1))
    ξr = bmin + (bmax-bmin)*(X-m.xmin)/(m.xmax-m.xmin)
    px, py, pz = ξr, ξr*Ȳ, ξr*Z̄
    R = m.Rot
    return (R[1]*px+R[4]*py+R[7]*pz, R[2]*px+R[5]*py+R[8]*pz, R[3]*px+R[6]*py+R[9]*pz)
end
function mapx(m::CubeMap, ξ, η, ζ)
    X = _aff(m.ξlo, m.ξhi, ξ); Y = _aff(m.ηlo, m.ηhi, η); Z = _aff(m.ζlo, m.ζhi, ζ)
    a = 1/sqrt(1 + X^2*Y^2 + X^2*Z^2 + Y^2*Z^2 - X^2*Y^2*Z^2)
    bmin = m.xmin*(1 + m.cmin*(a-1))
    return (bmin*X, bmin*Y, bmin*Z)
end
const _ROTS = (  # proper rotations taking (1,0,0) to ±x̂, ±ŷ, ±ẑ (column-major)
    (1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0),      # +x
    (-1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,1.0),    # −x  (rot_z π)
    (0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,0.0,1.0),     # +y  (rot_z π/2): x̂→ŷ
    (0.0,-1.0,0.0, 1.0,0.0,0.0, 0.0,0.0,1.0),     # −y
    (0.0,0.0,1.0, 0.0,1.0,0.0, -1.0,0.0,0.0),     # +z  (rot_y −π/2): x̂→ẑ
    (0.0,0.0,-1.0, 0.0,1.0,0.0, 1.0,0.0,0.0))     # −z

struct BallElem
    N::NTuple{3,Int}          # nodes per reference direction
    p::NTuple{3,Int}
    off::Int                  # flat index of node (1,1,1); node (i,j,k) ↦ off + (i-1) + N1(j-1) + N1N2(k-1)
    region::Symbol            # :cube | :center | :interior | :surface | :exterior
    map::Union{WedgeMap,CubeMap}
end
@inline nidx(e::BallElem, i, j, k) = e.off + (i-1) + e.N[1]*(j-1) + e.N[1]*e.N[2]*(k-1)
@inline nnodes(e::BallElem) = e.N[1]*e.N[2]*e.N[3]

struct BallFace
    elem::Int; side::Int                 # side: 1=−ξ 2=+ξ 3=−η 4=+η 5=−ζ 6=+ζ
    nodes::Vector{Int}                   # flat node indices on this face (own element)
    partner::Vector{Int}                 # matching node of the neighbour element (0 = outer boundary)
    A::Matrix{Float64}                   # 3×nf outward coordinate area vectors at the face nodes
    coef::Vector{Float64}                # 1/(J w_face) at the face nodes
end

"""
    ball_grid(R; nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3, surface=:quadratic, xmax_fac=3.0,
              n_ext=5) -> NamedTuple

Radial structure after Table II of Hébert–Kidder–Teukolsky (B1), rescaled from their
r_NS = 8.125 to the areal radius `R`, with the stellar surface on a shell boundary:
rounded cube to 1.3s (c=0.55), centre shells to 1.9s (c 0.85) and 2.5s (c 1), interior
spherical shells, `surface=:quadratic` → 5 shells of h=R/16 across [R−h, R+4h] with radial
degree 2 (the I2 configuration of the 1D study), `:linear` → 10 shells of h=R/32 across
[R−2h, R+8h] with degree 1 (their B1/I1), `n_ext` exterior shells to `xmax_fac`·R.
Each wedge is split `nt`×`nt` tangentially with degree `p_t`; the cube `nt`³ with degree `p_t`.
"""
function ball_grid(R::Float64; nt::Int=2, p_t::Int=3, p_int::Int=3, p_surf::Int=2, p_ext::Int=3,
                   surface::Symbol=:quadratic, xmax_fac::Float64=3.0, n_ext::Int=5)
    s = R/8.125
    if surface == :quadratic
        h = R/16; rin = R-h; rout = R+4h; nsurf = 5
    else
        h = R/32; rin = R-2h; rout = R+8h; nsurf = 10
    end
    center = [(1.3s, 1.9s, 0.55, 0.85), (1.9s, 2.5s, 0.85, 1.0)]
    ib = [2.5s, 3.0s, 3.6s, 4.33s, 5.2s, 6.24s, 0.5*(6.24s+rin), rin]
    interior = [(ib[i], ib[i+1]) for i in 1:length(ib)-1]
    surf = [(rin + (i-1)*h, rin + i*h) for i in 1:nsurf]
    eb = [rout + (xmax_fac*R - rout)*((i/n_ext)^1.3) for i in 0:n_ext]
    exterior = [(eb[i], eb[i+1]) for i in 1:n_ext]
    return (cube=(1.3s, 0.55), center=center, interior=interior, surface=surf, exterior=exterior,
            nt=nt, p_t=p_t, p_int=p_int, p_surf=p_surf, p_ext=p_ext)
end

function _build_elems(g)
    elems = BallElem[]; off = 1
    nt = g.nt; pt = g.p_t
    # rounded cube split nt³, degree p_t in every direction
    xmin_c, c_c = g.cube
    for k in 1:nt, j in 1:nt, i in 1:nt
        lo(i) = -1 + 2(i-1)/nt; hi(i) = -1 + 2i/nt
        m = CubeMap(xmin_c, c_c, lo(i), hi(i), lo(j), hi(j), lo(k), hi(k))
        push!(elems, BallElem((pt+1,pt+1,pt+1), (pt,pt,pt), off, :cube, m)); off += (pt+1)^3
    end
    # shells: (xmin, xmax, cmin, cmax, region, p_r)
    shells = Tuple{Float64,Float64,Float64,Float64,Symbol,Int}[]
    for (a,b,c1,c2) in g.center; push!(shells, (a,b,c1,c2,:center,g.p_int)); end
    for (a,b) in g.interior; push!(shells, (a,b,1.0,1.0,:interior,g.p_int)); end
    for (a,b) in g.surface;  push!(shells, (a,b,1.0,1.0,:surface,g.p_surf)); end
    for (a,b) in g.exterior; push!(shells, (a,b,1.0,1.0,:exterior,g.p_ext)); end
    for (xmin,xmax,cmin,cmax,reg,pr) in shells, w in 1:6, jt in 1:nt, it in 1:nt
        lo(i) = -1 + 2(i-1)/nt; hi(i) = -1 + 2i/nt
        m = WedgeMap(xmin, xmax, cmin, cmax, lo(it), hi(it), lo(jt), hi(jt), _ROTS[w])
        push!(elems, BallElem((pr+1,pt+1,pt+1), (pr,pt,pt), off, reg, m)); off += (pr+1)*(pt+1)^2
    end
    return elems, off-1
end

# ----------------------------------------------------------------------------------
# engine / state
# ----------------------------------------------------------------------------------
struct DGBall3DEngine
    elems::Vector{BallElem}
    faces::Vector{BallFace}
    bases::Dict{Int,LGLBasis}
    K::Int; Ntot::Int
    Γ::Float64; κ::Float64
    eos::BarotropicEOS
    R::Float64; M::Float64
    ρ_atm::Float64; ρ_cut::Float64; ε_atm::Float64; p_atm::Float64; εfac_max::Float64
    flux::Symbol; limiter::Symbol; limit_regions::Vector{Symbol}
    filter::Symbol; filt_vars::Symbol; filt_α::Float64; filt_s_center::Int; filt_s_shell::Int
    wellbalanced::Bool; entropy_floor::Bool
    cfl::Float64; dxmin::Float64
    # node geometry
    x::Vector{Float64}; y::Vector{Float64}; z::Vector{Float64}; r::Vector{Float64}
    J::Vector{Float64}                     # det ∂x/∂x̄
    Ji::Array{Float64,3}                   # Ji[b,a,i] = ∂x̄^b/∂x^a
    wq::Vector{Float64}                    # tensor LGL weight w_i w_j w_k (quadrature weight without J)
    α::Vector{Float64}; elam::Vector{Float64}; sqrtγ::Vector{Float64}
    nx::Vector{Float64}; ny::Vector{Float64}; nz::Vector{Float64}
    Φp::Vector{Float64}; λp::Vector{Float64}; gor::Vector{Float64}   # Φ', λ', (e^λ−1)/r
    origin::Vector{Int}                    # nodes at the origin
    # equilibrium and its static residual
    Deq::Vector{Float64}; Sxeq::Vector{Float64}; Syeq::Vector{Float64}; Szeq::Vector{Float64}; τeq::Vector{Float64}
    Seq::NTuple{5,Vector{Float64}}
end

mutable struct DGBall3DState
    D::Vector{Float64}; Sx::Vector{Float64}; Sy::Vector{Float64}; Sz::Vector{Float64}; τ::Vector{Float64}
    ρ::Vector{Float64}; ε::Vector{Float64}; p::Vector{Float64}
    vx::Vector{Float64}; vy::Vector{Float64}; vz::Vector{Float64}     # contravariant v^a
    nfix::Int
end
DGBall3DState(n) = DGBall3DState((zeros(n) for _ in 1:11)..., 0)

# ----------------------------------------------------------------------------------
# thermodynamics / conversions (Γ-law, covariant momentum, Cowling metric)
# ----------------------------------------------------------------------------------
@inline _p_ig(Γ,ρ,ε) = (Γ-1)*ρ*ε
@inline _εpoly(κ,Γ,ρ) = κ*ρ^(Γ-1)/(Γ-1)
@inline _cs2(Γ,ρ,ε) = (p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ; clamp(Γ*p/(ρ*h),0.0,1.0))
# lower an index: v_a = γ_ab v^b = v^a + (e^λ−1) n_a (n·v)
@inline function _lower(el,nx,ny,nz,vx,vy,vz)
    nv=nx*vx+ny*vy+nz*vz; q=el-1
    return vx+q*nx*nv, vy+q*ny*nv, vz+q*nz*nv
end
@inline function _raise(el,nx,ny,nz,wx,wy,wz)          # v^a = γ^{ab} w_b, γ^{ab} = δ − (1−e^{−λ}) n n
    nw=nx*wx+ny*wy+nz*wz; q=1/el-1
    return wx+q*nx*nw, wy+q*ny*nw, wz+q*nz*nw
end
# undensitized conserved from primitives
@inline function _p2c(Γ,ρ,ε,vx,vy,vz,el,nx,ny,nz)
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ
    vlx,vly,vlz=_lower(el,nx,ny,nz,vx,vy,vz)
    v2=clamp(vx*vlx+vy*vly+vz*vlz,0.0,1-1e-14); W=1/sqrt(1-v2)
    D̂=ρ*W; f=ρ*h*W^2
    return D̂, f*vlx, f*vly, f*vlz, f-p-D̂
end
# characteristic speeds along the Euclidean unit direction m (covector components m_a = m^a)
@inline function _speeds(Γ,ρ,ε,vx,vy,vz,el,nx,ny,nz,α,mx,my,mz)
    cs2=_cs2(Γ,ρ,ε)
    vlx,vly,vlz=_lower(el,nx,ny,nz,vx,vy,vz); v2=clamp(vx*vlx+vy*vly+vz*vlz,0.0,1-1e-14)
    vm=vx*mx+vy*my+vz*mz; nm=nx*mx+ny*my+nz*mz; γmm=1-(1-1/el)*nm^2
    a=α/(1-v2*cs2); disc=sqrt(max(cs2*(1-v2)*(γmm*(1-v2*cs2)-vm^2*(1-cs2)),0.0))
    return a*(vm*(1-cs2)-disc), a*(vm*(1-cs2)+disc)
end

# Γ-law cons2prim with Galeazzi-type fixing (see DGStarHP). Inputs undensitized covariant.
@inline function _c2p_resid(Γ,D̂,Ŝ,τ̂,p)
    v=Ŝ/(τ̂+D̂+p); W=1/sqrt(max(1-v^2,1e-16)); ρ=D̂/W
    ε=(τ̂+D̂*(1-W)+p*(1-W^2))/(D̂*W)
    return (Γ-1)*ρ*ε-p, ρ, ε, v, W
end
function _c2p(eng::DGBall3DEngine, D̂, Ŝx, Ŝy, Ŝz, τ̂in, el, nx, ny, nz, p0)
    Γ=eng.Γ; κ=eng.κ
    if !(D̂ > eng.ρ_cut) || !isfinite(D̂) || !isfinite(τ̂in) || !isfinite(Ŝx+Ŝy+Ŝz)
        return eng.ρ_atm, eng.ε_atm, 0.0,0.0,0.0, eng.p_atm, 1
    end
    code=0
    τ̂ = τ̂in < 0.0 ? 1e-300 : τ̂in; τ̂ != τ̂in && (code=2)
    ux,uy,uz=_raise(el,nx,ny,nz,Ŝx,Ŝy,Ŝz); S2=max(Ŝx*ux+Ŝy*uy+Ŝz*uz,0.0); Ŝ=sqrt(S2)
    Smax=sqrt(max(τ̂*(τ̂+2D̂),0.0))*(1-1e-10)
    scale=1.0
    if Ŝ > Smax; scale=Smax/Ŝ; Ŝ=Smax; code=2; end
    p=max(p0,1e-300); ok=false
    for it in 1:30
        f,_,_,_,_=_c2p_resid(Γ,D̂,Ŝ,τ̂,p); dp=max(1e-8*p,1e-300)
        f2,_,_,_,_=_c2p_resid(Γ,D̂,Ŝ,τ̂,p+dp); df=(f2-f)/dp
        df ≥ 0 && break
        pn=p-f/df; pn ≤ 0 && (pn=0.5*p)
        if abs(pn-p) ≤ 1e-13*max(pn,1e-300); p=pn; ok=true; break; end
        p=pn
    end
    if !ok
        plo=0.0; phi=max((Γ-1)*τ̂,1e-300)
        while _c2p_resid(Γ,D̂,Ŝ,τ̂,phi)[1] > 0 && phi < 1e10; phi*=2; end
        for _ in 1:200
            pm=0.5*(plo+phi); fm=_c2p_resid(Γ,D̂,Ŝ,τ̂,pm)[1]
            if fm > 0; plo=pm else phi=pm end
            (phi-plo) ≤ 1e-13*phi && break
        end
        p=0.5*(plo+phi)
    end
    _,ρ,ε,vmag,W=_c2p_resid(Γ,D̂,Ŝ,τ̂,p)
    if !(ρ > eng.ρ_atm)
        return eng.ρ_atm, eng.ε_atm, 0.0,0.0,0.0, eng.p_atm, 1
    end
    εmin= eng.entropy_floor ? _εpoly(κ,Γ,ρ) : 0.0; εmax=eng.εfac_max*_εpoly(κ,Γ,ρ)
    if ε < εmin; ε=εmin; code=2 elseif ε > εmax; ε=εmax; code=2 end
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ
    # covariant velocity v_a = Ŝ_a/(ρ h W²) (with the bound rescaling), then raise
    fac = S2 > 0 ? scale/(ρ*h*W^2) : 0.0
    vlx=Ŝx*fac; vly=Ŝy*fac; vlz=Ŝz*fac
    vx,vy,vz=_raise(el,nx,ny,nz,vlx,vly,vlz)
    return ρ,ε,vx,vy,vz,p,code
end

function _update_prims!(st::DGBall3DState, eng::DGBall3DEngine)
    Threads.@threads for i in 1:eng.Ntot
        @inbounds begin
            sg=eng.sqrtγ[i]; el=eng.elam[i]; nx=eng.nx[i]; ny=eng.ny[i]; nz=eng.nz[i]
            ρ,ε,vx,vy,vz,p,code=_c2p(eng, st.D[i]/sg, st.Sx[i]/sg, st.Sy[i]/sg, st.Sz[i]/sg, st.τ[i]/sg, el,nx,ny,nz, st.p[i])
            wasatm = st.ρ[i] ≤ eng.ρ_atm
            st.ρ[i]=ρ; st.ε[i]=ε; st.p[i]=p; st.vx[i]=vx; st.vy[i]=vy; st.vz[i]=vz
            if code != 0
                D̂,Ŝx,Ŝy,Ŝz,τ̂=_p2c(eng.Γ,ρ,ε,vx,vy,vz,el,nx,ny,nz)
                st.D[i]=sg*D̂; st.Sx[i]=sg*Ŝx; st.Sy[i]=sg*Ŝy; st.Sz[i]=sg*Ŝz; st.τ[i]=sg*τ̂
                (code==2 || !wasatm) && (st.nfix+=1)
            end
        end
    end
end

# ----------------------------------------------------------------------------------
# setup
# ----------------------------------------------------------------------------------
@inline function _lin(xs,ys,x)
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end

"""
    setup_dgball3d(eos, εc; grid=ball_grid(R), Γ=2.0, κ=eos.κ, atm_fac=1e-13, cut_fac=10.0,
                   εfac_max=100.0, flux=:hll, limiter=:wb, limit_regions=[:surface],
                   filter=:all, filt_vars=:momentum, filt_α=36.0, filt_s_center=6, filt_s_shell=12,
                   wellbalanced=true, entropy_floor=true, cfl=0.25,
                   h_tov=2e-4, grid_kwargs...) -> (engine, state)

Build the cubed-sphere DG star. `grid` is a `ball_grid` NamedTuple (or pass its keyword
arguments, e.g. `nt=2, p_int=3, surface=:quadratic`).

THE FILTER IS NOT OPTIONAL ON THIS GRID. Without it every seeded run develops an exponential
instability with e-folding time ≈ 250 M⊙ (err[D̃] and the central density drift, mass is
created, the star blows up at t ≈ 600–1900 M⊙ depending on the seed), the aliasing
instability of S̃_i that the paper reports "on O(100 M⊙) timescales" and cures the same way:
the exponential modal filter exp[−α(i/p)^s] on the MOMENTUM deviation only, α = 36, s = 6 in
the central cube/shells and s = 12 in the cubed-sphere shells, after every full step
(`filter=:all, filt_vars=:momentum`). Filtering all five variables (`filt_vars=:all`) destroys
the star instead. Measured (nt=2 grid, ℓ=2 seed 1e-3, t ∈ [600,800]): err[D̃] 6.0e-4 and
growing without the filter, 1.1e-5 and falling with it; the mode amplitude is 40% lower.
"""
function setup_dgball3d(eos::BarotropicEOS, εc::Float64; grid=nothing, Γ::Float64=2.0,
        κ::Float64=(eos isa ShumPolytrope ? eos.κ : 100.0), atm_fac::Float64=1e-13, cut_fac::Float64=10.0,
        εfac_max::Float64=100.0, flux::Symbol=:hll, limiter::Symbol=:wb, limit_regions::Vector{Symbol}=[:surface],
        filter::Symbol=:all, filt_vars::Symbol=:momentum, filt_α::Float64=36.0, filt_s_center::Int=6, filt_s_shell::Int=12,
        wellbalanced::Bool=true, entropy_floor::Bool=true,
        cfl::Float64=0.25, h_tov::Float64=2e-4, grid_kwargs...)
    star=solve_tov(eos,εc;h=h_tov); R,M=star.R,star.M
    g = grid === nothing ? ball_grid(R; grid_kwargs...) : grid
    elems,Ntot=_build_elems(g); K=length(elems)
    bases=Dict{Int,LGLBasis}()
    for e in elems, p in e.p; haskey(bases,p) || (bases[p]=build_lgl_basis(p)); end
    # ---- node coordinates
    x=zeros(Ntot); y=zeros(Ntot); z=zeros(Ntot)
    for e in elems
        b1=bases[e.p[1]]; b2=bases[e.p[2]]; b3=bases[e.p[3]]
        for k in 1:e.N[3], j in 1:e.N[2], i in 1:e.N[1]
            xx,yy,zz=mapx(e.map, b1.ξ[i], b2.ξ[j], b3.ξ[k]); n=nidx(e,i,j,k)
            x[n]=xx; y[n]=yy; z[n]=zz
        end
    end
    r=sqrt.(x.^2 .+ y.^2 .+ z.^2)
    # ---- discrete Jacobian ∂x^a/∂x̄^b (D along each reference line), inverse, determinant, weights
    Jm=zeros(3,3,Ntot); J=zeros(Ntot); Ji=zeros(3,3,Ntot); wq=zeros(Ntot); dxmin=Inf
    coords=(x,y,z)
    for e in elems
        bs=(bases[e.p[1]],bases[e.p[2]],bases[e.p[3]])
        for k in 1:e.N[3], j in 1:e.N[2], i in 1:e.N[1]
            n=nidx(e,i,j,k); wq[n]=bs[1].w[i]*bs[2].w[j]*bs[3].w[k]
            for a in 1:3
                # ∂x^a/∂ξ
                s=0.0; for l in 1:e.N[1]; s+=bs[1].D[i,l]*coords[a][nidx(e,l,j,k)]; end; Jm[a,1,n]=s
                s=0.0; for l in 1:e.N[2]; s+=bs[2].D[j,l]*coords[a][nidx(e,i,l,k)]; end; Jm[a,2,n]=s
                s=0.0; for l in 1:e.N[3]; s+=bs[3].D[k,l]*coords[a][nidx(e,i,j,l)]; end; Jm[a,3,n]=s
            end
            Jn=Jm[:,:,n]; J[n]=det(Jn); Ji[:,:,n]=inv(Jn)
            # smallest node spacing along the reference lines
            i<e.N[1] && (dxmin=min(dxmin, norm([coords[a][nidx(e,i+1,j,k)]-coords[a][n] for a in 1:3])))
            j<e.N[2] && (dxmin=min(dxmin, norm([coords[a][nidx(e,i,j+1,k)]-coords[a][n] for a in 1:3])))
            k<e.N[3] && (dxmin=min(dxmin, norm([coords[a][nidx(e,i,j,k+1)]-coords[a][n] for a in 1:3])))
        end
    end
    all(J .> 0) || error("negative Jacobian: $(count(J .≤ 0)) nodes (mapping orientation)")
    # ---- faces: area vectors, coefficients, geometric matching
    faces=_build_faces(elems,bases,x,y,z,Jm,J)
    # ---- Cowling metric at the nodes (areal Cartesian coordinates)
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p; et=star.ε
    m_of(rr)= rr<R ? _lin(rt,mt,rr) : M
    ν_of(rr)= rr<R ? _lin(rt,νt,rr) : log(1-2M/rr)
    p_of(rr)= rr<R ? max(_lin(rt,pt,rr),0.0) : 0.0
    e_of(rr)= rr<R ? max(_lin(rt,et,rr),0.0) : 0.0
    α=zeros(Ntot); elam=zeros(Ntot); sqrtγ=zeros(Ntot); nx=zeros(Ntot); ny=zeros(Ntot); nz=zeros(Ntot)
    Φp=zeros(Ntot); λp=zeros(Ntot); gor=zeros(Ntot); origin=Int[]
    for n in 1:Ntot
        rr=r[n]
        if rr < 1e-9
            α[n]=exp(0.5*ν_of(1e-6)); elam[n]=1.0; sqrtγ[n]=1.0; push!(origin,n); continue   # all source coefficients vanish at r=0
        end
        m=m_of(rr); f=max(1-2m/rr,1e-12)
        α[n]=exp(0.5*ν_of(rr)); elam[n]=1/f; sqrtγ[n]=sqrt(elam[n])
        nx[n]=x[n]/rr; ny[n]=y[n]/rr; nz[n]=z[n]/rr
        Φp[n]=(m+4π*rr^3*p_of(rr))/(rr*(rr-2m))
        mp=4π*rr^2*e_of(rr); λp[n]=(2mp/rr-2m/rr^2)/f
        gor[n]=(elam[n]-1)/rr
    end
    ρc=rho_from_p(eos,pressure(eos,εc)); ρ_atm=atm_fac*ρc; ε_atm=_εpoly(κ,Γ,ρ_atm); p_atm=_p_ig(Γ,ρ_atm,ε_atm)
    st=DGBall3DState(Ntot)
    for n in 1:Ntot
        rr=r[n]
        ρ = rr < R*(1-1e-10) ? max(rho_from_p(eos,p_of(rr)),ρ_atm) : ρ_atm
        ρ ≤ cut_fac*ρ_atm && (ρ=ρ_atm)
        ε=_εpoly(κ,Γ,ρ); st.ρ[n]=ρ; st.ε[n]=ε; st.p[n]=_p_ig(Γ,ρ,ε)
        D̂,Ŝx,Ŝy,Ŝz,τ̂=_p2c(Γ,ρ,ε,0.0,0.0,0.0,elam[n],nx[n],ny[n],nz[n]); sg=sqrtγ[n]
        st.D[n]=sg*D̂; st.Sx[n]=sg*Ŝx; st.Sy[n]=sg*Ŝy; st.Sz[n]=sg*Ŝz; st.τ[n]=sg*τ̂
    end
    Seq=(zeros(Ntot),zeros(Ntot),zeros(Ntot),zeros(Ntot),zeros(Ntot))
    eng=DGBall3DEngine(elems,faces,bases,K,Ntot,Γ,κ,eos,R,M,ρ_atm,cut_fac*ρ_atm,ε_atm,p_atm,εfac_max,
                       flux,limiter,limit_regions,filter,filt_vars,filt_α,filt_s_center,filt_s_shell,wellbalanced,entropy_floor,cfl,dxmin,
                       x,y,z,r,J,Ji,wq,α,elam,sqrtγ,nx,ny,nz,Φp,λp,gor,origin,
                       copy(st.D),copy(st.Sx),copy(st.Sy),copy(st.Sz),copy(st.τ),Seq)
    wellbalanced && _raw_rhs!(Seq..., st, eng)
    return eng, st
end

# face node lists per side, ordered; the outward area vector from the discrete Jacobian columns
function _face_nodes(e::BallElem, side::Int)
    N1,N2,N3=e.N; out=Int[]
    if side ≤ 2
        i = side==1 ? 1 : N1
        for k in 1:N3, j in 1:N2; push!(out, nidx(e,i,j,k)); end
    elseif side ≤ 4
        j = side==3 ? 1 : N2
        for k in 1:N3, i in 1:N1; push!(out, nidx(e,i,j,k)); end
    else
        k = side==5 ? 1 : N3
        for j in 1:N2, i in 1:N1; push!(out, nidx(e,i,j,k)); end
    end
    out
end
function _build_faces(elems, bases, x, y, z, Jm, J)
    faces=BallFace[]
    key(px,py,pz)=(round(Int64,px*1e8), round(Int64,py*1e8), round(Int64,pz*1e8))
    # face centroid → list of face ids
    cent=Dict{NTuple{3,Int64},Vector{Int}}()
    for (ke,e) in enumerate(elems), side in 1:6
        nodes=_face_nodes(e,side); nf=length(nodes)
        A=zeros(3,nf); coef=zeros(nf)
        dir=(side+1)÷2; sgn = isodd(side) ? -1.0 : 1.0
        t1,t2 = dir==1 ? (2,3) : (dir==2 ? (3,1) : (1,2))         # right-handed: t1×t2 ∥ +dir
        b=bases[e.p[dir]]; wf = isodd(side) ? b.w[1] : b.w[end]
        for (q,n) in enumerate(nodes)
            u=(Jm[1,t1,n],Jm[2,t1,n],Jm[3,t1,n]); v=(Jm[1,t2,n],Jm[2,t2,n],Jm[3,t2,n])
            c=(u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1])
            A[1,q]=sgn*c[1]; A[2,q]=sgn*c[2]; A[3,q]=sgn*c[3]
            coef[q]=1/(J[n]*wf)
        end
        cx=sum(x[nodes])/nf; cy=sum(y[nodes])/nf; cz=sum(z[nodes])/nf
        push!(faces, BallFace(ke, side, nodes, zeros(Int,nf), A, coef))
        push!(get!(cent, key(cx,cy,cz), Int[]), length(faces))
    end
    # partner faces: same centroid; then node ↔ node by position
    for (_,ids) in cent
        length(ids) ≤ 2 || error("face shared by $(length(ids)) elements — grid is not conforming")
        length(ids) == 1 && continue                     # outer boundary face
        fa,fb=faces[ids[1]],faces[ids[2]]
        pos=Dict{NTuple{3,Int64},Int}()
        for n in fb.nodes; pos[key(x[n],y[n],z[n])]=n; end
        for (q,n) in enumerate(fa.nodes)
            m=get(pos, key(x[n],y[n],z[n]), 0); m==0 && error("face node without partner at ($(x[n]),$(y[n]),$(z[n]))")
            fa.partner[q]=m
        end
        pos=Dict{NTuple{3,Int64},Int}()
        for n in fa.nodes; pos[key(x[n],y[n],z[n])]=n; end
        for (q,n) in enumerate(fb.nodes); fb.partner[q]=pos[key(x[n],y[n],z[n])]; end
    end
    return faces
end

# ----------------------------------------------------------------------------------
# RHS
# ----------------------------------------------------------------------------------
# densitized fluxes f^a (a=1..3) of the 5 conserved variables at node n: returns 3×5 tuple-of-tuples
@inline function _fluxes(eng::DGBall3DEngine, n::Int, ρ, ε, vx, vy, vz)
    Γ=eng.Γ; el=eng.elam[n]; nx=eng.nx[n]; ny=eng.ny[n]; nz=eng.nz[n]; A=eng.α[n]*eng.sqrtγ[n]
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ
    vlx,vly,vlz=_lower(el,nx,ny,nz,vx,vy,vz); v2=clamp(vx*vlx+vy*vly+vz*vlz,0.0,1-1e-14); W=1/sqrt(1-v2)
    D̂=ρ*W; f=ρ*h*W^2; Ŝx=f*vlx; Ŝy=f*vly; Ŝz=f*vlz; τ̂=f-p-D̂
    fx=(A*D̂*vx, A*(Ŝx*vx+p), A*Ŝy*vx, A*Ŝz*vx, A*(τ̂+p)*vx)
    fy=(A*D̂*vy, A*Ŝx*vy, A*(Ŝy*vy+p), A*Ŝz*vy, A*(τ̂+p)*vy)
    fz=(A*D̂*vz, A*Ŝx*vz, A*Ŝy*vz, A*(Ŝz*vz+p), A*(τ̂+p)*vz)
    q=(A/eng.α[n]*D̂, A/eng.α[n]*Ŝx, A/eng.α[n]*Ŝy, A/eng.α[n]*Ŝz, A/eng.α[n]*τ̂)   # densitized state √γ Û
    return fx,fy,fz,q
end

function _raw_rhs!(rD,rSx,rSy,rSz,rτ, st::DGBall3DState, eng::DGBall3DEngine)
    Ntot=eng.Ntot; Γ=eng.Γ
    # nodal fluxes (stored: 3 directions × 5 variables)
    F=eng_fluxbuf(eng)
    Threads.@threads for n in 1:Ntot
        @inbounds begin
            fx,fy,fz,_=_fluxes(eng,n,st.ρ[n],st.ε[n],st.vx[n],st.vy[n],st.vz[n])
            for c in 1:5; F[c,1,n]=fx[c]; F[c,2,n]=fy[c]; F[c,3,n]=fz[c]; end
        end
    end
    rhs=(rD,rSx,rSy,rSz,rτ)
    # volume term: −Σ_b Σ_a Ji[b,a] (D^b F^a)
    Threads.@threads for e in eng.elems
        @inbounds begin
            bs=(eng.bases[e.p[1]],eng.bases[e.p[2]],eng.bases[e.p[3]]); N1,N2,N3=e.N
            for k in 1:N3, j in 1:N2, i in 1:N1
                n=nidx(e,i,j,k)
                for c in 1:5
                    acc=0.0
                    for a in 1:3
                        d1=0.0; for l in 1:N1; d1+=bs[1].D[i,l]*F[c,a,nidx(e,l,j,k)]; end
                        d2=0.0; for l in 1:N2; d2+=bs[2].D[j,l]*F[c,a,nidx(e,i,l,k)]; end
                        d3=0.0; for l in 1:N3; d3+=bs[3].D[k,l]*F[c,a,nidx(e,i,j,l)]; end
                        acc+=eng.Ji[1,a,n]*d1+eng.Ji[2,a,n]*d2+eng.Ji[3,a,n]*d3
                    end
                    rhs[c][n]=-acc
                end
            end
        end
    end
    # face terms
    Threads.@threads for fc in eng.faces
        @inbounds for (q,n) in enumerate(fc.nodes)
            Ax=fc.A[1,q]; Ay=fc.A[2,q]; Az=fc.A[3,q]; Am=sqrt(Ax^2+Ay^2+Az^2); Am>0 || continue
            mx=Ax/Am; my=Ay/Am; mz=Az/Am
            m=fc.partner[q]
            fxL,fyL,fzL,qL=_fluxes(eng,n,st.ρ[n],st.ε[n],st.vx[n],st.vy[n],st.vz[n])
            λmL,λpL=_speeds(Γ,st.ρ[n],st.ε[n],st.vx[n],st.vy[n],st.vz[n],eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n],eng.α[n],mx,my,mz)
            if m==0      # outer boundary: atmosphere at rest (outflow)
                fxR,fyR,fzR,qR=_fluxes(eng,n,eng.ρ_atm,eng.ε_atm,0.0,0.0,0.0)
                λmR,λpR=_speeds(Γ,eng.ρ_atm,eng.ε_atm,0.0,0.0,0.0,eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n],eng.α[n],mx,my,mz)
            else
                fxR,fyR,fzR,qR=_fluxes(eng,n,st.ρ[m],st.ε[m],st.vx[m],st.vy[m],st.vz[m])   # same position: own metric
                λmR,λpR=_speeds(Γ,st.ρ[m],st.ε[m],st.vx[m],st.vy[m],st.vz[m],eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n],eng.α[n],mx,my,mz)
            end
            cmin=min(λmL,λmR,0.0); cmax=max(λpL,λpR,0.0)
            for c in 1:5
                fnL=fxL[c]*mx+fyL[c]*my+fzL[c]*mz; fnR=fxR[c]*mx+fyR[c]*my+fzR[c]*mz
                if eng.flux==:hll
                    fs = cmax-cmin ≤ 1e-300 ? 0.5*(fnL+fnR) : (cmax*fnL-cmin*fnR+cmin*cmax*(qR[c]-qL[c]))/(cmax-cmin)
                else
                    amax=max(abs(λmL),abs(λpL),abs(λmR),abs(λpR)); fs=0.5*(fnL+fnR)-0.5*amax*(qR[c]-qL[c])
                end
                rhs[c][n]-=fc.coef[q]*Am*(fs-fnL)
            end
        end
    end
    # sources
    Threads.@threads for n in 1:Ntot
        @inbounds begin
            ρ=st.ρ[n]; ε=st.ε[n]; p=st.p[n]; vx=st.vx[n]; vy=st.vy[n]; vz=st.vz[n]
            el=eng.elam[n]; nx=eng.nx[n]; ny=eng.ny[n]; nz=eng.nz[n]; α=eng.α[n]; sg=eng.sqrtγ[n]
            h=1+ε+p/ρ; vlx,vly,vlz=_lower(el,nx,ny,nz,vx,vy,vz); W2=1/(1-clamp(vx*vlx+vy*vly+vz*vlz,0.0,1-1e-14))
            f=ρ*h*W2; vn=vx*nx+vy*ny+vz*nz
            c1=f*el*eng.λp[n]*vn^2 + p*eng.λp[n] - 2*f*eng.gor[n]*vn^2      # coefficient of n_i in T^{ab}∂_iγ_ab
            c2=2*f*eng.gor[n]*vn                                          # coefficient of v^i
            g=(f-p)*α*eng.Φp[n]                                            # gravity: −√γ E ∂_iα = −√γ (ρhW²−p) α Φ' n_i
            rSx[n]+=sg*(0.5*α*(c1*nx+c2*vx) - g*nx)
            rSy[n]+=sg*(0.5*α*(c1*ny+c2*vy) - g*ny)
            rSz[n]+=sg*(0.5*α*(c1*nz+c2*vz) - g*nz)
            rτ[n]+= -sg*f*vn*α*eng.Φp[n]
        end
    end
    return nothing
end
# per-engine flux buffer (allocated once)
const _FLUXBUF = Dict{UInt,Array{Float64,3}}()
function eng_fluxbuf(eng::DGBall3DEngine)
    k=objectid(eng)
    haskey(_FLUXBUF,k) || (_FLUXBUF[k]=zeros(5,3,eng.Ntot))
    _FLUXBUF[k]
end

function _rhs!(rD,rSx,rSy,rSz,rτ, st::DGBall3DState, eng::DGBall3DEngine)
    _raw_rhs!(rD,rSx,rSy,rSz,rτ, st, eng)
    if eng.wellbalanced
        S=eng.Seq
        @inbounds Threads.@threads for n in 1:eng.Ntot
            rD[n]-=S[1][n]; rSx[n]-=S[2][n]; rSy[n]-=S[3][n]; rSz[n]-=S[4][n]; rτ[n]-=S[5][n]
        end
    end
end

# ----------------------------------------------------------------------------------
# limiters (equilibrium-preserving scaling, J-weighted means) and filter
# ----------------------------------------------------------------------------------
@inline _q3(D,Sx,Sy,Sz,τ,el,nx,ny,nz) = begin
    ux,uy,uz=_raise(el,nx,ny,nz,Sx,Sy,Sz); (τ+D) - sqrt(D^2 + max(Sx*ux+Sy*uy+Sz*uz,0.0))
end
function _elem_mean(eng::DGBall3DEngine, U, e::BallElem)
    s=0.0; w=0.0
    @inbounds for k in 1:e.N[3], j in 1:e.N[2], i in 1:e.N[1]
        n=nidx(e,i,j,k); ww=eng.wq[n]*eng.J[n]; s+=ww*U[n]; w+=ww
    end
    s/w
end
function _limit_wb_elem!(st::DGBall3DState, eng::DGBall3DEngine, e::BallElem)
    nn=nnodes(e); o=e.off-1
    # J-weighted means of the deviations and of the background
    mD=0.0; mSx=0.0; mSy=0.0; mSz=0.0; mτ=0.0; Deqm=0.0; τeqm=0.0; wSm=0.0; wt=0.0
    @inbounds for q in 1:nn
        n=o+q; ww=eng.wq[n]*eng.J[n]; wt+=ww
        mD+=ww*(st.D[n]-eng.Deq[n]); mSx+=ww*(st.Sx[n]-eng.Sxeq[n]); mSy+=ww*(st.Sy[n]-eng.Syeq[n]); mSz+=ww*(st.Sz[n]-eng.Szeq[n]); mτ+=ww*(st.τ[n]-eng.τeq[n])
        Deqm+=ww*eng.Deq[n]; τeqm+=ww*eng.τeq[n]; wSm+=ww*(eng.Deq[n]-eng.sqrtγ[n]*eng.ρ_atm)
    end
    mD/=wt; mSx/=wt; mSy/=wt; mSz/=wt; mτ/=wt; Deqm/=wt; τeqm/=wt; wSm/=wt
    RD=zeros(nn); RSx=zeros(nn); RSy=zeros(nn); RSz=zeros(nn); Rτ=zeros(nn); ΔK=zeros(nn); ΔKm=0.0
    @inbounds for q in 1:nn
        n=o+q; wD=eng.Deq[n]/Deqm; wsS = wSm>0 ? (eng.Deq[n]-eng.sqrtγ[n]*eng.ρ_atm)/wSm : wD
        RD[q]=eng.Deq[n]+mD*wD; RSx[q]=eng.Sxeq[n]+mSx*wsS; RSy[q]=eng.Syeq[n]+mSy*wsS; RSz[q]=eng.Szeq[n]+mSz*wsS
        el=eng.elam[n]; nx=eng.nx[n]; ny=eng.ny[n]; nz=eng.nz[n]
        ΔK[q]=_q3(eng.Deq[n],eng.Sxeq[n],eng.Syeq[n],eng.Szeq[n],0.0,el,nx,ny,nz)-_q3(RD[q],RSx[q],RSy[q],RSz[q],0.0,el,nx,ny,nz)
        ΔKm+=eng.wq[n]*eng.J[n]*ΔK[q]
    end
    ΔKm/=wt
    @inbounds for q in 1:nn; n=o+q; Rτ[q]=eng.τeq[n]+ΔK[q]+(mτ-ΔKm)*eng.τeq[n]/τeqm; end
    θ=1.0
    @inbounds for q in 1:nn
        n=o+q; Dfl=0.5*eng.sqrtγ[n]*eng.ρ_atm
        if st.D[n] < Dfl
            RD[q] > Dfl || return _limit_mean_elem!(st,eng,e)
            θ=min(θ, clamp((RD[q]-Dfl)/(RD[q]-st.D[n]+1e-300),0.0,1.0))
        end
    end
    if θ<1; @inbounds for q in 1:nn; n=o+q; st.D[n]=RD[q]+θ*(st.D[n]-RD[q]); end; end
    θq=1.0
    @inbounds for q in 1:nn
        n=o+q; el=eng.elam[n]; nx=eng.nx[n]; ny=eng.ny[n]; nz=eng.nz[n]
        qR=_q3(RD[q],RSx[q],RSy[q],RSz[q],Rτ[q],el,nx,ny,nz); qR > 0 || return _limit_mean_elem!(st,eng,e)
        qn=_q3(st.D[n],st.Sx[n],st.Sy[n],st.Sz[n],st.τ[n],el,nx,ny,nz); qε=1e-12*qR
        qn < qε && (θq=min(θq, clamp((qR-qε)/(qR-qn+1e-300),0.0,1.0)))
    end
    if θq<1
        @inbounds for q in 1:nn
            n=o+q
            st.D[n]=RD[q]+θq*(st.D[n]-RD[q]); st.Sx[n]=RSx[q]+θq*(st.Sx[n]-RSx[q]); st.Sy[n]=RSy[q]+θq*(st.Sy[n]-RSy[q])
            st.Sz[n]=RSz[q]+θq*(st.Sz[n]-RSz[q]); st.τ[n]=Rτ[q]+θq*(st.τ[n]-Rτ[q])
        end
    end
    return (θ<1 || θq<1) ? 1 : 0
end
function _limit_mean_elem!(st::DGBall3DState, eng::DGBall3DEngine, e::BallElem)
    D̄=_elem_mean(eng,st.D,e); S̄x=_elem_mean(eng,st.Sx,e); S̄y=_elem_mean(eng,st.Sy,e); S̄z=_elem_mean(eng,st.Sz,e); τ̄=_elem_mean(eng,st.τ,e)
    nn=nnodes(e); o=e.off-1; θ=1.0
    @inbounds for q in 1:nn
        n=o+q; Dfl=eng.sqrtγ[n]*eng.ρ_atm
        st.D[n] < Dfl && (θ=min(θ, D̄>Dfl ? clamp((D̄-Dfl)/(D̄-st.D[n]+1e-300),0.0,1.0) : 0.0))
    end
    n0=o+1; q̄=_q3(D̄,S̄x,S̄y,S̄z,τ̄,eng.elam[n0],eng.nx[n0],eng.ny[n0],eng.nz[n0]); qε=1e-12*max(q̄,0.0); θq=1.0
    @inbounds for q in 1:nn
        n=o+q; Dn=D̄+θ*(st.D[n]-D̄)
        qn=_q3(Dn,st.Sx[n],st.Sy[n],st.Sz[n],st.τ[n],eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n])
        qn < qε && (θq=min(θq, q̄>qε ? clamp((q̄-qε)/(q̄-qn+1e-300),0.0,1.0) : 0.0))
    end
    θ=min(θ,θq)
    @inbounds for q in 1:nn
        n=o+q
        st.D[n]=D̄+θ*(st.D[n]-D̄); st.Sx[n]=S̄x+θ*(st.Sx[n]-S̄x); st.Sy[n]=S̄y+θ*(st.Sy[n]-S̄y); st.Sz[n]=S̄z+θ*(st.Sz[n]-S̄z); st.τ[n]=τ̄+θ*(st.τ[n]-τ̄)
    end
    return 2
end
function _limit!(st::DGBall3DState, eng::DGBall3DEngine)
    eng.limiter == :none && return
    Threads.@threads for e in eng.elems
        e.region in eng.limit_regions || continue
        if eng.limiter == :wb; _limit_wb_elem!(st,eng,e)
        elseif eng.limiter == :mean; _limit_mean_elem!(st,eng,e)
        else error("unknown limiter $(eng.limiter)") end
    end
end

# exponential modal filter (tensor product) on the deviation from equilibrium, applied after
# every full step to the elements selected by `filter` (:all, or :center = cube + centre shells,
# the paper's "central portion"), with the paper's strengths: exp[−α(i/p)^s], s = filt_s_center in
# the cube/centre shells and filt_s_shell in the cubed-sphere shells; `filt_vars` = :momentum
# (their choice: S̃_i only) or :all.
function _filter!(st::DGBall3DState, eng::DGBall3DEngine)
    eng.filter == :none && return
    Threads.@threads for e in eng.elems
        central = e.region == :cube || e.region == :center
        (eng.filter == :all || central) || continue
        N1,N2,N3=e.N; bs=(eng.bases[e.p[1]],eng.bases[e.p[2]],eng.bases[e.p[3]])
        sexp = central ? eng.filt_s_center : eng.filt_s_shell
        σ=[ [exp(-eng.filt_α*((i-1)/max(e.p[d],1))^sexp) for i in 1:e.N[d]] for d in 1:3 ]
        buf=zeros(N1,N2,N3); tmp=zeros(N1,N2,N3)
        vars = eng.filt_vars == :all ? ((st.D,eng.Deq),(st.Sx,eng.Sxeq),(st.Sy,eng.Syeq),(st.Sz,eng.Szeq),(st.τ,eng.τeq)) :
                                        ((st.Sx,eng.Sxeq),(st.Sy,eng.Syeq),(st.Sz,eng.Szeq))
        for (U,Ueq) in vars
            for k in 1:N3, j in 1:N2, i in 1:N1; n=nidx(e,i,j,k); buf[i,j,k]=U[n]-Ueq[n]; end
            for k in 1:N3, j in 1:N2
                v=bs[1].invV*buf[:,j,k]; v.*=σ[1]; tmp[:,j,k]=bs[1].V*v
            end
            for k in 1:N3, i in 1:N1
                v=bs[2].invV*tmp[i,:,k]; v.*=σ[2]; buf[i,:,k]=bs[2].V*v
            end
            for j in 1:N2, i in 1:N1
                v=bs[3].invV*buf[i,j,:]; v.*=σ[3]; tmp[i,j,:]=bs[3].V*v
            end
            for k in 1:N3, j in 1:N2, i in 1:N1; n=nidx(e,i,j,k); U[n]=Ueq[n]+tmp[i,j,k]; end
        end
    end
end

# ----------------------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------------------
# Lagrange interpolant of a nodal field at reference coordinates (ξ,η,ζ) of element e
function _interp(eng::DGBall3DEngine, e::BallElem, U, ξ, η, ζ)
    bs=(eng.bases[e.p[1]],eng.bases[e.p[2]],eng.bases[e.p[3]]); ξs=(ξ,η,ζ)
    L=[[begin ℓ=1.0; for c in 1:e.N[d]; c==a && continue; ℓ*=(ξs[d]-bs[d].ξ[c])/(bs[d].ξ[a]-bs[d].ξ[c]); end; ℓ end for a in 1:e.N[d]] for d in 1:3]
    s=0.0
    @inbounds for k in 1:e.N[3], j in 1:e.N[2], i in 1:e.N[1]; s+=L[1][i]*L[2][j]*L[3][k]*U[nidx(e,i,j,k)]; end
    s
end
"""central rest-mass density: the nodal value(s) at the origin when it is a node (even `nt`),
otherwise the interpolant of the cube element containing the origin (odd `nt`)."""
function dgball3d_central_density(st::DGBall3DState, eng::DGBall3DEngine)
    isempty(eng.origin) || return sum(st.ρ[eng.origin])/length(eng.origin)
    for e in eng.elems
        e.region == :cube || continue
        m=e.map::CubeMap
        (m.ξlo < 0 < m.ξhi && m.ηlo < 0 < m.ηhi && m.ζlo < 0 < m.ζhi) || continue
        ξ=-1+2*(0-m.ξlo)/(m.ξhi-m.ξlo); η=-1+2*(0-m.ηlo)/(m.ηhi-m.ηlo); ζ=-1+2*(0-m.ζlo)/(m.ζhi-m.ζlo)
        return _interp(eng,e,st.ρ,ξ,η,ζ)
    end
    error("no cube element contains the origin")
end
dgball3d_errD(st::DGBall3DState, eng::DGBall3DEngine) = sqrt(sum(abs2, st.D .- eng.Deq)/sum(abs2, eng.Deq))
dgball3d_baryon_mass(st::DGBall3DState, eng::DGBall3DEngine) = sum(eng.wq .* eng.J .* st.D)
dgball3d_volume(eng::DGBall3DEngine) = sum(eng.wq .* eng.J)
"""real spherical-harmonic moment ∫ D̃ r^ℓ Y_ℓm d³x for (ℓ,m) ∈ {(0,0),(2,0),(2,2),(2,1)} (unnormalized Y)."""
function dgball3d_moment(st::DGBall3DState, eng::DGBall3DEngine, l::Int, m::Int)
    s=0.0
    @inbounds for n in 1:eng.Ntot
        r=eng.r[n]; r>0 || continue
        Y = (l==0) ? 1.0 : (l==2 && m==0) ? (3*eng.nz[n]^2-1)/2 : (l==2 && m==2) ? (eng.nx[n]^2-eng.ny[n]^2) : (l==2 && m==1) ? eng.nx[n]*eng.nz[n] : error("moment ($l,$m) not implemented")
        s+=eng.wq[n]*eng.J[n]*st.D[n]*r^l*Y
    end
    s
end
"""max static momentum residual / gravitational force per region, on nodes with ρ > 1e-6 ρ_c."""
function dgball3d_static_residual(st::DGBall3DState, eng::DGBall3DEngine)
    r=(zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot))
    _raw_rhs!(r..., st, eng)
    ρc=dgball3d_central_density(st,eng); out=Dict{Symbol,Float64}()
    for e in eng.elems
        w=get(out,e.region,0.0)
        for q in 1:nnodes(e)
            n=e.off+q-1; ρ=st.ρ[n]; ρ > 1e-6*ρc || continue
            grav=eng.α[n]*eng.sqrtγ[n]*(ρ*(1+st.ε[n])+st.p[n])*abs(eng.Φp[n])
            grav > 0 || continue
            w=max(w, sqrt(r[2][n]^2+r[3][n]^2+r[4][n]^2)/grav)
        end
        out[e.region]=w
    end
    out
end
"""max |rhs| over all nodes and variables (for the free-stream test)."""
function dgball3d_rhs_norm(st::DGBall3DState, eng::DGBall3DEngine)
    r=(zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot),zeros(eng.Ntot))
    _raw_rhs!(r..., st, eng)
    maximum(maximum(abs, v) for v in r)
end

# ----------------------------------------------------------------------------------
# evolution
# ----------------------------------------------------------------------------------
function _maxspeed(st::DGBall3DState, eng::DGBall3DEngine)
    a=0.0
    @inbounds for n in 1:eng.Ntot
        for (mx,my,mz) in ((1.0,0.0,0.0),(0.0,1.0,0.0),(0.0,0.0,1.0))
            λm,λp=_speeds(eng.Γ,st.ρ[n],st.ε[n],st.vx[n],st.vy[n],st.vz[n],eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n],eng.α[n],mx,my,mz)
            a=max(a,abs(λm),abs(λp))
        end
    end
    a
end

"""
    evolve_dgball3d!(st, eng; tmax, sample_dt=2.0, dt=0.0, verbose=false) -> rec

SSP-RK3 with the limiter after every substep and the filter after every step; Δt fixed from
cfl·Δx_min/a_max at t=0 unless given. Records t, ρc, errD, Mb, q00, q20, q22, vatm (max |v| for
r > R), nfix.
"""
function evolve_dgball3d!(st::DGBall3DState, eng::DGBall3DEngine; tmax::Float64, sample_dt::Float64=2.0,
                          dt::Float64=0.0, verbose::Bool=false)
    n=eng.Ntot
    r=(zeros(n),zeros(n),zeros(n),zeros(n),zeros(n)); U0=(zeros(n),zeros(n),zeros(n),zeros(n),zeros(n))
    U=(st.D,st.Sx,st.Sy,st.Sz,st.τ)
    _update_prims!(st,eng)
    dt = dt>0 ? dt : eng.cfl*eng.dxmin/max(_maxspeed(st,eng),1e-3)
    ts=Float64[]; ρc=Float64[]; errD=Float64[]; Mb=Float64[]; q00=Float64[]; q20=Float64[]; q22=Float64[]; vatm=Float64[]
    outside=eng.r .> eng.R
    record!(t) = begin
        push!(ts,t); push!(ρc,dgball3d_central_density(st,eng)); push!(errD,dgball3d_errD(st,eng)); push!(Mb,dgball3d_baryon_mass(st,eng))
        push!(q00,dgball3d_moment(st,eng,0,0)); push!(q20,dgball3d_moment(st,eng,2,0)); push!(q22,dgball3d_moment(st,eng,2,2))
        vm=0.0; @inbounds for i in 1:n; outside[i] && (vm=max(vm, sqrt(st.vx[i]^2+st.vy[i]^2+st.vz[i]^2))); end; push!(vatm,vm)
    end
    record!(0.0); t=0.0; last=0.0; ns=0; t0=time()
    stage!(c0,c1,h) = begin
        _rhs!(r..., st, eng)
        Threads.@threads for i in 1:n
            @inbounds for c in 1:5; U[c][i]=c0*U0[c][i]+c1*(U[c][i]+h*r[c][i]); end
        end
        _limit!(st,eng); _update_prims!(st,eng)
    end
    while t < tmax-1e-12
        h=min(dt,tmax-t)
        for c in 1:5; copyto!(U0[c],U[c]); end
        stage!(0.0,1.0,h); stage!(0.75,0.25,h); stage!(1/3,2/3,h)
        _filter!(st,eng); eng.filter != :none && _update_prims!(st,eng)
        t+=h; ns+=1
        if t-last ≥ sample_dt-1e-12
            record!(t); last=t
            verbose && @info "t=$(round(t,digits=1)) errD=$(errD[end]) ρc/ρc0-1=$(ρc[end]/ρc[1]-1) [$(round(time()-t0)) s]"
            (isfinite(ρc[end]) && ρc[end] < 100*ρc[1]) || break
        end
    end
    return (t=ts, ρc=ρc, errD=errD, Mb=Mb, q00=q00, q20=q20, q22=q22, vatm=vatm, nfix=st.nfix, dt=dt, nsteps=ns, seconds=time()-t0)
end

"""radial seed v^a = A sin(πr/R) n^a (vanishes at the surface)."""
function seed_dgball3d_radial!(st::DGBall3DState, eng::DGBall3DEngine; A::Float64=1e-3)
    @inbounds for n in 1:eng.Ntot
        rr=eng.r[n]; (rr < eng.R && st.ρ[n] > 10eng.ρ_atm) || continue
        v=A*sin(π*rr/eng.R); vx=v*eng.nx[n]; vy=v*eng.ny[n]; vz=v*eng.nz[n]
        st.vx[n]=vx; st.vy[n]=vy; st.vz[n]=vz
        D̂,Ŝx,Ŝy,Ŝz,τ̂=_p2c(eng.Γ,st.ρ[n],st.ε[n],vx,vy,vz,eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n]); sg=eng.sqrtγ[n]
        st.D[n]=sg*D̂; st.Sx[n]=sg*Ŝx; st.Sy[n]=sg*Ŝy; st.Sz[n]=sg*Ŝz; st.τ[n]=sg*τ̂
    end
end
"""ℓ=2, m=0 seed v^a = A (r/R) Y₂₀(n) sin(πr/R)^0 n^a  (Y₂₀ = (3n_z²−1)/2), as in DGCart3D but vanishing at R when `taper=true`."""
function seed_dgball3d_l2!(st::DGBall3DState, eng::DGBall3DEngine; A::Float64=1e-3, taper::Bool=true)
    @inbounds for n in 1:eng.Ntot
        rr=eng.r[n]; (rr < eng.R && st.ρ[n] > 10eng.ρ_atm) || continue
        Y=(3*eng.nz[n]^2-1)/2; shape = taper ? sin(π*rr/eng.R) : rr/eng.R
        v=A*shape*Y; vx=v*eng.nx[n]; vy=v*eng.ny[n]; vz=v*eng.nz[n]
        st.vx[n]=vx; st.vy[n]=vy; st.vz[n]=vz
        D̂,Ŝx,Ŝy,Ŝz,τ̂=_p2c(eng.Γ,st.ρ[n],st.ε[n],vx,vy,vz,eng.elam[n],eng.nx[n],eng.ny[n],eng.nz[n]); sg=eng.sqrtγ[n]
        st.D[n]=sg*D̂; st.Sx[n]=sg*Ŝx; st.Sy[n]=sg*Ŝy; st.Sz[n]=sg*Ŝz; st.τ[n]=sg*τ̂
    end
end

end # module DGBall3D
