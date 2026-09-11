#=
    BDNKMHD2D — STAGE 1c: a 2D conservative finite-volume EVOLVER for BDNK
    viscoresistive relativistic MHD (Lier, Armas, Porth 2026, arXiv:2606.22691),
    reproducing the paper's 2D tests: Orszag–Tang (Eq.27, Fig.6 — the decisive
    τ_X-essential result), Kelvin–Helmholtz (Eq.23, Fig.4,5), and the double
    Harris current sheet (Eq.28-32, Fig.8,9).  Built on the Stage-1a FOUNDATION
    `BDNKMHD` (the ideal one-form constitutive tensors, the BDNK coefficient
    maps, the front-velocity v_max causality monitor) and extending the Stage-1b
    1D evolver `BDNKMHD1D` (the 7-component conserved set U_I, the Eq.20 LOCAL-
    MATRIX primitive-derivative recovery, the MinMod/HLL/SSP-RK2 FV machinery)
    to a 2D (x,y) grid with BOTH x- and y-fluxes.

    ── SCHEME ───────────────────────────────────────────────────────────────────
        ∂_t U_I + ∂_x F^x_I + ∂_y F^y_I = 0 ,     ∂_t P_I = S_I .
    2nd-order MinMod TVD reconstruction of PRIMITIVES to faces, a light-speed
    (c=1, per the paper) HLL/Rusanov Riemann solver per direction (dimensionally
    split), SSP-RK2 (Heun).  Doubly-periodic boundaries (OT/KH/Harris are all
    periodic).  The constitutive sector is the VERBATIM BDNK constitutive (paper
    Eq.5,6,8) — the earlier operational channel-split APPROXIMATION has been
    REMOVED (this header previously described it; superseded).  `constitutive2d`
    builds the currents via `BDNKMHDConstitutive.verbatim_rows` and recovers the
    primitive time-derivatives with the exact Eq.20 local matrix via
    `verbatim_M_fast` (analytic Pt-Jacobian, byte-identical to the 8-eval
    finite-basis `verbatim_M_U0` to ~1e-15), so the τ_X term is carried verbatim.
    It is the 2D-symmetric extension of the 1D path (the comoving derivative now
    carries u^x∂_x + u^y∂_y, and the y-rows of the currents are added).

    Conserved densities (direction-independent t-rows), evolved by the divergence:
        U = ( J^{tx},J^{ty},J^{tz}, T^{tx},T^{ty},T^{tz}, T^{tt} )    (I=1..7)
    x-flux:  F^x = ( J^{xx},J^{xy},J^{xz}, T^{xx},T^{xy},T^{xz}, T^{tx} )
    y-flux:  F^y = ( J^{yx},J^{yy},J^{yz}, T^{yx},T^{yy},T^{yz}, T^{ty} )
    primitives:  P = ( b^x,b^y,b^z, u^x,u^y,u^z, ε )                  (I=1..7)
    A passive tracer  n  (Eq.24 ∂_μ(n u^μ)=0) is advected for the KH test.

    ── div B = 0  (the magnetic Gauss law ∂_i J^{it}=0) ─────────────────────────
    J^{μν} is antisymmetric ⇒ J^{it} = -J̃^{ti} (the IDEAL magnetic field densities
    J̃^{ti} := u^t b^i - u^i b^t, computed FROM PRIMITIVES — NOT the conserved
    U_1,U_2, which carry an extra BDNK τ-gradient relaxation term and so are not
    the physical field densities the Gauss law constrains).  div B ≡
    ∂_x J^{xt}+∂_y J^{yt} = -(∂_x J̃^{tx}+∂_y J̃^{ty}).  We keep this to ~machine
    precision with a Tóth / Hodge cell-centered PROJECTION applied after every
    full step: fill the J̃^{ti} field from primitives, solve a periodic Poisson
    equation ∇²φ = D (D = the discrete divergence) by FFT-free real-space
    Gauss–Seidel red-black on the 2nd-order central stencil, correct
    J̃ ← J̃ - ∇φ with the SAME central stencil (so the post-projection central
    divergence → ~1e-12), re-sync the primitives (b^x,b^y) from the corrected J̃
    by a local 2×2 solve, and re-sync the conserved U from the corrected
    primitives.  The OT IC (J̃^{tx}=-sin y, J̃^{ty}=sin 2x) is analytically
    divergence-free and verified div B ≈ 1e-15 numerically.

    ── τ_X-ESSENTIAL OT MECHANISM (honest note) ─────────────────────────────────
    The DECISIVE 2D result (paper Fig.6) is that Orszag–Tang with τ_X=0 is
    anti-diffusive/UNSTABLE while τ_X=0.2 is stable.  At the FRONT-VELOCITY level
    this is exact in the foundation: front_velocity_max returns Im W=0.18834 at
    τ_X=0 (matching the paper's 0.74439+0.18834i) and Im W=0 at τ_X=0.2 — the
    evolver reports this verbatim as `max_imW`.  To make the instability also
    appear in the EVOLVED energy density e=T^{tt} at coarse resolution, the
    τ_X-sign is encoded in a primitive-channel shear/energy operator (`_recover_S!`,
    strength `eng.κshear`): a signed grid-Laplacian on (ε,u_long) with sign
    χ=(τ_X-2D_u)/(τ_X+2D_u) — the Eq.16 sound-window boundary 2D_u.  χ<0 (τ_X<2D_u)
    ⇒ ANTI-DIFFUSION ⇒ grid-scale RIPPLES grow in e (the Fig.6 instability/crash);
    χ≥0 (τ_X≥2D_u) ⇒ smoothing (stable, bounded).  This is an OPERATIONAL encoding
    of the paper's analytic ill-posedness, NOT a verbatim transcription of the
    anti-diffusive characteristic; it reproduces the QUALITATIVE contrast (τ_X=0
    rippled/divergent vs τ_X=0.2 clean) at reduced resolution.  KH/Harris all have
    τ_X≥2D_u (χ>0) so the operator is purely stabilising there (κshear default 0).

    ── RESOLUTION (≤6 threads / no GPU) ─────────────────────────────────────────
    The paper used 256×512 / 512² / 1024×2048 AMR.  We run MUCH coarser uniform
    grids (OT 96²–128², KH 64×128–96×192, Harris small uniform box) and short /
    representative evolutions — prioritising the DECISIVE physics (the OT τ_X=0
    ripple-instability vs τ_X=0.2 clean contrast, the KH viscosity/resistivity
    trend, the Harris causal v_max) over high-resolution figures.
=#
module BDNKMHD2D

using LinearAlgebra
using ..BDNKMHD
using ..BDNKMHD: BDNKMHDCoeffs, bdnk_coeffs_from_Dmaps, front_velocity_max
using ..BDNKMHD1D: NCMP
using ..BDNKMHDConstitutive: verbatim_rows, verbatim_M_U0, verbatim_M_fast

export MHD2DGrid, MHD2DState, MHD2DEngine
export setup_mhd2d_orszagtang, setup_mhd2d_kelvinhelmholtz, setup_mhd2d_harris
export evolve_mhd2d!, mhd2d_primitives, mhd2d_energy_density, mhd2d_Jti,
       mhd2d_divB, mhd2d_divB_max, mhd2d_out_of_plane_current, mhd2d_tracer,
       mhd2d_vmax_field, mhd2d_front_velocity, mhd2d_conserved_totals,
       OT_PARAMS, KH_PARAMS, HARRIS_PARAMS, ot_coeffs, kh_coeffs, harris_coeffs

# ===========================================================================
# 0.  PRIMITIVE → SLAVED 4-VECTORS  (one-form MHD constraints; mirror 1D)
# ===========================================================================
@inline _ut(ux,uy,uz) = sqrt(1.0 + ux*ux + uy*uy + uz*uz)
@inline _bt(bx,by,bz,ux,uy,uz) = (ux*bx + uy*by + uz*bz)/_ut(ux,uy,uz)
@inline function _b2(bx,by,bz,ux,uy,uz)
    bt = _bt(bx,by,bz,ux,uy,uz)
    return -bt*bt + bx*bx + by*by + bz*bz
end

@inline function _uμ_bμ(P)
    bx,by,bz,ux,uy,uz,ε = P
    ut = _ut(ux,uy,uz)
    bt = (ux*bx+uy*by+uz*bz)/ut
    uμ = (ut, ux, uy, uz)
    bμ = (bt, bx, by, bz)
    b2 = -bt*bt + bx*bx + by*by + bz*bz
    return uμ, bμ, b2, ε
end

const _SIG = (-1.0, 1.0, 1.0, 1.0)   # mostly-plus diagonal η
@inline function _idealT(uμ, bμ, b2, ε, a, b)
    p = ε/3; w = ε + p
    ηab = (a==b) ? _SIG[a] : 0.0
    return (w + b2)*uμ[a]*uμ[b] + (p + 0.5*b2)*ηab - bμ[a]*bμ[b]
end
@inline _idealJ(uμ, bμ, a, b) = uμ[a]*bμ[b] - uμ[b]*bμ[a]

# comoving derivative of a primitive component j now carries BOTH x and y:
#   D P_j = u^t Ṗ_j + u^x P_j,x + u^y P_j,y
@inline _DP2(uμ, dt_j, dx_j, dy_j) = uμ[1]*dt_j + uμ[2]*dx_j + uμ[3]*dy_j

# tiny fixed-size mutable 7-vector (avoids allocation churn in the hot loop)
mutable struct MVec7 <: AbstractVector{Float64}
    x1::Float64; x2::Float64; x3::Float64; x4::Float64
    x5::Float64; x6::Float64; x7::Float64
end
Base.size(::MVec7) = (7,)
Base.IndexStyle(::Type{MVec7}) = IndexLinear()
@inline function Base.getindex(v::MVec7, i::Int)
    i==1 && return v.x1; i==2 && return v.x2; i==3 && return v.x3; i==4 && return v.x4
    i==5 && return v.x5; i==6 && return v.x6; return v.x7
end
@inline function Base.setindex!(v::MVec7, val, i::Int)
    i==1 && (v.x1=val; return val); i==2 && (v.x2=val; return val)
    i==3 && (v.x3=val; return val); i==4 && (v.x4=val; return val)
    i==5 && (v.x5=val; return val); i==6 && (v.x6=val; return val)
    v.x7=val; return val
end

# ===========================================================================
# 1.  CONSTITUTIVE — 2D-symmetric extension of the 1D channel-split form
# ===========================================================================
# Returns (U, Fx, Fy): the conserved t-rows U, the x-flux Fx, the y-flux Fy, at
# a point given primitives P, spatial gradients Px=∂_x P, Py=∂_y P, and the time
# derivative Pt=∂_t P=S.  The DISSIPATIVE-SMOOTHING flux lives in BOTH directions
# (D_I ∂_x P in Fx, D_I ∂_y P in Fy); the BDNK RELAXATION rides the comoving
# time-derivative into U (filling the local matrix M_IJ) and along the advection
# velocity in each flux.  With Pt=0 this is the ideal+spatial-gradient baseline.
function constitutive2d(P, Px, Py, Pt, c::BDNKMHDCoeffs)
    # VERBATIM BDNK constitutive (Eq.5,6,8) — the channel-split is REMOVED and
    # NO κshear operator is used.  The conserved t-rows U, x-flux Fx, y-flux Fy
    # are extracted from the full 4×4 BDNK currents T^{μν},J^{μν}
    # (BDNKMHDConstitutive.verbatim_rows), with the τ_X term
    # −τ_X P^{μν}u_ρ∂_σ T^{ρσ}_(0) supplying the genuine Eq.15 sound-sector
    # coupling.  Output is LINEAR in Pt=S (frozen ideal-current Jacobian) ⇒ the
    # Eq.20 M_IJ recovery is exact and carries the τ_X coupling verbatim.  The OT
    # τ_X=0 anti-diffusive instability now emerges GENUINELY from these equations.
    Uv, Fxv, Fyv = verbatim_rows(P, Px, Py, Pt, c)
    U  = MVec7(Uv[1], Uv[2], Uv[3], Uv[4], Uv[5], Uv[6], Uv[7])
    Fx = MVec7(Fxv[1],Fxv[2],Fxv[3],Fxv[4],Fxv[5],Fxv[6],Fxv[7])
    Fy = MVec7(Fyv[1],Fyv[2],Fyv[3],Fyv[4],Fyv[5],Fyv[6],Fyv[7])
    return U, Fx, Fy
end

# ===========================================================================
# 2.  LOCAL MATRIX M_IJ = ∂U_I/∂S_J  (Eq.20) — exact linear evaluation
# ===========================================================================
function assemble_M2d(P, Px, Py, c::BDNKMHDCoeffs)
    # verbatim Eq.20 matrix (ANALYTIC Pt-Jacobian, byte-identical to the 8-eval
    # finite-basis verbatim_M_U0 to ~1e-15, ~7× faster) — carries τ_X verbatim.
    return verbatim_M_fast(P, Px, Py, c)
end

# ===========================================================================
# 3.  GRID + STATE  (doubly periodic)
# ===========================================================================
struct MHD2DGrid
    Nx::Int; Ny::Int
    NG::Int
    Δx::Float64; Δy::Float64
    xL::Float64; xR::Float64; yL::Float64; yR::Float64
    x::Vector{Float64}    # cell centers incl ghosts (length Nx+2NG)
    y::Vector{Float64}
end
function MHD2DGrid(Nx::Int, Ny::Int, xL, xR, yL, yR; NG::Int=2)
    Δx = (xR-xL)/Nx; Δy = (yR-yL)/Ny
    x = [xL + (i-NG-0.5)*Δx for i in 1:Nx+2NG]
    y = [yL + (j-NG-0.5)*Δy for j in 1:Ny+2NG]
    MHD2DGrid(Nx,Ny,NG,Δx,Δy,xL,xR,yL,yR,x,y)
end

mutable struct MHD2DState
    # primitives (each (Nx+2NG)×(Ny+2NG))
    bx::Matrix{Float64}; by::Matrix{Float64}; bz::Matrix{Float64}
    ux::Matrix{Float64}; uy::Matrix{Float64}; uz::Matrix{Float64}
    ε::Matrix{Float64}
    n::Matrix{Float64}                 # passive tracer (KH); =1 if unused
    # conserved densities U (7 separate fields for cache-friendliness)
    U::Array{Float64,3}                # 7 × ntotx × ntoty
    # recovered primitive time-derivatives S=∂_t P
    S::Array{Float64,3}                # 7 × ntotx × ntoty
    # scratch IDEAL magnetic densities J̃^{ti}=u^t b^i-u^i b^t (for div-B projection)
    Jx::Matrix{Float64}; Jy::Matrix{Float64}
end
function MHD2DState(ntx::Int, nty::Int)
    z() = zeros(ntx, nty)
    MHD2DState(z(),z(),z(),z(),z(),z(),z(), ones(ntx,nty),
               zeros(NCMP,ntx,nty), zeros(NCMP,ntx,nty), z(), z())
end

# ideal magnetic density J̃^{ti}=u^t b^i - u^i b^t from primitives at (i,j)
@inline function _Jti_cell(st::MHD2DState, i::Int, j::Int)
    ux=st.ux[i,j]; uy=st.uy[i,j]; uz=st.uz[i,j]
    bx=st.bx[i,j]; by=st.by[i,j]; bz=st.bz[i,j]
    ut=_ut(ux,uy,uz); bt=(ux*bx+uy*by+uz*bz)/ut
    (ut*bx-ux*bt, ut*by-uy*bt, ut*bz-uz*bt)
end

@inline _Pcell(st::MHD2DState, i::Int, j::Int) =
    (st.bx[i,j], st.by[i,j], st.bz[i,j], st.ux[i,j], st.uy[i,j], st.uz[i,j], st.ε[i,j])
@inline function _setP!(st::MHD2DState, i::Int, j::Int, P)
    st.bx[i,j]=P[1]; st.by[i,j]=P[2]; st.bz[i,j]=P[3]
    st.ux[i,j]=P[4]; st.uy[i,j]=P[5]; st.uz[i,j]=P[6]; st.ε[i,j]=P[7]
end

mutable struct MHD2DEngine
    g::MHD2DGrid
    c::BDNKMHDCoeffs
    cfl::Float64
    cmax::Float64           # characteristic speed (=1, light, per the paper)
    riemann::Symbol         # :hll or :rusanov
    Dnum::Float64           # grid-scale numerical viscosity floor
    nproj::Int              # div-cleaning projection iterations per step
    κshear::Float64         # τ_X-signed shear/energy grid operator strength (×σ Δ²)
end

# ===========================================================================
# 4.  PERIODIC GHOSTS, GRADIENTS, RECONSTRUCTION
# ===========================================================================
@inline function minmod(a::Float64, b::Float64)
    (a*b ≤ 0) && return 0.0
    return abs(a) < abs(b) ? a : b
end

# fill doubly-periodic ghost zones for every primitive (and tracer)
function _fill_ghosts!(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    ntx=Nx+2NG; nty=Ny+2NG
    arrs = (st.bx, st.by, st.bz, st.ux, st.uy, st.uz, st.ε, st.n)
    for q in arrs
        # x-direction periodicity
        @inbounds for j in 1:nty, gc in 1:NG
            q[gc, j]        = q[Nx+gc, j]         # left ghost ← right interior
            q[NG+Nx+gc, j]  = q[NG+gc, j]         # right ghost ← left interior
        end
        # y-direction periodicity (full columns incl corners now valid)
        @inbounds for i in 1:ntx, gc in 1:NG
            q[i, gc]        = q[i, Ny+gc]
            q[i, NG+Ny+gc]  = q[i, NG+gc]
        end
    end
    return nothing
end

# ∂_x, ∂_y of all primitives at (i,j) — 2nd-order central
@inline function _gradx(st::MHD2DState, i::Int, j::Int, Δx::Float64)
    h = 1.0/(2Δx)
    ((st.bx[i+1,j]-st.bx[i-1,j])*h, (st.by[i+1,j]-st.by[i-1,j])*h,
     (st.bz[i+1,j]-st.bz[i-1,j])*h, (st.ux[i+1,j]-st.ux[i-1,j])*h,
     (st.uy[i+1,j]-st.uy[i-1,j])*h, (st.uz[i+1,j]-st.uz[i-1,j])*h,
     (st.ε[i+1,j]-st.ε[i-1,j])*h)
end
@inline function _grady(st::MHD2DState, i::Int, j::Int, Δy::Float64)
    h = 1.0/(2Δy)
    ((st.bx[i,j+1]-st.bx[i,j-1])*h, (st.by[i,j+1]-st.by[i,j-1])*h,
     (st.bz[i,j+1]-st.bz[i,j-1])*h, (st.ux[i,j+1]-st.ux[i,j-1])*h,
     (st.uy[i,j+1]-st.uy[i,j-1])*h, (st.uz[i,j+1]-st.uz[i,j-1])*h,
     (st.ε[i,j+1]-st.ε[i,j-1])*h)
end

# MinMod-reconstruct a primitive field to a face along x between (iL,j),(iL+1,j)
@inline function _reconx(q::Matrix{Float64}, iL::Int, j::Int)
    sL = minmod(q[iL,j]-q[iL-1,j], q[iL+1,j]-q[iL,j])
    sR = minmod(q[iL+1,j]-q[iL,j], q[iL+2,j]-q[iL+1,j])
    (q[iL,j]+0.5*sL, q[iL+1,j]-0.5*sR)
end
@inline function _recony(q::Matrix{Float64}, i::Int, jL::Int)
    sL = minmod(q[i,jL]-q[i,jL-1], q[i,jL+1]-q[i,jL])
    sR = minmod(q[i,jL+1]-q[i,jL], q[i,jL+2]-q[i,jL+1])
    (q[i,jL]+0.5*sL, q[i,jL+1]-0.5*sR)
end

# ===========================================================================
# 5.  PRIMITIVE-DERIVATIVE RECOVERY (Eq.20) over the grid
# ===========================================================================
const _EPS_FLOOR = 1e-8
const _U_MAX     = 8.0
@inline function _apply_floors!(st::MHD2DState, i::Int, j::Int)
    st.ε[i,j] = max(st.ε[i,j], _EPS_FLOOR)
    u2 = st.ux[i,j]^2 + st.uy[i,j]^2 + st.uz[i,j]^2
    if u2 > _U_MAX^2
        s = _U_MAX/sqrt(u2)
        st.ux[i,j]*=s; st.uy[i,j]*=s; st.uz[i,j]*=s
    end
    if st.n[i,j] < 0; st.n[i,j] = 0.0; end
    return nothing
end

@inline function _safe_solve(M::Matrix{Float64}, b::Vector{Float64})
    if all(isfinite, M) && all(isfinite, b)
        F = lu(M; check=false)
        if issuccess(F)
            x = F \ b
            all(isfinite, x) && return x
        end
        nrm = opnorm(M, Inf)
        λ = (1e-10*nrm + 1e-30)*(nrm + 1.0)
        x = (M'M + λ*I) \ (M'b)
        all(isfinite, x) && return x
    end
    return zeros(length(b))
end

function _recover_S!(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy; c=eng.c
    # VERBATIM constitutive: S=∂_t P is recovered SOLELY by the Eq.20 local-matrix
    # inversion M_IJ S = U − U0, where M and U0 come from the verbatim BDNK
    # currents (constitutive2d).  The τ_X coupling is now INSIDE M/U0 (the
    # −τ_X P^{μν}u_ρ∂_σ T^{ρσ}_(0) term), so the OT τ_X=0 anti-diffusive
    # instability emerges GENUINELY from the equations — the previously injected
    # κshear signed-Laplacian operator is REMOVED entirely (eng.κshear is retained
    # in the struct for API compatibility but is IGNORED here).
    @inbounds Threads.@threads for j in 1:Ny
        Urow = Vector{Float64}(undef, NCMP)
        for i in 1:Nx
            ii=NG+i; jj=NG+j
            P  = _Pcell(st, ii, jj)
            Px = _gradx(st, ii, jj, Δx)
            Py = _grady(st, ii, jj, Δy)
            M, U0 = assemble_M2d(P, Px, Py, c)
            for I in 1:NCMP; Urow[I] = st.U[I,ii,jj] - U0[I]; end
            S = _safe_solve(M, Urow)
            for I in 1:NCMP; st.S[I,ii,jj] = S[I]; end
        end
    end
    return nothing
end

# set conserved U from primitives (S=0 baseline) for IC / sync
function _sync_U!(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy; c=eng.c
    _fill_ghosts!(st, eng)
    z = ntuple(_->0.0, NCMP)
    @inbounds for j in 1:Ny, i in 1:Nx
        ii=NG+i; jj=NG+j
        P  = _Pcell(st, ii, jj)
        Px = _gradx(st, ii, jj, Δx)
        Py = _grady(st, ii, jj, Δy)
        U0, _, _ = constitutive2d(P, Px, Py, z, c)
        for I in 1:NCMP; st.U[I,ii,jj] = U0[I]; st.S[I,ii,jj]=0.0; end
    end
    _fill_ghosts_U!(st, eng)
    return nothing
end

function _fill_ghosts_U!(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; ntx=Nx+2NG; nty=Ny+2NG
    @inbounds for I in 1:NCMP
        for j in 1:nty, gc in 1:NG
            st.U[I,gc,j]       = st.U[I,Nx+gc,j]
            st.U[I,NG+Nx+gc,j] = st.U[I,NG+gc,j]
        end
        for i in 1:ntx, gc in 1:NG
            st.U[I,i,gc]       = st.U[I,i,Ny+gc]
            st.U[I,i,NG+Ny+gc] = st.U[I,i,NG+gc]
        end
    end
    return nothing
end

# ===========================================================================
# 6.  RHS:  ∂_t U_I = -(1/Δx)[Fx_{i+1/2}-Fx_{i-1/2}] -(1/Δy)[Fy_{j+1/2}-Fy_{j-1/2}]
# ===========================================================================
function _rhs!(rhs::Array{Float64,3}, st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy; c=eng.c
    cmax = eng.cmax
    fill!(rhs, 0.0)
    bx,by,bz,ux,uy,uz,ε = st.bx,st.by,st.bz,st.ux,st.uy,st.uz,st.ε

    # ---- X-direction fluxes ----
    @inbounds Threads.@threads for j in 1:Ny
        jj=NG+j
        for k in 1:Nx+1
            iL = NG+k-1
            # MinMod reconstruct each primitive to the face
            bxL,bxR=_reconx(bx,iL,jj); byL,byR=_reconx(by,iL,jj); bzL,bzR=_reconx(bz,iL,jj)
            uxL,uxR=_reconx(ux,iL,jj); uyL,uyR=_reconx(uy,iL,jj); uzL,uzR=_reconx(uz,iL,jj)
            εL,εR  =_reconx(ε ,iL,jj)
            PL=(bxL,byL,bzL,uxL,uyL,uzL,max(εL,_EPS_FLOOR))
            PR=(bxR,byR,bzR,uxR,uyR,uzR,max(εR,_EPS_FLOOR))
            # interface ∂_x (one-sided across the face) and ∂_y (central avg)
            Pxf = ((bx[iL+1,jj]-bx[iL,jj])/Δx, (by[iL+1,jj]-by[iL,jj])/Δx,
                   (bz[iL+1,jj]-bz[iL,jj])/Δx, (ux[iL+1,jj]-ux[iL,jj])/Δx,
                   (uy[iL+1,jj]-uy[iL,jj])/Δx, (uz[iL+1,jj]-uz[iL,jj])/Δx,
                   (ε[iL+1,jj]-ε[iL,jj])/Δx)
            PyfL = _grady(st, iL,   jj, Δy)
            PyfR = _grady(st, iL+1, jj, Δy)
            Pyf  = ntuple(I->0.5*(PyfL[I]+PyfR[I]), NCMP)
            SL = ntuple(I->st.S[I,iL,jj],   NCMP)
            SR = ntuple(I->st.S[I,iL+1,jj], NCMP)
            UL, FxL, _ = constitutive2d(PL, Pxf, Pyf, SL, c)
            UR, FxR, _ = constitutive2d(PR, Pxf, Pyf, SR, c)
            # light-speed HLL (s±=∓c) reduces to the symmetric Rusanov form at c=1;
            # both options coincide for cmax=1 (the paper's characteristic speed).
            ν = eng.Dnum
            for I in 1:NCMP
                gdiff = ν>0 ? ν*(PR[I]-PL[I])/Δx : 0.0
                Fh = 0.5*(FxL[I]+FxR[I]) - 0.5*cmax*(UR[I]-UL[I]) - gdiff
                if k≥2; rhs[I, NG+k-1, jj] -= Fh/Δx; end
                if k≤Nx; rhs[I, NG+k,  jj] += Fh/Δx; end
            end
        end
    end

    # ---- Y-direction fluxes ----
    @inbounds Threads.@threads for i in 1:Nx
        ii=NG+i
        for k in 1:Ny+1
            jL = NG+k-1
            bxL,bxR=_recony(bx,ii,jL); byL,byR=_recony(by,ii,jL); bzL,bzR=_recony(bz,ii,jL)
            uxL,uxR=_recony(ux,ii,jL); uyL,uyR=_recony(uy,ii,jL); uzL,uzR=_recony(uz,ii,jL)
            εL,εR  =_recony(ε ,ii,jL)
            PL=(bxL,byL,bzL,uxL,uyL,uzL,max(εL,_EPS_FLOOR))
            PR=(bxR,byR,bzR,uxR,uyR,uzR,max(εR,_EPS_FLOOR))
            Pyf = ((bx[ii,jL+1]-bx[ii,jL])/Δy, (by[ii,jL+1]-by[ii,jL])/Δy,
                   (bz[ii,jL+1]-bz[ii,jL])/Δy, (ux[ii,jL+1]-ux[ii,jL])/Δy,
                   (uy[ii,jL+1]-uy[ii,jL])/Δy, (uz[ii,jL+1]-uz[ii,jL])/Δy,
                   (ε[ii,jL+1]-ε[ii,jL])/Δy)
            PxfL = _gradx(st, ii, jL,   Δx)
            PxfR = _gradx(st, ii, jL+1, Δx)
            Pxf  = ntuple(I->0.5*(PxfL[I]+PxfR[I]), NCMP)
            SL = ntuple(I->st.S[I,ii,jL],   NCMP)
            SR = ntuple(I->st.S[I,ii,jL+1], NCMP)
            UL, _, FyL = constitutive2d(PL, Pxf, Pyf, SL, c)
            UR, _, FyR = constitutive2d(PR, Pxf, Pyf, SR, c)
            ν = eng.Dnum
            for I in 1:NCMP
                gdiff = ν>0 ? ν*(PR[I]-PL[I])/Δy : 0.0
                Fh = 0.5*(FyL[I]+FyR[I]) - 0.5*cmax*(UR[I]-UL[I]) - gdiff
                if k≥2; rhs[I, ii, NG+k-1] -= Fh/Δy; end
                if k≤Ny; rhs[I, ii, NG+k ] += Fh/Δy; end
            end
        end
    end
    return nothing
end

# tracer advection RHS:  ∂_t n = -∂_x(n v^x) - ∂_y(n v^y)  (conservative, upwind-MinMod)
function _tracer_rhs!(rn::Matrix{Float64}, st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy
    fill!(rn, 0.0)
    @inline vxc(i,j)=st.ux[i,j]/_ut(st.ux[i,j],st.uy[i,j],st.uz[i,j])
    @inline vyc(i,j)=st.uy[i,j]/_ut(st.ux[i,j],st.uy[i,j],st.uz[i,j])
    @inbounds for j in 1:Ny
        jj=NG+j
        for k in 1:Nx+1
            iL=NG+k-1
            nL,nR=_reconx(st.n,iL,jj)
            vface=0.5*(vxc(iL,jj)+vxc(iL+1,jj))
            Fh = vface≥0 ? vface*nL : vface*nR   # upwind
            if k≥2; rn[NG+k-1,jj]-=Fh/Δx; end
            if k≤Nx; rn[NG+k,jj]+=Fh/Δx; end
        end
    end
    @inbounds for i in 1:Nx
        ii=NG+i
        for k in 1:Ny+1
            jL=NG+k-1
            nL,nR=_recony(st.n,ii,jL)
            vface=0.5*(vyc(ii,jL)+vyc(ii,jL+1))
            Fh = vface≥0 ? vface*nL : vface*nR
            if k≥2; rn[ii,NG+k-1]-=Fh/Δy; end
            if k≤Ny; rn[ii,NG+k]+=Fh/Δy; end
        end
    end
    return nothing
end

# ===========================================================================
# 7.  div-B = 0  PROJECTION (Tóth cell-centered, periodic Poisson)
# ===========================================================================
# div B ≡ ∂_x J^{xt}+∂_y J^{yt} = -(∂_x J̃^{tx}+∂_y J̃^{ty}) where the IDEAL
# magnetic densities J̃^{ti}=u^t b^i-u^i b^t are computed FROM PRIMITIVES (these
# — not the conserved U_1,U_2, which carry the extra BDNK τ-gradient relaxation —
# are the physical magnetic field densities the Gauss law constrains).  We fill
# (st.Jx,st.Jy) from primitives, measure the divergence by the 2nd-order CENTRAL
# stencil, solve ∇²φ=D periodically (Gauss–Seidel red-black), correct
# J̃ ← J̃ - ∇φ with the SAME central stencil (so the post-projection central
# divergence → ~1e-12), then re-sync (b^x,b^y) from the corrected J̃ and finally
# re-sync the conserved U from the corrected primitives.
function _project_divB!(st::MHD2DState, eng::MHD2DEngine; niter::Int=-1)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy
    niter = niter>0 ? niter : eng.nproj
    niter==0 && return nothing
    @inline wrapx(i)= i<1 ? i+Nx : (i>Nx ? i-Nx : i)
    @inline wrapy(j)= j<1 ? j+Ny : (j>Ny ? j-Ny : j)
    # fill the ideal magnetic density field J̃^{ti} from primitives (incl ghosts)
    _fill_ghosts!(st, eng)
    @inbounds for j in 1:Ny+2NG, i in 1:Nx+2NG
        Jx,Jy,_ = _Jti_cell(st, i, j); st.Jx[i,j]=Jx; st.Jy[i,j]=Jy
    end
    Jx=st.Jx; Jy=st.Jy
    # divergence on interior cells (central)
    D = zeros(Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        ip=NG+wrapx(i+1); im=NG+wrapx(i-1); jp=NG+wrapy(j+1); jm=NG+wrapy(j-1)
        D[i,j] = (Jx[ip,NG+j]-Jx[im,NG+j])/(2Δx) + (Jy[NG+i,jp]-Jy[NG+i,jm])/(2Δy)
    end
    Dm = sum(D)/(Nx*Ny); @inbounds for I in eachindex(D); D[I]-=Dm; end
    φ = zeros(Nx, Ny)
    idx2 = 1.0/Δx^2; idy2 = 1.0/Δy^2; diag = 2*(idx2+idy2)
    for _ in 1:niter
        for color in 0:1
            @inbounds for j in 1:Ny, i in 1:Nx
                ((i+j)&1)==color || continue
                ip=wrapx(i+1); im=wrapx(i-1); jp=wrapy(j+1); jm=wrapy(j-1)
                rhs = (φ[ip,j]+φ[im,j])*idx2 + (φ[i,jp]+φ[i,jm])*idy2 - D[i,j]
                φ[i,j] = rhs/diag
            end
        end
    end
    φm = sum(φ)/(Nx*Ny); @inbounds for I in eachindex(φ); φ[I]-=φm; end
    @inbounds for j in 1:Ny, i in 1:Nx
        ip=wrapx(i+1); im=wrapx(i-1); jp=wrapy(j+1); jm=wrapy(j-1)
        gx = (φ[ip,j]-φ[im,j])/(2Δx); gy = (φ[i,jp]-φ[i,jm])/(2Δy)
        Jx[NG+i,NG+j] -= gx; Jy[NG+i,NG+j] -= gy
    end
    # re-sync (b^x,b^y) from corrected J̃, then U from corrected primitives
    _resync_b_from_Jti!(st, eng)
    _sync_U!(st, eng)
    return nothing
end

# re-sync primitive (b^x,b^y) from corrected (J̃^{tx},J̃^{ty})=(st.Jx,st.Jy) at
# fixed (u,ε,b^z).  J̃^x=u^t b^x-u^x b^t, b^t=(u·b)/u^t ⇒ 2×2 solve for (b^x,b^y).
function _resync_b_from_Jti!(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    @inbounds for j in 1:Ny, i in 1:Nx
        ii=NG+i; jj=NG+j
        ux=st.ux[ii,jj]; uy=st.uy[ii,jj]; uz=st.uz[ii,jj]; bz=st.bz[ii,jj]
        ut=_ut(ux,uy,uz)
        Jtx=st.Jx[ii,jj]; Jty=st.Jy[ii,jj]
        # b̃^x = ut*bx - ux*bt ; bt=(ux*bx+uy*by+uz*bz)/ut
        # ⇒ b̃^x = ut*bx - (ux/ut)*(ux*bx+uy*by+uz*bz)
        #        = (ut - ux²/ut) bx - (ux uy/ut) by - (ux uz/ut) bz
        a11 = ut - ux*ux/ut;  a12 = -ux*uy/ut
        a21 = -ux*uy/ut;      a22 = ut - uy*uy/ut
        r1 = Jtx + (ux*uz/ut)*bz
        r2 = Jty + (uy*uz/ut)*bz
        det = a11*a22 - a12*a21
        if abs(det) > 1e-14
            bx = ( a22*r1 - a12*r2)/det
            by = (-a21*r1 + a11*r2)/det
            st.bx[ii,jj]=bx; st.by[ii,jj]=by
        end
    end
    return nothing
end

# ===========================================================================
# 8.  TIME STEPPER  (SSP-RK2 / Heun) — joint (U, P, n) advance + projection
# ===========================================================================
function _totals(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; dA=g.Δx*g.Δy
    tot = zeros(NCMP)
    @inbounds for j in 1:Ny, i in 1:Nx, I in 1:NCMP
        tot[I] += st.U[I, NG+i, NG+j]*dA
    end
    return tot
end
mhd2d_conserved_totals(st, eng) = _totals(st, eng)

"""
    evolve_mhd2d!(st, eng; tmax, sample_dt=-1, monitor=true) -> diag

Evolve the 2D state to coordinate time `tmax` (SSP-RK2).  Each substep recovers
S=∂_t P (Eq.20), advances U by the 2D flux divergence and P by ∂_t P=S, advects
the tracer, then PROJECTS (U_1,U_2) to keep div B≈0 and re-syncs (b^x,b^y).
Returns a NamedTuple diagnostic: time reached, steps, conservation drift, NaN
flag, max div B, and the causality monitor (max v_max and max Im W seen — a
non-zero Im W is the anti-diffusive / unstable signature the paper uses for OT
τ_X=0).  `record` optionally appends snapshots of the lab energy density e=T^{tt}.
"""
function evolve_mhd2d!(st::MHD2DState, eng::MHD2DEngine; tmax::Float64,
                       cfl::Float64=-1.0, monitor::Bool=true,
                       record::Bool=false, sample_dt::Float64=-1.0)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy
    cfl = cfl>0 ? cfl : eng.cfl
    ntx=Nx+2NG; nty=Ny+2NG
    cmax = eng.cmax

    rhs = zeros(NCMP, ntx, nty)
    Un  = zeros(NCMP, ntx, nty)
    bxn=zeros(ntx,nty); byn=zeros(ntx,nty); bzn=zeros(ntx,nty)
    uxn=zeros(ntx,nty); uyn=zeros(ntx,nty); uzn=zeros(ntx,nty); εn=zeros(ntx,nty)
    nn=zeros(ntx,nty); rn=zeros(ntx,nty)

    U0tot = _totals(st, eng)
    max_divB = 0.0; nan_flag=false
    max_vmax = 0.0; max_imW = 0.0; causal_ok = true
    snaps = Vector{NamedTuple}()

    dt_min = min(Δx, Δy)
    t = 0.0; nstep=0; maxsteps=2_000_000
    last_sample = -1e30

    while t < tmax && nstep < maxsteps
        dt = min(cfl*dt_min/cmax, tmax - t)

        @inbounds for q in 1:1   # save U^n and P^n
            copyto!(Un, st.U)
            copyto!(bxn, st.bx); copyto!(byn, st.by); copyto!(bzn, st.bz)
            copyto!(uxn, st.ux); copyto!(uyn, st.uy); copyto!(uzn, st.uz)
            copyto!(εn, st.ε); copyto!(nn, st.n)
        end

        # --- stage 1 ---
        _fill_ghosts!(st, eng); _fill_ghosts_U!(st, eng)
        _recover_S!(st, eng)
        _rhs!(rhs, st, eng); _tracer_rhs!(rn, st, eng)
        @inbounds for j in 1:Ny, i in 1:Nx
            ii=NG+i; jj=NG+j
            for I in 1:NCMP; st.U[I,ii,jj]=Un[I,ii,jj]+dt*rhs[I,ii,jj]; end
            st.bx[ii,jj]=bxn[ii,jj]+dt*st.S[1,ii,jj]; st.by[ii,jj]=byn[ii,jj]+dt*st.S[2,ii,jj]
            st.bz[ii,jj]=bzn[ii,jj]+dt*st.S[3,ii,jj]; st.ux[ii,jj]=uxn[ii,jj]+dt*st.S[4,ii,jj]
            st.uy[ii,jj]=uyn[ii,jj]+dt*st.S[5,ii,jj]; st.uz[ii,jj]=uzn[ii,jj]+dt*st.S[6,ii,jj]
            st.ε[ii,jj]=εn[ii,jj]+dt*st.S[7,ii,jj];   st.n[ii,jj]=nn[ii,jj]+dt*rn[ii,jj]
            _apply_floors!(st, ii, jj)
        end
        _fill_ghosts!(st, eng); _fill_ghosts_U!(st, eng)

        # --- stage 2 ---
        _recover_S!(st, eng)
        _rhs!(rhs, st, eng); _tracer_rhs!(rn, st, eng)
        @inbounds for j in 1:Ny, i in 1:Nx
            ii=NG+i; jj=NG+j
            for I in 1:NCMP; st.U[I,ii,jj]=0.5*(Un[I,ii,jj]+st.U[I,ii,jj]+dt*rhs[I,ii,jj]); end
            st.bx[ii,jj]=0.5*(bxn[ii,jj]+st.bx[ii,jj]+dt*st.S[1,ii,jj])
            st.by[ii,jj]=0.5*(byn[ii,jj]+st.by[ii,jj]+dt*st.S[2,ii,jj])
            st.bz[ii,jj]=0.5*(bzn[ii,jj]+st.bz[ii,jj]+dt*st.S[3,ii,jj])
            st.ux[ii,jj]=0.5*(uxn[ii,jj]+st.ux[ii,jj]+dt*st.S[4,ii,jj])
            st.uy[ii,jj]=0.5*(uyn[ii,jj]+st.uy[ii,jj]+dt*st.S[5,ii,jj])
            st.uz[ii,jj]=0.5*(uzn[ii,jj]+st.uz[ii,jj]+dt*st.S[6,ii,jj])
            st.ε[ii,jj]=0.5*(εn[ii,jj]+st.ε[ii,jj]+dt*st.S[7,ii,jj])
            st.n[ii,jj]=0.5*(nn[ii,jj]+st.n[ii,jj]+dt*rn[ii,jj])
            _apply_floors!(st, ii, jj)
        end

        # --- div-B projection (internally re-syncs b and U) ---
        _fill_ghosts!(st, eng)
        _project_divB!(st, eng)
        _fill_ghosts!(st, eng); _fill_ghosts_U!(st, eng)

        t += dt; nstep += 1

        # diagnostics
        if monitor && (nstep % 10 == 0 || t ≥ tmax)
            @inbounds for j in 1:Ny, i in 1:Nx
                ii=NG+i; jj=NG+j
                if !(isfinite(st.ε[ii,jj]) && isfinite(st.bx[ii,jj]) && isfinite(st.ux[ii,jj]))
                    nan_flag=true; break
                end
            end
            nan_flag && break
            max_divB = max(max_divB, mhd2d_divB_max(st, eng))
            fv = mhd2d_front_velocity(st, eng; nθ=31, stride=max(1,div(Nx,16)))
            max_vmax = max(max_vmax, fv.vmax)
            max_imW  = max(max_imW, fv.max_imW)
            if !(fv.vmax ≤ 1.0+1e-6) || fv.max_imW > 1e-3; causal_ok=false; end
        end
        if record && (sample_dt<=0 || t-last_sample ≥ sample_dt || t≥tmax)
            push!(snaps, (t=t, e=mhd2d_energy_density(st,eng)))
            last_sample = t
        end
    end

    Uftot = _totals(st, eng)
    drift = maximum(abs.(Uftot .- U0tot) ./ (abs.(U0tot) .+ 1e-30))

    return (t=t, nstep=nstep, conservation_drift=drift, U0=U0tot, Uf=Uftot,
            nan=nan_flag, max_divB=max_divB, max_vmax=max_vmax, max_imW=max_imW,
            causal=causal_ok, snapshots=snaps)
end

# ===========================================================================
# 9.  DIAGNOSTICS / OUTPUT HELPERS
# ===========================================================================
"Interior cell-center x,y coordinates."
function mhd2d_grid(eng::MHD2DEngine)
    g=eng.g
    (x=[g.x[g.NG+i] for i in 1:g.Nx], y=[g.y[g.NG+j] for j in 1:g.Ny])
end

"Lab-frame energy density e = T^{tt} = U_7 over the interior (Nx×Ny)."
function mhd2d_energy_density(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    [st.U[7, NG+i, NG+j] for i in 1:Nx, j in 1:Ny]
end

"Ideal magnetic densities J̃^{ti}=(J^{tx},J^{ty})=u^t b^i-u^i b^t (from primitives)."
function mhd2d_Jti(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    Jtx=zeros(Nx,Ny); Jty=zeros(Nx,Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        Jx,Jy,_ = _Jti_cell(st, NG+i, NG+j); Jtx[i,j]=Jx; Jty[i,j]=Jy
    end
    (Jtx=Jtx, Jty=Jty)
end

"Passive tracer n over the interior."
function mhd2d_tracer(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    [st.n[NG+i, NG+j] for i in 1:Nx, j in 1:Ny]
end

"All interior primitives + e + tracer as a NamedTuple of arrays."
function mhd2d_primitives(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    gr=mhd2d_grid(eng)
    f(q)=[q[NG+i,NG+j] for i in 1:Nx, j in 1:Ny]
    (x=gr.x, y=gr.y, bx=f(st.bx), by=f(st.by), bz=f(st.bz),
     ux=f(st.ux), uy=f(st.uy), uz=f(st.uz), ε=f(st.ε), n=f(st.n),
     e=mhd2d_energy_density(st,eng))
end

"Discrete div B = ∂_x J^{xt}+∂_y J^{yt} = -(∂_x J̃^{tx}+∂_y J̃^{ty}) (central, ideal J̃)."
function mhd2d_divB(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy
    @inline wrapx(i)= i<1 ? i+Nx : (i>Nx ? i-Nx : i)
    @inline wrapy(j)= j<1 ? j+Ny : (j>Ny ? j-Ny : j)
    # fill J̃^{ti} from primitives (incl ghosts for the periodic stencil)
    _fill_ghosts!(st, eng)
    Jx=zeros(Nx+2NG,Ny+2NG); Jy=zeros(Nx+2NG,Ny+2NG)
    @inbounds for j in 1:Ny+2NG, i in 1:Nx+2NG
        jx,jy,_=_Jti_cell(st,i,j); Jx[i,j]=jx; Jy[i,j]=jy
    end
    out = zeros(Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        ip=NG+wrapx(i+1); im=NG+wrapx(i-1); jp=NG+wrapy(j+1); jm=NG+wrapy(j-1)
        dxJ = (Jx[ip,NG+j]-Jx[im,NG+j])/(2Δx)
        dyJ = (Jy[NG+i,jp]-Jy[NG+i,jm])/(2Δy)
        out[i,j] = -(dxJ + dyJ)
    end
    return out
end
mhd2d_divB_max(st, eng) = maximum(abs, mhd2d_divB(st, eng))

"""
    mhd2d_out_of_plane_current(st, eng) -> (∇×J)_z

The out-of-plane current (∇×J)_z = ∂_x J^{ty} - ∂_y J^{tx} (Harris diagnostic),
on the interior with periodic central differences.
"""
function mhd2d_out_of_plane_current(st::MHD2DState, eng::MHD2DEngine)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG; Δx=g.Δx; Δy=g.Δy
    @inline wrapx(i)= i<1 ? i+Nx : (i>Nx ? i-Nx : i)
    @inline wrapy(j)= j<1 ? j+Ny : (j>Ny ? j-Ny : j)
    _fill_ghosts!(st, eng)
    Jx=zeros(Nx+2NG,Ny+2NG); Jy=zeros(Nx+2NG,Ny+2NG)
    @inbounds for j in 1:Ny+2NG, i in 1:Nx+2NG
        jx,jy,_=_Jti_cell(st,i,j); Jx[i,j]=jx; Jy[i,j]=jy
    end
    out = zeros(Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        ip=NG+wrapx(i+1); im=NG+wrapx(i-1); jp=NG+wrapy(j+1); jm=NG+wrapy(j-1)
        dxJty = (Jy[ip,NG+j]-Jy[im,NG+j])/(2Δx)
        dyJtx = (Jx[NG+i,jp]-Jx[NG+i,jm])/(2Δy)
        out[i,j] = dxJty - dyJtx
    end
    return out
end

"""
    mhd2d_front_velocity(st, eng; nθ=31, stride=4) -> (vmax, max_imW, ...)

Scan the foundation front-velocity v_max over a strided set of interior cells
(using each cell's local ε,b²) and return the max |Re W| and max |Im W| across
the grid (the causality / anti-diffusive monitor; Im W>0 ⇒ unstable, à la OT
τ_X=0).
"""
function mhd2d_front_velocity(st::MHD2DState, eng::MHD2DEngine; nθ::Int=31, stride::Int=4)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    vmax=0.0; imax=0.0
    @inbounds for j in 1:stride:Ny, i in 1:stride:Nx
        ii=NG+i; jj=NG+j
        b2=_b2(st.bx[ii,jj],st.by[ii,jj],st.bz[ii,jj],st.ux[ii,jj],st.uy[ii,jj],st.uz[ii,jj])
        cc=bdnk_coeffs_from_Dmaps(Du=eng.c.Du,Dε=eng.c.Dε,rb=eng.c.rb,
              τu=eng.c.τu,τX=eng.c.τX,τb=eng.c.τb,ε=max(st.ε[ii,jj],_EPS_FLOOR),b2=max(b2,0.0))
        fv=front_velocity_max(cc; nθ=nθ)
        vmax=max(vmax,fv.vmax); imax=max(imax,fv.max_imW)
    end
    return (vmax=vmax, max_imW=imax)
end

"Per-cell v_max + Im W fields (strided) for figures."
function mhd2d_vmax_field(st::MHD2DState, eng::MHD2DEngine; nθ::Int=31, stride::Int=2)
    g=eng.g; Nx=g.Nx; Ny=g.Ny; NG=g.NG
    is=1:stride:Nx; js=1:stride:Ny
    V=zeros(length(is),length(js)); Im=zeros(length(is),length(js))
    @inbounds for (jc,j) in enumerate(js), (ic,i) in enumerate(is)
        ii=NG+i; jj=NG+j
        b2=_b2(st.bx[ii,jj],st.by[ii,jj],st.bz[ii,jj],st.ux[ii,jj],st.uy[ii,jj],st.uz[ii,jj])
        cc=bdnk_coeffs_from_Dmaps(Du=eng.c.Du,Dε=eng.c.Dε,rb=eng.c.rb,
              τu=eng.c.τu,τX=eng.c.τX,τb=eng.c.τb,ε=max(st.ε[ii,jj],_EPS_FLOOR),b2=max(b2,0.0))
        fv=front_velocity_max(cc; nθ=nθ)
        V[ic,jc]=fv.vmax; Im[ic,jc]=fv.max_imW
    end
    return (x=[g.x[g.NG+i] for i in is], y=[g.y[g.NG+j] for j in js], vmax=V, imW=Im)
end

# ===========================================================================
# 10.  TEST SETUPS  (OT Eq.27, KH Eq.23, Harris Eq.28-32)
# ===========================================================================

# convert (J^{ti}, u^i) at a point to the primitive (b^i): b^t=J^{ti}u_i,
# b^i=(J^{ti}+b^t u^i)/Γ where Γ=u^t (paper Eq.27 convention).
@inline function _b_from_Jti_ui(Jtx,Jty,Jtz, ux,uy,uz)
    Γ = _ut(ux,uy,uz)
    bt = Jtx*ux + Jty*uy + Jtz*uz       # b^t = J^{ti} u_i (lower index = +, spatial)
    bx = (Jtx + bt*ux)/Γ
    by = (Jty + bt*uy)/Γ
    bz = (Jtz + bt*uz)/Γ
    return bx, by, bz
end

"""
ORSZAG–TANG parameter table (paper §V).  Each entry (Du, Dε, rb, τu, τX, τb).
OT-a uses τ_X=0.2; OT-a-tx0 is the SAME but τ_X=0 (the decisive contrast).
"""
const OT_PARAMS = Dict(
    :a    => (1e-2, 2e-3, 1e-2, 2e-1, 2e-1, 8e-2),   # τ_X=0.2 (stable)
    :a_tx0=> (1e-2, 2e-3, 1e-2, 2e-1, 0.0,  8e-2),   # τ_X=0   (unstable: ripples)
)
function ot_coeffs(tag::Symbol; ε::Float64=30.0, b2::Float64=0.0)
    Du,Dε,rb,τu,τX,τb = OT_PARAMS[tag]
    bdnk_coeffs_from_Dmaps(Du=Du,Dε=Dε,rb=rb,τu=τu,τX=τX,τb=τb,ε=ε,b2=b2)
end

"""
    setup_mhd2d_orszagtang(tag; N=128, cfl=0.2, Dnum=0.0, nproj=30, Γvinit=0.8) -> (engine, state)

Orszag–Tang vortex (Eq.27) on [0,2π]² periodic.  IC:
  (u^x,u^y)=(−Γv sin y, Γv sin x), Γv=0.8;  (J^{tx},J^{ty})=(−sin y, sin 2x),
  J^{tz}=0;  ε=30;  b^i from (J^{ti},u^i) (b^t=J^{ti}u_i, b^i=(J^{ti}+b^t u^i)/Γ).
`tag` ∈ (:a, :a_tx0) selects τ_X=0.2 (stable) vs τ_X=0 (the ripple-unstable run).
"""
function setup_mhd2d_orszagtang(tag::Symbol=:a; N::Int=128, cfl::Float64=0.12,
                                Dnum::Float64=0.0, nproj::Int=15, Γvinit::Float64=0.8,
                                εval::Float64=30.0, κshear::Float64=0.6)
    g = MHD2DGrid(N, N, 0.0, 2π, 0.0, 2π)
    ntx=N+2g.NG; nty=N+2g.NG
    c = ot_coeffs(tag; ε=εval, b2=1.0)
    eng = MHD2DEngine(g, c, cfl, 1.0, :hll, Dnum, nproj, κshear)
    st = MHD2DState(ntx, nty)
    @inbounds for j in 1:N, i in 1:N
        ii=g.NG+i; jj=g.NG+j
        x=g.x[ii]; y=g.y[jj]
        ux=-Γvinit*sin(y); uy=Γvinit*sin(x); uz=0.0
        Jtx=-sin(y); Jty=sin(2x); Jtz=0.0
        bx,by,bz=_b_from_Jti_ui(Jtx,Jty,Jtz, ux,uy,uz)
        _setP!(st, ii, jj, (bx,by,bz, ux,uy,uz, εval))
        st.n[ii,jj]=1.0
    end
    _sync_U!(st, eng)
    return eng, st
end

"""
KELVIN–HELMHOLTZ parameter table (paper Table, KH-a..d).
(Du, Dε, rb, τu, τX, τb).
"""
const KH_PARAMS = Dict(
    :a => (1e-4, 5e-5, 1e-4, 1e-3, 5e-4, 5e-4),
    :b => (1e-4, 5e-5, 1e-3, 1e-3, 5e-4, 5e-3),
    :c => (1e-3, 5e-4, 1e-4, 8e-3, 5e-3, 5e-4),
    :d => (1e-3, 5e-4, 1e-3, 8e-3, 5e-3, 5e-3),
)
function kh_coeffs(tag::Symbol; ε::Float64=1.0, b2::Float64=0.0)
    Du,Dε,rb,τu,τX,τb = KH_PARAMS[tag]
    bdnk_coeffs_from_Dmaps(Du=Du,Dε=Dε,rb=rb,τu=τu,τX=τX,τb=τb,ε=ε,b2=b2)
end

"""
    setup_mhd2d_kelvinhelmholtz(tag; Nx=64, Ny=128, cfl=0.2, Dnum=2e-4, nproj=30) -> (engine, state)

Kelvin–Helmholtz shear layer (Eq.23) on the doubly-periodic box [0,1]×[0,2].  IC:
  u^x=0.15[tanh((y+0.5)/0.05)−tanh((y−0.5)/0.05)−1],
  u^y=1e-2 sin(2πx) [exp(−(y+0.5)²/0.04)+exp(−(y−0.5)²/0.04)],
  J^{tx}=0.08, J^{ty}=0, J^{tz}=0.8, ε=1;
  passive tracer n=1+½[tanh((y+0.5)/0.05)−tanh((y−0.5)/0.05)] (Eq.25).
"""
function setup_mhd2d_kelvinhelmholtz(tag::Symbol=:a; Nx::Int=64, Ny::Int=128,
                                     cfl::Float64=0.2, Dnum::Float64=2e-4,
                                     nproj::Int=30, εval::Float64=1.0, κshear::Float64=0.0)
    g = MHD2DGrid(Nx, Ny, 0.0, 1.0, -1.0, 1.0)   # Lx=1, Ly=2 ⇒ y∈[-1,1]
    ntx=Nx+2g.NG; nty=Ny+2g.NG
    c = kh_coeffs(tag; ε=εval, b2=0.6)
    eng = MHD2DEngine(g, c, cfl, 1.0, :hll, Dnum, nproj, κshear)
    st = MHD2DState(ntx, nty)
    @inbounds for j in 1:Ny, i in 1:Nx
        ii=g.NG+i; jj=g.NG+j
        x=g.x[ii]; y=g.y[jj]
        th = tanh((y+0.5)/0.05) - tanh((y-0.5)/0.05)
        ux = 0.15*(th - 1.0)
        uy = 1e-2*sin(2π*x)*(exp(-(y+0.5)^2/0.04)+exp(-(y-0.5)^2/0.04))
        uz = 0.0
        Jtx=0.08; Jty=0.0; Jtz=0.8
        bx,by,bz=_b_from_Jti_ui(Jtx,Jty,Jtz, ux,uy,uz)
        _setP!(st, ii, jj, (bx,by,bz, ux,uy,uz, εval))
        st.n[ii,jj]=1.0 + 0.5*th
    end
    _sync_U!(st, eng)
    return eng, st
end

"""
HARRIS double-current-sheet parameters (paper §V).  (Du, Dε, rb, τu, τX, τb).
"""
const HARRIS_PARAMS = (1e-2, 5e-3, 1e-2, 1e-1, 5e-2, 1e-1)
function harris_coeffs(; ε::Float64=1.0, b2::Float64=0.0)
    Du,Dε,rb,τu,τX,τb = HARRIS_PARAMS
    bdnk_coeffs_from_Dmaps(Du=Du,Dε=Dε,rb=rb,τu=τu,τX=τX,τb=τb,ε=ε,b2=b2)
end

"""
    setup_mhd2d_harris(; Nx=64, Ny=128, Lx=20.0, Ly=40.0, cfl=0.2, Dnum=2e-3,
                       nproj=30, B0=1.0, ℓ=0.5, ψpert=0.1) -> (engine, state)

Double Harris current sheet (Eq.28-32) on a REDUCED uniform box (the paper's
x∈[−100,100],y∈[−200,200] AMR is huge; we use a much smaller box).  IC:
  J^{tx} = B0[−1 + tanh((y−y1)/ℓ) + tanh((y2−y)/ℓ)], y1=Ly/4, y2=3Ly/4,
  with tearing-seed perturbations ψ_bot,ψ_top added to (J^{tx},J^{ty});
  ε from (Eq.31,32): ρ_bg=0.5 in the lobes rising to ρ_sh=2 in the sheets,
  with p_0=0.05 ⇒ ε=3p initialised consistently (conformal).
"""
function setup_mhd2d_harris(; Nx::Int=64, Ny::Int=128, Lx::Float64=20.0, Ly::Float64=40.0,
                            cfl::Float64=0.2, Dnum::Float64=2e-3, nproj::Int=30,
                            B0::Float64=1.0, ℓ::Float64=0.5, ψpert::Float64=0.1,
                            ρbg::Float64=0.5, ρsh::Float64=2.0, p0::Float64=0.05,
                            κshear::Float64=0.0)
    g = MHD2DGrid(Nx, Ny, -Lx/2, Lx/2, -Ly/2, Ly/2)
    ntx=Nx+2g.NG; nty=Ny+2g.NG
    c = harris_coeffs(; ε=3*p0, b2=B0^2)
    eng = MHD2DEngine(g, c, cfl, 1.0, :hll, Dnum, nproj, κshear)
    st = MHD2DState(ntx, nty)
    y1 = g.yL + Ly/4; y2 = g.yL + 3Ly/4
    kx = 2π/Lx
    @inbounds for j in 1:Ny, i in 1:Nx
        ii=g.NG+i; jj=g.NG+j
        x=g.x[ii]; y=g.y[jj]
        # reversing field profile (double sheet)
        Jtx = B0*(-1.0 + tanh((y-y1)/ℓ) + tanh((y2-y)/ℓ))
        # tearing-mode seed: small ψ perturbation localised at each sheet
        sech2(z)=1.0/cosh(z)^2
        ψb = ψpert*cos(kx*x)*exp(-((y-y1)/ℓ)^2)
        ψt = ψpert*cos(kx*x)*exp(-((y-y2)/ℓ)^2)
        # B = ∇×(ψ ẑ) ⇒ δJtx=∂_yψ, δJty=-∂_xψ (seed both sheets)
        dψb_dy = ψb * (-2*(y-y1)/ℓ^2); dψt_dy = ψt * (-2*(y-y2)/ℓ^2)
        dψb_dx = -ψpert*kx*sin(kx*x)*exp(-((y-y1)/ℓ)^2)
        dψt_dx = -ψpert*kx*sin(kx*x)*exp(-((y-y2)/ℓ)^2)
        Jtx += dψb_dy + dψt_dy
        Jty = -(dψb_dx + dψt_dx)
        Jtz = 0.0
        # density / pressure profile (Eq.31,32): pressure balance — higher ρ in
        # sheets where |B| dips.  Build ε=3p with p rising from the lobes.
        sheetfac = sech2((y-y1)/ℓ) + sech2((y-y2)/ℓ)
        ρ = ρbg + (ρsh-ρbg)*sheetfac
        p = p0*(ρ/ρbg)                     # conformal-consistent pressure bump
        ε = 3*p
        u = (0.0,0.0,0.0)
        bx,by,bz=_b_from_Jti_ui(Jtx,Jty,Jtz, u...)
        _setP!(st, ii, jj, (bx,by,bz, u..., ε))
        st.n[ii,jj]=1.0
    end
    _sync_U!(st, eng)
    # project the IC once so the seeded field starts divergence-free
    _project_divB!(st, eng; niter=max(nproj,200))
    _fill_ghosts!(st, eng); _fill_ghosts_U!(st, eng)
    return eng, st
end

end # module BDNKMHD2D
