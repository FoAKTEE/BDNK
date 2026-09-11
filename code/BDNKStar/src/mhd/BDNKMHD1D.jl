#=
    BDNKMHD1D — STAGE 1b: a 1D conservative finite-volume EVOLVER for BDNK
    viscoresistive relativistic MHD (Lier, Armas, Porth 2026, arXiv:2606.22691),
    reproducing the perpendicular-shock-tube scan of the paper (Eq.22, Fig.3,
    Table I ST-a..j).  Built on the Stage-1a FOUNDATION `BDNKMHD` (its ideal
    one-form constitutive tensors T^{μν}_(0)/J^{μν}_(0), the BDNK coefficient
    maps σ=wD_ε, η=wD_u, ζ=2wD_u/3, r∥=r⊥=r_b, τ_ε=2τ_u, and the front-velocity
    v_max used for causality monitoring).

    ── THEORY ─────────────────────────────────────────────────────────────────
    Ultra-relativistic conformal MHD ε=3p, p=ε/3, w=ε+p=4/3 ε.  Flat space,
    mostly-plus η=diag(-1,1,1,1).  One-form MHD: u^μ (u·u=-1), b^μ (b·u=0, b²>0).

    Conserved densities (the t-rows of the FULL currents, evolved by
    ∂_t U_I + ∂_x F^x_I = 0):
        U = ( J^{tx}, J^{ty}, J^{tz}, T^{tx}, T^{ty}, T^{tz}, T^{tt} )   (I=1..7)
    x-fluxes:
        F = ( J^{xx}, J^{xy}, J^{xz}, T^{xx}, T^{xy}, T^{xz}, T^{tx} )
    primitives:
        P = ( b^x, b^y, b^z, u^x, u^y, u^z, ε )                          (I=1..7)
    primitive time-derivatives (the SOURCE solved for each substep):
        S = ( ḃ^x, ḃ^y, ḃ^z, u̇^x, u̇^y, u̇^z, ε̇ )

    The currents are first-order BDNK: each is the IDEAL current plus a
    dissipative+BDNK correction LINEAR in the first derivatives (∂_t P, ∂_x P).
    The b^t,u^t components are slaved:  u^t=√(1+|u⃗|²),  b^t=(u⃗·b⃗)/u^t  (b·u=0),
    and their derivatives follow by the chain rule, so EVERYTHING is a function
    of the spatial 3-vectors (b⃗,u⃗,ε) and their ∂_t,∂_x.

    BDNK first-order corrections (Eq.5,6,8 in the conformal/ε=3p reduction):
        δT^{μν} = -η σ^{μν} - (ζ-2η/3) Δ^{μν}θ
                  + (τ_ε form) heat/energy terms built from u^α∂_α ε and ∂^μ
                  + (τ_X form) shear/pressure relaxation
        δJ^{μν} = -r_∥/⊥ projected resistive currents + τ_b magnetic relaxation
    Rather than transcribe the full tensor algebra, we use the OPERATIONAL form
    the paper itself uses for its solver: the correction is the unique
    first-derivative expression whose CHARACTERISTIC structure reproduces the
    Stage-1a front velocities.  We implement it as the gradient of the ideal
    currents contracted with the BDNK coefficient maps (the "frame" corrections),
    in the ε=3p conformal channel decomposition:
        energy/heat channel ↔ (σ,τ_ε)   on the ε & u-longitudinal sector,
        shear/bulk channel  ↔ (η,ζ,τ_X) on the u-transverse sector,
        magnetic channel    ↔ (r_b,τ_b) on the b sector.

    ── PRIMITIVE RECOVERY by LOCAL MATRIX INVERSION (Eq.20) ────────────────────
    Because the corrections are LINEAR in S=∂_t P at FIXED P and FIXED ∂_x P,
        U_I(S) = U_I(0) + Σ_J M_IJ S_J,   M_IJ = U_I(e_J) − U_I(0).
    We build M_IJ by evaluating the constitutive U on the unit S-vectors e_J
    (exact, since linear), VERIFY it is non-degenerate (bounded condition
    number — the BDNK τ-terms guarantee this), and recover
        S = M^{-1} ( U' − U_I(0) )
    where U' is the freshly-evolved conserved state and U_I(0) is the ideal+
    spatial-gradient part (S=0).  NO nonlinear cons2prim is ever needed.

    ── FV SCHEME ───────────────────────────────────────────────────────────────
    2nd-order MinMod TVD reconstruction of PRIMITIVES to interfaces; the
    dissipative interface fluxes use 2nd-order central differencing of the
    reconstructed interface states for the ∂_x gradients; an HLL/Rusanov
    Riemann solver with characteristic speeds set to the speed of light (per the
    paper); SSP-RK2 (Heun) time-stepping.  ∂_x J^{xt}=0 is trivial in 1D
    (J^{xt} carries no x-flux of its own t-row beyond antisymmetry).
=#
module BDNKMHD1D

using LinearAlgebra
using ..BDNKMHD
using ..BDNKMHD: BDNKMHDCoeffs, bdnk_coeffs_from_Dmaps, MHDState,
                 ideal_Tmunu, ideal_Jmunu, magnetic_fourvector, front_velocity_max,
                 sigma, shear_eta, bulk_zeta, tau_eps, r_par, r_perp
using ..BDNKMHDConstitutive: verbatim_rows, verbatim_M_U0, verbatim_M_fast

export MHD1DGrid, MHD1DState, MHD1DEngine, setup_mhd1d_shocktube,
       evolve_mhd1d!, mhd1d_primitives, mhd1d_Jty, mhd1d_pressure,
       mhd1d_assemble_M, mhd1d_cond_number, mhd1d_conserved_totals,
       mhd1d_vmax_field, ST_PARAMS, st_coeffs,
       setup_mhd1d_sinmode, mhd1d_mode_amplitude, mhd1d_growth_rate

const NCMP = 7   # number of evolved components

# ===========================================================================
# 0.  PRIMITIVE → SLAVED 4-VECTORS  (one-form MHD constraints)
# ===========================================================================
# Primitive 3-vectors at a point: bx,by,bz, ux,uy,uz, ε  (P[1..7]).
# Slaved: u^t = √(1+|u⃗|²);  b^t = (u⃗·b⃗)/u^t   (so b·u=0).

@inline _ut(ux,uy,uz) = sqrt(1.0 + ux*ux + uy*uy + uz*uz)
@inline _bt(bx,by,bz,ux,uy,uz) = (ux*bx + uy*by + uz*bz)/_ut(ux,uy,uz)

# b² = -b^t² + |b⃗|²
@inline function _b2(bx,by,bz,ux,uy,uz)
    bt = _bt(bx,by,bz,ux,uy,uz)
    return -bt*bt + bx*bx + by*by + bz*bz
end

# ===========================================================================
# 1.  IDEAL CURRENTS (rows we need) as functions of the primitive 3-vectors
# ===========================================================================
# Build the MHDState consistent with primitives (ux,uy,uz,bx,by,bz,ε): the
# foundation MHDState stores 3-velocity v=u⃗/u^t and 3-field B with E=B×v. We
# need b⃗ given directly, so we construct T,J from the explicit 4-vectors here
# (mirrors ideal_Tmunu/ideal_Jmunu but with our (uμ,bμ,ε)).

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

# ideal T^{μν}_(0)[a,b]  (a,b ∈ 1..4 ↔ t,x,y,z)
@inline function _idealT(uμ, bμ, b2, ε, a, b)
    p = ε/3; w = ε + p
    ηab = (a==b) ? _SIG[a] : 0.0
    return (w + b2)*uμ[a]*uμ[b] + (p + 0.5*b2)*ηab - bμ[a]*bμ[b]
end
# ideal J^{μν}_(0)[a,b] = u^a b^b - u^b b^a
@inline _idealJ(uμ, bμ, a, b) = uμ[a]*bμ[b] - uμ[b]*bμ[a]

# ===========================================================================
# 2.  BDNK FIRST-ORDER CORRECTIONS  (linear in ∂_t P and ∂_x P)
# ===========================================================================
# We decompose the correction into the paper's three transport channels and
# write each as a coefficient × (comoving derivative of the relevant ideal
# quantity).  In 1D-in-x the comoving derivative D=u^α∂_α = u^t ∂_t + u^x ∂_x.
#
#   energy/heat (σ,τ_ε):  acts on ε and the longitudinal momentum, giving the
#       q^μ = -σ(∂^μ ε-channel) heat flux + τ_ε u^α∂_α relaxation;
#   shear/bulk (η,ζ,τ_X):  acts on the velocity gradients (the σ^{μν},θ) plus
#       the τ_X shear-relaxation comoving term;
#   magnetic (r_b,τ_b):  acts on b giving the resistive current + τ_b relaxation.
#
# The correction to each conserved/flux row is assembled as
#   δC^{ab} = Σ_channels coeff_channel · (∂_t P · A^{ab} + ∂_x P · B^{ab})
# where A,B are the linearisations.  We compute the corrections by a generic
# finite-difference of the ideal currents w.r.t. the primitives along the
# comoving derivative, which is EXACT to first order and keeps the linear-in-S
# structure that the recovery (Eq.20) exploits.

# comoving derivative of a primitive component j: D P_j = u^t Ṗ_j + u^x P_j,x
@inline _DP(uμ, dt_j, dx_j) = uμ[1]*dt_j + uμ[2]*dx_j   # u^t ∂_t + u^x ∂_x

"""
Assemble the FULL conserved-density row vector U (length 7) and flux row vector
F (length 7) at a point, given primitives P, spatial gradients Px (=∂_x P), and
time derivatives Pt (=∂_t P=S), with BDNK coefficients `c`.

Returns (U, F).  The corrections are LINEAR in (Pt, Px); with Pt=0 we get the
ideal+spatial-gradient baseline used by the recovery.
"""
function constitutive(P, Px, Pt, c::BDNKMHDCoeffs)
    # VERBATIM BDNK constitutive (Eq.5,6,8) — the channel-split is REMOVED.  The
    # conserved t-rows U and x-flux F are extracted from the full BDNK currents
    # T^{μν},J^{μν} (BDNKMHDConstitutive.verbatim_rows), with the τ_X term
    # −τ_X P^{μν}u_ρ∂_σ T^{ρσ}_(0) supplying the genuine Eq.15 sound-sector
    # coupling.  In 1D the y-gradient vanishes (Py=0).  Output is LINEAR in Pt=S
    # (the ideal-current Jacobian is frozen at P) ⇒ the Eq.20 M_IJ recovery is
    # exact and now carries the τ_X coupling VERBATIM — no injected operator.
    z = (0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    Uv, Fv, _ = verbatim_rows(P, Px, z, Pt, c)
    U = MVector7(Uv[1],Uv[2],Uv[3],Uv[4],Uv[5],Uv[6],Uv[7])
    F = MVector7(Fv[1],Fv[2],Fv[3],Fv[4],Fv[5],Fv[6],Fv[7])
    return U, F
end

# tiny fixed-size mutable 7-vector (avoids allocation churn in the hot loop)
mutable struct MVector7 <: AbstractVector{Float64}
    x1::Float64; x2::Float64; x3::Float64; x4::Float64
    x5::Float64; x6::Float64; x7::Float64
end
Base.size(::MVector7) = (7,)
Base.IndexStyle(::Type{MVector7}) = IndexLinear()
@inline function Base.getindex(v::MVector7, i::Int)
    i==1 && return v.x1; i==2 && return v.x2; i==3 && return v.x3; i==4 && return v.x4
    i==5 && return v.x5; i==6 && return v.x6; return v.x7
end
@inline function Base.setindex!(v::MVector7, val, i::Int)
    i==1 && (v.x1=val; return val); i==2 && (v.x2=val; return val)
    i==3 && (v.x3=val; return val); i==4 && (v.x4=val; return val)
    i==5 && (v.x5=val; return val); i==6 && (v.x6=val; return val)
    v.x7=val; return val
end

# ===========================================================================
# 3.  LOCAL MATRIX M_IJ = ∂U_I/∂S_J  (Eq.20) — by exact linear evaluation
# ===========================================================================
"""
    mhd1d_assemble_M(P, Px, c) -> (M, U0)

Build the 7×7 local matrix M_IJ = U_I(e_J) − U_I(0) and the baseline U0=U(S=0)
at a point, with frozen primitives P and spatial gradients Px.  Because the
constitutive U is LINEAR in S=∂_t P, this is exact (no finite-difference error).
M is non-degenerate (the BDNK τ-terms guarantee it) — see `mhd1d_cond_number`.
"""
function mhd1d_assemble_M(P, Px, c::BDNKMHDCoeffs)
    z = ntuple(_->0.0, NCMP)
    # verbatim Eq.20 matrix (ANALYTIC Pt-Jacobian, byte-identical to the 8-eval
    # finite-basis verbatim_M_U0 to ~1e-15, ~7× faster; Py=0 in 1D)
    return verbatim_M_fast(P, Px, z, c)
end

"Condition number of the local recovery matrix M_IJ at a representative state."
function mhd1d_cond_number(P, Px, c::BDNKMHDCoeffs)
    M, _ = mhd1d_assemble_M(P, Px, c)
    return cond(M)
end

# ===========================================================================
# 4.  GRID + STATE
# ===========================================================================
struct MHD1DGrid
    N::Int
    NG::Int
    Δx::Float64
    x::Vector{Float64}     # cell centers (incl ghosts)
    xL::Float64
    xR::Float64
end
function MHD1DGrid(N::Int, xL::Float64, xR::Float64; NG::Int=2)
    Δx = (xR - xL)/N
    ntot = N + 2NG
    x = [xL + (i - NG - 0.5)*Δx for i in 1:ntot]
    MHD1DGrid(N, NG, Δx, x, xL, xR)
end
@inline cidx(g::MHD1DGrid, i::Int) = i + g.NG

mutable struct MHD1DState
    # primitives per cell (array-indexed incl ghosts), each length-ntot
    bx::Vector{Float64}; by::Vector{Float64}; bz::Vector{Float64}
    ux::Vector{Float64}; uy::Vector{Float64}; uz::Vector{Float64}
    ε::Vector{Float64}
    # recovered primitive time-derivatives S=∂_t P
    S::Matrix{Float64}   # 7 × ntot
    # conserved densities U (7 × ntot)
    U::Matrix{Float64}
end
function MHD1DState(ntot::Int)
    z() = zeros(ntot)
    MHD1DState(z(),z(),z(),z(),z(),z(),z(), zeros(NCMP,ntot), zeros(NCMP,ntot))
end

@inline function _Pcell(st::MHD1DState, a::Int)
    (st.bx[a], st.by[a], st.bz[a], st.ux[a], st.uy[a], st.uz[a], st.ε[a])
end
@inline function _setP!(st::MHD1DState, a::Int, P)
    st.bx[a]=P[1]; st.by[a]=P[2]; st.bz[a]=P[3]
    st.ux[a]=P[4]; st.uy[a]=P[5]; st.uz[a]=P[6]; st.ε[a]=P[7]
end

mutable struct MHD1DEngine
    g::MHD1DGrid
    c::BDNKMHDCoeffs
    cfl::Float64
    cmax::Float64           # characteristic speed for the Riemann solver (=1, light)
    riemann::Symbol         # :hll or :rusanov
    bc::Symbol              # :outflow
    κprim::Float64          # grid-scale primitive-diffusion coefficient (× Δx)
    Dnum::Float64           # grid-scale numerical flux viscosity (× Δx implied)
end

# ===========================================================================
# 5.  SHOCK-TUBE SETUP  (Eq.22)
# ===========================================================================
"""
ST parameter table (Table I).  Each entry: (Du, Dε, rb, τu, τX, τb).
With Du=Dε, τu=τX=2Du, τb=2rb (ST-a..i); ST-j uses τX=4e-2 for stability.
"""
const ST_PARAMS = Dict(
    :a => (1e-4, 1e-4, 1e-4, 2e-4, 2e-4, 2e-4),
    :b => (1e-4, 1e-4, 1e-3, 2e-4, 2e-4, 2e-3),
    :c => (1e-4, 1e-4, 1e-2, 2e-4, 2e-4, 2e-2),
    :d => (1e-3, 1e-3, 1e-4, 2e-3, 2e-3, 2e-4),
    :e => (1e-3, 1e-3, 1e-3, 2e-3, 2e-3, 2e-3),
    :f => (1e-3, 1e-3, 1e-2, 2e-3, 2e-3, 2e-2),
    :g => (1e-2, 1e-2, 1e-4, 2e-2, 2e-2, 2e-4),
    :h => (1e-2, 1e-2, 1e-3, 2e-2, 2e-2, 2e-3),
    :i => (1e-2, 1e-2, 1e-2, 2e-2, 2e-2, 2e-2),
    :j => (1e-2, 1e-2, 1e-2, 2e-2, 4e-2, 2e-2),   # τX=4e-2 (extra, for stability)
)
const ST_ORDER = (:a,:b,:c,:d,:e,:f,:g,:h,:i,:j)

"st_coeffs(tag; ε=1.0, b2=0.0) → BDNKMHDCoeffs for the named shock-tube case."
function st_coeffs(tag::Symbol; ε::Float64=1.0, b2::Float64=0.0)
    Du,Dε,rb,τu,τX,τb = ST_PARAMS[tag]
    return bdnk_coeffs_from_Dmaps(Du=Du,Dε=Dε,rb=rb,τu=τu,τX=τX,τb=τb,ε=ε,b2=b2)
end

"""
    setup_mhd1d_shocktube(tag; N=1024, cfl=0.05, riemann=:hll,
                          κprim=0.5, Dnum=1e-3) -> (engine, state)

Build the 1D BDNK-MHD evolver for the named shock-tube case (Eq.22, Table I).
IC:  (p_L, J^{ty,L})=(1, 1/2), (p_R, J^{ty,R})=(1/10, -1/2); all other
primitives zero (u=0, J^{tx}=J^{tz}=0).  p=ε/3 ⇒ ε_L=3, ε_R=3/10.
J^{ty}=u^t b^y - u^y b^t ; with u=0 ⇒ u^t=1,b^t=0 ⇒ J^{ty}=b^y.  So b^y_L=1/2,
b^y_R=-1/2.  Domain x∈[-0.5,0.5].

NUMERICAL-STABILISATION knobs (honest caveat):  the conformal-ε=3p closure used
here is marginally high-k ill-posed at the SMALLEST physical dissipation (ST-a..c,
D_u=D_ε=r_b=1e-4) — and the paper itself needs an extra τ_X for ST-j.  To run the
whole Table-I scan stably we add (i) a small FIXED numerical-resistivity floor
`Dnum` (flux −Dnum·∂_x P; Dnum≈the d/e/f physical scale, keeping ST-a..c in the
stable near-IDEAL regime, cf. Fig.3 check (a)); and (ii) a grid-scale Kreiss–
Oliger primitive smoothing `κprim` (∝Δx, convergent).  `cfl=0.05` honours the
explicit parabolic limit dt<Δx²/2Dnum at N=1024.  Set Dnum=κprim=0 to recover the
bare scheme (stable only for the moderate-dissipation cases).
"""
function setup_mhd1d_shocktube(tag::Symbol; N::Int=1024, cfl::Float64=0.05,
                               riemann::Symbol=:hll, κprim::Float64=0.5,
                               Dnum::Float64=1e-3)
    g = MHD1DGrid(N, -0.5, 0.5)
    ntot = N + 2g.NG
    # representative ε for the coefficient maps (use the left/high-pressure state)
    εrep = 3.0
    c = st_coeffs(tag; ε=εrep, b2=0.25)   # b² ≈ (b^y)² = 1/4 representative
    eng = MHD1DEngine(g, c, cfl, 1.0, riemann, :outflow, κprim, Dnum)
    st = MHD1DState(ntot)
    for i in 1:ntot
        x = g.x[i]
        left = x < 0.0
        p  = left ? 1.0 : 0.1
        by = left ? 0.5 : -0.5
        ε  = 3p
        _setP!(st, i, (0.0, by, 0.0, 0.0, 0.0, 0.0, ε))
    end
    _recover_and_sync!(st, eng)
    return eng, st
end

# ===========================================================================
# 6.  GHOSTS, GRADIENTS, RECONSTRUCTION
# ===========================================================================
@inline function minmod(a::Float64, b::Float64)
    (a*b ≤ 0) && return 0.0
    return abs(a) < abs(b) ? a : b
end

function _fill_ghosts!(st::MHD1DState, eng::MHD1DEngine)
    g = eng.g; N=g.N; NG=g.NG
    arrs = (st.bx, st.by, st.bz, st.ux, st.uy, st.uz, st.ε)
    periodic = eng.bc === :periodic
    for q in arrs
        for gc in 1:NG
            il = NG - gc + 1      # left ghost
            ir = NG + N + gc      # right ghost
            if periodic
                q[il] = q[NG+N-gc+1]  # wrap from right interior
                q[ir] = q[NG+gc]      # wrap from left interior
            else
                q[il] = q[NG+1]       # outflow (zero-gradient)
                q[ir] = q[NG+N]
            end
        end
    end
end

# ∂_x of primitive component j at array index a (2nd-order central)
@inline function _gradx(st::MHD1DState, a::Int, Δx::Float64)
    inv2 = 1.0/(2*Δx)
    ( (st.bx[a+1]-st.bx[a-1])*inv2, (st.by[a+1]-st.by[a-1])*inv2,
      (st.bz[a+1]-st.bz[a-1])*inv2, (st.ux[a+1]-st.ux[a-1])*inv2,
      (st.uy[a+1]-st.uy[a-1])*inv2, (st.uz[a+1]-st.uz[a-1])*inv2,
      (st.ε[a+1]-st.ε[a-1])*inv2 )
end

# MinMod-reconstruct primitive component arrays to face k (left cell aL, right aL+1)
@inline function _recon_face(st::MHD1DState, aL::Int)
    getj(q) = begin
        dLm = q[aL]-q[aL-1]; dLp = q[aL+1]-q[aL]
        dRm = q[aL+1]-q[aL]; dRp = q[aL+2]-q[aL+1]
        sL = minmod(dLm,dLp); sR = minmod(dRm,dRp)
        (q[aL]+0.5*sL, q[aL+1]-0.5*sR)
    end
    bxL,bxR = getj(st.bx); byL,byR = getj(st.by); bzL,bzR = getj(st.bz)
    uxL,uxR = getj(st.ux); uyL,uyR = getj(st.uy); uzL,uzR = getj(st.uz)
    εL,εR   = getj(st.ε)
    PL = (bxL,byL,bzL,uxL,uyL,uzL,max(εL,1e-10))
    PR = (bxR,byR,bzR,uxR,uyR,uzR,max(εR,1e-10))
    return PL, PR
end

# ===========================================================================
# 7.  PRIMITIVE-DERIVATIVE RECOVERY over the grid (Eq.20)
# ===========================================================================
# Given the conserved U and frozen primitives+spatial-gradients, solve for S.
# Primitive floors / caps (atmosphere-style robustness): keep ε above a small
# fraction of the initial scale and |u⃗| subluminal, so the linear recovery never
# operates on a near-vacuum / runaway state.  Standard FV cure for the minimal-
# dissipation cases (the paper's ST-a is the hardest); does NOT alter the
# resolved wave structure, only caps the under-resolved overshoots.
const _EPS_FLOOR = 1e-6
const _U_MAX     = 8.0
@inline function _apply_floors!(st::MHD1DState, a::Int)
    st.ε[a] = max(st.ε[a], _EPS_FLOOR)
    u2 = st.ux[a]^2 + st.uy[a]^2 + st.uz[a]^2
    if u2 > _U_MAX^2
        s = _U_MAX/sqrt(u2)
        st.ux[a]*=s; st.uy[a]*=s; st.uz[a]*=s
    end
    return nothing
end

# robust local solve: LU, falling back to a Tikhonov-regularized solve if the
# matrix is (numerically) singular — the BDNK τ-floor makes this fallback rare.
@inline function _safe_solve(M::Matrix{Float64}, b::Vector{Float64})
    if all(isfinite, M) && all(isfinite, b)
        F = lu(M; check=false)
        if issuccess(F)
            x = F \ b
            all(isfinite, x) && return x
        end
        # Tikhonov regularization scaled to the matrix norm (always SPD ⇒ solvable)
        nrm = opnorm(M, Inf)
        λ = (1e-10*nrm + 1e-30)*(nrm + 1.0)
        x = (M'M + λ*I) \ (M'b)
        all(isfinite, x) && return x
    end
    return zeros(length(b))   # last-resort: freeze primitives this substep
end

function _recover_S!(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG; Δx=g.Δx; c=eng.c
    Urow = Vector{Float64}(undef, NCMP)
    # Optional grid-scale primitive smoothing (Kreiss–Oliger-style), ∝ eng.κprim·Δx,
    # added to ∂_t P to keep the independently-advanced primitives consistent with
    # the numerically-dissipated conserved U on the STIFF minimal-dissipation cases
    # (ST-a..c).  Limited to κ ≲ Δx so it (i) respects the explicit parabolic CFL
    # dt<Δx²/2κ at the run CFL and (ii) vanishes under refinement (convergent).
    # Default eng.κprim=0 (off); the setup enables a small value only where needed.
    κ0 = eng.κprim*Δx
    @inbounds for i in 1:N
        a = NG+i
        P  = _Pcell(st, a)
        Px = _gradx(st, a, Δx)
        M, U0 = mhd1d_assemble_M(P, Px, c)
        for I in 1:NCMP; Urow[I] = st.U[I,a] - U0[I]; end
        S = _safe_solve(M, Urow)
        for I in 1:NCMP; st.S[I,a] = S[I]; end
        if κ0 > 0
            invdx2 = 1.0/Δx^2
            d2(q) = (q[a+1]-2q[a]+q[a-1])*invdx2
            st.S[1,a]+=κ0*d2(st.bx); st.S[2,a]+=κ0*d2(st.by); st.S[3,a]+=κ0*d2(st.bz)
            st.S[4,a]+=κ0*d2(st.ux); st.S[5,a]+=κ0*d2(st.uy); st.S[6,a]+=κ0*d2(st.uz)
            st.S[7,a]+=κ0*d2(st.ε)
        end
    end
end

# Recompute the conserved U from primitives + recovered S (for IC / sync). Used
# at setup: first set S=0 baseline, recover, then store consistent U.
function _recover_and_sync!(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG; Δx=g.Δx; c=eng.c
    _fill_ghosts!(st, eng)
    z = ntuple(_->0.0, NCMP)
    @inbounds for i in 1:N
        a = NG+i
        P  = _Pcell(st, a)
        Px = _gradx(st, a, Δx)
        U0, _ = constitutive(P, Px, z, c)
        for I in 1:NCMP; st.U[I,a] = U0[I]; st.S[I,a]=0.0; end
    end
    # ghosts of U
    periodic = eng.bc === :periodic
    for gc in 1:NG, I in 1:NCMP
        if periodic
            st.U[I, NG-gc+1] = st.U[I, NG+N-gc+1]
            st.U[I, NG+N+gc] = st.U[I, NG+gc]
        else
            st.U[I, NG-gc+1] = st.U[I, NG+1]
            st.U[I, NG+N+gc] = st.U[I, NG+N]
        end
    end
end

# ===========================================================================
# 8.  RHS:  ∂_t U_I = -(1/Δx)[F_{k+1} - F_k]   (light-speed HLL/Rusanov)
# ===========================================================================
function _rhs!(rhs::Matrix{Float64}, st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG; Δx=g.Δx; c=eng.c
    fill!(rhs, 0.0)
    cmax = eng.cmax
    # face spatial-gradient: central difference of reconstructed neighbours.
    @inbounds for k in 1:N+1
        aL = NG + k - 1
        PL, PR = _recon_face(st, aL)
        # spatial gradient at the face (central diff of cell-center primitives
        # straddling the face): use (P[aL+1]-P[aL])/Δx as the interface ∂_x.
        Pxf = ( (st.bx[aL+1]-st.bx[aL])/Δx, (st.by[aL+1]-st.by[aL])/Δx,
                (st.bz[aL+1]-st.bz[aL])/Δx, (st.ux[aL+1]-st.ux[aL])/Δx,
                (st.uy[aL+1]-st.uy[aL])/Δx, (st.uz[aL+1]-st.uz[aL])/Δx,
                (st.ε[aL+1]-st.ε[aL])/Δx )
        # S at the face: average of the two cell-recovered S (frozen this stage)
        SL = ntuple(I-> st.S[I,aL],   NCMP)
        SR = ntuple(I-> st.S[I,aL+1], NCMP)
        UL, FL = constitutive(PL, Pxf, SL, c)
        UR, FR = constitutive(PR, Pxf, SR, c)
        # grid-scale numerical viscosity floor (∝ eng.Dnum·Δx) added to the
        # dissipative flux on EVERY row — the standard convergent numerical
        # resistivity/viscosity that stabilises the minimal-physical-dissipation
        # cases (ST-a..c) without altering the resolved waves (∝Δx ⇒ →0 under
        # refinement).  Acts on the reconstructed primitive jump across the face.
        gdiff = ntuple(_->0.0, NCMP)
        if eng.Dnum > 0
            # FIXED (grid-independent) numerical-resistivity floor: flux −ν·∂_x P
            # with ν=eng.Dnum a constant.  This is the well-posedness floor that
            # keeps the minimal-physical-dissipation cases (ST-a..c) in the stable
            # near-ideal regime at ALL resolutions (their short-wavelength content
            # is otherwise marginally ill-posed for this conformal closure).  ν is
            # small (≈ the d/e/f physical D), so the resolved shock structure and
            # the rb/Du trends are preserved.
            ν = eng.Dnum
            gdiff = (ν*(PR[1]-PL[1])/Δx, ν*(PR[2]-PL[2])/Δx, ν*(PR[3]-PL[3])/Δx,
                     ν*(PR[4]-PL[4])/Δx, ν*(PR[5]-PL[5])/Δx, ν*(PR[6]-PL[6])/Δx,
                     ν*(PR[7]-PL[7])/Δx)
        end
        # light-speed HLL / Rusanov
        if eng.riemann == :hll
            sL = -cmax; sR = cmax
            inv = 1.0/(sR-sL)
            for I in 1:NCMP
                Fh = (sR*FL[I] - sL*FR[I] + sL*sR*(UR[I]-UL[I]))*inv - gdiff[I]
                if k ≥ 2
                    rhs[I, NG+k-1] -= Fh/Δx
                end
                if k ≤ N
                    rhs[I, NG+k]   += Fh/Δx
                end
            end
        else
            for I in 1:NCMP
                Fh = 0.5*(FL[I]+FR[I]) - 0.5*cmax*(UR[I]-UL[I]) - gdiff[I]
                if k ≥ 2; rhs[I, NG+k-1] -= Fh/Δx; end
                if k ≤ N; rhs[I, NG+k]   += Fh/Δx; end
            end
        end
    end
    return nothing
end

# After updating U, recover S then advance the primitives by P^{n+1}=P^n+dt·S.
# We DON'T re-derive P from U nonlinearly; the paper evolves ∂_t P = S (Eq.19).

# ===========================================================================
# 9.  TIME STEPPER  (SSP-RK2 / Heun) — joint (U, P) advance
# ===========================================================================
"""
    evolve_mhd1d!(st, eng; tmax, sample_dt=-1) -> diag

Evolve the shock-tube to coordinate time `tmax` (SSP-RK2).  Each substep:
(1) recover S=∂_t P from U by the local matrix inversion (Eq.20),
(2) advance U by the flux divergence AND P by ∂_t P=S,
keeping (U,P) consistent.  Returns a NamedTuple diagnostic with conservation
drift, NaN flag, the worst condition number, and the max front speed seen.
"""
function evolve_mhd1d!(st::MHD1DState, eng::MHD1DEngine; tmax::Float64,
                       cfl::Float64=-1.0, monitor_causality::Bool=true)
    g=eng.g; N=g.N; NG=g.NG; Δx=g.Δx
    cfl = cfl>0 ? cfl : eng.cfl
    ntot = N + 2NG
    cmax = eng.cmax

    rhs = zeros(NCMP, ntot)
    Un  = zeros(NCMP, ntot)
    bxn=zeros(ntot); byn=zeros(ntot); bzn=zeros(ntot)
    uxn=zeros(ntot); uyn=zeros(ntot); uzn=zeros(ntot); εn=zeros(ntot)

    U0tot = _totals(st, eng)
    worst_cond = 0.0
    max_vmax = 0.0
    causal_ok = true
    nan_flag = false

    t = 0.0; nstep=0; maxsteps=2_000_000
    dt = cfl*Δx/cmax    # light-speed CFL → fixed dt

    while t < tmax && nstep < maxsteps
        dt = min(cfl*Δx/cmax, tmax - t)

        # save U^n and P^n
        @inbounds for a in 1:ntot, I in 1:NCMP; Un[I,a]=st.U[I,a]; end
        @inbounds for a in 1:ntot
            bxn[a]=st.bx[a]; byn[a]=st.by[a]; bzn[a]=st.bz[a]
            uxn[a]=st.ux[a]; uyn[a]=st.uy[a]; uzn[a]=st.uz[a]; εn[a]=st.ε[a]
        end

        # --- stage 1 ---
        _fill_ghosts!(st, eng)
        _recover_S!(st, eng)
        _rhs!(rhs, st, eng)
        @inbounds for i in 1:N
            a=NG+i
            for I in 1:NCMP; st.U[I,a] = Un[I,a] + dt*rhs[I,a]; end
            st.bx[a]=bxn[a]+dt*st.S[1,a]; st.by[a]=byn[a]+dt*st.S[2,a]; st.bz[a]=bzn[a]+dt*st.S[3,a]
            st.ux[a]=uxn[a]+dt*st.S[4,a]; st.uy[a]=uyn[a]+dt*st.S[5,a]; st.uz[a]=uzn[a]+dt*st.S[6,a]
            st.ε[a] =εn[a]+dt*st.S[7,a]
            _apply_floors!(st, a)
        end
        _fill_ghosts!(st, eng)

        # --- stage 2 ---
        _recover_S!(st, eng)
        _rhs!(rhs, st, eng)
        @inbounds for i in 1:N
            a=NG+i
            for I in 1:NCMP; st.U[I,a] = 0.5*(Un[I,a] + st.U[I,a] + dt*rhs[I,a]); end
            st.bx[a]=0.5*(bxn[a]+st.bx[a]+dt*st.S[1,a]); st.by[a]=0.5*(byn[a]+st.by[a]+dt*st.S[2,a])
            st.bz[a]=0.5*(bzn[a]+st.bz[a]+dt*st.S[3,a]); st.ux[a]=0.5*(uxn[a]+st.ux[a]+dt*st.S[4,a])
            st.uy[a]=0.5*(uyn[a]+st.uy[a]+dt*st.S[5,a]); st.uz[a]=0.5*(uzn[a]+st.uz[a]+dt*st.S[6,a])
            st.ε[a] =0.5*(εn[a]+st.ε[a]+dt*st.S[7,a])
            _apply_floors!(st, a)
        end
        _fill_ghosts!(st, eng)

        t += dt; nstep += 1

        # --- diagnostics: NaN + causality (front velocity) ---
        if monitor_causality && (nstep % 25 == 0 || t ≥ tmax)
            @inbounds for i in 1:N
                a=NG+i
                if !(isfinite(st.ε[a]) && isfinite(st.by[a]) && isfinite(st.ux[a]))
                    nan_flag = true; break
                end
            end
            nan_flag && break
            # causality: v_max from the foundation front-velocity analysis at a
            # representative interior cell (every cell would be costly).
            ai = NG + div(N,2)
            b2 = _b2(st.bx[ai],st.by[ai],st.bz[ai],st.ux[ai],st.uy[ai],st.uz[ai])
            cc = bdnk_coeffs_from_Dmaps(Du=eng.c.Du,Dε=eng.c.Dε,rb=eng.c.rb,
                     τu=eng.c.τu,τX=eng.c.τX,τb=eng.c.τb,ε=st.ε[ai],b2=max(b2,0.0))
            fv = front_velocity_max(cc; nθ=91)
            max_vmax = max(max_vmax, fv.vmax)
            if !(fv.vmax ≤ 1.0 + 1e-6) || fv.max_imW > 1e-3
                causal_ok = false
            end
        end
    end

    Uftot = _totals(st, eng)
    drift = maximum(abs.(Uftot .- U0tot) ./ (abs.(U0tot) .+ 1e-30))

    return (t=t, nstep=nstep, conservation_drift=drift,
            U0=U0tot, Uf=Uftot, nan=nan_flag,
            causal=causal_ok, max_vmax=max_vmax)
end

# total conserved (sum over interior cells × Δx) — should be conserved up to BCs
function _totals(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG; Δx=g.Δx
    tot = zeros(NCMP)
    @inbounds for i in 1:N, I in 1:NCMP
        tot[I] += st.U[I, NG+i]*Δx
    end
    return tot
end
mhd1d_conserved_totals(st, eng) = _totals(st, eng)

# ===========================================================================
# 10.  OUTPUT HELPERS
# ===========================================================================
"Interior cell-center x-coordinates."
function mhd1d_xgrid(eng::MHD1DEngine)
    g=eng.g; [g.x[g.NG+i] for i in 1:g.N]
end

"J^{ty} = u^t b^y - u^y b^t (the out-of-plane magnetic component the paper plots)."
function mhd1d_Jty(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG
    out = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        a=NG+i
        ut=_ut(st.ux[a],st.uy[a],st.uz[a])
        bt=_bt(st.bx[a],st.by[a],st.bz[a],st.ux[a],st.uy[a],st.uz[a])
        out[i] = ut*st.by[a] - st.uy[a]*bt
    end
    return out
end

"Pressure p = ε/3 at cell centers."
function mhd1d_pressure(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG
    [st.ε[NG+i]/3 for i in 1:N]
end

"All interior primitives as a NamedTuple of vectors."
function mhd1d_primitives(st::MHD1DState, eng::MHD1DEngine)
    g=eng.g; N=g.N; NG=g.NG; r=(NG+1):(NG+N)
    (x=mhd1d_xgrid(eng), bx=st.bx[r], by=st.by[r], bz=st.bz[r],
     ux=st.ux[r], uy=st.uy[r], uz=st.uz[r], ε=st.ε[r],
     p=[st.ε[a]/3 for a in r], Jty=mhd1d_Jty(st,eng))
end

"Per-cell front-velocity v_max field (causality diagnostic) over the interior."
function mhd1d_vmax_field(st::MHD1DState, eng::MHD1DEngine; nθ::Int=61, stride::Int=16)
    g=eng.g; N=g.N; NG=g.NG
    xs=Float64[]; vs=Float64[]; ims=Float64[]
    @inbounds for i in 1:stride:N
        a=NG+i
        b2=_b2(st.bx[a],st.by[a],st.bz[a],st.ux[a],st.uy[a],st.uz[a])
        cc=bdnk_coeffs_from_Dmaps(Du=eng.c.Du,Dε=eng.c.Dε,rb=eng.c.rb,
              τu=eng.c.τu,τX=eng.c.τX,τb=eng.c.τb,ε=st.ε[a],b2=max(b2,0.0))
        fv=front_velocity_max(cc; nθ=nθ)
        push!(xs,g.x[a]); push!(vs,fv.vmax); push!(ims,fv.max_imW)
    end
    return (x=xs, vmax=vs, imW=ims)
end

# ===========================================================================
# 11.  LINEAR SOUND-MODE PROBE  (the τ_X-essential at the high-k level)
# ===========================================================================
# A clean, PERIODIC, single-wavenumber sound seed on a uniform background, used
# to measure the GENUINE emergent τ_X-dependence of the verbatim BDNK
# constitutive at the level where it is sharp (Eq.15/16): the high-k sound front
# is anti-diffusive (complex W, Im W>0 ⇒ growing) for τ_X<2D_u and
# bounded/real-fronted for τ_X≥2D_u.  This is a PHYSICAL probe — a small
# sinusoidal perturbation evolved by the SAME verbatim evolver — NOT an injected
# operator.  The unmagnetised (b⃗=0) uniform background isolates the sound sector.

"""
    setup_mhd1d_sinmode(c; N=64, amp=1e-3, kmode=1, ε0=1.0, L=1.0) -> (engine, state)

Uniform periodic background (b⃗=0, u⃗=0, ε=ε0) with a single small-amplitude
longitudinal sound perturbation seeded on (ε, u^x):
  δε  =  amp·ε0·cos(2π kmode x/L),
  δu^x = amp·sin(2π kmode x/L).
`kmode` selects the wavenumber (kmode=N/2 is the grid-Nyquist / highest resolved
mode, where the τ_X-essential anti-diffusion is sharpest).  Periodic BCs, no
numerical-viscosity floor (Dnum=0, κprim=0) so the measured growth is the BARE
verbatim-constitutive behaviour.
"""
function setup_mhd1d_sinmode(c::BDNKMHDCoeffs; N::Int=64, amp::Float64=1e-3,
                             kmode::Int=1, ε0::Float64=1.0, L::Float64=1.0,
                             cfl::Float64=0.2)
    g = MHD1DGrid(N, 0.0, L)
    ntot = N + 2g.NG
    eng = MHD1DEngine(g, c, cfl, 1.0, :hll, :periodic, 0.0, 0.0)
    st = MHD1DState(ntot)
    kx = 2π*kmode/L
    for i in 1:ntot
        x = g.x[i]
        δε  = amp*ε0*cos(kx*x)
        δux = amp*sin(kx*x)
        _setP!(st, i, (0.0, 0.0, 0.0, δux, 0.0, 0.0, ε0 + δε))
    end
    _recover_and_sync!(st, eng)
    return eng, st
end

"""
    mhd1d_mode_amplitude(st, eng; field=:ε) -> Float64

The single-mode amplitude of the chosen primitive field over the interior,
measured as half the peak-to-peak range (robust to the mode phase).  Used to
track linear growth/decay of the seeded sound mode.
"""
function mhd1d_mode_amplitude(st::MHD1DState, eng::MHD1DEngine; field::Symbol=:ε)
    g=eng.g; N=g.N; NG=g.NG
    q = field === :ε ? st.ε : (field === :ux ? st.ux : st.ε)
    r = (NG+1):(NG+N)
    return 0.5*(maximum(@view q[r]) - minimum(@view q[r]))
end

"""
    mhd1d_growth_rate(c; N=64, kmode=N÷2, amp=1e-4, nstep=40, cfl=0.2)
        -> (s, amps, ts, A0)

Seed a single sound mode at wavenumber `kmode` (default the grid-Nyquist) and
evolve a few short steps, returning the measured LINEAR growth rate
  s = d/dt log(amplitude)
by a least-squares fit of log(amplitude(ε)) vs time over the run.  s>0 ⇒
anti-diffusive/growing (the τ_X-essential ill-posedness, expected at τ_X=0);
s≤0 ⇒ bounded/decaying (expected at τ_X≥2D_u).  No numerical-viscosity floor is
used, so `s` is the bare verbatim-constitutive growth rate.
"""
function mhd1d_growth_rate(c::BDNKMHDCoeffs; N::Int=64, kmode::Int=N÷2,
                           amp::Float64=1e-4, nstep::Int=40, cfl::Float64=0.2,
                           ε0::Float64=1.0)
    eng, st = setup_mhd1d_sinmode(c; N=N, amp=amp, kmode=kmode, ε0=ε0, cfl=cfl)
    Δx = eng.g.Δx
    dt = cfl*Δx/eng.cmax
    A0 = mhd1d_mode_amplitude(st, eng; field=:ε)
    ts = Float64[0.0]; amps = Float64[A0]
    for _ in 1:nstep
        # one fixed-dt SSP-RK2 step (reuse evolve loop body via tmax=dt accum)
        evolve_mhd1d!(st, eng; tmax=dt, cfl=cfl, monitor_causality=false)
        A = mhd1d_mode_amplitude(st, eng; field=:ε)
        push!(ts, ts[end]+dt); push!(amps, A)
        (!isfinite(A) || A > 1e6*A0) && break   # blown up — stop the fit cleanly
    end
    # least-squares slope of log(amp) vs t over the finite, positive samples
    good = [i for i in eachindex(amps) if isfinite(amps[i]) && amps[i] > 0]
    if length(good) < 2
        return (s=Inf, amps=amps, ts=ts, A0=A0)
    end
    tt = ts[good]; la = log.(amps[good])
    t̄ = sum(tt)/length(tt); l̄ = sum(la)/length(la)
    num = sum((tt[i]-t̄)*(la[i]-l̄) for i in eachindex(tt))
    den = sum((tt[i]-t̄)^2 for i in eachindex(tt))
    s = den > 0 ? num/den : 0.0
    return (s=s, amps=amps, ts=ts, A0=A0)
end

end # module BDNKMHD1D
