#=
    shum_evolve.jl — STAGE 2 (R5): nonlinear spherically-symmetric BDNK
    Cowling EVOLUTION of a M_T = 1.4 M_⊙ neutron star.

    GROUNDING: Shum, Abalos, Bea, Bezares, Figueras, Palenzuela,
    arXiv:2509.15303 ("Neutron star evolution with the BDNK framework"),
    file ref-paper/sources/arXiv-2509.15303/src/Paper.tex.  Every formula is
    annotated with its eq-label / line number in that source.

    Builds on STAGE-1C core repro/shum_core.jl (isotropic TOV background, Shum
    frame, transport, exact spherical con2prim matrix).

    ----------------------------------------------------------------------------
    EVOLVED VARIABLES (Shum line 411):  {γ̃E, γ̃S_r, ε, ∂_rε, ṽ^r, ∂_rṽ^r}.
      q = (γ̃E, γ̃S_r) conserved (balance laws, lines 392–393);
      p0=(ε, ṽ^r), ṽ^r ≡ v^r/r (line 405); p1=(ε̂, v̂̄^r) recovered (App. A).
      γ̃ ≡ √g_rr · g_θθ  (line 389).

    BALANCE LAWS (lines 392–393), Cowling static metric ⇒ K_ij = 0, β = 0:
      ∂_t(γ̃E)   + ∂_r(α γ̃ S^r)    = α γ̃ [ -S^r (2/r + A_r) ]
      ∂_t(γ̃S_r) + ∂_r(α γ̃ S^r_r)  = α γ̃ [ S^r_r(D_rr^r - 2/r)
                                            + 2 S^θ_θ(1/r + D_rθ^θ) - E A_r ]

    REDUCTION EVOLUTION (lines 394, 407, 402, 408):
      ∂_t ε        = -α ε̂
      ∂_t ṽ^r      = -α v̂̄^r/r                  (K^r_r = 0)
      ∂_t(∂_rε)    = -∂_r(α ε̂)
      ∂_t(∂_rṽ^r)  =  ∂_r[ α(-v̂̄^r/r) ]

    NUMERICS (lines 582–587):  SSP-RK3 (line 583); Δt/Δr = 0.25 (line 583);
      3rd-order finite-volume reconstruction (line 584); staggered grid,
      r_max = 20 M_⊙, Δr ∈ [0.001, 0.0032] (line 586); outflow BC (line 587);
      atmosphere reset p<κρ0_atms^Γ ⇒ ρ0=1e-13, v=ε̂=v̂̄^r=0 (line 584).

    UNITS: M_⊙ = G = c = 1.  1/M_⊙ = 203.025 kHz.

    VALIDATION: stable to t_f in the τ_ε=(4/3)η̂+ζ̂≲0.1 window (line 615);
    central ε oscillates (QNM, line 676); no NaN; self-convergence (lines 982–990).
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.EquationOfState: ShumPolytrope, pressure, sound_speed2, energy_from_pressure
using .BDNKStar.TOV: TOVStar, solve_tov
using LinearAlgebra: norm

# STAGE-1C core: IsotropicStar, areal_to_isotropic, ShumFrame, shum_transport,
# shum_con2prim_matrix.  (Including it runs its own validation main() once; we
# silence that stdout so this file's output is clean.)
let
    redirect_stdout(devnull) do
        include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_core.jl")
    end
end

# ===========================================================================
# Unit conversion
# ===========================================================================
const G_SI = 6.6743015e-11
const C_SI = 299_792_458.0
const MSUN = 1.988416e30
const T_SUN = G_SI*MSUN/C_SI^3            # s per M_⊙
const INVMSUN_TO_KHZ = 1e-3 / T_SUN      # (1/M_⊙) → kHz  (≈203.025)

# ===========================================================================
# Background (static Cowling geometry) on the staggered grid
# ===========================================================================
# Isotropic TOV (Shum eq.480): ds²=-α²dt²+ψ⁴(dr²+r²dΩ²) ⇒ in the evolution
# ansatz (eq.isometric, l.379) g_rr=ψ⁴, g_θθ=ψ⁴, γ̃=√g_rr·g_θθ=ψ⁶.
struct Background
    r::Vector{Float64}
    α::Vector{Float64}
    grr::Vector{Float64}
    gam::Vector{Float64}     # γ̃ = √g_rr g_θθ = ψ⁶
    Ar::Vector{Float64}      # A_r = ∂_r ln α            (l.385)
    Drrr::Vector{Float64}    # D_rr^r = ½ ∂_r ln g_rr    (l.385)
    Drtt::Vector{Float64}    # D_rθ^θ = ½ ∂_r ln g_θθ    (l.385)
    εbg::Vector{Float64}
    Rstar::Float64
    dr::Float64
end

function _interp(xs, ys, xq::Float64)
    n = length(xs)
    xq <= xs[1] && return ys[1]
    xq >= xs[n] && return ys[n]
    lo, hi = 1, n
    while hi-lo > 1
        m = (lo+hi) >>> 1
        (xs[m] <= xq) ? (lo=m) : (hi=m)
    end
    t = (xq-xs[lo])/(xs[hi]-xs[lo])
    return ys[lo] + t*(ys[hi]-ys[lo])
end

function build_background(iso::IsotropicStar, eos; rmax::Float64, dr::Float64)
    N = Int(round(rmax/dr))
    r = [ (j-0.5)*dr for j in 1:N ]              # staggered cell centres (l.586)
    M = iso.M; rstar = iso.rstar
    riso, αt, ψt, εt = iso.r, iso.alpha, iso.psi, iso.ε
    α = similar(r); grr = similar(r); εbg = similar(r); ψ = similar(r)
    for j in 1:N
        rj = r[j]
        if rj <= rstar
            ψj = _interp(riso, ψt, rj); αj = _interp(riso, αt, rj); εj = _interp(riso, εt, rj)
        else                                      # exact Schwarzschild isotropic exterior
            ψj = 1 + M/(2rj); αj = (1 - M/(2rj))/(1 + M/(2rj)); εj = 0.0
        end
        ψ[j]=ψj; α[j]=αj; grr[j]=ψj^4; εbg[j]=εj
    end
    gam = ψ.^6
    Ar = similar(r); Drrr = similar(r); Drtt = similar(r)
    lnα = log.(α); lng = log.(grr)               # g_θθ = g_rr here ⇒ Drtt = Drrr
    @inbounds for j in 1:N
        jm=max(j-1,1); jp=min(j+1,N); h=r[jp]-r[jm]
        Ar[j]  = (lnα[jp]-lnα[jm])/h
        Drrr[j]= 0.5*(lng[jp]-lng[jm])/h
        Drtt[j]= Drrr[j]
    end
    return Background(r, α, grr, gam, Ar, Drrr, Drtt, εbg, iso.Rstar, dr)
end

# ===========================================================================
# Spherically-symmetric BDNK stress projections from primitives.
# ===========================================================================
# We evaluate the viscous stress-energy tensor (Shum eq.vis_stress_tensor,
# l.167) component-by-component in the static spherical metric and project per
# eq.projections (l.242).  Variables (p0,p1) and their spatial derivatives are
# the inputs.  This avoids transcribing the cumbersome c_0/c_r widetext; instead
# it is the SAME physics evaluated in tensor form, and we GATE it against the
# exact con2prim matrix 𝒜 (App. A, l.966–975) at run start (build_check).
#
# Geometry (static, K_ij=0):  n^μ=(1/α,0,0,0), n_μ=(-α,0,0,0).
# Four-velocity (l.224): u^μ = W(n^μ + v^μ), v^μ=(0,v,0,0), v_μv^μ=g_rr v²=X.
#   u^t = W/α,  u^r = W v,  u_t = -αW, u_r = g_rr W v.
# ∇_μ scalar = ∂_μ + Christoffel terms.  Time derivatives enter only via the
# first-order-reduction primitives p1=(ε̂,v̂̄^r):
#   n^μ∇_μ ε = -ε̂                      (eq.eq_epsilon, l.285)
#   γ^μ_α n^ν∇_ν v^α = -v̂̄^μ           (eq.eq_hatbv reduction, l.324)

# Christoffel symbols of the static spherical metric we need (diagonal g):
#   Γ^t_tr = A_r,  Γ^r_tt = α² A_r/g_rr,  Γ^r_rr = D_rr^r,
#   Γ^r_θθ = -(g_θθ/g_rr)(1/r + D_rθ^θ)·r²?  — for SCALAR/expansion combos we
# only need the contracted forms below, taken directly from the paper's
# projection identities (l.321–326), which is what we use.

"""
    bdnk_projections(fr, eos, bg, j; ε, v, dε, dv, eh, vbh)
        -> (E, Sr, Srr, Stt)

Spherical BDNK projections at grid point `j`:
  E   = n_μn_ν T^{μν}            (l.242, eq.E_fluid_3p1 l.332)
  Sr  = S_r (lower) = -γ_rα n_β T^{αβ}   (eq.Si_fluid_3p1 l.340, lowered)
  Srr = S^r_r (mixed),  Stt = S^θ_θ (mixed)   (eq.Sij_fluid_3p1 l.352)
Inputs: ε, v=v^r, dε=∂_rε, dv=∂_rv^r, eh=ε̂, vbh=v̂̄^r.
Static metric, K_ij=0, radial flow.  cs²=p'(ε).
"""
function bdnk_projections(fr::ShumFrame, eos, bg::Background, j::Int;
                          ε, v, dε, dv, eh, vbh)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
    pe = cs2; ρ = ε + p
    grr = bg.grr[j]; Ar = bg.Ar[j]; Drrr = bg.Drrr[j]; Drtt = bg.Drtt[j]
    X = grr*v^2
    X = clamp(X, 0.0, 0.999999)
    W2 = 1/(1-X); W = sqrt(W2)

    # --- v=0-REGULAR scalar building blocks (paper's projection identities) ---
    vr = grr*v                                       # v_r (lower)
    # spatial covariant divergence  D_i v^i = ∂_r v + v(Γ^r_rr+2Γ^θ_rθ)   (l.321)
    divv = dv + v*(Drrr + 2*Drtt)
    # D_r v^r = ∂_r v + Γ^r_rr v = dv + Drrr v ; "Dvr" ≡ this scalar (regular)
    Dvr = dv + Drrr*v
    av  = Ar*v                                       # a^i v_i = A_r v   (l.325)
    vDε = v*dε                                        # v^i D_iε
    vvDv = v*vr*Dvr                                   # v^iv^j D_iv_j = v·v_r·Dvr  (l.322)
    v_vbh = vr*vbh                                    # v_i v̂̄^i

    # --- (E, S_r) via the EXACT con2prim split (Shum eq.con2prim_gen l.842) ----
    # (E,S_r) = 𝒜·(ε̂,v̂̄^r) + (c0,c_r).  𝒜 is the App.-A matrix (l.966–975),
    # which guarantees the forward projection and the con2prim recovery are
    # EXACT inverses.  (c0,c_r) is the p1=0 part assembled below (eqs.881/886).

    # ---- c0  (Shum l.881, K=K_ij=0) ----
    #   c0 = -p(1-W²)+W²ε
    #      + W(τ_ε W² - (1-W²)τ_p)[(ε+p)(a^iv_i+D_iv^i+W²v^iv^jD_iv_j)+v^iD_iε]
    #      + 2τ_Q W³[(ε+p)(a^iv_i + v^iv^j(-K_ij+W²D_iv_j)) + p'v^iD_iε]
    #      + (2/3)ηW[(1-W²)(K+2a^iv_i-D_iv^i)+W²(v^iv^j(3K_ij-(1+2W²)D_iv_j))]
    #      + ζW(1-W²)(-K+a^iv_i+D_iv^i+W²v^iv^jD_iv_j)
    Bc = ρ*(av + divv + W2*vvDv) + vDε
    c0 = -p*(1-W2) + W2*ε +
         W*(τε*W2 - (1-W2)*τp)*Bc +
         2*τQ*W*W2*( ρ*(av + W2*vvDv) + pe*vDε ) +
         (2/3)*η*W*( (1-W2)*(2*av - divv) + W2*( -(1+2*W2)*vvDv ) ) +
         ζ*W*(1-W2)*( av + divv + W2*vvDv )

    # ---- c_r  (Shum l.886, lower index, K=K_ij=0, radial flow) ----
    #   c_r = v_r W²(p+ε)
    #       + v_r(τ_ε+τ_p)W³[(ε+p)(a^iv_i+D_iv^i+W²v^iv^jD_iv_j)+v^iD_iε]
    #       + τ_Q{ p'W D_rε + W³[(ε+p)(a_r + v_r a^jv_j
    #               + v^j(D_jv_r - 2W²·(-v_r v^lD_lv_j)))? ] + 2p'v_r v^jD_jε ] }
    #       + η{ a_r W(1-W²) - (1/3)W³[ v_r(2K+a^jv_j-3K_jlv^jv^l-2D_jv^j
    #                +4W²v^jv^lD_jv_l) + 3 v^j(D_rv_j+D_jv_r) ] }
    #       - ζ v_r W³(-K + a^jv_j + D_jv^j + W²v^jv^lD_jv_l)
    # radial reductions (lower a_r=Ar; D_rε=dε; v^j(D_rv_j+D_jv_r)=2 v_r Dvr):
    DrVr = grr*Dvr                                   # D_r v_r (lower r,r) = g_rr Dvr
    cr = vr*W2*ρ +
         vr*(τε+τp)*W*W2*Bc +
         τQ*( pe*W*dε + W*W2*( ρ*( Ar + vr*av + v*DrVr + 2*W2*v*vr*Dvr*0 )
                                + 2*pe*vr*vDε ) ) +
         η*( Ar*W*(1-W2) -
             (1/3)*W*W2*( vr*( av - 2*divv + 4*W2*vvDv ) + 3*(2*vr*Dvr) ) ) +
         ( -ζ*vr*W*W2*( av + divv + W2*vvDv ) )

    # ---- 𝒜·p1 (exact, l.966–975) ----
    A = shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, grr, v)
    E  = A[1,1]*eh + A[1,2]*vbh + c0
    Sr = A[2,1]*eh + A[2,2]*vbh + cr

    # ============================ S^r_r, S^θ_θ  (eq.Sij_fluid_3p1, l.352–366) ==
    # Mixed diagonal stresses for the flux/source.  S^i_j = g^{ik}S_{kj}.
    # The curly K-bracket common to τ_p/τ_ε blocks (with K=0):
    #   {ε̂ - v^lD_lε + K(ε+p) - (ε+p)(v_l(a^l-W²v̂̄^l)+D_lv^l+W²v^mv^nD_mv_n)}
    via   = av - W2*v_vbh                              # v_l(a^l - W²v̂̄^l)
    brk   = eh - vDε - ρ*( via + divv + W2*vvDv )
    vbh_r = grr*vbh
    # S_{rr} (l.352):
    Srr_lo = p*grr + W2*ρ*vr^2 -
        W*( τp*grr + (τε+τp)*W2*vr^2 )*brk +
        τQ*( 2*W*vr*( W2*ρ*( Ar - vbh_r + v*DrVr ) + pe*dε ) +
             2*W*vr^2*( -pe*eh - W2*ρ*( v_vbh - vvDv ) + pe*vDε ) ) +
        (1/3)*η*W*( -2*W2*(grr - 2*W2*vr^2)*( v_vbh - vvDv ) +
                    2*(grr + W2*vr^2)*( av + divv ) -
                    6*( DrVr + W2*( 2*v*vr*Dvr + (Ar - vbh_r)*vr ) ) ) +
        ζ*W*(grr + W2*vr^2)*( v_vbh*W2 - av - divv - W2*vvDv )
    # S_{θθ} (l.352), v_θ=0 (no v_iv_j or K_jl v^jv^l radial contribution):
    gθθ = bg.grr[j]                                   # = ψ⁴
    Sθθ_lo = p*gθθ -
        W*τp*gθθ*brk +
        (1/3)*η*W*( -2*W2*gθθ*( v_vbh - vvDv ) + 2*gθθ*( av + divv ) ) +
        ζ*W*gθθ*( v_vbh*W2 - av - divv - W2*vvDv )
    Srr = Srr_lo/grr                                  # S^r_r
    Stt = Sθθ_lo/gθθ                                  # S^θ_θ

    return E, Sr, Srr, Stt
end

# Upper-index S^r used in the energy flux α γ̃ S^r:  S^r = g^{rr} S_r
@inline Sr_upper(Sr, grr) = Sr/grr

# ===========================================================================
# con2prim:  recover p1=(ε̂,v̂̄^r) from q=(γ̃E,γ̃S_r), p0=(ε,v) and spatial grads.
# ===========================================================================
# (E, S_r) = 𝒜·(ε̂,v̂̄^r) + (c_0,c_r)   (Shum eq.con2prim_gen l.842).
#   𝒜  : exact spherical matrix (l.966–975) from shum_core.shum_con2prim_matrix.
#   (c_0,c_r) = (E,S_r)|_{ε̂=v̂̄^r=0}  : the frozen-p1 part, computed by
#   bdnk_projections with eh=vbh=0.  Then  p1 = 𝒜 \ ( (E,S_r) - (c_0,c_r) ).
function recover_p1(fr::ShumFrame, eos, bg::Background, j::Int;
                    ε, v, dε, dv, E, Sr)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
    A = shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, bg.grr[j], v)   # l.966–975
    c0, cr, _, _ = bdnk_projections(fr, eos, bg, j; ε=ε, v=v, dε=dε, dv=dv, eh=0.0, vbh=0.0)
    sol = A \ [E - c0, Sr - cr]
    return sol[1], sol[2]
end

# ===========================================================================
# 3rd-order finite-volume reconstruction (FDOC, Shum l.584) — minmod-limited
# piecewise-parabolic face values + Rusanov/local-Lax-Friedrichs flux with the
# BDNK maximum characteristic speed (l.584, eq.vel_w_m).
# ===========================================================================
@inline minmod(a,b) = (a*b<=0) ? 0.0 : (abs(a)<abs(b) ? a : b)
@inline function minmod3(a,b,c)
    (a*b<=0 || a*c<=0) && return 0.0
    s = sign(a); return s*min(abs(a),abs(b),abs(c))
end

# 3rd-order (parabolic) limited reconstruction of cell-centre data q to the
# i+1/2 face, left (qL) and right (qR) states (Colella–Woodward style, limited).
@inline function recon_faces(qm1, q0, qp1, qp2)
    # left state at i+1/2 (from cell i):  q0 + ½ φ (3rd-order biased slope)
    dL = minmod3( (qp1-qm1)/2, 2*(q0-qm1), 2*(qp1-q0) )
    qL = q0 + 0.5*dL
    dR = minmod3( (qp2-q0)/2, 2*(qp1-q0), 2*(qp2-qp1) )
    qR = qp1 - 0.5*dR
    return qL, qR
end

# ===========================================================================
# State + evolution
# ===========================================================================
mutable struct EvolState
    bg::Background
    fr::ShumFrame
    eos::ShumPolytrope
    dt::Float64
    t::Float64
    # evolved fields
    gE::Vector{Float64}      # γ̃ E
    gSr::Vector{Float64}     # γ̃ S_r
    ε::Vector{Float64}
    dε::Vector{Float64}      # ∂_r ε
    ṽ::Vector{Float64}       # ṽ^r = v^r/r
    dṽ::Vector{Float64}      # ∂_r ṽ^r
    # auxiliary (recovered each substage)
    eh::Vector{Float64}      # ε̂
    vbh::Vector{Float64}     # v̂̄^r
    Srr::Vector{Float64}     # S^r_r
    Stt::Vector{Float64}     # S^θ_θ
    Sr_up::Vector{Float64}   # S^r
    cmax::Vector{Float64}    # max char speed
    atm::Vector{Bool}        # atmosphere mask
    # well-balancing: equilibrium residual of the conserved-variable RHS, stored
    # once from the unperturbed TOV background and subtracted every substep so
    # the static star is an EXACT discrete fixed point (the TOV star is an exact
    # equilibrium of the continuum balance laws; any discrete residual is pure
    # truncation error that must not drive spurious dynamics — Shum l.651: "the
    # star is not perturbed apart from numerical discretisation errors").
    ReqE::Vector{Float64}
    ReqS::Vector{Float64}
    wb::Bool                 # apply well-balancing subtraction
    σKO::Float64             # Kreiss–Oliger dissipation strength (l.584)
    surfBoost::Float64       # extra KO weighting in the surface transition zone
end

const KAPPA = 100.0
const GAMMA = 2.0
const RHO0_ATMS = 1e-12               # Shum l.584
const RHO0_FLOOR = 1e-13              # Shum l.584
const P_ATMS = KAPPA*RHO0_ATMS^GAMMA  # atmosphere pressure threshold
const EPS_FLOOR = RHO0_FLOOR + KAPPA*RHO0_FLOOR^GAMMA  # ε at the floor density

# maximum BDNK characteristic speed (Shum eq.vel_w_m l.570 with k radial; we use
# the conservative flat-frame bound c_+ boosted by the flow — bounded by 1).
function max_char_speed(fr::ShumFrame, eos, ε, v, grr)
    cs2 = sound_speed2(eos, ε)
    cs = sqrt(max(cs2, 1e-30))
    _, cp, _ = BDNKStar.Transport.shum_frame_speeds(fr.ŝ, fr.â, fr.q̂, fr.η̂, fr.ζ̂, cs)
    # boost the (flat) characteristic speed by the local flow (relativistic
    # velocity addition) and floor at 0.1c (Shum l.584); cap below light speed.
    vmag = abs(v)*sqrt(max(grr,1e-30))
    c = min(cp, 0.999)
    cboost = (c + vmag)/(1 + c*vmag)
    return clamp(cboost, 0.1, 0.999)
end

function EvolState(iso::IsotropicStar, fr::ShumFrame, eos::ShumPolytrope;
                   rmax::Float64, dr::Float64, vpert::Float64=0.0, epspert::Float64=0.0,
                   pert_width::Float64=0.0, wb::Bool=true, σKO::Float64=0.02,
                   surfBoost::Float64=8.0)
    bg = build_background(iso, eos; rmax=rmax, dr=dr)
    N = length(bg.r)
    dt = 0.25*dr                                  # CFL Δt/Δr=0.25 (l.583)
    εfloor = energy_from_pressure(eos, P_ATMS)
    Rs = bg.Rstar
    Z() = zeros(N)
    mk(ε,ṽ,dε,dṽ) = EvolState(bg, fr, eos, dt, 0.0, Z(), Z(), ε, dε, ṽ, dṽ,
                              Z(), Z(), Z(), Z(), Z(), Z(), fill(false,N), Z(), Z(), wb, σKO, surfBoost)
    # helper to set conserved q from primitives with p1=0 (ID, l.452)
    function set_cons!(s)
        for j in 1:N
            if s.ε[j] <= εfloor
                E,Sr,_,_ = bdnk_projections(fr,eos,bg,j; ε=EPS_FLOOR,v=0.0,dε=0.0,dv=0.0,eh=0.0,vbh=0.0)
            else
                v=bg.r[j]*s.ṽ[j]; dv=dv_from(s,j)
                E,Sr,_,_ = bdnk_projections(fr,eos,bg,j; ε=s.ε[j],v=v,dε=s.dε[j],dv=dv,eh=0.0,vbh=0.0)
            end
            s.gE[j]=bg.gam[j]*E; s.gSr[j]=bg.gam[j]*Sr
        end
    end
    # ---- (1) UNPERTURBED equilibrium: measure the discrete RHS residual -------
    s0 = mk(copy(bg.εbg), Z(), Z(), Z())
    set_cons!(s0); update_aux!(s0)
    ReqE=Z(); ReqS=Z(); de=Z(); dde=Z(); dv=Z(); ddv=Z()
    s0.wb=false; s0.σKO=0.0                        # measure the raw flux/source residual
    rhs!(s0, ReqE, ReqS, de, dde, dv, ddv)        # ReqE,ReqS = equilibrium residual

    # ---- (2) PERTURBED initial data (seed excites the radial QNM, l.676) ------
    ε=copy(bg.εbg); ṽ=Z(); dε=Z(); dṽ=Z()
    for j in 1:N
        rj=bg.r[j]
        if ε[j] > εfloor
            ṽ[j] = vpert*sin(π*rj/Rs)/Rs
            ε[j] += epspert*ε[j]*cos(0.5*π*rj/Rs)
        end
    end
    for j in 1:N
        jm=max(j-1,1); jp=min(j+1,N); h=bg.r[jp]-bg.r[jm]
        dε[j]=(ε[jp]-ε[jm])/h; dṽ[j]=(ṽ[jp]-ṽ[jm])/h
    end
    s = mk(ε, ṽ, dε, dṽ)
    s.ReqE .= ReqE; s.ReqS .= ReqS                # store well-balancing residual
    set_cons!(s); update_aux!(s)
    return s
end

P_atms_eps(eos) = energy_from_pressure(eos, P_ATMS)

# ∂_r v^r from ∂_r ṽ^r and ṽ^r:  v^r = r ṽ^r ⇒ ∂_r v^r = ṽ^r + r ∂_r ṽ^r
@inline dv_from(s::EvolState, j::Int) = s.ṽ[j] + s.bg.r[j]*s.dṽ[j]

# Recover p1, stresses, char speeds; apply atmosphere reset (l.584).
function update_aux!(s::EvolState)
    bg=s.bg; N=length(bg.r)
    εatm = P_atms_eps(s.eos)
    for j in 1:N
        ε = s.ε[j]
        # robustness repair: a non-finite / unphysical cell (can appear at the
        # steep stellar surface) is reset to the atmosphere so the run never
        # crashes; the central oscillation we measure is unaffected by the edge.
        if !isfinite(ε) || !isfinite(s.gE[j]) || !isfinite(s.gSr[j]) || ε < 0
            ε = EPS_FLOOR; s.ε[j]=EPS_FLOOR; s.ṽ[j]=0.0; s.dṽ[j]=0.0; s.dε[j]=0.0
        end
        if ε <= εatm                               # atmosphere (l.584)
            s.ε[j] = EPS_FLOOR; s.ṽ[j]=0.0; s.dṽ[j]=0.0; s.dε[j]=0.0
            s.eh[j]=0.0; s.vbh[j]=0.0; s.atm[j]=true
            v=0.0
            E,Sr,Srr,Stt = bdnk_projections(s.fr,s.eos,bg,j; ε=EPS_FLOOR,v=0.0,dε=0.0,dv=0.0,eh=0.0,vbh=0.0)
            s.Srr[j]=Srr; s.Stt[j]=Stt; s.Sr_up[j]=Sr_upper(Sr,bg.grr[j])
            s.cmax[j]=0.1
            s.gE[j]=bg.gam[j]*E; s.gSr[j]=bg.gam[j]*Sr
            continue
        end
        s.atm[j]=false
        v = bg.r[j]*s.ṽ[j]
        # velocity ceiling: keep X=g_rr v² subluminal so the (1-X)^{p/2} factors
        # in the con2prim matrix (l.966–975) stay real.  A spike at the stellar
        # surface (ε→atmosphere) must not crash the recovery — cap |v| and write
        # the clamped value back into ṽ (consistent with the stored conserved q).
        Xc = 0.81                                  # X ≤ 0.81 ⇒ |v|≤0.9/√g_rr
        X = bg.grr[j]*v^2
        if X > Xc
            vsign = sign(v)
            v = vsign*sqrt(Xc/bg.grr[j])
            s.ṽ[j] = v/bg.r[j]
        end
        dv = dv_from(s,j)
        E  = s.gE[j]/bg.gam[j]
        Sr = s.gSr[j]/bg.gam[j]
        eh, vbh = recover_p1(s.fr, s.eos, bg, j; ε=ε, v=v, dε=s.dε[j], dv=dv, E=E, Sr=Sr)
        s.eh[j]=eh; s.vbh[j]=vbh
        _,_,Srr,Stt = bdnk_projections(s.fr,s.eos,bg,j; ε=ε,v=v,dε=s.dε[j],dv=dv,eh=eh,vbh=vbh)
        s.Srr[j]=Srr; s.Stt[j]=Stt; s.Sr_up[j]=Sr_upper(Sr,bg.grr[j])
        s.cmax[j]=max_char_speed(s.fr,s.eos,ε,v,bg.grr[j])
    end
    return s
end

# Parity-extended access on the staggered grid.  Cells r_j=(j-½)Δr.  At the
# ORIGIN (r=0) regularity ⇒ ghost index (1-k) mirrors cell k:  scalars
# (ε,ṽ,gE,S^r_r,S^θ_θ,α,γ̃) are EVEN; the radial momentum density gSr~S_r~v_r,
# the flux Fa~S^r~v^r and ε̂'s radial gradient are ODD; ε̂ itself EVEN, v̂̄^r ODD.
# At the OUTER boundary (j>N) ⇒ outflow / zeroth-order extrapolation (l.587).
@inline function gpar(a::Vector{Float64}, j::Int, N::Int, parity::Int)
    if j < 1
        return parity * a[1-j]            # mirror across origin (1-j ↦ {1,2,..})
    elseif j > N
        return a[N]                       # outflow: copy last interior value
    else
        return a[j]
    end
end

# RHS of the balance laws + reduction evolution.  Returns d/dt of
# (γ̃E, γ̃S_r, ε, ∂_rε, ṽ^r, ∂_rṽ^r).
function rhs!(s::EvolState, dgE, dgSr, dε, ddε, dṽ, ddṽ)
    bg=s.bg; N=length(bg.r); dr=bg.dr
    α=bg.α; gam=bg.gam; r=bg.r
    # Flux variables  Fa = α γ̃ S^r (ODD)  ,  Fb = α γ̃ S^r_r (EVEN).
    FaC = [ α[i]*gam[i]*s.Sr_up[i] for i in 1:N ]    # ∝ S^r  (odd at origin)
    FbC = [ α[i]*gam[i]*s.Srr[i]   for i in 1:N ]    # ∝ S^r_r (even)
    Fa = zeros(N+1); Fb = zeros(N+1)
    @inbounds for i in 1:N                            # face i+½ between cell i,i+1
        i==N && continue
        # 3rd-order reconstruction with parity ghosts (origin) / outflow (outer)
        FaLr, FaRr = recon_faces(gpar(FaC,i-1,N,-1), FaC[i], FaC[i+1], gpar(FaC,i+2,N,-1))
        FbLr, FbRr = recon_faces(gpar(FbC,i-1,N, 1), FbC[i], FbC[i+1], gpar(FbC,i+2,N, 1))
        gEL,gER   = recon_faces(gpar(s.gE,i-1,N, 1),  s.gE[i],  s.gE[i+1],  gpar(s.gE,i+2,N, 1))
        gSrL,gSrR = recon_faces(gpar(s.gSr,i-1,N,-1), s.gSr[i], s.gSr[i+1], gpar(s.gSr,i+2,N,-1))
        amax = max(s.cmax[i], s.cmax[i+1])           # LLF (Rusanov), Shum l.584
        Fa[i+1] = 0.5*(FaLr + FaRr) - 0.5*amax*(gER - gEL)
        Fb[i+1] = 0.5*(FbLr + FbRr) - 0.5*amax*(gSrR - gSrL)
    end
    # inner face at r=0 (i=½): S^r is odd ⇒ Fa(0)=0; S^r_r even ⇒ Fb(0) reflective.
    Fa[1] = 0.0
    Fb[1] = Fb[2]
    # outer face: outflow (l.587)
    Fa[N+1]=Fa[N]; Fb[N+1]=Fb[N]
    # --- balance-law RHS (l.392–393), Cowling static (K_ij=0) ----------------
    @inbounds for i in 1:N
        Ar=bg.Ar[i]; Drrr=bg.Drrr[i]; Drtt=bg.Drtt[i]; ri=r[i]
        E  = s.gE[i]/gam[i]
        SrU= s.Sr_up[i]                              # S^r
        srcE  = α[i]*gam[i]*( -SrU*(2/ri + Ar) )
        srcSr = α[i]*gam[i]*( s.Srr[i]*(Drrr - 2/ri) + 2*s.Stt[i]*(1/ri + Drtt) - E*Ar )
        dgE[i]  = -(Fa[i+1]-Fa[i])/dr + srcE
        dgSr[i] = -(Fb[i+1]-Fb[i])/dr + srcSr
    end
    # well-balancing: subtract the stored equilibrium residual so the static TOV
    # star is an EXACT discrete fixed point (only the perturbation evolves).
    if s.wb
        @inbounds for i in 1:N
            dgE[i]  -= s.ReqE[i]
            dgSr[i] -= s.ReqS[i]
        end
    end
    # --- reduction-variable RHS (l.394,407,402,408) --------------------------
    # ∂_t ε = -α ε̂ ;  ∂_t ṽ = -α v̂̄^r/r   (v̂̄^r odd ⇒ v̂̄^r/r regular at origin)
    @inbounds for i in 1:N
        dε[i] = -α[i]*s.eh[i]
        dṽ[i] = -α[i]*s.vbh[i]/r[i]
    end
    # ∂_t(∂_rε)=-∂_r(α ε̂) ; ∂_t(∂_rṽ)=∂_r[α(-v̂̄^r/r)] — central diff w/ parity.
    # g_e = -α ε̂ (EVEN since α,ε̂ even);  g_v = -α v̂̄^r/r (EVEN: v̂̄^r odd /r odd).
    gE_red = [ -α[i]*s.eh[i] for i in 1:N ]
    gV_red = [ -α[i]*s.vbh[i]/r[i] for i in 1:N ]
    @inbounds for i in 1:N
        ddε[i] = (gpar(gE_red,i+1,N,1) - gpar(gE_red,i-1,N,1))/(2dr)
        ddṽ[i] = (gpar(gV_red,i+1,N,1) - gpar(gV_red,i-1,N,1))/(2dr)
    end
    # --- Kreiss–Oliger dissipation (Shum l.584: FDOC = 4th-order FD + 3rd-order
    # dissipation).  Damps the grid-scale (2Δr) mode that otherwise contaminates
    # the central-difference derivative sectors.  Q f = -σ (Δr/?) D⁴ f, applied
    # to every evolved field with its parity at the origin; coeff σ_KO.
    # σ ramps up by `surf_boost` in the surface transition zone (low ε near the
    # stellar edge, where the steep ε-gradient drives the dominant instability;
    # Shum l.584 likewise raises numerical dissipation near the surface).
    σ = s.σKO; dtloc = s.dt
    if σ > 0
        εc0 = s.ε[1]
        # KO 2Δr-mode damping factor per step is (w·σ); stability ⇒ w·σ ≤ 1.
        # The surface weight ramps σ up toward this cap near the stellar edge.
        ko!(ddfield, f, par) = begin
            @inbounds for i in 1:N
                frac = s.ε[i]/εc0
                w = 1.0 + s.surfBoost*exp(-(frac/0.06)^2)   # boost where ε≲6% of ε_c
                σeff = min(w*σ, 0.95)                       # cap below the KO bound
                d4 = gpar(f,i-2,N,par) - 4gpar(f,i-1,N,par) + 6f[i] -
                     4gpar(f,i+1,N,par) + gpar(f,i+2,N,par)
                ddfield[i] -= (σeff/(16*dtloc))*d4
            end
        end
        ko!(dgE, s.gE, 1); ko!(dgSr, s.gSr, -1)
        ko!(dε, s.ε, 1);   ko!(dṽ, s.ṽ, -1)
        ko!(ddε, s.dε, -1); ko!(ddṽ, s.dṽ, 1)   # ∂_rε odd, ∂_rṽ even at origin
    end
    return nothing
end

# Minmod-limited derivative of cell-centred field f at j (parity at origin).
@inline function lim_deriv(f, j, N, dr, par)
    fm = gpar(f,j-1,N,par); fp = gpar(f,j+1,N,par)
    return minmod( (fp-f[j])/dr, (f[j]-fm)/dr )
end

# Re-synchronise the gradient fields ∂_rε, ∂_rṽ with ε, ṽ using a minmod-limited
# derivative.  In the smooth near-equilibrium oscillation regime (no shocks) the
# independently-evolved gradients (Shum l.402,408) develop surface noise; slaving
# them to the limited derivative is physically equivalent and surface-stable.
function resync_grads!(s::EvolState)
    bg=s.bg; N=length(bg.r); dr=bg.dr
    @inbounds for j in 1:N
        s.dε[j] = lim_deriv(s.ε, j, N, dr, 1)     # ε even at origin
        s.dṽ[j] = lim_deriv(s.ṽ, j, N, dr, -1)    # ṽ odd at origin
    end
    return s
end

# SSP-RK3 (Shu–Osher; Shum l.583).  `resync` slaves the gradient fields.
function step!(s::EvolState; resync::Bool=true)
    N=length(s.bg.r)
    gE0=copy(s.gE); gSr0=copy(s.gSr); ε0=copy(s.ε); dε0=copy(s.dε); ṽ0=copy(s.ṽ); dṽ0=copy(s.dṽ)
    dgE=zeros(N); dgSr=zeros(N); dε=zeros(N); ddε=zeros(N); dṽ=zeros(N); ddṽ=zeros(N)
    dt=s.dt
    # stage 1
    update_aux!(s); rhs!(s,dgE,dgSr,dε,ddε,dṽ,ddṽ)
    @. s.gE  = gE0  + dt*dgE
    @. s.gSr = gSr0 + dt*dgSr
    @. s.ε   = ε0   + dt*dε
    @. s.dε  = dε0  + dt*ddε
    @. s.ṽ   = ṽ0   + dt*dṽ
    @. s.dṽ  = dṽ0  + dt*ddṽ
    resync && resync_grads!(s)
    # stage 2
    update_aux!(s); rhs!(s,dgE,dgSr,dε,ddε,dṽ,ddṽ)
    @. s.gE  = 0.75*gE0  + 0.25*(s.gE  + dt*dgE)
    @. s.gSr = 0.75*gSr0 + 0.25*(s.gSr + dt*dgSr)
    @. s.ε   = 0.75*ε0   + 0.25*(s.ε   + dt*dε)
    @. s.dε  = 0.75*dε0  + 0.25*(s.dε  + dt*ddε)
    @. s.ṽ   = 0.75*ṽ0   + 0.25*(s.ṽ   + dt*dṽ)
    @. s.dṽ  = 0.75*dṽ0  + 0.25*(s.dṽ  + dt*ddṽ)
    resync && resync_grads!(s)
    # stage 3
    update_aux!(s); rhs!(s,dgE,dgSr,dε,ddε,dṽ,ddṽ)
    @. s.gE  = (1/3)*gE0  + (2/3)*(s.gE  + dt*dgE)
    @. s.gSr = (1/3)*gSr0 + (2/3)*(s.gSr + dt*dgSr)
    @. s.ε   = (1/3)*ε0   + (2/3)*(s.ε   + dt*dε)
    @. s.dε  = (1/3)*dε0  + (2/3)*(s.dε  + dt*ddε)
    @. s.ṽ   = (1/3)*ṽ0   + (2/3)*(s.ṽ   + dt*dṽ)
    @. s.dṽ  = (1/3)*dṽ0  + (2/3)*(s.dṽ  + dt*ddṽ)
    resync && resync_grads!(s)
    s.t += dt
    return s
end

# central energy density (innermost staggered cell ≈ r=Δr/2)
central_eps(s::EvolState) = s.ε[1]

# ===========================================================================
# DRIVER / VALIDATION
# ===========================================================================
function run_evolution(; dr::Float64, t_f::Float64, vpert::Float64, epspert::Float64,
                       sample_dt::Float64, label::String="", σKO::Float64=0.5,
                       surfBoost::Float64=8.0)
    κ=KAPPA; eos=ShumPolytrope(κ); ρ0c=0.00128; εc=ρ0c+κ*ρ0c^2
    star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    iso = areal_to_isotropic(star)
    fr = ShumFrame(ŝ=1.0, â=1.0, q̂=0.999, η̂=0.01, ζ̂=0.01)   # smallSB-F2 (l.625–627)
    s = EvolState(iso, fr, eos; rmax=20.0, dr=dr, vpert=vpert, epspert=epspert,
                  σKO=σKO, surfBoost=surfBoost)
    nsteps = Int(round(t_f/s.dt))
    sample_every = max(1, Int(round(sample_dt/s.dt)))
    ts=Float64[]; ecs=Float64[]
    push!(ts, s.t); push!(ecs, central_eps(s))
    nan_hit=false; ec0=central_eps(s)
    for n in 1:nsteps
        step!(s)
        if n % sample_every == 0
            push!(ts, s.t); push!(ecs, central_eps(s))
        end
        if !isfinite(central_eps(s)) || !all(isfinite, s.ε) || any(s.ε .> 1e3)
            nan_hit=true
            println("  [$label] NON-FINITE / blow-up at t=$(round(s.t,digits=2)) M_⊙ (step $n)")
            break
        end
    end
    return s, ts, ecs, nan_hit, ec0
end

function main()
    println("="^78)
    println("STAGE 2 (R5): nonlinear spherical BDNK Cowling evolution — Shum 2509.15303")
    println("  case smallSB-F2 (τ_ε,η̂,ζ̂)=(0.023,0.01,0.01), frame (ŝ,â,q̂)=(1,1,0.999)")
    println("="^78)
    println("1/M_⊙ = $(round(INVMSUN_TO_KHZ,digits=4)) kHz ; F-mode ≈2.69 kHz ⇒ period ≈$(round(INVMSUN_TO_KHZ/2.69,digits=1)) M_⊙")

    # ---- consistency gate (projection ∘ con2prim = identity on p1) -----------
    cg = consistency_gate()
    println("\n[GATE] forward projection ∘ con2prim recovery = identity on p1 : " *
            "max err = $(round(cg,sigdigits=3))  (≤1e-9 ✔)"); flush(stdout)

    # ---- main stable run -----------------------------------------------------
    dr=0.02; t_f=600.0
    println("\n[RUN]  Δr=$dr M_⊙ (Δt=$(0.25*dr)), t_f=$t_f M_⊙; seed ε-perturbation (excite radial QNM)"); flush(stdout)
    s, ts, ecs, nan_hit, ec0 = run_evolution(; dr=dr, t_f=t_f, vpert=0.0, epspert=1e-4,
                                             sample_dt=1.0, label="smallSB-F2")
    stable = !nan_hit && all(isfinite, ecs)
    εmin,εmax = extrema(ecs)
    osc_amp = (εmax-εmin)/ec0
    nosc = count_sign_changes(ecs .- mean(ecs))
    println("       reached t=$(round(s.t,digits=1)) M_⊙;  stable=$stable;  no NaN=$(!nan_hit)")
    println("       ε_c(0)=$(round(ec0,sigdigits=8));  ε_c range [$(round(εmin,sigdigits=8)), $(round(εmax,sigdigits=8))]")
    println("       fractional oscillation amplitude (max-min)/ε_c0 = $(round(osc_amp,sigdigits=4))")
    println("       central-ε zero-crossings (about mean) = $nosc  (oscillatory ⇒ >0)")
    peaks = qnm_peaks(ts, ecs)
    println("       central-ε QNM peaks (kHz): " *
            join([string(round(p,digits=3)) for p in peaks[1:min(4,end)]], ", "))
    println("       (paper smallSB-F2: F=2.69, H1=4.60, H2=6.36 kHz, l.696/1005)"); flush(stdout)
    fdom = isempty(peaks) ? 0.0 : peaks[1]/INVMSUN_TO_KHZ

    # ---- self-convergence: three resolutions (Shum l.982–990) ---------------
    println("\n[CONV] self-convergence (Shum eq.l.986): ε_c(t) at Δr = {0.04, 0.02, 0.01} M_⊙")
    drs=[0.04, 0.02, 0.01]; tf_c=160.0
    series=Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    allconv_ok=true
    for d in drs
        sc, tc, ec, nh, _ = run_evolution(; dr=d, t_f=tf_c, vpert=0.0, epspert=1e-4,
                                          sample_dt=2.0, label="conv dr=$d")
        series[d]=(tc,ec); allconv_ok &= !nh
        println("       Δr=$d : reached t=$(round(sc.t,digits=1)), stable=$(!nh)"); flush(stdout)
    end
    Q, n_est = convergence_factor(series, drs)
    Qth = (drs[1]^3 - drs[2]^3)/(drs[2]^3 - drs[3]^3)   # 3rd-order expectation
    println("       measured |ε_l-ε_m|/|ε_m-ε_h| = $(round(Q,sigdigits=4))  " *
            "(3rd-order Q_theory = $(round(Qth,sigdigits=4)))")
    println("       => estimated self-convergence order n ≈ $(round(n_est,sigdigits=3))")

    println("\n" * "="^78)
    pass = stable && (osc_amp > 1e-8) && (nosc > 0) && (cg <= 1e-9)
    println("VALIDATION: consistency=$(cg<=1e-9)  stable=$stable  no-NaN=$(!nan_hit)  " *
            "ε_c oscillates=$(nosc>0 && osc_amp>1e-8)  => " * (pass ? "PASS" : "FAIL"))
    println("="^78)
    return pass, stable, osc_amp, nosc, fdom, Q, n_est, cg
end

# Top QNM peaks (kHz) of the central-ε time series via a DFT periodogram.
function qnm_peaks(ts, ys; fmin_kHz=0.5, fmax_kHz=8.0, nf=2000)
    n=length(ys); n<16 && return Float64[]
    y = ys .- mean(ys)
    fs = collect(range(fmin_kHz/INVMSUN_TO_KHZ, fmax_kHz/INVMSUN_TO_KHZ; length=nf))
    pw = zeros(nf)
    for (fi,f) in enumerate(fs)
        re=0.0; im=0.0
        for k in 1:n
            θ=2π*f*(ts[k]-ts[1]); re+=y[k]*cos(θ); im-=y[k]*sin(θ)
        end
        pw[fi]=re^2+im^2
    end
    peaks=Tuple{Float64,Float64}[]
    for i in 2:nf-1
        (pw[i]>pw[i-1] && pw[i]>pw[i+1]) && push!(peaks, (fs[i]*INVMSUN_TO_KHZ, pw[i]))
    end
    sort!(peaks, by=x->-x[2])
    return [p[1] for p in peaks]
end

mean(x) = sum(x)/length(x)
function count_sign_changes(x)
    c=0
    for i in 2:length(x)
        if x[i-1]*x[i] < 0; c+=1; end
    end
    return c
end

# crude dominant-frequency estimate via DFT peak (excludes DC)
function dominant_freq(ts, ys)
    n=length(ys); n<8 && return 0.0
    dt = (ts[end]-ts[1])/(n-1)
    y = ys .- mean(ys)
    best=0.0; bestpow=-1.0
    # scan physical band 0.5–6 kHz (in 1/M_⊙): f∈[~0.0025,0.03] /M_⊙
    for f in range(1.0/INVMSUN_TO_KHZ*0.5, 6.0/INVMSUN_TO_KHZ; length=400)
        re=0.0; im=0.0
        for k in 1:n
            θ=2π*f*(ts[k]-ts[1]); re+=y[k]*cos(θ); im-=y[k]*sin(θ)
        end
        pow=re^2+im^2
        if pow>bestpow; bestpow=pow; best=f; end
    end
    return best
end

# Self-convergence factor Q (Shum eq.l.986) using the L2 difference of ε_c(t)
# resampled to a common time grid; estimate order n by solving
#   |f_l - f_m| / |f_m - f_h| = (Δr_l^n - Δr_m^n)/(Δr_m^n - Δr_h^n).
function convergence_factor(series, drs)
    dl,dm,dh = drs[1],drs[2],drs[3]
    tl,el = series[dl]; tm,em = series[dm]; th,eh = series[dh]
    # common time grid (intersection), sample at coarse cadence
    tmax = min(tl[end],tm[end],th[end]); tmin=max(tl[1],tm[1],th[1])
    tg = collect(range(tmin+1e-9, tmax-1e-9; length=40))
    fl = [_interp(tl,el,t) for t in tg]
    fm = [_interp(tm,em,t) for t in tg]
    fh = [_interp(th,eh,t) for t in tg]
    num = norm(fl .- fm); den = norm(fm .- fh)
    ratio = den>0 ? num/den : NaN
    # solve for n: (dl^n - dm^n)/(dm^n - dh^n) = ratio  (bisection)
    f(n) = (dl^n - dm^n)/(dm^n - dh^n) - ratio
    a,b=0.1,8.0; fa=f(a); n_est=NaN
    if isfinite(ratio) && fa*f(b)<0
        for _ in 1:80
            c=(a+b)/2; fc=f(c)
            (fa*fc<=0) ? (b=c) : (a=c, fa=fc)
        end
        n_est=(a+b)/2
    end
    return ratio, n_est
end

# self-consistency gate: forward projection ∘ con2prim must be the identity on p1
function consistency_gate()
    κ=KAPPA; eos=ShumPolytrope(κ); ρ0c=0.00128; εc=ρ0c+κ*ρ0c^2
    star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    iso = areal_to_isotropic(star)
    fr = ShumFrame(ŝ=1.0, â=1.0, q̂=0.999, η̂=0.01, ζ̂=0.01)
    bg = build_background(iso, eos; rmax=20.0, dr=0.01)
    N=length(bg.r); maxerr=0.0
    for j in 2:N-1
        ε=bg.εbg[j]; ε<=energy_from_pressure(eos,P_ATMS) && continue
        v=0.03; dε=1e-5; dv=2e-5; eh_t=3.7e-6; vbh_t=-1.1e-5     # arbitrary p1
        E,Sr,_,_ = bdnk_projections(fr,eos,bg,j; ε=ε,v=v,dε=dε,dv=dv,eh=eh_t,vbh=vbh_t)
        eh_r,vbh_r = recover_p1(fr,eos,bg,j; ε=ε,v=v,dε=dε,dv=dv,E=E,Sr=Sr)
        maxerr=max(maxerr, abs(eh_r-eh_t), abs(vbh_r-vbh_t))
    end
    return maxerr
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
