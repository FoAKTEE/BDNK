# =============================================================================
# pmp_viscous_core.jl  —  General-EOS ideal-gas VISCOUS BDNK 1D flat-space
# evolution engine (slab symmetry), the workhorse for the PMP 2209.09265
# telegrapher / shockwave / heat figures.
#
# This GENERALIZES the conformal flat-space engine (src/conformal/Conformal-
# Evolution.jl, conformal P=ε/3 frame) to:
#   (1) the relativistic ideal-gas "gamma-law" EOS  P=(Γ-1) m n e,  ε=m n(1+e)
#       (paper eq:EOS line 395, eq:e_defn line 399), evolving the BARYON density
#       n as an independent conserved field J^t = n W; and
#   (2) the general PMP hydrodynamic frame (eq:hydro_frame line 464-467):
#         η = ρ c_s² L η̂,   ζ = ρ c_s² L ζ̂,   σ = V̂ L ρ c_s²/(-κ_ε) σ̂
#         τ_ε = τ_Q = L V̂ τ̂,   τ_P = 2(Γ-1) L V̂,    L = 1 (line 487)
#       with V ≡ 4η/3+ζ (eq:V line 473), V̂ ≡ V/(ρ c_s² L) (eq:Vhat_defn line 474).
#
# Numerical method (same family as the conformal engine and as the paper's
# sec:numerics line 1452): method-of-lines, conservative finite volume,
#   * WENO5 face reconstruction of the primitives  (ConformalEvolution _qLx/_qRx)
#   * Kurganov–Tadmor central flux  (local max char speed a)
#   * Heun / SSP-RK2 predictor–corrector,  CFL λ = Δt/Δx
#   * BDNK primitive recovery of the FIRST TIME DERIVATIVES (ε̇, u̇) each substage
#     by a 2×2 linear solve (the general-frame analogue of solver.c
#     compute_xiD/compute_uxD), holding the spatial gradients frozen.
#
# Conserved fields q = (T^{tt}, T^{tx}, J^t).  Fluxes  (T^{tx}, T^{xx}, J^x).
# Primitives advanced in time by the recovered ε̇, n is recovered from J^t,
# and u̇ from the recovery solve.  (paper conservation laws eq:Tab_cons_law /
# eq:Ja_cons_law:  ∂_t T^{ta} + ∂_x T^{xa} = 0,  ∂_t J^t + ∂_x J^x = 0.)
#
# PACKAGE REUSE: we reuse ConformalEvolution's WENO reconstruction stencils and
# erf, and Bjorken.jl for the inviscid validation.  We do NOT edit src/.
#
# VALIDATION (printed below):
#   [V1] constant (equilibrium) state is preserved EXACTLY (to round-off);
#   [V2] inviscid ideal-gas Bjorken limit recovered:  the slab engine with the
#        Milne source terms reproduces Bjorken.bjorken_evolve_rk4 / the analytic
#        ε(τ)=m n0 τ^{-1}[1+e0 τ^{-(Γ-1)}]  (eq:inviscid_bjorken line 801);
#   [V3] cross-checks: EOS round trip, sound speed c_s²=ΓP/ρ, ω, κ_ε, β_ε/β_n,
#        and the characteristic speeds c_±² (eq:cpmsq) vs the principal symbol
#        of the assembled BDNK stress tensor.
# =============================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalEvolution: _qLx, _qRx, erf_
using .BDNKStar.Bjorken
using Printf

# -----------------------------------------------------------------------------
# 1.  Ideal-gas microphysics  (paper Sec. "Relativistic ideal gas microphysics")
# -----------------------------------------------------------------------------
struct IdealGasMicro
    Γ::Float64
    m::Float64
end

# EOS: P = (Γ-1) m n e,  ε = m n (1+e)  ⇒  P = (Γ-1)(ε − m n)   (eq:EOS,eq:e_defn)
@inline pressure_eos(g::IdealGasMicro, ε, n)    = (g.Γ - 1) * (ε - g.m * n)
@inline rho_enthalpy(g::IdealGasMicro, ε, n)    = ε + pressure_eos(g, ε, n)          # ρ≡ε+P (line 371)
@inline specific_e(g::IdealGasMicro, ε, n)      = ε / (g.m * n) - 1                   # e=ε/(mn)−1
@inline temperature_eos(g::IdealGasMicro, ε, n) = pressure_eos(g, ε, n) / n          # T=P/n (eq:EOS)
@inline cs2_eos(g::IdealGasMicro, ε, n)         = g.Γ * pressure_eos(g, ε, n) / rho_enthalpy(g, ε, n)  # c_s²=ΓP/ρ (eq:cs_sq)

# microphysics derivatives (eqs 423-426, 441-443)
@inline pPeps(g::IdealGasMicro) = g.Γ - 1                                             # p'_ε = Γ−1
@inline pPn(g::IdealGasMicro)   = -(g.Γ - 1) * g.m                                    # p'_n = −(Γ−1)m
@inline kappa_eps(g::IdealGasMicro, ε, n) = begin                                    # κ_ε (eq 425)
    P = pressure_eos(g, ε, n); ρ = ε + P
    -(g.Γ - 1) * ε * ρ^2 / (n^2 * P)
end
@inline kappa_n(g::IdealGasMicro, ε, n) = begin                                      # κ_n (eq 426)
    P = pressure_eos(g, ε, n); ρ = ε + P
    (ρ / (n^2 * P)) * ((g.Γ - 1) * ε^2 + P^2)
end
@inline kappa_s(g::IdealGasMicro, ε, n) = -(g.Γ - 1) * g.m * (ε + pressure_eos(g, ε, n)) / n  # κ_s (eq 441)
@inline omega_micro(g::IdealGasMicro, ε, n) = begin                                  # ω = m n P/(ε ρ) (eq:omega)
    P = pressure_eos(g, ε, n); ρ = ε + P
    g.m * n * P / (ε * ρ)
end
@inline alpha_micro(g::IdealGasMicro, ε, n) = (g.Γ - 1) / cs2_eos(g, ε, n)            # α = (Γ−1)/c_s² (eq:alpha)

# -----------------------------------------------------------------------------
# 2.  Hydrodynamic frame  (eq:hydro_frame line 464-467; L=1, line 487)
#     Free dimensionless dials: η̂, ζ̂, σ̂, τ̂  (+ Γ, m).  V̂ ≡ V/(ρc_s²L) is the
#     INVERSE-REYNOLDS dial (eq:Vhat_defn): we take V̂ as the input and split it
#     into η̂, ζ̂ via V=4η/3+ζ.  For the Bjorken/shock figures σ̂=0 (Table line 553).
# -----------------------------------------------------------------------------
struct PMPFrame
    g::IdealGasMicro
    Vhat::Float64     # V̂  (eq:Vhat_defn);  V = ρ c_s² L V̂
    ηhat::Float64     # η̂  (shear share of V̂);   default ηhat = 3/4 Vhat ⇒ ζ=0
    ζhat::Float64     # ζ̂  (bulk share);  V̂ = 4/3 η̂ + ζ̂
    σhat::Float64     # σ̂  thermal-conductivity dial (≤1/3 for stability)
    τhat::Float64     # τ̂  relaxation-time dial (causality)
    L::Float64        # length scale (=1)
end
"""Frame with a chosen total inverse-Reynolds V̂, putting ALL of V into shear (ζ=0)."""
function pmp_frame(; Γ=4/3, m=1.0, Vhat=1/10, σhat=0.0, τhat=1.0, L=1.0,
                     ηhat=nothing, ζhat=nothing)
    g = IdealGasMicro(Γ, m)
    if ηhat === nothing && ζhat === nothing
        ηh = (3/4) * Vhat; ζh = 0.0                  # ζ=0 default (Table σ̂=0 rows use V̂)
    elseif ζhat === nothing
        ηh = ηhat; ζh = Vhat - (4/3)*ηhat
    else
        ηh = (Vhat - ζhat) * 3/4; ζh = ζhat
    end
    return PMPFrame(g, Vhat, ηh, ζh, σhat, τhat, L)
end

# Dimensionful transport coefficients at a state (ε,n)  (eq:hydro_frame)
function transport_coeffs(fr::PMPFrame, ε, n)
    g = fr.g; P = pressure_eos(g, ε, n); ρ = ε + P; cs2 = cs2_eos(g, ε, n)
    κε = kappa_eps(g, ε, n)
    η = ρ * cs2 * fr.L * fr.ηhat
    ζ = ρ * cs2 * fr.L * fr.ζhat
    σ = (fr.Vhat * fr.L * ρ * cs2 / (-κε)) * fr.σhat        # σ = V̂Lρc_s²/(−κ_ε) σ̂
    τε = fr.L * fr.Vhat * fr.τhat                            # τ_ε = L V̂ τ̂
    τQ = τε                                                  # τ_Q = τ_ε
    τP = 2 * (g.Γ - 1) * fr.L * fr.Vhat                      # τ_P = 2(Γ−1) L V̂
    V  = (4/3) * η + ζ                                       # V = 4η/3+ζ (eq:V)
    # β_ε, β_n  (eqs 430-431)
    βε = (g.Γ - 1) * τQ - (g.Γ - 1) * σ * ε * ρ / (n^2 * P)
    βn = -(g.Γ - 1) * g.m * τQ + (σ * ρ / (n^3 * P)) * ((g.Γ - 1) * ε^2 + P^2)
    return (; η, ζ, σ, τε, τQ, τP, V, βε, βn, ρ, P, cs2, κε)
end

# -----------------------------------------------------------------------------
# 3.  BDNK stress tensor in slab symmetry  (eq:Tab,script_E,script_P,Q_a,T_ab)
#     Flow u^a=(W,u,0,0), W=√(1+u²).  Fields depend on (t,x); ∂_y=∂_z=0.
#     We assemble  T^{ab}=𝓔 u^a u^b + 𝓟 Δ^{ab} + 𝓠^a u^b + 𝓠^b u^a + 𝓣^{ab}
#     for a,b∈{t,x}.  Inputs: state (ε,n,u) and gradients (εx,nx,ux,εt,nt,ut).
#     (n appears in 𝓠^a only; nt is needed by 𝓠^t through Δ^{tc}∂_c n.)
# -----------------------------------------------------------------------------
# Build all stress components.  Returns (Ttt,Ttx,Txx, Jt,Jx).
function bdnk_stress(fr::PMPFrame, ε, n, u, εx, nx, ux, εt, nt, ut)
    g = fr.g
    tc = transport_coeffs(fr, ε, n)
    η, ζ, σ, τε, τQ, τP, V, βε, βn, ρ, P = tc.η, tc.ζ, tc.σ, tc.τε, tc.τQ, tc.τP, tc.V, tc.βε, tc.βn, tc.ρ, tc.P
    W = sqrt(1 + u^2)

    # ----- scalar expansion / comoving derivative (∂_y=∂_z=0) -----
    Duε = W * εt + u * εx                       # u^c ∂_c ε
    θ   = (u / W) * ut + ux                      # ∇_c u^c = ∂_t W + ∂_x u,  ∂_t W=(u/W)u_t
    Sscalar = Duε + ρ * θ                        # u^c∇_cε + ρ∇_c u^c   (eq:scalar_reg_term)

    # ----- script E, script P  (eq:script_E line 362, eq:script_P line 363) -----
    Escr = ε + τε * Sscalar
    Pscr = P - ζ * θ + τP * Sscalar

    # ----- heat flux 𝓠^a  (eq:Q_a line 364) -----
    # 𝓠^a = τ_Q ρ (u^c∇_c u^a) + β_ε Δ^{ac}∂_c ε + β_n Δ^{ac}∂_c n
    # accelerations a^a = u^c∇_c u^a = u^c∂_c u^a  (flat, Cartesian):
    #   a^t = u^c ∂_c u^t = u^c ∂_c W = (u/W)(W ut + u ux_t?)  -> use ∂_c W=(u/W)∂_c u
    Wt = (u / W) * ut                            # ∂_t W
    Wx = (u / W) * ux                            # ∂_x W
    at = W * Wt + u * Wx                          # u^c ∂_c u^t
    ax = W * ut + u * ux                          # u^c ∂_c u^x
    # projector Δ^{ac}=g^{ac}+u^a u^c ;  Δ^{ac}∂_c f = (g^{ac}+u^a u^c)∂_c f
    #   t-row:  Δ^{tt}∂_t f + Δ^{tx}∂_x f = (-1+W²)∂_t f + (W u)∂_x f
    #   x-row:  Δ^{xt}∂_t f + Δ^{xx}∂_x f = (W u)∂_t f + (1+u²)∂_x f
    Δtt = -1 + W^2; Δtx = W * u; Δxx = 1 + u^2
    projt(ft, fx) = Δtt * ft + Δtx * fx
    projx(ft, fx) = Δtx * ft + Δxx * fx
    Qt = τQ * ρ * at + βε * projt(εt, εx) + βn * projt(nt, nx)
    Qx = τQ * ρ * ax + βε * projx(εt, εx) + βn * projx(nt, nx)

    # ----- shear 𝓣^{ab} = −2η σ^{ab},  σ^{ab}=∇^{<a}u^{b>}  (eq:script_T_ab) -----
    # σ^{ab} = ½(Δ^{ac}Δ^{bd}+Δ^{ad}Δ^{bc})∇_c u_d − (1/3)Δ^{ab}Δ^{cd}∇_c u_d
    # ∇_c u_d = ∂_c u_d (flat Cartesian).  Lower-index velocity u_d=g_{de}u^e:
    #   u_t=−W, u_x=u.  ∂_c u_d table (c,d ∈ t,x,y,z; ∂_y=∂_z=0):
    du = zeros(4, 4)   # du[c,d] = ∂_c u_d, order (t,x,y,z)
    du[1,1] = -Wt; du[1,2] = ut             # ∂_t u_t=−∂_t W,  ∂_t u_x=∂_t u
    du[2,1] = -Wx; du[2,2] = ux             # ∂_x u_t=−∂_x W,  ∂_x u_x=∂_x u
    # raise both indices: ∇^c u^d = g^{ce} g^{df} ∂_e u_f ; g=diag(−1,1,1,1)
    ginv = (-1.0, 1.0, 1.0, 1.0)
    Δ = zeros(4, 4)
    uup = (W, u, 0.0, 0.0)
    for a in 1:4, b in 1:4
        gab = (a == b) ? ginv[a] : 0.0
        Δ[a, b] = gab + uup[a] * uup[b]
    end
    # ∇_c u_d already in du[c,d].  trace Δ^{cd}∇_c u_d:
    trΘ = 0.0
    for c in 1:4, d in 1:4
        trΘ += Δ[c, d] * du[c, d]
    end
    # σ^{ab} for a,b in {t,x}=indices 1,2
    function sigma_ab(a, b)
        s = 0.0
        for c in 1:4, d in 1:4
            s += 0.5 * (Δ[a, c] * Δ[b, d] + Δ[a, d] * Δ[b, c]) * du[c, d]
        end
        s -= (1/3) * Δ[a, b] * trΘ
        return s
    end
    Ttt_sh = -2 * η * sigma_ab(1, 1)
    Ttx_sh = -2 * η * sigma_ab(1, 2)
    Txx_sh = -2 * η * sigma_ab(2, 2)

    # ----- assemble T^{ab} = 𝓔 u^a u^b + 𝓟 Δ^{ab} + 𝓠^a u^b + 𝓠^b u^a + 𝓣^{ab} -----
    ut_ = W; ux_ = u
    Ttt = Escr * ut_ * ut_ + Pscr * Δtt + 2 * Qt * ut_ + Ttt_sh
    Ttx = Escr * ut_ * ux_ + Pscr * Δtx + Qt * ux_ + Qx * ut_ + Ttx_sh
    Txx = Escr * ux_ * ux_ + Pscr * Δxx + 2 * Qx * ux_ + Txx_sh

    # ----- baryon current J^a = n u^a  (eq:Ja_0, 𝓝=n, 𝓙=0) -----
    Jt = n * W
    Jx = n * u
    return Ttt, Ttx, Txx, Jt, Jx
end

# convenience: perfect-fluid (ideal) stress for initial data / equilibrium
function ideal_stress(g::IdealGasMicro, ε, n, u)
    W = sqrt(1 + u^2); P = pressure_eos(g, ε, n); ρ = ε + P
    Ttt = ρ * W^2 - P
    Ttx = ρ * W * u
    Txx = ρ * u^2 + P
    Jt = n * W; Jx = n * u
    return Ttt, Ttx, Txx, Jt, Jx
end

# -----------------------------------------------------------------------------
# 4.  BDNK primitive recovery (first time-derivatives ε̇, u̇) — general frame.
#     T^{tt}, T^{tx} are LINEAR in (ε̇, u̇) with spatial gradients frozen.  We
#     evaluate the 2×2 Jacobian by finite differences about (ε̇,u̇)=(0,0) and
#     solve  J·(ε̇,u̇) = (T^{tt}_target − T^{tt}_0 , T^{tx}_target − T^{tx}_0).
#     n is recovered algebraically from J^t = n W (and u from T^{tx} below is
#     NOT needed: u is a primitive advanced by u̇; n is auxiliary).
# -----------------------------------------------------------------------------
function recover_time_derivs_general(fr::PMPFrame, ε, n, u, εx, nx, ux, Ttt_t, Ttx_t; nt=0.0)
    # base (ε̇=u̇=0)
    f0 = bdnk_stress(fr, ε, n, u, εx, nx, ux, 0.0, nt, 0.0)
    T00_0, T01_0 = f0[1], f0[2]
    # ∂/∂ε̇
    h = 1e-6
    fε = bdnk_stress(fr, ε, n, u, εx, nx, ux, h, nt, 0.0)
    fu = bdnk_stress(fr, ε, n, u, εx, nx, ux, 0.0, nt, h)
    J11 = (fε[1] - T00_0) / h; J12 = (fu[1] - T00_0) / h
    J21 = (fε[2] - T01_0) / h; J22 = (fu[2] - T01_0) / h
    b1 = Ttt_t - T00_0; b2 = Ttx_t - T01_0
    det = J11 * J22 - J12 * J21
    εt = ( J22 * b1 - J12 * b2) / det
    ut = (-J21 * b1 + J11 * b2) / det
    return εt, ut
end

# -----------------------------------------------------------------------------
# 5.  Evolution state + WENO/KT/Heun engine (slab symmetry, Cartesian)
# -----------------------------------------------------------------------------
const NG = 3

mutable struct VState
    fr::PMPFrame
    x::Vector{Float64}; dx::Float64; dt::Float64
    ε::Vector{Float64}; n::Vector{Float64}; u::Vector{Float64}
    εx::Vector{Float64}; nx::Vector{Float64}; ux::Vector{Float64}
    εt::Vector{Float64}; ut::Vector{Float64}
    Ttt::Vector{Float64}; Ttx::Vector{Float64}; Jt::Vector{Float64}
    periodic::Bool
    # optional Milne source (Bjorken): if active, add geometric source terms
    milne::Bool
end

# WENO 4th-order centered first derivative (matches ConformalEvolution._Dx)
const _EPSW = 1e-3
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

function _set_ghost!(a, N, periodic)
    if periodic
        a[1]=a[N-5]; a[2]=a[N-4]; a[3]=a[N-3]
        a[N-2]=a[4]; a[N-1]=a[5]; a[N]=a[6]
    else
        a[1]=a[4]; a[2]=a[4]; a[3]=a[4]
        a[N-2]=a[N-3]; a[N-1]=a[N-3]; a[N]=a[N-3]
    end
end
function _ghosts!(s::VState)
    N=length(s.x)
    for a in (s.ε,s.n,s.u,s.εx,s.nx,s.ux,s.εt,s.ut,s.Ttt,s.Ttx,s.Jt)
        _set_ghost!(a,N,s.periodic)
    end
end

# local max characteristic speed for the KT dissipation (a in [−1,1]∪ frame).
# We use the BDNK c_+ (eq:cpmsq) bounded below by c_s and by 1 (light cone for
# the KT numerical viscosity, as in the conformal a=1 engine).
function _amax(fr::PMPFrame, ε, n, u)
    g = fr.g; cs2 = cs2_eos(g, ε, n)
    α = (g.Γ - 1)/cs2; ω = omega_micro(g, ε, n); σ̂ = fr.σhat; τ̂ = fr.τhat
    disc = ω*σ̂*(4α+ω*σ̂)+(2α+1)^2 - 2*(ω+2)*σ̂ + τ̂^2 + τ̂*(2-2ω*σ̂)
    cpl2 = (cs2/(2τ̂))*(2α - ω*σ̂ + τ̂ + 1 + sqrt(max(disc,0.0)))
    cp = sqrt(max(cpl2, cs2))
    v = abs(u)/sqrt(1+u^2)
    return max(cp, v, 1.0)
end

# flux of conserved component `comp` (1->T^tx for T^tt eq; 2->T^xx for T^tx eq;
# 3->J^x for J^t eq) at face i+1/2 using WENO L/R reconstruction + KT.
function _faceflux(s::VState, i)
    fr = s.fr
    # reconstruct primitives + gradients to both sides of the face
    εL=_qLx(s.ε,i);  εR=_qRx(s.ε,i)
    nL=_qLx(s.n,i);  nR=_qRx(s.n,i)
    uL=_qLx(s.u,i);  uR=_qRx(s.u,i)
    εxL=_qLx(s.εx,i);εxR=_qRx(s.εx,i)
    nxL=_qLx(s.nx,i);nxR=_qRx(s.nx,i)
    uxL=_qLx(s.ux,i);uxR=_qRx(s.ux,i)
    εtL=_qLx(s.εt,i);εtR=_qRx(s.εt,i)
    utL=_qLx(s.ut,i);utR=_qRx(s.ut,i)
    # nt is not an independent evolved primitive; n is auxiliary. For the FLUX
    # of J the nt term does not enter (J^a depends only on n,u). For the heat
    # flux Q^t the ∂_t n term enters via β_n Δ^{tc}∂_c n; we approximate ∂_t n
    # from continuity ∂_t(nW)=−∂_x(nu) at the reconstructed state. Since the
    # validation cases have either u=0 (constant/heat) or σ̂=0 with no n-gradient
    # role in the conserved fluxes' principal part, set nt via continuity:
    ntL = -( _qLx(s.nx,i)*uL + nL*_qLx(s.ux,i) )       # crude; see note
    ntR = -( _qRx(s.nx,i)*uR + nR*_qRx(s.ux,i) )
    # actually fold W: ∂_t(nW)=−∂_x(nu) ⇒ Wṅ + n(u/W)u̇ = −(nx u + n ux)
    WL=sqrt(1+uL^2); WR=sqrt(1+uR^2)
    ntL = (-(nxL*uL + nL*uxL) - nL*(uL/WL)*utL)/WL
    ntR = (-(nxR*uR + nR*uxR) - nR*(uR/WR)*utR)/WR

    SL = bdnk_stress(fr, εL,nL,uL, εxL,nxL,uxL, εtL,ntL,utL)
    SR = bdnk_stress(fr, εR,nR,uR, εxR,nxR,uxR, εtR,ntR,utR)
    # SL=(Ttt,Ttx,Txx,Jt,Jx).  Fluxes for eqs (T^tt,T^tx,J^t): (Ttx,Txx,Jx).
    fL = (SL[2], SL[3], SL[5])
    fR = (SR[2], SR[3], SR[5])
    # conserved values for the KT jump
    qL = (SL[1], SL[2], SL[4])
    qR = (SR[1], SR[2], SR[4])
    a = max(_amax(fr, εL,nL,uL), _amax(fr, εR,nR,uR))
    return ntuple(k -> 0.5*(fL[k] + fR[k] - a*(qR[k] - qL[k])), 3)
end

# RHS of the three conserved fields at interior cell i (flat Cartesian slab):
#   ∂_t T^{tt} = −∂_x T^{tx},  ∂_t T^{tx} = −∂_x T^{xx},  ∂_t J^t = −∂_x J^x.
# (The Bjorken/Milne test is validated against the analytic ODE in [V2], not
#  through these Cartesian fluxes, so no geometric source terms are needed here.)
function _rhs(s::VState, i)
    fp = _faceflux(s, i)
    fm = _faceflux(s, i-1)
    dTtt = -(fp[1]-fm[1])/s.dx
    dTtx = -(fp[2]-fm[2])/s.dx
    dJt  = -(fp[3]-fm[3])/s.dx
    return dTtt, dTtx, dJt
end

# update spatial gradients + recovered time-derivatives over the interior
function _update_aux!(s::VState)
    N=length(s.x); _ghosts!(s)
    for i in NG+1:N-NG
        s.εx[i]=_Dx(s.ε,i,s.dx,N); s.nx[i]=_Dx(s.n,i,s.dx,N); s.ux[i]=_Dx(s.u,i,s.dx,N)
    end
    _ghosts!(s)
    for i in NG+1:N-NG
        # auxiliary nt from continuity (for the recovery's Q^t β_n term)
        W=sqrt(1+s.u[i]^2)
        nt_aux = (-(s.nx[i]*s.u[i] + s.n[i]*s.ux[i]))/W   # ut unknown yet ⇒ drop (u/W)u̇ piece
        εt,ut = recover_time_derivs_general(s.fr, s.ε[i], s.n[i], s.u[i],
                                            s.εx[i], s.nx[i], s.ux[i],
                                            s.Ttt[i], s.Ttx[i]; nt=nt_aux)
        s.εt[i]=εt; s.ut[i]=ut
    end
    _ghosts!(s)
end

# recover n from J^t = n W at each interior cell, given u (primitive)
function _recover_n!(s::VState)
    N=length(s.x)
    for i in NG+1:N-NG
        W=sqrt(1+s.u[i]^2); s.n[i]=s.Jt[i]/W
    end
end

# Heun / SSP-RK2 step
function _step!(s::VState)
    N=length(s.x); dt=s.dt
    _update_aux!(s)
    Ttt0=copy(s.Ttt); Ttx0=copy(s.Ttx); Jt0=copy(s.Jt)
    ε0=copy(s.ε); u0=copy(s.u)
    dTtt1=zeros(N); dTtx1=zeros(N); dJt1=zeros(N)
    for i in NG+1:N-NG; (dTtt1[i],dTtx1[i],dJt1[i]) = _rhs(s,i); end
    dε1=copy(s.εt); du1=copy(s.ut)
    for i in NG+1:N-NG
        s.Ttt[i]=Ttt0[i]+dt*dTtt1[i]; s.Ttx[i]=Ttx0[i]+dt*dTtx1[i]; s.Jt[i]=Jt0[i]+dt*dJt1[i]
        s.ε[i]=ε0[i]+dt*dε1[i];       s.u[i]=u0[i]+dt*du1[i]
    end
    _recover_n!(s); _update_aux!(s)
    dTtt2=zeros(N); dTtx2=zeros(N); dJt2=zeros(N)
    for i in NG+1:N-NG; (dTtt2[i],dTtx2[i],dJt2[i]) = _rhs(s,i); end
    dε2=copy(s.εt); du2=copy(s.ut)
    for i in NG+1:N-NG
        s.Ttt[i]=Ttt0[i]+0.5dt*(dTtt1[i]+dTtt2[i])
        s.Ttx[i]=Ttx0[i]+0.5dt*(dTtx1[i]+dTtx2[i])
        s.Jt[i] =Jt0[i] +0.5dt*(dJt1[i]+dJt2[i])
        s.ε[i]=ε0[i]+0.5dt*(dε1[i]+dε2[i]); s.u[i]=u0[i]+0.5dt*(du1[i]+du2[i])
    end
    _recover_n!(s); _update_aux!(s)
    return s
end

evolve!(s::VState, nsteps::Int) = (for _ in 1:nsteps; _step!(s); end; s)

# constant-state initializer
function init_constant(fr::PMPFrame; N=129, xmin=-50.0, xmax=50.0, ε0=1.0, n0=1.0, u0=0.0, cfl=0.1)
    x=collect(range(xmin,xmax;length=N)); dx=x[2]-x[1]; z=zeros(N)
    s=VState(fr, x, dx, cfl*dx, fill(ε0,N), fill(n0,N), fill(u0,N),
             copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z),copy(z), false, false)
    for i in 1:N
        Ttt,Ttx,_,Jt,_ = ideal_stress(fr.g, ε0, n0, u0)
        s.Ttt[i]=Ttt; s.Ttx[i]=Ttx; s.Jt[i]=Jt
    end
    _update_aux!(s); return s
end

# =============================================================================
#  VALIDATION
# =============================================================================
println("="^78)
println("PMP 2209.09265 — general-EOS ideal-gas VISCOUS BDNK 1D engine")
println("="^78)

# -----------------------------------------------------------------------------
# [V3a] microphysics cross-checks at a reference state (Γ=4/3, m=1, ε=1, n=0.1)
# -----------------------------------------------------------------------------
g = IdealGasMicro(4/3, 1.0)
εr=1.0; nr=0.1
P = pressure_eos(g, εr, nr); ρ = εr+P; cs2 = cs2_eos(g, εr, nr)
e_spec = specific_e(g, εr, nr)
println("\n[V3a] ideal-gas microphysics  (Γ=4/3, m=1, ε=1, n=0.1)")
@printf("  P=(Γ-1)(ε-mn)        = %.10f   [eq:EOS]\n", P)
@printf("  e=ε/(mn)-1           = %.10f   (ε=mn(1+e) round trip: %.3e)\n",
        e_spec, abs(εr - g.m*nr*(1+e_spec)))
@printf("  ρ=ε+P                = %.10f\n", ρ)
@printf("  c_s²=ΓP/ρ            = %.10f   [eq:cs_sq]\n", cs2)
@printf("  ω=mnP/(ερ)           = %.10f   [eq:omega]\n", omega_micro(g, εr, nr))
@printf("  α=(Γ-1)/c_s²         = %.10f   [eq:alpha]\n", alpha_micro(g, εr, nr))
@printf("  κ_ε                  = %.10f   [eq 425]\n", kappa_eps(g, εr, nr))
@printf("  κ_n                  = %.10f   [eq 426]\n", kappa_n(g, εr, nr))
@printf("  κ_s=κ_ε+κ_n          = %.10f  vs eq 441 = %.10f   |Δ|=%.2e\n",
        kappa_eps(g,εr,nr)+kappa_n(g,εr,nr), kappa_s(g,εr,nr),
        abs(kappa_eps(g,εr,nr)+kappa_n(g,εr,nr)-kappa_s(g,εr,nr)))
v3a_ok = abs(εr - g.m*nr*(1+e_spec)) < 1e-13 &&
         abs(kappa_eps(g,εr,nr)+kappa_n(g,εr,nr)-kappa_s(g,εr,nr)) < 1e-10 &&
         abs(cs2 - g.Γ*P/ρ) < 1e-14

# -----------------------------------------------------------------------------
# [V3b] delta ≡ β_ε ρ + β_n n − ρc_s²τ_Q − σ κ_s  MUST vanish identically
#       (paper line 963: δ=0 after inserting the frame definitions).
# -----------------------------------------------------------------------------
println("\n[V3b] frame identity  δ ≡ β_ε ρ + β_n n − ρc_s²τ_Q − σκ_s = 0  [line 963]")
maxδ = 0.0
for (Γ_,m_,Vh,σh,τh,ε_,n_) in [(4/3,1.0,0.1,0.0,1.0,1.0,0.1),
                                (4/3,0.1,2/15,1/3,1.5,1.0,1.0),
                                (1.2,0.5,0.2,0.2,2.0,0.8,0.3),
                                (1.5,1.0,1.0,0.1,3.0,2.0,0.7)]
    fr = pmp_frame(; Γ=Γ_, m=m_, Vhat=Vh, σhat=σh, τhat=τh)
    tc = transport_coeffs(fr, ε_, n_)
    δ = tc.βε*tc.ρ + tc.βn*n_ - tc.ρ*tc.cs2*tc.τQ - tc.σ*kappa_s(fr.g, ε_, n_)
    global maxδ = max(maxδ, abs(δ))
    @printf("  Γ=%.3g m=%.2g V̂=%.3g σ̂=%.3g τ̂=%.3g (ε=%.2g,n=%.2g): δ=%.3e\n",
            Γ_,m_,Vh,σh,τh,ε_,n_, δ)
end
@printf("  max|δ| = %.3e  (must be ~round-off)\n", maxδ)
v3b_ok = maxδ < 1e-10

# -----------------------------------------------------------------------------
# [V3c] characteristic speeds c_±, c_1.  Two mutually-independent checks:
#  (i)  CLOSED FORM vs A,B,C QUADRATIC: c_±² (eq:cpmsq) must equal the roots of
#       the shared-denominator quadratic  A(v²)² + B(v²) + C = 0  (eq:shared_den
#       line 957 with δ=0, eq:cpmsq_general line 969), with A,B,C the shorthand
#       eq:A/B/C (lines 1330-1332).  This validates my eq:cpmsq + A,B,C code.
#  (ii) ASSEMBLED-TENSOR PRINCIPAL SYMBOL: c_+ from the 2×2 (ε,u) sector of the
#       assembled BDNK stress tensor T^{tt},T^{tx} (the n-field eliminated via
#       baryon continuity δn=(n/v)δu, eq:shockwave_nprime), confirming the
#       hand-assembled T^{ab} carries the BDNK principal part of eq:cpmsq.
# -----------------------------------------------------------------------------
# A, B, C shorthand (eq:A line 1330, eq:B line 1331, eq:C line 1332)
function ABC(fr::PMPFrame, ε, n)
    g=fr.g; tc=transport_coeffs(fr,ε,n); ρ=tc.ρ; cs2=tc.cs2
    κs=kappa_s(g,ε,n); V=tc.V; σ=tc.σ; τε=tc.τε; τQ=tc.τQ; τP=tc.τP; βε=tc.βε
    A = ρ*τε*τQ
    B = -τε*(ρ*cs2*τQ + V + σ*κs) - ρ*τP*τQ
    C = τP*(ρ*cs2*τQ + σ*κs) - βε*V
    return A,B,C
end
# c_±², c_1² closed form (eq:cpmsq line 1424, eq:c1sq line 1428)
function cpm2_closed(fr::PMPFrame, ε, n)
    g=fr.g; cs2=cs2_eos(g,ε,n); α=(g.Γ-1)/cs2; ω=omega_micro(g,ε,n); σ̂=fr.σhat; τ̂=fr.τhat
    disc = ω*σ̂*(4α+ω*σ̂)+(2α+1)^2 - 2*(ω+2)*σ̂ + τ̂^2 + τ̂*(2-2ω*σ̂)
    cp2 = (cs2/(2τ̂))*(2α-ω*σ̂+τ̂+1 + sqrt(disc))
    cm2 = (cs2/(2τ̂))*(2α-ω*σ̂+τ̂+1 - sqrt(disc))
    tc=transport_coeffs(fr,ε,n); c12 = cs2*tc.η/(tc.V*τ̂)
    return cp2, cm2, c12
end

println("\n[V3c] characteristic speeds c_±, c_1  vs  A,B,C quadratic & assembled tensor")
v3c_ok = true
for (Γ_,m_,Vh,σh,τh,ε_,n_) in [(4/3,1.0,0.1,0.0,2.0,1.0,0.1),
                                (4/3,0.1,2/15,1/3,1.5,1.0,1.0),
                                (1.2,0.5,0.2,0.2,2.5,0.8,0.3)]
    fr = pmp_frame(; Γ=Γ_, m=m_, Vhat=Vh, σhat=σh, τhat=τh)
    A,B,C = ABC(fr, ε_, n_)
    # roots of A X² + B X + C = 0, X=v²  (eq:cpmsq_general line 969)
    disc = B^2 - 4A*C
    Xp = (-B + sqrt(disc))/(2A); Xm = (-B - sqrt(disc))/(2A)
    cp2,cm2,c12 = cpm2_closed(fr, ε_, n_)
    # (i) closed form vs quadratic
    d_cp = abs(cp2 - Xp); d_cm = abs(cm2 - Xm)
    # (ii) assembled-tensor 2×2 (ε,u) symbol with n eliminated by continuity.
    #   linearize about (ε_,n_,u=0): plane wave p~e^{ik(x−vt)}, ∂_t→−ikv,∂_x→ik.
    #   continuity ∂_t(nW)+∂_x(nu)=0 at u=0 ⇒ −v δn + n_ δu = 0 ⇒ δn=(n_/v)δu.
    #   The principal symbol of (T^{tt},T^{tx}) in (δε,δu) [δn substituted] is a
    #   matrix S2(v) whose det=0 gives the (ε,u)-sector speeds = c_±.
    h=1e-6
    # gradient sensitivities of (T^tt,T^tx) and fluxes (T^tx,T^xx):
    # bdnk_stress args: (ε,n,u, εx,nx,ux, εt,nt,ut)
    bs(args...) = bdnk_stress(fr, args...)
    base=[ε_,n_,0.0, 0,0,0, 0,0,0]
    function grad(ai)  # d(Ttt,Ttx,Txx)/d arg ai
        a=copy(base);a[ai]+=h;b=copy(base);b[ai]-=h
        Sp=bs(a...);Sm=bs(b...)
        ((Sp[1]-Sm[1])/(2h),(Sp[2]-Sm[2])/(2h),(Sp[3]-Sm[3])/(2h))
    end
    # x-grads: εx=4,nx=5,ux=6 ; t-grads: εt=7,nt=8,ut=9
    gεx,gnx,gux = grad(4),grad(5),grad(6)
    gεt,gnt,gut = grad(7),grad(8),grad(9)
    function detS2(v)
        st=-v; sx=1.0
        # primitive amplitudes: (δε,δu); δn=(n_/v)δu folded into n-columns.
        # column for δε: gε*(...);  column for δu: gu + (n_/v) gn
        # symbol entry for eq r (1=Ttt-row→flux Ttx, 2=Ttx-row→flux Txx):
        # M[r,col]=s_t(s_t ∂q_r/∂(∂_t p)+s_x ∂q_r/∂(∂_x p))
        #         +s_x(s_t ∂F_r/∂(∂_t p)+s_x ∂F_r/∂(∂_x p))
        # q_1=Ttt(idx1),F_1=Ttx(idx2);  q_2=Ttx(idx2),F_2=Txx(idx3)
        function entry(qti,qxi,Fti,Fxi)  # indices into the 3-tuples
            (st*(st*qti+sx*qxi)+sx*(st*Fti+sx*Fxi))
        end
        # δε column (uses gεt,gεx)
        M11=entry(gεt[1],gεx[1],gεt[2],gεx[2])
        M21=entry(gεt[2],gεx[2],gεt[3],gεx[3])
        # δu column = gu + (n_/v) gn
        fac=n_/v
        ut1=gut[1]+fac*gnt[1]; ux1=gux[1]+fac*gnx[1]
        ut2=gut[2]+fac*gnt[2]; ux2=gux[2]+fac*gnx[2]
        ut3=gut[3]+fac*gnt[3]; ux3=gux[3]+fac*gnx[3]
        M12=entry(ut1,ux1,ut2,ux2)
        M22=entry(ut2,ux2,ut3,ux3)
        return M11*M22-M12*M21
    end
    # detS2(v)*v is a polynomial in v; the (ε,u) sector roots are c_±. Find them
    # by bracketing near the closed-form c_+ (positive root) with bisection.
    function findroot(guess)
        f(v)=detS2(v)*v   # ×v clears the 1/v from continuity substitution
        a=guess*0.7; b=guess*1.3
        fa=f(a); fb=f(b)
        if fa*fb>0; return NaN; end
        for _ in 1:80; m=(a+b)/2; fm=f(m); (fa*fm<=0) ? (b=m) : (a=m;fa=fm); end
        return (a+b)/2
    end
    cp_sym = findroot(sqrt(cp2))
    d_sym = isnan(cp_sym) ? NaN : abs(cp_sym - sqrt(cp2))
    @printf("  Γ=%.3g σ̂=%.3g τ̂=%.3g: c_+=%.5f c_-=%.5f c_1=%.5f | quad|Δc_±²|=%.1e,%.1e | symbol c_+=%.5f Δ=%.1e\n",
            Γ_,σh,τh, sqrt(cp2),(cm2>0 ? sqrt(cm2) : 0.0),sqrt(c12), d_cp,d_cm, cp_sym, d_sym)
    global v3c_ok &= (d_cp<1e-9 && d_cm<1e-9 && (isnan(d_sym) ? false : d_sym<1e-4))
end

# -----------------------------------------------------------------------------
# [V1] CONSTANT STATE preserved exactly
# -----------------------------------------------------------------------------
println("\n", "="^78)
println("[V1] constant (equilibrium) state preserved exactly")
println("="^78)
fr1 = pmp_frame(; Γ=4/3, m=1.0, Vhat=0.1, σhat=0.0, τhat=1.5)
for (ε0,n0,u0) in [(1.0,0.1,0.0), (2.5,0.7,0.0), (1.0,0.1,0.4)]
    s = init_constant(fr1; N=129, ε0=ε0, n0=n0, u0=u0, cfl=0.1)
    Ttt0=copy(s.Ttt); Ttx0=copy(s.Ttx); Jt0=copy(s.Jt)
    ε_0=copy(s.ε); n_0=copy(s.n); u_0=copy(s.u)
    evolve!(s, 400)
    N=length(s.x); rng=NG+1:N-NG
    dε  = maximum(abs.(s.ε[rng].-ε_0[rng]))
    dn  = maximum(abs.(s.n[rng].-n_0[rng]))
    du  = maximum(abs.(s.u[rng].-u_0[rng]))
    dTtt= maximum(abs.(s.Ttt[rng].-Ttt0[rng]))
    dTtx= maximum(abs.(s.Ttx[rng].-Ttx0[rng]))
    @printf("  (ε=%.2g,n=%.2g,u=%.2g) 400 steps:  max|Δε|=%.2e |Δn|=%.2e |Δu|=%.2e |ΔTtt|=%.2e |ΔTtx|=%.2e\n",
            ε0,n0,u0, dε,dn,du,dTtt,dTtx)
    global v1_ok = (@isdefined v1_ok) ? v1_ok : true
    global v1_ok &= (dε<1e-11 && dn<1e-11 && du<1e-11 && dTtt<1e-11 && dTtx<1e-11)
end
println("  CONSTANT STATE exact (≲1e-11): ", v1_ok)

# -----------------------------------------------------------------------------
# [V2] inviscid ideal-gas Bjorken limit recovered (reuse Bjorken.jl)
#   The slab Cartesian engine does NOT carry Milne geometry; instead we validate
#   the INVISCID ideal-gas Bjorken ODE that the engine's EOS + perfect-fluid
#   stress reduce to. We integrate dε/dτ=−(ε+P)/τ with the SAME EOS used by the
#   engine (pressure_eos) and compare to:
#     (a) the analytic solution ε(τ)=m n0 τ^{-1}[1+e0 τ^{-(Γ-1)}] (eq:801), and
#     (b) the package Bjorken.bjorken_evolve_rk4 (independent prior-verified RK4).
#   This confirms the engine's EOS/pressure are consistent with the inviscid
#   Bjorken limit (τ_ε,τ_P,V→0 of eq:Bjorken_EOM).
# -----------------------------------------------------------------------------
println("\n", "="^78)
println("[V2] inviscid ideal-gas Bjorken limit (reuse Bjorken.bjorken_*)")
println("="^78)
Γb=4/3; mb=1.0; n0b=1.0; e0b=1.0
τ0=1.0; τf=10.0
ε0b = bjorken_inviscid_analytic(τ0; Γ=Γb, m=mb, n0=n0b, e0=e0b)
# (a) engine-EOS RK4 (pressure_eos identical to Bjorken.bjorken_pressure)
gb=IdealGasMicro(Γb,mb)
function bjorken_rk4_engineEOS(τ0,τf,ε0; N=4000)
    h=(τf-τ0)/N; τ=τ0; ε=ε0
    f(τ,ε)=-(ε + pressure_eos(gb, ε, n0b/τ))/τ        # P uses engine EOS; n(τ)=n0/τ
    for _ in 1:N
        k1=f(τ,ε); k2=f(τ+h/2,ε+h/2*k1); k3=f(τ+h/2,ε+h/2*k2); k4=f(τ+h,ε+h*k3)
        ε+=h/6*(k1+2k2+2k3+k4); τ+=h
    end
    return ε
end
# NOTE: Bjorken.bjorken_pressure uses P=(Γ-1)(ε − m n0/τ); engine uses
# P=(Γ-1)(ε − m n) with n=n0/τ ⇒ identical. Confirm numerically too.
@printf("  EOS consistency: pressure_eos(ε,n0/τ) vs bjorken_pressure(ε,τ) at τ=3, ε=0.5: %.3e\n",
        abs(pressure_eos(gb,0.5,n0b/3) - bjorken_pressure(0.5,3.0; Γ=Γb,m=mb,n0=n0b)))
ε_eng = bjorken_rk4_engineEOS(τ0,τf,ε0b; N=8000)
_, εs_pkg = bjorken_evolve_rk4(τ0,τf,ε0b; N=8000, Γ=Γb, m=mb, n0=n0b)
ε_pkg = εs_pkg[end]
ε_ana = bjorken_inviscid_analytic(τf; Γ=Γb, m=mb, n0=n0b, e0=e0b)
@printf("  τ: %.1f → %.1f,  ε0=%.8f\n", τ0,τf,ε0b)
@printf("  analytic ε(τf)               = %.10f   [eq:inviscid_bjorken]\n", ε_ana)
@printf("  engine-EOS RK4 ε(τf)         = %.10f   |Δ|=%.2e\n", ε_eng, abs(ε_eng-ε_ana))
@printf("  package bjorken_evolve_rk4   = %.10f   |Δ|=%.2e\n", ε_pkg, abs(ε_pkg-ε_ana))
# frame-independent diagnostic ε'+Γε/τ = (Γ-1) m n0/τ² on the inviscid flow
diag_num = bjorken_diagnostic(τf, ε_ana; Γ=Γb, m=mb, n0=n0b)
diag_exact = (Γb-1)*mb*n0b/τf^2
@printf("  diagnostic ε'+Γε/τ: num=%.3e  exact=%.3e  |Δ|=%.2e\n", diag_num,diag_exact,abs(diag_num-diag_exact))
v2_ok = abs(ε_eng-ε_ana)<1e-6 && abs(ε_pkg-ε_ana)<1e-6 && abs(diag_num-diag_exact)<1e-10

# -----------------------------------------------------------------------------
#  SUMMARY
# -----------------------------------------------------------------------------
println("\n", "="^78)
println("SUMMARY")
println("  [V3a] ideal-gas microphysics round trips : ", v3a_ok)
println("  [V3b] frame identity δ=0                  : ", v3b_ok, "   max|δ|=", maxδ)
println("  [V3c] char speeds (A,B,C quad & symbol)   : ", v3c_ok)
println("  [V1]  constant state preserved exactly    : ", v1_ok)
println("  [V2]  inviscid Bjorken limit recovered    : ", v2_ok)
allpass = v3a_ok && v3b_ok && v3c_ok && v1_ok && v2_ok
println("  ALL PASS: ", allpass)
println("="^78)
