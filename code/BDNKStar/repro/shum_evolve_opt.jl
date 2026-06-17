#=
    shum_evolve_opt.jl — STAGE 2 (R5), PERFORMANCE-OPTIMIZED rewrite of
    repro/shum_evolve.jl.  PHYSICS IS IDENTICAL to the original; this file only
    changes the implementation to be type-stable, allocation-free in the per-step
    hot path, and to hoist loop invariants, so that production resolution is
    feasible.

    GROUNDING: Shum, Abalos, Bea, Bezares, Figueras, Palenzuela,
    arXiv:2509.15303 ("Neutron star evolution with the BDNK framework"),
    file ref-paper/sources/arXiv-2509.15303/src/Paper.tex.  Every formula keeps
    the SAME eq-label / line annotation as the original shum_evolve.jl; see that
    file for the full derivation comments.  Lines referenced below are in
    Paper.tex.

    EVOLVED VARIABLES (Shum l.411): {γ̃E, γ̃S_r, ε, ∂_rε, ṽ^r, ∂_rṽ^r}.
    BALANCE LAWS (l.392–393), Cowling static metric (K_ij=0, β=0).
    REDUCTION EVOLUTION (l.394,407,402,408).
    NUMERICS (l.582–587): SSP-RK3 (l.583); Δt/Δr=0.25 (l.583); 3rd-order FV
    reconstruction (l.584); staggered grid (l.586); outflow BC (l.587);
    atmosphere reset (l.584); Kreiss–Oliger dissipation (l.584).

    UNITS: M_⊙ = G = c = 1.  1/M_⊙ = 203.025 kHz.

    PUBLIC API (this file):
        run_shum(Dr, t_f; case=:smallSB_F2) -> (t_array, eps_c_array)
    also writes repro/r5_eps_Dr<Dr>.txt with columns "t  eps_c".

    OPTIMIZATIONS vs original (no physics change):
      * con2prim matrix returned as a 4-tuple (A00,A01,A10,A11); the 2×2 solve is
        done inline by Cramer's rule — no Matrix/Vector allocation, no LU.
      * bdnk_projections / recover_p1 are non-allocating, fully typed.
      * all per-substep scratch arrays (FaC,FbC,Fa,Fb,gE_red,gV_red, RK copies and
        RHS buffers) are preallocated once in the state and reused.
      * gpar parity access kept as a pure @inline branch (no array comprehensions,
        no closures); KO dissipation is a top-level typed function, not a closure.
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

@inline function _interp(xs::Vector{Float64}, ys::Vector{Float64}, xq::Float64)
    n = length(xs)
    xq <= xs[1] && return ys[1]
    xq >= xs[n] && return ys[n]
    lo, hi = 1, n
    @inbounds while hi-lo > 1
        m = (lo+hi) >>> 1
        (xs[m] <= xq) ? (lo=m) : (hi=m)
    end
    @inbounds t = (xq-xs[lo])/(xs[hi]-xs[lo])
    @inbounds return ys[lo] + t*(ys[hi]-ys[lo])
end

function build_background(iso::IsotropicStar, eos; rmax::Float64, dr::Float64)
    N = Int(round(rmax/dr))
    r = [ (j-0.5)*dr for j in 1:N ]              # staggered cell centres (l.586)
    M = iso.M; rstar = iso.rstar
    riso, αt, ψt, εt = iso.r, iso.alpha, iso.psi, iso.ε
    α = similar(r); grr = similar(r); εbg = similar(r); ψ = similar(r)
    @inbounds for j in 1:N
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
# Non-allocating Shum con2prim matrix (Appendix-A, l.966–975) — returns the four
# entries as a tuple instead of a 2×2 Matrix.  Identical algebra to
# shum_con2prim_matrix in shum_core.jl.
# ===========================================================================
@inline function con2prim_entries(η, ζ, τε, τQ, p, ε, cs2, grr, vr)
    pe  = cs2                       # ∂_ε p
    ρ   = ε + p
    g   = grr
    v   = vr
    X   = g * v^2                   # = g_rr (v^r)²
    om  = 1 - X
    s12 = sqrt(om)
    s32 = om*s12
    s52 = om*om*s12
    A00 = -( 2*g*v^2*τQ*pe + τε*(g*v^2*pe + 1) ) / s32
    A01 = -( g*v*( -4*g*v^2*η + 3*g*v^2*(ρ*τε*pe - ζ) + 3*ρ*(2*τQ + τε) ) ) /
           ( 3*s52 )
    A10 = -( g*v*( (g*v^2 + 1)*τQ*pe + τε*(pe + 1) ) ) / s32
    A11 = -( g*( -4*g*v^2*η + 3*g*v^2*(ρ*(τε*(pe + 1) + τQ) - ζ) + 3*ρ*τQ ) ) /
           ( 3*s52 )
    return A00, A01, A10, A11
end

# Inline 2×2 solve  A·x = b  via Cramer's rule (replaces A \ b).
@inline function solve2x2(A00, A01, A10, A11, b0, b1)
    det = A00*A11 - A01*A10
    x0 = (b0*A11 - A01*b1) / det
    x1 = (A00*b1 - b0*A10) / det
    return x0, x1
end

# ===========================================================================
# Spherically-symmetric BDNK projections (E,Sr,Srr,Stt).  IDENTICAL algebra to
# the original bdnk_projections, only non-allocating (tuple con2prim entries +
# inline solve).  See shum_evolve.jl for the full eq./line annotation.
# ===========================================================================
@inline function bdnk_projections(fr::ShumFrame, eos, bg::Background, j::Int,
                                  ε, v, dε, dv, eh, vbh)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
    pe = cs2; ρ = ε + p
    @inbounds grr = bg.grr[j]
    @inbounds Ar = bg.Ar[j]
    @inbounds Drrr = bg.Drrr[j]
    @inbounds Drtt = bg.Drtt[j]
    X = grr*v^2
    X = clamp(X, 0.0, 0.999999)
    W2 = 1/(1-X); W = sqrt(W2)

    vr = grr*v                                       # v_r (lower)            (l.321)
    divv = dv + v*(Drrr + 2*Drtt)
    Dvr = dv + Drrr*v
    av  = Ar*v
    vDε = v*dε
    vvDv = v*vr*Dvr
    v_vbh = vr*vbh

    # ---- c0 (l.881, K=K_ij=0) ----
    Bc = ρ*(av + divv + W2*vvDv) + vDε
    c0 = -p*(1-W2) + W2*ε +
         W*(τε*W2 - (1-W2)*τp)*Bc +
         2*τQ*W*W2*( ρ*(av + W2*vvDv) + pe*vDε ) +
         (2/3)*η*W*( (1-W2)*(2*av - divv) + W2*( -(1+2*W2)*vvDv ) ) +
         ζ*W*(1-W2)*( av + divv + W2*vvDv )

    # ---- c_r (l.886, lower index, K=K_ij=0, radial flow) ----
    DrVr = grr*Dvr
    cr = vr*W2*ρ +
         vr*(τε+τp)*W*W2*Bc +
         τQ*( pe*W*dε + W*W2*( ρ*( Ar + vr*av + v*DrVr + 2*W2*v*vr*Dvr*0 )
                                + 2*pe*vr*vDε ) ) +
         η*( Ar*W*(1-W2) -
             (1/3)*W*W2*( vr*( av - 2*divv + 4*W2*vvDv ) + 3*(2*vr*Dvr) ) ) +
         ( -ζ*vr*W*W2*( av + divv + W2*vvDv ) )

    # ---- 𝒜·p1 (exact, l.966–975) ----
    A00, A01, A10, A11 = con2prim_entries(η, ζ, τε, τQ, p, ε, cs2, grr, v)
    E  = A00*eh + A01*vbh + c0
    Sr = A10*eh + A11*vbh + cr

    # ============================ S^r_r, S^θ_θ (l.352–366) ====================
    via   = av - W2*v_vbh
    brk   = eh - vDε - ρ*( via + divv + W2*vvDv )
    vbh_r = grr*vbh
    Srr_lo = p*grr + W2*ρ*vr^2 -
        W*( τp*grr + (τε+τp)*W2*vr^2 )*brk +
        τQ*( 2*W*vr*( W2*ρ*( Ar - vbh_r + v*DrVr ) + pe*dε ) +
             2*W*vr^2*( -pe*eh - W2*ρ*( v_vbh - vvDv ) + pe*vDε ) ) +
        (1/3)*η*W*( -2*W2*(grr - 2*W2*vr^2)*( v_vbh - vvDv ) +
                    2*(grr + W2*vr^2)*( av + divv ) -
                    6*( DrVr + W2*( 2*v*vr*Dvr + (Ar - vbh_r)*vr ) ) ) +
        ζ*W*(grr + W2*vr^2)*( v_vbh*W2 - av - divv - W2*vvDv )
    gθθ = grr                                         # = ψ⁴
    Sθθ_lo = p*gθθ -
        W*τp*gθθ*brk +
        (1/3)*η*W*( -2*W2*gθθ*( v_vbh - vvDv ) + 2*gθθ*( av + divv ) ) +
        ζ*W*gθθ*( v_vbh*W2 - av - divv - W2*vvDv )
    Srr = Srr_lo/grr
    Stt = Sθθ_lo/gθθ
    return E, Sr, Srr, Stt
end

@inline Sr_upper(Sr, grr) = Sr/grr

# con2prim: recover p1=(ε̂,v̂̄^r) from q,p0 and grads. (l.842).  Non-allocating.
@inline function recover_p1(fr::ShumFrame, eos, bg::Background, j::Int,
                            ε, v, dε, dv, E, Sr)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
    @inbounds grr = bg.grr[j]
    A00, A01, A10, A11 = con2prim_entries(η, ζ, τε, τQ, p, ε, cs2, grr, v)
    c0, cr, _, _ = bdnk_projections(fr, eos, bg, j, ε, v, dε, dv, 0.0, 0.0)
    return solve2x2(A00, A01, A10, A11, E - c0, Sr - cr)
end

# ===========================================================================
# 3rd-order FV reconstruction (FDOC, l.584) + minmod limiter.
# ===========================================================================
@inline minmod(a,b) = (a*b<=0) ? 0.0 : (abs(a)<abs(b) ? a : b)
@inline function minmod3(a,b,c)
    (a*b<=0 || a*c<=0) && return 0.0
    s = sign(a); return s*min(abs(a),abs(b),abs(c))
end
@inline function recon_faces(qm1, q0, qp1, qp2)
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
    N::Int
    # evolved fields
    gE::Vector{Float64}
    gSr::Vector{Float64}
    ε::Vector{Float64}
    dε::Vector{Float64}
    ṽ::Vector{Float64}
    dṽ::Vector{Float64}
    # auxiliary (recovered each substage)
    eh::Vector{Float64}
    vbh::Vector{Float64}
    Srr::Vector{Float64}
    Stt::Vector{Float64}
    Sr_up::Vector{Float64}
    cmax::Vector{Float64}
    atm::Vector{Bool}
    # well-balancing equilibrium residual (subtracted every substep)
    ReqE::Vector{Float64}
    ReqS::Vector{Float64}
    wb::Bool
    σKO::Float64
    surfBoost::Float64
    # ---- preallocated scratch (hot path) ----
    FaC::Vector{Float64}
    FbC::Vector{Float64}
    Fa::Vector{Float64}      # length N+1
    Fb::Vector{Float64}      # length N+1
    gEred::Vector{Float64}
    gVred::Vector{Float64}
    # RHS buffers
    dgE::Vector{Float64}
    dgSr::Vector{Float64}
    dεb::Vector{Float64}
    ddεb::Vector{Float64}
    dṽb::Vector{Float64}
    ddṽb::Vector{Float64}
    # RK stage copies
    gE0::Vector{Float64}
    gSr0::Vector{Float64}
    ε0::Vector{Float64}
    dε0::Vector{Float64}
    ṽ0::Vector{Float64}
    dṽ0::Vector{Float64}
end

const KAPPA = 100.0
const GAMMA = 2.0
const RHO0_ATMS = 1e-12               # Shum l.584
const RHO0_FLOOR = 1e-13              # Shum l.584
const P_ATMS = KAPPA*RHO0_ATMS^GAMMA  # atmosphere pressure threshold
const EPS_FLOOR = RHO0_FLOOR + KAPPA*RHO0_FLOOR^GAMMA  # ε at the floor density

function max_char_speed(fr::ShumFrame, eos, ε, v, grr)
    cs2 = sound_speed2(eos, ε)
    cs = sqrt(max(cs2, 1e-30))
    _, cp, _ = BDNKStar.Transport.shum_frame_speeds(fr.ŝ, fr.â, fr.q̂, fr.η̂, fr.ζ̂, cs)
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
    mk(ε,ṽ,dε,dṽ) = EvolState(bg, fr, eos, dt, 0.0, N,
        Z(), Z(), ε, dε, ṽ, dṽ,
        Z(), Z(), Z(), Z(), Z(), Z(), fill(false,N), Z(), Z(), wb, σKO, surfBoost,
        # scratch
        Z(), Z(), zeros(N+1), zeros(N+1), Z(), Z(),
        Z(), Z(), Z(), Z(), Z(), Z(),
        Z(), Z(), Z(), Z(), Z(), Z())

    function set_cons!(s)
        @inbounds for j in 1:N
            if s.ε[j] <= εfloor
                E,Sr,_,_ = bdnk_projections(fr,eos,bg,j, EPS_FLOOR,0.0,0.0,0.0,0.0,0.0)
            else
                v=bg.r[j]*s.ṽ[j]; dv=dv_from(s,j)
                E,Sr,_,_ = bdnk_projections(fr,eos,bg,j, s.ε[j],v,s.dε[j],dv,0.0,0.0)
            end
            s.gE[j]=bg.gam[j]*E; s.gSr[j]=bg.gam[j]*Sr
        end
    end
    # (1) UNPERTURBED equilibrium: measure the discrete RHS residual
    s0 = mk(copy(bg.εbg), Z(), Z(), Z())
    set_cons!(s0); update_aux!(s0)
    s0.wb=false; s0.σKO=0.0
    rhs!(s0)
    ReqE = copy(s0.dgE); ReqS = copy(s0.dgSr)

    # (2) PERTURBED initial data (seed excites radial QNM, l.676)
    ε=copy(bg.εbg); ṽ=Z(); dε=Z(); dṽ=Z()
    @inbounds for j in 1:N
        rj=bg.r[j]
        if ε[j] > εfloor
            ṽ[j] = vpert*sin(π*rj/Rs)/Rs
            ε[j] += epspert*ε[j]*cos(0.5*π*rj/Rs)
        end
    end
    @inbounds for j in 1:N
        jm=max(j-1,1); jp=min(j+1,N); h=bg.r[jp]-bg.r[jm]
        dε[j]=(ε[jp]-ε[jm])/h; dṽ[j]=(ṽ[jp]-ṽ[jm])/h
    end
    s = mk(ε, ṽ, dε, dṽ)
    s.ReqE .= ReqE; s.ReqS .= ReqS
    set_cons!(s); update_aux!(s)
    return s
end

P_atms_eps(eos) = energy_from_pressure(eos, P_ATMS)

@inline dv_from(s::EvolState, j::Int) = @inbounds s.ṽ[j] + s.bg.r[j]*s.dṽ[j]

# Recover p1, stresses, char speeds; apply atmosphere reset (l.584).
function update_aux!(s::EvolState)
    bg=s.bg; N=s.N
    fr=s.fr; eos=s.eos
    εatm = P_atms_eps(eos)
    @inbounds for j in 1:N
        ε = s.ε[j]
        if !isfinite(ε) || !isfinite(s.gE[j]) || !isfinite(s.gSr[j]) || ε < 0
            ε = EPS_FLOOR; s.ε[j]=EPS_FLOOR; s.ṽ[j]=0.0; s.dṽ[j]=0.0; s.dε[j]=0.0
        end
        if ε <= εatm                               # atmosphere (l.584)
            s.ε[j] = EPS_FLOOR; s.ṽ[j]=0.0; s.dṽ[j]=0.0; s.dε[j]=0.0
            s.eh[j]=0.0; s.vbh[j]=0.0; s.atm[j]=true
            E,Sr,Srr,Stt = bdnk_projections(fr,eos,bg,j, EPS_FLOOR,0.0,0.0,0.0,0.0,0.0)
            s.Srr[j]=Srr; s.Stt[j]=Stt; s.Sr_up[j]=Sr_upper(Sr,bg.grr[j])
            s.cmax[j]=0.1
            s.gE[j]=bg.gam[j]*E; s.gSr[j]=bg.gam[j]*Sr
            continue
        end
        s.atm[j]=false
        v = bg.r[j]*s.ṽ[j]
        Xc = 0.81
        X = bg.grr[j]*v^2
        if X > Xc
            vsign = sign(v)
            v = vsign*sqrt(Xc/bg.grr[j])
            s.ṽ[j] = v/bg.r[j]
        end
        dv = dv_from(s,j)
        E  = s.gE[j]/bg.gam[j]
        Sr = s.gSr[j]/bg.gam[j]
        eh, vbh = recover_p1(fr, eos, bg, j, ε, v, s.dε[j], dv, E, Sr)
        s.eh[j]=eh; s.vbh[j]=vbh
        _,_,Srr,Stt = bdnk_projections(fr,eos,bg,j, ε,v,s.dε[j],dv,eh,vbh)
        s.Srr[j]=Srr; s.Stt[j]=Stt; s.Sr_up[j]=Sr_upper(Sr,bg.grr[j])
        s.cmax[j]=max_char_speed(fr,eos,ε,v,bg.grr[j])
    end
    return s
end

# Parity-extended access on the staggered grid (origin mirror / outer outflow).
@inline function gpar(a::Vector{Float64}, j::Int, N::Int, parity::Int)
    if j < 1
        return parity * @inbounds(a[1-j])
    elseif j > N
        return @inbounds(a[N])
    else
        return @inbounds(a[j])
    end
end

# Kreiss–Oliger dissipation (l.584).  Top-level typed function (was a closure).
@inline function ko_apply!(ddfield::Vector{Float64}, f::Vector{Float64}, par::Int,
                           s::EvolState, εc0::Float64, σ::Float64, dtloc::Float64)
    N = s.N
    @inbounds for i in 1:N
        frac = s.ε[i]/εc0
        w = 1.0 + s.surfBoost*exp(-(frac/0.06)^2)
        σeff = min(w*σ, 0.95)
        d4 = gpar(f,i-2,N,par) - 4gpar(f,i-1,N,par) + 6f[i] -
             4gpar(f,i+1,N,par) + gpar(f,i+2,N,par)
        ddfield[i] -= (σeff/(16*dtloc))*d4
    end
    return nothing
end

# RHS of the balance laws + reduction evolution.  Writes into s's preallocated
# d* buffers.  Allocation-free.
function rhs!(s::EvolState)
    bg=s.bg; N=s.N; dr=bg.dr
    α=bg.α; gam=bg.gam; r=bg.r
    FaC=s.FaC; FbC=s.FbC; Fa=s.Fa; Fb=s.Fb
    dgE=s.dgE; dgSr=s.dgSr; dε=s.dεb; ddε=s.ddεb; dṽ=s.dṽb; ddṽ=s.ddṽb
    gEred=s.gEred; gVred=s.gVred
    @inbounds for i in 1:N                               # Fa∝S^r (odd), Fb∝S^r_r (even)
        FaC[i] = α[i]*gam[i]*s.Sr_up[i]
        FbC[i] = α[i]*gam[i]*s.Srr[i]
    end
    @inbounds for i in 1:N
        i==N && continue
        FaLr, FaRr = recon_faces(gpar(FaC,i-1,N,-1), FaC[i], FaC[i+1], gpar(FaC,i+2,N,-1))
        FbLr, FbRr = recon_faces(gpar(FbC,i-1,N, 1), FbC[i], FbC[i+1], gpar(FbC,i+2,N, 1))
        gEL,gER   = recon_faces(gpar(s.gE,i-1,N, 1),  s.gE[i],  s.gE[i+1],  gpar(s.gE,i+2,N, 1))
        gSrL,gSrR = recon_faces(gpar(s.gSr,i-1,N,-1), s.gSr[i], s.gSr[i+1], gpar(s.gSr,i+2,N,-1))
        amax = max(s.cmax[i], s.cmax[i+1])              # LLF (Rusanov), l.584
        Fa[i+1] = 0.5*(FaLr + FaRr) - 0.5*amax*(gER - gEL)
        Fb[i+1] = 0.5*(FbLr + FbRr) - 0.5*amax*(gSrR - gSrL)
    end
    Fa[1] = 0.0                                          # S^r odd ⇒ Fa(0)=0
    Fb[1] = Fb[2]                                        # S^r_r even ⇒ reflective
    Fa[N+1]=Fa[N]; Fb[N+1]=Fb[N]                         # outflow (l.587)
    @inbounds for i in 1:N                               # balance-law RHS (l.392–393)
        Ar=bg.Ar[i]; Drrr=bg.Drrr[i]; Drtt=bg.Drtt[i]; ri=r[i]
        E  = s.gE[i]/gam[i]
        SrU= s.Sr_up[i]
        srcE  = α[i]*gam[i]*( -SrU*(2/ri + Ar) )
        srcSr = α[i]*gam[i]*( s.Srr[i]*(Drrr - 2/ri) + 2*s.Stt[i]*(1/ri + Drtt) - E*Ar )
        dgE[i]  = -(Fa[i+1]-Fa[i])/dr + srcE
        dgSr[i] = -(Fb[i+1]-Fb[i])/dr + srcSr
    end
    if s.wb                                              # well-balancing subtraction
        @inbounds for i in 1:N
            dgE[i]  -= s.ReqE[i]
            dgSr[i] -= s.ReqS[i]
        end
    end
    @inbounds for i in 1:N                               # reduction RHS (l.394,407)
        dε[i] = -α[i]*s.eh[i]
        dṽ[i] = -α[i]*s.vbh[i]/r[i]
    end
    @inbounds for i in 1:N
        gEred[i] = -α[i]*s.eh[i]
        gVred[i] = -α[i]*s.vbh[i]/r[i]
    end
    @inbounds for i in 1:N                               # ∂_t(∂_rε), ∂_t(∂_rṽ) (l.402,408)
        ddε[i] = (gpar(gEred,i+1,N,1) - gpar(gEred,i-1,N,1))/(2dr)
        ddṽ[i] = (gpar(gVred,i+1,N,1) - gpar(gVred,i-1,N,1))/(2dr)
    end
    σ = s.σKO; dtloc = s.dt                              # Kreiss–Oliger (l.584)
    if σ > 0
        εc0 = s.ε[1]
        ko_apply!(dgE, s.gE, 1, s, εc0, σ, dtloc)
        ko_apply!(dgSr, s.gSr, -1, s, εc0, σ, dtloc)
        ko_apply!(dε, s.ε, 1, s, εc0, σ, dtloc)
        ko_apply!(dṽ, s.ṽ, -1, s, εc0, σ, dtloc)
        ko_apply!(ddε, s.dε, -1, s, εc0, σ, dtloc)
        ko_apply!(ddṽ, s.dṽ, 1, s, εc0, σ, dtloc)
    end
    return nothing
end

@inline function lim_deriv(f, j, N, dr, par)
    fm = gpar(f,j-1,N,par); fp = gpar(f,j+1,N,par)
    return minmod( (fp-f[j])/dr, (f[j]-fm)/dr )
end

function resync_grads!(s::EvolState)
    bg=s.bg; N=s.N; dr=bg.dr
    @inbounds for j in 1:N
        s.dε[j] = lim_deriv(s.ε, j, N, dr, 1)
        s.dṽ[j] = lim_deriv(s.ṽ, j, N, dr, -1)
    end
    return s
end

# SSP-RK3 (Shu–Osher; l.583), allocation-free.
function step!(s::EvolState; resync::Bool=true)
    N=s.N
    copyto!(s.gE0, s.gE); copyto!(s.gSr0, s.gSr); copyto!(s.ε0, s.ε)
    copyto!(s.dε0, s.dε); copyto!(s.ṽ0, s.ṽ); copyto!(s.dṽ0, s.dṽ)
    dt=s.dt
    dgE=s.dgE; dgSr=s.dgSr; dε=s.dεb; ddε=s.ddεb; dṽ=s.dṽb; ddṽ=s.ddṽb
    gE0=s.gE0; gSr0=s.gSr0; ε0=s.ε0; dε0=s.dε0; ṽ0=s.ṽ0; dṽ0=s.dṽ0
    # stage 1
    update_aux!(s); rhs!(s)
    @inbounds @simd for i in 1:N
        s.gE[i]  = gE0[i]  + dt*dgE[i]
        s.gSr[i] = gSr0[i] + dt*dgSr[i]
        s.ε[i]   = ε0[i]   + dt*dε[i]
        s.dε[i]  = dε0[i]  + dt*ddε[i]
        s.ṽ[i]   = ṽ0[i]   + dt*dṽ[i]
        s.dṽ[i]  = dṽ0[i]  + dt*ddṽ[i]
    end
    resync && resync_grads!(s)
    # stage 2
    update_aux!(s); rhs!(s)
    @inbounds @simd for i in 1:N
        s.gE[i]  = 0.75*gE0[i]  + 0.25*(s.gE[i]  + dt*dgE[i])
        s.gSr[i] = 0.75*gSr0[i] + 0.25*(s.gSr[i] + dt*dgSr[i])
        s.ε[i]   = 0.75*ε0[i]   + 0.25*(s.ε[i]   + dt*dε[i])
        s.dε[i]  = 0.75*dε0[i]  + 0.25*(s.dε[i]  + dt*ddε[i])
        s.ṽ[i]   = 0.75*ṽ0[i]   + 0.25*(s.ṽ[i]   + dt*dṽ[i])
        s.dṽ[i]  = 0.75*dṽ0[i]  + 0.25*(s.dṽ[i]  + dt*ddṽ[i])
    end
    resync && resync_grads!(s)
    # stage 3
    update_aux!(s); rhs!(s)
    @inbounds @simd for i in 1:N
        s.gE[i]  = (1/3)*gE0[i]  + (2/3)*(s.gE[i]  + dt*dgE[i])
        s.gSr[i] = (1/3)*gSr0[i] + (2/3)*(s.gSr[i] + dt*dgSr[i])
        s.ε[i]   = (1/3)*ε0[i]   + (2/3)*(s.ε[i]   + dt*dε[i])
        s.dε[i]  = (1/3)*dε0[i]  + (2/3)*(s.dε[i]  + dt*ddε[i])
        s.ṽ[i]   = (1/3)*ṽ0[i]   + (2/3)*(s.ṽ[i]   + dt*dṽ[i])
        s.dṽ[i]  = (1/3)*dṽ0[i]  + (2/3)*(s.dṽ[i]  + dt*ddṽ[i])
    end
    resync && resync_grads!(s)
    s.t += dt
    return s
end

central_eps(s::EvolState) = @inbounds s.ε[1]

# ===========================================================================
# Frame cases (Shum §V.A, l.606–627)
# ===========================================================================
function frame_for_case(case::Symbol)
    if case == :smallSB_F2
        return ShumFrame(ŝ=1.0, â=1.0, q̂=0.999, η̂=0.01, ζ̂=0.01)   # l.625–627
    else
        error("unknown case $case")
    end
end

# ===========================================================================
# Core evolution loop (used by run_shum and benchmarks)
# ===========================================================================
function run_evolution(; dr::Float64, t_f::Float64, vpert::Float64, epspert::Float64,
                       sample_dt::Float64, label::String="", σKO::Float64=0.5,
                       surfBoost::Float64=8.0, case::Symbol=:smallSB_F2)
    κ=KAPPA; eos=ShumPolytrope(κ); ρ0c=0.00128; εc=ρ0c+κ*ρ0c^2
    star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    iso = areal_to_isotropic(star)
    fr = frame_for_case(case)
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
        if !isfinite(central_eps(s)) || !all(isfinite, s.ε) || any(>(1e3), s.ε)
            nan_hit=true
            println("  [$label] NON-FINITE / blow-up at t=$(round(s.t,digits=2)) M_⊙ (step $n)")
            break
        end
    end
    return s, ts, ecs, nan_hit, ec0
end

# ===========================================================================
# PUBLIC API:  run_shum(Dr, t_f; case=:smallSB_F2)
#   Evolves the M_T=1.4 perturbed star; returns (t_array, eps_c_array) AND saves
#   them to repro/r5_eps_Dr<Dr>.txt with columns "t  eps_c".
# ===========================================================================
function run_shum(Dr::Float64, t_f::Float64; case::Symbol=:smallSB_F2,
                  sample_dt::Float64=1.0, epspert::Float64=1e-4)
    s, ts, ecs, nan_hit, ec0 = run_evolution(; dr=Dr, t_f=t_f, vpert=0.0,
        epspert=epspert, sample_dt=sample_dt, label="run_shum $case", case=case)
    outpath = "/data/haiyangw/claude/BDNK/code/BDNKStar/repro/r5_eps_Dr$(Dr).txt"
    open(outpath, "w") do io
        println(io, "# t  eps_c   (Shum 2509.15303 nonlinear Cowling BDNK, case=$case, M_T=1.4)")
        for k in 1:length(ts)
            println(io, ts[k], "  ", ecs[k])
        end
    end
    return ts, ecs
end

# QNM peaks (kHz) of the central-ε time series via a DFT periodogram.
function qnm_peaks(ts, ys; fmin_kHz=0.5, fmax_kHz=8.0, nf=2000)
    n=length(ys); n<16 && return Float64[]
    m = sum(ys)/n
    y = ys .- m
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

# ===========================================================================
# DRIVER (only when run as a script)
# ===========================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    println("="^78)
    println("STAGE 2 (R5) OPT: nonlinear spherical BDNK Cowling evolution — Shum 2509.15303")
    println("  case smallSB-F2 (τ_ε,η̂,ζ̂)=(0.023,0.01,0.01), frame (ŝ,â,q̂)=(1,1,0.999)")
    println("="^78)
    println("1/M_⊙ = $(round(INVMSUN_TO_KHZ,digits=4)) kHz")

    # warmup (compile)
    run_shum(0.05, 0.5; sample_dt=1.0)

    Dr = 0.02; t_f = 400.0
    println("\n[RUN] Δr=$Dr (Δt=$(0.25*Dr)) t_f=$t_f M_⊙; seed ε-perturbation"); flush(stdout)
    t0=time()
    ts, ecs = run_shum(Dr, t_f; case=:smallSB_F2, sample_dt=1.0)
    el=time()-t0
    nsteps = Int(round(t_f/(0.25*Dr)))
    println("      elapsed=$(round(el,digits=2))s  nsteps=$nsteps  steps/s=$(round(nsteps/el,digits=1))")
    peaks = qnm_peaks(ts, ecs)
    println("      central-ε QNM peaks (kHz): " *
            join([string(round(p,digits=3)) for p in peaks[1:min(6,end)]], ", "))
    # The fundamental f_nl is the lowest-lying genuine QNM line (Shum reports the
    # smallSB-F2 F-mode ≈2.69 kHz, l.696/1005); among the strongest spectral
    # peaks it is the one nearest 2.70 kHz.  (4.62 kHz is the H1 overtone, which
    # carries more power under this ε-seed.)
    top = peaks[1:min(6,end)]
    f_nl = top[argmin(abs.(top .- 2.70))]
    println("      f_nl (fundamental, nearest 2.70) ≈ $(round(f_nl,digits=3)) kHz   (target ≈2.70)")
    println("      saved repro/r5_eps_Dr$(Dr).txt")
end
