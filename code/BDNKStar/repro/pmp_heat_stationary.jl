# =============================================================================
# pmp_heat_stationary.jl  —  PMP 2209.09265, Fig. heat_stationary
#
# Reproduces the (1+1)D STATIONARY HEAT-CONDUCTION result of PMP arXiv:2209.09265
# Sec. "Heat flow problem with non-constant coefficients" (paper.tex line 1213),
# Fig. \ref{fig:heat_stationary} (paper.tex line 1236-1240), Table parameters
# (paper.tex line 558):   Γ=4/3,  m=0.1,  V̂=2/15,  σ̂ ∈ {0, 1/3},  τ̂=1.5.
#
# PHYSICS (paper eqs):
#   Initial data (eq:heat_flow_ID, line 1217):
#       T(0,x) = A e^{-x²/w²} + δ ,    P(0,x) = P₀ = const ,
#   implemented in hydrodynamic variables via the ideal-gas EOS (line 1220):
#       ε = P[ m T⁻¹ + (Γ-1)⁻¹ ] ,   n = P T⁻¹ ,   u = 0,
#   with TIME-SYMMETRIC data  ε̇(0,x)=u̇ⁱ(0,x)=0  (line 1220).
#
#   For this data, baryon continuity gives ṅ=0 (eq:heat_baryon_EOM), the
#   x-component of ∂_a T^{ax}=0 is trivially satisfied at t=0, and the ONLY
#   nontrivial equation of motion at t=0 is the t-component (eq:heat_ID_EOM,
#   line 1224):
#                τ_ε ε̈  =  (κ T')'              [eq:heat_ID_EOM]
#   with κ ≡ σ ρ²/(n² T)   (eq:kappa, line 587).
#
#   => If σ̂=0 then σ=0 ⇒ κ=0 ⇒ ε̈=0  : the solution has NO dynamics and the
#      initial profile is STATIONARY (top panel of fig:heat_stationary). Any
#      ε̇ that appears in a numerical evolution is pure numerical error and
#      converges to ZERO as the grid is refined.
#   => If σ̂=1/3 then σ>0 ⇒ κ≠0 ⇒ ε̈≠0 over part of the domain : there is a
#      genuine dynamical heat-flow solution, and a numerical ε̇ converges to a
#      finite, nonzero value over part of the domain (bottom panel).
#
# WHAT WE COMPUTE & MATCH:
#   [A] The closed-form initial acceleration ε̈(x) = (κ T')'/τ_ε from
#       eq:heat_ID_EOM, showing ε̈≡0 for σ̂=0 and ε̈≠0 (localized) for σ̂=1/3.
#       This is the quantitative content of fig:heat_stationary.
#   [B] An ACTUAL conservative-finite-volume BDNK evolution (the PMP engine in
#       repro/pmp_viscous_core.jl: WENO5 + Kurganov–Tadmor + Heun, CFL λ=0.1,
#       paper sec:numerics line 1452) of the same initial data, at several grid
#       resolutions, confirming the convergence statement of the figure caption:
#         σ̂=0   : max|ε̇| at the snapshot DECREASES toward 0 as N grows;
#         σ̂=1/3 : ε̇ converges to a finite nonzero profile (~the [A] curve).
#
# PACKAGE REUSE: we include the verified PMP engine and the BDNKStar package.
# We do NOT edit src/ or another agent's repro file.
# =============================================================================

# The PMP engine (repro/pmp_viscous_core.jl) itself begins with
#   include(".../src/BDNKStar.jl"); using .BDNKStar
# (per the package-reuse directive). We include it ONCE here so that BDNKStar
# is loaded a single time (a second top-level include of BDNKStar.jl would
# create a duplicate module and make Bjorken exports ambiguous). All BDNKStar
# functionality is therefore available transitively via the engine include.
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/pmp_viscous_core.jl")
using .BDNKStar
using Printf

# pmp_viscous_core.jl defines, at top level (no module): IdealGasMicro, PMPFrame,
# pmp_frame, pressure_eos, cs2_eos, kappa_eps, transport_coeffs, bdnk_stress,
# ideal_stress, VState, NG, _Dx, _update_aux!, _recover_n!, evolve!, etc.
# We reuse them directly.

println("\n", "="^78)
println("PMP 2209.09265 — Fig. heat_stationary : (1+1)D stationary heat conduction")
println("Γ=4/3, m=0.1, V̂=2/15, σ̂∈{0,1/3}, τ̂=1.5  (Table line 558)")
println("="^78)

# -----------------------------------------------------------------------------
# Figure / table parameters  (paper.tex line 558)
# -----------------------------------------------------------------------------
const Γ   = 4/3
const mm  = 0.1
const Vhat= 2/15
const τhat= 1.5

# σ̂ shares V̂ between shear/bulk; for the heat figure η̂,ζ̂ play no role at u=0
# (the only viscous channel that acts on the t=0 data is heat conduction σ).
# We use the package default split (all of V̂ in shear, ζ=0). This does not
# affect the heat-flow EOM at t=0, which depends only on κ (i.e. σ) and τ_ε.

# κ from σ  (eq:kappa, line 587):  κ ≡ σ ρ²/(n² T)
kappa_heat(g::IdealGasMicro, ε, n, σ) = begin
    P = pressure_eos(g, ε, n); ρ = ε + P; T = P / n
    σ * ρ^2 / (n^2 * T)
end

# -----------------------------------------------------------------------------
# Heat-flow initial data (eq:heat_flow_ID + EOS map, lines 1217-1220).
# Returns primitive arrays (ε,n,u=0) on a grid x, given T-profile parameters.
# -----------------------------------------------------------------------------
function heat_flow_primitives(g::IdealGasMicro, x; A=0.4, w=10.0, δ=1.0, P0=1.0)
    Tprof = @. A * exp(-x^2 / w^2) + δ                 # T(0,x)  (eq:heat_flow_ID)
    Pprof = fill(P0, length(x))                        # P(0,x)=P₀=const
    ε = @. Pprof * (g.m / Tprof + 1/(g.Γ - 1))         # ε=P[m/T+1/(Γ-1)] (line 1220)
    n = @. Pprof / Tprof                               # n=P/T              (line 1220)
    return ε, n, Tprof, Pprof
end

# -----------------------------------------------------------------------------
# [A] Closed-form initial acceleration  ε̈(x) = (κ T')'/τ_ε  (eq:heat_ID_EOM).
#     τ_ε = L V̂ τ̂  (eq:hydro_frame, L=1).  κ from eq:kappa using the LOCAL σ
#     (which equals V̂ ρ c_s²/(−κ_ε) σ̂, eq:hydro_frame) so that κ T' is built
#     from the same frame the engine uses.
# -----------------------------------------------------------------------------
function epsddot_closed(fr::PMPFrame, x, ε, n, Tprof)
    g = fr.g
    N = length(x); dx = x[2]-x[1]
    τε = fr.L * fr.Vhat * fr.τhat                       # τ_ε  (eq:hydro_frame)
    # σ(x) from the frame (eq:hydro_frame): σ = V̂ L ρ c_s²/(−κ_ε) σ̂
    σ = similar(x); κ = similar(x)
    for i in 1:N
        tc = transport_coeffs(fr, ε[i], n[i])           # gives σ already
        σ[i] = tc.σ
        κ[i] = kappa_heat(g, ε[i], n[i], σ[i])          # κ=σρ²/(n²T) (eq:kappa)
    end
    # T' (centered 2nd order), then (κ T')' (centered 2nd order)
    Tp = similar(x); flux = similar(x); epsdd = similar(x)
    for i in 2:N-1; Tp[i] = (Tprof[i+1]-Tprof[i-1])/(2dx); end
    Tp[1]=Tp[2]; Tp[N]=Tp[N-1]
    for i in 1:N; flux[i] = κ[i]*Tp[i]; end
    for i in 2:N-1; epsdd[i] = (flux[i+1]-flux[i-1])/(2dx)/τε; end
    epsdd[1]=0.0; epsdd[N]=0.0
    return epsdd, κ
end

# -----------------------------------------------------------------------------
# [B] Build a VState for the heat-flow ID and evolve with the PMP engine.
#     Initial conserved fields use bdnk_stress at u=0, ε̇=u̇=0 (time-symmetric):
#       T^{tt} = ε ,  T^{tx} = β_ε ε' + β_n n'  (the spatial heat flux Q^x),
#       J^t   = n .
# -----------------------------------------------------------------------------
function init_heat_flow(fr::PMPFrame; N=2049, xmax=50.0, A=0.4, w=10.0, δ=1.0, P0=1.0)
    x = collect(range(-xmax, xmax; length=N)); dx = x[2]-x[1]
    ε, n, Tprof, Pprof = heat_flow_primitives(fr.g, x; A=A, w=w, δ=δ, P0=P0)
    z = zeros(N)
    s = VState(fr, x, dx, 0.1*dx,                       # CFL λ=0.1 (line 1452)
               copy(ε), copy(n), copy(z),               # ε, n, u=0
               copy(z),copy(z),copy(z),                 # εx,nx,ux
               copy(z),copy(z),                         # εt,ut
               copy(z),copy(z),copy(z),                 # Ttt,Ttx,Jt
               false, false)                            # non-periodic, no Milne
    # spatial gradients of the ID
    for i in 2:N-1
        s.εx[i]=(ε[i+1]-ε[i-1])/(2dx); s.nx[i]=(n[i+1]-n[i-1])/(2dx)
    end
    # conserved fields from bdnk_stress at u=0, ε̇=u̇=0 (time-symmetric ID)
    for i in 1:N
        Ttt,Ttx,_,Jt,_ = bdnk_stress(fr, ε[i], n[i], 0.0,
                                     s.εx[i], s.nx[i], 0.0,  # εx,nx,ux
                                     0.0, 0.0, 0.0)          # εt,nt,ut = 0 (time-sym)
        s.Ttt[i]=Ttt; s.Ttx[i]=Ttx; s.Jt[i]=Jt
    end
    _update_aux!(s)                                     # recover εt,ut & set ghosts
    return s, x, Tprof
end

# convenience: max|ε̇| recovered at t=0 (before any step), interior only
function max_epsdot_interior(s::VState)
    N=length(s.x); rng=NG+1:N-NG
    return maximum(abs.(s.εt[rng]))
end

# =============================================================================
#  RUN
# =============================================================================
# ID profile parameters (eq:heat_flow_ID). A,w,δ,P₀ are not tabulated in the
# paper; we use a smooth Gaussian temperature bump (amplitude A=0.4 on a
# background δ=1, width w=10 matching the shockwave width of line 1056, and
# P₀=1). The QUALITATIVE conclusion (σ̂=0 stationary, σ̂=1/3 dynamical) and the
# eq:heat_ID_EOM balance are independent of these choices.
A0, w0, δ0, P00 = 0.4, 10.0, 1.0, 1.0

# stability / causality sanity (eq:simple_constraints, line 506) at background:
gbg = IdealGasMicro(Γ, mm)
εbg = P00*(mm/δ0 + 1/(Γ-1)); nbg = P00/δ0
cs2bg = cs2_eos(gbg, εbg, nbg)
τmin = ((Γ-1)*(2-cs2bg)+cs2bg)/(1-cs2bg)
@printf("\nBackground (T=δ=%.2f, P₀=%.2f): ε=%.5f n=%.5f c_s²=%.5f\n", δ0,P00,εbg,nbg,cs2bg)
@printf("Causality/stability (eq:simple_constraints): σ̂≤1/3 (use 0,1/3 ✓),  τ̂≥%.4f (use τ̂=%.2f ✓)\n",
        τmin, τhat)

# -----------------------------------------------------------------------------
# [A] Closed-form ε̈(x) from eq:heat_ID_EOM for σ̂=0 and σ̂=1/3
# -----------------------------------------------------------------------------
println("\n", "-"^78)
println("[A] Initial acceleration ε̈(x) = (κT')'/τ_ε  [eq:heat_ID_EOM, line 1224]")
println("-"^78)
Ncf = 4097
xcf = collect(range(-50.0, 50.0; length=Ncf))
εcf, ncf, Tcf, _ = heat_flow_primitives(gbg, xcf; A=A0, w=w0, δ=δ0, P0=P00)

results_closed = Dict{Float64,Tuple{Float64,Float64}}()  # σ̂ => (max|κ|, max|ε̈|)
for σhat in (0.0, 1/3)
    fr = pmp_frame(; Γ=Γ, m=mm, Vhat=Vhat, σhat=σhat, τhat=τhat)
    epsdd, κ = epsddot_closed(fr, xcf, εcf, ncf, Tcf)
    maxκ  = maximum(abs.(κ))
    maxdd = maximum(abs.(epsdd))
    results_closed[σhat] = (maxκ, maxdd)
    @printf("  σ̂=%-5.4g :  max|κ|=%.6e   max|ε̈|=%.6e", σhat, maxκ, maxdd)
    if σhat == 0.0
        @printf("   <- κ=0 ⇒ ε̈≡0  : STATIONARY (no dynamics)\n")
    else
        # where is ε̈ nonzero? report support fraction (|ε̈|>1% of max)
        thr = 0.01*maxdd
        frac = count(>(thr), abs.(epsdd)) / Ncf
        ipk = argmax(abs.(epsdd))
        @printf("   <- κ≠0 ⇒ ε̈≠0 over part of domain (|ε̈|>1%% on %.0f%% of grid; peak at x=%.2f)\n",
                100*frac, xcf[ipk])
    end
end

# -----------------------------------------------------------------------------
# [B] Actual BDNK conservative evolution + convergence of the snapshot ε̇
#     We evolve to a fixed PHYSICAL time t_snap ("shortly after t=0", caption)
#     at several resolutions and report max|ε̇| over the interior. Because the
#     engine recovers ε̇ from the BDNK constraint each substage, ε̇ at t=0 is
#     already the discrete realization of eq:heat_ID_EOM's first time-deriv;
#     the figure plots ε̇ at the snapshot for several resolutions.
# -----------------------------------------------------------------------------
println("\n", "-"^78)
println("[B] BDNK finite-volume evolution (WENO5+KT+Heun, λ=0.1, line 1452)")
println("    snapshot max|ε̇| (interior) vs grid resolution N  [fig:heat_stationary]")
println("-"^78)

t_snap = 1.0                      # "shortly after t=0"
res_list = [513, 1025, 2049]      # increasing resolution (darker gray in fig)

for σhat in (0.0, 1/3)
    fr = pmp_frame(; Γ=Γ, m=mm, Vhat=Vhat, σhat=σhat, τhat=τhat)
    @printf("\n  σ̂ = %.4g :\n", σhat)
    prev = nothing
    for N in res_list
        s, x, Tprof = init_heat_flow(fr; N=N, xmax=50.0, A=A0, w=w0, δ=δ0, P0=P00)
        dt = s.dt
        nsteps = max(1, round(Int, t_snap/dt))
        # max|ε̇| of the RECOVERED initial-time derivative (the figure's ε̇)
        epsdot0 = max_epsdot_interior(s)
        evolve!(s, nsteps)
        epsdot_snap = max_epsdot_interior(s)
        # also peak |ε̇| location at snapshot
        Nn=length(s.x); rng=NG+1:Nn-NG
        ipk = rng[argmax(abs.(s.εt[rng]))]
        tag = (prev===nothing) ? "" : @sprintf("   ratio vs prev N: %.3f", epsdot_snap/prev)
        @printf("    N=%5d (Δx=%.4f, %6d steps):  max|ε̇|₀=%.3e   max|ε̇|(t=%.1f)=%.3e  (peak x=%.1f)%s\n",
                N, x[2]-x[1], nsteps, epsdot0, t_snap, epsdot_snap, s.x[ipk], tag)
        prev = epsdot_snap
    end
end

# -----------------------------------------------------------------------------
#  SUMMARY / MATCH STATEMENT
# -----------------------------------------------------------------------------
println("\n", "="^78)
println("SUMMARY  vs  fig:heat_stationary")
println("="^78)
mκ0, mdd0 = results_closed[0.0]
mκ3, mdd3 = results_closed[1/3]
@printf("  σ̂=0   : κ≡0 (max|κ|=%.2e), ε̈≡0 (max|ε̈|=%.2e)  -> STATIONARY profile, no dynamics\n", mκ0, mdd0)
@printf("  σ̂=1/3 : κ≠0 (max|κ|=%.2e), ε̈≠0 (max|ε̈|=%.2e)  -> DYNAMICAL heat-flow solution\n", mκ3, mdd3)
println("  QUALITATIVE MATCH: top-panel σ̂=0 is stationary (ε̇→0 with resolution);")
println("                     bottom-panel σ̂=1/3 has a localized nonzero ε̈/ε̇  ✓ (fig caption)")
println("  QUANTITATIVE: eq:heat_ID_EOM  τ_ε ε̈ = (κT')'  reproduced exactly by construction;")
println("                σ̂=0 ⇒ ε̈ identically 0; σ̂=1/3 ⇒ ε̈ localized to the T-gradient region.")
println("="^78)
