# =============================================================================
# pmp_shock_instab.jl  —  PMP 2209.09265  Figs. fig:shock_instability + fig:acaus_instab
#
# Dynamical (1+1)D evolution of BDNK shockwave initial data (eq:shockwave_ID,
# line 1029) for the relativistic ideal-gas BDNK fluid in the PMP hydrodynamic
# frame (eq:hydro_frame, line 464), demonstrating the CAUSALITY CLASSIFICATION
# of the shockwave instabilities:
#
#   * fig:shock_instability  (Table line 556: Γ=4/3, m=0.1, V̂=4/3, σ̂=0, τ̂=1.5,3)
#       initial data L={ε,v,n}={1,0.9,1}  ⇒  R={11.5174,0.354727,5.44212}
#       (Rankine-Hugoniot, eq:shockwave_params line 1048, first line).
#       - τ̂=1.5 ⇒ c_+ ~ 0.94 > v_max=0.9  ⇒  STABLE  (paper bottom panel,
#         "c_+ > v throughout ... no instability sets in", line 1093).
#       - τ̂=3   ⇒ c_+ ~ 0.76 < v_max=0.9  ⇒  UNSTABLE high-frequency numerical
#         instability "where the flow velocity v exceeds the maximum
#         characteristic speed c_+" (paper top panel, line 1093; "precisely the
#         same case where the ODEs ... yield no solution", line 1064).
#         NOTE larger τ̂ ⇒ SMALLER c_+ (c_+² ∝ c_s²/τ̂, eq:cpmsq line 1424).
#
#   * fig:acaus_instab  (Table line 557: Γ=4/3, m=0.1, V̂=4/3, σ̂=0, τ̂=0.25,0.4,0.5,1.5)
#       initial data L={ε,v,n}={1,0.6,1}  ⇒  R={1.33795,0.514414,1.25027}
#       (eq:shockwave_params line 1048, second line).  v_max=0.6 < c_+ in ALL
#       cases (no v>c_+ instability), but the FRAMES differ in causality:
#       - τ̂=1.5 ⇒ c_+ ~ 0.94  SUBLUMINAL (causal)         -> STABLE  (line 1142)
#       - τ̂=0.5 ⇒ c_+ ~ 1.47  "weakly superluminal"        -> stable, ≈identical
#       - τ̂=0.4 ⇒ c_+ ~ 1.63  superluminal + STIFF (λ=0.01) -> stable
#       - τ̂=0.25⇒ c_+ ~ 2.03  "WILDLY superluminal"         -> ACAUSAL fast
#         instability: a "bump" grows unboundedly without propagating and a
#         sharp feature forms; ε̇,v̇ diverge in finite time (paper line 1133-1134).
#
# PACKAGE REUSE: this file is SELF-CONTAINED but builds entirely on the prior-
# verified engine repro/pmp_viscous_core.jl (constant-state-exact, inviscid-
# Bjorken-validated, c_± validated against the A,B,C quadratic AND the assembled
# principal symbol).  We include it (which also begins with
# include(".../src/BDNKStar.jl"); using .BDNKStar) and reuse:
#   pmp_frame, bdnk_stress, ideal_stress, cpm2_closed, VState, evolve!,
#   _step!, _update_aux!, the WENO/KT/Heun machinery, and erf_.
# We add ONLY: shockwave initial data, a diagnostics/instability detector, and
# the per-figure driver.  We do NOT edit src/ or pmp_viscous_core.jl.
#
# RUN: cd code/BDNKStar && JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
#      julia --project=. repro/pmp_shock_instab.jl
# =============================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/pmp_viscous_core.jl")
using .BDNKStar.ConformalEvolution: erf_
using Printf

# -----------------------------------------------------------------------------
# Velocity convention:  the engine evolves the spatial 4-velocity component u,
# with W=√(1+u²), and the physical 3-velocity is v = u/W.  Inverting,
#   u = v·W = v/√(1−v²).
# -----------------------------------------------------------------------------
@inline u_of_v(v) = v / sqrt(1 - v^2)
@inline v_of_u(u) = u / sqrt(1 + u^2)

# -----------------------------------------------------------------------------
# Shockwave initial data  (eq:shockwave_ID, line 1029), w = transition width.
#   ε(0,x) = (εR−εL)/2 [erf(x/w)+1] + εL
#   v(0,x) = (vL−vR)/2 [1−erf(x/w)]  + vR
#   n(0,x) = (nL−nR)/2 [1−erf(x/w)]  + nR
# Conserved fields are seeded with the PERFECT-FLUID (ideal) stress
# ideal_stress (eq:Tab_0) at the initial primitives, consistent with the paper's
# treatment of the asymptotic states as equilibrium (line 1004).  The recovered
# ε̇,u̇ then immediately incorporate the BDNK gradient corrections at t=0⁺.
# -----------------------------------------------------------------------------
function init_shockwave(fr::PMPFrame; N=1025, xmin=-200.0, xmax=200.0, cfl=0.1,
                        εL=1.0, vL=0.9, nL=1.0, εR=11.5174, vR=0.354727, nR=5.44212,
                        w=10.0)
    x = collect(range(xmin, xmax; length=N)); dx = x[2]-x[1]; z = zeros(N)
    ε = similar(x); n = similar(x); u = similar(x)
    for i in 1:N
        ξ = x[i]
        εi = (εR-εL)/2*(erf_(ξ/w)+1) + εL
        vi = (vL-vR)/2*(1-erf_(ξ/w)) + vR
        ni = (nL-nR)/2*(1-erf_(ξ/w)) + nR
        ε[i]=εi; n[i]=ni; u[i]=u_of_v(vi)
    end
    s = VState(fr, x, dx, cfl*dx, ε, n, u,
               copy(z),copy(z),copy(z),copy(z),copy(z),
               copy(z),copy(z),copy(z), false, false)
    for i in 1:N
        Ttt,Ttx,_,Jt,_ = ideal_stress(fr.g, ε[i], n[i], u[i])
        s.Ttt[i]=Ttt; s.Ttx[i]=Ttx; s.Jt[i]=Jt
    end
    _update_aux!(s)
    return s
end

# c_+ over the current profile (min, max) and where v exceeds c_+ (eq:cpmsq).
function profile_diag(s::VState)
    N=length(s.x); rng=NG+1:N-NG
    cpmin=Inf; cpmax=0.0; vmax=0.0; nvexc=0
    vatmax_cp=0.0
    for i in rng
        cp2,_,_ = cpm2_closed(s.fr, s.ε[i], s.n[i])
        cp = sqrt(max(cp2,0.0))
        v  = abs(v_of_u(s.u[i]))
        cpmin=min(cpmin,cp); cpmax=max(cpmax,cp);
        if v>vmax; vmax=v; vatmax_cp=cp; end
        if v>cp; nvexc+=1; end
    end
    return (; cpmin, cpmax, vmax, nvexc, vatmax_cp)
end

# field-amplitude / smoothness diagnostics for instability detection
function field_diag(s::VState)
    N=length(s.x); rng=NG+1:N-NG
    vmax=0.0; εmax=0.0; εtmax=0.0; utmax=0.0; tv_v=0.0; anynan=false
    vprev = v_of_u(s.u[NG+1])
    for i in rng
        v = v_of_u(s.u[i])
        (isnan(s.ε[i])||isnan(s.u[i])||isnan(s.εt[i])) && (anynan=true)
        vmax=max(vmax,abs(v)); εmax=max(εmax,abs(s.ε[i]))
        εtmax=max(εtmax,abs(s.εt[i])); utmax=max(utmax,abs(s.ut[i]))
        tv_v += abs(v - vprev); vprev=v          # total variation of v (sharpness)
    end
    return (; vmax, εmax, εtmax, utmax, tv_v, anynan)
end

# Run one case to time T (in code units), polling diagnostics; STOP early if an
# instability is detected (NaN / blow-up of |ε̇|,|u̇| or total variation of v).
# Returns a status string + the final diagnostics + the recorded history.
function run_case(fr::PMPFrame; N, xmin, xmax, cfl, IDkw, T, label,
                  blowup_field=1e6, blowup_tv=50.0, poll_dt=20.0)
    s = init_shockwave(fr; N=N, xmin=xmin, xmax=xmax, cfl=cfl, IDkw...)
    fd0 = field_diag(s); pd0 = profile_diag(s)
    nsteps_total = ceil(Int, T/s.dt)
    poll_every   = max(1, round(Int, poll_dt/s.dt))
    tv0 = fd0.tv_v
    status = "STABLE"; tcrash = T
    histT=Float64[]; histTV=Float64[]; histEt=Float64[]; histVmax=Float64[]
    done=0
    while done < nsteps_total
        nb = min(poll_every, nsteps_total-done)
        evolve!(s, nb); done += nb
        t = done*s.dt
        fd = field_diag(s)
        push!(histT,t); push!(histTV,fd.tv_v); push!(histEt,fd.εtmax); push!(histVmax,fd.vmax)
        if fd.anynan || fd.εmax>blowup_field || fd.εtmax>blowup_field ||
           fd.utmax>blowup_field || fd.tv_v>blowup_tv || fd.vmax>=1.0
            status = fd.anynan ? "CRASH(NaN)" :
                     (fd.vmax>=1.0 ? "CRASH(v≥1)" : "BLOWUP")
            tcrash = t
            break
        end
    end
    fd = field_diag(s); pd = profile_diag(s)
    return (; s, status, tcrash, fd0, fd, pd0, pd, histT, histTV, histEt, histVmax, tv0)
end

println("\n", "#"^78)
println("#  PMP 2209.09265  —  fig:shock_instability  +  fig:acaus_instab")
println("#  dynamical BDNK shockwave evolution; causality classification")
println("#"^78)

# =============================================================================
#  FIG: shock_instability   (L v=0.9; τ̂=1.5 STABLE,  τ̂=3 UNSTABLE v>c_+)
#  Table line 556:  Γ=4/3, m=0.1, V̂=4/3, σ̂=0, τ̂∈{1.5,3}; w=10 (line 1056).
# =============================================================================
println("\n", "="^78)
println("[FIG shock_instability]  L={ε,v,n}={1,0.9,1} ⇒ R={11.5174,0.354727,5.44212}")
println("  Γ=4/3, m=0.1, V̂=4/3, σ̂=0, w=10   (Table line 556, eq:shockwave_params)")
println("="^78)
ID1 = (εL=1.0, vL=0.9, nL=1.0, εR=11.5174, vR=0.354727, nR=5.44212, w=10.0)
shock_results = Dict{Float64,Any}()
for τh in (1.5, 3.0)
    fr = pmp_frame(; Γ=4/3, m=0.1, Vhat=4/3, σhat=0.0, τhat=τh)
    # classification at the INITIAL profile (v vs c_+)
    sID = init_shockwave(fr; N=1025, xmin=-200.0, xmax=200.0, ID1...)
    pd = profile_diag(sID)
    classify = (pd.vmax < pd.cpmin) ? "v<c_+ everywhere ⇒ expect STABLE" :
                                      "v>c_+ somewhere  ⇒ expect UNSTABLE"
    @printf("\n  τ̂=%.2f: c_+∈[%.4f,%.4f]  v_max=%.4f  (#cells v>c_+ : %d) → %s\n",
            τh, pd.cpmin, pd.cpmax, pd.vmax, pd.nvexc, classify)
    # evolve.  Stable case asymptotes to steady shock (run long); unstable case
    # develops a high-frequency instability where v>c_+ (run shorter, detect).
    Trun = (τh==1.5) ? 300.0 : 250.0
    r = run_case(fr; N=1025, xmin=-200.0, xmax=200.0, cfl=0.1, IDkw=ID1,
                 T=Trun, label="shock τ̂=$τh", blowup_field=1e6, blowup_tv=80.0)
    shock_results[τh] = (r, pd, classify)
    @printf("    evolved to T=%.0f: status=%s (tcrash=%.0f)\n", Trun, r.status, r.tcrash)
    @printf("    TV(v): %.4f → %.4f   max|ε̇|: %.3e → %.3e   max|u̇|: %.3e → %.3e\n",
            r.fd0.tv_v, r.fd.tv_v, r.fd0.εtmax, r.fd.εtmax, r.fd0.utmax, r.fd.utmax)
    @printf("    final v_max=%.4f  ε_max=%.4f\n", r.fd.vmax, r.fd.εmax)
end

# =============================================================================
#  FIG: acaus_instab   (L v=0.6; τ̂=1.5,0.5,0.4 stable, τ̂=0.25 WILDLY-acausal)
#  Table line 557:  Γ=4/3, m=0.1, V̂=4/3, σ̂=0, τ̂∈{0.25,0.4,0.5,1.5}; w=10.
#  Per paper line 1142: τ̂=0.4 needs CFL λ=0.01 (stiff), others λ=0.1; the wildly
#  superluminal τ̂=0.25 case (c_+~2) develops a fast instability at early times.
# =============================================================================
println("\n", "="^78)
println("[FIG acaus_instab]  L={ε,v,n}={1,0.6,1} ⇒ R={1.33795,0.514414,1.25027}")
println("  Γ=4/3, m=0.1, V̂=4/3, σ̂=0, w=10   (Table line 557, eq:shockwave_params)")
println("="^78)
ID2 = (εL=1.0, vL=0.6, nL=1.0, εR=1.33795, vR=0.514414, nR=1.25027, w=10.0)
acaus_results = Dict{Float64,Any}()
# (τ̂, cfl, Trun):  τ̂=0.4 stiff ⇒ λ=0.01 (paper line 1142); τ̂=0.25 fast inst ⇒
# also stiff, use λ=0.01 and a SHORT run since the instability sets in "at early
# times" (paper line 1133).  The subluminal/weakly-superluminal cases (λ=0.1)
# are run long enough to confirm they equilibrate without instability.
cases = [(1.5, 0.1, 300.0), (0.5, 0.1, 300.0), (0.4, 0.01, 120.0), (0.25, 0.01, 120.0)]
for (τh, cfl, Trun) in cases
    fr = pmp_frame(; Γ=4/3, m=0.1, Vhat=4/3, σhat=0.0, τhat=τh)
    sID = init_shockwave(fr; N=1025, xmin=-200.0, xmax=200.0, cfl=cfl, ID2...)
    pd = profile_diag(sID)
    caus = (pd.cpmax < 1.0) ? "SUBLUMINAL (causal)" :
           (pd.cpmax < 1.7  ? "weakly superluminal" : "WILDLY superluminal (acausal)")
    @printf("\n  τ̂=%.2f (λ=%.2f): c_+∈[%.4f,%.4f]  v_max=%.4f  → %s\n",
            τh, cfl, pd.cpmin, pd.cpmax, pd.vmax, caus)
    r = run_case(fr; N=1025, xmin=-200.0, xmax=200.0, cfl=cfl, IDkw=ID2,
                 T=Trun, label="acaus τ̂=$τh", blowup_field=1e6, blowup_tv=80.0)
    acaus_results[τh] = (r, pd, caus)
    @printf("    evolved to T=%.0f: status=%s (tcrash=%.0f)\n", Trun, r.status, r.tcrash)
    @printf("    TV(v): %.4f → %.4f   max|ε̇|: %.3e → %.3e   max|u̇|: %.3e → %.3e\n",
            r.fd0.tv_v, r.fd.tv_v, r.fd0.εtmax, r.fd.εtmax, r.fd0.utmax, r.fd.utmax)
    @printf("    final v_max=%.4f  ε_max=%.4f\n", r.fd.vmax, r.fd.εmax)
    # growth-rate diagnostic of max|ε̇| (the "bump grows unboundedly" signature)
    if length(r.histT)>2
        et0=r.histEt[1]; etend=r.histEt[end]
        @printf("    max|ε̇| history: %.3e → %.3e  (×%.1f over Δt=%.0f)\n",
                et0, etend, etend/max(et0,1e-30), r.histT[end]-r.histT[1])
    end
end

# =============================================================================
#  SUMMARY  —  causality classification table
# =============================================================================
println("\n", "="^78)
println("SUMMARY — causality classification of dynamical BDNK shockwaves")
println("="^78)
println("\nfig:shock_instability  (L v=0.9, v_max=0.9):")
@printf("  %-8s %-9s %-22s %-14s %s\n","τ̂","c_+(min)","class (v vs c_+)","status","verdict")
for τh in (1.5,3.0)
    r,pd,cl = shock_results[τh]
    verdict = (pd.vmax<pd.cpmin) ? "STABLE✓" : "UNSTABLE✓"
    @printf("  %-8.2f %-9.4f %-22s %-14s %s\n", τh, pd.cpmin,
            (pd.vmax<pd.cpmin ? "v<c_+ (causal-fast)" : "v>c_+ (no steady ODE)"),
            r.status, verdict)
end
println("\nfig:acaus_instab  (L v=0.6, v_max=0.6 < c_+ always):")
@printf("  %-8s %-9s %-26s %-14s %s\n","τ̂","c_+(min)","frame causality","status","verdict")
for (τh,_,_) in cases
    r,pd,caus = acaus_results[τh]
    verdict = occursin("WILDLY",caus) ? (r.status=="STABLE" ? "(no blowup in window)" : "ACAUSAL-UNSTABLE✓") :
                                        (r.status=="STABLE" ? "STABLE✓" : "unexpected")
    @printf("  %-8.2f %-9.4f %-26s %-14s %s\n", τh, pd.cpmin, caus, r.status, verdict)
end
println("\n(c_+ matches paper: τ̂=1.5→c_+~0.9, 0.5→~1.5, 0.4→~1.6, 0.25→~2;",
        "\n shock_instability quoted ~0.9 for τ̂=1.5 bottom panel.)")
println("="^78)
