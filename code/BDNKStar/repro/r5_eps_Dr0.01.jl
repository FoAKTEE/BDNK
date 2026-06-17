#=
    r5_eps_Dr0.01.jl — STAGE 2 (R5) production run of the Shum smallSB-F2
    nonlinear spherical BDNK Cowling evolution at Δr = 0.01 M_⊙.

    GROUNDING: Shum, Abalos, Bea, Bezares, Figueras, Palenzuela,
    arXiv:2509.15303 ("Neutron star evolution with the BDNK framework"),
    ref-paper/sources/arXiv-2509.15303/src/Paper.tex.
      - evolved variables {γ̃E, γ̃S_r, ε, ∂_rε, ṽ^r, ∂_rṽ^r}        (l.411)
      - balance laws (l.392–393), Cowling static metric (K_ij=0, β=0)
      - reduction evolution (l.394,407,402,408)
      - SSP-RK3 (l.583); Δt/Δr=0.25 (l.583); 3rd-order FV recon (l.584);
        staggered grid (l.586); outflow BC (l.587); atmosphere reset (l.584);
        Kreiss–Oliger dissipation (l.584)
      - smallSB-F2 frame (ŝ,â,q̂,η̂,ζ̂)=(1,1,0.999,0.01,0.01)        (l.625–627)

    This driver REUSES the verified, optimized engine repro/shum_evolve_opt.jl
    (which itself includes src/BDNKStar.jl and repro/shum_core.jl).  It does NOT
    re-derive anything; it only (a) compiles via a tiny warmup, (b) picks the
    largest t_f that completes within the wallclock budget by measuring the
    realized steps/s at Δr=0.01, and (c) runs run_shum(0.01, t_f), which writes
    repro/r5_eps_Dr0.01.txt (columns "t  eps_c").
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_evolve_opt.jl")

const DR   = 0.01
const DT   = 0.25*DR
const TARGET_TF = 8000.0      # task target
const LONG_TF   = 4000.0      # PF/long fallback
# Wallclock budget for the production run (s).  Keep well under the harness cap.
const BUDGET_S  = 7800.0

println("="^78)
println("STAGE 2 (R5): nonlinear spherical BDNK Cowling evolution — Shum 2509.15303")
println("  case smallSB-F2 (ŝ,â,q̂,η̂,ζ̂)=(1,1,0.999,0.01,0.01); Δr=$DR (Δt=$DT)")
println("  target t_f=$TARGET_TF M_⊙  (PF/long fallback t_f=$LONG_TF)")
println("="^78)
println("1/M_⊙ = $(round(INVMSUN_TO_KHZ,digits=4)) kHz")
flush(stdout)

# ---- (1) compile warmup (cheap; output suppressed) -------------------------
print("compiling (warmup run_shum(0.05, 0.5)) ... "); flush(stdout)
let
    redirect_stdout(devnull) do
        run_shum(0.05, 0.5; sample_dt=1.0)
    end
end
println("done"); flush(stdout)

# ---- (2) measure realized steps/s at Δr=0.01 -------------------------------
println("\n[probe] measuring steps/s at Δr=$DR ..."); flush(stdout)
t_probe = 5.0
tp0 = time()
let
    redirect_stdout(devnull) do
        run_shum(DR, t_probe; case=:smallSB_F2, sample_dt=1.0)
    end
end
probe_el = time() - tp0
probe_steps = Int(round(t_probe/DT))
sps = probe_steps/probe_el
println("  probe: t_f=$t_probe  steps=$probe_steps  elapsed=$(round(probe_el,digits=2))s  steps/s=$(round(sps,digits=1))")
flush(stdout)

# ---- (3) choose largest feasible t_f --------------------------------------
# steps for a given t_f is t_f/DT; pick t_f so that t_f/DT/sps <= BUDGET_S.
tf_feasible = sps*BUDGET_S*DT
t_f = if tf_feasible >= TARGET_TF
    TARGET_TF
elseif tf_feasible >= LONG_TF
    LONG_TF
else
    # round DOWN to a clean multiple of 100 that fits the budget
    floor(tf_feasible/100)*100
end
t_f = max(t_f, 100.0)
println("\n[plan] feasible t_f≈$(round(tf_feasible,digits=1)) within $(BUDGET_S)s budget ⇒ running t_f=$t_f M_⊙")
nsteps_plan = Int(round(t_f/DT))
println("       planned nsteps=$nsteps_plan"); flush(stdout)

# ---- (4) PRODUCTION run ----------------------------------------------------
println("\n[RUN] run_shum($DR, $t_f) ..."); flush(stdout)
t0 = time()
ts, ecs = run_shum(DR, t_f; case=:smallSB_F2, sample_dt=1.0)
el = time() - t0
nsteps = Int(round(t_f/DT))

# ---- (5) stability / oscillation diagnostics -------------------------------
stable = all(isfinite, ecs) && all(x->0.0 < x < 1e3, ecs)
ec0 = ecs[1]; ecmin = minimum(ecs); ecmax = maximum(ecs)
# count interior local extrema of central-ε to confirm oscillation
n_extrema = 0
for i in 2:length(ecs)-1
    if (ecs[i]-ecs[i-1])*(ecs[i+1]-ecs[i]) < 0
        global n_extrema += 1
    end
end
peaks = qnm_peaks(ts, ecs)

println("\n" * "="^78)
println("RESULT  (Shum smallSB-F2, Δr=$DR)")
println("  Dr (Δr)         = $DR  M_⊙")
println("  t_f reached     = $(ts[end])  M_⊙   (requested $t_f)")
println("  #steps          = $nsteps")
println("  wallclock       = $(round(el,digits=2)) s   ($(round(el/60,digits=2)) min)")
println("  steps/s         = $(round(nsteps/el,digits=1))")
println("  central ε(0)    = $ec0")
println("  central ε range = [$ecmin, $ecmax]   (Δ=$(ecmax-ecmin))")
println("  central-ε interior extrema (oscillation count) = $n_extrema")
println("  central-ε oscillates = $(n_extrema >= 2)")
println("  stable (no NaN/blowup) = $stable")
println("  central-ε QNM peaks (kHz): " *
        join([string(round(p,digits=3)) for p in peaks[1:min(6,end)]], ", "))
println("  saved: /data/haiyangw/claude/BDNK/code/BDNKStar/repro/r5_eps_Dr0.01.txt")
println("="^78)
flush(stdout)
