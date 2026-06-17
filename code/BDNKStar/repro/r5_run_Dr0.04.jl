#=
    r5_run_Dr0.04.jl — STAGE 2 (R5) production run of the Shum smallSB-F2
    nonlinear spherical BDNK Cowling evolution at Δr=0.04 M_⊙.

    Calls run_shum(0.04, t_f) from repro/shum_evolve_opt.jl (PHYSICS unchanged).
    GROUNDING: Shum et al., arXiv:2509.15303, ref-paper/sources/arXiv-2509.15303/src/Paper.tex.
      EVOLVED VARS l.411; BALANCE LAWS l.392-393; REDUCTION l.394,407,402,408;
      NUMERICS l.582-587 (SSP-RK3, Δt/Δr=0.25, 3rd-order FV, staggered grid,
      outflow BC, atmosphere reset, KO dissipation).  case smallSB-F2 l.625-627.

    Writes repro/r5_eps_Dr0.04.txt with columns "t  eps_c".
    Reports: Dr, t_f reached, #steps, wallclock, stable(no NaN), central-eps oscillation.
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar

# Bring in the optimized evolution driver (defines run_shum, run_evolution, etc.).
# It includes shum_core.jl internally; silence any include-time stdout.
let
    redirect_stdout(devnull) do
        include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_evolve_opt.jl")
    end
end

const Dr = 0.04

# --- choose t_f from a short benchmark so the run completes in the time budget ---
println("="^78)
println("[R5 run Dr=0.04] Shum smallSB-F2 nonlinear BDNK Cowling (2509.15303)")
println("  Δr=$Dr  Δt=$(0.25*Dr)  N=$(Int(round(20.0/Dr))) cells")
println("="^78)
flush(stdout)

# warmup / compile (also exercises full code path once)
print("[warmup] compiling..."); flush(stdout)
twc = @elapsed run_shum(0.05, 0.5; sample_dt=1.0)
println(" done ($(round(twc,digits=1))s)"); flush(stdout)

# short timing probe at Δr=0.04 to estimate steps/s
t_probe = 4.0
nb = Int(round(t_probe/(0.25*Dr)))
println("[probe] t_f=$t_probe ($nb steps) to estimate throughput..."); flush(stdout)
tb = @elapsed run_shum(Dr, t_probe; sample_dt=1.0)
sps = nb/tb
println("[probe] elapsed=$(round(tb,digits=2))s  steps/s=$(round(sps,digits=1))"); flush(stdout)

# Time budget ~ 35 min for the production run; pick the largest of the
# preferred ladder {8000,4000,2000,1000} that fits, but at least 1000.
budget_s = 35*60.0
function steps_for(tf); Int(round(tf/(0.25*Dr))); end
tf_target = 8000.0
for cand in (8000.0, 4000.0, 2000.0, 1000.0)
    if steps_for(cand)/sps <= budget_s
        global tf_target = cand
        break
    end
    global tf_target = cand   # fallback to smallest if none fit
end
tf = tf_target
nsteps = steps_for(tf)
println("[choose] t_f=$tf  ($nsteps steps, est $(round(nsteps/sps/60,digits=1)) min)"); flush(stdout)

# --- production run ---
println("[RUN] run_shum($Dr, $tf) ..."); flush(stdout)
t0 = time()
ts, ecs = run_shum(Dr, tf; case=:smallSB_F2, sample_dt=1.0)
wall = time() - t0

t_reached = ts[end]
ec0 = ecs[1]
stable = all(isfinite, ecs) && all(<(1e3), ecs)
ecmin, ecmax = minimum(ecs), maximum(ecs)
# count sign changes of (eps_c - mean) as an oscillation diagnostic
m = sum(ecs)/length(ecs)
sc = 0
for k in 2:length(ecs)
    ((ecs[k]-m)*(ecs[k-1]-m) < 0) && (global sc += 1)
end
osc = sc >= 2
rel_amp = (ecmax-ecmin)/ec0

println("="^78)
println("[RESULT] Dr=$Dr")
println("  t_f reached      = $(round(t_reached,digits=2)) M_sun  (target t_f=$tf, ladder target 8000)")
println("  #steps           = $nsteps")
println("  wallclock        = $(round(wall,digits=2)) s  ($(round(wall/60,digits=2)) min)")
println("  steps/s          = $(round(nsteps/wall,digits=1))")
println("  stable (no NaN)  = $stable")
println("  eps_c[0]         = $ec0")
println("  eps_c min/max    = $ecmin / $ecmax")
println("  rel amplitude    = $(round(rel_amp,sigdigits=4))")
println("  oscillates       = $osc  (mean-crossings=$sc)")
peaks = qnm_peaks(ts, ecs)
println("  central-eps QNM peaks (kHz): " *
        join([string(round(p,digits=3)) for p in peaks[1:min(6,end)]], ", "))
println("  saved: /data/haiyangw/claude/BDNK/code/BDNKStar/repro/r5_eps_Dr0.04.txt")
println("="^78)
flush(stdout)
