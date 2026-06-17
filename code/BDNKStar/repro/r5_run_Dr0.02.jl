#=
    r5_run_Dr0.02.jl — driver for [R5/run Dr=0.02].
    Uses repro/shum_evolve_opt.jl run_shum(Dr, t_f) for the Shum smallSB-F2
    nonlinear Cowling BDNK evolution (arXiv:2509.15303,
    ref-paper/sources/arXiv-2509.15303/src/Paper.tex, l.392-393 balance laws,
    l.582-587 numerics).  Δr=0.02, Δt=0.25·Δr=0.005.

    Output: repro/r5_eps_Dr0.02.txt (columns: t  eps_c).
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar

# Bring in run_shum and helpers (silence its own driver: PROGRAM_FILE guard
# only fires when shum_evolve_opt.jl is the script, so plain include is safe).
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_evolve_opt.jl")

const Dr   = 0.02
const T_F  = 8000.0          # target; PF/long runs may fall back to 4000

println("="^72)
println("[R5/run Dr=0.02] Shum smallSB-F2 nonlinear Cowling BDNK evolution")
println("  Δr=$Dr  Δt=$(0.25*Dr)  target t_f=$T_F M_⊙")
println("="^72)
flush(stdout)

# warmup / compile (short)
run_shum(0.05, 0.5; sample_dt=1.0)

t0 = time()
ts, ecs = run_shum(Dr, T_F; case=:smallSB_F2, sample_dt=1.0)
el = time() - t0

t_reached = ts[end]
nsteps    = Int(round(t_reached/(0.25*Dr)))
stable    = all(isfinite, ecs)
ec0       = ecs[1]
# central-eps oscillation diagnostics
ecmin = minimum(ecs); ecmax = maximum(ecs)
amp   = (ecmax - ecmin)
# count sign changes of (ec - mean) -> oscillation evidence
m = sum(ecs)/length(ecs)
sgn = 0; crossings = 0
for e in ecs
    s = sign(e - m)
    if s != 0
        if sgn != 0 && s != sgn
            crossings += 1
        end
        sgn = s
    end
end

println("\n[RESULT]")
println("  Dr            = $Dr")
println("  t_f reached   = $t_reached M_⊙  (target $T_F)")
println("  nsteps        = $nsteps")
println("  wallclock     = $(round(el,digits=2)) s  ($(round(el/60,digits=2)) min)")
println("  steps/s       = $(round(nsteps/el,digits=1))")
println("  stable (finite)= $stable")
println("  eps_c[0]      = $ec0")
println("  eps_c min/max = $ecmin / $ecmax")
println("  osc amplitude = $amp   (rel $(round(amp/ec0,digits=6)))")
println("  mean-crossings= $crossings  (oscillation: $(crossings>=4))")
println("  samples       = $(length(ts))")
println("  output file   = repro/r5_eps_Dr0.02.txt")
flush(stdout)
println("DONE")
