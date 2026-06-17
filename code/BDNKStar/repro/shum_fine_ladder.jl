#=
    Shum fine-Δr ladder for error_fit / decay-continuum.  The default σKO=0.5
    is unstable at Δr≲0.0032 (blows up ~t=4); σKO=2.0 stabilises (KO damps the
    grid-scale instability, not the smooth f-mode).  Runs ONE Δr (from ARGS) at
    σKO=2.0 and saves r5_eps_Dr<Δr>.txt.

    Launch (one per Δr, parallel):  julia shum_fine_ladder.jl 0.0032
=#
include(joinpath(@__DIR__, "shum_evolve_opt.jl"))

Dr  = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.0032
t_f = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2000.0
σKO = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 2.0

println("[fine ladder] Δr=$Dr  t_f=$t_f  σKO=$σKO  (N≈$(round(Int,20/Dr)))"); flush(stdout)
s, ts, ecs, nan_hit, ec0 = run_evolution(; dr=Dr, t_f=t_f, vpert=0.0, epspert=1e-4,
                                          sample_dt=1.0, σKO=σKO, case=:smallSB_F2)
outpath = joinpath(@__DIR__, "r5_eps_Dr$(Dr).txt")
open(outpath, "w") do io
    println(io, "# t  eps_c   (Shum smallSB-F2, fine Δr=$Dr, σKO=$σKO, M_T=1.4)")
    for k in 1:length(ts); println(io, ts[k], "  ", ecs[k]); end
end
println("[fine ladder] Δr=$Dr DONE: nan_hit=$nan_hit  rows=$(length(ts))  finite=$(all(isfinite,ecs))  ",
        "εc range=$(round(minimum(ecs),sigdigits=5))-$(round(maximum(ecs),sigdigits=5))  -> $outpath")
