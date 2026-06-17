#=
    Pandya 2201.12317 — steady planar shock (SMOOTH_SHOCK ID), Julia engine vs the
    reference C code.  Same setup as parameters.h with ID_TYPE=SMOOTH_SHOCK:
      conformal PMP-luminal frame eta0=10^{1/4}/(3π); εL=1, vL=0.8 → εR≈4.4074
      x∈[−200,200], N=129, CFL=0.1, dt=0.3125, save every 10 steps.

    Run: julia code/BDNKStar/repro/conf_shock_overlay.jl
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame

eta0 = 10.0^0.25/(3π)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)
s    = init_smooth_shock(fr; N=129, xmin=-200.0, xmax=200.0, εL=1.0, vL=0.8, cfl=0.1)

rows  = [1, 16, 31, 51, 102]
nstep = [10*(r-1) for r in rows]
x = copy(s.x); cols = Vector{Vector{Float64}}(); done = 0
push!(cols, copy(energy_density(s)))
println("  t=0  εL=", round(minimum(cols[1]),sigdigits=5), "  εR=", round(maximum(cols[1]),sigdigits=6))
for k in 2:length(nstep)
    evolve!(s, nstep[k]-done); global done = nstep[k]
    push!(cols, copy(energy_density(s)))
    println("  evolved to t=$(round(done*s.dt,digits=3))  εmax=", round(maximum(cols[end]),sigdigits=6),
            " nan=", any(!isfinite, cols[end]))
end
open(joinpath(@__DIR__,"conf_shock_overlay.txt"),"w") do io
    println(io, "# x  t0 t46.875 t93.75 t156.25 t315.625  (Julia smooth-shock vs C rows ", rows, ")")
    for i in eachindex(x)
        print(io, x[i]); for c in cols; print(io, "  ", c[i]); end; println(io)
    end
end
println("SAVED conf_shock_overlay.txt")
