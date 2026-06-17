#=
    Pandya 2201.12317 Conv_plot — self-convergence of the conformal-BDNK Gaussian
    clump.  Evolve the Julia engine at N=129,257,513 (RES=1,2,4; same final time
    t=318.75, dt=CFL·dx) and save the final ε(x); the viz computes the Richardson
    self-convergence order and overlays it on the reference C code's order.

    Run: julia code/BDNKStar/repro/conf_convergence.jl
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame

eta0 = 10.0^0.25/(3π)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)

for (N, nsteps) in [(129,1020),(257,2040),(513,4080)]
    s = init_gaussian(fr; N=N, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)
    evolve!(s, nsteps)
    ε = energy_density(s)
    open(joinpath(@__DIR__,"conf_conv_N$(N).txt"),"w") do io
        println(io, "# x  eps   (N=$N, t=$(round(nsteps*s.dt,digits=3)), εmax=$(round(maximum(ε),sigdigits=6)))")
        for i in eachindex(s.x); println(io, s.x[i], "  ", ε[i]); end
    end
    println("SAVED conf_conv_N$(N).txt  εmax=", round(maximum(ε),sigdigits=6),
            " nan=", any(!isfinite, ε))
end
println("JULIA_CONVERGENCE_DONE")
