#=
    Conformal-BDNK engine — convergence order on a WELL-RESOLVED SMOOTH datum
    (wide Gaussian ε=exp(−x²/2500)+0.1, half-width ~50, ≫ dx).  The earlier
    conf_convergence (narrow clump, w=25) gave p≈1.35 because the steep feature
    was under-resolved at those N; on a smooth well-resolved field the true scheme
    order should emerge.  Uses N=513,1025,2049 (aligned), N=2049 as the reference,
    locked dt (isolate spatial order).

    Run: julia code/BDNKStar/repro/conf_conv_smooth.jl
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame

eta0 = 10.0^0.25/(3π)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)

# evolve a smooth wide Gaussian at three aligned resolutions to t=40, locked dt
function run(N, nsteps)
    s = init_gaussian(fr; N=N, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=2500.0, c=0.1, cfl=0.1)
    evolve!(s, nsteps); return energy_density(s)
end
u513  = run(513, 2040)
u1025 = run(1025, 4080)
u2049 = run(2049, 8160)

# dx-weighted L1 error vs the N=2049 reference
errN(uN, uR, N) = (step = 2048 ÷ (N-1); sum(abs.(uN .- uR[1:step:end])) * 400.0/(N-1))
e513  = errN(u513,  u2049, 513)
e1025 = errN(u1025, u2049, 1025)
p = log2(e513/e1025)
open(joinpath(@__DIR__,"conf_conv_smooth.txt"),"w") do io
    println(io, "# N  err_vs_2049")
    println(io, "513  ", e513)
    println(io, "1025 ", e1025)
    println(io, "# order p = ", p)
end
println("SMOOTH conformal convergence: err(513)=", round(e513,sigdigits=4),
        " err(1025)=", round(e1025,sigdigits=4), " -> ORDER p=", round(p,digits=3))
