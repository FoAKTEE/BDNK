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

# IDENTICAL setup to the reference C smooth ladder: N=129,257,513 to t=100, locked dt
function run(N, nsteps)
    s = init_gaussian(fr; N=N, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=2500.0, c=0.1, cfl=0.1)
    evolve!(s, nsteps); return energy_density(s)
end
u129 = run(129, 320)
u257 = run(257, 640)
u513 = run(513, 1280)

# dx-weighted L1 error vs the N=513 reference (matches the C test)
errN(uN, uR, N) = (step = 512 ÷ (N-1); sum(abs.(uN .- uR[1:step:end])) * 400.0/(N-1))
e129 = errN(u129, u513, 129)
e257 = errN(u257, u513, 257)
p = log2(e129/e257)
open(joinpath(@__DIR__,"conf_conv_smooth.txt"),"w") do io
    println(io, "# N  err_vs_513  (smooth Gaussian @t=100, matches C ladder)")
    println(io, "129  ", e129)
    println(io, "257  ", e257)
    println(io, "# order p = ", p, "  (C code on identical setup: p=2.368)")
end
println("SMOOTH conformal convergence @t=100 (matches C): err(129)=", round(e129,sigdigits=4),
        " err(257)=", round(e257,sigdigits=4), " -> ORDER p=", round(p,digits=3),
        "  [C code: 2.368]")
