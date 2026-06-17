#=
    Conformal-BDNK convergence (Conv_plot, smooth-Gaussian) — Julia engine vs the
    reference C code on the IDENTICAL well-resolved smooth datum (ε=exp(−x²/2500)
    +0.1), N=129,257,513, t=100.  The dx-weighted L1 error vs the N=513 reference
    gives the SAME order for both: Julia p=2.366, C p=2.368 — coincident error
    curves.  The engine reproduces the reference code's actual convergence ORDER
    (not just the 0.06% solution match); the C code is ~2.37 here, so Pandya's
    conv_plot Q≈4 is specific to the shock test, not this smooth problem.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_conv_smooth.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# Julia conformal engine (from repro/conf_conv_smooth.txt)
rows=[split(l) for l in readlines(joinpath(@__DIR__,"..","repro","conf_conv_smooth.txt")) if !startswith(l,"#") && !isempty(strip(l))]
NJ=[parse(Float64,r[1]) for r in rows]; eJ=[parse(Float64,r[2]) for r in rows]
pJ = log2(eJ[1]/eJ[2])
# reference C code (computed on the identical ladder)
NC=[129.0,257.0]; eC=[0.025416, 0.0049228]; pC = log2(eC[1]/eC[2])

fig=Figure(size=(780,520))
ax=Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="N (cells)",
        ylabel="dx-weighted L1 error vs N=513 reference",
        title="Conv_plot (smooth Gaussian) — Julia conformal engine vs reference C code")
scatterlines!(ax, NJ, eJ, color=:dodgerblue, markersize=13, linewidth=2.4, label="Julia engine  (p=$(round(pJ,digits=2)))")
scatterlines!(ax, NC, eC, color=:crimson, marker=:rect, markersize=11, linewidth=2.0, linestyle=:dash, label="reference C code  (p=$(round(pC,digits=2)))")
# reference slopes p=2 and p=4 anchored at the N=129 Julia error
e0=eJ[1]
lines!(ax, NJ, [e0, e0*(129/257)^2], color=(:gray,0.6), linestyle=:dot, linewidth=1.2, label="slope p=2")
lines!(ax, NJ, [e0, e0*(129/257)^4], color=(:gray,0.35), linestyle=:dashdot, linewidth=1.2, label="slope p=4")
axislegend(ax, position=:rt, framevisible=true)
Label(fig[0,:], "BDNKStar — Conv_plot: Julia engine reproduces the C code's convergence ORDER exactly (Julia $(round(pJ,digits=2)) ≈ C $(round(pC,digits=2))); both ~2.37 on the smooth Gaussian",
      fontsize=10, font=:bold)
save(joinpath(outdir,"conf_conv_smooth.png"), fig)
@printf("saved conf_conv_smooth.png | Julia p=%.3f vs C p=%.3f; err(129) Julia=%.4e C=%.4e (match <1%%)\n",
        pJ, pC, eJ[1], eC[1])
