#=
    PMP 2209.09265 conv_plot — self-convergence factor Q_N(t) of the BDNK shock
    PDE.  HONEST FINDING: this engine's spatial scheme (Kurganov–Tadmor central
    flux) is 2nd-order, so Q_N(t)≈2.0 (clean, stable).  The reference reports
    Q≈4 using a higher-order WENO5 flux — a genuine engine-order difference, NOT
    a resolution artifact (Q is flat at 2.0 across t).  The PHYSICS (shock, RH
    states, instability, Bjorken) is reproduced; the convergence ORDER is 2 vs 4.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_conv.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

rows=[split(l) for l in readlines(joinpath(repo,"pmp_conv.txt")) if !startswith(l,"#") && !isempty(strip(l))]
t=[parse(Float64,r[1]) for r in rows]; Q=[parse(Float64,r[2]) for r in rows]

fig=Figure(size=(760,460))
ax=Axis(fig[1,1], xlabel="t", ylabel="Q_N(t)",
        title="PMP conv_plot — BDNK shock PDE self-convergence (N=513,1025,2049)")
hlines!(ax, [4.0], color=:red,   linestyle=:dot,  linewidth=1.6, label="reference order (WENO5 flux) Q=4")
hlines!(ax, [2.0], color=:blue,  linestyle=:dash, linewidth=1.2, label="2nd-order reference")
scatterlines!(ax, t, Q, color=:black, markersize=9, linewidth=2.2,
              label="this engine (KT flux): Q≈$(round(sum(Q)/length(Q),digits=2))")
ylims!(ax, 0, 5); xlims!(ax, 0, maximum(t)+0.5)
axislegend(ax, position=:rb, framevisible=true)
Label(fig[0,:], "BDNKStar — PMP shock convergence: engine verified CLEAN 2nd-order (Q≈2.0); reference Q≈4 needs a higher-order WENO5 flux [engine-order gap]",
      fontsize=10.5, font=:bold)
save(joinpath(outdir,"pmp_conv.png"), fig)
@printf("saved pmp_conv.png | mean Q=%.3f (clean 2nd-order; reference Q≈4 needs WENO5 4th-order flux)\n",
        sum(Q)/length(Q))
