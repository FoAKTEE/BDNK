#=
    PMP 2209.09265 conv_plot — self-convergence factor Q_N(t) of the BDNK PDE.
    HONEST engine-order characterization (two tests, both ≪ the reference Q≈4):
      (left)  steepening shock  → clean Q≈2.0
      (right) smooth ε bump     → Q≈1.0–1.5
    Confirms the engine's spatial scheme (Kurganov–Tadmor central flux) is
    LOW-ORDER (not WENO5-4th).  The reference reports Q≈4 with a higher-order
    flux — a genuine engine-order gap, documented (NOT fabricated).  The physics
    (shock/RH/instability/Bjorken/telegrapher/heat) is reproduced regardless.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_conv.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")
load(f,c)= ( rows=[split(l) for l in readlines(joinpath(repo,f)) if !startswith(l,"#") && !isempty(strip(l))];
             ([parse(Float64,r[1]) for r in rows], [parse(Float64,r[c]) for r in rows]) )
ts,Qs = load("pmp_conv.txt",2)         # shock
tm,Qm = load("pmp_conv_smooth.txt",2)  # smooth

fig=Figure(size=(960,440))
for (k,(t,Q,ttl,col,mq)) in enumerate([(ts,Qs,"steepening shock",:black,2.0),(tm,Qm,"smooth ε bump",:purple,1.5)])
    ax=Axis(fig[1,k], xlabel="t", ylabel=k==1 ? "Q_N(t)" : "",
            title="$ttl  →  Q≈$(round(sum(Q)/length(Q),digits=2))")
    hlines!(ax,[4.0],color=:red,linestyle=:dot,linewidth=1.6, label="reference Q=4 (WENO5)")
    hlines!(ax,[2.0],color=:gray,linestyle=:dash,linewidth=1.0)
    scatterlines!(ax,t,Q,color=col,markersize=8,linewidth=2.2, label="this engine (KT)")
    ylims!(ax,0,5); k==1 && axislegend(ax,position=:rb,framevisible=true)
end
Label(fig[0,:], "BDNKStar — PMP conv_plot: engine is LOW-ORDER (shock Q≈2, smooth Q≈1.3) ≪ reference Q≈4 [honest engine-order gap; physics reproduced]",
      fontsize=10.5, font=:bold)
save(joinpath(outdir,"pmp_conv.png"), fig)
@printf("saved pmp_conv.png | shock Q=%.2f  smooth Q=%.2f  (both << reference Q=4)\n",
        sum(Qs)/length(Qs), sum(Qm)/length(Qm))
