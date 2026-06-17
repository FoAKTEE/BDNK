#=
    PMP 2209.09265 conv_plot — PDE self-convergence FACTOR Q_N(t) of the BDNK shock.
    Q_N = ||u_513−u_1025|| / ||u_1025−u_2049||  (the convergence FACTOR, = 2^order).
    The paper's conv_plot red line is Q_N=4 = 2² → order 2 (shock-limited WENO5);
    the ODE/RK4 tests report Q_N≈16=2⁴ (cf. Bjorken Q→16, reproduced separately).
    Our shock gives Q_N≈4.0 (order 2.0) — REPRODUCES the PMP conv_plot.  Earlier I
    plotted the order (2.0) and misread Q_N=4 as 4th order; Q_N=4 is the factor.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_conv.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

# shock: factor = d_lo/d_hi (cols 3,4 of pmp_conv.txt)
rs=[split(l) for l in readlines(joinpath(repo,"pmp_conv.txt")) if !startswith(l,"#") && !isempty(strip(l))]
ts=[parse(Float64,r[1]) for r in rs]; Qs=[parse(Float64,r[3])/parse(Float64,r[4]) for r in rs]
# smooth (for context): factor = 2^order
rm=[split(l) for l in readlines(joinpath(repo,"pmp_conv_smooth.txt")) if !startswith(l,"#") && !isempty(strip(l))]
tm=[parse(Float64,r[1]) for r in rm]; Qm=[2.0^parse(Float64,r[2]) for r in rm]

fig=Figure(size=(820,500))
ax=Axis(fig[1,1], xlabel="t", ylabel="Q_N(t)  (convergence factor)",
        title="PMP conv_plot — BDNK shock PDE convergence FACTOR (N=513,1025,2049)")
hlines!(ax,[4.0],color=:red,linestyle=:dot,linewidth=1.8, label="paper Q_N=4 (=2², order 2)")
hlines!(ax,[16.0],color=:gray,linestyle=:dash,linewidth=1.0, label="RK4/ODE Q_N=16 (cf. Bjorken Q→16)")
scatterlines!(ax,ts,Qs,color=:black,markersize=8,linewidth=2.2, label="this engine: shock Q_N≈$(round(sum(Qs)/length(Qs),digits=2))")
scatterlines!(ax,tm,Qm,color=:purple,markersize=6,linewidth=1.6, linestyle=:dot, label="this engine: smooth Q_N≈$(round(sum(Qm)/length(Qm),digits=2))")
ylims!(ax,0,18); xlims!(ax,0,maximum(ts)+0.5)
axislegend(ax,position=:rc,framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce PMP conv_plot: shock convergence FACTOR Q_N≈4.0 matches the paper (Q_N=4 = 2² = order 2, shock-limited WENO5)",
      fontsize=10.5, font=:bold)
save(joinpath(outdir,"pmp_conv.png"), fig)
@printf("saved pmp_conv.png | shock Q_N=%.3f (paper 4.0); smooth Q_N=%.2f; [Q_N=2^order]\n",
        sum(Qs)/length(Qs), sum(Qm)/length(Qm))
