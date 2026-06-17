#=
    PMP 2209.09265 shockwave_plot (fig:shockwave_profile) — steady planar BDNK
    shock structure ε (solid), v (dash-dot), n (dotted) vs x.  Left {1,0.8,0.1},
    Γ=4/3, m=0.1, V̂=2/15, σ̂=0, τ̂=1.5; relaxed steady state (black) vs the initial
    smooth shock (green), showing the BDNK shock sharpens to its steady width.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_shockwave.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

rows=[split(l) for l in readlines(joinpath(repo,"pmp_shockwave.txt")) if !startswith(l,"#") && !isempty(strip(l))]
x =[parse(Float64,r[1]) for r in rows]
e0=[parse(Float64,r[2]) for r in rows]; v0=[parse(Float64,r[3]) for r in rows]; n0=[parse(Float64,r[4]) for r in rows]
ef=[parse(Float64,r[5]) for r in rows]; vf=[parse(Float64,r[6]) for r in rows]; nf=[parse(Float64,r[7]) for r in rows]

fig=Figure(size=(720,500))
ax=Axis(fig[1,1], xlabel="x", ylabel="f(x)",
        title="PMP shockwave_plot — steady BDNK shock {1,0.8,0.1}_L (Γ=4/3, V̂=2/15, τ̂=1.5)")
# relaxed steady (black) and initial smooth (green)
lines!(ax, x, ef, color=:black,     linewidth=2.0, label="f=ε")
lines!(ax, x, vf, color=:black,     linewidth=2.0, linestyle=:dashdot, label="f=v")
lines!(ax, x, nf, color=:black,     linewidth=2.0, linestyle=:dot,     label="f=n")
lines!(ax, x, e0, color=:seagreen,  linewidth=1.6)
lines!(ax, x, v0, color=:seagreen,  linewidth=1.6, linestyle=:dashdot)
lines!(ax, x, n0, color=:seagreen,  linewidth=1.6, linestyle=:dot)
xlims!(ax, -6, 6)
axislegend(ax, position=:lt, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce PMP shockwave_plot: steady BDNK shock ε/v/n (black=relaxed, green=initial); εL=1→εR=$(round(ef[argmax(ef)],digits=2)), vL=0.8→vR=$(round(vf[end],digits=3))",
      fontsize=11, font=:bold)
save(joinpath(outdir,"pmp_shockwave.png"), fig)
@printf("saved pmp_shockwave.png | εL=%.3f εR=%.3f  vL=%.3f vR=%.3f  nL=%.3f nR=%.3f\n",
        ef[1], maximum(ef), vf[1], vf[end], nf[1], nf[end])
