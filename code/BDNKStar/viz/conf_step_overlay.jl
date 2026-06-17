#=
    Pandya 2201.12317 — STEP (Riemann) ID, Julia engine vs reference C code.
    Sharp discontinuity εL=1 / εR=0.1, u=0; conformal PMP-luminal frame; N=129.
    Completes the code-validation triad (Gaussian, SMOOTH_SHOCK, STEP).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_step_overlay.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

eta0 = 10.0^0.25/(3π)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)
s    = init_step(fr; N=129, xmin=-200.0, xmax=200.0, εL=1.0, εR=0.1, cfl=0.1)
rows = [1,16,31,51,102]; nstep=[10*(r-1) for r in rows]
xj = copy(s.x); juls = Vector{Vector{Float64}}(); done=0
push!(juls, copy(energy_density(s)))
for k in 2:length(nstep); evolve!(s, nstep[k]-done); global done=nstep[k]; push!(juls, copy(energy_density(s))); end

clines=[l for l in readlines(joinpath(repo,"conf_cref_step_eps.txt")) if !isempty(strip(l))]
dx=400.0/128; xc=[-200.0+i*dx for i in 0:128]
cref=[[parse(Float64,t) for t in split(clines[r])] for r in rows]
times=[0.0,46.875,93.75,156.25,315.625]; cols=[:black,:crimson,:seagreen,:dodgerblue,:darkorange]

fig=Figure(size=(900,560))
ax=Axis(fig[1,1], xlabel="x", ylabel="ε(x)",
        title="Pandya STEP/Riemann (εL=1, εR=0.1): Julia engine (lines) vs reference C code (dots)")
maxdev=0.0
for k in 1:5
    lines!(ax, xj, juls[k], color=cols[k], linewidth=2.0, label="t=$(times[k])")
    scatter!(ax, xc[1:3:end], cref[k][1:3:end], color=cols[k], markersize=6)
    m=abs.(xj).<80
    dev=maximum(abs.(juls[k][m].-cref[k][m])./max.(cref[k][m],1e-3)); global maxdev=max(maxdev,dev)
    @printf("  t=%-9.3f Julia[min,max]=[%.4f,%.4f] C[min,max]=[%.4f,%.4f] dev=%.3g\n",
            times[k], minimum(juls[k]), maximum(juls[k]), minimum(cref[k]), maximum(cref[k]), dev)
end
xlims!(ax,-80,80)
axislegend(ax, "snapshot time", position=:rc, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Pandya STEP/Riemann vs reference C code (max rel-dev = $(round(maxdev*100,sigdigits=2))%)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"conf_step_overlay.png"), fig)
println("saved conf_step_overlay.png | overall max rel-dev (|x|<80) = ", round(maxdev*100,sigdigits=3), "%")
