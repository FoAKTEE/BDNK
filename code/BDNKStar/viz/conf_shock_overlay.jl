#=
    Pandya 2201.12317 steady planar shock — Julia engine vs reference C code.
    ε(x) at snapshot rows {1,16,31,51,102} for the SMOOTH_SHOCK ID (εL=1, vL=0.8).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_shock_overlay.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

jl = [split(l) for l in readlines(joinpath(repo,"conf_shock_overlay.txt")) if !startswith(l,"#") && !isempty(strip(l))]
xj   = [parse(Float64,r[1]) for r in jl]
juls = [[parse(Float64,r[k]) for r in jl] for k in 2:6]

clines = [l for l in readlines(joinpath(repo,"conf_cref_shock_eps.txt")) if !isempty(strip(l))]
rows = [1,16,31,51,102]; dx = 400.0/128; xc = [-200.0 + i*dx for i in 0:128]
cref_cols = [[parse(Float64,t) for t in split(clines[r])] for r in rows]
times = [0.0, 46.875, 93.75, 156.25, 315.625]
cols  = [:black, :crimson, :seagreen, :dodgerblue, :darkorange]

fig = Figure(size=(900, 560))
ax = Axis(fig[1,1], xlabel="x", ylabel="ε(x)",
          title="Pandya steady shock (εL=1, vL=0.8): Julia engine (lines) vs reference C code (dots)")
maxdev = 0.0
for k in 1:5
    lines!(ax, xj, juls[k], color=cols[k], linewidth=2.0, label="t=$(times[k])")
    scatter!(ax, xc[1:3:end], cref_cols[k][1:3:end], color=cols[k], markersize=6)
    m = abs.(xj) .< 80
    dev = maximum(abs.(juls[k][m] .- cref_cols[k][m]) ./ max.(cref_cols[k][m], 1e-3))
    global maxdev = max(maxdev, dev)
    @printf("  t=%-9.3f  Julia εmax=%.5f  C εmax=%.5f  max-rel-dev(|x|<80)=%.3g\n",
            times[k], maximum(juls[k]), maximum(cref_cols[k]), dev)
end
xlims!(ax, -80, 80)
axislegend(ax, "snapshot time", position=:rc, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Pandya steady shock vs reference C code (εR=4.4074; max rel-dev = $(round(maxdev*100,sigdigits=2))%)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"conf_shock_overlay.png"), fig)
println("saved conf_shock_overlay.png | overall max rel-dev (|x|<80) = ", round(maxdev*100,sigdigits=3), "%")
