#=
    Pandya 2201.12317 — Julia conformal-BDNK engine vs the published reference C
    code (1D_conformal_bdnk).  Gaussian clump ε=exp(−x²/25)+0.1 evolved by BOTH
    codes; overlay ε(x) at the reference snapshot rows {1,16,31,51,102}.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_overlay.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")
cref = joinpath(repo, "conf_cref_eps.txt")   # committed copy of the reference C output

# Julia engine output
jl = [split(l) for l in readlines(joinpath(repo,"conf_gaussian_overlay.txt")) if !startswith(l,"#") && !isempty(strip(l))]
xj   = [parse(Float64,r[1]) for r in jl]
juls = [[parse(Float64,r[k]) for r in jl] for k in 2:6]    # 5 time columns

# C reference: each line = ε(x) at one saved time; pick rows 1,16,31,51,102
clines = [l for l in readlines(cref) if !isempty(strip(l))]
rows = [1,16,31,51,102]
dx = 400.0/128; xc = [-200.0 + i*dx for i in 0:128]
cref_cols = [[parse(Float64,t) for t in split(clines[r])] for r in rows]
times = [0.0, 46.875, 93.75, 156.25, 315.625]
cols  = [:black, :crimson, :seagreen, :dodgerblue, :darkorange]

fig = Figure(size=(900, 560))
ax = Axis(fig[1,1], xlabel="x", ylabel="ε(x)",
          title="Pandya conformal-BDNK: Julia engine (lines) vs reference C code (dots)")
maxdev = 0.0
for k in 1:5
    lines!(ax, xj, juls[k], color=cols[k], linewidth=2.0, label="t=$(times[k])")
    scatter!(ax, xc[1:4:end], cref_cols[k][1:4:end], color=cols[k], markersize=6)
    # relative deviation near the pulse (|x|<60, where ε is well above the 0.1 floor)
    m = abs.(xj) .< 60
    dev = maximum(abs.(juls[k][m] .- cref_cols[k][m]) ./ max.(cref_cols[k][m], 1e-3))
    global maxdev = max(maxdev, dev)
    @printf("  t=%-9.3f  Julia εmax=%.5f  C εmax=%.5f  max-rel-dev(|x|<60)=%.3g\n",
            times[k], maximum(juls[k]), maximum(cref_cols[k]), dev)
end
xlims!(ax, -120, 120)
axislegend(ax, "snapshot time", position=:rt, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Pandya 2201.12317 conformal evolution vs reference C code (max rel-dev = $(round(maxdev*100,sigdigits=2))%)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"conf_overlay.png"), fig)
println("saved conf_overlay.png | overall max rel-dev (|x|<60) = ", round(maxdev*100,sigdigits=3), "%")
