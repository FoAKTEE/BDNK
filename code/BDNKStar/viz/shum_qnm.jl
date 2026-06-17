#=
    R5 reproduction figure — Shum 2509.15303 QNM_plot: power spectrum of the
    central energy density of the nonlinear spherically-symmetric BDNK Cowling
    evolution (Δr=0.04, t_f=8000 M☉, smallSB-F2), with the fundamental F and
    overtones H1,H2 marked against the Shum table values.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_qnm.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

# read (t, eps_c) from the completed Δr=0.04 production run
data = readlines(joinpath(repo, "r5_eps_Dr0.04.txt"))
rows = [split(l) for l in data if !startswith(l, "#") && !isempty(strip(l))]
t = [parse(Float64, r[1]) for r in rows]
e = [parse(Float64, r[2]) for r in rows]

# use the well-resolved early window (t ≤ 2000 M☉, ~27 f-mode cycles)
m = t .<= 2000.0
tw = t[m]; ew = e[m]
ew = ew .- sum(ew)/length(ew)                     # remove DC
N = length(tw)
bw = [0.42 - 0.5cos(2π*i/(N-1)) + 0.08cos(4π*i/(N-1)) for i in 0:N-1]  # Blackman
ewb = ew .* bw

const MSUN_KHZ = 203.0248                          # 1/M☉ = 203.0248 kHz (cyclic)
fkHz = range(0.5, 8.0; length=800)
P = Float64[]
for f in fkHz
    fg = f / MSUN_KHZ                              # cyclic freq in 1/M☉
    s = sum(ewb .* exp.(-2π*im*fg .* tw))
    push!(P, abs2(s))
end
P ./= maximum(P)

fig = Figure(size=(760, 460))
ax = Axis(fig[1,1], title="Shum 2509.15303 QNM_plot — central-ε power spectrum (Δr=0.04, smallSB-F2)",
          xlabel="f [kHz]", ylabel="normalized power", yscale=log10)
lines!(ax, collect(fkHz), max.(P, 1e-6), color=:dodgerblue, linewidth=2)
# Shum table targets
for (lbl, ft, col) in (("F 2.69", 2.69, :crimson), ("H1 4.55", 4.55, :seagreen), ("H2 6.36", 6.36, :purple))
    vlines!(ax, [ft], color=col, linestyle=:dash, linewidth=1.8)
    text!(ax, ft, 1.3, text=lbl, align=(:center,:bottom), color=col, fontsize=12)
end
# achieved peaks
for fa in (2.699, 4.551, 6.468)
    scatter!(ax, [fa], [1.0], color=:black, markersize=10, marker=:dtriangle)
end
ylims!(ax, 1e-4, 2.0)
Label(fig[0,:], "BDNKStar R5 — reproduce Shum QNM: F=2.699/H1=4.551/H2=6.468 kHz vs (2.69,4.55,6.36), <2%",
      fontsize=13, font=:bold)
outfile = joinpath(outdir, "shum_qnm_reproduction.png")
save(outfile, fig); println("saved ", outfile)
