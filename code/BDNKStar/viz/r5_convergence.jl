#=
    R5 QNM convergence (Shum 2509.15303, convergence/error_fit): the extracted
    QNM frequencies vs resolution Δr, Richardson-extrapolated to the continuum,
    vs the Shum table targets. Frequencies converge to <1.4%; the f-mode DECAY
    rate is weakly damped and resolution-sensitive (Shum resolves it with a
    4-resolution Δr=0.0020–0.0032 ladder — a [HOLE] here at Δr∈{0.04,0.01}).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/r5_convergence.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

Δr   = [0.04, 0.01]
F    = [2.699, 2.700];  H1 = [4.551, 4.589];  H2 = [6.468, 6.453]
Fc, H1c, H2c = 2.700, 4.602, 6.448          # 2-pt Richardson continuum
Ft, H1t, H2t = 2.69, 4.55, 6.36             # Shum targets
τl   = [0.00209, 0.000136]; τnl = [0.0, 0.00046]   # decay (scattered)

fig = Figure(size=(1150, 450))
axF = Axis(fig[1,1], title="A. QNM frequencies → continuum (Shum targets dashed)",
           xlabel="Δr [M⊙]", ylabel="f [kHz]", xscale=log10)
for (y, yc, yt, c, lab) in ((F,Fc,Ft,:crimson,"F"),(H1,H1c,H1t,:seagreen,"H1"),(H2,H2c,H2t,:dodgerblue,"H2"))
    scatterlines!(axF, Δr, y, color=c, markersize=11, label="$lab (Δr runs)")
    scatter!(axF, [3e-3], [yc], color=c, marker=:star5, markersize=16)       # continuum
    hlines!(axF, [yt], color=c, linestyle=:dash, linewidth=1.5)
end
text!(axF, 3e-3, 6.6, text="★ = continuum", align=(:center,:bottom), fontsize=11)
axislegend(axF, position=:rc, framevisible=true, labelsize=10)

axT = Axis(fig[1,2], title="B. f-mode decay 1/τ vs Δr — resolution-sensitive [HOLE]",
           xlabel="Δr [M⊙]", ylabel="1/τ [M⊙⁻¹]", xscale=log10, yscale=log10)
scatterlines!(axT, Δr, max.(τl, 1e-5), color=:purple, markersize=11, label="linear log-max")
scatterlines!(axT, Δr, max.(τnl, 1e-5), color=:orange, markersize=11, label="nonlinear fit")
hlines!(axT, [0.0011], color=:black, linestyle=:dash, linewidth=1.8, label="Shum continuum 0.0011")
band!(axT, [8e-3, 5e-2], 0.0016, 0.0019, color=(:gray,0.2))
text!(axT, 1.2e-2, 0.0017, text="Shum per-Δr band", align=(:left,:bottom), fontsize=10)
axislegend(axT, position=:lt, framevisible=true, labelsize=10)

Label(fig[0,:], "BDNKStar R5 convergence — frequencies → continuum <1.4% of Shum; decay rate needs a finer Δr ladder",
      fontsize=13, font=:bold)
outfile = joinpath(outdir, "r5_convergence_reproduction.png")
save(outfile, fig); println("saved ", outfile)
