#=
    Shum R5 decay-rate extraction sensitivity — why error_fit / decay-continuum are
    blocked.  The f-mode FREQUENCY ω is robust (resolution-robust eigenvalue, f<2%),
    but the DECAY RATE 1/τ is method/window-sensitive: the damped-sinusoid fit on
    the Δr=0.04 series gives 1/τ ∈ [0.0004,0.0015] across fit windows, and the
    log-of-maxima gives 0.00209 — together BRACKETING the paper's continuum (0.0011)
    and per-Δr band (0.0016–0.0019) but not robustly pinned.  A clean error_fit
    needs both the fine Δr (which blow up) AND a robust extractor.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_decay_sensitivity.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# damped-sinusoid 1/τ across fit windows (measured), + the log-of-maxima value
windows = ["[2000,5000]","[2500,6000]","[3000,7000]","[2000,7000]","[1500,5500]","[3000,8000]"]
itau_ds = [0.0006516,0.001091,0.001464,0.0007687,0.0004268,0.001482]
itau_lm = 0.00209   # log-of-maxima (per-Δr method)

fig=Figure(size=(820,500))
ax=Axis(fig[1,1], xlabel="extracted 1/τ  [1/M_⊙]", ylabel="",
        yticks=(1:length(windows), ["damped-fit "*w for w in windows]),
        title="Shum R5 decay-rate extraction sensitivity (Δr=0.04) — frequency robust, decay rate is not")
# paper reference bands
vspan!(ax, 0.0016, 0.0019, color=(:orange,0.25))
vlines!(ax, [0.0011], color=:purple, linestyle=:dash, linewidth=2, label="paper continuum 0.0011")
text!(ax, 0.00175, 6.4, text="paper per-Δr\n0.0016–0.0019", fontsize=10, color=:darkorange)
scatter!(ax, itau_ds, 1:length(windows), color=:steelblue, markersize=13, label="damped-sinusoid fit (window-sensitive)")
scatter!(ax, [itau_lm], [length(windows)+0.7], color=:black, marker=:star5, markersize=16, label="log-of-maxima (per-Δr) 0.00209")
xlims!(ax, 0, 0.0024)
axislegend(ax, position=:rb, framevisible=true)
Label(fig[0,:], "BDNKStar — Shum decay-rate extraction sensitivity: 1/τ spans 0.0004–0.0021 (method/window) bracketing paper values; error_fit needs fine Δr + robust extractor",
      fontsize=9.5, font=:bold)
save(joinpath(outdir,"shum_decay_sensitivity.png"), fig)
println("saved shum_decay_sensitivity.png | damped-fit 1/τ range [",
        minimum(itau_ds),",",maximum(itau_ds),"]; log-max 0.00209; paper continuum 0.0011, per-Δr 0.0016-0.0019")
