#=
    Shum 2509.15303 stable_evol_comparing_tau — the equilibrium energy-density
    profile ε·M_⊙² vs r/M_⊙ (isotropic) of the M_T=1.4 BDNK star.  In the paper
    the four frame cases (smallSB-F2, medS-F2, highB-F9, medSB-F9) all evolve from
    — and remain on — THIS equilibrium (stability); they differ only at the <0.1%
    level shown in the insets.  All four share the same TOV background, so the
    dominant curve is the equilibrium ε(r) computed here directly from the Shum
    TOV solution (ρ0c=0.00128 → M_T=1.4, εc=ρ0c+κρ0c²=0.00144384).

    SCOPE: this reproduces the t=0 / main equilibrium curve (which the BDNK
    evolution preserves).  The sub-0.1% per-frame inset splitting requires the
    four high-precision evolutions; flagged PRELIMINARY for that part.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_profile.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "shum_core.jl"))   # ShumPolytrope, solve_tov, areal_to_isotropic
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

κ=100.0; eos=ShumPolytrope(κ); ρ0c=0.00128; εc=ρ0c+κ*ρ0c^2
star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
iso  = areal_to_isotropic(star)
r = iso.r; ε = iso.ε                          # ε·M_⊙² (numerically ε, M_⊙=1)
# extend with the vacuum atmosphere out to r=18 (ε=0) to match the paper's x-range
rfull = vcat(r, [iso.rstar, iso.rstar+1e-9, 18.0]); εfull = vcat(ε, [0.0,0.0,0.0])

fig = Figure(size=(820, 560))
ax = Axis(fig[1,1], xlabel="r / M_⊙", ylabel="ε M_⊙²",
          title="Shum stable_evol_comparing_tau — equilibrium ε(r), M_T=$(round(star.M,digits=3)) M_⊙")
scatter!(ax, rfull, εfull, color=:black, markersize=5, label="t=0 equilibrium (4 frames overlap)")
lines!(ax, rfull, εfull, color=(:orange,0.8), linewidth=1.6)
xlims!(ax, 0, 18); ylims!(ax, -0.00005, 0.00150)
axislegend(ax, position=(0.98, 0.60), framevisible=true)

# center-zoom inset (top) — matches original's 0.0014428–0.0014434 near r~0.05–0.10
axc = Axis(fig[1,1], width=Relative(0.34), height=Relative(0.28), halign=0.40, valign=0.95,
           xticklabelsize=9, yticklabelsize=8, title="center zoom", titlesize=9)
mask = r .<= 0.14
scatter!(axc, r[mask], ε[mask], color=:black, markersize=5)
lines!(axc, r[mask], ε[mask], color=(:orange,0.8), linewidth=1.5)
# surface-zoom inset (bottom-right) — ε→0 near r≈r_star
axs = Axis(fig[1,1], width=Relative(0.30), height=Relative(0.26), halign=0.97, valign=0.32,
           xticklabelsize=9, yticklabelsize=8, title="surface zoom (×1e-6)", titlesize=9)
sm = (r .>= iso.rstar-0.4) .& (r .<= iso.rstar)
scatter!(axs, r[sm], ε[sm].*1e6, color=:black, markersize=5)
lines!(axs, r[sm], ε[sm].*1e6, color=(:orange,0.8), linewidth=1.5)

Label(fig[0,:], "BDNKStar — reproduce Shum stable_evol_comparing_tau equilibrium profile (εc=$(round(εc,sigdigits=6)), r_star=$(round(iso.rstar,digits=2)) M_⊙)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"shum_profile.png"), fig)
println("saved shum_profile.png | M_T=", round(star.M,sigdigits=6),
        "  εc=", round(εc,sigdigits=6), " (paper ~0.00144)  r_star(iso)=", round(iso.rstar,sigdigits=4),
        " M_⊙ (paper ~8)")
