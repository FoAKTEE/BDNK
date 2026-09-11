#=
    IS-vs-BDNK bulk-viscous contrast figure (STAGE 2, node s2.is_contrast) —
    Chabanov–Rezzolla arXiv:2311.13027.

      A  the causal limiter: effective viscous sound speed c_s'^2 vs ζ/τ_Π,
         naive (eq:limit_sound) vs limited (eq:causal_limit clamps to c_max²).
         The naive branch goes superluminal; the limiter holds it at c_max².
      B  inverse-Reynolds-number band R⁻¹∈[σp/(e+p),(e−p)/(e+p)] vs ρ for the
         ideal EOS p=κρ², with the paper's ρ_c=1.28e-3 point and its verbatim
         targets R⁻¹_min≃−0.1, R⁻¹_max≃0.8 marked.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/is_contrast.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# ---- Panel A data: causal limiter -----------------------------------------
cs2_eq, ζ, ρc, hp, cmax2 = 0.2, 1.0e-3, 0.0011, 1.23, 0.9
ratios = 10 .^ range(-5, 1; length=400)             # ζ/τ_Π
cs2_naive = [cs2_viscous(cs2_eq, ζ, ζ/r, ρc, hp) for r in ratios]
cs2_lim   = [apply_causality_fix(cs2_eq, ζ, ζ/r, ρc, hp; cmax2=cmax2)[2] for r in ratios]

# ---- Panel B data: inverse-Reynolds band ----------------------------------
Γ, κ, σ = 2.0, 100.0, -0.9
ρs   = 10 .^ range(-5, log10(8e-3); length=400)
pB(ρ) = κ * ρ^Γ
eB(ρ) = ρ + pB(ρ) / (Γ - 1)
Rmin = [reynolds_inv_min(σ, pB(ρ), eB(ρ)) for ρ in ρs]
Rmax = [reynolds_inv_max(pB(ρ), eB(ρ))    for ρ in ρs]
ρ_paper = 1.28e-3
Rmin_p  = reynolds_inv_min(σ, pB(ρ_paper), eB(ρ_paper))
Rmax_p  = reynolds_inv_max(pB(ρ_paper), eB(ρ_paper))

fig = Figure(size=(960, 420))

axA = Axis(fig[1,1], xscale=log10, xlabel="ζ / τ_Π", ylabel="c_s'^2",
           title="A  IS causal limiter (eq:limit_sound → eq:causal_limit)")
lines!(axA, ratios, cs2_naive, color=:crimson, label="naive c_s'^2")
lines!(axA, ratios, cs2_lim,   color=:seagreen, linewidth=2.5, label="limited (≤ c_max²)")
hlines!(axA, [1.0], color=:black, linestyle=:dash, label="light cone c²=1")
hlines!(axA, [cmax2], color=:gray, linestyle=:dot, label="c_max² = 0.9")
ylims!(axA, 0, 2.0)
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xscale=log10, xlabel="ρ  [M_⊙^{-2}]", ylabel="R⁻¹ = Π/(p+e)",
           title="B  inverse-Reynolds band (Chabanov–Rezzolla §3)")
band!(axB, ρs, Rmin, Rmax, color=(:steelblue, 0.25))
lines!(axB, ρs, Rmax, color=:steelblue, label="R⁻¹_max = (e−p)/(e+p)")
lines!(axB, ρs, Rmin, color=:orange,    label="R⁻¹_min = σp/(e+p)")
scatter!(axB, [ρ_paper], [Rmax_p], color=:steelblue, markersize=12)
scatter!(axB, [ρ_paper], [Rmin_p], color=:orange,    markersize=12)
# right-aligned, anchored just left of the ρ_c marker so the paper-target
# parentheticals stay inside the panel frame
text!(axB, ρ_paper*0.85, Rmax_p; text="ρ_c: R⁻¹_max=$(round(Rmax_p,digits=3)) (≃0.8)",
      align=(:right,:center), fontsize=11)
text!(axB, ρ_paper*0.85, Rmin_p; text="R⁻¹_min=$(round(Rmin_p,digits=3)) (≃−0.1)",
      align=(:right,:center), fontsize=11)
axislegend(axB, position=:lb, framevisible=true)

Label(fig[0, :], "BDNKStar — IS-vs-BDNK bulk-viscous contrast (s2.is_contrast): " *
      "MIS causal limiter + reproduced Reynolds limits (2311.13027)",
      fontsize=14, font=:bold)

save(joinpath(outdir, "is_contrast.png"), fig)
println("saved is_contrast.png | Rmin_paper=", Rmin_p, " Rmax_paper=", Rmax_p)
