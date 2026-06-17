#=
    Caballero–Yunes heat-conduction stability criterion figure (2506.09149):
      adiabatic c_s² vs fixed-baryon c_n², and Δ = c_s²−c_n². Stability needs
      Δ ≥ 0. A cold barotrope is marginal (Δ=0); the Γ-law ideal gas violates it
      (Δ = −(Γ−1)/(1+Γϵ) < 0) ⇒ conditional heat-conduction instability.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/heat_criterion.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
ρ = 1e-3
ϵs = range(0.02, 3.0; length=200)

fig = Figure(size=(1150, 460))
axA = Axis(fig[1,1], title="A. Sound speeds (ideal gas)", xlabel="ϵ (specific internal energy)", ylabel="c²")
for (Γ, col) in ((4/3,:dodgerblue),(5/3,:seagreen),(2.0,:crimson))
    gas = IdealGas(Γ)
    lines!(axA, ϵs, [sound_speed2(gas,ρ,ϵ) for ϵ in ϵs], color=col, linewidth=2.5, label="c_s² Γ=$(round(Γ,digits=2))")
    lines!(axA, ϵs, [cn2(gas,ρ,ϵ) for ϵ in ϵs], color=col, linewidth=2, linestyle=:dash)
end
axislegend(axA, position=:rb, framevisible=true, labelsize=10)
text!(axA, 0.05, 0.02, text="dashed = c_n² (fixed baryon)", align=(:left,:bottom), fontsize=11)

axB = Axis(fig[1,2], title="B. CY criterion  Δ = c_s² − c_n²  (stable ⇔ Δ≥0)",
           xlabel="ϵ", ylabel="c_s² − c_n²")
hlines!(axB, [0.0], color=:black, linewidth=2)
band!(axB, [0.02,3.0], -0.6, 0.0, color=(:red,0.07))
text!(axB, 1.5, -0.05, text="Δ<0: heat-conduction UNSTABLE", align=(:center,:top), color=:red, fontsize=12)
for (Γ, col) in ((4/3,:dodgerblue),(5/3,:seagreen),(2.0,:crimson))
    gas = IdealGas(Γ)
    lines!(axB, ϵs, [sound_speed2(gas,ρ,ϵ)-cn2(gas,ρ,ϵ) for ϵ in ϵs], color=col, linewidth=2.5,
           label="ideal Γ=$(round(Γ,digits=2)): −(Γ−1)/(1+Γϵ)")
end
# cold barotrope marginal: c_s²=c_n² ⇒ Δ=0
lines!(axB, ϵs, zeros(length(ϵs)), color=:gray, linestyle=:dot, linewidth=2, label="cold barotrope (Δ=0, marginal)")
axislegend(axB, position=:rb, framevisible=true, labelsize=10)
ylims!(axB, -0.5, 0.1)

Label(fig[0, :], "BDNKStar — Caballero–Yunes heat-conduction stability criterion c_s²−c_n²≥0 (ideal gas violates it)",
      fontsize=14, font=:bold)
outfile = joinpath(outdir, "heat_criterion.png")
save(outfile, fig); println("saved ", outfile)
