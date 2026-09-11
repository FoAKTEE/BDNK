#=
    shum_frame_analysis.jl — FLUID-FRAME (Shum-style) characteristic-structure
    figure validating the BDNKStar code against Shum et al. arXiv:2509.15303
    eqs.(67–71) (fluid-rest-frame BDNK characteristic speeds + well-posedness).

      Panel A  c0, c±/cs vs q̂  (ŝ=â=1, η̂=ζ̂=0.01): the frame-space characteristic
               structure, with the production frame q̂=0.999 and the well-posed
               boundary q̂=ŝ=1 (Shum eq.71) marked.
      Panel B  c±(r) and c0(r) along the Shum M=1.4 star (production frame), with
               the luminal line c=1 — c₊=√3·cs stays subluminal throughout.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_frame_analysis.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
cs_c = sqrt(sound_speed2(eos, εc))

# production frame
ŝ, â, η̂, ζ̂ = 1.0, 1.0, 0.01, 0.01
q̂_prod = 0.999

# ---- Panel A data: scan q̂ ∈ (0,1] at ŝ=â=1, cs = stellar centre ----------
q̂s = collect(range(0.01, 1.0, length=200))
c0A = Float64[]; cpA = Float64[]; cmA = Float64[]
for q̂ in q̂s
    c0, cp, cm = shum_frame_speeds(ŝ, â, q̂, η̂, ζ̂, cs_c)
    push!(c0A, c0/cs_c); push!(cpA, cp/cs_c); push!(cmA, cm/cs_c)
end

# ---- Panel B data: c0,c±(r) along the star (production frame) -------------
rr = Float64[]; c0r = Float64[]; cpr = Float64[]; cmr = Float64[]
for i in eachindex(star.r)
    ε = star.ε[i]; ε <= 0 && continue
    cs = sqrt(sound_speed2(eos, ε))
    c0, cp, cm = shum_frame_speeds(ŝ, â, q̂_prod, η̂, ζ̂, cs)
    push!(rr, star.r[i]/star.R); push!(c0r, c0); push!(cpr, cp); push!(cmr, cm)
end
cp_max = maximum(cpr)

# ===========================================================================
fig = Figure(size=(1000, 440))
cP, cM, c0c = :crimson, :navy, :teal

# Panel A: frame-space characteristic structure
axA = Axis(fig[1,1], xlabel="frame parameter  q̂   (ŝ=â=1)",
           ylabel="characteristic speed  c / cs",
           title="A  fluid-frame characteristic structure (Shum eq.67-71)")
lines!(axA, q̂s, cpA, color=cP,  linewidth=2.4, label="c₊/cs  (modified sound)")
lines!(axA, q̂s, cmA, color=cM,  linewidth=2.4, label="c₋/cs  (slow mode)")
lines!(axA, q̂s, c0A, color=c0c, linewidth=2.4, label="c0/cs  (diffusive)")
hlines!(axA, [sqrt(3.0)], color=:gray, linestyle=:dot, label="√3")
vlines!(axA, [1.0], color=:black, linestyle=:dash, label="well-posed bd  q̂=ŝ")
vlines!(axA, [q̂_prod], color=:darkorange, linestyle=:dashdot, label="production q̂=0.999")
text!(axA, 0.06, sqrt(3.0)+0.03, text="c₊→√3 cs", color=cP, fontsize=11)
axislegend(axA, position=:cb, framevisible=true, nbanks=2, labelsize=9)

# Panel B: speeds along the star
axB = Axis(fig[1,2], xlabel="r / R   (areal)", ylabel="characteristic speed  c",
           title="B  production-frame c±(r) along the Shum M=1.4 star")
lines!(axB, rr, cpr, color=cP,  linewidth=2.4, label="c₊ = √3 cs(r)")
lines!(axB, rr, cmr, color=cM,  linewidth=2.4, label="c₋(r)")
lines!(axB, rr, c0r, color=c0c, linewidth=2.4, label="c0(r)")
hlines!(axB, [1.0], color=:black, linestyle=:dash, label="luminal c=1")
text!(axB, 0.30, cp_max+0.03,
      text=@sprintf("max c₊ = %.3f < 1  (causal)", cp_max), color=cP, fontsize=11)
ylims!(axB, 0.0, 1.05)
axislegend(axB, position=:rt, framevisible=true, labelsize=9)

Label(fig[0, :],
      "BDNKStar — Shum fluid-frame analysis: c₊=√3 cs, c₋≈0.0183 cs at production " *
      "frame (1,1,0.999); well-posed 0<q̂<ŝ; c₊=√3cs subluminal throughout the M=1.4 star",
      fontsize=12, font=:bold)

save(joinpath(outdir, "shum_frame_analysis.png"), fig)
@printf("saved shum_frame_analysis.png   (centre cs=%.4f, max c₊=%.4f)\n", cs_c, cp_max)
