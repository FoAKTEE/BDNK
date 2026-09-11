#=
    STAGE-3 result: full-frame (causal) BDNK shear viscosity dissipates the stellar
    mode energy on the boundary-conforming (r,θ) FULL star (surface on r=R — the
    geometry that makes the velocity recovery stable, where the staircased Cartesian
    surface blew up).

      A  quadratic mode energy E(t)/E₀ for a range of η̂: the inviscid scheme has a
         slow grid-scale GROWTH; shear viscosity ν_mom∇²δS (tied to η̂) turns it into
         monotonically stronger DECAY — genuine dissipation.
      B  E_end/E₀ vs η̂: monotone, crossing 1 (growth→decay). The faithful viscous
         diagnostic is the ENERGY, not amplitude/quadrupole decay (which is masked by
         Kreiss–Oliger dissipation and by dephasing of the non-eigenmode seed).

    KEY POINT (hard-won): of the three BDNK frame channels, only the shear-viscous
    MOMENTUM force ν_mom∇²δS dissipates; the conduction/heat term κ_Q c_s²∂ε is
    REACTIVE (sets frequency/stability, not damping), and the relaxation τ_R is the
    causal-frame stiffness knob.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/sphbdnk_viscous.jl
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
s = build_sphstar(eos, εc; Nr=64, Nθ=12)

function run(η̂)
    e  = setup_sphbdnk(s; η̂=η̂, σ_ko=0.02)          # physical frame: all transport ∝ η̂
    st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, e, 3; A=1e-3)
    ts, _, en = evolve_sphbdnk!(st, e; dt=0.005, nsteps=round(Int, 600/0.005), sample=40)
    ts, en ./ en[1]
end

ηs = [0.0, 0.02, 0.04, 0.06, 0.08]
series = [run(η̂) for η̂ in ηs]
ratios = [ser[2][end] for ser in series]
println("η̂:        ", ηs)
println("E_end/E0:  ", round.(ratios, digits=4))

fig = Figure(size=(980, 430))
cols = [:black, :navy, :teal, :darkorange, :crimson]

axA = Axis(fig[1,1], xlabel="t  [M⊙ geometric]", ylabel="mode energy  E(t)/E₀",
           yscale=log10, title="A  BDNK shear viscosity dissipates the mode energy")
for (k, ser) in enumerate(series)
    lines!(axA, ser[1], max.(ser[2], 1e-3), color=cols[k], label=@sprintf("η̂=%.2f", ηs[k]))
end
hlines!(axA, [1.0], color=:gray, linestyle=:dash)
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xlabel="shear-viscosity knob  η̂", ylabel="E_end/E₀",
           title="B  energy dissipation increases monotonically with η̂")
scatter!(axB, ηs, ratios, color=:crimson, markersize=13)
lines!(axB, ηs, ratios, color=:crimson)
hlines!(axB, [1.0], color=:gray, linestyle=:dash, label="growth ↑ / decay ↓")
axislegend(axB, position=:rt, framevisible=true)

Label(fig[0, :], "BDNKStar — STAGE 3: full-frame causal BDNK on the boundary-conforming (r,θ) full star — " *
      "shear viscosity ν_mom∇²δS dissipates the mode energy ∝ η̂ (energy diagnostic)",
      fontsize=12, font=:bold)

save(joinpath(outdir, "sphbdnk_viscous.png"), fig)
println("saved sphbdnk_viscous.png")
