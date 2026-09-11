#=
    STAGE-3 viscous result: BDNK shear viscosity damps the non-radial f-mode.

      A  excess mode damping rate 1/τ (above the inviscid numerical floor) vs the
         dimensionless shear-viscosity knob η̂, from full 3+1D Cowling evolutions —
         a clean linear law 1/τ ∝ η̂ (the BDNK/Navier–Stokes viscous scaling).
      B  the ℓ=2 quadrupole time-series: inviscid (long-lived) vs viscous (damped).

    The HMNS-viscosity effect — viscosity damping the stellar oscillation — in 3D.

    SCOPE: this evolves the LEADING-order (Navier–Stokes-limit) BDNK shear+bulk
    stress as a flat divergence ∂_jπ^{ij}, WITHOUT the frame relaxation times
    (τ_ε,τ_P,τ_Q) or the energy/heat sectors — valid in the long-wavelength regime
    where BDNK≈NS. The BDNK frame's role is UV causality (checked separately by
    the Causality monitor, not enforced by this parabolic evolution). The slope
    is a geometric, convention-dependent rate (η₀=η̂·ε₀), not a calibrated η/s.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bdnk_viscous_damping.jl
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
s = build_star3d(eos, εc; N=32, Lfac=1.2)

function run(η̂)
    e = setup_evo3d(s; σ_ko=0.004, η̂=η̂)
    st = EvolState(s.grid.N); seed_l2!(st, e; A=1e-3)
    ts, q2, _ = evolve3d!(st, e; dt=0.18*s.grid.dx, nsteps=5200, sample=4)
    ts, q2, damping_rate(ts, q2)
end

ηs = [0.0, 0.04, 0.08, 0.12, 0.16]
series = [run(η̂) for η̂ in ηs]
d0 = series[1][3]
excess = [s3[3] - d0 for s3 in series]
# `damping_rate` returns NaN when a run holds <3 peaks (too short, or overdamped
# and no longer ringing). All five points here are measurable at these settings
# (11–13 peaks), but guard anyway: a single NaN would poison both the scatter and
# the through-origin fit below without any visible error.
meas = isfinite.(excess)
all(meas) || println("  NOTE: rate unmeasurable at η̂ = ", ηs[.!meas], " — omitted")
println("η̂: ", ηs[meas]); println("excess 1/τ: ", excess[meas])

fig = Figure(size=(960, 430))

axA = Axis(fig[1,1], xlabel="shear-viscosity knob  η̂", ylabel="excess damping  1/τ  [geom]",
           title="A  BDNK viscous damping of the f-mode: 1/τ ∝ η̂")
scatter!(axA, ηs[meas], excess[meas], color=:crimson, markersize=12)
# linear-through-origin fit slope (measurable points only)
slope = sum(ηs[meas] .* excess[meas]) / sum(ηs[meas] .^ 2)
lines!(axA, [0, maximum(ηs)], [0, slope*maximum(ηs)], color=:gray, linestyle=:dash,
       label=@sprintf("linear fit, 1/τ = %.3f η̂", slope))
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xlabel="t  [M⊙ geometric]", ylabel="ℓ=2 quadrupole of δε",
           title="B  ℓ=2 oscillation: inviscid (η̂=0) vs viscous (η̂=0.16)")
lines!(axB, series[1][1], series[1][2], color=:black, label="η̂=0 (inviscid)")
lines!(axB, series[end][1], series[end][2], color=:crimson, label="η̂=0.16 (viscous)")
axislegend(axB, position=:rt, framevisible=true)

Label(fig[0, :], "BDNKStar — STAGE 3: leading-order (NS-limit) BDNK shear viscosity damps the f-mode " *
      "(3+1D Cowling; 1/τ ∝ η̂; frame ⇒ UV causality, checked separately)", fontsize=12, font=:bold)

save(joinpath(outdir, "bdnk_viscous_damping.png"), fig)
println("saved bdnk_viscous_damping.png  (slope 1/τ per η̂ = ", round(slope, digits=4), ")")
