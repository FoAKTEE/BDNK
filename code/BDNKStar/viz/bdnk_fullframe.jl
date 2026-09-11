#=
    STAGE-3 capstone: FULL-FRAME (causal, hyperbolic) BDNK vs the parabolic
    Navier–Stokes-limit, in 3+1D Cowling.

      A  ℓ=2 amplitude |q(t)| at large shear viscosity (η̂=2), same grid and same
         dt=0.16dx for both engines. The parabolic NS-limit (CowlingEvolve3D) is
         outside its diffusive CFL here and OVERFLOWS (|q|→1.2e302, then NaN,
         within 1400 steps). The full-frame causal BDNK (CowlingBDNK3D, frame
         times τ_ε,τ_P,τ_Q + time-derivative recovery) stays FINITE at the same
         dt — that is the hyperbolic advantage, and it is real.
      B  BDNK excess damping 1/τ vs η̂, sampled over the STABLE range η̂ ∈ [0,0.5]
         (outside it the runs are unbounded — see below). Linear in η̂ up to
         η̂≈0.25, then rolling over as the mode approaches the overdamped regime.

    TWO CLAIMS THIS FIGURE USED TO MAKE, BOTH MEASURED AND BOTH WRONG:
      • "BDNK stays BOUNDED at η̂=2". It does not. From a 7.5e-4 seed it reaches
        9.1 by n=1400, 5.4e+09 by n=2600, 3.4e+19 by n=4000. It is finite, and it
        grows far more slowly than the NS overflow, but it is NOT bounded. The
        advantage demonstrated here is over the TIMESTEP (dt∼dx vs dt∼dx²), not
        over stability at large η̂.
      • "damping stable across the range, saturating". The apparent saturation was
        artefact: at η̂≥1 the runs are unstable, and damping_rate either could not
        measure them (<3 peaks ⇒ NaN) or fitted a NEGATIVE rate to a growing
        signal. Panel B now keeps only bounded, ringing runs — see the filter.
    Boundedness in this engine holds for η̂ ≲ 1 at dt=0.08dx (measured at N=22:
    max|q| ratio 1.000 up to η̂=1, 1.6e3 at η̂=1.5); the divergence is pre-existing
    and is not removed by any excision setting.

    Caveat: the full-frame engine uses a core excision to tame surface-recovery
    stiffness, which shifts the absolute mode frequency; this figure is about the
    causal STABILITY structure, not a precise frequency reproduction.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bdnk_fullframe.jl
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
s = build_star3d(eos, εc; N=28, Lfac=1.2)
# dt chosen in the parabolic-CFL-unstable window for the NS-limit at η̂=2; the
# hyperbolic BDNK (dt∼dx, not dx²) stays FINITE here — same grid, same dt.
# (Finite, not bounded: it still grows slowly at η̂=2. See the header.)
dt = 0.16 * s.grid.dx

e_ns = setup_evo3d(s; σ_ko=0.02, η̂=2.0)
st_ns = EvolState(s.grid.N); seed_l2!(st_ns, e_ns; A=1e-3)
ts_ns, q_ns, _ = evolve3d!(st_ns, e_ns; dt=dt, nsteps=1400, sample=4)

e_bd = setup_bdnk3d(s; η̂=2.0, ζ̂=0.0, τ̂=8.0, σ_ko=0.03)
st_bd = BDNKState(s.grid.N); seed_bdnk_l2!(st_bd, e_bd; A=1e-3)
ts_bd, q_bd, _ = evolve_bdnk3d!(st_bd, e_bd; dt=dt, nsteps=1400, sample=4)

# clamp to [1e-6,1e3] so the NS runaway (→1e302, then NONFINITE) shows as hitting
# the cap rather than blowing the axis to 300 decades
absfin(t,q) = (tt=Float64[]; aa=Float64[]; for i in eachindex(q); isfinite(q[i]) || break; a=clamp(abs(q[i]),1e-6,1e3); push!(tt,t[i]); push!(aa,a); end; (tt,aa))
tn,an = absfin(ts_ns,q_ns); tb,ab = absfin(ts_bd,q_bd)
println(@sprintf("NS-limit finite=%s |q|max=%.2e   |   full-frame BDNK finite=%s |q|max=%.2e",
        all(isfinite,q_ns), maximum(abs,filter(isfinite,q_ns)), all(isfinite,q_bd), maximum(abs,filter(isfinite,q_bd))))

# Panel B: BDNK damping vs η̂
ηs = [0.0, 0.125, 0.25, 0.375, 0.5]; dr = Float64[]; qmax = Float64[]
for η̂ in ηs
    e = setup_bdnk3d(s; η̂=η̂, ζ̂=0.0, τ̂=8.0, σ_ko=0.03)
    st = BDNKState(s.grid.N); seed_bdnk_l2!(st, e; A=1e-3)
    ts, q2, _ = evolve_bdnk3d!(st, e; dt=dt, nsteps=2600, sample=4)
    push!(dr, damping_rate(ts, q2)); push!(qmax, maximum(abs, q2))
end
d0 = dr[1]; excess = dr .- d0
# THE η̂ GRID IS THE STABLE RANGE, NOT THE OLD [0,0.5,1,2,4]. A damping rate is
# only meaningful for a bounded, ringing run, and the old grid failed that in two
# independent ways (measured, nsteps=2600, dt=0.16dx, N=28):
#     η̂        0        0.5       1.0       2.0       4.0
#     peaks    8        6         2         3         0
#     rate    +2.37e-2 +2.98e-2   NaN      -4.53e-2   NaN
#     |q|max   7.5e-4   7.5e-4    7.7e-1    5.4e+09   7.8e+21
#  • η̂=1,4 were UNMEASURABLE — <3 peaks, so damping_rate returns NaN.
#  • η̂=2 returned a perfectly FINITE rate that is NEGATIVE: a least-squares fit to
#    a GROWING signal (|q| up 13 orders from the seed). An isfinite check alone
#    would pass it straight onto the plot as a large "damping".
# Hence the filter below tests BOTH finiteness and boundedness, and is kept even
# though every point on the present grid passes it.
# ON THE PRESENT GRID (all bounded at the 7.466e-4 seed, 6–10 peaks):
#     η̂        0       0.125    0.25     0.375    0.5
#     peaks    8       8        9        10       6
#     excess   0       4.97e-3  9.95e-3  1.14e-2  6.08e-3
#     excess/η̂ —       0.0398   0.0398   0.0303   0.0122
# So 1/τ is LINEAR in η̂ only up to η̂≈0.25 (0.0398 twice), then rolls over and
# falls. The rollover coincides with the mode losing its ring (peaks 10→6 at
# η̂=0.5), i.e. the approach to the overdamped regime where the envelope fit is
# measuring less and less of an oscillation — read the last point with that in mind.
bounded = qmax .< 5 * qmax[1]
keep = isfinite.(excess) .& bounded
if !all(keep)
    @printf("  NOTE: omitted from panel B — unmeasurable (<3 peaks) at η̂ = %s; unbounded at η̂ = %s\n",
            string(ηs[.!isfinite.(excess)]), string(ηs[.!bounded]))
end

fig = Figure(size=(1180, 450))
axA = Axis(fig[1,1], xlabel="t  [M⊙ geometric]", ylabel="ℓ=2 amplitude |q|", yscale=log10,
           title="A  η̂=2, same dt: NS-limit overflows; BDNK finite but growing")
lines!(axA, tn, an, color=:crimson, label="NS-limit (parabolic) — overflows to 1e302, then NaN")
lines!(axA, tb, ab, color=:seagreen, label="full-frame BDNK (causal) — finite, slow growth")
ylims!(axA, 1e-6, 3e3)
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xlabel="shear-viscosity knob η̂", ylabel="excess damping 1/τ [geom]",
           title="B  BDNK excess damping: linear to η̂≈0.25, then rolls over")
scatterlines!(axB, ηs[keep], excess[keep], color=:seagreen, markersize=12)
# guide to the eye: the linear-in-η̂ regime, fitted through the origin on η̂ ≤ 0.25
lin = keep .& (ηs .<= 0.25)
if count(lin) ≥ 2
    slope = sum(ηs[lin] .* excess[lin]) / sum(ηs[lin] .^ 2)
    lines!(axB, [0, maximum(ηs[keep])], [0, slope*maximum(ηs[keep])],
           color=:gray, linestyle=:dash, label=@sprintf("linear regime: 1/τ = %.3f η̂", slope))
    axislegend(axB, position=:lt, framevisible=true)
end

Label(fig[0, :], "BDNKStar — STAGE 3: at dt=0.16dx the parabolic NS-limit overflows while the causal " *
      "BDNK stays finite — but at η̂=2 neither is bounded", fontsize=12, font=:bold)
save(joinpath(outdir, "bdnk_fullframe.png"), fig)
println("saved bdnk_fullframe.png  (BDNK damping vs η̂: ",
        round.(excess[keep], digits=4), " at η̂ = ", ηs[keep], ")")
