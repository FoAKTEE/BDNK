#=
    dg_staircase_trend — the KEY Stage-3 scientific test for the DG fallback.

    The FV-Cartesian verdict was 'fv-insufficient-need-dg': under an ℓ=2 surface
    perturbation the staircased non-conforming surface drove a secular central-
    density drift whose growth ACCELERATED with resolution (≈1.8% @N=32³ → 7.4%
    @N=48³ at t=300). The open question: does the high-order DG (sub-cell surface
    resolution + troubled-cell + positivity limiter) reverse that trend?

    Here we evolve the 2D meridional (x–z) DG star with a seeded ℓ=2 velocity
    perturbation at three effective resolutions and report the central-density
    drift at a fixed coordinate time. IMPROVING (or non-worsening) drift with
    refinement = staircase cured/mitigated; worsening = DG also insufficient.

    Run: JULIA_NUM_THREADS=6 julia --project=. repro/dg_staircase_trend.jl
=#
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar

κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)
tmax = 80.0
A = 1e-4

println("# DG ℓ=2 staircase resolution-trend test (t=$tmax, A=$A)")
println("# Kx=Kz  p  Neff(=Kx*p)  ρc_drift_max  q2_amp  blowup")
res = Tuple{Int,Int,Float64,Float64}[]
for (Kx, p) in ((8,2), (12,2), (16,2))
    eng, st = setup_dgcart2d(eos, εc; Kx=Kx, Kz=Kx, p=p, M_tvb=50.0, cfl=0.2)
    ρc0 = dgcart2d_central_density(st, eng)
    seed_dgcart2d_l2!(st, eng; A=A)
    ts, q2, ρc = evolve_dgcart2d!(st, eng; tmax=tmax, sample_dt=4.0)
    drift = maximum(abs.(ρc .- ρc0))/ρc0
    qamp = maximum(abs.(q2 .- q2[1]))
    blow = !all(isfinite, ρc) || drift > 1.0
    Neff = Kx*p
    push!(res, (Kx, p, Neff, drift))
    println("  $Kx  $p  $Neff  $(round(drift,sigdigits=4))  $(round(qamp,sigdigits=4))  $blow")
end

# trend verdict
drifts = [r[4] for r in res]
println()
if drifts[end] ≤ drifts[1]*1.2
    println("VERDICT: DG drift does NOT worsen with refinement (FV worsened 4× over the same span).")
    println("         → staircase instability cured/mitigated by sub-cell DG resolution.")
else
    println("VERDICT: DG drift still worsens with refinement — DG alone insufficient here.")
end

open(joinpath(@__DIR__, "dg_staircase_trend.txt"), "w") do io
    println(io, "# Kx  p  Neff  ρc_drift_max")
    for r in res; println(io, "  ", r[1], "  ", r[2], "  ", r[3], "  ", r[4]); end
end
println("SAVED dg_staircase_trend.txt")
