#=
    GAUSSIAN-NOISE MODE EXCITATION (STAGE 3, boundary-conforming causal SphBDNK).

    A generic causal random "pluck" — band-limited Gaussian-noise initial data
    (seed_sphbdnk_noise!) seeded into the scalar density δρ — excites the star's
    physical (Cowling) mode spectrum. Four panels:

      A  ℓ=2 power spectrum from a random pluck, with the eigensolver f/p₁/p₂
         frequencies marked — the random peaks land on them.
      B  ℓ=0/2/4 projected-moment spectra stacked: each angular channel rings at
         its own f/p ladder.
      C  ideal vs viscous (η̂) ℓ=2 spectrum — shear viscosity DAMPS/broadens the
         peaks (lower, wider).
      D  seed-independence: the f-peak frequency coincides across RNG seeds; only
         the amplitude is a random-draw variable.

    Reusable ID: seed_sphbdnk_noise! builds A·ε₀·[Σ aₙ sin(nπr/R)]·[Σ bₗ Pₗ(cosθ)]
    with aₙ,bₗ~N(0,1); the sin/even-Legendre basis is constraint-valid (ε₀+δε>0)
    and parity-regular at centre+poles by construction. The production frame is
    causal (c₊=√3·cs ≤ 1 at every radius; 0<q̂<ŝ well-posed).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/noise_excitation.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar: Units
using CairoMakie
using Printf
using Statistics
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
Nr, Nθ = 64, 16
s = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ)
kHz2cyc(f) = f*Units.Msun_to_km*Units.kHz_to_km
println(@sprintf("Shum star: R=%.3f Msun (%.2f km)  M=%.4f Msun  Nr=%d Nθ=%d",
        s.R, s.R*Units.Msun_to_km, s.M, Nr, Nθ))

# eigensolver references (independent benchmark; peaks should land here)
benchℓ = Dict(l => nonradial_cowling_spectrum(eos, εc; l=l, nmodes=3, N=4000, nscan=800)[1]
              for l in (2, 4))
println("Cowling ℓ=2 (f,p1,p2) = ", round.(benchℓ[2], digits=3), " kHz")
println("Cowling ℓ=4 (f,p1,p2) = ", round.(benchℓ[4], digits=3), " kHz")

# Hann-windowed periodogram of a moment time series, evaluated in kHz.
function pgram(ts, q, νkHz)
    n=length(q); q̄=q .- mean(q)
    w=[0.5-0.5*cos(2π*(k-1)/(n-1)) for k in 1:n]; q̄ .*= w
    P=zeros(length(νkHz))
    @inbounds for (m,ν) in enumerate(νkHz)
        ω=2π*ν*kHz2cyc(1.0); re=0.0; im=0.0
        for k in 1:n; sc=sincos(ω*ts[k]); re+=q̄[k]*sc[2]; im-=q̄[k]*sc[1]; end
        P[m]=re^2+im^2
    end
    P
end

# random pluck → time series of ℓ=0/2/4 projected moments + mode energy
function pluck(seed; η̂=0.0, T=2600.0, A=1e-4)
    e  = setup_sphbdnk(s; η̂=η̂)
    st = SphBDNKState(Nr, Nθ)
    seed_sphbdnk_noise!(st, e; A=A, seed=seed)
    dt = 0.2*s.grid.dr; nst = round(Int, T/dt)
    ks = [SphBDNKState(Nr, Nθ) for _ in 1:5]
    ts=Float64[]; M=Dict(l=>Float64[] for l in (0,2,4)); en=Float64[]
    for n in 0:nst
        if n % 4 == 0
            push!(ts, n*dt)
            for l in (0,2,4); push!(M[l], lℓ_quad_sphbdnk(st,e,l)); end
            push!(en, mode_energy(st,e))
        end
        n == nst && break
        BDNKStar.SphBDNK._rk4!(st, e, dt, ks...)
    end
    ts, M, en
end

# ── data ──────────────────────────────────────────────────────────────────────
ν2 = collect(range(0.8, 7.0; length=1600))   # ℓ=2 spectroscopy band
ν4 = collect(range(0.8, 8.0; length=1600))

println("running random pluck (seed 42)...")
ts, M, en = pluck(42; η̂=0.0)
P2 = pgram(ts, M[2], ν2)
P0 = pgram(ts, M[0], ν2)
P4 = pgram(ts, M[4], ν4)

println("running viscous comparison (η̂=0.08)...")
tsv, Mv, env = pluck(42; η̂=0.08)
P2v = pgram(tsv, Mv[2], ν2)
println(@sprintf("  ideal E_end/E0=%.3f   viscous E_end/E0=%.3f   f-peak power %.3f×",
        en[end]/en[1], env[end]/env[1], maximum(P2v)/maximum(P2)))

println("running seed sweep...")
seeds = (42, 7, 20260620)
Pseeds = Dict{Int,Vector{Float64}}()
for sd in seeds
    tss, Ms, _ = pluck(sd; η̂=0.0)
    Pseeds[sd] = pgram(tss, Ms[2], ν2)
    println(@sprintf("  seed %d: f-peak = %.3f kHz", sd, ν2[argmax(Pseeds[sd][ν2 .< 2.6])]))
end

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size=(1180, 880))
modecols = [:firebrick, :darkorange, :goldenrod]   # f, p1, p2

# A — ℓ=2 spectrum with eigenfrequencies marked
axA = Axis(fig[1,1], xlabel="frequency  [kHz]", ylabel="power  (ℓ=2 moment)",
           title="A  random pluck → ℓ=2 spectrum; peaks land on the eigenfrequencies")
lines!(axA, ν2, P2 ./ maximum(P2), color=:navy)
for (k, f) in enumerate(benchℓ[2])
    vlines!(axA, [f], color=modecols[k], linestyle=:dash,
            label=@sprintf("%s = %.2f kHz", k==1 ? "f" : "p$(k-1)", f))
end
xlims!(axA, 0.8, 7.0); axislegend(axA, position=:rt, framevisible=true)

# B — ℓ=0/2/4 stacked
axB = Axis(fig[1,2], xlabel="frequency  [kHz]", ylabel="normalised power  (offset)",
           title="B  ℓ=0 / 2 / 4 projected-moment spectra")
lines!(axB, ν2, P0 ./ maximum(P0) .+ 0.0, color=:seagreen,   label="ℓ=0 (radial)")
lines!(axB, ν2, P2 ./ maximum(P2) .+ 1.2, color=:navy,       label="ℓ=2")
lines!(axB, ν4, P4 ./ maximum(P4) .+ 2.4, color=:purple,     label="ℓ=4")
for (k,f) in enumerate(benchℓ[2]); vlines!(axB, [f], ymin=0.4, ymax=0.73, color=:navy,   linestyle=:dot); end
for (k,f) in enumerate(benchℓ[4]); vlines!(axB, [f], ymin=0.73, ymax=1.0, color=:purple, linestyle=:dot); end
xlims!(axB, 0.8, 8.0); axislegend(axB, position=:rt, framevisible=true)

# C — ideal vs viscous (peaks broaden / damp)
axC = Axis(fig[2,1], xlabel="frequency  [kHz]", ylabel="power  (ℓ=2 moment)",
           title="C  shear viscosity damps & broadens the peaks")
lines!(axC, ν2, P2  ./ maximum(P2), color=:navy,    label="ideal  η̂=0")
lines!(axC, ν2, P2v ./ maximum(P2), color=:crimson, label="viscous  η̂=0.08")
vlines!(axC, [benchℓ[2][1]], color=:gray, linestyle=:dash, label="f-mode")
xlims!(axC, 1.2, 5.0); axislegend(axC, position=:rt, framevisible=true)

# D — seed-independence (f-peak frequency coincides; amplitude varies)
axD = Axis(fig[2,2], xlabel="frequency  [kHz]", ylabel="normalised power  (ℓ=2)",
           title="D  seed-independence: f-peak frequency fixed, amplitude varies")
scols = [:black, :dodgerblue, :darkorange]
for (k, sd) in enumerate(seeds)
    lines!(axD, ν2, Pseeds[sd] ./ maximum(Pseeds[sd]), color=scols[k], label="seed $sd")
end
vlines!(axD, [benchℓ[2][1]], color=:crimson, linestyle=:dash,
        label=@sprintf("eigensolver f=%.2f", benchℓ[2][1]))
xlims!(axD, 1.4, 2.6); axislegend(axD, position=:rt, framevisible=true)

Label(fig[0, :], "BDNKStar — STAGE 3: a generic causal Gaussian-noise pluck excites the star's physical " *
      "mode spectrum (peaks match the Cowling eigensolver), seed-independent in frequency, damped by viscosity",
      fontsize=13, font=:bold)

save(joinpath(outdir, "noise_excitation.png"), fig)
println("saved figures/noise_excitation.png")
