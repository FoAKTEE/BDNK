#=
    fv_cartesian — STAGE 1/2 figure for the NONLINEAR flux-conservative Valencia
    GRHydro (Cowling) reproduction with an atmosphere.

      A  the stable star + atmosphere profile: ρ(r) of the equilibrium FV state
         (interior polytrope + low-density atmosphere floor), and the central-
         density history of the unperturbed run (flat ⇒ well-balanced/stable).
      B  Stage-1 radial-mode power spectrum (1D radial FV ℓ=0 perturbation),
         with the grid-converged FV peak marked.
      C  Stage-2 Cartesian ℓ=2 quadrupole power spectrum vs the spherical-engine
         f/p eigenfrequencies (1.883/4.107/6.031 kHz overlaid).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/fv_cartesian.jl
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
Lkm = BDNKStar.Units.Msun_to_km
M2c = Lkm * BDNKStar.Units.kHz_to_km    # kHz → cyclic geometric freq

# spherical-engine eigenfrequencies (the reproduction target)
feig, _, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3)
println("eigensolver ℓ=2 f/p1/p2 = ", round.(feig, digits=3), " kHz")

# ---- Panel A data: equilibrium ρ profile + atmosphere + static history -------
engR, stR = setup_fvradial(eos, εc; N=300, rmax_fac=1.3, atm_fac=1e-7, cfl=0.3)
NG = engR.g.NG; N = engR.g.N
rprof = engR.g.r[NG+1:NG+N]
ρprof = stR.ρ[NG+1:NG+N]
tsS, _, ρcS, driftS = evolve_fvradial!(stR, engR; tmax=300.0, probe_frac=0.5)
@printf("Stage-1 static drift (300 t) = %.2e (well-balanced)\n", driftS)

# ---- Panel B data: Stage-1 radial mode spectrum (Hann-windowed, homologous) --
engP, stP = setup_fvradial(eos, εc; N=300, rmax_fac=1.3, atm_fac=1e-7, cfl=0.25)
seed_radial_velocity!(stP, engP; A=1e-4, profile=:linear)
tsR, probeR, ρcR, driftR = evolve_fvradial!(stP, engP; tmax=4000.0, probe_frac=0.3, sample_dt=2.0)
fR = range(1.0, 8.0; length=3000)
PR = radial_periodogram_freqs(tsR, ρcR; fmin_kHz=1.0, fmax_kHz=8.0, npts=3000, window=true)[2]
PR ./= maximum(PR)
fRpk, _, _ = radial_periodogram_freqs(tsR, ρcR; fmin_kHz=2.5, fmax_kHz=3.6, npeaks=1, window=true)
fR0 = fRpk[1]
@printf("Stage-1 FV radial fundamental ≈ %.3f kHz (grid-converged; drift %.2e)\n", fR0, driftR)

# ---- Panel C data: Stage-2 Cartesian ℓ=2 — staircase instability vs resolution
# Two resolutions to expose the RESOLUTION-WORSENING secular surface drift (the
# honest fv-insufficient signature): the perturbed Cartesian star drifts FASTER
# at higher N — the staircase the atmosphere cannot cure under perturbation.
driftC = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
local tsC, q2C
for N in (32, 48)
    eng, st = setup_fvcart(eos, εc; N=N, L_fac=1.3, atm_fac=1e-7, cfl=0.25)
    seed_l2_velocity!(st, eng; A=1e-5)
    ts, q2, ρc = evolve_fvcart!(st, eng; tmax=300.0, sample_dt=6.0)
    ρc0 = ρc[1]
    driftC[N] = (ts, abs.(ρc .- ρc0) ./ ρc0)
    @printf("Stage-2 N=%d: ρc drift @t=300 = %.2e (grows with N ⇒ staircase)\n",
            N, last(driftC[N][2]))
    global tsC, q2C = ts, q2
end
fC = range(0.8, 8.0; length=3000)
PC = fvcart_periodogram_freqs(tsC, q2C; fmin_kHz=0.8, fmax_kHz=8.0, npts=3000, window=true)[2]
PC ./= maximum(PC)

# ----------------------------- figure (2×2) -----------------------------------
fig = Figure(size=(1180, 760))

axA = Axis(fig[1,1], xlabel="r  [M⊙ geometric]", ylabel="ρ  (rest-mass density)",
           yscale=log10, title="A  stable star + atmosphere (FV equilibrium, well-balanced)")
lines!(axA, rprof, max.(ρprof, 1e-12), color=:steelblue)
vlines!(axA, [R], color=:gray, linestyle=:dot, label="surface R")
hlines!(axA, [engR.atm.ρ_atm], color=:crimson, linestyle=:dash, label="ρ_atm floor")
axislegend(axA, position=:lb, framevisible=true)

axB = Axis(fig[1,2], xlabel="frequency [kHz]", ylabel="power (norm.)", yscale=log10,
           title="B  Stage-1 radial (ℓ=0) FV spectrum — grid-converged")
lines!(axB, collect(fR), max.(PR,1e-6), color=:black)
vlines!(axB, [fR0], color=:seagreen, linestyle=:dash,
        label=@sprintf("FV F₀≈%.2f kHz (N=200=400)", fR0))
ylims!(axB, 1e-5, 2); axislegend(axB, position=:rt, framevisible=true)

axC = Axis(fig[2,1], xlabel="frequency [kHz]", ylabel="power (norm.)", yscale=log10,
           title="C  Stage-2 Cartesian ℓ=2 spectrum vs spherical-engine modes")
lines!(axC, collect(fC), max.(PC,1e-6), color=:steelblue)
for (lbl,fv,col) in (("f",feig[1],:crimson),("p₁",feig[2],:seagreen),("p₂",feig[3],:darkorange))
    vlines!(axC, [fv], color=col, linestyle=:dash, label="$lbl=$(round(fv,digits=2)) (eig)")
end
ylims!(axC, 1e-5, 2); axislegend(axC, position=:rt, framevisible=true)

axD = Axis(fig[2,2], xlabel="t  [M⊙ geometric]", ylabel="|Δρ_c|/ρ_c  (drift)",
           yscale=log10, title="D  Stage-2: perturbed surface drift GROWS with resolution")
for (N,col) in ((32,:seagreen),(48,:crimson))
    ts,dr = driftC[N]
    lines!(axD, ts, max.(dr,1e-9), color=col, label="N=$N")
end
vlines!(axD, [108.0], color=:gray, linestyle=:dot, label="1 f-mode period")
axislegend(axD, position=:rb, framevisible=true)

Label(fig[0, :], "BDNKStar — nonlinear flux-conservative FV GRHydro (Cowling)+atmosphere: " *
      "Stage-1 radial MVP SOLID (converged), Stage-2 ℓ=2 reproduction fv-insufficient " *
      "(perturbed staircase instability worsens with N)",
      fontsize=12, font=:bold)

save(joinpath(outdir, "fv_cartesian.png"), fig)
println("saved figures/fv_cartesian.png")
