#=
    STAGE-3 Phase-2 validation figure: the 3+1D Cowling time-domain evolution
    reproduces the non-radial mode spectrum of the Phase-1 eigensolver.

      A  power spectrum of the ℓ=2 quadrupole from a full 3D Cartesian evolution
         (perturb ℓ=2 → evolve → FFT), with the Phase-1 eigensolver f/p1/p2 mode
         frequencies overlaid as dashed lines.
      B  the ℓ=2 quadrupole time-series that is Fourier-analysed.

    Two independent methods (frequency-domain shooting vs full 3D grid evolution)
    agreeing at the few-percent level — the cross-validation of Stage 3.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/cowling3d_fmode.jl
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
M2c = BDNKStar.Units.Msun_to_km * BDNKStar.Units.kHz_to_km   # kHz → cyclic geometric

# --- Phase-1 eigensolver reference lines ---
feig, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3, N=4000, nscan=600)

# --- Phase-2: full 3D evolution ---
N = 44
s = build_star3d(eos, εc; N=N, Lfac=1.2)
e = setup_evo3d(s; σ_ko=0.004)
st = EvolState(s.grid.N); seed_l2!(st, e; A=1e-3)
ts, q2, qc = evolve3d!(st, e; dt=0.30*s.grid.dx, nsteps=15000, sample=6)

# power spectrum over a kHz band
fkhz = range(0.5, 7.0; length=3000)
P = periodogram(ts, q2, fkhz .* M2c); P ./= maximum(P)
# refined evolved f-mode peak in the f band
fband = range(1.4, 2.4; length=2000); Pf = periodogram(ts, q2, fband .* M2c)
mi = argmax(Pf); fevo = fband[mi]
println(@sprintf("evolved f-mode = %.4f kHz  vs eigensolver %.4f kHz  (%.2f%%)",
        fevo, feig[1], 100*abs(fevo-feig[1])/feig[1]))

fig = Figure(size=(960, 430))
axA = Axis(fig[1,1], xlabel="frequency [kHz]", ylabel="power (normalised)", yscale=log10,
           title="A  3D-evolution ℓ=2 power spectrum vs eigensolver modes")
lines!(axA, collect(fkhz), max.(P, 1e-6), color=:steelblue)
for (lbl,fv,col) in (("f",feig[1],:crimson),("p₁",feig[2],:seagreen),("p₂",feig[3],:darkorange))
    vlines!(axA, [fv], color=col, linestyle=:dash, label="$lbl=$(round(fv,digits=3)) kHz (eig)")
end
ylims!(axA, 1e-5, 2); axislegend(axA, position=:rt, framevisible=true)

axB = Axis(fig[1,2], xlabel="t  [M⊙ geometric]", ylabel="ℓ=2 quadrupole of δε",
           title="B  evolved ℓ=2 quadrupole time-series  (N=$(N)³, $(round(Int,ts[end]/108)) f-periods)")
lines!(axB, ts, q2, color=:black)

Label(fig[0, :], "BDNKStar — STAGE 3 Phase 2: 3+1D Cowling evolution reproduces the " *
      "non-radial f-mode (evolved $(round(fevo,digits=3)) vs eigensolver $(round(feig[1],digits=3)) kHz)",
      fontsize=13, font=:bold)

save(joinpath(outdir, "cowling3d_fmode.png"), fig)
println("saved cowling3d_fmode.png")
