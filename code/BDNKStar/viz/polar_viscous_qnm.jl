#=
    STAGE-3 result: POLAR (even-parity) VISCOUS quasi-normal modes — the even-parity
    counterpart of the repo's axial viscous-QNM track. Frequency-domain matrix
    eigenproblem (PolarViscousModes): eigenvalues λ of the 1D ℓ-reduced linearised polar
    BDNK–Cowling operator give freq=|Im λ| and damping γ=−Re λ.

      A  f- and p₁-mode damping rate γ vs the shear-viscosity knob η̂: γ>0 rising ∝η̂
         (viscosity DAMPS the polar modes), with p₁ damping faster than f (higher k ⇒
         more shear). The ideal (η̂→0) limit is non-dissipative (γ→0).
      B  the complex spectrum at η̂=0.02: the f/p modes are WEAKLY damped (high-Q, low γ);
         the BDNK frame/relaxation modes are HEAVILY damped (overdamped) — these are the
         polar non-hydrodynamic modes, the only "viscosity-driven" addition (no long-lived
         η-mode analog, unlike the axial sector).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/polar_viscous_qnm.jl
=#
using Pkg; Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie, LinearAlgebra, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2

# Panel A: γ(η̂) for f and p₁
ηs = collect(0.0:0.005:0.04)
γf = Float64[]; γp = Float64[]; ff = Float64[]
for η̂ in ηs
    m = polar_qnm(eos, εc; l=2, η̂=η̂, Nr=140, nmodes=2)
    push!(γf, qnm_damping(m[1])); push!(γp, qnm_damping(m[2])); push!(ff, qnm_freq_kHz(m[1]))
end
sf = sum(ηs.*γf)/sum(ηs.^2); sp = sum(ηs.*γp)/sum(ηs.^2)
@printf("f-mode γ slope dγ/dη̂ = %.3f ; p₁-mode = %.3f  (p₁ damps faster)\n", sf, sp)

# Panel B: full complex spectrum at η̂=0.02
s = build_sphstar(eos, εc; Nr=140, Nθ=2)
L = polar_bdnk_operator(s, 2; η̂=0.02)        # default frame = shear-viscous isolation
ev = eigvals(L)
osc = [e for e in ev if imag(e)>1e-6]
fr  = [qnm_freq_kHz(e) for e in osc]; γa = [max(qnm_damping(e),1e-5) for e in osc]
mν  = polar_qnm(eos, εc; l=2, η̂=0.02, Nr=140, nmodes=2)

fig = Figure(size=(980, 430))
axA = Axis(fig[1,1], xlabel="shear-viscosity knob η̂", ylabel="damping rate γ = −Re λ [geom]",
           title="A  viscosity damps the polar f/p modes (γ ∝ η̂)")
scatter!(axA, ηs, γf, color=:crimson, markersize=11, label="f-mode")
lines!(axA, ηs, sf.*ηs, color=:crimson, linestyle=:dash)
scatter!(axA, ηs, γp, color=:navy, markersize=11, label="p₁-mode")
lines!(axA, ηs, sp.*ηs, color=:navy, linestyle=:dash)
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xlabel="frequency |Im λ| [kHz]", ylabel="damping γ = −Re λ [geom]",
           yscale=log10, title="B  spectrum at η̂=0.02: high-Q f/p vs overdamped frame modes")
scatter!(axB, fr, γa, color=(:gray,0.5), markersize=6, label="all modes")
scatter!(axB, [qnm_freq_kHz(mν[1])], [qnm_damping(mν[1])], color=:crimson, markersize=13, label="f")
scatter!(axB, [qnm_freq_kHz(mν[2])], [qnm_damping(mν[2])], color=:navy, markersize=13, label="p₁")
axislegend(axB, position=:rb, framevisible=true)

Label(fig[0,:], "BDNKStar — STAGE 3: polar (even-parity) viscous quasi-normal modes of the M=1.4 M☉ star " *
      "(frequency-domain BDNK–Cowling eigenproblem)", fontsize=12, font=:bold)
save(joinpath(outdir,"polar_viscous_qnm.png"), fig)
println("saved polar_viscous_qnm.png")
