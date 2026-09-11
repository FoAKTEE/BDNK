#=
    STAGE-3 result: GRAVITY (g) modes are a genuine new mode family that appears
    ONLY when the star is stratified. The relativistic Cowling polar eigensolver
    (GravityModes.jl) extends the validated barotropic f/p (W,V) system with the
    buoyancy / Schwarzschild-discriminant term:

        N²(r) = g² (1/c_e² − 1/c_s²) e^{ν−λ}   (relativistic Brunt–Väisälä)

    with c_s² the ADIABATIC sound speed (perturbation index Γ, frozen composition)
    and c_e²=dp/dε the EQUILIBRIUM sound speed of the hydrostatic structure (index
    Γ_struct). Γ>Γ_struct ⇒ N²>0 ⇒ a tower of STABLE, real g-modes BELOW the f-mode.

      A  N²(r): ≈0 everywhere for the barotropic (Γ=Γ_struct) star — NO buoyancy,
         NO g-modes; strictly POSITIVE in the interior for the stratified
         (Γ>Γ_struct) star — the buoyancy that traps the g-modes.
      B  the spectrum: barotropic = f/p only; stratified = a g-tower (g₁>g₂>…)
         BELOW the f-mode, with the f and p-modes essentially unchanged.

    Eqs.: Jaikumar et al. 2021 / Shirke et al. arXiv:2506.18892 Eq.12-13, with the
    Schwarzschild-discriminant sign of Gaertig & Kokkotas arXiv:0905.0821 Eq.12.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/gravity_modes.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# Two stars sharing K, ρc — differ ONLY in the structure index Γ_struct.
K = 30.0; ρc = 3e-3; Γ = 2.0
star_baro = solve_tov_idealgas(Γ=Γ, Γ_struct=Γ,   K=K, ρc=ρc, h=2e-4)   # N²≡0
star_strat= solve_tov_idealgas(Γ=Γ, Γ_struct=5/3, K=K, ρc=ρc, h=2e-4)   # N²>0

rb, N2b, _, _ = brunt_vaisala(star_baro)
rs, N2s, _, _ = brunt_vaisala(star_strat)

# convert N² (geometric km⁻²) to a kHz "Brunt frequency" for an intuitive y-axis
N2_to_kHz(N2) = sign(N2) * freq_kHz_from_omega2_g(abs(N2))

sp_b = gmode_spectrum(star_baro;  l=2, N=9000,  nscan=2000, ω2lo=1e-5, ω2hi=0.3)
sp_s = gmode_spectrum(star_strat; l=2, N=12000, nscan=3000, ω2lo=1e-5, ω2hi=0.3)

@printf("BAROTROPIC  M=%.3f R=%.2f km  N2max=%.2e  has_g=%s  f=%.3f kHz\n",
        mass_solar(star_baro), star_baro.R, maximum(N2b), sp_b.has_gmodes, sp_b.f_freq)
@printf("STRATIFIED  M=%.3f R=%.2f km  N2max=%.2e  has_g=%s  f=%.3f kHz\n",
        mass_solar(star_strat), star_strat.R, maximum(N2s), sp_s.has_gmodes, sp_s.f_freq)
println("  g-modes [kHz]: ", round.(sp_s.g_freqs_kHz, digits=4))
println("  p-modes [kHz]: ", round.(sp_s.p_freqs, digits=4))

fig = Figure(size=(1040, 440))

# Panel A: N²(r) profile
axA = Axis(fig[1,1], xlabel="r / R", ylabel="N²  [km⁻²]",
           title="A  Brunt–Väisälä N²(r): buoyancy only when stratified")
lines!(axA, rb ./ star_baro.R,  N2b, color=:navy,   linewidth=2.5,
       label=@sprintf("barotropic Γ_struct=%.3f  (N²≈0)", star_baro.Γ_struct))
lines!(axA, rs ./ star_strat.R, N2s, color=:crimson, linewidth=2.5,
       label="stratified Γ_struct=5/3  (N²>0)")
hlines!(axA, [0.0], color=:gray, linestyle=:dash)
axislegend(axA, position=:lt, framevisible=true)

# Panel B: spectrum — g below f, p above f
axB = Axis(fig[1,2], xlabel="mode family", ylabel="frequency  [kHz]",
           title="B  g-modes are a NEW family below the f-mode (stratified only)",
           xticks=([1,2], ["barotropic", "stratified Γ>Γ_struct"]))
# barotropic column (x=1): f + p only
scatter!(axB, fill(1.0, 1), [sp_b.f_freq], color=:black, markersize=15, marker=:star5, label="f-mode")
scatter!(axB, fill(1.0, length(sp_b.p_freqs)), sp_b.p_freqs, color=:dodgerblue, markersize=11, label="p-modes")
# stratified column (x=2): g (red, below) + f (star) + p (blue, above)
scatter!(axB, fill(2.0, length(sp_s.g_freqs_kHz)), sp_s.g_freqs_kHz, color=:crimson, markersize=11, marker=:diamond, label="g-modes")
scatter!(axB, fill(2.0, 1), [sp_s.f_freq], color=:black, markersize=15, marker=:star5)
scatter!(axB, fill(2.0, length(sp_s.p_freqs)), sp_s.p_freqs, color=:dodgerblue, markersize=11)
hlines!(axB, [sp_s.f_freq], color=:gray, linestyle=:dot)
text!(axB, 2.05, sp_s.f_freq; text="f", align=(:left,:center), fontsize=12)
axislegend(axB, position=:rc, framevisible=true)
xlims!(axB, 0.5, 2.8)

Label(fig[0, :], "BDNKStar — STAGE 3: relativistic Cowling g-mode eigensolver — buoyancy (N²>0) " *
      "produces a NEW g-mode tower below the f-mode only in the stratified star",
      fontsize=12, font=:bold)

save(joinpath(outdir, "gravity_modes.png"), fig)
println("saved gravity_modes.png")
