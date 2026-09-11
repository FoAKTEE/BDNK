#=
    doubly_diffusive_modes.jl — LOCAL (WKB) doubly-diffusive instabilities of a
    thermally-stratified, viscous + heat-conducting BDNK (causal first-order)
    neutron-star fluid element. Two KNOWN families re-derived in the BDNK context
    (DoublyDiffusive.jl), wired to the toolkit's own N² (brunt_vaisala) and
    transport sector (conformal_frame_PMP → η, κ_Q, τ):

      A  GSF / rotational doubly-diffusive — the unstable wavevector cone vs the
         shear (κ_epi²) and the Prandtl ratio κ_th/ν: thermal diffusion erodes the
         stabilizing N² and WIDENS the cone beyond the adiabatic Solberg–Høiland
         limit. (Goldreich–Schubert 1967; Acheson 1978; Menou-Balbus-Spruit 2004.)
      B  Thermohaline / semiconvective overstable growth Re s(k) vs the Prandtl
         number Pr=ν/χ on a convectively STABLE (N²>0) stratification — overstable
         only for Pr<1. (Stern 1960; Baines–Gill 1969.)
      C  Heat-coupled g-mode: the oscillatory g-root acquires Re s ∝ κ_Q (DRIVEN
         for N²_T<0, Pr<1) — the diffusive continuation of the thermohaline branch.
      D  BDNK causal high-k cutoff: the NS diffusive rate |Re s|~Dk² grows without
         bound; the causal telegrapher D k²/(1+sτ) SATURATES it near the relaxation
         pole 1/τ.

    HONEST: the instabilities are textbook; the new content is their BDNK
    realization + the causal D/(1+sτ) modification (UV cutoff; τ-invariant onset).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/doubly_diffusive_modes.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# ── background: stratified, SH-stable ideal-gas star (Γ>Γ_struct ⇒ N²>0) ───── #
star = solve_tov_idealgas(; Γ=2.0, Γ_struct=1.9, K=100.0, ρc=1.28e-3)
r, N2, _, _ = brunt_vaisala(star)
i0   = argmin(abs.(r .- 0.5*star.R))
N2_0 = N2[i0]
tc   = conformal_frame_PMP(star.ε[i0])
w0   = star.ε[i0] + star.p[i0]

@printf("Background: M=%.3f Msun  R=%.2f km  N²=%.3e km⁻² (>0 ⇒ SH stable)\n",
        mass_solar(star), star.R, N2_0)
@printf("Transport frame: η=%.3e κ_Q=%.3e τ_ε=%.3e τ_Q=%.3e  Pr⁻¹=%.2f\n",
        tc.η, tc.κQ, tc.τε, tc.τQ, tc.κQ/tc.η)

# WKB scaling (matches repro/gsf_doubly_diffusive.jl): pin νk²≈0.3|N|
k_test = 50.0 / star.R
ν_kin  = 0.3 * sqrt(abs(N2_0)) / k_test^2
τη = (tc.τε/tc.η) * ν_kin

# rotation Ω(ϖ)=Ω0(ϖ/R)^q at ϖ/R=0.5
Ω0 = 0.05 * sqrt(star.ε[i0])
Ωloc(q) = Ω0 * 0.5^q

fig = Figure(size=(1180, 900))

# ─────────────────────────────────────────────────────────────────────────── #
# Panel A: GSF unstable wavevector cone vs shear q and Prandtl ratio κ_th/ν
# ─────────────────────────────────────────────────────────────────────────── #
axA = Axis(fig[1,1], xlabel="shear  q = dlnΩ/dlnϖ", ylabel="unstable cone half-angle  ψ_max  [deg]",
           title="A  GSF: thermal diffusion WIDENS the unstable cone (erodes N²>0)")
qgrid = range(-3.5, -1.5; length=40)
for (Prinv, col, lbl) in [(0.0, :black, "adiabatic (Solberg–Høiland)"),
                          (3.0, :goldenrod, "κ_th/ν=3"),
                          (30.0, :crimson, "κ_th/ν=30 (GSF regime)")]
    cones = Float64[]
    for q in qgrid
        κ2 = gsf_epicyclic(Ωloc(q), q)
        if Prinv == 0
            c, _, _ = gsf_unstable_cone(κ2, N2_0, 0.0, 0.0, k_test)
        else
            c, _, _ = gsf_unstable_cone(κ2, N2_0, ν_kin, Prinv*ν_kin, k_test)
        end
        push!(cones, c)
    end
    lines!(axA, collect(qgrid), cones, color=col, linewidth=2.6, label=lbl)
end
vlines!(axA, [-2.0], color=:gray, linestyle=:dash)   # Rayleigh boundary κ_epi²=0
text!(axA, -2.0, 5; text="Rayleigh\nκ_epi²=0", align=(:center,:bottom), fontsize=9, color=:gray)
axislegend(axA, position=:lt, framevisible=true, labelsize=10)

# ─────────────────────────────────────────────────────────────────────────── #
# Panel B: thermohaline overstable growth Re s(k) vs Prandtl number Pr=ν/χ
# ─────────────────────────────────────────────────────────────────────────── #
axB = Axis(fig[1,2], xlabel="wavenumber  k  [km⁻¹]", ylabel="overstable growth  Re s  [km⁻¹]",
           xscale=log10, title="B  Thermohaline overstability vs Pr=ν/χ (net N²>0, stable)")
N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0)
χ0  = 1.0e-3
k_res = sqrt(sqrt(abs(N2T+N2mu)) / χ0)
kgrid = 10 .^ range(log10(k_res)-2.0, log10(k_res)+1.5; length=240)
for (Pr, col) in [(0.05, :navy), (0.2, :seagreen), (0.5, :goldenrod), (1.5, :crimson)]
    res = Float64[]
    for k in kgrid
        rs = dd_cubic_roots(k; ν=Pr*χ0, χ=χ0, N2mu=N2mu, N2T=N2T)
        push!(res, max(first(most_unstable(rs)), 0.0))
    end
    lines!(axB, kgrid, res, color=col, linewidth=2.4,
           label=@sprintf("Pr=%.2f%s", Pr, Pr<1 ? "" : " (stable)"))
end
hlines!(axB, [0.0], color=:gray, linestyle=:dot)
axislegend(axB, position=:rt, framevisible=true, labelsize=10)

# ─────────────────────────────────────────────────────────────────────────── #
# Panel C: heat-coupled g-mode — Re s ∝ κ_Q (driven), Im s ≈ √N² (oscillatory)
# ─────────────────────────────────────────────────────────────────────────── #
axC = Axis(fig[2,1], xlabel="thermal diffusion  χ k²/N   (∝ κ_Q, weak-diffusion)",
           ylabel="g-mode  Re s  [km⁻¹]",
           title="C  Heat-coupled g-mode: DRIVEN ∝ κ_Q  (N²_T<0, Pr<1)")
# Weak-diffusion regime χk²≲N keeps the g-root oscillatory; its Re s ∝ κ_Q is the
# heat-coupled drive. We sweep χ so the dimensionless χk²/N goes 0→~0.5; ν is held
# at the BDNK frame Prandtl (Pr=ν/χ=η/κ_Q<1, the driven regime).
Pr_frame = tc.η / tc.κQ               # =ν/χ from the frame (<1 ⇒ driven)
Nbar = sqrt(abs(N2T+N2mu))
k_gm = sqrt(k_res)                    # a moderate wavenumber for the g-mode
xrat = range(0.0, 0.5; length=40)     # χk²/N
res_g = Float64[]; ims_g = Float64[]
for x in xrat
    χg = x * Nbar / k_gm^2
    s = gmode_root(k_gm, χg; ν=Pr_frame*χg, N2mu=N2mu, N2T=N2T)
    push!(res_g, real(s)); push!(ims_g, abs(imag(s)))
end
facs = xrat
lines!(axC, collect(facs), res_g, color=:purple, linewidth=2.8, label="Re s (growth ∝ κ_Q)")
hlines!(axC, [0.0], color=:gray, linestyle=:dot)
axC2 = Axis(fig[2,1], yaxisposition=:right, ylabel="Im s ≈ √N²  [km⁻¹]",
            yticklabelcolor=:teal, ygridvisible=false)
hidespines!(axC2); hidexdecorations!(axC2)
lines!(axC2, collect(facs), ims_g, color=:teal, linewidth=2.0, linestyle=:dash)
linkxaxes!(axC, axC2)
axislegend(axC, position=:lt, framevisible=true, labelsize=10)
text!(axC2, 0.45, ims_g[end]; text="Im s (osc.)", align=(:right,:bottom),
      color=:teal, fontsize=10)

# ─────────────────────────────────────────────────────────────────────────── #
# Panel D: BDNK causal high-k cutoff — NS rate ~Dk² unbounded vs BDNK saturation
# ─────────────────────────────────────────────────────────────────────────── #
axD = Axis(fig[2,2], xlabel="wavenumber  k  [km⁻¹]", ylabel="fastest |Re s|  [km⁻¹]",
           xscale=log10, yscale=log10,
           title="D  BDNK causal cutoff: NS ~Dk² unbounded; BDNK saturates ~1/τ")
χh = tc.κQ/w0; νh = 0.2*χh; τπ = tc.τε; τQ = tc.τQ
kk = 10 .^ range(log10(k_res)-1, log10(k_res)+4.5; length=260)
fast(rs) = maximum(abs.(real.(rs)))
ns = [fast(dd_cubic_roots(k; ν=νh, χ=χh, N2mu=N2mu, N2T=N2T)) for k in kk]
bd = [fast(dd_bdnk_roots(k; ν=νh, χ=χh, N2mu=N2mu, N2T=N2T, τπ=τπ, τQ=τQ)) for k in kk]
lines!(axD, kk, ns, color=:crimson, linewidth=2.6, label="Navier–Stokes (∝ Dk², acausal)")
lines!(axD, kk, bd, color=:navy,    linewidth=2.6, label="BDNK causal D/(1+sτ)")
hlines!(axD, [1/τQ], color=:teal, linestyle=:dash)
text!(axD, kk[end], 1/τQ; text="1/τ_Q", align=(:right,:bottom), color=:teal, fontsize=10)
hlines!(axD, [1/τπ], color=:orange, linestyle=:dash)
text!(axD, kk[1], 1/τπ; text="1/τ_π", align=(:left,:bottom), color=:orange, fontsize=10)
axislegend(axD, position=:lt, framevisible=true, labelsize=10)

Label(fig[0, :],
      "BDNKStar — LOCAL doubly-diffusive instabilities in causal (BDNK) hydro: GSF / rotational (A) + " *
      "thermohaline overstability (B) + heat-coupled g-mode (C) + BDNK high-k cutoff (D)\n" *
      "KNOWN families (GSF 1967; thermohaline/Stern 1960) re-derived on a BDNK star — new: the causal D/(1+sτ) signature",
      fontsize=12, font=:bold)

save(joinpath(outdir, "doubly_diffusive_modes.png"), fig)
println("saved doubly_diffusive_modes.png")
