#=
    STAGE-3 result: FULL GENERAL-RELATIVISTIC POLAR (even-parity) quasi-normal
    modes — the Cowling (frozen-metric) approximation DROPPED.  The metric
    perturbations (Lindblom–Detweiler {H1,K,W,X} system, DL85/KGK) are evolved
    coupled to the fluid and matched to the exterior ZERILLI equation with the
    outgoing-wave (QNM) boundary condition (analytic Chandrasekhar–Detweiler
    transform of the Regge–Wheeler Leaver continued fraction).  The mode acquires
    a finite gravitational-wave damping time τ (Im ω ≠ 0) — impossible in the
    Cowling limit (τ=∞).  Module: PolarGRModes.

      A  the converged full-GR polar QNM in the complex-ω plane (a GW-DAMPED
         mode, Imω≠0), contrasted with the REAL Cowling f-mode (Imω=0, τ=∞).
      B  the interior LD metric response |H1(R)| vs real frequency — the fluid
         resonance structure — with the Cowling f-mode marked.

    The full-GR f-MODE (EOS1 star) is f=2.86 kHz, τ_GW=0.115 s — BELOW the Cowling
    f-mode (3.30 kHz; the Cowling approx overestimates f) with a finite GW damping
    time, matching the Andersson–Kokkotas (1998) universal relations to ≈8%/10%.
    A polar w-mode (f≈9.2 kHz, τ≈64 μs) is also recovered.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/polar_gr_qnm.jl
=#
using Pkg; Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const PG = BDNKStar.PolarGRModes
const MSUN = BDNKStar.Units.Msun_to_km

eos, star = build_polar_star()
M, R = star.M, star.R
G = 6.6743015e-11; c = 299_792_458.0
εc = (3e15 * 1e3) * G / c^2 * 1e6

# full-GR polar f-mode (GW-damped) and the polar w-mode
qnm = polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(2.8, 0.15), star=star)
@printf("full-GR f-mode:    f=%.4f kHz  τ=%.4e s  (Imω=%.3e ≠ 0 ⇒ GW-damped)\n",
        qnm.f_kHz, qnm.tau_s, imag(qnm.omega))
wm = polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(9.0, 5e-5), star=star)
@printf("full-GR w-mode:    f=%.4f kHz  τ=%.4e s\n", wm.f_kHz, wm.tau_s)

# Cowling f-mode (real, τ=∞) for the same star
fc, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=2, N=6000, nscan=800,
                                      Lunit_km=1.0, ω2lo=1e-4, ω2hi=0.2)
f_cowling = fc[1]
@printf("Cowling f-mode:    f=%.4f kHz  τ=∞ (metric frozen)\n", f_cowling)

# AK1998 universal-relation f-mode prediction (context)
xak = sqrt((M/(1.4*MSUN))/(R/10)^3); f_ak = 0.78 + 1.635*xak
@printf("AK1998 f-mode:     f=%.4f kHz (universal relation)\n", f_ak)

# interior metric response |H1(R)| across the f/p band (real ω)
fs = collect(1.4:0.01:4.2)
resp = [interior_metric_response(star, eos, PG.polar_ftau_to_omega(f, 1e6), 2) for f in fs]

# ---- figure ----
fig = Figure(size=(1040, 430))

# Panel A — complex-ω plane (kHz vs 1/τ)
axA = Axis(fig[1,1], xlabel="Re f = Re ω/2π  [kHz]",
           ylabel="damping rate 1/τ  [s⁻¹]",
           title="A  full-GR polar QNM in the complex plane (GW-damped, Imω≠0)")
# Cowling f-mode: real (1/τ = 0)
scatter!(axA, [f_cowling], [0.0], color=:gray, markersize=15, marker=:diamond,
         label=@sprintf("Cowling f-mode (τ=∞): %.2f kHz", f_cowling))
# AK prediction (real axis, context)
vlines!(axA, [f_ak], color=:seagreen, linestyle=:dot, label=@sprintf("AK1998 f: %.2f kHz", f_ak))
# full-GR f-mode: finite GW damping, BELOW Cowling
scatter!(axA, [qnm.f_kHz], [1/qnm.tau_s], color=:crimson, markersize=15,
         label=@sprintf("full-GR f-mode: %.2f kHz, τ=%.0f ms", qnm.f_kHz, qnm.tau_s*1e3))
axislegend(axA, position=:lt, framevisible=true, labelsize=10)
ylims!(axA, -2, 14)

# Panel B — interior LD metric response (fluid resonance structure)
axB = Axis(fig[1,2], xlabel="frequency f  [kHz]", ylabel="interior |H₁(R)|  (log)",
           yscale=log10, title="B  Lindblom–Detweiler interior metric response")
lines!(axB, fs, resp, color=:navy, linewidth=1.8)
vlines!(axB, [f_cowling], color=:gray, linestyle=:dash, label="Cowling f-mode")
vlines!(axB, [qnm.f_kHz], color=:crimson, linestyle=:dash, label="full-GR QNM")
axislegend(axB, position=:lb, framevisible=true, labelsize=10)

title_str = @sprintf("BDNKStar — STAGE 3: FULL-GR polar (even-parity) QNMs of the M=%.2f M⊙, R=%.2f km EOS1 star — Cowling DROPPED, mode now GW-damped (Lindblom–Detweiler interior + Zerilli exterior)", M/MSUN, R)
Label(fig[0,:], title_str, fontsize=11, font=:bold)

save(joinpath(outdir, "polar_gr_qnm.png"), fig)
println("saved figures/polar_gr_qnm.png")
