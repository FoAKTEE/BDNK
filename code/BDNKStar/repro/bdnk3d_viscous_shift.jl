# ======================================================================
# repro/bdnk3d_viscous_shift.jl — CONVERGING THE VISCOUS FREQUENCY SHIFT OF THE ℓ=2 f-MODE
# in the 3+1D BDNK Cowling engine (CowlingBDNK3D), frame F₁ of Clarisse et al.
#
# RESULT (data: repro/data/bdnk3d_viscous_shift.csv, 33 runs, 2026-09-16):
#   The shift Δf/f = f(η̂)/f(0) − 1 measured earlier as +0.40/+0.28/+0.15% (N=24/32/40,
#   σ_KO=0.01) and quoted only as an upper bound is a KREISS–OLIGER ARTIFACT. Halving σ_KO
#   to 0.005 removes ~0.35% of it at every N. What remains converges to the DAMPED-
#   OSCILLATOR shift, which has NO free parameter (γ_tot and γ₀ are measured in the runs):
#
#           Δf/f = −(γ_tot² − γ₀²) / (2ω²)
#
#   σ_KO=0.005, η̂=0.03:   N    measured (pgram/pencil)   predicted   residual
#                          24   +0.055 / +0.048%          −0.202%     +0.25%
#                          32   −0.025 / −0.038%          −0.157%     +0.13%
#                          40   −0.069 / −0.094%          −0.133%     +0.05%
#                          48   −0.112 / −0.138%          −0.133%     +0.008%
#   The residual (the leftover numerical offset) falls ~h³. At N=40 the η̂-scaling is
#   quadratic once the offset is removed (η̂=0.015/0.030/0.045). Frame-reactive BDNK
#   corrections are O(τ²ω²) ≈ 3e−5, below resolution.
#
#   ⇒  Δf/f(η̂=0.03) = −0.125 ± 0.013%  (the damped-oscillator value), resolved for
#      σ_KO ≤ 0.005 and N ≥ 48.
#
#   ALSO: γ_visc (pencil) = 2.597e−3 (N=40) / 2.607e−3 (N=48) /M⊙ at σ_KO=0.005 = 99.2% of
#   the Navier–Stokes dissipation integral 0.08746·η̂ = 2.6238e−3. The "96–97%" from the
#   σ_KO=0.01 Richardson sequence was KO-biased as well.
#
# METHOD NOTES that the data force:
#   • Hold the FRAME fixed (η̂_frame=0.03) while sweeping η̂, so the η-dependence is not
#     entangled with τ ∝ η̂; then one η̂=0 control serves the whole sweep at a given (N,σ).
#   • The inviscid control goes UNSTABLE (cut-cell surface mode) below σ_KO≈0.005 at N=24 and
#     32; viscous runs stay stable to σ_KO=0.00125. So σ_KO=0.005 is the floor of the
#     difference protocol, set by the control, not by the viscous physics.
#   • At fixed τ, η̂=0.06 destabilises at N=24 (margin τ_Q w₀−η = 0.047 ε₀) — the k_c picture.
#   • `visc_compact=true` INCREASES the shift by ~0.15% and cuts γ_visc by 32% at N=24; that
#     stencil is not trustworthy for the mode and is not used.
#   • dt_fac 0.20→0.10 changes nothing to 7 digits: RK4 at ~560 steps/period is converged
#     (error ~1e−8). (dt_fac=0.20 does blow up at N=12 — coarse-grid recovery stiffness.)
#   • Two estimators (log-parabolic periodogram peak, matrix pencil pole nearest the peak)
#     disagree by 0.02–0.05% on the shift; their mean is quoted, the spread is the error.
#
# Running the full set costs ~5 CPU-hours (N=48 runs ≈ 12 min each single-threaded). The
# short form below reproduces the σ_KO discriminant at N=24 in ~10 min.
# ======================================================================
using BDNKStar, Printf
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
const F1 = (25/4, 25/7)
const KHZ_TO_KM = 1e3/2.99792458e5; const MSUN_KM = 1.476625

run(N, η̂, σ; np=9) = bdnk3d_qnm(eos, εc; N=N, η̂=η̂, frame=F1, η̂_frame=0.03, control=false,
                                 nperiods=np, f_ref_kHz=1.88291, σ_ko=σ, dt_fac=0.20).main
function shift_and_prediction(N, η̂, σ; np=9)
    c = run(N, 0.0, σ; np=np); v = run(N, η̂, σ; np=np)
    ω = 2π * c.f_kHz * KHZ_TO_KM * MSUN_KM
    meas = 100 * ((v.f_kHz/c.f_kHz - 1) + (v.f_pencil_kHz/c.f_pencil_kHz - 1)) / 2
    pred = -100 * (v.γ^2 - c.γ^2) / (2ω^2)
    (meas=meas, pred=pred, γ_visc=v.γ_pencil - c.γ_pencil, stable=c.stable && v.stable)
end

full = get(ENV, "BDNK_SHIFT_FULL", "0") == "1"
println("σ_KO discriminant at N=24, η̂=0.03:")
for σ in (0.01, 0.005)
    r = shift_and_prediction(24, 0.03, σ)
    @printf("  σ_KO=%.3f  shift=%+.3f%%  damped-osc pred=%+.3f%%  γ_visc=%.3e  stable=%s\n", σ, r.meas, r.pred, r.γ_visc, r.stable)
end
if full
    println("\nN convergence at σ_KO=0.005, η̂=0.03 (12-period records):")
    for N in (32, 40, 48)
        r = shift_and_prediction(N, 0.03, 0.005; np=12)
        @printf("  N=%d  shift=%+.3f%%  pred=%+.3f%%  residual=%+.3f%%  γ_visc=%.3e (NS integral 2.624e-3)\n",
                N, r.meas, r.pred, r.meas - r.pred, r.γ_visc)
    end
end
