#=
    SECONDARY TEST — linearity / amplitude onset of the SphBDNK mode machinery.

    The mode work (seed_sphbdnk_n!, mode_energy, evolve_sphbdnk!) is LINEAR perturbation
    theory. We verify the regime is genuinely linear and locate the nonlinearity onset by
    scanning the seed amplitude A over [1e-5 .. 1e-1] for a fixed radial-overtone ℓ=2 seed.

    Diagnostics per A (all from the SAME source: SphBDNK):
      (1) mode_energy E(t).  Quadratic in the perturbation ⇒ E ∝ A² in the linear regime.
          Report E_peak/A² and E_end/A²; a CONSTANT (A-independent) ratio ⇒ linear.
      (2) Oscillation FREQUENCY from the ℓ=2 quadrupole moment q2(t) (LINEAR in A):
          peak of the periodogram. A-independent ⇒ linear; a SHIFT ⇒ nonlinear detuning.
      (3) Nonlinearity onset: where E/A² departs from constant, where the frequency shifts,
          and where 2nd-harmonic power (P(2f)/P(f) of q2/A) grows (harmonic generation is the
          hallmark of quadratic nonlinearity; pure linear theory cannot make harmonics).

    NOTE: the seed sets δρ only (velocities start at 0); E builds up over the first
    oscillations, so we report the PEAK energy and the late-time energy, both ÷A².
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using Printf

# Same star as viz/sphbdnk_viscous.jl (Shum Γ=2 polytrope), inviscid (η̂=0) so the
# only A-dependence we probe is the perturbation nonlinearity, not transport scaling.
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
s = build_sphstar(eos, εc; Nr=64, Nθ=12)

const NSEED = 3                       # fixed radial overtone (ℓ=2, n=3) as requested
const DT    = 0.005
const TEND  = 600.0
const NSTEP = round(Int, TEND/DT)
const SAMP  = 40

# frequency grid (geometric cyclic freq, length⁻¹). f-mode band for this star is O(1e-2).
νs = range(1e-4, 0.08; length=4000)

function run_amp(A)
    e  = setup_sphbdnk(s; η̂=0.0, σ_ko=0.02)
    st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, e, NSEED; A=A)
    E0 = mode_energy(st, e)                              # seed energy (δρ only)
    ts, q2, en = evolve_sphbdnk!(st, e; dt=DT, nsteps=NSTEP, sample=SAMP)
    # peak / end energy
    Epk = maximum(en); Eend = en[end]
    # fundamental frequency from q2 (linear observable)
    P  = periodogram(ts, q2, νs)
    kpk = argmax(P); f1 = νs[kpk]; Pf1 = P[kpk]
    # 2nd-harmonic power: nearest periodogram bin to 2*f1
    k2  = argmin(abs.(collect(νs) .- 2f1)); Pf2 = P[k2]
    harm = Pf2 / Pf1
    return (; A, E0, Epk, Eend, f1, Pf1, harm,
              q2amp = maximum(abs.(q2 .- sum(q2)/length(q2))))
end

As = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
res = [run_amp(A) for A in As]

# Reference normalisation from the SMALLEST amplitude (deepest in linear regime).
r0 = res[1]
f_ref      = r0.f1
EpkA2_ref  = r0.Epk  / r0.A^2
EendA2_ref = r0.Eend / r0.A^2

println("="^104)
@printf("SphBDNK amplitude-linearity scan — ℓ=2 n=%d seed, η̂=0, Nr=%d Nθ=%d, T=%g, dt=%g\n",
        NSEED, s.grid.Nr, s.grid.Nθ, TEND, DT)
println("="^104)
@printf("%-9s %-13s %-13s %-12s %-11s %-11s %-11s\n",
        "A", "Epeak/A^2", "Eend/A^2", "f1[geom]", "f1/f_ref", "P2/P1", "(E/A2)/ref")
println("-"^104)
for r in res
    EpkA2  = r.Epk  / r.A^2
    EendA2 = r.Eend / r.A^2
    @printf("%-9.0e %-13.5e %-13.5e %-12.6f %-11.6f %-11.3e %-11.5f\n",
            r.A, EpkA2, EendA2, r.f1, r.f1/f_ref, r.harm, EpkA2/EpkA2_ref)
end
println("-"^104)

# kHz conversion for the reference fundamental (for physical sanity).
@printf("reference fundamental: f_ref = %.6f geom  = %.4f kHz\n",
        f_ref, freq_kHz_cyclic(f_ref))
println()

# ---- automated verdicts -------------------------------------------------------
println("LINEARITY CHECKS (relative to A=1e-5 reference):")
for r in res
    dE   = abs(r.Epk/r.A^2 / EpkA2_ref - 1)        # fractional departure of E/A^2
    dfr  = abs(r.f1/f_ref - 1)                     # fractional frequency shift
    flag_E = dE  > 0.05 ? "  <-- E/A^2 departs >5%"   : ""
    flag_f = dfr > 0.01 ? "  <-- freq shift >1%"      : ""
    flag_h = r.harm > 1e-2 ? "  <-- 2nd harmonic >1%" : ""
    @printf("  A=%.0e : dE/A2=%7.2f%%  df=%7.3f%%  P2/P1=%8.2e %s%s%s\n",
            r.A, 100dE, 100dfr, r.harm, flag_E, flag_f, flag_h)
end
println()
println("scan complete.")
