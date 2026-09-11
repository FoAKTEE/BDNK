using Test
using BDNKStar
using BDNKStar: Units

# ─────────────────────────────────────────────────────────────────────────────
# CONVERGENCE + NUMERICAL-vs-PHYSICAL VISCOSITY SYNTHESIS
#
# This suite distils the convergence studies of four BDNKStar solver tracks into
# fast, ROBUST regression gates.  We assert qualitative signatures (monotonicity,
# Richardson-consistency, sign, physical-vs-numerical discrimination), NOT brittle
# exact eigenvalues — because some of these operators are only semi-quantitative
# or marginally stable, and the honest finding is about the TREND, not the number.
#
# Tracks & verdicts (from the analysis agents):
#   • RModes  (LOM98 Lane–Emden integrals)        → CONVERGED, clean 2nd order
#   • AxialViscousModes (w-/η-mode shooting QNM)   → CONVERGED, grid-saturated
#   • PolarViscousModes (f/p Cowling eigenproblem) → STABLE only in Nr≲200 window
#   • SphBDNK (time-domain viscous dissipation)    → MIXED: ν_mom is a genuine sink,
#       but the measured Γ at small η̂ is contaminated by a resolution-dependent,
#       σ_ko-controlled grid-scale GROWTH of the η̂=0 baseline.  Physical sign/trend
#       only above η̂≈0.06; quantitative Γ(η̂≲0.05) is NOT trustworthy.
# ─────────────────────────────────────────────────────────────────────────────

eos_for_sph() = ShumPolytrope(100.0)
εc_for_sph()  = 0.00128 + 100 * 0.00128^2

"Least-squares exponential rate Γ of an energy time series: E(t)≈E0·exp(−Γt).
 Γ>0 ⇒ decay (physical dissipation); Γ<0 ⇒ growth (numerical floor)."
function _energy_rate(ts::AbstractVector, en::AbstractVector)
    m = (en .> 0) .& isfinite.(en)
    t = collect(ts[m]); y = log.(en[m])
    n = length(t); n < 3 && return NaN
    t̄ = sum(t)/n; ȳ = sum(y)/n
    Sxx = sum((t .- t̄).^2); Sxy = sum((t .- t̄).*(y .- ȳ))
    return -(Sxy/Sxx)            # slope of log E = −Γ, so Γ = −slope
end

"Run a short SphBDNK quadrupole-overtone evolution and return its energy rate Γ."
function _sph_gamma(s; η̂, σ_ko, T=120.0, dt=0.01)
    ev = setup_sphbdnk(s; η̂=η̂, σ_ko=σ_ko)
    st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, ev, 3; A=1e-3)
    nst = round(Int, T/dt)
    ts, _, en = evolve_sphbdnk!(st, ev; dt=dt, nsteps=nst, sample=max(1, nst÷200))
    return _energy_rate(ts, en), all(isfinite, en), en[end]/en[1]
end

@testset "Convergence + physical-vs-numerical viscosity synthesis" begin

    # ── (a) RModes: LOM98 Lane–Emden integrals converge (clean 2nd order) ──────
    # The structure integrals J̃, Ĩ and the timescales τ_GW, τ_sv settle MONOTONELY
    # toward the LOM98 Table-I values as Lane–Emden resolution N grows, with the
    # Richardson rate consistent with the trapezoid order (p≈2).
    @testset "(a) RMode LOM98 integrals: monotone, Richardson-consistent 2nd order" begin
        v1 = rmode_validate_lom98(N=1000)
        v2 = rmode_validate_lom98(N=2000)
        v4 = rmode_validate_lom98(N=4000)
        v8 = rmode_validate_lom98(N=8000)

        # benchmark agreement at the production N (LOM98 Table I)
        @test isapprox(v4.J̃,  1.635e-2; rtol=2e-3)
        @test isapprox(v4.Ĩ,  0.2614;  rtol=2e-3)
        @test isapprox(v4.τ_GW, -3.257; rtol=2e-3)
        @test isapprox(v4.τ_sv, 2.516e8; rtol=2e-3)

        # MONOTONE in N (sequence tends to a limit, here from below for J̃,Ĩ)
        @test v1.J̃ < v2.J̃ < v4.J̃ < v8.J̃
        @test v1.Ĩ < v2.Ĩ < v4.Ĩ < v8.Ĩ

        # Richardson rate p = log2(|d12| / |d24|) ≈ 2 for a 2nd-order scheme.
        for q in (:J̃, :Ĩ)
            d12 = getfield(v2,q) - getfield(v1,q)
            d24 = getfield(v4,q) - getfield(v2,q)
            p = log2(abs(d12)/abs(d24))
            @test 1.6 < p < 2.4                      # clean 2nd order (theory p=2)
        end
        # successive increments SHRINK (converging, not diverging)
        @test abs(v8.J̃ - v4.J̃) < abs(v4.J̃ - v2.J̃) < abs(v2.J̃ - v1.J̃)
    end

    # ── (b) PolarViscousModes: eigenvalue stable across the converged Nr window ─
    # The f/p₁ Cowling QNMs are NOT strictly convergent (marginally-stable operator),
    # but they ARE stable inside the documented Nr≲200 window: same sign of damping,
    # bounded drift, p₁ damps more than f, γ monotone in η̂.  We assert that window.
    @testset "(b) Polar-viscous eigenvalue: stable in the Nr≲200 window" begin
        eos = ShumPolytrope(100.0); εc = 0.00128 + 100 * 0.00128^2

        # ideal (η̂=0) f-mode frequency is resolution-STABLE across 96↔160 (≲5%)
        f96  = qnm_freq_kHz(polar_qnm(eos, εc; l=2, η̂=0.0, Nr=96,  nmodes=1)[1])
        f160 = qnm_freq_kHz(polar_qnm(eos, εc; l=2, η̂=0.0, Nr=160, nmodes=1)[1])
        @test 1.7 < f96 < 2.4 && 1.7 < f160 < 2.4
        @test abs(f160 - f96)/f96 < 0.08             # drift bounded inside the window

        # viscous (η̂=0.03) f-mode: positive damping at both resolutions, factor<2 apart
        γ96  = qnm_damping(polar_qnm(eos, εc; l=2, η̂=0.03, Nr=96,  nmodes=1, nstep=6)[1])
        γ160 = qnm_damping(polar_qnm(eos, εc; l=2, η̂=0.03, Nr=160, nmodes=1, nstep=6)[1])
        @test γ96 > 5e-3 && γ160 > 5e-3              # clearly damped (same sign)
        @test 0.5 < γ160/γ96 < 2.0                   # stable (documented factor≲2)

        # γ increases monotonically with η̂ (the physically robust differential fact)
        γη = [qnm_damping(polar_qnm(eos, εc; l=2, η̂=η, Nr=96, nmodes=1, nstep=6)[1])
              for η in (0.0, 0.01, 0.02, 0.03)]
        @test issorted(γη)
    end

    # ── (c) THE CRUX — SphBDNK: physical-vs-numerical viscous discrimination ────
    # Honest "mixed" verdict, encoded:
    #   (c1) the ν_mom shear channel is a GENUINE energy sink: at η̂ above the
    #        physical threshold the energy ratio E_end/E0 DECREASES with η̂.
    #   (c2) the η̂=0 baseline is a NUMERICAL grid-scale rate controlled by σ_ko
    #        (sloped Γ vs σ_ko), whereas adding viscosity REDUCES that σ_ko slope
    #        (the rate becomes more physical / flatter, though not perfectly flat).
    #   (c3) (Γ−Γ₀) GROWS with η̂ — viscosity adds damping over the inviscid floor.
    # Short, coarse runs: we assert SIGNS and TRENDS, never exact Γ values.
    @testset "(c) SphBDNK: viscous dissipation is physical where discrimination showed" begin
        s = build_sphstar(eos_for_sph(), εc_for_sph(); Nr=48, Nθ=10)

        # (c1) energy ratio falls with η̂ across the physically meaningful range
        function eratio(η̂)
            ev = setup_sphbdnk(s; η̂=η̂, σ_ko=0.02)
            st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, ev, 3; A=1e-3)
            _, _, en = evolve_sphbdnk!(st, ev; dt=0.01, nsteps=round(Int,120/0.01), sample=40)
            (en[end]/en[1], all(isfinite, en))
        end
        r0, ok0 = eratio(0.0)
        r6, ok6 = eratio(0.06)
        r10, ok10 = eratio(0.10)
        @test ok0 && ok6 && ok10                     # all runs finite/bounded
        @test r10 < r6 < r0                          # MORE viscosity ⇒ MORE dissipation
        @test r10 < 1.0                              # net DECAY above the η̂≈0.06 threshold

        # (c2) numerical-channel signature: Γ(η̂=0) is σ_ko-DEPENDENT (sloped);
        #      adding viscosity REDUCES that σ_ko-sensitivity (closer to flat/physical).
        Γ0_lo, _, _ = _sph_gamma(s; η̂=0.0,  σ_ko=0.005)
        Γ0_hi, _, _ = _sph_gamma(s; η̂=0.0,  σ_ko=0.04)
        Γv_lo, _, _ = _sph_gamma(s; η̂=0.08, σ_ko=0.005)
        Γv_hi, _, _ = _sph_gamma(s; η̂=0.08, σ_ko=0.04)
        slope0 = abs(Γ0_hi - Γ0_lo)                  # σ_ko-sensitivity of the bare floor
        slopev = abs(Γv_hi - Γv_lo)                  # σ_ko-sensitivity once viscous
        @test slopev < slope0                        # viscosity flattens the σ_ko dependence
        @test Γ0_lo < 0 && Γ0_hi < 0                 # η̂=0 is GROWTH (numerical floor), not decay

        # (c3) (Γ−Γ₀) grows with η̂: viscosity adds dissipation above the inviscid floor.
        Γ0,  _, _ = _sph_gamma(s; η̂=0.0,  σ_ko=0.02)
        Γ6,  _, _ = _sph_gamma(s; η̂=0.06, σ_ko=0.02)
        Γ10, _, _ = _sph_gamma(s; η̂=0.10, σ_ko=0.02)
        @test (Γ6  - Γ0) > 0                         # viscosity damps above the floor
        @test (Γ10 - Γ0) > (Γ6 - Γ0)                 # and more so at larger η̂ (monotone)
    end

    # ── (d) The numerical η̂=0 floor SHRINKS with resolution, but never vanishes ─
    # Honest caveat: the inviscid grid-scale GROWTH does weaken at higher Nr at fixed
    # σ_ko — so the contamination is real but resolution-dependent, NOT scheme-order
    # clean.  We assert the magnitude of the floor at σ_ko=0.005 decreases 48→64.
    @testset "(d) SphBDNK η̂=0 numerical floor weakens with resolution" begin
        s48 = build_sphstar(eos_for_sph(), εc_for_sph(); Nr=48, Nθ=10)
        s64 = build_sphstar(eos_for_sph(), εc_for_sph(); Nr=64, Nθ=12)
        Γ48, _, _ = _sph_gamma(s48; η̂=0.0, σ_ko=0.005)
        Γ64, _, _ = _sph_gamma(s64; η̂=0.0, σ_ko=0.005)
        @test Γ48 < 0 && Γ64 < 0                     # growth at both resolutions
        @test abs(Γ64) < abs(Γ48)                    # but the floor weakens with Nr
    end
end
