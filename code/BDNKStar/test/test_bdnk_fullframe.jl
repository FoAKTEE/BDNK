using Test
using BDNKStar

# STAGE 3 — FULL-FRAME (causal, hyperbolic) BDNK: the 8-field conserved+primitive
# scheme with frame relaxation times τ_ε,τ_P,τ_Q and the diagonal time-derivative
# recovery. (Quantitative spectrum lives in viz/bdnk_fullframe.jl; the absolute
# f-mode is shifted by the core excision used to tame surface-recovery stiffness.)
#
# THIS TEST WAS REWRITTEN AFTER IT WAS FOUND TO PASS/FAIL FOR THE WRONG REASONS.
# Three separate defects, all measured, none of them in the physics of the scheme:
#
# (1) SENTINEL COLLISION. `damping_rate` used to return 0.0 when it found fewer
#     than 3 peaks in |q₂| — indistinguishable from a genuine "no damping". At the
#     original nsteps=1200 even the η̂=0 run yields only 2 peaks, so the damping
#     comparison was decided by a sentinel, not by physics. Fixed on both sides:
#     `damping_rate` now returns NaN when it cannot measure, and this test both
#     lengthens the record until each run rings ≥3 times and asserts `isfinite`.
#
# (2) THE TREND WAS READ WHERE THERE IS NO MODE. At η̂=2 (τ̂=8) the quadrupole is
#     OVERDAMPED: it stops ringing (1 peak) and leaves a slowly relaxing tail, so
#     its envelope decays SLOWER than at η̂=0 and "the damping rate of the f-mode"
#     is not a well-posed observable there. The trend is now read at η̂=0.5, where
#     both runs ring cleanly (7 and 9 peaks) and separate by ~45%. η̂=0.25 is NOT
#     usable — its rate sits below the η̂=0 one, inside the Kreiss–Oliger floor
#     that dominates both.
#
# (3) THE BOUNDEDNESS CLAIM AT η̂=2 IS FALSE; THE OLD WINDOW JUST STOPPED EARLY.
#     max|q₂| at N=22, τ̂=8, dt=0.08dx:
#         η̂     n=1200     n=2000     n=3000     n=4000
#         0.0   7.50e-04   7.50e-04   7.50e-04   7.50e-04
#         2.0   7.50e-04   1.13e-02   7.48e+00   2.12e+03
#     The divergence is PRE-EXISTING (it is worse on the pre-lapse-fix code:
#     1.6e+05 at n=4000) and is not removed by any excision setting. Scanning η̂
#     at n=4000 gives max|q₂|/max|q₂|(η̂=0) = 1.000 for η̂ ≤ 1.0 and 1.6e+03 at
#     η̂=1.5 — so boundedness is asserted at η̂=1.0, where it genuinely holds, and
#     the η̂=2 case is recorded as @test_broken rather than hidden.
#
# SCOPE NOTE ON "THE HYPERBOLIC ADVANTAGE": the old comment here claimed the
# parabolic NS-limit engine goes nonfinite at η̂=2. AT THIS FILE'S dt IT DOES NOT —
# dt=0.08dx is well inside the NS diffusive CFL, and that engine stays finite to at
# least η̂=4 (its damping rises 2.87e-2 → 5.55e-2 from η̂=0 → 0.5, so its viscosity
# is genuinely active). The advantage is real but it is about the TIMESTEP, and it
# only shows at a dt the parabolic scheme cannot take: at dt=0.16dx (N=28) the NS
# engine overflows to 1.2e302 within 1400 steps while BDNK stays finite. That is
# what viz/bdnk_fullframe.jl exercises. This unit test runs at the smaller dt, so
# it does NOT test the hyperbolic advantage and no longer claims to.
@testset "Full-frame causal BDNK (3+1D Cowling): hyperbolic stability + viscous damping" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100 * 0.00128^2
    s = build_star3d(eos, εc; N=22, Lfac=1.2)

    # CAUSALITY (BDN condition a): τ_Q w₀ − η > 0 over the evolved core
    eη = setup_bdnk3d(s; η̂=2.0, ζ̂=0.0, τ̂=8.0, σ_ko=0.03)
    @test bdnk_causal_denominator(eη) > 0

    function damp(η̂)
        e = setup_bdnk3d(s; η̂=η̂, ζ̂=0.0, τ̂=8.0, σ_ko=0.03)
        st = BDNKState(s.grid.N); seed_bdnk_l2!(st, e; A=1e-3)
        ts, q2, _ = evolve_bdnk3d!(st, e; dt=0.08*s.grid.dx, nsteps=4000, sample=4)
        (damping_rate(ts, q2), all(isfinite, q2), maximum(abs, q2))
    end
    d0, ok0, m0 = damp(0.0)
    dv, okv, _  = damp(0.5)            # oscillatory ⇒ the damping rate is meaningful
    _,  ok1, m1 = damp(1.0)            # largest η̂ that stays bounded
    _,  ok2, m2 = damp(2.0)            # diverges — see (3) above

    @test ok0 && okv && ok1 && ok2     # finite (the η̂=2 growth is bounded, not NaN)

    # VISCOUS DAMPING, measured where the mode still oscillates
    @test isfinite(d0) && isfinite(dv) # GUARD: real fits (damping_rate → NaN if <3 peaks)
    @test dv > d0

    # BOUNDEDNESS: holds up to η̂=1, fails at η̂≥1.5 (known, pre-existing)
    @test m1 < 5 * m0
    @test_broken m2 < 5 * m0
end
