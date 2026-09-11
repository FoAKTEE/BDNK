using Test
using BDNKStar
using .BDNKStar: TOV, Transport, Causality

# STAGE 3, viscous extension: the BDNK first-order shear+bulk dissipative stress
# (−2ησ + ζθ) added to the 3+1D Cowling evolution. Validates that (a) the BDNK
# frame is causal at stellar conditions, and (b) viscosity damps the ℓ=2 mode
# (the HMNS-viscosity effect). Quantitative 1/τ∝η̂ scaling lives in
# viz/bdnk_viscous_damping.jl.
@testset "BDNK viscous non-radial modes (Cowling 3+1D + first-order viscosity)" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100 * 0.00128^2

    # (1) CAUSALITY: BDNK characteristic speeds real, non-negative, subluminal
    stov = TOV.solve_tov(eos, εc; h=2e-4); ec = stov.ε[1]
    tc = Transport.TransportCoefficients(η=0.16*ec, ζ=0.0, κQ=0.0,
                                         τε=2.0, τP=2.0, τQ=1.5, L=stov.R)
    fl = Causality.causality_flag(pressure(eos,ec), ec, sound_speed2(eos,ec), tc)
    @test fl.real_speeds && fl.nonneg && fl.subluminal && fl.causal

    # (2) VISCOSITY DAMPS the ℓ=2 mode: a finite η̂ run damps faster than inviscid
    s = build_star3d(eos, εc; N=22, Lfac=1.2)
    function damp(η̂)
        e = setup_evo3d(s; σ_ko=0.006, η̂=η̂)
        st = EvolState(s.grid.N); seed_l2!(st, e; A=1e-3)
        ts, q2, _ = evolve3d!(st, e; dt=0.16*s.grid.dx, nsteps=1600, sample=4)
        (damping_rate(ts, q2), all(isfinite, q2), maximum(abs, q2))
    end
    d0, ok0, m0 = damp(0.0)            # inviscid: numerical-floor damping only
    dv, okv, mv = damp(0.16)          # viscous: extra physical damping
    @test ok0 && okv                  # both stable / finite
    @test mv < 10 * m0                # bounded
    # GUARD: damping_rate returns NaN when the record holds <3 peaks (mode too
    # short-lived or overdamped to measure). Assert both are real fits, else the
    # comparison below would silently be decided by an unmeasurable run.
    @test isfinite(d0) && isfinite(dv)
    @test dv > d0                     # viscosity increases the mode damping rate
    @test d0 ≥ 0                      # inviscid engine intact (viscous block inert at η̂=0)
end
