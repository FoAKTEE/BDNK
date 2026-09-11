using Test
using BDNKStar

# STAGE 3 — POLAR (even-parity) VISCOUS quasi-normal modes, the even-parity viscous
# counterpart of the repo's axial viscous-QNM track. Frequency-domain matrix eigenvalue
# problem: eigenvalues λ of the 1D ℓ-reduced linearised polar BDNK–Cowling operator give
# the QNMs (freq = |Im λ|, damping γ = −Re λ); modes are tracked by complex continuation
# from the resolution-clean η̂=0 ideal modes. Key physics: the polar fluid already carries
# the f/p modes, so viscosity DAMPS them (no η-mode-like viscosity-driven family); γ rises
# with η̂ and p₁ damps more than f. This solver is SEMI-QUANTITATIVE (precision ~tens-of-%,
# valid Nr≲200) — we therefore test the ROBUST qualitative facts, not precise rates.
@testset "Polar viscous QNM (even-parity Cowling): viscosity damps the f/p modes" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100 * 0.00128^2
    bench,_,_ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=3, N=3000, nscan=700)

    # IDEAL limit (η̂=0): the two tracked modes are f, p₁ — near-non-dissipative (γ≈0)
    m0 = polar_qnm(eos, εc; l=2, η̂=0.0, Nr=96, nmodes=2)
    f0, p10 = qnm_freq_kHz(m0[1]), qnm_freq_kHz(m0[2])
    @test 1.7 < f0 < 2.4                                   # f-mode (conservative frame, ~+8% vs shooting)
    @test 3.3 < p10 < 4.6                                  # p₁
    @test abs(f0 - bench[1])/bench[1] < 0.12               # within ~12% of the shooting eigensolver
    @test abs(qnm_damping(m0[1])) < 5e-3                   # η̂=0 ⇒ essentially non-dissipative
    @test abs(qnm_damping(m0[2])) < 5e-3

    # VISCOUS (η̂=0.03): positive damping, larger than ideal; p₁ damps more than f
    mν = polar_qnm(eos, εc; l=2, η̂=0.03, Nr=96, nmodes=2, nstep=6)
    @test qnm_damping(mν[1]) > qnm_damping(m0[1]) + 3e-3   # f-mode acquires damping
    @test qnm_damping(mν[1]) > 5e-3                        # clearly damped
    @test qnm_damping(mν[2]) > qnm_damping(mν[1])          # p₁ damps MORE than f (higher k)

    # monotonic γ(η̂) for the f-mode (continuation makes this resolution-robust in-range)
    γ = [qnm_damping(polar_qnm(eos, εc; l=2, η̂=η̂, Nr=96, nmodes=1, nstep=6)[1]) for η̂ in (0.0,0.01,0.02,0.03)]
    @test issorted(γ)                                      # damping increases with viscosity

    # resolution-stability WITHIN the validated range (Nr 96↔140): γ_f same sign, factor≲2
    γ140 = qnm_damping(polar_qnm(eos, εc; l=2, η̂=0.03, Nr=140, nmodes=1, nstep=6)[1])
    @test γ140 > 5e-3 && 0.5 < γ140/qnm_damping(mν[1]) < 2.0
end
