using Test
using BDNKStar

# STAGE 2 node s2.is_contrast — Israel–Stewart causal bulk-viscous closure,
# the 1D-reproducible physics of Chabanov–Rezzolla (arXiv:2311.13027).
# Verbatim targets from progress/reproduction/2311.13027.md.
@testset "IS bulk viscosity (Chabanov–Rezzolla 2311.13027): causal limiter + Reynolds limits" begin

    # --- §3 inverse-Reynolds-number limits (ideal EOS p = κ ρ^Γ) -------------
    # Verbatim conditions: Γ_id = 2, κ = 100, ρ = 0.00128, σ = −0.9
    # Verbatim targets:    R⁻¹_min ≃ −0.1 ,  R⁻¹_max ≃ 0.8
    # NB: this validates the analytic ideal-EOS *estimate* against the paper's
    # stated numbers, not a dynamical IS evolution measurement.
    Γ, κ, ρ, σ = 2.0, 100.0, 0.00128, -0.9
    p = κ * ρ^Γ
    e = ρ + p / (Γ - 1)                      # e = ρ(1+ε), ε = p/((Γ−1)ρ)
    Rmin = reynolds_inv_min(σ, p, e)
    Rmax = reynolds_inv_max(p, e)
    @info "Chabanov–Rezzolla §3 Reynolds limits" Rmin Rmax target_min=-0.1 target_max=0.8
    @test isapprox(Rmin, -0.1; atol=0.01)    # computed −0.0917 vs paper ≃ −0.1
    @test isapprox(Rmax,  0.8; atol=0.01)    # computed  0.7962 vs paper ≃  0.8
    @test -1.0 ≤ Rmin ≤ 0.0 ≤ Rmax ≤ 1.0     # both within the physical band

    # --- causal limiter (iii) self-consistency: Eq.31 → Eq.30 gives c_max² ----
    cs2_eq, ζ, ρc, hp, cmax2 = 0.2, 1.0e-3, 0.0011, 1.23, 0.9
    τbad = ζ / 5                              # ζ/τ_Π = 5 ⇒ grossly superluminal naive
    cs2_naive = cs2_viscous(cs2_eq, ζ, τbad, ρc, hp)
    @test cs2_naive > 1.0                     # naive branch IS superluminal
    (τf, cs2f, fixed) = apply_causality_fix(cs2_eq, ζ, τbad, ρc, hp; cmax2=cmax2)
    @test fixed
    @test abs(cs2f - cmax2) < 1e-14           # FP-exact: limiter IS the algebraic inverse
    # direct route matches the wrapper
    τf2 = tauPi_causal(cmax2, cs2_eq, ζ, ρc, hp)
    @test cs2_viscous(cs2_eq, ζ, τf2, ρc, hp) ≈ cmax2 atol=1e-14
    @test τf ≈ τf2
    # PHYSICAL properties NOT guaranteed by the algebraic inverse (the checks that matter):
    @test cs2f < cs2_naive                    # the fix LOWERS the sound speed
    @test 0 ≤ cs2f < 1                        # ...into the physical subluminal band
    @test τf > τbad                           # limiter "locally increases τ_Π" (eq:causal_limit)
    @test cs2_viscous(cs2_eq, ζ, 2τbad, ρc, hp) <
          cs2_viscous(cs2_eq, ζ, τbad, ρc, hp) # c_s'^2 strictly decreasing in τ_Π (ζ>0)

    # --- limiter is a no-op when the naive sound speed is already subluminal --
    τbig = 5.0                                # ζ/τ_Π = 2e-4 ⇒ naive ≈ 0.35 < c_max²
    (τu, cs2u, fixed2) = apply_causality_fix(cs2_eq, ζ, τbig, ρc, hp; cmax2=cmax2)
    @test !fixed2
    @test τu == τbig
    @test cs2u < cmax2

    # --- limiter guards: ill-posed (c_s²≥c_max²), equal boundary, superluminal target
    @test_throws ErrorException tauPi_causal(0.1, 0.2, ζ, ρc, hp)   # c_max² < c_s²_eq
    @test_throws ErrorException tauPi_causal(0.2, 0.2, ζ, ρc, hp)   # equal boundary (strict >)
    @test_throws ErrorException tauPi_causal(1.5, 0.2, ζ, ρc, hp)   # superluminal target c_max²≥1

    # --- negative-Π regime (limiter (i): Π=σp): h' must stay positive ---------
    # (branches (i) Π=σp and (ii) Π=e−p are out of scope for this node; only the
    #  causal limiter (iii) is implemented — we just guard the sign-flip risk.)
    hpm = enthalpy_prime(1.23, σ * p, ρ)      # h' = h + Π/ρ with physical Π=σp<0
    @test hpm > 0
    @test isfinite(cs2_viscous(cs2_eq, ζ, 1.0, ρ, hpm))

    # --- constitutive / conserved-variable identities + properties ------------
    @test bulk_pressure_NS(1e-3, 0.05) == -1e-3 * 0.05            # Π_NS = −ζΘ
    @test enthalpy_prime(1.2, -0.01, 1e-3) == 1.2 + (-0.01) / 1e-3 # h' = h + Π/ρ
    @test reynolds_inv(0.5, 2.0, 3.0) == 0.5 / 5.0                # R⁻¹ = Π/(p+e)
    @test reynolds_inv(σ * p, p, e) ≈ reynolds_inv_min(σ, p, e)   # the two code paths agree
    @test DPi_conserved(2.0, 1.5, -0.3) == 2.0 * 1.5 * (-0.3)     # DΠ = ρWΠ
    @test DPi_conserved(ρ, 1.005, 0.0) == 0.0                     # Π_atmo=0 ⇒ DΠ=0
    @test sign(DPi_conserved(ρ, 1.005, -3e-5)) == -1             # DΠ tracks sign of Π
end
