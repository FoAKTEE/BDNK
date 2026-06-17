using Test
using BDNKStar

# STEP 0 GATE: causality monitor operative. We validate the *machinery* (the
# characteristic-speed solver and the pointwise predicate). The physical
# classification of specific named frames is cross-checked against the source
# papers in a later substage; here we assert the biquadratic is solved exactly
# and the predicate fields behave.

@testset "Causality: characteristic speeds solve the biquadratic (Vieta residual)" begin
    poly = PolytropeEnergy(100.0, 1.0)
    for e in (5e-4, 1e-3, 2e-3)
        p = pressure(poly, e); cs2 = sound_speed2(poly, e); cs = sqrt(cs2)
        tc = TransportCoefficients(η=1e-3, ζ=5e-4, κQ=1e-3, τε=2.0, τP=2.0, τQ=1.5, L=1.0)
        c2m, c2p, disc = characteristic_speeds(p, e, cs, tc)
        # roots of  Λ2 x² - 2Λ1 x + Λ0 = 0  ⇔ x = (Λ1 ± √disc)/Λ2
        λ0 = BDNKStar.Causality.Λ0(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
        λ1 = BDNKStar.Causality.Λ1(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
        λ2 = BDNKStar.Causality.Λ2(p, e, cs, tc.η, tc.ζ, tc.τε, tc.τP, tc.τQ, tc.L)
        scale = abs(λ2*c2p^2) + abs(2λ1*c2p) + abs(λ0) + eps()
        @test abs(λ2*c2m^2 - 2λ1*c2m + λ0) / scale < 1e-10
        @test abs(λ2*c2p^2 - 2λ1*c2p + λ0) / scale < 1e-10
        # Vieta sum (well-conditioned) and product (ill-conditioned: λ0/λ2 is a
        # tiny difference of large terms, so a looser tol is appropriate; the
        # residual plug-back above is the authoritative correctness check).
        @test isapprox(c2m + c2p, 2λ1/λ2; rtol=1e-10)
        @test isapprox(c2m * c2p, λ0/λ2; rtol=1e-6)
    end
end

@testset "Causality: Shum frame characteristic speeds (2509.15303 eq.67-71)" begin
    # production frame (ŝ,â,q̂)=(1,1,0.999); ratios c±/cs are frame-only
    c0, cp, cm = shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, 1.0)
    @test isapprox(cp, sqrt(3.0); atol=1e-3)     # c₊ = √3 cs
    @test isapprox(cm, 0.0183; atol=1e-3)        # c₋ = 0.0183 cs
    @test shum_frame_wellposed(1.0, 1.0, 0.999)  # 0 < q̂ < ŝ
    @test !shum_frame_wellposed(1.0, 1.0, 1.5)   # q̂ > ŝ violates well-posedness
end

@testset "Causality: predicate fields and subluminal detection" begin
    poly = PolytropeEnergy(100.0, 1.0)
    e = 1e-3; p = pressure(poly, e); cs2 = sound_speed2(poly, e)
    tc = TransportCoefficients(η=1e-3, ζ=5e-4, κQ=1e-3, τε=2.0, τP=2.0, τQ=1.5, L=1.0)
    f = causality_flag(p, e, cs2, tc)
    @test f isa NamedTuple
    @test isfinite(f.c2_plus) && isfinite(f.c2_minus)
    @test (f.subluminal == (f.c2_plus ≤ 1 + 1e-12))
    @test (f.causal == (f.real_speeds && f.nonneg && f.subluminal))
    # conformal reference frame evaluates without error
    tcc = conformal_frame_PMP(e)
    fc = causality_flag(p, e, cs2, tcc)
    @test isfinite(fc.disc)
end
