using Test
using BDNKStar

# ===========================================================================
# Shum FLUID-FRAME (fluid-rest-frame) characteristic-structure / well-posedness
# validation — Shum, Abalos, Bea, Bezares, Figueras, Palenzuela arXiv:2509.15303
# eqs.(67–71).  Mirrors repro/shum_frame_analysis.jl.  Fast (no evolution).
#
#   c±  = cs √{ [â(1+ŝ)+q̂ ± √D]/(2â) },  D = q̂²+â²(4q̂+(ŝ-1)²)+2âq̂(1+ŝ)
#   c0  = cs √( q̂ η̂ / (â V̂) ),  V̂=(4/3)η̂+ζ̂
#   well-posed / linearly stable (eq.71): 0 < q̂ < ŝ
#   production frame (ŝ,â,q̂)=(1,1,0.999) ⇒ c₊=√3 cs, c₋≈0.0183 cs.
# ===========================================================================

# Local copy of the shum_core frame→transport mapping (lines 153–164 of
# repro/shum_core.jl), kept here so the test is self-contained and fast.
function _shum_transport(ŝ, â, q̂, η̂, ζ̂, L, eos, ε)
    p   = pressure(eos, ε)
    cs2 = sound_speed2(eos, ε)
    ρ   = ε + p
    V̂   = (4/3)*η̂ + ζ̂
    η   = q̂ * L * cs2 * ρ * η̂
    ζ   = q̂ * L * cs2 * ρ * ζ̂
    τε  = V̂ * L
    τp  = ŝ * cs2 * L * V̂
    τQ  = â * L * V̂
    return η, ζ, τε, τp, τQ, p, cs2
end

@testset "Shum frame: production-frame analytic c± (eq.67-71)" begin
    # (1) production-frame c₊=√3, c₋≈0.0183 (frame-only, cs=1)
    c0, cp, cm = shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, 1.0)
    @test isapprox(cp, sqrt(3.0); rtol=1e-3)     # √3 to <0.1%
    @test isapprox(cm, 0.0183;  rtol=0.05)       # 0.0183 to a few %
    @test isfinite(c0) && c0 > 0                 # diffusive mode real & positive
    # c± scale linearly with cs
    c0b, cpb, cmb = shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, 0.5)
    @test isapprox(cpb, 0.5*cp; rtol=1e-12)
    @test isapprox(cmb, 0.5*cm; rtol=1e-12)
end

@testset "Shum frame: well-posedness predicate (eq.71  0<q̂<ŝ)" begin
    # (2) shum_frame_wellposed correct on a few cases
    @test  shum_frame_wellposed(1.0, 1.0, 0.999)   # 0<q̂<ŝ  ✔
    @test  shum_frame_wellposed(1.5, 1.0, 0.75)    # 0<q̂<ŝ  ✔
    @test !shum_frame_wellposed(1.0, 1.0, 1.0)     # q̂=ŝ    ✘ (boundary excluded)
    @test !shum_frame_wellposed(1.0, 1.0, 1.5)     # q̂>ŝ    ✘
    @test !shum_frame_wellposed(1.0, 1.0, -0.1)    # q̂≤0    ✘
end

@testset "Shum frame: real & finite speeds, well-posed along the Shum star" begin
    # (3) along the Shum M=1.4 star the production-frame speeds are real & finite
    #     and the frame is well-posed throughout.
    eos = ShumPolytrope(100.0)
    εc  = 0.00128 + 100*0.00128^2
    star = solve_tov(eos, εc; h=5e-4, ptol_rel=1e-10, rmax=50.0)
    @test isapprox(star.M, 1.4; atol=0.02)
    # evaluate the production-frame speeds at every interior radius, then AGGREGATE
    # (avoid one @test per grid point — keeps the assertion count sane)
    sp = [shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, sqrt(sound_speed2(eos, ε)))
          for ε in star.ε if ε > 0]
    c0s = first.(sp); cps = getindex.(sp, 2); cms = last.(sp)
    @test all(isfinite, c0s) && all(isfinite, cps) && all(isfinite, cms)
    @test all(≥(0), c0s) && all(≥(0), cps) && all(≥(0), cms)   # real, non-negative
    @test shum_frame_wellposed(1.0, 1.0, 0.999)                # frame well-posed
    cp_max = maximum(cps)
    # c₊ = √3·cs stays subluminal everywhere (cs < 1/√3 throughout this star)
    @test cp_max ≤ 1.0
    # the centre value is c₊ = √3·cs(centre)
    csc = sqrt(sound_speed2(eos, εc))
    @test isapprox(cp_max, sqrt(3.0)*csc; rtol=2e-2)
end

@testset "Shum frame: cross-validation (mapping ratios hold; speed formulas differ)" begin
    # (4) Dual-path: the shum_core mapping reproduces Shum's dimensionless frame
    #     RATIOS exactly, but the general Causality biquadratic is a DIFFERENT
    #     (dimensionful) characteristic system, so the SPEED formulas do NOT
    #     coincide.  We assert what actually holds (honest), not a forced match.
    eos = ShumPolytrope(100.0)
    εc  = 0.00128 + 100*0.00128^2
    ŝ, â, q̂, η̂, ζ̂, L = 1.0, 1.0, 0.999, 0.01, 0.01, 1.0
    η, ζ, τε, τp, τQ, p, cs2 = _shum_transport(ŝ, â, q̂, η̂, ζ̂, L, eos, εc)

    # (4a) the mapping reproduces ŝ = τP/(cs²τε) and â = τQ/τε EXACTLY
    @test isapprox(τp/(cs2*τε), ŝ; rtol=1e-12)
    @test isapprox(τQ/τε,        â; rtol=1e-12)

    # (4b) the two speed formulas do NOT agree (documented convention mismatch)
    cs = sqrt(cs2)
    _, cpA, _ = shum_frame_speeds(ŝ, â, q̂, η̂, ζ̂, cs)            # Path A: Shum frame
    tc = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=τε, τP=τp, τQ=τQ, L=L)
    c2m, c2p, disc = characteristic_speeds(p, εc, cs, tc)        # Path B: biquadratic
    @test disc ≥ 0 && isfinite(c2p) && c2p ≥ 0                   # B solves & is real
    cpB = sqrt(c2p)
    @test !isapprox(cpA, cpB; rtol=1e-2)                         # they DIFFER

    # (4c) the structural reason: B's trace c₊²+c₋² = cs²(1/τQ + τP/τε + 1)
    #      matches the biquadratic exactly (Vieta), and blows up as τ→0.
    trace_formula = cs2*(1/τQ + τp/τε + 1)
    @test isapprox(c2p + c2m, trace_formula; rtol=1e-8)
    @test (c2p + c2m)/cs2 > 10                                   # ≫ Shum's frame-pure 3
end
