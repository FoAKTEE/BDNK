using Test
using BDNKStar
using BDNKStar.Transport: shum_frame_speeds, shum_frame_wellposed, TransportCoefficients
using BDNKStar.Causality: characteristic_speeds, causality_flag
using LinearAlgebra: eigvals

# ===========================================================================
# CAUSALITY PHASE-DIAGRAM + RECONCILIATION test.
# Mirrors repro/bdnk_phase_diagram.jl and repro/causality_reconcile.jl.
#
#   (a) RECONCILIATION.  The Causality.jl biquadratic (Keeble–Redondo-Yuste
#       port) is a DIFFERENT non-dimensionalisation than the Kovtun large-k
#       limit; the documented relationship is asserted, NOT a forced numeric
#       match.  The genuine characteristic-speed = k→∞ slope identity is the
#       Kovtun eq.4.18 route, which we verify against the full sound quartic.
#   (b) PRODUCTION SHUM frame (ŝ,â,q̂)=(1,1,0.999) is inside the causal+stable
#       region: 0 < q̂ < ŝ (well-posed) and c₊ ≤ 1 (causal) at stellar cs.
#   (c) The M=1.4 production-frame Shum star is causal+well-posed at every
#       radius (0 superluminal radii, speeds real, 0<q̂<ŝ throughout).
# Fast: a single small TOV solve, all else analytic.
# ===========================================================================

const ŝP, âP, q̂P = 1.0, 1.0, 0.999
const η̂P, ζ̂P     = 0.01, 0.01

# --- Kovtun sound-channel large-k machinery (companion-matrix roots) --------
# The full rest-frame sound quartic c4 ω⁴+c3 ω³+c2 ω²+c1 ω+c0 (Kovtun
# 1907.08191).  Its large-k slopes ω/k → c∞ ARE the characteristic speeds.
function _kovtun_largek_speed2(; cs, ε1, ε2, θ, π1, γs=1.0, w0=1.0, k=1e7)
    cs2 = cs^2; k2 = k^2
    c4 = cs2*ε1*θ
    c3 = im*w0*(cs2*ε1 + θ)
    c2 = -(w0^2 + k2*cs2*(cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1))
    c1 = -im*k2*w0*(γs + cs2^2*ε1 + cs2*θ)
    c0 = k2*cs2*(w0^2 + k2*θ*(cs2*(ε2+π1 - cs2*ε1) - γs))
    a = ComplexF64[c0, c1, c2, c3, c4] ./ c4
    C = zeros(ComplexF64, 4, 4); C[2,1]=1; C[3,2]=1; C[4,3]=1
    for i in 1:4; C[i,4] = -a[i]; end
    rts = eigvals(C)
    return maximum((real(z)/k)^2 for z in rts)        # max c∞²
end

# Kovtun eq.4.18 ANALYTIC large-k biquadratic (dimensionless bars). Its largest
# root c2p is the analytic c∞².
function _kovtun_eq418_c2p(; cs, ε1, ε2, θ, π1, γs=1.0)
    cs2 = cs^2
    ε̄1 = cs2*ε1/γs; ε̄2 = ε2/γs; θ̄ = θ/γs; π̄1 = π1/γs
    a = ε̄1*θ̄/cs2
    b = ε̄1*(ε̄2+π̄1-ε̄1-1/cs2) - θ̄*(ε̄2+π̄1) - ε̄2*π̄1
    c = θ̄*(cs2*(ε̄2+π̄1-ε̄1) - 1)
    disc = b^2 - 4a*c
    return (-b + sqrt(disc))/(2a)
end

@testset "Causality map (a): RECONCILIATION — char-speed = Kovtun large-k limit" begin
    # The genuine identity (characteristic speed = k→∞ slope of ω(k)) is carried
    # by the Kovtun eq.4.18 route: the analytic large-k biquadratic equals the
    # numeric large-k slopes of the full sound quartic, to high precision, across
    # several frames. (This is the "reconciled" relationship the note documents.)
    cs = 0.5; γs = 1.0; w0 = 1.0
    # Kovtun frames (ε1,ε2,θ,π1) built from Shum-style ratios ŝ=cs²ε1/γs, â=θ/γs.
    frames = [(ε1=1*γs/cs^2, ε2=0.0, θ=1*γs, π1=8*γs),   # prod-like ŝ=â=1, p̂=8
              (ε1=3*γs/cs^2, ε2=0.0, θ=4*γs, π1=8*γs),   # ŝ=3 â=4 p̂=8
              (ε1=3*γs/cs^2, ε2=0.0, θ=15*γs, π1=8*γs)]  # ŝ=3 â=15 p̂=8
    for fc in frames
        c2_num = _kovtun_largek_speed2(; cs=cs, ε1=fc.ε1, ε2=fc.ε2, θ=fc.θ,
                                        π1=fc.π1, γs=γs, w0=w0)
        c2_an  = _kovtun_eq418_c2p(; cs=cs, ε1=fc.ε1, ε2=fc.ε2, θ=fc.θ,
                                    π1=fc.π1, γs=γs)
        @test isapprox(c2_num, c2_an; rtol=1e-4)         # numeric ≡ analytic
    end

    # The DOCUMENTED SPLIT: feeding dimensionful relaxation times into the
    # Keeble–Redondo-Yuste Causality.jl biquadratic does NOT reproduce Kovtun's
    # large-k limit (different physical char-system + non-dimensionalisation).
    # We assert the split is real: its trace c₊²+c₋² = cs²(1/τQ+τP/τε+1), which
    # is NOT the Kovtun large-k sum, and is NOT τ-rescale invariant (so the
    # speeds are not k→∞ slopes as fed) — exactly the documented behaviour.
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    p = pressure(eos, εc); cs2 = sound_speed2(eos, εc); csc = sqrt(cs2); ρ = εc + p
    V̂ = (4/3)*η̂P + ζ̂P; L = 1.0
    η = q̂P*L*cs2*ρ*η̂P; ζ = q̂P*L*cs2*ρ*ζ̂P
    τε = V̂*L; τp = ŝP*cs2*L*V̂; τQ = âP*L*V̂
    tc  = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=τε, τP=τp, τQ=τQ, L=L)
    tc2 = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=2τε, τP=2τp, τQ=2τQ, L=L)
    c2m,  c2p,  _ = characteristic_speeds(p, εc, csc, tc)
    c2m2, c2p2, _ = characteristic_speeds(p, εc, csc, tc2)
    @test isapprox(c2p + c2m, cs2*(1/τQ + τp/τε + 1); rtol=1e-8)   # trace identity
    @test !isapprox(c2p, c2p2; rtol=1e-3)                          # NOT τ-rescale invariant
end

@testset "Causality map (b): production frame inside causal+stable region" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    csc = sqrt(sound_speed2(eos, εc))
    # well-posed / linearly stable (Shum eq.71): 0 < q̂ < ŝ
    @test shum_frame_wellposed(ŝP, âP, q̂P)
    @test 0 < q̂P < ŝP
    # causal: c₊ ≤ 1 at the stellar-centre sound speed
    c0, cp, cm = shum_frame_speeds(ŝP, âP, q̂P, η̂P, ζ̂P, csc)
    @test cp ≤ 1.0
    @test isapprox(cp, sqrt(3.0)*csc; rtol=1e-3)        # c₊ = √3 cs
    @test isfinite(c0) && c0 ≥ 0 && cm ≥ 0              # all modes real, non-neg
    # the WP boundary q̂=ŝ is the binding constraint (margin ~1e-3 by design)
    @test (ŝP - q̂P) > 0 && (ŝP - q̂P) < 1e-2
end

@testset "Causality map (c): M=1.4 Shum star causal+well-posed at every radius" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    star = solve_tov(eos, εc; h=5e-4, ptol_rel=1e-10, rmax=50.0)
    @test isapprox(star.M, 1.4; atol=0.02)
    n_superlum = 0; allreal = true
    for ε in star.ε
        ε > 0 || continue
        cs = sqrt(max(sound_speed2(eos, ε), 0.0))
        c0, cp, cm = shum_frame_speeds(ŝP, âP, q̂P, η̂P, ζ̂P, cs)
        allreal &= isfinite(c0) && isfinite(cp) && isfinite(cm) && c0≥0 && cp≥0 && cm≥0
        cp > 1 && (n_superlum += 1)
    end
    @test allreal                                       # speeds real ∀r
    @test n_superlum == 0                               # 0 superluminal radii
    @test shum_frame_wellposed(ŝP, âP, q̂P)              # well-posed ∀r (EOS-indep.)
    # caveat probe: a STIFFER EOS (cs > 1/√3) WOULD violate subluminality; here
    # max cs stays below 1/√3, which is exactly why this star is causal.
    cs_max = maximum(sqrt(max(sound_speed2(eos, ε), 0.0)) for ε in star.ε if ε > 0)
    @test cs_max < 1/sqrt(3)
end
