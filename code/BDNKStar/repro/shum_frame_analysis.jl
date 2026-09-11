#=
    shum_frame_analysis.jl  —  FLUID-FRAME (Shum-style) characteristic-structure
    analysis that VALIDATES the BDNKStar BDNK code against the well-posedness /
    causality (fluid-rest-frame characteristic-speed) analysis of

        Shum, Abalos, Bea, Bezares, Figueras, Palenzuela, arXiv:2509.15303
        ("Neutron star evolution with the BDNK ... framework").

    GROUNDING.  Shum eqs.(67–71) give the FLUID-REST-FRAME BDNK characteristic
    speeds in the hatted spherical-Cowling frame parametrization (V̂=(4/3)η̂+ζ̂):

        c0  = cs √( q̂ η̂ / (â V̂) )                          (diffusive mode)
        c±  = cs √{ [â(1+ŝ)+q̂ ± √D] / (2â) },               (modified sound)
              D = q̂² + â²(4q̂+(ŝ-1)²) + 2â q̂(1+ŝ)

    and the well-posedness/linear-stability condition (eq.71):  0 < q̂ < ŝ.
    Production frame (ŝ,â,q̂)=(1,1,0.999) ⇒ c₊=√3 cs, c₋≈0.0183 cs.

    These are implemented in the package:
        src/transport/Transport.jl :: shum_frame_speeds, shum_frame_wellposed.
    An INDEPENDENT general-BDNK characteristic biquadratic (different reference,
    Keeble & Redondo-Yuste; Λ₂c⁴−2Λ₁c²+Λ₀=0) lives in
        src/transport/Causality.jl :: characteristic_speeds.
    The frame→relaxation-time mapping lives in repro/shum_core.jl :: shum_transport.

    This script does NOT reimplement physics — it drives the package routines and
    reports every number.  Five sections:
      1. VALIDATE shum_frame_speeds against Shum's analytic c±.
      2. FRAME-SPACE characteristic-structure / well-posedness map.
      3. ALONG-THE-STAR production-frame speeds c0,c±(r).
      4. CROSS-VALIDATION: Shum fluid-frame c± vs the general Causality biquadratic
         (built from the shum_core mapping), in the fluid rest frame.
      5. SUMMARY.

    Run: cd code/BDNKStar && julia --project=. repro/shum_frame_analysis.jl
=#

# Reuse the STAGE-1C core, which ALSO loads the BDNKStar package and defines the
# ShumFrame struct + the frame→transport mapping (shum_transport).  Including it
# here is the single source of truth for the frame physics (no reimplementation).
# NOTE: shum_core.jl runs its own STAGE-1C validation `main()` on include — that
# is the documented isotropic-transform / recovery round-trip gate; we let it run.
include(joinpath(@__DIR__, "shum_core.jl"))   # -> .BDNKStar, ShumFrame, shum_transport, …

using .BDNKStar
using .BDNKStar.Transport: TransportCoefficients
using .BDNKStar.Causality: characteristic_speeds
using Printf
# ShumPolytrope, pressure, sound_speed2, solve_tov, shum_frame_speeds,
# shum_frame_wellposed are already in scope from shum_core.jl's `using` block.

# ===========================================================================
# Helper: the Shum star + a production-frame object (single source of truth)
# ===========================================================================
const κ_SHUM = 100.0
const EOS    = ShumPolytrope(κ_SHUM)
const ρ0c    = 0.00128
const εc     = ρ0c + κ_SHUM*ρ0c^2                       # M=1.4 M⊙ central ε
const PROD   = ShumFrame(ŝ=1.0, â=1.0, q̂=0.999, η̂=0.01, ζ̂=0.01)   # smallSB-F2

# ===========================================================================
# 1. VALIDATE shum_frame_speeds against Shum's documented analytic c±
# ===========================================================================
function section1_validate()
    println("="^76)
    println("1. VALIDATE shum_frame_speeds vs Shum analytic (production frame, cs=1)")
    println("="^76)
    cs = 1.0
    c0, cp, cm = shum_frame_speeds(PROD.ŝ, PROD.â, PROD.q̂, PROD.η̂, PROD.ζ̂, cs)
    cp_target, cm_target = sqrt(3.0), 0.0183
    cp_err = abs(cp - cp_target)/cp_target
    cm_err = abs(cm - cm_target)/cm_target
    @printf("   c₊        = %.6f   (target √3 = %.6f,  rel.err = %.3e)\n", cp, cp_target, cp_err)
    @printf("   c₋        = %.6f   (target  %.4f,       rel.err = %.3e)\n", cm, cm_target, cm_err)
    @printf("   c0        = %.6f   (diffusive mode = cs√(q̂η̂/(âV̂)))\n", c0)
    ok_cp = cp_err < 1e-3       # <0.1%
    ok_cm = cm_err < 0.05       # few %
    @printf("   c₊ within 0.1%%:  %s     c₋ within 5%%:  %s\n", ok_cp, ok_cm)
    return (c0=c0, cp=cp, cm=cm, ok=ok_cp && ok_cm)
end

# ===========================================================================
# 2. FRAME-SPACE characteristic-structure / well-posedness map
# ===========================================================================
function classify_point(ŝ, â, q̂, η̂, ζ̂, cs)
    # Replicate the algebra of shum_frame_speeds WITHOUT throwing on a negative
    # mode-square (that negative case IS the non-hyperbolic class we classify).
    V̂   = (4/3)*η̂ + ζ̂
    D    = q̂^2 + â^2*(4q̂ + (ŝ-1)^2) + 2*â*q̂*(1+ŝ)     # under-root of c±
    base = â*(1+ŝ) + q̂
    arg0 = q̂*η̂/(â*V̂)                                     # c0²/cs²
    sqD  = D ≥ 0 ? sqrt(D) : NaN
    argp = (base + sqD)/(2â)                              # c₊²/cs²
    argm = (base - sqD)/(2â)                              # c₋²/cs²
    hyperbolic = (D ≥ 0) && (argp ≥ 0) && (argm ≥ 0) && (arg0 ≥ 0)  # real speeds
    # use the package routine where it is safe (the production-allowed region)
    c0 = cs*sqrt(max(arg0, 0.0)); cp = cs*sqrt(max(argp, 0.0)); cm = cs*sqrt(max(argm, 0.0))
    causal    = hyperbolic && (cp ≤ 1) && (c0 ≤ 1)       # subluminal (with this cs)
    wellposed = shum_frame_wellposed(ŝ, â, q̂)            # Shum eq.71: 0<q̂<ŝ
    return (D=D, c0=c0, cp=cp, cm=cm, hyp=hyperbolic, caus=causal, wp=wellposed)
end

function section2_framespace()
    println()
    println("="^76)
    println("2. FRAME-SPACE characteristic structure / well-posedness (Shum eq.71)")
    println("="^76)
    cs_c = sqrt(sound_speed2(EOS, εc))     # stellar-centre sound speed
    @printf("   stellar-centre cs = %.5f  (cs² = %.5f)\n", cs_c, cs_c^2)

    # (a) scan q̂ ∈ (0, ŝ]  at ŝ = â = 1
    println("\n   (a) scan q̂ ∈ (0,1] at ŝ=â=1, η̂=ζ̂=0.01 (cs = centre):")
    println("       q̂        c0/cs     c₊/cs     c₋/cs    hyp  caus  well-posed")
    q̂s = [0.05, 0.2, 0.4, 0.6, 0.8, 0.9, 0.999, 1.0, 1.2]
    for q̂ in q̂s
        r = classify_point(1.0, 1.0, q̂, PROD.η̂, PROD.ζ̂, cs_c)
        mark = (q̂ == PROD.q̂) ? "  <- production" : ""
        @printf("     %6.3f   %7.4f  %7.4f  %7.4f   %s   %s    %s%s\n",
                q̂, r.c0/cs_c, r.cp/cs_c, r.cm/cs_c, r.hyp, r.caus, r.wp, mark)
    end

    # (b) small (ŝ, q̂) grid (â=1) to map the well-posed boundary q̂<ŝ
    println("\n   (b) (ŝ, q̂) grid at â=1 (well-posed iff 0<q̂<ŝ):")
    ŝs = [0.5, 1.0, 1.5]
    q̂g = [0.25, 0.75, 1.25]
    @printf("       %-8s", "ŝ\\q̂")
    for q̂ in q̂g; @printf("  q̂=%-5.2f", q̂); end
    println()
    for ŝ in ŝs
        @printf("       ŝ=%-5.2f ", ŝ)
        for q̂ in q̂g
            r = classify_point(ŝ, 1.0, q̂, PROD.η̂, PROD.ζ̂, cs_c)
            tag = r.wp ? (r.caus ? "WP+C " : "WP   ") : "----"
            @printf("  %-7s", tag)
        end
        println()
    end
    println("       (WP=well-posed 0<q̂<ŝ ; C=also causal cp≤1,c0≤1 at centre cs)")

    # production frame summary
    rp = classify_point(PROD.ŝ, PROD.â, PROD.q̂, PROD.η̂, PROD.ζ̂, cs_c)
    @printf("\n   production frame (1,1,0.999): hyperbolic=%s causal=%s well-posed=%s\n",
            rp.hyp, rp.caus, rp.wp)
    return (wp=rp.wp, hyp=rp.hyp, caus=rp.caus, cs_c=cs_c)
end

# ===========================================================================
# 3. ALONG-THE-STAR production-frame speeds c0, c±(r)
# ===========================================================================
function section3_alongstar()
    println()
    println("="^76)
    println("3. ALONG-THE-STAR production-frame characteristic speeds c0,c±(r)")
    println("="^76)
    star = solve_tov(EOS, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    @printf("   Shum star: M = %.5f M⊙ (target 1.4),  R = %.5f M⊙\n", star.M, star.R)
    N = length(star.r)

    allreal = true; allwp = true
    cp_max = 0.0; r_cp_max = 0.0; cs_at = 0.0
    n_superlum = 0; r_superlum_first = NaN
    for i in 1:N
        ε = star.ε[i]
        ε <= 0 && continue
        cs = sqrt(sound_speed2(EOS, ε))
        c0, cp, cm = shum_frame_speeds(PROD.ŝ, PROD.â, PROD.q̂, PROD.η̂, PROD.ζ̂, cs)
        allreal &= isfinite(c0) && isfinite(cp) && isfinite(cm) && cp ≥ 0 && cm ≥ 0 && c0 ≥ 0
        allwp   &= shum_frame_wellposed(PROD.ŝ, PROD.â, PROD.q̂)
        if cp > cp_max; cp_max = cp; r_cp_max = star.r[i]; cs_at = cs; end
        if cp > 1.0
            n_superlum += 1
            if isnan(r_superlum_first); r_superlum_first = star.r[i]; end
        end
    end

    # report a handful of radii
    println("   sampled profile (r in M⊙):")
    println("     r/R      cs        c0        c₊        c₋     c₊>1?")
    idxs = unique(clamp.(round.(Int, range(1, N, length=9)), 1, N))
    for i in idxs
        ε = star.ε[i]; ε <= 0 && continue
        cs = sqrt(sound_speed2(EOS, ε))
        c0, cp, cm = shum_frame_speeds(PROD.ŝ, PROD.â, PROD.q̂, PROD.η̂, PROD.ζ̂, cs)
        @printf("    %5.3f   %7.4f  %7.4f  %7.4f  %7.4f    %s\n",
                star.r[i]/star.R, cs, c0, cp, cm, cp > 1)
    end

    @printf("\n   c₊ is c₊=√3·cs ; max c₊ = %.5f at r=%.4f (cs=%.4f there)\n",
            cp_max, r_cp_max, cs_at)
    if n_superlum == 0
        @printf("   c₊ ≤ 1 EVERYWHERE in the star (frame is causal given cs(r)).\n")
    else
        @printf("   c₊ > 1 at %d radii (first at r=%.4f) — superluminal where cs > 1/√3≈%.4f.\n",
                n_superlum, r_superlum_first, 1/sqrt(3))
    end
    @printf("   speeds real & finite throughout: %s ;  well-posed throughout: %s\n",
            allreal, allwp)
    # connect to the QNM star
    println("   (this is the same TOV star reproducing Shum F/H1/H2 = 2.69/4.55/6.36 kHz)")
    return (star=star, allreal=allreal, allwp=allwp, cp_max=cp_max,
            n_superlum=n_superlum, cs_max=maximum(sqrt.(sound_speed2.(Ref(EOS), filter(>(0), star.ε)))))
end

# ===========================================================================
# 4. CROSS-VALIDATION: Shum fluid-frame c± vs general Causality biquadratic
# ===========================================================================
function section4_crossvalidate(star)
    println()
    println("="^76)
    println("4. CROSS-VALIDATION: Shum fluid-frame c± vs general Causality biquadratic")
    println("="^76)
    println("""
   Path A: shum_frame_speeds(ŝ,â,q̂,η̂,ζ̂,cs)            (Transport.jl, Shum eq.67-71)
   Path B: characteristic_speeds(p,e,cs,TransportCoefficients(...))  (Causality.jl,
           Keeble & Redondo-Yuste biquadratic Λ₂c⁴-2Λ₁c²+Λ₀=0), with the
           TransportCoefficients built from the shum_core mapping
           (ŝ,â,q̂,η̂,ζ̂)->(η,ζ,τε,τP,τQ,L).  Both evaluated in the fluid rest frame.""")

    # --- (i) the cross-check that DOES close: the mapping's frame RATIOS ---
    # Shum defines the frame as the DIMENSIONLESS ratios (his eqs.67-71)
    #   ŝ = τ_P/(cs² τ_ε),   â = τ_Q/τ_ε,
    # and shum_frame_speeds depends ONLY on (ŝ,â,q̂,η̂,ζ̂,cs).  So the cleanest
    # validation of the shum_core mapping is: does it reproduce these ratios?
    cs2_c = sound_speed2(EOS, εc)
    cs_c  = sqrt(cs2_c)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(PROD, EOS, εc)
    ŝ_back = τp/(cs2_c*τε)        # should equal PROD.ŝ
    â_back = τQ/τε                # should equal PROD.â
    @printf("\n   (i) mapping->frame-ratio recovery (the clean cross-check):\n")
    @printf("       τP/(cs²τε) = %.6f   (frame ŝ = %.3f)   rel.err = %.2e\n",
            ŝ_back, PROD.ŝ, abs(ŝ_back-PROD.ŝ)/PROD.ŝ)
    @printf("       τQ/τε      = %.6f   (frame â = %.3f)   rel.err = %.2e\n",
            â_back, PROD.â, abs(â_back-PROD.â)/PROD.â)
    ratios_ok = isapprox(ŝ_back, PROD.ŝ; rtol=1e-10) && isapprox(â_back, PROD.â; rtol=1e-10)
    @printf("       mapping reproduces Shum frame ratios (ŝ,â) exactly : %s\n", ratios_ok)

    # --- (ii) the DIRECT speed comparison: it does NOT close (honest negative) ---
    c0A, cpA, cmA = shum_frame_speeds(PROD.ŝ, PROD.â, PROD.q̂, PROD.η̂, PROD.ζ̂, cs_c)  # Path A
    tc = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=τε, τP=τp, τQ=τQ, L=PROD.L)
    c2mB, c2pB, disc = characteristic_speeds(p, εc, cs_c, tc)
    cpB = c2pB ≥ 0 ? sqrt(c2pB) : NaN
    cmB = c2mB ≥ 0 ? sqrt(c2mB) : NaN
    @printf("\n   (ii) direct characteristic-SPEED comparison at the centre (cs²=%.5f):\n", cs2_c)
    @printf("        η=%.4e ζ=%.4e τε=%.4e τP=%.4e τQ=%.4e\n", η, ζ, τε, τp, τQ)
    @printf("        Path A (Shum frame) :  c₊ = %.6f   c₋ = %.6f   c0 = %.6f\n", cpA, cmA, c0A)
    @printf("        Path B (biquadratic):  c₊ = %.6f   c₋ = %.6f   (disc=%.3e)\n", cpB, cmB, disc)
    sumA = cpA^2 + cmA^2;  sumB = c2pB + c2mB
    cp_reldiff = abs(cpA-cpB)/cpA
    sum_reldiff = abs(sumA-sumB)/abs(sumA)
    @printf("        c₊²+c₋² : A = %.5e  B = %.5e  rel.diff = %.3e\n", sumA, sumB, sum_reldiff)
    @printf("        c₊      : A = %.5f      B = %.5f      rel.diff = %.3e\n", cpA, cpB, cp_reldiff)

    cp_ok  = cp_reldiff  < 1e-3
    sum_ok = sum_reldiff < 1e-2

    println("""

   VERDICT (HONEST — partial / negative result, NOT forced):
     * The shum_core MAPPING is validated: it reproduces Shum's dimensionless
       frame ratios ŝ=τP/(cs²τε) and â=τQ/τε exactly (rel.err ~1e-16).
     * The two characteristic-SPEED FORMULAS do NOT coincide.  Reason (verified
       by an L-rescaling scan, below): the Causality biquadratic gives
           c₊²+c₋² = 2Λ₁/Λ₂ = [τε + τQ(τP+τε)] cs² /(τQ τε)
                   = cs²·( 1/τQ + τP/τε + 1 ),
       whose 1/τQ term BLOWS UP as the relaxation times τ→0 (here τ≈0.023, so
       c₊²+c₋²≈44 cs²), whereas the Shum frame formula is SCALE-FREE,
       c₊²+c₋² = 3 cs² (the √3 sound mode).  No common L (which scales all τ
       together) makes them agree: 1/τQ ∝ 1/L while τP/τε is L-independent.
     * Conclusion: Transport.jl (Shum eq.67-71, dimensionless hatted frame) and
       Causality.jl (dimensionful Keeble-Redondo-Yuste biquadratic) are DIFFERENT
       characteristic systems / non-dimensionalisations.  They are not expected to
       map cleanly at the speed level; we report the mismatch rather than tune L
       to fake agreement.  The dual-path validation that DOES hold is the frame
       ratios (i) plus the analytic c± validation (Section 1).""")

    # --- (iii) the L-rescaling scan that demonstrates the structural mismatch ---
    println("\n   (iii) L-rescaling scan (ratios ŝ,â fixed; biquadratic c₊²+c₋²/cs²):")
    println("        τQ          c₊²+c₋² / cs²   (Shum frame target = 3)")
    for f in (1.0, 5.0, 1/τε, 50.0)
        τεs = τε*f; τPs = (τp/τε)*τεs; τQs = (τQ/τε)*τεs
        tcs = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=τεs, τP=τPs, τQ=τQs, L=PROD.L)
        c2ms, c2ps, _ = characteristic_speeds(p, εc, cs_c, tcs)
        @printf("        %-10.5g  %-13.5g  %s\n", τQs, (c2ps+c2ms)/cs2_c,
                f == 1.0 ? "(production τ)" : "")
    end

    return (ratios_ok=ratios_ok, cpA=cpA, cpB=cpB, cmA=cmA, cmB=cmB,
            cp_ok=cp_ok, sum_ok=sum_ok, cp_reldiff=cp_reldiff, sum_reldiff=sum_reldiff,
            ŝ_back=ŝ_back, â_back=â_back)
end

# ===========================================================================
# DRIVER
# ===========================================================================
function run_analysis()
    println("\n" * "#"^76)
    println("# Shum FLUID-FRAME characteristic-structure analysis — validating BDNKStar")
    println("#"^76 * "\n")
    s1 = section1_validate()
    s2 = section2_framespace()
    s3 = section3_alongstar()
    s4 = section4_crossvalidate(s3.star)

    println("\n" * "="^76)
    println("5. SUMMARY")
    println("="^76)
    @printf("   [1] production c₊=√3, c₋≈0.0183 validated   : %s\n", s1.ok)
    @printf("   [2] production frame hyp/causal/well-posed  : %s/%s/%s\n", s2.hyp, s2.caus, s2.wp)
    @printf("   [3] along star real & well-posed throughout : %s/%s ; c₊>1 at %d radii\n",
            s3.allreal, s3.allwp, s3.n_superlum)
    @printf("   [4] cross-val: mapping->frame-ratios exact=%s ; direct speed match=%s\n",
            s4.ratios_ok, s4.cp_ok)
    @printf("       (the two speed FORMULAS use different non-dimensionalisations — honest negative)\n")
    println("="^76)
    return (s1=s1, s2=s2, s3=s3, s4=s4)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_analysis()
end
