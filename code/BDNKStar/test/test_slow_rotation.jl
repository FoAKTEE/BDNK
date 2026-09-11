using Test
using BDNKStar
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# SLOW-ROTATION (Hartle 1967) frame dragging + moment of inertia on the static
# TOV background. Validation gate against PUBLISHED universal relations:
#
#  PRIMARY  Breu & Rezzolla (2016), MNRAS 459, 646 [arXiv:1601.06083], Eq. (20)
#           + Table 2 SLOW-ROTATION row:
#             Ī = I/M³ = 0.8134 C⁻¹ + 0.2101 C⁻² + 3.175e-3 C⁻³ − 2.717e-4 C⁻⁴
#           (χ²red = 0.1184; ⟨L₁⟩ ≈ 1.1%, L∞ ≈ 3%).  C = M/R.
#  CROSS    Lattimer & Schutz (2005), ApJ 629, 979:
#             Ĩ = I/(MR²) = 0.237 + 0.674 C + 4.48 C⁴.
#  ANCHOR   PSR J0737-3039A I_{1.338} = 1.425×10⁴⁵ g cm² (Reed et al. 2022,
#           arXiv:2204.09000) — concrete observational/ab-initio I scale.
#
# The four Read et al. (2009) realistic EOS (SLy/APR4/H4/MS1) all land within
# ~1.6% of Breu-Rezzolla at 1.4 M⊙. The ShumPolytrope(100) Γ=2 TOY EOS sits at
# the soft edge of the universal band (~10% below Breu) and is NOT gated — it is
# outside the realistic-EOS set the relation was fit to (documented caveat).
# ─────────────────────────────────────────────────────────────────────────────

const _GC = BDNKStar.Units.gram_per_cm3_to_km_minus2   # g/cm³ → km⁻²

# central-density window (km⁻²) bracketing a 1.4 M⊙ star for each EOS
const _ε14WIN = Dict(
    :SLy  => (6e14*_GC, 1.6e15*_GC),
    :APR4 => (6e14*_GC, 1.6e15*_GC),
    :H4   => (4e14*_GC, 1.2e15*_GC),
    :MS1  => (2.5e14*_GC, 1.0e15*_GC),
)

"Find the model in [εlo,εhi] whose mass is closest to `Mtarget` (coarse scan)."
function _near_mass(eos, εlo, εhi, Mtarget; n=30, h=0.01)
    best = nothing; bd = Inf
    for εc in exp.(range(log(εlo), log(εhi); length=n))
        r = moment_of_inertia(eos, εc; h=h)
        d = abs(r.M_Msun - Mtarget)
        if d < bd; bd = d; best = r; end
    end
    return best
end

@testset "Slow rotation (Hartle 1967): frame dragging + moment of inertia" begin

    # ── 1. Basic sanity: I>0, drag ∈ (0,1), exterior matching ────────────────
    @testset "positivity, drag fraction, profile" begin
        eos = piecewise_polytrope(:SLy)
        r = moment_of_inertia(eos, 1e15*_GC; h=0.005)
        @test r.I_geo_km3 > 0
        @test r.I_cgs > 0
        @test r.Ibar > 0
        @test 0 < r.drag_ratio < 1                    # ω(0)/Ω in (0,1)
        @test r.Ω > 0
        @test r.J > 0
        # ϖ(r) = Ω - ω(r) is positive and monotonically INCREASING outward
        # (ω = drag decreases outward), so ω(r)/Ω is largest at the center:
        @test all(r.ϖ .> 0)
        @test issorted(r.ϖ)                            # ϖ rises center→surface
        @test r.ω_over_Ω[1] ≈ r.drag_ratio rtol=1e-12
        @test r.ω_over_Ω[1] > r.ω_over_Ω[end]          # dragging rises toward center
        @test r.ω_over_Ω[end] > 0                       # still some drag at the surface
    end

    # ── 2. Normalization independence (linear-homogeneous ODE) ───────────────
    @testset "I independent of central normalization ϖc" begin
        eos = piecewise_polytrope(:APR4)
        r1 = moment_of_inertia(eos, 1.1e15*_GC; ϖc=1.0,  h=0.005)
        r2 = moment_of_inertia(eos, 1.1e15*_GC; ϖc=37.0, h=0.005)
        @test isapprox(r1.I_geo_km3, r2.I_geo_km3; rtol=1e-10)
        @test isapprox(r1.Ibar,      r2.Ibar;      rtol=1e-10)
        @test isapprox(r1.drag_ratio,r2.drag_ratio;rtol=1e-10)
    end

    # ── 3. cgs conversion factor derived from Units (≈1.3466e43 g cm²/km³) ────
    @testset "km³ → g cm² conversion (Units-derived)" begin
        @test isapprox(KM3_TO_G_CM2, 1.3466e43; rtol=1e-3)
        eos = piecewise_polytrope(:SLy)
        r = moment_of_inertia(eos, 1e15*_GC; h=0.005)
        @test isapprox(r.I_cgs, r.I_geo_km3 * KM3_TO_G_CM2; rtol=1e-14)
    end

    # ── 4. VALIDATION GATE: 1.4 M⊙ I vs Breu-Rezzolla & per-EOS scale ─────────
    @testset "1.4 M⊙ I vs published relations (Breu 2016 / Lattimer-Schutz)" begin
        ngate = 0
        for sym in (:SLy, :APR4, :H4, :MS1)
            eos = piecewise_polytrope(sym)
            εlo, εhi = _ε14WIN[sym]
            r = _near_mass(eos, εlo, εhi, 1.4)
            @test isapprox(r.M_Msun, 1.4; atol=0.05)         # actually near 1.4
            predB = IBAR_FROM_C_BREU(r.C)
            predLS = MR2_FROM_C_LS(r.C)
            pctB  = 100*(r.Ibar - predB)/predB
            pctLS = 100*(r.I_MR2 - predLS)/predLS
            @info "Slow-rotation I @1.4 M⊙" EOS=sym M=r.M_Msun R_km=r.R_km C=r.C I_cgs=r.I_cgs Ibar=r.Ibar Breu=predB pct_Breu=pctB I_MR2=r.I_MR2 LS=predLS pct_LS=pctLS drag=r.drag_ratio
            # gate: within the relation's quoted scatter (~5%, L∞≈3% for Breu)
            @test abs(pctB) ≤ 6.0
            @test abs(pctLS) ≤ 6.0
            # sanity: 1.4 M⊙ I ~ 1–2×10⁴⁵ g cm². The very stiff MS1 (R≈14.9 km)
            # genuinely sits slightly above 2×10⁴⁵ (I_1.4≈2.1×10⁴⁵) — stiffer EOS
            # ⇒ larger I; allow 2.2×10⁴⁵ for it (honest EOS-dependent band).
            Ihi = sym === :MS1 ? 2.2e45 : 2.0e45
            @test 1.0e45 ≤ r.I_cgs ≤ Ihi
            ngate += 1
        end
        @test ngate ≥ 3                                     # ≥3 EOS gated (have 4)
    end

    # ── 5. Ī monotonically DECREASING with compactness (mass sequence) ───────
    # Breu Ī ~ 1/C² ⇒ more compact (larger C) ⇒ smaller Ī. Check on an SLy seq.
    @testset "Ī monotone in compactness (SLy mass sequence)" begin
        eos = piecewise_polytrope(:SLy)
        Cs = Float64[]; Ibars = Float64[]
        for εc in exp.(range(log(6e14*_GC), log(2.0e15*_GC); length=10))
            r = moment_of_inertia(eos, εc; h=0.008)
            push!(Cs, r.C); push!(Ibars, r.Ibar)
        end
        p = sortperm(Cs)
        @test issorted(Ibars[p]; rev=true)                 # Ī falls as C rises
        # and each point tracks Breu within scatter over the realistic range
        for (C, Ib) in zip(Cs, Ibars)
            @test abs(100*(Ib - IBAR_FROM_C_BREU(C))/IBAR_FROM_C_BREU(C)) ≤ 8.0
        end
    end

    # ── 6. ShumPolytrope toy EOS: runs, positive, drag-rising (NOT gated) ─────
    @testset "ShumPolytrope(100) toy EOS — runs, soft-edge caveat" begin
        sh = ShumPolytrope(100.0)
        r = moment_of_inertia(sh, 2e-3; h=0.008)
        @test r.I_geo_km3 > 0
        @test 0 < r.drag_ratio < 1
        @test r.ω_over_Ω[1] > r.ω_over_Ω[end]
        # documented: Γ=2 toy sits ~5-15% below the realistic-EOS Breu band
        pct = 100*(r.Ibar - IBAR_FROM_C_BREU(r.C))/IBAR_FROM_C_BREU(r.C)
        @info "ShumPolytrope(100) (toy, ungated)" M=r.M_Msun C=r.C Ibar=r.Ibar Breu=IBAR_FROM_C_BREU(r.C) pct_Breu=pct
        @test pct < 0                                       # below the band, as expected
    end
end
