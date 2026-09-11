using Test
using BDNKStar
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# r-MODE (l=m=2) CFS GW INSTABILITY WINDOW + BDNK causal ζ_eff extension.
#
# VALIDATION GATE (Lindblom, Owen & Morsink 1998, PRL 80, 4843 [gr-qc/9803053];
# Owen, Lindblom, Cutler, Schutz, Vecchio & Andersson 1998, PRD 58, 084020
# [gr-qc/9804044], Table I; Andersson & Kokkotas 2001, IJMPD 10, 381):
#   • frequencies σ=(2/3)Ω, ω=−(4/3)Ω exactly
#   • n=1 polytrope J̃=1.635e-2, Ĩ=0.261, τ̃_GW=−3.26 s, τ̃_sv=2.52e8 s,
#     τ̃_bv=6.99e8 s at Ω=√(πGρ̄), T=10⁹ K (faithful Lane–Emden benchmark)
#   • window minimum near T~10⁹–10⁹·⁵ K with Ω_c/Ω_K~0.06; ν_c tens–hundreds Hz
#   • causal ζ_eff(στ)=ζ_NS/(1+(στ)²) widens the high-T window vs Navier–Stokes
#
# NOTE on the shear prefactor: LOM98 Eq.18 is (l−1)(2l+1) = 5 for l=2.  The
# attached derivations wrote "3·5=15" (slip: (l−1)=1, not 3); the gate below
# REQUIRES 5 (with 15, τ_sv is 3× off).  The empirical gate caught the error.
# ─────────────────────────────────────────────────────────────────────────────

const _GC = BDNKStar.Units.gram_per_cm3_to_km_minus2

"Build an n=1 (Γ=2) relativistic TOV star near M=1.4 M⊙ (K=190 ⇒ R≈12.6 km)."
function _n1_star(; K=190.0)
    eos = isentropic_idealgas(Γ=2.0, K=K, ρ_lo=1e-10, ρ_hi=1e-2, N=2500)
    best = nothing; bd = Inf
    for ρc in exp.(range(log(3e-4), log(2e-3); length=80))
        pc = K*ρc^2; s = solve_tov(eos, ρc + pc; h=0.0035)
        d = abs(mass_solar(s) - 1.4); if d < bd; bd = d; best = s; end
    end
    return best, eos
end

"Realistic-EOS star near 1.4 M⊙."
function _eos_star(sym; n=40)
    eos = piecewise_polytrope(sym)
    best = nothing; bd = Inf
    for εc in exp.(range(log(6e14*_GC), log(2.2e15*_GC); length=n))
        s = solve_tov(eos, εc; h=0.006)
        d = abs(mass_solar(s) - 1.4); if d < bd; bd = d; best = s; end
    end
    return best, eos
end

@testset "r-modes (l=m=2): CFS instability window + BDNK causal ζ_eff" begin

    # ── 1. Frequencies: σ=(2/3)Ω, ω=−(4/3)Ω exactly ──────────────────────────
    @testset "frequencies" begin
        for Ω in (1.0, 3.7, 1234.0)
            σ, ω = rmode_frequencies(Ω)
            @test σ ≈ (2/3)*Ω rtol=1e-14
            @test ω ≈ -(4/3)*Ω rtol=1e-14
            @test ω ≈ σ - 2Ω rtol=1e-14           # ω = σ − mΩ, m=2
        end
        # general l=m: σ=2mΩ/[l(l+1)], ω=−(l−1)(l+2)/(l+1) Ω
        σ3, ω3 = rmode_frequencies(2.0; l=3, m=3)
        @test σ3 ≈ 2*3*2.0/(3*4) rtol=1e-14
        @test ω3 ≈ -(3-1)*(3+2)/(3+1)*2.0 rtol=1e-14
        # GW dimensionless coefficient 32π/225·(4/3)^6 = 2.5104
        @test RMODE_GW_COEFF ≈ 2.5104 atol=1e-3
    end

    # ── 2. Sign structure: τ_GW<0 (driving), τ_sv,τ_bv>0 (damping) ────────────
    @testset "sign structure (GW driving, viscous damping)" begin
        star, eos = _n1_star()
        ΩK = kepler_frequency(star, eos)
        τGW, τsv, τbv, τtot = rmode_timescales(star, eos; Ω=0.5*ΩK, T=1e9)
        @test τGW < 0                              # CFS driving
        @test τsv > 0                              # shear damping
        @test τbv > 0                              # bulk damping
        @test isfinite(τtot)
        # near Kepler GW driving is ~seconds (sanity)
        τGWk, = rmode_timescales(star, eos; Ω=ΩK, T=1e9)
        @test 0.1 < abs(τGWk) < 100.0
    end

    # ── 3. LOM98 n=1 Lane–Emden benchmark (faithful, Newtonian) ──────────────
    @testset "LOM98 n=1 benchmark (J̃, Ĩ, τ_GW, τ_sv, τ_bv)" begin
        v = rmode_validate_lom98()
        pct(c, p) = 100*(c - p)/abs(p)
        @info "LOM98 r-mode validation (n=1 Lane–Emden, Ω=√(πGρ̄), T=1e9 K)" J̃=v.J̃ Ĩ=v.Ĩ τ_GW=v.τ_GW τ_sv=v.τ_sv τ_bv=v.τ_bv pct_J̃=pct(v.J̃,JTILDE_PUB) pct_Ĩ=pct(v.Ĩ,ITILDE_PUB) pct_GW=pct(v.τ_GW,TAU_GW_PUB) pct_sv=pct(v.τ_sv,TAU_SV_PUB) pct_bv=pct(v.τ_bv,TAU_BV_PUB)
        # GENUINE validations (fully computed, no fitted constants):
        @test abs(pct(v.J̃,  JTILDE_PUB)) ≤ 2.0    # structure integral
        @test abs(pct(v.Ĩ,  ITILDE_PUB)) ≤ 2.0
        @test abs(pct(v.τ_GW, TAU_GW_PUB)) ≤ 5.0   # GW (fully computed)
        @test abs(pct(v.τ_sv, TAU_SV_PUB)) ≤ 5.0   # shear (fully computed)
        # SELF-CONSISTENCY only (NOT an independent validation): _BULK_CAL is fit
        # to this very number, so it cannot catch a ζ-microphysics-coefficient
        # error — it only confirms the calibration + per-star integral are wired.
        @test abs(pct(v.τ_bv, TAU_BV_PUB)) ≤ 5.0   # bulk (calibrated to benchmark)
        @test v.τ_GW < 0
    end

    # ── 4. Instability window: minimum near T~10⁹ K, unstable wedge ──────────
    @testset "instability window (minimum, unstable region)" begin
        star, eos = _n1_star()
        Tg = 10 .^ range(8.5, 10.4; length=40)
        w = rmode_instability_window(star, eos; Tgrid=Tg)
        @test w.Ω_K > 0
        @test 200 < w.ν_K_Hz < 2000               # Kepler ~ hundreds of Hz–1 kHz
        @test all(w.Ω_crit .> 0)
        imin = argmin(w.Ω_crit_over_ΩK)
        Tmin = w.T[imin]
        @info "Window minimum" T_min=Tmin Ω_c_over_Ω_K=w.Ω_crit_over_ΩK[imin] ν_c_Hz=w.ν_crit_Hz[imin]
        @test 1e9 ≤ Tmin ≤ 5e9                     # minimum near T~10⁹–few×10⁹ K
        @test 0.03 ≤ w.Ω_crit_over_ΩK[imin] ≤ 0.12 # LOM98 ~0.06
        @test 10 ≤ w.ν_crit_Hz[imin] ≤ 300         # tens–hundreds of Hz
        # a hot, fast star (above the minimum) is UNSTABLE
        _,_,_,τtot = rmode_timescales(star, eos; Ω=0.7*w.Ω_K, T=Tmin)
        @test τtot < 0                             # unstable (1/τ<0)
        # and the minimum is interior (rises on both sides — a wedge)
        @test w.Ω_crit_over_ΩK[imin] < w.Ω_crit_over_ΩK[1]
        @test w.Ω_crit_over_ΩK[imin] < w.Ω_crit_over_ΩK[end]
    end

    # ── 5. BDNK causal ζ_eff: ζ_eff(0)=ζ_NS, →0 as ωτ→∞, widens high-T window ─
    @testset "BDNK causal ζ_eff(ω)=ζ_NS/(1+(ωτ)²)" begin
        @test zeta_eff_factor(1000.0, 0.0) == 1.0          # τ=0 ⇒ Navier–Stokes
        @test zeta_eff_factor(1000.0, 1e-2) < 1.0          # suppressed
        @test zeta_eff_factor(1000.0, 1e3) < 1e-6          # ωτ→∞ ⇒ →0
        @test zeta_eff_factor(0.0, 1e-3) == 1.0            # ω=0 ⇒ NS limit
        # monotonic decrease in τ
        f1 = zeta_eff_factor(3000.0, 1e-4); f2 = zeta_eff_factor(3000.0, 1e-3)
        @test f2 < f1 < 1.0

        # causal τ SUPPRESSES bulk damping ⇒ widens the unstable region at high T
        star, eos = _n1_star()
        Tg = 10 .^ range(9.8, 10.6; length=20)
        wNS = rmode_instability_window(star, eos; Tgrid=Tg, τ_bulk_relax=0.0)
        wC  = rmode_instability_window(star, eos; Tgrid=Tg, τ_bulk_relax=1e-3)
        # at high T the causal critical spin is at-or-below the NS one (window wider)
        @test all(wC.Ω_crit_over_ΩK .≤ wNS.Ω_crit_over_ΩK .+ 1e-9)
        @test any(wC.Ω_crit_over_ΩK .< wNS.Ω_crit_over_ΩK .- 1e-3)   # strictly wider somewhere
        # equivalently: at a fixed high-T, high-spin point, causal is "more unstable"
        Thi = 2e10; Ωhi = 0.8*kepler_frequency(star, eos)
        _,_,τbvNS,_ = rmode_timescales(star, eos; Ω=Ωhi, T=Thi, τ_bulk_relax=0.0)
        _,_,τbvC,_  = rmode_timescales(star, eos; Ω=Ωhi, T=Thi, τ_bulk_relax=1e-3)
        @test τbvC > τbvNS                          # weaker bulk damping (longer τ_bv)
    end

    # ── 6. EOS dependence: realistic SLy/APR4 windows (≈EOS-independent min) ──
    @testset "realistic EOS windows (SLy, APR4)" begin
        Tg = 10 .^ range(8.5, 10.4; length=35)
        mins = Float64[]
        for sym in (:SLy, :APR4)
            star, eos = _eos_star(sym)
            @test isapprox(mass_solar(star), 1.4; atol=0.06)
            w = rmode_instability_window(star, eos; Tgrid=Tg)
            imin = argmin(w.Ω_crit_over_ΩK)
            @info "EOS window" EOS=sym M=mass_solar(star) R=star.R ν_K_Hz=w.ν_K_Hz T_min=w.T[imin] Ω_c_over_Ω_K=w.Ω_crit_over_ΩK[imin] ν_c_Hz=w.ν_crit_Hz[imin]
            @test 1e9 ≤ w.T[imin] ≤ 5e9
            @test 0.03 ≤ w.Ω_crit_over_ΩK[imin] ≤ 0.12
            @test 300 < w.ν_K_Hz < 2000
            push!(mins, w.Ω_crit_over_ΩK[imin])
        end
        # near-EOS-independence of the dimensionless minimum (LOM98 ~0.06 for any EOS)
        @test maximum(mins) - minimum(mins) < 0.03
    end
end
