using Test
using BDNKStar
using Printf

# =============================================================================
# ℓ=2 TIDAL DEFORMABILITY — static Love number k₂/Λ, the I-Love-Q completion,
# the GW170817 confrontation, and the DYNAMICAL (frequency-dependent) Λ_eff(ω).
#
# STATIC method (source-verified): Postnikov–Prakash–Lattimer 2010 (PRD 82,
# 024016; arXiv:1004.5098) y(r)=rH'/H Riccati ODE on the TOV background, k₂ via
# Hinderer 2008 (ApJ 677,1216) + arXiv:0711.2420 erratum, Λ=(2/3)k₂C^{-5}.
# I-Love: Yagi–Yunes 2013 (Science 341,365; arXiv:1302.4499) Ī–Λ̄ fit, with Ī
# from SlowRotation.moment_of_inertia.  GW170817: Abbott et al. 2018 (PRL 121,
# 161101) Λ_1.4 ≲ 800 (90% upper).
#
# DYNAMICAL method: Hinderer et al. 2016 (PRL 116,181101) / Steinhoff et al.
# 2016 (PRD 94,104028) effective-Love-number form Λ_eff(ω)=Σ_n Λ_n ω_n²/(ω_n²−ω²),
# Σ Λ_n = Λ_static; f-mode dominated.
#
# Kept fast: coarse central-density bisection on solve_tov, h=1.5e-3 km y-ODE.
# =============================================================================

# central energy density (km⁻²) tuned to a target gravitational mass (M⊙)
function _εc_for_mass(eos, Mt; lo=3e-4, hi=4e-3)
    f(εc) = mass_solar(solve_tov(eos, εc; h=2.5e-3)) - Mt
    a, b = lo, hi; fa, fb = f(a), f(b); it = 0
    while fa*fb > 0 && it < 40; b *= 1.15; fb = f(b); it += 1; end
    for _ in 1:70
        c = 0.5*(a+b); fc = f(c)
        fa*fc ≤ 0 ? (b = c; fb = fc) : (a = c; fa = fc)
        abs(b-a) < 1e-7 && break
    end
    return 0.5*(a+b)
end

# Published Λ(1.4) references (CITED): Read 2009 / Hinderer 2010 / Lackey&Wade
# 2015 standard PP-EOS values. NOTE the MS1 preset here is the genuine (stiffest)
# MS1 (R_1.4≈14.9 km), NOT MS1b (R_1.4=14.5, Λ≈1220); its Λ≈1387 is correct.
const _PUB_L14 = Dict(:SLy=>300.0, :APR4=>245.0, :H4=>900.0)

@testset "ℓ=2 tidal deformability (k₂, Λ) + I-Love-Q + dynamical Λ_eff(ω)" begin

    # cache the 1.4 M⊙ stars / Λ once for the soft EOS (reused across subtests)
    L14 = Dict{Symbol,NamedTuple}()
    for sym in (:SLy, :APR4)
        eos = piecewise_polytrope(sym; N=4000)
        εc  = _εc_for_mass(eos, 1.4)
        star = solve_tov(eos, εc; h=1.5e-3)
        k2, Λ, yR = tidal_love_number(eos, star)
        C = star.M/star.R
        L14[sym] = (eos=eos, εc=εc, star=star, k2=k2, Λ=Λ, yR=yR, C=C,
                    M=mass_solar(star), R=star.R)
    end

    # =====================================================================
    # (a) Λ(1.4) for >=2 EOS within tol of published (SLy, APR4)
    # =====================================================================
    @testset "a. Λ(1.4 M⊙) vs published (SLy, APR4)" begin
        for sym in (:SLy, :APR4)
            r = L14[sym]; pub = _PUB_L14[sym]
            dev = abs(r.Λ - pub)/pub
            @info "Λ(1.4)" EOS=sym M=r.M R=r.R C=r.C k2=r.k2 Λ=r.Λ Λ_pub=pub dev_pct=100*dev
            @test isapprox(r.M, 1.4; atol=0.02)
            @test 0.05 < r.k2 < 0.15            # standard k₂(1.4) band
            @test dev < 0.05                    # within 5% of published Λ(1.4)
        end
        # ordering: APR4 (softer/smaller R) is less deformable than SLy
        @test L14[:APR4].Λ < L14[:SLy].Λ
    end

    # =====================================================================
    # (b) I-Love-Q: our (Ī, Λ) lie on the Yagi-Yunes universal fit, EOS-blind
    # =====================================================================
    @testset "b. I-Love universal relation (Yagi-Yunes 2013)" begin
        maxdev = 0.0; npts = 0
        for sym in (:SLy, :APR4, :H4)
            eos = piecewise_polytrope(sym; N=4000)
            for Mt in (1.2, 1.4, 1.6)
                εc = _εc_for_mass(eos, Mt)
                star = solve_tov(eos, εc; h=1.5e-3)
                _, Λ, _ = tidal_love_number(eos, star)
                Ibar = moment_of_inertia(star).Ibar
                Iyy  = ibar_yagi_yunes(Λ)
                dev  = abs(Ibar - Iyy)/Iyy
                maxdev = max(maxdev, dev); npts += 1
                @test dev < 0.03                # on the universal curve to <3%
            end
        end
        @info "I-Love EOS-independence" npts=npts max_dev_pct=100*maxdev
        @test npts ≥ 9
        @test maxdev < 0.03                     # EOS-independent (universal)
    end

    # =====================================================================
    # (c) GW170817: soft SLy/APR4 satisfy Λ_1.4 < 800; stiff MS1 disfavored
    # =====================================================================
    @testset "c. GW170817 confrontation (Λ_1.4 ≲ 800)" begin
        for sym in (:SLy, :APR4)
            Λ = L14[sym].Λ
            @info "GW170817 consistent" EOS=sym Λ=Λ
            @test Λ < 800                       # consistent with the 90% upper bound
        end
        # stiffest EOS (MS1): Λ_1.4 ≫ 800 ⇒ DISFAVORED by GW170817 (LVC 1805.11581)
        eos = piecewise_polytrope(:MS1; N=4000)
        εc  = _εc_for_mass(eos, 1.4)
        star = solve_tov(eos, εc; h=1.5e-3)
        _, Λms1, _ = tidal_love_number(eos, star)
        @info "GW170817 disfavored" EOS=:MS1 Λ=Λms1 R=star.R
        @test Λms1 > 800                        # exceeds the bound ⇒ disfavored
        @test star.R > 14.0                     # MS1 is the stiff/large-R EOS
    end

    # =====================================================================
    # (d) DYNAMICAL: Λ_eff(ω→0)=Λ_static and Λ_eff amplifies toward resonance
    # =====================================================================
    @testset "d. dynamical Λ_eff(ω): static limit + f-mode amplification" begin
        Λs = L14[:SLy].Λ
        f_f = 2.0                                # representative ℓ=2 f-mode [kHz]
        g_freqs = [0.66, 0.46, 0.36]            # weak g-mode tower
        m = lambda_eff_modes(Λs, f_f; g_freqs_kHz=g_freqs, g_frac=1e-3)

        # static sum rule Σ Λ_n = Λ_static (by construction)
        @test isapprox(sum(m.mode_Λ), Λs; rtol=1e-12)

        # (d.i) static limit Λ_eff(ω→0) = Λ_static
        Λ0 = lambda_eff(m, 1e-7)
        @info "static limit" Λ_eff_0=Λ0 Λ_static=Λs Δrel=abs(Λ0-Λs)/Λs
        @test isapprox(Λ0, Λs; rtol=1e-6)

        # (d.ii) amplification: Λ_eff grows monotonically toward the f-mode
        # resonance f_GW = ω_f/π = 2 f_f, staying > Λ_static and diverging.
        f_res = 2*f_f                            # f_GW resonance frequency [kHz]
        fGWs = [0.5, 1.0, 1.5, 2.0]            # below resonance (f_res = 4 kHz)
        vals = [lambda_eff(m; f_GW_kHz=fg) for fg in fGWs]
        @info "f-mode amplification" f_res_kHz=f_res Λeff_over_Λs=vals./Λs
        @test all(v -> v > Λs, vals)            # enhanced above the static value
        @test issorted(vals)                    # grows toward resonance
        @test vals[end]/Λs > 1.2                # sizeable enhancement near 2 kHz
        @test lambda_eff(m; f_GW_kHz=0.95*f_res)/Λs > 5  # strong near resonance

        # pole sign: just above resonance Λ_eff flips negative (driven oscillator)
        @test lambda_eff(m; f_GW_kHz=1.01*f_res) < 0
    end
end
