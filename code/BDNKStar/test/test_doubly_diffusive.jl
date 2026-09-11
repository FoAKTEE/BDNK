using Test
using BDNKStar

# LOCAL (WKB) DOUBLY-DIFFUSIVE dispersion analysis (DoublyDiffusive.jl): two KNOWN
# instability families re-derived in the BDNK causal framework on the toolkit's own
# stratified ideal-gas star (N² from brunt_vaisala) and transport frame (η, κ_Q, τ
# from conformal_frame_PMP):
#   (1) GSF / rotational doubly-diffusive (Goldreich–Schubert 1967; Fricke 1968;
#       local Acheson 1978 / Menou-Balbus-Spruit 2004),
#   (2) thermohaline / semiconvective overstability + heat-coupled g-mode
#       (Stern 1960; Baines–Gill 1969).
# These are NOT new instabilities; the new content is their BDNK realization plus
# the causal telegrapher modification D k² → D k²/(1+sτ).
#
# Asserted here (all LOCAL/algebraic, fast):
#  (a) IDEAL-limit validation — GSF reduces to Solberg–Høiland/Rayleigh; thermohaline
#      reduces to the Brunt–Väisälä g-mode (ν=χ=0 ⇒ s=±i√N²) and to monotonic
#      Schwarzschild/Ledoux convection when N²<0.
#  (b) the doubly-diffusive instability appears (Re s>0) under the right conditions
#      (SH-stable-but-Rayleigh-unstable GSF cone widening; Pr<1 thermohaline overstable).
#  (c) growth scales with the diffusivity and VANISHES as it →0 in the doubly-diffusive
#      band (the hallmark: the mode exists ONLY because of differential diffusion).
#  (d) the BDNK relaxation-time signature: thermohaline high-k diffusive rate SATURATES
#      at ~1/τ (UV cutoff) while NS grows ∝Dk²; GSF causal rates reduce to NS as τ→0
#      and suppress the growth rate while leaving the marginal boundary τ-invariant.

@testset "Doubly-diffusive (GSF / thermohaline) local BDNK dispersion" begin

    # ─ stratified, SH-stable background (Γ>Γ_struct ⇒ N²>0) ─────────────────── #
    star = solve_tov_idealgas(; Γ=2.0, Γ_struct=1.9, K=100.0, ρc=1.28e-3)
    r, N2, ce2, cs2 = brunt_vaisala(star)
    i0  = argmin(abs.(r .- 0.5*star.R))
    N2_0 = N2[i0]
    @test N2_0 > 0                                   # convectively/SH STABLE

    # transport frame (η, κ_Q, τ) at this point
    tc = conformal_frame_PMP(star.ε[i0])
    @test tc.η > 0 && tc.κQ > 0

    # WKB diffusivity scaling (as in repro/gsf_doubly_diffusive.jl): pin νk²≈0.3|N|
    k_test = 50.0 / star.R
    ν_kin  = 0.3 * sqrt(abs(N2_0)) / k_test^2
    Pr_inv = 30.0                                    # κ_th/ν (radiative/GSF regime)
    κth    = Pr_inv * ν_kin
    τη = (tc.τε/tc.η) * ν_kin
    τQ = (tc.τQ/tc.κQ) * κth

    # rotation: power-law Ω(ϖ)=Ω0(ϖ/R)^q, Ω0 ≈ 5% Ω_dyn, evaluated at ϖ/R=0.5.
    # The LOCAL Ω (not Ω0) sets the epicyclic frequency κ_epi²=4Ω²(1+q/2).
    Ω0    = 0.05 * sqrt(star.ε[i0])
    ϖfrac = 0.5
    Ωloc(q) = Ω0 * ϖfrac^q
    κepi(q) = gsf_epicyclic(Ωloc(q), q)

    # ====================================================================== #
    # (a) IDEAL-LIMIT VALIDATION
    # ====================================================================== #
    @testset "(a) ideal limit ⇒ Solberg–Høiland/Rayleigh and Brunt g-mode" begin
        # GSF ideal (ν=κth=0): s²=−(κ²cos²ψ+N²sin²ψ). Stable ∀ψ iff κ²≥0 ∧ N²≥0.
        κ2_stable   = κepi(-0.3)                # mild shear ⇒ κ²>0
        κ2_rayleigh = κepi(-2.5)                # steep ⇒ κ²<0 (Rayleigh-unstable)
        @test κ2_stable   > 0
        @test κ2_rayleigh < 0
        s_SH  = maximum(gsf_max_growth_NS(κ2_stable,   N2_0, 0.0, 0.0, k_test, ψ)
                        for ψ in range(0, π/2; length=181))
        s_Ray = maximum(gsf_max_growth_NS(κ2_rayleigh, N2_0, 0.0, 0.0, k_test, ψ)
                        for ψ in range(0, π/2; length=181))
        @test s_SH  ≤ 1e-12                           # SH criterion satisfied ⇒ stable
        @test s_Ray > 1e-6                            # Rayleigh violated ⇒ unstable
        # the ideal Rayleigh growth peaks at the axial wavevector (ψ=0)
        @test gsf_max_growth_NS(κ2_rayleigh, N2_0, 0.0, 0.0, k_test, 0.0) > 1e-6

        # Thermohaline ideal (ν=χ=0): s=±i√N² (Brunt–Väisälä g-mode), exactly.
        N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0); N2net = N2T + N2mu
        rs0 = dd_cubic_roots(k_test; ν=0.0, χ=0.0, N2mu=N2mu, N2T=N2T)
        s0  = rs0[argmax(abs.(imag.(rs0)))]
        @test isapprox(abs(imag(s0)), sqrt(N2net); rtol=1e-6)
        @test abs(real(s0)) < 1e-9                    # purely oscillatory (marginal)

        # Thermohaline N²<0 ⇒ direct MONOTONIC convection (Schwarzschild/Ledoux)
        rsneg = dd_cubic_roots(k_test; ν=ν_kin, χ=κth, N2mu=-abs(N2mu), N2T=-abs(N2T))
        re, sneg = most_unstable(rsneg)
        @test re > 0                                  # unstable
        @test abs(imag(sneg)) < 1e-9                  # monotonic (Im s=0)
    end

    # ====================================================================== #
    # (b) THE DOUBLY-DIFFUSIVE INSTABILITY APPEARS (Re s>0)
    # ====================================================================== #
    @testset "(b) GSF cone widening + thermohaline overstability appear" begin
        # GSF: κ²<0 (Rayleigh-unstable on cylinders) but N²>0. Thermal diffusion
        # erodes the stabilizing N² ⇒ the unstable wavevector cone WIDENS.
        κ2 = κepi(-3.0)
        @test κ2 < 0
        cone_ideal, peak_ideal, _ = gsf_unstable_cone(κ2, N2_0, 0.0,   0.0, k_test)
        cone_NS,    peak_NS,    _ = gsf_unstable_cone(κ2, N2_0, ν_kin, κth, k_test)
        @test peak_NS    > 1e-6                        # GSF UNSTABLE
        @test cone_NS    > cone_ideal                  # diffusion WIDENS the cone
        @test cone_ideal > 0                           # adiabatic cone is non-empty

        # genuinely doubly-diffusive band: a ψ that is adiabatically STABLE but
        # diffusively UNSTABLE.
        ψ_dd = deg2rad(cone_ideal + 15)
        @test gsf_max_growth_NS(κ2, N2_0, 0.0,   0.0, k_test, ψ_dd) ≤ 1e-12  # ideal: stable
        @test gsf_max_growth_NS(κ2, N2_0, ν_kin, κth, k_test, ψ_dd) > 1e-6   # diffusive: UNSTABLE

        # Thermohaline: net N²>0 (convectively stable) but Pr=ν/χ<1 ⇒ OVERSTABLE.
        N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0)
        χ = 1.0e-3; ν_th = 0.16*χ                       # Pr=0.16<1 (frame value)
        k_res  = sqrt(sqrt(abs(N2T+N2mu)) / χ)
        k_scan = 10 .^ range(log10(k_res)-3, log10(k_res)+3; length=600)
        g, kstar = dd_max_growth(k_scan; ν=ν_th, χ=χ, N2mu=N2mu, N2T=N2T)
        @test g > 1e-9                                  # overstable growth
        rstar = dd_cubic_roots(kstar; ν=ν_th, χ=χ, N2mu=N2mu, N2T=N2T)
        @test dd_is_overstable(rstar)                   # Re s>0 AND Im s≠0 (oscillatory)

        # Pr>1 (momentum diffuses faster than heat) ⇒ NO overstability.
        χ2 = ν_th / 2.0                                 # Pr=ν/χ=2>1
        g_hiPr, _ = dd_max_growth(k_scan; ν=ν_th, χ=χ2, N2mu=N2mu, N2T=N2T)
        @test g_hiPr ≤ 1e-9
    end

    # ====================================================================== #
    # (c) GROWTH SCALES WITH DIFFUSIVITY AND VANISHES AS IT →0
    # ====================================================================== #
    @testset "(c) doubly-diffusive growth vanishes as diffusivity→0" begin
        # GSF: in the diffusively-opened band the growth is driven ENTIRELY by the
        # diffusion (the band is adiabatically stable). Scaling both diffusivities
        # by f→0 (at fixed Pr) the growth → 0 — the doubly-diffusive hallmark. (The
        # dependence is non-monotone at intermediate f — growth peaks at moderate
        # diffusivity — but it strictly vanishes as f→0 and is small at small f.)
        κ2 = κepi(-3.0)
        cone_ideal, _, _ = gsf_unstable_cone(κ2, N2_0, 0.0, 0.0, k_test)
        ψ_dd = deg2rad(cone_ideal + 15)
        gf(f) = gsf_max_growth_NS(κ2, N2_0, f*ν_kin, f*κth, k_test, ψ_dd)
        @test gf(0.25) > 1e-6                            # unstable at finite diffusivity
        @test gf(0.0)  ≤ 1e-12                           # VANISHES at zero diffusivity
        @test gf(0.01) < gf(0.25)                        # small diffusivity ⇒ small growth
        @test gf(0.01) > gf(0.0)                         # and approaches 0 from above

        # GSF Prandtl threshold: the opened band needs κth/ν above a threshold
        # (heat must outrun momentum diffusion to win against N²).
        s_lowPr  = gsf_max_growth_NS(κ2, N2_0, ν_kin,   1.0*ν_kin, k_test, ψ_dd)  # κth/ν=1
        s_highPr = gsf_max_growth_NS(κ2, N2_0, ν_kin, 100.0*ν_kin, k_test, ψ_dd)  # κth/ν=100
        @test s_highPr > s_lowPr                        # stronger heat diffusion ⇒ more unstable

        # Thermohaline: at a FIXED wavenumber (the resonance) the overstable growth
        # → 0 as the diffusivities → 0 at fixed Pr (the doubly-diffusive hallmark:
        # the drive is the differential diffusion itself). (NB the k-MAXIMISED growth
        # is ~N·f(Pr) and only the resonant k∝1/√χ shifts; so we probe at fixed k.)
        N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0)
        χ0 = 1.0e-3
        k_res = sqrt(sqrt(abs(N2T+N2mu)) / χ0)
        gx(f) = (rs = dd_cubic_roots(k_res; ν=0.16*f*χ0, χ=f*χ0, N2mu=N2mu, N2T=N2T);
                 first(most_unstable(rs)))
        @test gx(0.25) > 1e-6                            # overstable at finite χ
        @test gx(0.10) < gx(0.25)                       # growth scales DOWN with χ
        @test gx(0.01) < gx(0.10)                       # monotone toward 0
        @test gx(0.0)  ≤ 1e-12                           # VANISHES as χ→0

        # Thermohaline: NO destabilising thermal source (N²_T=0, all buoyancy
        # stable) ⇒ overstability VANISHES for any Pr (necessary condition).
        k_scan = 10 .^ range(log10(k_res)-3, log10(k_res)+3; length=600)
        gns, _ = dd_max_growth(k_scan; ν=0.16*χ0, χ=χ0, N2mu=abs(N2_0), N2T=0.0)
        @test gns ≤ 1e-9
    end

    # ====================================================================== #
    # (d) BDNK RELAXATION-TIME SIGNATURE
    # ====================================================================== #
    @testset "(d) BDNK causal signature: high-k cutoff + τ→0 reduction" begin
        N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0)
        ν_th = 0.2*tc.κQ/(star.ε[i0]+star.p[i0])
        χ_th = tc.κQ/(star.ε[i0]+star.p[i0])
        τπ = tc.τε; τQ_th = tc.τQ
        k_res = sqrt(sqrt(abs(N2T+N2mu)) / χ_th)

        # (i) THERMOHALINE UV CUTOFF: the NS diffusive rate ~Dk² grows without
        #     bound at high k; BDNK SATURATES near the relaxation pole 1/τ.
        fastest(rs) = maximum(abs.(real.(rs)))
        klo = k_res
        khi = k_res * 1.0e4
        ns_lo = fastest(dd_cubic_roots(klo; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T))
        ns_hi = fastest(dd_cubic_roots(khi; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T))
        bd_lo = fastest(dd_bdnk_roots(klo; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T, τπ=τπ, τQ=τQ_th))
        bd_hi = fastest(dd_bdnk_roots(khi; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T, τπ=τπ, τQ=τQ_th))
        @test ns_hi > 1e3 * ns_lo                       # NS rate explodes ∝k²
        @test bd_hi < 5.0 * bd_lo                       # BDNK rate SATURATES
        # the BDNK saturation sits near the relaxation poles 1/τ
        @test bd_hi < 10.0 * max(1/τπ, 1/τQ_th)

        # (ii) τ→0 reduction: BDNK → NS for the physical (low-k) overstable branch.
        k_scan = 10 .^ range(log10(k_res)-3, log10(k_res)+3; length=600)
        gNS, _  = dd_max_growth(k_scan; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T)
        gBDe, _ = dd_max_growth(k_scan; ν=ν_th, χ=χ_th, N2mu=N2mu, N2T=N2T,
                                bdnk=true, τπ=1e-12*τπ, τQ=1e-12*τQ_th)
        @test isapprox(gNS, gBDe; rtol=1e-4)            # τ→0 ⇒ BDNK ≡ NS

        # (iii) GSF causal reduction + growth-rate suppression with τ-invariant onset.
        κ2 = κepi(-3.0)
        cone_ideal, _, _ = gsf_unstable_cone(κ2, N2_0, 0.0, 0.0, k_test)
        ψ_dd = deg2rad(cone_ideal + 15)
        sNS  = gsf_max_growth_NS(κ2, N2_0, ν_kin, κth, k_test, ψ_dd)
        sBD0 = gsf_max_growth_BDNK(κ2, N2_0, ν_kin, κth, k_test, ψ_dd, 1e-12*τη, 1e-12*τQ)
        @test isapprox(sNS, sBD0; rtol=1e-5)            # τ→0 ⇒ BDNK ≡ NS (GSF)
        # finite (amplified) τ SUPPRESSES the secular GSF growth rate (does not raise it)
        sBDτ = gsf_max_growth_BDNK(κ2, N2_0, ν_kin, κth, k_test, ψ_dd, 5000*τη, 5000*τQ)
        @test sBDτ < sNS                                # causal slows already-unstable modes
        @test sBDτ > 0                                  # but it is still unstable (onset unmoved)
    end

    # ====================================================================== #
    # heat-coupled g-mode: Re s ∝ κ_Q, driven (sign set by N²_T<0, Pr<1)
    # ====================================================================== #
    @testset "(e) heat-coupled g-mode driven ∝ κ_Q (thermohaline continuation)" begin
        N2T = -0.6*abs(N2_0); N2mu = +1.0*abs(N2_0)
        χ_th = tc.κQ/(star.ε[i0]+star.p[i0])
        ν_th = (tc.η/tc.κQ) * χ_th                       # keep frame Pr
        k_res = sqrt(sqrt(abs(N2T+N2mu)) / χ_th)
        k_gm  = 0.2 * k_res                              # weak-diffusion (χk²≲N)
        # the g-root stays oscillatory and acquires Re s>0 ∝ κ_Q (driven)
        s_gm = gmode_root(k_gm, 0.1*χ_th; ν=ν_th, N2mu=N2mu, N2T=N2T)
        @test abs(imag(s_gm)) > 1e-9                     # still oscillatory (g-mode)
        slope = gmode_growth_slope(k_gm; ν=ν_th, N2mu=N2mu, N2T=N2T,
                                   χa=0.02*χ_th, χb=0.08*χ_th)
        @test slope > 0                                  # DRIVEN/overstable g-mode ∝ κ_Q
    end
end
