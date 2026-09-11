using Test
using BDNKStar

# STAGE 3 — GRAVITY (g) modes: relativistic Cowling polar (ℓ≥2) eigensolver that
# INCLUDES the buoyancy / Schwarzschild-discriminant term, so it yields f, p AND g
# modes. The g-modes are a genuine NEW family that appears ONLY when the star is
# stratified (cs²≠ce², i.e. adiabatic index Γ ≠ structure index Γ_struct).
#
# Equations: the (U,V) relativistic-Cowling system of Jaikumar et al. 2021 /
# Shirke et al. (arXiv:2506.18892 Eq.12-13), with the Schwarzschild-discriminant
# sign of Gaertig & Kokkotas (arXiv:0905.0821 Eq.12). N² = g²(1/ce²−1/cs²)e^{ν−λ},
# g = −(dp/dr)/(ε+p) the local gravity, cs² the ADIABATIC (frozen-composition,
# index Γ) sound speed, ce²=dp/dε the EQUILIBRIUM (structure, index Γ_struct) one.
# Validated barotropic limit: when Γ=Γ_struct ⇒ N²≡0 ⇒ the solver reduces to the
# f/p tower and matches NonRadialModes (the validated barotropic anchor) to ~1%.
#
# Physics (the easy-to-flip sign, asserted numerically below):
#   Γ > Γ_struct  ⇒  cs²>ce²  ⇒  N²>0  ⇒  STABLE real g-modes BELOW the f-mode.
#   Γ = Γ_struct  ⇒  N²≡0     ⇒  NO g-modes (only f/p).

const _K = 30.0
const _ρc = 3e-3

@testset "Gravity (g) modes: relativistic Cowling, buoyancy from Γ vs Γ_struct" begin

    # ---------------------------------------------------------------------- #
    # (a) BAROTROPIC star (Γ_struct = Γ): N²≈0 everywhere, NO g-modes, and the
    #     f/p tower matches the validated NonRadialModes solver to a few %.
    # ---------------------------------------------------------------------- #
    @testset "(a) barotropic ⇒ N²≈0, no g-modes, f/p match NonRadialModes" begin
        Γ = 2.0
        sb = solve_tov_idealgas(Γ=Γ, Γ_struct=Γ, K=_K, ρc=_ρc, h=2e-4)
        r, N2, ce2, cs2 = brunt_vaisala(sb)
        @test maximum(abs.(N2)) < 1e-8           # N²≈0 (machine zero) everywhere
        @test all(isapprox.(ce2, cs2; rtol=1e-6))# ce²=cs² (no stratification)

        spb = gmode_spectrum(sb; l=2)
        @test spb.N2_max < 1e-8
        @test spb.has_gmodes == false            # NO buoyancy tower
        @test isempty(spb.g_freqs_kHz)
        @test spb.f_freq > 0 && length(spb.p_freqs) ≥ 2
        @test spb.f_freq < spb.p_freqs[1]        # f below p₁ (f/p ordering)
        @test issorted(spb.p_freqs)

        # f/p must agree with the validated barotropic eigensolver (same star, km units)
        ig = isentropic_idealgas(Γ=Γ, K=_K, ρ_lo=1e-9, ρ_hi=2*_ρc, N=1200)
        εc = _ρc + _K*_ρc^Γ/(Γ-1)
        fp,_,_ = nonradial_cowling_spectrum(ig, εc; l=2, nmodes=3, N=8000,
                                            Lunit_km=1.0, ω2lo=1e-4, ω2hi=0.3, nscan=1500)
        @info "(a) barotropic f/p [kHz]" gmode=round.([spb.f_freq; spb.p_freqs[1:2]],digits=3) nonradial=round.(fp,digits=3)
        @test isapprox(spb.f_freq,    fp[1]; rtol=0.03)   # f-mode within 3%
        @test isapprox(spb.p_freqs[1], fp[2]; rtol=0.03)  # p₁ within 3%
    end

    # ---------------------------------------------------------------------- #
    # (b)+(c) STRATIFIED, convectively STABLE star (Γ=2 > Γ_struct=5/3):
    #     N²>0 in the interior, a REAL g-mode tower (g₁>g₂>…) BELOW the f-mode,
    #     with f and the p-modes still present and ordered f < p₁ < p₂.
    # ---------------------------------------------------------------------- #
    @testset "(b)+(c) stratified Γ>Γ_struct ⇒ N²>0, real g-tower below f, f/p intact" begin
        ss = solve_tov_idealgas(Γ=2.0, Γ_struct=5/3, K=_K, ρc=_ρc, h=2e-4)
        r, N2, ce2, cs2 = brunt_vaisala(ss)
        @test maximum(N2) > 0                     # buoyancy present
        @test all(N2 .> -1e-10)                   # N²>0 in the interior (no convective instability)
        @test all(cs2 .> ce2)                     # adiabatic stiffer than equilibrium ⇒ stable

        sps = gmode_spectrum(ss; l=2, N=5000, nscan=2000)
        @info "(b) stratified spectrum [kHz]" g=round.(sps.g_freqs_kHz,digits=3) f=round(sps.f_freq,digits=3) p=round.(sps.p_freqs[1:min(3,end)],digits=3) N2max=round(sps.N2_max,sigdigits=3)
        @test sps.has_gmodes == true              # the NEW family appears
        @test length(sps.g_freqs_kHz) ≥ 2         # at least two g-modes resolved
        @test sps.g_freqs_kHz[1] > sps.g_freqs_kHz[2]          # g₁ > g₂ (tower)
        @test issorted(sps.g_freqs_kHz; rev=true)             # g₁>g₂>… descending
        @test all(sps.g_freqs_kHz .< sps.f_freq)              # ALL g-modes BELOW the f-mode
        @test all(sps.g_freqs_kHz .> 0)                       # real (no instability)
        @test length(sps.p_freqs) ≥ 1
        @test sps.f_freq < sps.p_freqs[1]                     # f below p₁ (ordering)
        @test issorted(sps.p_freqs)
    end

    # ---------------------------------------------------------------------- #
    # (d) RESOLUTION CONVERGENCE of at least g₁ and the f-mode.
    # ---------------------------------------------------------------------- #
    @testset "(d) resolution convergence of g₁ and the f-mode" begin
        ss = solve_tov_idealgas(Γ=2.0, Γ_struct=5/3, K=_K, ρc=_ρc, h=2e-4)
        sp_lo  = gmode_spectrum(ss; l=2, N=3000, nscan=2000)
        sp_hi  = gmode_spectrum(ss; l=2, N=7000, nscan=2000)
        @info "(d) convergence [kHz]" f_lo=round(sp_lo.f_freq,digits=4) f_hi=round(sp_hi.f_freq,digits=4) g1_lo=round(sp_lo.g_freqs_kHz[1],digits=4) g1_hi=round(sp_hi.g_freqs_kHz[1],digits=4)
        @test isapprox(sp_lo.f_freq, sp_hi.f_freq; rtol=5e-3)             # f-mode converged <0.5%
        @test isapprox(sp_lo.g_freqs_kHz[1], sp_hi.g_freqs_kHz[1]; rtol=5e-3)  # g₁ converged <0.5%
        @test isapprox(sp_lo.g_freqs_kHz[2], sp_hi.g_freqs_kHz[2]; rtol=1e-2)  # g₂ converged <1%
    end

    # ---------------------------------------------------------------------- #
    # (e) SANITY: the g-mode frequency rises with the stratification (hence N²).
    #     CLEAN knob: fix Γ_struct ⇒ IDENTICAL background star (the hydrostatic
    #     structure depends only on Γ_struct, K, ρc) and vary ONLY the adiabatic
    #     index Γ. Larger Γ ⇒ larger cs² ⇒ larger N² ⇒ higher g-modes, with the
    #     stellar M,R unchanged (so the rise is buoyancy, not a structure artefact).
    # ---------------------------------------------------------------------- #
    @testset "(e) g-mode frequency increases with stratification (N²)" begin
        weak   = solve_tov_idealgas(Γ=1.8, Γ_struct=5/3, K=_K, ρc=_ρc, h=2e-4)
        strong = solve_tov_idealgas(Γ=2.2, Γ_struct=5/3, K=_K, ρc=_ρc, h=2e-4)
        @test isapprox(weak.M, strong.M; rtol=1e-10)   # SAME background star (M,R)
        @test isapprox(weak.R, strong.R; rtol=1e-10)
        _, N2w, _, _ = brunt_vaisala(weak)
        _, N2s, _, _ = brunt_vaisala(strong)
        @test maximum(N2s) > maximum(N2w)              # stiffer adiabat ⇒ larger N²
        spw = gmode_spectrum(weak;   l=2)
        sps = gmode_spectrum(strong; l=2)
        @info "(e) stratification vs g₁" Γ_weak=1.8 Γ_strong=2.2 g1_weak=round(spw.g_freqs_kHz[1],digits=3) g1_strong=round(sps.g_freqs_kHz[1],digits=3) N2max_weak=round(maximum(N2w),sigdigits=3) N2max_strong=round(maximum(N2s),sigdigits=3)
        @test spw.has_gmodes && sps.has_gmodes
        @test sps.g_freqs_kHz[1] > spw.g_freqs_kHz[1]  # stronger buoyancy ⇒ higher g₁
    end
end
