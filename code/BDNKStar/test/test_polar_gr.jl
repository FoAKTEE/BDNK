using Test
using BDNKStar
using Printf
const PG = BDNKStar.PolarGRModes
const AX = BDNKStar.AxialViscousModes
const MSUN_KM = BDNKStar.Units.Msun_to_km

# ============================================================================
# STAGE 3 — FULL-GR POLAR (even-parity) quasi-normal modes (Cowling DROPPED).
#
# Formulation: the Lindblom–Detweiler interior metric+fluid system in the
# numerically-stable {H1,K,W,X} variables (Detweiler & Lindblom 1985, ApJ 292,12;
# compact form Krüger–Gaertig–Kokkotas 2011, arXiv:1106.2736 Appendix A),
# matched to the exterior ZERILLI equation with the OUTGOING-wave (QNM) boundary
# condition.  The outgoing BC is imposed analytically via the Chandrasekhar–
# Detweiler (Darboux) transform of the VALIDATED Regge–Wheeler Leaver continued
# fraction (AxialViscousModes.up_logderiv).  Convention e^{+iωt}, ω=2πf+i/τ.
#
# Benchmark star: EOS1 polytrope (κ=100, n=1, ρc=3e15 g/cc) ⇒ M≈1.27 M⊙,
# R≈8.86 km — the SAME star the axial w-mode solver used.
#
# VALIDATED (full f-mode now recovered): (1) Zerilli potential → Schwarzschild form;
#  (2) the surface (H1,K)→Zerilli map; (3) the Darboux+Leaver outgoing solution vs
#  direct integration; (4) the full-GR f-MODE: f=2.86 kHz (BELOW the Cowling 3.30 —
#  the GR overestimate-correction, ratio 0.87) with a finite GW damping τ=0.115 s,
#  robust across seeds, agreeing with the Andersson–Kokkotas (1998) universal relations
#  to ≈8% (frequency) and ≈10% (damping); and a polar w-mode (f≈9.2 kHz, τ≈64 μs).
#  This required correcting the LD X′ equation against DL85 eq (11): the H0 term sign,
#  the X-term e^{+ν/2} in eq (5), and the ½ on the V coefficient −½ℓ(ℓ+1)ν′/r².
# ============================================================================
@testset "Polar full-GR QNM (Lindblom–Detweiler interior + Zerilli exterior)" begin
    eos, star = build_polar_star()
    M, R = star.M, star.R
    G = 6.6743015e-11; c = 299_792_458.0
    εc = (3e15 * 1e3) * G / c^2 * 1e6

    # ---- the benchmark star (same as the axial solver) ----
    @test isapprox(M/MSUN_KM, 1.27; atol=0.02)
    @test isapprox(R, 8.86;     atol=0.05)
    @test 0.20 < M/R < 0.22                              # compact EOS1 star

    # ---- unit conversion round-trip (e^{+iωt}: ω = 2πf + i/τ) ----
    ω = PG.polar_ftau_to_omega(4.0, 0.02)
    f_kHz, τ_s = PG.polar_omega_to_ftau(ω)
    @test isapprox(f_kHz, 4.0;  rtol=1e-12)
    @test isapprox(τ_s,   0.02; rtol=1e-12)
    @test imag(ω) > 0                                    # damped ⇒ Imω>0

    # ---- (1) Zerilli potential → standard Schwarzschild black-hole form ----
    # At large r, V_Z → ℓ(ℓ+1)/r² (centrifugal); positive and decaying.
    Mbh = 1.0
    VZ_far = zerilli_potential(Mbh, 1e6, 2)
    @test isapprox(VZ_far, 2*3/(1e6)^2; rtol=1e-3)       # ℓ(ℓ+1)=6
    @test zerilli_potential(Mbh, 10.0, 2) > 0           # positive barrier
    # n = (ℓ-1)(ℓ+2)/2 = 2 for ℓ=2 enters the denominator (nr+3M)²:
    @test zerilli_potential(Mbh, 3.0, 2) > zerilli_potential(Mbh, 8.0, 2)

    # ---- (2) surface (H1,K) → Zerilli map is ODE-self-consistent ----
    ωt = PG.polar_ftau_to_omega(3.0, 0.1)
    H1, K, rs = PG.interior_surface(star, eos, ωt, 2; rmin=1e-3, nint=6000,
                                    surf_frac=1-1e-4)
    Z0, dZ0 = PG.surface_to_zerilli(M, rs, 2, ωt, H1, K)
    # step the Zerilli ODE a hair and finite-difference dZ/dr* vs the analytic seed
    h = 1e-3
    Zp, _ = PG.integrate_zerilli(M, ωt, 2, Z0, dZ0; r0=rs, r1=rs+h, nsteps=20)
    fR = 1 - 2M/rs
    dZ_num = (Zp - Z0) / (h/fR)                          # numerical dZ/dr*
    @test isapprox(dZ_num, dZ0; rtol=1e-3)              # analytic seed consistent
    @test abs(Z0) > 0                                    # nontrivial surface data

    # ---- (3) analytic Darboux+Leaver outgoing solution == direct integration ----
    # (validates the Chandrasekhar transform + the e^{+iωt}↔e^{−iωt} convention map)
    a = 20*R
    Zu, dZu = PG.outgoing_zerilli_leaver(M, ωt, 2, a; Ncf=900)
    Lz_leaver = dZu / Zu
    Zs, dZs = PG.outgoing_seed(M, ωt, 2, a + 3000.0)
    Zi, dZi = PG.integrate_zerilli(M, ωt, 2, Zs, dZs; r0=a+3000.0, r1=a, nsteps=30000)
    Lz_integ = dZi / Zi
    @test isapprox(real(Lz_leaver), real(Lz_integ); atol=1e-3)
    @test isapprox(imag(Lz_leaver), imag(Lz_integ); rtol=2e-2)   # outgoing phase −iω
    @test imag(Lz_leaver) < 0                            # outgoing e^{−iωr*} ⇒ Im<0

    # ---- (4) HEADLINE: the full-GR f-MODE — GR-lowered frequency + GW damping ----
    # (recovered after correcting eq-11 H0 sign, eq-5 X-term e^{+ν/2}, eq-11 V-term ½.)
    res = polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(2.8, 0.15), star=star)
    @info "polar full-GR f-mode" f_kHz=res.f_kHz tau_s=res.tau_s omega=res.omega residual=res.residual
    @test res.converged
    @test res.residual < 1e-6
    @test imag(res.omega) > 0                            # DAMPED (Imω≠0): coupled to GR
    @test 0.05 < res.tau_s < 0.5                         # finite GW damping time (f-mode ballpark)
    @test 2.4 < res.f_kHz < 3.2                          # f-mode of this M=1.27 M⊙ star

    # robustness: every nearby seed converges to the SAME f-mode root
    res2 = polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(2.4, 0.10), star=star)
    @test res2.converged
    @test isapprox(res2.f_kHz, res.f_kHz; rtol=2e-3)
    @test isapprox(res2.tau_s, res.tau_s; rtol=2e-2)

    # ---- the GR signature: f_GR < f_Cowling (Cowling overestimates the f-mode) ----
    fc, _, _ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=1, N=6000,
                                          nscan=800, Lunit_km=1.0, ω2lo=1e-4, ω2hi=0.2)
    f_cowling = fc[1]
    @test res.f_kHz < f_cowling                          # GR LOWERS the f-mode (τ_Cowling=∞)
    @test 0.80 < res.f_kHz/f_cowling < 0.95             # ~13% below Cowling (typical)

    # ---- validation vs the Andersson–Kokkotas (1998) universal relations ----
    # (R is in km here; M_⊙ via Msun_to_km).  AK has ~10–20% intrinsic scatter.
    Mb = (M/MSUN_KM)/1.4; Rb = R/10.0
    f_AK   = 0.78 + 1.635*sqrt(Mb/Rb^3)
    tau_AK = 1/((Mb^3/Rb^4)*(22.85 - 14.65*(Mb/Rb)))
    @info "f-mode vs AK1998" f_GR=res.f_kHz f_AK=f_AK tau_GR=res.tau_s tau_AK=tau_AK
    @test abs(res.f_kHz - f_AK)/f_AK   < 0.20           # frequency within AK scatter
    @test abs(res.tau_s - tau_AK)/tau_AK < 0.30         # GW damping time within AK scatter

    # ---- a polar w-mode also exists (high f, short τ) ----
    resw = polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(9.0, 5e-5), star=star)
    @test resw.converged && resw.f_kHz > 6.0 && resw.tau_s < 1e-3
end
