#=
    radial_shoot_crosscheck.jl

    INDEPENDENT cross-check of the full-GR fundamental radial mode of the
    BDNKStar benchmark star, by a SHOOTING method — a DIFFERENT route than the
    Sturm–Liouville matrix eigensolver `chandrasekhar_radial_omega2` and the
    time-domain dynamical-GR engine DynGR1D (F_dyn ≈ 2.03 kHz).

    FORMULATION (relativistic radial pulsation, first-order (ξ, Δp) pair)
    --------------------------------------------------------------------
        ξ  ≡ Δr/r   (relative Lagrangian displacement)
        Δp           (Lagrangian pressure perturbation)
    Gondek, Haensel & Zdunik (1997) A&A 325, 217 (verbatim coefficients V,W,X,Y as
    reproduced in Pereira, Coelho & Rueda 2018, arXiv:1706.09371). SAME metric
    convention as TOV.jl / the engine:  ds² = -e^{ν} dt² + e^{λ} dr² + r² dΩ²
    with ν the FULL potential (g_tt=-e^ν), Φ ≡ ν/2, e^{λ}=1/(1-2m/r):

      dξ/dr  = V ξ + W Δp,   dΔp/dr = X ξ + Y Δp,   with
        V = -3/r - (dp/dr)/(p+ε)
        W = -(1/r)·1/(Γ₁ p)
        X = ω² e^{λ-ν}(p+ε) r - 4(dp/dr) + (dp/dr)² r/(p+ε) - 8π e^{λ}(p+ε) p r
        Y = (dp/dr)/(p+ε) - 4π (p+ε) r e^{λ}                       [the Δp term]

    Γ₁ = (p+ε)/p · c_s²  is the adiabatic index; ν stored FULL by TOV.jl, so the
    inertial term ω² e^{λ-ν} uses the SAME ν the dynamical engine evolves.
    [EARLIER-BUG NOTE: a previous attempt used Y = -ν' (i.e. -2(m+4πr³p)/(r(r-2m)));
     the CORRECT Y is (dp/dr)/(p+ε) - 4π(p+ε)r e^λ = -ν'/2 - 4π(p+ε)r e^λ. The
     wrong Y pushed the fundamental up to ~3.6 kHz; the correct Y lands on ~2 kHz.]

    BACKGROUND: SAME TOV star / EOS / εc the dynamical engine used
    (ShumPolytrope(κ=100), εc = 0.0015).

    BOUNDARY CONDITIONS
      centre r→0:  ξ'(0) finite ⇒  Δp = -3 Γ₁ p ξ   (regularity).
      surface R :  Δp(R) = 0  — the free-surface eigenvalue condition.
    Shoot from the centre (ξ=1, Δp from regularity), RK4 to R, root-find ω² on
    Δp(R)=0; the FUNDAMENTAL is the lowest ω² (eigenfunction ξ has NO interior
    nodes), confirmed by node-counting.
=#

using BDNKStar
using BDNKStar: solve_tov, ShumPolytrope, pressure, sound_speed2,
                energy_from_pressure,
                chandrasekhar_radial_omega2, cowling_radial_omega2,
                radial_cowling_spectrum
using BDNKStar.Units: kHz_to_km
using Printf

const CONV_kHz = 1.0 / kHz_to_km
fkHz(ω2) = sqrt(max(ω2, 0.0)) / (2π) * CONV_kHz

# ---- benchmark star (identical to test_dyngr.jl) --------------------------
const EOS  = ShumPolytrope(100.0)
const EPSC = 0.0015
const HTOV = 2e-5

star = solve_tov(EOS, EPSC; h=HTOV)
const R = star.R
const M = star.M
@printf("# benchmark star: ShumPolytrope(kappa=100), eps_c=%.5g km^-2\n", EPSC)
@printf("# R=%.5f km  M=%.5f km  2M/R=%.5f (compactness)\n", R, M, 2M/R)

# ---------------------------------------------------------------------------
# CO-INTEGRATED shooting: the TOV background (m,p,ν) and the perturbation (ξ,Δp)
# are advanced on ONE RK4 grid (the SAME RK4 stepper TOV.jl uses). This makes the
# background C∞-consistent with the integrator — interpolating a pre-tabulated
# TOV onto the perturbation grid injects C0 kinks that the stiff Δp/(Γ₁p)
# coupling amplifies into hundreds of spurious nodes; co-integration removes that.
# State y = (m, p, ν, ξ, Δp). ν is FULL (g_tt=-e^ν); we work with e^{-ν} relative
# to the centre (ν(0)=0) since the inertial term ω²e^{λ-ν} only needs ν up to an
# additive constant that the eigenvalue scale absorbs — but to be exact we apply
# the Schwarzschild shift e^{ν(R)}=1-2M/R by rescaling ω² post-hoc (see below).
# ---------------------------------------------------------------------------
@inline function deriv(r, y, ω2)
    m, p, ν, ξ, Δp = y
    p ≤ 0 && return (0.0, 0.0, 0.0, 0.0, 0.0)
    ε   = energy_from_pressure(EOS, p)
    den = r * (r - 2m); fac = m + 4π * r^3 * p
    dm  = 4π * r^2 * ε
    dν  = 2 * fac / den                                   # ν' (full)
    pp  = -(ε + p) * fac / den                            # dp/dr (TOV)
    eλ  = 1.0 / max(1 - 2m / r, 1e-12)
    eνm = exp(-ν)                                          # e^{-ν} (ν(0)=0 here)
    cs2 = clamp(sound_speed2(EOS, ε), 1e-12, 1.0)
    Γ1  = cs2 * (ε + p) / p
    # GHZ (1997) coefficients V,W,X,Y  (Pereira–Coelho–Rueda 2018 transcription):
    V = -3.0 / r - pp / (ε + p)
    W = -(1.0 / r) / (Γ1 * p)
    X = ω2 * eλ * eνm * (ε + p) * r - 4.0 * pp + pp^2 * r / (ε + p) -
        8π * eλ * (ε + p) * p * r
    Y = pp / (ε + p) - 4π * (ε + p) * r * eλ
    dξ  = V * ξ + W * Δp
    dΔp = X * ξ + Y * Δp
    return (dm, pp, dν, dξ, dΔp)
end

# RK4 from a Taylor-seeded regular centre to the surface (p→0). Returns
# (Δp(R), ξ(R), interior node count of ξ, ν(R), R_int). ν is referenced to ν(0)=0.
function shoot(ω2; h::Float64=5e-4, ptol_rel::Float64=1e-8)
    pc = pressure(EOS, EPSC); ptol = ptol_rel * pc
    r0 = h
    m = (4π/3) * EPSC * r0^3
    p = pc - 2π*(EPSC + pc)*(EPSC/3 + pc) * r0^2
    ν = 0.0
    ε0 = energy_from_pressure(EOS, p)
    Γ10 = clamp(sound_speed2(EOS, ε0), 1e-12, 1.0) * (ε0 + p) / p
    ξ = 1.0; Δp = -3.0 * Γ10 * p * ξ                       # centre regularity
    y = (m, p, ν, ξ, Δp); r = r0; nodes = 0; prevξ = ξ
    Δp_s = Δp; ξ_s = ξ; ν_s = ν; R_int = r0
    while p > ptol && r < 100.0
        k1 = deriv(r,     y, ω2)
        y2 = ntuple(i -> y[i] + h/2*k1[i], 5); k2 = deriv(r+h/2, y2, ω2)
        y3 = ntuple(i -> y[i] + h/2*k2[i], 5); k3 = deriv(r+h/2, y3, ω2)
        y4 = ntuple(i -> y[i] + h*k3[i],   5); k4 = deriv(r+h,   y4, ω2)
        yn = ntuple(i -> y[i] + h/6*(k1[i] + 2k2[i] + 2k3[i] + k4[i]), 5)
        r += h; pn = yn[2]
        if pn ≤ ptol
            frac = y[2] / (y[2] - pn)
            ξ_s  = y[4] + frac*(yn[4] - y[4])
            Δp_s = y[5] + frac*(yn[5] - y[5])
            ν_s  = y[3] + frac*(yn[3] - y[3])
            R_int = (r - h) + frac*h
            break
        end
        # count GENUINE interior nodes only (exclude the surface boundary layer
        # p < 1e-4 p_c where ξ diverges as p→0 — a coordinate artifact, not a node)
        if pn > 1e-4*pc && sign(yn[4]) != sign(prevξ)
            nodes += 1
        end
        prevξ = yn[4]
        p = pn; y = yn
    end
    return Δp_s, ξ_s, nodes, ν_s, R_int
end

# The inertial term ω² e^{λ-ν} uses ν referenced to ν(0)=0. The PHYSICAL ν is
# shifted so e^{ν(R)} = 1-2M/R (Schwarzschild match, exactly TOV.jl's convention).
# Shifting ν → ν + Δν multiplies e^{-ν} by e^{-Δν}, i.e. it rescales ω²_raw by
# e^{-Δν}. So the physical eigenvalue is ω²_phys = ω²_raw · e^{Δν}, with
# Δν = ν_phys(R) - ν_raw(R) = ln(1-2M/R) - ν_raw(R). We solve in ω²_raw then map.
const ΔΝ_SHIFT = let (_, _, _, νR, Rint) = shoot(1e-3)
    log(1 - 2M/Rint) - νR
end
ω2_phys(ω2_raw) = ω2_raw * exp(ΔΝ_SHIFT)

surfΔp(ω2_raw) = shoot(ω2_raw)[1]

# ---- scan for the LOWEST sign change of Δp(R) (fundamental: ξ nodeless), bisect
function find_fundamental(; ω2lo=1e-4, ω2hi=1.5e-2, nscan=400)
    grid = collect(range(ω2lo, ω2hi; length=nscan))
    prevv = surfΔp(grid[1]); prevω = grid[1]; bracket = nothing; ndb = 0
    for ω2 in grid[2:end]
        cur = surfΔp(ω2)
        if sign(cur) != sign(prevv) && bracket === nothing
            bracket = (prevω, ω2); ndb = shoot(0.5*(prevω+ω2))[3]
        end
        prevv = cur; prevω = ω2
    end
    bracket === nothing && error("no Δp(R) sign change in [$ω2lo,$ω2hi]")
    a, c = bracket; fa = surfΔp(a)
    for _ in 1:200
        mid = 0.5*(a+c); fm = surfΔp(mid)
        sign(fm) == sign(fa) ? (a = mid; fa = fm) : (c = mid)
        (c - a) < 1e-13 && break
    end
    return 0.5*(a+c), ndb
end

println("\n=== SHOOTING (full-GR radial pulsation, co-integrated GHZ (xi,Dp)) ===")
ω2_raw, nd = find_fundamental()
ω2_shoot = ω2_phys(ω2_raw)
@printf("fundamental: omega^2 = %.6e km^-2   F = %.4f kHz   (interior nodes of xi = %d)\n",
        ω2_shoot, fkHz(ω2_shoot), nd)
@printf("   (raw omega^2=%.6e, Schwarzschild-nu shift factor e^Dnu=%.5f)\n",
        ω2_raw, exp(ΔΝ_SHIFT))

println("\n=== SL matrix eigensolver (chandrasekhar_radial_omega2) ===")
ω2_sl = chandrasekhar_radial_omega2(EOS, EPSC; nmodes=2, N=3000, h_tov=HTOV)
@printf("omega^2_SL[1] = %.6e   F_SL = %.4f kHz   (overtone F=%.4f kHz)\n",
        ω2_sl[1], fkHz(ω2_sl[1]), fkHz(ω2_sl[2]))

println("\n=== relativistic COWLING (frozen metric) ===")
fcow, ω2cow, _ = radial_cowling_spectrum(EOS, EPSC; N=2000, h_tov=5e-5, nmodes=2)
@printf("omega^2_Cowling[1] = %.6e   F_Cowling = %.4f kHz\n", ω2cow[1], fcow[1])

println("\n=== SUMMARY: three independent routes ===")
const FDYN = 2.03
@printf("  F_dyn   (DynGR1D time-domain) = %.4f kHz  [engine benchmark]\n", FDYN)
@printf("  F_SL    (SL matrix eig)       = %.4f kHz\n", fkHz(ω2_sl[1]))
@printf("  F_shoot (this script)         = %.4f kHz\n", fkHz(ω2_shoot))
@printf("  F_Cowling (frozen metric)     = %.4f kHz\n", fcow[1])
@printf("\n  shoot vs SL     : %+.2f%%\n", 100*(fkHz(ω2_shoot)-fkHz(ω2_sl[1]))/fkHz(ω2_sl[1]))
@printf("  shoot vs engine : %+.2f%%\n", 100*(fkHz(ω2_shoot)-FDYN)/FDYN)
@printf("  omega^2_GR=%.6e  vs  omega^2_Cowling=%.6e\n", ω2_shoot, ω2cow[1])
@printf("  SIGN: omega^2_GR %s omega^2_Cowling  -> %s\n",
        ω2_shoot < ω2cow[1] ? "<" : ">",
        ω2_shoot < ω2cow[1] ?
            "GR BELOW Cowling: spacetime response softens restoring force (CORRECT)" :
            "GR ABOVE Cowling: the earlier sign bug")
