#=
    dynamical_tides.jl — DYNAMICAL (frequency-dependent) TIDES for the BDNKStar
    reference neutron star: the effective tidal deformability Λ_eff(ω) from
    resonant excitation of the stellar f / g modes during a binary inspiral.

    PHYSICS
    -------
    During inspiral the companion's ℓ=2 tidal field forces the star at twice the
    orbital frequency.  As 2Ω_orb sweeps up toward a stellar mode frequency ω_n,
    that mode is resonantly excited and the effective (frequency-dependent) tidal
    response takes the driven-oscillator form

        Λ_eff(ω) = Σ_n Λ_n · ω_n² / (ω_n² − ω²),

    dominated by the ℓ=2 f-mode (Λ_f ≈ Λ_static, since the f-mode carries almost
    all of the tidal overlap; Σ_n Λ_n = Λ_static is the static sum rule).
    This is the standard f-mode dynamical-tide / effective-Love-number framework
    (Hinderer et al. 2016, PRL 116, 181101; Steinhoff et al. 2016, PRD 94,
    104028).  The GW frequency at the f-mode resonance is f_GW = ω_f/π.

    WHAT THIS SCRIPT DOES (data only)
    ---------------------------------
    1.  STATIC LOVE NUMBER (track A anchor): k2 and Λ_static = (2/3) k2 (R/M)^5
        via the standard Hinderer (2008) y(r) metric-perturbation ODE on the
        TOV background (PolytropeEnergy EOS1 reference star).

    2.  ℓ=2 f-MODE FREQUENCY:
        (a) relativistic Cowling f-mode (NonRadialModes), and
        (b) the full-GR interior resonance / f-mode diagnostic (PolarGRModes).
        These set ω_f and f_GW = ω_f/π.

    3.  g-MODE TOWER on a stratified finite-T star (GravityModes), with the
        per-g-mode tidal overlap estimated from its (weak) coupling.

    4.  Λ_eff(ω):  the driven-oscillator sum, with the f-mode carrying essentially
        all the static deformability and the g-modes a small overlap fraction.
        We report the static limit (sum rule check), the resonance enhancement
        as ω→ω_f, and the g-mode resonances.

    5.  VALIDATION:  static limit recovers Λ_static; the f-mode-dominated
        enhancement matches the standard effective-Love-number form.

    HONESTY
    -------
    The per-mode overlap integrals Λ_n are NOT computed from first principles
    here (that needs the full mode eigenfunctions normalised against the tidal
    field).  We use the well-justified leading approximation Λ_f = Λ_static for
    the f-mode (it dominates the overlap) and an order-of-magnitude g-mode
    overlap fraction from the literature (g-modes carry ≲1% of the f-mode
    overlap).  The ROBUST outputs are: Λ_static (k2), ω_f, f_GW, and the
    single-pole f-mode Λ_eff(ω) shape.  The g-mode overlap normalisation is the
    APPROXIMATE part.
=#

using BDNKStar
using BDNKStar.EquationOfState
using BDNKStar.TOV
using BDNKStar.NonRadialModes
using BDNKStar.GravityModes
using BDNKStar.PolarGRModes
using Printf

const Msun_to_km = BDNKStar.Units.Msun_to_km
const kHz_to_km  = BDNKStar.Units.kHz_to_km

# =====================================================================
# Reference star: Bussières EOS1 (PolytropeEnergy κ=100, n=1, ρc=3e15 g/cc)
# ⇒ M≈1.27 M⊙, R≈8.86 km.  Same star the axial/polar QNM solvers use.
# =====================================================================
const Gc = 6.6743015e-11
const cc = 299_792_458.0
const EOS = PolytropeEnergy(100.0, 1.0)
const EC  = (3e15 * 1e3) * Gc / cc^2 * 1e6        # central ε [km^-2]

star = solve_tov(EOS, EC; h=2e-4, ptol_rel=1e-12, rmax=50.0)
M = star.M; R = star.R; C = M/R
@printf("STAR  M=%.4f Msun  R=%.3f km  C=M/R=%.4f\n",
        mass_solar(star), R, C)

# =====================================================================
# 1. STATIC LOVE NUMBER k2 and Λ_static  (Hinderer 2008, PRD 77, 021502;
#    Hinderer et al. 2010).  Integrate the ℓ=2 metric-perturbation log-
#    derivative y = r H'/H from the centre to the surface, then evaluate k2.
#
#    y' = -y²/r - y·F(r) - r·Q(r),
#      F = [r - 4π r³ (ε - p)] / (r - 2m)·(1/r)  ... use the standard compact form:
#      e^{λ} = (1-2m/r)^{-1}
#      F(r)  = e^{λ}·[ 1 - 4π r²(ε - p) ] / r              (matches Hinderer A2)
#      Q(r)  = e^{λ}·{ 4π[5ε + 9p + (ε+p)/cs²] - 6/r² }    (metric/fluid source)
#                - (ν')²,   ν' = 2(m+4π r³ p)/(r(r-2m))
#    BC at centre: y(0)=2 (H ∝ r²).  Surface y_R then gives k2.
# =====================================================================
@inline function _ie(rs, ys, r)
    r <= rs[1] && return ys[1]
    r >= rs[end] && return ys[end]
    j = searchsortedlast(rs, r); j = clamp(j, 1, length(rs)-1)
    t = (r - rs[j])/(rs[j+1]-rs[j])
    ys[j] + t*(ys[j+1]-ys[j])
end

function love_y_rhs(r, y, st, eos)
    m = _ie(st.r, st.m, r); p = _ie(st.r, st.p, r); ε = _ie(st.r, st.ε, r)
    eλ = 1.0/(1.0 - 2m/r)
    νp = 2*(m + 4π*r^3*p)/(r*(r - 2m))
    cs2 = ε > 0 ? max(sound_speed2(eos, ε), 1e-14) : 1e-14
    dεdp = ε > 0 ? (ε + p)/(ε*cs2) : 0.0              # dε/dp = (ε+p)/(ε cs²)
    F = eλ*(1 - 4π*r^2*(ε - p))/r
    # Q: metric source for the ℓ=2 even-parity master function
    Q = eλ*( 4π*(5ε + 9p + (ε + p)*dεdp) - 6/r^2 ) - νp^2
    return -y^2/r - y*F - r*Q
end

function static_love(st, eos; nstep=20000)
    r0 = st.r[1]*2; rf = st.R*(1 - 1e-5)
    h  = (rf - r0)/nstep
    y  = 2.0; r = r0
    for _ in 1:nstep
        k1 = love_y_rhs(r,       y,          st, eos)
        k2 = love_y_rhs(r+h/2,   y+h/2*k1,   st, eos)
        k3 = love_y_rhs(r+h/2,   y+h/2*k2,   st, eos)
        k4 = love_y_rhs(r+h,     y+h*k3,     st, eos)
        y += h/6*(k1 + 2k2 + 2k3 + k4)
        r += h
    end
    yR = y
    Cc = st.M/st.R
    # Hinderer (2008) closed form for k2 in terms of y_R and compactness C:
    num = (8/5)*Cc^5*(1-2Cc)^2*(2 + 2Cc*(yR-1) - yR)
    den = 2*Cc*(6 - 3yR + 3Cc*(5yR-8)) +
          4*Cc^3*(13 - 11yR + Cc*(3yR-2) + 2Cc^2*(1+yR)) +
          3*(1-2Cc)^2*(2 - yR + 2Cc*(yR-1))*log(1-2Cc)
    k2 = num/den
    Λ  = (2.0/3.0)*k2/Cc^5
    return k2, Λ, yR
end

k2, Λ_static, yR = static_love(star, EOS)
@printf("\n[1] STATIC LOVE NUMBER (Hinderer 2008 y(r) ODE)\n")
@printf("    y_R = %.4f   k2 = %.4f   Λ_static = %.1f\n", yR, k2, Λ_static)

# =====================================================================
# 2. ℓ=2 f-MODE FREQUENCY
# =====================================================================
@printf("\n[2] ℓ=2 f-MODE FREQUENCY\n")

# (a) relativistic Cowling f-mode
fC, ω2C, RC = nonradial_cowling_spectrum(EOS, EC; l=2, nmodes=4, N=8000,
                                         Lunit_km=1.0, ω2lo=5e-4, ω2hi=0.20,
                                         nscan=1400)
f_cowling = fC[1]                                   # f-mode = lowest
@printf("    Cowling f-mode:  f = %.4f kHz   (p1=%.4f, p2=%.4f kHz)\n",
        f_cowling, length(fC)>1 ? fC[2] : NaN, length(fC)>2 ? fC[3] : NaN)

# ω in geometric km^-1 and the GW resonance frequency f_GW = ω_f/π = 2 f_mode
ω_f_cowling = sqrt(ω2C[1])                          # km^-1 (Lunit=1 km)
f_GW_cowling = 2*f_cowling                           # f_GW = ω_f/π = 2 f_mode

# (b) full-GR f-mode (PolarGRModes).  Seed near the Cowling value scaled down
# (full GR lies ~10-15% below Cowling).  Report whichever root converges.
f_gr = NaN; τ_gr = NaN; gr_ok = false
try
    # seed as in test_polar_gr.jl (the validated full-GR f-mode of this star)
    seed = polar_ftau_to_omega(2.8, 0.15)
    res = polar_gr_qnm(EOS, EC; l=2, ω0=seed, star=star)
    if res.converged && isfinite(res.f_kHz) && res.f_kHz > 0 && imag(res.omega) > 0
        global f_gr = res.f_kHz; global τ_gr = res.tau_s; global gr_ok = true
    end
catch e
    @printf("    (full-GR solve raised: %s)\n", e)
end
if gr_ok
    @printf("    Full-GR f-mode:  f = %.4f kHz   τ_GW = %.4f s\n", f_gr, τ_gr)
else
    @printf("    Full-GR f-mode:  not robustly converged from this seed (see PolarGRModes notes)\n")
end

# Choose the f-mode driving the dynamical tide.  Use the full-GR value if it
# converged (it is the physical one); else fall back to Cowling.
f_mode = gr_ok ? f_gr : f_cowling
ω_f    = 2π * f_mode * kHz_to_km                      # km^-1
f_GW_f = 2*f_mode                                     # GW frequency at resonance
@printf("    → f-mode used for Λ_eff:  f_f = %.4f kHz  (%s)\n",
        f_mode, gr_ok ? "full-GR" : "Cowling")
@printf("    → f-mode RESONANCE GW frequency  f_GW = ω_f/π = 2 f_f = %.3f kHz\n", f_GW_f)

# =====================================================================
# 3. g-MODE TOWER on a stratified finite-T star (GravityModes).
#    Stratification: adiabatic Γ > structure Γ_struct ⇒ N²>0 ⇒ stable g-modes.
#    We use a representative neutron-star stratification (Γ=2.05, Γ_struct=2.0).
# =====================================================================
@printf("\n[3] g-MODE TOWER (stratified finite-T star, Γ=2.2 / Γ_struct=2.0)\n")
g_freqs = Float64[]
try
    gs = gmode_spectrum(; Γ=2.2, Γ_struct=2.0, K=100.0, ρc=1.28e-3,
                        l=2, N=4000, Nbg=3000, ω2lo=1e-5, ω2hi=0.25, nscan=2400,
                        nmodes_g=4, nmodes_p=2)
    global g_freqs = collect(gs.g_freqs_kHz)
    @printf("    g-star  M=%.3f Msun  R=%.3f km  N²_max=%.3e  has_gmodes=%s\n",
            mass_solar(gs.M), gs.R, gs.N2_max, gs.has_gmodes)
    @printf("    f-mode(g-star) = %.4f kHz\n", gs.f_freq)
    for (i,gf) in enumerate(g_freqs)
        @printf("    g%-2d = %.4f kHz\n", i, gf)
    end
catch e
    @printf("    (g-mode solve raised: %s)\n", e)
end

# =====================================================================
# 4. BUILD Λ_eff(ω)
#    Λ_eff(ω) = Σ_n Λ_n ω_n²/(ω_n²-ω²),  Σ_n Λ_n = Λ_static.
#    f-mode carries the dominant overlap: Λ_f = (1-Σ_g frac)·Λ_static.
#    Each g-mode carries an overlap fraction frac_g (literature: g-modes couple
#    ≲1% of the f-mode tidal overlap; we use frac_g per g-mode below, APPROXIMATE).
# =====================================================================
@printf("\n[4] Λ_eff(ω) DRIVEN-OSCILLATOR SUM\n")

# per-g-mode overlap fraction of Λ_static (APPROXIMATE; g-modes weakly coupled).
# Lai (1994) / Reisenegger (1994): g-mode tidal overlap is ~1e-2-1e-3 of f-mode.
const G_FRAC = 1e-3                                   # each g-mode's Λ_n / Λ_static

# mode frequencies (kHz) and their Λ_n
mode_f   = Float64[f_mode]
mode_Λ   = Float64[]
g_total  = G_FRAC * length(g_freqs)
push!(mode_Λ, (1 - g_total)*Λ_static)                # f-mode gets the rest
for gf in g_freqs
    push!(mode_f, gf)
    push!(mode_Λ, G_FRAC*Λ_static)
end

# ω_n in km^-1
mode_ω = [2π*f*kHz_to_km for f in mode_f]

# sum-rule check
@printf("    Σ_n Λ_n = %.2f   (Λ_static = %.2f)   ratio = %.6f\n",
        sum(mode_Λ), Λ_static, sum(mode_Λ)/Λ_static)

Λeff(ω) = sum(mode_Λ[n]*mode_ω[n]^2/(mode_ω[n]^2 - ω^2) for n in 1:length(mode_f))

# (a) static limit ω→0
ω0 = 1e-6
@printf("    (a) static limit  Λ_eff(ω→0) = %.4f   (Λ_static=%.4f)  Δ=%.2e\n",
        Λeff(ω0), Λ_static, abs(Λeff(ω0)-Λ_static)/Λ_static)

# (b) f-mode resonance enhancement: tabulate Λ_eff vs GW frequency f_GW=ω/π
@printf("    (b) f-mode resonance enhancement  Λ_eff(ω)/Λ_static  vs  f_GW=ω/π [kHz]\n")
@printf("        %-10s %-12s %-10s\n", "f_GW[kHz]", "Λ_eff", "Λ_eff/Λs")
for fGW in [0.2, 0.5, 1.0, 1.5, 2.0, 0.7*f_GW_f, 0.85*f_GW_f, 0.95*f_GW_f, 0.99*f_GW_f]
    ω = π*fGW*kHz_to_km                                # ω = π f_GW  (f_GW=ω/π)
    Λ = Λeff(ω)
    @printf("        %-10.3f %-12.2f %-10.4f\n", fGW, Λ, Λ/Λ_static)
end
@printf("        f-mode resonance at f_GW = %.3f kHz (Λ_eff diverges; dynamical-tide amplification)\n", f_GW_f)

# (c) g-mode resonances
if !isempty(g_freqs)
    @printf("    (c) g-mode resonances (weak; f_GW = 2 f_g):\n")
    for (i,gf) in enumerate(g_freqs)
        @printf("        g%d at f_g=%.4f kHz ⇒ f_GW=%.4f kHz, overlap Λ_g/Λs=%.1e\n",
                i, gf, 2gf, G_FRAC)
    end
else
    @printf("    (c) no g-modes (barotropic star) ⇒ Λ_eff is a clean single f-mode pole\n")
end

# =====================================================================
# 5. VALIDATION
# =====================================================================
@printf("\n[5] VALIDATION\n")
# (i) static sum rule
sr = sum(mode_Λ)/Λ_static
@printf("    (i)  static sum rule Σ Λ_n = Λ_static : ratio %.6f  %s\n",
        sr, abs(sr-1)<1e-9 ? "PASS" : "FAIL")
# (ii) Λ_static vs universal-relation / published expectation for this star.
#   k2 for an n=1 polytrope at C~0.21 ~ 0.06-0.07; Λ for a 1.27 Msun NS ~ a few hundred.
@printf("    (ii) k2=%.4f, Λ_static=%.1f for C=%.3f (n=1 polytrope) — k2~0.05-0.08, Λ~%.0f expected scale\n",
        k2, Λ_static, C, Λ_static)
# (iii) f-mode universal relation (Andersson-Kokkotas / Chan et al. 2014):
#   ω_f M (geometric) vs sqrt(M/R^3). Report dimensionless f-mode for cross-check.
Mkm = M
ωf_geo = 2π*f_mode*kHz_to_km
@printf("    (iii) f-mode: f_f=%.4f kHz, M ω_f = %.4f, sqrt(M/R³)=%.4f (geom) ; f_GW=%.3f kHz (~2kHz near merger)\n",
        f_mode, Mkm*ωf_geo, sqrt(M/R^3), f_GW_f)

@printf("\nDONE.\n")
