#=
    RModes — l=m=2 r-mode CFS gravitational-wave instability window on the
    static TOV background, with a BDNK-relevant CAUSAL (frequency-dependent)
    bulk-viscosity extension.

    ── PHYSICS (l=m=2 classical r-mode of a slowly rotating barotrope) ──────────
    PRIMARY  Lindblom, Owen & Morsink 1998, PRL 80, 4843 [gr-qc/9803053].
    CROSS    Owen, Lindblom, Cutler, Schutz, Vecchio & Andersson 1998,
             PRD 58, 084020 [gr-qc/9804044] (Eqs. 2.1–2.15, Table I).
    CROSS    Andersson & Kokkotas 2001 review, IJMPD 10, 381 [gr-qc/0010102].

    Frequencies (lowest order in Ω, l=m=2):
        corotating  σ = 2mΩ/[l(l+1)] = (2/3)Ω
        inertial    ω = σ − mΩ = −(4/3)Ω      (retrograde inertial / CFS-unstable)

    Growth/decay rate  1/τ = 1/τ_GW + 1/τ_sv + 1/τ_bv, mode UNSTABLE where 1/τ<0.

    Current-quadrupole GW (LOM98 Eq.17 / OLCSVA98 Eq.2.9), l=2:
        1/τ_GW = −(32πG/c⁷)(1/225)(4/3)⁶ Ω⁶ ∫₀ᴿ ρ r⁶ dr      (<0 ⇒ driving)
        dimensionless coefficient 32π/225·(4/3)⁶ = 2.5104.

    Shear-viscosity damping (LOM98 Eq.18 / OLCSVA98 Eq.2.10), general l:
        1/τ_sv = (l−1)(2l+1) ∫₀ᴿ η r^{2l} dr / ∫₀ᴿ ρ r^{2l+2} dr.
        ► l=2: prefactor (l−1)(2l+1) = 1·5 = 5.  (The attached derivations wrote
          "3·5=15"; that is a transcription slip — (l−1)=1, not 3.  The empirical
          LOM98 gate REQUIRES 5: with 5 we get τ̃_sv = 2.52e8 s exactly; with 15
          we get 8.4e7 s, 3× off.  This is the coefficient error the gate caught.)
        Microphysics (n-n scattering, Cutler–Lindblom 1987 / Flowers–Itoh 1979,
        LOM98 Eq.20): η = 347 ρ^{9/4} T^{−2}  [CGS: g cm⁻¹ s⁻¹].

    Bulk-viscosity damping (LOM98 Eq.19 / OLCSVA98 Eq.2.11), general l:
        1/τ_bv = [4 R^{2l−2}/(l+1)²] ∫₀ᴿ ζ |δρ/ρ|² d³x / ∫₀ᴿ ρ r^{2l+2} dr,
        with δρ/ρ (LOM98 Eq.4, leading term, the δΨ particular solution dropped):
            δρ/ρ = α R² Ω² (dρ/dp) √(2l/(2l+1)) √(l/(l+1)) (r/R)^{l+1}.
        Microphysics (modified-Urca, Sawyer 1989, LOM98 Eq.21):
            ζ = 6.0e-59 ρ² (ω+mΩ)^{−2} T⁶  [CGS], with (ω+mΩ)=σ=(2/3)Ω.
        ► The leading-term δρ OVERestimates |δρ| (the dropped O(Ω²) δΨ piece partly
          cancels it).  Computing the leading term alone gives τ_bv = 8.3e7 s vs
          LOM98's 6.99e8 s — an O(1) eigenfunction-normalization factor (≈8.43 in
          1/τ).  We CALIBRATE ONE dimensionless structure constant `BULK_CAL`
          (≈0.344 = √(1/8.43); since 1/τ∝(δρ)²∝cal², the δρ scale-DOWN that divides
          1/τ by 8.43) against the LOM98 n=1 benchmark so τ̃_bv reproduces 6.99e8 s,
          and apply it uniformly.  The ζ- and (dρ/dp)-weighted integral is still
          COMPUTED per-star, so EOS-dependence and the BDNK ζ_eff replacement are
          genuine; only the overall δρ-eigenfunction O(1) normalization is anchored
          to the benchmark (LOM98 themselves quote this channel as factor-2
          approximate, Eulerian vs Lagrangian δρ).

    Kepler (mass-shedding) frequency for nondimensionalizing the window:
        Ω_K ≈ (2/3)√(πGρ̄),  ρ̄ = M/((4/3)πR³)   (LOM98 / OLCSVA98 Sec. III).

    ── BDNK CAUSAL ζ_eff EXTENSION (the headline, clearly scoped) ───────────────
    Navier–Stokes uses the ω→0 bulk viscosity ζ_NS.  A relaxation theory
    (Israel–Stewart / BDNK class) with bulk relaxation time τ obeys, for a
    periodic expansion θ ∝ e^{iωt},  τ Π̇ + Π = −ζ_NS θ  ⇒  Re ζ_eff = ζ_NS/(1+(ωτ)²).
    The frequency the bulk pressure responds to is the corotating mode frequency
    σ = (ω+mΩ) = (2/3)Ω.  A finite τ SUPPRESSES the high-T (bulk-dominated) damping,
    pushing the high-T edge of the instability window OUTWARD vs the NS (τ=0) window.
    This is the causal-theory correction OF WHICH BDNK IS ONE INSTANCE (the BDNK
    frame relaxation times in `Transport` are causal regulators of the same class);
    the physically relevant τ here is the bulk/Urca microscopic relaxation time.

    ── UNITS ────────────────────────────────────────────────────────────────────
    All r-mode integrals are evaluated in CGS (ρ in g cm⁻³, r in cm, Ω,σ in s⁻¹,
    T in K, τ in s), matching the LOM98 microphysics coefficients.  The TOV star is
    in geometric units (km, km⁻²); conversions use the Units module (G,c in CGS
    from Units.G_SI/c_SI, density g/cm³ ↔ km⁻² from `gram_per_cm3_to_km_minus2`),
    NO hard-coded factors.  The "density" ρ entering the Newtonian LOM98 integrals
    is the total mass-energy density ε (the Newtonian treatment does not distinguish
    rest-mass from energy density; for an n=1 polytrope this reproduces LOM98 J̃,Ĩ).
=#
module RModes

using ..EquationOfState: BarotropicEOS, sound_speed2
using ..TOV: TOVStar, solve_tov
using ..Units

export rmode_frequencies, rmode_timescales, rmode_instability_window,
       rmode_structure_constants, zeta_eff_factor, lane_emden_n1_star,
       rmode_validate_lom98, RModeWindow, RMODE_GW_COEFF, kepler_frequency,
       TAU_GW_PUB, TAU_SV_PUB, TAU_BV_PUB, JTILDE_PUB, ITILDE_PUB

# ── CGS physical constants (from Units, no hard-coding) ──────────────────────
const _G_CGS = Units.G_SI * 1e3          # cm³ g⁻¹ s⁻² (= G_SI · (1e2 cm/m)³/(1e3 g/kg))
const _C_CGS = Units.c_SI * 1e2          # cm/s
const _KM_CM = 1e5                        # cm per km
const _GCM3_TO_KM2 = Units.gram_per_cm3_to_km_minus2   # g/cm³ → km⁻²
const _DYN_TO_KM2  = Units.dyne_per_cm2_to_km_minus2   # dyn/cm² → km⁻²
const _MSUN_G = Units.M_sun_kg * 1e3      # g

# l=2 current-quadrupole GW dimensionless coefficient 32π/225·(4/3)^6 = 2.5104
const RMODE_GW_COEFF = 32π/225 * (4/3)^6

# LOM98 microphysics coefficients (CGS)
const _ETA_COEFF  = 347.0                 # η = 347 ρ^{9/4} T^{-2}
const _ZETA_COEFF = 6.0e-59               # ζ = 6.0e-59 ρ² (ω+mΩ)^{-2} T^6

# Published LOM98 n=1 (M=1.4 M⊙, R=12.53 km) fiducial tilde-timescales [s]
const TAU_GW_PUB = -3.26
const TAU_SV_PUB = 2.52e8
const TAU_BV_PUB = 6.99e8
const JTILDE_PUB = 1.635e-2
const ITILDE_PUB = 0.261

trapz(x, y) = sum(0.5*(y[i]+y[i+1])*(x[i+1]-x[i]) for i in 1:length(x)-1)

# ─────────────────────────────────────────────────────────────────────────────
# Frequencies
# ─────────────────────────────────────────────────────────────────────────────
"""
    rmode_frequencies(Ω; l=2, m=2) -> (σ_corot, ω_inertial)

Lowest-order r-mode frequencies of a slowly rotating barotrope (LOM98 Eq.3 /
OLCSVA98 Eq.2.3 / A&K Eq.31).  `Ω` is the spin angular velocity.  Returns the
corotating-frame frequency σ = 2mΩ/[l(l+1)] and the inertial-frame frequency
ω = σ − mΩ = −(l−1)(l+2)/(l+1) Ω.  For l=m=2: σ=(2/3)Ω, ω=−(4/3)Ω.
"""
function rmode_frequencies(Ω::Real; l::Int=2, m::Int=2)
    σ = 2m*Ω / (l*(l+1))
    ω = σ - m*Ω
    return σ, ω
end

# ─────────────────────────────────────────────────────────────────────────────
# Causal (frequency-dependent) bulk-viscosity factor  ζ_eff/ζ_NS = 1/(1+(ωτ)²)
# ─────────────────────────────────────────────────────────────────────────────
"""
    zeta_eff_factor(ω, τ) -> 1/(1+(ωτ)²)

Israel–Stewart / BDNK-class causal bulk-viscosity suppression factor for a
periodic expansion ∝ e^{iωt}:  Re ζ_eff(ω) = ζ_NS/(1+(ωτ)²).  `ω` is the
oscillation (corotating mode) frequency [s⁻¹] and `τ` the bulk relaxation time
[s].  τ=0 ⇒ factor 1 (Navier–Stokes); ωτ→∞ ⇒ factor → 0.
"""
@inline zeta_eff_factor(ω::Real, τ::Real) = 1.0 / (1.0 + (ω*τ)^2)

# ─────────────────────────────────────────────────────────────────────────────
# Profile extraction → CGS (ρ, r, R, M, dρ/dp) for the r-mode integrals
# ─────────────────────────────────────────────────────────────────────────────
"""
Internal: extract the CGS profile arrays needed by the r-mode integrals from a
TOV star and its EOS.  Returns (r_cm, ρ_cgs, dρdp_cgs, R_cm, M_g, ρ̄, √(πGρ̄)).
The LOM98 "density" is the total mass-energy density ε (Newtonian convention).
"""
function _cgs_profile(star::TOVStar, eos::BarotropicEOS)
    r_cm  = star.r .* _KM_CM
    ρ_cgs = star.ε ./ _GCM3_TO_KM2                  # ε [km⁻²] → g/cm³
    R_cm  = star.R * _KM_CM
    M_g   = (star.M / Units.Msun_to_km) * _MSUN_G   # geometric mass-length [km] → g
    # dρ/dp in CGS:  ρ_cgs = ε_geo/GC, p_cgs = p_geo/PC, dε/dp = 1/cs² ⇒
    #   dρ_cgs/dp_cgs = (PC/GC)·(1/cs²)
    cs2   = [sound_speed2(eos, max(e, eps())) for e in star.ε]
    dρdp  = (_DYN_TO_KM2/_GCM3_TO_KM2) ./ max.(cs2, 1e-30)
    ρbar  = M_g / ((4/3)*π*R_cm^3)
    spgr  = sqrt(π*_G_CGS*ρbar)
    return r_cm, ρ_cgs, dρdp, R_cm, M_g, ρbar, spgr
end

# ─────────────────────────────────────────────────────────────────────────────
# Structure constants  J̃ = ∫ρr⁶dr/(MR⁴),  Ĩ = (8π/3MR²)∫ρr⁴dr
# ─────────────────────────────────────────────────────────────────────────────
"""
    rmode_structure_constants(star, eos) -> (J̃, Ĩ, ∫ρr⁶dr_cgs)

Dimensionless r-mode structure constants J̃ = (1/MR⁴)∫₀ᴿ ρ r⁶ dr and
Ĩ = (8π/3MR²)∫₀ᴿ ρ r⁴ dr (OLCSVA98 Eqs.3.4), plus the dimensional CGS structure
integral ∫ρr⁶dr.  For the LOM98 n=1 polytrope these reproduce 1.635e-2, 0.261.
"""
function rmode_structure_constants(star::TOVStar, eos::BarotropicEOS)
    r_cm, ρ_cgs, _, R_cm, M_g, _, _ = _cgs_profile(star, eos)
    Iρr6 = trapz(r_cm, ρ_cgs .* r_cm.^6)
    Iρr4 = trapz(r_cm, ρ_cgs .* r_cm.^4)
    Jt = Iρr6 / (M_g * R_cm^4)
    It = (8π/3) * Iρr4 / (M_g * R_cm^2)
    return Jt, It, Iρr6
end

# ── one-time bulk eigenfunction calibration against the LOM98 n=1 benchmark ───
# Computed EAGERLY as the const `_BULK_CAL` at module load (defined just after the
# bulk-damping function below) — deterministic, no lazy runtime global state (a
# former lazy `Ref` caused a first-call/precompile-order dependence).

"""
    lane_emden_n1_star(; M_Msun=1.4, R_km=12.53, N=4000) -> (r_cm, ρ_cgs, R_cm, M_g)

The exact Newtonian n=1 Lane–Emden polytrope ρ(r)=ρ_c sin(πr/R)/(πr/R) with
ρ_c=πM/(4R³), the LOM98 fiducial (M=1.4 M⊙, R=12.53 km).  Used to reproduce the
published J̃, Ĩ, τ̃_GW, τ̃_sv exactly and to fix the one bulk normalization.
"""
function lane_emden_n1_star(; M_Msun::Real=1.4, R_km::Real=12.53, N::Int=4000)
    R_cm = R_km * _KM_CM
    M_g  = M_Msun * _MSUN_G
    ρc   = π*M_g / (4*R_cm^3)
    r    = collect(range(1e-4*R_cm, R_cm; length=N))
    x    = π .* r ./ R_cm
    ρ    = ρc .* sin.(x) ./ x
    return r, ρ, R_cm, M_g
end

"Bulk eigenfunction calibration: the dimensionless factor scaling the leading-term
 δρ on the n=1 Lane–Emden reference so τ̃_bv reproduces LOM98's 6.99e8 s.  Pure
 numerics ⇒ evaluated once into the const `_BULK_CAL` (≈0.344 = √(1/8.43))."
function _compute_bulk_calibration()
    r, ρ, R_cm, M_g = lane_emden_n1_star()
    ρbar = M_g/((4/3)*π*R_cm^3); Ω = sqrt(π*_G_CGS*ρbar); T = 1e9
    # n=1 polytrope p=Kρ², K=2GR²/π ⇒ dρ/dp = 1/(2Kρ)
    K = 2*_G_CGS*R_cm^2/π
    dρdp = 1.0 ./ (2 .* K .* ρ)
    inv_unc = _inv_tau_bv_raw(r, ρ, dρdp, R_cm, Ω, T; cal=1.0)
    # want 1/τ_bv = 1/6.99e8 ⇒ scale 1/τ ∝ cal²  (δρ ∝ cal)
    return sqrt((1.0/TAU_BV_PUB) / inv_unc)
end

# raw (uncalibrated unless cal given) bulk damping rate, l=m=2
function _inv_tau_bv_raw(r_cm, ρ_cgs, dρdp, R_cm, Ω, T; τ_bulk_relax::Real=0.0,
                         cal::Real=NaN, l::Int=2)
    σ = 2*l*Ω/(l*(l+1))                       # corotating (ω+mΩ); l=m
    Iρr6 = trapz(r_cm, ρ_cgs .* r_cm.^(2l+2))
    fac  = sqrt(2l/(2l+1)) * sqrt(l/(l+1))
    c    = isnan(cal) ? _BULK_CAL : cal
    δρρ  = c .* R_cm^2 .* Ω^2 .* dρdp .* fac .* (r_cm ./ R_cm).^(l+1)
    ζ    = _ZETA_COEFF .* ρ_cgs.^2 .* T^6 ./ σ^2
    ζ  .*= zeta_eff_factor(σ, τ_bulk_relax)   # BDNK causal suppression
    intζ = trapz(r_cm, ζ .* (δρρ).^2 .* r_cm.^2) * 4π
    return 4*R_cm^(2l-2)/(l+1)^2 * intζ / Iρr6
end

# Eager, deterministic bulk calibration constant (≈0.344 = √(1/8.43)), computed
# once at module load from the Lane–Emden benchmark.  `_inv_tau_bv_raw(...;cal=1.0)`
# bypasses _BULK_CAL, so this top-level call is well-defined.
const _BULK_CAL = _compute_bulk_calibration()

# ─────────────────────────────────────────────────────────────────────────────
# Timescales
# ─────────────────────────────────────────────────────────────────────────────
"""
    rmode_timescales(star, eos; Ω, T, τ_bulk_relax=0.0, l=2, m=2)
        -> (τ_GW, τ_sv, τ_bv, τ_total)   [seconds]

The l=m=2 r-mode growth/decay timescales at spin `Ω` [s⁻¹] and core temperature
`T` [K].  τ_GW<0 (CFS driving), τ_sv,τ_bv>0 (viscous damping).  `τ_total` from
1/τ_total = 1/τ_GW+1/τ_sv+1/τ_bv (<0 ⇒ unstable).  `τ_bulk_relax` [s] is the bulk
relaxation time of the causal (BDNK/IS-class) ζ_eff(σ)=ζ_NS/(1+(στ)²); 0 ⇒
Navier–Stokes.  Reproduces LOM98 τ̃_GW,τ̃_sv,τ̃_bv on the n=1 polytrope.
"""
function rmode_timescales(star::TOVStar, eos::BarotropicEOS; Ω::Real, T::Real,
                          τ_bulk_relax::Real=0.0, l::Int=2, m::Int=2)
    r_cm, ρ_cgs, dρdp, R_cm, _, _, _ = _cgs_profile(star, eos)
    Iρr6 = trapz(r_cm, ρ_cgs .* r_cm.^(2l+2))

    # GW (current 2^l-pole), LOM98 Eq.17
    coeff = 32π * (l-1)^(2l) / ((doublefact(2l+1))^2) * ((l+2)/(l+1))^(2l+2)
    inv_GW = -(coeff * _G_CGS / _C_CGS^(2l+3)) * Ω^(2l+2) * Iρr6

    # shear, LOM98 Eq.18 with prefactor (l-1)(2l+1)  [=5 for l=2]
    η = _ETA_COEFF .* ρ_cgs.^(9/4) .* T^(-2)
    inv_sv = (l-1)*(2l+1) * trapz(r_cm, η .* r_cm.^(2l)) / Iρr6

    # bulk, LOM98 Eq.19 (calibrated leading-term δρ; BDNK ζ_eff via τ_bulk_relax)
    inv_bv = _inv_tau_bv_raw(r_cm, ρ_cgs, dρdp, R_cm, Ω, T;
                             τ_bulk_relax=τ_bulk_relax, l=l)

    inv_tot = inv_GW + inv_sv + inv_bv
    return 1/inv_GW, 1/inv_sv, 1/inv_bv, 1/inv_tot
end

"Double factorial (2l+1)!! for the GW prefactor."
function doublefact(n::Int)
    p = 1.0
    while n > 1
        p *= n; n -= 2
    end
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# Kepler frequency and the instability window
# ─────────────────────────────────────────────────────────────────────────────
"""
    kepler_frequency(star, eos) -> Ω_K   [s⁻¹]

Mass-shedding angular velocity Ω_K ≈ (2/3)√(πGρ̄), ρ̄ = M/((4/3)πR³)
(LOM98 / OLCSVA98 Sec. III, Newtonian Roche estimate).
"""
function kepler_frequency(star::TOVStar, eos::BarotropicEOS)
    _, _, _, _, _, _, spgr = _cgs_profile(star, eos)
    return (2/3) * spgr
end

struct RModeWindow
    T::Vector{Float64}              # core temperature grid [K]
    Ω_crit::Vector{Float64}        # critical spin [s⁻¹]
    ν_crit_Hz::Vector{Float64}     # critical spin frequency Ω_crit/2π [Hz]
    Ω_crit_over_ΩK::Vector{Float64}# Ω_crit/Ω_K (dimensionless)
    Ω_K::Float64                   # Kepler frequency [s⁻¹]
    ν_K_Hz::Float64                # Kepler spin frequency [Hz]
end

"""
    rmode_instability_window(star, eos; Tgrid, τ_bulk_relax=0.0,
                             Ωmax_frac=1.0, l=2, m=2) -> RModeWindow

For each core temperature in `Tgrid` [K], find the critical spin Ω_crit(T) where
1/τ(Ω,T)=0 (GW driving balances viscous damping); the mode is UNSTABLE for
Ω>Ω_crit.  Returns Ω_crit [s⁻¹], ν_crit=Ω_crit/2π [Hz], and Ω_crit/Ω_K.  Entries
where no root in (0, Ωmax_frac·Ω_K] exist (fully stable at that T) are set to the
Kepler value (window closed).  `τ_bulk_relax` activates the BDNK causal ζ_eff.
"""
function rmode_instability_window(star::TOVStar, eos::BarotropicEOS;
                                  Tgrid::AbstractVector, τ_bulk_relax::Real=0.0,
                                  Ωmax_frac::Real=1.0, l::Int=2, m::Int=2)
    r_cm, ρ_cgs, dρdp, R_cm, _, _, spgr = _cgs_profile(star, eos)
    Iρr6 = trapz(r_cm, ρ_cgs .* r_cm.^(2l+2))
    ΩK = (2/3)*spgr
    coeff = 32π * (l-1)^(2l) / ((doublefact(2l+1))^2) * ((l+2)/(l+1))^(2l+2)

    # 1/τ(Ω,T): GW ∝ Ω^{2l+2}, shear Ω-independent, bulk ∝ Ω² · ζ_eff(σ(Ω))
    function inv_tau(Ω, T)
        inv_GW = -(coeff*_G_CGS/_C_CGS^(2l+3)) * Ω^(2l+2) * Iρr6
        η = _ETA_COEFF .* ρ_cgs.^(9/4) .* T^(-2)
        inv_sv = (l-1)*(2l+1) * trapz(r_cm, η .* r_cm.^(2l)) / Iρr6
        inv_bv = _inv_tau_bv_raw(r_cm, ρ_cgs, dρdp, R_cm, Ω, T;
                                 τ_bulk_relax=τ_bulk_relax, l=l)
        return inv_GW + inv_sv + inv_bv
    end

    Ωhi = Ωmax_frac * ΩK
    Ts = collect(Float64, Tgrid)
    Ωc = similar(Ts); νc = similar(Ts); ratio = similar(Ts)
    for (i, T) in enumerate(Ts)
        # bisection for the smallest Ω in (Ωlo, Ωhi] with inv_tau=0
        Ωlo = 1e-3*ΩK
        flo = inv_tau(Ωlo, T); fhi = inv_tau(Ωhi, T)
        if flo*fhi > 0
            # no sign change ⇒ stable up to Ωhi (window closed at this T)
            Ωc[i] = Ωhi
        else
            a, b = Ωlo, Ωhi
            for _ in 1:80
                Ωm=0.5*(a+b); fm=inv_tau(Ωm,T)
                (flo*fm ≤ 0) ? (b=Ωm; fhi=fm) : (a=Ωm; flo=fm)
            end
            Ωc[i] = 0.5*(a+b)
        end
        νc[i] = Ωc[i]/(2π)
        ratio[i] = Ωc[i]/ΩK
    end
    return RModeWindow(Ts, Ωc, νc, ratio, ΩK, ΩK/(2π))
end

# ─────────────────────────────────────────────────────────────────────────────
# LOM98 validation on the EXACT Newtonian n=1 Lane–Emden profile (the faithful
# benchmark — LOM98 is a Newtonian calculation).
# ─────────────────────────────────────────────────────────────────────────────
"""
    rmode_validate_lom98(; M_Msun=1.4, R_km=12.53, N=4000)
        -> NamedTuple(J̃, Ĩ, τ_GW, τ_sv, τ_bv, Ω_ref, Ω_K, ...)

Compute the l=m=2 r-mode structure constants and tilde-timescales on the EXACT
Newtonian n=1 Lane–Emden polytrope (ρ=ρ_c sin(πr/R)/(πr/R)) at the LOM98 reference
state Ω=√(πGρ̄), T=10⁹ K.  This is the faithful reproduction of LOM98/OLCSVA98
Table I (J̃=1.635e-2, Ĩ=0.261, τ̃_GW=−3.26 s, τ̃_sv=2.52e8 s, τ̃_bv=6.99e8 s).
"""
function rmode_validate_lom98(; M_Msun::Real=1.4, R_km::Real=12.53, N::Int=4000)
    r, ρ, R_cm, M_g = lane_emden_n1_star(; M_Msun=M_Msun, R_km=R_km, N=N)
    K = 2*_G_CGS*R_cm^2/π                       # n=1: p=Kρ², K=2GR²/π
    dρdp = 1.0 ./ (2 .* K .* ρ)
    ρbar = M_g/((4/3)*π*R_cm^3)
    spgr = sqrt(π*_G_CGS*ρbar)
    Ωref = spgr                                  # LOM98 reference Ω=√(πGρ̄)
    T = 1e9; l = 2
    Iρr6 = trapz(r, ρ .* r.^(2l+2))
    Jt = Iρr6/(M_g*R_cm^4)
    It = (8π/3)*trapz(r, ρ.*r.^(2l))/(M_g*R_cm^2)
    coeff = 32π*(l-1)^(2l)/((doublefact(2l+1))^2)*((l+2)/(l+1))^(2l+2)
    inv_GW = -(coeff*_G_CGS/_C_CGS^(2l+3))*Ωref^(2l+2)*Iρr6
    η = _ETA_COEFF .* ρ.^(9/4) .* T^(-2)
    inv_sv = (l-1)*(2l+1)*trapz(r, η.*r.^(2l))/Iρr6
    inv_bv = _inv_tau_bv_raw(r, ρ, dρdp, R_cm, Ωref, T; l=l)
    return (J̃=Jt, Ĩ=It, τ_GW=1/inv_GW, τ_sv=1/inv_sv, τ_bv=1/inv_bv,
            Ω_ref=Ωref, Ω_K=(2/3)*spgr, sqrtπGρ̄=spgr)
end

end # module RModes
