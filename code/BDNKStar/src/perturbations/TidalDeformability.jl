#=
    TidalDeformability — ℓ=2 tidal Love number k₂, dimensionless tidal
    deformability Λ, and the frequency-dependent (dynamical) effective tidal
    deformability Λ_eff(ω).  Completes the I-Love-Q triangle (the "Love" leg,
    the "I" already in SlowRotation) and supplies the inviscid-GR tidal baseline
    that the BDNK viscous QNM machinery extends.

    STATIC  (Postnikov–Prakash–Lattimer 2010, PRD 82 024016 [arXiv:1004.5098],
    Eqs. eqH/eq:Q/eq:met/eq:y; k₂ closed form Hinderer 2008 ApJ 677,1216 + the
    arXiv:0711.2420 erratum).  Integrate the first-order Riccati ODE for
    y(r)=rH'/H from the centre y(0)=2:

        r y' + y² + y·F(r) + r²·Q(r) = 0,

    with  e^λ = 1/(1−2m/r),
        F(r) = e^λ [1 + 4π r²(p−ε)],
        Q(r) = 4π e^λ (5ε + 9p + (ε+p)/c_s²) − 6 e^λ/r² − (ν')²,
        ν'  = 2 e^λ (m + 4π r³ p)/r²,     c_s² = dp/dε.
    Surface y_R = y(R).  With C = M/R,

        k₂ = (8C⁵/5)(1−2C)²[2 + 2C(y_R−1) − y_R]
             / { 2C[6 − 3y_R + 3C(5y_R−8)]
               + 4C³[13 − 11y_R + C(3y_R−2) + 2C²(1+y_R)]
               + 3(1−2C)²[2 − y_R + 2C(y_R−1)] ln(1−2C) },
        Λ = (2/3) k₂ C^{-5}.

    DYNAMICAL  (Hinderer et al. 2016 PRL 116,181101 / Steinhoff et al. 2016 PRD
    94,104028 effective-Love-number form).  As the tidal driving frequency ω
    sweeps toward a stellar mode ω_n, the effective tidal response takes the
    driven-oscillator form

        Λ_eff(ω) = Σ_n Λ_n · ω_n² / (ω_n² − ω²),     Σ_n Λ_n = Λ_static,

    dominated by the ℓ=2 f-mode (Λ_f ≈ Λ_static).  Λ_eff(ω→0)=Λ_static and
    amplifies toward the f-mode resonance.
=#
module TidalDeformability

using ..EquationOfState
using ..TOV
using ..Units: kHz_to_km

export tidal_love_number, lambda_from_k2, lambda_eff, lambda_eff_modes,
       ibar_yagi_yunes, LambdaEffModes, YY_ILOVE_COEFFS, fmode_resonance_fGW_kHz

# ── k₂ closed form (Hinderer 2008 + erratum; matches PPL2010 eq:k2C) ──────────
"k₂ given compactness C=M/R and the surface log-derivative y_R = R H'(R)/H(R)."
function k2_of(C::Float64, y::Float64)
    num = (8C^5/5)*(1-2C)^2*(2 + 2C*(y-1) - y)
    den = 2C*(6 - 3y + 3C*(5y-8)) +
          4C^3*(13 - 11y + C*(3y-2) + 2C^2*(1+y)) +
          3*(1-2C)^2*(2 - y + 2C*(y-1))*log(1-2C)
    return num/den
end

"Λ = (2/3) k₂ C^{-5}."
lambda_from_k2(k2::Float64, C::Float64) = (2.0/3.0)*k2/C^5

# ── y(r) integration on the TOV grid (RK4, midpoint EOS lookups) ─────────────
# RHS:  y' = -[ y² + y F + r² Q ] / r
function _yR(eos, star)
    r = star.r; m = star.m
    @inline function FQ(rr, mm, p, ε, cs2)
        eλ = 1/(1 - 2mm/rr)
        F  = eλ*(1 + 4π*rr^2*(p - ε))
        νp = 2*eλ*(mm + 4π*rr^3*p)/rr^2
        Q  = 4π*eλ*(5ε + 9p + (ε+p)/cs2) - 6*eλ/rr^2 - νp^2
        return F, Q
    end
    @inline function dydr(rr, yy, mm, p, ε, cs2)
        F, Q = FQ(rr, mm, p, ε, cs2)
        return -(yy^2 + yy*F + rr^2*Q)/rr
    end
    @inline function state(i)
        εi = star.ε[i]; pi = star.p[i]
        cs2 = εi > 0 ? max(sound_speed2(eos, εi), 1e-14) : 1.0
        return pi, εi, cs2
    end
    @inline function state_mid(i)
        εi = 0.5*(star.ε[i] + star.ε[i+1])
        pmid = εi > 0 ? pressure(eos, εi) : 0.0
        cs2  = εi > 0 ? max(sound_speed2(eos, εi), 1e-14) : 1.0
        mmid = 0.5*(m[i] + m[i+1])
        return mmid, pmid, εi, cs2
    end

    y = 2.0
    n = length(r)
    # integrate up to the last node with ε>0, then carry y to R across the
    # (tiny) crust edge with the vacuum equation (p=ε=0 ⇒ F=1, Q=-(ν')²).
    ilast = n
    while ilast > 1 && star.ε[ilast] ≤ 0
        ilast -= 1
    end
    @inbounds for i in 1:ilast-1
        hh = r[i+1] - r[i]
        p0,ε0,c0 = state(i)
        k1 = dydr(r[i], y, m[i], p0, ε0, c0)
        mm,pm,εm,cm = state_mid(i)
        rm = r[i] + hh/2
        k2 = dydr(rm, y + hh/2*k1, mm, pm, εm, cm)
        k3 = dydr(rm, y + hh/2*k2, mm, pm, εm, cm)
        p1,ε1,c1 = state(i+1)
        k4 = dydr(r[i+1], y + hh*k3, m[i+1], p1, ε1, c1)
        y += hh/6*(k1 + 2k2 + 2k3 + k4)
    end
    @inbounds for i in ilast:n-1
        hh = r[i+1] - r[i]
        eλ = 1/(1 - 2m[i]/r[i])
        νp = 2*eλ*m[i]/r[i]^2
        y += hh*(-(y^2 + y*eλ)/r[i] + r[i]*νp^2)
    end
    return y
end

"""
    tidal_love_number(eos, εc; h=1e-3, kwargs...) -> (k2, Λ, y_R)

ℓ=2 tidal Love number `k2`, dimensionless tidal deformability `Λ=(2/3)k2 C^{-5}`
(C=M/R), and the surface log-derivative `y_R = R H'(R)/H(R)` for the barotropic
`eos` at central energy density `εc` [km⁻²].  Integrates the Postnikov–Prakash–
Lattimer (2010) y(r) Riccati ODE on the TOV background by RK4.

`tidal_love_number(eos, star)` reuses a precomputed `TOVStar`.
"""
function tidal_love_number(eos::BarotropicEOS, εc::Float64; h::Float64=1e-3,
                           ptol_rel::Float64=1e-10, rmax::Float64=100.0)
    star = solve_tov(eos, εc; h=h, ptol_rel=ptol_rel, rmax=rmax)
    return tidal_love_number(eos, star)
end

function tidal_love_number(eos::BarotropicEOS, star::TOVStar)
    C  = star.M/star.R
    yR = _yR(eos, star)
    k2 = k2_of(C, yR)
    Λ  = lambda_from_k2(k2, C)
    return (k2, Λ, yR)
end

# ── DYNAMICAL: driven-oscillator effective tidal deformability Λ_eff(ω) ──────
"""
    LambdaEffModes(mode_f_kHz, mode_Λ)

A mode decomposition for the dynamical tidal response: cyclic mode frequencies
`mode_f_kHz` [kHz] (the f-mode first, then any g-modes) and the per-mode tidal
overlaps `mode_Λ` (Σ = Λ_static by the static sum rule).  Build with
`lambda_eff_modes`; evaluate with `lambda_eff(modes, ω)`.
"""
struct LambdaEffModes
    mode_f_kHz::Vector{Float64}   # cyclic frequencies f_n [kHz]
    mode_Λ::Vector{Float64}       # per-mode overlaps Λ_n  (Σ = Λ_static)
    mode_ω::Vector{Float64}       # ω_n = 2π f_n [km⁻¹]
end

"""
    lambda_eff_modes(Λ_static, f_f_kHz; g_freqs_kHz=Float64[], g_frac=1e-3)

Assemble the dynamical-tide mode decomposition.  The ℓ=2 f-mode `f_f_kHz`
carries the dominant tidal overlap Λ_f = (1 − N_g·g_frac)·Λ_static; each g-mode
carries the (small, literature-order) overlap fraction `g_frac`·Λ_static
(Lai 1994 / Reisenegger 1994: g-modes couple ≲1% of the f-mode overlap).  The
sum rule Σ_n Λ_n = Λ_static holds by construction.
"""
function lambda_eff_modes(Λ_static::Float64, f_f_kHz::Float64;
                          g_freqs_kHz::AbstractVector{<:Real}=Float64[],
                          g_frac::Float64=1e-3)
    mode_f = Float64[f_f_kHz]
    g_total = g_frac*length(g_freqs_kHz)
    mode_Λ = Float64[(1 - g_total)*Λ_static]
    for gf in g_freqs_kHz
        push!(mode_f, float(gf))
        push!(mode_Λ, g_frac*Λ_static)
    end
    mode_ω = [2π*f*kHz_to_km for f in mode_f]
    return LambdaEffModes(mode_f, mode_Λ, mode_ω)
end

"""
    lambda_eff(modes, ω) -> Λ_eff

Driven-oscillator effective tidal deformability
`Λ_eff(ω) = Σ_n Λ_n ω_n²/(ω_n² − ω²)` at driving frequency `ω` [km⁻¹].
`lambda_eff(modes; f_GW_kHz=...)` evaluates at the GW frequency f_GW=ω/π.
"""
lambda_eff(m::LambdaEffModes, ω::Real) =
    sum(m.mode_Λ[n]*m.mode_ω[n]^2/(m.mode_ω[n]^2 - ω^2) for n in 1:length(m.mode_f_kHz))

lambda_eff(m::LambdaEffModes; f_GW_kHz::Real) =
    lambda_eff(m, π*f_GW_kHz*kHz_to_km)

# convenience: the f-mode resonance GW frequency f_GW = ω_f/π = 2 f_f
fmode_resonance_fGW_kHz(m::LambdaEffModes) = 2*m.mode_f_kHz[1]

# ── Yagi–Yunes (2013) I–Love fit:  ln Ī = Σ coeff (ln Λ)^k  (Table I, Ī–Λ̄) ──
const YY_ILOVE_COEFFS = (1.47, 0.0817, 0.0149, 2.87e-4, -3.64e-5)

"Yagi–Yunes (2013, Science 341,365) universal I–Love fit: Ī(Λ)."
function ibar_yagi_yunes(Λ::Real)
    x = log(Λ)
    a,b,c,d,e = YY_ILOVE_COEFFS
    return exp(a + b*x + c*x^2 + d*x^3 + e*x^4)
end

end # module TidalDeformability
