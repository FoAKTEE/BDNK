#=
    SlowRotation — Hartle (1967) slow-rotation frame dragging + moment of inertia
    on the STATIC TOV background.

    PRIMARY: J. B. Hartle, "Slowly Rotating Relativistic Stars. I," ApJ 150, 1005
    (1967): metric eq. (30), j(r)=exp[-(ν+λ)/2] eq. (40), frame-drag ODE
    (divergence form) eqs. (43)/(46), center regularity eq. (44), exterior
    ϖ=Ω-2J/r³ eq. (47), J=IΩ eq. (48).
    CROSS-CHECK: Fattoyev & Piekarewicz, Phys. Rev. C 82, 025810 (2010)
    [arXiv:1006.3758]: I integral eq. (2), ODE eq. (7), j eq. (8), BCs eq. (9),
    consistency ϖ̃'(R)=6GI/R⁴ eq. (10).

    ── METRIC CONVENTION (read from src/tov/TOV.jl, NOT assumed) ───────────────
    TOVStar.ν is HARTLE's RAW exponent: g_tt = -e^{ν}, fixed by the surface match
    `νsurf = log(1 - 2M/R)` in TOV.jl (an e^{2Φ} convention would use 0.5·log).
    Hence the physical metric potential is Φ = ν/2, and

        j(r) = e^{-(ν+λ)/2} = e^{-ν/2} · √(1 - 2m/r)        (Hartle eq. 40)

    with e^{-λ} = 1 - 2m/r. In the canonical Φ-form j = e^{-Φ}√(1-2m/r), identical
    since Φ = ν/2. We therefore use star.ν directly with a factor 1/2.

    ── CANONICAL ODE (implemented as two first-order ODEs, sign-safe) ──────────
        d/dr[ r⁴ j(r) dϖ/dr ] + 4 r³ (dj/dr) ϖ = 0,     ϖ ≡ Ω - ω.
    Let g(r) = r⁴ j(r) ϖ'(r). Then
        dg/dr  = -4 r³ (dj/dr) ϖ,
        dϖ/dr  = g / (r⁴ j).
    Center BCs (regular branch): ϖ(0)=ϖc (arbitrary, =1), ϖ'(0)=0 ⇒ g(0)=0.

    ── SURFACE / I ────────────────────────────────────────────────────────────
        J = (R⁴/6) ϖ'(R),   Ω = ϖ(R) + (R/3) ϖ'(R),   I = J/Ω.
    Linear+homogeneous ⇒ I independent of ϖc and of Ω. Central drag fraction
    ω(0)/Ω = 1 - ϖ(0)/Ω.

    ── UNITS ──────────────────────────────────────────────────────────────────
    I_geo in km³ (G=c=1). cgs: I[g cm²] = I[km³]·(mass-per-km in g)·(cm/km)².
    The mass-per-km factor is derived from the Units module (Units.kg_to_m =
    G/c²), NOT hard-coded; numerically ≈1.3466×10⁴³ g cm² per km³. Ī = I/M³ (M km).
=#
module SlowRotation

using ..EquationOfState: BarotropicEOS
using ..TOV: TOVStar, solve_tov, mass_solar
using ..Units

export SlowRotResult, moment_of_inertia, IBAR_FROM_C_BREU, MR2_FROM_C_LS,
       KM3_TO_G_CM2

# km³ of moment of inertia → g cm². 1 km of geometric mass = 1e6/kg_to_m grams;
# (1 km)² = 1e10 cm². Both derived from Units (no hard-coded c²/G).
const KM3_TO_G_CM2 = (1e6 / Units.kg_to_m) * 1e10   # ≈ 1.3466e43 g cm² / km³

# ── Published universal relations (CITED; see test/viz headers) ──────────────
"""
    IBAR_FROM_C_BREU(C) -> Ī = I/M³

Breu & Rezzolla (2016), MNRAS 459, 646 [arXiv:1601.06083], Eq. (20), Table 2,
SLOW-ROTATION row:  Ī = ā₁C⁻¹ + ā₂C⁻² + ā₃C⁻³ + ā₄C⁻⁴ with
ā₁=8.134e-1, ā₂=2.101e-1, ā₃=3.175e-3, ā₄=-2.717e-4 (χ²red=0.1184; ⟨L₁⟩≈1.1%).
"""
IBAR_FROM_C_BREU(C) = 0.8134/C + 0.2101/C^2 + 3.175e-3/C^3 - 2.717e-4/C^4

"""
    MR2_FROM_C_LS(C) -> Ĩ = I/(M R²)

Lattimer & Schutz (2005), ApJ 629, 979: Ĩ = ã₀ + ã₁C + ã₄C⁴ with
ã₀=0.237, ã₁=0.674, ã₄=4.48 (equivalently 0.237(1+2.844C+18.91C⁴)).
(These compactness-form coefficients are as quoted by Breu & Rezzolla 2016;
the LS original uses an equivalent differently-normalized dimensional form.)
"""
MR2_FROM_C_LS(C) = 0.237 + 0.674*C + 4.48*C^4

struct SlowRotResult
    M_Msun::Float64
    M_km::Float64
    R_km::Float64
    C::Float64               # compactness M/R
    I_geo_km3::Float64       # moment of inertia, geometric (km³)
    I_cgs::Float64           # g cm²
    Ibar::Float64            # I/M³ (dimensionless)
    I_MR2::Float64           # I/(M R²) (dimensionless)
    J::Float64               # angular momentum (geometric, ϖc-normalized cancels in I)
    Ω::Float64               # spin (geometric, ϖc-normalized)
    drag_ratio::Float64      # ω(0)/Ω central frame-drag fraction
    r::Vector{Float64}       # radial grid [km]
    ϖ::Vector{Float64}       # ϖ(r) = Ω-ω(r) (normalized ϖ(0)=ϖc)
    ω_over_Ω::Vector{Float64}# local dragging ω(r)/Ω
end

# j(r) and dj/dr on the TOV grid.  j = e^{-ν/2}√(1-2m/r).
# dj/dr computed analytically from TOV: ν' = 2(m+4πr³p)/(r(r-2m)),
# m' = 4πr²ε, so d/dr[√(1-2m/r)] = -(m'/r - m/r²)/√(1-2m/r) handled below.
function _j_and_dj(star::TOVStar)
    n = length(star.r)
    j  = Vector{Float64}(undef, n)
    dj = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        r = star.r[i]; m = star.m[i]; p = star.p[i]; ε = star.ε[i]; ν = star.ν[i]
        f   = 1 - 2m/r                      # = e^{-λ}
        sf  = sqrt(max(f, 0.0))
        eΦm = exp(-ν/2)                     # e^{-ν/2} = e^{-Φ}
        j[i] = eΦm * sf
        # ν' (Hartle raw exponent); Φ' = ν'/2
        denom = r*(r - 2m)
        νp = 2*(m + 4π*r^3*p)/denom
        mp = 4π*r^2*ε
        # d/dr √f = (1/(2√f)) d f/dr, df/dr = -2(m'/r - m/r²) = -2(mp*r - m)/r²
        dfdr = -2*(mp*r - m)/r^2
        dsf  = sf > 0 ? dfdr/(2*sf) : 0.0
        dj[i] = -0.5*νp*eΦm*sf + eΦm*dsf    # d/dr[e^{-ν/2}]·sf + e^{-ν/2}·dsf
    end
    return j, dj
end

"""
    moment_of_inertia(eos, εc; ϖc=1.0, kwargs...) -> SlowRotResult

Solve Hartle's slow-rotation frame-dragging ODE on the static TOV background with
central energy density `εc` (km⁻²) and return the moment of inertia and dragging
profile. Extra keywords (`h`, `rmax`, `ptol_rel`) are forwarded to `solve_tov`.

The result is independent of `ϖc` (linear-homogeneous ODE); the keyword exists
only to make the normalization-independence testable.
"""
function moment_of_inertia(eos::BarotropicEOS, εc::Float64; ϖc::Float64=1.0,
                           h::Float64=0.005, rmax::Float64=60.0,
                           ptol_rel::Float64=1e-10)
    star = solve_tov(eos, εc; h=h, rmax=rmax, ptol_rel=ptol_rel)
    return moment_of_inertia(star; ϖc=ϖc)
end

"""
    moment_of_inertia(star::TOVStar; ϖc=1.0) -> SlowRotResult

Slow-rotation moment of inertia from a precomputed TOV background.
"""
function moment_of_inertia(star::TOVStar; ϖc::Float64=1.0)
    r = star.r
    n = length(r)
    j, dj = _j_and_dj(star)

    ϖ = Vector{Float64}(undef, n)
    g = Vector{Float64}(undef, n)   # g = r⁴ j ϖ'
    ϖ[1] = ϖc
    g[1] = 0.0                       # ϖ'(0)=0 ⇒ g(0)=0 (regular branch)

    # RHS of the two first-order ODEs at (r, ϖ, g)
    @inline function rhs(rr, ϖϖ, gg, jr, djr)
        dϖ = gg / (rr^4 * jr)
        dg = -4 * rr^3 * djr * ϖϖ
        return dϖ, dg
    end

    # RK4 with midpoint interpolation of (j, dj) on the nonuniform TOV grid.
    @inbounds for i in 1:n-1
        hh = r[i+1] - r[i]
        jm  = 0.5*(j[i] + j[i+1])
        djm = 0.5*(dj[i] + dj[i+1])
        rm  = r[i] + hh/2

        k1ϖ, k1g = rhs(r[i],   ϖ[i],          g[i],          j[i],  dj[i])
        k2ϖ, k2g = rhs(rm,     ϖ[i]+hh/2*k1ϖ, g[i]+hh/2*k1g, jm,    djm)
        k3ϖ, k3g = rhs(rm,     ϖ[i]+hh/2*k2ϖ, g[i]+hh/2*k2g, jm,    djm)
        k4ϖ, k4g = rhs(r[i+1], ϖ[i]+hh*k3ϖ,   g[i]+hh*k3g,   j[i+1],dj[i+1])

        ϖ[i+1] = ϖ[i] + hh/6*(k1ϖ + 2k2ϖ + 2k3ϖ + k4ϖ)
        g[i+1] = g[i] + hh/6*(k1g + 2k2g + 2k3g + k4g)
    end

    R = star.R; M = star.M
    # ϖ'(R) from g(R) = R⁴ j(R) ϖ'(R).  At the surface j(R)=e^{-ν(R)/2}√(1-2M/R)=1:
    # the TOV ν-shift sets e^{-ν(R)/2}=1/√(1-2M/R), which cancels the √(1-2M/R) factor
    # (Schwarzschild-exterior value).  DO NOT "simplify" to √(1-2M/R) — that is a 20%
    # error.  j[n] already holds the correct value (≈1); we divide by it explicitly.
    ϖpR = g[n] / (R^4 * j[n])
    J = (R^4 / 6) * ϖpR
    Ω = ϖ[n] + (R/3) * ϖpR
    I_geo = J / Ω

    drag = 1 - ϖ[1]/Ω                       # ω(0)/Ω = 1 - ϖ(0)/Ω
    ω_over_Ω = 1 .- ϖ ./ Ω

    C = M / R
    return SlowRotResult(
        mass_solar(star), M, R, C,
        I_geo, I_geo*KM3_TO_G_CM2, I_geo/M^3, I_geo/(M*R^2),
        J, Ω, drag, copy(r), ϖ, ω_over_Ω)
end

end # module SlowRotation
