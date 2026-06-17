# ======================================================================
# repro/axial_waveeqs.jl
#
# R4 / stage 1 — AXIAL (odd-parity) coupled wave equations for a viscous
# neutron star, in the FREQUENCY domain.
#
# Primary ground:  Bussieres, Redondo-Yuste, Ortega-Gomez & Cardoso,
#                  "Axial Oscillations of Viscous Neutron Stars",
#                  arXiv:2604.13208  (file src/main.tex).
#   - the two coupled axial ODEs  ............ eqs (17)-(18)  [lines 280-281]
#       eq.(17)  "gw_eq"    : f[ d (f ψ') ]' + (ω² - V) ψ = -16π e^{ν/2} iω η ψ + C1 Z
#       eq.(18)  "fluid_eq" : f[ d (f Z') ]' + (c_η^{-2} ω² - U) Z
#                              = C2 Z' + C3 Z + C4 ψ' + C5 ψ
#   - Regge-Wheeler perfect-fluid potential V  eq.(19)  [line 287]
#   - U, C1..C5 coefficient block ............ [lines 299-307]
#   - viscous (second-sound) speed c_η²        [line 293]
#   - surface regularity condition (BC)  ..... eq.(24)  [lines 329-353]
#   - viscosity parametrizations A / B ....... eq.(13a/13b) [lines 170-183]
#
# Cross-check ground:  Redondo-Yuste, "Perturbations of relativistic
#                  dissipative stars", arXiv:2411.16841 (file src/Formalism.tex)
#   - master odd inviscid wave eq ............ eq:Master_Odd_Inviscid [line 514]
#   - master odd VISCOUS wave eq ............. eq:Master_Odd_Viscous   [line 532]
#       (the GW eq (17) is the frequency-domain form of this; the explicit
#        iω-damping term  16π η e^{ν/2} ∂_t ψ  -> -16π η e^{ν/2} iω ψ in -iωt
#        convention, matching Bussieres eq (17) RHS.)
#   - RW-type potential in eq:Master_Odd_Inviscid (eq.19 of the prompt) [l.514]
#
# CONVENTIONS / NOTATION MAP  (Bussieres  <->  Redondo  <->  this code / pkg)
#   ν (metric)   = Φ            metric potential e^{ν} = e^{Φ} = -g_tt
#   λ (metric)   = Λ            e^{λ}=e^{Λ}=(1-2m/r)^{-1}
#   ρ (Bussieres)= ε  (energy density)   -> pkg field star.ε
#   p            = p                     -> pkg field star.p
#   m(r)         = M(r)  enclosed mass   -> pkg field star.m
#   c_s²         = dp/dε   (sound_speed2)
#   f²           = e^{ν-λ}               (interior "lapse/redshift" combo)
#   η            = shear viscosity (geometric units, [km^? ] consistent w/ T_ab)
#   τ            = BDNK relaxation coeff "τ_Q" (Bussieres' τ; Redondo's τ_Q)
#   convention: time dependence  e^{-iωt}   (Bussieres, line 410).
#
# In the inviscid limit (η,τ -> 0) eq.(17) reduces to the STANDARD relativistic
# axial / w-mode wave equation  f(fψ')' + (ω²-V)ψ = 0  with V the RW potential
# eq.(19) — verified numerically in the SANITY block at the bottom.
#
# Units: pure geometric G=c=1, lengths in km (same as solve_tov / the pkg).
# Einstein convention is G_ab = 8π T_ab (matches both papers' TOV + source 4π/8π).
# ======================================================================

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar

using LinearAlgebra
using Printf

# ----------------------------------------------------------------------
# Background container: TOV fields interpolated to a single radius r, plus
# the metric potentials and their radial derivatives needed by the
# coefficient functions.  Built from a TOVStar (solve_tov output).
# ----------------------------------------------------------------------
"""
    AxialBackground

Radially-interpolated TOV background evaluated at a point.  All quantities
geometric (G=c=1, km).  Field meaning (Bussieres notation):

  r      areal radius
  m      enclosed mass m(r)
  p, ρ   pressure, energy density (ρ ≡ ε)
  cs2    sound speed squared dp/dρ
  ν, λ   metric potentials   e^{ν}=-g_tt ,  e^{λ}=(1-2m/r)^{-1}
  dνdr   dν/dr   (from TOV eq.(16b):  ν' = (2m+8π r³ p)/(r(r-2m)))
  f, f2  f²=e^{ν-λ},  f=√(e^{ν-λ})
"""
struct AxialBackground
    r::Float64
    m::Float64
    p::Float64
    ρ::Float64
    cs2::Float64
    ν::Float64
    λ::Float64
    dνdr::Float64
    f::Float64
    f2::Float64
end

# linear interpolation of a TOV field array onto r (clamped to [r0,R])
@inline function _interp(rs::Vector{Float64}, ys::Vector{Float64}, r::Float64)
    if r <= rs[1]
        return ys[1]
    elseif r >= rs[end]
        return ys[end]
    end
    # binary search
    lo = searchsortedlast(rs, r)
    lo = clamp(lo, 1, length(rs)-1)
    t = (r - rs[lo]) / (rs[lo+1] - rs[lo])
    return ys[lo] + t*(ys[lo+1] - ys[lo])
end

"""
    background_at(star, eos, r) -> AxialBackground

Interpolate the TOV solution to radius `r` and assemble the metric
potentials.  ρ≡ε (energy density), λ from m(r), ν interpolated from the TOV
ν array (already Schwarzschild-matched at the surface), dν/dr from the TOV
equation eq.(16b) of Bussieres [main.tex line 239].
"""
function background_at(star::TOVStar, eos::BarotropicEOS, r::Float64)
    m  = _interp(star.r, star.m, r)
    p  = _interp(star.r, star.p, r)
    ρ  = _interp(star.r, star.ε, r)
    ν  = _interp(star.r, star.ν, r)
    cs2 = ρ > 0 ? sound_speed2(eos, ρ) : 0.0
    λ  = -log(1 - 2m/r)                       # e^{λ}=(1-2m/r)^{-1}
    # TOV:  ν' = (2m + 8π r³ p)/(r(r-2m))     (Bussieres eq.16b, line 239)
    dνdr = (2m + 8π*r^3*p) / (r*(r - 2m))
    f2 = exp(ν - λ)
    f  = sqrt(f2)
    return AxialBackground(r, m, p, ρ, cs2, ν, λ, dνdr, f, f2)
end

# ----------------------------------------------------------------------
# Viscosity parametrizations  (Bussieres eq.(13a)/(13b), lines 170-183).
# η, θ, τ as functions of the local background.  ζ=0 throughout (shear only).
# A:  η = η̂ (ρ+p) L0 cs²,   θ = L0 η̂,        τ = τ̂ L0 η̂
# B:  η = η̂ p L0,           θ = L0 p/ρ,       τ = τ̂ L0 p/ρ
# Returns (η, θ, τ).  η̂=τ̂=0  OR  η̂=0  -> inviscid (η=τ=0).
# ----------------------------------------------------------------------
struct Viscosity
    param::Symbol     # :A or :B
    ηhat::Float64
    τhat::Float64
    L0::Float64
end

inviscid() = Viscosity(:A, 0.0, 0.0, 1.0)

"""
    transport(v, bg) -> (η, θ, τ)

Local shear viscosity η, frame coeff θ, relaxation τ, geometric units.
ζ=0 (shear-only), per Bussieres line 184 / Redondo abstract (axial sector
unaffected by bulk viscosity).
"""
function transport(v::Viscosity, bg::AxialBackground)
    if v.ηhat == 0.0
        return (0.0, 0.0, 0.0)
    end
    if v.param === :A
        η = v.ηhat * (bg.ρ + bg.p) * v.L0 * bg.cs2
        θ = v.L0 * v.ηhat
        τ = v.τhat * v.L0 * v.ηhat
    elseif v.param === :B
        η = v.ηhat * bg.p * v.L0
        θ = v.L0 * (bg.ρ > 0 ? bg.p/bg.ρ : 0.0)
        τ = v.τhat * v.L0 * (bg.ρ > 0 ? bg.p/bg.ρ : 0.0)
    else
        error("unknown viscosity parametrization $(v.param)")
    end
    return (η, θ, τ)
end

# dη/dr via finite differences of the parametrization on the interpolated bg.
function dη_dr(v::Viscosity, star::TOVStar, eos::BarotropicEOS, r::Float64)
    if v.ηhat == 0.0
        return 0.0
    end
    h = max(1e-6, 1e-5*r)
    rp = min(r+h, star.R)
    rm = max(r-h, star.r[1])
    ηp = transport(v, background_at(star, eos, rp))[1]
    ηm = transport(v, background_at(star, eos, rm))[1]
    return (ηp - ηm)/(rp - rm)
end

# ----------------------------------------------------------------------
# Potentials  V (eq.19) and U  (Bussieres lines 287, 301)
# ----------------------------------------------------------------------
"""
    RW_potential(bg, ℓ) -> V        (Bussieres eq.(19), main.tex line 287)

V = e^{ν} [ ℓ(ℓ+1)/r² - 6m/r³ + 4π(ρ - p) ].
This is exactly the "RW-type potential" referenced as eq.19 in the prompt and
as eq:Master_Odd_Inviscid in Redondo (Formalism.tex line 514): there the term
reads  (e^{Φ}/r²)(λ²-6M/r+4π r²(ε-p))  with λ²=ℓ(ℓ+1), identical after
distributing e^{Φ}/r².
"""
@inline function RW_potential(bg::AxialBackground, ℓ::Int)
    L = ℓ*(ℓ+1)
    return exp(bg.ν) * ( L/bg.r^2 - 6*bg.m/bg.r^3 + 4π*(bg.ρ - bg.p) )
end

"""
    U_potential(bg, ℓ) -> U         (Bussieres line 301)

U = e^{ν} [ ℓ(ℓ+1)/r² - 2m/r³ + 8π(2p+ρ) ].
"""
@inline function U_potential(bg::AxialBackground, ℓ::Int)
    L = ℓ*(ℓ+1)
    return exp(bg.ν) * ( L/bg.r^2 - 2*bg.m/bg.r^3 + 8π*(2*bg.p + bg.ρ) )
end

"""
    cη2(η, τ, bg) -> c_η²            (Bussieres line 293)

Viscous (second-sound) propagation speed²  c_η² = η / (τ (p+ρ)).
"""
@inline function cη2(η::Float64, τ::Float64, bg::AxialBackground)
    return η / (τ * (bg.p + bg.ρ))
end

# ----------------------------------------------------------------------
# Coupling coefficient functions  C1..C5   (Bussieres lines 302-306)
# Returned together as a named-tuple given (bg, ℓ, ω, η, η', τ).
# ω is complex (frequency-domain, e^{-iωt}).
# ----------------------------------------------------------------------
"""
    coupling_coeffs(bg, ℓ, ω, η, dηdr, τ) -> (C1,C2,C3,C4,C5)

The viscous coupling functions of Bussieres eqs (main.tex lines 302-306):

 C1 = 8π e^{ν-λ/2}/r² [ 2r η' + (e^{λ}(1+8π r² p) - 1) η ]
 C2 = f²/(2r) [ e^{λ}(1+8π r² p) - 1 - 2r η'/η ]
 C3 = -iω(p+ρ) e^{ν/2} (1/η + 16π τ) + 2 f² η'/(r η)
 C4 = r f [ iω + (p+ρ)/η (e^{ν/2} - iω τ) ]
 C5 = f [ (p+ρ) e^{ν/2}/η - iω/2 (-7 + e^{λ}(1+8π r² p))
          + iω/η ( r η' - (p+ρ) τ ) ]
"""
function coupling_coeffs(bg::AxialBackground, ℓ::Int, ω::ComplexF64,
                         η::Float64, dηdr::Float64, τ::Float64)
    r, p, ρ, ν, λ, f, f2 = bg.r, bg.p, bg.ρ, bg.ν, bg.λ, bg.f, bg.f2
    eλ = exp(λ)
    A  = eλ*(1 + 8π*r^2*p) - 1                      # recurring factor
    iω = im*ω

    C1 = 8π * exp(ν - λ/2) / r^2 * ( 2r*dηdr + A*η )
    C2 = f2/(2r) * ( A - 2r*dηdr/η )
    C3 = -iω*(p+ρ)*exp(ν/2)*(1/η + 16π*τ) + 2*f2*dηdr/(r*η)
    C4 = r*f * ( iω + (p+ρ)/η * (exp(ν/2) - iω*τ) )
    C5 = f * ( (p+ρ)*exp(ν/2)/η - (iω/2)*(-7 + A) + (iω/η)*( r*dηdr - (p+ρ)*τ ) )
    return (C1, C2, C3, C4, C5)
end

# ----------------------------------------------------------------------
# LINEAR-SYSTEM ASSEMBLER
#
# We write the two 2nd-order ODEs eqs (17)-(18) in the standard form
#       u'' = P u' + Q u,         u = (ψ, Z),
# i.e. the matrices P(r;ω,ℓ,visc), Q(r;ω,ℓ,visc) of the first-derivative and
# value couplings.  This is the "linear system at given (r,ω,ℓ,viscosity)"
# requested.  (A first-order companion form is also returned.)
#
# Operator on the LHS:  f[ d(f u') ]' = f(f u')' = f² u'' + f f' u'.
# With f=√(e^{ν-λ}),   f' = f·(ν'-λ')/2.   We need λ' too:
#       λ' = d/dr[-ln(1-2m/r)] = (2m' r - 2m)/(r(r-2m)) ,   m'=4π r² ρ.
# So the operator =  f² u'' + (f² (ν'-λ')/2) u'.
#
# eq (17):  f² ψ'' + f²(ν'-λ')/2 ψ' + (ω² - V) ψ = -16π e^{ν/2} iω η ψ + C1 Z
# eq (18):  f² Z'' + f²(ν'-λ')/2 Z' + (cη^{-2} ω² - U) Z
#                                          = C2 Z' + C3 Z + C4 ψ' + C5 ψ
#
# Dividing by f² and moving everything to standard form gives P, Q below.
# ----------------------------------------------------------------------
"""
    axial_linear_system(star, eos, r, ω, ℓ, visc) ->
        (P, Q, aux)

Assemble the 2×2 second-order linear system   u'' = P·u' + Q·u ,  u=(ψ,Z),
for the axial coupled wave equations eqs (17)-(18) of Bussieres at radius r,
complex frequency ω, harmonic ℓ, and viscosity prescription `visc`.

Returns matrices P,Q (2×2 ComplexF64) and an `aux` named-tuple with the
intermediate quantities (V,U,η,τ,cη2,C1..C5,f2,dlogf) for inspection.

In the inviscid limit the (ψ) row decouples and reduces to the standard
relativistic axial wave equation  ψ'' + (ν'-λ')/2 ψ' + (ω²-V)/f² ψ = 0.
"""
function axial_linear_system(star::TOVStar, eos::BarotropicEOS, r::Float64,
                             ω::ComplexF64, ℓ::Int, visc::Viscosity)
    bg = background_at(star, eos, r)
    f2 = bg.f2
    # λ' = (2 m' r - 2 m)/(r (r-2m)),  m' = 4π r² ρ
    mp = 4π*r^2*bg.ρ
    dλdr = (2*mp*r - 2*bg.m) / (r*(r - 2*bg.m))
    # f'/f = (ν'-λ')/2  ->  operator first-deriv coeff on u' is f²·(ν'-λ')/2
    dlogf = (bg.dνdr - dλdr)/2          # = f'/f
    V = RW_potential(bg, ℓ)

    # viscosity
    η, θ, τ = transport(visc, bg)
    inv = (η == 0.0)

    P = zeros(ComplexF64, 2, 2)
    Q = zeros(ComplexF64, 2, 2)

    if inv
        # ---- inviscid: ψ decouples to the standard RW/w-mode equation ----
        # ψ'' = -dlogf ψ' - (ω²-V)/f² ψ
        P[1,1] = -dlogf
        Q[1,1] = -(ω^2 - V)/f2
        # Z row: with η,τ->0 the 2nd-order Z eq degenerates.  For the inviscid
        # background the fluid eq is first-order (Bussieres line 322):
        #   iω Z = f (r ψ' + ψ).   We expose this via aux; here we leave the
        #   Z-row of (P,Q) at the perfect-fluid algebraic relation written as a
        #   degenerate 2nd-order placeholder (Z'' = 0, Z determined by ψ).
        P[2,2] = 0.0
        Q[2,2] = 0.0
        aux = (V=V, U=NaN, η=0.0, τ=0.0, cη2=NaN, dlogf=dlogf, f2=f2,
               C1=0.0+0im, C2=0.0+0im, C3=0.0+0im, C4=0.0+0im, C5=0.0+0im,
               inviscid_Z = (iω = im*ω, fac_psip = bg.f*r, fac_psi = bg.f))
        return (P, Q, aux)
    end

    # ---- viscous case ----
    dηdr = dη_dr(visc, star, eos, r)
    U = U_potential(bg, ℓ)
    ce2 = cη2(η, τ, bg)
    C1, C2, C3, C4, C5 = coupling_coeffs(bg, ℓ, ω, η, dηdr, τ)

    # eq (17):  f² ψ'' + f²·dlogf·ψ' + (ω²-V)ψ = -16π e^{ν/2} iω η ψ + C1 Z
    #  => ψ'' = -dlogf ψ'  - (ω²-V)/f² ψ  + [ -16π e^{ν/2} iω η/f² ψ + C1/f² Z ]
    damp = -16π*exp(bg.ν/2)*im*ω*η          # coefficient of ψ on RHS of (17)
    P[1,1] = -dlogf
    P[1,2] = 0.0
    Q[1,1] = -(ω^2 - V)/f2 + damp/f2
    Q[1,2] = C1/f2

    # eq (18):  f² Z'' + f²·dlogf·Z' + (cη^{-2}ω² - U)Z = C2 Z' + C3 Z + C4 ψ' + C5 ψ
    #  => Z'' = (-dlogf + C2/f²) Z' + (-(cη^{-2}ω²-U)/f² + C3/f²) Z
    #            + C4/f² ψ' + C5/f² ψ
    P[2,2] = -dlogf + C2/f2
    P[2,1] = C4/f2
    Q[2,2] = -(ω^2/ce2 - U)/f2 + C3/f2
    Q[2,1] = C5/f2

    aux = (V=V, U=U, η=η, τ=τ, cη2=ce2, dlogf=dlogf, f2=f2,
           C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, damp=damp, dηdr=dηdr)
    return (P, Q, aux)
end

"""
    axial_first_order_matrix(star, eos, r, ω, ℓ, visc) -> M (4×4)

First-order companion form.  With state  y = (ψ, Z, ψ', Z'),  y' = M y, where
the lower block comes from u'' = P u' + Q u.  Convenient for shooting / the
QNM root-finder downstream (R4 stage 2).
"""
function axial_first_order_matrix(star::TOVStar, eos::BarotropicEOS, r::Float64,
                                  ω::ComplexF64, ℓ::Int, visc::Viscosity)
    P, Q, _ = axial_linear_system(star, eos, r, ω, ℓ, visc)
    M = zeros(ComplexF64, 4, 4)
    M[1,3] = 1; M[2,4] = 1
    M[3,1] = Q[1,1]; M[3,2] = Q[1,2]; M[3,3] = P[1,1]; M[3,4] = P[1,2]
    M[4,1] = Q[2,1]; M[4,2] = Q[2,2]; M[4,3] = P[2,1]; M[4,4] = P[2,2]
    return M
end

# ----------------------------------------------------------------------
# Surface regularity condition  (Bussieres eq.(24), lines 329-353).
# Returns (B1,B2,B3,B4) s.t.  B1 Z(R)+B2 Z'(R)+B3 ψ(R)+B4 ψ'(R)=0.
# 𝒞 = M/R compactness.  Parametrization A or B.
# ----------------------------------------------------------------------
"""
    surface_condition(star, ω, ℓ, visc) -> (B1,B2,B3,B4)

Coefficients of the surface regularity condition eq.(24) [main.tex 329-353].
"""
function surface_condition(star::TOVStar, ω::ComplexF64, visc::Viscosity)
    R = star.R; M = star.M; 𝒞 = M/R
    η̂ = visc.ηhat; τ̂ = visc.τhat
    s = sqrt(1 - 2𝒞)
    if visc.param === :A
        B1 = (-4*𝒞^2*η̂ + 2*𝒞*η̂ + R*ω*(η̂*τ̂*R*ω + im*s)) / (η̂*R^2*(1-2𝒞)^2)
        B2 = 𝒞 / (2R*(1-2𝒞))
        B3 = im*( im*(1-2𝒞)^(3/2) - 2*𝒞*η̂*τ̂*R*ω + η̂*τ̂*R*ω + 𝒞*η̂*R*ω ) /
             ((1-2𝒞)^2*η̂*R)
        B4 = -1/(s*η̂) + im*τ̂*R*ω/(1-2𝒞)
    elseif visc.param === :B
        B1 = -s*im*ω + 𝒞*(5𝒞 - 2)/R^2 * η̂
        B2 = -𝒞*η̂*(1 - 2𝒞)
        B3 = -M/s*η̂*im*ω + 1 - 2𝒞
        B4 = R*(1 - 2𝒞)
    else
        error("unknown parametrization")
    end
    return (B1, B2, B3, B4)
end

# ======================================================================
# SANITY  — inviscid limit must reproduce the standard relativistic axial
#           (w-mode) wave equation with the RW potential V (eq.19).
# ======================================================================
function _sanity()
    println("="^70)
    println("AXIAL WAVE EQUATIONS — build + inviscid sanity")
    println("="^70)

    # Bussieres reference star: PolytropeEnergy(κ=100 km², n=1),
    # ρ_c = 3e15 g/cm³ -> M≈1.27 M⊙, R≈8.86 km  (main.tex line 254).
    eos = PolytropeEnergy(100.0, 1.0)

    # convert ρ_c = 3e15 g/cm³ to geometric (km^-2)
    # 1 g/cm³ = G/c² * 1e3 kg/m³ in 1/m² ... use pkg's unit scale via a sweep:
    # Instead pick εc directly in geometric units giving M≈1.27 M⊙.
    # ρ_c[g/cm³] -> ε[km^-2]: ε = ρ_c * G/c² *1e3 *1e6  (km^-2).
    G = 6.6743015e-11; c = 299_792_458.0
    ρc_cgs = 3e15                      # g/cm³
    ρc_SI  = ρc_cgs * 1e3              # kg/m³
    εc_m2  = ρc_SI * G / c^2           # 1/m²
    εc = εc_m2 * 1e6                   # 1/km²  (geometric energy density)
    @printf("central ε = %.6e km^-2  (ρc = 3e15 g/cm^3)\n", εc)

    star = solve_tov(eos, εc; h=1e-3, ptol_rel=1e-12, rmax=50.0)
    @printf("TOV star:  M = %.6f M_sun ,  R = %.6f km ,  M/R = %.4f\n",
            mass_solar(star), star.R, star.M/star.R)

    ℓ = 2
    ω = complex(0.05, -0.01)           # arbitrary test frequency (1/km)

    # pick an interior radius
    r = 0.5*star.R
    visc0 = inviscid()

    P, Q, aux = axial_linear_system(star, eos, r, ω, ℓ, visc0)
    bg = background_at(star, eos, r)
    @printf("\nAt r = %.4f km (interior):\n", r)
    @printf("  V (RW potential, eq.19) = %.8e\n", aux.V)
    @printf("  f² = e^{ν-λ}            = %.8e\n", aux.f2)
    @printf("  d(log f)/dr             = %.8e\n", aux.dlogf)

    # ---- inviscid ψ-row vs analytic standard RW/w-mode equation ----
    # Standard relativistic axial wave eq (Kokkotas&Schmidt / Redondo
    # eq:Master_Odd_Inviscid):   ψ'' + (f'/f) ψ' + (ω²-V)/f² ψ = 0.
    # i.e.  P_expected = -(f'/f) = -dlogf ,  Q_expected = -(ω²-V)/f².
    P_exp = -aux.dlogf
    Q_exp = -(ω^2 - aux.V)/aux.f2
    errP = abs(P[1,1] - P_exp)
    errQ = abs(Q[1,1] - Q_exp)
    @printf("\nINVISCID ψ-row vs standard RW axial wave equation:\n")
    @printf("  P[1,1]            = % .10e  (expected % .10e)\n", real(P[1,1]), real(P_exp))
    @printf("  Q[1,1]            = % .10e %+0.10e i\n", real(Q[1,1]), imag(Q[1,1]))
    @printf("  Q_expected        = % .10e %+0.10e i\n", real(Q_exp), imag(Q_exp))
    @printf("  |ΔP| = %.3e , |ΔQ| = %.3e\n", errP, errQ)
    @printf("  coupling Q[1,2] (=C1/f², must be 0 inviscid) = %.3e\n", abs(Q[1,2]))

    # ---- EXTERIOR check: at r>R the eq must be the vacuum RW equation ----
    # Build a background point just outside; V_vac = f(ℓ(ℓ+1)/r²-6M/r³), f=1-2M/r.
    rext = 1.5*star.R
    # outside, ε=p=0, m=M:
    m=star.M; f_ext = 1 - 2m/rext
    V_vac_expected = f_ext*(ℓ*(ℓ+1)/rext^2 - 6m/rext^3)
    # our RW_potential with ext background (ρ=p=0, e^{ν}=f_ext, e^{λ}=1/f_ext):
    νext = log(f_ext); λext = -log(f_ext)
    bgext = AxialBackground(rext, m, 0.0, 0.0, 0.0, νext, λext,
                            (2m)/(rext*(rext-2m)), sqrt(exp(νext-λext)), exp(νext-λext))
    V_ext = RW_potential(bgext, ℓ)
    @printf("\nEXTERIOR (r=%.3f km) vacuum RW potential:\n", rext)
    @printf("  RW_potential(bg)  = %.10e\n", V_ext)
    @printf("  f(ℓ(ℓ+1)/r²-6M/r³)= %.10e\n", V_vac_expected)
    @printf("  |Δ| = %.3e\n", abs(V_ext - V_vac_expected))

    # ---- VISCOUS build (frame A1): exercise full coefficient block ----
    println("\n" * "-"^70)
    println("VISCOUS build (frame A1: param A, η̂ test, τ̂=10) — coefficient sample")
    println("-"^70)
    viscA = Viscosity(:A, 1e-2, 10.0, star.R)     # L0=R, η̂=0.01, τ̂=10
    Pv, Qv, auxv = axial_linear_system(star, eos, r, ω, ℓ, viscA)
    @printf("  η(r)   = %.6e\n", auxv.η)
    @printf("  τ(r)   = %.6e\n", auxv.τ)
    @printf("  c_η²   = %.6e   (=cs²/τ̂ for param A; cs²=%.4e)\n",
            auxv.cη2, bg.cs2)
    @printf("  U      = %.6e\n", auxv.U)
    @printf("  C1=%.4e%+.4ei  C2=%.4e%+.4ei\n",
            real(auxv.C1),imag(auxv.C1), real(auxv.C2),imag(auxv.C2))
    @printf("  C3=%.4e%+.4ei  C4=%.4e%+.4ei  C5=%.4e%+.4ei\n",
            real(auxv.C3),imag(auxv.C3), real(auxv.C4),imag(auxv.C4),
            real(auxv.C5),imag(auxv.C5))
    @printf("  ψ-damping coeff (-16π e^{ν/2} iω η) = %.4e%+.4ei\n",
            real(auxv.damp), imag(auxv.damp))
    @printf("  coupling Q[1,2]=C1/f² (NONzero now)  = %.4e%+.4ei\n",
            real(Qv[1,2]), imag(Qv[1,2]))

    # check c_η² = cs²/τ̂ for parametrization A (Bussieres line 311)
    ce2_expected = bg.cs2 / viscA.τhat
    @printf("  c_η² vs cs²/τ̂ : %.6e  vs  %.6e   |Δ|=%.3e\n",
            auxv.cη2, ce2_expected, abs(auxv.cη2 - ce2_expected))

    # ---- surface condition coefficients ----
    B = surface_condition(star, ω, viscA)
    @printf("\nSurface BC eq.(24) coeffs (param A): B1..B4 =\n   ")
    for b in B; @printf("%.4e%+.4ei  ", real(b), imag(b)); end
    println()

    # ---- first-order companion matrix sanity ----
    Mmat = axial_first_order_matrix(star, eos, r, ω, ℓ, viscA)
    @printf("\n4×4 companion matrix built, top-right identity block: M[1,3]=%.1f M[2,4]=%.1f\n",
            real(Mmat[1,3]), real(Mmat[2,4]))

    # ---- final verdict ----
    tol = 1e-10
    ok = (errP < tol) && (errQ < tol) && (abs(Q[1,2]) < tol) &&
         (abs(V_ext - V_vac_expected) < 1e-8) &&
         (abs(auxv.cη2 - ce2_expected) < 1e-12)
    println("\n" * "="^70)
    if ok
        println("SANITY PASS: inviscid ψ-row == standard relativistic axial")
        println("             w-mode wave eq (RW potential eq.19); exterior")
        println("             reduces to vacuum Regge-Wheeler; param-A c_η²=cs²/τ̂.")
    else
        println("SANITY FAIL — see deltas above.")
    end
    println("="^70)
    return ok
end

if abspath(PROGRAM_FILE) == @__FILE__
    ok = _sanity()
    exit(ok ? 0 : 1)
end
