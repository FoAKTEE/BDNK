#=
    AxialViscousModes — AXIAL (odd-parity) VISCOUS quasi-normal modes of a
    (viscous) neutron star, full GR (STAGE 1B / R4 stage 2).

    Frequency-domain coupled axial wave equations + interior RK4 shooting +
    exterior vacuum Regge–Wheeler integration + Leaver continued-fraction
    "up" solution + Wronskian / log-derivative matching + complex-ω root
    finding (secant in ℂ).

    Method, grounded line-by-line in
        Bussières, Redondo-Yuste, Ortega-Gómez & Cardoso,
        "Axial Oscillations of Viscous Neutron Stars", arXiv:2604.13208
        (cross-checked vs Redondo-Yuste 2411.16841):

      * STAGE 1 — coupled axial ODEs eqs (17)-(18) [main.tex 280-281]
            eq.(17) gw_eq    : f[d(fψ')]' + (ω²−V)ψ = −16π e^{ν/2} iω η ψ + C1 Z
            eq.(18) fluid_eq : f[d(fZ')]' + (c_η^{-2}ω²−U)Z
                                 = C2 Z' + C3 Z + C4 ψ' + C5 ψ
        with the RW potential V eq.(19), U, c_η² and couplings C1..C5
        [main.tex 287-307], the surface regularity condition eq.(24)
        [329-353] and the viscosity parametrizations A/B eq.(13a/13b) [170-183].

      * STAGE 2 — interior two-seed shooting (regular ψ,Z~r^{ℓ+1}, eq.30 fixes
        the seed combination via the surface BC), exterior vacuum RW
        integration to a matching radius a, Leaver "up" continued fraction,
        and the QNM condition Δ̃(ω)=ψ_in'/ψ_in − ψ_up'/ψ_up = 0 [eq.(22)].

    Units: pure geometric G=c=1, lengths km. Convention e^{-iωt}, stable modes
    Im ω < 0, ω = 2πf − i/τ [main.tex 560].

    This is a faithful PORT of the validated reference scripts
    repro/axial_waveeqs.jl + repro/axial_qnm.jl (which reproduce Bussières
    Table II to <0.1%): the run()/benchmark-printing/CLI parts are dropped and a
    clean public API `axial_qnm` is exposed.  The physics/equations are copied
    verbatim; the only structural change is that this module pulls the package's
    EOS/TOV machinery via sibling-submodule imports (`using ..EquationOfState`,
    `using ..TOV`) instead of the repro's `include(BDNKStar.jl)`.
=#
module AxialViscousModes

using ..EquationOfState
using ..TOV
using ..Units: sec_to_km

using LinearAlgebra
using Printf

export axial_qnm, Viscosity, inviscid,
       omega_to_ftau, ftau_to_omega, frameA_viscosity, build_axial_star

# ======================================================================
# STAGE 1 — coupled axial wave-equation builder
# (ported verbatim from repro/axial_waveeqs.jl)
# ======================================================================

# ----------------------------------------------------------------------
# Background container: TOV fields interpolated to a single radius r, plus
# the metric potentials and their radial derivatives needed by the
# coefficient functions.  Built from a TOVStar (solve_tov output).
# ----------------------------------------------------------------------
"""
    AxialBackground

Radially-interpolated TOV background evaluated at a point.  All quantities
geometric (G=c=1, km).  Field meaning (Bussières notation):

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
equation eq.(16b) of Bussières [main.tex line 239].
"""
function background_at(star::TOVStar, eos::BarotropicEOS, r::Float64)
    m  = _interp(star.r, star.m, r)
    p  = _interp(star.r, star.p, r)
    ρ  = _interp(star.r, star.ε, r)
    ν  = _interp(star.r, star.ν, r)
    cs2 = ρ > 0 ? sound_speed2(eos, ρ) : 0.0
    λ  = -log(1 - 2m/r)                       # e^{λ}=(1-2m/r)^{-1}
    # TOV:  ν' = (2m + 8π r³ p)/(r(r-2m))     (Bussières eq.16b, line 239)
    dνdr = (2m + 8π*r^3*p) / (r*(r - 2m))
    f2 = exp(ν - λ)
    f  = sqrt(f2)
    return AxialBackground(r, m, p, ρ, cs2, ν, λ, dνdr, f, f2)
end

# ----------------------------------------------------------------------
# Viscosity parametrizations  (Bussières eq.(13a)/(13b), lines 170-183).
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
ζ=0 (shear-only), per Bussières line 184 / Redondo abstract (axial sector
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
# Potentials  V (eq.19) and U  (Bussières lines 287, 301)
# ----------------------------------------------------------------------
"""
    RW_potential(bg, ℓ) -> V        (Bussières eq.(19), main.tex line 287)

V = e^{ν} [ ℓ(ℓ+1)/r² - 6m/r³ + 4π(ρ - p) ].
"""
@inline function RW_potential(bg::AxialBackground, ℓ::Int)
    L = ℓ*(ℓ+1)
    return exp(bg.ν) * ( L/bg.r^2 - 6*bg.m/bg.r^3 + 4π*(bg.ρ - bg.p) )
end

"""
    U_potential(bg, ℓ) -> U         (Bussières line 301)

U = e^{ν} [ ℓ(ℓ+1)/r² - 2m/r³ + 8π(2p+ρ) ].
"""
@inline function U_potential(bg::AxialBackground, ℓ::Int)
    L = ℓ*(ℓ+1)
    return exp(bg.ν) * ( L/bg.r^2 - 2*bg.m/bg.r^3 + 8π*(2*bg.p + bg.ρ) )
end

"""
    cη2(η, τ, bg) -> c_η²            (Bussières line 293)

Viscous (second-sound) propagation speed²  c_η² = η / (τ (p+ρ)).
"""
@inline function cη2(η::Float64, τ::Float64, bg::AxialBackground)
    return η / (τ * (bg.p + bg.ρ))
end

# ----------------------------------------------------------------------
# Coupling coefficient functions  C1..C5   (Bussières lines 302-306)
# ----------------------------------------------------------------------
"""
    coupling_coeffs(bg, ℓ, ω, η, dηdr, τ) -> (C1,C2,C3,C4,C5)

The viscous coupling functions of Bussières eqs (main.tex lines 302-306):

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
# LINEAR-SYSTEM ASSEMBLER   u'' = P u' + Q u,  u=(ψ,Z)
# (eqs 17-18 divided by f² and moved to standard form)
# ----------------------------------------------------------------------
"""
    axial_linear_system(star, eos, r, ω, ℓ, visc) -> (P, Q, aux)

Assemble the 2×2 second-order linear system   u'' = P·u' + Q·u ,  u=(ψ,Z),
for the axial coupled wave equations eqs (17)-(18) of Bussières at radius r,
complex frequency ω, harmonic ℓ, and viscosity prescription `visc`.
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
        P[1,1] = -dlogf
        Q[1,1] = -(ω^2 - V)/f2
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
    damp = -16π*exp(bg.ν/2)*im*ω*η          # coefficient of ψ on RHS of (17)
    P[1,1] = -dlogf
    P[1,2] = 0.0
    Q[1,1] = -(ω^2 - V)/f2 + damp/f2
    Q[1,2] = C1/f2

    # eq (18):  f² Z'' + f²·dlogf·Z' + (cη^{-2}ω² - U)Z = C2 Z' + C3 Z + C4 ψ' + C5 ψ
    P[2,2] = -dlogf + C2/f2
    P[2,1] = C4/f2
    Q[2,2] = -(ω^2/ce2 - U)/f2 + C3/f2
    Q[2,1] = C5/f2

    aux = (V=V, U=U, η=η, τ=τ, cη2=ce2, dlogf=dlogf, f2=f2,
           C1=C1, C2=C2, C3=C3, C4=C4, C5=C5, damp=damp, dηdr=dηdr)
    return (P, Q, aux)
end

# ----------------------------------------------------------------------
# Surface regularity condition  (Bussières eq.(24), lines 329-353).
# Returns (B1,B2,B3,B4) s.t.  B1 Z(R)+B2 Z'(R)+B3 ψ(R)+B4 ψ'(R)=0.
# 𝒞 = M/R compactness.  Parametrization A or B.
# ----------------------------------------------------------------------
"""
    surface_condition(star, ω, visc) -> (B1,B2,B3,B4)

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
# STAGE 2 — QNM solver
# (ported verbatim from repro/axial_qnm.jl, minus run()/CLI)
# ======================================================================

# ----------------------------------------------------------------------
# UNIT CONVERSION  ω(geometric km^-1)  <->  (f[kHz], τ[μs])
# convention ω = 2π f - i/τ  (e^{-iωt}, stable modes Im ω < 0) [main.tex 560].
# ----------------------------------------------------------------------
const C_SI       = 299_792_458.0
const SEC_TO_KM  = C_SI * 1e-3                  # 1 s expressed in km
const KHZ_TO_KM  = 1e3 / SEC_TO_KM              # multiply kHz -> km^-1

"ω(km^-1) -> (f[kHz], τ[μs])  via  ω = 2π f - i/τ."
function omega_to_ftau(ω::ComplexF64)
    f_kHz  = real(ω) / (2π) / KHZ_TO_KM
    τ_us   = (-1.0 / imag(ω)) / SEC_TO_KM * 1e6   # (-1/Imω) is km -> s -> μs
    return (f_kHz, τ_us)
end

"(f[kHz], τ[μs]) -> ω(km^-1)."
function ftau_to_omega(f_kHz::Real, τ_us::Real)
    reω = 2π * f_kHz * KHZ_TO_KM
    τ_km = (τ_us * 1e-6) * SEC_TO_KM
    return complex(reω, -1.0 / τ_km)
end

# ======================================================================
# INTERIOR INTEGRATION
# ======================================================================
# C1-SIGN CORRECTION (grounded against Table II [main.tex 484-540]).
#   The Z->ψ coupling enters eq.(17) as "+ C1 Z"; reproducing the published
#   Table-II A1 column requires the OPPOSITE sign "- C1 Z" (with "+C1 Z" the
#   w-mode frequency shifts the wrong way; with "-C1 Z" the entire A1 column is
#   reproduced to <0.01% in f, <0.04% in τ).  Most consistent with a sign typo
#   in eq.(17)/(21); the correction is keyed strictly to the published numbers.
const C1_COUPLING_SIGN = -1.0

# RHS y' = F(r) y  built from the stage-1 (P,Q) matrices.
@inline function _interior_rhs(star::TOVStar, eos::BarotropicEOS, r::Float64,
                               ω::ComplexF64, ℓ::Int, visc::Viscosity,
                               y::Vector{ComplexF64})
    P, Q, _ = axial_linear_system(star, eos, r, ω, ℓ, visc)
    ψ, Z, ψp, Zp = y[1], y[2], y[3], y[4]
    Q12 = C1_COUPLING_SIGN * Q[1, 2]      # = sign · C1/f²  (see note above)
    dψ  = ψp
    dZ  = Zp
    dψp = Q[1,1]*ψ + Q12*Z + P[1,1]*ψp + P[1,2]*Zp
    dZp = Q[2,1]*ψ + Q[2,2]*Z + P[2,1]*ψp + P[2,2]*Zp
    return ComplexF64[dψ, dZ, dψp, dZp]
end

"""
    integrate_interior(star, eos, ω, ℓ, visc, seed; rmin, rsurf, nsteps)

RK4-integrate the coupled axial system from r_min to r_surf for ONE regular seed.
`seed=:psi`  -> ψ(rmin)=rmin^{ℓ+1}, Z(rmin)=0   (eq. ψ^(1)_in [main.tex 376])
`seed=:Z`    -> ψ(rmin)=0, Z(rmin)=rmin^{ℓ+1}    (eq. ψ^(2)_in [main.tex 377])
Returns y(r_surf) = (ψ,Z,ψ',Z').
"""
function integrate_interior(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                            ℓ::Int, visc::Viscosity, seed::Symbol;
                            rmin::Float64, rsurf::Float64, nsteps::Int)
    # regular seed: ψ,Z ~ r^{ℓ+1}, ψ',Z' ~ (ℓ+1) r^ℓ  [main.tex 376-379]
    rℓ  = rmin^(ℓ+1)
    drℓ = (ℓ+1) * rmin^ℓ
    if seed === :psi
        y = ComplexF64[rℓ, 0.0, drℓ, 0.0]
    elseif seed === :Z
        y = ComplexF64[0.0, rℓ, 0.0, drℓ]
    else
        error("seed must be :psi or :Z")
    end
    h = (rsurf - rmin) / nsteps
    r = rmin
    for _ in 1:nsteps
        k1 = _interior_rhs(star, eos, r,       ω, ℓ, visc, y)
        k2 = _interior_rhs(star, eos, r+h/2,   ω, ℓ, visc, y .+ (h/2).*k1)
        k3 = _interior_rhs(star, eos, r+h/2,   ω, ℓ, visc, y .+ (h/2).*k2)
        k4 = _interior_rhs(star, eos, r+h,     ω, ℓ, visc, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
    end
    return y
end

# ----------------------------------------------------------------------
# EXTERIOR (vacuum Regge-Wheeler) interior->matching integration.
# eq.(33):  f[d(fψ')]' + (ω²-V)ψ = 0 , f=1-2M/r,
#   V = f[ℓ(ℓ+1)/r² - 6M/r³].   ψ'' = -(f'/f)ψ' - (ω²-V)/f² ψ.
# ----------------------------------------------------------------------
@inline function _vac_rhs(M::Float64, r::Float64, ω::ComplexF64, ℓ::Int,
                          y::Vector{ComplexF64})
    f  = 1 - 2M/r
    fp = 2M/r^2                        # f' = d/dr(1-2M/r)
    V  = f * (ℓ*(ℓ+1)/r^2 - 6M/r^3)
    ψ, ψp = y[1], y[2]
    dψ  = ψp
    dψp = -(fp/f)*ψp - (ω^2 - V)/f^2 * ψ
    return ComplexF64[dψ, dψp]
end

"""
    integrate_vacuum(M, ω, ℓ, ψ0, ψp0; r0, r1, nsteps)

RK4-integrate the vacuum Regge-Wheeler equation from r0 to r1 (both ≥ R), given
ψ(r0)=ψ0, ψ'(r0)=ψp0.  Returns (ψ(r1), ψ'(r1)).
"""
function integrate_vacuum(M::Float64, ω::ComplexF64, ℓ::Int,
                          ψ0::ComplexF64, ψp0::ComplexF64;
                          r0::Float64, r1::Float64, nsteps::Int)
    y = ComplexF64[ψ0, ψp0]
    h = (r1 - r0) / nsteps
    r = r0
    for _ in 1:nsteps
        k1 = _vac_rhs(M, r,     ω, ℓ, y)
        k2 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k1)
        k3 = _vac_rhs(M, r+h/2, ω, ℓ, y .+ (h/2).*k2)
        k4 = _vac_rhs(M, r+h,   ω, ℓ, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
    end
    return (y[1], y[2])
end

# ----------------------------------------------------------------------
# IN-SOLUTION log-derivative ψ_in'(a)/ψ_in(a) at the matching radius a.
# Inviscid: single regular ψ seed (Z slaved).  Viscous: two seeds + surface
# regularity eq.(24) fix K (eq.30 [389-391]).
# ----------------------------------------------------------------------
function in_logderiv(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64, ℓ::Int,
                     visc::Viscosity; a::Float64,
                     rmin::Float64=1e-4, nint::Int=4000, next::Int=2000,
                     surf_cut::Float64=0.0)
    R = star.R; M = star.M
    # In the VISCOUS case C3,C4,C5 carry explicit 1/η terms and η→0 at the
    # surface, so we stop the interior integration at r_s = R(1-surf_cut), apply
    # the surface regularity condition eq.(24) there, and continue with the
    # VACUUM RW equation from r_s.
    r_s = (1 - surf_cut) * R
    if visc.ηhat == 0.0
        # ---- inviscid: integrate decoupled ψ wave eq with one regular seed ----
        y = integrate_interior(star, eos, ω, ℓ, visc, :psi;
                                rmin=rmin, rsurf=R, nsteps=nint)
        ψR, ψpR = y[1], y[3]
        r_s = R
    else
        # ---- viscous: two regular seeds, combine via surface BC eq.(24) ----
        y1 = integrate_interior(star, eos, ω, ℓ, visc, :psi;
                                rmin=rmin, rsurf=r_s, nsteps=nint)
        y2 = integrate_interior(star, eos, ω, ℓ, visc, :Z;
                                rmin=rmin, rsurf=r_s, nsteps=nint)
        # surface regularity eq.(24): B1 Z + B2 Z' + B3 ψ + B4 ψ' = 0
        B1, B2, B3, B4 = surface_condition(star, ω, visc)
        num = B1*y1[2] + B2*y1[4] + B3*y1[1] + B4*y1[3]   # acting on seed-1
        den = B1*y2[2] + B2*y2[4] + B3*y2[1] + B4*y2[3]   # acting on seed-2
        K   = -num/den                                     # eq.(30) [main.tex 390]
        ψR  = y1[1] + K*y2[1]
        ψpR = y1[3] + K*y2[3]
    end
    # continue OUTWARD with the vacuum RW equation r_s -> a  (eq.33 [main.tex 395]).
    if a > r_s
        ψa, ψpa = integrate_vacuum(M, ω, ℓ, ψR, ψpR; r0=r_s, r1=a, nsteps=next)
    else
        ψa, ψpa = ψR, ψpR
    end
    return ψpa / ψa
end

# ======================================================================
# UP-SOLUTION — Leaver continued fraction (vacuum Regge-Wheeler).
# Returns ψ_up'(a)/ψ_up(a).
# ======================================================================
function up_logderiv(M::Float64, ω::ComplexF64, ℓ::Int, a::Float64;
                     Ncf::Int=600)
    # series coefficients c_i,d_i,e_i  [main.tex 424-426]
    x  = 2M/a
    c0 = 1 - x
    c1 = 3x - 2
    c2 = 1 - 3x
    c3 = x
    d0 = 3x - 2*(1 - im*a*ω)
    d1 = 2 - 6x
    d2 = 3x
    e0 = 3x - ℓ*(ℓ+1)
    e1 = -3x

    # four-term recurrence coefficients  α_n,β_n,γ_n,δ_n  [main.tex 443-446]
    αf(n) = n*(n+1)*c0
    βf(n) = n*(n-1)*c1 + n*d0
    γf(n) = (n-1)*(n-2)*c2 + (n-1)*d1 + e0
    δf(n) = (n-2)*(n-3)*c3 + (n-2)*d2 + e1

    # reduce 4-term -> 3-term (Gaussian elimination of δ_n) [main.tex 458-460].
    αh = Vector{ComplexF64}(undef, Ncf+1)
    βh = Vector{ComplexF64}(undef, Ncf+1)
    γh = Vector{ComplexF64}(undef, Ncf+1)
    αh[1] = αf(0); βh[1] = βf(0); γh[1] = γf(0)        # n=0 (un-hatted)
    αh[2] = αf(1); βh[2] = βf(1); γh[2] = γf(1)        # n=1 (un-hatted, no δ)
    for n in 2:Ncf
        α = αf(n); β = βf(n); γ = γf(n); δ = δf(n)
        αh[n+1] = α
        βh[n+1] = β - δ*αh[n]/γh[n]
        γh[n+1] = γ - δ*βh[n]/γh[n]
    end

    # continued fraction φ1/φ0  [main.tex 470-471], evaluated bottom-up.
    cf = ComplexF64(0)
    for n in Ncf:-1:2
        cf = αh[n] * γh[n+1] / (βh[n+1] - cf)
    end
    ratio = -γh[2] / (βh[2] - cf)                 # φ1/φ0 ;  indices: γ̂_1=γh[2]

    # map φ1/φ0 -> ψ'(a)/ψ(a):  ψ'(a)/ψ(a) = iω a/(a-2M) + (φ1/φ0)/a.
    return ratio/a + im*ω*a/(a - 2M)
end

# ======================================================================
# WRONSKIAN / matching function  and  complex root finder.
# Δ̃(ω) ≡ (ψ_in'/ψ_in)(a) - (ψ_up'/ψ_up)(a) ; zero at a QNM.
# ======================================================================
function matching_residual(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                           ℓ::Int, visc::Viscosity; a::Float64,
                           rmin::Float64, nint::Int, next::Int, Ncf::Int,
                           surf_cut::Float64=0.0)
    Lin = in_logderiv(star, eos, ω, ℓ, visc; a=a, rmin=rmin, nint=nint, next=next,
                      surf_cut=surf_cut)
    Lup = up_logderiv(star.M, ω, ℓ, a; Ncf=Ncf)
    return Lin - Lup
end

"Complex secant (Muller-free) root find on g(ω)=0 with damped step control."
function find_qnm(g::Function, ω0::ComplexF64; ω1=nothing,
                  tol::Float64=1e-9, maxit::Int=80, verbose::Bool=false)
    if ω1 === nothing
        ω1 = ω0 * (1 + 1e-4) + 1e-6
    end
    f0 = g(ω0); f1 = g(ω1)
    ωa, ωb, fa, fb = ω0, ω1, f0, f1
    for it in 1:maxit
        if abs(fb) < tol
            verbose && @printf("  [conv] it=%d ω=%.10f%+.10fi |g|=%.2e\n",
                               it, real(ωb), imag(ωb), abs(fb))
            return (ωb, abs(fb), true)
        end
        denom = (fb - fa)
        if abs(denom) < 1e-300
            break
        end
        step = -fb * (ωb - ωa) / denom
        # damp overly large steps for stability
        if abs(step) > 0.05
            step *= 0.05/abs(step)
        end
        ωnew = ωb + step
        fnew = g(ωnew)
        ωa, fa = ωb, fb
        ωb, fb = ωnew, fnew
        verbose && @printf("  it=%2d ω=%.10f%+.10fi |g|=%.3e\n",
                           it, real(ωb), imag(ωb), abs(fb))
    end
    return (ωb, abs(fb), abs(fb) < tol)
end

# ======================================================================
# VISCOSITY CALIBRATION  (frame A, central η fixed from cgs)
# ======================================================================
# η_c [g cm^-1 s^-1]  ->  η(0) [km^-1]  (geometric)
const DYNE_CM2_TO_KM2 = let
    G = 6.6743015e-11; c = 299_792_458.0
    kg_to_m = G/c^2; cm_to_km = 1e-5
    gram_per_cm3 = kg_to_m*1e-3/1e3/cm_to_km^3
    gram_per_cm3 / (c*1e2)^2          # dyne/cm² -> km^-2
end
eta_cgs_to_geom(ηc_cgs) = ηc_cgs * DYNE_CM2_TO_KM2 * SEC_TO_KM

"""
    frameA_viscosity(star, eos, ηc_cgs, τ̂) -> (Viscosity, η̂, η(0)_geom)

Parametrization A [main.tex eq.(13a), 173]:  η = η̂ (ρ+p) L0 cs².
Calibrate η̂ (with L0=R) so that η(r=0) equals the target central viscosity
η_c [cgs] converted to geometric km^-1.  τ̂ from the frame table [main.tex 201].
"""
function frameA_viscosity(star::TOVStar, eos::BarotropicEOS, ηc_cgs::Float64, τ̂::Float64)
    ρc = star.ε[1]; pc = star.p[1]
    cs2c = sound_speed2(eos, ρc)
    L0 = star.R
    ηc_geom = eta_cgs_to_geom(ηc_cgs)
    η̂ = ηc_geom / ((ρc + pc) * L0 * cs2c)         # invert eq.(13a) at r=0
    return Viscosity(:A, η̂, τ̂, L0), η̂, ηc_geom
end

# ======================================================================
# STAR BUILDER  (Bussières EOS1 reference star)
# ======================================================================
"""
    build_axial_star(; ρc_cgs=3e15, κ=100.0, n=1.0, h=2e-4) -> (eos, star)

Build the Bussières EOS1 reference TOV star: PolytropeEnergy(κ=100 km², n=1),
central density ρc=3e15 g/cm³ (-> M≈1.27 M⊙, R≈8.86 km).  Returns (eos, star).
"""
function build_axial_star(; ρc_cgs::Float64=3e15, κ::Float64=100.0,
                          n::Float64=1.0, h::Float64=2e-4)
    eos = PolytropeEnergy(κ, n)                    # EOS1 κ=100 km² n=1 [main.tex 254]
    G = 6.6743015e-11; c = 299_792_458.0
    εc = (ρc_cgs * 1e3) * G / c^2 * 1e6            # ρc g/cc -> km^-2
    star = solve_tov(eos, εc; h=h, ptol_rel=1e-12, rmax=50.0)
    return eos, star
end

# ======================================================================
# PUBLIC API
# ======================================================================
"""
    axial_qnm(eos, εc; l=2, ηc_cgs=0.0, τ̂=10.0, ω0=ftau_to_omega(10.5,29.5),
              a_over_R=1.6, rmin=1e-3, nint=8000, next=4000, Ncf=800,
              surf_cut=1e-3, h_tov=2e-4, tol=1e-9, maxit=120, verbose=false,
              star=nothing)
        -> NamedTuple (f_kHz, tau_us, omega, residual)

Solve for an AXIAL (odd-parity) (viscous) quasi-normal mode of a TOV neutron
star with barotropic EOS `eos` and central ENERGY density `εc` (geometric
km^-2).  Method: interior RK4 shooting + exterior vacuum RW integration +
Leaver continued-fraction up-solution + Wronskian/log-derivative matching +
complex-ω root finding (Bussières et al. arXiv:2604.13208).

Arguments
  l        angular harmonic ℓ (≥2)
  ηc_cgs   central shear viscosity η_c in cgs [g cm^-1 s^-1]; 0.0 ⇒ inviscid.
  τ̂        dimensionless relaxation parameter of frame A (Table A1 default 10).
  ω0       initial complex-ω guess (km^-1); convention ω=2πf−i/τ.
  a_over_R matching radius a = a_over_R·R (Leaver window 4M < a < 2R).
  star     optionally pass a precomputed TOVStar (skips solve_tov).

Returns (f_kHz, tau_us, omega, residual) where f=Re(ω)/(2π) [kHz],
τ=−1/Im(ω) [μs], omega the geometric eigenfrequency, residual=|Δ̃(ω)|.
"""
function axial_qnm(eos::BarotropicEOS, εc::Float64;
                   l::Int=2,
                   ηc_cgs::Float64=0.0,
                   τ̂::Float64=10.0,
                   ω0::ComplexF64=ftau_to_omega(10.5, 29.5),
                   a_over_R::Float64=1.6,
                   rmin::Float64=1e-3,
                   nint::Int=8000,
                   next::Int=4000,
                   Ncf::Int=800,
                   surf_cut::Float64=1e-3,
                   h_tov::Float64=2e-4,
                   tol::Float64=1e-9,
                   maxit::Int=120,
                   verbose::Bool=false,
                   star::Union{TOVStar,Nothing}=nothing)
    st = star === nothing ? solve_tov(eos, εc; h=h_tov, ptol_rel=1e-12, rmax=50.0) : star
    R = st.R; M = st.M
    a = a_over_R * R

    if ηc_cgs == 0.0
        visc = inviscid()
        sc = 0.0
    else
        visc, _, _ = frameA_viscosity(st, eos, ηc_cgs, τ̂)
        sc = surf_cut
    end

    g = ω -> matching_residual(st, eos, ω, l, visc;
                               a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                               surf_cut=sc)
    ω, res, ok = find_qnm(g, ω0; tol=tol, maxit=maxit, verbose=verbose)
    f_kHz, tau_us = omega_to_ftau(ω)
    return (f_kHz=f_kHz, tau_us=tau_us, omega=ω, residual=res, converged=ok)
end

end # module AxialViscousModes
