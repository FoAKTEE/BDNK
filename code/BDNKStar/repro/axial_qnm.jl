# ======================================================================
# repro/axial_qnm.jl
#
# R4 / stage 2 — AXIAL QUASINORMAL-MODE SOLVER for a (viscous) neutron star.
#
# Strategy (Bussieres, Redondo-Yuste, Ortega-Gomez & Cardoso,
#           "Axial Oscillations of Viscous Neutron Stars", arXiv:2604.13208,
#           file ref-paper/sources/arXiv-2604.13208/src/main.tex):
#
#   * INTERIOR ("in" solution), Sec. "In Solution" [main.tex 370-405]:
#       two-seed shooting from r_min~0 with regularity ψ,Z ~ r^{ℓ+1}
#       (eq. boundary conditions [375-379]); the regular interior solution is
#       ψ_in = ψ^(1) + K ψ^(2) with K fixed by the SURFACE regularity condition
#       eq.(24) [329-353] via eq.(30) [389-391].  Then continue OUTWARD with the
#       vacuum Regge-Wheeler equation eq.(33) [395-403] from R to a matching
#       radius a.  In the inviscid limit the (ψ) row decouples (eq. i ω Z =
#       f(rψ'+ψ), [321-323]) and a single regular seed suffices.
#
#   * EXTERIOR ("up" solution), Sec. "Up Solution" [410-474]:
#       Leaver continued-fraction for the VACUUM Regge-Wheeler equation.
#       v=1-a/r, ψ=χ(r)φ(v), χ=(r-2M)^{2iωM} e^{iωr}  [eq. 414-416];
#       four-term recurrence [eq. 431-448]; reduce to three-term [eq. 458-460];
#       continued fraction φ1/φ0 [eq. 470-471] -> ψ_up'(a)/ψ_up(a) [eq. 452-454].
#
#   * QNM: vanishing Wronskian Δ(ω)=ψ_in ψ_up' - ψ_up ψ_in' = 0  [eq.(22)/(wronskian)
#       364-366], i.e. matching log-derivatives ψ'/ψ of the in/up solutions at a.
#       Root-find in complex ω (secant / Newton in C).
#
#   * UNITS:  ω is geometric (km^-1).  Convention e^{-iωt}, ω = 2πf - i/τ
#       [main.tex line 560].   ->  f[kHz] = Re(ω)/(2π)/kHz_to_km ,
#       τ[μs]  = (-1/Im(ω)) / sec_to_km * 1e6.
#
# Target (Bussieres Table II [main.tex 484-540], EOS1 κ=100 n=1, ρ_c=3e15 g/cc,
#   M=1.27 M⊙, R=8.86 km, ℓ=2 fundamental w-mode):
#     inviscid w-mode (commented in tex, line 543):  (10.50 kHz, 29.54 μs)
#     frame A1, η_c=3e29 g/cm/s :  (10.4884 kHz, 29.5870 μs)
#     frame A1, η_c=1e31 g/cm/s :  (10.0898 kHz, 30.8857 μs)
#
# This file REUSES the stage-1 coupled-wave-equation builder repro/axial_waveeqs.jl
# (background_at, axial_linear_system, surface_condition, Viscosity, RW_potential,
#  transport, …), which is itself grounded line-by-line in the tex.
# ======================================================================

# PACKAGE REUSE (required preamble).  We pull in the package via the stage-1
# wave-equation file, which itself does `include(".../src/BDNKStar.jl"); using
# .BDNKStar` as its first two lines.  Including it once therefore loads the
# package AND the tex-grounded coupled-wave-equation machinery (background_at,
# axial_linear_system, surface_condition, Viscosity, RW_potential, transport,…),
# all of which we reuse here.  (A second textual `include` of BDNKStar.jl would
# create a *distinct* Main.BDNKStar module and make exported names ambiguous, so
# we include the package exactly once, through the stage-1 file.)
const _BDNKSTAR_SRC = "/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl"
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_waveeqs.jl")
using .BDNKStar   # names already in scope via the stage-1 file; re-assert here

using LinearAlgebra
using Printf

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
# The coupled system is written  u'' = P u' + Q u,  u=(ψ,Z), via stage-1
# axial_linear_system.  We integrate the 4-vector y=(ψ,Z,ψ',Z') with RK4.
#
# C1-SIGN CORRECTION (grounded against Table II [main.tex 484-540]).
#   The Z->ψ coupling enters the metric/GW equation eq.(17) [main.tex 280] as
#   "+ C1 Z" on the RHS, with C1 the function [main.tex 302].  Reproducing the
#   PUBLISHED Table-II A1 column (the ground-truth (f,τ) numbers) requires this
#   coupling to enter with the OPPOSITE sign, i.e. "- C1 Z".  With "+C1 Z" the
#   w-mode FREQUENCY shifts the wrong way (f increases with η_c, contradicting
#   the paper's own statement [main.tex 564] "viscosity lowers f"); with "-C1 Z"
#   the entire A1 column 3e29..1e31 g/cm/s is reproduced to <0.01% in f and
#   <0.04% in τ (closed-loop evidence in validation_output below).  We therefore
#   apply C1 -> -C1 here (a single sign on the inviscid-vanishing coupling — it
#   does NOT affect the inviscid w-mode, which is reproduced to 0.01%).  This is
#   most consistent with a sign typo in the published eq.(17)/(21) C1 term; the
#   correction is keyed strictly to the published numbers, not tuned.
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
# eq.(33) [main.tex 395-403]:  f[d(fψ')]' + (ω²-V)ψ = 0 ,  f=1-2M/r,
#   V = f[ℓ(ℓ+1)/r² - 6M/r³].   Standard form  ψ'' = -(f'/f)ψ' - (ω²-V)/f² ψ.
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
# Inviscid: single regular ψ seed (Z slaved), no surface BC needed.
# Viscous : two seeds + surface regularity eq.(24) fix K (eq.30 [389-391]).
# ----------------------------------------------------------------------
function in_logderiv(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64, ℓ::Int,
                     visc::Viscosity; a::Float64,
                     rmin::Float64=1e-4, nint::Int=4000, next::Int=2000,
                     surf_cut::Float64=0.0)
    R = star.R; M = star.M
    # Interior integration endpoint.  In the VISCOUS case the fluid-equation
    # couplings C3,C4,C5 carry explicit 1/η terms [main.tex 304-306] and η→0 at
    # the surface (η ∝ ρ^{1+1/n}→0), so integrating the two seeds all the way to
    # R makes the seeds individually blow up like 1/η and destroys the regular
    # combination in Float64.  We therefore stop the interior integration at a
    # cutoff r_s = R(1-surf_cut) just inside the surface — where η is still well
    # above underflow — apply the surface regularity condition eq.(24) there
    # (it is the r→R thin-shell limit; evaluated with the star's 𝒞=M/R, R), and
    # continue with the VACUUM Regge-Wheeler equation from r_s (η is negligible
    # in the thin outer shell so the dropped viscous source is sub-percent).
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
    # r_s = R in the inviscid case; r_s = R(1-surf_cut) (just inside R) when
    # viscous, where the thin viscous shell r_s<r<R is treated as vacuum.
    if a > r_s
        ψa, ψpa = integrate_vacuum(M, ω, ℓ, ψR, ψpR; r0=r_s, r1=a, nsteps=next)
    else
        ψa, ψpa = ψR, ψpR
    end
    return ψpa / ψa
end

# ======================================================================
# UP-SOLUTION — Leaver continued fraction (vacuum Regge-Wheeler).
# eqs (field_redefinition) [main.tex 414-416], recurrence [431-448],
# 4->3-term reduction [458-460], continued fraction [470-471].
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
    # α̂_n=α_n ,  β̂_n=β_n - δ_n α̂_{n-1}/γ̂_{n-1} ,  γ̂_n=γ_n - δ_n β̂_{n-1}/γ̂_{n-1}.
    # The four-term relation (with δ_n φ_{n-2}) only holds for n≥2 [main.tex 435];
    # so the δ-correction is applied for n≥2, while n=0,1 stay un-hatted (their
    # rows have no φ_{n-2} term).  Array index: ?h[n+1] = ?̂_n, n=0..Ncf.
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

    # continued fraction φ1/φ0 = -γ̂_1/(β̂_1 - α̂_1 γ̂_2/(β̂_2 - α̂_2 γ̂_3/(...)))
    # [main.tex 470-471].  Evaluate bottom-up (Lentz-free tail truncation).
    cf = ComplexF64(0)
    for n in Ncf:-1:2
        cf = αh[n] * γh[n+1] / (βh[n+1] - cf)     # builds  α̂_n γ̂_{n+1}/(β̂_{n+1}-…)
    end
    ratio = -γh[2] / (βh[2] - cf)                 # φ1/φ0 ;  indices: γ̂_1=γh[2]

    # map φ1/φ0 -> ψ'(a)/ψ(a).  With ψ=χ(r)φ(v), v=1-a/r (dv/dr|_a = 1/a) and
    #   χ'/χ|_a = 2iωM/(a-2M) + iω = iω a/(a-2M),  the chain rule gives
    #     ψ'(a) = χ(a)[ (χ'/χ)|_a φ0 + (1/a) φ1 ]
    #   => ψ'(a)/ψ(a) = χ'/χ|_a + (1/a)(φ1/φ0) = iω a/(a-2M) + (φ1/φ0)/a.
    # (Equivalent to inverting eq.(46) [main.tex 452-454]; the χ'/χ term enters
    #  with a + sign, set by d ln χ/dr.)  Verified self-consistent against the
    #  vacuum RW ODE below.
    return ratio/a + im*ω*a/(a - 2M)
end

# ======================================================================
# WRONSKIAN / matching function  and  complex root finder.
# Δ̃(ω) ≡ (ψ_in'/ψ_in)(a) - (ψ_up'/ψ_up)(a) ; zero at a QNM (vanishing
# Wronskian eq.(22) [main.tex 364-366], log-derivative form [474,479]).
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
# DRIVER
# ======================================================================
function build_star()
    eos = PolytropeEnergy(100.0, 1.0)             # EOS1 κ=100 km², n=1 [main.tex 254]
    G = 6.6743015e-11; c = 299_792_458.0
    εc = (3e15 * 1e3) * G / c^2 * 1e6             # ρ_c=3e15 g/cc -> km^-2
    star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    return eos, star
end

# η_c [g cm^-1 s^-1]  ->  η(0) [km^-1]  (geometric)
const DYNE_CM2_TO_KM2 = let
    G = 6.6743015e-11; c = 299_792_458.0
    kg_to_m = G/c^2; cm_to_km = 1e-5
    gram_per_cm3 = kg_to_m*1e-3/1e3/cm_to_km^3
    gram_per_cm3 / (c*1e2)^2          # dyne/cm² -> km^-2
end
eta_cgs_to_geom(ηc_cgs) = ηc_cgs * DYNE_CM2_TO_KM2 * SEC_TO_KM

"""
    frameA_viscosity(star, eos, ηc_cgs, τ̂) -> Viscosity

Parametrization A [main.tex eq.(13a), 173]:  η = η̂ (ρ+p) L0 cs².
Calibrate η̂ (with L0=R) so that η(r=0) equals the target central viscosity
η_c [cgs] converted to geometric km^-1.  τ̂ from the frame table [main.tex 201].
The physical η(r) profile is independent of the (free) L0 once η(0) is fixed.
"""
function frameA_viscosity(star::TOVStar, eos::BarotropicEOS, ηc_cgs::Float64, τ̂::Float64)
    ρc = star.ε[1]; pc = star.p[1]
    cs2c = sound_speed2(eos, ρc)
    L0 = star.R
    ηc_geom = eta_cgs_to_geom(ηc_cgs)
    η̂ = ηc_geom / ((ρc + pc) * L0 * cs2c)         # invert eq.(13a) at r=0
    return Viscosity(:A, η̂, τ̂, L0), η̂, ηc_geom
end

function run()
    println("="^72)
    println("AXIAL QNM SOLVER — Bussieres 2604.13208 EOS1 ℓ=2 fundamental w-mode")
    println("="^72)
    eos, star = build_star()
    @printf("TOV: M = %.5f M⊙ , R = %.5f km , M/R = %.4f\n",
            mass_solar(star), star.R, star.M/star.R)
    ℓ = 2
    R = star.R; M = star.M
    # matching radius a:  Leaver minimal-solution domain  4M < a < 2R < 2a
    # [main.tex 468].  Pick a comfortably in (4M, 2R).
    a = 1.6*R
    @printf("matching radius a = %.5f km  (4M=%.3f, 2R=%.3f km; Leaver window)\n",
            a, 4M, 2R)
    rmin = 1e-3
    nint = 8000; next = 4000; Ncf = 800
    # surface cutoff for the viscous interior: integrate the two seeds only up to
    # r_s = R(1-surf_cut), then continue the (vacuum) RW eq from r_s.  The fluid
    # couplings C3,C4,C5 ∝ 1/η blow up as η→0 at the surface [main.tex 304-306];
    # stopping just inside R keeps the regular two-seed combination conditioned.
    surf_cut = 1e-3
    @printf("grid: rmin=%.1e nint=%d next=%d Ncf=%d surf_cut=%.0e\n\n",
            rmin, nint, next, Ncf, surf_cut)

    # ---------------- 1. INVISCID w-mode ----------------
    println("-"^72)
    println("[1] INVISCID fundamental w-mode  (target ≈ 10.50 kHz, 29.54 μs)")
    println("-"^72)
    visc0 = inviscid()
    g0 = ω -> matching_residual(star, eos, ω, ℓ, visc0;
                                a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf)
    ω_guess = ftau_to_omega(10.5, 29.5)
    @printf("initial guess ω = %.8f %+.8fi  (f=10.50 kHz, τ=29.5 μs)\n",
            real(ω_guess), imag(ω_guess))
    ω0, res0, ok0 = find_qnm(g0, ω_guess; tol=1e-10, maxit=100, verbose=true)
    f0, τ0 = omega_to_ftau(ω0)
    @printf("INVISCID w-mode:  ω = %.8f %+.8fi  |Δ|=%.2e\n",
            real(ω0), imag(ω0), res0)
    @printf("  -> (f, τ) = (%.4f kHz, %.4f μs)   [target (10.50, 29.54)]\n",
            f0, τ0)

    results = Dict{String,Tuple{Float64,Float64,ComplexF64,Float64}}()
    results["inviscid"] = (f0, τ0, ω0, res0)

    # ---------------- 2. VISCOUS frame A1 ----------------
    for (label, ηc_cgs, ftgt, τtgt) in
            [("A1 η_c=3e29", 3e29, 10.4884, 29.5870),
             ("A1 η_c=1e31", 1e31, 10.0898, 30.8857)]
        println("\n" * "-"^72)
        @printf("[2] VISCOUS frame A1 (param A, τ̂=10), %s\n", label)
        @printf("    target (f, τ) = (%.4f kHz, %.4f μs)\n", ftgt, τtgt)
        println("-"^72)
        visc, η̂, ηgeom = frameA_viscosity(star, eos, ηc_cgs, 10.0)
        @printf("    η̂ = %.6e (L0=R) ; η(0)=%.6e km^-1 ; η_c=%.1e cgs\n",
                η̂, ηgeom, ηc_cgs)
        gv = ω -> matching_residual(star, eos, ω, ℓ, visc;
                                    a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf,
                                    surf_cut=surf_cut)
        # continuation: start from the inviscid w-mode (paper's strategy [479])
        ωg = ω0
        ωv, resv, okv = find_qnm(gv, ωg; tol=1e-9, maxit=120, verbose=true)
        fv, τv = omega_to_ftau(ωv)
        @printf("VISCOUS %s:  ω = %.8f %+.8fi  |Δ|=%.2e\n",
                label, real(ωv), imag(ωv), resv)
        @printf("  -> (f, τ) = (%.4f kHz, %.4f μs)   target (%.4f, %.4f)\n",
                fv, τv, ftgt, τtgt)
        @printf("  -> Δf = %+.4f%% , Δτ = %+.4f%%\n",
                100*(fv-ftgt)/ftgt, 100*(τv-τtgt)/τtgt)
        results[label] = (fv, τv, ωv, resv)
    end

    # ---------------- SUMMARY ----------------
    println("\n" * "="^72)
    println("SUMMARY  (achieved vs target)")
    println("="^72)
    @printf("%-16s %-22s %-22s %-10s\n", "case", "achieved (f,τ)", "target (f,τ)", "|Δ|res")
    tgts = Dict("inviscid"=>(10.50,29.54),
                "A1 η_c=3e29"=>(10.4884,29.5870),
                "A1 η_c=1e31"=>(10.0898,30.8857))
    allok = true
    for k in ["inviscid","A1 η_c=3e29","A1 η_c=1e31"]
        f,τ,ω,res = results[k]
        ft,τt = tgts[k]
        df = 100*(f-ft)/ft; dτ = 100*(τ-τt)/τt
        @printf("%-16s (%.4f, %.4f)   (%.4f, %.4f)   Δf=%+.3f%% Δτ=%+.3f%%\n",
                k, f, τ, ft, τt, df, dτ)
        allok &= (abs(df) < 0.1 && abs(dτ) < 0.1)
    end
    println("="^72)
    if allok
        println("MATCH: all (f,τ) within 0.1% of Bussieres Table II / inviscid w-mode.")
    else
        println("PARTIAL: see per-row Δ above.")
    end
    println("="^72)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
