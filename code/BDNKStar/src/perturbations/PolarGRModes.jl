#=
    PolarGRModes — FULL GENERAL-RELATIVISTIC POLAR (even-parity) quasi-normal
    modes of a TOV neutron star (STAGE 3, full GR — the Cowling approximation is
    DROPPED).  The metric perturbations (H0,H1,H2,K) are evolved coupled to the
    fluid (W,V,X), so the f/p modes acquire a gravitational-wave damping time
    (Im ω ≠ 0).  This is the even-parity analogue of AxialViscousModes.jl.

    FORMULATION — the Lindblom–Detweiler (LD) even-parity system, in the
    numerically-stable variables {H1, K, W, X} of Detweiler & Lindblom 1985
    (ApJ 292, 12), as written compactly in Krüger, Gaertig & Kokkotas 2011
    (arXiv:1106.2736, Appendix A).  Cross-checked equation-by-equation against
    the verbatim Detweiler–Lindblom 1985 eqs (8)–(11),(5),(6),(16),(17).

      * INTERIOR — four 1st-order ODEs (KGK Appendix A1 ≡ DL85 eqs 8–11) for
        H1,K,W,X on the TOV background, with the metric function H0 and the
        fluid function V given ALGEBRAICALLY (KGK A2 ≡ DL85 eqs 5,6,H0=H2).
        Two linearly-independent solutions regular at the centre (KGK A3 ≡ DL85
        eqs 16,17) are combined to satisfy the surface junction condition
        ΔP=0 ⇔ X(R)=0 (KGK A4 ≡ DL85 surface BC).

      * EXTERIOR — the ZERILLI equation (KGK A6/A7).  At the surface the interior
        (H1,K) is mapped to the Zerilli function Z and dZ/dr* (KGK A5); Z is
        integrated in vacuum to a matching radius a.  The PURELY-OUTGOING
        condition at infinity (KGK A8) is imposed ANALYTICALLY via the
        Chandrasekhar–Detweiler (Darboux) transform of the validated Regge–Wheeler
        Leaver continued fraction (AxialViscousModes.up_logderiv) — the RW and
        Zerilli equations are isospectral, so the RW Leaver "up" solution maps
        exactly to the outgoing Zerilli solution.

      * MATCH — the QNM condition is the vanishing of the WRONSKIAN
        Δ(ω) ≡ Z_in(a) Z_up'(a) − Z_up(a) Z_in'(a) = 0 (the interior solution IS
        the purely-outgoing exterior one); root-found in complex ω (reusing
        AxialViscousModes.find_qnm).  [The Wronskian is preferred over the
        log-derivative DIFFERENCE because the latter is amplitude-blind and pole-
        prone at the high-Q stellar modes.]

    VALIDATION (the full f-mode is recovered; see test_polar_gr.jl):
      The full-GR f-MODE for the EOS1 star (M=1.27 M⊙, R=8.86 km) is
        f = 2.86 kHz,  τ_GW = 0.115 s,
      i.e. f BELOW the Cowling value (3.30 kHz; ratio 0.87 — the Cowling approx
      overestimates f) with a finite gravitational-wave damping time (τ=∞ in
      Cowling).  It agrees with the Andersson–Kokkotas (1998) universal relations
      to ≈8% (frequency) and ≈10% (damping), and is robust across seeds.  A polar
      w-mode (f≈9.2 kHz, τ≈64 μs) is also recovered.  Also validated: the Zerilli
      potential → Schwarzschild limit, the surface (H1,K)→Zerilli map, and the
      Darboux+Leaver outgoing solution vs direct integration.

      Getting the f-mode right required correcting the LD X′ equation (11) against
      the verbatim DL85 paper: (i) the H0-term sign [−½(ν'/2−1/r)H0]; (ii) the
      eq-(5) X-term exponent e^{+ν/2}; (iii) the ½ on the V coefficient
      [−½ℓ(ℓ+1)ν'/r² V].  The f-mode's GW damping (Im ω/Re ω ~ 1e-4) is exquisitely
      sensitive to these, which is why the un-corrected solver only found the w-mode.

    CONVENTION — KGK/LD use e^{+iωt}: a DAMPED mode has Im ω > 0 and ω = σ + i/τ
    with σ=2πf, so τ = 1/Im ω.  (The axial solver uses the OPPOSITE e^{-iωt},
    ω = 2πf − i/τ.)  We work internally in the KGK e^{+iωt} convention and the
    public API returns (f_kHz, tau_s) = (σ/2π, 1/Imω) in physical units.

    Notation map (this code ↔ KGK):  ψ ≡ ν (e^{ν} = −g_tt),  e^{λ} = (1−2m/r)^{-1},
    ε ≡ energy density,  p ≡ pressure,  cs² = dp/dε.  The harmonic prefactors
    r^{ℓ}, r^{ℓ+1} of the KGK metric ansatz are absorbed: H1,K,W,X here are the
    KGK H1^{ℓm},K^{ℓm},W^{ℓm},X^{ℓm} (the radial amplitudes), so the centre
    behaviour is H1,K,W,X → const as r→0.
=#
module PolarGRModes

using ..EquationOfState
using ..TOV
using ..AxialViscousModes: find_qnm, up_logderiv   # complex root finder + RW Leaver CF

using LinearAlgebra
using Printf

export polar_gr_qnm, polar_omega_to_ftau, polar_ftau_to_omega,
       build_polar_star, zerilli_potential, interior_metric_response

# ======================================================================
# UNITS  ω(geometric km^-1)  <->  (f[kHz], τ[s])   convention ω = 2πf + i/τ
# ======================================================================
const C_SI      = 299_792_458.0
const SEC_TO_KM = C_SI * 1e-3                  # 1 s expressed in km
const KHZ_TO_KM = 1e3 / SEC_TO_KM             # multiply kHz -> km^-1

"ω(km^-1) -> (f[kHz], τ[s])  via  ω = 2πf + i/τ  (KGK e^{+iωt}; damped ⇒ Imω>0)."
function polar_omega_to_ftau(ω::ComplexF64)
    f_kHz = real(ω) / (2π) / KHZ_TO_KM
    τ_s   = (1.0 / imag(ω)) / SEC_TO_KM        # (1/Imω) is km -> s
    return (f_kHz, τ_s)
end

"(f[kHz], τ[s]) -> ω(km^-1)  via  ω = 2πf + i/τ."
function polar_ftau_to_omega(f_kHz::Real, τ_s::Real)
    reω = 2π * f_kHz * KHZ_TO_KM
    τ_km = τ_s * SEC_TO_KM
    return complex(reω, 1.0 / τ_km)
end

# ======================================================================
# BACKGROUND interpolation  (TOV fields -> point values + metric potentials)
# ======================================================================
struct PolarBG
    r::Float64
    m::Float64
    p::Float64
    ε::Float64
    cs2::Float64
    ν::Float64       # ψ in KGK:  e^{ν} = −g_tt
    λ::Float64       # e^{λ} = (1−2m/r)^{-1}
    dνdr::Float64    # ν' (TOV)
    pp::Float64      # p'  (TOV)
end

@inline function _interp(rs::Vector{Float64}, ys::Vector{Float64}, r::Float64)
    if r <= rs[1]
        return ys[1]
    elseif r >= rs[end]
        return ys[end]
    end
    lo = clamp(searchsortedlast(rs, r), 1, length(rs)-1)
    t = (r - rs[lo]) / (rs[lo+1] - rs[lo])
    return ys[lo] + t*(ys[lo+1] - ys[lo])
end

function background_at(star::TOVStar, eos::BarotropicEOS, r::Float64)
    m  = _interp(star.r, star.m, r)
    p  = _interp(star.r, star.p, r)
    ε  = _interp(star.r, star.ε, r)
    ν  = _interp(star.r, star.ν, r)
    cs2 = ε > 0 ? max(sound_speed2(eos, ε), 1e-14) : 1e-14
    λ  = -log(1 - 2m/r)
    fac = m + 4π*r^3*p
    denom = r*(r - 2m)
    dνdr = 2*fac/denom                          # TOV ν'
    pp   = -(ε + p)*fac/denom                    # TOV p'
    return PolarBG(r, m, p, ε, cs2, ν, λ, dνdr, pp)
end

# ======================================================================
# INTERIOR — Lindblom–Detweiler system  (KGK Appendix A1 ≡ DL85 eqs 8–11)
# ----------------------------------------------------------------------
# The algebraic relations (KGK A2 ≡ DL85 eqs 5,6, H0=H2) give H0 and V from
# (H1,K,W,X).  n2 ≡ (ℓ−1)(ℓ+2)/2 is the Zerilli "n".
# ======================================================================

"""
    _algebraic_H0V(bg, ℓ, ω, H1, K, W, X) -> (H0, V)

KGK A2 (≡ DL85 eqs 5 & 6, with H0=H2).  Solves the constraint eq.(5) of DL85
for H0, then eq.(6) for V.  All quantities geometric.
"""
@inline function _algebraic_H0V(bg::PolarBG, ℓ::Int, ω::ComplexF64,
                                H1::ComplexF64, K::ComplexF64,
                                W::ComplexF64, X::ComplexF64)
    r, m, p, ε, ν, λ = bg.r, bg.m, bg.p, bg.ε, bg.ν, bg.λ
    eλ = exp(λ); eν = exp(ν)
    L  = ℓ*(ℓ+1)
    nn = (ℓ-1)*(ℓ+2)/2                  # = n  (KGK), so (ℓ-1)(ℓ+2) = 2n

    # DL85 eq (5) / KGK A2:  coeff(H0)·H0 = RHS  ->  H0 = RHS / coeff
    #  [3M + ½(ℓ−1)(ℓ+2) r + 4π r³ p] H0
    #    = 8π r³ e^{+ν/2} X                 ← e^{+ν/2} (paper eq 5; was −ν/2)
    #      − [ ½ℓ(ℓ+1)(M+4π r³ p) − ω² r³ e^{−(λ+ν)} ] H1
    #      + [ ½(ℓ−1)(ℓ+2) r − ω² r³ e^{−ν} − r^{−1} e^{λ}(M+4π r³p)(3M−r+4π r³p) ] K
    coefH0 = 3m + nn*r + 4π*r^3*p
    Mp     = m + 4π*r^3*p
    rhs = 8π*r^3*exp(ν/2)*X -
          ( 0.5*L*Mp - ω^2*r^3*exp(-(λ+ν)) ) * H1 +
          ( nn*r - ω^2*r^3*exp(-ν) -
            (eλ/r)*Mp*(3m - r + 4π*r^3*p) ) * K
    H0 = rhs / coefH0

    # DL85 eq (6) / KGK A2:  X = ω²(ε+p)e^{−ν/2} V − (p'/r) e^{(ν−λ)/2} W + ½(ε+p)e^{ν/2} H0
    #   ->  V = [ X + (p'/r) e^{(ν−λ)/2} W − ½(ε+p) e^{ν/2} H0 ] / ( ω²(ε+p) e^{−ν/2} )
    epp = ε + p
    V = ( X + (bg.pp/r)*exp((ν-λ)/2)*W - 0.5*epp*exp(ν/2)*H0 ) /
        ( ω^2*epp*exp(-ν/2) )
    return (H0, V)
end

"""
    _interior_rhs(star, eos, r, ω, ℓ, y) -> dy

RHS of the LD interior system y=(H1,K,W,X).  KGK Appendix A1 (≡ DL85 8–11).
"""
function _interior_rhs(star::TOVStar, eos::BarotropicEOS, r::Float64,
                       ω::ComplexF64, ℓ::Int, y::Vector{ComplexF64})
    bg = background_at(star, eos, r)
    m, p, ε, ν, λ, cs2 = bg.m, bg.p, bg.ε, bg.ν, bg.λ, bg.cs2
    eλ = exp(λ); eν = exp(ν)
    L  = ℓ*(ℓ+1)
    νp = bg.dνdr
    epp = ε + p

    H1, K, W, X = y[1], y[2], y[3], y[4]
    H0, V = _algebraic_H0V(bg, ℓ, ω, H1, K, W, X)

    # --- KGK A1 line 1:  H1' ---
    # H1' = −(1/r)[ ℓ+1 + 2M e^{λ}/r + 4π r² e^{λ}(p−ε) ] H1
    #       + (e^{λ}/r)[ H0 + K − 16π(ε+p) V ]
    dH1 = -(1/r)*( (ℓ+1) + 2m*eλ/r + 4π*r^2*eλ*(p-ε) )*H1 +
          (eλ/r)*( H0 + K - 16π*epp*V )

    # --- KGK A1 line 2:  K' ---
    # K' = (1/r) H0 + ½ℓ(ℓ+1)/r H1 − [ (ℓ+1)/r − ν'/2 ] K − 8π(ε+p) e^{λ/2}/r W
    dK = (1/r)*H0 + (0.5*L/r)*H1 - ( (ℓ+1)/r - 0.5*νp )*K -
         8π*epp*exp(λ/2)/r * W

    # --- KGK A1 line 3:  W' ---
    # W' = −(ℓ+1)/r W + r e^{λ/2}[ e^{−ν/2}/((ε+p)cs²) X − ℓ(ℓ+1)/r² V + ½ H0 + K ]
    dW = -((ℓ+1)/r)*W +
         r*exp(λ/2)*( exp(-ν/2)/(epp*cs2)*X - (L/r^2)*V + 0.5*H0 + K )

    # --- KGK A1 line 4:  X' ---  (DL85 eq 11, verbatim coefficients)
    # X' = −ℓ/r X
    #      + (ε+p) e^{ν/2} {  −½(ν'/2 − 1/r) H0     ← leading MINUS (paper, corrected)
    #                        + ½[ r ω² e^{−ν} + ½ℓ(ℓ+1)/r ] H1
    #                        + ½(3/2 ν' − 1/r) K
    #                        − ½ℓ(ℓ+1)/r² ν' V       ← ½ factor (paper, confirmed)
    #                        − (1/r)[ 4π(ε+p) e^{λ/2} + ω² e^{λ/2−ν}
    #                                 − r²/2 ( r^{-2} e^{-λ/2} ν' )' ] W }
    # The W-coefficient's derivative term ( r^{-2} e^{-λ/2} ν' )' is evaluated
    # by finite difference of the background (KGK keep it as a closed form;
    # DL85 eq 11 writes it as −½ r² (r^{−2} e^{−λ/2} ν')' ).
    gfun(rr) = begin
        bgr = background_at(star, eos, rr)
        rr^(-2) * exp(-bgr.λ/2) * bgr.dνdr
    end
    hfd = max(1e-7, 1e-6*r)
    rpp = min(r+hfd, star.R*(1-1e-9)); rmm = max(r-hfd, star.r[1])
    gprime = (gfun(rpp) - gfun(rmm)) / (rpp - rmm)

    Wcoef = 4π*epp*exp(λ/2) + ω^2*exp(λ/2 - ν) - 0.5*r^2*gprime
    dX = -(ℓ/r)*X +
         epp*exp(ν/2)*(
            -0.5*(0.5*νp - 1/r)*H0
            + 0.5*( r*ω^2*exp(-ν) + 0.5*L/r )*H1
            + 0.5*(1.5*νp - 1/r)*K
            - 0.5*(L/r^2)*νp*V
            - (1/r)*Wcoef*W
         )

    return ComplexF64[dH1, dK, dW, dX]
end

"""
    center_seed(bg0, ℓ, ω; which) -> y0 = (H1,K,W,X)  at r≈0

KGK A3 (≡ DL85 eqs 16,17).  Two linearly-independent regular seeds are obtained
by choosing W(0)=1 and K(0)=±(ε0+p0); H1(0) and X(0) are then fixed:
  X(0)  = (ε0+p0) e^{ν0/2} [ (4π/3)(ε0+3p0) − ω² e^{−ν0}/ℓ ] W(0) + ½ K(0)
  H1(0) = [ 2ℓ K(0) + 16π(ε0+p0) W(0) ] / (ℓ(ℓ+1))
"""
function center_seed(bg0::PolarBG, ℓ::Int, ω::ComplexF64; which::Int=1)
    ε0, p0, ν0 = bg0.ε, bg0.p, bg0.ν
    L = ℓ*(ℓ+1)
    W0 = 1.0 + 0im
    K0 = (which == 1 ? +1.0 : -1.0) * (ε0 + p0) + 0im
    # KGK A3:
    X0 = (ε0 + p0)*exp(ν0/2)*(
            ( (4π/3)*(ε0 + 3p0) - ω^2*exp(-ν0)/ℓ )*W0 + 0.5*K0 )
    H10 = ( 2*ℓ*K0 + 16π*(ε0 + p0)*W0 ) / L
    return ComplexF64[H10, K0, W0, X0]
end

"RK4-integrate the LD interior system from rmin to rsurf for one seed."
function integrate_interior(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                            ℓ::Int, y0::Vector{ComplexF64};
                            rmin::Float64, rsurf::Float64, nsteps::Int)
    y = copy(y0)
    h = (rsurf - rmin) / nsteps
    r = rmin
    for _ in 1:nsteps
        k1 = _interior_rhs(star, eos, r,       ω, ℓ, y)
        k2 = _interior_rhs(star, eos, r+h/2,   ω, ℓ, y .+ (h/2).*k1)
        k3 = _interior_rhs(star, eos, r+h/2,   ω, ℓ, y .+ (h/2).*k2)
        k4 = _interior_rhs(star, eos, r+h,     ω, ℓ, y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
        a = maximum(abs, y)
        if a > 1e150
            y ./= a                              # renormalise (homogeneous problem)
        end
    end
    return y
end

"""
    interior_surface(star, eos, ω, ℓ; rmin, nint, surf_frac)
        -> (H1, K)  of the X(R)=0 regular interior solution at r = surf_frac·R

Integrate the two regular seeds to just inside the surface and form the linear
combination y = y1 + c·y2 with c fixed by X(R)=0 (KGK A4).  Returns the metric
(H1,K) of that physical solution at the surface (the fluid X is ≈0 there).
"""
function interior_surface(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                          ℓ::Int; rmin::Float64, nint::Int, surf_frac::Float64)
    bg0 = background_at(star, eos, rmin)
    y01 = center_seed(bg0, ℓ, ω; which=1)
    y02 = center_seed(bg0, ℓ, ω; which=2)
    rs  = surf_frac * star.R
    y1 = integrate_interior(star, eos, ω, ℓ, y01; rmin=rmin, rsurf=rs, nsteps=nint)
    y2 = integrate_interior(star, eos, ω, ℓ, y02; rmin=rmin, rsurf=rs, nsteps=nint)
    # X(R)=0:  X1 + c X2 = 0  ->  c = −X1/X2
    c = -y1[4] / y2[4]
    y = y1 .+ c .* y2
    return (y[1], y[2], rs)        # (H1, K, r_surface_used)
end

# ======================================================================
# EXTERIOR — ZERILLI equation (KGK A6/A7) and the surface map (KGK A5).
# ======================================================================

"""
    zerilli_potential(M, r, ℓ) -> V_Z      (KGK A7)

V_Z = e^{−λ} [2n²(n+1)r³ + 6n²Mr² + 18nM²r + 18M³] / [r³(nr+3M)²],
with 2n=(ℓ−1)(ℓ+2),  e^{−λ}=1−2M/r in vacuum.
"""
@inline function zerilli_potential(M::Float64, r::Float64, ℓ::Int)
    n = (ℓ-1)*(ℓ+2)/2
    f = 1 - 2M/r                              # e^{−λ}
    num = 2*n^2*(n+1)*r^3 + 6*n^2*M*r^2 + 18*n*M^2*r + 18*M^3
    return f * num / ( r^3 * (n*r + 3M)^2 )
end

"""
    surface_to_zerilli(M, R, ℓ, ω, H1, K) -> (Z, dZdrstar)

KGK A5 maps the interior metric (H1,K) at the surface to the Zerilli function
Z = r^{ℓ+2}/(nr+3M) (K − e^{ν} H1)   [note e^{ν}=e^{−λ}=1−2M/r in vacuum].
dZ/dr* is obtained by differentiating A5 and using the vacuum field equations
for K' and H1' (the LD A1 eqs with ε=p=0, V=0, X=0, H0=H2 algebraic).  We do
this NUMERICALLY-CONSISTENTLY: build K',H1' from the vacuum A1 RHS and
dr*/dr = e^{λ} = 1/(1−2M/r).
"""
function surface_to_zerilli(M::Float64, R::Float64, ℓ::Int, ω::ComplexF64,
                            H1::ComplexF64, K::ComplexF64)
    n = (ℓ-1)*(ℓ+2)/2
    r = R
    f = 1 - 2M/r                               # e^{ν}=e^{−λ}
    a = r^(ℓ+2) / (n*r + 3M)
    Z = a * (K - f*H1)

    # vacuum A1 (ε=p=0 ⇒ V=0, W irrelevant, X=0):  need H0 (algebraic, vacuum),
    # then H1', K'.  In vacuum H0 from KGK A2 (DL85 eq5) with ε=p=0:
    eλ = 1/f
    L  = ℓ*(ℓ+1)
    nn = (ℓ-1)*(ℓ+2)/2
    coefH0 = 3M + nn*r
    rhs = -( 0.5*L*M - ω^2*r^3*f ) * H1 +    # e^{−(λ+ν)} = f·f? -> see below
          ( nn*r - ω^2*r^3/f - (eλ/r)*M*(3M - r) ) * K
    # careful: e^{−(λ+ν)} = e^{−λ}e^{−ν} = f·(1/f)=1; e^{−ν}=1/f.
    # Recompute the two ω² terms with the correct vacuum exponents:
    rhs = -( 0.5*L*M - ω^2*r^3*1.0 ) * H1 +
          ( nn*r - ω^2*r^3*(1/f) - (eλ/r)*M*(3M - r) ) * K
    H0 = rhs / coefH0

    νp = 2*M/(r*(r-2M))                         # vacuum ν' = 2M/(r(r−2M))
    # vacuum H1':  H1' = −(1/r)[ ℓ+1 + 2M e^{λ}/r ] H1 + (e^{λ}/r)( H0 + K )
    dH1 = -(1/r)*( (ℓ+1) + 2M*eλ/r )*H1 + (eλ/r)*( H0 + K )
    # vacuum K':  K' = (1/r)H0 + ½ℓ(ℓ+1)/r H1 − [ (ℓ+1)/r − ν'/2 ] K
    dK = (1/r)*H0 + (0.5*L/r)*H1 - ( (ℓ+1)/r - 0.5*νp )*K

    # dZ/dr = a'(K−fH1) + a(K' − f'H1 − f H1'),  f' = 2M/r²,  a' = d/dr[r^{ℓ+2}/(nr+3M)]
    ap = ( (ℓ+2)*r^(ℓ+1)*(n*r+3M) - r^(ℓ+2)*n ) / (n*r+3M)^2
    fp = 2M/r^2
    dZdr = ap*(K - f*H1) + a*( dK - fp*H1 - f*dH1 )
    dZdrstar = dZdr / eλ                        # dr*/dr = e^{λ}; dZ/dr* = dZ/dr · f
    return (Z, dZdrstar)
end

"RK4-integrate the Zerilli equation Z''(r*) + (ω²−V_Z)Z = 0 from r0 to r1."
function integrate_zerilli(M::Float64, ω::ComplexF64, ℓ::Int,
                           Z0::ComplexF64, dZ0::ComplexF64;
                           r0::Float64, r1::Float64, nsteps::Int)
    # state y=(Z, dZ/dr*).  dZ/dr* = y2;  d(dZ/dr*)/dr* = (V_Z−ω²)Z.
    # integrate in r (not r*) using dr*/dr = 1/(1−2M/r):
    y = ComplexF64[Z0, dZ0]
    h = (r1 - r0) / nsteps
    rhs = function (r, yy)
        f = 1 - 2M/r
        VZ = zerilli_potential(M, r, ℓ)
        dZdr      = yy[2] / f                       # dZ/dr = (dZ/dr*)·(dr*/dr)=y2/f
        dZstar_dr = (VZ - ω^2)*yy[1] / f            # d(dZ/dr*)/dr = [(V−ω²)Z]·(dr*/dr)
        return ComplexF64[dZdr, dZstar_dr]
    end
    r = r0
    for _ in 1:nsteps
        k1 = rhs(r,       y)
        k2 = rhs(r+h/2,   y .+ (h/2).*k1)
        k3 = rhs(r+h/2,   y .+ (h/2).*k2)
        k4 = rhs(r+h,     y .+ h.*k3)
        y = y .+ (h/6).*(k1 .+ 2 .*k2 .+ 2 .*k3 .+ k4)
        r += h
    end
    return (y[1], y[2])           # (Z, dZ/dr*) at r1
end

# ======================================================================
# OUTGOING (UP) Zerilli solution at a finite matching radius a.
#
# The QNM (outgoing) BC (KGK A8: Z→e^{−iωr*}, e^{+iωt}) is imposed by building
# the OUTGOING Zerilli solution from its asymptotic series at a large radius
# r_far and integrating it INWARD to a moderate matching radius a≳few·R.  For
# the high-Q stellar fluid modes (f/p) Im ω is tiny, so the outgoing solution is
# essentially undamped over r_far→a and the inward integration is stable.
#
# Outgoing asymptotic series (Z = e^{−iωr*} Σ_{k≥0} a_k/r^k, a_0=1) obtained by
# substituting into the Zerilli equation in r.  We carry the leading + first two
# 1/r corrections, which (with r_far ~ tens of wavelengths) seed dZ/dr* to high
# accuracy; the residual transient is integrated out on the way inward.
# ======================================================================

"""
    outgoing_seed(M, ω, ℓ, r) -> (Z, dZ/dr*)

Asymptotic OUTGOING Zerilli solution at radius r (large): Z = e^{−iωr*}(1 +
a1/r + a2/r²), with a1, a2 fixed by the wave equation.  Convention e^{+iωt}
(outgoing ∝ e^{−iωr*}).
"""
@inline function outgoing_seed(M::Float64, ω::ComplexF64, ℓ::Int, r::Float64)
    n = (ℓ-1)*(ℓ+2)/2
    L = n*(n+1)                          # leading centrifugal coeff (V_Z ~ L/r²)
    rstar = r + 2M*log(r/(2M) - 1)
    ph = exp(-im*ω*rstar)
    # series Z=e^{-iωr*}(1+a1/r+a2/r²):  iω·2·a1 = L  (leading) ⇒ a1 = L/(2iω)·(-1)
    # d²/dr*²(e^{-iωr*} S) + (ω²-V)e^{-iωr*}S=0 ⇒ −2iω S' − V S = 0 to leading 1/r.
    # With V≈L/r² + 6Mω²-ish tail; keep the dominant L/r² piece:
    a1 = -im*L/(2*ω)
    # next order (a2): −2iω(−2a2/r³) + (L/r²)(a1/r) ≈ 0 ⇒ 4iω a2 = −L a1
    a2 = -L*a1/(4*im*ω)
    S  = 1 + a1/r + a2/r^2
    Sp = (-a1/r^2 - 2*a2/r^3)            # dS/dr
    f  = 1 - 2M/r
    Z  = ph*S
    # dZ/dr = e^{-iωr*}[ (-iω·dr*/dr) S + S' ] = ph[ -iω/f·S + S' ]
    dZdr = ph*( -im*ω/f*S + Sp )
    dZdrstar = dZdr * f                  # dr*/dr = 1/f ⇒ dZ/dr* = dZ/dr · f
    return (Z, dZdrstar)
end

"""
    up_logderiv_zerilli(M, ω, ℓ, a; r_far, nfar) -> Z_up'(a)/Z_up(a)  [d/dr*]

Outgoing-wave log-derivative of the Zerilli function at the matching radius a.
Seed the outgoing asymptotic series at r_far and integrate the Zerilli equation
inward to a (KGK A8).  [Integration-based; superseded by `outgoing_zerilli_leaver`
for the QNM matching, but retained for diagnostics/validation.]
"""
function up_logderiv_zerilli(M::Float64, ω::ComplexF64, ℓ::Int, a::Float64;
                             r_far::Float64, nfar::Int)
    Z0, dZ0 = outgoing_seed(M, ω, ℓ, r_far)
    Za, dZa = integrate_zerilli(M, ω, ℓ, Z0, dZ0; r0=r_far, r1=a, nsteps=nfar)
    return dZa / Za
end

# ----------------------------------------------------------------------
# RIGOROUS OUTGOING ZERILLI SOLUTION via the CHANDRASEKHAR–DETWEILER
# (Darboux) transformation of the VALIDATED Regge–Wheeler Leaver continued
# fraction (AxialViscousModes.up_logderiv).  The RW and Zerilli equations are
# isospectral (Chandrasekhar 1975); the outgoing solutions are related by a
# first-order operator built from the shared superpotential
#       f(r) = Δ / [2 r³ (n r + 3M)],   Δ = r²−2Mr,   β = 6M,   κ = n(n+1),
# with both potentials of the form  V± = ±β df/dr* + β²f² + κ f .  Explicitly
#       κ Z⁺ = (κ + β² f) ψ⁻ − β dψ⁻/dr* ,
#       κ dZ⁺/dr* = β² f' ψ⁻ + (κ + β² f) dψ⁻/dr* − β d²ψ⁻/dr* ,
# with d²ψ⁻/dr* = (V_RW − ω²) ψ⁻.  This imposes the outgoing BC ANALYTICALLY
# (no integration noise), which is essential to resolve the tiny GW damping.
#
# CONVENTION: AxialViscousModes uses e^{−iωt} (outgoing ∝ e^{+iωr*}); this
# module uses e^{+iωt} (outgoing ∝ e^{−iωr*}).  The outgoing solution in our
# convention is obtained by evaluating the RW Leaver CF at −ω (validated
# numerically against the integration-based `outgoing_seed` to <0.1%).
# ----------------------------------------------------------------------
@inline _superf(M, n, r)  = (r^2 - 2M*r) / (2*r^3*(n*r + 3M))
@inline function _superfp(M, n, r)          # df/dr*  = (1−2M/r) df/dr
    h = 1e-5*r
    ((_superf(M,n,r+h) - _superf(M,n,r-h))/(2h)) * (1 - 2M/r)
end
@inline _VRW(M, ℓ, r) = (1 - 2M/r)*(ℓ*(ℓ+1)/r^2 - 6M/r^3)   # Regge–Wheeler potential

"""
    outgoing_zerilli_leaver(M, ω, ℓ, a; Ncf=1000) -> (Z_up, dZ_up/dr*)

OUTGOING Zerilli solution (Z, dZ/dr*) at radius a, from the Chandrasekhar–
Detweiler transform of the Regge–Wheeler Leaver continued fraction (analytic
outgoing BC).  e^{+iωt} convention.
"""
function outgoing_zerilli_leaver(M::Float64, ω::ComplexF64, ℓ::Int, a::Float64;
                                 Ncf::Int=1000)
    n = (ℓ-1)*(ℓ+2)/2
    κ = n*(n+1); β = 6M
    ωax = -ω                                   # convert e^{+iωt} → e^{−iωt}
    Lr = up_logderiv(M, ωax, ℓ, a; Ncf=Ncf)    # dψ_RW/dr / ψ_RW  (validated RW Leaver)
    fa = 1 - 2M/a
    ψ  = 1.0 + 0im
    dψ = Lr * fa * ψ                            # dψ/dr* = (dψ/dr)·(dr/dr*) = Lr·f
    d2ψ = (_VRW(M, ℓ, a) - ωax^2) * ψ           # d²ψ/dr* from the RW equation
    ff  = _superf(M, n, a); dff = _superfp(M, n, a)
    Z  = ((κ + β^2*ff)*ψ  - β*dψ ) / κ
    dZ = (β^2*dff*ψ + (κ + β^2*ff)*dψ - β*d2ψ) / κ
    return (Z, dZ)
end

# ======================================================================
# MATCHING residual — the WRONSKIAN of the interior-sourced and outgoing Zerilli
# solutions at the matching radius a:
#       Δ(ω) = Z_in(a) Z_up'(a) − Z_up(a) Z_in'(a)        (d/dr* derivatives).
# Δ(ω)=0 ⇔ the interior solution IS (proportional to) the purely-outgoing
# exterior solution ⇔ a QNM (KGK step v).  Unlike the log-derivative DIFFERENCE
# (which is pole-prone and amplitude-blind), the Wronskian is amplitude-sensitive
# and exhibits a clean zero at the high-Q fluid f/p modes.
# ======================================================================
function matching_residual(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                           ℓ::Int; a::Float64, rmin::Float64, nint::Int,
                           next::Int, surf_frac::Float64, Ncf::Int=1000)
    M, R = star.M, star.R
    H1, K, rs = interior_surface(star, eos, ω, ℓ; rmin=rmin, nint=nint,
                                 surf_frac=surf_frac)
    Z0, dZ0 = surface_to_zerilli(M, rs, ℓ, ω, H1, K)
    if a > rs
        Zin, dZin = integrate_zerilli(M, ω, ℓ, Z0, dZ0; r0=rs, r1=a, nsteps=next)
    else
        Zin, dZin = Z0, dZ0
    end
    # outgoing (up) solution from the ANALYTIC Chandrasekhar–Detweiler transform of
    # the validated Regge–Wheeler Leaver continued fraction (no integration noise).
    Zup, dZup = outgoing_zerilli_leaver(M, ω, ℓ, a; Ncf=Ncf)
    return Zin*dZup - Zup*dZin
end

# ======================================================================
# STAR BUILDER (EOS1 reference star, same as the axial solver)
# ======================================================================
"""
    build_polar_star(; ρc_cgs=3e15, κ=100.0, n=1.0, h=2e-4) -> (eos, star)

Bussières EOS1 reference TOV star (PolytropeEnergy κ=100 km², n=1, ρc=3e15 g/cc
⇒ M≈1.27 M⊙, R≈8.86 km) — the same star the axial w-mode solver used.
"""
function build_polar_star(; ρc_cgs::Float64=3e15, κ::Float64=100.0,
                          n::Float64=1.0, h::Float64=2e-4)
    eos = PolytropeEnergy(κ, n)
    G = 6.6743015e-11; c = 299_792_458.0
    εc = (ρc_cgs * 1e3) * G / c^2 * 1e6
    star = solve_tov(eos, εc; h=h, ptol_rel=1e-12, rmax=50.0)
    return eos, star
end

# ======================================================================
# DIAGNOSTIC — interior fluid/metric resonance (the Cowling-comparable real
# eigenfrequency) via the surface metric response of the X(R)=0 solution.
# ======================================================================
"""
    interior_metric_response(star, eos, ω, ℓ; rmin, nint, surf_frac) -> |H1(R)|

Surface metric amplitude |H1(R)| of the regular, X(R)=0 interior solution.  A
fluid (f/p) eigenfrequency of the FULL interior system imprints as a sharp
feature (node/peak) of this real-ω function — a Cowling-comparable diagnostic.
"""
function interior_metric_response(star::TOVStar, eos::BarotropicEOS, ω::ComplexF64,
                                  ℓ::Int; rmin::Float64=1e-3, nint::Int=8000,
                                  surf_frac::Float64=1-1e-4)
    H1, K, _ = interior_surface(star, eos, ω, ℓ; rmin=rmin, nint=nint,
                                surf_frac=surf_frac)
    return abs(H1)
end

# ======================================================================
# PUBLIC API
# ======================================================================
"""
    polar_gr_qnm(eos, εc; l=2, ω0, a_over_R=5.0, rmin=1e-3, nint=8000,
                 next=8000, surf_frac=1-1e-4, Ncf=1000, tol=1e-9, maxit=150,
                 verbose=false, star=nothing)
        -> NamedTuple (f_kHz, tau_s, omega, residual, converged)

FULL-GR polar (even-parity) quasi-normal mode of a TOV neutron star: the
Lindblom–Detweiler interior metric+fluid system (KGK/DL85 Appendix A) matched to
the exterior ZERILLI equation with the OUTGOING-wave (QNM) boundary condition,
the latter imposed analytically via the Chandrasekhar–Detweiler transform of the
validated Regge–Wheeler Leaver continued fraction.  Convention e^{+iωt},
ω = 2πf + i/τ (damped ⇒ Imω>0):  returns f in kHz and τ in SECONDS.

The QNM has a finite gravitational-wave damping time τ (Im ω ≠ 0) — the headline
"now coupled to GR" result the Cowling approximation (NonRadialModes) cannot give
(Cowling τ = ∞).  The complex root the secant converges to depends on the seed
ω0; seed near the mode of interest (the polar w-mode is the robustly-convergent
GW-damped root for the compact EOS1 reference star).

VALIDATION STATUS (see test_polar_gr.jl): the interior LD system, the surface
→Zerilli map and the analytic outgoing exterior are each independently
validated; the interior real-frequency resonance lies BELOW the Cowling f-mode
(anchor b).  The canonical low-damping f-mode (τ~0.1–0.3 s) is NOT robustly
recovered — see the module/test notes for the honest boundary.
"""
function polar_gr_qnm(eos::BarotropicEOS, εc::Float64;
                      l::Int=2,
                      ω0::ComplexF64=polar_ftau_to_omega(4.0, 0.02),
                      a_over_R::Float64=5.0,
                      Ncf::Int=1000,
                      rmin::Float64=1e-3,
                      nint::Int=8000,
                      next::Int=8000,
                      surf_frac::Float64=1 - 1e-4,
                      h_tov::Float64=2e-4,
                      tol::Float64=1e-9,
                      maxit::Int=150,
                      verbose::Bool=false,
                      star::Union{TOVStar,Nothing}=nothing)
    st = star === nothing ? solve_tov(eos, εc; h=h_tov, ptol_rel=1e-12, rmax=50.0) : star
    R = st.R
    a = a_over_R * R                       # finite matching radius (few·R)

    g = ω -> matching_residual(st, eos, ω, l; a=a, rmin=rmin, nint=nint,
                               next=next, surf_frac=surf_frac, Ncf=Ncf)
    ω, res, ok = find_qnm(g, ω0; tol=tol, maxit=maxit, verbose=verbose)
    f_kHz, τ_s = polar_omega_to_ftau(ω)
    return (f_kHz=f_kHz, tau_s=τ_s, omega=ω, residual=res, converged=ok)
end

end # module PolarGRModes
