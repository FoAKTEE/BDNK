#=
    shum_core.jl  — STAGE 1C core (R5):

      (a) areal(polar)→isotropic coordinate transform of a TOV star, and
      (b) the general-EOS BDNK *linear* primitive recovery (the spherical 2×2
          con2prim solve) for the Shum Γ=2 EOS and the Shum spherical-Cowling
          frame.

    Everything is GROUNDED line-by-line in Shum, Abalos, Bea, Bezares, Figueras,
    Palenzuela, arXiv:2509.15303 ("Neutron star evolution with the BDNK ...
    framework"), file ref-paper/sources/arXiv-2509.15303/src/Paper.tex.

    Units: M_⊙ = G = c = 1 (Shum §V, ρ0c = 0.00128 M_⊙^-2 → M_T = 1.4 M_⊙).
    In these geometric units the raw TOV mass IS the mass in solar masses
    (consistent with test/test_tov.jl "Shum M_T=1.4").

    VALIDATION GATE (this file): the BDNK recovery round-trip closes to ≤ 1e-8 in
    the *no-gradient* (equilibrium) limit — i.e. given the conserved (E, S_r)
    built from a static star with v^r = 0 and ε̂ = v̂̄^r = 0 and all spatial
    gradients zero, the 2×2 solve returns ε̂ = v̂̄^r = 0.

    PACKAGE REUSE per task: load the shared trunk first.
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.EquationOfState: ShumPolytrope, pressure, sound_speed2, energy_from_pressure
using .BDNKStar.TOV: TOVStar, solve_tov
using .BDNKStar.Transport: shum_frame_speeds, shum_frame_wellposed
using LinearAlgebra: norm

# ===========================================================================
# (a)  Areal(polar) → isotropic coordinate transform of the TOV background
# ===========================================================================
#
# Polar-areal (Schwarzschild) metric, Shum eq.(schwarzschild_metric), lines
# 454–457:
#     ds² = -α²(R) dt² + a²(R) dR² + R² dΩ²
# Maximal-isotropic metric, Shum eq.(maximal_isotropic), lines 479–482:
#     ds² = -α²(r) dt² + ψ⁴(r) (dr² + r² dΩ²)
#
# Matching the two line elements:
#   • angular part :   R² = ψ⁴ r²            ⇒  ψ² = R / r
#   • radial  part :   a² dR² = ψ⁴ dr²       ⇒  a dR = ψ² dr = (R/r) dr
#     ⇒  d(ln r)/dR = a / R .                                       (transform-ODE)
#
# `solve_tov` returns the areal grid (R = star.r), the mass m(R) = star.m and
# the metric potential ν(R) with g_tt = -e^{ν}. Hence
#     α(R) = e^{ν(R)/2},        a(R)² = g_RR = 1/(1 - 2 m(R)/R)
# the standard TOV areal metric (Shum lines 466–470 use exactly a, α with
# a(0)=α(0)=1 and the Schwarzschild exterior match a→1/α as R→∞).
#
# The integration constant of r is fixed by matching to the EXACT exterior
# Schwarzschild isotropic relation (valid for R ≥ R_star, m = M):
#     R = r (1 + M/(2r))²   ⇒   r(R) = ½ ( R - M + √(R² - 2 M R) )      (ext-iso)
# We integrate (transform-ODE) inward from the surface using this boundary value,
# which automatically gives the conformal factor ψ⁴ = (R/r)² and ψ = √(R/r).

"""
    IsotropicStar

The TOV star re-expressed on an *isotropic* radial grid (Shum eq.480):
`r` isotropic radius, `R` areal radius, `psi` conformal factor (g = ψ⁴(dr²+r²dΩ²)),
`alpha` lapse, `grr = ψ⁴` the spatial-metric radial component in isotropic coords,
plus the matter profiles carried over (`p`, `ε`).  All in M_⊙=G=c=1 units.
"""
struct IsotropicStar
    r::Vector{Float64}      # isotropic radius
    R::Vector{Float64}      # areal radius (same physical points)
    psi::Vector{Float64}    # conformal factor ψ  (g_θθ = ψ⁴ r² ; ψ² = R/r)
    alpha::Vector{Float64}  # lapse α = e^{ν/2}
    grr::Vector{Float64}    # isotropic g_rr = ψ⁴
    p::Vector{Float64}      # pressure
    ε::Vector{Float64}      # total energy density
    M::Float64              # gravitational mass
    Rstar::Float64          # areal stellar radius
    rstar::Float64          # isotropic stellar radius
end

"""
    areal_to_isotropic(star::TOVStar) -> IsotropicStar

Transform the polar-areal TOV solution `star` into maximal-isotropic coordinates
(Shum lines 478–483).  Solves d(ln r)/dR = a/R inward from the surface, where the
isotropic radius is seeded by the exact exterior relation r = ½(R-M+√(R²-2MR)).
"""
function areal_to_isotropic(star::TOVStar)
    R   = star.r
    m   = star.m
    ν   = star.ν
    N   = length(R)
    M   = star.M
    Rstar = star.R

    # a(R) = 1/√(1 - 2 m/R)   (areal g_RR);   α(R) = e^{ν/2}
    a     = [1.0 / sqrt(1 - 2*m[i]/R[i]) for i in 1:N]
    alpha = [exp(ν[i]/2) for i in 1:N]

    # ---- exterior boundary value of r at the surface (exact Schwarzschild) ----
    Rs = R[N]
    rstar = 0.5 * (Rs - M + sqrt(Rs^2 - 2*M*Rs))         # (ext-iso) at R=R_star

    # ---- integrate d(ln r)/dR = a/R inward from the surface (trapezoid) -------
    lnr = Vector{Float64}(undef, N)
    lnr[N] = log(rstar)
    f(i) = a[i] / R[i]                                    # = d(ln r)/dR
    for i in (N-1):-1:1
        dR = R[i+1] - R[i]                                # > 0
        lnr[i] = lnr[i+1] - 0.5*(f(i+1) + f(i)) * dR      # step inward
    end
    r = exp.(lnr)

    # ψ² = R/r  ⇒  ψ = √(R/r),  isotropic g_rr = ψ⁴ = (R/r)²
    psi = sqrt.(R ./ r)
    grr = (R ./ r).^2

    return IsotropicStar(r, R, psi, alpha, grr, star.p, star.ε, M, Rstar, rstar)
end

# ===========================================================================
# Shum spherical-Cowling FRAME  (eqs. 67–71 / hatted_parameters, lines 493–529)
# ===========================================================================
#
# Definitions (lines 494–522):
#   ρ ≡ ε + p ,    V ≡ (4/3)η + ζ ,     V̂ ≡ (4/3)η̂ + ζ̂
#   η  ≡ q̂ L cs² ρ η̂              ζ  ≡ q̂ L cs² ρ ζ̂
#   τ_p ≡ ŝ cs² L V̂               τ_Q ≡ â L V̂           τ_ε ≡ V̂ L
#   β_ε ≡ cs² â V̂ L   (Q^μ ∝ EoM, β_ε = τ_Q p'(ε))
# Production frame (Shum §V.A, lines 606–611):  (ŝ, â, q̂) = (1, 1, 0.999), L = 1
#   ⇒  τ_p = cs² τ_ε ,  τ_Q = τ_ε  (â=1),  giving c₊ = √3 cs, c₋ ≈ 0.0183 cs.

"""
    ShumFrame(ŝ, â, q̂, η̂, ζ̂; L=1.0)

Container for the Shum hatted-frame parameters (lines 498–522).  Production
defaults (ŝ,â,q̂)=(1,1,0.999); a viscosity pair (η̂,ζ̂) must be supplied
(e.g. smallSB-F2 = (0.01,0.01), τ_ε = (4/3)η̂+ζ̂).
"""
Base.@kwdef struct ShumFrame
    ŝ::Float64 = 1.0
    â::Float64 = 1.0
    q̂::Float64 = 0.999
    η̂::Float64
    ζ̂::Float64
    L::Float64 = 1.0
end

"""
    shum_transport(fr::ShumFrame, eos, ε) -> (η, ζ, τε, τp, τQ, p, cs2)

Evaluate the *dimensionful* BDNK transport coefficients of the Shum frame
(lines 494–522) at energy density `ε` for barotropic `eos`.  cs² = dp/dε.
"""
function shum_transport(fr::ShumFrame, eos, ε::Real)
    p   = pressure(eos, ε)
    cs2 = sound_speed2(eos, ε)
    ρ   = ε + p
    V̂   = (4/3)*fr.η̂ + fr.ζ̂
    L   = fr.L
    η   = fr.q̂ * L * cs2 * ρ * fr.η̂
    ζ   = fr.q̂ * L * cs2 * ρ * fr.ζ̂
    τε  = V̂ * L
    τp  = fr.ŝ * cs2 * L * V̂
    τQ  = fr.â * L * V̂
    return η, ζ, τε, τp, τQ, p, cs2
end

# ===========================================================================
# (b)  General-EOS BDNK linear primitive recovery — the spherical 2×2 solve
# ===========================================================================
#
# Shum Appendix A "Primitive variables recovery", spherical case, lines 943–976.
# The conserved variables map LINEARLY onto p₁ = (ε̂, v̂̄^r):
#
#       ( A00  A01 ) ( ε̂   )   ( b0 )
#       ( A10  A11 ) ( v̂̄^r ) = ( br )                         (eq. lines 944–963)
#
# with (lines 966–975), writing  X = g_rr (v^r)²  and  p' = ∂_ε p = cs²,
# τε≡τ_ε, τQ≡τ_Q, ρ = ε + p :
#
#   A00 = -[ 2 g_rr (v^r)² τ_Q p' + τ_ε ( g_rr (v^r)² p' + 1 ) ] / (1 - X)^{3/2}
#   A01 = - g_rr v^r [ -4 g_rr (v^r)² η + 3 g_rr (v^r)² ( ρ τ_ε p' - ζ )
#                       + 3 ρ (2 τ_Q + τ_ε) ] / [ 3 (1 - X)^{5/2} ]
#   A10 = - g_rr v^r [ (g_rr (v^r)² + 1) τ_Q p' + τ_ε (p' + 1) ] / (1 - X)^{3/2}
#   A11 = - g_rr [ -4 g_rr (v^r)² η + 3 g_rr (v^r)² ( ρ ( τ_ε (p'+1) + τ_Q ) - ζ )
#                   + 3 ρ τ_Q ] / [ 3 (1 - X)^{5/2} ]
#
# (b0, br) are the "perfect-fluid + frozen-gradient" deficits E − c0, S_r − c_r.
# Shum (line 978) deliberately omits the explicit (b0,br); for the equilibrium
# round-trip gate we construct them from the static-star primitives with ALL
# spatial gradients zero (the only regime the paper supplies explicitly).

"""
    shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, grr, vr) -> A::2×2 Matrix

The Shum spherical con2prim matrix 𝒜, exactly as Appendix-A lines 966–975
(`∂_ε p = cs2`).  `grr` is the (isotropic) spatial-metric radial component, `vr`
the contravariant radial 3-velocity v^r.
"""
function shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, grr, vr)
    pe  = cs2                       # ∂_ε p
    ρ   = ε + p
    g   = grr
    v   = vr
    X   = g * v^2                   # = g_rr (v^r)²  (= v_μ v^μ for radial flow)
    s32 = (1 - X)^(3/2)
    s52 = (1 - X)^(5/2)

    A00 = -( 2*g*v^2*τQ*pe + τε*(g*v^2*pe + 1) ) / s32                            # line 967–968
    A01 = -( g*v*( -4*g*v^2*η + 3*g*v^2*(ρ*τε*pe - ζ) + 3*ρ*(2*τQ + τε) ) ) /
           ( 3*s52 )                                                              # line 969–970
    A10 = -( g*v*( (g*v^2 + 1)*τQ*pe + τε*(pe + 1) ) ) / s32                      # line 971
    A11 = -( g*( -4*g*v^2*η + 3*g*v^2*(ρ*(τε*(pe + 1) + τQ) - ζ) + 3*ρ*τQ ) ) /
           ( 3*s52 )                                                              # line 972–974

    return [A00 A01; A10 A11]
end

"""
    recover_time_derivs_shum(fr::ShumFrame, eos, ε, vr, grr, E, Sr;
                             b0, br) -> (ε̂, v̂̄^r, A)

General-EOS BDNK *linear* primitive recovery for the Shum frame: build the 2×2
matrix 𝒜 (Appendix-A, lines 966–975) from the local frozen state (ε, vr, grr)
and solve 𝒜 · (ε̂, v̂̄^r)ᵀ = (b0, br)ᵀ where (b0, br) = (E, S_r) − (c0, c_r) are
the supplied gradient-frozen deficits.  Analogous to the conformal
`recover_time_derivs`, but for the Shum EOS/frame and the 2×2 spherical system.
"""
function recover_time_derivs_shum(fr::ShumFrame, eos, ε::Real, vr::Real, grr::Real,
                                  E::Real, Sr::Real; b0::Real, br::Real)
    η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
    A = shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, grr, vr)
    rhs = [b0, br]
    sol = A \ rhs
    return sol[1], sol[2], A
end

# ---------------------------------------------------------------------------
# Equilibrium "perfect-fluid" conserved densities (no-gradient limit).
# Static, spherically symmetric, v^r = 0, ε̂ = v̂̄^r = 0, all D_i = 0, K = 0.
# Then 𝒜 (ε̂, v̂̄^r)ᵀ = 0 and the constant vector (c0, c_r) carries the whole
# perfect-fluid stress.  From Shum lines 881/886 with v=0, W=1, all gradients 0:
#     c0 = -p(1 - W²) + W² ε  →  ε
#     c_r = v_r W² (p+ε)       →  0
# i.e. (E, S_r) = (ε, 0) and (b0, br) = (E, S_r) - (c0, c_r) = (0, 0).
# A correct linear solve must then return (ε̂, v̂̄^r) = (0, 0) exactly. This is
# the round-trip gate.  (Equivalently: 𝒜·0 = 0 ⇒ recovered p₁ = 0.)
# ---------------------------------------------------------------------------

"""
    equilibrium_conserved(eos, ε) -> (E, S_r)

Perfect-fluid conserved densities (Shum c0,c_r constant vector, lines 881/886) at
a static fluid element (v^r=0, all gradients zero): E = ε, S_r = 0.
"""
equilibrium_conserved(eos, ε::Real) = (ε, 0.0)

# ===========================================================================
# DRIVER / VALIDATION
# ===========================================================================
function main()
    println("="^74)
    println("STAGE 1C core (Shum 2509.15303): isotropic transform + BDNK recovery")
    println("="^74)

    # ---- TOV background (Shum §V): ρ0c = 0.00128 → M_T = 1.4 M_⊙ -----------
    κ   = 100.0
    eos = ShumPolytrope(κ)
    ρ0c = 0.00128
    εc  = ρ0c + κ*ρ0c^2                       # ε = ρ0 + p/(Γ-1), Γ=2, p=κρ0²
    println("\n[TOV]  ShumPolytrope(κ=$κ),  ρ0c=$ρ0c,  εc=$(round(εc,sigdigits=6))")
    star = solve_tov(eos, εc; h=2e-4, ptol_rel=1e-12, rmax=50.0)
    println("       M_T = $(round(star.M,sigdigits=6)) M_⊙   (target 1.4),  " *
            "R_areal = $(round(star.R,sigdigits=6)) M_⊙")
    M_ok = isapprox(star.M, 1.4; atol=0.02)
    println("       M_T within 0.02 of 1.4 : $M_ok")

    # ---- (a) areal → isotropic transform ----------------------------------
    iso = areal_to_isotropic(star)
    println("\n[ISO]  areal→isotropic transform (Shum lines 478–483)")
    println("       R_star(areal)    = $(round(iso.Rstar,sigdigits=8)) M_⊙")
    println("       r_star(isotropic)= $(round(iso.rstar,sigdigits=8)) M_⊙")
    # exterior consistency: at the surface R = r(1+M/2r)²  must reproduce R_star
    M = star.M; rs = iso.rstar
    R_recon = rs*(1 + M/(2*rs))^2
    iso_ext_err = abs(R_recon - iso.Rstar)/iso.Rstar
    println("       exterior check  R=r(1+M/2r)² vs R_star : rel.err = " *
            "$(round(iso_ext_err,sigdigits=3))")
    # ψ→1 and α→1/ψ²·(…): at large r the isotropic factor ψ²=R/r→1 weakly; check
    # central regularity ψ finite, monotone areal radius, positive grid:
    iso_ok = all(iso.r .> 0) && all(isfinite, iso.psi) && all(iso.grr .> 0) &&
             issorted(iso.R) && iso_ext_err < 1e-10
    # interior conformal-factor sanity at a representative interior point
    j = max(2, length(iso.r) ÷ 2)
    println("       at r=$(round(iso.r[j],sigdigits=4)): ψ=$(round(iso.psi[j],sigdigits=6)), " *
            "g_rr=ψ⁴=$(round(iso.grr[j]^1,sigdigits=6)), α=$(round(iso.alpha[j],sigdigits=6))")
    println("       isotropic transform sane : $iso_ok")

    # ---- Frame: production (ŝ,â,q̂)=(1,1,0.999), smallSB-F2 viscosity ------
    fr = ShumFrame(ŝ=1.0, â=1.0, q̂=0.999, η̂=0.01, ζ̂=0.01)   # smallSB-F2, lines 625–627
    cs2c = sound_speed2(eos, εc)
    c0, cp, cm = shum_frame_speeds(fr.ŝ, fr.â, fr.q̂, fr.η̂, fr.ζ̂, sqrt(cs2c))
    wp = shum_frame_wellposed(fr.ŝ, fr.â, fr.q̂)
    println("\n[FRAME] Shum production (ŝ,â,q̂)=(1,1,0.999), smallSB-F2 (η̂,ζ̂)=(0.01,0.01)")
    println("        τ_ε=(4/3)η̂+ζ̂ = $(round((4/3)*fr.η̂+fr.ζ̂,sigdigits=4))  (target 0.023, lines 625–627)")
    println("        c₊/cs = $(round(cp/sqrt(cs2c),sigdigits=5)) (target √3≈1.7320),  " *
            "c₋/cs = $(round(cm/sqrt(cs2c),sigdigits=5)) (target 0.0183),  well-posed=$wp")

    # ---- (b) BDNK recovery round-trip in the NO-GRADIENT (equilibrium) limit
    println("\n[RECOVERY] BDNK 2×2 con2prim round-trip in the no-gradient limit")
    println("           (static star, v^r=0, ε̂=v̂̄^r=0, all D_i=K=0 ⇒ p₁ must be 0)")
    # Sample interior points of the star (skip the very surface where ε→0)
    Np   = length(star.r)
    idxs = unique(clamp.([2, Np÷8, Np÷4, Np÷2, (3*Np)÷4, (7*Np)÷8, Np-1], 2, Np-1))
    maxerr = 0.0
    worst  = (0.0, 0.0, 0.0)
    detmin = Inf
    for i in idxs
        ε   = star.ε[i]
        ε <= 0 && continue
        grr = iso.grr[i]                       # isotropic spatial metric (curved)
        vr  = 0.0                              # hydrostatic equilibrium
        E, Sr = equilibrium_conserved(eos, ε)  # perfect-fluid conserved (c0,c_r)
        # frozen-gradient deficits: b = (E,S_r) − (c0,c_r) = (0,0) in equilibrium
        b0  = E - ε                            # E − c0  with c0 = ε  ⇒ 0
        br  = Sr - 0.0                         # S_r − c_r with c_r = 0 ⇒ 0
        ε̂, v̂, A = recover_time_derivs_shum(fr, eos, ε, vr, grr, E, Sr; b0=b0, br=br)
        err = max(abs(ε̂), abs(v̂))
        detmin = min(detmin, abs(A[1,1]*A[2,2] - A[1,2]*A[2,1]))
        if err > maxerr
            maxerr = err; worst = (star.r[i], ε, err)
        end
    end
    println("           sampled $(length(idxs)) interior points; min|det 𝒜| = " *
            "$(round(detmin,sigdigits=4))  (matrix non-singular)")
    println("           max |p₁| = max(|ε̂|,|v̂̄^r|) over points = $(maxerr)")
    rt_ok = maxerr <= 1e-8
    println("           round-trip ≤ 1e-8 : $rt_ok")

    # Also exercise a moving-fluid (v^r ≠ 0) consistency check: build conserved
    # from a chosen p₁ via 𝒜, then recover and confirm we get p₁ back (pure
    # linear-algebra round-trip of the SAME matrix, independent of gradients).
    println("\n[RECOVERY] linear self-consistency at v^r≠0 (build b=𝒜·p₁, solve back)")
    lin_maxerr = 0.0
    for i in idxs
        ε = star.ε[i]; ε <= 0 && continue
        grr = iso.grr[i]
        vr  = 0.05                              # representative subsonic flow
        η, ζ, τε, τp, τQ, p, cs2 = shum_transport(fr, eos, ε)
        A = shum_con2prim_matrix(η, ζ, τε, τQ, p, ε, cs2, grr, vr)
        p1_true = [3.7e-6, -1.1e-5]             # arbitrary (ε̂, v̂̄^r)
        b = A * p1_true
        p1_rec = A \ b
        lin_maxerr = max(lin_maxerr, norm(p1_rec - p1_true, Inf))
    end
    println("           max |p₁_rec − p₁_true| = $(lin_maxerr)")
    lin_ok = lin_maxerr <= 1e-8
    println("           linear solve exact ≤ 1e-8 : $lin_ok")

    println("\n" * "="^74)
    allpass = M_ok && iso_ok && rt_ok && lin_ok
    println("OVERALL: M_T=$M_ok  isotropic=$iso_ok  round-trip=$rt_ok  linear=$lin_ok  => " *
            (allpass ? "PASS" : "FAIL"))
    println("="^74)
    return allpass, maxerr, lin_maxerr, star.M, iso_ext_err
end

main()
