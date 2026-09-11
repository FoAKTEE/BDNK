#=
    DoublyDiffusive.jl — LOCAL (WKB) doubly-diffusive dispersion analysis of a
    thermally-stratified, viscous + heat-conducting BDNK (causal first-order)
    fluid element. Two known instability families, re-derived in the BDNK
    framework with the toolkit's own buoyancy (GravityModes.brunt_vaisala) and
    transport sector (Transport.conformal_frame_PMP → η, κ_Q, τ):

      (1) GSF / rotational doubly-diffusive (Goldreich–Schubert 1967; Fricke
          1968; local form Acheson 1978; Menou, Balbus & Spruit 2004). A
          differentially-rotating, thermally-stratified element. Restoring
          agents: epicyclic κ_epi² (rotation), buoyancy N² (diffuses with the
          THERMAL diffusivity κ_th), momentum (diffuses with ν). Meridional
          plane-wave ξ∝exp(i k·r + s t), k=(k_R,k_z)=k(sinψ,cosψ).

            (s+d_ν)²(s+d_κ) + (s+d_κ)κ_epi²cos²ψ + (s+d_ν)N²sin²ψ = 0   (NS)

          d_ν=νk², d_κ=κ_th k². IDEAL limit (d_ν,d_κ→0): s²=−(κ_epi²cos²ψ
          +N²sin²ψ) ⇒ stable ∀ψ iff κ_epi²≥0 ∧ N²≥0 — the Solberg–Høiland /
          Rayleigh criterion.

      (2) Thermohaline / semiconvective doubly-diffusive overstability (Stern
          1960; Kato; Baines–Gill; Acheson). A non-rotating element with TWO
          opposed buoyancy components — a destabilising THERMAL part N²_T<0
          (conduction-coupled, diffuses with χ=κ_th) and a stabilising
          compositional/adiabatic part N²_μ>0 (non-diffusive) — net N²>0
          (convectively stable). Doubly-diffusive cubic (Acheson/Baines–Gill):

            (s+νk²)·s·(s+χk²) + N²_T·s + N²_μ·(s+χk²) = 0                (NS)

          IDEAL limit (ν=χ=0): s=±i√N² (Brunt–Väisälä g-mode). N²<0 ⇒ direct
          monotonic Schwarzschild/Ledoux convection. Overstable (Re s>0,
          Im s≠0) when N²_μ>0, N²_T<0 and Pr=ν/χ<1 (heat leaks faster than
          momentum). The heat-coupled g-mode is the weak-diffusion (χk²≲N)
          continuation of this branch: the oscillatory g-root acquires a
          Re s ∝ κ_Q whose sign is set by N²_T.

    BDNK CAUSAL MODIFICATION (the headline, common to both families). Navier–
    Stokes uses INSTANTANEOUS diffusion D k². BDNK / IS-class causal transport
    gives the flux a finite relaxation time τ: for a mode ∝ e^{st} the effective
    diffusivity is the telegrapher response

        D k²  →  D k² / (1 + s τ)      (τ_η for momentum, τ_Q for heat),

    one extra relaxation pole per channel (cubic → quintic). Clearing the (1+sτ)
    denominators and solving the resulting polynomial by companion-matrix
    eigenvalues exposes the BDNK signature:
      • thermohaline branch: the NS diffusive rate |Re s|~Dk² grows WITHOUT
        bound at high k (acausal); BDNK SATURATES it near the relaxation pole
        ~1/τ — a finite-k cutoff — while leaving the low-k overstable band intact.
      • GSF branch: at the MARGINAL boundary s→0, D/(1+sτ)→Dk², so BDNK≡NS there;
        finite τ does NOT move the onset but SUPPRESSES the secular growth rate
        of already-unstable modes by O(sτ).

    HONEST PROVENANCE. The instabilities themselves are TEXTBOOK (GSF 1967;
    thermohaline / Stern 1960) — NOT new families. What is new is their first
    realization on a BDNK causal-hydro neutron-star background wired to this
    toolkit's N² and transport frame, plus the explicit BDNK causal modification
    (telegrapher D/(1+sτ)) and its honest consequences (UV cutoff for
    thermohaline; growth-rate suppression with τ-invariant onset for GSF).

    REFERENCES: Goldreich & Schubert 1967 ApJ 150 571; Fricke 1968 ZAp 68 317;
    Acheson 1978 Phil.Trans.R.Soc.A 289 459; Knobloch & Spruit 1982 A&A 113 261;
    Menou, Balbus & Spruit 2004 ApJ 607 564; Stern 1960 Tellus 12 172;
    Baines & Gill 1969 JFM 37 289. Causal regulator: telegrapher/Cattaneo (BDNK).
=#
module DoublyDiffusive

using LinearAlgebra: eigvals

export gsf_roots_NS, gsf_roots_BDNK, gsf_max_growth_NS, gsf_max_growth_BDNK,
       gsf_epicyclic, gsf_unstable_cone,
       dd_cubic_roots, dd_bdnk_roots, dd_max_growth, dd_is_overstable,
       gmode_root, gmode_growth_slope,
       most_unstable, fastest_rate, real_roots_poly

# --------------------------------------------------------------------------- #
#  Polynomial root helpers — companion-matrix eigenvalues, no external deps.
#  Coefficients are ASCENDING: c[1] + c[2] s + c[3] s² + …
# --------------------------------------------------------------------------- #

"All (complex) roots of a polynomial with ASCENDING coefficients, via the
 companion matrix. Trailing ~0 leading coefficients are trimmed."
function roots_companion(coeffs::AbstractVector{<:Real})
    c = collect(float.(coeffs))
    sc = maximum(abs.(c))
    while length(c) > 1 && abs(c[end]) < 1e-30 * max(sc, 1.0)
        pop!(c)
    end
    n = length(c) - 1
    n ≤ 0 && return ComplexF64[]
    a = c ./ c[end]                  # monic, ascending
    C = zeros(ComplexF64, n, n)
    for i in 1:n-1
        C[i+1, i] = 1.0
    end
    for i in 1:n
        C[i, n] = -a[i]
    end
    return eigvals(C)
end

"Real roots of a polynomial given ASCENDING coefficients (imag part ~0)."
function real_roots_poly(coeffs::AbstractVector{<:Real}; tol::Float64=1e-7)
    rts = roots_companion(coeffs)
    return Float64[real(z) for z in rts if abs(imag(z)) < tol * max(1.0, abs(real(z)))]
end

"Multiply two ASCENDING-coefficient polynomials."
function poly_mul(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    c = zeros(length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i+j-1] += a[i] * b[j]
    end
    return c
end

"Sum a list of ASCENDING-coefficient polynomials (zero-padded to common length)."
function poly_add(ps...)
    L = maximum(length, ps)
    c = zeros(L)
    for p in ps, i in eachindex(p)
        c[i] += p[i]
    end
    return c
end

# --------------------------------------------------------------------------- #
#  (1) GSF / rotational doubly-diffusive dispersion
# --------------------------------------------------------------------------- #

"""
    gsf_epicyclic(Ω, q_shear)

Radial epicyclic frequency squared for a power-law differential rotation
Ω(ϖ)=Ω0(ϖ/R)^q,  κ_epi² = 4Ω²(1 + q/2) = (1/ϖ³)d(ϖ²Ω)²/dϖ. κ_epi²<0 ⇒ the
angular-momentum gradient is Rayleigh-unstable on cylinders.
"""
gsf_epicyclic(Ω::Real, q_shear::Real) = 4Ω^2 * (1 + q_shear/2)

"""
    gsf_roots_NS(κ2, N2, ν, κth, k, ψ) -> Vector{Float64}

Real roots s of the Navier–Stokes (instantaneous-diffusion) GSF dispersion
cubic  (s+d_ν)²(s+d_κ) + (s+d_κ)κ²cos²ψ + (s+d_ν)N²sin²ψ = 0,
with d_ν=νk², d_κ=κth·k². ψ = angle of k to the rotation axis (cosψ=k_z/k).
Setting ν=κth=0 gives the IDEAL Solberg–Høiland limit s²=−(κ²cos²ψ+N²sin²ψ).
"""
function gsf_roots_NS(κ2::Real, N2::Real, ν::Real, κth::Real, k::Real, ψ::Real)
    dν = ν   * k^2
    dκ = κth * k^2
    c2 = cos(ψ)^2; s2 = sin(ψ)^2
    a3 = 1.0
    a2 = 2dν + dκ
    a1 = dν^2 + 2dν*dκ + κ2*c2 + N2*s2
    a0 = dν^2*dκ + dκ*κ2*c2 + dν*N2*s2
    return real_roots_poly([a0, a1, a2, a3])
end

"""
    gsf_roots_BDNK(κ2, N2, ν, κth, k, ψ, τη, τQ) -> Vector{Float64}

Real roots of the BDNK CAUSAL GSF dispersion: the diffusive decay rates carry
finite relaxation times d_ν→νk²/(1+sτη), d_κ→κth·k²/(1+sτQ). Clearing the
(1+sτ) denominators raises the cubic to degree ≤5; solved by companion matrix.
As τη,τQ→0 this reduces exactly to `gsf_roots_NS`.
"""
function gsf_roots_BDNK(κ2::Real, N2::Real, ν::Real, κth::Real, k::Real, ψ::Real,
                        τη::Real, τQ::Real)
    c2 = cos(ψ)^2; s2 = sin(ψ)^2
    Dν = ν*k^2; Dκ = κth*k^2
    # (s + d_ν(s))·(1+sτη) = Dν + s + τη s²  ≡ Pν   (numerator of s+d_ν)
    Pν  = [Dν, 1.0, τη]
    Pκ  = [Dκ, 1.0, τQ]
    Dnν = [1.0, τη]                  # 1 + s τη
    DnQ = [1.0, τQ]                  # 1 + s τQ
    # common denominator Dnν² DnQ; sum of numerators = 0:
    Pν2  = poly_mul(Pν, Pν)
    Dnν2 = poly_mul(Dnν, Dnν)
    num1 = poly_mul(Pν2, Pκ)                              # (s+dν)²(s+dκ)
    num2 = poly_mul(poly_mul(Pκ, Dnν2), [κ2*c2])         # (s+dκ)κ²c2
    num3 = poly_mul(poly_mul(Pν, poly_mul(Dnν, DnQ)), [N2*s2])  # (s+dν)N²s2
    return real_roots_poly(poly_add(num1, num2, num3))
end

"Max real growth rate over a root list (−Inf if empty)."
fastest_rate(roots) = isempty(roots) ? -Inf : maximum(roots)

gsf_max_growth_NS(κ2, N2, ν, κth, k, ψ) =
    fastest_rate(gsf_roots_NS(κ2, N2, ν, κth, k, ψ))
gsf_max_growth_BDNK(κ2, N2, ν, κth, k, ψ, τη, τQ) =
    fastest_rate(gsf_roots_BDNK(κ2, N2, ν, κth, k, ψ, τη, τQ))

"""
    gsf_unstable_cone(κ2, N2, ν, κth, k; τη=0, τQ=0, nψ=181) -> (ψmax_deg, peak_s, ψpeak_deg)

Scan ψ∈[0,π/2] and return the half-opening angle (deg) of the unstable
wavevector cone, the peak growth rate, and the ψ where it peaks. τ>0 selects
the BDNK causal dispersion; τ=0 the NS one.
"""
function gsf_unstable_cone(κ2::Real, N2::Real, ν::Real, κth::Real, k::Real;
                           τη::Real=0.0, τQ::Real=0.0, nψ::Int=181, tol::Float64=1e-12)
    cone = 0.0; peak = -Inf; ψpk = 0.0
    for ψ in range(0, π/2; length=nψ)
        s = (τη == 0 && τQ == 0) ? gsf_max_growth_NS(κ2, N2, ν, κth, k, ψ) :
                                   gsf_max_growth_BDNK(κ2, N2, ν, κth, k, ψ, τη, τQ)
        if s > tol
            cone = max(cone, rad2deg(ψ))
        end
        if s > peak
            peak = s; ψpk = ψ
        end
    end
    return (cone, peak, rad2deg(ψpk))
end

# --------------------------------------------------------------------------- #
#  (2) Thermohaline / semiconvective doubly-diffusive dispersion
# --------------------------------------------------------------------------- #

"""
    dd_cubic_roots(k; ν, χ, N2mu, N2T) -> Vector{ComplexF64}

Roots s of the Navier–Stokes thermohaline doubly-diffusive cubic
  (s+νk²)·s·(s+χk²) + N²_T·s + N²_μ·(s+χk²) = 0,
ν momentum diffusivity, χ thermal diffusivity, N²_T destabilising thermal
buoyancy (conduction-coupled, <0), N²_μ stabilising compositional buoyancy
(non-diffusive, >0). ν=χ=0 ⇒ s=±i√(N²_T+N²_μ) (ideal Brunt–Väisälä g-mode).
"""
function dd_cubic_roots(k::Real; ν::Real, χ::Real, N2mu::Real, N2T::Real)
    Dν = ν*k^2; Dχ = χ*k^2
    # (s+Dν)·s·(s+Dχ) = s³ + (Dν+Dχ)s² + Dν Dχ s
    # + N²_T s
    # + N²_μ(s+Dχ) = N²_μ Dχ + N²_μ s
    a0 = N2mu*Dχ
    a1 = Dν*Dχ + N2T + N2mu
    a2 = Dν + Dχ
    a3 = 1.0
    return roots_companion([a0, a1, a2, a3])
end

"""
    dd_bdnk_roots(k; ν, χ, N2mu, N2T, τπ, τQ) -> Vector{ComplexF64}

BDNK CAUSAL thermohaline dispersion: each diffusivity Dk² → Dk²/(1+τs)
(τπ momentum, τQ heat). Clearing the two relaxation denominators raises the
cubic to degree 5. As τπ,τQ→0 reduces to `dd_cubic_roots`. The high-k
diffusive rate saturates near 1/τ (the causal UV cutoff).
"""
function dd_bdnk_roots(k::Real; ν::Real, χ::Real, N2mu::Real, N2T::Real,
                       τπ::Real, τQ::Real)
    Dν = ν*k^2; Dχ = χ*k^2
    Mν = [1.0, τπ]                   # 1 + τπ s
    Mχ = [1.0, τQ]                   # 1 + τQ s
    sP = [0.0, 1.0]                  # s
    Lν = poly_add(poly_mul(sP, Mν), [Dν])   # s·Mν + Dν  = (s+Dν/Mν)·Mν
    Lχ = poly_add(poly_mul(sP, Mχ), [Dχ])   # s·Mχ + Dχ
    # Lν·s·Lχ + N²_T·s·Mν·Mχ + N²_μ·Lχ·Mν = 0  (degree 5)
    t1 = poly_mul(poly_mul(Lν, sP), Lχ)
    t2 = poly_mul([N2T], poly_mul(sP, poly_mul(Mν, Mχ)))
    t3 = poly_mul([N2mu], poly_mul(Lχ, Mν))
    return roots_companion(poly_add(t1, t2, t3))
end

"(max Re s, the root achieving it) over a root list."
function most_unstable(rs)
    isempty(rs) && return (-Inf, 0.0 + 0.0im)
    j = argmax(real.(rs))
    return (real(rs[j]), rs[j])
end

"Is the most-unstable root genuinely OVERSTABLE (Re s>0 AND Im s≠0)?"
function dd_is_overstable(rs; tol::Float64=1e-12)
    re, s = most_unstable(rs)
    return re > tol && abs(imag(s)) > tol
end

"""
    dd_max_growth(k_scan; ν, χ, N2mu, N2T, bdnk=false, τπ=0, τQ=0) -> (max Re s, k*)

Max growth rate of the thermohaline dispersion over a wavenumber list.
`bdnk=true` uses the causal dispersion with (τπ,τQ).
"""
function dd_max_growth(k_scan; ν, χ, N2mu, N2T, bdnk::Bool=false, τπ=0.0, τQ=0.0)
    best = -Inf; kbest = first(k_scan)
    for k in k_scan
        rs = bdnk ? dd_bdnk_roots(k; ν=ν, χ=χ, N2mu=N2mu, N2T=N2T, τπ=τπ, τQ=τQ) :
                    dd_cubic_roots(k; ν=ν, χ=χ, N2mu=N2mu, N2T=N2T)
        re, _ = most_unstable(rs)
        if re > best
            best = re; kbest = k
        end
    end
    return (best, kbest)
end

# --------------------------------------------------------------------------- #
#  Heat-coupled g-mode (weak-diffusion continuation of the thermohaline branch)
# --------------------------------------------------------------------------- #

"""
    gmode_root(k, χ; ν, N2mu, N2T) -> ComplexF64

The oscillatory g-mode root (largest |Im s|) of the thermohaline cubic, tracked
by continuation from the χ=0 Brunt oscillation. In the weak-diffusion regime
χk²≲N its Re s is the heat-coupled g-mode growth(>0)/decay(<0) rate ∝ κ_Q.
"""
function gmode_root(k::Real, χ::Real; ν::Real, N2mu::Real, N2T::Real)
    rs = dd_cubic_roots(k; ν=ν, χ=χ, N2mu=N2mu, N2T=N2T)
    return rs[argmax(abs.(imag.(rs)))]
end

"""
    gmode_growth_slope(k; ν, N2mu, N2T, χa, χb) -> Float64

Leading slope d(Re s_gmode)/dχ of the heat-coupled g-mode between two small
thermal diffusivities χa<χb (weak-diffusion). Sign>0 ⇒ DRIVEN/overstable
g-mode (set by N²_T<0 and Pr<1); <0 ⇒ damped.
"""
function gmode_growth_slope(k::Real; ν::Real, N2mu::Real, N2T::Real,
                            χa::Real, χb::Real)
    sa = real(gmode_root(k, χa; ν=ν, N2mu=N2mu, N2T=N2T))
    sb = real(gmode_root(k, χb; ν=ν, N2mu=N2mu, N2T=N2T))
    return (sb - sa) / (χb - χa)
end

end # module DoublyDiffusive
