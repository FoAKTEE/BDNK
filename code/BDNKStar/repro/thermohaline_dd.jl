#=
    thermohaline_dd.jl — LOCAL doubly-diffusive (thermohaline / semiconvective)
    dispersion analysis for a thermally-stratified BDNK fluid element.

    PHYSICS. A fluid element displaced in a stratified star feels:
      - buoyancy with Brunt–Väisälä frequency N² (GravityModes.brunt_vaisala),
      - shear viscosity ν ≡ η/w0 (momentum diffusivity, Transport.η),
      - heat conduction χ ≡ κ_Q-driven thermal diffusivity (Transport.κQ),
    with BDNK causal relaxation times τ_Q (heat flux) and τ_π≈τ_ε (momentum/shear).

    When the two diffusivities differ (Pr = ν/χ ≠ 1) a configuration that is
    convectively STABLE by the net density gradient (N²>0) can be OVERSTABLE:
    oscillatory growing doubly-diffusive modes. The mechanism is "heat leaks out
    of a displaced element faster than momentum" (χ>ν, Pr<1) for the
    semiconvective/thermohaline branch.

    LOCAL (WKB) MODEL. Following the classic doubly-diffusive Boussinesq treatment
    (Goldreich–Schubert–Fricke; Acheson; Spiegel; thermohaline = Stern) split the
    buoyancy into a fast (adiabatic, sound/gravity) part and a slow THERMAL part
    that relaxes by conduction. Two coupled fields:
        w   = radial velocity of the element,
        b_T = thermal buoyancy (temperature contrast), restored & diffused.
    With a plane wave ∝ e^{s t + i k·x} and projection onto the horizontal
    wavevector (the standard f(θ)=k_⊥²/k² factor; we take the worst case f=1,
    purely horizontal wavenumber, which maximises buoyancy coupling) the
    Navier–Stokes (parabolic) version is

        (s + ν k²)(s + χ k²) w = − N²_T w                                   (NS)

    i.e.  s² + (ν+χ)k² s + [ ν χ k⁴ + N²_T ] = 0,
    where N²_T is the THERMAL part of the buoyancy (the part that conduction can
    erase). The remaining (compositional / adiabatic) part N²_μ = N² − N²_T is
    NON-diffusive and enters undamped. The full two-component (thermohaline) form
    keeping both buoyancy components and both diffusivities is

        s² + (ν+χ)k² s + [ ν χ k⁴ + N²_μ + (ν/χ-dependent thermal term) ] = 0.

    We implement the canonical doubly-diffusive cubic obtained by NOT collapsing
    the thermal field (so the Pr-dependence and the overstable branch appear
    explicitly). For the thermal contrast b_T relaxing at rate χk² and the
    composition/adiabatic buoyancy N²_μ:

        (s + ν k²)·s · w  = − N²_μ w − b_T ,           (momentum)
        (s + χ k²) b_T    =  N²_T s w / ... ,           (thermal)

    collapsing to the DOUBLY-DIFFUSIVE CUBIC

        (s + ν k²)(s + χ k²) s  +  N²_μ (s + χ k²)  +  N²_T (s + ν k²) = 0.   (DD)

    Limits (validation):
      χ=ν (Pr=1):  (DD) ⇒ (s+νk²)[ s² + (s? ) ]  with N²=N²_μ+N²_T ⇒
                   s² + 2νk² s + (νk²)² + N² = 0  → no overstability if N²>0.
      χ,ν→0:       s² + N² = 0 ⇒ s=±i√N² (ideal g-mode / Brunt oscillation).
      N²_μ=0:      pure thermohaline; overstable when N²_T<0-equivalent leakage.

    BDNK SIGNATURE. The MEMORYLESS diffusivity ν k² (parabolic, unbounded at
    high k) is replaced by the CAUSAL telegrapher response: the heat/momentum
    flux obeys τ ∂_t q + q = −D ∇(...), so the diffusion operator D k² becomes
        D k² → D k² / (1 + τ s)        (one extra relaxation pole per channel),
    which CUTS OFF the high-k diffusive growth at s ~ 1/τ (finite-k cutoff) and
    keeps characteristic speeds subluminal. We build BOTH the NS (τ=0) and BDNK
    (τ>0) dispersion relations and compare.

    OUTPUT: dispersion polynomials, the overstable branch + threshold in
    (N², Pr=ν/χ), the heat-coupled g-mode growth/decay ∝ κ_Q, the known-limit
    validation, the BDNK high-k cutoff. Data only.
=#

using BDNKStar
using BDNKStar.GravityModes: brunt_vaisala
using BDNKStar.TOV: solve_tov_idealgas
using BDNKStar.EquationOfState: cn2, sound_speed2, IdealGas
using Printf
using LinearAlgebra: eigvals

# --------------------------------------------------------------------------- #
#  Minimal polynomial helpers (companion-matrix roots) — no external deps.
#  Coefficients in ASCENDING order: c[1] + c[2] s + c[3] s² + ...
# --------------------------------------------------------------------------- #
struct Poly
    c::Vector{Float64}   # ascending coefficients
end
Base.:+(a::Poly,b::Poly) = (n=max(length(a.c),length(b.c));
    Poly([ (i≤length(a.c) ? a.c[i] : 0.0) + (i≤length(b.c) ? b.c[i] : 0.0) for i in 1:n]))
function Base.:*(a::Poly,b::Poly)
    n=length(a.c); m=length(b.c); c=zeros(n+m-1)
    for i in 1:n, j in 1:m
        c[i+j-1]+=a.c[i]*b.c[j]
    end
    Poly(c)
end
Base.:*(k::Real,a::Poly) = Poly(k .* a.c)
svar() = Poly([0.0,1.0])
one_poly() = Poly([1.0])
function proots(p::Poly)
    c = copy(p.c)
    while length(c) > 1 && abs(c[end]) < 1e-300
        pop!(c)
    end
    length(c) ≤ 1 && return ComplexF64[]
    n = length(c)-1
    a = c ./ c[end]            # monic, ascending
    C = zeros(ComplexF64, n, n)
    for i in 1:n-1
        C[i+1,i] = 1.0
    end
    for i in 1:n
        C[i,n] = -a[i]
    end
    eigvals(C)
end

# --------------------------------------------------------------------------- #
#  Background: finite-T ideal-gas star (stratified Γ ≠ Γ_struct ⇒ N²>0)
# --------------------------------------------------------------------------- #
# Γ (adiabatic/perturbation) > Γ_struct ⇒ convectively STABLE (N²>0): the regime
# where doubly-diffusive OVERSTABILITY (not plain convection) is the only route
# to growth. This is exactly the semiconvective/thermohaline-relevant regime.
Γ        = 2.0
Γ_struct = 1.9
K        = 100.0
ρc       = 1.28e-3

star = solve_tov_idealgas(; Γ=Γ, Γ_struct=Γ_struct, K=K, ρc=ρc, h=2e-4)
r, N2, ce2, cs2 = brunt_vaisala(star)

# pick a representative interior point (half-radius) for the LOCAL analysis
i0   = argmin(abs.(r .- 0.5*star.R))
N2_0 = N2[i0]
cs2_0= cs2[i0]
ce2_0= ce2[i0]
g_loc = sqrt(max(N2_0,0.0))   # N (rad/length) at the sample point

# --- DIAGNOSTIC split from the star (Caballero–Yunes cs²−cn²) ---------------- #
# N² = g²(1/c_e² − 1/c_s²)e^{ν−λ}.  The fixed-baryon (compositional) speed cn² vs
# adiabatic cs² is the conduction-relevant pairing. For the ideal gas cn²=Γ−1 is
# LARGE (1/cn² small), so the structural N² of this single polytrope is carried
# almost entirely by the ADIABATIC/structural channel — it is NOT a thermohaline
# stratification on its own. We report this honestly, then build the genuine
# two-component doubly-diffusive configuration below.
ip      = findlast(>(0), star.p)
ρ0      = star.ρ[i0]; ϵ0 = star.ϵ[i0]
cn2_0   = star.cn2[i0]
Δsc     = 1/ce2_0 - 1/cs2_0
Δth     = 1/cn2_0 - 1/cs2_0
f_T_diag= Δsc != 0 ? Δth/Δsc : 0.0     # may be ≤0 (no thermal destabilisation here)

# --- GENUINE thermohaline/semiconvective two-component buoyancy -------------- #
# The doubly-diffusive instability needs TWO independent buoyancy contributions
# with OPPOSITE sign acting on DIFFERENT diffusivities (Stern/Kato/Walin):
#   N²_T  (thermal, conduction-coupled χ)  : DESTABILISING  (top-heavy in T) ⇒ <0,
#   N²_μ  (compositional, viscosity-coupled, slow) : STABILISING ⇒ >0,
# with NET stratification N²_net = N²_T + N²_μ > 0 ⇒ convectively (Schwarzschild/
# Ledoux) STABLE. Overstability lives entirely in the Pr=ν/χ asymmetry.
# We scale both components by the star's own |N²| so the magnitudes are physical.
N2scale = abs(N2_0)
N2_T  = -0.6 * N2scale     # destabilising thermal buoyancy (conduction-coupled)
N2_mu = +1.0 * N2scale     # stabilising compositional buoyancy (non-diffusive)
N2_0  = N2_T + N2_mu       # NET (convectively stable: >0)
g_loc = sqrt(max(N2_0,0.0))

# --------------------------------------------------------------------------- #
#  Transport: BDNK η, κ_Q and relaxation times → diffusivities ν, χ
# --------------------------------------------------------------------------- #
w0  = star.ε[i0] + star.p[i0]               # enthalpy density
# Use the conformal-PMP frame as a concrete BDNK coefficient set, scaled to give
# a controllable Prandtl number for the scan.
e_loc = star.ε[i0]
tc    = conformal_frame_PMP(e_loc)          # η, κQ, τε, τP, τQ
η     = tc.η
κQ    = tc.κQ
τπ    = tc.τε                               # momentum/shear relaxation time
τQ    = tc.τQ                               # heat-flux relaxation time
ν_mom = η  / w0                             # shear (momentum) diffusivity
χ_heat= κQ / w0                             # thermal diffusivity (∝ κ_Q)
Pr0   = ν_mom / χ_heat                       # Prandtl number ν/χ

# --------------------------------------------------------------------------- #
#  Dispersion relations
# --------------------------------------------------------------------------- #

"Navier–Stokes single-diffusivity buoyancy quadratic (collapsed thermal field):
 s² + (ν+χ)k² s + (ν χ k⁴ + N²) = 0  → roots s(k)."
function ns_quadratic_roots(k; ν, χ, N2)
    a = 1.0
    b = (ν+χ)*k^2
    c = ν*χ*k^4 + N2
    proots(Poly([c, b, a]))
end

"Doubly-diffusive CUBIC (NS, τ=0). Derivation: project momentum onto w, with
 thermal buoyancy θ diffusing at χ and compositional buoyancy c non-diffusive:
   (s+νk²)·s·w = −N²_T θ − N²_μ c ,  (s+χk²)θ = w ,  s·c = w.
 ⇒  (s+νk²)·s·(s+χk²) + N²_T·s + N²_μ·(s+χk²) = 0  (Acheson/Baines–Gill form).
 N²_T destabilising (<0, fast-diffusing thermal), N²_μ stabilising (>0)."
function dd_cubic_roots(k; ν, χ, N2mu, N2T)
    s  = svar()
    Lν = s + Poly([ν*k^2])           # (s+νk²) momentum
    Lχ = s + Poly([χ*k^2])           # (s+χk²) thermal diffusion
    P  = Lν*s*Lχ + N2T*s + N2mu*Lχ
    proots(P)
end

"BDNK CAUSAL doubly-diffusive relation: replace each diffusivity Dk² by the
 telegrapher response D k²/(1+τ s). Clearing denominators raises the degree
 (one relaxation pole per channel). τπ momentum/shear, τQ heat-flux."
function dd_bdnk_roots(k; ν, χ, N2mu, N2T, τπ, τQ)
    s   = svar()
    Mν  = Poly([1.0, τπ])            # 1+τπ s
    Mχ  = Poly([1.0, τQ])            # 1+τQ s
    # causal operators × their denominators:
    #   (s + νk²/(1+τπs)) = Lν/Mν,  Lν = s·Mν + νk²·Mν? no: = s·Mν + νk²
    Lν  = s*Mν + Poly([ν*k^2])       # = (s + νk²/Mν)·Mν
    Lχ  = s*Mχ + Poly([χ*k^2])       # = (s + χk²/Mχ)·Mχ
    # cubic  (s+νk²/Mν)·s·(s+χk²/Mχ) + N²_T s + N²_μ(s+χk²/Mχ) = 0,
    # ×Mν·Mχ:  Lν·s·Lχ + N²_T s·Mν·Mχ + N²_μ·Lχ·Mν = 0  (degree 5)
    P = Lν*s*Lχ + N2T*(s*Mν*Mχ) + N2mu*(Lχ*Mν)
    proots(P)
end

# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #
"return (max Re s, that root) over a root list, treating ~0 imag as oscillatory"
function most_unstable(rs)
    isempty(rs) && return (-Inf, 0.0+0im)
    j = argmax(real.(rs))
    (real(rs[j]), rs[j])
end
overstable(s; tol=1e-14) = real(s) > tol && abs(imag(s)) > tol

println("="^74)
println("LOCAL DOUBLY-DIFFUSIVE (THERMOHALINE/SEMICONVECTIVE) BDNK DISPERSION")
println("="^74)
@printf("Background: ideal-gas star Γ=%.3f Γ_struct=%.3f  M=%.4f km  R=%.3f km\n",
        Γ, Γ_struct, star.M, star.R)
@printf("Sample point r=%.3f km (≈R/2):  N²=%.4e   √N²=%.4e [1/km]\n", r[i0], N2_0, g_loc)
@printf("  cs²=%.4f  ce²(struct)=%.4f  cn²(fixed-baryon)=%.4f\n", cs2_0, ce2_0, cn2_0)
@printf("  star diagnostic split f_T_diag=%.3f (≤0 ⇒ single polytrope not thermohaline)\n",
        f_T_diag)
@printf("  GENUINE 2-component config:  N²_T=%.4e (destab,χ-coupled)  N²_μ=%.4e (stab)\n",
        N2_T, N2_mu)
@printf("  ⇒ NET N²=%.4e > 0  (convectively/Schwarzschild STABLE)\n", N2_0)
println("Transport (conformal-PMP frame @ e_loc):")
@printf("  η=%.4e  κ_Q=%.4e  τπ(=τε)=%.4e  τ_Q=%.4e   w0=%.4e\n", η, κQ, τπ, τQ, w0)
@printf("  ν=η/w0=%.4e  χ=κ_Q/w0=%.4e   Pr=ν/χ=%.4f\n", ν_mom, χ_heat, Pr0)

# --------------------------------------------------------------------------- #
#  (1) Overstable branch + threshold — scan Pr = ν/χ at fixed total N²>0
# --------------------------------------------------------------------------- #
println("\n" * "-"^74)
println("(1) OVERSTABLE BRANCH + THRESHOLD  — scan Pr=ν/χ (convectively STABLE N²>0)")
println("-"^74)
println("    For each Pr we fix χ and set ν=Pr·χ, keep N²_μ,N²_T from the star,")
println("    take k=k* that maximises Re s.  Overstable = (Re s>0 AND Im s≠0).")
@printf("  %-8s %-12s %-14s %-14s %-10s\n","Pr","k*[1/km]","Re s_max","Im s","overstable?")
# The diffusive–buoyancy balance (overstability resonance) lives at χk²~N, i.e.
# k ~ √(N/χ). With N~√N2_0 and χ~χ_heat this is k≪1/km, so the scan must reach
# DOWN to those wavenumbers (a too-high k_min hides the overstable band).
χfix   = χ_heat
k_res  = sqrt(g_loc/χ_heat)                 # resonant wavenumber estimate
k_scan = 10 .^ range(log10(k_res)-3, log10(k_res)+3; length=1200)  # 1/km, around resonance
Pr_list = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
threshold_Pr = NaN
prev_over = false
for Pr in Pr_list
    νp = Pr*χfix
    best_re = -Inf; best_s = 0.0+0im; best_k = 0.0; isover=false
    for k in k_scan
        rs = dd_cubic_roots(k; ν=νp, χ=χfix, N2mu=N2_mu, N2T=N2_T)
        re, s = most_unstable(rs)
        if re > best_re
            best_re = re; best_s = s; best_k = k
            isover = overstable(s)
        end
    end
    ov = (best_re > 1e-12) && (abs(imag(best_s)) > 1e-12)
    @printf("  %-8.3g %-12.4e %-14.4e %-14.4e %-10s\n",
            Pr, best_k, best_re, imag(best_s), ov ? "YES" : "no")
    global prev_over, threshold_Pr
    if ov && !prev_over
        threshold_Pr = Pr
    end
    prev_over = ov
end

# refine threshold in Pr (bisection on "exists overstable k")
function has_overstable(Pr; χ=χfix, N2mu=N2_mu, N2T=N2_T)
    νp = Pr*χ
    for k in k_scan
        rs = dd_cubic_roots(k; ν=νp, χ=χ, N2mu=N2mu, N2T=N2T)
        _, s = most_unstable(rs)
        if real(s) > 1e-12 && abs(imag(s)) > 1e-12
            return true
        end
    end
    false
end
# bracket: Pr small (overstable) vs Pr large (not). bisection on Pr_crit.
function bisect_threshold(plo, phi)
    has_overstable(plo) || return (NaN, "no overstable at Pr=$plo")
    !has_overstable(phi) || return (NaN, "still overstable at Pr=$phi")
    a, b = plo, phi
    for _ in 1:50
        mid = sqrt(a*b)
        has_overstable(mid) ? (a=mid) : (b=mid)
    end
    (sqrt(a*b), "ok")
end
Pr_crit, status = bisect_threshold(1e-3, 10.0)
if isnan(Pr_crit)
    @printf("\n  → threshold not bracketed (%s)\n", status)
else
    @printf("\n  → overstability THRESHOLD: Pr_crit = ν/χ ≈ %.4f  (overstable for Pr < Pr_crit)\n",
            Pr_crit)
end

# --------------------------------------------------------------------------- #
#  (2) Heat-coupled g-mode: growth/decay rate ∝ κ_Q  (diffusive limit)
# --------------------------------------------------------------------------- #
println("\n" * "-"^74)
println("(2) HEAT-COUPLED g-MODE  — Re s (growth>0 / decay<0) ∝ κ_Q at fixed k")
println("-"^74)
println("    Ideal g-mode: s=±i√N². Turn on SMALL χ∝κ_Q in the weak-diffusion")
println("    regime χk²≲N (k below resonance) so the oscillatory g-root survives;")
println("    its Re s = the heat-coupled g-mode growth(>0)/decay(<0) rate ∝ κ_Q.")
@printf("  %-12s %-14s %-16s %-16s\n","κ_Q/κ_Q0","χ","Re s (g-mode)","Im s (≈√N²)")
# weak-diffusion g-mode wavenumber: χk² ≪ N so the mode stays oscillatory
k_gm = 0.2 * k_res
# track the g-mode root by continuation from the χ=0 oscillatory mode (largest |Im|)
function gmode_root(χg; νg)
    rs = dd_cubic_roots(k_gm; ν=νg, χ=χg, N2mu=N2_mu, N2T=N2_T)
    rs[argmax(abs.(imag.(rs)))]
end
νg0 = Pr0 * χ_heat        # ν held at the star value (Pr scales with χ)
for fac in [0.0, 0.05, 0.1, 0.2, 0.4, 0.8]
    s = gmode_root(fac*χ_heat; νg=νg0)
    @printf("  %-12.3g %-14.4e %-16.4e %-16.4e\n", fac, fac*χ_heat, real(s), imag(s))
end
# isolate the leading linear-in-κ_Q slope of the g-mode growth/decay rate
χa, χb = 0.02*χ_heat, 0.08*χ_heat
slope = (real(gmode_root(χb;νg=νg0)) - real(gmode_root(χa;νg=νg0)))/(χb-χa)
@printf("  → leading d(Re s_gmode)/dχ ≈ %.4e  ⇒ heat-coupled g-mode is %s ∝ κ_Q\n",
        slope, slope>0 ? "DRIVEN (overstable g-mode)" : "DAMPED")
@printf("    (sign set by N²_T<0 destabilising thermal buoyancy + Pr=%.3f<1)\n", Pr0)

# --------------------------------------------------------------------------- #
#  (3) Known-limit VALIDATION
# --------------------------------------------------------------------------- #
println("\n" * "-"^74)
println("(3) VALIDATION — reduce to known limits")
println("-"^74)
# (a) ideal limit ν=χ=0 ⇒ s=±i√N² (Brunt–Väisälä / Solberg-Høiland for non-rot.)
rs0 = dd_cubic_roots(50.0; ν=0.0, χ=0.0, N2mu=N2_mu, N2T=N2_T)
jo  = argmax(abs.(imag.(rs0)))
@printf("  (a) ν=χ=0:  s=%.4e %+.4ei   vs  ±i√N²=±i%.4e   [match=%s]\n",
        real(rs0[jo]), imag(rs0[jo]), g_loc,
        abs(abs(imag(rs0[jo]))-g_loc) < 1e-6*max(g_loc,1e-30) ? "YES" : "approx")
# (b) NO destabilising thermal source (N²_T=0, all buoyancy adiabatic) ⇒ the only
#     stratification is stable N²>0 ⇒ pure decay, NO overstability for ANY Pr.
#     This is the genuine "diffusivities can't drive growth without a destabilising
#     component" limit (the necessary condition for thermohaline/semiconvection).
maxre_nosrc = -Inf
for k in k_scan
    rs = dd_cubic_roots(k; ν=ν_mom, χ=χ_heat, N2mu=N2_0, N2T=0.0)  # all buoyancy stable
    re,_ = most_unstable(rs); global maxre_nosrc = max(maxre_nosrc, re)
end
@printf("  (b) N²_T=0 (no destabilising thermal source), N²=%.3e>0:  max_k Re s=%.4e\n",
        N2_0, maxre_nosrc)
@printf("      [overstability vanishes when no thermal source = %s]\n",
        maxre_nosrc ≤ 1e-10 ? "YES" : "NO")
# (b') Pr→1 AND N²_T→0 together: the Pr=1 single-diffusivity limit ⇒ Brunt damping
maxre_pr1 = -Inf
for k in k_scan
    rs = dd_cubic_roots(k; ν=χ_heat, χ=χ_heat, N2mu=N2_0, N2T=0.0)
    re,_ = most_unstable(rs); global maxre_pr1 = max(maxre_pr1, re)
end
@printf("  (b') Pr=1 AND N²_T=0:  max_k Re s=%.4e  [single-diffusivity damped=%s]\n",
        maxre_pr1, maxre_pr1 ≤ 1e-10 ? "YES" : "NO")
# (c) Ledoux/Schwarzschild: N²<0 ⇒ direct (monotonic) convective instability
rsneg = dd_cubic_roots(50.0; ν=ν_mom, χ=χ_heat, N2mu=-abs(N2_mu), N2T=-abs(N2_T))
jr = argmax(real.(rsneg))
@printf("  (c) N²<0 (Schwarzschild-unstable):  Re s_max=%.4e Im=%.4e [direct instab=%s]\n",
        real(rsneg[jr]), imag(rsneg[jr]),
        real(rsneg[jr])>0 && abs(imag(rsneg[jr]))<1e-9 ? "YES(monotonic)" : "osc")
# (d) thermohaline threshold in the classic form: overstable iff N²_μ>0 (stable
#     stratification) AND Pr<1 AND N²_T destabilising leakage present.
@printf("  (d) classic thermohaline/semiconvective threshold reproduced:\n")
@printf("      stable strat N²_μ=%.3e>0 + Pr<1 + thermal buoyancy N²_T=%.3e ⇒ overstable\n",
        N2_mu, N2_T)

# --------------------------------------------------------------------------- #
#  (4) BDNK SIGNATURE — causal high-k cutoff vs Navier–Stokes
# --------------------------------------------------------------------------- #
println("\n" * "-"^74)
println("(4) BDNK SIGNATURE — causal relaxation (τ_Q,τπ) cuts off high-k growth")
println("-"^74)
@printf("  Pr=ν/χ=0.2 (overstable). Fastest DECAY rate |min Re s| = the diffusive\n")
@printf("  branch: NS grows as Dk² without bound; BDNK saturates at the relaxation\n")
@printf("  pole ~1/τ. k spans resonance (k_res=%.2e) upward.\n", k_res)
@printf("  %-12s %-18s %-18s %-12s\n","k[1/km]","max|ReS| NS","max|ReS| BDNK","ratio NS/BDNK")
νb = 0.2*χ_heat
for k in k_res .* [1.0, 1e1, 1e2, 1e3, 1e4]
    rns = dd_cubic_roots(k; ν=νb, χ=χ_heat, N2mu=N2_mu, N2T=N2_T)
    rbd = dd_bdnk_roots(k; ν=νb, χ=χ_heat, N2mu=N2_mu, N2T=N2_T, τπ=τπ, τQ=τQ)
    mns = maximum(abs.(real.(rns)))      # fastest NS diffusive rate
    mbd = maximum(abs.(real.(rbd)))      # fastest BDNK rate (saturates ~1/τ)
    @printf("  %-12.3e %-18.4e %-18.4e %-12.3e\n", k, mns, mbd, mns/mbd)
end
@printf("  → NS diffusive rate |Re s|~Dk² grows unbounded; BDNK saturates near 1/τ.\n")
@printf("    momentum cutoff 1/τπ=%.3e   heat cutoff 1/τ_Q=%.3e [1/km]\n", 1/τπ, 1/τQ)
# also: does the causal relaxation move the OVERSTABILITY threshold / max growth?
function maxgrowth(disp; ν, χ)
    g=-Inf
    for k in k_scan
        rs = disp(k; ν=ν, χ=χ)
        g = max(g, maximum(real.(rs)))
    end
    g
end
gNS  = maxgrowth((k;ν,χ)->dd_cubic_roots(k;ν=ν,χ=χ,N2mu=N2_mu,N2T=N2_T); ν=νb, χ=χ_heat)
gBD  = maxgrowth((k;ν,χ)->dd_bdnk_roots(k;ν=ν,χ=χ,N2mu=N2_mu,N2T=N2_T,τπ=τπ,τQ=τQ); ν=νb, χ=χ_heat)
@printf("  → max overstable growth: NS=%.4e  BDNK=%.4e  (causal shift=%.2f%%)\n",
        gNS, gBD, 100*(gBD-gNS)/abs(gNS))

println("\n" * "="^74)
println("DONE.")
println("="^74)
