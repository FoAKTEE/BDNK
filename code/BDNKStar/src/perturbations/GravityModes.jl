#=
    GravityModes — relativistic Cowling POLAR (ℓ≥2) pulsation eigensolver that
    INCLUDES the buoyancy / Schwarzschild-discriminant term, so it yields the
    g (gravity / buoyancy) modes in ADDITION to the f and p modes.

    This is the STRATIFIED counterpart of NonRadialModes (the validated
    barotropic f/p (W,V) eigensolver). The only physical addition is buoyancy:
    a fluid element perturbed adiabatically (index Γ, "frozen composition")
    against a background whose hydrostatic structure follows a DIFFERENT index
    Γ_struct ("equilibrium") feels a restoring buoyant force set by the
    Brunt–Väisälä frequency N². When Γ = Γ_struct (isentropic/barotropic) N²≡0
    and the spectrum collapses back to f/p only.

    EQUATIONS (relativistic Cowling, NON-rotating). The (U,V) form used here is
    the McDermott–Van Horn–Scholl / Reisenegger–Goldreich / Kantor–Gusakov
    system as written by Jaikumar et al. 2021 and Shirke et al. (DM g-modes,
    arXiv:2506.18892, their Eq. 12-13), cross-checked against the Schwarzschild
    discriminant of Gaertig & Kokkotas 2009 (arXiv:0905.0821, their Eq. 12):

        U = r² e^{λ/2} ξ_r ,   V = δP/(ε+P)               (fluid variables)

        dU/dr = (g/c_s²) U + e^{λ/2}[ ℓ(ℓ+1) e^{ν}/ω² − r²/c_s² ] V          (12a)
        dV/dr = e^{λ/2−ν}[ (ω² − N²)/r² ] U + g Δ(c⁻²) V                      (12b)

    with
        g       = −(dp/dr)/(ε+p) = ν'/2        local gravity (>0; from TOV)
        Δ(c⁻²)  = 1/c_e² − 1/c_s²              the Schwarzschild-discriminant factor
        N²      = g² Δ(c⁻²) e^{ν−λ}            relativistic Brunt–Väisälä frequency (13)
        c_s²    = adiabatic sound speed (perturbation index Γ, frozen composition)
        c_e²    = dp/dε ALONG the star = equilibrium sound speed (structure index Γ_struct)

    Metric ds² = −e^{ν}dt² + e^{λ}dr² + r²dΩ²,  e^{λ}=(1−2m/r)^{-1},
    ν' = 2(m+4πr³p)/(r(r−2m)) (TOV, ν increasing outward) — IDENTICAL convention
    to NonRadialModes. (Refs. write the line element −e^{−ν}; the load-bearing
    object is g=ν'/2 and the *relative* ν,λ factors in (12), which we pin to the
    barotropic f-mode anchor — see validation note below.)

    SIGN of the physics (the easy-to-flip claim, validated numerically):
        Δ(c⁻²) > 0  ⇔  c_s² > c_e²  ⇔  Γ > Γ_struct  ⇒  N² > 0
        ⇒ STABLE, real, oscillatory g-modes (a tower BELOW the f-mode, g₁>g₂>…→0).
        Δ(c⁻²) < 0 (Γ < Γ_struct) ⇒ N² < 0 ⇒ convectively UNSTABLE (imaginary g).
        Γ = Γ_struct ⇒ N²≡0 ⇒ NO g-modes (recovers NonRadialModes f/p).

    METHOD: RK4 shoot from the regular centre, Brent root-find ω² on the surface
    Lagrangian Δp=0 condition (V finite ⇒ U(R)=0 is the cleaner surface BC here;
    we use the surface mismatch D(ω²)=U(rf)). Reuses the NonRadialModes/Numerics
    template. The buoyancy term makes the W-eq source ω-DEPENDENT in a way that
    creates the low-frequency g-branch (the (ω²−N²)/r² factor flips sign across
    N², trapping g-modes where ω²<N²), on top of the f/p branch.

    EOS-GENERAL: works on any IdealGasStar (with its cs²,cn² and the analytic
    equilibrium c_e² of the structure polytrope) or any TOV+barotrope where N²≈0.
=#
module GravityModes

using ..EquationOfState
using ..TOV: IdealGasStar, solve_tov_idealgas, TOVStar, solve_tov
using ..Numerics: brent
using ..Units: Msun_to_km, kHz_to_km

export gmode_spectrum, brunt_vaisala, freq_kHz_from_omega2_g

# ω² (geometric km⁻²) → kHz.  The ideal-gas star is built in km units (ρc in
# km⁻²), so the geometric length unit IS 1 km ⇒ Lunit_km=1.0.
freq_kHz_from_omega2_g(ω2::Real; Lunit_km::Real=1.0) =
    sqrt(abs(ω2)) / (Lunit_km * 2π * kHz_to_km)

@inline function _lin(xs, ys, x)              # linear interp on ascending xs
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j]) / (xs[j+1] - xs[j])
    ys[j] + t * (ys[j+1] - ys[j])
end

# O(1) linear interpolation on a UNIFORM grid x = x0 .+ (0:n-1)*dx (ascending).
@inline function _lin_uniform(ys, x0, dx, n, x)
    s = (x - x0) / dx
    s ≤ 0 && return ys[1]
    s ≥ n - 1 && return ys[n]
    j = unsafe_trunc(Int, s) + 1
    t = s - (j - 1)
    ys[j] + t * (ys[j+1] - ys[j])
end

# --------------------------------------------------------------------------- #
#   Background profiles needed by the g-mode system, packed for interpolation
# --------------------------------------------------------------------------- #
"""
Internal background bundle on the interior grid (ascending r):
fields r, ν, λ, g (=ν'/2), cs2 (adiabatic), ce2 (equilibrium), N2, plus R.
`ν` is shifted so e^{ν(R)}=1−2M/R (Schwarzschild match), exactly as in TOV.
"""
struct GBackground
    r::Vector{Float64}
    ν::Vector{Float64}
    λ::Vector{Float64}     # e^{λ}=(1−2m/r)^{-1} ⇒ λ=−log(1−2m/r)
    g::Vector{Float64}     # local gravity ν'/2
    cs2::Vector{Float64}   # adiabatic sound speed²
    ce2::Vector{Float64}   # equilibrium sound speed² (structure dp/dε)
    N2::Vector{Float64}    # Brunt–Väisälä N²
    R::Float64
    M::Float64
end

# Equilibrium (structure-polytrope) sound speed²: c_e² = dp/dε along p=Kρ^Γs,
# ε=ρ+p/(Γs−1).  Analytic:  dp/dρ = ΓsKρ^{Γs−1}=:x ;  dε/dρ = 1 + x/(Γs−1) ;
#   c_e² = x / (1 + x/(Γs−1)).
@inline function _ce2_struct(ρ::Real, K::Real, Γs::Real)
    x = Γs * K * ρ^(Γs - 1)
    return x / (1 + x / (Γs - 1))
end

"""
    brunt_vaisala(star::IdealGasStar) -> (r, N2, ce2, cs2)

Relativistic Brunt–Väisälä frequency N²(r)=g²(1/c_e²−1/c_s²)e^{ν−λ} on the
finite-temperature ideal-gas star, with g=−(dp/dr)/(ε+p)=ν'/2 the local gravity,
c_s² the ADIABATIC sound speed (index Γ, frozen composition) and c_e²=dp/dε the
EQUILIBRIUM sound speed of the hydrostatic structure (index Γ_struct). Returns
the interior arrays (surface point dropped where p→0). N²>0 ⇔ c_s²>c_e² ⇔
Γ>Γ_struct ⇔ convectively STABLE g-modes; N²≡0 for the barotropic Γ=Γ_struct star.
"""
function brunt_vaisala(star::IdealGasStar)
    ni = findlast(>(0), star.p)
    ni === nothing && error("star has no interior pressure")
    r = star.r[1:ni]; m = star.m[1:ni]; p = star.p[1:ni]; ε = star.ε[1:ni]
    ν = star.ν[1:ni]; ρ = star.ρ[1:ni]; cs2 = star.cs2[1:ni]
    K = star.K; Γs = star.Γ_struct
    λ  = @. -log(1 - 2m / r)
    νp = @. 2 * (m + 4π * r^3 * p) / (r * (r - 2m))   # ν' (TOV)
    g  = @. νp / 2                                     # local gravity = −p'/(ε+p)
    ce2 = [_ce2_struct(ρi, K, Γs) for ρi in ρ]
    Δinv = @. 1.0 / ce2 - 1.0 / cs2                   # Δ(c⁻²) = 1/c_e²−1/c_s²
    N2 = @. g^2 * Δinv * exp(ν - λ)
    return r, N2, ce2, cs2
end

# --------------------------------------------------------------------------- #
#   Uniform-grid background bundle for the FAST eigensolver inner loop.
#   All fields resampled onto x = r0 .+ (0:Nbg-1)*dx (ascending). O(1) lookups.
# --------------------------------------------------------------------------- #
struct GBgUniform
    r0::Float64
    dx::Float64
    n::Int
    ν::Vector{Float64}
    λ::Vector{Float64}
    g::Vector{Float64}
    cs2::Vector{Float64}
    ce2::Vector{Float64}
    N2::Vector{Float64}
end

@inline function _bg_at(b::GBgUniform, r)
    ν   = _lin_uniform(b.ν,   b.r0, b.dx, b.n, r)
    λ   = _lin_uniform(b.λ,   b.r0, b.dx, b.n, r)
    g   = _lin_uniform(b.g,   b.r0, b.dx, b.n, r)
    cs2 = max(_lin_uniform(b.cs2, b.r0, b.dx, b.n, r), 1e-12)
    ce2 = max(_lin_uniform(b.ce2, b.r0, b.dx, b.n, r), 1e-12)
    N2  = _lin_uniform(b.N2,  b.r0, b.dx, b.n, r)
    return ν, λ, g, cs2, ce2, N2
end

# RHS of the (U,V) g-mode system (Jaikumar/Shirke Eq. 12) at radius r.
@inline function _uv_rhs(b::GBgUniform, r, U, V, l, ω2)
    ν, λ, g, cs2, ce2, N2 = _bg_at(b, r)
    eλh = exp(λ / 2)
    enu = exp(ν)
    Δinv = 1.0 / ce2 - 1.0 / cs2
    dU = (g / cs2) * U + eλh * (l * (l + 1) * enu / ω2 - r^2 / cs2) * V
    dV = eλh / enu * (ω2 - N2) / r^2 * U + g * Δinv * V
    return dU, dV
end

# regular-centre seed for U,V at r0 (ξ_r∼r^{l−1} ⇒ U=r²e^{λ/2}ξ_r∼r^{l+1}; V from
# the dU-eq dominant balance so the centre is regular).
@inline function _seed_uv(b::GBgUniform, l, ω2, r0)
    ν0, λ0, g0, cs20, _, _ = _bg_at(b, r0)
    U = r0^(l + 1)
    eλh0 = exp(λ0 / 2)
    brack = eλh0 * (l * (l + 1) * exp(ν0) / ω2 - r0^2 / cs20)
    V = brack != 0 ? ((l + 1) * U / r0 - (g0 / cs20) * U) / brack : 0.0
    return U, V
end

# integrate centre→rf with RK4; return surface value U(rf) (the surface mismatch
# D(ω²) — Lagrangian Δp finite ⇒ at the free surface ξ_r regular ⇒ U(R) drives it).
function _shoot_uv(b::GBgUniform, l, ω2, r0, rf, nstep)
    h = (rf - r0) / nstep
    U, V = _seed_uv(b, l, ω2, r0)
    r = r0
    for _ in 1:nstep
        k1U, k1V = _uv_rhs(b, r,        U,           V,           l, ω2)
        k2U, k2V = _uv_rhs(b, r + h/2,  U + h/2*k1U, V + h/2*k1V, l, ω2)
        k3U, k3V = _uv_rhs(b, r + h/2,  U + h/2*k2U, V + h/2*k2V, l, ω2)
        k4U, k4V = _uv_rhs(b, r + h,    U + h*k3U,   V + h*k3V,   l, ω2)
        U += h/6 * (k1U + 2k2U + 2k3U + k4U)
        V += h/6 * (k1V + 2k2V + 2k3V + k4V)
        r += h
        a = abs(U) + abs(V)
        if a > 1e150
            U /= a; V /= a
        end
    end
    return U
end

# count radial nodes (sign changes) of ξ_r ∝ U/(r²e^{λ/2}) for eigenvalue ω². The
# node count classifies the mode family: f-mode has 0 nodes; p_n (ω>ω_f) and g_n
# (ω<ω_f) each have n nodes — the standard Cowling/Sturm node theorem. Nodes are
# counted only in the BULK interior [edge, 1−edge]·(rf−r0): U→0 at both the regular
# centre and the free surface, where round-off introduces SPURIOUS sign flips that
# would corrupt the (sensitive) f-mode 0-node count. A genuine radial node sits in
# the bulk, so trimming the thin boundary layers makes the classifier N-robust.
function _node_count(b::GBgUniform, l, ω2, r0, rf, nstep; edge::Float64=0.03)
    h = (rf - r0) / nstep
    U, V = _seed_uv(b, l, ω2, r0)
    r = r0; nodes = 0; Uprev = U
    rlo = r0 + edge * (rf - r0)
    rhi = rf - edge * (rf - r0)
    for _ in 1:nstep
        k1U, k1V = _uv_rhs(b, r,        U,           V,           l, ω2)
        k2U, k2V = _uv_rhs(b, r + h/2,  U + h/2*k1U, V + h/2*k1V, l, ω2)
        k3U, k3V = _uv_rhs(b, r + h/2,  U + h/2*k2U, V + h/2*k2V, l, ω2)
        k4U, k4V = _uv_rhs(b, r + h,    U + h*k3U,   V + h*k3V,   l, ω2)
        U += h/6 * (k1U + 2k2U + 2k3U + k4U)
        V += h/6 * (k1V + 2k2V + 2k3V + k4V)
        r += h
        a = abs(U) + abs(V)
        if a > 1e150
            U /= a; V /= a; Uprev /= a
        end
        if rlo ≤ r ≤ rhi && U * Uprev < 0
            nodes += 1
        end
        Uprev = U
    end
    return nodes
end

# resample the (non-uniform-safe) interior background onto a uniform grid in r.
function _build_uniform_bg(rt, νt, λt, gt, cs2t, ce2t, N2t, r0, rf, Nbg)
    dx = (rf - r0) / (Nbg - 1)
    ν   = Vector{Float64}(undef, Nbg); λ   = similar(ν); g   = similar(ν)
    cs2 = similar(ν); ce2 = similar(ν); N2  = similar(ν)
    for i in 1:Nbg
        r = r0 + (i - 1) * dx
        ν[i]   = _lin(rt, νt, r);   λ[i]   = _lin(rt, λt, r);  g[i]   = _lin(rt, gt, r)
        cs2[i] = _lin(rt, cs2t, r); ce2[i] = _lin(rt, ce2t, r); N2[i] = _lin(rt, N2t, r)
    end
    return GBgUniform(r0, dx, Nbg, ν, λ, g, cs2, ce2, N2)
end

"""
    gmode_spectrum(star::IdealGasStar; l=2, N=4000, Nbg=3000, Lunit_km=1.0,
                   ω2lo=1e-4, ω2hi=0.25, nscan=1600, nmodes_g=4, nmodes_p=4)
        -> NamedTuple

Relativistic Cowling polar (ℓ≥2) spectrum of a finite-temperature ideal-gas TOV
star INCLUDING buoyancy. Solves the (U,V) system by RK4 shooting + Brent root-find
on the surface mismatch U(rf)=0. The smooth background is resampled onto a uniform
grid of `Nbg` nodes for O(1) interpolation, and the eigensolver takes `N` RK4 steps
(both decoupled from the much-finer TOV grid). Splits the eigenvalues at the f-mode:

    g-modes  : ω² below the f-mode (the buoyancy tower, present iff N²>0 somewhere)
    f-mode   : the lowest mode with NO radial node interior to it (fundamental)
    p-modes  : ω² above the f-mode

Returns `(; g_freqs_kHz, f_freq, p_freqs, all_freqs_kHz, ω2_all, N2_max, R, M, has_gmodes)`.
Frequencies in kHz. Because the ideal-gas star is in km units, `Lunit_km=1.0`.
"""
function gmode_spectrum(star::IdealGasStar;
        l::Int=2, N::Int=4000, Nbg::Int=3000, Lunit_km::Float64=1.0,
        ω2lo::Float64=1e-4, ω2hi::Float64=0.25, nscan::Int=1600,
        nmodes_g::Int=4, nmodes_p::Int=4)

    ni = findlast(>(0), star.p)
    R = star.R; M = star.M
    r, N2, ce2, cs2 = brunt_vaisala(star)
    rt = star.r[1:ni]; mt = star.m[1:ni]; pt = star.p[1:ni]
    νt = star.ν[1:ni]
    λt = @. -log(1 - 2mt / rt)
    νp = @. 2 * (mt + 4π * rt^3 * pt) / (rt * (rt - 2mt))
    gt = @. νp / 2

    r0 = rt[1] + (R - rt[1]) * 1e-4
    rf = R * (1 - 1e-3)
    bg = _build_uniform_bg(rt, νt, λt, gt, cs2, ce2, N2, r0, rf, Nbg)
    D(ω2) = _shoot_uv(bg, l, ω2, r0, rf, N)

    grid = exp.(range(log(ω2lo), log(ω2hi); length=nscan))  # log scan: resolves g-tower
    ω2s = Float64[]
    Dprev = D(first(grid)); ω2prev = first(grid)
    for ω2 in Iterators.drop(grid, 1)
        Dc = D(ω2)
        if isfinite(Dc) && isfinite(Dprev) && Dc * Dprev < 0
            res = brent(D, ω2prev, ω2; xtol=1e-13)
            res.converged && push!(ω2s, res.root)
        end
        Dprev = Dc; ω2prev = ω2
    end
    sort!(ω2s)
    # dedupe near-identical roots (log scan can double-bracket)
    uniq = Float64[]
    for w in ω2s
        if isempty(uniq) || abs(w - uniq[end]) > 1e-9 * max(w, 1e-9)
            push!(uniq, w)
        end
    end
    ω2s = uniq

    N2max = maximum(N2)
    freqs = [freq_kHz_from_omega2_g(w; Lunit_km=Lunit_km) for w in ω2s]

    # classify by RADIAL NODE COUNT (Sturm/Cowling node theorem): the f-mode is the
    # mode with the FEWEST radial nodes; modes ABOVE it (ω²>ω²_f) are p-modes (node
    # count rising upward); modes BELOW it (ω²<ω²_f) are g-modes (node count rising
    # DOWNward, accumulating at ω²→0). Buoyancy (N²>0) is what populates the g-branch.
    # (argmin is robust to a constant offset in the discrete node count; ties resolve
    # to the lowest-frequency member, which is the fundamental.)
    nodes = [_node_count(bg, l, w, r0, rf, N) for w in ω2s]
    fpos = argmin(nodes)
    ω2_f = ω2s[fpos]
    f_freq = freqs[fpos]
    gidx = findall(w -> w < ω2_f, ω2s)              # below the f-mode ⇒ g-branch
    pidx = findall(w -> w > ω2_f, ω2s)              # above the f-mode ⇒ p-branch
    g_freqs = sort(freqs[gidx]; rev=true)          # g₁ > g₂ > … (descending toward 0)
    p_freqs = sort(freqs[pidx])                      # p₁ < p₂ < … (ascending)
    has_g = N2max > 1e-10 && !isempty(gidx)

    return (; g_freqs_kHz = g_freqs[1:min(nmodes_g, length(g_freqs))],
              f_freq = f_freq,
              p_freqs = p_freqs[1:min(nmodes_p, length(p_freqs))],
              all_freqs_kHz = freqs,
              ω2_all = ω2s,
              nodes = nodes,
              N2_max = N2max,
              R = R, M = M,
              has_gmodes = has_g)
end

"""
    gmode_spectrum(; Γ, Γ_struct, K, ρc, kwargs...) -> NamedTuple

Convenience wrapper: build the ideal-gas star with `solve_tov_idealgas` then run
`gmode_spectrum` on it. `Γ` is the adiabatic (perturbation) index, `Γ_struct` the
hydrostatic structure index; Γ>Γ_struct ⇒ N²>0 ⇒ stable g-modes.
"""
function gmode_spectrum(; Γ::Real, Γ_struct::Real=Γ, K::Real=100.0, ρc::Real,
                        h::Float64=2e-4, kwargs...)
    star = solve_tov_idealgas(; Γ=Γ, Γ_struct=Γ_struct, K=K, ρc=ρc, h=h)
    return gmode_spectrum(star; kwargs...)
end

end # module GravityModes
