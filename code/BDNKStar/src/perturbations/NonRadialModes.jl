#=
    NonRadialModes — relativistic Cowling NON-radial (polar, ℓ≥2) stellar
    pulsations: the f / p modes (STAGE 3, the eigensolver foundation for the
    3+1D non-radial spectrum).

    Frozen-metric (Cowling) fluid perturbations on a TOV background, McDermott–
    Van Horn–Scholl / Sotani W(r),V(r) form. Grounded in arXiv:2412.18569
    Eqs. (11)-(14) (cross-checked vs Counsell et al. 2409.20178 and Yoshida &
    Lee 2002, astro-ph/0210591). Metric ds² = −e^{ν}dt² + e^{λ}dr² + r²dΩ² with
    ν,λ FULL exponents (e^{ν}=−g_tt, e^{λ}=(1−2m/r)^{-1}); the source paper's
    Φ,λ_half map as ν=2Φ, so Φ'=ν'/2.

        dW/dr = (1/c_s²)[ ω² r² e^{λ−ν} V + (ν'/2) W ] − ℓ(ℓ+1) e^{λ/?}... V   (raw)
    implemented with e^{λ}=(1−2m/r)^{-1/2} (proper-length factor) as:
        dW/dr = (1/c_s²)[ ω² r² (e^{λ}/e^{ν}) V + (ν'/2) W ] − ℓ(ℓ+1) e^{λ} V    (11)
        dV/dr =  ν' V − e^{λ} W / r²                                            (12)
    BCs: centre  W∼r^{ℓ+1}, V∼−r^ℓ/ℓ  (13);  surface Δp=0:
        ω² R² (e^{λ}/e^{ν}) V(R) + (ν'/2) W(R) = 0                              (14)
    Here ν' = 2(m+4πr³p)/(r(r−2m)) (TOV), c_s²=dp/dε. ω² appears in the W-eq and
    the surface condition ⇒ the eigenvalue problem is solved by SHOOTING: RK4
    from the regular centre, root-find ω² on the surface mismatch (Numerics.brent).

    Cowling f/p-modes are real (no metric back-reaction ⇒ no gravitational
    damping); the imaginary w-modes need the full GR system (a later stage).
=#
module NonRadialModes

using ..EquationOfState
using ..TOV
using ..Numerics
using ..Units: Msun_to_km, kHz_to_km

export nonradial_cowling_spectrum, freq_kHz_from_omega2

# ω² (geometric, length-unit⁻²) → kHz.  Lunit_km = km-length of the geometric
# unit (Msun_to_km for M⊙=G=c=1 EOS scaling, 1.0 for km units).
freq_kHz_from_omega2(ω2::Real; Lunit_km::Real=Msun_to_km) =
    sqrt(ω2) / (Lunit_km * 2π * kHz_to_km)

@inline function _lin(xs, ys, x)              # linear interp on ascending xs
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j]) / (xs[j+1] - xs[j])
    ys[j] + t * (ys[j+1] - ys[j])
end

# RHS of the (W,V) Cowling system at radius r, trial ω².  `bg(r)→(m,ν,ε)`.
@inline function _wv_rhs(eos, bg, r, W, V, l, ω2)
    m, ν, ε = bg(r)
    p   = pressure(eos, ε)
    cs2 = max(sound_speed2(eos, ε), 1e-14)
    elam = 1.0 / sqrt(1.0 - 2m / r)                       # e^{λ} = (1−2m/r)^{-1/2}
    enu  = exp(ν)                                          # e^{ν}=e^{2Φ}=−g_tt
    νp   = 2.0 * (m + 4π * r^3 * p) / (r * (r - 2m))       # ν' = 2Φ'
    B    = ω2 * r^2 * (elam / enu) * V + 0.5 * νp * W      # Δp bracket (eq 11/14)
    dW   = B / cs2 - l * (l + 1) * elam * V
    dV   = νp * V - elam * W / r^2
    return dW, dV, B
end

# integrate centre→r_f with RK4; return surface mismatch D(ω²)=B(r_f) and W(r_f)
function _shoot(eos, bg, l, ω2, r0, rf, nstep)
    h = (rf - r0) / nstep
    W = r0^(l + 1); V = -r0^l / l                         # regular-centre BC (eq 13)
    Blast = 0.0
    r = r0
    for _ in 1:nstep
        k1W, k1V, _  = _wv_rhs(eos, bg, r,         W,            V,            l, ω2)
        k2W, k2V, _  = _wv_rhs(eos, bg, r + h/2,   W + h/2*k1W,  V + h/2*k1V,  l, ω2)
        k3W, k3V, _  = _wv_rhs(eos, bg, r + h/2,   W + h/2*k2W,  V + h/2*k2V,  l, ω2)
        k4W, k4V, _  = _wv_rhs(eos, bg, r + h,     W + h*k3W,    V + h*k3V,    l, ω2)
        W += h/6 * (k1W + 2k2W + 2k3W + k4W)
        V += h/6 * (k1V + 2k2V + 2k3V + k4V)
        r += h
        # renormalise to avoid overflow for badly-off trial ω² (homogeneous problem)
        a = abs(W) + abs(V)
        if a > 1e150
            W /= a; V /= a
        end
    end
    _, _, Blast = _wv_rhs(eos, bg, rf, W, V, l, ω2)
    return Blast, W
end

"""
    nonradial_cowling_spectrum(eos, εc; l=2, nmodes=5, N=8000, h_tov=2e-4,
                               Lunit_km=Msun_to_km,
                               ω2lo=5e-4, ω2hi=0.12, nscan=1200)
        -> (freqs_kHz, ω²_unit, R)

ℓ≥2 relativistic Cowling f/p-mode spectrum of a TOV star (barotropic EOS).
`εc` = central ENERGY density in the EOS's geometric units. Solves the (W,V)
system by RK4 shooting + Brent root-finding on the surface Δp=0 condition.
Returns the lowest `nmodes` mode frequencies in kHz (f-mode first), the squared
eigenfrequencies, and the stellar radius.
"""
function nonradial_cowling_spectrum(eos::BarotropicEOS, εc::Float64;
        l::Int=2, nmodes::Int=5, N::Int=8000, h_tov::Float64=2e-4,
        Lunit_km::Float64=Msun_to_km,
        ω2lo::Float64=5e-4, ω2hi::Float64=0.12, nscan::Int=1200)

    star = solve_tov(eos, εc; h=h_tov)
    R = star.R
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε
    εf = εc * 1e-9
    bg(r) = (_lin(rt, mt, r), _lin(rt, νt, r), max(_lin(rt, et, r), εf))

    r0 = R * 1e-4                    # regular-centre start
    rf = R * (1 - 1e-3)              # just inside surface (c_s²→0 exactly at R)
    D(ω2) = _shoot(eos, bg, l, ω2, r0, rf, N)[1]

    # scan ω² for sign changes of the surface mismatch, then Brent-refine each
    grid = range(ω2lo, ω2hi; length=nscan)
    ω2s = Float64[]
    Dprev = D(first(grid)); ω2prev = first(grid)
    for ω2 in Iterators.drop(grid, 1)
        Dc = D(ω2)
        if isfinite(Dc) && isfinite(Dprev) && Dc * Dprev < 0
            res = brent(D, ω2prev, ω2; xtol=1e-12)
            res.converged && push!(ω2s, res.root)
        end
        Dprev = Dc; ω2prev = ω2
    end
    sort!(ω2s)
    k = min(nmodes, length(ω2s))
    freqs = [freq_kHz_from_omega2(w; Lunit_km=Lunit_km) for w in ω2s[1:k]]
    return freqs, ω2s[1:k], R
end

end # module NonRadialModes
