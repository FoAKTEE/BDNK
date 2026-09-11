#=
    PiecewisePolytrope — GENERAL / REALISTIC cold neutron-star EOS via the
    standard Read–Lackey–Owen–Friedman (2009) PIECEWISE-POLYTROPE
    parametrization [arXiv:0812.2163, Phys. Rev. D 79, 124032 (2009)].

    A fixed low-density CRUST (a 4-piece polytropic fit to the SLy/SLy4 crust,
    Read et al. Table II) is joined to N = 3 high-density polytropic segments

        p = K_i ρ^{Γ_i}   on   [ρ_{i-1}, ρ_i]

    with the standard FIXED dividing densities ρ1 = 10^14.7 g/cm³ and
    ρ2 = 10^15.0 g/cm³ (Read et al. §III). Each candidate EOS is then specified
    by FOUR numbers: (log10 p1, Γ1, Γ2, Γ3), where p1 = p(ρ1) in dyn/cm².
    Continuity of p fixes the K_i and the matching density ρ0 between the crust
    and the first high-density piece; the total energy density follows the first
    law,  ε_i = (1+a_i) ρ + K_i ρ^{Γ_i}/(Γ_i−1),  with a_i fixed by continuity
    of ε across segments (Read et al. Eqs. (5)-(7)):

        a_i = a_{i-1} + (n_{i-1} − n_i) p(ρ_i)/ρ_i,   n_i = 1/(Γ_i − 1),  a_0 = 0.

    The whole construction is done in GEOMETRIC units (G=c=1, lengths km, ε & p
    in km⁻²) using the Units module conversions — the CGS table parameters are
    converted, NO conversion factors are hard-coded here. Following the existing
    `isentropic_idealgas` pattern, the analytic ρ→(p,ε,c_s²) curve is TABULATED
    into a `TabulatedBarotrope`, so the realistic EOS plugs straight into
    `solve_tov` and every perturbation/QNM solver.

    ── CITED reference values ──────────────────────────────────────────────────
    Crust (Read et al. Table II, SLy4 4-piece fit), taken VERBATIM from the
    reference implementation in LALSimulation (LALSimNeutronStarEOSPiecewise-
    Polytrope.c), which states "In the paper PRD 79, 124032 (2009) which used
    cgs, the k_i values in Table II should be multiplied by c²":
        ρ boundaries [g/cm³]: 0, 2.44033979e7, 3.78358138e11, 2.62780487e12
        K (CGS, p[dyn/cm²]=Kρ^Γ): 6.11252036792443e12, 9.54352947022931e14,
                                  4.787640050002652e22, 3.593885515256112e13
        Γ: 1.58424999, 1.28732904, 0.62223344, 1.35692395

    High-density parameters (Read et al. 2009 Table III, least-squares fits):
        EOS    log10 p1   Γ1     Γ2     Γ3      M_max[M⊙]  R_1.4[km]
        SLy    34.384     3.005  2.988  2.851   2.049      11.736
        APR4   34.269     2.830  3.445  3.348   ~2.20      ~11.4
        H4     34.669     2.909  2.246  2.144   2.032      13.774
        MS1    34.858     3.224  3.033  1.325   2.767      14.918
    (SLy, H4, MS1 parameters + M_max/R_1.4 reproduced from Lackey & Wade 2015,
     arXiv:1410.8866 Table I, "values that minimize the least-squares residual
     defined in Read et al. (2009)"; APR4 parameters from Read et al. Table III
     and the standard PP literature; APR4 M_max/R_1.4 are the well-known APR
     values M_max≈2.2 M⊙, R_1.4≈11.4 km. Crust + dividing densities cross-checked
     against LALSimulation.)
=#
module PiecewisePolytrope

using ..Units: gram_per_cm3_to_km_minus2, dyne_per_cm2_to_km_minus2
using ..EquationOfState: TabulatedBarotrope

export piecewise_polytrope

# ── CGS crust (Read et al. Table II, SLy4 fit; see module docstring) ──────────
const _CRUST_RHO_CGS = (0.0, 2.44033979e7, 3.78358138e11, 2.62780487e12)  # g/cm³
const _CRUST_K_CGS   = (6.11252036792443e12, 9.54352947022931e14,
                        4.787640050002652e22, 3.593885515256112e13)        # CGS
const _CRUST_GAMMA   = (1.58424999, 1.28732904, 0.62223344, 1.35692395)

# Fixed high-density dividing densities (Read et al. §III)
const _RHO1_CGS = 10.0^14.7   # g/cm³  (reference density for p1)
const _RHO2_CGS = 10.0^15.0   # g/cm³

# Published Read et al. 2009 Table III parameters (log10 p1 [dyn/cm²], Γ1,Γ2,Γ3)
const _PRESETS = Dict(
    :SLy  => (logp1 = 34.384, Γ1 = 3.005, Γ2 = 2.988, Γ3 = 2.851),
    :APR4 => (logp1 = 34.269, Γ1 = 2.830, Γ2 = 3.445, Γ3 = 3.348),
    :H4   => (logp1 = 34.669, Γ1 = 2.909, Γ2 = 2.246, Γ3 = 2.144),
    :MS1  => (logp1 = 34.858, Γ1 = 3.224, Γ2 = 3.033, Γ3 = 1.325),
)

# Published reference M_max [M⊙] and R_1.4 [km] for validation (see docstring).
const PUBLISHED = Dict(
    :SLy  => (M_max = 2.049, R_14 = 11.736),
    :APR4 => (M_max = 2.20,  R_14 = 11.40),
    :H4   => (M_max = 2.032, R_14 = 13.774),
    :MS1  => (M_max = 2.767, R_14 = 14.918),
)

"""
    _segments(logp1, Γ1, Γ2, Γ3) -> (ρ_geo, K_geo, Γ, a, n)

Assemble the full set of polytropic segments (crust + 3 high-density pieces) in
GEOMETRIC units. Returns, for the M segments: the starting density `ρ_geo[i]`
[km⁻²], polytropic constant `K_geo[i]` (geometric, p=Kρ^Γ), index `Γ[i]`, and
the first-law constants `a[i]` and `n[i]=1/(Γ[i]−1)` enforcing ε-continuity.

If the crust–core matching density ρ0 falls outside the last crust / first core
piece, a single bridging polytrope is inserted between ρ=5e15 and 1e16 g/cm³
(exactly the LALSimulation fallback), so the construction is robust for stiff/
soft p1 the same way the reference is.
"""
function _segments(logp1::Float64, Γ1::Float64, Γ2::Float64, Γ3::Float64)
    # crust in geometric units: convert ρ (g/cm³) and the polytrope so that
    # p[km⁻²] = K_geo ρ_geo^Γ. With p_cgs = K_cgs ρ_cgs^Γ and
    # p_geo = p_cgs·(dyn/cm²→km⁻²), ρ_geo = ρ_cgs·(g/cm³→km⁻²):
    #   K_geo = K_cgs · (dyn→km⁻²) / (g→km⁻²)^Γ
    gcgs = gram_per_cm3_to_km_minus2
    pcgs = dyne_per_cm2_to_km_minus2
    Kgeo_cgs(Kc, Γ) = Kc * pcgs / gcgs^Γ

    ρ  = Float64[ρc * gcgs for ρc in _CRUST_RHO_CGS]              # 4 crust ρ (geo)
    K  = Float64[Kgeo_cgs(_CRUST_K_CGS[i], _CRUST_GAMMA[i]) for i in 1:4]
    Γ  = collect(Float64, _CRUST_GAMMA)

    # high-density pieces
    ρ1 = _RHO1_CGS * gcgs
    ρ2 = _RHO2_CGS * gcgs
    p1 = 10.0^logp1 * pcgs                                        # p(ρ1), geometric
    k1 = p1 / ρ1^Γ1
    k2 = p1 / ρ1^Γ2
    k3 = k2 * ρ2^(Γ2 - Γ3)

    # crust–core matching density ρ0 (continuity p_crust = p_core at ρ0)
    Kc_last, Γc_last = K[4], Γ[4]
    ρ0 = (Kc_last / k1)^(1.0 / (Γ1 - Γc_last))

    if ρ[4] < ρ0 < ρ1
        # 7 segments: 4 crust + (Γ1 from ρ0) + (Γ2 from ρ1) + (Γ3 from ρ2)
        append!(ρ, (ρ0, ρ1, ρ2))
        append!(K, (k1, k2, k3))
        append!(Γ, (Γ1, Γ2, Γ3))
    else
        # insert a bridging polytrope between 5e15 and 1e16 g/cm³ (LAL fallback)
        ρj1 = 5.0e15 * gcgs
        ρj2 = 1.0e16 * gcgs
        pj1 = Kc_last * ρj1^Γc_last
        pj2 = k1 * ρj2^Γ1
        Γj  = log(pj2 / pj1) / log(ρj2 / ρj1)
        kj  = pj1 / ρj1^Γj
        append!(ρ, (ρj1, ρj2, ρ1, ρ2))
        append!(K, (kj, k1, k2, k3))
        append!(Γ, (Γj, Γ1, Γ2, Γ3))
    end

    # first-law constants a_i (ε-continuity) and n_i
    M = length(ρ)
    n = [1.0 / (Γ[i] - 1.0) for i in 1:M]
    a = zeros(M)
    for i in 2:M
        p_i = K[i] * ρ[i]^Γ[i]
        a[i] = a[i-1] + (n[i-1] - n[i]) * p_i / ρ[i]
    end
    return ρ, K, Γ, a, n
end

"""
    piecewise_polytrope(; logp1, Γ1, Γ2, Γ3, name="",
                        ρ_lo=1e3, ρ_hi=1e16, N=2000) -> TabulatedBarotrope

Build a Read et al. (2009) piecewise-polytrope `BarotropicEOS` from the four
high-density parameters (`logp1` = log10 of p(ρ1) in dyn/cm², `Γ1,Γ2,Γ3`),
joined to the fixed SLy crust. Densities `ρ_lo`/`ρ_hi` are the tabulation range
in g/cm³ (default 10³–10¹⁶, comfortably bracketing any NS); `N` is the number
of log-spaced rest-mass-density nodes.

The analytic ρ→(p,ε) curve along the segments is tabulated as a barotrope:
the table stores p(ε) and c_s² = dp/dε = Γ K ρ^Γ / (ε + p) at each node, exactly
as the reference EOS interface expects. Returns a `TabulatedBarotrope` in
geometric units (ε, p in km⁻²).
"""
function piecewise_polytrope(; logp1::Real, Γ1::Real, Γ2::Real, Γ3::Real,
                             name::AbstractString="",
                             ρ_lo::Real=1e3, ρ_hi::Real=1e16, N::Int=2000)
    ρseg, Kseg, Γseg, aseg, nseg = _segments(float(logp1), float(Γ1), float(Γ2), float(Γ3))
    gcgs = gram_per_cm3_to_km_minus2

    # locate the segment index for a geometric rest-mass density
    @inline function _seg(ρ)
        i = 1
        @inbounds for j in eachindex(ρseg)
            ρseg[j] ≤ ρ && (i = j)
        end
        return i
    end

    ρ_lo_g = ρ_lo * gcgs
    ρ_hi_g = ρ_hi * gcgs
    ρgrid  = exp.(range(log(ρ_lo_g), log(ρ_hi_g); length=N))      # ascending ⇒ ε ascending

    p   = Vector{Float64}(undef, N)
    e   = Vector{Float64}(undef, N)
    cs2 = Vector{Float64}(undef, N)
    @inbounds for j in 1:N
        ρ = ρgrid[j]
        i = _seg(ρ)
        Γ = Γseg[i]; K = Kseg[i]; a = aseg[i]; n = nseg[i]
        pj = K * ρ^Γ
        ej = (1.0 + a) * ρ + n * pj                              # ε = (1+a)ρ + p/(Γ-1)
        p[j]   = pj
        e[j]   = ej
        cs2[j] = (Γ * pj) / (ej + pj)                            # dp/dε = Γp/(ε+p)
    end
    return TabulatedBarotrope(log.(e), p, cs2)
end

"""
    piecewise_polytrope(preset::Symbol; kwargs...) -> TabulatedBarotrope

Named Read et al. (2009) Table-III preset: `:SLy`, `:APR4`, `:H4`, `:MS1`.
"""
function piecewise_polytrope(preset::Symbol; kwargs...)
    haskey(_PRESETS, preset) ||
        error("unknown piecewise-polytrope preset $preset; choose from $(keys(_PRESETS))")
    p = _PRESETS[preset]
    return piecewise_polytrope(; logp1=p.logp1, Γ1=p.Γ1, Γ2=p.Γ2, Γ3=p.Γ3,
                               name=String(preset), kwargs...)
end

end # module PiecewisePolytrope
