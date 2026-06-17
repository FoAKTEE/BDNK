#=
    Recovery — the conservative <-> primitive map (the recurring "wall" on the
    trunk, cf. bdnk_hmns_plan.tex §Plan).

    Special-relativistic / flat-space form (the Cowling-core building block;
    curved-space lapse/shift factors enter at STAGE 2). Two ideal-fluid
    inversions, both reduced to a single 1-D root find in the pressure:

      barotropic (e, v)   <-> (E = T^{tt},  S = T^{tx})
      general    (ρ,v,ϵ)  <-> (D = ρW,  S = ρhW²v,  τ = ρhW²-p-D)

    BDNK twist: the conserved densities carry the first-derivative dissipative
    corrections, so the inversion is the *ideal* inversion shifted by frozen
    source terms (δE, δS) evaluated from grid data (Pandya–Most–Pretorius
    2209.09265). In the ideal limit δ→0 the BDNK recovery reduces *exactly* to
    the ideal inversion, so the round-trip closes to machine precision — the
    STEP 0 gate (≤ 1e-10 in smooth regions).
=#
module Recovery

using ..Numerics
using ..EquationOfState

export prim2cons_barotropic, cons2prim_barotropic,
       prim2cons_general, cons2prim_general,
       cons2prim_bdnk_barotropic, ConsBarotropic, PrimBarotropic,
       lorentz_W

@inline lorentz_W(v::Real) = 1 / sqrt(1 - v^2)

# ---- barotropic ideal inversion -------------------------------------------
"""
    prim2cons_barotropic(eos, e, v) -> (E, S, p)

Perfect-fluid conserved densities from barotropic primitives. T^{μν} =
(e+p)u^μu^ν + p η^{μν}, η = diag(-1,1,…). E = T^{tt}, S = T^{tx}.
"""
function prim2cons_barotropic(eos::BarotropicEOS, e::Real, v::Real)
    p  = pressure(eos, e)
    W2 = 1 / (1 - v^2)
    hh = (e + p) * W2                     # (e+p) W²
    E  = hh - p
    S  = hh * v
    return E, S, p
end

"""
    cons2prim_barotropic(eos, E, S) -> (e, v, p, info)

Invert (E,S) -> (e,v). Eliminate v = S/(E+p), e = E - S²/(E+p); solve the scalar
residual g(p) = pressure(eos, e(p)) - p = 0 by Brent on a guaranteed bracket.
For the conformal limit p = e/3 this reproduces the closed form
e = -E + √(4E² - 3S²) (cf. 1D_conformal_bdnk/solver.c).
"""
function cons2prim_barotropic(eos::BarotropicEOS, E::Real, S::Real)
    # bracket: need E + p > |S| (so v² < 1) and p ≥ 0.
    p_lo = max(0.0, abs(S) - E) + 1e-300
    e_hi = E                                  # p>0 ⇒ e = E - S²/(E+p) < E
    p_hi = pressure(eos, e_hi) + abs(S)       # generous upper bracket
    g(p) = pressure(eos, E - S^2/(E + p)) - p
    # widen p_hi until sign change (monotone problem, terminates fast)
    while g(p_lo) * g(p_hi) > 0 && p_hi < 1e300
        p_hi *= 2
    end
    r = brent(g, p_lo, p_hi; xtol=1e-15)
    p = r.root
    v = S / (E + p)
    e = E - S^2 / (E + p)
    return e, v, p, r
end

# ---- general (finite-T) ideal inversion -----------------------------------
"""
    prim2cons_general(eos, ρ, v, ϵ) -> (D, S, τ, p)

Perfect-fluid conserved densities from (ρ, v, ϵ). T^{μν} = ρh u^μu^ν + p η^{μν},
h = 1 + ϵ + p/ρ. D = ρW, S = ρhW²v, τ = ρhW² - p - D.
"""
function prim2cons_general(eos::GeneralEOS, ρ::Real, v::Real, ϵ::Real)
    p  = pressure(eos, ρ, ϵ)
    W  = lorentz_W(v)
    h  = 1 + ϵ + p/ρ
    ρhW2 = ρ * h * W^2
    D = ρ * W
    S = ρhW2 * v
    τ = ρhW2 - p - D
    return D, S, τ, p
end

"""
    cons2prim_general(eos, D, S, τ) -> (ρ, v, ϵ, p, info)

Invert (D,S,τ) -> (ρ,v,ϵ). With E = τ + D: v = S/(E+p), W = 1/√(1-v²),
ρ = D/W, h = (E+p)/(ρW²), ϵ = h - 1 - p/ρ; solve f(p) = pressure(eos,ρ,ϵ) - p.
"""
function cons2prim_general(eos::GeneralEOS, D::Real, S::Real, τ::Real)
    E = τ + D
    function prim_of_p(p)
        v = S / (E + p)
        W = lorentz_W(v)
        ρ = D / W
        h = (E + p) / (ρ * W^2)
        ϵ = h - 1 - p/ρ
        return ρ, v, ϵ, W
    end
    function f(p)
        ρ, v, ϵ, _ = prim_of_p(p)
        return pressure(eos, ρ, ϵ) - p
    end
    p_lo = max(0.0, abs(S) - E) + 1e-300
    p_hi = max(E, D) + abs(S) + 1.0
    while f(p_lo) * f(p_hi) > 0 && p_hi < 1e300
        p_hi *= 2
    end
    r = brent(f, p_lo, p_hi; xtol=1e-15)
    p = r.root
    ρ, v, ϵ, _ = prim_of_p(p)
    return ρ, v, ϵ, p, r
end

# ---- BDNK gradient-frozen recovery ----------------------------------------
"""
    cons2prim_bdnk_barotropic(eos, E, S, δE, δS) -> (e, v, p, info)

BDNK recovery: the conserved densities (E,S) = (E_PF + δE, S_PF + δS) carry the
first-order dissipative corrections (δE, δS) which depend on gradients evaluated
and *frozen* from grid data. Subtract them and run the ideal inversion. In the
ideal/equilibrium limit (δE,δS)→0 this is exactly `cons2prim_barotropic`, so the
round-trip closes to machine precision.
"""
function cons2prim_bdnk_barotropic(eos::BarotropicEOS, E::Real, S::Real,
                                   δE::Real, δS::Real)
    return cons2prim_barotropic(eos, E - δE, S - δS)
end

# convenience aggregates for tests / drivers
struct PrimBarotropic; e::Float64; v::Float64; p::Float64; end
struct ConsBarotropic; E::Float64; S::Float64; end

end # module Recovery
