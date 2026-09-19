#=
    DGStarFD — the DG / finite-difference HYBRID for the radial (Cowling) TOV star: the scheme
    Deppe et al. (PRD 105, 123031; arXiv:2109.12033) endorse after finding that every classical
    DG limiter they tried fails on a neutron star ("the only limiting strategy we can endorse is
    a discontinuous Galerkin–finite-difference hybrid method"), implemented on the hp grids of
    DGStarHP so that it can be compared limiter-for-limiter with the pure-DG schemes of
    VALIDATION.md §7.9 on identical grids.

    THE IDEA (Dumbser et al. 2014; SpECTRE `src/Evolution/DgSubcell/`). Every element carries two
    equivalent representations of the solution: the DG polynomial on N Legendre–Gauss–Lobatto
    nodes, and M = 2N−1 finite-volume cell averages on a subgrid of the same element. A
    troubled-cell indicator decides, at the start of every step, which one is evolved:

      • untroubled → the high-order DG operator (spectral accuracy where the flow is smooth);
      • troubled  → a second-order finite-volume scheme on the subcells (monotonised-central
        reconstruction of the primitives to the subcell faces, HLL flux), which is robust at a
        discontinuity and positivity-friendly at the stellar surface.

    The two representations are exchanged by the conservative pair of DGSubcell: P (DG → FD, the
    exact cell averages of the DG polynomial) and R (FD → DG, the constrained least squares that
    preserves the element integral), with R∘P = identity. Elements exchange data through the
    SAME common numerical flux at their shared face whichever grid each side is using, so the
    scheme is conservative across DG/FD interfaces; an FD element reconstructs its face state
    from its outermost subcells plus one ghost cell taken from the neighbour (projected to the
    subgrid if the neighbour is running DG), which is SpECTRE's ghost-zone exchange.

    TROUBLED-CELL INDICATOR (a priori, evaluated on the current solution at the start of each
    step, with hysteresis on the way back):
      • physical admissibility — any node/cell at or below the atmosphere cutoff, or with a
        non-positive pressure proxy, forces FD. On this star that alone keeps the surface
        elements permanently on the subgrid, which is the point of the method;
      • Persson spectral-decay indicator on D̃ (α = 4, top mode), SpECTRE's `PerssonTci`;
      • an element on FD returns to DG only if its RECONSTRUCTED DG solution passes both checks
        for `min_clear` consecutive steps (SpECTRE's `FdToDgTci`, `MinimumClearTcis`).
    SpECTRE additionally runs the relaxed discrete maximum principle as an a-posteriori check
    with a step rollback; here the indicator is a priori (its `AlwaysUseSubcells`-style mode),
    and a global step redo at half Δt is the only a-posteriori safety net. DGSubcell provides
    `rdmp_troubled` for the rollback variant.

    WHAT IS SHARED WITH DGStarHP: the grid builder (`hp_grid`, `_build_elems`), the Γ-law
    thermodynamics, the primitive recovery with the Galeazzi-type fixing, the atmosphere and the
    parity projection of the doubled domain. WHAT IS NOT: there is no slope limiter and no
    exponential filter — robustness comes from the subcell scheme alone — and the geometric
    source is written as 2α e^{λ/2} x p rather than α√γ (2p/x), which is algebraically the same
    but regular at x = 0, because a subcell centre does sit on the origin (M is odd).

    The well-balanced residual subtraction of DGStarHP is deliberately absent: the paper's
    scheme lets the star settle to its own numerical equilibrium, and err[D̃] then measures the
    dynamics, exactly as in their Figs. 9 and 12.
=#
module DGStarFD

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon: rho_from_p
using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis
using ..DGSubcell
using ..DGSubcell: SubcellOps, subcell_ops, persson_ratio, persson_threshold, mc_slope
using ..DGStarHP: HPElem, hp_grid, _build_elems, _p_ig, _εpoly, _cs2, _p2c, _λpm, _lin, _c2p_resid
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram

export DGStarFDEngine, DGStarFDState, setup_dgstarfd, evolve_dgstarfd!,
       seed_dgstarfd_radial!, dgstarfd_central_density, dgstarfd_errD,
       dgstarfd_baryon_mass, dgstarfd_fd_fraction, dgstarfd_active_map, dgstarfd_spectrum

# ----------------------------------------------------------------------------------
# engine / state
# ----------------------------------------------------------------------------------
struct DGStarFDEngine
    elems::Vector{HPElem}
    bases::Dict{Int,LGLBasis}
    ops::Dict{Int,SubcellOps}
    K::Int; Ntot::Int; Mtot::Int; kc::Int
    soff::Vector{Int}; sM::Vector{Int}; shx::Vector{Float64}   # subcell offset, count, width per element
    sfoff::Vector{Int}                                          # subcell-face offset per element (M+1 faces)
    Γ::Float64; κ::Float64
    eos::BarotropicEOS
    R::Float64; M::Float64
    ρ_atm::Float64; ρ_cut::Float64; ε_atm::Float64; p_atm::Float64; εfac_max::Float64
    flux::Symbol; entropy_floor::Bool; symmetrize::Bool
    tci_α::Float64; tci_nmodes::Int; min_clear::Int; always_fd::Bool
    cfl::Float64; dxmin::Float64
    # DG node geometry (Ntot)
    x::Vector{Float64}; r::Vector{Float64}; α::Vector{Float64}; elam::Vector{Float64}
    sqrtγ::Vector{Float64}; Φx::Vector{Float64}; λx::Vector{Float64}
    # subcell-centre geometry (Mtot)
    sx::Vector{Float64}; sr::Vector{Float64}; sα::Vector{Float64}; selam::Vector{Float64}
    ssqrtγ::Vector{Float64}; sΦx::Vector{Float64}; sλx::Vector{Float64}
    # subcell-face geometry (Σ(M+1))
    sfx::Vector{Float64}; sfα::Vector{Float64}; sfelam::Vector{Float64}; sfsqrtγ::Vector{Float64}
    # equilibrium (projected TOV) on both grids, for the err[D̃] diagnostic
    Deq::Vector{Float64}; sDeq::Vector{Float64}
    # which grid each element is on: 0 = DG, 1 = FD; and the FD→DG clear counter
    active::Vector{Int}; clear::Vector{Int}
end

mutable struct DGStarFDState
    D::Vector{Float64}; S::Vector{Float64}; τ::Vector{Float64}
    ρ::Vector{Float64}; ε::Vector{Float64}; p::Vector{Float64}; v::Vector{Float64}
    sD::Vector{Float64}; sS::Vector{Float64}; sτ::Vector{Float64}
    sρ::Vector{Float64}; sε::Vector{Float64}; sp::Vector{Float64}; sv::Vector{Float64}
    nfix::Int
end
DGStarFDState(n, m) = DGStarFDState((zeros(n) for _ in 1:7)..., (zeros(m) for _ in 1:7)..., 0)

# ----------------------------------------------------------------------------------
# thermodynamics helpers with explicit geometry (the subcell grid has its own metric values)
# ----------------------------------------------------------------------------------
# densitized state, flux and coordinate characteristic speeds at a point
@inline function _qfg(Γ::Float64, α::Float64, el::Float64, sg::Float64, ρ, ε, v̂)
    p = _p_ig(Γ, ρ, ε); h = 1 + ε + p/ρ
    W = 1/sqrt(1 - clamp(v̂^2, 0.0, 1 - 1e-14)); el2 = sqrt(el); vx = v̂/el2
    D̃ = sg*ρ*W; S̃ = sg*ρ*h*W^2*el2*v̂; τ̃ = sg*(ρ*h*W^2 - p) - D̃
    q = (D̃, S̃, τ̃); f = (α*D̃*vx, α*(S̃*vx + sg*p), α*(τ̃*vx + sg*p*vx))
    cs2 = _cs2(Γ, ρ, ε); λm, λp = _λpm(cs2, v̂); c = α/el2
    return q, f, c*λm, c*λp
end

# geometric source, written without the 1/x that is singular at a subcell centre on the origin:
#   α√γ (2p/x) = α e^{λ/2} x² (2p/x) = 2 α e^{λ/2} x p
@inline function _src(Γ::Float64, α, el, sg, xx, Φ, λp, ρ, ε, p, v̂)
    h = 1 + ε + p/ρ; W2 = 1/(1 - clamp(v̂^2, 0.0, 1 - 1e-14)); el2 = sqrt(el); vx = v̂/el2
    sS = α*sg*((ρ*h*W2*v̂^2 + p)*0.5*λp) + 2*α*el2*xx*p - α*sg*(ρ*h*W2 - p)*Φ
    sτ = -α*sg*ρ*h*W2*vx*Φ
    return sS, sτ
end

@inline function _numflux(flux::Symbol, qL, fL, smL, spL, qR, fR, smR, spR)
    if flux == :hll
        cmin = min(smL, smR, 0.0); cmax = max(spL, spR, 0.0)
        if cmax - cmin ≤ 1e-300
            return (0.5*(fL[1]+fR[1]), 0.5*(fL[2]+fR[2]), 0.5*(fL[3]+fR[3]))
        end
        iv = 1/(cmax - cmin)
        return ((cmax*fL[1] - cmin*fR[1] + cmin*cmax*(qR[1]-qL[1]))*iv,
                (cmax*fL[2] - cmin*fR[2] + cmin*cmax*(qR[2]-qL[2]))*iv,
                (cmax*fL[3] - cmin*fR[3] + cmin*cmax*(qR[3]-qL[3]))*iv)
    else
        a = max(abs(smL), abs(spL), abs(smR), abs(spR))
        return (0.5*(fL[1]+fR[1]) - 0.5*a*(qR[1]-qL[1]),
                0.5*(fL[2]+fR[2]) - 0.5*a*(qR[2]-qL[2]),
                0.5*(fL[3]+fR[3]) - 0.5*a*(qR[3]-qL[3]))
    end
end

# Γ-law primitive recovery with the Galeazzi-type fixing (as DGStarHP; code 0 clean, 1 atmosphere
# reset, 2 conserved variables altered)
function _c2p(eng::DGStarFDEngine, D̂::Float64, Ŝx::Float64, τ̂in::Float64, elam::Float64, p0::Float64)
    Γ = eng.Γ; κ = eng.κ
    if !(D̂ > eng.ρ_cut) || !isfinite(D̂) || !isfinite(Ŝx) || !isfinite(τ̂in)
        return eng.ρ_atm, eng.ε_atm, 0.0, eng.p_atm, 1
    end
    code = 0
    Ŝ = abs(Ŝx)/sqrt(elam)
    τ̂ = τ̂in < 0.0 ? 1e-300 : τ̂in; τ̂ != τ̂in && (code = 2)
    Smax = sqrt(max(τ̂*(τ̂ + 2D̂), 0.0))*(1 - 1e-10)
    if Ŝ > Smax; Ŝ = Smax; code = 2; end
    p = max(p0, 1e-300); ok = false
    for _ in 1:30
        f, _, _, _, _ = _c2p_resid(Γ, D̂, Ŝ, τ̂, p); dp = max(1e-8*p, 1e-300)
        f2, _, _, _, _ = _c2p_resid(Γ, D̂, Ŝ, τ̂, p + dp); df = (f2 - f)/dp
        df ≥ 0 && break
        pn = p - f/df; pn ≤ 0 && (pn = 0.5*p)
        if abs(pn - p) ≤ 1e-13*max(pn, 1e-300); p = pn; ok = true; break; end
        p = pn
    end
    if !ok
        plo = 0.0; phi = max((Γ-1)*τ̂, 1e-300)
        while _c2p_resid(Γ, D̂, Ŝ, τ̂, phi)[1] > 0 && phi < 1e10; phi *= 2; end
        for _ in 1:200
            pm = 0.5*(plo + phi); fm = _c2p_resid(Γ, D̂, Ŝ, τ̂, pm)[1]
            if fm > 0; plo = pm else phi = pm end
            (phi - plo) ≤ 1e-13*phi && break
        end
        p = 0.5*(plo + phi)
    end
    _, ρ, ε, v̂, _ = _c2p_resid(Γ, D̂, Ŝ, τ̂, p)
    if !(ρ > eng.ρ_atm)
        return eng.ρ_atm, eng.ε_atm, 0.0, eng.p_atm, 1
    end
    εmin = eng.entropy_floor ? _εpoly(κ, Γ, ρ) : 0.0; εmax = eng.εfac_max*_εpoly(κ, Γ, ρ)
    if ε < εmin; ε = εmin; code = 2 elseif ε > εmax; ε = εmax; code = 2 end
    p = _p_ig(Γ, ρ, ε)
    return ρ, ε, copysign(v̂, Ŝx), p, code
end

# primitives on the DG nodes of one element (and rewrite the conserved if they were fixed)
function _prims_dg_elem!(st::DGStarFDState, eng::DGStarFDEngine, k::Int)
    e = eng.elems[k]
    @inbounds for a in 1:e.N
        i = e.off + a - 1; sg = eng.sqrtγ[i]; el = eng.elam[i]
        ρ, ε, v̂, p, code = _c2p(eng, st.D[i]/sg, st.S[i]/sg, st.τ[i]/sg, el, st.p[i])
        wasatm = st.ρ[i] ≤ eng.ρ_atm
        st.ρ[i] = ρ; st.ε[i] = ε; st.v[i] = v̂; st.p[i] = p
        if code != 0
            D̂, Ŝ, τ̂ = _p2c(eng.Γ, ρ, ε, v̂, el)
            st.D[i] = sg*D̂; st.S[i] = sg*Ŝ; st.τ[i] = sg*τ̂
            (code == 2 || !wasatm) && (st.nfix += 1)
        end
    end
end
# primitives on the subcells of one element
function _prims_fd_elem!(st::DGStarFDState, eng::DGStarFDEngine, k::Int)
    o = eng.soff[k]
    @inbounds for c in 1:eng.sM[k]
        i = o + c - 1; sg = eng.ssqrtγ[i]; el = eng.selam[i]
        ρ, ε, v̂, p, code = _c2p(eng, st.sD[i]/sg, st.sS[i]/sg, st.sτ[i]/sg, el, st.sp[i])
        wasatm = st.sρ[i] ≤ eng.ρ_atm
        st.sρ[i] = ρ; st.sε[i] = ε; st.sv[i] = v̂; st.sp[i] = p
        if code != 0
            D̂, Ŝ, τ̂ = _p2c(eng.Γ, ρ, ε, v̂, el)
            st.sD[i] = sg*D̂; st.sS[i] = sg*Ŝ; st.sτ[i] = sg*τ̂
            (code == 2 || !wasatm) && (st.nfix += 1)
        end
    end
end

# ----------------------------------------------------------------------------------
# conservative transfer between the two representations of one element
# ----------------------------------------------------------------------------------
"project the DG solution of element k onto its subcells (conservative)"
function _project_elem!(st::DGStarFDState, eng::DGStarFDEngine, k::Int)
    e = eng.elems[k]; P = eng.ops[e.p].P; o = eng.soff[k]; M = eng.sM[k]
    @inbounds for c in 1:M
        dD = 0.0; dS = 0.0; dτ = 0.0
        for a in 1:e.N
            j = e.off + a - 1; w = P[c, a]
            dD += w*st.D[j]; dS += w*st.S[j]; dτ += w*st.τ[j]
        end
        i = o + c - 1; st.sD[i] = dD; st.sS[i] = dS; st.sτ[i] = dτ
    end
end
"reconstruct the DG solution of element k from its subcells (mean-preserving least squares)"
function _reconstruct_elem!(st::DGStarFDState, eng::DGStarFDEngine, k::Int)
    e = eng.elems[k]; R = eng.ops[e.p].R; o = eng.soff[k]; M = eng.sM[k]
    @inbounds for a in 1:e.N
        dD = 0.0; dS = 0.0; dτ = 0.0
        for c in 1:M
            i = o + c - 1; w = R[a, c]
            dD += w*st.sD[i]; dS += w*st.sS[i]; dτ += w*st.sτ[i]
        end
        j = e.off + a - 1; st.D[j] = dD; st.S[j] = dS; st.τ[j] = dτ
    end
end

# ----------------------------------------------------------------------------------
# setup
# ----------------------------------------------------------------------------------
"""
    setup_dgstarfd(eos, εc; grid=:I1, Γ=2.0, κ=eos.κ, atm_fac=1e-13, cut_fac=10.0,
                   εfac_max=100.0, flux=:hll, entropy_floor=true, symmetrize=true,
                   tci_α=4.0, tci_nmodes=1, min_clear=5, always_fd=false,
                   p_center=7, center_fac=3, cfl=0.25, xmax_fac=3.0, h_tov=2e-4)
        -> (engine, state)

Build the DG–FD hybrid radial star on the TOV star of central energy density `εc`. `grid` is a
preset of [`hp_grid`](@ref) or an explicit region list, so the hybrid runs on exactly the grids
of the pure-DG study (VALIDATION.md §7.9). `always_fd=true` puts every element on the subgrid
(SpECTRE's `AlwaysUseSubcells`: a pure finite-volume run, the reference for how much the DG
elements buy). `min_clear` is how many consecutive steps an FD element's reconstructed DG
solution must pass the indicator before it returns to DG.
"""
function setup_dgstarfd(eos::BarotropicEOS, εc::Float64; grid=:I1, Γ::Float64=2.0,
        κ::Float64=(eos isa ShumPolytrope ? eos.κ : 100.0), atm_fac::Float64=1e-13,
        cut_fac::Float64=10.0, εfac_max::Float64=100.0, flux::Symbol=:hll,
        entropy_floor::Bool=true, symmetrize::Bool=true, tci_α::Float64=4.0,
        tci_nmodes::Int=1, min_clear::Int=5, always_fd::Bool=false,
        p_center::Int=7, center_fac::Int=3, cfl::Float64=0.25, xmax_fac::Float64=3.0,
        h_tov::Float64=2e-4)
    star = solve_tov(eos, εc; h=h_tov); R, M = star.R, star.M
    regions = grid isa Symbol ? hp_grid(R; preset=grid, xmax_fac=xmax_fac) : grid
    bases = Dict{Int,LGLBasis}()
    elems, Ntot = _build_elems(regions, bases; p_center=p_center, center_fac=center_fac)
    K = length(elems); kc = findfirst(e -> e.region == :center, elems)
    ops = Dict{Int,SubcellOps}(p => subcell_ops(b) for (p, b) in bases)
    rt = star.r; mt = star.m; νt = star.ν; pt = star.p; et = star.ε
    m_of(r) = r < R ? _lin(rt, mt, r) : M
    ν_of(r) = r < R ? _lin(rt, νt, r) : log(1 - 2M/r)
    p_of(r) = r < R ? max(_lin(rt, pt, r), 0.0) : 0.0
    e_of(r) = r < R ? max(_lin(rt, et, r), 0.0) : 0.0
    # metric coefficients at an arbitrary x (regular at x = 0)
    function geom(xx)
        rr = abs(xx)
        if rr < 1e-12
            return exp(0.5*ν_of(0.0)), 1.0, 0.0, 0.0, 0.0     # α, e^λ, √γ, Φ'_x, λ'_x
        end
        mm = m_of(rr); ff = max(1 - 2mm/rr, 1e-12)
        αv = exp(0.5*ν_of(rr)); elv = 1/ff; sgv = sqrt(elv)*rr^2
        Φpv = (mm + 4π*rr^3*p_of(rr))/(rr*(rr - 2mm))
        mpv = 4π*rr^2*e_of(rr); λpv = (2mpv/rr - 2mm/rr^2)/ff
        return αv, elv, sgv, sign(xx)*Φpv, sign(xx)*λpv
    end
    # --- DG node geometry
    x = zeros(Ntot); r = zeros(Ntot); α = zeros(Ntot); elam = zeros(Ntot)
    sqrtγ = zeros(Ntot); Φx = zeros(Ntot); λx = zeros(Ntot); dxmin = Inf
    for e in elems
        b = bases[e.p]
        for a in 1:e.N
            i = e.off + a - 1
            xx = 0.5*(e.xlo + e.xhi) + e.J*b.ξ[a]
            x[i] = xx; r[i] = abs(xx)
            α[i], elam[i], sqrtγ[i], Φx[i], λx[i] = geom(xx)
            a < e.N && (dxmin = min(dxmin, e.J*(b.ξ[a+1] - b.ξ[a])))
        end
    end
    # --- subcell geometry (centres and faces)
    sM = [ops[e.p].M for e in elems]
    soff = zeros(Int, K); sfoff = zeros(Int, K); off = 1; foff = 1
    for k in 1:K; soff[k] = off; off += sM[k]; sfoff[k] = foff; foff += sM[k] + 1; end
    Mtot = off - 1; Ftot = foff - 1
    shx = [(elems[k].xhi - elems[k].xlo)/sM[k] for k in 1:K]
    sx = zeros(Mtot); sr = zeros(Mtot); sα = zeros(Mtot); selam = zeros(Mtot)
    ssqrtγ = zeros(Mtot); sΦx = zeros(Mtot); sλx = zeros(Mtot)
    sfx = zeros(Ftot); sfα = zeros(Ftot); sfelam = zeros(Ftot); sfsqrtγ = zeros(Ftot)
    xq, wq = DGSubcell.gauss_legendre(5)
    for k in 1:K
        e = elems[k]; h = shx[k]; dxmin = min(dxmin, h)
        for c in 1:sM[k]
            i = soff[k] + c - 1
            xx = e.xlo + (c - 0.5)*h
            sx[i] = xx; sr[i] = abs(xx)
            sα[i], selam[i], _, sΦx[i], sλx[i] = geom(xx)
            # √γ = e^{λ/2}x² must be the CELL AVERAGE, not the value at the centre: the subcell
            # grid has an odd number of cells, so the central element's middle cell is centred on
            # x = 0 where the point value vanishes and the densitized state would carry no
            # information at all (the recovery then returns atmosphere and the star empties from
            # the centre — observed before this was fixed). The average is also the consistent
            # finite-volume reading of q̄ = (1/h)∫√γ ρW dx, and is O(h²) from the point value
            # everywhere else.
            sgav = 0.0
            for q in eachindex(xq)
                _, _, sgq, _, _ = geom(xx + 0.5*h*xq[q]); sgav += wq[q]*sgq
            end
            ssqrtγ[i] = 0.5*sgav
        end
        for c in 1:sM[k]+1
            i = sfoff[k] + c - 1
            xx = e.xlo + (c - 1)*h
            sfx[i] = xx
            sfα[i], sfelam[i], sfsqrtγ[i], _, _ = geom(xx)
        end
    end
    # --- initial data (the projected TOV star) on both grids
    ρc = rho_from_p(eos, pressure(eos, εc))
    ρ_atm = atm_fac*ρc; ε_atm = _εpoly(κ, Γ, ρ_atm); p_atm = _p_ig(Γ, ρ_atm, ε_atm)
    st = DGStarFDState(Ntot, Mtot)
    ρ_of(xx) = begin
        rr = abs(xx)
        ρv = rr < R*(1 - 1e-10) ? max(rho_from_p(eos, p_of(rr)), ρ_atm) : ρ_atm
        ρv ≤ cut_fac*ρ_atm ? ρ_atm : ρv
    end
    for i in 1:Ntot
        ρ = ρ_of(x[i]); ε = _εpoly(κ, Γ, ρ)
        st.ρ[i] = ρ; st.ε[i] = ε; st.p[i] = _p_ig(Γ, ρ, ε); st.v[i] = 0.0
        D̂, Ŝ, τ̂ = _p2c(Γ, ρ, ε, 0.0, elam[i]); sg = sqrtγ[i]
        st.D[i] = sg*D̂; st.S[i] = sg*Ŝ; st.τ[i] = sg*τ̂
    end
    for i in 1:Mtot
        ρ = ρ_of(sx[i]); ε = _εpoly(κ, Γ, ρ)
        st.sρ[i] = ρ; st.sε[i] = ε; st.sp[i] = _p_ig(Γ, ρ, ε); st.sv[i] = 0.0
        D̂, Ŝ, τ̂ = _p2c(Γ, ρ, ε, 0.0, selam[i]); sg = ssqrtγ[i]
        st.sD[i] = sg*D̂; st.sS[i] = sg*Ŝ; st.sτ[i] = sg*τ̂
    end
    active = fill(always_fd ? 1 : 0, K); clear = zeros(Int, K)
    eng = DGStarFDEngine(elems, bases, ops, K, Ntot, Mtot, kc, soff, sM, shx, sfoff,
                         Γ, κ, eos, R, M, ρ_atm, cut_fac*ρ_atm, ε_atm, p_atm, εfac_max,
                         flux, entropy_floor, symmetrize, tci_α, tci_nmodes, min_clear, always_fd,
                         cfl, dxmin, x, r, α, elam, sqrtγ, Φx, λx,
                         sx, sr, sα, selam, ssqrtγ, sΦx, sλx, sfx, sfα, sfelam, sfsqrtγ,
                         copy(st.D), copy(st.sD), active, clear)
    always_fd || _tci!(st, eng; initial=true)
    # the initial indicator projects the troubled elements onto their subgrids, so refresh the
    # equilibrium reference afterwards: err[D̃] must start from zero on whichever grid is active
    copyto!(eng.Deq, st.D); copyto!(eng.sDeq, st.sD)
    return eng, st
end

# ----------------------------------------------------------------------------------
# troubled-cell indicator
# ----------------------------------------------------------------------------------
# Is the DG representation of element k admissible and smooth? The order of the checks follows
# SpECTRE's `TciOptions`: an element lying ENTIRELY in the atmosphere is uniform and therefore
# perfectly smooth — putting it on the subgrid would only waste subcells — while an element that
# straddles the density cutoff is the surface, which is exactly what the subcells are for. The
# spectral-decay test is applied to the undensitized D̂ = ρW, i.e. with the √γ = e^{λ/2}r² factor
# of the geometry divided out, so that it measures the fluid and not the coordinate volume.
function _dg_ok(st::DGStarFDState, eng::DGStarFDEngine, k::Int)
    e = eng.elems[k]; b = eng.bases[e.p]; N = e.N
    D̂ = Vector{Float64}(undef, N); dmax = 0.0; straddles = false
    @inbounds for a in 1:N
        i = e.off + a - 1; sg = eng.sqrtγ[i]
        d = st.D[i]/sg; τ̂ = st.τ[i]/sg; Ŝ2 = (st.S[i]/sg)^2/eng.elam[i]
        (isfinite(d) && isfinite(τ̂) && d > 0.0 && τ̂ > 0.0 && Ŝ2 < τ̂*(τ̂ + 2d)) || return false
        D̂[a] = d; dmax = max(dmax, d); d ≤ eng.ρ_cut && (straddles = true)
    end
    dmax ≤ eng.ρ_cut && return true      # entirely atmosphere: uniform, no subcells needed
    straddles && return false             # spans the stellar surface: subcells
    return persson_ratio(D̂, b; nmodes=eng.tci_nmodes) ≤
           persson_threshold(b; α=eng.tci_α, nmodes=eng.tci_nmodes)
end

"""
    _tci!(st, eng; initial=false)

Choose the grid of every element for the coming step. A DG element whose solution is
inadmissible or spectrally rough drops to the subgrid (its DG state is projected); an FD element
returns to DG only after its reconstructed DG solution has been acceptable `min_clear`
consecutive times.
"""
function _tci!(st::DGStarFDState, eng::DGStarFDEngine; initial::Bool=false)
    eng.always_fd && return (0, 0)
    ndrop = 0; nback = 0
    for k in 1:eng.K
        if eng.active[k] == 0
            if !_dg_ok(st, eng, k)
                _project_elem!(st, eng, k); _prims_fd_elem!(st, eng, k)
                eng.active[k] = 1; eng.clear[k] = 0; ndrop += 1
            end
        else
            # trial reconstruction: keep the subcell state, test the DG candidate
            Dsave = [st.D[eng.elems[k].off + a - 1] for a in 1:eng.elems[k].N]
            Ssave = [st.S[eng.elems[k].off + a - 1] for a in 1:eng.elems[k].N]
            τsave = [st.τ[eng.elems[k].off + a - 1] for a in 1:eng.elems[k].N]
            _reconstruct_elem!(st, eng, k)
            ok = _dg_ok(st, eng, k)
            eng.clear[k] = ok ? eng.clear[k] + 1 : 0
            if ok && eng.clear[k] ≥ eng.min_clear && !initial
                _prims_dg_elem!(st, eng, k); eng.active[k] = 0; nback += 1
            else
                e = eng.elems[k]
                @inbounds for a in 1:e.N
                    j = e.off + a - 1; st.D[j] = Dsave[a]; st.S[j] = Ssave[a]; st.τ[j] = τsave[a]
                end
            end
        end
    end
    return ndrop, nback
end

# ----------------------------------------------------------------------------------
# hybrid right-hand side
# ----------------------------------------------------------------------------------
# face state (ρ, ε, v̂) that element k presents at its left (side=-1) or right (side=+1) edge
@inline function _edge_state(st::DGStarFDState, eng::DGStarFDEngine, k::Int, side::Int,
                             ghost::NTuple{4,Float64})
    e = eng.elems[k]
    if eng.active[k] == 0
        i = side < 0 ? e.off : e.off + e.N - 1
        return st.ρ[i], st.ε[i], st.v[i]
    end
    o = eng.soff[k]; M = eng.sM[k]; h = eng.shx[k]
    gρ, gε, gv, gh = ghost                      # neighbour cell values and its width
    fac = h/(0.5*(h + gh))                      # express the one-sided difference per own width
    if side < 0
        i0 = o; i1 = o + 1
        # the ghost value expressed at the element's own uniform spacing
        sρ = mc_slope(st.sρ[i0] - (st.sρ[i0] - gρ)*fac, st.sρ[i0], st.sρ[i1])
        sε = mc_slope(st.sε[i0] - (st.sε[i0] - gε)*fac, st.sε[i0], st.sε[i1])
        sv = mc_slope(st.sv[i0] - (st.sv[i0] - gv)*fac, st.sv[i0], st.sv[i1])
        return st.sρ[i0] - 0.5*sρ, st.sε[i0] - 0.5*sε, st.sv[i0] - 0.5*sv
    else
        i0 = o + M - 1; i1 = o + M - 2
        sρ = mc_slope(st.sρ[i1], st.sρ[i0], st.sρ[i0] + (gρ - st.sρ[i0])*fac)
        sε = mc_slope(st.sε[i1], st.sε[i0], st.sε[i0] + (gε - st.sε[i0])*fac)
        sv = mc_slope(st.sv[i1], st.sv[i0], st.sv[i0] + (gv - st.sv[i0])*fac)
        return st.sρ[i0] + 0.5*sρ, st.sε[i0] + 0.5*sε, st.sv[i0] + 0.5*sv
    end
end

# the ghost values element k sees beyond its left/right edge (neighbour's outermost subcell)
@inline function _ghost(st::DGStarFDState, eng::DGStarFDEngine, k::Int, side::Int)
    kn = side < 0 ? k - 1 : k + 1
    if kn < 1 || kn > eng.K
        return (eng.ρ_atm, eng.ε_atm, 0.0, eng.shx[k])       # outflow into the atmosphere
    end
    o = eng.soff[kn]; M = eng.sM[kn]
    i = side < 0 ? o + M - 1 : o
    return (st.sρ[i], st.sε[i], st.sv[i], eng.shx[kn])
end

function _rhs!(rD, rS, rτ, srD, srS, srτ, st::DGStarFDState, eng::DGStarFDEngine)
    Γ = eng.Γ; K = eng.K
    fill!(rD, 0.0); fill!(rS, 0.0); fill!(rτ, 0.0)
    fill!(srD, 0.0); fill!(srS, 0.0); fill!(srτ, 0.0)
    # --- primitives on the active grid, and subcell primitives wherever a ghost is needed
    for k in 1:K
        eng.active[k] == 0 ? _prims_dg_elem!(st, eng, k) : _prims_fd_elem!(st, eng, k)
    end
    for k in 1:K
        eng.active[k] == 0 || continue
        need = (k > 1 && eng.active[k-1] == 1) || (k < K && eng.active[k+1] == 1)
        need || continue
        _project_elem!(st, eng, k); _prims_fd_elem!(st, eng, k)
    end
    # --- common flux at every element face kf = 0..K
    F̂ = Vector{NTuple{3,Float64}}(undef, K + 1)
    @inbounds for kf in 0:K
        if kf == 0
            kR = 1; xf = eng.elems[1].xlo
            αf, elf, sgf = eng.sfα[eng.sfoff[1]], eng.sfelam[eng.sfoff[1]], eng.sfsqrtγ[eng.sfoff[1]]
            ρR, εR, vR = _edge_state(st, eng, kR, -1, _ghost(st, eng, kR, -1))
            ρL, εL, vL = eng.ρ_atm, eng.ε_atm, 0.0
        elseif kf == K
            kL = K; i = eng.sfoff[K] + eng.sM[K]
            αf, elf, sgf = eng.sfα[i], eng.sfelam[i], eng.sfsqrtγ[i]
            ρL, εL, vL = _edge_state(st, eng, kL, +1, _ghost(st, eng, kL, +1))
            ρR, εR, vR = eng.ρ_atm, eng.ε_atm, 0.0
        else
            kL = kf; kR = kf + 1; i = eng.sfoff[kR]
            αf, elf, sgf = eng.sfα[i], eng.sfelam[i], eng.sfsqrtγ[i]
            ρL, εL, vL = _edge_state(st, eng, kL, +1, _ghost(st, eng, kL, +1))
            ρR, εR, vR = _edge_state(st, eng, kR, -1, _ghost(st, eng, kR, -1))
        end
        qL, fL, smL, spL = _qfg(Γ, αf, elf, sgf, ρL, εL, vL)
        qR, fR, smR, spR = _qfg(Γ, αf, elf, sgf, ρR, εR, vR)
        F̂[kf+1] = _numflux(eng.flux, qL, fL, smL, spL, qR, fR, smR, spR)
    end
    # --- DG elements: strong-form volume term + boundary correction with the common flux
    @inbounds for k in 1:K
        eng.active[k] == 0 || continue
        e = eng.elems[k]; b = eng.bases[e.p]; Dm = b.D; N = e.N; o = e.off - 1
        fD = Vector{Float64}(undef, N); fS = similar(fD); fτ = similar(fD)
        for a in 1:N
            i = o + a
            _, f, _, _ = _qfg(Γ, eng.α[i], eng.elam[i], eng.sqrtγ[i], st.ρ[i], st.ε[i], st.v[i])
            fD[a] = f[1]; fS[a] = f[2]; fτ[a] = f[3]
        end
        for a in 1:N
            sD = 0.0; sS = 0.0; sτ = 0.0
            for c in 1:N
                d = Dm[a, c]; sD += d*fD[c]; sS += d*fS[c]; sτ += d*fτ[c]
            end
            rD[o+a] = -sD/e.J; rS[o+a] = -sS/e.J; rτ[o+a] = -sτ/e.J
        end
        # boundary corrections (strong form: F̂ − F_interior at the two endpoints)
        FL = F̂[k]; FR = F̂[k+1]
        w1 = b.w[1]; wN = b.w[N]
        rD[o+1] += (FL[1] - fD[1])/(e.J*w1); rS[o+1] += (FL[2] - fS[1])/(e.J*w1); rτ[o+1] += (FL[3] - fτ[1])/(e.J*w1)
        rD[o+N] -= (FR[1] - fD[N])/(e.J*wN); rS[o+N] -= (FR[2] - fS[N])/(e.J*wN); rτ[o+N] -= (FR[3] - fτ[N])/(e.J*wN)
        for a in 1:N
            i = o + a
            sS, sτv = _src(Γ, eng.α[i], eng.elam[i], eng.sqrtγ[i], eng.x[i], eng.Φx[i], eng.λx[i],
                           st.ρ[i], st.ε[i], st.p[i], st.v[i])
            rS[i] += sS; rτ[i] += sτv
        end
    end
    # --- FD elements: finite-volume divergence on the subcells
    @inbounds for k in 1:K
        eng.active[k] == 1 || continue
        o = eng.soff[k]; M = eng.sM[k]; h = eng.shx[k]; fo = eng.sfoff[k]
        gL = _ghost(st, eng, k, -1); gR = _ghost(st, eng, k, +1)
        facL = h/(0.5*(h + gL[4])); facR = h/(0.5*(h + gR[4]))
        # cell-centred primitives with one ghost on each side, in local indices 0..M+1
        getρ(c) = c == 0 ? st.sρ[o] - (st.sρ[o] - gL[1])*facL :
                          (c == M+1 ? st.sρ[o+M-1] + (gR[1] - st.sρ[o+M-1])*facR : st.sρ[o+c-1])
        getε(c) = c == 0 ? st.sε[o] - (st.sε[o] - gL[2])*facL :
                          (c == M+1 ? st.sε[o+M-1] + (gR[2] - st.sε[o+M-1])*facR : st.sε[o+c-1])
        getv(c) = c == 0 ? st.sv[o] - (st.sv[o] - gL[3])*facL :
                          (c == M+1 ? st.sv[o+M-1] + (gR[3] - st.sv[o+M-1])*facR : st.sv[o+c-1])
        Fl = Vector{NTuple{3,Float64}}(undef, M + 1)
        Fl[1] = F̂[k]; Fl[M+1] = F̂[k+1]
        for j in 2:M                      # interior subcell face between cells j−1 and j
            i = fo + j - 1
            αf = eng.sfα[i]; elf = eng.sfelam[i]; sgf = eng.sfsqrtγ[i]
            ρm = getρ(j-2); ρ0 = getρ(j-1); ρp = getρ(j); ρq = getρ(j+1)
            εm = getε(j-2); ε0 = getε(j-1); εp = getε(j); εq = getε(j+1)
            vm = getv(j-2); v0 = getv(j-1); vp = getv(j); vq = getv(j+1)
            ρL = ρ0 + 0.5*mc_slope(ρm, ρ0, ρp); ρR = ρp - 0.5*mc_slope(ρ0, ρp, ρq)
            εL = ε0 + 0.5*mc_slope(εm, ε0, εp); εR = εp - 0.5*mc_slope(ε0, εp, εq)
            vL = v0 + 0.5*mc_slope(vm, v0, vp); vR = vp - 0.5*mc_slope(v0, vp, vq)
            ρL = max(ρL, eng.ρ_atm); ρR = max(ρR, eng.ρ_atm)
            εL = max(εL, 0.0); εR = max(εR, 0.0)
            qL, fL, smL, spL = _qfg(Γ, αf, elf, sgf, ρL, εL, vL)
            qR, fR, smR, spR = _qfg(Γ, αf, elf, sgf, ρR, εR, vR)
            Fl[j] = _numflux(eng.flux, qL, fL, smL, spL, qR, fR, smR, spR)
        end
        for c in 1:M
            i = o + c - 1
            srD[i] = -(Fl[c+1][1] - Fl[c][1])/h
            srS[i] = -(Fl[c+1][2] - Fl[c][2])/h
            srτ[i] = -(Fl[c+1][3] - Fl[c][3])/h
            sS, sτv = _src(Γ, eng.sα[i], eng.selam[i], eng.ssqrtγ[i], eng.sx[i], eng.sΦx[i],
                           eng.sλx[i], st.sρ[i], st.sε[i], st.sp[i], st.sv[i])
            srS[i] += sS; srτ[i] += sτv
        end
    end
    return nothing
end

# ----------------------------------------------------------------------------------
# parity projection of the doubled domain (D̃, τ̃ even in x; S̃ odd)
# ----------------------------------------------------------------------------------
function _symmetrize!(st::DGStarFDState, eng::DGStarFDEngine)
    K = eng.K
    @inbounds for k in 1:(K÷2)
        km = K + 1 - k
        if eng.active[k] == 0 && eng.active[km] == 0
            e = eng.elems[k]; m = eng.elems[km]; N = e.N
            for a in 1:N
                i = e.off + a - 1; j = m.off + N - a
                d = 0.5*(st.D[i] + st.D[j]); st.D[i] = d; st.D[j] = d
                t = 0.5*(st.τ[i] + st.τ[j]); st.τ[i] = t; st.τ[j] = t
                s = 0.5*(st.S[i] - st.S[j]); st.S[i] = s; st.S[j] = -s
            end
        elseif eng.active[k] == 1 && eng.active[km] == 1
            o = eng.soff[k]; om = eng.soff[km]; M = eng.sM[k]
            for c in 1:M
                i = o + c - 1; j = om + M - c
                d = 0.5*(st.sD[i] + st.sD[j]); st.sD[i] = d; st.sD[j] = d
                t = 0.5*(st.sτ[i] + st.sτ[j]); st.sτ[i] = t; st.sτ[j] = t
                s = 0.5*(st.sS[i] - st.sS[j]); st.sS[i] = s; st.sS[j] = -s
            end
        end
    end
    kc = eng.kc
    if eng.active[kc] == 0
        e = eng.elems[kc]; N = e.N
        @inbounds for a in 1:(N÷2)
            i = e.off + a - 1; j = e.off + N - a
            d = 0.5*(st.D[i] + st.D[j]); st.D[i] = d; st.D[j] = d
            t = 0.5*(st.τ[i] + st.τ[j]); st.τ[i] = t; st.τ[j] = t
            s = 0.5*(st.S[i] - st.S[j]); st.S[i] = s; st.S[j] = -s
        end
    else
        o = eng.soff[kc]; M = eng.sM[kc]
        @inbounds for c in 1:(M÷2)
            i = o + c - 1; j = o + M - c
            d = 0.5*(st.sD[i] + st.sD[j]); st.sD[i] = d; st.sD[j] = d
            t = 0.5*(st.sτ[i] + st.sτ[j]); st.sτ[i] = t; st.sτ[j] = t
            s = 0.5*(st.sS[i] - st.sS[j]); st.sS[i] = s; st.sS[j] = -s
        end
        M % 2 == 1 && (st.sS[o + (M+1)÷2 - 1] = 0.0)      # the cell centred on x = 0
    end
end

# ----------------------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------------------
"""central rest-mass density: the subcell centred on x=0 when the central element is on the
subgrid, otherwise the DG interpolant at x=0."""
function dgstarfd_central_density(st::DGStarFDState, eng::DGStarFDEngine)
    kc = eng.kc
    if eng.active[kc] == 1
        M = eng.sM[kc]
        return st.sρ[eng.soff[kc] + (M + 1)÷2 - 1]
    end
    e = eng.elems[kc]; b = eng.bases[e.p]; N = e.N; s = 0.0
    @inbounds for a in 1:N
        ℓ = 1.0
        for c in 1:N; c == a && continue; ℓ *= (0.0 - b.ξ[c])/(b.ξ[a] - b.ξ[c]); end
        s += ℓ*st.ρ[e.off + a - 1]
    end
    s
end

"""normalized L2 error of D̃ against the projected TOV star, on whichever grid each element uses."""
function dgstarfd_errD(st::DGStarFDState, eng::DGStarFDEngine)
    num = 0.0; den = 0.0
    @inbounds for k in 1:eng.K
        if eng.active[k] == 0
            e = eng.elems[k]
            for a in 1:e.N
                i = e.off + a - 1; num += (st.D[i] - eng.Deq[i])^2; den += eng.Deq[i]^2
            end
        else
            o = eng.soff[k]
            for c in 1:eng.sM[k]
                i = o + c - 1; num += (st.sD[i] - eng.sDeq[i])^2; den += eng.sDeq[i]^2
            end
        end
    end
    sqrt(num/max(den, 1e-300))
end

"""baryon mass M_b = 2π∫D̃ dx over the symmetric domain, each element on its active grid."""
function dgstarfd_baryon_mass(st::DGStarFDState, eng::DGStarFDEngine)
    s = 0.0
    @inbounds for k in 1:eng.K
        if eng.active[k] == 0
            e = eng.elems[k]; w = eng.bases[e.p].w
            for a in 1:e.N; s += w[a]*e.J*st.D[e.off + a - 1]; end
        else
            o = eng.soff[k]; h = eng.shx[k]
            for c in 1:eng.sM[k]; s += h*st.sD[o + c - 1]; end
        end
    end
    2π*s
end

"fraction of elements currently evolved on the subgrid"
dgstarfd_fd_fraction(eng::DGStarFDEngine) = count(==(1), eng.active)/eng.K
"radial extent (x range) of the elements currently on the subgrid"
function dgstarfd_active_map(eng::DGStarFDEngine)
    fd = [k for k in 1:eng.K if eng.active[k] == 1]
    isempty(fd) && return (n=0, xlo=NaN, xhi=NaN, r_over_R=Float64[])
    return (n=length(fd), xlo=minimum(eng.elems[k].xlo for k in fd),
            xhi=maximum(eng.elems[k].xhi for k in fd),
            r_over_R=sort(unique(round.([0.5*(abs(eng.elems[k].xlo) + abs(eng.elems[k].xhi))/eng.R for k in fd], digits=3))))
end

"""Hanning-windowed periodogram of ρ_c(t); peaks sorted by frequency (kHz), times in M⊙."""
function dgstarfd_spectrum(ts::Vector{Float64}, ρc::Vector{Float64}; fmin_kHz=1.0, fmax_kHz=15.0,
                           npts::Int=6000, npeaks::Int=6, tmax::Float64=Inf, prominence::Float64=3.0)
    m = ts .≤ tmax; t = ts[m]; q = copy(ρc[m]); n = length(q); μ = sum(q)/n
    @inbounds for i in 1:n; q[i] = (q[i] - μ)*(0.5 - 0.5*cos(2π*(i-1)/(n-1))); end
    conv = kHz_to_km*Msun_to_km
    νgrid = collect(range(fmin_kHz*conv, fmax_kHz*conv; length=npts))
    P = periodogram(t, q, νgrid)
    peaks = Tuple{Float64,Float64}[]
    win = max(5, npts ÷ 60)
    for i in 2:npts-1
        (P[i] > P[i-1] && P[i] ≥ P[i+1]) || continue
        lo = max(1, i-win); hi = min(npts, i+win); med = sort(P[lo:hi])[(hi-lo)÷2+1]
        P[i] > prominence*med && push!(peaks, (P[i], νgrid[i]/conv))
    end
    sort!(peaks, by=x->-x[1]); peaks = peaks[1:min(npeaks, length(peaks))]
    sort!(peaks, by=x->x[2])
    return [pk[2] for pk in peaks], νgrid./conv, P
end

# ----------------------------------------------------------------------------------
# evolution
# ----------------------------------------------------------------------------------
function _maxspeed(st::DGStarFDState, eng::DGStarFDEngine)
    a = 0.0
    @inbounds for k in 1:eng.K
        if eng.active[k] == 0
            e = eng.elems[k]
            for q in 1:e.N
                i = e.off + q - 1
                cs2 = _cs2(eng.Γ, st.ρ[i], st.ε[i]); λm, λp = _λpm(cs2, st.v[i])
                c = eng.α[i]/sqrt(eng.elam[i]); a = max(a, c*abs(λm), c*abs(λp))
            end
        else
            o = eng.soff[k]
            for c in 1:eng.sM[k]
                i = o + c - 1
                cs2 = _cs2(eng.Γ, st.sρ[i], st.sε[i]); λm, λp = _λpm(cs2, st.sv[i])
                cc = eng.sα[i]/sqrt(eng.selam[i]); a = max(a, cc*abs(λm), cc*abs(λp))
            end
        end
    end
    a
end

"""
    evolve_dgstarfd!(st, eng; tmax, sample_dt=1.0, dt=0.0, verbose=false) -> rec

SSP-RK3 evolution of the hybrid. The troubled-cell indicator runs once per step (a priori); the
active grid can change between steps but is fixed within one step's three stages. Returns
vectors `t, ρc, errD, Mb, fd_fraction, vatm` plus `nfix`, `dt`, `nsteps`, `ndrop`, `nback`.
"""
function evolve_dgstarfd!(st::DGStarFDState, eng::DGStarFDEngine; tmax::Float64,
                          sample_dt::Float64=1.0, dt::Float64=0.0, verbose::Bool=false)
    n = eng.Ntot; m = eng.Mtot
    rD = zeros(n); rS = zeros(n); rτ = zeros(n)
    srD = zeros(m); srS = zeros(m); srτ = zeros(m)
    D0 = zeros(n); S0 = zeros(n); τ0 = zeros(n)
    sD0 = zeros(m); sS0 = zeros(m); sτ0 = zeros(m)
    for k in 1:eng.K
        eng.active[k] == 0 ? _prims_dg_elem!(st, eng, k) : _prims_fd_elem!(st, eng, k)
    end
    dt = dt > 0 ? dt : eng.cfl*eng.dxmin/max(_maxspeed(st, eng), 1e-3)
    ts = Float64[]; ρc = Float64[]; errD = Float64[]; Mb = Float64[]; fdf = Float64[]; vatm = Float64[]
    ndrop = 0; nback = 0
    outside = eng.r .> eng.R; soutside = eng.sr .> eng.R
    record!(t) = begin
        push!(ts, t); push!(ρc, dgstarfd_central_density(st, eng)); push!(errD, dgstarfd_errD(st, eng))
        push!(Mb, dgstarfd_baryon_mass(st, eng)); push!(fdf, dgstarfd_fd_fraction(eng))
        vm = 0.0
        @inbounds for k in 1:eng.K
            if eng.active[k] == 0
                e = eng.elems[k]
                for a in 1:e.N; i = e.off + a - 1; outside[i] && (vm = max(vm, abs(st.v[i]))); end
            else
                o = eng.soff[k]
                for c in 1:eng.sM[k]; i = o + c - 1; soutside[i] && (vm = max(vm, abs(st.sv[i]))); end
            end
        end
        push!(vatm, vm)
    end
    record!(0.0)
    t = 0.0; last = 0.0; ns = 0; t0 = time()
    stage!(c0, c1, h) = begin
        _rhs!(rD, rS, rτ, srD, srS, srτ, st, eng)
        @inbounds for k in 1:eng.K
            if eng.active[k] == 0
                e = eng.elems[k]
                for a in 1:e.N
                    i = e.off + a - 1
                    st.D[i] = c0*D0[i] + c1*(st.D[i] + h*rD[i])
                    st.S[i] = c0*S0[i] + c1*(st.S[i] + h*rS[i])
                    st.τ[i] = c0*τ0[i] + c1*(st.τ[i] + h*rτ[i])
                end
            else
                o = eng.soff[k]
                for c in 1:eng.sM[k]
                    i = o + c - 1
                    st.sD[i] = c0*sD0[i] + c1*(st.sD[i] + h*srD[i])
                    st.sS[i] = c0*sS0[i] + c1*(st.sS[i] + h*srS[i])
                    st.sτ[i] = c0*sτ0[i] + c1*(st.sτ[i] + h*srτ[i])
                end
            end
        end
        eng.symmetrize && _symmetrize!(st, eng)
        for k in 1:eng.K
            eng.active[k] == 0 ? _prims_dg_elem!(st, eng, k) : _prims_fd_elem!(st, eng, k)
        end
    end
    while t < tmax - 1e-12
        d, b = _tci!(st, eng); ndrop += d; nback += b
        h = min(dt, tmax - t)
        copyto!(D0, st.D); copyto!(S0, st.S); copyto!(τ0, st.τ)
        copyto!(sD0, st.sD); copyto!(sS0, st.sS); copyto!(sτ0, st.sτ)
        stage!(0.0, 1.0, h); stage!(0.75, 0.25, h); stage!(1/3, 2/3, h)
        t += h; ns += 1
        if t - last ≥ sample_dt - 1e-12
            record!(t); last = t
            verbose && @info "t=$(round(t,digits=1)) errD=$(errD[end]) ρc/ρc0-1=$(ρc[end]/ρc[1]-1) FD=$(round(100*fdf[end]))% [$(round(time()-t0)) s]"
            (isfinite(ρc[end]) && ρc[end] < 100*ρc[1]) || break
        end
    end
    return (t=ts, ρc=ρc, errD=errD, Mb=Mb, fd_fraction=fdf, vatm=vatm,
            nfix=st.nfix, dt=dt, nsteps=ns, ndrop=ndrop, nback=nback, seconds=time()-t0)
end

"""seed a radial velocity v̂ = A sin(πr/R) inside the star (odd in x), on both representations."""
function seed_dgstarfd_radial!(st::DGStarFDState, eng::DGStarFDEngine; A::Float64=1e-3)
    @inbounds for i in 1:eng.Ntot
        rr = eng.r[i]; (rr < eng.R && st.ρ[i] > 10eng.ρ_atm) || continue
        v̂ = A*sin(π*rr/eng.R)*sign(eng.x[i]); st.v[i] = v̂
        D̂, Ŝ, τ̂ = _p2c(eng.Γ, st.ρ[i], st.ε[i], v̂, eng.elam[i]); sg = eng.sqrtγ[i]
        st.D[i] = sg*D̂; st.S[i] = sg*Ŝ; st.τ[i] = sg*τ̂
    end
    @inbounds for i in 1:eng.Mtot
        rr = eng.sr[i]; (rr < eng.R && st.sρ[i] > 10eng.ρ_atm) || continue
        v̂ = A*sin(π*rr/eng.R)*sign(eng.sx[i]); st.sv[i] = v̂
        D̂, Ŝ, τ̂ = _p2c(eng.Γ, st.sρ[i], st.sε[i], v̂, eng.selam[i]); sg = eng.ssqrtγ[i]
        st.sD[i] = sg*D̂; st.sS[i] = sg*Ŝ; st.sτ[i] = sg*τ̂
    end
end

end # module DGStarFD
