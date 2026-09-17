#=
    DGStarHP — hp-adapted nodal RKDG for the radial (Cowling) TOV star, after
    Hébert, Kidder & Teukolsky, PRD 98, 044041 (2018) [arXiv:1804.02003], Sec. VI.A.

    What is taken from the paper
      • a SYMMETRIC, STAGGERED domain x ∈ [−x_max, x_max], r = |x|, with the central
        element straddling the origin and odd polynomial degree there, so that no
        node sits on the r=0 coordinate singularity (their "staggered grid");
      • regions of different element SIZE and polynomial ORDER: high order (p=3)
        inside and outside the star, thin low-order (p=1) elements across the
        surface — the hp-adapted "I1/I2/I1R" grids of their Table I, rescaled to
        our star.  We additionally put the stellar surface EXACTLY on an element
        boundary (theirs sits inside a thin element);
      • Valencia GRHD on the frozen TOV metric in areal coordinates, with the
        COVARIANT momentum, the Γ-law ideal-gas EOS for the evolution (the
        polytrope is initial data only, so shock heating is allowed), and the
        atmosphere fixing of Galeazzi et al. 2013 (density cut, atmosphere reset,
        entropy bounds κρ^(Γ−1)/(Γ−1) ≤ ε ≤ 100× that);
      • HLL flux with Davis speed estimates (or LLF), SSP-RK3;
      • the minmod (ΛΠ¹) slope limiter applied ONLY to the low-order surface
        elements, followed by the physical-state check (D>ρ_atm, τ>0,
        S²<τ(τ+2D)) with slope reduction;
      • the exponential modal filter exp[−α(i/p)^s] on the high-order elements
        (available; not needed once the reflection parity is enforced, see below).

    What differs
      • WELL-BALANCED: the static DG residual of the projected TOV star is stored
        and subtracted (lake-at-rest), so the exact projection is a fixed point
        on every grid.  The paper's star instead "settles to its numerical
        equilibrium" with err[D̃] ~ 1e-4–1e-3; here that number measures the
        DYNAMICS only.
      • the limiter acts on the DENSITIZED conserved variables, so the minmod
        conserves ∫D̃ dx = baryon mass exactly (the paper reports a slow M_b
        growth because their limiter ignores the 4πr² volume element);
      • the filter acts on the deviation from equilibrium by default (filtering
        the equilibrium itself would erode its modal content every step);
      • the state is projected onto the physical reflection parity after every
        stage (symmetrize=true): the doubled domain has an unphysical antisymmetric
        (translation) mode that grows at the dynamical rate from round-off;
      • three limiters can be compared on the same grid:
          :minmod  — the paper's choice (slope limiter, surface elements only);
          :mean    — mean-based positivity scaling U ← Ū + θ(U−Ū), the
                     Zhang–Shu / MRS-type rescaling the paper found to let the
                     star expand past its surface, and the mechanism that eroded
                     DGCart3D (VALIDATION.md §7.8);
          :wb      — the same scaling about the EQUILIBRIUM-PRESERVING reference
                     of DGCart3D (deviation shared ∝ background, kinetic energy
                     carried), the fix of §7.8.

    Equations (areal radius r, x-coordinate on the symmetric domain, β=0, K_ij=0):
        √γ = e^{λ/2} x²,  α = e^{Φ},  Φ'_x = sign(x) Φ'(|x|),  λ'_x = sign(x) λ'(|x|)
        q = (D̃, S̃_x, τ̃) = √γ (ρW, ρhW² v_x, ρhW² − p − ρW),   v_x = e^{λ/2} v̂,
        f = α (D̃ v^x, S̃_x v^x + √γ p, τ̃ v^x + √γ p v^x),        v^x = e^{−λ/2} v̂,
        s_S = α√γ[(ρhW² v̂² + p) λ'_x/2 + 2p/x] − α√γ (ρhW² − p) Φ'_x,
        s_τ = −α√γ ρhW² v^x Φ'_x,
    which reduce to TOV (p' = −(ε+p)Φ') for the static star.  (The older DGStar
    module used an orthonormal-momentum shortcut whose static residual is not
    TOV; it is rescued only by its well-balanced subtraction.  Do not reuse it
    for dynamics.)

    Diagnostics reproduce the paper's Figs. 9, 12, 13: err[D̃] (their eqs. 35–36),
    ‖S̃‖, ρ_c(t), M_b(t), and the Hanning-windowed spectrum of ρ_c against the
    linear radial Cowling spectrum (RadialModes.radial_cowling_spectrum).

    RESULTS (repro/dgstarhp_hkt.jl, VALIDATION.md §7.9; F=2.686, H1=4.550, H2=6.341,
    H3=8.108 kHz):
      • the paper's I1 + minmod, no subtraction: err[D̃] 1.1e-3 at 1e4 M⊙ (theirs 7e-4),
        ρ_c −3e-4 (theirs −5e-4), M_b conserved to 1e-11, F and H3 in the spectrum of
        the settling transient. Their I1R improvement is NOT reproduced (5e-3, slowly
        growing); any seeded oscillation on linear elements + minmod blows a wind.
      • I2 (quadratic surface) + :wb limiter: no subtraction → settles to 1.0e-3 (I2)
        and 1.6e-5 (I2R) with a quiet atmosphere; with subtraction → fixed point to
        1e-12; seeded on I2R → F, H1, H2, H3 all within 0.7% of linear theory.
=#
module DGStarHP

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, energy_from_pressure
using ..TOV
using ..TOV: solve_tov
using ..FVCommon: rho_from_p
using ..DGCommon
using ..DGCommon: LGLBasis, build_lgl_basis
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram

export DGStarHPEngine, DGStarHPState, hp_grid, setup_dgstarhp, evolve_dgstarhp!,
       seed_dgstarhp_radial!, dgstarhp_central_density, dgstarhp_errD,
       dgstarhp_baryon_mass, dgstarhp_spectrum, dgstarhp_surface_state

# ----------------------------------------------------------------------------------
# grid description
# ----------------------------------------------------------------------------------
struct HPElem
    p::Int; N::Int; off::Int          # degree, nodes, offset of first node in flat arrays
    xlo::Float64; xhi::Float64; J::Float64; h::Float64
    region::Symbol                     # :center | :interior | :surface | :exterior
end

"""
    hp_grid(R; preset=:I1, xmax_fac=3.0) -> Vector{NamedTuple}

Region list `(xlo, xhi, K, p, region)` for x ≥ 0 (mirrored automatically). The presets
rescale Table I of Hébert–Kidder–Teukolsky from their r_NS = 8.125 to the star's areal
radius `R`, and put the surface on an element boundary:

* `:I1`  — interior 12½ elements p=3 (h≈0.72 for R=9.59), 10 linear elements of
           h=R/32 across [R−2h, R+8h], exterior 7 elements p=3 to 3R.
* `:I2`  — as I1 but 5 quadratic (p=2) surface elements of width 2h.
* `:I1R` — interior 50½ elements p=3, 20 linear elements of h=R/160 across
           [R−2h, R+18h], exterior 30 elements p=3.
* `:I2R` — as I1R with 10 quadratic surface elements of h=R/80 across [R−h, R+9h]
           (the refined form of the configuration recommended below).
* `:coarse` — a cheap test grid: 6½ interior p=3, 6 linear surface, 4 exterior p=3.

Custom grids: pass any vector of `(xlo, xhi, K, p, region)` with `xlo` of the first
region equal to 0 (the first region is split so that its central element straddles 0).
"""
function hp_grid(R::Float64; preset::Symbol=:I1, xmax_fac::Float64=3.0)
    xmax = xmax_fac*R
    if preset == :I1
        h = R/32; ri = R-2h; ro = R+8h
        return [(xlo=0.0, xhi=ri, K=12, p=3, region=:interior), (xlo=ri, xhi=ro, K=10, p=1, region=:surface),
                (xlo=ro, xhi=xmax, K=7, p=3, region=:exterior)]
    elseif preset == :I2
        h = R/16; ri = R-h; ro = R+4h
        return [(xlo=0.0, xhi=ri, K=12, p=3, region=:interior), (xlo=ri, xhi=ro, K=5, p=2, region=:surface),
                (xlo=ro, xhi=xmax, K=7, p=3, region=:exterior)]
    elseif preset == :I1R
        h = R/160; ri = R-2h; ro = R+18h
        return [(xlo=0.0, xhi=ri, K=50, p=3, region=:interior), (xlo=ri, xhi=ro, K=20, p=1, region=:surface),
                (xlo=ro, xhi=xmax, K=30, p=3, region=:exterior)]
    elseif preset == :I2R
        h = R/80; ri = R-h; ro = R+9h
        return [(xlo=0.0, xhi=ri, K=50, p=3, region=:interior), (xlo=ri, xhi=ro, K=10, p=2, region=:surface),
                (xlo=ro, xhi=xmax, K=30, p=3, region=:exterior)]
    elseif preset == :coarse
        h = R/16; ri = R-2h; ro = R+4h
        return [(xlo=0.0, xhi=ri, K=6, p=3, region=:interior), (xlo=ri, xhi=ro, K=6, p=1, region=:surface),
                (xlo=ro, xhi=xmax, K=4, p=3, region=:exterior)]
    else
        error("unknown preset $preset")
    end
end

# build the symmetric element list from the x≥0 region list; the first region contributes
# K full elements per side plus ONE central element of the same width straddling 0.
function _build_elems(regions, bases::Dict{Int,LGLBasis}; p_center::Int=5, center_fac::Int=3)
    r1 = regions[1]; r1.xlo == 0.0 || error("first region must start at x=0")
    isodd(p_center) || error("the central element needs odd degree so that no node sits at x=0")
    isodd(center_fac) || error("center_fac must be odd (the central element spans center_fac nominal widths)")
    h1 = (r1.xhi - r1.xlo)/(r1.K + 0.5)                    # nominal interior spacing
    m = (center_fac-1)÷2                                   # the central element absorbs m elements per side
    pos = Tuple{Float64,Float64,Int,Symbol}[]              # (xlo,xhi,p,region) for x>0 (excluding center)
    for j in (m+1):r1.K
        push!(pos, (h1/2 + (j-1)*h1, h1/2 + j*h1, r1.p, r1.region))
    end
    for reg in regions[2:end]
        h = (reg.xhi - reg.xlo)/reg.K
        for j in 1:reg.K
            push!(pos, (reg.xlo + (j-1)*h, reg.xlo + j*h, reg.p, reg.region))
        end
    end
    elems = HPElem[]; off = 1
    for (xlo,xhi,p,reg) in reverse(pos)                    # negative side, outermost first
        haskey(bases,p) || (bases[p] = build_lgl_basis(p))
        N = p+1; push!(elems, HPElem(p,N,off,-xhi,-xlo,(xhi-xlo)/2,xhi-xlo,reg)); off += N
    end
    haskey(bases,p_center) || (bases[p_center] = build_lgl_basis(p_center))
    hc = center_fac*h1
    push!(elems, HPElem(p_center, p_center+1, off, -hc/2, hc/2, hc/2, hc, :center)); off += p_center+1
    for (xlo,xhi,p,reg) in pos
        N = p+1; push!(elems, HPElem(p,N,off,xlo,xhi,(xhi-xlo)/2,xhi-xlo,reg)); off += N
    end
    return elems, off-1
end

# ----------------------------------------------------------------------------------
# engine / state
# ----------------------------------------------------------------------------------
struct DGStarHPEngine
    elems::Vector{HPElem}
    bases::Dict{Int,LGLBasis}
    K::Int; Ntot::Int; kc::Int                 # kc = index of the central element
    Γ::Float64; κ::Float64
    eos::BarotropicEOS
    R::Float64; M::Float64
    ρ_atm::Float64; ρ_cut::Float64; ε_atm::Float64; p_atm::Float64; εfac_max::Float64
    flux::Symbol; limiter::Symbol; p_lim::Int
    filter::Symbol; p_filt::Int; filt_α::Float64; filt_s::Int; filter_deviation::Bool
    wellbalanced::Bool; symmetrize::Bool; entropy_floor::Bool
    cfl::Float64; dxmin::Float64
    # node geometry (flat, length Ntot)
    x::Vector{Float64}; r::Vector{Float64}; α::Vector{Float64}; elam::Vector{Float64}
    sqrtγ::Vector{Float64}; Φx::Vector{Float64}; λx::Vector{Float64}
    # equilibrium (projected TOV) and its static DG residual
    Deq::Vector{Float64}; Seq::Vector{Float64}; τeq::Vector{Float64}
    Seq_D::Vector{Float64}; Seq_S::Vector{Float64}; Seq_τ::Vector{Float64}
end

mutable struct DGStarHPState
    D::Vector{Float64}; S::Vector{Float64}; τ::Vector{Float64}       # densitized conserved
    ρ::Vector{Float64}; ε::Vector{Float64}; p::Vector{Float64}; v::Vector{Float64}   # v = v̂ (orthonormal, along +x)
    nfix::Int                                                        # atmosphere/entropy fixes applied so far
end
DGStarHPState(n) = DGStarHPState((zeros(n) for _ in 1:7)..., 0)

# ----------------------------------------------------------------------------------
# thermodynamics (Γ-law) and conversions
# ----------------------------------------------------------------------------------
@inline _p_ig(Γ,ρ,ε) = (Γ-1)*ρ*ε
@inline _εpoly(κ,Γ,ρ) = κ*ρ^(Γ-1)/(Γ-1)
@inline function _cs2(Γ,ρ,ε)
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ
    return clamp(Γ*p/(ρ*h), 0.0, 1.0)
end
# primitives (ρ, ε, v̂) → undensitized (D̂, Ŝ_x covariant, τ̂); elam = e^λ
@inline function _p2c(Γ,ρ,ε,v̂,elam)
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ; W=1/sqrt(1-clamp(v̂^2,0.0,1-1e-14))
    D̂=ρ*W; Ŝ=ρ*h*W^2*sqrt(elam)*v̂; τ̂=ρ*h*W^2-p-D̂
    return D̂,Ŝ,τ̂
end
# orthonormal wave speeds (v̂ ± c_s relativistic composition)
@inline function _λpm(cs2,v̂)
    cs=sqrt(cs2); v2=clamp(v̂^2,0.0,1-1e-14)
    a=1/(1-v2*cs2); disc=sqrt(max(cs2*(1-v2)*(1-v2*cs2-v2*(1-cs2)),0.0))
    return a*(v̂*(1-cs2)-disc), a*(v̂*(1-cs2)+disc)
end

# Γ-law cons2prim with the Galeazzi-type fixing. Inputs are UNdensitized; p0 is a starting
# guess (the node's previous pressure). Returns (ρ, ε, v̂, p, code): code 0 = clean, 1 = reset
# to atmosphere, 2 = conserved variables altered (τ<0 / S bound / entropy bounds).
@inline function _c2p_resid(Γ,D̂,Ŝ,τ̂,p)
    v=Ŝ/(τ̂+D̂+p); W=1/sqrt(max(1-v^2,1e-16)); ρ=D̂/W
    ε=(τ̂+D̂*(1-W)+p*(1-W^2))/(D̂*W)
    return (Γ-1)*ρ*ε-p, ρ, ε, v, W
end
function _c2p(eng::DGStarHPEngine, D̂::Float64, Ŝx::Float64, τ̂in::Float64, elam::Float64, p0::Float64)
    Γ=eng.Γ; κ=eng.κ
    if !(D̂ > eng.ρ_cut) || !isfinite(D̂) || !isfinite(Ŝx) || !isfinite(τ̂in)
        return eng.ρ_atm, eng.ε_atm, 0.0, eng.p_atm, 1
    end
    code=0
    Ŝ=abs(Ŝx)/sqrt(elam)                                  # orthonormal momentum magnitude
    τ̂ = τ̂in < 0.0 ? 1e-300 : τ̂in; τ̂ != τ̂in && (code=2)
    Smax=sqrt(max(τ̂*(τ̂+2D̂),0.0))*(1-1e-10)               # S² < τ(τ+2D): the cold-fluid bound
    if Ŝ > Smax; Ŝ=Smax; code=2; end
    # Newton on f(p) = (Γ−1)ρε − p from p0 (derivative by finite difference), bisection fallback
    p=max(p0, 1e-300); ok=false
    for it in 1:30
        f,_,_,_,_=_c2p_resid(Γ,D̂,Ŝ,τ̂,p)
        dp=max(1e-8*p, 1e-300)
        f2,_,_,_,_=_c2p_resid(Γ,D̂,Ŝ,τ̂,p+dp)
        df=(f2-f)/dp
        df ≥ 0 && break
        pn=p-f/df
        pn ≤ 0 && (pn=0.5*p)
        if abs(pn-p) ≤ 1e-13*max(pn,1e-300); p=pn; ok=true; break; end
        p=pn
    end
    if !ok
        plo=0.0; phi=max((Γ-1)*τ̂, 1e-300)
        while _c2p_resid(Γ,D̂,Ŝ,τ̂,phi)[1] > 0 && phi < 1e10; phi *= 2; end
        for _ in 1:200
            pm=0.5*(plo+phi); fm=_c2p_resid(Γ,D̂,Ŝ,τ̂,pm)[1]
            if fm > 0; plo=pm else phi=pm end
            (phi-plo) ≤ 1e-13*phi && break
        end
        p=0.5*(plo+phi)
    end
    _,ρ,ε,v̂,W=_c2p_resid(Γ,D̂,Ŝ,τ̂,p)
    if !(ρ > eng.ρ_atm)
        return eng.ρ_atm, eng.ε_atm, 0.0, eng.p_atm, 1
    end
    εmin=eng.entropy_floor ? _εpoly(κ,Γ,ρ) : 0.0; εmax=eng.εfac_max*_εpoly(κ,Γ,ρ)   # entropy bounds
    if ε < εmin; ε=εmin; code=2
    elseif ε > εmax; ε=εmax; code=2 end
    p=_p_ig(Γ,ρ,ε)
    return ρ,ε,copysign(v̂,Ŝx),p,code
end

function _update_prims!(st::DGStarHPState, eng::DGStarHPEngine)
    @inbounds for i in 1:eng.Ntot
        sg=eng.sqrtγ[i]; el=eng.elam[i]
        ρ,ε,v̂,p,code=_c2p(eng, st.D[i]/sg, st.S[i]/sg, st.τ[i]/sg, el, st.p[i])
        wasatm = st.ρ[i] ≤ eng.ρ_atm
        st.ρ[i]=ρ; st.ε[i]=ε; st.v[i]=v̂; st.p[i]=p
        if code != 0
            D̂,Ŝ,τ̂=_p2c(eng.Γ,ρ,ε,v̂,el)
            st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂
            (code == 2 || !wasatm) && (st.nfix+=1)          # count only genuine fixes
        end
    end
end

# ----------------------------------------------------------------------------------
# setup
# ----------------------------------------------------------------------------------
"""
    setup_dgstarhp(eos, εc; grid=:I1, Γ=2.0, κ=eos.κ, atm_fac=1e-13, cut_fac=10.0,
                   εfac_max=100.0, flux=:hll, limiter=:minmod, p_lim=2,
                   filter=:none, p_filt=3, filt_α=36.0, filt_s=6, filter_deviation=true,
                   wellbalanced=true, symmetrize=true, entropy_floor=true, p_center=7, center_fac=3, cfl=0.25, xmax_fac=3.0, h_tov=2e-4) -> (engine, state)

Build the hp-adapted radial DG star on the TOV star of central energy density `εc`
(G=c=M⊙=1). `wellbalanced=false` drops the static-residual subtraction (the paper's scheme:
the star then settles to its numerical equilibrium). NOTE: the subtraction only pays off with a limiter that leaves the
projected equilibrium invariant (:wb, :mean, :none): with :minmod, which rewrites the surface
elements' slopes whenever a neighbour estimate is smaller, the star settles to err[D̃] ≈ 1e-3
with or without it (I1 grid, t=2000: 1.0e-3 vs 0.9e-3). The paper's scheme is
wellbalanced=false.

`p_center` (odd, default 7) is the degree of the central element: with p=3 the flux x²p(x)
of the static star has an unrepresentable x⁴ term whose aliased derivative is as large as the
gravitational force at the innermost nodes (static momentum residual 3 × gravity, on every
grid); the p=3 element next to it has the same defect at its inner node (0.6), so the central
element spans `center_fac` (odd, default 3) nominal interior widths. Measured static momentum
residual / gravity (I1 grid): centre 3.0 (p=3) → 3.8e-2 (p=5, fac 3) → 1.7e-4 (p=7, fac 3);
first interior element 0.57 → 1.8e-2 with fac 3.
Linear (p=1) surface elements likewise carry a static residual
of 0.5 of gravity because the quadratic pressure p∝(R−r)² is not representable; quadratic
surface elements (:I2) hold it exactly.

`symmetrize=true` (default) projects the state onto the physical reflection parity (D̃, τ̃ even,
S̃ odd) after every stage. The doubled symmetric domain admits an UNPHYSICAL antisymmetric
mode (odd δρ, even v: a translation of the star) that the spherical problem does not have; on
every grid and for every flux, CFL and limiter tried it grows from round-off at 0.049/M⊙
(e-folding 20 M⊙, the dynamical time), first visible at the central element, and destroys the
star by t≈600 M⊙; its measured antisymmetric fraction is >1 (pure antisymmetry gives √2).
Enforcing the parity removes it exactly. The exponential filter (α=36, s=6 on ALL high-order
elements, the paper's central-cube setting) also suppresses it, but only by damping the linear
mode of every element by 5% per step, and s≥12 does not. `entropy_floor=false` lifts the lower
entropy bound (a test switch; it is not the cause of the instability). `grid` is a preset symbol for [`hp_grid`](@ref) or an explicit region list.
The evolution EOS is the Γ-law p=(Γ−1)ρε with the polytrope (κ, Γ) as initial adiabat and
as the lower entropy bound. `atm_fac` sets ρ_atm = atm_fac·ρ_c (the paper: ~1e-13),
`cut_fac` sets ρ_cut = cut_fac·ρ_atm.
"""
function setup_dgstarhp(eos::BarotropicEOS, εc::Float64; grid=:I1, Γ::Float64=2.0,
        κ::Float64=(eos isa ShumPolytrope ? eos.κ : 100.0), atm_fac::Float64=1e-13, cut_fac::Float64=10.0,
        εfac_max::Float64=100.0, flux::Symbol=:hll, limiter::Symbol=:minmod, p_lim::Int=2,
        filter::Symbol=:none, p_filt::Int=3, filt_α::Float64=36.0, filt_s::Int=6, filter_deviation::Bool=true,
        wellbalanced::Bool=true, symmetrize::Bool=true, entropy_floor::Bool=true, p_center::Int=7, center_fac::Int=3, cfl::Float64=0.25, xmax_fac::Float64=3.0, h_tov::Float64=2e-4)
    star=solve_tov(eos,εc;h=h_tov); R,M=star.R,star.M
    regions = grid isa Symbol ? hp_grid(R; preset=grid, xmax_fac=xmax_fac) : grid
    bases=Dict{Int,LGLBasis}()
    elems,Ntot=_build_elems(regions,bases; p_center=p_center, center_fac=center_fac)
    K=length(elems); kc=findfirst(e->e.region==:center, elems)
    rt=star.r; mt=star.m; νt=star.ν; pt=star.p; et=star.ε
    m_of(r)= r<R ? _lin(rt,mt,r) : M
    ν_of(r)= r<R ? _lin(rt,νt,r) : log(1-2M/r)
    p_of(r)= r<R ? max(_lin(rt,pt,r),0.0) : 0.0
    e_of(r)= r<R ? max(_lin(rt,et,r),0.0) : 0.0
    x=zeros(Ntot); r=zeros(Ntot); α=zeros(Ntot); elam=zeros(Ntot); sqrtγ=zeros(Ntot); Φx=zeros(Ntot); λx=zeros(Ntot)
    dxmin=Inf
    for e in elems
        b=bases[e.p]
        for a in 1:e.N
            i=e.off+a-1
            xx=0.5*(e.xlo+e.xhi)+e.J*b.ξ[a]; x[i]=xx; rr=abs(xx); r[i]=rr
            m=m_of(rr); f=max(1-2m/rr,1e-12)
            α[i]=exp(0.5*ν_of(rr)); elam[i]=1/f; sqrtγ[i]=sqrt(elam[i])*rr^2
            Φp=(m+4π*rr^3*p_of(rr))/(rr*(rr-2m))
            mp=4π*rr^2*e_of(rr)
            λp=(2mp/rr-2m/rr^2)/f
            Φx[i]=sign(xx)*Φp; λx[i]=sign(xx)*λp
            a<e.N && (dxmin=min(dxmin, e.J*(b.ξ[a+1]-b.ξ[a])))
        end
    end
    ρc=rho_from_p(eos,pressure(eos,εc)); ρ_atm=atm_fac*ρc; ε_atm=_εpoly(κ,Γ,ρ_atm); p_atm=_p_ig(Γ,ρ_atm,ε_atm)
    st=DGStarHPState(Ntot)
    for i in 1:Ntot
        rr=r[i]
        # nodes at (or within round-off of) the surface, and any node whose projected density
        # falls below the atmosphere cut, are atmosphere: a "fluid" node with ρ ~ 1e-10 ρ_c has
        # c_s ~ 1e-6 and is shock-heated to the entropy cap by round-off velocities, and that
        # is enough to blow the surface element (seen: e-folding ~100 M⊙ from the r=R node).
        ρ = rr < R*(1-1e-10) ? max(rho_from_p(eos,p_of(rr)), ρ_atm) : ρ_atm
        ρ ≤ cut_fac*ρ_atm && (ρ=ρ_atm)
        ε=_εpoly(κ,Γ,ρ)
        st.ρ[i]=ρ; st.ε[i]=ε; st.p[i]=_p_ig(Γ,ρ,ε); st.v[i]=0.0
        D̂,Ŝ,τ̂=_p2c(Γ,ρ,ε,0.0,elam[i]); sg=sqrtγ[i]
        st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂
    end
    eng=DGStarHPEngine(elems,bases,K,Ntot,kc,Γ,κ,eos,R,M,ρ_atm,cut_fac*ρ_atm,ε_atm,p_atm,εfac_max,
                       flux,limiter,p_lim,filter,p_filt,filt_α,filt_s,filter_deviation,wellbalanced,symmetrize,entropy_floor,cfl,dxmin,
                       x,r,α,elam,sqrtγ,Φx,λx, copy(st.D),copy(st.S),copy(st.τ), zeros(Ntot),zeros(Ntot),zeros(Ntot))
    wellbalanced && _raw_rhs!(eng.Seq_D,eng.Seq_S,eng.Seq_τ, st, eng)
    return eng, st
end

@inline function _lin(xs,ys,x)
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end

# ----------------------------------------------------------------------------------
# RHS
# ----------------------------------------------------------------------------------
# densitized state and α-weighted flux at a node from its primitives
@inline function _qf(eng::DGStarHPEngine, i::Int, ρ, ε, v̂)
    Γ=eng.Γ; el=eng.elam[i]; sg=eng.sqrtγ[i]; α=eng.α[i]
    p=_p_ig(Γ,ρ,ε); h=1+ε+p/ρ; W=1/sqrt(1-clamp(v̂^2,0.0,1-1e-14))
    vx=v̂/sqrt(el)                                         # coordinate velocity v^x
    D̃=sg*ρ*W; S̃=sg*ρ*h*W^2*sqrt(el)*v̂; τ̃=sg*(ρ*h*W^2-p)-D̃
    q=(D̃,S̃,τ̃); f=(α*D̃*vx, α*(S̃*vx+sg*p), α*(τ̃*vx+sg*p*vx))
    cs2=_cs2(Γ,ρ,ε); λm,λp=_λpm(cs2,v̂); c=α/sqrt(el)      # coordinate speeds
    return q,f,c*λm,c*λp
end

function _raw_rhs!(rD,rS,rτ, st::DGStarHPState, eng::DGStarHPEngine)
    Γ=eng.Γ; Ntot=eng.Ntot
    fD=Vector{Float64}(undef,Ntot); fS=similar(fD); fτ=similar(fD)
    @inbounds for i in 1:Ntot
        _,f,_,_=_qf(eng,i,st.ρ[i],st.ε[i],st.v[i])
        fD[i]=f[1]; fS[i]=f[2]; fτ[i]=f[3]
    end
    # volume term: −(1/J) D f
    @inbounds for e in eng.elems
        Dm=eng.bases[e.p].D; N=e.N; o=e.off-1
        for a in 1:N
            sD=0.0; sS=0.0; sτ=0.0
            for c in 1:N
                d=Dm[a,c]; sD+=d*fD[o+c]; sS+=d*fS[o+c]; sτ+=d*fτ[o+c]
            end
            rD[o+a]=-sD/e.J; rS[o+a]=-sS/e.J; rτ[o+a]=-sτ/e.J
        end
    end
    # faces kf = 0..K (kf=0: left of element 1; kf=K: right of element K)
    @inbounds for kf in 0:eng.K
        if kf==0
            eR=eng.elems[1]; iR=eR.off
            qR,fR,smR,spR=_qf(eng,iR,st.ρ[iR],st.ε[iR],st.v[iR])
            qL,fL,smL,spL=_qf(eng,iR,eng.ρ_atm,eng.ε_atm,0.0)          # outflow: atmosphere at rest
        elseif kf==eng.K
            eL=eng.elems[eng.K]; iL=eL.off+eL.N-1
            qL,fL,smL,spL=_qf(eng,iL,st.ρ[iL],st.ε[iL],st.v[iL])
            qR,fR,smR,spR=_qf(eng,iL,eng.ρ_atm,eng.ε_atm,0.0)
        else
            eL=eng.elems[kf]; eR=eng.elems[kf+1]; iL=eL.off+eL.N-1; iR=eR.off
            qL,fL,smL,spL=_qf(eng,iL,st.ρ[iL],st.ε[iL],st.v[iL])
            qR,fR,smR,spR=_qf(eng,iR,st.ρ[iR],st.ε[iR],st.v[iR])
        end
        if eng.flux==:hll
            cmin=min(smL,smR,0.0); cmax=max(spL,spR,0.0)
            if cmax-cmin ≤ 1e-300
                F̂=(0.5*(fL[1]+fR[1]),0.5*(fL[2]+fR[2]),0.5*(fL[3]+fR[3]))
            else
                inv=1/(cmax-cmin)
                F̂=((cmax*fL[1]-cmin*fR[1]+cmin*cmax*(qR[1]-qL[1]))*inv, (cmax*fL[2]-cmin*fR[2]+cmin*cmax*(qR[2]-qL[2]))*inv,
                   (cmax*fL[3]-cmin*fR[3]+cmin*cmax*(qR[3]-qL[3]))*inv)
            end
        else
            amax=max(abs(smL),abs(spL),abs(smR),abs(spR))
            F̂=(0.5*(fL[1]+fR[1])-0.5*amax*(qR[1]-qL[1]), 0.5*(fL[2]+fR[2])-0.5*amax*(qR[2]-qL[2]), 0.5*(fL[3]+fR[3])-0.5*amax*(qR[3]-qL[3]))
        end
        if kf≥1
            eL=eng.elems[kf]; iL=eL.off+eL.N-1; wN=eng.bases[eL.p].w[eL.N]
            rD[iL]-=(F̂[1]-fL[1])/(eL.J*wN); rS[iL]-=(F̂[2]-fL[2])/(eL.J*wN); rτ[iL]-=(F̂[3]-fL[3])/(eL.J*wN)
        end
        if kf≤eng.K-1
            eR=eng.elems[kf+1]; iR=eR.off; w1=eng.bases[eR.p].w[1]
            rD[iR]+=(F̂[1]-fR[1])/(eR.J*w1); rS[iR]+=(F̂[2]-fR[2])/(eR.J*w1); rτ[iR]+=(F̂[3]-fR[3])/(eR.J*w1)
        end
    end
    # sources
    @inbounds for i in 1:Ntot
        ρ=st.ρ[i]; ε=st.ε[i]; v̂=st.v[i]; p=st.p[i]
        h=1+ε+p/ρ; W2=1/(1-clamp(v̂^2,0.0,1-1e-14)); el=eng.elam[i]
        α=eng.α[i]; sg=eng.sqrtγ[i]; xx=eng.x[i]
        vx=v̂/sqrt(el)
        rS[i]+= α*sg*((ρ*h*W2*v̂^2+p)*0.5*eng.λx[i] + 2p/xx) - α*sg*(ρ*h*W2-p)*eng.Φx[i]
        rτ[i]+= -α*sg*ρ*h*W2*vx*eng.Φx[i]
    end
    return nothing
end

function _rhs!(rD,rS,rτ, st::DGStarHPState, eng::DGStarHPEngine)
    _raw_rhs!(rD,rS,rτ, st, eng)
    @inbounds for i in 1:eng.Ntot
        rD[i]-=eng.Seq_D[i]; rS[i]-=eng.Seq_S[i]; rτ[i]-=eng.Seq_τ[i]
    end
end

# ----------------------------------------------------------------------------------
# limiters
# ----------------------------------------------------------------------------------
@inline function _mean(eng::DGStarHPEngine, U, e::HPElem)
    w=eng.bases[e.p].w; s=0.0
    @inbounds for a in 1:e.N; s+=w[a]*U[e.off+a-1]; end
    0.5*s
end
# L2 slope du/dξ of the linear part, from the orthonormal-Legendre modal coefficient û₁
# (ψ₁ = √(3/2) ξ). NOT the LGL-quadrature formula 1.5 Σ w ξ u: the N-point LGL rule integrates
# ξ·u (degree p+1) exactly only for p ≥ 2, and for p=1 it overestimates the slope by 3×, which
# made the minmod limiter replace the surface elements' slopes by the neighbour estimates on
# every stage (a bug that destroyed the star on the I1 grid).
@inline function _slope(eng::DGStarHPEngine, U, e::HPElem)
    b=eng.bases[e.p]; s=0.0
    @inbounds for a in 1:e.N; s+=b.invV[2,a]*U[e.off+a-1]; end
    s*sqrt(1.5)
end
@inline _minmod3(a,b,c) = (sign(a)==sign(b)==sign(c)) ? sign(a)*min(abs(a),abs(b),abs(c)) : 0.0

# physical-state check of an element's nodes (undensitized): D>ρ_atm, τ>0, S²<τ(τ+2D)
function _physical(st::DGStarHPState, eng::DGStarHPEngine, e::HPElem)
    @inbounds for a in 1:e.N
        i=e.off+a-1; sg=eng.sqrtγ[i]
        D̂=st.D[i]/sg; τ̂=st.τ[i]/sg; Ŝ2=(st.S[i]/sg)^2/eng.elam[i]
        (D̂ ≥ (1-1e-8)*eng.ρ_atm && τ̂ > 0.0 && Ŝ2 < τ̂*(τ̂+2D̂)) || return false
    end
    return true
end

# the paper's minmod ΛΠ¹: slope estimates from the element's own linear mode and the
# neighbouring means; if the smallest is not the element's own, replace the solution by
# mean + limited slope (all higher modes dropped). Then the physical-state check with slope
# reduction. Applied to elements with p ≤ p_lim only.
function _limit_minmod!(st::DGStarHPState, eng::DGStarHPEngine)
    K=eng.K
    for k in 1:K
        e=eng.elems[k]; e.p ≤ eng.p_lim || continue
        b=eng.bases[e.p]
        em = k>1 ? eng.elems[k-1] : e; ep = k<K ? eng.elems[k+1] : e
        for U in (st.D, st.S, st.τ)
            ū=_mean(eng,U,e); a1=_slope(eng,U,e)/e.J                    # du/dx of the linear mode
            ūm = k>1 ? _mean(eng,U,em) : ū; ūp = k<K ? _mean(eng,U,ep) : ū
            a2=(ūp-ū)/(0.5*e.h); a3=(ū-ūm)/(0.5*e.h)
            a=_minmod3(a1,a2,a3)
            if a != a1                       # limiter activates: mean + limited slope, higher modes dropped
                @inbounds for i in 1:e.N; U[e.off+i-1]=ū+a*e.J*b.ξ[i]; end
            end
        end
        # physical-state check: halve the slopes (about the means) until every node is admissible
        if !_physical(st,eng,e)
            for it in 1:40
                for U in (st.D, st.S, st.τ)
                    ū=_mean(eng,U,e)
                    @inbounds for i in 1:e.N; j=e.off+i-1; U[j]=ū+0.5*(U[j]-ū); end
                end
                _physical(st,eng,e) && break
            end
        end
    end
end

# mean-based positivity scaling (the Zhang–Shu / MRS-type rescaling), U ← Ū + θ(U−Ū),
# θ from D ≥ √γρ_atm and the pressure proxy q = τ + D − √(D² + S²/e^λ) > 0.
@inline _q(D,S,τ,el) = (τ+D) - sqrt(D^2 + S^2/el)
function _limit_mean!(st::DGStarHPState, eng::DGStarHPEngine)
    for e in eng.elems
        e.p ≤ eng.p_lim || continue
        D̄=_mean(eng,st.D,e); S̄=_mean(eng,st.S,e); τ̄=_mean(eng,st.τ,e)
        θ=1.0
        @inbounds for a in 1:e.N
            i=e.off+a-1; Dfl=eng.sqrtγ[i]*eng.ρ_atm
            if st.D[i] < Dfl
                θ=min(θ, D̄>Dfl ? clamp((D̄-Dfl)/(D̄-st.D[i]+1e-300),0.0,1.0) : 0.0)
            end
        end
        if θ<1; @inbounds for a in 1:e.N; i=e.off+a-1; st.D[i]=D̄+θ*(st.D[i]-D̄); end; end
        θq=1.0; sgm=_mean(eng,eng.sqrtγ,e)
        q̄=_q(D̄,S̄,τ̄,eng.elam[e.off]); qε=1e-12*max(q̄,0.0)
        @inbounds for a in 1:e.N
            i=e.off+a-1; qn=_q(st.D[i],st.S[i],st.τ[i],eng.elam[i])
            if qn < qε
                θq=min(θq, q̄>qε ? clamp((q̄-qε)/(q̄-qn+1e-300),0.0,1.0) : 0.0)
            end
        end
        if θq<1
            @inbounds for a in 1:e.N
                i=e.off+a-1
                st.D[i]=D̄+θq*(st.D[i]-D̄); st.S[i]=S̄+θq*(st.S[i]-S̄); st.τ[i]=τ̄+θq*(st.τ[i]-τ̄)
            end
        end
    end
end

# the same scaling about the EQUILIBRIUM-PRESERVING reference (DGCart3D, VALIDATION §7.8):
# R_D = D_eq + δD̄·w_D, R_S = S_eq + δS̄·w_S, R_τ = τ_eq + ΔK + (δτ̄ − ΔK̄)·w_τ.
function _limit_wb!(st::DGStarHPState, eng::DGStarHPEngine)
    for e in eng.elems
        e.p ≤ eng.p_lim || continue
        N=e.N; o=e.off-1; w=eng.bases[e.p].w
        mD=0.0; mS=0.0; mτ=0.0; Deqm=0.0; τeqm=0.0; wSm=0.0
        @inbounds for a in 1:N
            i=o+a
            mD+=0.5*w[a]*(st.D[i]-eng.Deq[i]); mS+=0.5*w[a]*(st.S[i]-eng.Seq[i]); mτ+=0.5*w[a]*(st.τ[i]-eng.τeq[i])
            Deqm+=0.5*w[a]*eng.Deq[i]; τeqm+=0.5*w[a]*eng.τeq[i]; wSm+=0.5*w[a]*(eng.Deq[i]-eng.sqrtγ[i]*eng.ρ_atm)
        end
        RD=zeros(N); RS=zeros(N); Rτ=zeros(N); ΔK=zeros(N); ΔKm=0.0
        @inbounds for a in 1:N
            i=o+a; wD=eng.Deq[i]/Deqm; wsS = wSm>0 ? (eng.Deq[i]-eng.sqrtγ[i]*eng.ρ_atm)/wSm : wD
            RD[a]=eng.Deq[i]+mD*wD; RS[a]=eng.Seq[i]+mS*wsS
            ΔK[a]=_q(eng.Deq[i],eng.Seq[i],0.0,eng.elam[i]) - _q(RD[a],RS[a],0.0,eng.elam[i])
            ΔKm+=0.5*w[a]*ΔK[a]
        end
        @inbounds for a in 1:N
            i=o+a; Rτ[a]=eng.τeq[i]+ΔK[a]+(mτ-ΔKm)*eng.τeq[i]/τeqm
        end
        # step 1: D ≥ ½√γρ_atm
        θ=1.0; feasible=true
        @inbounds for a in 1:N
            i=o+a; Dfl=0.5*eng.sqrtγ[i]*eng.ρ_atm
            if st.D[i] < Dfl
                RD[a] > Dfl || (feasible=false; break)
                θ=min(θ, clamp((RD[a]-Dfl)/(RD[a]-st.D[i]+1e-300),0.0,1.0))
            end
        end
        if !feasible; _limit_mean_elem!(st,eng,e); continue; end
        if θ<1; @inbounds for a in 1:N; i=o+a; st.D[i]=RD[a]+θ*(st.D[i]-RD[a]); end; end
        # step 2: pressure proxy about the same reference
        θq=1.0
        @inbounds for a in 1:N
            i=o+a; qR=_q(RD[a],RS[a],Rτ[a],eng.elam[i])
            qR > 0 || (feasible=false; break)
            qn=_q(st.D[i],st.S[i],st.τ[i],eng.elam[i]); qε=1e-12*qR
            qn < qε && (θq=min(θq, clamp((qR-qε)/(qR-qn+1e-300),0.0,1.0)))
        end
        if !feasible; _limit_mean_elem!(st,eng,e); continue; end
        if θq<1
            @inbounds for a in 1:N
                i=o+a
                st.D[i]=RD[a]+θq*(st.D[i]-RD[a]); st.S[i]=RS[a]+θq*(st.S[i]-RS[a]); st.τ[i]=Rτ[a]+θq*(st.τ[i]-Rτ[a])
            end
        end
    end
end
function _limit_mean_elem!(st::DGStarHPState, eng::DGStarHPEngine, e::HPElem)
    D̄=_mean(eng,st.D,e); S̄=_mean(eng,st.S,e); τ̄=_mean(eng,st.τ,e); θ=1.0
    @inbounds for a in 1:e.N
        i=e.off+a-1; Dfl=eng.sqrtγ[i]*eng.ρ_atm
        st.D[i] < Dfl && (θ=min(θ, D̄>Dfl ? clamp((D̄-Dfl)/(D̄-st.D[i]+1e-300),0.0,1.0) : 0.0))
    end
    θq=1.0; q̄=_q(D̄,S̄,τ̄,eng.elam[e.off]); qε=1e-12*max(q̄,0.0)
    @inbounds for a in 1:e.N
        i=e.off+a-1; Dn=D̄+θ*(st.D[i]-D̄); qn=_q(Dn,st.S[i],st.τ[i],eng.elam[i])
        qn < qε && (θq=min(θq, q̄>qε ? clamp((q̄-qε)/(q̄-qn+1e-300),0.0,1.0) : 0.0))
    end
    θ=min(θ,θq)
    @inbounds for a in 1:e.N
        i=e.off+a-1
        st.D[i]=D̄+θ*(st.D[i]-D̄); st.S[i]=S̄+θ*(st.S[i]-S̄); st.τ[i]=τ̄+θ*(st.τ[i]-τ̄)
    end
end

# project onto the physical reflection parity: D̃, τ̃ even in x, S̃_x odd. Element k mirrors
# K+1−k and node a mirrors N+1−a. This removes the unphysical antisymmetric modes that the
# doubled domain admits but the spherical problem does not.
function _symmetrize!(st::DGStarHPState, eng::DGStarHPEngine)
    K=eng.K
    @inbounds for k in 1:(K÷2)
        e=eng.elems[k]; m=eng.elems[K+1-k]; N=e.N
        for a in 1:N
            i=e.off+a-1; j=m.off+N-a
            d=0.5*(st.D[i]+st.D[j]); st.D[i]=d; st.D[j]=d
            t=0.5*(st.τ[i]+st.τ[j]); st.τ[i]=t; st.τ[j]=t
            sv=0.5*(st.S[i]-st.S[j]); st.S[i]=sv; st.S[j]=-sv
        end
    end
    e=eng.elems[eng.kc]; N=e.N                          # central element: mirror within itself
    @inbounds for a in 1:(N÷2)
        i=e.off+a-1; j=e.off+N-a
        d=0.5*(st.D[i]+st.D[j]); st.D[i]=d; st.D[j]=d
        t=0.5*(st.τ[i]+st.τ[j]); st.τ[i]=t; st.τ[j]=t
        sv=0.5*(st.S[i]-st.S[j]); st.S[i]=sv; st.S[j]=-sv
    end
end

function _limit!(st::DGStarHPState, eng::DGStarHPEngine)
    eng.symmetrize && _symmetrize!(st,eng)
    if eng.limiter==:minmod; _limit_minmod!(st,eng)
    elseif eng.limiter==:mean; _limit_mean!(st,eng)
    elseif eng.limiter==:wb; _limit_wb!(st,eng)
    elseif eng.limiter==:none; nothing
    else error("unknown limiter $(eng.limiter)") end
end

# exponential modal filter on elements with p ≥ p_filt (on the deviation from equilibrium
# unless filter_deviation=false)
function _filter!(st::DGStarHPState, eng::DGStarHPEngine)
    eng.filter == :none && return
    for e in eng.elems
        e.p ≥ eng.p_filt || continue
        eng.filter == :center && e.region != :center && continue
        b=eng.bases[e.p]; N=e.N; o=e.off-1
        σ=[exp(-eng.filt_α*((i-1)/e.p)^eng.filt_s) for i in 1:N]
        for (U,Ueq) in ((st.D,eng.Deq),(st.S,eng.Seq),(st.τ,eng.τeq))
            u=[eng.filter_deviation ? U[o+a]-Ueq[o+a] : U[o+a] for a in 1:N]
            û=b.invV*u; û.*=σ; u=b.V*û
            @inbounds for a in 1:N; U[o+a]=(eng.filter_deviation ? Ueq[o+a] : 0.0)+u[a]; end
        end
    end
end

# ----------------------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------------------
"""central rest-mass density: the interpolant of ρ at x=0 in the central element."""
function dgstarhp_central_density(st::DGStarHPState, eng::DGStarHPEngine)
    e=eng.elems[eng.kc]; b=eng.bases[e.p]; N=e.N; s=0.0
    @inbounds for a in 1:N                                   # Lagrange basis at ξ=0
        ℓ=1.0
        for c in 1:N; c==a && continue; ℓ*=(0.0-b.ξ[c])/(b.ξ[a]-b.ξ[c]); end
        s+=ℓ*st.ρ[e.off+a-1]
    end
    s
end
"""normalized L2 error of D̃ against the projected TOV star over all nodes (their eqs. 35–36)."""
dgstarhp_errD(st::DGStarHPState, eng::DGStarHPEngine) =
    sqrt(sum(abs2, st.D .- eng.Deq)/sum(abs2, eng.Deq))
"""baryon mass M_b = 4π∫D̃ dr = 2π∫D̃ dx over the symmetric domain (LGL quadrature)."""
function dgstarhp_baryon_mass(st::DGStarHPState, eng::DGStarHPEngine)
    s=0.0
    for e in eng.elems
        w=eng.bases[e.p].w
        @inbounds for a in 1:e.N; s+=w[a]*e.J*st.D[e.off+a-1]; end
    end
    2π*s
end
"""(max|v̂| in the atmosphere region r>R, max r with ρ>10ρ_atm) — the paper's 'star extends
beyond its surface' diagnostic (their Fig. 10 discussion)."""
function dgstarhp_surface_state(st::DGStarHPState, eng::DGStarHPEngine)
    vmax=0.0; rmax=0.0
    @inbounds for i in 1:eng.Ntot
        eng.r[i] > eng.R && (vmax=max(vmax,abs(st.v[i])))
        st.ρ[i] > 10eng.ρ_atm && (rmax=max(rmax,eng.r[i]))
    end
    return vmax, rmax
end

# ----------------------------------------------------------------------------------
# evolution
# ----------------------------------------------------------------------------------
@inline function _maxspeed(st::DGStarHPState, eng::DGStarHPEngine)
    a=0.0
    @inbounds for i in 1:eng.Ntot
        cs2=_cs2(eng.Γ,st.ρ[i],st.ε[i]); λm,λp=_λpm(cs2,st.v[i]); c=eng.α[i]/sqrt(eng.elam[i])
        a=max(a,c*abs(λm),c*abs(λp))
    end
    a
end

"""
    evolve_dgstarhp!(st, eng; tmax, sample_dt=1.0, dt=0.0) -> rec

SSP-RK3 evolution. `dt=0` picks Δt = cfl·Δx_min/a_max at t=0 and keeps it fixed (the
paper's practice); the limiter runs after every substep, the filter after every full step.
Returns a NamedTuple of vectors sampled every `sample_dt`: `t, ρc, errD, Snorm, Mb, vatm,
rstar` (max atmosphere speed and radius of the ρ>10ρ_atm region) plus `nfix`.
"""
function evolve_dgstarhp!(st::DGStarHPState, eng::DGStarHPEngine; tmax::Float64, sample_dt::Float64=1.0,
                          dt::Float64=0.0)
    n=eng.Ntot
    rD=zeros(n); rS=zeros(n); rτ=zeros(n); D0=zeros(n); S0=zeros(n); τ0=zeros(n)
    _update_prims!(st,eng)
    dt = dt>0 ? dt : eng.cfl*eng.dxmin/max(_maxspeed(st,eng),1e-3)
    ts=Float64[]; ρc=Float64[]; errD=Float64[]; Sn=Float64[]; Mb=Float64[]; vatm=Float64[]; rstar=Float64[]
    record!(t) = begin
        push!(ts,t); push!(ρc,dgstarhp_central_density(st,eng)); push!(errD,dgstarhp_errD(st,eng))
        push!(Sn,sqrt(sum(abs2,st.S)/n)); push!(Mb,dgstarhp_baryon_mass(st,eng))
        va,rs=dgstarhp_surface_state(st,eng); push!(vatm,va); push!(rstar,rs)
    end
    record!(0.0)
    t=0.0; last=0.0; ns=0
    stage!(c0,c1,h) = begin            # U ← c0·U0 + c1·(U + h·rhs(U))
        _rhs!(rD,rS,rτ,st,eng)
        @inbounds for i in 1:n
            st.D[i]=c0*D0[i]+c1*(st.D[i]+h*rD[i]); st.S[i]=c0*S0[i]+c1*(st.S[i]+h*rS[i]); st.τ[i]=c0*τ0[i]+c1*(st.τ[i]+h*rτ[i])
        end
        _limit!(st,eng); _update_prims!(st,eng)
    end
    while t < tmax-1e-12
        h=min(dt,tmax-t)
        copyto!(D0,st.D); copyto!(S0,st.S); copyto!(τ0,st.τ)
        stage!(0.0,1.0,h); stage!(0.75,0.25,h); stage!(1/3,2/3,h)
        _filter!(st,eng); eng.filter != :none && _update_prims!(st,eng)
        t+=h; ns+=1
        if t-last ≥ sample_dt-1e-12
            record!(t); last=t
            (isfinite(ρc[end]) && ρc[end] < 100*ρc[1]) || break
        end
    end
    return (t=ts, ρc=ρc, errD=errD, Snorm=Sn, Mb=Mb, vatm=vatm, rstar=rstar, nfix=st.nfix, dt=dt, nsteps=ns)
end

"""
    seed_dgstarhp_radial!(st, eng; A=1e-3, profile=:sine)

Seed a radial velocity inside the star (odd in x): `profile=:sine` gives v̂ = A sin(π r/R),
which vanishes at the surface (the thin surface elements are not kicked directly);
`:linear` gives the homologous v̂ = A r/R.
"""
function seed_dgstarhp_radial!(st::DGStarHPState, eng::DGStarHPEngine; A::Float64=1e-3, profile::Symbol=:sine)
    @inbounds for i in 1:eng.Ntot
        rr=eng.r[i]; rr < eng.R || continue
        st.ρ[i] > 10eng.ρ_atm || continue
        shape = profile==:sine ? sin(π*rr/eng.R) : rr/eng.R
        v̂=A*shape*sign(eng.x[i]); st.v[i]=v̂
        D̂,Ŝ,τ̂=_p2c(eng.Γ,st.ρ[i],st.ε[i],v̂,eng.elam[i]); sg=eng.sqrtγ[i]
        st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂
    end
end

"""
    dgstarhp_spectrum(ts, ρc; fmin_kHz=1.0, fmax_kHz=15.0, npts=6000, npeaks=6, tmax=Inf)
        -> (peaks_kHz, νgrid_kHz, power)

Hanning-windowed periodogram of ρ_c(t) (the paper's procedure, Fig. 13) over t ≤ tmax;
peaks sorted by frequency. Times in M⊙.
"""
function dgstarhp_spectrum(ts::Vector{Float64}, ρc::Vector{Float64}; fmin_kHz=1.0, fmax_kHz=15.0,
                           npts::Int=6000, npeaks::Int=6, tmax::Float64=Inf, prominence::Float64=3.0)
    m=ts .≤ tmax; t=ts[m]; q=ρc[m]; n=length(q); μ=sum(q)/n
    @inbounds for i in 1:n; q[i]=(q[i]-μ)*(0.5-0.5*cos(2π*(i-1)/(n-1))); end
    conv=kHz_to_km*Msun_to_km                                # kHz → 1/M⊙
    νgrid=collect(range(fmin_kHz*conv, fmax_kHz*conv; length=npts))
    P=periodogram(t,q,νgrid)
    # local maxima that stand `prominence` above the running median of their neighbourhood
    peaks=Tuple{Float64,Float64}[]
    win=max(5, npts ÷ 60)
    for i in 2:npts-1
        (P[i]>P[i-1] && P[i]≥P[i+1]) || continue
        lo=max(1,i-win); hi=min(npts,i+win); med=sort(P[lo:hi])[(hi-lo)÷2+1]
        P[i] > prominence*med && push!(peaks,(P[i],νgrid[i]/conv))
    end
    sort!(peaks, by=x->-x[1]); peaks=peaks[1:min(npeaks,length(peaks))]
    sort!(peaks, by=x->x[2])
    return [pk[2] for pk in peaks], νgrid./conv, P
end

end # module DGStarHP
