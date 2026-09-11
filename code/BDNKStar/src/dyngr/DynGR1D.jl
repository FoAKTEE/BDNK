#=
    DynGR1D — STAGE 2: a 1+1D spherically-symmetric DYNAMICAL-GR hydro engine.

    THE FIRST engine in this toolkit that DROPS the Cowling approximation and
    EVOLVES THE SPACETIME coupled to the matter. Every other time-domain engine
    (SphBDNK, FVCartesian, DGCart2/3D, cowling3d, DGStar) freezes the metric at
    the TOV background. Here the metric m(r,t), α(r,t) is RE-SOLVED from the
    current matter every substep — the genuine dynamical-GR coupling.

    FORMULATION (radial/areal gauge + polar slicing, constrained evolution):
        ds² = −α(r,t)² dt² + X(r,t)² dr² + r² dΩ²,   X² = 1/(1−2m/r).

    Matter: radial Valencia GRHydro (D, S, τ) in the orthonormal radial frame,
    REUSING FVCommon's cons2prim_barotrope / prim2cons_barotrope / HLL / MinMod.
    We densitize by √γ̃ = X r² (NOT the frozen e^{λ/2} r² — X is the LIVE metric).

    Constrained metric solve each substep (the hallmark of dynamical GR):
      • Hamiltonian constraint for the mass function:
            ∂m/∂r = 4π r² E,     E = ρ h W² − p   (Eulerian energy density),
        m(0)=0, integrated outward (trapezoid). Then X² = 1/(1−2m/r).
      • Polar-slicing lapse ODE:
            ∂(ln α)/∂r = X² [ m/r² + 4π r S_rr ],   S_rr = ρ h W² v² + p,
        integrated outward and matched to the Schwarzschild exterior at the
        outer boundary:  α(r_out) = √(1 − 2 m(r_out)/r_out).

    These two ODEs make the metric RESPOND to the matter motion. With no
    velocity the constraint reproduces the TOV m(r), α(r) (the static check).

    Conserved evolution (well-balanced finite-volume balance law):
        ∂_t(√γ̃ u) + (1/Δr)[ A_{k+1}F_{k+1} − A_k F_k ] = √γ̃ s,
    A = α X r² the face area weight, F the orthonormal radial fluxes. The
    geometric/gravity source uses the LIVE α, X, and ∂_r(ln α). We subtract the
    initial-equilibrium RHS (well-balancing) so the static TOV star starts at
    machine-zero RHS, then the metric solve keeps it balanced as it evolves.

    VALIDATION (each test a Cowling engine CANNOT do):
      1. static TOV stationarity over many dynamical times;
      2. full-GR radial mode F matching the Chandrasekhar GR eigenvalue and
         DIFFERING from the Cowling value;
      3. F²→0 at the maximum-mass turning point;
      4. collapse of a supercritical star: α_c→0, 2m/r→1, ρ_c↑.
=#
module DynGR1D

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2,
                         energy_from_pressure
using ..TOV
using ..TOV: TOVStar, solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, cons2prim_barotrope, prim2cons_barotrope,
                  minmod, hll_flux, lax_friedrichs_flux, rho_from_p, eos_cs2
using ..Numerics: brent
using ..RadialModes
using LinearAlgebra: eigen, Symmetric

export DynGRGrid, DynGRState, DynGREngine, setup_dyngr, evolve_dyngr!,
       seed_dyngr_velocity!, seed_dyngr_toroidal!, deplete_pressure!, dyngr_metric!,
       dyngr_radial_freq, chandrasekhar_radial_omega2, cowling_radial_omega2,
       dyngr_central_lapse, dyngr_max_2mor, dyngr_central_density

# ---------------------------------------------------------------------------
# Grid (cell-centered, NG ghosts). Identical layout to FVRadial.
# ---------------------------------------------------------------------------
struct DynGRGrid
    N::Int
    NG::Int
    Δr::Float64
    r::Vector{Float64}        # cell centers (incl ghosts)
    rL::Vector{Float64}       # interior faces 0..rmax (length N+1)
end

function DynGRGrid(N::Int, rmax::Float64; NG::Int=2)
    Δr = rmax / N
    ntot = N + 2NG
    r  = [ (i - NG - 0.5)*Δr for i in 1:ntot ]
    rL = [ (i-1)*Δr for i in 1:N+1 ]
    DynGRGrid(N, NG, Δr, r, rL)
end
@inline cidx(g::DynGRGrid, i::Int) = i + g.NG

# ---------------------------------------------------------------------------
# State: undensitized primitives + densitized conserved + LIVE metric arrays.
# ---------------------------------------------------------------------------
mutable struct DynGRState
    D::Vector{Float64}        # densitized conserved (by √γ̃ = X r²)
    S::Vector{Float64}
    τ::Vector{Float64}
    ρ::Vector{Float64}
    p::Vector{Float64}
    ε::Vector{Float64}
    v::Vector{Float64}        # orthonormal radial velocity
    # MHD: toroidal magnetic field B^φ̂ (orthonormal Eulerian component) — the only
    # div-B-free magnetic configuration in spherical symmetry (poloidal B breaks it).
    # b² = B²/W² (comoving), v·B=0 ⇒ b^t̂=0; magnetic energy/stress ½B²(1+v²) source
    # the metric.  B≡0 (default) ⇒ pure-fluid path is byte-identical to the original.
    B::Vector{Float64}        # Eulerian toroidal field B^φ̂  (orthonormal)
    B̃::Vector{Float64}        # densitized conserved magnetic var √γ̃·B (advected)
    # LIVE metric (cell centers)
    m::Vector{Float64}        # mass function m(r,t)
    X::Vector{Float64}        # X = 1/√(1−2m/r)
    α::Vector{Float64}        # lapse
    dlnα::Vector{Float64}     # ∂_r ln α (cell centers)
    sqrtg::Vector{Float64}    # √γ̃ = X r²
    # LIVE metric (faces, length N+1)
    Xf::Vector{Float64}
    αf::Vector{Float64}
    rf2::Vector{Float64}
end

function DynGRState(ntot::Int, Nf::Int)
    z() = zeros(ntot); zf() = zeros(Nf)
    DynGRState(z(),z(),z(), z(),z(),z(),z(), z(),z(),
               z(),z(),z(),z(),z(), zf(),zf(),zf())
end

mutable struct DynGREngine
    g::DynGRGrid
    eos::BarotropicEOS
    atm::AtmospherePars
    cfl::Float64
    riemann::Symbol
    R::Float64
    M::Float64
    εc::Float64
    freeze_metric::Bool       # if true, metric is NOT re-solved (Cowling control)
    mhd::Bool                 # if true, evolve the toroidal-field GRMHD sector
    r_resist::Float64         # BDNK resistivity r_b (0 ⇒ ideal; >0 ⇒ field diffuses/decays)
    ν_visc::Float64           # kinematic shear/bulk viscosity η̂=η/(ε+p) (0 ⇒ inviscid;
                              # >0 ⇒ radial momentum diffuses ⇒ stellar oscillations damp)
    # Velocity BC seen by the VISCOUS flux at the star/atmosphere interface.
    #   :zero    — atmosphere keeps v=0 (the cons2prim reset). Because v is PINNED just
    #              outside the star, the viscous operator feels a NO-SLIP RIGID WALL and
    #              develops a wall stress ⇒ γ∝√ν (the rigid-wall Stokes law). Default,
    #              i.e. the historical behaviour — keeps every existing test unchanged.
    #   :outflow — zero-gradient: an atmosphere cell inherits its interior neighbour's v,
    #              so no artificial wall stress is applied and the damping is bulk (γ∝ν).
    # This knob exists to TEST the 2026-07-19 audit's no-slip diagnosis; see
    # repro/viscous_damping.jl and VALIDATION.md §6.
    atm_vbc::Symbol
    # well-balancing equilibrium RHS (subtracted; static star → zero RHS)
    Seq_D::Vector{Float64}
    Seq_S::Vector{Float64}
    Seq_τ::Vector{Float64}
end

@inline _lininterp(xs, ys, x) = begin
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j])/(xs[j+1]-xs[j])
    ys[j] + t*(ys[j+1]-ys[j])
end

@inline function _p_of_rho(eos::BarotropicEOS, ρ::Float64)
    eos isa ShumPolytrope && return eos.κ*ρ^2
    return pressure(eos, ρ)
end

# ---------------------------------------------------------------------------
# Eulerian source densities from primitives: E (energy), S_rr (radial stress).
# E = ρhW² − p,  S_rr = ρhW² v² + p   (orthonormal radial frame).
# ---------------------------------------------------------------------------
@inline function _euler_densities(eos, ρ, p, v)
    ε = energy_from_pressure(eos, p)
    h = ρ > 0 ? (ε + p)/ρ : 1.0
    v2 = clamp(v^2, 0.0, 1.0-1e-12)
    W2 = 1.0/(1.0 - v2)
    ρhW2 = ρ*h*W2
    E   = ρhW2 - p
    Srr = ρhW2*v2 + p
    return E, Srr
end

# ---------------------------------------------------------------------------
# THE DYNAMICAL-GR METRIC SOLVE. Integrate the Hamiltonian constraint for m(r)
# and the polar-slicing ODE for α(r) from the CURRENT matter. Writes the live
# metric arrays (centers + faces). This is what makes the spacetime evolve.
# ---------------------------------------------------------------------------
function dyngr_metric!(st::DynGRState, eng::DynGREngine)
    eng.freeze_metric && return eng.M    # Cowling control: metric stays frozen
    return _solve_metric!(st, eng)
end

# unconditional constraint solve (used at setup even when freeze_metric=true)
function _solve_metric!(st::DynGRState, eng::DynGREngine)
    g = eng.g; eos = eng.eos; N = g.N; NG = g.NG; Δr = g.Δr
    rL = g.rL

    # --- 1. mass function m(r) from ∂m/∂r = 4π r² E (Hamiltonian constraint).
    # Evaluate E at cell centers; integrate to FACES by trapezoid (m at faces),
    # then m at centers = average of bracketing faces.
    # Eulerian energy density at interior cell centers.
    Ec = Vector{Float64}(undef, N)        # E at interior centers (i=1..N)
    @inbounds for i in 1:N
        ai = NG + i
        E, _ = _euler_densities(eos, st.ρ[ai], st.p[ai], st.v[ai])
        # toroidal-field magnetic energy ½B²(1+v²) (Eulerian; =0 if unmagnetized)
        Emag = 0.5*st.B[ai]^2*(1 + st.v[ai]^2)
        Ec[i] = E + Emag
    end
    # face integrand q(r) = 4π r² E. Use center E to estimate the face integrand
    # by extrapolation: at face k (radius rL[k]) take the average of the two
    # bracketing center values of (4π r² E). m at faces accumulated by trapezoid.
    qcenter = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        r = g.r[NG+i]
        qcenter[i] = 4π * r^2 * Ec[i]
    end
    mface = Vector{Float64}(undef, N+1)
    mface[1] = 0.0                          # m(0)=0
    @inbounds for k in 2:N+1
        # integral over cell (k-1): width Δr, integrand ≈ qcenter[k-1]
        mface[k] = mface[k-1] + Δr * qcenter[k-1]
    end
    # outside the matter, m stays at its surface value (atmosphere E≈0 already).

    # m at interior centers = midpoint of bracketing faces.
    @inbounds for i in 1:N
        ai = NG + i
        st.m[ai] = 0.5*(mface[i] + mface[i+1])
        r = g.r[ai]
        f = 1 - 2*st.m[ai]/r
        f = max(f, 1e-10)
        st.X[ai] = 1.0/sqrt(f)
        st.sqrtg[ai] = st.X[ai]*r^2
    end
    # faces
    @inbounds for k in 1:N+1
        r = rL[k]
        if r ≤ 0
            st.Xf[k] = 1.0; st.rf2[k] = 0.0
        else
            f = max(1 - 2*mface[k]/r, 1e-10)
            st.Xf[k] = 1.0/sqrt(f); st.rf2[k] = r^2
        end
    end
    Mtot = mface[N+1]
    eng.M = Mtot

    # --- 2. lapse α(r) from polar slicing: ∂(ln α)/∂r = X²[ m/r² + 4π r S_rr ].
    # Integrate the RHS outward to faces, then exp, then rescale so that
    # α(r_out) = √(1 − 2M/r_out)  (Schwarzschild exterior match).
    rhsc = Vector{Float64}(undef, N)        # ∂_r ln α at interior centers
    @inbounds for i in 1:N
        ai = NG + i; r = g.r[ai]
        _, Srr = _euler_densities(eos, st.ρ[ai], st.p[ai], st.v[ai])
        # toroidal magnetic radial stress ½B²(1+v²) (=0 if unmagnetized)
        Srrmag = 0.5*st.B[ai]^2*(1 + st.v[ai]^2)
        X2 = st.X[ai]^2
        rhsc[i] = X2*( st.m[ai]/r^2 + 4π*r*(Srr + Srrmag) )
    end
    # ln α accumulated at faces (unnormalized); start lnα_face[1]=0 at center.
    lnαface = Vector{Float64}(undef, N+1)
    lnαface[1] = 0.0
    @inbounds for k in 2:N+1
        lnαface[k] = lnαface[k-1] + Δr*rhsc[k-1]
    end
    # ln α at interior centers (midpoint of faces)
    lnαc = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        lnαc[i] = 0.5*(lnαface[i] + lnαface[i+1])
    end
    # match: α(r_out_face) = √(1−2M/r_out). Use the last interior face.
    rout = rL[N+1]
    αout_target = sqrt(max(1 - 2*Mtot/rout, 1e-10))
    shift = log(αout_target) - lnαface[N+1]
    @inbounds for i in 1:N
        ai = NG + i
        st.α[ai] = exp(lnαc[i] + shift)
        st.dlnα[ai] = rhsc[i]
    end
    @inbounds for k in 1:N+1
        st.αf[k] = exp(lnαface[k] + shift)
    end
    # ghost-cell metric (reflect inner, Schwarzschild outer)
    @inbounds for gc in 1:NG
        ii = NG - gc + 1; jj = NG + gc
        st.m[ii]=st.m[jj]; st.X[ii]=st.X[jj]; st.α[ii]=st.α[jj]
        st.dlnα[ii]=-st.dlnα[jj]; st.sqrtg[ii]=st.sqrtg[jj]
        io = NG + N + gc; jo = NG + N
        st.m[io]=Mtot; r=g.r[io]; f=max(1-2*Mtot/r,1e-10)
        st.X[io]=1/sqrt(f); st.sqrtg[io]=st.X[io]*r^2
        st.α[io]=sqrt(f); st.dlnα[io]=st.dlnα[jo]
    end
    return Mtot
end

# ---------------------------------------------------------------------------
# TOROIDAL-FIELD GRMHD prim<->cons (orthonormal radial frame).  For B=B^φ̂ with
# v·B=0: b^t̂=0, b²=B²/W²; E_f≡ρhW²=(ε+p)W²;
#   D = ρW,   S = (E_f + B²) v,   τ = E_f + B² − p − B²/(2W²) − D = τ_fluid + ½B²(1+v²).
# ---------------------------------------------------------------------------
@inline function prim2cons_toroidal(eos::BarotropicEOS, ρ, p, v, B)
    W  = 1.0/sqrt(1.0 - clamp(v^2,0.0,1.0-1e-12))
    ε  = energy_from_pressure(eos, p)
    Ef = (ε + p)*W^2
    B2 = B^2
    D  = ρ*W
    S  = (Ef + B2)*v
    τ  = Ef + B2 - p - B2/(2*W^2) - D
    return D, S, τ
end

# S_pred(v)=(E_f(v)+B²)v for the toroidal cons2prim root-find (module-level ⇒
# specialised per concrete EOS, NO closure boxing in the per-cell hot loop).
@inline function _S_of_v_tor(eos::BarotropicEOS, D̂::Float64, B2::Float64, v::Float64)
    W = 1.0/sqrt(1.0 - v^2)
    ρ = D̂/W
    p = _p_of_rho(eos, ρ)
    ε = energy_from_pressure(eos, p)
    return ((ε + p)*W^2 + B2)*v
end

# Recover (ρ,p,ε,|v|,W) from (D̂,Ŝ,τ̂,B̂).  S_pred(v)=(E_f(v)+B²)v is monotone in v
# (E_f∝D·W for a polytrope) ⇒ bisection on v∈[0,vmax].  Reduces to the fluid
# cons2prim when B=0.
function cons2prim_toroidal(eos::BarotropicEOS, D̂::Float64, Ŝ::Float64, τ̂::Float64,
                            B̂::Float64, atm::AtmospherePars)
    if D̂ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0, 1.0, true)
    end
    B2 = B̂^2; Sabs = abs(Ŝ)
    # bracket [0, vmax]; S_pred(0)=0 ≤ Sabs, S_pred(vmax) large
    vlo, vhi = 0.0, atm.vmax
    if _S_of_v_tor(eos, D̂, B2, vhi) < Sabs
        v = vhi                       # cap at vmax (mildly superluminal cons)
    else
        for _ in 1:64
            vm = 0.5*(vlo+vhi)
            (_S_of_v_tor(eos, D̂, B2, vm) < Sabs) ? (vlo = vm) : (vhi = vm)
        end
        v = 0.5*(vlo+vhi)
    end
    W = 1.0/sqrt(1.0 - v^2)
    ρ = D̂/W
    p = max(_p_of_rho(eos, ρ), atm.p_atm)
    ε = energy_from_pressure(eos, p)
    if ρ ≤ atm.ρ_cut
        return (atm.ρ_atm, atm.p_atm, atm.ε_atm, 0.0, 1.0, true)
    end
    return (ρ, p, max(ε,0.0), v, W, false)
end

# ---------------------------------------------------------------------------
# cons2prim over the grid (LIVE √γ̃). Re-densitizes after atmosphere resets.
# ---------------------------------------------------------------------------
function _update_primitives!(st::DynGRState, eng::DynGREngine)
    eos = eng.eos; atm = eng.atm; mhd = eng.mhd; g = eng.g
    @inbounds for i in eachindex(st.D)
        sg = st.sqrtg[i]
        if sg ≤ 0
            st.ρ[i]=atm.ρ_atm; st.p[i]=atm.p_atm; st.ε[i]=atm.ε_atm; st.v[i]=0.0
            mhd && (st.B[i]=0.0)
            continue
        end
        D̂ = st.D[i]/sg; Ŝ = st.S[i]/sg; τ̂ = st.τ[i]/sg
        if mhd
            # recover the Eulerian toroidal field from the conserved Φ_B = X r B
            Xr = st.X[i]*g.r[i]
            B = Xr > 0 ? st.B̃[i]/Xr : 0.0
            st.B[i] = B
            ρ,p,ε,vmag,W,isatm = cons2prim_toroidal(eos, D̂, abs(Ŝ), τ̂, B, atm)
            v = isatm ? 0.0 : sign(Ŝ)*vmag
            st.ρ[i]=ρ; st.p[i]=p; st.ε[i]=ε; st.v[i]=v
            if isatm
                D̂2,Ŝ2,τ̂2 = prim2cons_toroidal(eos, ρ, p, 0.0, B)
                st.D[i]=sg*D̂2; st.S[i]=sg*Ŝ2; st.τ[i]=sg*τ̂2
            end
        else
            ρ,p,ε,vmag,W,isatm = cons2prim_barotrope(eos, D̂, abs(Ŝ), τ̂, atm)
            v = isatm ? 0.0 : sign(Ŝ)*vmag
            st.ρ[i]=ρ; st.p[i]=p; st.ε[i]=ε; st.v[i]=v
            if isatm
                D̂2,Ŝ2,τ̂2 = prim2cons_barotrope(eos, ρ, p, 0.0)
                st.D[i]=sg*D̂2; st.S[i]=sg*Ŝ2; st.τ[i]=sg*τ̂2
            end
        end
    end
end

function _fill_ghosts!(st::DynGRState, eng::DynGREngine)
    g=eng.g; NG=g.NG; N=g.N; eos=eng.eos; atm=eng.atm
    for gc in 1:NG
        ii = NG - gc + 1; jj = NG + gc
        st.ρ[ii]=st.ρ[jj]; st.p[ii]=st.p[jj]; st.ε[ii]=st.ε[jj]; st.v[ii]=-st.v[jj]
        sg = st.sqrtg[ii]
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, st.ρ[ii], st.p[ii], abs(st.v[ii]))
        st.D[ii]=sg*D̂; st.S[ii]=sg*sign(st.v[ii])*Ŝ; st.τ[ii]=sg*τ̂
        if eng.mhd
            st.B[ii]=st.B[jj]; st.B̃[ii]=st.B̃[jj]   # toroidal field even at centre
        end
    end
    for gc in 1:NG
        ii = NG + N + gc
        st.ρ[ii]=atm.ρ_atm; st.p[ii]=atm.p_atm; st.ε[ii]=atm.ε_atm; st.v[ii]=0.0
        sg=st.sqrtg[ii]
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, atm.ρ_atm, atm.p_atm, 0.0)
        st.D[ii]=sg*D̂; st.S[ii]=sg*Ŝ; st.τ[ii]=sg*τ̂
        if eng.mhd
            st.B[ii]=0.0; st.B̃[ii]=0.0           # atmosphere unmagnetised
        end
    end
end

@inline function _reconstruct(q, aL::Int)
    sL = minmod(q[aL]-q[aL-1], q[aL+1]-q[aL])
    sR = minmod(q[aL+1]-q[aL], q[aL+2]-q[aL+1])
    return q[aL]+0.5*sL, q[aL+1]-0.5*sR
end

@inline function _phys_flux(eos, ρ, p, v)
    ε = energy_from_pressure(eos, p)
    h = ρ>0 ? (ε+p)/ρ : 1.0
    W = 1.0/sqrt(1.0-clamp(v^2,0.0,1.0-1e-12))
    D̂ = ρ*W; Ŝ = ρ*h*W^2*v; τ̂ = ρ*h*W^2 - p - ρ*W
    return (D̂,Ŝ,τ̂), (D̂*v, Ŝ*v+p, (τ̂+p)*v)
end

@inline function _wavespeeds(eos, ρ, p, ε, v)
    cs = sqrt(clamp(eos_cs2(eos, ε), 0.0, 1.0))
    v2 = clamp(v^2,0.0,1.0-1e-12)
    a = 1.0/(1.0 - v2*cs^2)
    disc = sqrt(max(cs^2*(1-v2)*(1 - v2*cs^2 - v^2*(1-cs^2)),0.0))
    return a*(v*(1-cs^2) - disc), a*(v*(1-cs^2) + disc)
end

# ---------------------------------------------------------------------------
# RAW RHS using the LIVE metric (st.α, st.X, st.dlnα). Same orthonormal-frame
# balance law as FVRadial but with α,X,Φ' taken from the current metric solve.
# ---------------------------------------------------------------------------
function _raw_rhs!(rD, rS, rτ, rB, st::DynGRState, eng::DynGREngine)
    g=eng.g; eos=eng.eos; N=g.N; NG=g.NG; Δr=g.Δr
    fill!(rD,0.0); fill!(rS,0.0); fill!(rτ,0.0); fill!(rB,0.0)
    mhd = eng.mhd
    @inbounds for k in 1:N+1
        aL = NG + k - 1
        ρL,ρR = _reconstruct(st.ρ, aL)
        pL,pR = _reconstruct(st.p, aL)
        vL,vR = _reconstruct(st.v, aL)
        ρL=max(ρL,eng.atm.ρ_atm); ρR=max(ρR,eng.atm.ρ_atm)
        pL=max(pL,eng.atm.p_atm); pR=max(pR,eng.atm.p_atm)
        vL=clamp(vL,-0.999,0.999); vR=clamp(vR,-0.999,0.999)
        εL=energy_from_pressure(eos,pL); εR=energy_from_pressure(eos,pR)
        UL,FL = _phys_flux(eos, ρL, pL, vL)
        UR,FR = _phys_flux(eos, ρR, pR, vR)
        λmL,λpL = _wavespeeds(eos, ρL, pL, εL, vL)
        λmR,λpR = _wavespeeds(eos, ρR, pR, εR, vR)
        if eng.riemann == :hll
            sLs = min(λmL, λmR, 0.0); sRs = max(λpL, λpR, 0.0)
            Fn = hll_flux(UL, FL, sLs, UR, FR, sRs)
        else
            amax = max(abs(λmL),abs(λpL),abs(λmR),abs(λpR))
            Fn = lax_friedrichs_flux(UL, FL, UR, FR, amax)
        end
        A = st.αf[k]*st.Xf[k]*st.rf2[k]      # LIVE face area weight α X r²
        FD=A*Fn[1]; FS=A*Fn[2]; Fτ=A*Fn[3]; FB=0.0
        if mhd
            BL,BR = _reconstruct(st.B, aL)
            # toroidal magnetic flux on the S,τ rows: U_mag=(B²v, ½B²(1+v²)),
            # F_mag=(½B²(1+v²), B²v), Lax–Friedrichs with the light-speed bound (1).
            UmagS_L=BL^2*vL;          UmagS_R=BR^2*vR
            UmagT_L=0.5*BL^2*(1+vL^2); UmagT_R=0.5*BR^2*(1+vR^2)
            FSmag = 0.5*(UmagT_L+UmagT_R) - 0.5*(UmagS_R-UmagS_L)   # F_S,mag=½B²(1+v²)
            FTmag = 0.5*(UmagS_L+UmagS_R) - 0.5*(UmagT_R-UmagT_L)   # F_τ,mag=B²v
            FS += A*FSmag; Fτ += A*FTmag
            # induction flux F=α v r B (upwind);  conserved Φ_B = X r B  (advected)
            rf = sqrt(st.rf2[k]); vf = 0.5*(vL+vR)
            Bup = vf ≥ 0 ? BL : BR
            FB = st.αf[k]*vf*rf*Bup
            if eng.r_resist > 0
                # BDNK resistive diffusion of the toroidal field: F_resist = −r⊥ ∂_l B,
                # r⊥ = r_b w/(w+b²) (general-EOS maps), proper radial gradient (1/X)∂_r B.
                pf=0.5*(pL+pR); εf=energy_from_pressure(eos,pf); wf=εf+pf
                Wf2=1.0/(1.0-clamp(vf^2,0.0,1.0-1e-12)); Bf=0.5*(BL+BR); b2f=Bf^2/Wf2
                rperp = eng.r_resist*wf/(wf+b2f)
                FB -= rperp*st.αf[k]*rf/max(st.Xf[k],1e-10)*(BR-BL)/Δr
            end
        end
        if eng.ν_visc > 0
            # Newtonian shear/bulk radial-momentum diffusion (parabolic). Viscous
            # radial stress τ_v = −μ (1/X)∂_r v, dynamic μ = ν_visc·w (w=ε+p ⇒ kinematic
            # diffusivity = ν_visc). Conserved S-flux = (αXr²)·τ_v = −α r² μ ∂_r v (the X
            # cancels via the proper gradient). Energy flux = v·τ_v (viscous work). Use the
            # RAW cell-center gradient (v[aL+1]−v[aL])/Δr — the minmod-reconstructed jump
            # vanishes for smooth modes. v≡0 equilibrium ⇒ zero flux ⇒ well-balancing intact.
            pf_v=0.5*(pL+pR); wf_v=energy_from_pressure(eos,pf_v)+pf_v
            μf=eng.ν_visc*wf_v
            vLc=st.v[aL]; vRc=st.v[aL+1]
            if eng.atm_vbc === :outflow
                # Zero-gradient at the star/atmosphere interface: an atmosphere cell
                # inherits its interior neighbour's velocity, so the viscous operator
                # applies NO wall stress there. (With the default :zero the atmosphere
                # v is pinned to 0 by cons2prim, which acts as a NO-SLIP RIGID WALL —
                # the artifact identified by the 2026-07-19 audit.)
                aL_atm = st.ρ[aL]   ≤ eng.atm.ρ_cut
                aR_atm = st.ρ[aL+1] ≤ eng.atm.ρ_cut
                if aL_atm && !aR_atm
                    vLc = vRc
                elseif aR_atm && !aL_atm
                    vRc = vLc
                elseif aL_atm && aR_atm
                    vLc = 0.0; vRc = 0.0
                end
            end
            dvdl=(vRc-vLc)/Δr
            FSv=-st.αf[k]*st.rf2[k]*μf*dvdl
            vmid=0.5*(vLc+vRc)
            FS += FSv; Fτ += vmid*FSv
        end
        if k ≥ 2
            aLc = NG + (k-1); rD[aLc]-=FD/Δr; rS[aLc]-=FS/Δr; rτ[aLc]-=Fτ/Δr; rB[aLc]-=FB/Δr
        end
        if k ≤ N
            aRc = NG + k; rD[aRc]+=FD/Δr; rS[aRc]+=FS/Δr; rτ[aRc]+=Fτ/Δr; rB[aRc]+=FB/Δr
        end
    end
    # geometric/gravity source with LIVE metric. Φ' ≡ ∂_r ln α (polar slicing).
    @inbounds for i in 1:N
        ai = NG + i; r = g.r[ai]
        α = st.α[ai]; X = st.X[ai]; Φp = st.dlnα[ai]
        ρ=st.ρ[ai]; p=st.p[ai]; ε=st.ε[ai]; v=st.v[ai]
        W2 = 1.0/(1.0-clamp(v^2,0.0,1.0-1e-12))
        sg = st.sqrtg[ai]
        rS[ai] += α*sg*( -(ε+p)*W2*Φp ) + 2.0*α*X*r*p
        rτ[ai] += α*sg*( -(ε+p)*W2*v*Φp )
        if mhd
            # toroidal magnetic gravity source (frame energy b²W²=B²; angular
            # magnetic stresses cancel ⇒ no 2αXr term).  Φ_B advection has no source.
            B2 = st.B[ai]^2
            rS[ai] += α*sg*( -B2*Φp )
            rτ[ai] += α*sg*( -B2*v*Φp )
        end
    end
    return nothing
end

function _rhs!(rD, rS, rτ, rB, st::DynGRState, eng::DynGREngine)
    _raw_rhs!(rD, rS, rτ, rB, st, eng)
    @inbounds for i in eachindex(rD)
        rD[i]-=eng.Seq_D[i]; rS[i]-=eng.Seq_S[i]; rτ[i]-=eng.Seq_τ[i]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Setup: TOV initial data, initialize LIVE metric from constraint solve, store
# the well-balancing equilibrium RHS.
# ---------------------------------------------------------------------------
"""
    setup_dyngr(eos, εc; N=400, rmax_fac=1.3, atm_fac=1e-7, cfl=0.3,
                riemann=:hll, h_tov=2e-4, wellbalanced=true) -> (engine, state)

Build the 1+1D dynamical-GR engine on the TOV star of central energy density
`εc`. The metric is initialized from TOV but then re-solved from the matter
every substep (dynamical GR). If `wellbalanced`, the static TOV star starts at
machine-zero RHS.
"""
function setup_dyngr(eos::BarotropicEOS, εc::Float64; N::Int=400,
                     rmax_fac::Float64=1.3, atm_fac::Float64=1e-7,
                     cfl::Float64=0.3, riemann::Symbol=:hll, h_tov::Float64=2e-4,
                     wellbalanced::Bool=true, freeze_metric::Bool=false,
                     ν_visc::Float64=0.0, atm_vbc::Symbol=:zero)
    star = solve_tov(eos, εc; h=h_tov)
    R, M = star.R, star.M
    rmax = rmax_fac*R
    g = DynGRGrid(N, rmax)
    ntot = N + 2g.NG; Nf = N+1
    rt=star.r; pt=star.p
    p_of(r) = r ≤ R ? max(_lininterp(rt, pt, r), 0.0) : 0.0

    ρc = rho_from_p(eos, pressure(eos, εc))
    ρ_atm = atm_fac*ρc
    p_atm = _p_of_rho(eos, ρ_atm)
    ε_atm = ρ_atm > 0 ? energy_from_pressure(eos, p_atm) : 0.0
    atm = AtmospherePars(ρ_atm, p_atm, ε_atm, 5*ρ_atm, 0.999)

    st = DynGRState(ntot, Nf)
    # init primitives (equilibrium, v=0)
    @inbounds for i in 1:ntot
        r = g.r[i]
        p = r ≤ R ? max(p_of(r), p_atm) : p_atm
        ρ = max(rho_from_p(eos, p), ρ_atm)
        ε = energy_from_pressure(eos, p)
        st.ρ[i]=ρ; st.p[i]=p; st.ε[i]=ε; st.v[i]=0.0
    end
    eng = DynGREngine(g, eos, atm, cfl, riemann, R, M, εc, freeze_metric, false, 0.0,
                      ν_visc, atm_vbc, zeros(ntot), zeros(ntot), zeros(ntot))
    # initial metric solve from the matter (reproduces TOV m,α)
    _solve_metric!(st, eng)
    # densitize conserved with the live √γ̃
    @inbounds for i in 1:ntot
        sg = st.sqrtg[i]
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, st.ρ[i], st.p[i], 0.0)
        st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂
    end
    _update_primitives!(st, eng); _fill_ghosts!(st, eng)
    _solve_metric!(st, eng)
    if wellbalanced
        _raw_rhs!(eng.Seq_D, eng.Seq_S, eng.Seq_τ, zeros(length(eng.Seq_D)), st, eng)
    end
    return eng, st
end

# ---------------------------------------------------------------------------
# Time stepping: SSP-RK2, with a METRIC RE-SOLVE at every substage (the metric
# is slaved to the matter — constrained evolution, no free metric evolution).
# ---------------------------------------------------------------------------
@inline function _max_speed(st::DynGRState, eng::DynGREngine)
    g=eng.g; eos=eng.eos; N=g.N; NG=g.NG; a=0.0
    @inbounds for i in 1:N
        ai=NG+i
        # coordinate characteristic speed ~ α/X × orthonormal speed
        λm,λp=_wavespeeds(eos, st.ρ[ai], st.p[ai], st.ε[ai], st.v[ai])
        fac = st.α[ai]/max(st.X[ai],1e-3)
        a=max(a, fac*max(abs(λm), abs(λp)))
    end
    return a
end

"""
    evolve_dyngr!(st, eng; tmax, sample_dt=-1, probe_frac=0.5, cfl=-1,
                  collapse=false) -> NamedTuple

Evolve to coordinate time `tmax`, re-solving the metric every substage. Records
time series: probe velocity, central density ρ_c, central lapse α_c, max(2m/r).
`drift` = max relative central-density change (static-stability diagnostic).
Returns (ts, probe, ρc, αc, max2mor, drift, collapsed).
"""
function evolve_dyngr!(st::DynGRState, eng::DynGREngine; tmax::Float64,
                       sample_dt::Float64=-1.0, probe_frac::Float64=0.5,
                       cfl::Float64=-1.0)
    g=eng.g; N=g.N; NG=g.NG; Δr=g.Δr
    cfl = cfl>0 ? cfl : eng.cfl
    ntot=N+2NG
    sample_dt = sample_dt>0 ? sample_dt : 5*Δr
    ai_probe = NG + clamp(round(Int, probe_frac*eng.R/Δr + 0.5), 1, N)
    ic = NG + 1
    ρc0 = st.ρ[ic]

    rD=zeros(ntot); rS=zeros(ntot); rτ=zeros(ntot); rB=zeros(ntot)
    Dn=zeros(ntot); Sn=zeros(ntot); τn=zeros(ntot); Bn=zeros(ntot)
    mhd = eng.mhd

    ts=Float64[]; probe=Float64[]; ρc_h=Float64[]; αc_h=Float64[]; m2_h=Float64[]
    t=0.0; last=-1e30; drift=0.0; collapsed=false
    nstep=0; maxsteps=5_000_000
    _update_primitives!(st, eng); _fill_ghosts!(st, eng); dyngr_metric!(st, eng)

    while t < tmax && nstep < maxsteps
        amax = max(_max_speed(st, eng), 1e-3)
        dt = min(cfl*Δr/amax, tmax-t)
        # explicit parabolic limit for the BDNK resistive diffusion (r⊥≈r_resist):
        # dt < 0.4·Δr²/(2 r⊥) keeps the field-diffusion stable (else IMEX needed).
        if eng.r_resist > 0
            dt = min(dt, 0.4*Δr^2/(2*eng.r_resist))
        end
        # explicit parabolic limit for the viscous momentum diffusion (kinematic ν_visc).
        if eng.ν_visc > 0
            dt = min(dt, 0.4*Δr^2/(2*eng.ν_visc))
        end

        @inbounds for i in 1:ntot
            Dn[i]=st.D[i]; Sn[i]=st.S[i]; τn[i]=st.τ[i]; mhd && (Bn[i]=st.B̃[i])
        end
        # stage 1
        _rhs!(rD,rS,rτ, rB, st, eng)
        @inbounds for i in 1:ntot
            st.D[i]=Dn[i]+dt*rD[i]; st.S[i]=Sn[i]+dt*rS[i]; st.τ[i]=τn[i]+dt*rτ[i]
            mhd && (st.B̃[i]=Bn[i]+dt*rB[i])
        end
        _update_primitives!(st, eng); _fill_ghosts!(st, eng); dyngr_metric!(st, eng)
        # stage 2
        _rhs!(rD,rS,rτ, rB, st, eng)
        @inbounds for i in 1:ntot
            st.D[i]=0.5*(Dn[i]+st.D[i]+dt*rD[i])
            st.S[i]=0.5*(Sn[i]+st.S[i]+dt*rS[i])
            st.τ[i]=0.5*(τn[i]+st.τ[i]+dt*rτ[i])
            mhd && (st.B̃[i]=0.5*(Bn[i]+st.B̃[i]+dt*rB[i]))
        end
        _update_primitives!(st, eng); _fill_ghosts!(st, eng); dyngr_metric!(st, eng)
        t += dt; nstep += 1

        # diagnostics
        max2mor = 0.0
        @inbounds for i in 1:N
            ai=NG+i; r=g.r[ai]
            v = 2*st.m[ai]/r
            v > max2mor && (max2mor = v)
        end
        if t - last ≥ sample_dt
            push!(ts,t); push!(probe,st.v[ai_probe]); push!(ρc_h,st.ρ[ic])
            push!(αc_h,st.α[ic]); push!(m2_h,max2mor)
            last=t; drift=max(drift, abs(st.ρ[ic]-ρc0)/ρc0)
        end
        # collapse / blow-up detection
        if !isfinite(st.ρ[ic]) || max2mor > 0.98 || st.α[ic] < 1e-3
            collapsed = true
            push!(ts,t); push!(probe,st.v[ai_probe]); push!(ρc_h,st.ρ[ic])
            push!(αc_h,st.α[ic]); push!(m2_h,max2mor)
            break
        end
    end
    return (ts=ts, probe=probe, ρc=ρc_h, αc=αc_h, max2mor=m2_h,
            drift=drift, collapsed=collapsed)
end

# ---------------------------------------------------------------------------
# Perturbation seeds.
# ---------------------------------------------------------------------------
"""
    seed_dyngr_velocity!(st, eng; A=1e-3, profile=:linear)

Seed a small ℓ=0 radial velocity (homologous `:linear` v=A r/R, or `:sin`).
Negative `A` is an INWARD kick (drives collapse of a marginally-stable star).
"""
function seed_dyngr_velocity!(st::DynGRState, eng::DynGREngine; A::Float64=1e-3,
                              profile::Symbol=:linear)
    g=eng.g; eos=eng.eos; N=g.N; NG=g.NG
    @inbounds for i in 1:N
        ai=NG+i; r=g.r[ai]
        if r < eng.R
            v = profile===:linear ? A*(r/eng.R) : A*sin(π*r/eng.R)
            st.v[ai]=v
            sg=st.sqrtg[ai]
            D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, st.ρ[ai], st.p[ai], abs(v))
            st.D[ai]=sg*D̂; st.S[ai]=sg*sign(v)*Ŝ; st.τ[ai]=sg*τ̂
        end
    end
    _fill_ghosts!(st, eng); dyngr_metric!(st, eng)
end

"""
    seed_dyngr_toroidal!(st, eng; B0=0.005, profile=:pressure) -> M

Enable the toroidal-field GRMHD sector (`eng.mhd=true`) and seed an Eulerian
toroidal field B^φ̂(r)=B0·prof(r) inside the star (`:pressure` ⇒ ∝p(r), `:sin` ⇒
sin(πr/R)).  Re-solves the metric WITH the magnetic energy/stress, then rebuilds
the conserved (D,S,τ) via the toroidal magnetic prim2cons and the densitized
induction variable Φ_B=X·r·B.  Call AFTER any velocity seed (uses the current v).
The induction is exactly conservative (∫Φ_B conserved); the field flux-freezes
(B∝ρrW) and amplifies under compression; the ADM mass is conserved under evolution.
"""
function seed_dyngr_toroidal!(st::DynGRState, eng::DynGREngine; B0::Float64=0.005,
                              profile::Symbol=:pressure, r_resist::Float64=0.0,
                              rebalance::Bool=false)
    eng.mhd = true; eng.r_resist = r_resist
    g=eng.g; N=g.N; NG=g.NG; R=eng.R; eos=eng.eos
    pc = maximum(@view st.p[NG+1:NG+N])
    @inbounds for i in eachindex(st.B)
        r = g.r[i]
        prof = (0 < r ≤ R) ? (profile===:pressure ? max(st.p[i],0.0)/pc : sin(π*r/R)) : 0.0
        st.B[i] = B0*max(prof, 0.0)
    end
    dyngr_metric!(st, eng)                         # re-solve WITH magnetic source
    @inbounds for i in eachindex(st.D)
        sg = st.sqrtg[i]; Xr = st.X[i]*g.r[i]
        D̂,Ŝ,τ̂ = prim2cons_toroidal(eos, st.ρ[i], st.p[i], st.v[i], st.B[i])
        st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂; st.B̃[i]=Xr*st.B[i]
    end
    _update_primitives!(st, eng); _fill_ghosts!(st, eng); dyngr_metric!(st, eng)
    if rebalance
        # well-balance the MAGNETISED state: subtract its initial RHS so the magnetised
        # star is (approximately) stationary — defines the equilibrium for clean MODE
        # perturbations (the f-mode then carries the magnetic-pressure restoring force).
        _raw_rhs!(eng.Seq_D, eng.Seq_S, eng.Seq_τ, zeros(length(eng.Seq_D)), st, eng)
    end
    return dyngr_metric!(st, eng)
end

"""
    deplete_pressure!(st, eng; factor=0.9)

Multiply the pressure (and consistently ρ, ε) by `factor<1` to trigger collapse
of an otherwise-stable star by pressure depletion. Re-densitizes & re-solves
the metric.
"""
function deplete_pressure!(st::DynGRState, eng::DynGREngine; factor::Float64=0.9)
    g=eng.g; eos=eng.eos; N=g.N; NG=g.NG; atm=eng.atm
    @inbounds for i in 1:N
        ai=NG+i; r=g.r[ai]
        if r < eng.R
            p = max(factor*st.p[ai], atm.p_atm)
            ρ = max(rho_from_p(eos,p), atm.ρ_atm)
            ε = energy_from_pressure(eos,p)
            st.p[ai]=p; st.ρ[ai]=ρ; st.ε[ai]=ε
            sg=st.sqrtg[ai]
            D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, ρ, p, abs(st.v[ai]))
            st.D[ai]=sg*D̂; st.S[ai]=sg*sign(st.v[ai])*Ŝ; st.τ[ai]=sg*τ̂
        end
    end
    _fill_ghosts!(st, eng); dyngr_metric!(st, eng)
end

# diagnostics
dyngr_central_lapse(st::DynGRState, eng::DynGREngine)   = st.α[eng.g.NG+1]
dyngr_central_density(st::DynGRState, eng::DynGREngine) = st.ρ[eng.g.NG+1]
function dyngr_max_2mor(st::DynGRState, eng::DynGREngine)
    g=eng.g; N=g.N; NG=g.NG; mx=0.0
    @inbounds for i in 1:N
        ai=NG+i; v=2*st.m[ai]/g.r[ai]; v>mx && (mx=v)
    end
    return mx
end

# ---------------------------------------------------------------------------
# Radial-mode frequency extraction (windowed periodogram of ρ_c(t)).
# ---------------------------------------------------------------------------
"""
    dyngr_radial_freq(ts, ρc; nf=4000) -> (f_geom, P, fgrid)

Dominant oscillation frequency (geometric units, km⁻¹) of the supplied time
series via a Hann-windowed periodogram. Convert to kHz with f[kHz] =
f[km⁻¹]/kHz_to_km (εc is in km⁻²; NO M⊙ rescaling). `fmin`,`fmax` are in km⁻¹.
"""
function dyngr_radial_freq(ts::Vector{Float64}, ρc::Vector{Float64}; nf::Int=4000,
                           fmin::Float64=0.0, fmax::Float64=-1.0)
    n=length(ts)
    n < 8 && return (NaN, Float64[], Float64[])
    T = ts[end]-ts[1]
    dtmin = minimum(diff(ts))
    fmax = fmax>0 ? fmax : 0.5/dtmin
    fmin = fmin>0 ? fmin : 1.0/T
    fgrid = collect(range(fmin, fmax; length=nf))
    q = copy(ρc); μ=sum(q)/n
    @inbounds for i in 1:n
        q[i]=(q[i]-μ)*(0.5-0.5*cos(2π*(i-1)/(n-1)))
    end
    P = Vector{Float64}(undef, nf)
    @inbounds for k in 1:nf
        ω=2π*fgrid[k]; sr=0.0; si=0.0
        for i in 1:n
            sr+=q[i]*cos(ω*ts[i]); si+=q[i]*sin(ω*ts[i])
        end
        P[k]=sr^2+si^2
    end
    # dominant interior peak
    fpk=fgrid[1]; Pmx=-1.0
    for k in 2:nf-1
        if P[k]>P[k-1] && P[k]≥P[k+1] && P[k]>Pmx
            Pmx=P[k]; fpk=fgrid[k]
        end
    end
    return (fpk, P, fgrid)
end

# ---------------------------------------------------------------------------
# REFERENCE radial-pulsation eigenvalues (full-GR Sturm–Liouville eigensolver).
#
# FULL-GR (Chandrasekhar 1964) relativistic radial pulsation, written as the
# self-adjoint Sturm–Liouville LAWE for the renormalized displacement
# ζ = r² e^{−ν/2} ξ  (ξ ≡ Δr/r):
#
#     d/dr( P dζ/dr ) + ( Q + ω² W ) ζ = 0 .
#
# We take the weight functions P, W, Q VERBATIM from Kokkotas & Ruoff (2001)
# A&A 366, 565 [gr-qc/0011093] eqs.(14)-(17), translated into our metric
# convention g_tt=−e^{ν} (Φ=ν/2, e^{λ}=1/(1−2m/r)); see the full convention
# derivation inside `chandrasekhar_radial_omega2`. In our convention:
#     P = Γ₁ p e^{(λ+3ν)/2} / r²
#     W = (ε+p) e^{(3λ+ν)/2} / r²
#     Q = e^{(λ+3ν)/2}(ε+p)/r² · [ ν'²/4 + 2ν'/r − 8π e^{λ} p ] .
# Discretized as a symmetric tridiagonal generalized eigenproblem A ζ = ω² B ζ
# (face-centered P, B = diag(W)) and solved with LAPACK generalized eigen.
#
# CORRECTNESS: the corrected exponents place the full-GR fundamental BELOW the
# Cowling value (ω²_GR ≈ 1.981e-3 < ω²_Cow ≈ 7.063e-3 km⁻², F_GR≈2.124 kHz) —
# the spacetime response softens the restoring force, matching the dynamical
# engine. Independently cross-checked by GHZ(1997) shooting in
# repro/radial_shoot_crosscheck.jl (F_shoot = 2.1236 kHz, 0.00% vs SL). The
# earlier 'GR above Cowling' bug was an exponent SWAP (P,Q ↔ W exponents).
# ---------------------------------------------------------------------------
"""
    chandrasekhar_radial_omega2(eos, εc; N=2000, h_tov=2e-4, nmodes=1) -> ω²[kHz-geom]

Fundamental full-GR radial-pulsation eigenvalue ω² (geometric, km⁻²) from the
Bardeen–Thorne–Meltzer (1966) / Chandrasekhar (1964) Sturm–Liouville radial
equation, discretized as a generalized matrix eigenproblem with the regular-
center and zero-Δp surface conditions. Returns the lowest `nmodes` ω² (ascending;
negative ⇒ radially unstable).
"""
function chandrasekhar_radial_omega2(eos::BarotropicEOS, εc::Float64; N::Int=1500,
                                     h_tov::Float64=2e-4, nmodes::Int=1)
    star = solve_tov(eos, εc; h=h_tov)
    R=star.R
    rt=star.r; mt=star.m; νt=star.ν; εt=star.ε; pt=star.p
    bg(r)=(max(_lininterp(rt,mt,r),0.0), _lininterp(rt,νt,r),
           max(_lininterp(rt,εt,r),1e-20), max(_lininterp(rt,pt,r),1e-30))
    # Self-adjoint LAWE in the SOURCE convention of Kokkotas & Ruoff (2001),
    # A&A 366, 565 [gr-qc/0011093], eqs. (14)-(17), with their metric
    #   ds² = −e^{2ν_KR} dt² + e^{2λ_KR} dr² + r² dΩ²   (their eq. 2)
    # and renormalized displacement  ζ = r² e^{−ν_KR} ξ  (their eq. 13):
    #   d/dr( P dζ/dr ) + ( Q + ω² W ) ζ = 0                         (KR eq. 14)
    #   r² W = (ε+p) e^{3λ_KR + ν_KR}                                (KR eq. 15)
    #   r² P = Γ₁ p e^{λ_KR + 3ν_KR}                                 (KR eq. 16)
    #   r² Q = e^{λ_KR + 3ν_KR}(ε+p)[ (ν_KR')² + (4/r)ν_KR'
    #                                  − 8π e^{2λ_KR} p ]            (KR eq. 17)
    #
    # CONVENTION TRANSLATION (the prior factor-2 / exponent-swap bug). OUR code
    # stores g_tt=−e^{ν} (Φ=ν/2) and e^{λ}=1/(1−2m/r), i.e.
    #   ν_KR = ν/2,   λ_KR = λ/2,   ν_KR' = ν'/2.
    # Hence the KR exponents become, in OUR raw (ν,λ):
    #   3λ_KR+ν_KR = (3λ+ν)/2   → W
    #   λ_KR+3ν_KR = (λ+3ν)/2   → P, Q
    #   e^{2λ_KR}  = e^{λ}       (= 1/(1−2m/r))
    #   (ν_KR')²   = ν'²/4 ,  (4/r)ν_KR' = 2ν'/r .
    # So in OUR convention (ζ = r² e^{−ν/2} ξ unchanged):
    #   P = Γ₁ p e^{(λ+3ν)/2} / r²
    #   W = (ε+p) e^{(3λ+ν)/2} / r²
    #   Q = e^{(λ+3ν)/2}(ε+p)/r² [ ν'²/4 + 2ν'/r − 8π e^{λ} p ].
    # (The previous version applied (3λ+ν)/2 to P,Q and (λ+3ν)/2 to W — the two
    #  exponents were SWAPPED, which over-stiffened ω² and inverted the GR/Cowling
    #  ordering. Dropping the −8π e^{λ}p term in Q recovers the Cowling restoring
    #  force; with it, ω²_GR sits below ω²_Cowling, as the dynamical engine shows.)
    n=N
    r=collect(range(R/n, R*(1-1e-6); length=n))
    dr=r[2]-r[1]
    P=zeros(n); Wt=zeros(n); Q=zeros(n)
    for i in 1:n
        m,ν,ε,p=bg(r[i])
        eλ = 1/max(1-2m/r[i],1e-12)
        λ  = log(eλ)
        cs2=clamp(sound_speed2(eos,ε),1e-12,1.0)
        Γ1 = cs2*(ε+p)/p
        νp = 2*(m+4π*r[i]^3*p)/(r[i]*(r[i]-2m))      # dν/dr (our raw ν)
        ePQ = exp((λ+3ν)/2)                          # P,Q exponent (λ+3ν)/2
        eW  = exp((3λ+ν)/2)                          # W   exponent (3λ+ν)/2
        P[i]  = Γ1*p*ePQ/r[i]^2
        Wt[i] = (ε+p)*eW/r[i]^2
        Q[i]  = ePQ*(ε+p)/r[i]^2 * ( νp^2/4 + 2νp/r[i] - 8π*eλ*p )
    end
    # assemble tridiagonal generalized eigenproblem A ζ = ω² B ζ,
    # A = −d/dr(P d/dr) − Q  (FD), B = diag(Wt). Dirichlet-like: ζ(0)=0 (regular
    # center: ξ finite ⇒ ζ∝r³→0), Δp(R)=0 enforced as natural (free) surface via
    # one-sided P→0. Use a simple symmetric FD for −(Pζ')'.
    A=zeros(n,n); B=zeros(n,n)
    Pf=zeros(n+1)   # P at faces i+1/2
    for i in 1:n-1
        Pf[i+1]=0.5*(P[i]+P[i+1])
    end
    Pf[1]=P[1]; Pf[n+1]=0.0   # surface face: P→0 (Δp=0 natural BC)
    for i in 1:n
        aw = Pf[i]/dr^2
        ae = Pf[i+1]/dr^2
        A[i,i] = aw+ae - Q[i]
        i>1 && (A[i,i-1] = -aw)
        i<n && (A[i,i+1] = -ae)
        B[i,i] = Wt[i]
    end
    # generalized eigenvalues
    F = eigen(Symmetric(A), Symmetric(B))
    ω2 = sort(real.(F.values))
    k=min(nmodes,length(ω2))
    return ω2[1:k]
end

"""
    cowling_radial_omega2(eos, εc; N=1500, h_tov=5e-5, nmodes=1) -> ω²[geom]

Relativistic COWLING (frozen-metric) radial-pulsation eigenvalue ω² (geometric,
km⁻²) — the reference for the FROZEN-metric mode. This delegates to the package's
validated Stage-1A eigensolver `radial_cowling_spectrum` (verbatim NSO operator),
so it is the trusted Cowling ground truth. The full-GR `chandrasekhar_radial_omega2`
sits slightly BELOW this (the spacetime response softens the restoring force),
which is the radial-mode signature distinguishing dynamical GR from Cowling.
"""
function cowling_radial_omega2(eos::BarotropicEOS, εc::Float64; N::Int=1500,
                               h_tov::Float64=5e-5, nmodes::Int=1)
    _, ω2, _ = RadialModes.radial_cowling_spectrum(eos, εc; N=N, h_tov=h_tov,
                                                   nmodes=nmodes)
    return ω2
end

end # module DynGR1D
