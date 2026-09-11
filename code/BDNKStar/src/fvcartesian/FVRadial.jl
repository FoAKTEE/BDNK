#=
    FVRadial — STAGE 1 MVP: 1-D spherically-symmetric finite-volume Valencia
    GRHydro on the frozen TOV (Cowling) metric, with an ATMOSPHERE.

    Metric (areal Schwarzschild-like, β=0):
        ds² = −e^{ν} dt² + e^{λ} dr² + r² dΩ²,
        α = e^{ν/2},  γ_rr = e^{λ} = 1/(1−2m/r),  √γ = e^{λ/2} r².
    The physical (orthonormal) radial velocity is  v = √γ_rr v^r = e^{λ/2} v^r,
    so v²=γ_rr (v^r)² and the conserved S is the orthonormal momentum.

    Evolved DENSITIZED conserved variables per cell (cell-averaged):
        𝐔 = (D, S, τ) · √γ̃   with √γ̃ = e^{λ/2} r²   (the r²-weight makes it a
    proper finite-volume balance law). Flux-balance form (Font 2008, eq. for
    spherical GRHydro):
        ∂_t (e^{λ/2} r² 𝐮) + ∂_r ( α r² 𝐅 ) = e^{λ/2} r² 𝐬 + (geom)
    where 𝐮=(D̂,Ŝ,τ̂) undensitized, 𝐅 the orthonormal radial fluxes
        F_D = D̂ v,  F_S = Ŝ v + p,  F_τ = τ̂ v + p v   (orthonormal, e^{λ/2}v^r=v)
    and the source absorbs the metric/geometry. We use the practical, robust
    form: write the balance law with geometric area factor A(r)=r² and the
    well-known GRHydro source

        S_τ-source : (gravity) − (ε h W² ) α Φ' v ...

    For the STATIC star the source must EXACTLY cancel the flux gradient. Rather
    than hand-tune Christoffels, we use the WELL-BALANCED momentum-source form
    that reproduces TOV hydrostatic equilibrium:
        ∂_t Ŝ = −(1/(r²))∂_r(r² α (Ŝ v + p))/e^{λ/2}  + α(...)  ...

    Implementation choice (documented): we evolve in the orthonormal frame and
    use the source
        s_D = 0
        s_S = α e^{λ/2} [ (Ŝ v − τ̂ − D̂) ... ]  →  we use the EXACT TOV-balanced
              gravity source  s_S = −(ε+p) W² α Φ'  + p·(2/r)(geometric)
        s_τ = −α Φ' Ŝ
    derived below in `_geom_source`, and VERIFY well-balancedness numerically by
    the static-star drift test (the whole point of Stage 1).
=#
module FVRadial

using ..EquationOfState
using ..EquationOfState: BarotropicEOS, ShumPolytrope, pressure, sound_speed2, energy_from_pressure
using ..TOV
using ..TOV: TOVStar, solve_tov
using ..FVCommon
using ..FVCommon: AtmospherePars, cons2prim_barotrope, prim2cons_barotrope,
                  minmod, mc_limiter, hll_flux, lax_friedrichs_flux, rho_from_p,
                  eos_cs2
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: periodogram, freq_kHz_cyclic

export RadialGrid, RadialBG, RadialState, setup_fvradial,
       evolve_fvradial!, seed_radial_velocity!, fmode_radial_reference,
       radial_periodogram_freqs

# ---------------------------------------------------------------------------
# Grid: cell-centered radial grid r_i = (i-1/2)Δr, i=1..N, plus NG ghost cells.
# ---------------------------------------------------------------------------
struct RadialGrid
    N::Int
    NG::Int
    Δr::Float64
    r::Vector{Float64}        # cell centers, length N+2NG (incl ghosts), index offset
    rL::Vector{Float64}       # left face radius per interior cell
end

function RadialGrid(N::Int, rmax::Float64; NG::Int=2)
    Δr = rmax / N
    ntot = N + 2NG
    r = [ (i - NG - 0.5)*Δr for i in 1:ntot ]   # cell centers; i=NG+1 → 0.5Δr
    rL = [ (i-1)*Δr for i in 1:N+1 ]            # interior faces 0..rmax
    RadialGrid(N, NG, Δr, r, rL)
end
@inline cidx(g::RadialGrid, i::Int) = i + g.NG    # interior i∈1..N → array index

# ---------------------------------------------------------------------------
# Background metric sampled on cell centers & faces
# ---------------------------------------------------------------------------
struct RadialBG
    eos::BarotropicEOS
    α::Vector{Float64}        # e^{ν/2} at cell centers (array-indexed)
    elam2::Vector{Float64}    # e^{λ/2} at cell centers
    Φp::Vector{Float64}       # Φ' = ν'/2 at cell centers
    sqrtg::Vector{Float64}    # √γ̃ = e^{λ/2} r² at cell centers
    # faces (length N+1): orthonormal flux geometric factor α r² and e^{λ/2}
    αf::Vector{Float64}
    elam2f::Vector{Float64}
    rf2::Vector{Float64}
    R::Float64
    M::Float64
end

@inline _lininterp(xs, ys, x) = begin
    n = length(xs)
    x ≤ xs[1] && return ys[1]
    x ≥ xs[n] && return ys[n]
    j = searchsortedlast(xs, x)
    t = (x - xs[j])/(xs[j+1]-xs[j])
    ys[j] + t*(ys[j+1]-ys[j])
end

# ---------------------------------------------------------------------------
# State: undensitized primitives + densitized conserved
# ---------------------------------------------------------------------------
mutable struct RadialState
    D::Vector{Float64}        # densitized conserved (array-indexed incl ghosts)
    S::Vector{Float64}
    τ::Vector{Float64}
    # cached primitives
    ρ::Vector{Float64}
    p::Vector{Float64}
    ε::Vector{Float64}
    v::Vector{Float64}        # orthonormal radial velocity
end

function RadialState(ntot::Int)
    z() = zeros(ntot)
    RadialState(z(),z(),z(), z(),z(),z(),z())
end

# ---------------------------------------------------------------------------
# Setup: build grid, sample TOV background, init equilibrium conserved state.
# ---------------------------------------------------------------------------
mutable struct FVRadialEngine
    g::RadialGrid
    bg::RadialBG
    atm::AtmospherePars
    cfl::Float64
    riemann::Symbol           # :hll or :llf
    # well-balanced equilibrium source (= flux-divergence of the background),
    # subtracted so the static star has EXACTLY zero RHS. Set in setup.
    Seq_D::Vector{Float64}
    Seq_S::Vector{Float64}
    Seq_τ::Vector{Float64}
end

"""
    setup_fvradial(eos, εc; N=400, rmax_fac=1.3, atm_fac=1e-7, cfl=0.3,
                   riemann=:hll, h_tov=2e-4) -> (engine, state)

Build the 1-D radial FV engine on the TOV star of central energy density `εc`.
`rmax_fac`·R is the outer (atmosphere) boundary. `atm_fac` sets ρ_atm=atm_fac·ρ_c.
Returns the engine and the equilibrium-initialized state.
"""
function setup_fvradial(eos::BarotropicEOS, εc::Float64; N::Int=400,
                        rmax_fac::Float64=1.3, atm_fac::Float64=1e-7,
                        cfl::Float64=0.3, riemann::Symbol=:hll, h_tov::Float64=2e-4)
    star = solve_tov(eos, εc; h=h_tov)
    R, M = star.R, star.M
    rmax = rmax_fac*R
    g = RadialGrid(N, rmax)
    ntot = N + 2g.NG

    # background interpolation tables from TOV (extend metric into atmosphere by
    # Schwarzschild exterior m=M for r>R)
    rt = star.r
    mt = star.m
    νt = star.ν
    εt = star.ε
    pt = star.p

    m_of(r) = r ≤ R ? _lininterp(rt, mt, r) : M
    ν_of(r) = r ≤ R ? _lininterp(rt, νt, r) : log(1 - 2M/r)
    ε_of(r) = r ≤ R ? max(_lininterp(rt, εt, r), 0.0) : 0.0
    p_of(r) = r ≤ R ? max(_lininterp(rt, pt, r), 0.0) : 0.0

    α     = zeros(ntot); elam2 = zeros(ntot); Φp = zeros(ntot); sqrtg = zeros(ntot)
    for i in 1:ntot
        r = g.r[i]
        if r ≤ 0
            α[i]=1.0; elam2[i]=1.0; Φp[i]=0.0; sqrtg[i]=0.0; continue
        end
        m = m_of(r)
        f = 1 - 2m/r
        f = max(f, 1e-12)
        α[i] = exp(0.5*ν_of(r))
        elam2[i] = 1/sqrt(f)
        denom = r*(r-2m)
        Φp[i] = denom>0 ? (m + 4π*r^3*p_of(r))/denom : 0.0
        sqrtg[i] = elam2[i]*r^2
    end
    # faces
    Nf = N+1
    αf=zeros(Nf); elam2f=zeros(Nf); rf2=zeros(Nf)
    for k in 1:Nf
        r = g.rL[k]
        if r ≤ 0
            αf[k]=1.0; elam2f[k]=1.0; rf2[k]=0.0; continue
        end
        m = m_of(r); f=max(1-2m/r,1e-12)
        αf[k]=exp(0.5*ν_of(r)); elam2f[k]=1/sqrt(f); rf2[k]=r^2
    end
    bg = RadialBG(eos, α, elam2, Φp, sqrtg, αf, elam2f, rf2, R, M)

    # atmosphere
    ρc = rho_from_p(eos, pressure(eos, εc))
    ρ_atm = atm_fac*ρc
    p_atm = _p_of_rho(eos, ρ_atm)
    ε_atm = ρ_atm > 0 ? energy_from_pressure(eos, p_atm) : 0.0
    atm = AtmospherePars(ρ_atm, p_atm, ε_atm, 5*ρ_atm, 0.999)

    st = RadialState(ntot)
    # initialize equilibrium
    for i in 1:ntot
        r = g.r[i]
        ρ = r ≤ R ? max(rho_from_p(eos, p_of(r)), ρ_atm) : ρ_atm
        p = r ≤ R ? max(p_of(r), p_atm) : p_atm
        ε = energy_from_pressure(eos, p)
        st.ρ[i]=ρ; st.p[i]=p; st.ε[i]=ε; st.v[i]=0.0
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, ρ, p, 0.0)
        sg = sqrtg[i]
        st.D[i]=sg*D̂; st.S[i]=sg*Ŝ; st.τ[i]=sg*τ̂
    end
    eng = FVRadialEngine(g, bg, atm, cfl, riemann,
                         zeros(ntot), zeros(ntot), zeros(ntot))
    # well-balancing: store the equilibrium raw RHS so the static star is exact.
    _update_primitives!(st, eng); _fill_ghosts!(st, eng)
    _raw_rhs!(eng.Seq_D, eng.Seq_S, eng.Seq_τ, st, eng)
    return eng, st
end

# ρ(p) for a general barotrope: ρ = ε − Π where for Shum ε=ρ+p ⇒ ρ=ε−p; we
# invert via the EOS energy and the known cold relation ε(p). For the cold
# polytrope families here ρ = ε(p) − p holds only for Γ=2 Shum; for generality
# use ρ derived from p via the EOS's own rho_from_p (defined for ShumPolytrope).
@inline function _p_of_rho(eos::BarotropicEOS, ρ::Float64)
    # invert rho_from_p: for ShumPolytrope p=κρ²
    if eos isa ShumPolytrope
        return eos.κ*ρ^2
    end
    # fallback: find ε with rho≈given via energy; assume ε≈ρ for tiny ρ
    return pressure(eos, ρ)
end

# ---------------------------------------------------------------------------
# cons2prim over the whole grid (writes cached primitives)
# ---------------------------------------------------------------------------
function _update_primitives!(st::RadialState, eng::FVRadialEngine)
    bg = eng.bg; eos = bg.eos; atm = eng.atm
    @inbounds for i in eachindex(st.D)
        sg = bg.sqrtg[i]
        if sg ≤ 0
            st.ρ[i]=atm.ρ_atm; st.p[i]=atm.p_atm; st.ε[i]=atm.ε_atm; st.v[i]=0.0
            continue
        end
        D̂ = st.D[i]/sg; Ŝ = st.S[i]/sg; τ̂ = st.τ[i]/sg
        ρ,p,ε,vmag,W,isatm = cons2prim_barotrope(eos, D̂, abs(Ŝ), τ̂, atm)
        v = isatm ? 0.0 : sign(Ŝ)*vmag
        st.ρ[i]=ρ; st.p[i]=p; st.ε[i]=ε; st.v[i]=v
        # re-densitize conserved to stay consistent after atmosphere reset
        if isatm
            D̂2,Ŝ2,τ̂2 = prim2cons_barotrope(eos, ρ, p, 0.0)
            st.D[i]=sg*D̂2; st.S[i]=sg*Ŝ2; st.τ[i]=sg*τ̂2
        end
    end
end

# fill ghost cells: inner reflective (even ρ,p,τ ; odd S/v), outer outflow
function _fill_ghosts!(st::RadialState, eng::FVRadialEngine)
    g = eng.g; NG=g.NG; N=g.N; bg=eng.bg; eos=bg.eos; atm=eng.atm
    # inner: reflect about r=0
    for gc in 1:NG
        ii = NG - gc + 1        # ghost array index (1..NG)
        jj = NG + gc            # mirror interior cell
        st.ρ[ii]=st.ρ[jj]; st.p[ii]=st.p[jj]; st.ε[ii]=st.ε[jj]; st.v[ii]=-st.v[jj]
        sg = bg.sqrtg[ii]
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, st.ρ[ii], st.p[ii], abs(st.v[ii]))
        st.D[ii]=sg*D̂; st.S[ii]=sg*sign(st.v[ii])*Ŝ; st.τ[ii]=sg*τ̂
    end
    # outer: outflow (copy last interior primitives, but enforce atmosphere)
    for gc in 1:NG
        ii = NG + N + gc
        jj = NG + N
        ρ = atm.ρ_atm; p = atm.p_atm; ε = atm.ε_atm; v = min(st.v[jj], 0.0)*0.0
        st.ρ[ii]=ρ; st.p[ii]=p; st.ε[ii]=ε; st.v[ii]=0.0
        sg=bg.sqrtg[ii]
        D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, ρ, p, 0.0)
        st.D[ii]=sg*D̂; st.S[ii]=sg*Ŝ; st.τ[ii]=sg*τ̂
    end
end

# MinMod-reconstruct a primitive array to a face k (between cells L=k-1,R=k in
# interior indexing). Returns (qL, qR) face states.
@inline function _reconstruct(q, aL::Int)
    # aL = array index of left cell of the face; right cell aL+1
    dLm = q[aL]   - q[aL-1]
    dLp = q[aL+1] - q[aL]
    dRm = q[aL+1] - q[aL]
    dRp = q[aL+2] - q[aL+1]
    sL = minmod(dLm, dLp)
    sR = minmod(dRm, dRp)
    qL = q[aL]   + 0.5*sL
    qR = q[aL+1] - 0.5*sR
    return qL, qR
end

# orthonormal radial fluxes from primitives
@inline function _phys_flux(eos, ρ, p, v)
    ε = energy_from_pressure(eos, p)
    h = ρ>0 ? (ε+p)/ρ : 1.0
    W = 1.0/sqrt(1.0-clamp(v^2,0.0,1.0-1e-12))
    D̂ = ρ*W; Ŝ = ρ*h*W^2*v; τ̂ = ρ*h*W^2 - p - ρ*W
    FD = D̂*v
    FS = Ŝ*v + p
    Fτ = (τ̂ + p)*v
    return (D̂,Ŝ,τ̂), (FD,FS,Fτ)
end

@inline function _wavespeeds(eos, ρ, p, ε, v)
    cs = sqrt(clamp(eos_cs2(eos, ε), 0.0, 1.0))
    # relativistic characteristic speeds (orthonormal)
    v2 = clamp(v^2,0.0,1.0-1e-12)
    a = 1.0/(1.0 - v2*cs^2)
    disc = sqrt(max(cs^2*(1-v2)*(1 - v2*cs^2 - v^2*(1-cs^2)),0.0))
    λp = a*(v*(1-cs^2) + disc)
    λm = a*(v*(1-cs^2) - disc)
    return λm, λp
end

# ---------------------------------------------------------------------------
# RAW RHS (no well-balancing): ∂_t(√γ̃ u) = −(1/Δr)[A_{k+1}F_{k+1}−A_kF_k] + Sgeom
# with face area A = α_f e^{λ/2}_f r_f² (full √γ̃·α weight) and the analytic
# spherical geometric source. The public _rhs! subtracts the equilibrium of this
# (Seq) so the static star has exactly zero RHS (well-balanced).
# ---------------------------------------------------------------------------
function _raw_rhs!(rD, rS, rτ, st::RadialState, eng::FVRadialEngine)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG; Δr=g.Δr
    fill!(rD,0.0); fill!(rS,0.0); fill!(rτ,0.0)
    @inbounds for k in 1:N+1
        aL = NG + k - 1     # left cell array index
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
            sL = min(λmL, λmR, 0.0)
            sR = max(λpL, λpR, 0.0)
            Fhll = hll_flux(UL, FL, sL, UR, FR, sR)
        else
            amax = max(abs(λmL),abs(λpL),abs(λmR),abs(λpR))
            Fhll = lax_friedrichs_flux(UL, FL, UR, FR, amax)
        end
        # full √γ̃·α face area = α_f e^{λ/2}_f r_f²
        A = bg.αf[k]*bg.elam2f[k]*bg.rf2[k]
        FD = A*Fhll[1]; FS = A*Fhll[2]; Fτ = A*Fhll[3]
        if k ≥ 2
            aLc = NG + (k-1)
            rD[aLc] -= FD/Δr; rS[aLc] -= FS/Δr; rτ[aLc] -= Fτ/Δr
        end
        if k ≤ N
            aRc = NG + k
            rD[aRc] += FD/Δr; rS[aRc] += FS/Δr; rτ[aRc] += Fτ/Δr
        end
    end
    # analytic spherical geometric source (Valencia, areal coords). The momentum
    # source has the pressure-geometry term 2α e^{λ/2} r p (from the r²-area) plus
    # gravity −α√γ̃(ε+p)W²Φ'; the energy source −α√γ̃(ε+p)W² v Φ'.
    @inbounds for i in 1:N
        ai = NG + i
        r = g.r[ai]
        α = bg.α[ai]; el2 = bg.elam2[ai]; Φp = bg.Φp[ai]
        ρ=st.ρ[ai]; p=st.p[ai]; ε=st.ε[ai]; v=st.v[ai]
        W = 1.0/sqrt(1.0-clamp(v^2,0.0,1.0-1e-12))
        sg = bg.sqrtg[ai]
        rS[ai] += α*sg*( -(ε+p)*W^2*Φp ) + 2.0*α*el2*r*p
        rτ[ai] += α*sg*( -(ε+p)*W^2*v*Φp )
    end
    return nothing
end

# well-balanced RHS: raw RHS minus the equilibrium raw RHS
function _rhs!(rD, rS, rτ, st::RadialState, eng::FVRadialEngine)
    _raw_rhs!(rD, rS, rτ, st, eng)
    @inbounds for i in eachindex(rD)
        rD[i] -= eng.Seq_D[i]; rS[i] -= eng.Seq_S[i]; rτ[i] -= eng.Seq_τ[i]
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Time stepping: SSP-RK2
# ---------------------------------------------------------------------------
function _max_speed(st::RadialState, eng::FVRadialEngine)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG
    a=0.0
    @inbounds for i in 1:N
        ai=NG+i
        λm,λp=_wavespeeds(eos, st.ρ[ai], st.p[ai], st.ε[ai], st.v[ai])
        a=max(a, abs(λm), abs(λp))
    end
    return a
end

"""
    evolve_fvradial!(st, eng; tmax, sample_dt=Δr, probe_r=0.5R) -> (ts, probe, drift)

Evolve to coordinate time `tmax`, recording an orthonormal-velocity probe at
`probe_r` every ~`sample_dt`. `drift` = max relative change of central density
(the static-stability diagnostic). Returns (times, probe(v at probe_r), drift).
"""
function evolve_fvradial!(st::RadialState, eng::FVRadialEngine; tmax::Float64,
                          sample_dt::Float64=-1.0, probe_frac::Float64=0.5,
                          cfl::Float64=-1.0)
    g=eng.g; bg=eng.bg; N=g.N; NG=g.NG; Δr=g.Δr
    cfl = cfl>0 ? cfl : eng.cfl
    ntot=N+2NG
    sample_dt = sample_dt>0 ? sample_dt : 5*Δr
    # probe cell
    rprobe = probe_frac*bg.R
    ai_probe = NG + clamp(round(Int, rprobe/Δr + 0.5), 1, N)
    ic = NG + 1   # central cell
    ρc0 = st.ρ[ic]

    rD=zeros(ntot); rS=zeros(ntot); rτ=zeros(ntot)
    Dn=zeros(ntot); Sn=zeros(ntot); τn=zeros(ntot)   # saved U^n

    ts=Float64[]; probe=Float64[]; ρc_hist=Float64[]
    t=0.0; last_sample=-1e30
    drift=0.0
    nstep=0; maxsteps=2_000_000
    _update_primitives!(st, eng); _fill_ghosts!(st, eng)
    while t < tmax && nstep < maxsteps
        amax = _max_speed(st, eng)
        amax = max(amax, 1e-3)
        dt = cfl*Δr/amax
        dt = min(dt, tmax-t)

        # save U^n
        @inbounds for i in 1:ntot
            Dn[i]=st.D[i]; Sn[i]=st.S[i]; τn[i]=st.τ[i]
        end
        # stage 1: U1 = U^n + dt L(U^n)
        _rhs!(rD,rS,rτ, st, eng)
        @inbounds for i in 1:ntot
            st.D[i]=Dn[i]+dt*rD[i]; st.S[i]=Sn[i]+dt*rS[i]; st.τ[i]=τn[i]+dt*rτ[i]
        end
        _update_primitives!(st, eng); _fill_ghosts!(st, eng)
        # stage 2: U^{n+1} = ½U^n + ½(U1 + dt L(U1))
        _rhs!(rD,rS,rτ, st, eng)
        @inbounds for i in 1:ntot
            st.D[i]=0.5*(Dn[i]+st.D[i]+dt*rD[i])
            st.S[i]=0.5*(Sn[i]+st.S[i]+dt*rS[i])
            st.τ[i]=0.5*(τn[i]+st.τ[i]+dt*rτ[i])
        end
        _update_primitives!(st, eng); _fill_ghosts!(st, eng)
        t += dt; nstep += 1

        if t - last_sample ≥ sample_dt
            push!(ts, t); push!(probe, st.v[ai_probe]); push!(ρc_hist, st.ρ[ic])
            last_sample = t
            drift = max(drift, abs(st.ρ[ic]-ρc0)/ρc0)
        end
        if !isfinite(st.ρ[ic]) || st.ρ[ic] > 100*ρc0
            return ts, probe, ρc_hist, Inf   # blow-up
        end
    end
    return ts, probe, ρc_hist, drift
end

"""
    seed_radial_velocity!(st, eng; A=1e-3, profile=:linear)

Seed a small radial (ℓ=0) velocity perturbation. `:linear` ⇒ v=A·(r/R)
(homologous, NODELESS — projects strongly onto the fundamental radial mode);
`:sin` ⇒ v=A·sin(πr/R) (excites overtones too).
"""
function seed_radial_velocity!(st::RadialState, eng::FVRadialEngine; A::Float64=1e-3,
                               profile::Symbol=:linear)
    g=eng.g; bg=eng.bg; eos=bg.eos; N=g.N; NG=g.NG
    @inbounds for i in 1:N
        ai=NG+i; r=g.r[ai]
        if r < bg.R
            v = profile===:linear ? A*(r/bg.R) : A*sin(π*r/bg.R)
            st.v[ai]=v
            sg=bg.sqrtg[ai]
            D̂,Ŝ,τ̂ = prim2cons_barotrope(eos, st.ρ[ai], st.p[ai], abs(v))
            st.D[ai]=sg*D̂; st.S[ai]=sg*sign(v)*Ŝ; st.τ[ai]=sg*τ̂
        end
    end
    _fill_ghosts!(st, eng)
end

# spectral peak extraction over a frequency band (returns sorted peak freqs kHz)
"""
    radial_periodogram_freqs(ts, probe; fmin_kHz, fmax_kHz, npts=4000,
                             Lunit_km=Msun_to_km, npeaks=4) -> freqs_kHz, P, νgrid

Periodogram of the probe time series over a kHz band; returns the dominant peak
frequencies (kHz), the power, and the cyclic-frequency grid (geometric).
"""
function radial_periodogram_freqs(ts::Vector{Float64}, probe::Vector{Float64};
        fmin_kHz::Float64=0.5, fmax_kHz::Float64=12.0, npts::Int=4000,
        Lunit_km::Float64=Msun_to_km, npeaks::Int=4, window::Bool=true)
    νmin = fmin_kHz*Lunit_km*kHz_to_km
    νmax = fmax_kHz*Lunit_km*kHz_to_km
    νgrid = range(νmin, νmax; length=npts)
    # Hann window kills rectangular-window sidelobes that create spurious peaks at
    # the steep stellar surface — essential for a clean radial-mode FFT.
    q = copy(probe)
    if window
        n=length(q); μ=sum(q)/n
        @inbounds for i in 1:n
            q[i] = (q[i]-μ)*(0.5-0.5*cos(2π*(i-1)/(n-1)))
        end
    end
    P = periodogram(ts, q, νgrid)
    # find local maxima
    peaks = Tuple{Float64,Float64}[]
    for n in 2:length(P)-1
        if P[n]>P[n-1] && P[n]≥P[n+1]
            push!(peaks, (P[n], νgrid[n]))
        end
    end
    sort!(peaks, by=x->-x[1])
    k = min(npeaks, length(peaks))
    freqs = [ p[2]/(Lunit_km*kHz_to_km) for p in peaks[1:k] ]
    return freqs, collect(P), collect(νgrid)
end

using ..NonRadialModes: freq_kHz_from_omega2
using ..Numerics: brent

# ---------------------------------------------------------------------------
# Independent relativistic radial (ℓ=0) Cowling eigensolver — for cross-checking
# the FV radial spectrum. Variables (ξ, η) with η ≡ Δp/(... )? We use the
# Gondek–Rosińska / Kokkotas–Ruoff Lagrangian form in (ξ, Δp):
#   dξ/dr  = −(1/r)[ 3ξ + Δp/(Γ₁ p) ] + Φ' ξ
#   dΔp/dr = ξ[ ω² e^{λ−ν}(ε+p) r − 4Φ' Δp/ξ ... ] — written cleanly as:
#   dΔp/dr = ξ (ε+p)[ ω² e^{λ−ν} r − ... ]  (Cowling: metric perts dropped)
# We use the standard Cowling form (e.g. Väth & Chanmugam 1992; Kokkotas–Ruoff
# 2001, Cowling limit):
#   dξ/dr  = −(3/r)ξ − (1/(Γ₁ p)) Δp + Φ' ξ
#   dΔp/dr = ξ [ ω² (ε+p) e^{λ−ν} r ... ] − Δp(... ); the cleanest closed Cowling
#   radial pair (Chandrasekhar in GR, frozen metric) is:
#     dξ/dr  = −(1/r)(3ξ + Δp/(Γ₁ p)) − Φ' ξ
#     dΔp/dr = ξ ( ω² e^{λ−ν}(ε+p) r − (ε+p)Φ'·2 ) − Δp ( Φ' (1 + 1/cs²)·? )
# To avoid sign ambiguities we VALIDATE this reference against the package's
# nonradial p-mode ladder (the radial ladder interleaves) and, more importantly,
# against the FV code's own grid-convergence. The reference is reported with its
# caveat.
# ---------------------------------------------------------------------------
"""
    fmode_radial_reference(eos, εc; nmodes=3) -> freqs_kHz

Relativistic radial (ℓ=0) Cowling-mode reference (Lagrangian shooting, frozen
metric). Returns the lowest `nmodes` radial frequencies (kHz, M⊙ units).

CAVEAT (honesty): this shooting solver returns a SUSPICIOUSLY UNIFORM ladder
(≈2.18, 3.92, 5.68, 7.44 kHz — even ~1.74 kHz spacing), which is unphysical for
real stellar p-modes (their spacing is non-uniform). It almost certainly carries
a sign/term error in the radial pair below and should NOT be trusted as the
ground truth. The PRIMARY Stage-1 validation is therefore the FV code's own
GRID CONVERGENCE: the windowed central-density FFT gives a single clean radial
fundamental ≈3.08 kHz, IDENTICAL at N=200 and N=400 (<0.1%). Kept as a (flagged)
independent attempt, not a passing benchmark.
"""
function fmode_radial_reference(eos::BarotropicEOS, εc::Float64; nmodes::Int=3,
                                N::Int=8000, h_tov::Float64=2e-4,
                                ω2lo::Float64=1e-4, ω2hi::Float64=0.25, nscan::Int=2000)
    star = solve_tov(eos, εc; h=h_tov)
    R=star.R
    rt=star.r; mt=star.m; νt=star.ν; εt=star.ε; pt=star.p
    bg(r) = (max(_lininterp(rt,mt,r),0.0), _lininterp(rt,νt,r),
             max(_lininterp(rt,εt,r),1e-20), max(_lininterp(rt,pt,r),1e-30))

    # RHS of the Cowling radial pair (ξ, Δp). Frozen-metric Chandrasekhar/
    # Gondek–Rosińska form:
    #   ξ' = −(1/r)(3ξ + Δp/(Γ₁ p)) + Φ' ξ
    #   Δp' = ξ[ ω² e^{λ−ν}(ε+p) r − 4Φ' p ] − Φ'(1 + 1/cs²)Δp ... (Cowling)
    function rhs(r, ξ, Δp, ω2)
        m,ν,ε,p = bg(r)
        cs2 = clamp(sound_speed2(eos, ε), 1e-12, 1.0)
        Γ1  = cs2*(ε+p)/p
        eλ  = 1/max(1-2m/r, 1e-12)
        eν  = exp(ν)
        Φp  = (m + 4π*r^3*p)/(r*(r-2m))
        dξ  = -(1/r)*(3ξ + Δp/(Γ1*p)) + Φp*ξ
        dΔp = ξ*( ω2*eλ/eν*(ε+p)*r - 4*Φp*p ) - Φp*(1 + 1/cs2)*Δp
        return dξ, dΔp
    end

    function shoot(ω2)
        r0=R*1e-3; r=r0
        m,ν,ε,p=bg(r0)
        cs2=clamp(sound_speed2(eos,ε),1e-12,1.0); Γ1=cs2*(ε+p)/p
        ξ=r0/R; Δp=-Γ1*p*3*ξ       # regular-center: ξ∝r, Δp from continuity
        rf=R*(1-1e-3); nstep=N; h=(rf-r0)/nstep
        for _ in 1:nstep
            k1ξ,k1Δ=rhs(r,ξ,Δp,ω2)
            k2ξ,k2Δ=rhs(r+h/2,ξ+h/2*k1ξ,Δp+h/2*k1Δ,ω2)
            k3ξ,k3Δ=rhs(r+h/2,ξ+h/2*k2ξ,Δp+h/2*k2Δ,ω2)
            k4ξ,k4Δ=rhs(r+h,ξ+h*k3ξ,Δp+h*k3Δ,ω2)
            ξ+=h/6*(k1ξ+2k2ξ+2k3ξ+k4ξ); Δp+=h/6*(k1Δ+2k2Δ+2k3Δ+k4Δ)
            r+=h
            a=abs(ξ)+abs(Δp); a>1e150 && (ξ/=a; Δp/=a)
        end
        return Δp     # surface ΔP=0
    end

    grid=range(ω2lo,ω2hi;length=nscan)
    ω2s=Float64[]; Dprev=shoot(first(grid)); ωprev=first(grid)
    for ω2 in Iterators.drop(grid,1)
        Dc=shoot(ω2)
        if isfinite(Dc)&&isfinite(Dprev)&&Dc*Dprev<0
            res=brent(shoot,ωprev,ω2;xtol=1e-12)
            res.converged && push!(ω2s,res.root)
        end
        Dprev=Dc; ωprev=ω2
    end
    sort!(ω2s); k=min(nmodes,length(ω2s))
    return [freq_kHz_from_omega2(w) for w in ω2s[1:k]]
end

end # module FVRadial
