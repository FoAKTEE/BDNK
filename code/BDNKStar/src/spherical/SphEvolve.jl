#=
    SphEvolve — IDEAL linearized Cowling pulsations on the boundary-conforming
    (r,θ) spherical background (STAGE-3 rewrite, Gate 2a). The surface sits on the
    coordinate line r=R, so there is no staircased boundary — the failure mode of
    the Cartesian engine.

    Linearized ideal-fluid Cowling equations (axisymmetric, evolving δε, δS_r, δS_θ;
    geometric connection terms cancel between the √γ-divergence and the Valencia
    metric-derivative source, leaving the clean curved-space form):
        ∂_t δε   = −(1/√γ)[∂_r(α√γ δS^r) + ∂_θ(α√γ δS^θ)] − αΦ' δS^r
        ∂_t δS_r = −∂_r δp − Φ'(δp + δε)
        ∂_t δS_θ = −∂_θ δp
    with √γ=e^{Λ/2}r²sinθ, δS^r=e^{−Λ}δS_r, δS^θ=δS_θ/r², δp=c_s²δε, Φ'=ν'/2.

    Method of lines: 2nd-order central FD + Kreiss–Oliger, RK4. Ghost cells carry
    the boundary conditions: centre (r→0) regularity [δε,δS_θ even, δS_r odd],
    surface (r=R) free surface [δε,δS→0 just outside, c_s²→0 makes δp→0 there],
    poles (θ=0,π) regularity for an ℓ=2,m=0 (even-parity) perturbation
    [δε,δS_r even, δS_θ odd].
=#
module SphEvolve

using ..SphBackground: SphStar
using ..CowlingEvolve3D: periodogram, freq_kHz_cyclic, damping_rate
using ..Units: Msun_to_km, kHz_to_km

export SphState, SphEvo, setup_sphevo, seed_sph_l2!, evolve_sph!, l2_quad_sph,
       periodogram, freq_kHz_cyclic, damping_rate

# 2D fields with one ghost layer on each side: index 0..Nr+1, 0..Nθ+1
mutable struct SphState
    δε::Matrix{Float64}; δSr::Matrix{Float64}; δSθ::Matrix{Float64}
end
SphState(Nr, Nθ) = SphState(zeros(Nr+2, Nθ+2), zeros(Nr+2, Nθ+2), zeros(Nr+2, Nθ+2))

struct SphEvo
    s::SphStar
    σ_ko::Float64
end
setup_sphevo(s::SphStar; σ_ko::Float64=0.02) = SphEvo(s, σ_ko)

# ghost index helper: physical cell (i,j) → array index (i+1,j+1)
@inline _a(i,j) = (i+1, j+1)

# fill ghost cells from the BC parity rules
function _fill_ghosts!(st::SphState, Nr, Nθ)
    δε, δSr, δSθ = st.δε, st.δSr, st.δSθ
    @inbounds for jj in 1:Nθ
        j = jj+1
        # centre r→0 (i=0 ghost ↔ i=1): δε,δSθ even, δSr odd
        δε[1,j]  =  δε[2,j];   δSθ[1,j] =  δSθ[2,j];  δSr[1,j] = -δSr[2,j]
        # surface r=R (i=Nr+1 ghost): FREE surface — zero-gradient (outflow); Δp=0 is
        # enforced physically by c_s²→0, not by an over-stiff Dirichlet δε=0
        δε[Nr+2,j] = δε[Nr+1,j]; δSr[Nr+2,j] = δSr[Nr+1,j]; δSθ[Nr+2,j] = δSθ[Nr+1,j]
    end
    @inbounds for ii in 1:Nr
        i = ii+1
        # pole θ=0 (j=0 ghost ↔ j=1): δε,δSr even, δSθ odd
        δε[i,1]  =  δε[i,2];   δSr[i,1]  =  δSr[i,2];  δSθ[i,1]  = -δSθ[i,2]
        # pole θ=π (j=Nθ+1 ghost): same parity
        δε[i,Nθ+2] = δε[i,Nθ+1]; δSr[i,Nθ+2] = δSr[i,Nθ+1]; δSθ[i,Nθ+2] = -δSθ[i,Nθ+1]
    end
    return nothing
end

@inline _ko1(A,i,j) = (A[i-2,j]-4A[i-1,j]+6A[i,j]-4A[i+1,j]+A[i+2,j]) +
                      (A[i,j-2]-4A[i,j-1]+6A[i,j]-4A[i,j+1]+A[i,j+2])

# background helpers indexed by GHOST array index a (physical radial index = a-1)
@inline _h0(s,a) = (s.ε0[a-1]+s.p0[a-1])/s.ρ0[a-1]        # specific enthalpy
@inline _w0(s,a) = s.ε0[a-1]+s.p0[a-1]
@inline _δpG(s,q,a,b) = s.cs2[a-1]*_h0(s,a)*q[a,b]         # δp = c_s² h₀ δρ at ghost idx
@inline _FrG(s,q,Sr,a,b) = s.α[a-1]*(s.sqγr[a-1]*s.grid.sinθ[b-1])*s.ρ0[a-1]*(Sr[a,b]/(s.eΛ[a-1]*_w0(s,a)))
@inline _FθG(s,q,Sθ,a,b) = s.α[a-1]*(s.sqγr[a-1]*s.grid.sinθ[b-1])*s.ρ0[a-1]*(Sθ[a,b]/(s.grid.r[a-1]^2*_w0(s,a)))

function _rhs!(d::SphState, st::SphState, e::SphEvo)
    s = e.s; Nr = s.grid.Nr; Nθ = s.grid.Nθ; dr = s.grid.dr; dθ = s.grid.dθ; σ = e.σ_ko
    _fill_ghosts!(st, Nr, Nθ)
    q, δSr, δSθ = st.δε, st.δSr, st.δSθ            # q ≡ δρ (rest-mass perturbation)
    @inbounds for jj in 2:Nθ-1, ii in 2:Nr-1
        i = ii+1; j = jj+1
        sqγ = s.sqγr[ii]*s.grid.sinθ[jj]
        α = s.α[ii]; Φp = s.νp[ii]/2; h0 = _h0(s,i)
        # CONTINUITY (rest mass, α in flux, no source): ∂_t δρ = −(1/√γ)∂_i(√γ α ρ₀ δv^i)
        divF = (_FrG(s,q,δSr,i+1,j) - _FrG(s,q,δSr,i-1,j))/(2dr) +
               (_FθG(s,q,δSθ,i,j+1) - _FθG(s,q,δSθ,i,j-1))/(2dθ)
        d.δε[i,j] = -divF/sqγ - σ*_ko1(q,i,j)
        # MOMENTUM (overall lapse α): ∂_t δS_r = −α[∂_rδp + Φ'(δp+δε)], ∂_t δS_θ = −α∂_θδp
        dpr = (_δpG(s,q,i+1,j) - _δpG(s,q,i-1,j))/(2dr)
        dpθ = (_δpG(s,q,i,j+1) - _δpG(s,q,i,j-1))/(2dθ)
        δp_here = s.cs2[ii]*h0*q[i,j]; δε_here = h0*q[i,j]
        d.δSr[i,j] = -α*(dpr + Φp*(δp_here + δε_here)) - σ*_ko1(δSr,i,j)
        d.δSθ[i,j] = -α*dpθ - σ*_ko1(δSθ,i,j)
    end
    return nothing
end

function _rk4!(st::SphState, e::SphEvo, dt, k1,k2,k3,k4, tmp)
    F=(st.δε,st.δSr,st.δSθ); T=(tmp.δε,tmp.δSr,tmp.δSθ)
    _rhs!(k1,st,e); K=(k1.δε,k1.δSr,k1.δSθ); for n in 1:3; @. T[n]=F[n]+dt/2*K[n]; end
    _rhs!(k2,tmp,e); K=(k2.δε,k2.δSr,k2.δSθ); for n in 1:3; @. T[n]=F[n]+dt/2*K[n]; end
    _rhs!(k3,tmp,e); K=(k3.δε,k3.δSr,k3.δSθ); for n in 1:3; @. T[n]=F[n]+dt*K[n]; end
    _rhs!(k4,tmp,e)
    K1=(k1.δε,k1.δSr,k1.δSθ);K2=(k2.δε,k2.δSr,k2.δSθ);K3=(k3.δε,k3.δSr,k3.δSθ);K4=(k4.δε,k4.δSr,k4.δSθ)
    for n in 1:3; @. F[n] += dt/6*(K1[n]+2K2[n]+2K3[n]+K4[n]); end
end

"""seed an ℓ=2, m=0 density perturbation δε = A ε₀ (r/R)(3cos²θ−1)."""
function seed_sph_l2!(st::SphState, e::SphEvo; A::Float64=1e-3)
    s = e.s
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        st.δε[ii+1,jj+1] = A*s.ε0[ii]*(s.grid.r[ii]/s.R)*(3*s.grid.cosθ[jj]^2 - 1)
    end
end

"""√γ-weighted ℓ=2 quadrupole moment of δε (the non-radial diagnostic)."""
function l2_quad_sph(st::SphState, e::SphEvo)
    s = e.s; acc = 0.0
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        acc += st.δε[ii+1,jj+1]*(3*s.grid.cosθ[jj]^2 - 1)*s.sqγr[ii]*s.grid.sinθ[jj]
    end
    acc * s.grid.dr * s.grid.dθ
end

"""evolve nsteps; return (ts, l2-quadrupole, central δε)."""
function evolve_sph!(st::SphState, e::SphEvo; dt::Float64, nsteps::Int, sample::Int=1)
    Nr=e.s.grid.Nr; Nθ=e.s.grid.Nθ
    k1=SphState(Nr,Nθ);k2=SphState(Nr,Nθ);k3=SphState(Nr,Nθ);k4=SphState(Nr,Nθ);tmp=SphState(Nr,Nθ)
    ts=Float64[]; q2=Float64[]; qc=Float64[]
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt); push!(q2, l2_quad_sph(st,e)); push!(qc, st.δε[2, Nθ÷2+1])
        end
        n == nsteps && break
        _rk4!(st, e, dt, k1,k2,k3,k4, tmp)
    end
    return ts, q2, qc
end

end # module SphEvolve
