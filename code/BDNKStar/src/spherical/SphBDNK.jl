#=
    SphBDNK — full-frame (causal) BDNK viscous pulsations on the boundary-conforming
    (r,θ) spherical background (STAGE-3 rewrite). Two results live here.

    (1) STABILITY. Putting the surface on the coordinate line r=R makes the BDNK
    velocity RECOVERY stable on the FULL star, curing the blow-up the staircased
    Cartesian surface caused. The velocity is a PRIMITIVE evolved through the frame
    recovery (the formerly unstable part):
        ∂_t δv^i = (δS^i − w₀ δv^i − κ_Q c_s² ∂^iδε)/den,   den = max(τ_R w₀ − η, floor)
    Continuity/momentum are the clean curved-space form:
        ∂_t δρ   = −(1/√γ) ∂_i(√γ α ρ₀ δv^i)
        ∂_t δS_i = −α[∂_i δp + Φ'(δp+δε)δ_i^r] + ν_mom ∇²δS_i,   δp = c_s²δε

    (2) VISCOUS DAMPING — where it actually lives (hard-won, with two false leads).
    There are THREE transport channels, all physically ∝η̂ but with very different roles:
      • τ_R  (relaxation)  — sets the recovery stiffness den; a causal-frame knob.
      • κ_Q  (conduction)  — the heat term −κ_Q c_s²∂ε; REACTIVE. It shifts the mode
                             frequency and stabilizes, but does NOT dissipate. Proven
                             by a decoupled scan: at fixed relaxation, modal "decay"
                             FALLS as κ_Q rises (it just detunes from the numerical sink).
      • ν_mom (shear visc) — the force ν_mom∇²δS_i in the MOMENTUM. THIS is the genuine
                             modal damper: as a friction on the velocity-like momentum it
                             makes the (δρ,δS) oscillator damped, decay = ν_mom k²/2 > 0,
                             and it dissipates the quadratic mode energy monotonically.

    CRUCIAL METHOD NOTE: read the damping from the QUADRATIC ENERGY
        E = ½∫√γ [ w₀(e^Λ(δv^r)² + r²(δv^θ)²) + c_s² δε²/w₀ ]   (see `mode_energy`)
    NOT from time-domain amplitude / quadrupole-moment decay. The latter is dominated by
    Kreiss–Oliger dissipation and by DEPHASING of a non-eigenmode seed, both of which a
    smoothing operator suppresses — so every dissipative term spuriously looks like it
    REDUCES "decay." In the energy, ν_mom turns the inviscid scheme's slow grid-scale
    growth into monotonically stronger decay as η̂ rises (e.g. n=3 quadrupole seed,
    E_end/E0: 14.6 (η̂=0, growing) → 0.80 → 0.51 → 0.42 → 0.37 for η̂=0.02..0.08).

    ⚠ CONVERGENCE / NUMERICAL-vs-PHYSICAL CAVEAT (discrimination battery, 2026-06;
    see test/test_convergence.jl, figures/convergence_viscous.png):
      • QUALITATIVE result is ROBUST/PHYSICAL: ν_mom∇²δS genuinely removes mode energy
        (E_end/E0 falls monotonically with η̂), κ_Q is REACTIVE (energy unchanged), and
        net decay strengthens with resolution at low Kreiss–Oliger. Signs/trends survive
        every refinement.
      • QUANTITATIVE damping RATE Γ is NOT trustworthy at the reference η̂ (~0.05): it is
        σ_ko-dependent (spread > |Γ|, sign-flipping), does NOT converge under Nr
        refinement, and is ~100× below the analytic ν_mom·k²/2 (linearity R²=0.21). The
        inviscid scheme has a large grid-scale numerical GROWTH (~10×) masking it below
        η̂≈0.06. dt-independence is clean → contamination is SPATIAL (KO + unconverged
        grid-scale growth), not RK4 time-integration.
      • USE: trust the SIGN/TREND here; for QUANTITATIVE viscous damping use the converged
        frequency-domain solvers (AxialViscousModes <0.05% vs Bussières; PolarViscousModes
        within Nr≲200). Physical signal robust only for η̂≳0.06–0.08, low σ_ko, high Nr.

    Method of lines RK4 + Kreiss–Oliger, ghost-cell BCs (centre/surface/poles).
=#
module SphBDNK

using ..SphBackground: SphStar
using ..CowlingEvolve3D: periodogram, freq_kHz_cyclic, damping_rate

export SphBDNKState, SphBDNKEvo, setup_sphbdnk, seed_sphbdnk_l2!, seed_sphbdnk_n!,
       seed_sphbdnk_noise!,
       evolve_sphbdnk!, l2_quad_sphbdnk, lℓ_quad_sphbdnk, mode_energy,
       ell_reduced_operator,
       periodogram, freq_kHz_cyclic, damping_rate

# Self-contained reproducible Gaussian generator (keeps the core package STDLIB-only:
# Random is NOT a declared dependency). splitmix64 → uniform → Box–Muller N(0,1).
mutable struct _RNG; s::UInt64; end
@inline function _u01(r::_RNG)
    r.s += 0x9e3779b97f4a7c15
    z = r.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    z = z ⊻ (z >> 31)
    (z >> 11) * (1.0/9007199254740992.0)         # 53-bit mantissa in [0,1)
end
@inline function _gauss(r::_RNG)
    u1 = _u01(r); u1 = u1 < 1e-300 ? 1e-300 : u1
    u2 = _u01(r)
    sqrt(-2*log(u1))*cos(2π*u2)
end
_gaussvec(r::_RNG, n::Int) = [_gauss(r) for _ in 1:n]

mutable struct SphBDNKState
    δρ::Matrix{Float64}; δSr::Matrix{Float64}; δSθ::Matrix{Float64}
    δvr::Matrix{Float64}; δvθ::Matrix{Float64}                      # contravariant velocity
    ξ::Matrix{Float64}          # radial Lagrangian displacement (STRATIFICATION only;
end                             # ∂_tξ = α δv^r. Decouples identically when B≡0.)
SphBDNKState(Nr, Nθ) = SphBDNKState((zeros(Nr+2, Nθ+2) for _ in 1:6)...)
_flds(st::SphBDNKState) = (st.δρ, st.δSr, st.δSθ, st.δvr, st.δvθ, st.ξ)

struct SphBDNKEvo
    s::SphStar
    η̂::Float64; τ̂::Float64; cτ::Float64; κ̂::Float64; cκ::Float64
    ν̂::Float64; cν::Float64; den_floor::Float64; σ_ko::Float64
    # TRUE BDNK SHEAR STRESS (η_sh>0 switches it on; 0 ⇒ exactly the legacy scheme).
    # η = η_sh·ε₀, the SAME normalisation as CowlingBDNK3D and as the dissipation
    # integral of test_cross_method — so damping rates are directly comparable.
    η_sh::Float64
    dlneL::Vector{Float64}          # d(ln g_rr)/dr, background
    # √γ-weighted MIXED stress π^i_j, staged so the divergence is a plain difference
    Frr::Matrix{Float64}; Ftr::Matrix{Float64}
    Frt::Matrix{Float64}; Ftt::Matrix{Float64}
    Ppp::Matrix{Float64}            # π^φ_φ (unweighted; only ever needed locally)
    rc_frac::Float64                # centre floor for 1/r² in the θ recovery
    sbc::Vector{Float64}            # free-surface coeff  p₀′/(c_s²h₀)  (finite as c_s²→0)
    fs::Bool                        # impose the Lagrangian free-surface BC at ii=Nr
    # STRATIFICATION (g-modes). cs2ad = ADIABATIC (frozen-composition) sound speed;
    # Bstrat = p₀′/(cs2ad·h₀) − ρ₀′ = ρ₀Φ′Δ(c⁻²) the buoyancy coefficient. With
    # cs2ad ≡ s.cs2 we get Bstrat ≡ 0 and the ξ block decouples exactly, so the
    # barotropic scheme is bit-identical. Field 1 then MEANS ρ̃ = ρ̂ − Bξ, chosen so
    # δp = cs2ad·h₀·ρ̃ EXACTLY and the free-surface coefficient sbc is unchanged.
    cs2ad::Vector{Float64}
    Bstrat::Vector{Float64}
end
# Three general-frame channels, deliberately SEPARATE (independent in the general BDNK
# frame, tied only by the causality inequalities). Each = baseline + (∝η̂) part:
#   • relaxation τ_R = τ̂ + cτ·η̂  → den = τ_R w₀ − η  (recovery stiffness; cτ≥1.5 ⇒ den>0);
#   • conduction κ_Q = κ̂ + cκ·η̂  → heat term −κ_Q c_s²∂ε  (REACTIVE: frequency/stability);
#   • shear visc ν_mom = ν̂ + cν·η̂ → momentum force ν_mom∇²δS  (the GENUINE energy damper).
# Decoupling these is what let us prove κ_Q is reactive and ν_mom is the dissipator
# (scan one at a time, read the quadratic ENERGY). Physical frame: all parts ∝η̂.
function setup_sphbdnk(s::SphStar; η̂=0.0, τ̂=0.03, cτ=2.0, κ̂=0.03, cκ=1.0,
                       ν̂=0.0, cν=1.0, den_frac=0.02, σ_ko=0.02, η_sh=0.0, free_surface=false, rc_frac=0.0,
                       cs2_ad=nothing)
    wmax = maximum(s.ε0 .+ s.p0)
    Nr, Nθ, dr = s.grid.Nr, s.grid.Nθ, s.grid.dr
    dl = similar(s.eΛ)
    @inbounds for i in 1:Nr
        dl[i] = i == 1  ? (log(s.eΛ[2]) - log(s.eΛ[1]))/dr :
                i == Nr ? (log(s.eΛ[Nr]) - log(s.eΛ[Nr-1]))/dr :
                          (log(s.eΛ[i+1]) - log(s.eΛ[i-1]))/(2dr)
    end
    # FREE-SURFACE COEFFICIENT. The physical condition at r=R is Lagrangian,
    # Δp = δp + ξ^r ∂_r p₀ = 0. Differentiating (p₀ static, ∂_tξ^r = δv^r) gives the
    # time-domain form ∂_tδp = −δv^r ∂_r p₀, i.e. ∂_tδρ = −δv^r·[∂_r p₀/(c_s²h₀)].
    # Both ∂_r p₀ and c_s² vanish at the surface (p ∼ (R−r)^{n+1}) but their RATIO is
    # FINITE: measured −1.221e−4 → −1.089e−4 over the outer 5%, smoothly converging.
    sbc = similar(s.p0)
    @inbounds for i in 1:Nr
        dp = i == 1  ? (s.p0[2]-s.p0[1])/dr :
             i == Nr ? (3s.p0[Nr]-4s.p0[Nr-1]+s.p0[Nr-2])/(2dr) :
                       (s.p0[i+1]-s.p0[i-1])/(2dr)
        h0i = (s.ε0[i]+s.p0[i])/s.ρ0[i]
        sbc[i] = dp/(max(s.cs2[i],1e-30)*h0i)
    end
    cs2ad = cs2_ad === nothing ? copy(s.cs2) : collect(cs2_ad)
    Bstrat = zeros(Nr)
    if cs2_ad !== nothing
        @inbounds for i in 1:Nr
            dp  = i==1 ? (s.p0[2]-s.p0[1])/dr : i==Nr ? (s.p0[Nr]-s.p0[Nr-1])/dr :
                         (s.p0[i+1]-s.p0[i-1])/(2dr)
            dρ0 = i==1 ? (s.ρ0[2]-s.ρ0[1])/dr : i==Nr ? (s.ρ0[Nr]-s.ρ0[Nr-1])/dr :
                         (s.ρ0[i+1]-s.ρ0[i-1])/(2dr)
            h0i = (s.ε0[i]+s.p0[i])/s.ρ0[i]
            Bstrat[i] = dp/(max(cs2ad[i],1e-30)*h0i) - dρ0
        end
    end
    z() = zeros(Nr+2, Nθ+2)
    SphBDNKEvo(s, η̂, τ̂, cτ, κ̂, cκ, ν̂, cν, den_frac*τ̂*wmax, σ_ko,
               η_sh, dl, z(), z(), z(), z(), z(), rc_frac, sbc, free_surface,
               cs2ad, Bstrat)
end

@inline _h0(s,a) = (s.ε0[a-1]+s.p0[a-1])/s.ρ0[a-1]
@inline _w0i(s,a) = s.ε0[a-1]+s.p0[a-1]

function _fill_ghosts!(st::SphBDNKState, Nr, Nθ)
    @inbounds for jj in 1:Nθ
        j = jj+1
        for (F,par) in ((st.δρ,1),(st.δSr,-1),(st.δSθ,1),(st.δvr,-1),(st.δvθ,1),(st.ξ,-1))  # centre: δSr,δvr odd
            F[1,j] = par*F[2,j]; F[Nr+2,j] = F[Nr+1,j]                            # surface: outflow
        end
    end
    @inbounds for ii in 1:Nr
        i = ii+1
        for (F,par) in ((st.δρ,1),(st.δSr,1),(st.δSθ,-1),(st.δvr,1),(st.δvθ,-1),(st.ξ,1))  # poles: δSθ,δvθ odd
            F[i,1] = par*F[i,2]; F[i,Nθ+2] = par*F[i,Nθ+1]
        end
    end
end

@inline _ci(a,N) = clamp(a-1, 1, N)
@inline _Fr(s,vr,a,b) = (ar=_ci(a,s.grid.Nr); ab=_ci(b,s.grid.Nθ); s.α[ar]*(s.sqγr[ar]*s.grid.sinθ[ab])*s.ρ0[ar]*vr[a,b])
@inline _Fθ(s,vθ,a,b) = (ar=_ci(a,s.grid.Nr); ab=_ci(b,s.grid.Nθ); s.α[ar]*(s.sqγr[ar]*s.grid.sinθ[ab])*s.ρ0[ar]*vθ[a,b])
@inline _ko1(A,i,j) = (A[i-2,j]-4A[i-1,j]+6A[i,j]-4A[i+1,j]+A[i+2,j]) +
                      (A[i,j-2]-4A[i,j-1]+6A[i,j]-4A[i,j+1]+A[i,j+2])
# centre-regularized Laplacian e^{−Λ}∂²_r + (1/r²_eff)∂²_θ; r²_eff floors the 1/r²
# angular blow-up at the coordinate centre.
@inline _lap(A,i,j,eΛ,r2eff,dr2,dθ2) = (A[i+1,j]-2A[i,j]+A[i-1,j])/(eΛ*dr2) +
                                        (A[i,j+1]-2A[i,j]+A[i,j-1])/(r2eff*dθ2)

function _rhs!(d::SphBDNKState, st::SphBDNKState, e::SphBDNKEvo)
    s=e.s; Nr=s.grid.Nr; Nθ=s.grid.Nθ; dr=s.grid.dr; dθ=s.grid.dθ; σ=e.σ_ko; dr2=dr^2; dθ2=dθ^2
    _fill_ghosts!(st, Nr, Nθ)
    δρ,δSr,δSθ,δvr,δvθ,ξ = st.δρ,st.δSr,st.δSθ,st.δvr,st.δvθ,st.ξ
    # ---- TRUE BDNK SHEAR STRESS π^{ij} = −2η σ^{ij} (staged; see below) ----------
    # σ from t189_sphshear.m, EXACT for the axisymmetric polar sector on
    # γ=diag(g_rr, r², r²sin²θ). Both indices are projected with Δ^{μν}, which adds
    # the Φ′ pieces that a purely spatial gradient misses (t188). Controls: the
    # derivation has γ_ij σ^{ij}=0 identically and reduces to the textbook
    # Navier–Stokes form at g_rr=α=1 (both residuals exactly 0).
    #     θ    = ∂_r v^r + v^r(2/r + Φ′ + ½∂_r ln g_rr) + ∂_θ v^θ + cotθ v^θ
    #     σ^rr = (1/3g_rr)[2∂_r v^r + v^r(∂_r ln g_rr + 2Φ′ − 2/r) − ∂_θ v^θ − cotθ v^θ]
    #     σ^rθ = ½[Φ′v^θ/g_rr + ∂_θ v^r/r² + ∂_r v^θ/g_rr]
    #     σ^θθ = (1/3r²)[(v^r/r)(1 − rΦ′ − ½r∂_r ln g_rr) − cotθ v^θ + 2∂_θ v^θ − ∂_r v^r]
    #     σ^φφ = (csc²θ/3r²)[(v^r/r)(1 − rΦ′ − ½r∂_r ln g_rr) + 2cotθ v^θ − ∂_θ v^θ − ∂_r v^r]
    # Stored √γ-weighted and index-lowered so the divergence below is a plain
    # difference. Filled on ALL physical cells; the momentum loop only reads
    # neighbours inside that set, so π needs no ghosts of its own.
    if e.η_sh > 0
        Frr,Ftr,Frt,Ftt,Ppp = e.Frr,e.Ftr,e.Frt,e.Ftt,e.Ppp
        @inbounds for jj in 1:Nθ, ii in 1:Nr
            i=ii+1; j=jj+1
            r=s.grid.r[ii]; grr=s.eΛ[ii]; Φp=s.νp[ii]/2; dl=e.dlneL[ii]
            sn=s.grid.sinθ[jj]; ct=s.grid.cosθ[jj]/sn; sq=s.sqγr[ii]*sn
            vr=δvr[i,j]; vt=δvθ[i,j]
            dvr_r=(δvr[i+1,j]-δvr[i-1,j])/(2dr); dvr_t=(δvr[i,j+1]-δvr[i,j-1])/(2dθ)
            dvt_r=(δvθ[i+1,j]-δvθ[i-1,j])/(2dr); dvt_t=(δvθ[i,j+1]-δvθ[i,j-1])/(2dθ)
            gfac = (vr/r)*(1 - r*Φp - 0.5*r*dl)
            srr = (2dvr_r + vr*(dl + 2Φp - 2/r) - dvt_t - ct*vt)/(3grr)
            srt = 0.5*(Φp*vt/grr + dvr_t/r^2 + dvt_r/grr)
            stt = (gfac - ct*vt + 2dvt_t - dvr_r)/(3r^2)
            spp = (gfac + 2ct*vt - dvt_t - dvr_r)/(3r^2*sn^2)
            m2η = -2*e.η_sh*s.ε0[ii]
            Frr[i,j]=sq*m2η*srr*grr;  Ftr[i,j]=sq*m2η*srt*grr      # √γ π^r_r, √γ π^θ_r
            Frt[i,j]=sq*m2η*srt*r^2;  Ftt[i,j]=sq*m2η*stt*r^2      # √γ π^r_θ, √γ π^θ_θ
            Ppp[i,j]=m2η*spp*(r*sn)^2                              # π^φ_φ
        end
    end
    @inbounds for jj in (e.fs ? (1:Nθ) : (2:Nθ-1)), ii in (e.fs ? (1:Nr) : (2:Nr-1))
        i=ii+1; j=jj+1
        koOK = (2 <= ii <= Nr-1) && (2 <= jj <= Nθ-1)   # _ko1 needs i±2, one ghost only
        sqγ=s.sqγr[ii]*s.grid.sinθ[jj]; α=s.α[ii]; Φp=s.νp[ii]/2; cs2=s.cs2[ii]
        h0=_h0(s,i); w0=_w0i(s,i); eΛ=s.eΛ[ii]; r2=s.grid.r[ii]^2
        # STRATIFICATION: field 1 is ρ̃ = ρ̂ − Bξ, so δp = cs2ad·h₀·ρ̃ EXACTLY while the
        # ENERGY perturbation keeps the ξ piece: δε = h₀(ρ̃ + Bξ). Both reduce to the
        # barotropic forms when Bstrat ≡ 0.
        cs2a=e.cs2ad[ii]; Bst=e.Bstrat[ii]
        δε=h0*(δρ[i,j] + Bst*ξ[i,j]); δp=cs2a*h0*δρ[i,j]
        # BDN causality (a) is τ_Q w₀ − η > 0 for the ACTUAL shear viscosity, so the
        # true-stress channel η_sh must enter here too — otherwise η_sh can be raised
        # past the causality bound with `den` none the wiser, and the run goes acausal
        # while still looking well-posed. (Found by a convergence study that refused
        # to converge: τ̂=0.03 with η_sh=0.04 gives τ_Q w₀ − η = −9.5e−6.)
        η=(e.η̂ + e.η_sh)*s.ε0[ii]
        # frame relaxation τ_R (stiffness), conduction κ_Q (reactive heat term), and the
        # shear-viscous momentum friction ν_mom (the genuine modal damper), independently
        # parametrized — see the struct note.
        τR=e.τ̂ + e.cτ*e.η̂; κQ=e.κ̂ + e.cκ*e.η̂; νm=e.ν̂ + e.cν*e.η̂
        den=max(τR*w0-η, e.den_floor); r2eff=max(r2,(0.15*s.R)^2)
        # CONTINUITY (rest mass)
        divF=(_Fr(s,δvr,i+1,j)-_Fr(s,δvr,i-1,j))/(2dr)+(_Fθ(s,δvθ,i,j+1)-_Fθ(s,δvθ,i,j-1))/(2dθ)
        # FREE SURFACE: at ii=Nr replace continuity by ∂_tδρ = −δv^r·sbc, which is
        # ∂_t(Δp)=0. The old scheme never wrote ii=Nr at all, so for ANY oscillatory
        # eigenmode that cell is identically zero (λx=0 ⇒ x=0) — a Dirichlet-zero
        # wall for δρ at the surface, producing a spurious ∂_rδp in the last evolved
        # cell, hence a δS_r jump, hence the δv^r spike (δv^r ≈ δS_r/(g_rr w₀)).
        d.δρ[i,j] = (e.fs && ii==Nr) ? -δvr[i,j]*e.sbc[ii] - (koOK ? σ*_ko1(δρ,i,j) : 0.0) :
                                       -divF/sqγ - α*Bst*δvr[i,j] - (koOK ? σ*_ko1(δρ,i,j) : 0.0)
        # MOMENTUM = ideal pressure stress + radial gravity + shear-viscous friction
        # ν_mom∇²δS_i. As a friction on the (velocity-like) momentum this is genuinely
        # DISSIPATIVE: it makes the (δρ,δS) oscillator damped, decay = ν_mom k²/2 > 0.
        ip=min(ii+1,Nr); im=max(ii-1,1)
        dpr=(e.cs2ad[ip]*_h0(s,ip+1)*δρ[i+1,j]-e.cs2ad[im]*_h0(s,im+1)*δρ[i-1,j])/(2dr)
        dpθ=cs2a*h0*(δρ[i,j+1]-δρ[i,j-1])/(2dθ)
        # NB the buoyancy restoring force is ALREADY carried by δε here: with
        # δε = h₀(ρ̃+Bξ), Φp*(δp+δε) = Φp*h₀[(1+cs2)ρ̃ + Bξ] — exactly the derived form.
        # Adding a separate Φp*h₀*B*ξ would DOUBLE-COUNT it (it shifted f by +11%).
        d.δSr[i,j]=-α*(dpr + Φp*(δp+δε)) + νm*_lap(δSr,i,j,eΛ,r2eff,dr2,dθ2) - (koOK ? σ*_ko1(δSr,i,j) : 0.0)
        d.δSθ[i,j]=-α*dpθ + νm*_lap(δSθ,i,j,eΛ,r2eff,dr2,dθ2) - (koOK ? σ*_ko1(δSθ,i,j) : 0.0)
        # TRUE BDNK SHEAR STRESS divergence (t188): ∂_tδS_j|visc = −α[∇_iπ^i_j + Φ′π^r_j],
        # using ∇_iπ^i_j = (1/√γ)∂_i(√γ π^i_j) − Γ^m_ij π^i_m with the spherical
        # Christoffels Γ^r_rr=½∂_r ln g_rr, Γ^θ_rθ=Γ^φ_rφ=1/r, Γ^r_θθ=−r/g_rr,
        # Γ^φ_θφ=cotθ. This REPLACES nothing — it adds alongside ν_mom, which stays
        # available and inert at its default 0, so every existing test is untouched.
        if e.η_sh > 0
            r=s.grid.r[ii]; ct=s.grid.cosθ[jj]/s.grid.sinθ[jj]
            prr=e.Frr[i,j]/sqγ; ptr=e.Ftr[i,j]/sqγ
            prt=e.Frt[i,j]/sqγ; ptt=e.Ftt[i,j]/sqγ; ppp=e.Ppp[i,j]
            divr = ((e.Frr[i+1,j]-e.Frr[i-1,j])/(2dr) + (e.Ftr[i,j+1]-e.Ftr[i,j-1])/(2dθ))/sqγ -
                   (0.5*e.dlneL[ii]*prr + (ptt + ppp)/r)
            divt = ((e.Frt[i+1,j]-e.Frt[i-1,j])/(2dr) + (e.Ftt[i,j+1]-e.Ftt[i,j-1])/(2dθ))/sqγ -
                   (prt/r - r*ptr/eΛ + ct*ppp)
            d.δSr[i,j] -= α*(divr + Φp*prr)
            d.δSθ[i,j] -= α*(divt + Φp*prt)
        end
        # VELOCITY RECOVERY (BDNK frame) — stable on the boundary-conforming grid. The
        # −κ_Q c_s²∂ε heat term is REACTIVE (sets frequency/stability, not the damping).
        dεr=h0*((δρ[i+1,j]+e.Bstrat[ip]*ξ[i+1,j])-(δρ[i-1,j]+e.Bstrat[im]*ξ[i-1,j]))/(2dr)
        dεθ=h0*((δρ[i,j+1]+Bst*ξ[i,j+1])-(δρ[i,j-1]+Bst*ξ[i,j-1]))/(2dθ)
        # CENTRE REGULARISATION of the θ recovery (rc_frac>0; default 0 ⇒ inert).
        # δS_θ/r² and ∂_θδε/r² blow up as r→dr/2 — the same 1/r² the Laplacian
        # already floors via r2eff. Without it the ℓ-reduced operator carries an
        # unstable mode that is: centre-localised (peak always at r≈dr/2), GRID-SCALE
        # (f ∝ Nr: 184→1324 kHz at Nr=100→800, vs 2–10 kHz physical), with
        # maxRe·dr ≈ const so it worsens with refinement, and 99.96% δv^θ.
        # MECHANISM: δS_θ is tiny (2.4e−6 of the vector) but NOT zero, and
        # 1/(r²·den) ≈ 176 at the innermost cell amplifies it past the −1/τ_R = −16.7
        # decay δv^θ would otherwise have — the amplitude rides on δv^θ while the
        # SMALL couplings feed it (same structure as the 3+1D instability).
        # This floor HALVES maxRe (1.329→0.68 at Nr=400) and then saturates; it does
        # NOT cure the mode. Full stabilisation needs σ_ko≈0.5 (25× default), which
        # costs +2.6% on the f-mode and inflates the numerical floor 3.5×. Larger τ̂
        # makes it WORSE. The f-mode is INSENSITIVE to rc_frac (identical to 5
        # digits), so the floor is free; and in the frequency domain the mode is
        # harmless anyway (300–700 kHz, removed by the band/Q filter). It is only
        # blocking for TIME-DOMAIN use — hence free_surface/rc_frac stay opt-in.
        r2r = e.rc_frac > 0 ? max(r2, (e.rc_frac*s.R)^2) : r2
        recr=(δSr[i,j]/eΛ - w0*δvr[i,j] - κQ*cs2*(dεr/eΛ))/den
        recθ=(δSθ[i,j]/r2r - w0*δvθ[i,j] - κQ*cs2*(dεθ/r2r))/den
        d.δvr[i,j]=recr - (koOK ? σ*_ko1(δvr,i,j) : 0.0)
        d.δvθ[i,j]=recθ - (koOK ? σ*_ko1(δvθ,i,j) : 0.0)
        d.ξ[i,j]  = α*δvr[i,j]                       # ∂_tξ^r = α δv^r
    end
    return nothing
end

function _rk4!(st, e, dt, k1,k2,k3,k4, tmp)
    F=_flds(st); T=_flds(tmp)
    _rhs!(k1,st,e);  K=_flds(k1); for n in 1:6; @. T[n]=F[n]+dt/2*K[n]; end
    _rhs!(k2,tmp,e); K=_flds(k2); for n in 1:6; @. T[n]=F[n]+dt/2*K[n]; end
    _rhs!(k3,tmp,e); K=_flds(k3); for n in 1:6; @. T[n]=F[n]+dt*K[n]; end
    _rhs!(k4,tmp,e)
    K1=_flds(k1);K2=_flds(k2);K3=_flds(k3);K4=_flds(k4)
    for n in 1:6; @. F[n] += dt/6*(K1[n]+2K2[n]+2K3[n]+K4[n]); end
end

function seed_sphbdnk_l2!(st::SphBDNKState, e::SphBDNKEvo; A=1e-3)
    s=e.s
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        st.δρ[ii+1,jj+1] = A*s.ε0[ii]*(s.grid.r[ii]/s.R)*(3*s.grid.cosθ[jj]^2 - 1)
    end
end

function l2_quad_sphbdnk(st::SphBDNKState, e::SphBDNKEvo)
    s=e.s; acc=0.0
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        acc += st.δρ[ii+1,jj+1]*(3*s.grid.cosθ[jj]^2-1)*s.sqγr[ii]*s.grid.sinθ[jj]
    end
    acc*s.grid.dr*s.grid.dθ
end

"""seed an ℓ=2 perturbation with radial overtone n: δρ = A ε₀ sin(nπr/R)(3cos²θ−1)."""
function seed_sphbdnk_n!(st::SphBDNKState, e::SphBDNKEvo, n::Int; A=1e-3)
    s=e.s
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        st.δρ[ii+1,jj+1] = A*s.ε0[ii]*sin(n*π*s.grid.r[ii]/s.R)*(3*s.grid.cosθ[jj]^2 - 1)
    end
end

@inline function _legendre(l::Int, x::Float64)
    l == 0 ? 1.0 :
    l == 1 ? x :
    l == 2 ? (3x^2 - 1)/2 :
    l == 3 ? (5x^3 - 3x)/2 :
    l == 4 ? (35x^4 - 30x^2 + 3)/8 :
    l == 5 ? (63x^5 - 70x^3 + 15x)/8 :
             (231x^6 - 315x^4 + 105x^2 - 5)/16   # l=6
end

"""
    seed_sphbdnk_noise!(st, e; A=1e-4, Kr=6, Lset=2:2:4, seed=1) -> (minfac, maxδρrel)

Seed CONSTRAINT-VALID, BAND-LIMITED GAUSSIAN-NOISE initial data for the scalar
density perturbation δρ. A generic causal random "pluck" used for mode
spectroscopy: it excites the star's whole physical mode spectrum, but only on the
RESOLVED scales (band-limited), so the power lands on the low-order f/p modes and
not on the grid scale.

  δρ(r,θ) = A·ε₀(r) · [Σ_{n=1..Kr} aₙ sin(nπr/R)] · [Σ_{ℓ∈Lset} bₗ Pₗ(cosθ)]

with aₙ,bₗ ~ N(0,1) drawn from MersenneTwister(`seed`). The radial sin-basis
vanishes at r=0 (regular centre) and r=R (Δp≈0 surface); using EVEN Legendre
orders (Lset=2:2:4) makes the angular factor even at both poles — exactly the
parity `_fill_ghosts!` expects for the scalar δρ (even at centre AND poles), so
the ghost-cell regularity BCs hold by construction.

Positivity (ε₀+δε>0) is enforced by a global validity-rescale: if any cell would
go non-positive the whole field is scaled down by `minfac`. Velocities start at 0
(trivially subluminal) and are recovered by the evolution. Returns the applied
positivity factor and the realised max|δρ|/ε₀.
"""
function seed_sphbdnk_noise!(st::SphBDNKState, e::SphBDNKEvo;
                             A=1e-4, Kr::Int=6, Lset=2:2:4, seed::Int=1)
    s = e.s; Nr=s.grid.Nr; Nθ=s.grid.Nθ; R=s.R
    rng = _RNG(UInt64(seed) ⊻ 0xda3e39cb94b95bdb)
    Ls = collect(Lset)
    ar = _gaussvec(rng, Kr); bl = _gaussvec(rng, length(Ls))
    # raw field
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        r = s.grid.r[ii]; c = s.grid.cosθ[jj]
        radial = 0.0
        for n in 1:Kr; radial += ar[n]*sin(n*π*r/R); end
        ang = 0.0
        for (m,l) in enumerate(Ls); ang += bl[m]*_legendre(l, c); end
        st.δρ[ii+1,jj+1] = A*s.ε0[ii]*radial*ang
    end
    # positivity-rescale: ε₀ + δε > 0 with δε = h₀ δρ
    minfac = 1.0
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        h0 = _h0(s, ii+1); δε = h0*st.δρ[ii+1,jj+1]
        if s.ε0[ii] + δε ≤ 0
            f = 0.9*s.ε0[ii]/abs(δε)        # scale this cell back below positivity bound
            minfac = min(minfac, f)
        end
    end
    if minfac < 1.0
        @inbounds for jj in 1:Nθ, ii in 1:Nr; st.δρ[ii+1,jj+1] *= minfac; end
    end
    maxrel = 0.0
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        maxrel = max(maxrel, abs(st.δρ[ii+1,jj+1])/s.ε0[ii])
    end
    (minfac, maxrel)
end

"""ℓ-th Legendre projected density moment  Mₗ = ∫ δρ Pₗ(cosθ) √γ sinθ dr dθ."""
function lℓ_quad_sphbdnk(st::SphBDNKState, e::SphBDNKEvo, l::Int)
    s=e.s; acc=0.0
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        acc += st.δρ[ii+1,jj+1]*_legendre(l, s.grid.cosθ[jj])*s.sqγr[ii]*s.grid.sinθ[jj]
    end
    acc*s.grid.dr*s.grid.dθ
end

"""Quadratic mode energy E = ½∫√γ[ w₀(e^Λ(δv^r)² + r²(δv^θ)²) + c_s² δε²/w₀ ].
This — NOT amplitude/quadrupole decay — is the faithful viscous-dissipation diagnostic."""
function mode_energy(st::SphBDNKState, e::SphBDNKEvo)
    s=e.s; E=0.0
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        sg=s.sqγr[ii]*s.grid.sinθ[jj]; w0=s.ε0[ii]+s.p0[ii]; h0=w0/s.ρ0[ii]
        δε=h0*st.δρ[ii+1,jj+1]
        E += sg*( w0*(s.eΛ[ii]*st.δvr[ii+1,jj+1]^2 + s.grid.r[ii]^2*st.δvθ[ii+1,jj+1]^2)
                  + s.cs2[ii]*δε^2/w0 )
    end
    0.5*E*s.grid.dr*s.grid.dθ
end

"""
    ell_reduced_operator(e, l) -> Matrix

The ℓ-reduced linear operator, 5·(Nr−2) square, built by PROJECTING the 2D
`_rhs!` rather than re-deriving anything. `_rhs!` is linear in the state, so
applying it to (radial delta)×(angular harmonic) gives exact columns; for a
spherically symmetric background the operator is block-diagonal in ℓ, so the
projection isolates that block (discretely, with O(dθ²) leakage).

WHY THIS EXISTS. On the full (r,θ) grid the spectrum mixes every ℓ, every radial
overtone and the frame modes — 574 oscillatory modes at Nr=48, all with Q=1–6 —
and the f-mode simply cannot be picked out: nearest-eigenvalue continuation in
η_sh returns anything from −69% to +632% of the expected damping depending on
which eigenvalue is followed. Reduced to ℓ=2 the tower is clean and identifiable:
f=1.8896 kHz against the validated `nonradial_cowling_spectrum` value 1.8829
(0.36%), then p1..p4 at 4.16, 6.12, 8.01, 9.84 kHz.

Basis: δρ, δS_r, δv^r ~ P_ℓ(cosθ);  δS_θ, δv^θ ~ ∂_θP_ℓ. The projection uses the
same truncated quadrature (θ cells 2:Nθ−1) as `_rhs!` writes, in numerator and
denominator alike, so the normalisation is consistent.
"""
function ell_reduced_operator(e::SphBDNKEvo, l::Int=2)
    s=e.s; Nr=s.grid.Nr; Nθ=s.grid.Nθ
    l == 2 || throw(ArgumentError("only ℓ=2 harmonics are tabulated here"))
    P(c)=0.5*(3c^2-1); dP(θ)=-3*cos(θ)*sin(θ)
    ang=(Float64[P(s.grid.cosθ[j]) for j in 1:Nθ], Float64[P(s.grid.cosθ[j]) for j in 1:Nθ],
         Float64[dP(s.grid.θ[j])   for j in 1:Nθ], Float64[P(s.grid.cosθ[j]) for j in 1:Nθ],
         Float64[dP(s.grid.θ[j])   for j in 1:Nθ], Float64[P(s.grid.cosθ[j]) for j in 1:Nθ])
    # 6th entry: ξ is a RADIAL displacement, so it carries P_l like δρ and δv^r.
    jr = e.fs ? (1:Nθ) : (2:Nθ-1); ir = e.fs ? (1:Nr) : (2:Nr-1); Nrr=length(ir)
    nrm=[sum(ang[k][j]^2*s.grid.sinθ[j] for j in jr) for k in 1:length(ang)]
    # 6 fields now (ξ added for stratification); with Bstrat≡0 the ξ block is
    # decoupled and contributes Nrr zero eigenvalues, which the band/Q filter drops.
    st=SphBDNKState(Nr,Nθ); d=SphBDNKState(Nr,Nθ); nf=length(_flds(st)); L=zeros(nf*Nrr,nf*Nrr)
    for (kf,F) in enumerate(_flds(st)), (kc,ii) in enumerate(ir)
        for A in _flds(st); fill!(A,0.0); end
        for j in 1:Nθ; F[ii+1,j+1]=ang[kf][j]; end
        _rhs!(d,st,e); c=(kf-1)*Nrr+kc
        for (lf,D) in enumerate(_flds(d)), (lc,kk) in enumerate(ir)
            L[(lf-1)*Nrr+lc,c]=sum(D[kk+1,j+1]*ang[lf][j]*s.grid.sinθ[j] for j in jr)/nrm[lf]
        end
    end
    L
end

"""evolve nsteps; return (ts, ℓ=2 quadrupole moment, mode energy)."""
function evolve_sphbdnk!(st::SphBDNKState, e::SphBDNKEvo; dt, nsteps, sample=1)
    Nr=e.s.grid.Nr; Nθ=e.s.grid.Nθ
    k1=SphBDNKState(Nr,Nθ);k2=SphBDNKState(Nr,Nθ);k3=SphBDNKState(Nr,Nθ);k4=SphBDNKState(Nr,Nθ);tmp=SphBDNKState(Nr,Nθ)
    ts=Float64[]; q2=Float64[]; en=Float64[]
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt); push!(q2, l2_quad_sphbdnk(st,e)); push!(en, mode_energy(st,e))
        end
        n == nsteps && break
        _rk4!(st, e, dt, k1,k2,k3,k4, tmp)
    end
    return ts, q2, en
end

end # module SphBDNK
