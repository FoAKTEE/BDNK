#=
    PolarViscousModes — POLAR (even-parity) VISCOUS quasi-normal modes of a neutron
    star, in the relativistic Cowling approximation (STAGE 3). The even-parity, viscous
    counterpart of the repo's axial viscous-QNM track (repro/axial_*.jl): where the
    axial sector hosts viscosity-DRIVEN η-modes (no perfect-fluid counterpart), the
    polar sector already carries the f/p modes, so viscosity DAMPS them.

    Method (frequency domain as a matrix eigenvalue problem). The linearised polar
    BDNK–Cowling system (the same operator the time-domain `SphBDNK` engine evolves)
    is reduced to 1D for a fixed ℓ — the angular θ-derivatives collapse to ℓ(ℓ+1)
    factors for a P_ℓ(cosθ) harmonic — giving ∂_t U = L_ℓ U with the 5 radial fields
        U = (ρ̂, R̂≡δŜ_r, T̂≡δŜ_θ, û^r, û^θ).
    The QNMs are the eigenvalues λ of L_ℓ:  a mode ∝ e^{λt} = e^{-iωt} has
        frequency  f = |Im λ| / (2π·M⊙_km·kHz_km),   damping rate  γ = −Re λ  (γ>0 ⇒ decays).
    The full complex spectrum exposes BOTH the (weakly-damped, high-Q) f/p modes and
    the heavily-damped non-hydrodynamic frame/relaxation modes.

    Transport channels (each = baseline + ∝η̂, as in SphBDNK):
        relaxation τ_R=τ̂+cτη̂ (den=τ_R w₀−η),  conduction κ_Q=κ̂+cκη̂ (REACTIVE heat term),
        shear viscosity ν_mom=ν̂+cνη̂  (momentum force ν_mom[e^{−Λ}∂²_r − ℓ(ℓ+1)/r²]δŜ —
        the genuine DAMPER).
    Validated: η̂→0 reproduces the f/p tower (NonRadialModes, ~5% conservative-frame
    systematic) with γ→0; η̂>0 gives γ>0 rising ∝η̂, p-modes damping more than f
    (higher k ⇒ more shear), consistent with the SphBDNK time-domain energy dissipation.

    NOTE: ν_mom carries a reactive (frequency-shifting) part as well as the dissipative
    one — track modes by complex continuation in η̂ (small steps). Absolute frequencies
    inherit the conservative-form/free-surface systematic; the damping is differential
    and robust.
=#
module PolarViscousModes

using ..SphBackground: SphStar, build_sphstar
using ..EquationOfState: BarotropicEOS
using ..Units: Msun_to_km, kHz_to_km
using LinearAlgebra: eigvals

export polar_bdnk_operator, polar_qnm, qnm_freq_kHz, qnm_damping

"""QNM angular eigenvalue λ (of ∂_tU=LU) → cyclic frequency in kHz."""
qnm_freq_kHz(λ) = abs(imag(λ)) / (2π * Msun_to_km * kHz_to_km)
"""QNM damping rate γ = −Re λ (geometric units; γ>0 ⇒ the mode decays)."""
qnm_damping(λ) = -real(λ)

"""
    polar_bdnk_operator(s, ℓ; η̂=0, τ̂=0.12, cτ=0, κ̂=0, cκ=0, ν̂=0, cν=1, den_frac=0.02)

Dense matrix L of the 1D ℓ-reduced linearised polar BDNK–Cowling operator on the
boundary-conforming radial grid of `s` (a `SphStar`). Fields stacked as
[ρ̂; R̂; T̂; û^r; û^θ], each length `s.grid.Nr`. Eigenvalues of L are the QNMs.
Defaults isolate the SHEAR-viscous damping (no reactive heat term, fixed relaxation).
"""
function polar_bdnk_operator(s::SphStar, ℓ::Int; η̂=0.0, τ̂=0.12, cτ=0.0, κ̂=0.0,
                             cκ=0.0, ν̂=0.0, cν=1.0, den_frac=0.02, free_surface=false,
                             Γratio=1.0, sbc_lapse=false, cs2_ad=nothing)
    Nr=s.grid.Nr; dr=s.grid.dr; dr2=dr^2
    α=s.α; ρ0=s.ρ0; ε0=s.ε0; w0=s.ε0 .+ s.p0; h0=w0 ./ ρ0; ce2=s.cs2
    eΛ=s.eΛ; Φp=s.νp ./ 2; r=s.grid.r; sγr=s.sqγr
    # ---- STRATIFICATION (g-modes) --------------------------------------------------
    # The background is BAROTROPIC, so its own sound speed ce2 = (dp/dε)_struct governs
    # the structure and the Brunt–Väisälä frequency vanishes identically ⇒ NO g-modes.
    # Γratio = Γ₁/Γ_struct > 1 makes the ADIABATIC sound speed of the perturbations
    # exceed the equilibrium one, giving N² > 0 and a genuine g-branch:
    #     cs2 = Γratio·ce2 ,  Δ(c⁻²) = 1/ce2 − 1/cs2 = (1−1/Γratio)/ce2 ,
    #     B   = ρ₀Φ′Δ(c⁻²) = p₀′/(cs2·h₀) − ρ₀′ ,  N² = Φ′²Δ(c⁻²)α²/e^Λ.
    # B is the ONLY new background coefficient; B>0 ⟺ N²>0 ⟺ convectively stable.
    # The state gains ONE field, the radial Lagrangian displacement ξ, with ∂_tξ = α u^r,
    # and the first field becomes the pressure-carrying density ρ̃ = ρ̂ − Bξ so that
    # δp = cs2·h₀·ρ̃ EXACTLY, while δε = h₀(ρ̃ + Bξ).
    # BAROTROPIC LIMIT IS EXACT: Γratio=1 ⇒ B≡0 ⇒ every new term vanishes and ξ decouples,
    # so we drop it entirely and the operator is bit-identical to the 5-field version.
    # (Independently derived twice — Lagrangian-displacement and two-sound-speed routes —
    #  which agree via Ap = −cs2·h₀·B. The BDNK viscous/relaxation sector is UNTOUCHED:
    #  buoyancy enters only through the ideal pressure closure and the gravity term.)
    # cs2_ad: the ADIABATIC (frozen-composition / index-Γ) sound-speed profile of a
    # genuinely stratified star, e.g. from solve_tov_idealgas + brunt_vaisala. Preferred
    # over the Γratio knob, which is an artificial uniform rescaling of a barotropic star.
    strat = cs2_ad !== nothing || Γratio != 1.0
    cs2 = cs2_ad !== nothing ? collect(cs2_ad) : (Γratio != 1.0 ? (Γratio .* ce2) : ce2)
    B = zeros(Nr)
    if strat
        @inbounds for i in 1:Nr
            dp  = i==1 ? (s.p0[2]-s.p0[1])/dr : i==Nr ? (s.p0[Nr]-s.p0[Nr-1])/dr :
                         (s.p0[i+1]-s.p0[i-1])/(2dr)
            dρ0 = i==1 ? (ρ0[2]-ρ0[1])/dr    : i==Nr ? (ρ0[Nr]-ρ0[Nr-1])/dr :
                         (ρ0[i+1]-ρ0[i-1])/(2dr)
            B[i] = dp/(max(cs2[i],1e-30)*h0[i]) - dρ0
        end
    end
    τR=τ̂+cτ*η̂; κQ=κ̂+cκ*η̂; νm=ν̂+cν*η̂
    den=max.(τR.*w0 .- η̂.*ε0, den_frac*τ̂*maximum(w0))
    c = sγr .* α .* ρ0; G = cs2 .* h0; r2eff = max.(r.^2, (0.15*s.R)^2)
    # FREE-SURFACE COEFFICIENT, ported from SphBDNK (where it is already validated).
    # The physical condition at r=R is LAGRANGIAN, Δp = δp + ξ^r∂_r p₀ = 0; differentiating
    # (p₀ static, ∂_tξ^r = δv^r) gives ∂_tρ̃ = −δv^r·sbc with sbc = ∂_r p₀/(c_s² h₀).
    # Both ∂_r p₀ and c_s² vanish at the surface but their RATIO is finite. Written in ρ̃
    # (which carries δp) the coefficient is unchanged by stratification.
    # WITHOUT this the outer BC is a bare zero-gradient ghost at r=R−dr/2, which mimics a
    # NO-SLIP WALL and injects an η̂-independent offset into γ — the defect that made the
    # f-mode rate and its apparent η̂-scaling exponent unreliable (and produced a since-
    # retracted γ∝√η̂ claim).
    sbc = zeros(Nr)
    if free_surface
        @inbounds for i in 1:Nr
            dp = i==1  ? (s.p0[2]-s.p0[1])/dr :
                 i==Nr ? (3s.p0[Nr]-4s.p0[Nr-1]+s.p0[Nr-2])/(2dr) :
                         (s.p0[i+1]-s.p0[i-1])/(2dr)
            sbc[i] = dp/(max(cs2[i],1e-30)*h0[i])
            # OPEN QUESTION (settled empirically below): one derivation argues ∂_tξ^r = α u^r
            # implies the surface row should carry a factor α. SphBDNK (validated) has none.
            sbc_lapse && (sbc[i] *= α[i])
        end
    end
    eg(x)=[x[1]; x; x[Nr]]; cg=eg(c); Gg=eg(G); h0g=eg(h0); Bg=eg(B)
    function applyL(ρ,R,T,ur,uθ,ξ)
        ρg=[ρ[1];ρ;ρ[Nr]]; Rg=[-R[1];R;R[Nr]]; Tg=[T[1];T;T[Nr]]   # centre parity / surface zero-grad
        urg=[-ur[1];ur;ur[Nr]]; ξg=[-ξ[1];ξ;ξ[Nr]]
        dρ=zeros(Nr);dR=zeros(Nr);dT=zeros(Nr);dur=zeros(Nr);duθ=zeros(Nr);dξ=zeros(Nr)
        @inbounds for i in 1:Nr
            ip=i+1
            if free_surface && i == Nr
                dρ[i] = -ur[i]*sbc[i]                       # ∂_t(Δp) = 0 at the surface
            else
                dρ[i] = -(cg[ip+1]*urg[ip+1]-cg[ip-1]*urg[ip-1])/(2dr*sγr[i]) +
                        ℓ*(ℓ+1)*α[i]*ρ0[i]*uθ[i] - α[i]*B[i]*ur[i]
            end
            lapR=(Rg[ip+1]-2R[i]+Rg[ip-1])/(eΛ[i]*dr2) - ℓ*(ℓ+1)*R[i]/r2eff[i]
            dR[i] = -α[i]*((Gg[ip+1]*ρg[ip+1]-Gg[ip-1]*ρg[ip-1])/(2dr) +
                           Φp[i]*h0[i]*(1+cs2[i])*ρ[i] +
                           Φp[i]*h0[i]*B[i]*ξ[i]) + νm*lapR      # ← buoyancy restoring force
            lapT=(Tg[ip+1]-2T[i]+Tg[ip-1])/(eΛ[i]*dr2) - ℓ*(ℓ+1)*T[i]/r2eff[i]
            dT[i] = -α[i]*cs2[i]*h0[i]*ρ[i] + νm*lapT
            # δε = h₀(ρ̃ + Bξ) enters the REACTIVE heat term only (κQ; zero by default)
            dεr=((h0g[ip+1]*(ρg[ip+1]+Bg[ip+1]*ξg[ip+1])) -
                 (h0g[ip-1]*(ρg[ip-1]+Bg[ip-1]*ξg[ip-1])))/(2dr)
            dur[i]=(R[i]/eΛ[i] - w0[i]*ur[i] - κQ*cs2[i]*dεr/eΛ[i])/den[i]
            duθ[i]=(T[i]/r[i]^2 - w0[i]*uθ[i] - κQ*cs2[i]*h0[i]*(ρ[i]+B[i]*ξ[i])/r[i]^2)/den[i]
            dξ[i]= α[i]*ur[i]                                    # ∂_tξ^r = α u^r
        end
        dρ,dR,dT,dur,duθ,dξ
    end
    nf = strat ? 6 : 5
    n=nf*Nr; L=zeros(n,n); e=zeros(n); zc=zeros(Nr)
    @inbounds for k in 1:n
        fill!(e,0.0); e[k]=1.0
        ξv = strat ? view(e,5Nr+1:6Nr) : zc
        a,b,d,f,g,x=applyL(view(e,1:Nr),view(e,Nr+1:2Nr),view(e,2Nr+1:3Nr),
                           view(e,3Nr+1:4Nr),view(e,4Nr+1:5Nr),ξv)
        L[1:Nr,k]=a; L[Nr+1:2Nr,k]=b; L[2Nr+1:3Nr,k]=d; L[3Nr+1:4Nr,k]=f; L[4Nr+1:5Nr,k]=g
        strat && (L[5Nr+1:6Nr,k]=x)
    end
    L
end

"""
    polar_qnm(eos, εc; l=2, η̂=0.0, Nr=120, nmodes=2, nstep=8, frame...) -> Vector{ComplexF64}

The lowest `nmodes` polar QNM eigenvalues λ (f, p₁) for harmonic ℓ=`l` at viscosity `η̂`,
obtained by **complex continuation**: the η̂=0 ideal modes (resolution-clean, the lowest
two `|Im|` candidates) are tracked to finite η̂ in `nstep` substeps by nearest-complex
matching. Use `qnm_freq_kHz`/`qnm_damping` on each.

SCOPE / CAVEATS (adversarially reviewed). This 5-field conservative+recovery operator
RELIABLY tracks only the **f and p₁** modes — the recovery/frame degrees of freedom add
extra eigenvalues that mix with the physical tower above p₁ (the matrix 3rd mode does NOT
match the shooting p₂; use NonRadialModes for the higher inviscid tower). The operator is
only marginally stable: spurious low-freq modes appear and PROLIFERATE with Nr, which is
why the naive "lowest |Im|" pick is not resolution-robust. Continuation from the converged
ideal modes avoids them only for **Nr ≲ 200** — at higher Nr the spurious modes are dense
enough to still break the track (use Nr≈96–160). This solver is therefore SEMI-QUANTITATIVE.
!!! danger "Outer-BC defect (2026-07-19 audit) — what NOT to trust"
    The outer boundary is a bare zero-gradient (Neumann) ghost imposed at a point that is not
    even the surface (`SphBackground.jl:71` puts the outermost cell centre at R−dr/2). On a
    truncated free surface that **mimics a no-slip wall** and injects a one-cell dissipation
    artifact — an **η̂-independent OFFSET added to γ**. Consequences, measured:
      * the damping **SIGN/TREND** is robust (γ>0, p₁ damps more than f, monotone in η̂) — SAFE;
      * the **p₁ SLOPE** dγ/dη̂ over η̂∈[0.01,0.04], Nr≈110 matches the linear dissipation
        integral to ~15 % (`test_cross_method.jl`) — SAFE, and this is the anchor to use;
      * the **f-mode absolute rate and its apparent η̂-SCALING EXPONENT are NOT reliable** — the
        offset makes an affine γ=offset+slope·η̂ look like a spurious sub-linear power law
        (this generated a since-RETRACTED "γ∝√η̂" claim; see VALIDATION.md §6). Fit
        AFFINE-vs-power-law before ever quoting an exponent.
    The ~5–8 % frequency systematic at η̂→0 and the "+20–30 % reactive pull" by η̂≈0.03 are
    consequences of THIS boundary condition, not of frame choice. Also require η̂ < τ̂ (≈0.12)
    or the BDNK recovery denominator goes negative and is floored over the whole star.
Report the p₁ damping SLOPE and sign/trend as physical; treat the f-mode rate, all absolute
frequencies, and any scaling EXPONENT as unreliable until the free-surface BC is derived and
the constant `νm` is replaced by η(r)=η̂(ε₀+p₀). For PRECISION QNMs use the full-GR shooting
solver (the axial counterpart `AxialViscousModes` reproduces benchmarks to ~1 %). Extra
keywords pass to `polar_bdnk_operator` (τ̂, cτ, κ̂, cκ, ν̂, cν).
"""
function polar_qnm(eos::BarotropicEOS, εc::Float64; l::Int=2, η̂=0.0, Nr::Int=120,
                   nmodes::Int=2, nstep::Int=8, fmin=0.4, fmax=15.0, qfac=0.6, gtol=0.1,
                   warn=true, frame...)
    s = build_sphstar(eos, εc; Nr=Nr, Nθ=2)
    # high-Q, non-growing, in-band candidates, sorted by frequency
    filt(ev) = sort([e for e in ev if imag(e) > 1e-5 && fmin < qnm_freq_kHz(e) < fmax &&
                       abs(real(e)) < qfac*abs(imag(e)) && real(e) < gtol*abs(imag(e))],
                    by=e->abs(imag(e)))
    ref = filt(eigvals(polar_bdnk_operator(s, l; η̂=0.0, frame...)))   # ideal, resolution-clean
    k = min(nmodes, length(ref)); tracked = ref[1:k]
    if η̂ != 0 && k > 0                                                # continuation 0 → η̂
        for ηs in range(0.0, η̂; length=nstep+1)[2:end]
            cv = filt(eigvals(polar_bdnk_operator(s, l; η̂=ηs, frame...)))
            isempty(cv) && break
            tracked = [cv[argmin(abs.(cv .- t))] for t in tracked]    # nearest-complex
        end
    end
    if warn && k < nmodes
        @warn "polar_qnm: $k clean modes (< nmodes=$nmodes); this 5-field operator reliably tracks only the low (f,p₁) modes — use NonRadialModes for the higher inviscid tower."
    end
    tracked
end

end # module PolarViscousModes
