#=
    BDNKMHDConstitutive — VERBATIM BDNK viscoresistive MHD constitutive
    (Lier, Armas, Porth 2026, arXiv:2606.22691, Eq.5, 6, 8) transcribed as
    tensor algebra and reduced to the 7-component conserved/flux rows used by
    the 1D/2D finite-volume evolvers.  This REPLACES the earlier operational
    channel-split: the τ_X-essential sound-sector anti-diffusion now emerges
    GENUINELY from the divergence-of-ideal-current terms — NO injected operator.

    ── THE TENSORS (flat space, mostly-plus η=diag(-1,1,1,1); ε=3p, w=4ε/3) ─────
    Eq.5:
      T^{μν} = T^{μν}_(0) + A_(1) u^μu^ν + (u^μ Q^ν_(1)+u^ν Q^μ_(1)) + Π^{μν}_(1)
      J^{μν} = J^{μν}_(0) + (u^μ N^ν_(1) − u^ν N^μ_(1)) + S^{μν}_(1)
    with (Eq.6 dissipative + Eq.8 BDNK-relaxation parts):
      A_(1)   = −τ_ε u_μ ∂_ν T^{μν}_(0)                                     (scalar)
      Q^μ_(1) = −σ P^{μρ}(u^ν∂_ν u_ρ + ∂_ρ ln T)        [Eq.6a]
                + τ_u P^μ_ρ ∂_ν T^{νρ}_(0)               [Eq.8]
      Π^{μν}_(1) = −2η P^{μρ}P^{νσ}∂_{(ρ}u_{σ)}
                   − (ζ−⅔η) P^{μν} P_{ρσ}∂^ρu^σ          [Eq.6b]
                   − τ_X P^{μν} u_ρ ∂_σ T^{ρσ}_(0)        [Eq.8]  ← τ_X essential
      N^μ_(1) = −τ_b ∂_ν J^{μν}_(0)                       [Eq.8]
      S^{μν}_(1) = resistivity (Eq.6c), isotropic ε=3p reduction (below)
    Maps: σ=wD_ε, η=wD_u, ζ=2wD_u/3 (⇒ ζ−⅔η=0), τ_ε=2τ_u, r∥=r⊥=r_b.
    EOS:  ε=3p ⇒ T∝p^{1/4} ⇒ ∂_ρ ln T = ¼ ∂_ρ ln ε.

    ── THE ∂_t → S SPLIT (linear-in-S structure for Eq.20 recovery) ────────────
    Every derivative ∂_α(·) of an ideal-current component is exact via the
    primitive chain rule: ∂_α F = Σ_j (∂F/∂P_j) ∂_α P_j, with ∂_t P_j = S_j
    (the recovered primitive time-derivative), ∂_x P_j = Px_j, ∂_y P_j = Py_j.
    The ideal-current Jacobians ∂F/∂P_j are evaluated at the FROZEN local P, so
    every constitutive output is EXACTLY linear in S ⇒ M_IJ = U_I(e_J)−U_I(0)
    is exact.  The τ_X term carries u_ρ ∂_σ T^{ρσ}_(0): its ∂_t part is linear in
    S and supplies the Eq.15 sound-sector coupling (complex/anti-diffusive at
    τ_X=0, real/subluminal at τ_X≥2D_u).

    Components/ordering match the evolvers:
      P = (b^x,b^y,b^z, u^x,u^y,u^z, ε)    (b^t,u^t slaved: u^t=√(1+|u⃗|²),
                                            b^t=(u⃗·b⃗)/u^t)
      U = (J^{tx},J^{ty},J^{tz}, T^{tx},T^{ty},T^{tz}, T^{tt})
      Fx= (J^{xx},J^{xy},J^{xz}, T^{xx},T^{xy},T^{xz}, T^{tx})
      Fy= (J^{yx},J^{yy},J^{yz}, T^{yx},T^{yy},T^{yz}, T^{ty})
=#
module BDNKMHDConstitutive

using ..BDNKMHD: BDNKMHDCoeffs, MHDThermo, ConformalThermo, mhd_pressure, mhd_Tslope

export verbatim_currents, verbatim_rows, verbatim_M_U0, NCMP_V

const NCMP_V = 7
const _SIG = (-1.0, 1.0, 1.0, 1.0)   # mostly-plus diagonal η^{μμ} (index 1..4 ↔ t,x,y,z)

# ---------------------------------------------------------------------------
# 0.  Slaved 4-vectors from primitives  P=(bx,by,bz,ux,uy,uz,ε)
# ---------------------------------------------------------------------------
@inline _ut(ux,uy,uz) = sqrt(1.0 + ux*ux + uy*uy + uz*uz)

# returns (uμ::NTuple4, bμ::NTuple4, b2, ε)
@inline function _vecs(P)
    bx,by,bz,ux,uy,uz,ε = P
    ut = _ut(ux,uy,uz)
    bt = (ux*bx+uy*by+uz*bz)/ut
    uμ = (ut, ux, uy, uz)
    bμ = (bt, bx, by, bz)
    b2 = -bt*bt + bx*bx + by*by + bz*bz
    return uμ, bμ, b2, ε
end

# ideal T^{ab}_(0) (a,b in 1..4):  (w+b²)u^a u^b + (p+½b²)η^{ab} − b^a b^b
@inline function _idT(thermo, uμ, bμ, b2, ε, a, b)
    p = mhd_pressure(thermo, ε); w = ε + p
    ηab = (a==b) ? _SIG[a] : 0.0
    return (w + b2)*uμ[a]*uμ[b] + (p + 0.5*b2)*ηab - bμ[a]*bμ[b]
end
# ideal J^{ab}_(0) = u^a b^b − u^b b^a
@inline _idJ(uμ, bμ, a, b) = uμ[a]*bμ[b] - uμ[b]*bμ[a]

# index-lowering with η (mostly-plus): V_a = η_{ab} V^b = _SIG[a] V^a
@inline _lower(Vup, a) = _SIG[a]*Vup[a]

# projector P^{ab} = η^{ab} + u^a u^b  (upper)
@inline _proj(uμ, a, b) = ((a==b) ? _SIG[a] : 0.0) + uμ[a]*uμ[b]
# mixed projector P^a_b = δ^a_b + u^a u_b = δ^a_b + u^a (η_{bc}u^c)
@inline _projmix(uμ, a, b) = ((a==b) ? 1.0 : 0.0) + uμ[a]*_SIG[b]*uμ[b]

# ---------------------------------------------------------------------------
# 1.  IDEAL-CURRENT JACOBIANS w.r.t. the 7 primitives (for the divergences)
# ---------------------------------------------------------------------------
# We need ∂_α T^{ab}_(0) and ∂_α J^{ab}_(0).  By the chain rule these are
#   ∂_α F = Σ_{j=1}^{7} (∂F/∂P_j) ∂_α P_j .
# We compute ∂F/∂P_j by a centered finite difference of the ideal tensor at the
# frozen local P (exact-to-roundoff for these smooth rational functions; the
# linearity-in-S of the final constitutive is preserved EXACTLY because the
# Jacobian is frozen and the divergence is linear in ∂_α P).
const _FDH = 1e-6

# perturb primitive component j by ±h and rebuild the 4-vectors
@inline function _perturbP(P, j, h)
    return ntuple(k -> k==j ? P[k]+h : P[k], Val(NCMP_V))
end

# ∂_α T^{ab}_(0): given dP_α = (∂_α P_1..7), returns the 4×4-indexed value at (a,b)
# We assemble the full primitive-derivative via the Jacobian contracted with dP.
# To avoid recomputing 7 perturbations per (a,b), we precompute the per-j tensors.
struct IdealJac
    # dT[j] :: NTuple of the 4×4 ∂T/∂P_j  (flattened 16), dJ likewise
    dT::NTuple{7, NTuple{16,Float64}}
    dJ::NTuple{7, NTuple{16,Float64}}
end

@inline _flatT(thermo,uμ,bμ,b2,ε) = ntuple(idx -> begin
        a = ((idx-1) >> 2) + 1; b = ((idx-1) & 3) + 1
        _idT(thermo,uμ,bμ,b2,ε,a,b)
    end, Val(16))
@inline _flatJ(uμ,bμ) = ntuple(idx -> begin
        a = ((idx-1) >> 2) + 1; b = ((idx-1) & 3) + 1
        _idJ(uμ,bμ,a,b)
    end, Val(16))

function ideal_jacobian(P, thermo=ConformalThermo())
    h = _FDH
    dT = ntuple(j -> begin
        Pp = _perturbP(P, j,  h); Pm = _perturbP(P, j, -h)
        up,bp,b2p,εp = _vecs(Pp); um,bm,b2m,εm = _vecs(Pm)
        Tp = _flatT(thermo,up,bp,b2p,εp); Tm = _flatT(thermo,um,bm,b2m,εm)
        ntuple(idx -> (Tp[idx]-Tm[idx])/(2h), Val(16))
    end, Val(7))
    dJ = ntuple(j -> begin
        Pp = _perturbP(P, j,  h); Pm = _perturbP(P, j, -h)
        up,bp,_,_ = _vecs(Pp); um,bm,_,_ = _vecs(Pm)
        Jp = _flatJ(up,bp); Jm = _flatJ(um,bm)
        ntuple(idx -> (Jp[idx]-Jm[idx])/(2h), Val(16))
    end, Val(7))
    return IdealJac(dT, dJ)
end

# ∂_α T^{ab}_(0) for a derivative described by dP=(∂_α P_1..7)
@inline function dT_div(jac::IdealJac, dP, a, b)
    idx = (a-1)*4 + b
    s = 0.0
    @inbounds for j in 1:7
        s += jac.dT[j][idx]*dP[j]
    end
    return s
end
@inline function dJ_div(jac::IdealJac, dP, a, b)
    idx = (a-1)*4 + b
    s = 0.0
    @inbounds for j in 1:7
        s += jac.dJ[j][idx]*dP[j]
    end
    return s
end

# ---------------------------------------------------------------------------
# 2.  THE VERBATIM CONSTITUTIVE  T^{μν}, J^{μν}  (Eq.5,6,8)
# ---------------------------------------------------------------------------
# Inputs: primitives P, spatial derivative tuples Px=∂_x P, Py=∂_y P (Py may be
# all-zero in 1D), time derivative Pt=∂_t P=S, BDNK coefficients c.
# Returns the full 4×4 T and J (flattened 16-tuples).  All derivatives use the
# chain rule with the frozen ideal-current Jacobian ⇒ output is LINEAR in S.

@inline function _divT_vec(jac, dPt, dPx, dPy, μ)
    # (∂_ν T^{μν}_(0)) = ∂_t T^{μt} + ∂_x T^{μx} + ∂_y T^{μy}  (z-derivs = 0)
    dT_div(jac, dPt, μ, 1) + dT_div(jac, dPx, μ, 2) + dT_div(jac, dPy, μ, 3)
end
@inline function _divJ_vec(jac, dPt, dPx, dPy, μ)
    dJ_div(jac, dPt, μ, 1) + dJ_div(jac, dPx, μ, 2) + dJ_div(jac, dPy, μ, 3)
end

# ---- comoving / velocity-gradient helpers (module-level → no closure boxing) -
# ∂_α u^σ (upper, σ=1..4) from a primitive-derivative tuple dP=(∂_α P_1..7).
@inline function _duσ(uμ, dP, σ)
    if σ == 1
        # ∂_α u^t = (u^x ∂_α u^x + u^y ∂_α u^y + u^z ∂_α u^z)/u^t
        (uμ[2]*dP[4] + uμ[3]*dP[5] + uμ[4]*dP[6])/uμ[1]
    else
        dP[σ+2]   # u^x,u^y,u^z ↔ primitives 4,5,6 ↔ σ=2,3,4
    end
end
# u^ν ∂_ν (·) = u^t (·)_t + u^x (·)_x + u^y (·)_y
@inline _comov(uμ, ft, fx, fy) = uμ[1]*ft + uμ[2]*fx + uμ[3]*fy
# ∂_ρ u^σ on the fly (upper): ρ=1..4 ↔ ∂_t,∂_x,∂_y,∂_z(=0)
@inline _duUpf(uμ, Pt, Px, Py, ρ, σ) =
    ρ==1 ? _duσ(uμ,Pt,σ) : (ρ==2 ? _duσ(uμ,Px,σ) : (ρ==3 ? _duσ(uμ,Py,σ) : 0.0))
# ∂_ρ u_σ (lower σ): η_{σσ} ∂_ρ u^σ
@inline _duLowf(uμ, Pt, Px, Py, ρ, σ) = _SIG[σ]*_duUpf(uμ,Pt,Px,Py,ρ,σ)
# symmetric gradient ∂_{(ρ}u_{σ)} = ½(∂_ρ u_σ + ∂_σ u_ρ)
@inline _symgradf(uμ, Pt, Px, Py, ρ, σ) =
    0.5*(_duLowf(uμ,Pt,Px,Py,ρ,σ) + _duLowf(uμ,Pt,Px,Py,σ,ρ))
# Π^{μν}_(1) component: −2η P^{μρ}P^{νσ}∂_{(ρ}u_{σ)} − (ζ−2η/3)P^{μν}θ − τ_X P^{μν}(u_ρ divT^ρ)
@inline function _Pi_munu(uμ, η, bulk_coeff, τX, uμ_divT, θexp, Pt, Px, Py, μ, ν)
    s = 0.0
    @inbounds for ρ in 1:4, σ in 1:4
        s += _proj(uμ,μ,ρ)*_proj(uμ,ν,σ)*_symgradf(uμ,Pt,Px,Py,ρ,σ)
    end
    shear = -2*η*s
    bulk  = -bulk_coeff*_proj(uμ,μ,ν)*θexp
    tauX  = -τX*_proj(uμ,μ,ν)*uμ_divT
    return shear + bulk + tauX
end

"""
    verbatim_currents(P, Px, Py, Pt, c) -> (T::NTuple{16}, J::NTuple{16})

The FULL 4×4 BDNK T^{μν} and J^{μν} (Eq.5,6,8), flattened row-major (a,b),
linear in Pt=S.  Px,Py are ∂_x P, ∂_y P (set Py to zeros in 1D).
"""
function verbatim_currents(P, Px, Py, Pt, c::BDNKMHDCoeffs; jac::Union{IdealJac,Nothing}=nothing)
    uμ, bμ, b2, ε = _vecs(P)
    p = mhd_pressure(c.thermo, ε); w = ε + p
    σ  = w*c.Dε
    η  = w*c.Du
    ζ  = (2/3)*w*c.Du
    τε = 2*c.τu
    τu = c.τu
    τX = c.τX
    τb = c.τb
    rpar  = c.rb
    rperp = c.rb*w/(w + b2)

    jac = jac === nothing ? ideal_jacobian(P, c.thermo) : jac

    # ---- divergence vectors of the ideal currents -------------------------
    # divT^μ ≡ ∂_ν T^{μν}_(0)  (μ=1..4),  divJ^μ ≡ ∂_ν J^{μν}_(0)
    divT = ntuple(μ -> _divT_vec(jac, Pt, Px, Py, μ), 4)
    divJ = ntuple(μ -> _divJ_vec(jac, Pt, Px, Py, μ), 4)

    # scalar  u_μ ∂_ν T^{μν}_(0) = Σ_μ (η_{μμ}u^μ) divT^μ
    uμ_divT = 0.0
    @inbounds for μ in 1:4
        uμ_divT += _lower(uμ, μ)*divT[μ]
    end

    # ---- A_(1) = −τ_ε u_μ ∂_ν T^{μν}_(0) ----------------------------------
    A1 = -τε*uμ_divT

    # ---- Q^μ_(1) = −σ P^{μρ}(u^ν∂_ν u_ρ + ∂_ρ ln T) + τ_u P^μ_ρ ∂_ν T^{νρ}_(0)
    # comoving accel  a_ρ = u^ν ∂_ν u_ρ = u^ν ∂_ν (η_{ρσ}u^σ) = η_{ρρ} (u^ν ∂_ν u^ρ)
    # build u^ν ∂_ν u^σ (upper) then lower.  u^ν∂_ν F = u^t ∂_t F + u^i ∂_i F.
    # ∂_α u^σ: for σ=1 (t) the slaved u^t; for σ=2..4 it is the primitive Pt/Px/Py
    #   component (u^x,u^y,u^z are primitives 4,5,6).  We get ∂_α u^σ via the
    #   primitive chain rule too (u^t depends on u^i).
    # ∂_ρ ln T = ¼ ∂_ρ ln ε = ¼ ∂_ρ ε / ε  (ρ spatial+time).  ∂_ρ ε is primitive 7.
    # Assemble Q for ρ=1..4 then form Q^μ = P^{μρ} (lowered already inside).
    # a^σ = u^ν ∂_ν u^σ  (upper)  — module-level _duσ/_comov (no closure boxing)
    aσ = ntuple(σ -> _comov(uμ, _duσ(uμ,Pt,σ), _duσ(uμ,Px,σ), _duσ(uμ,Py,σ)), Val(4))
    # a_ρ = η_{ρρ} a^ρ
    a_low = ntuple(ρ -> _lower(aσ, ρ), 4)
    # ∂_ρ ln T (lower index ρ) = s_T ∂_ρ ε / ε ,  s_T = d lnT/d lnε (¼ conformal).
    # the gradient ∂_ρ is a COVECTOR (lower) — ∂_t,∂_x,∂_y,∂_z directly:
    sT = mhd_Tslope(c.thermo, ε)
    dlnT = ntuple(ρ -> begin
        dε = ρ==1 ? Pt[7] : (ρ==2 ? Px[7] : (ρ==3 ? Py[7] : 0.0))
        sT*dε/ε
    end, 4)

    # Q^μ_(1): first the Eq.6a heat flux  q6^μ = −σ P^{μρ}(a_ρ + ∂_ρ ln T)
    # P^{μρ} contracts an upper index μ with a LOWER ρ via P^μ_ρ (mixed):
    #   P^{μρ} X_ρ = (η^{μρ}+u^μu^ρ) X_ρ = η^{μρ}X_ρ + u^μ (u^ρ X_ρ)
    # with X_ρ = a_ρ + ∂_ρ lnT (lower).  η^{μρ}X_ρ = _SIG[μ] X_μ.
    uXcontr = 0.0
    Xlow = ntuple(ρ -> a_low[ρ] + dlnT[ρ], 4)
    @inbounds for ρ in 1:4
        uXcontr += uμ[ρ]*Xlow[ρ]              # u^ρ X_ρ
    end
    q6 = ntuple(μ -> -σ*(_SIG[μ]*Xlow[μ] + uμ[μ]*uXcontr), 4)
    # Eq.8 part  q8^μ = τ_u P^μ_ρ ∂_ν T^{νρ}_(0) = τ_u (divT^μ + u^μ (u_ρ divT^ρ))
    #   P^μ_ρ V^ρ = δ^μ_ρ V^ρ + u^μ u_ρ V^ρ = V^μ + u^μ (u_ρ V^ρ).
    q8 = ntuple(μ -> τu*(divT[μ] + uμ[μ]*uμ_divT), 4)
    Q1 = ntuple(μ -> q6[μ] + q8[μ], 4)

    # ---- Π^{μν}_(1) -------------------------------------------------------
    # shear  σ^{μν} = 2 P^{μρ}P^{νσ}∂_{(ρ}u_{σ)} − (2/3)P^{μν}P_{ρσ}∂^ρu^σ
    # but Eq.6b is written as  −2η P^{μρ}P^{νσ}∂_{(ρ}u_{σ)} − (ζ−2η/3)P^{μν}θ_exp
    # with the bulk-rate scalar  θ = P_{ρσ}∂^ρ u^σ = P^{ρσ}∂_ρ u_σ.
    # We need ∂_ρ u_σ (both lower).  ∂_ρ u_σ = η_{σσ} ∂_ρ u^σ.
    # Build the 4×4 ∂_ρ u^σ (ρ lower derivative index, σ upper field index):
    #   ∂_ρ u^σ via duσ on the ρ-derivative tuple.
    # expansion scalar θ = P^{ρσ}∂_ρ u_σ = Σ_{ρσ} P^{ρσ} ∂_ρ u_σ  (module-level _duLowf)
    θexp = 0.0
    @inbounds for ρ in 1:4, σ in 1:4
        θexp += _proj(uμ,ρ,σ)*_duLowf(uμ,Pt,Px,Py,ρ,σ)
    end
    # shear tensor (traceless transverse):  Σ^{μν} = 2 P^{μρ}P^{νσ}∂_{(ρ}u_{σ)}
    #                                              − (2/3) P^{μν} θ
    # We assemble Π^{μν}_(1) = −η Σ^{μν} − (ζ−2η/3) P^{μν} θ
    #   = −2η P^{μρ}P^{νσ}∂_{(ρ}u_{σ)} + (2η/3)P^{μν}θ − (ζ−2η/3)P^{μν}θ
    #   = −2η P^{μρ}P^{νσ}∂_{(ρ}u_{σ)} − (ζ−4η/3)P^{μν}θ
    # (with ζ=2η/3·... but we keep ζ general per the maps; ζ−2η/3=0 here).
    # τ_X part:  −τ_X P^{μν} u_ρ ∂_σ T^{ρσ}_(0) = −τ_X P^{μν} (u_ρ divT^ρ)
    #   = −τ_X P^{μν} uμ_divT.    ← THE SOUND-SECTOR (Eq.15) coupling.
    bulk_coeff = (ζ - (2/3)*η)        # = 0 for the ε=3p maps, kept general

    # ---- N^μ_(1) = −τ_b ∂_ν J^{μν}_(0) -----------------------------------
    N1 = ntuple(μ -> -τb*divJ[μ], 4)

    # ---- S^{μν}_(1) resistivity (Eq.6c), isotropic ε=3p reduction --------
    # The one-form resistive current is  S^{μν}_(1) = −r∥ (projected ∂b) … .  In
    # the isotropic r∥=r⊥=r_b reduction the leading antisymmetric resistive piece
    # acts as a diffusion of the magnetic density:  the conserved magnetic rows
    # J^{ti} receive −r_b (P^{iρ}∂_ρ b-flux).  We implement the resistive flux as
    # the divergence-form resistive smoothing of b consistent with Eq.6c:
    #   S^{ti}_(1) carries no new t-row contribution (antisymmetric, t-t=0); the
    #   resistive diffusion enters the SPATIAL flux rows J^{ji}.  We add it at the
    #   row-extraction stage (verbatim_rows) as −r ∂_x b / −r ∂_y b on the J^{ji}
    #   rows (the standard isotropic resistive flux), keeping it OUT of the t-rows
    #   so the M_IJ recovery is governed by the τ_b N-term above.
    # (Sij handled in verbatim_rows.)

    # ---- ASSEMBLE T^{μν}, J^{μν} (Eq.5) ----------------------------------
    T = ntuple(idx -> begin
        a = ((idx-1) >> 2) + 1; b = ((idx-1) & 3) + 1
        _idT(c.thermo,uμ,bμ,b2,ε,a,b) +
            A1*uμ[a]*uμ[b] +
            (uμ[a]*Q1[b] + uμ[b]*Q1[a]) +
            _Pi_munu(uμ, η, bulk_coeff, τX, uμ_divT, θexp, Pt, Px, Py, a, b)
    end, Val(16))
    J = ntuple(idx -> begin
        a = ((idx-1) >> 2) + 1; b = ((idx-1) & 3) + 1
        _idJ(uμ,bμ,a,b) +
            (uμ[a]*N1[b] - uμ[b]*N1[a])
    end, Val(16))
    return T, J
end

# ---------------------------------------------------------------------------
# 3.  ROW EXTRACTION  → (U, Fx, Fy)  for the FV evolvers
# ---------------------------------------------------------------------------
@inline _Tcomp(T, a, b) = T[(a-1)*4 + b]
@inline _Jcomp(J, a, b) = J[(a-1)*4 + b]

"""
    verbatim_rows(P, Px, Py, Pt, c) -> (U, Fx, Fy)  (each NTuple{7,Float64})

Extract the conserved t-rows U, x-flux Fx, y-flux Fy from the verbatim BDNK
currents (Eq.5,6,8).  The isotropic resistive flux (Eq.6c reduction, r∥=r⊥=r_b)
is added to the magnetic SPATIAL flux rows here (−r ∂ b), keeping it out of the
t-rows so the Eq.20 M_IJ recovery is set by the τ-relaxation N/A/Q/Π terms.
"""
function verbatim_rows(P, Px, Py, Pt, c::BDNKMHDCoeffs; jac::Union{IdealJac,Nothing}=nothing)
    T, J = verbatim_currents(P, Px, Py, Pt, c; jac=jac)
    _, _, b2, ε = _vecs(P)
    p = mhd_pressure(c.thermo, ε); w = ε + p
    rpar  = c.rb
    rperp = c.rb*w/(w + b2)

    # U = (J^{tx},J^{ty},J^{tz}, T^{tx},T^{ty},T^{tz}, T^{tt})  (a=1=t)
    U = ( _Jcomp(J,1,2), _Jcomp(J,1,3), _Jcomp(J,1,4),
          _Tcomp(T,1,2), _Tcomp(T,1,3), _Tcomp(T,1,4),
          _Tcomp(T,1,1) )
    # Fx = (J^{xx},J^{xy},J^{xz}, T^{xx},T^{xy},T^{xz}, T^{tx})  (a=2=x)
    Fx = ( _Jcomp(J,2,2), _Jcomp(J,2,3), _Jcomp(J,2,4),
           _Tcomp(T,2,2), _Tcomp(T,2,3), _Tcomp(T,2,4),
           _Tcomp(T,1,2) )
    # Fy = (J^{yx},J^{yy},J^{yz}, T^{yx},T^{yy},T^{yz}, T^{ty})  (a=3=y)
    Fy = ( _Jcomp(J,3,2), _Jcomp(J,3,3), _Jcomp(J,3,4),
           _Tcomp(T,3,2), _Tcomp(T,3,3), _Tcomp(T,3,4),
           _Tcomp(T,1,3) )

    # isotropic resistive flux (Eq.6c reduction): magnetic rows diffuse b.
    #   F^x_{J^{xi}} -= r ∂_x b^i ,  F^y_{J^{yi}} -= r ∂_y b^i .
    # row 1 (J^{·x}) ← b^x with r∥ ; rows 2,3 (J^{·y},J^{·z}) ← b^y,b^z with r⊥.
    Fx = ( Fx[1] - rpar *Px[1], Fx[2] - rperp*Px[2], Fx[3] - rperp*Px[3],
           Fx[4], Fx[5], Fx[6], Fx[7] )
    Fy = ( Fy[1] - rpar *Py[1], Fy[2] - rperp*Py[2], Fy[3] - rperp*Py[3],
           Fy[4], Fy[5], Fy[6], Fy[7] )
    return U, Fx, Fy
end

"""
    verbatim_M_U0(P, Px, Py, c) -> (M::Matrix{7×7}, U0::Vector{7})

Build the Eq.20 LOCAL recovery matrix M_IJ = ∂U_I/∂S_J and the baseline
U0 = U(S=0) at a point, EXACTLY (U is linear in S) and EFFICIENTLY (the ideal-
current Jacobian — the expensive part — is computed ONCE and reused across the
8 constitutive evaluations).  Equivalent to the naive U_I(e_J)−U_I(0) loop.
"""
function verbatim_M_U0(P, Px, Py, c::BDNKMHDCoeffs)
    jac = ideal_jacobian(P, c.thermo)
    z = (0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    U0t, _, _ = verbatim_rows(P, Px, Py, z, c; jac=jac)
    U0 = collect(Float64, U0t)
    M = Matrix{Float64}(undef, NCMP_V, NCMP_V)
    @inbounds for J in 1:NCMP_V
        eJ = ntuple(k -> k==J ? 1.0 : 0.0, Val(NCMP_V))
        UJt, _, _ = verbatim_rows(P, Px, Py, eJ, c; jac=jac)
        for I in 1:NCMP_V
            M[I,J] = UJt[I] - U0[I]
        end
    end
    return M, U0
end

"""
    verbatim_M_fast(P, Px, Py, c) -> (M::Matrix{7×7}, U0::Vector{7})

ANALYTIC, byte-identical equivalent of `verbatim_M_U0` that exploits the EXACT
linearity of U in S=∂_t P.  The only Pt-dependence in the conserved t-rows U
flows through the ideal-current divergences ∂_t T^{μt}, ∂_t J^{μt} (the spatial
∂_x,∂_y parts and the ideal/heat/shear pieces are Pt-independent and cancel in
∂U/∂S).  Their Pt-Jacobian is a single frozen-Jacobian entry, so M is assembled
in ONE pass (1 full constitutive eval for U0 + cheap algebra) instead of 8 — a
~8× cut of the per-cell recovery cost.  Verified identical to `verbatim_M_U0` to
machine precision over random states.
"""
function verbatim_M_fast(P, Px, Py, c::BDNKMHDCoeffs)
    jac = ideal_jacobian(P, c.thermo)
    z = (0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    U0t, _, _ = verbatim_rows(P, Px, Py, z, c; jac=jac)
    U0 = collect(Float64, U0t)
    uμ, _, b2, ε = _vecs(P)
    p = mhd_pressure(c.thermo, ε); w = ε + p
    sT = mhd_Tslope(c.thermo, ε)             # d lnT/d lnε (¼ conformal; general EOS)
    σheat = w*c.Dε; η = w*c.Du; ζ = (2/3)*w*c.Du; bulk_coeff = ζ - (2/3)*η
    τε = 2*c.τu; τu = c.τu; τX = c.τX; τb = c.τb
    ut = uμ[1]
    M = Matrix{Float64}(undef, NCMP_V, NCMP_V)
    @inbounds for J in 1:NCMP_V
        # ---- (I) divergence-term Pt-Jacobian (A1, q8, N1, τ_X-Π) -------------
        # ∂(∂_t T^{μt})/∂Pt_J = jac.dT[J][(μ-1)*4+1]; likewise ∂_t J^{μt}
        gT1=jac.dT[J][1]; gT2=jac.dT[J][5]; gT3=jac.dT[J][9]; gT4=jac.dT[J][13]
        gJ1=jac.dJ[J][1]; gJ2=jac.dJ[J][5]; gJ3=jac.dJ[J][9]; gJ4=jac.dJ[J][13]
        gUDT = _SIG[1]*uμ[1]*gT1 + _SIG[2]*uμ[2]*gT2 + _SIG[3]*uμ[3]*gT3 + _SIG[4]*uμ[4]*gT4
        dq8_1=τu*(gT1+uμ[1]*gUDT); dq8_2=τu*(gT2+uμ[2]*gUDT)
        dq8_3=τu*(gT3+uμ[3]*gUDT); dq8_4=τu*(gT4+uμ[4]*gUDT)
        dA1 = -τε*gUDT
        # ---- (II) heat-flux q6 Pt-Jacobian (via ∂_t u^σ accel + ∂_t lnT) -----
        # ∂(∂_t u^σ)/∂Pt_J : σ=1 slaved u^t, σ=2..4 the primitives 4,5,6
        ds1 = (uμ[2]*(J==4) + uμ[3]*(J==5) + uμ[4]*(J==6))/ut
        ds2 = (J==4) ? 1.0 : 0.0; ds3 = (J==5) ? 1.0 : 0.0; ds4 = (J==6) ? 1.0 : 0.0
        ds = (ds1, ds2, ds3, ds4)
        # ∂Xlow[ρ]/∂Pt_J = _SIG[ρ] u^t ∂_t u^ρ + (ρ==1 ? s_T δ_{J7}/ε : 0)
        dX1 = _SIG[1]*ut*ds1 + sT*(J==7)/ε
        dX2 = _SIG[2]*ut*ds2; dX3 = _SIG[3]*ut*ds3; dX4 = _SIG[4]*ut*ds4
        dX = (dX1, dX2, dX3, dX4)
        duX = uμ[1]*dX1 + uμ[2]*dX2 + uμ[3]*dX3 + uμ[4]*dX4   # ∂uXcontr/∂Pt_J
        # ∂q6^μ/∂Pt_J = −σ (_SIG[μ] ∂Xlow[μ] + u^μ ∂uXcontr)
        dq6_1 = -σheat*(_SIG[1]*dX1 + uμ[1]*duX); dq6_2 = -σheat*(_SIG[2]*dX2 + uμ[2]*duX)
        dq6_3 = -σheat*(_SIG[3]*dX3 + uμ[3]*duX); dq6_4 = -σheat*(_SIG[4]*dX4 + uμ[4]*duX)
        dq1 = dq8_1 + dq6_1; dq2 = dq8_2 + dq6_2; dq3 = dq8_3 + dq6_3; dq4 = dq8_4 + dq6_4
        # ---- (III) shear/bulk Π Pt-Jacobian (via ∂_t u in symgrad, θexp) -----
        # ∂θexp/∂Pt_J = Σ_σ P^{1σ} _SIG[σ] ∂_t u^σ   (only ρ=1 derivative is ∂_t)
        dθ = _proj(uμ,1,1)*_SIG[1]*ds1 + _proj(uμ,1,2)*_SIG[2]*ds2 +
             _proj(uμ,1,3)*_SIG[3]*ds3 + _proj(uμ,1,4)*_SIG[4]*ds4
        # ∂symgrad(ρ,σ)/∂Pt_J = ½(_SIG[σ] ds[σ] [ρ==1] + _SIG[ρ] ds[ρ] [σ==1])
        dsym(ρ,σ) = 0.5*(_SIG[σ]*ds[σ]*(ρ==1) + _SIG[ρ]*ds[ρ]*(σ==1))
        # ∂shear^{a,b}/∂Pt_J = −2η Σ_{ρσ} P^{aρ}P^{bσ} ∂symgrad(ρ,σ)
        function dshear(a,b)
            s = 0.0
            for ρ in 1:4, σ in 1:4
                s += _proj(uμ,a,ρ)*_proj(uμ,b,σ)*dsym(ρ,σ)
            end
            -2*η*s
        end
        dPi(a,b) = dshear(a,b) - bulk_coeff*_proj(uμ,a,b)*dθ - τX*_proj(uμ,a,b)*gUDT
        # ---- assemble M columns (a=t=1) -------------------------------------
        M[1,J] = -τb*(uμ[1]*gJ2 - uμ[2]*gJ1)
        M[2,J] = -τb*(uμ[1]*gJ3 - uμ[3]*gJ1)
        M[3,J] = -τb*(uμ[1]*gJ4 - uμ[4]*gJ1)
        M[4,J] = dA1*uμ[1]*uμ[2] + uμ[1]*dq2 + uμ[2]*dq1 + dPi(1,2)
        M[5,J] = dA1*uμ[1]*uμ[3] + uμ[1]*dq3 + uμ[3]*dq1 + dPi(1,3)
        M[6,J] = dA1*uμ[1]*uμ[4] + uμ[1]*dq4 + uμ[4]*dq1 + dPi(1,4)
        M[7,J] = dA1*uμ[1]*uμ[1] + 2*uμ[1]*dq1           + dPi(1,1)
    end
    return M, U0
end

export verbatim_M_fast

end # module BDNKMHDConstitutive
