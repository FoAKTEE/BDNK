#=
    CowlingBDNK3D — FULL-FRAME (causal, hyperbolic) first-order BDNK viscous
    hydrodynamics, linearized on the frozen TOV background, 3+1D Cowling.

    The upgrade over CowlingEvolve3D's leading (Navier–Stokes-limit) viscous
    force: here the BDNK frame relaxation times (τ_ε,τ_P,τ_Q) are kept, so the
    dissipative corrections enter the CONSERVED densities and the system is the
    genuine hyperbolic BDNK theory (Bemfica–Disconzi–Noronha, PRX 12 021044),
    not parabolic Navier–Stokes. This removes the diffusive CFL blow-up: the NS
    limit is unstable for η̂≳1, the causal BDNK system is not.

    General-frame constitutive relations (BDN eq. finaltheory, single-component σ=0):
      T^{μν} = (ε+A)u^μu^ν + (P+Π)Δ^{μν} − 2ησ^{μν} + u^μQ^ν + u^νQ^μ
      A   = τ_ε[Dε + (ε+P)θ]
      Π   = −ζθ + τ_P[Dε + (ε+P)θ]
      Q^ν = τ_Q[(ε+P)Du^ν + c_s² Δ^{νl}∂_lε]          (D=u·∇, θ=∇·u)

    Linearized about the static star (u₀ at rest, leading order), the conserved
    densities carry ∂_t of the primitives, giving a DIAGONAL time-derivative
    recovery. D and θ are COVARIANT, and on the static background u^μ=(1/α,0,0,0)
    gives Dδε=(1/α)∂_tδε and (Du^i)_lin=(1/α)∂_tδvⁱ, so each recovered derivative
    is a PROPER-time one and converting to coordinate time multiplies by α
    (derived exactly in exact_solutions/mathematica/t188_recovery.m):
      δE  = δε + τ_ε((1/α)∂_tδε + w₀θ)      ⇒ ∂_tδε = α[(δE−δε)/τ_ε − w₀θ]
      δSⁱ = w₀δvⁱ + (τ_Q w₀−η)(1/α)∂_tδvⁱ + τ_Q c_s²γ^{il}∂_lδε
                                     ⇒ ∂_tδvⁱ = α(δSⁱ−w₀δvⁱ−τ_Q c_s²γ^{il}∂_lδε)/(τ_Q w₀−η)
    with θ = ∇_μu^μ = (1/(α√γ))∂_i(α√γ δvⁱ) = ∂_iδvⁱ + (Φ′+A′/A)(n·δv), NOT the
    flat divergence (t187).
    with (τ_Q w₀ − η) > 0 the BDN causality condition (a) ρτ_Q>η. The conserved
    densities evolve by ∇_μδT^{μν}=0; the primitives are advanced by the recovered
    ∂_t. Scheme: method-of-lines RK4 + Kreiss–Oliger, same metric/gravity
    structure as CowlingEvolve3D (leading weak-field divergences).

    SCOPE: linearized BDNK; the frame is real (causal, hyperbolic) and the
    dissipation is the full shear+bulk+heat first-order BDNK stress. The three
    problems that once blocked this engine are all closed — (1),(2) by construction,
    (3) by the stability analysis below — and the engine has a first-class QNM
    API: `setup_bdnk3d(...; frame=(a₁,a₂))` installs a Clarisse hydrodynamic frame,
    `seed_bdnk_ylm!` seeds any (ℓ,m) on the full cube, `evolve_bdnk3d_moments!`
    records arbitrary moments, and `bdnk3d_qnm` runs the gated two-estimator
    read-out with the η̂=0 difference protocol. See test/test_bdnk3d_qnm.jl.
    (1) THE FULL-STAR SURFACE TREATMENT — `cut=true` (cut cells, `setup_bdnk3d`)
        runs the FULL star r≤R with NO excision. See the status note there; the
        long list of failed regularizers that used to live in this file was, in
        hindsight, six different ways of working around a stair-stepped RIGID WALL
        where the physics needs a FREE SURFACE.
    (2) The curved-space covariant (α√γ-weighted) divergence — the viscous
        momentum sector (t184/t186), the expansion and shear (t187/t188), and the
        lapse in both recovery equations are exact on the curved background.
    VALIDATION (cut cells, full star, no excision, matrix-pencil estimator, T=500,
    τ̂=2, σ_ko=0.004; slope from the η̂=0→0.05 interval so the KO floor is subtracted):
        N        32        44        64
        f err  −0.699%   −0.343%   −0.192%   vs the validated eigensolver
        dγ/dη̂  0.09287   0.09170   0.08985   = 106.2%, 104.8%, 102.7% of the
                                              independent dissipation integral
        KO γ₀  5.696e−4  3.198e−4  1.324e−4  (falls like dx^2.4)
    The frequency converges at observed order ~1.5–2.2 — the first clean 2nd-order
    frequency convergence in this engine (it was ~0.5 with the mask) — and the
    damping ratio descends monotonically toward 100%, which argues the NEWTONIAN-form
    integral is the right reference (an earlier reading, that the ×0.95 relativistic
    correction put the truth near 108%, is not supported by the sequence).
    (3) THE STABILITY CEILING — UNDERSTOOD (2026-09). The linear stability condition
        of the 4×4 longitudinal sound sector this engine integrates is
                       τ_Q w₀ c_s²  >  (4/3) η          (NOT τ_Q w₀ > η),
        reproduced to ratio 1.0000 in 6/6 (c_s², η̂) combinations by a frozen-
        coefficient flat-space analysis. Under the Clarisse map τ = a·η/w₀ it becomes
        the η-INDEPENDENT a₂ c_s² > 4/3. Its violation is SUFFICIENT-NOT-NECESSARY on a
        grid: the unstable band lies ABOVE a critical wavenumber k_c = C/η̂
        (C = 0.4599/0.3019/0.1753 for F₁/F₂/F₃, minimum over the star at r/R≈0.994),
        and the 2nd-order centred stencil caps k_eff = sin(kh)/h ≤ 1/h = N/23.0, so the
        band is only RESOLVED for N > 353/232/135 (anchor star, η̂=0.03). That is why
        N=32/44 looked clean and N≥64 did not with a constant τ̂, and why F₁ — which
        violates the bound in 100% of cells — evolves stably for 46 periods through 12
        decades of decay at N≤40. The criterion "unstable iff k_c < 1/h", moved to the
        cheap η̂ axis at N=24, predicted 6/6 blow-ups in sign and 24–28% in growth rate,
        and F₃ (fewer violating cells, k_c 2.6× smaller) destabilises FIRST. The
        violation FRACTION does not order the instability; the critical WAVENUMBER does.
        `τ_stab` remains available as a safety factor for constant-τ̂ runs.
    THREE-FRAME RESULT (η̂=0.03, cut, σ_ko=0.01, N=24/32/40, difference protocol):
        γ_visc(F₁) = 2.35482e−3 → 2.44941e−3 → 2.48251e−3 /M⊙, observed order ≈3;
        frame spread of γ_visc 0.150% → 0.130% → 0.045% (shrinking: frame-invariant);
        numerical floor γ(η̂=0) = 3.27e−3 → 1.96e−3 → 1.30e−3, ∝ dx^1.8, i.e. 52% of
        γ_visc at N=40 — the difference protocol is MANDATORY, not optional;
        a₁ is inert (a₁=6.25/25/100 at fixed a₂: f spread 0.006%, γ spread 0.22%).
        At σ_ko=0.005 the pencil γ_visc is 2.597e−3 (N=40) / 2.607e−3 (N=48) = 99.2% of
        the NS dissipation integral 0.08746·η̂ = 2.6238e−3; the "96–97%" implied by the
        σ_ko=0.01 Richardson sequence was itself KO-biased.
    THE VISCOUS FREQUENCY SHIFT — CONVERGED (repro/bdnk3d_viscous_shift.jl, 33 runs).
        The +0.40/+0.28/+0.15% (N=24/32/40) at σ_ko=0.01 was a KREISS–OLIGER ARTIFACT:
        halving σ_ko removes ~0.35% of it at every N. What remains converges to the
        damped-oscillator shift, which has no free parameter,
                 Δf/f = −(γ_tot² − γ₀²)/(2ω²)
        (σ_ko=0.005, η̂=0.03; measured mean of the two estimators / predicted):
             N=24 +0.05% / −0.20%   N=32 −0.03% / −0.16%   N=40 −0.08% / −0.13%
             N=48 −0.125% / −0.133%   — residual +0.25→+0.13→+0.05→+0.01%, ~h³.
        η̂-scaling at N=40 is quadratic once the offset is removed; frame-reactive
        corrections are O(τ²ω²) ≈ 3e−5. So Δf/f(η̂=0.03) = −0.125 ± 0.013%.
        Protocol facts: hold the frame FIXED (η̂_frame) while sweeping η̂; the inviscid
        control is unstable below σ_ko≈0.005 (N=24, 32) while viscous runs survive to
        0.00125, so σ_ko=0.005 is the floor of the difference protocol; `visc_compact`
        shifts f by +0.15% and cuts γ_visc 32% at N=24 — do not use it for the mode.
    STILL OPEN: realistic crusted EOS excite a growing surface mode for N≥40 (surface
        index n>1.8 in ε∝(R−r)ⁿ; mechanism open); σ_ko is a ±0.25% systematic on
        absolute frequencies and THE dominant systematic of any viscous frequency shift;
        time-domain resolution Δf=1/T is 4.4% of f at the production T — every
        sub-percent number is a two-estimator parametric estimate, never a resolved line.
=#
module CowlingBDNK3D

using ..Background3D: Star3D, build_star3d
using ..EquationOfState: BarotropicEOS
using ..Units: Msun_to_km, kHz_to_km
using ..CowlingEvolve3D: Scratch, _cutcell_geometry, l2_quadrupole, central_deps, periodogram,
                         freq_kHz_cyclic, damping_rate, ylm_real, ylm_moment, analyze_qnm

export BDNKState, BDNK3D, setup_bdnk3d, seed_bdnk_l2!, evolve_bdnk3d!,
       bdnk_causal_denominator, bdnk_bound_violation,
       seed_bdnk_ylm!, evolve_bdnk3d_moments!, bdnk3d_qnm

# 8-field state: conserved (δE, δS_i) + primitive (δε, δv^i)
mutable struct BDNKState
    δE::Array{Float64,3}
    δSx::Array{Float64,3}; δSy::Array{Float64,3}; δSz::Array{Float64,3}
    δε::Array{Float64,3}
    δvx::Array{Float64,3}; δvy::Array{Float64,3}; δvz::Array{Float64,3}
end
BDNKState(N::Int) = BDNKState((zeros(N,N,N) for _ in 1:8)...)

struct BDNK3D
    s::Star3D
    w0::Array{Float64,3}; cs2::Array{Float64,3}
    η0::Array{Float64,3}; ζ0::Array{Float64,3}
    τε::Array{Float64,3}; τP::Array{Float64,3}; τQ::Array{Float64,3}
    Φp::Array{Float64,3}; qinv::Array{Float64,3}
    # Christoffels of γ_ij = δ_ij+(A²−1)n_in_j are closed-form in n-structures:
    #     Γ^a_bc = n^a[ gΓ1 n_b n_c + gΓ2 (δ_bc − n_b n_c) ],  gΓ1=A′/A, gΓ2=(A²−1)/(rA²)
    # (exact in all 27 components, t186_christ2.m).  Note gΓ^i_im = gΓ1 n_m, so the
    # covariant expansion is θ = ∂_iv^i + gΓ1 (n·v) = (1/√γ)∂_i(√γ v^i).
    gΓ1::Array{Float64,3}; gΓ2::Array{Float64,3}
    wfloor::Float64        # excise the outermost shell w₀<wfloor (graded surface excision)
    den_floor::Float64     # ABSOLUTE floor on the recovery denominator τ_Q w₀−η (legacy)
    den_rel::Float64       # RELATIVE floor, as a fraction of τ_Q w₀; >0 replaces den_floor.
                           # ∂_tδv = α(δS−w₀δv−…)/den relaxes δv at rate α·w₀/den, so an
                           # ABSOLUTE floor destroys that rate wherever w₀ is small — at
                           # N=64 the worst active cell had α·w₀/den = 9.5e−8 against the
                           # true 1/τ_Q = 0.5, i.e. δv became an undamped integrator driven
                           # by δS, in exactly the outer shells where the unstable
                           # eigenvector lives. Flooring against τ_Q w₀ instead keeps the
                           # rate in [α/τ_Q, α/(den_rel τ_Q)] everywhere, so δv always
                           # relaxes toward the conserved-implied δS/w₀ rather than
                           # integrating freely. No absolute floor is needed for
                           # well-posedness: η₀∝ε₀ ⇒ τ_Q w₀−η = ε₀[τ̂(1+p₀/ε₀)−η̂] > 0.
                           # NOTE it is nearly always INACTIVE and that is intended: for
                           # τ̂>η̂ one has τ_Q w₀−η ≈ (τ̂−η̂)ε₀ > den_rel·τ̂·ε₀, so it binds in
                           # 0% of cells and simply disables the absolute floor. It is a
                           # guard, not the fix — `wcut` below is the fix.
    wcut::Float64          # NO FLUID ⇒ NO VELOCITY. Below this w₀ the velocity primitive
                           # is not recovered at all (∂_tδv≡0, so δv stays 0). Rewriting
                           # the recovery as ∂_tδv = (α w₀/den)(δS/w₀ − δv − src/w₀) shows
                           # the driving TARGET is δS/w₀, which is unbounded as w₀→0: the
                           # active set (κ>0) includes an ATMOSPHERE-FLOOR shell at
                           # w₀/max≈9e−10 (14.1% of active cells at N=32, 5.9% at N=88),
                           # where δS/w₀ amplifies noise by ~1e9. The w₀ histogram has a
                           # clean GAP — the count is identical for thresholds 1e−9…1e−3 —
                           # so the cutoff is unambiguous. Same treatment, and same 1e−6
                           # default, as the velocity reconstruction in CowlingEvolve3D.
    hrsc_a::Float64        # Kurganov–Tadmor signal-speed floor; >0 ⇒ conservative HRSC dissipation
    hyb_lo::Float64        # hybrid surface blend: w₀ below this ⇒ pure NS-limit δv→δS/w₀
    hyb_hi::Float64        # w₀ above this ⇒ pure BDNK recovery; smoothstep between
    κ_c::Float64           # constraint-damping rate: slave primitives → conserved-implied
    cut::Bool              # cut-cell (aperture-weighted) surface treatment
    κmin::Float64          # small-cell floor on the volume fraction
    κ::Array{Float64,3}
    axm::Array{Float64,3}; axp::Array{Float64,3}
    aym::Array{Float64,3}; ayp::Array{Float64,3}
    azm::Array{Float64,3}; azp::Array{Float64,3}
    bc_w::Float64          # EXCISION BC: shell w₀<bc_w where primitives are slaved exactly
    κ_bc::Float64          # BC relaxation rate for any pre-existing primitive/conserved offset
    σ_ko::Float64
    visc_compact::Bool     # COMPACT viscous stencil (see _ko4x/_ko4y/_ko4z above): removes
                           # the odd/even decoupling of the 2-cell Laplacian that the
                           # centred-stress-then-centred-divergence composition produces
end

"""
    setup_bdnk3d(s; η̂=0.0, ζ̂=0.0, τ̂=2.0, σ_ko=0.01, den_frac=0.004)

Build the full-frame BDNK coefficients. η₀=η̂·ε₀, ζ₀=ζ̂·ε₀, single frame time
τ_ε=τ_P=τ_Q=τ̂ (geometric units; τ̂·w₀>η₀ for causality, see
`bdnk_causal_denominator`).

SURFACE TREATMENT (two regularizers of the velocity recovery, which is stiff
where τ_Q w₀−η → 0 at the tenuous surface):
  • `den_frac` — smooth floor on the recovery denominator (fraction of central);
  • `excise_frac` — graded excision: evolve only where w₀ ≥ excise_frac·max(w₀).

`hrsc_a>0` enables conservative Kurganov–Tadmor/Rusanov numerical viscosity on the
conserved fluxes (signal speed max(c_s,hrsc_a)). `hyb_lo_frac<hyb_hi_frac` enables
the HYBRID surface blend (BDNK recovery in the bulk → NS δv→δS/w₀ slaving in the
surface shell, smoothstep in w₀).

SURFACE TREATMENT — SOLVED BY CUT CELLS (`cut=true`, default false).
Set `cut=true` and `excise_frac=0.0` to evolve the FULL star with no excision. The
domain is bounded at the TRUE r=R: cell volume fractions κ and the six face
apertures come from the level set (both reduce to 1D quadratures of one
disc–rectangle area, verified against the analytic sphere — volume to ~1e−7, and
the closure condition reproduces 4πR² at 2nd order). All three divergences are
aperture-weighted — the continuity flux, the isotropic pressure gradient, and the
viscous divergence of √γπ^{ik} — and every boundary-face term vanishes for a free
surface (zero mass flux, δp_b=0, zero viscous traction), so the BC is imposed BY
OMISSION rather than as an extra constraint. `κmin` floors κ in small cells
(κ-clipping; not conservative — cell merging would be the proper refinement, worth
~±0.2% on the ideal engine).
MEASURED (τ̂=2, σ_ko=0.004, matrix pencil; target = dissipation integral 0.08746;
slope from the η̂=0→0.05 interval, which subtracts the KO floor):
    N=32  cut, no excision   f=−0.699%  dγ/dη̂=0.09287 = 106.2%   floor 5.696e−4
    N=44  cut, no excision   f=−0.343%  dγ/dη̂=0.09170 = 104.8%   floor 3.198e−4
    N=64  cut, no excision   f=−0.192%  dγ/dη̂=0.08985 = 102.7%   floor 1.324e−4
    N=44  mask + excise      f=+8.18%   dγ/dη̂=0.03418 =  39.1%   floor 1.1e−2
Cut cells converge (frequency at observed order ~1.5–2.2, floor like dx^2.4, damping
ratio descending monotonically toward 100%); the masked scheme DIVERGES (86.4%→39.1%).
CAVEATS: η̂=0.20 returns NaN, so the sweep uses η̂∈{0,0.05,0.10}; and the target is the
NEWTONIAN-FORM integral, while this engine's stress is the full curved-space one — but
the monotone descent 106.2→104.8→102.7% now argues the Newtonian form IS the right
reference, superseding an earlier reading that the measured γ_rel/γ_newt=0.95 put the
truth near 108%.

*** OPEN: STABILITY CEILING AT N≈64 — THE SURFACE VELOCITY-DRIFT MODE IS BACK. ***
Envelope verdict (last window / min window), T=500, dt=0.20dx; rows N, columns η̂:
    N          0.0           0.05            0.10
     32      clean 1.2     clean 1.0       clean 1.0
     44      clean 1.7     clean 1.0       clean 1.0
     64      clean 1.2     GROWTH 2.1e1    BLOW-UP 2.1e14
     88      GROWTH 8.4    BLOW-UP 6.3e8   BLOW-UP 1.5e39
By N=88 even the INVISCID run grows. The N=32/44 agreement was a window before the
mode emerged — do NOT read it as convergence. The damping numbers above survive
because the matrix pencil separates the growing mode from the f-mode (at N=64, η̂=0.05
the damping is 4.656e−3 on the clean first 55% vs 4.625e−3 on the full record, 0.7%).
IT IS A SPATIAL-OPERATOR EIGENVALUE, NOT A TIME-STEP CONSTRAINT (N=64, η̂=0.10, T=250):
    cfl      0.20      0.10      0.05
    rate   0.08881   0.08878   0.08878
identical envelopes window-by-window over a 4× range in dt — so CFL, stiffness, any
parabolic viscous limit, and RK4 itself are all excluded, and refining dt cannot help.
LOCALISED by matrix-free power iteration on the RK4 step map (rate 0.0885/0.0901 from
two seeds, vs 0.08878 in the time domain). The eigenvector is equipartitioned PURE
VELOCITY — δvx,δvy,δvz = 0.3363, 0.3354, 0.3283 of the power, with δE=3.9e−6,
δε=8.4e−7, δS≈1.8e−9 small but NOT zero — and 97.6% of it sits in the outer two
shells (65.1% at r/R 0.92–1.00, 32.5% at 0.83–0.92), per-cell density climbing
monotonically outward by EIGHT orders. Cut cells are 19.5% of the active set but
carry 52.1% of the power.
*** ROOT CAUSE — IT IS THE EQUATIONS, NOT THE DISCRETISATION. *** Grouped Jacobian
block ablation (N=16, η̂=1.0, dim 10912, rhs(0)=0 exactly, max Re λ=+2.899e−2, purely
real, eigenvector an exactly equipartitioned δv triplet) shows FOUR couplings each
stabilise when removed — δε←δE (−1.73e−3), δv←δε (−1.72e−3), δS←δv (−5.33e−4),
δE←δS (−1.66e−4) — i.e. the closed cycle
    δv --(viscous stress)--> δS --(continuity)--> δE --(relaxation)--> δε --(HEAT FLUX)--> δv
(distinct from the historical Π-closure cycle, which bypassed δE). Solving the 4×4
sound sector for max Re s over all κ then settles it: THE BDN STABILITY CONDITION FOR
THIS FRAME IS
        τ_Q w₀ c_s²  >  (4/3) η        NOT  τ_Q w₀ > η
(measured τ̂_min/(η̂/c_s²) = 1.33–1.37 over c_s²=0.1335…0.001, η̂=0.02…0.30; the 4/3 is
ζ+4η/3). `bdnk_causal_denominator` checks only the latter, weaker by 1/c_s². Since
c_s²→0 at the surface while η₀ ∝ ε₀, it FAILS for ANY constant τ̂ and any η̂>0 — at
η̂=0.10 the constant τ̂=2 violates it in 63.5% OF ACTIVE CELLS — and Re s grows without
bound with κ, so a finer grid merely resolves more of the unstable band (Nyquist
κ = 4.4/8.7/12.0 at N=32/64/88). That is why N=32/44 looked clean, why four separate
DISCRETISATION fixes all failed, and why the COMPACT viscous stencil made it WORSE:
κ₂=2sin(kh/2)/h reaches 2/h at Nyquist where the 2-cell κ₁ falls back to 0, so the
inaccurate stencil was accidentally PROTECTIVE.
FIX: `τ_stab` (position-dependent relaxation, default 0.0 = off). With τ_stab=2.0 the
condition holds in EVERY cell with a live δv DOF (the residual violations are exactly
the atmosphere cells, where `wcut` already freezes δv — the two compose). Measured at
the DEFAULT σ_ko=0.004: N=64 η̂=0.10 growth 1.00 (was 2.1e14); N=88 η̂=0.05 growth 1.02
(was 6.3e8); and it is RESOLUTION-CONVERGED — η̂=0.05 gives f=+2.760%/+2.885% and
γ=3.862e−3/3.980e−3 at N=64/88.
BUT IT CHANGES THE PHYSICS, because τ is a FRAME choice and now spans ~7 orders of
magnitude across the star:
    η̂          0.0        0.05       0.10
    f err    −0.180%    +2.760%    +9.601%
    dγ/dη̂ = 0.07441 (85.1%) first-interval, 0.06091 (69.6%) 3-point, per-interval
    0.0744 → 0.0474 — strongly NON-LINEAR, where before they bracketed 0.08746.
TWO CAVEATS ON READING THAT. (i) The CONSTANT τ̂ was the artificial choice: τ ~
η/(w c_s²) is the natural relaxation scaling (essentially Israel–Stewart τ_π), and
constant τ̂ made the theory ill-posed at high k, so the earlier 102.7% was obtained
where the EQUATIONS were unstable. (ii) 0.08746 is the NEWTONIAN NS-LIMIT integral,
which BDNK approaches as τ→0; here τ is large, and the plane-wave work already
measured the causal correction pushing damping BELOW NS (14% low at k=1.2), so 85% may
be the correct CAUSAL answer — NOT established, it needs a BDNK-consistent reference.
ALTERNATIVE NOT TRIED: taper η ∝ ε₀c_s² so η/(w₀c_s²) stays bounded at constant τ;
that changes the viscosity profile instead of the frame, but moves the reference too.

HISTORICAL — WHAT FAILED BEFORE CUT CELLS (kept: it documents what does NOT work,
and every entry was working around the wrong boundary condition):
  (i)   denominator floor alone .................... unstable;
  (ii)  conservative HRSC dissipation (hrsc_a≤0.8) .. unstable;
  (iii) heavy 2nd-order dissipation on ALL fields ... unstable (⇒ NOT a high-k mode);
  (iv)  hybrid NS-surface slaving (full star) ....... unstable BUT reads f=1.886 kHz,
        i.e. the eigensolver value 1.883 to 0.2%;
  (v)   hybrid + dissipation / + minimal excision ... still unstable;
  (vi)  constraint damping κ_c (slave primitives → conserved) ... unstable (κ_c only
        slows growth and pulls f toward the NS value).
SUPERSEDED CONCLUSION (was: "the instability is ROBUST … needs a fundamentally
different method"). The second half of that was right — a boundary-conforming
treatment WAS required — but the diagnosis of the cause was not. It is not
ill-posedness as the causality margin τ_Q w₀−η → 0; it is that the staircased mask
imposes δv=0 as well as δp=0, i.e. a RIGID WALL, whose boundary POSITION is wrong
by O(dx) per direction. (The stair-step VOLUME error is irrelevant — 0.006% at
N=72 while the frequency error was 5.7%; it is the position, not the volume.)
Point (iv) was the tell and was read correctly at the time: the full-star physics
is fine, only the surface numerics fail.

EXCISION AND THE ABSOLUTE FREQUENCY (LEGACY PATH, cut=false). Excising the outer shell replaces the free
surface (δp=0) by a rigid wall (δv=0), which RAISES the fundamental; the error is
set by the cut radius and is essentially independent of τ̂ (2.0 vs 0.5 moves it
<1 point), so it is a cavity effect, not a frame effect. At η=0, N=24, against the
validated eigensolver (1.8829 kHz):
    excise_frac  0.08      0.04      0.02
    r_cut/R      0.887     0.942     0.963
    error       +39.1%    +12.0%     +2.6%
DEFAULT excise_frac=0.02: +2.6%, comparable to the ideal engine's own −5.7%
surface error, so at η=0 this engine reduces to validated ideal hydro. The older
0.08 was chosen when the run was still unstable (the Π-closure bug) and cost a
factor ~1.4 in the absolute frequency.
VERIFIED STABLE at the new default — dense Jacobian max Re λ = −8.79e−3 (N=14,
η̂=0), −1.17e−2 (N=14, η̂=0.1), −5.27e−3 (N=16, η̂=0), −8.01e−3 (N=16, η̂=0.1);
and 30000-step evolutions at N=24 stay bounded and decay for η̂=0 and 0.05,
reading +3.1% / +3.4%.
CAVEAT, N≥20 ONLY: the update loop runs 3:N-2, so an active cell with an index in
the 2-cell grid halo is seeded but never advanced (frozen Dirichlet data inside
the star). With Lfac=1.2 the star's surface sits ~0.083(N−1) cells from the grid
edge, so this bites only on very coarse grids: at N=14 24/720 cells (3.3%) are
frozen; at N=20,24,32,44 it is exactly 0 (Lfac 1.2 and 1.5 alike). The wider
excise_frac=0.08 masked this by cutting those cells anyway.
"""
function setup_bdnk3d(s::Star3D; η̂::Float64=0.0, ζ̂::Float64=0.0, τ̂::Float64=2.0,
                      frame::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
                      η̂_frame::Union{Nothing,Real}=nothing,
                      σ_ko::Float64=0.01, den_frac::Float64=0.004, den_rel::Float64=0.0,
                      wcut_frac::Float64=1e-6, excise_frac::Float64=0.0,
                      hrsc_a::Float64=0.0, hyb_lo_frac::Float64=0.0, hyb_hi_frac::Float64=0.0,
                      κ_c::Float64=0.0, bc_frac::Float64=0.0, κ_bc::Float64=1.0,
                      cut::Bool=true, κmin::Float64=0.5, visc_compact::Bool=false,
                      τ_stab::Float64=0.0, cs2_floor::Float64=1e-8)
    w0 = s.ε0 .+ s.p0
    Φp = similar(w0); qinv = similar(w0); gΓ1 = similar(w0); gΓ2 = similar(w0)
    @inbounds for I in eachindex(w0)
        r = s.r[I]; e2λ = s.e2λ[I]; m = r*(1 - 1/e2λ)/2
        denom = r*(r - 2m)
        Φp[I] = denom > 0 ? (m + 4π*r^3*s.p0[I])/denom : 0.0
        qinv[I] = 1/e2λ - 1                       # = −(A²−1)/A², so γ^{ij}=δ^{ij}+qinv nⁱnʲ
        # A = √e2λ = (1−2m/r)^{−1/2} with m′=4πr²ε₀  ⇒  A′ = A³(4πr³ε₀ − m)/r².
        # Both coefficients vanish at the centre (m∼r³ ⇒ A′∼r, (A²−1)/r ∼ r), but the
        # background m(r) is LINEARLY interpolated, so at an on-grid origin (odd N,
        # where r is clamped to 1e−12) m/r² is spuriously ∼1/r and gΓ1,gΓ2 blow up to
        # ∼1e11. Impose the analytic r→0 limit inside the first half-cell.
        if r > 0.25*s.grid.dx
            gΓ1[I] = e2λ*(4π*r^3*s.ε0[I] - m)/r^2                  # A′/A
            gΓ2[I] = -qinv[I]/r                                    # (A²−1)/(rA²) = 2m/r²
        else
            gΓ1[I] = 0.0; gΓ2[I] = 0.0
        end
    end
    η0 = η̂ .* s.ε0; ζ0 = ζ̂ .* s.ε0
    # POSITION-DEPENDENT RELAXATION TIME (τ_stab>0). The BDN stability condition for
    # this frame is NOT the τ_Q w₀ > η that `bdnk_causal_denominator` checks — solving
    # the 4×4 sound sector (δE,δS,δε,δv) for max Re s over all κ gives
    #        τ_Q w₀ c_s² > (4/3) η
    # (measured τ̂_min/(η̂/c_s²) = 1.33–1.37 over c_s²=0.1335…0.001, η̂=0.02…0.30; the
    # 4/3 is the longitudinal combination ζ+4η/3). That is weaker than the checked
    # condition by a factor 1/c_s², and since c_s²→0 at the surface while η₀ ∝ ε₀ it
    # FAILS in some surface layer for ANY constant τ̂ and any η̂>0 — the growth rate
    # rises without bound with κ, so a finer grid simply resolves more of the unstable
    # band (Nyquist κ = 4.4/8.7/12.0 at N=32/64/88). This is why N=32/44 looked clean
    # and N≥64 did not, and why four separate DISCRETISATION fixes all failed.
    # τ_stab is the safety factor C in τ = max(τ̂, C·η₀/(w₀c_s²)); C ≥ 4/3 is required
    # and C ≈ 2 leaves margin. As τ→∞ the surface recovers IDEAL behaviour
    # (∂_tδv → −αc_s²∇δε/w₀, ∂_tδε → −αw₀θ), while the viscous stress −2ησ is
    # independent of τ, so the dissipation itself is retained.
    τ = fill(τ̂, size(w0))
    if τ_stab > 0
        @inbounds for I in eachindex(τ)
            c2 = max(s.cs2[I], cs2_floor)
            w  = max(w0[I], eps(Float64))
            τ[I] = max(τ̂, τ_stab * η0[I] / (w * c2))
        end
    end
    # BDNK frame: A = dE - de is the energy-frame correction, so the matching
    # PRESSURE correction is fixed by the EOS, Pi = c_s^2 A, i.e. tau_P = c_s^2 tau_e.
    # Then  iso = c_s^2 de + (tau_P/tau_e)(dE - de) = c_s^2 dE  -- the isotropic
    # stress depends on the CONSERVED energy only, and d(iso)/d(de) = 0.
    # With tau_P = tau_e instead one gets iso = dE + (c_s^2-1) de, whose NEGATIVE
    # coefficient (c_s^2-1) closes a primitive-only feedback cycle
    #     dv --(-w0 theta)--> de --(stress)--> dS --(dS/den)--> dv
    # that the exact-Jacobian block ablation identifies as the instability:
    # zeroing ANY of those three blocks stabilises, all to the same -1.02e-2.
    # ── HYDRODYNAMIC FRAME (Clarisse et al. arXiv:2510.16603 Table 1) ────────────────
    # τ_ε = a₁ (1/T)(η/s), τ_Q = a₂ (1/T)(η/s); at zero chemical potential w = sT so
    # (1/T)(η/s) = η/w and the stellar map is τ_ε = a₁η/w₀(x), τ_Q = a₂η/w₀(x), installed
    # cell by cell. τ_P = c_s²τ_ε is NOT a free slot (see below). Because η₀ ∝ ε₀ the
    # times are nearly uniform across the star: for the anchor star F₁ gives
    # τ_ε ∈ [0.168, 0.188] M⊙, F₂ [0.337, 0.375], F₃ [0.674, 0.750].
    # The frame must be built from a POSITIVE viscosity even for the η̂=0 control run of
    # the difference protocol (γ_visc = γ(η̂) − γ(0) at identical frame) — hence η̂_frame.
    # Under this map the module's stability bound τ_Q w₀ c_s² > (4/3)η becomes the
    # η-independent a₂c_s² > 4/3; its violation is SUFFICIENT-NOT-NECESSARY (see header).
    τQ = copy(τ)
    if frame !== nothing
        a1, a2 = float(frame[1]), float(frame[2])
        ηf = float(η̂_frame === nothing ? η̂ : η̂_frame)
        ηf > 0 || throw(ArgumentError("frame=(a₁,a₂) builds τ = a·η/w₀ and needs a positive " *
                                      "viscosity: pass η̂_frame>0 for an η̂=0 control run"))
        a2 > 1 || @warn "frame: a₂ ≤ 1 makes the recovery denominator τ_Q w₀−η = (a₂−1)η ≤ 0 (ill-posed)"
        @inbounds for I in eachindex(τ)
            ηI = ηf * s.ε0[I]; w = max(w0[I], eps(Float64))
            τ[I]  = a1 * ηI / w
            τQ[I] = a2 * ηI / w
        end
    end
    τP = s.cs2 .* τ
    den_floor = den_frac * (frame === nothing ? τ̂ * maximum(w0) : maximum(τQ .* w0))
    wfloor = excise_frac * maximum(w0)         # excise the outermost unstable shell
    wmax = maximum(w0)
    N3 = s.grid.N
    z3()=zeros(N3,N3,N3)
    κc,axm,axp,aym,ayp,azm,azp = cut ? _cutcell_geometry(s) :
                                 (z3(),z3(),z3(),z3(),z3(),z3(),z3())
    BDNK3D(s, w0, copy(s.cs2), η0, ζ0, τ, τP, τQ, Φp, qinv, gΓ1, gΓ2, wfloor, den_floor,
           den_rel, wcut_frac*maximum(w0),
           hrsc_a, hyb_lo_frac*wmax, hyb_hi_frac*wmax, κ_c,
           cut, κmin, κc, axm, axp, aym, ayp, azm, azp,
           bc_frac*wmax, κ_bc, σ_ko, visc_compact)
end

"""min over the star of (τ_Q w₀ − η₀); >0 ⇒ the BDN causality condition (a) holds."""
function bdnk_causal_denominator(e::BDNK3D)
    m = Inf
    @inbounds for I in eachindex(e.w0)
        e.s.interior[I] || continue
        m = min(m, e.τQ[I]*e.w0[I] - e.η0[I])
    end
    m
end

"""
    bdnk_bound_violation(e; η̂_ref=nothing) -> (cells, volume)

Fraction of interior cells (and of proper volume √γ dV) violating the sound-sector
stability bound τ_Q w₀ c_s² > (4/3) η. With `η̂_ref` the bound is evaluated at that
reference viscosity instead of the engine's own η₀ (so an η̂=0 control engine can report
the fraction of the frame it carries). Under the Clarisse map the bound is a₂c_s² > 4/3,
independent of η. THIS FRACTION DOES NOT ORDER THE INSTABILITY: the unstable band lies
above k_c = C/η̂ and is only resolved for N > 353/232/135 (F₁/F₂/F₃, anchor star); the
critical wavenumber does.
"""
function bdnk_bound_violation(e::BDNK3D; η̂_ref::Union{Nothing,Real}=nothing)
    nv = 0; ntot = 0; vv = 0.0; vtot = 0.0
    @inbounds for I in eachindex(e.w0)
        e.s.interior[I] || continue
        η = η̂_ref === nothing ? e.η0[I] : float(η̂_ref)*e.s.ε0[I]
        viol = e.τQ[I]*e.w0[I]*e.cs2[I] ≤ (4/3)*η
        ntot += 1; vtot += e.s.sqrtγ[I]
        if viol; nv += 1; vv += e.s.sqrtγ[I]; end
    end
    (cells = ntot > 0 ? nv/ntot : NaN, volume = vtot > 0 ? vv/vtot : NaN)
end

@inline _dx(A,i,j,k,h)=(A[i+1,j,k]-A[i-1,j,k])/(2h)
@inline _dy(A,i,j,k,h)=(A[i,j+1,k]-A[i,j-1,k])/(2h)
@inline _dz(A,i,j,k,h)=(A[i,j,k+1]-A[i,j,k-1])/(2h)
@inline function _ko(A,i,j,k)
    (A[i-2,j,k]-4A[i-1,j,k]+6A[i,j,k]-4A[i+1,j,k]+A[i+2,j,k]) +
    (A[i,j-2,k]-4A[i,j-1,k]+6A[i,j,k]-4A[i,j+1,k]+A[i,j+2,k]) +
    (A[i,j,k-2]-4A[i,j,k-1]+6A[i,j,k]-4A[i,j,k+1]+A[i,j,k+2])
end
# PER-DIRECTION 4th differences, for the COMPACT VISCOUS STENCIL (`visc_compact`).
# The stress is built from CENTRED cell-centred velocity differences and its
# divergence telescopes to a centred difference too, so ∂_d(η∂_d v) is discretised as
# the 2-CELL Laplacian L₂ = (v[i+2]−2v[i]+v[i−2])/4h², symbol −sin²(kh)/h². That
# leaves the odd/even sublattices UNCOUPLED (it vanishes at kh=π) and is 47% too weak
# at the kh≈1.34 of the observed unstable mode. The compact L₁ = (v[i+1]−2v[i]+v[i−1])/h²
# has symbol −4sin²(kh/2)/h², and exactly
#     L₁ − L₂ = −(1/4h²)·KO4_d
# so adding −(coeff)·(L₁−L₂) per direction converts the operator to the compact one.
# Collecting the same-direction coefficients in ∂_iπ^{ij} (−η for i≠j, −4η/3 for i=j):
#     correction^j = (η/4h²)[ Σ_d KO4_d(v^j) + ⅓ KO4_j(v^j) ]
# CHECK (shear branch): −αη[sin²(kh) + 4sin⁴(kh/2)]/h² = −4αη sin²(kh/2)/h² = −αη κ₂²
# with κ₂ = 2sin(kh/2)/h — i.e. precisely the compact symbol. The correction is O(h²)
# consistent (a difference of two consistent discretisations) and ∝ η, so the inviscid
# limit is untouched.
# VARIABLE η IS ESSENTIAL. The constant-η form of the correction is (η/4h²)·KO4, but
# η₀ = η̂·ε₀ varies by ORDERS OF MAGNITUDE across the surface layer — exactly where the
# mode lives — and there the two discretisations differ at O(η/h²), i.e. the same order
# as the viscous term itself, so the constant-η correction is not a small O(h²) term at
# all. (Measured: using it blew N=64, η̂=0.05 up to 1.6e10.) These return the difference
# (compact − 2-cell) with FACE-AVERAGED η, reducing to (η/4h²)·KO4_d when η is constant.
@inline function _cmp_x(u,ηa,i,j,k,h2)
    ηp=0.5*(ηa[i,j,k]+ηa[i+1,j,k]); ηm=0.5*(ηa[i,j,k]+ηa[i-1,j,k])
    (ηp*(u[i+1,j,k]-u[i,j,k]) - ηm*(u[i,j,k]-u[i-1,j,k]))/h2 -
    (ηa[i+1,j,k]*(u[i+2,j,k]-u[i,j,k]) - ηa[i-1,j,k]*(u[i,j,k]-u[i-2,j,k]))/(4h2)
end
@inline function _cmp_y(u,ηa,i,j,k,h2)
    ηp=0.5*(ηa[i,j,k]+ηa[i,j+1,k]); ηm=0.5*(ηa[i,j,k]+ηa[i,j-1,k])
    (ηp*(u[i,j+1,k]-u[i,j,k]) - ηm*(u[i,j,k]-u[i,j-1,k]))/h2 -
    (ηa[i,j+1,k]*(u[i,j+2,k]-u[i,j,k]) - ηa[i,j-1,k]*(u[i,j,k]-u[i,j-2,k]))/(4h2)
end
@inline function _cmp_z(u,ηa,i,j,k,h2)
    ηp=0.5*(ηa[i,j,k]+ηa[i,j,k+1]); ηm=0.5*(ηa[i,j,k]+ηa[i,j,k-1])
    (ηp*(u[i,j,k+1]-u[i,j,k]) - ηm*(u[i,j,k]-u[i,j,k-1]))/h2 -
    (ηa[i,j,k+1]*(u[i,j,k+2]-u[i,j,k]) - ηa[i,j,k-1]*(u[i,j,k]-u[i,j,k-2]))/(4h2)
end
# 2nd-order Laplacian (Σ D²): the Kurganov–Tadmor/Rusanov numerical viscosity kernel
@inline function _lap(A,i,j,k)
    (A[i+1,j,k]-2A[i,j,k]+A[i-1,j,k]) +
    (A[i,j+1,k]-2A[i,j,k]+A[i,j-1,k]) +
    (A[i,j,k+1]-2A[i,j,k]+A[i,j,k-1])
end

# RHS for all 8 fields. scratch arrays: δp, contravariant δS^i (Sux..), δT^{ij} (πxx..)
function _rhs!(d::BDNKState, st::BDNKState, e::BDNK3D, scr::Scratch)
    s = e.s; N = s.grid.N; h = s.grid.dx; σ = e.σ_ko
    Iso = scr.δp; Sux,Suy,Suz = scr.Sux,scr.Suy,scr.Suz   # Iso ≡ isotropic stress c_s²δε+Π
    vx,vy,vz = scr.vx,scr.vy,scr.vz
    Txx,Tyy,Tzz,Txy,Txz,Tyz = scr.πxx,scr.πyy,scr.πzz,scr.πxy,scr.πxz,scr.πyz
    # contravariant conserved momentum δS^i and δp; primitive velocity already in st
    @inbounds for I in eachindex(Iso)
        nS = s.nx[I]*st.δSx[I] + s.ny[I]*st.δSy[I] + s.nz[I]*st.δSz[I]
        Sux[I] = st.δSx[I] + e.qinv[I]*s.nx[I]*nS
        Suy[I] = st.δSy[I] + e.qinv[I]*s.ny[I]*nS
        Suz[I] = st.δSz[I] + e.qinv[I]*s.nz[I]*nS
        vx[I] = st.δvx[I]; vy[I] = st.δvy[I]; vz[I] = st.δvz[I]
    end
    # ISOTROPIC stress  Iso = c_s²δε + Π  (Π = −ζθ + τ_P(δE−δε)/τ_ε) is stored on the
    # grid, and the T arrays now hold the ANISOTROPIC stress PRE-MULTIPLIED by √γ,
    # i.e. Txx = √γ π^{xx}, so the momentum loop can form (1/√γ)∂_i(√γ π^{ik}) directly.
    # The two sectors are treated differently because they are exact differently:
    # t175/t177 verified the isotropic momentum equation IS the bare-lapse flat form
    # ∂_tδS_j = −α[∂_j Iso + Φ′ n_j(Iso+δE)] on the full curved background (every metric
    # factor cancels); the anisotropic one is NOT — see the momentum loop.
    #
    # CURVED-SPACE SHEAR.  σ^{ik} = ½(∇^iv^k+∇^kv^i) − ⅓γ^{ik}θ with the closed-form
    # Christoffels Γ^a_bc = n^a[gΓ1 n_bn_c + gΓ2(δ_bc−n_bn_c)] (t186):
    #     ∇_m v^k = ∂_m v^k + n^k C_m ,   C_m = gΓ1 n_m(n·v) + gΓ2(v_m − n_m(n·v))
    #     ∇^i v^k = ∇_i v^k + qinv nⁱ (n^m∇_m v^k) ,  n^m C_m = gΓ1(n·v)
    # THE EXPANSION IS THE 4D ONE and carries the LAPSE as well as √γ (t187):
    #     θ = ∇_μu^μ = (1/(α√γ))∂_i(α√γ v^i) = ∂_iv^i + (Φ′+gΓ1)(n·v).
    # σ must be projected on BOTH indices, σ^{μν}=Δ^{μα}Δ^{νβ}∇_(α u_β) − ⅓Δ^{μν}θ.
    # Because Δ^{it}=v^i/α is FIRST ORDER and ∇_(k u_t) has the O(1) piece ½α′n_k,
    # that adds a term the purely spatial gradient misses (t188, exact, 9/9):
    #     σ^{ij} = ½(∇^iv^j+∇^jv^i) + (Φ′/2A²)(v^in^j+v^jn^i) − ⅓γ^{ij}θ.
    # CONTROL: γ_ij σ^{ij} = [∂_iv^i+gΓ1(n·v)] + Φ′(n·v) − θ = 0 — the Φ′ term is
    # exactly what makes σ trace-free against the TRUE expansion. 1/A² = qinv+1.
    @inbounds for k in 2:N-1, j in 2:N-1, i in 2:N-1
        nx=s.nx[i,j,k]; ny=s.ny[i,j,k]; nz=s.nz[i,j,k]
        g1=e.gΓ1[i,j,k]; g2=e.gΓ2[i,j,k]; qi=e.qinv[i,j,k]; Φ=e.Φp[i,j,k]
        dxvx=_dx(vx,i,j,k,h); dyvy=_dy(vy,i,j,k,h); dzvz=_dz(vz,i,j,k,h)
        dyvx=_dy(vx,i,j,k,h); dzvx=_dz(vx,i,j,k,h)
        dxvy=_dx(vy,i,j,k,h); dzvy=_dz(vy,i,j,k,h)
        dxvz=_dx(vz,i,j,k,h); dyvz=_dy(vz,i,j,k,h)
        ux=vx[i,j,k]; uy=vy[i,j,k]; uz=vz[i,j,k]
        nv = nx*ux + ny*uy + nz*uz
        θ  = dxvx + dyvy + dzvz + (Φ + g1)*nv           # ∇_μ u^μ, the true expansion
        Π  = -e.ζ0[i,j,k]*θ + e.τP[i,j,k]*(st.δE[i,j,k]-st.δε[i,j,k])/e.τε[i,j,k]
        Iso[i,j,k] = e.cs2[i,j,k]*st.δε[i,j,k] + Π
        η  = e.η0[i,j,k]; sg = s.sqrtγ[i,j,k]
        Cx = g1*nx*nv + g2*(ux - nx*nv)
        Cy = g1*ny*nv + g2*(uy - ny*nv)
        Cz = g1*nz*nv + g2*(uz - nz*nv)
        ndx = nx*dxvx + ny*dyvx + nz*dzvx + g1*nv*nx    # n^m ∇_m v^x
        ndy = nx*dxvy + ny*dyvy + nz*dzvy + g1*nv*ny
        ndz = nx*dxvz + ny*dyvz + nz*dzvz + g1*nv*nz
        Dxx = dxvx + nx*Cx + qi*nx*ndx                  # D^{ik} = ∇^i v^k
        Dxy = dxvy + ny*Cx + qi*nx*ndy
        Dxz = dxvz + nz*Cx + qi*nx*ndz
        Dyx = dyvx + nx*Cy + qi*ny*ndx
        Dyy = dyvy + ny*Cy + qi*ny*ndy
        Dyz = dyvz + nz*Cy + qi*ny*ndz
        Dzx = dzvx + nx*Cz + qi*nz*ndx
        Dzy = dzvy + ny*Cz + qi*nz*ndy
        Dzz = dzvz + nz*Cz + qi*nz*ndz
        # π^{ik} = −η(D^{ik}+D^{ki}) + (2η/3)θ γ^{ik} − (ηΦ′/A²)(v^in^k+v^kn^i)
        t23 = (2η/3)*θ; hΦ = η*Φ*(qi+1)                 # qi+1 = 1/A²
        Txx[i,j,k]=sg*(-2η*Dxx + t23*(1+qi*nx*nx) - 2hΦ*ux*nx)
        Tyy[i,j,k]=sg*(-2η*Dyy + t23*(1+qi*ny*ny) - 2hΦ*uy*ny)
        Tzz[i,j,k]=sg*(-2η*Dzz + t23*(1+qi*nz*nz) - 2hΦ*uz*nz)
        Txy[i,j,k]=sg*(-η*(Dxy+Dyx) + t23*qi*nx*ny - hΦ*(ux*ny+uy*nx))
        Txz[i,j,k]=sg*(-η*(Dxz+Dzx) + t23*qi*nx*nz - hΦ*(ux*nz+uz*nx))
        Tyz[i,j,k]=sg*(-η*(Dyz+Dzy) + t23*qi*ny*nz - hΦ*(uy*nz+uz*ny))
    end
    fill!(d.δE,0.0); fill!(d.δSx,0.0); fill!(d.δSy,0.0); fill!(d.δSz,0.0)
    fill!(d.δε,0.0); fill!(d.δvx,0.0); fill!(d.δvy,0.0); fill!(d.δvz,0.0)
    @inbounds for k in 3:N-2, j in 3:N-2, i in 3:N-2
        # CUT CELLS: evolve wherever the cell has volume inside r=R. Otherwise the
        # legacy mask (interior minus the graded excision shell).
        if e.cut
            e.κ[i,j,k] > 0.0 || continue
        else
            (s.interior[i,j,k] && e.w0[i,j,k] ≥ e.wfloor) || continue
        end
        κcc = e.cut ? max(e.κ[i,j,k], e.κmin) : 1.0
        Axm=e.axm[i,j,k]; Axp=e.axp[i,j,k]; Aym=e.aym[i,j,k]
        Ayp=e.ayp[i,j,k]; Azm=e.azm[i,j,k]; Azp=e.azp[i,j,k]
        α=s.α[i,j,k]; sg=s.sqrtγ[i,j,k]; Φ=e.Φp[i,j,k]
        nx=s.nx[i,j,k]; ny=s.ny[i,j,k]; nz=s.nz[i,j,k]
        w=e.w0[i,j,k]; cs2=e.cs2[i,j,k]; τε=e.τε[i,j,k]; τQ=e.τQ[i,j,k]; η=e.η0[i,j,k]
        # --- DIAGONAL BDNK recovery of the primitive time-derivatives ---------
        # The constitutive relations are COVARIANT: A = τ_ε[Dε + wθ] and
        # Q^ν = τ_Q[w Du^ν + c_s²Δ^{νl}∂_lε] with D = u^μ∂_μ, θ = ∇_μu^μ. On the
        # static background u^μ=(1/α,0,0,0), so the recovered derivative is a PROPER-
        # time one and converting to coordinate time multiplies by α (t188, exact):
        #     Dδε = (1/α)∂_tδε ,   (Du^i)_lin = (1/α)∂_tδv^i ,  Δ^{il}∂_lδε = γ^{il}∂_lδε.
        # This is the same missing-lapse error as the momentum equation (t175), here
        # in the recovery. θ is the full expansion ∂_iv^i + (Φ′+gΓ1)(n·v).
        nv = nx*vx[i,j,k]+ny*vy[i,j,k]+nz*vz[i,j,k]
        θ = _dx(vx,i,j,k,h)+_dy(vy,i,j,k,h)+_dz(vz,i,j,k,h) + (Φ+e.gΓ1[i,j,k])*nv
        dtε = α*((st.δE[i,j,k]-st.δε[i,j,k])/τε - w*θ)
        # HYBRID surface blend: W=1 in the bulk ⇒ full causal BDNK recovery;
        # W→0 in the surface shell ⇒ den→w₀, τ_Q-source→0, so ∂_tδvⁱ→(δSⁱ−w₀δvⁱ)/w₀
        # = δSⁱ/w₀−δvⁱ, i.e. δvⁱ relaxes to the regular NS velocity δS^i/w₀ (no stiff recovery).
        W = e.hyb_hi > e.hyb_lo ?
            (let t=clamp((w-e.hyb_lo)/(e.hyb_hi-e.hyb_lo),0.0,1.0); t*t*(3-2t) end) : 1.0
        # NO FLUID ⇒ NO VELOCITY. Written as ∂_tδv = (α w/den)(δS/w − δv − src/w), the
        # recovery drives δv toward δS/w — a TARGET that is unbounded as w→0. With cut
        # cells and excise_frac=0 the active set reaches the atmosphere-floor shell
        # (w/max≈9e−10), where δS/w amplifies noise by ~1e9 and δv runs away. Below
        # `wcut` the velocity primitive is simply not recovered (see the field comment).
        # For the legacy masked path this is inert: there w ≥ excise_frac·max(w₀).
        hasfluid = w ≥ e.wcut
        dtvx = 0.0; dtvy = 0.0; dtvz = 0.0
        if hasfluid
            # RELATIVE floor (den_rel>0) keeps the δv relaxation rate α·w/den bounded
            # below by α/τ_Q; the legacy ABSOLUTE floor drives that rate to ~0 wherever
            # w₀ is small, which is what let δv integrate freely in the outer shells.
            denraw = W*(τQ*w - η) + (1-W)*max(w, e.den_floor)
            den = e.den_rel > 0 ? max(denraw, e.den_rel*τQ*w, 1e-300) :
                                  max(denraw, e.den_floor)
            # heat-flux gradient raised with γ, not δ: γ^{il}∂_lδε = ∂_iδε + qinv nⁱ(n·∂δε)
            ex=_dx(st.δε,i,j,k,h); ey=_dy(st.δε,i,j,k,h); ez=_dz(st.δε,i,j,k,h)
            qi = e.qinv[i,j,k]; nde = qi*(nx*ex+ny*ey+nz*ez)
            dtvx = α*(Sux[i,j,k]-w*st.δvx[i,j,k]-W*τQ*cs2*(ex+nde*nx))/den
            dtvy = α*(Suy[i,j,k]-w*st.δvy[i,j,k]-W*τQ*cs2*(ey+nde*ny))/den
            dtvz = α*(Suz[i,j,k]-w*st.δvz[i,j,k]-W*τQ*cs2*(ez+nde*nz))/den
        end
        # CONSTRAINT DAMPING: slave the primitives to the conserved-implied values,
        # suppressing the drift between the separately-integrated 8 fields.
        if e.κ_c > 0
            ws = max(w, e.den_floor)
            dtε  += e.κ_c*(st.δE[i,j,k] - st.δε[i,j,k])
            dtvx += e.κ_c*(Sux[i,j,k]/ws - st.δvx[i,j,k])
            dtvy += e.κ_c*(Suy[i,j,k]/ws - st.δvy[i,j,k])
            dtvz += e.κ_c*(Suz[i,j,k]/ws - st.δvz[i,j,k])
        end
        # --- conserved-density evolution  ∂_tδE=−∇·δS^i,  ∂_tδS_j=−∂_iδT^{ij} ---
        Fx(a,b,c)=s.α[a,b,c]*s.sqrtγ[a,b,c]*Sux[a,b,c]
        Fy(a,b,c)=s.α[a,b,c]*s.sqrtγ[a,b,c]*Suy[a,b,c]
        Fz(a,b,c)=s.α[a,b,c]*s.sqrtγ[a,b,c]*Suz[a,b,c]
        # aperture-weighted FV divergence; the boundary face carries ZERO mass flux
        # (free surface moves with the fluid) so it never enters the sum.
        divS = e.cut ?
            ( Axp*0.5*(Fx(i,j,k)+Fx(i+1,j,k)) - Axm*0.5*(Fx(i-1,j,k)+Fx(i,j,k))
            + Ayp*0.5*(Fy(i,j,k)+Fy(i,j+1,k)) - Aym*0.5*(Fy(i,j-1,k)+Fy(i,j,k))
            + Azp*0.5*(Fz(i,j,k)+Fz(i,j,k+1)) - Azm*0.5*(Fz(i,j,k-1)+Fz(i,j,k)) )/(κcc*h) :
            (Fx(i+1,j,k)-Fx(i-1,j,k)+Fy(i,j+1,k)-Fy(i,j-1,k)+Fz(i,j,k+1)-Fz(i,j,k-1))/(2h)
        ndotSup=nx*Sux[i,j,k]+ny*Suy[i,j,k]+nz*Suz[i,j,k]
        d.δE[i,j,k] = -divS/sg - α*Φ*ndotSup - σ*_ko(st.δE,i,j,k)
        # --- ISOTROPIC sector: exact as it stands (t175/t177) ------------------
        iso = Iso[i,j,k]
        grav = Φ*(iso + st.δE[i,j,k])
        # --- ANISOTROPIC (viscous) sector, FULL CURVED-SPACE FORM --------------
        # The momentum flux is −(1/√γ)∂_i(α√γ γ_jk δT^{ik}) + α·½δT^{ik}∂_jγ_ik.
        # Metric compatibility gives ½π^{ik}∂_jγ_ik = Γ^m_ji π^i_m exactly, so the
        # whole thing collapses to the manifestly covariant
        #     ∂_tδS_j|visc = −α ∇_i π^i_j − π^i_j ∂_iα ,   ∂_iα = αΦ′n_i.
        # With Γ^k_im π^{im} = n^k[gΓ1(nnπ) + gΓ2(tr_δπ − nnπ)] (t186) every piece is
        # a local contraction, so no extra grid arrays are needed:
        #     ∇_iπ^{ik} = (1/√γ)∂_i(√γ π^{ik}) + Γ^k_im π^{im}       [T holds √γ π]
        #     γ_jk ∇_iπ^{ik} = (∇·π)_j + (A²−1) n_j n·(∇·π)
        #     π^i_j n_i      = q_j + (A²−1) n_j (nnπ),   q^k ≡ n_iπ^{ik}
        # A→1 reduces this EXACTLY to the flat form −α∂_iπ^{ij} − αΦ′(π^{ij}n_i).
        sgi = 1/sg; Bm = s.e2λ[i,j,k] - 1                      # Bm = A²−1
        g1=e.gΓ1[i,j,k]; g2=e.gΓ2[i,j,k]
        pxx=Txx[i,j,k]*sgi; pyy=Tyy[i,j,k]*sgi; pzz=Tzz[i,j,k]*sgi
        pxy=Txy[i,j,k]*sgi; pxz=Txz[i,j,k]*sgi; pyz=Tyz[i,j,k]*sgi
        qx = nx*pxx + ny*pxy + nz*pxz
        qy = nx*pxy + ny*pyy + nz*pyz
        qz = nx*pxz + ny*pyz + nz*pzz
        nnπ = nx*qx + ny*qy + nz*qz
        gcon = g1*nnπ + g2*((pxx+pyy+pzz) - nnπ)               # Γ^k_im π^{im} = n^k gcon
        # same aperture weighting for the viscous stress divergence; the boundary
        # face contributes π^{ij}n_j A_b = 0 (no viscous traction on a free surface).
        fvx = e.cut ?
            ( Axp*0.5*(Txx[i,j,k]+Txx[i+1,j,k]) - Axm*0.5*(Txx[i-1,j,k]+Txx[i,j,k])
            + Ayp*0.5*(Txy[i,j,k]+Txy[i,j+1,k]) - Aym*0.5*(Txy[i,j-1,k]+Txy[i,j,k])
            + Azp*0.5*(Txz[i,j,k]+Txz[i,j,k+1]) - Azm*0.5*(Txz[i,j,k-1]+Txz[i,j,k]) )/(κcc*h) :
            (_dx(Txx,i,j,k,h)+_dy(Txy,i,j,k,h)+_dz(Txz,i,j,k,h))
        fvy = e.cut ?
            ( Axp*0.5*(Txy[i,j,k]+Txy[i+1,j,k]) - Axm*0.5*(Txy[i-1,j,k]+Txy[i,j,k])
            + Ayp*0.5*(Tyy[i,j,k]+Tyy[i,j+1,k]) - Aym*0.5*(Tyy[i,j-1,k]+Tyy[i,j,k])
            + Azp*0.5*(Tyz[i,j,k]+Tyz[i,j,k+1]) - Azm*0.5*(Tyz[i,j,k-1]+Tyz[i,j,k]) )/(κcc*h) :
            (_dx(Txy,i,j,k,h)+_dy(Tyy,i,j,k,h)+_dz(Tyz,i,j,k,h))
        fvz = e.cut ?
            ( Axp*0.5*(Txz[i,j,k]+Txz[i+1,j,k]) - Axm*0.5*(Txz[i-1,j,k]+Txz[i,j,k])
            + Ayp*0.5*(Tyz[i,j,k]+Tyz[i,j+1,k]) - Aym*0.5*(Tyz[i,j-1,k]+Tyz[i,j,k])
            + Azp*0.5*(Tzz[i,j,k]+Tzz[i,j,k+1]) - Azm*0.5*(Tzz[i,j,k-1]+Tzz[i,j,k]) )/(κcc*h) :
            (_dx(Txz,i,j,k,h)+_dy(Tyz,i,j,k,h)+_dz(Tzz,i,j,k,h))
        Dvx = fvx*sgi + nx*gcon
        Dvy = fvy*sgi + ny*gcon
        Dvz = fvz*sgi + nz*gcon
        if e.visc_compact
            # correction^j = (η/4h²)[ Σ_d KO4_d(v^j) + ⅓ KO4_j(v^j) ]  (derivation above).
            # Skipped where the ±2 stencil would reach a cell outside the domain, so this
            # cannot introduce a new boundary artifact.
            h2 = h*h; ηa = e.η0
            okx = !e.cut || (e.κ[i-2,j,k]>0 && e.κ[i-1,j,k]>0 && e.κ[i+1,j,k]>0 && e.κ[i+2,j,k]>0)
            oky = !e.cut || (e.κ[i,j-2,k]>0 && e.κ[i,j-1,k]>0 && e.κ[i,j+1,k]>0 && e.κ[i,j+2,k]>0)
            okz = !e.cut || (e.κ[i,j,k-2]>0 && e.κ[i,j,k-1]>0 && e.κ[i,j,k+1]>0 && e.κ[i,j,k+2]>0)
            # correction^j = −[ Σ_d c_dj (Lc−L2)_d(v^j) ],  c = 4/3 when d=j else 1
            if okx
                Dvx -= (4/3)*_cmp_x(vx,ηa,i,j,k,h2)
                Dvy -=       _cmp_x(vy,ηa,i,j,k,h2)
                Dvz -=       _cmp_x(vz,ηa,i,j,k,h2)
            end
            if oky
                Dvx -=       _cmp_y(vx,ηa,i,j,k,h2)
                Dvy -= (4/3)*_cmp_y(vy,ηa,i,j,k,h2)
                Dvz -=       _cmp_y(vz,ηa,i,j,k,h2)
            end
            if okz
                Dvx -=       _cmp_z(vx,ηa,i,j,k,h2)
                Dvy -=       _cmp_z(vy,ηa,i,j,k,h2)
                Dvz -= (4/3)*_cmp_z(vz,ηa,i,j,k,h2)
            end
        end
        nD = nx*Dvx + ny*Dvy + nz*Dvz
        vsx = -α*((Dvx + Bm*nx*nD) + Φ*(qx + Bm*nx*nnπ))
        vsy = -α*((Dvy + Bm*ny*nD) + Φ*(qy + Bm*ny*nnπ))
        vsz = -α*((Dvz + Bm*nz*nD) + Φ*(qz + Bm*nz*nnπ))
        # INT d_j Iso dV = OINT Iso n_j dA; the boundary face contributes Iso_b = 0
        # (free surface: Δp=0 with p₀′(R)=0 ⇒ δp(R)=0), so it drops out.
        gIx = e.cut ? (Axp*0.5*(Iso[i,j,k]+Iso[i+1,j,k]) - Axm*0.5*(Iso[i-1,j,k]+Iso[i,j,k]))/(κcc*h) : _dx(Iso,i,j,k,h)
        gIy = e.cut ? (Ayp*0.5*(Iso[i,j,k]+Iso[i,j+1,k]) - Aym*0.5*(Iso[i,j-1,k]+Iso[i,j,k]))/(κcc*h) : _dy(Iso,i,j,k,h)
        gIz = e.cut ? (Azp*0.5*(Iso[i,j,k]+Iso[i,j,k+1]) - Azm*0.5*(Iso[i,j,k-1]+Iso[i,j,k]))/(κcc*h) : _dz(Iso,i,j,k,h)
        d.δSx[i,j,k] = α*(-gIx - grav*nx) + vsx - σ*_ko(st.δSx,i,j,k)
        d.δSy[i,j,k] = α*(-gIy - grav*ny) + vsy - σ*_ko(st.δSy,i,j,k)
        d.δSz[i,j,k] = α*(-gIz - grav*nz) + vsz - σ*_ko(st.δSz,i,j,k)
        # HRSC: Kurganov–Tadmor/Rusanov numerical viscosity on the CONSERVED
        # fluxes, +½a·ΣD²(q)/h, signal speed a=max(c_s,hrsc_a) floored so the
        # surface (c_s→0) still gets dissipation — damps grid/surface modes (∝k²)
        # far more than the long-wavelength f-mode.
        if e.hrsc_a > 0
            visc = 0.5 * max(sqrt(cs2), e.hrsc_a) / h
            d.δE[i,j,k]  += visc*_lap(st.δE,i,j,k)
            d.δSx[i,j,k] += visc*_lap(st.δSx,i,j,k)
            d.δSy[i,j,k] += visc*_lap(st.δSy,i,j,k)
            d.δSz[i,j,k] += visc*_lap(st.δSz,i,j,k)
        end
        # --- EXCISION BOUNDARY CONDITION ---------------------------------------
        # The excision at w₀=wfloor is an ARTIFICIAL internal boundary, and the
        # primitives (δε,δvⁱ) have no condition there: they are free to drift from
        # the conserved fields. The exact Jacobian shows the resulting unstable
        # eigenvector is 100% δvⁱ with 0.00% in every conserved field, 91% of its
        # power in the two shells at that boundary, and growing with resolution
        # (+5.9e-3 at N=14 → +3.98e-2 at N=36). κ_c damps it in the BULK and is
        # NOT enough at production N (still +8e-3 at κ_c=5).
        # In a shell w₀<bc_w we therefore REMOVE the free primitive degree of
        # freedom instead of damping it: the primitives are advanced by the
        # conserved-implied derivatives (δε↦δE, δvⁱ↦δSⁱ/w₀), plus a relaxation
        # that eats any pre-existing offset. Reduces to the untouched scheme when
        # bc_frac=0.
        if e.bc_w > 0 && w < e.bc_w
            wb = max(w, e.den_floor)
            dtε  = d.δE[i,j,k]      + e.κ_bc*(st.δE[i,j,k]  - st.δε[i,j,k])
            dtvx = d.δSx[i,j,k]/wb  + e.κ_bc*(Sux[i,j,k]/wb - st.δvx[i,j,k])
            dtvy = d.δSy[i,j,k]/wb  + e.κ_bc*(Suy[i,j,k]/wb - st.δvy[i,j,k])
            dtvz = d.δSz[i,j,k]/wb  + e.κ_bc*(Suz[i,j,k]/wb - st.δvz[i,j,k])
        end
        # --- primitives advanced by the recovered time-derivatives ------------
        d.δε[i,j,k]  = dtε  - σ*_ko(st.δε,i,j,k)
        # KO is bypassed in the no-fluid shell too: its stencil reaches into the fluid,
        # so leaving it on would re-seed δv there from the neighbours it is meant to
        # stay decoupled from, and δv would not remain identically zero.
        d.δvx[i,j,k] = hasfluid ? dtvx - σ*_ko(st.δvx,i,j,k) : 0.0
        d.δvy[i,j,k] = hasfluid ? dtvy - σ*_ko(st.δvy,i,j,k) : 0.0
        d.δvz[i,j,k] = hasfluid ? dtvz - σ*_ko(st.δvz,i,j,k) : 0.0
    end
    return nothing
end

_fields(st::BDNKState) = (st.δE,st.δSx,st.δSy,st.δSz,st.δε,st.δvx,st.δvy,st.δvz)

function _rk4!(st::BDNKState, e::BDNK3D, dt, k1,k2,k3,k4, tmp, scr)
    F = _fields(st); T = _fields(tmp)
    _rhs!(k1, st, e, scr)
    K=_fields(k1); for n in 1:8; @. T[n] = F[n] + dt/2*K[n]; end
    _rhs!(k2, tmp, e, scr); K=_fields(k2); for n in 1:8; @. T[n] = F[n] + dt/2*K[n]; end
    _rhs!(k3, tmp, e, scr); K=_fields(k3); for n in 1:8; @. T[n] = F[n] + dt*K[n]; end
    _rhs!(k4, tmp, e, scr)
    K1=_fields(k1);K2=_fields(k2);K3=_fields(k3);K4=_fields(k4)
    for n in 1:8; @. F[n] += dt/6*(K1[n]+2K2[n]+2K3[n]+K4[n]); end
end

"""seed ℓ=2 density perturbation in BOTH δε and the conserved δE (δE≈δε at t=0)."""
function seed_bdnk_l2!(st::BDNKState, e::BDNK3D; A::Float64=1e-3)
    s = e.s; xs = s.grid.x; N = s.grid.N
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        (s.interior[i,j,k] && e.w0[i,j,k] ≥ e.wfloor) || continue   # seed evolved region
        X,Y,Z=xs[i],xs[j],xs[k]; r=s.r[i,j,k]
        v = A*s.ε0[i,j,k]*(r/s.R)*(2Z^2-X^2-Y^2)/r^2
        st.δε[i,j,k] = v; st.δE[i,j,k] = v        # at rest δE = δε + τ_ε w₀θ ≈ δε (θ=0)
    end
end

"""
    evolve_bdnk3d!(st, e; dt, nsteps, sample=1) -> (ts, q_l2, q_c)
"""
function evolve_bdnk3d!(st::BDNKState, e::BDNK3D; dt::Float64, nsteps::Int, sample::Int=1)
    N = e.s.grid.N
    k1=BDNKState(N); k2=BDNKState(N); k3=BDNKState(N); k4=BDNKState(N); tmp=BDNKState(N)
    scr = Scratch(N)
    ts=Float64[]; q2=Float64[]; qc=Float64[]
    xs=e.s.grid.x
    quad(st) = begin acc=0.0
        @inbounds for kk in 1:N, jj in 1:N, ii in 1:N
            e.s.interior[ii,jj,kk] || continue
            X,Y,Z=xs[ii],xs[jj],xs[kk]; r=e.s.r[ii,jj,kk]
            acc += st.δε[ii,jj,kk]*(2Z^2-X^2-Y^2)/r^2*e.s.sqrtγ[ii,jj,kk]
        end; acc*e.s.grid.dx^3
    end
    c = N÷2 + 1
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt); push!(q2, quad(st)); push!(qc, st.δε[c,c,c])
        end
        n == nsteps && break
        _rk4!(st, e, dt, k1,k2,k3,k4, tmp, scr)
    end
    return ts, q2, qc
end

# ═══════════════════════════════════════════════════════════════════════════════════
#  GENERAL (ℓ,m) SEEDING, MOMENT RECORDING, AND THE QNM DRIVER
# ═══════════════════════════════════════════════════════════════════════════════════
@inline _active_b(e::BDNK3D, i, j, k) =
    (e.cut ? (e.κ[i,j,k] > 0.0) : e.s.interior[i,j,k]) && e.w0[i,j,k] ≥ e.wfloor

"""
    seed_bdnk_ylm!(st, e; l=2, m=0, A=1e-3)

Seed δε = δE = A ε₀ (r/R)^ℓ Y_ℓm(n̂) (at rest δE = δε + τ_ε w₀θ = δε since θ=0), all
momenta and velocities zero. Any ℓ∈{2,3,4}, |m|≤ℓ; the full cube carries every m.
"""
function seed_bdnk_ylm!(st::BDNKState, e::BDNK3D; l::Int=2, m::Int=0, A::Float64=1e-3)
    s = e.s; N = s.grid.N
    @inbounds for k in 1:N, j in 1:N, i in 1:N
        _active_b(e, i, j, k) || continue
        r = s.r[i,j,k]
        v = A * s.ε0[i,j,k] * (r/s.R)^l * ylm_real(l, m, s.nx[i,j,k], s.ny[i,j,k], s.nz[i,j,k])
        st.δε[i,j,k] = v; st.δE[i,j,k] = v
    end
    for A_ in (st.δSx, st.δSy, st.δSz, st.δvx, st.δvy, st.δvz); fill!(A_, 0.0); end
    return nothing
end

"""
    evolve_bdnk3d_moments!(st, e; dt, nsteps, sample=1, moments=[(2,0)]) -> (ts, Q)

As `evolve_bdnk3d!`, recording the listed (ℓ,m) moments of δε: `Q[n,j]` is
`moments[j]` at sample `n`.
"""
function evolve_bdnk3d_moments!(st::BDNKState, e::BDNK3D; dt::Float64, nsteps::Int, sample::Int=1,
                                moments::Vector{Tuple{Int,Int}}=[(2,0)])
    N = e.s.grid.N
    k1=BDNKState(N); k2=BDNKState(N); k3=BDNKState(N); k4=BDNKState(N); tmp=BDNKState(N)
    scr = Scratch(N)
    ts = Float64[]; rows = Vector{Vector{Float64}}()
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt)
            push!(rows, [ylm_moment(st.δε, e.s; l=l, m=m) for (l,m) in moments])
        end
        n == nsteps && break
        _rk4!(st, e, dt, k1,k2,k3,k4, tmp, scr)
    end
    Q = Matrix{Float64}(undef, length(ts), length(moments))
    for (n, r) in enumerate(rows); Q[n, :] .= r; end
    return ts, Q
end

"""
    bdnk3d_qnm(eos, εc; N=32, l=2, m=0, η̂=0.0, ζ̂=0.0, τ̂=2.0, frame=nothing, η̂_frame=nothing,
               σ_ko=0.01, cut=true, κmin=0.5, A=1e-3, dt_fac=0.20, Lfac=1.2,
               f_ref_kHz=1.883, nperiods=12, samples_per_period=50,
               Lunit_km=Msun_to_km, control=(η̂>0), kw...) -> NamedTuple

The 3+1D BDNK quasi-normal-mode driver: build the star and engine, seed (ℓ,m), evolve for
`nperiods` of the reference frequency, and read the mode out with the gated two-estimator
analysis (`analyze_qnm`). With `control=true` the run is repeated at η̂=0 on the SAME grid
and frame and the viscous damping is the DIFFERENCE, γ_visc = γ(η̂) − γ(0): the numerical
floor γ(0) (Kreiss–Oliger + surface, ≈1.3e−3 M⊙⁻¹ at N=40) is 52% of γ_visc there, so
quoting γ(η̂) raw overstates the damping by 1.5×. Returns
  `main`, `ctrl`     — `analyze_qnm` results (ctrl = nothing without control)
  `f_kHz`, `γ`, `γ0`, `γ_visc`, `Q_visc` (= ω/2γ_visc), `τ_visc_ms`
  `df_kHz`           — Rayleigh resolution 1/T: nothing finer is RESOLVED
  `stable`           — both records finite and decaying
  `bound_violation`  — fraction of cells/volume violating τ_Q w₀ c_s² > (4/3)η at η̂
  `causal_den_min`   — min(τ_Q w₀ − η) over the star (must be > 0)
  `M`, `R_km`, `N`, `T`, `dt`
`frame=(a₁,a₂)` selects a Clarisse frame (τ built from η̂_frame, default η̂); otherwise the
scalar `τ̂` is used. `f_ref_kHz` only sets the record length and the search band — use
the 1D shooting value when available. Cost ∝ N⁴; N=24 with 12 periods is ~1 min,
N=48 with 23 periods ~15 min single-threaded.
"""
function bdnk3d_qnm(eos::BarotropicEOS, εc::Real; N::Int=32, l::Int=2, m::Int=0,
                    η̂::Real=0.0, ζ̂::Real=0.0, τ̂::Real=2.0,
                    frame::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
                    η̂_frame::Union{Nothing,Real}=nothing,
                    σ_ko::Real=0.01, cut::Bool=true, κmin::Real=0.5, A::Real=1e-3,
                    dt_fac::Real=0.20, Lfac::Real=1.2,
                    f_ref_kHz::Real=1.883, nperiods::Real=12, samples_per_period::Int=50,
                    Lunit_km::Real=Msun_to_km, control::Bool=(η̂ > 0), kw...)
    s = build_star3d(eos, float(εc); N=N, Lfac=float(Lfac))
    ν_ref = float(f_ref_kHz) * kHz_to_km * Lunit_km          # cyclic, code units
    P = 1/ν_ref; T = nperiods*P
    dt = dt_fac*s.grid.dx; nsteps = ceil(Int, T/dt)
    sample = max(1, floor(Int, P/(samples_per_period*dt)))
    ηf = η̂_frame === nothing ? (η̂ > 0 ? η̂ : nothing) : η̂_frame
    (frame !== nothing && ηf === nothing) &&
        throw(ArgumentError("bdnk3d_qnm: frame=(a₁,a₂) with η̂=0 needs η̂_frame>0"))
    # NB: distinct names inside the closure. A Julia closure that assigns a name which is
    # ALSO assigned in the enclosing function writes to the enclosing variable — the
    # control run would silently replace the viscous engine (caught: it reported a 0%
    # bound violation for F₁, which violates the bound in 100% of cells).
    function run(ηh)
        eng = setup_bdnk3d(s; η̂=float(ηh), ζ̂=float(ζ̂), τ̂=float(τ̂), frame=frame, η̂_frame=ηf,
                           σ_ko=float(σ_ko), cut=cut, κmin=float(κmin), kw...)
        st = BDNKState(N); seed_bdnk_ylm!(st, eng; l=l, m=m, A=float(A))
        tt, QQ = evolve_bdnk3d_moments!(st, eng; dt=dt, nsteps=nsteps, sample=sample, moments=[(l,m)])
        (analyze_qnm(tt, QQ[:,1]; ν_ref=ν_ref, Lunit_km=Lunit_km), eng, tt, QQ[:,1])
    end
    main, e, ts, q = run(η̂)
    ctrl = control ? run(0.0)[1] : nothing
    γ0 = ctrl === nothing ? NaN : ctrl.γ
    γv = ctrl === nothing ? NaN : main.γ - γ0
    ω  = 2π * main.f_kHz * kHz_to_km * Lunit_km                # geometric angular frequency
    # damping time in ms: τ[code] = 1/γ ; code length unit L_unit km ⇒ seconds = L/c
    τ_ms = isfinite(γv) && γv > 0 ? (1/γv) * Lunit_km / 299792.458 * 1e3 : NaN
    (main = main, ctrl = ctrl,
     f_kHz = main.f_kHz, f_pencil_kHz = main.f_pencil_kHz, γ = main.γ, γ0 = γ0, γ_visc = γv,
     Q_visc = isfinite(γv) && γv > 0 ? ω/(2γv) : NaN, τ_visc_ms = τ_ms,
     df_kHz = main.df_kHz, nperiods = main.nperiods,
     stable = main.stable && (ctrl === nothing || ctrl.stable),
     bound_violation = bdnk_bound_violation(e), causal_den_min = bdnk_causal_denominator(e),
     M = s.M, R_km = s.R * Lunit_km, N = N, T = T, dt = dt, ts = ts, q = q)
end

end # module CowlingBDNK3D
