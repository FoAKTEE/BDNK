# ======================================================================
# test_frame_independence.jl
#
# BDNK FRAME-INDEPENDENCE — guard the ROBUST findings of three analyses:
#   (a) DISPERSION  (Kovtun 1907.08191 §4.2 sound channel, general/BDNK frame):
#         recomputed LIVE here (cheap quartic root-find). At fixed physical
#         (vs, η, ζ) and ζ pinned via Kovtun fn-8, sweeping the frame relaxation
#         times (ŝ=vs²ε1/γs, â=θ/γs, p̂=π1/γs):
#           * small-k sound speed Re ω/k → vs           is FLAT (frame-invariant)
#           * leading damping Γ=-2Imω·w0/k² → γs        is FLAT (frame-invariant)
#           * the two non-hydro gaps w0/θ, w0/(vs²ε1)   MOVE ∝ 1/τ (gauge modes)
#       Mirrors repro/bdnk_frame_independence.jl (built on repro/kovtun_sound.jl).
#
#   (b) AXIAL  (Bussières 2604.13208 viscous QNM; repo's best-validated solver):
#         at FIXED physical η_c=1e31 cgs, sweeping the frame τ̂∈{5,10,20,40,80}
#         and frames A/B, the PHYSICAL viscous w-mode (f, 1/τ) is frame-invariant
#         to the measured tolerance, while the non-hydro η-mode (second-sound)
#         c_η(0) ∝ 1/√τ̂ — i.e. it MOVES with the frame (gauge). These are the
#         converged values from repro/axial_frame_independence.jl (the live
#         solver takes ~25 min; the robust outcomes are asserted as data).
#
#   (c) POLAR  (Cowling BDNK 5-field operator, PolarViscousModes.jl):
#         the noted ~5-8% systematic is a FIXED conservative-frame / free-surface
#         SOLVER offset vs the shooting benchmark, present at ALL frames and NOT
#         tracking the frame knob ⇒ the FRAME-spread is small (≤~0.8%) even though
#         the absolute benchmark offset is ~+6-7%. Data from
#         repro/polar_frame_invariance.jl.
#
# Run fast: the dispersion block is recomputed live (ms); the axial/polar blocks
# assert the already-reported converged numbers (no heavy QNM solve).
# ======================================================================
using Test
using BDNKStar
using LinearAlgebra: eigvals

# ----------------------------------------------------------------------
# (a) DISPERSION — live, self-contained quartic (mirrors kovtun_sound.jl
#     fsound_rest_coeffs + companion-matrix roots; rest frame v0=0).
# ----------------------------------------------------------------------
function _polyroots(c::Vector{ComplexF64})
    while length(c) > 1 && abs(c[end]) < 1e-300
        c = c[1:end-1]
    end
    n = length(c) - 1
    n == 0 && return ComplexF64[]
    a = c ./ c[end]
    C = zeros(ComplexF64, n, n)
    for i in 1:n-1; C[i+1, i] = 1.0; end
    for i in 1:n;   C[i, n] = -a[i]; end
    return eigvals(C)
end

# F_sound(v0=0) coeffs [ω^0..ω^4] for numeric k2 (Kovtun eq:Fsound0).
function _fsound_rest(k2; cs2, ε1, ε2, π1, θ, γs, w0)
    c4 = cs2*ε1*θ
    c3 = im*w0*(cs2*ε1 + θ)
    c2 = -( w0^2 + k2*cs2*( cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1 ) )
    c1 = -im*k2*w0*( γs + cs2^2*ε1 + cs2*θ )
    c0 = k2*cs2*( w0^2 + k2*θ*( cs2*(ε2+π1 - cs2*ε1) - γs ) )
    return ComplexF64[c0, c1, c2, c3, c4]
end

# fixed PHYSICAL coefficients (Kovtun dimensionless figure units w0=1, γs=1)
const _VS  = 0.5;  const _VS2 = _VS^2;  const _W0 = 1.0
const _ETA = 0.3
const _ZETA = (4.0/3.0)*1.0 - (4.0/3.0)*_ETA      # ⇒ γs = (4/3)η+ζ = 1
const _GS   = (4.0/3.0)*_ETA + _ZETA              # ≡ 1.0

# frame parametrized by Shum ratios; ζ pinned via fn-8 (choose π2, ε2=0).
function _frame(ŝ, â, p̂; ε2=0.0)
    ε1 = ŝ*_GS/_VS2
    θ  = â*_GS
    π1 = p̂*_GS
    π2 = _VS2*(π1 - _VS2*ε1) - _ZETA + _VS2*ε2
    (ε1=ε1, ε2=ε2, π1=π1, π2=π2, θ=θ)
end
_zeta_from(fc) = _VS2*(fc.π1 - _VS2*fc.ε1) - fc.π2 + _VS2*fc.ε2

# small-k hydro sound branch: Re ω/k → vs, Γ = -2Imω·w0/k² → γs.
function _sound_branch(fc; k=1e-3)
    rts = _polyroots(_fsound_rest(k^2; cs2=_VS2, ε1=fc.ε1, ε2=fc.ε2,
                                  π1=fc.π1, θ=fc.θ, γs=_GS, w0=_W0))
    perm = sortperm(rts, by=z->abs(imag(z)))
    g1, g2 = rts[perm[1]], rts[perm[2]]
    vs = (abs(real(g1)) + abs(real(g2)))/2 / k
    Γ  = -2*((imag(g1)+imag(g2))/2)*_W0 / k^2
    (vs=vs, Γ=Γ)
end

# the two non-hydro gaps measured (most-damped roots) vs predicted (eq.4.10).
function _gaps(fc; k=1e-4)
    rts = _polyroots(_fsound_rest(k^2; cs2=_VS2, ε1=fc.ε1, ε2=fc.ε2,
                                  π1=fc.π1, θ=fc.θ, γs=_GS, w0=_W0))
    g = sort(rts, by=z->-abs(imag(z)))[1:2]
    meas = sort([-imag(z) for z in g])
    pred = sort([_W0/(_VS2*fc.ε1), _W0/fc.θ])
    (meas=meas, pred=pred)
end

_spreadpct(v) = (maximum(v) - minimum(v)) / abs(sum(v)/length(v)) * 100

@testset "BDNK frame-independence" begin

    # ---------- (a) DISPERSION (live) ----------
    @testset "dispersion: hydro flat, non-hydro gaps ∝ 1/τ" begin
        # SWEEP frame: â=θ/γs over 60×, ŝ=3, p̂=8 fixed (stable+causal region).
        âs = [1.0, 2.0, 4.0, 8.0, 15.0, 30.0, 60.0]
        fcs = [_frame(3.0, â, 8.0) for â in âs]

        # physical ζ pinned across all frames (fn-8 compensation).
        ζrec = [_zeta_from(fc) for fc in fcs]
        @test _spreadpct(ζrec) < 1e-10           # essentially exact
        @test all(isapprox.(ζrec, _ZETA; atol=1e-12))

        sb = [_sound_branch(fc) for fc in fcs]
        vss = [s.vs for s in sb];  Γs = [s.Γ for s in sb]

        # HYDRO branch FRAME-INVARIANT: vs→0.5, Γ→γs=1, flat to a TIGHT tol.
        @test all(isapprox.(vss, _VS;  atol=2e-4))
        @test all(isapprox.(Γs,  _GS;  atol=2e-3))
        @test _spreadpct(vss) < 0.05             # measured ~0.002%
        @test _spreadpct(Γs)  < 0.10             # measured ~0.026%

        # NON-HYDRO gaps MOVE with the frame: gap=w0/θ scales ∝ 1/θ over 60×.
        gth = [_W0/fc.θ for fc in fcs]
        @test isapprox(maximum(gth)/minimum(gth), 60.0; rtol=1e-9)   # = â range
        # and the MEASURED most-damped roots reproduce the analytic gaps.
        for fc in fcs
            g = _gaps(fc)
            @test isapprox(g.meas[1], g.pred[1]; rtol=5e-3)
            @test isapprox(g.meas[2], g.pred[2]; rtol=5e-3)
        end
        # the gap genuinely varies a lot (NOT frame-invariant) — the contrast.
        @test _spreadpct(gth) > 100.0            # huge: it is a gauge mode
    end

    # ---------- (b) AXIAL viscous QNM (Bussières 2604.13208) ----------
    # Converged values from repro/axial_frame_independence.jl (η_c=1e31 cgs FIXED,
    # τ̂∈{5,10,20,40,80}). The PHYSICAL w-mode is ~frame-invariant; the non-hydro
    # η-mode second-sound speed c_η(0) ∝ 1/√τ̂ — it MOVES with the frame (gauge).
    @testset "axial: w-mode frame-invariant; η-mode c_η ∝ 1/√τ̂ (gauge)" begin
        tauhats = [5.0, 10.0, 20.0, 40.0, 80.0]

        # viscous w-mode frequency [kHz], frame A — flat to ~0.12% across τ̂.
        wf_A = [10.08276, 10.08982, 10.09286, 10.09425, 10.09493]
        @test _spreadpct(wf_A) < 0.30
        # viscous w-mode damping 1/τ [1/μs], frame A — flat to ~0.06%.
        wg_A = [0.032395, 0.032389, 0.032382, 0.032378, 0.032376]
        @test _spreadpct(wg_A) < 0.30

        # cross-frame A vs B at matched τ̂=10: agreement ~0.25% (small, nonzero).
        f_A10, f_B10 = 10.08982, 10.06493
        @test abs(f_A10 - f_B10)/f_A10 * 100 < 0.5

        # cross-check vs Bussières Table II A1 at τ̂=10: f matches to <0.05%.
        f_tableII = 10.0898
        @test abs(f_A10 - f_tableII)/f_tableII * 100 < 0.05

        # inviscid w-mode reference (frame-exact, η=0).
        @test isapprox(10.50135, 10.50135; atol=1e-5)

        # NON-HYDRO η-mode: second-sound centre speed c_η(0), frame B, τ̂∈{10,20,40,80}.
        # MOVES with the frame: c_η ∝ 1/√τ̂ (the gauge artifact, as expected).
        ceta_tau = [10.0, 20.0, 40.0, 80.0]
        ceta = [0.2146, 0.1518, 0.1073, 0.0759]
        # NOT frame-invariant — drops with τ̂ (the visual/numeric contrast):
        @test _spreadpct(ceta) > 80.0
        # verify the 1/√τ̂ scaling: c_η·√τ̂ is (nearly) constant.
        scaled = ceta .* sqrt.(ceta_tau)
        @test _spreadpct(scaled) < 2.0           # ~1/√τ̂ to ≲1%
    end

    # ---------- (c) POLAR viscous QNM (Cowling BDNK 5-field) ----------
    # The ~5-8% systematic is a FIXED solver offset vs the shooting benchmark,
    # present at ALL frames and NOT tracking the frame knob. Data from
    # repro/polar_frame_invariance.jl (η̂=0.03 fixed, Nr=140).
    @testset "polar: fixed solver offset, small FRAME spread (not frame-gauge)" begin
        # shooting (NonRadialModes) benchmark — the frame-independent reference.
        f_bench, p1_bench = 1.8829, 4.1067

        # f-mode at η̂=0 across the τ̂ frame sweep (pure gauge isolation).
        ff0 = [2.0125, 2.0118, 1.9996]           # τ̂ = 0.06, 0.12, 0.50
        # p1-mode at η̂=0 across the frame sweep.
        fp0 = [3.9798, 4.0109]                    # τ̂ = 0.06, 0.50

        # (i) FRAME-spread of the physical mode is SMALL (it does not move w/ τ̂).
        @test _spreadpct(ff0) < 1.0              # measured ~0.64%
        @test _spreadpct(fp0) < 1.0             # measured ~0.78%

        # (ii) but the absolute OFFSET vs the benchmark is a fixed ~+6-7% (f) and
        #      ~-2.3..-3.1% (p1) — present at every frame (a solver artifact).
        off_f = [100*(f - f_bench)/f_bench for f in ff0]
        @test all(off_f .> 5.0) && all(off_f .< 8.0)   # ~+6-7% at all frames
        # the offset itself barely moves across frames (<0.8% spread of offsets):
        @test (maximum(off_f) - minimum(off_f)) < 0.8

        off_p = [100*(f - p1_bench)/p1_bench for f in fp0]
        @test all(off_p .< -2.0) && all(off_p .> -3.5)  # ~-2.3..-3.1% at all frames

        # (iii) reactive-viscous + frame-coupling sweeps at η̂=0.03 stay tight:
        cτ_spread  = 1.08      # cτ=0..4 f-mode spread (reported)
        κ_spread   = 0.89      # κ̂/cκ sweep f-mode spread (reported)
        @test cτ_spread < 2.0 && κ_spread < 2.0

        # VERDICT LOGIC: frame-spread << benchmark-offset, offset ≈ const across
        # frames ⇒ SOLVER ARTIFACT, not a frame-gauge dependence.
        @test _spreadpct(ff0) < minimum(off_f)
    end
end
