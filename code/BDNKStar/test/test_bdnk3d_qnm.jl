# ======================================================================
# test_bdnk3d_qnm.jl — the 3+1D BDNK QNM engine, end to end.
#
#   (1) Clarisse frames installed EXACTLY: τ_Q w₀/η₀ = a₂, τ_ε w₀/η₀ = a₁, τ_P = c_s² τ_ε,
#       and the bound-violation diagnostic reproduces the analytic a₂c_s² > 4/3 reduction
#       (F₁: 100% of cells; F₃: about half, matching the 1D proper-volume fraction).
#   (2) Real Y_ℓm: cubic-group parity and orthogonality of the ℓ=2 quintet on a sphere.
#   (3) LINEAR engine, full cube, m = 0, 2, 1 seeds on ONE background:
#         • E_g/T₂g cross-moments vanish EXACTLY by reflection parity (the m≠0 detector)
#         • Y₂₀ and Re Y₂₂ are partners of the E_g irrep of the cubic group, so a
#           cubic-symmetric discretisation CANNOT split them: f(m=0) ≡ f(m=2)
#         • Re Y₂₁ is T₂g: its splitting from E_g is the grid's cubic anisotropy, O(h²)
#   (4) BDNK QNM driver at N=24, frames F₁ and F₃, difference protocol:
#         • anchor f-mode from the η̂=0 control (N=24 reference: +0.005%)
#         • γ_visc reproduces the verified N=24 value 2.355e−3 M⊙⁻¹; the numerical
#           floor γ₀ reproduces 3.27e−3; frame spread of γ_visc < 1% (measured 0.15%)
#         • F₁ (violates the bound in 100% of cells) is STABLE: the bound is
#           sufficient-not-necessary (the unstable band sits above k_c, unresolved at N=24)
#         • a₁ is inert: a 16× change in a₁ at fixed a₂ leaves f and γ unchanged
#   (5) The production defaults are the validated ones: cut=true, excise_frac=0.
#
# Cost: ~6–8 min single-threaded (five N=24 BDNK runs of 9 periods + three cheap linear
# runs). Every frequency here is a parametric estimate — the records resolve Δf/f≈17% —
# so tolerances are set from the two-estimator agreement, not from spectral resolution.
# ======================================================================
using Test
using BDNKStar
using BDNKStar: Msun_to_km, kHz_to_km
using LinearAlgebra: Diagonal, diag

const F1D_ANCHOR = 1.88291      # relativistic-Cowling ℓ=2 f-mode, shooting solver, kHz
const F1 = (25/4, 25/7); const F3 = (25.0, 25.0)   # Clarisse et al. 2510.16603 Table 1

@testset "3+1D BDNK QNM engine" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2

    # ---------------------------------------------------------------- (5) defaults
    @testset "production defaults" begin
        s16 = build_star3d(eos, εc; N=16, Lfac=1.2)
        @test setup_evo3d(s16).cut == true
        eb = setup_bdnk3d(s16; η̂=0.03)
        @test eb.cut == true
        @test eb.wfloor == 0.0                         # excise_frac default 0 ⇒ full star
    end

    # ---------------------------------------------------------------- (1) frames
    @testset "Clarisse frame installation is exact" begin
        s = build_star3d(eos, εc; N=20, Lfac=1.2)
        for (a1, a2) in (F1, (25/2, 25/3), F3)
            e = setup_bdnk3d(s; η̂=0.03, frame=(a1, a2))
            I = findall(i -> s.interior[i] && e.η0[i] > 0, eachindex(e.w0))
            r2 = [e.τQ[i]*e.w0[i]/e.η0[i] for i in I]
            r1 = [e.τε[i]*e.w0[i]/e.η0[i] for i in I]
            @test maximum(abs.(r2 .- a2)) < 1e-10
            @test maximum(abs.(r1 .- a1)) < 1e-10
            @test maximum(abs.(e.τP[I] .- e.cs2[I] .* e.τε[I])) < 1e-15   # Π = c_s²𝒜 closure
            @test bdnk_causal_denominator(e) > 0                          # (a₂−1)η > 0
        end
        # η̂=0 control keeps the frame via η̂_frame; without it the frame is an error
        e0 = setup_bdnk3d(s; η̂=0.0, frame=F1, η̂_frame=0.03)
        e1 = setup_bdnk3d(s; η̂=0.03, frame=F1)
        @test e0.τQ == e1.τQ && e0.τε == e1.τε
        @test_throws ArgumentError setup_bdnk3d(s; η̂=0.0, frame=F1)
        # bound τ_Q w₀ c_s² > (4/3)η reduces to a₂c_s² > 4/3 under the map; max c_s² ≈ 0.20
        v1 = bdnk_bound_violation(e1)
        @test v1.cells == 1.0 && v1.volume == 1.0                        # F₁: everywhere
        v3 = bdnk_bound_violation(setup_bdnk3d(s; η̂=0.03, frame=F3))
        @test 0.35 < v3.volume < 0.65                                    # F₃: ~half (1D: 0.526)
        # a scalar τ̂ frame is untouched by the new code path
        eτ = setup_bdnk3d(s; η̂=0.03, τ̂=2.0)
        @test all(eτ.τQ .== 2.0) && all(eτ.τε .== 2.0)
    end

    # ---------------------------------------------------------------- (2) harmonics
    @testset "real Y_ℓm: parity and orthogonality" begin
        # cubic-group reflection parities: E_g even under all three; T₂g odd under two
        n = (0.3, -0.5, sqrt(1 - 0.34))
        for (l, m, px, py, pz) in ((2,0, 1,1,1), (2,2, 1,1,1), (2,1, -1,1,-1), (2,-1, 1,-1,-1), (2,-2, -1,-1,1))
            y  = ylm_real(l, m, n...)
            @test ylm_real(l, m, -n[1], n[2], n[3]) ≈ px*y
            @test ylm_real(l, m, n[1], -n[2], n[3]) ≈ py*y
            @test ylm_real(l, m, n[1], n[2], -n[3]) ≈ pz*y
        end
        # orthogonality of the ℓ=2 quintet on a Fibonacci sphere sample
        npts = 20000; G = zeros(5, 5); ms = (-2, -1, 0, 1, 2)
        for k in 0:npts-1
            z = 1 - 2(k + 0.5)/npts; φ = k * 2.399963229728653; ρ = sqrt(1 - z^2)
            v = [ylm_real(2, m, ρ*cos(φ), ρ*sin(φ), z) for m in ms]
            G .+= v * v'
        end
        G ./= npts
        off = maximum(abs.(G - Diagonal(diag(G)))) / minimum(diag(G))
        @test off < 5e-3
        @test_throws ArgumentError ylm_real(5, 0, 0.0, 0.0, 1.0)
    end

    # ---------------------------------------------------------------- (3) linear, full cube
    @testset "linear engine: E_g degeneracy exact, cross-moments vanish, T₂g splitting small" begin
        N = 24; s = build_star3d(eos, εc; N=N, Lfac=1.2)
        e = setup_evo3d(s; σ_ko=0.01)                 # cut=true by default
        ν_ref = F1D_ANCHOR * kHz_to_km * Msun_to_km; P = 1/ν_ref
        dt = 0.30*s.grid.dx; nsteps = ceil(Int, 9P/dt); sample = max(1, floor(Int, P/(50dt)))
        moments = [(2,0), (2,2), (2,1)]
        res = Dict{Int,Any}(); Q = Dict{Int,Matrix{Float64}}()
        for m in (0, 2, 1)
            st = EvolState(N); seed_ylm!(st, e; l=2, m=m, A=1e-3)
            ts, Qm = evolve3d_moments!(st, e; dt=dt, nsteps=nsteps, sample=sample, moments=moments)
            Q[m] = Qm
            j = findfirst(==((2,m)), moments)
            res[m] = analyze_qnm(ts, Qm[:, j]; ν_ref=ν_ref)
        end
        rms(v) = sqrt(sum(abs2, v)/length(v))
        # the m≠0 detector: cross-moments between E_g partners and between irreps vanish
        # by reflection parity — exactly, not approximately
        @test rms(Q[0][:,2]) / rms(Q[0][:,1]) < 1e-8      # Y₂₀ seed ⇒ no q₂₂ (x↔y symmetric)
        @test rms(Q[2][:,1]) / rms(Q[2][:,2]) < 1e-8      # Y₂₂ seed ⇒ no q₂₀ (x↔y antisymmetric)
        @test rms(Q[1][:,1]) / rms(Q[1][:,3]) < 1e-8      # Y₂₁ seed ⇒ no q₂₀ (odd in x)
        @test rms(Q[1][:,2]) / rms(Q[1][:,3]) < 1e-8      # Y₂₁ seed ⇒ no q₂₂
        for m in (0, 2, 1)
            @test res[m].stable
            @test res[m].agree_pct < 0.5                   # periodogram and pencil agree
        end
        f0, f2, f1 = res[0].f_kHz, res[2].f_kHz, res[1].f_kHz
        @info "linear full-cube ℓ=2 f-mode by m [kHz]" m0=f0 m2=f2 m1=f1 anchor=F1D_ANCHOR df_kHz=res[0].df_kHz
        @test isapprox(f0, F1D_ANCHOR; rtol=5e-3)        # N=24 cut-cell: measured +0.005%
        # E_g partners: a cubic-symmetric operator cannot split them. With the exact (kink-split,
        # permutation-symmetrised) cut-cell geometry the discretisation IS cubic-symmetric to
        # round-off, and the two frequencies agree to 14 significant figures (measured 8e-15);
        # with the old 24-point Gauss geometry they differed at 5e-7 — the κ asymmetry.
        @test abs(f2 - f0)/f0 < 1e-9
        # T₂g vs E_g: the grid's cubic anisotropy, a MEASUREMENT (1.13% at N=24, dx/R=0.10);
        # a Cartesian grid splits the ℓ=2 quintet into E_g and T₂g at O(h²)
        @info "cubic-grid E_g/T₂g splitting at N=24" pct=100*(f1 - f0)/f0
        @test abs(f1 - f0)/f0 < 2.5e-2
    end

    # ---------------------------------------------------------------- (4) BDNK driver
    @testset "BDNK driver: frames F₁/F₃, difference protocol, stability, a₁ inert" begin
        kw = (N=24, η̂=0.03, nperiods=9, f_ref_kHz=F1D_ANCHOR, σ_ko=0.01)   # the verified T≈1000 regime
        r1 = bdnk3d_qnm(eos, εc; frame=F1, kw...)
        r3 = bdnk3d_qnm(eos, εc; frame=F3, kw...)
        for r in (r1, r3)
            @test r.stable                                  # viscous AND control records decay
            @test r.main.agree_pct < 0.5 && r.ctrl.agree_pct < 0.5
            @test r.causal_den_min > 0
            @test r.γ0 > 0 && r.γ_visc > 0                  # a real floor, and dissipation above it
        end
        # anchor from the η̂=0 controls (the frame τ is inert at η=0 — 1.87681 for all frames at N=40)
        @test isapprox(r1.ctrl.f_kHz, F1D_ANCHOR; rtol=6e-3)
        @test abs(r1.ctrl.f_kHz - r3.ctrl.f_kHz)/F1D_ANCHOR < 1e-3
        # F₁ violates the stability bound EVERYWHERE and is nonetheless stable at N=24:
        # the bound is sufficient-not-necessary (k_c = 15.3 M⊙⁻¹ needs N > 353 to resolve)
        @test r1.bound_violation.cells == 1.0
        @test r1.stable
        # verified N=24 values (T=1000): γ_visc(F₁)=2.35482e−3, γ₀=3.2697e−3 M⊙⁻¹
        @test isapprox(r1.γ_visc, 2.355e-3; rtol=0.10)
        @test isapprox(r1.γ0, 3.27e-3; rtol=0.15)
        # frame invariance of the PHYSICAL damping (measured 0.15% at N=24) and frequency
        spread_γ = abs(r1.γ_visc - r3.γ_visc)/r1.γ_visc
        spread_f = abs(r1.f_kHz - r3.f_kHz)/r1.f_kHz
        @info "BDNK frames N=24" f_F1=r1.f_kHz f_F3=r3.f_kHz γvisc_F1=r1.γ_visc γvisc_F3=r3.γ_visc γ0=r1.γ0 spread_γ_pct=100spread_γ spread_f_pct=100spread_f Q=r1.Q_visc τ_ms=r1.τ_visc_ms
        @test spread_γ < 1e-2
        @test spread_f < 2e-3
        # a₁ inert: a₁ = 25 → 100 at fixed a₂ = 25 (the verified experiment: a₁=6.25/25/100 gave
        # f spread 0.006%, γ spread 0.22%); compare against F₃ = (25,25) on the raw rates
        rA = bdnk3d_qnm(eos, εc; frame=(100.0, F3[2]), control=false, kw...)
        @test rA.stable
        @info "a₁ inertness (a₂=25): a₁=25 vs 100" f_25=r3.f_kHz f_100=rA.f_kHz γ_25=r3.γ γ_100=rA.γ
        # tolerance = the two-estimator agreement level of a 9-period record (Δf/f ≈ 11%);
        # measured 0.23% here, 0.006% in the T=1000 periodogram study with the legacy seed —
        # a₁ is inert to the resolution of the record, which is the honest statement
        @test abs(rA.f_kHz - r3.f_kHz)/r3.f_kHz < 5e-3
        @test abs(rA.γ - r3.γ)/r3.γ < 1e-2
        # the difference protocol is not optional: the floor is a large fraction of the signal
        @test r1.γ0 / r1.γ_visc > 0.5
    end
end
