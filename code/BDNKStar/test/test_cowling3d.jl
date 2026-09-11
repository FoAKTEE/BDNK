using Test
using BDNKStar

# STAGE 3 Phase 2: 3+1D Cowling evolution engine. Fast gate — verifies the
# frozen-background construction and the numerical STABILITY + oscillation of the
# linear evolution. The quantitative f-mode↔eigensolver validation (minutes of
# evolution + FFT) lives in viz/cowling3d_fmode.jl, not in the unit suite.
@testset "3+1D Cowling evolution (Phase 2): frozen background + stable ℓ=2 oscillation" begin
    eos = ShumPolytrope(100.0)
    εc  = 0.00128 + 100 * 0.00128^2                       # M=1.4 benchmark star
    s = build_star3d(eos, εc; N=24, Lfac=1.2)

    @test isapprox(s.M, 1.4; atol=0.05)
    # metric self-consistency: √γ = e^{λ} = √(g_rr) EXACTLY by construction
    @test all(isapprox.(s.sqrtγ, sqrt.(s.e2λ); rtol=1e-12))
    # reconstruct enclosed mass m(r) from g_rr=(1−2m/r)^{-1} at a mid-star cell
    c = s.grid.N ÷ 2
    i = argmin(abs.(s.grid.x .- 0.5s.R))
    r = s.r[i,c,c]; mr = r * (1 - 1/s.e2λ[i,c,c]) / 2
    @test 0 < mr < s.M                                    # physical enclosed mass
    @test any(s.interior) && !all(s.interior)             # star + atmosphere both present

    # stable linear evolution: seed ℓ=2, evolve — must stay bounded and oscillate
    e  = setup_evo3d(s; σ_ko=0.01)
    st = EvolState(s.grid.N); seed_l2!(st, e; A=1e-3)
    m0 = maximum(abs, st.δε)
    ts, q2, qc = evolve3d!(st, e; dt=0.25*s.grid.dx, nsteps=800, sample=4)
    @test all(isfinite, q2)
    @test maximum(abs, st.δε) < 10m0                      # bounded (no blow-up)
    @test sum(diff(sign.(q2)) .!= 0) ≥ 2                  # ℓ=2 quadrupole oscillates
end
