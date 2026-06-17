using Test
using BDNKStar

@testset "Conformal evolution: uniform state is exactly preserved" begin
    fr = pmp_luminal_frame(10.0)
    s = init_gaussian(fr; N=129, A=0.0, x0=0.0, w=25.0, c=0.5, cfl=0.1)  # ε≡0.5
    ε0 = copy(energy_density(s))
    evolve!(s, 100)
    εT = energy_density(s)
    @test maximum(abs.(εT[5:end-4] .- ε0[5:end-4])) == 0.0   # zero flux divergence
end

@testset "Conformal evolution: Gaussian stable + physical diffusion" begin
    fr = pmp_luminal_frame(10.0)
    s = init_gaussian(fr; N=257, A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)
    peak0 = maximum(energy_density(s))
    evolve!(s, 200)
    εT = energy_density(s)
    @test all(isfinite, εT)                # no NaN/blow-up
    @test all(εT .> 0)                      # positivity of energy density
    @test maximum(εT) < peak0               # viscous spreading lowers the peak
    @test minimum(εT) ≥ 0.09                # relaxes toward the background 0.1
end

@testset "Conformal evolution: steady shock — RH asymptotics + quasi-stationary" begin
    fr = pmp_luminal_frame(10.0)
    s = init_smooth_shock(fr; N=257, εL=1.0, vL=0.8, cfl=0.1)
    mid = (1.0 + 4.40741)/2
    center(st) = (ε=energy_density(st); st.x[findfirst(k->ε[k] ≥ mid, eachindex(ε))])
    c0 = center(s)
    evolve!(s, 200)
    εT = energy_density(s)
    @test all(isfinite, εT)
    @test isapprox(εT[10], 1.0; atol=0.02)        # left asymptotic state
    @test isapprox(εT[end-9], 4.40741; atol=0.05) # right asymptotic state (RH)
    @test abs(center(s) - c0) ≤ 2*s.dx            # shock stays put to ≤2 cells / 200 steps
end
