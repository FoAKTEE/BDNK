using Test
using BDNKStar

@testset "Bjorken flow (PMP 2209.09265): RK4 Q→16 + analytic + diagnostic" begin
    Γ, m, n0, e0 = 4/3, 1.0, 1.0, 1.0
    ε0 = bjorken_inviscid_analytic(1.0; Γ=Γ, m=m, n0=n0, e0=e0)
    errs = Float64[]
    for N in (250, 500, 1000, 2000)
        τs, εs = bjorken_evolve_rk4(1.0, 20.0, ε0; N=N, Γ=Γ, m=m, n0=n0)
        εa = [bjorken_inviscid_analytic(τ; Γ=Γ, m=m, n0=n0, e0=e0) for τ in τs]
        push!(errs, maximum(abs.(εs .- εa)))
    end
    Qs = [errs[i]/errs[i+1] for i in 1:length(errs)-1]
    @info "Bjorken RK4 convergence" errs Q=Qs
    @test errs[end] < 1e-9                          # RK4 matches analytic
    @test all(q -> 15.0 < q < 17.0, Qs)            # 4th-order convergence Q→16
    # frame-independent diagnostic identity (analytic)
    τs, εs = bjorken_evolve_rk4(1.0, 20.0, ε0; N=4000, Γ=Γ, m=m, n0=n0)
    dev = maximum(abs(bjorken_diagnostic(τs[i], εs[i]; Γ=Γ, m=m, n0=n0) - (Γ-1)*m*n0/τs[i]^2)
                  for i in 2:length(τs))
    @test dev < 1e-12
end
