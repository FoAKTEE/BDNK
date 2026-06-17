using Test
using BDNKStar

@testset "EOS: thermodynamic consistency" begin
    # --- cold energy polytrope p = κ e^{1+1/n} ---
    poly = PolytropeEnergy(100.0, 1.0)
    for e in (1e-4, 5e-4, 1e-3, 2e-3)
        p  = pressure(poly, e)
        cs2 = sound_speed2(poly, e)
        # analytic vs finite-difference dp/de
        de = 1e-8 * e
        cs2_fd = (pressure(poly, e+de) - pressure(poly, e-de)) / (2de)
        @test isapprox(cs2, cs2_fd; rtol=1e-6)
        # inverse: energy_from_pressure ∘ pressure = identity
        @test isapprox(energy_from_pressure(poly, p), e; rtol=1e-12)
        @test cs2 > 0
    end

    # --- Γ-law ideal gas p = (Γ-1) ρ ϵ ---
    for Γ in (5/3, 2.0)
        gas = IdealGas(Γ)
        for ρ in (1e-4, 1e-3, 1e-2), ϵ in (0.05, 0.5, 1.5)
            p = pressure(gas, ρ, ϵ)
            # ∂p/∂ρ|_ϵ and ∂p/∂ϵ|_ρ vs finite difference
            dρ = 1e-8*ρ; dϵ = 1e-8*ϵ
            dpdρ_fd = (pressure(gas, ρ+dρ, ϵ) - pressure(gas, ρ-dρ, ϵ)) / (2dρ)
            dpdϵ_fd = (pressure(gas, ρ, ϵ+dϵ) - pressure(gas, ρ, ϵ-dϵ)) / (2dϵ)
            @test isapprox(dpdrho_eps(gas, ρ, ϵ), dpdρ_fd; rtol=1e-6)
            @test isapprox(dpdeps_rho(gas, ρ, ϵ), dpdϵ_fd; rtol=1e-6)
            # sound speed = (∂p/∂ρ + (p/ρ²)∂p/∂ϵ)/h, subluminal in range
            cs2 = sound_speed2(gas, ρ, ϵ)
            @test 0 < cs2 < 1
            @test temperature(gas, ρ, ϵ) ≈ p/ρ
        end
    end
end

@testset "EOS: tabulated interpolation converges to base" begin
    base = PolytropeEnergy(100.0, 1.0)
    e_lo, e_hi = 1e-4, 2e-3
    es = range(2e-4, 1.8e-3; length=37)
    errs = Float64[]
    for N in (100, 200, 400)
        tab = tabulate(base, e_lo, e_hi, N)
        emax = maximum(abs(pressure(tab, e) - pressure(base, e)) / pressure(base, e) for e in es)
        push!(errs, emax)
    end
    @test errs[end] < errs[1]            # refining the table reduces error
    @test errs[end] < 1e-6               # 4th-order Hermite, easily met
end
