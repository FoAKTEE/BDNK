using Test
using BDNKStar

@testset "Caballero-Yunes heat-conduction criterion c_s²−c_n²" begin
    # Ideal gas: analytic c_s²−c_n² = −(Γ−1)/(1+Γϵ) < 0 ⇒ heat-conduction UNSTABLE
    for Γ in (4/3, 5/3, 2.0), ρ in (1e-4, 1e-2), ϵ in (0.05, 0.5, 2.0)
        gas = IdealGas(Γ)
        cs2 = sound_speed2(gas, ρ, ϵ)
        cn2v = cn2(gas, ρ, ϵ)
        @test isapprox(cn2v, Γ-1; rtol=1e-12)                      # (∂p/∂ε)_n = Γ-1
        @test isapprox(cs2 - cn2v, -(Γ-1)/(1+Γ*ϵ); rtol=1e-10)     # CY difference, exact
        @test cs2 - cn2v < 0                                        # violates criterion
        @test heat_conduction_stable(gas, ρ, ϵ) == false           # ⇒ unstable
    end
end
