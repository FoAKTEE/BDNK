using Test
using BDNKStar

const G2KM = BDNKStar.Units.gram_per_cm3_to_km_minus2

@testset "TOV: reproduce Bussières EOS1 (M=1.27 M☉, R=8.86 km)" begin
    # Bussières et al. 2604.13208 EOS1: energy polytrope p = κ ε^{1+1/n},
    # κ=100 km², n=1, central density ρ_c = 3e15 g/cm³.
    eos = PolytropeEnergy(100.0, 1.0)
    εc = 3e15 * G2KM
    st = solve_tov(eos, εc; h=5e-4)
    M = mass_solar(st); R = st.R
    @info "Bussières EOS1 reproduction" M R compactness=st.M/st.R
    @test isapprox(M, 1.27; atol=0.01)     # paper value 1.27 M☉ (2 d.p.)
    @test isapprox(R, 8.86; atol=0.02)     # paper value 8.86 km
end

@testset "TOV: reproduce Shum M_T=1.4 M☉ (M☉=G=c=1 units)" begin
    # Shum et al. 2509.15303: cold Γ=2 EOS p=κρ², κ=100, ρ0c=0.00128 M☉⁻²
    # In geometric M☉=G=c=1 units the raw TOV mass IS the mass in solar masses.
    eos = ShumPolytrope(100.0)
    ρ0c = 0.00128; εc = ρ0c + 100*ρ0c^2
    st = solve_tov(eos, εc; h=2e-4)
    @info "Shum reproduction" M_Msun=st.M R_Msun=st.R
    @test isapprox(st.M, 1.4; atol=0.02)        # paper value M_T=1.4 M☉
end

@testset "TOV: RK4 convergence + physical sanity" begin
    eos = PolytropeEnergy(100.0, 1.0)
    εc = 3e15 * G2KM
    s1 = solve_tov(eos, εc; h=1e-3)
    s2 = solve_tov(eos, εc; h=5e-4)
    # M, R converged between resolutions
    @test isapprox(s1.M, s2.M; rtol=1e-4)
    @test isapprox(s1.R, s2.R; rtol=1e-4)
    # monotone pressure profile, Buchdahl bound, surface pressure ~ 0
    @test issorted(s2.p; rev=true)
    @test 2*s2.M/s2.R < 8/9
    @test s2.p[end] == 0.0
    # Schwarzschild exterior match enforced: e^{ν(R)} = 1 - 2M/R
    @test isapprox(exp(s2.ν[end]), 1 - 2*s2.M/s2.R; rtol=1e-12)
end

@testset "TOV: mass-radius sequence is physical (M_max exists)" begin
    eos = PolytropeEnergy(100.0, 1.0)
    εcs = [k * G2KM for k in (1e15, 2e15, 3e15, 5e15, 8e15, 1.2e16, 2e16)]
    Ms = [mass_solar(solve_tov(eos, εc; h=1e-3)) for εc in εcs]
    # mass rises then turns over (a maximum mass in the sequence)
    @test maximum(Ms) > 1.3
    @test argmax(Ms) < length(Ms)          # turnover before the densest model
end
