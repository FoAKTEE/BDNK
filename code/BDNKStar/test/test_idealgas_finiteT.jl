using Test
using BDNKStar

# Finite-temperature ideal-gas TOV star + the BDNK heat-conduction sector ON it.
#
# `solve_tov_idealgas` integrates the hydrostatic structure with the ideal-gas-on-
# an-adiabat polytrope p = K ρ^{Γ_struct} (a barotrope — a static star must be),
# then reconstructs the genuine finite-T ideal-gas profiles ρ, ϵ, T = p/ρ, the
# adiabatic sound speed cs² and the fixed-baryon (equilibrium) sound speed
# cn² = Γ-1 via the micro-EOS p = (Γ-1)ρϵ with the adiabatic index Γ.
#
# This exercises the heat-conduction sector: a real temperature gradient dT/dr<0,
# the Caballero–Yunes criterion cs²−cn² (here violated, as the ideal gas always
# does), and BDNK causal heat conduction (real, subluminal characteristic speeds
# with κQ > 0) on the star. Γ=5/3 ⇒ cn²=2/3 strictly subluminal.
@testset "Finite-T ideal-gas star: heat-conduction sector exercised" begin
    Γ = 5/3; K = 15.0; ρc = 1e-3
    st = solve_tov_idealgas(Γ=Γ, K=K, ρc=ρc, h=2e-4)

    # --- the star builds with a sensible M, R ---
    M = mass_solar(st); R = st.R
    @test 0.8 < M < 2.0                         # NS-like gravitational mass
    @test 8.0 < R < 20.0                        # NS-like areal radius
    @test 2*st.M/st.R < 8/9                     # Buchdahl bound
    @test st.M ≈ st.m[end] rtol=1e-12
    @test st.Γ == Γ && st.Γ_struct == Γ         # isentropic by default

    # interior indices (p > 0)
    ni = findlast(>(0), st.p)

    # --- temperature profile: positive, with a NEGATIVE gradient ---
    @test all(st.T[1:ni] .> 0)                  # T(r) > 0 in the interior
    @test st.T[1] > st.T[ni]                    # T_center > T_surface
    @test issorted(st.T[1:ni]; rev=true)        # dT/dr < 0 everywhere (monotone)
    # finite-difference gradient is strictly negative away from the center
    dTdr = (st.T[ni] - st.T[ni-5]) / (st.r[ni] - st.r[ni-5])
    @test dTdr < 0

    # --- sound speeds in (0,1): adiabatic cs² and equilibrium cn² ---
    @test all(0 .< st.cs2[1:ni] .< 1)
    @test all(0 .< st.cn2[1:ni] .< 1)
    @test all(st.cn2[1:ni] .≈ Γ - 1)            # cn² = Γ-1 (fixed-baryon)

    # --- Caballero–Yunes heat-conduction criterion cs²−cn² along the star ---
    # the ideal gas violates it everywhere (cs²−cn² < 0): the sector is NON-trivial
    diff = st.cs2[1:ni] .- st.cn2[1:ni]
    @test all(diff .< 0)
    # exact analytic value cs²−cn² = −(Γ−1)/(1+Γϵ) at the center
    ϵc = st.ϵ[1]
    @test isapprox(diff[1], -(Γ-1)/(1+Γ*ϵc); rtol=1e-8)

    # --- BDNK characteristic speeds with a heat-flux coefficient κQ ≠ 0 ---
    # causal frame (τP > 1); confirm real, non-negative, subluminal heat conduction
    κQ = st.ε[1]^0.25 / (3π)                     # representative κ_Q > 0
    @test κQ > 0
    tc = TransportCoefficients(η=1e-3, ζ=1e-3, κQ=κQ,
                               τε=2.0, τP=2.0, τQ=1.0, L=1.0)
    fl = causality_flag(st.p[1], st.ε[1], st.cs2[1], tc)
    @test fl.real_speeds                         # disc ≥ 0  (real)
    @test fl.nonneg                              # c²∓ ≥ 0
    @test fl.subluminal                          # c²₊ ≤ 1
    @test fl.causal                              # CAUSAL heat conduction
    @test 0 ≤ fl.c2_minus && fl.c2_plus ≤ 1      # speeds-squared in [0,1]

    # causal along the WHOLE interior, not just the center
    allcausal = true
    for i in 1:ni
        κQi = st.ε[i]^0.25 / (3π)
        tci = TransportCoefficients(η=1e-3, ζ=1e-3, κQ=κQi,
                                    τε=2.0, τP=2.0, τQ=1.0, L=1.0)
        allcausal &= is_causal(st.p[i], st.ε[i], st.cs2[i], tci)
    end
    @test allcausal
end

@testset "Finite-T ideal-gas star: stratified Γ_struct ≠ Γ closure" begin
    # composition/entropy stratification: structure index ≠ adiabatic index.
    Γ = 5/3; Γs = 2.0
    st = solve_tov_idealgas(Γ=Γ, Γ_struct=Γs, K=100.0, ρc=1.3e-3, h=2e-4)
    ni = findlast(>(0), st.p)
    @test st.Γ_struct == Γs && st.Γ == Γ
    @test mass_solar(st) > 0 && st.R > 0
    @test st.T[1] > st.T[ni] > 0                 # still dT/dr < 0
    @test all(0 .< st.cs2[1:ni] .< 1)
    @test all(0 .< st.cn2[1:ni] .< 1)
    @test all((st.cs2[1:ni] .- st.cn2[1:ni]) .< 0)   # CY violated
end
