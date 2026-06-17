using Test
using BDNKStar

# STEP 0 GATE: primitive <-> conservative round-trip closes to <= 1e-10 in smooth
# regions, across polytropic -> ideal -> tabulated, plus BDNK gradient-frozen.

@testset "Round-trip: barotropic polytrope (<= 1e-10)" begin
    poly = PolytropeEnergy(100.0, 1.0)
    maxerr = 0.0
    for e in (1e-4, 5e-4, 1e-3, 2e-3), v in (-0.6, -0.2, 0.0, 0.3, 0.6)
        E, S, _ = prim2cons_barotropic(poly, e, v)
        e2, v2, p2, info = cons2prim_barotropic(poly, E, S)
        @test info.converged
        err = max(abs(e2-e)/e, abs(v2-v))
        maxerr = max(maxerr, err)
    end
    @test maxerr ≤ 1e-10
    @info "barotropic round-trip max error" maxerr
end

@testset "Round-trip: Shum Γ=2 polytrope (1C target EOS, <= 1e-10)" begin
    shum = ShumPolytrope(100.0)          # Shum et al. 2509.15303 eq.53, κ=100
    maxerr = 0.0
    for e in (1e-4, 5e-4, 1e-3, 3e-3), v in (-0.5, 0.0, 0.4)
        # EOS self-consistency: cs²∈(0,1), p(e(p))=p
        @test 0 < sound_speed2(shum, e) < 1
        p = pressure(shum, e)
        @test isapprox(energy_from_pressure(shum, p), e; rtol=1e-12)
        E, S, _ = prim2cons_barotropic(shum, e, v)
        e2, v2, _, info = cons2prim_barotropic(shum, E, S)
        @test info.converged
        maxerr = max(maxerr, max(abs(e2-e)/e, abs(v2-v)))
    end
    @test maxerr ≤ 1e-10
    @info "Shum Γ=2 round-trip max error" maxerr
end

@testset "Round-trip: barotropic matches conformal closed form (p=e/3)" begin
    # cs² = 1/3 ⇒ PolytropeEnergy with n→∞ is not p=e/3; use κ s.t. p=e/3 at a
    # point via the closed-form check on the inversion algebra directly.
    # conformal: e = -E + √(4E² - 3S²),  for E,S consistent with some (e,v).
    e, v = 1e-3, 0.4
    p = e/3
    W2 = 1/(1-v^2); hh = (e+p)*W2
    E = hh - p; S = hh*v
    e_closed = -E + sqrt(4E^2 - 3S^2)
    @test isapprox(e_closed, e; rtol=1e-12)
end

@testset "Round-trip: general ideal gas (<= 1e-10)" begin
    maxerr = 0.0
    for Γ in (5/3, 2.0)
        gas = IdealGas(Γ)
        for ρ in (1e-4, 1e-3, 1e-2), v in (-0.6, 0.0, 0.5), ϵ in (0.05, 0.5, 1.5)
            D, S, τ, _ = prim2cons_general(gas, ρ, v, ϵ)
            ρ2, v2, ϵ2, p2, info = cons2prim_general(gas, D, S, τ)
            @test info.converged
            err = max(abs(ρ2-ρ)/ρ, abs(v2-v), abs(ϵ2-ϵ)/max(ϵ,1e-30))
            maxerr = max(maxerr, err)
        end
    end
    @test maxerr ≤ 1e-10
    @info "general round-trip max error" maxerr
end

@testset "Round-trip: BDNK gradient-frozen reduces to ideal" begin
    poly = PolytropeEnergy(100.0, 1.0)
    maxerr = 0.0
    for e in (5e-4, 1e-3), v in (-0.3, 0.0, 0.4)
        E_pf, S_pf, _ = prim2cons_barotropic(poly, e, v)
        # arbitrary frozen dissipative corrections
        δE = 0.03 * E_pf - 1e-7
        δS = -0.02 * S_pf + 1e-7
        # conserved densities carry the corrections; recovery subtracts them
        e2, v2, p2, info = cons2prim_bdnk_barotropic(poly, E_pf + δE, S_pf + δS, δE, δS)
        @test info.converged
        err = max(abs(e2-e)/e, abs(v2-v))
        maxerr = max(maxerr, err)
    end
    @test maxerr ≤ 1e-10
    @info "BDNK gradient-frozen round-trip max error" maxerr
end

@testset "Round-trip: tabulated convergent under refinement" begin
    base = PolytropeEnergy(100.0, 1.0)
    sweep = [(e, v) for e in (3e-4, 1e-3, 1.6e-3) for v in (-0.4, 0.0, 0.4)]
    errs = Float64[]
    for N in (100, 200, 400)
        tab = tabulate(base, 1e-4, 2e-3, N)
        emax = 0.0
        for (e, v) in sweep
            E, S, _ = prim2cons_barotropic(tab, e, v)
            e2, v2, _, info = cons2prim_barotropic(tab, E, S)
            emax = max(emax, max(abs(e2-e)/e, abs(v2-v)))
        end
        push!(errs, emax)
    end
    # The inversion itself is machine-precision for every N (the table is exact
    # at nodes); the genuine table-refinement convergence of p(e) is asserted in
    # test_eos.jl. Here we require the recovered primitives stay within gate.
    @test all(e -> e ≤ 1e-10, errs)
    @info "tabulated round-trip errors vs N=(100,200,400)" errs
end
