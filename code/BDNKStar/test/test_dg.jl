#=
    test_dg — nodal Discontinuous-Galerkin (RKDG) relativistic-hydro FALLBACK
    scheme for the FV-Cartesian failure (verdict 'fv-insufficient-need-dg').

    Stage 1 (the MVP — DG CAPTURES SHOCKS, the user's explicit ask): a 1D special-
      relativistic shock tube vs the EXACT relativistic Riemann solution. Asserts
      the limiter works (no spurious oscillations / overshoot), positivity holds
      (ρ,p>0), and the L1 error CONVERGES under refinement.
    Stage 2 (DG on the stellar surface): the 1D radial TOV (Cowling) star is
      statically well-balanced (machine-zero drift) once the r=0 coordinate
      singularity is excised; DG resolves the surface sub-cell.
    Stage 3 (does DG cure the staircase?): 2D meridional Cartesian DG with the
      non-conforming surface is statically well-balanced and short-time bounded
      under an ℓ=2 perturbation.
=#
using Test
using BDNKStar

@testset "DG (RKDG) relativistic-hydro fallback" begin
    Γ = 5/3

    @testset "Stage 1: exact relativistic Riemann reference is correct" begin
        # Martí–Müller test 1: ρL=10,pL=13.33,ρR=1,pR≈0 ⇒ p*≈1.45, v*≈0.72.
        sol = exact_riemann_sr(10.0,0.0,13.33, 1.0,0.0,1e-6; Γ=Γ, x0=0.5)
        @test 1.40 < sol.pstar < 1.50
        @test 0.70 < sol.vstar < 0.74
    end

    @testset "Stage 1 (MVP): DG shock tube — capture, positivity, convergence" begin
        t = 0.35
        sol = exact_riemann_sr(10.0,0.0,13.33, 1.0,0.0,1e-6; Γ=Γ, x0=0.5)
        errs = Float64[]
        for K in (100, 200)
            dg = setup_srdg(; K=K, xL=0.0, xR=1.0, p=3, Γ=Γ, cfl=0.2, M_tvb=1.0)
            shocktube_initial!(dg; x0=0.5, ρL=10.0,vL=0.0,pL=13.33, ρR=1.0,vR=0.0,pR=1e-6)
            evolve_srdg!(dg; tmax=t)
            xm, ρm, vm, pm = srdg_cell_means(dg)
            # positivity held at every cell mean (no negative ρ or p)
            @test minimum(ρm) > 0.0
            @test minimum(pm) > 0.0
            # subluminal everywhere (no nonphysical overshoot to |v|≥1)
            @test maximum(abs.(vm)) < 1.0
            # no large overshoot beyond the exact extrema (limiter monotone):
            ρe = [sample_riemann_sr(sol, x, t)[1] for x in xm]
            @test maximum(ρm) < 1.05*maximum(ρe)      # density overshoot < 5%
            L1 = sum(abs.(ρm .- ρe)) * (xm[2]-xm[1])
            push!(errs, L1)
        end
        # error CONVERGES under refinement (shock-capturing works)
        @test errs[2] < 0.7*errs[1]
        @test errs[2] < 0.15
    end

    @testset "Stage 1: smooth flow stays high-order/bounded (no false limiting)" begin
        # a small smooth density pulse at rest must not blow up or clip
        dg = setup_srdg(; K=96, xL=0.0, xR=1.0, p=3, Γ=Γ, cfl=0.2, M_tvb=200.0)
        set_initial!(dg, x->1.0+0.2*exp(-100*(x-0.5)^2), x->0.0, x->1.0)
        evolve_srdg!(dg; tmax=0.05)
        _, ρm, _, pm = srdg_cell_means(dg)
        @test all(isfinite, ρm)
        @test minimum(ρm) > 0.0 && minimum(pm) > 0.0
        @test maximum(ρm) < 1.25            # near the exact peak ≈1.2, no overshoot
    end

    @testset "Stage 2: radial DG star is statically well-balanced (surface cured)" begin
        κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)
        eng, st = setup_dgstar(eos, εc; K=48, p=3, M_tvb=50.0, cfl=0.2)
        ts, probe, ρc, drift = evolve_dgstar!(st, eng; tmax=120.0, sample_dt=10.0)
        @test drift < 1e-8                  # machine-zero static drift (well-balanced)
        @test all(isfinite, ρc)
    end

    @testset "Stage 3: 2D Cartesian DG static well-balancing (non-conforming surface)" begin
        κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)
        eng, st = setup_dgcart2d(eos, εc; Kx=10, Kz=10, p=2, M_tvb=50.0, cfl=0.2)
        ρc0 = dgcart2d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart2d!(st, eng; tmax=24.0, sample_dt=4.0)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 1e-10   # well-balanced static (machine zero)
        @test all(isfinite, ρc)
    end

    @testset "Stage 3: 2D DG ℓ=2 oscillates and is short-time bounded" begin
        κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)
        eng, st = setup_dgcart2d(eos, εc; Kx=10, Kz=10, p=2, M_tvb=50.0, cfl=0.2)
        seed_dgcart2d_l2!(st, eng; A=1e-5)
        ρc0 = dgcart2d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart2d!(st, eng; tmax=80.0, sample_dt=4.0)
        @test all(isfinite, q2)
        # genuine oscillation: the quadrupole signal reverses direction (has both
        # a rising and a falling segment), i.e. it is not monotone drift.
        dq = diff(q2)
        @test any(dq .> 0) && any(dq .< 0)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 0.1       # bounded ~one dynamical time
    end
end
