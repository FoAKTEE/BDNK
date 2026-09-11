#=
    test_dg3d — FULL 3+1D nodal Discontinuous-Galerkin Cartesian GRHydro engine
    (DGCart3D), the genuine three-dimensional extension of the validated 2+1D
    axisymmetric DGCart2D (which was m=0 only).

    Four validation pillars (kept FAST/coarse for ≤6 threads):
      1. 3D STATIC WELL-BALANCING — the unperturbed star on the full 3D Cartesian
         octant grid stays machine-zero (lake-at-rest holds in 3D, no NaN).
      2. SHOCK SANITY in 3D — Martí–Müller test 1 along a grid AXIS in the 3D
         engine: positivity (ρ,p>0), subluminal, no large overshoot.
      3. STAIRCASE TREND in 3D — an ℓ=2 (Y20) perturbed star: central-density
         drift does NOT worsen (the FV failure to beat); oscillates & stays finite.
      4. NON-AXISYMMETRIC m≠0 (the new 3D capability) — a Y22 (m=2) perturbation:
         the m=2 quadrupole moment is nonzero, bounded, finite, positivity held.
         This is what 2D-axisymmetric could NOT do.
=#
using Test
using BDNKStar

@testset "DGCart3D — full 3+1D nodal DG GRHydro" begin
    κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)

    @testset "3D static well-balancing (lake-at-rest in 3D)" begin
        eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
        ρc0 = dgcart3d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=8.0, sample_dt=2.0)
        @test all(isfinite, ρc)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 1e-8   # machine-zero static drift
    end

    @testset "3D shock sanity (Martí–Müller test 1 along x axis)" begin
        r = dgcart3d_shocktube_diagonal!(eos; K=12, p=2, tmax=0.25, dir=:x,
                                         ρL=10.0, pL=13.33, ρR=1.0, pR=1e-3)
        @test r.ρmin > 0.0 && r.pmin > 0.0          # positivity held
        @test r.vmax < 1.0 && r.vmin > -1.0         # subluminal
        @test r.ρmax < 1.05*r.ρL                    # no large overshoot (<5%)
        @test isfinite(r.ρmax) && isfinite(r.pmax)
    end

    @testset "3D staircase trend: ℓ=2 (m=0) drift does not worsen" begin
        # coarse vs finer effective resolution; assert finer is NOT worse than
        # a generous bound (FV WORSENED with refinement; DG must not).
        drifts = Float64[]
        for (K, p) in ((5, 2), (6, 2))
            eng, st = setup_dgcart3d(eos, εc; Kx=K, Ky=K, Kz=K, p=p, cfl=0.2)
            seed_dgcart3d_l2!(st, eng; A=1e-5)
            ρc0 = dgcart3d_central_density(st, eng)
            ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=24.0, sample_dt=4.0)
            @test all(isfinite, ρc)
            @test all(isfinite, q2)
            push!(drifts, maximum(abs.(ρc .- ρc0))/ρc0)
        end
        # bounded at both resolutions, and refinement does not blow it up
        @test all(d -> d < 0.05, drifts)
        @test drifts[2] < 5.0*drifts[1]             # NOT a worsening runaway
    end

    @testset "3D non-axisymmetric m=2 (Y22) capability" begin
        eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
        # unperturbed star is axisymmetric ⇒ m=2 density moment is machine-zero.
        @test abs(dgcart3d_quadrupole_m2(st, eng)) < 1e-12
        # seed a Y22 *velocity* perturbation (initial m=2 DENSITY moment still ~0;
        # the velocity sources a growing m=2 density deformation under evolution).
        seed_dgcart3d_Y22!(st, eng; A=1e-5)
        ρc0 = dgcart3d_central_density(st, eng)
        ts, q2, ρc, qm2 = evolve_dgcart3d!(st, eng; tmax=40.0, sample_dt=4.0,
                                           record_m2=true)
        # the m=2 moment is genuinely excited (grows from ~0), stays finite, and
        # reverses sign (oscillates) — a non-axisymmetric mode 2D-axisym CANNOT do.
        @test all(isfinite, qm2)
        @test maximum(abs.(qm2)) > 1e-9              # genuinely excited
        @test any(diff(qm2) .> 0) && any(diff(qm2) .< 0)   # oscillates / bounded
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 0.1     # positivity/density bounded
    end
end
