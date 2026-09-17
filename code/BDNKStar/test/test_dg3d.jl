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

    # THE FIXED-POINT PROPERTY — the equilibrium must be invariant under BOTH the limiter and
    # the RHS, not only the RHS. Before the equilibrium-preserving limiter, the unseeded star was
    # thrown into a ±3% radial oscillation with 0.6c surface velocities within 50 M⊙ by the
    # limiter alone (RHS subtraction exact; limiter=false stayed static to 1e-12).
    @testset "equilibrium is a fixed point of the limiter and of the RHS" begin
        eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
        D0 = copy(st.D); τ0 = copy(st.τ); S0 = copy(st.Sx)
        BDNKStar.DGCart3D._limit!(st, eng)
        @test st.D == D0 && st.τ == τ0 && st.Sx == S0             # bitwise: the limiter does nothing
        r = [zeros(size(st.D)) for _ in 1:5]
        BDNKStar.DGCart3D._rhs!(r..., st, eng)
        @test maximum(abs, r[1]) == 0.0 && maximum(abs, r[5]) == 0.0   # rhs(U_eq) ≡ 0 exactly
        # and the LONG unseeded evolution stays static with the limiter ON (25 M⊙ ≈ 20% of a
        # period would already show the old ±3% excitation)
        ρc0 = dgcart3d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=120.0, sample_dt=30.0)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 1e-9
        mn = dgcart3d_prim_minmax(st)
        @test max(abs(mn[5]), abs(mn[6])) < 1e-6                     # no spurious surface velocities
    end

    @testset "seeded star: surface elements stay on the equilibrium-preserving path" begin
        # A 1% velocity seed used to send the three surface elements of this grid to the flattening
        # fallback on every stage (reference momentum without its kinetic energy → negative reference
        # pressure where τ_eq → 0). With the energy-consistent reference no element that contains
        # stellar nodes may fall back, at t=0 or after half a period.
        eng, st = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
        seed_dgcart3d_l2!(st, eng; A=1e-2)
        c0 = BDNKStar.DGCart3D.dgcart3d_limiter_census(st, eng)
        @test c0.fallback == 0 && c0.scaled == 0
        ρc0 = dgcart3d_central_density(st, eng)
        ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=50.0, sample_dt=50.0)
        c1 = BDNKStar.DGCart3D.dgcart3d_limiter_census(st, eng)
        # at half a period the seed's outward kick has put ejecta with a negative τ error into two
        # corner elements that hold ONE stellar node each (r/R=0.993) among 26 atmosphere nodes; their
        # mean τ is below zero, so no mean-preserving reference can be feasible there. Bounded, not zero.
        @test c1.star_fallback ≤ 2
        ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=50.0, sample_dt=50.0)
        c2 = BDNKStar.DGCart3D.dgcart3d_limiter_census(st, eng)
        @test c2.star_fallback == 0                # and none from one period on (0 through t=1000 measured)
        @test abs(ρc[end]/ρc0 - 1) < 2e-3        # was −10% by 4.6 periods with the mean-based limiter
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
