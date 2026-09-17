#=
    test_dgstarhp — the 1D hp-adapted DG star after Hébert–Kidder–Teukolsky (2018), Sec. VI.A
    (src/dg/DGStarHP.jl). Anchor star ShumPolytrope(100), ε_c = 0.00128 + 100·0.00128²
    (R = 9.586 M⊙, the star they evolve). Linear Cowling radial modes for it:
    F = 2.686, H1 = 4.550, H2 = 6.341 kHz (RadialModes.radial_cowling_spectrum).
    Numbers quoted in the comments are from repro/dgstarhp_hkt.jl (2026-09-17).
=#
using Test
using BDNKStar

@testset "DGStarHP — hp-adapted radial DG star (Hébert–Kidder–Teukolsky grids)" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2

    @testset "grids: surface on an element boundary, staggered centre, per-region order" begin
        for preset in (:I1, :I2, :I1R, :I2R, :coarse)
            eng, st = setup_dgstarhp(eos, εc; grid=preset)
            xb = vcat([e.xlo for e in eng.elems], eng.elems[end].xhi)
            @test minimum(abs.(xb .- eng.R)) < 1e-12 * eng.R          # R is an element boundary
            @test minimum(abs.(xb .+ eng.R)) < 1e-12 * eng.R          # and so is −R
            @test all(abs.(eng.x) .> 1e-6)                            # no node at the origin
            c = eng.elems[eng.kc]
            @test c.region == :center && isodd(c.p) && c.xlo == -c.xhi
            @test all(e -> e.p == 1 || e.p == 2, filter(e -> e.region == :surface, eng.elems))
        end
        @test_throws ErrorException setup_dgstarhp(eos, εc; grid=:I2, p_center=4)
    end

    @testset "well-balanced: the projected TOV star is a fixed point (I2, :wb limiter)" begin
        eng, st = setup_dgstarhp(eos, εc; grid=:I2, limiter=:wb)
        rec = evolve_dgstarhp!(st, eng; tmax=300.0, sample_dt=100.0)
        @test maximum(rec.errD) < 1e-9                    # measured 2e-10 … 9e-11 over 4000 M⊙
        @test maximum(abs.(rec.ρc ./ rec.ρc[1] .- 1)) < 1e-9
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-9
        @test maximum(rec.vatm) == 0.0                    # nothing moves in the atmosphere
    end

    @testset "the doubled domain's antisymmetric mode: grows without the parity projection" begin
        # unphysical translation mode, growth 0.049/M⊙ from round-off on every grid and flux
        eng, st = setup_dgstarhp(eos, εc; grid=:I2, limiter=:none, symmetrize=false)
        rec_a = evolve_dgstarhp!(st, eng; tmax=400.0, sample_dt=100.0)
        eng, st = setup_dgstarhp(eos, εc; grid=:I2, limiter=:none, symmetrize=true)
        rec_s = evolve_dgstarhp!(st, eng; tmax=400.0, sample_dt=100.0)
        @test rec_a.errD[end] > 1e3 * rec_s.errD[end]     # measured: 4e-6 vs 2e-13 at t=600
        @test rec_s.errD[end] < 1e-10
        rate = (log(rec_a.errD[end]) - log(rec_a.errD[2])) / (rec_a.t[end] - rec_a.t[2])
        @test 0.02 < rate < 0.1                           # measured 0.049/M⊙
    end

    @testset "seeded radial oscillation: the fundamental mode of linear theory" begin
        F_lin = radial_cowling_spectrum(eos, εc; nmodes=2)[1] ./ BDNKStar.Units.Msun_to_km  # kHz
        @test abs(F_lin[1] - 2.686) < 0.01
        eng, st = setup_dgstarhp(eos, εc; grid=:I2, limiter=:wb)
        seed_dgstarhp_radial!(st, eng; A=1e-3)            # v̂ = A sin(πr/R)
        rec = evolve_dgstarhp!(st, eng; tmax=2000.0, sample_dt=0.5)
        pk, _, _ = dgstarhp_spectrum(rec.t, rec.ρc; npeaks=3)
        @test any(f -> abs(f/F_lin[1] - 1) < 0.02, pk)    # F recovered (4000 M⊙ record: −0.2%)
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-7       # the :wb limiter is conservative
        @test rec.errD[end] < 5e-3
        @test maximum(abs.(st.v)) < 0.5                   # no relativistic atmosphere
    end

    @testset "the paper's own scheme: I1 grid, minmod on the linear surface elements, no subtraction" begin
        # settles to a numerical equilibrium instead of staying on the projection (their Fig. 9:
        # err[D̃] 1e-4 → 7e-4 by 1e4 M⊙; ours 1.0e-3 at 4000); mass is conserved because the
        # limiter acts on the densitized D̃ (theirs loses M_b at the 1e-4 level)
        eng, st = setup_dgstarhp(eos, εc; grid=:I1, limiter=:minmod, wellbalanced=false)
        rec = evolve_dgstarhp!(st, eng; tmax=1000.0, sample_dt=100.0)
        @test all(isfinite, rec.errD) && rec.errD[end] < 5e-3
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-8             # measured 3e-11
        @test abs(rec.ρc[end]/rec.ρc[1] - 1) < 2e-3             # measured −2.3e-4 at 4000
        @test maximum(abs.(st.v[eng.r .> 1.3*eng.R])) < 1e-3    # nothing blown far out
    end

    @testset "minmod slope estimate is the exact linear mode on every element order" begin
        # the LGL-quadrature formula 1.5 Σ w ξ u is 3× too large on p=1 elements; that bug made the
        # limiter rewrite the surface slopes every stage and destroyed the star on the I1 grid
        eng, st = setup_dgstarhp(eos, εc; grid=:I2)
        for e in (eng.elems[eng.kc], eng.elems[eng.kc+1], eng.elems[eng.kc+12])
            U = zeros(eng.Ntot); b = eng.bases[e.p]
            for a in 1:e.N; U[e.off+a-1] = 2.0 + 0.7*b.ξ[a] + (e.p ≥ 2 ? 0.3*(3b.ξ[a]^2-1)/2 : 0.0); end
            @test abs(BDNKStar.DGStarHP._slope(eng, U, e) - 0.7) < 1e-12
        end
        eng, st = setup_dgstarhp(eos, εc; grid=:I1)
        e = eng.elems[findfirst(el -> el.p == 1, eng.elems)]
        U = zeros(eng.Ntot); U[e.off] = 1.3; U[e.off+1] = 2.7
        @test abs(BDNKStar.DGStarHP._slope(eng, U, e) - 0.7) < 1e-12
    end
end
