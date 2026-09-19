#=
    test_dgstarfd — the DG/finite-difference subcell machinery (src/dg/DGSubcell.jl) and the
    hybrid radial Cowling star (src/dg/DGStarFD.jl), the scheme of Deppe et al. 2021 /
    SpECTRE `src/Evolution/DgSubcell/`. Anchor star ShumPolytrope(100), ε_c = 0.00128+100·0.00128²
    (R = 9.586 M⊙); linear radial Cowling modes F = 2.686, H1 = 4.550, H2 = 6.341 kHz.
    Numbers in the comments are from repro/dgstarfd_hybrid.jl (2026-09-19).
=#
using Test
using BDNKStar
using BDNKStar.DGCommon: build_lgl_basis
using BDNKStar.DGSubcell
using LinearAlgebra: opnorm, dot, I

@testset "DGSubcell — conservative DG ↔ subcell transfer and the troubled-cell indicators" begin
    @testset "projection/reconstruction pair: R∘P = identity, P conservative" begin
        for p in (1, 2, 3, 4, 5, 7)
            b = build_lgl_basis(p); ops = subcell_ops(b)
            @test ops.M == 2*ops.N - 1                       # SpECTRE's 2N−1 subcells
            @test opnorm(ops.R*ops.P - I) < 1e-12            # measured ≤ 9e-16
            @test all(abs.(sum(ops.P, dims=2) .- 1) .< 1e-13) # a constant projects to itself
            # the projection carries the element integral exactly, for any polynomial the DG
            # basis can represent
            u = [sum(b.ξ[j]^m for m in 0:p) for j in 1:b.N]
            @test abs(ops.h*sum(ops.P*u) - dot(b.w, u)) < 1e-12*max(abs(dot(b.w, u)), 1.0)
            @test maximum(abs, ops.R*(ops.P*u) - u) < 1e-12
        end
    end

    @testset "Persson indicator separates smooth from discontinuous data" begin
        b = build_lgl_basis(3)
        smooth = [exp(0.3*x) for x in b.ξ]
        jump = [x < 0 ? 1.0 : 2.0 for x in b.ξ]
        @test persson_threshold(b) ≈ 3.0^(-4.0)
        @test !persson_troubled(smooth, b)                    # measured ratio 1.3e-3
        @test persson_troubled(jump, b)                       # measured ratio 1.5e-1
        @test persson_ratio(fill(2.0, b.N), b) < 1e-14        # a constant has no high modes
    end

    @testset "relaxed discrete maximum principle and the MC slope" begin
        @test !rdmp_troubled(0.9, 1.1, 0.9, 1.1)
        @test rdmp_troubled(0.5, 1.1, 0.9, 1.1)               # candidate dips below the range
        @test !rdmp_troubled(0.9 - 5e-5, 1.1, 0.9, 1.1)       # ... but δ₀ = 1e-4 tolerates this
        @test mc_slope(1.0, 2.0, 3.0) ≈ 1.0                   # linear data: exact slope
        @test mc_slope(1.0, 2.0, 1.0) == 0.0                  # extremum: flattened
        @test mc_slope(1.0, 1.0, 5.0) == 0.0                  # one-sided jump: no overshoot
    end
end

@testset "DGStarFD — DG/FD hybrid radial star" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2

    @testset "the indicator selects the stellar surface, not the atmosphere" begin
        eng, st = setup_dgstarfd(eos, εc; grid=:I2)
        @test eng.Mtot > eng.Ntot                              # 2N−1 subcells per element
        @test dgstarfd_errD(st, eng) == 0.0                    # equilibrium reference is the initial state
        fd = [k for k in 1:eng.K if eng.active[k] == 1]
        @test !isempty(fd)
        # every element on the subgrid straddles the stellar surface, and every element that
        # straddles it is on the subgrid
        for k in 1:eng.K
            e = eng.elems[k]
            straddles = min(abs(e.xlo), abs(e.xhi)) < eng.R < max(abs(e.xlo), abs(e.xhi))
            straddles && @test eng.active[k] == 1
        end
        for k in fd
            e = eng.elems[k]
            @test min(abs(e.xlo), abs(e.xhi)) ≤ eng.R*1.05     # they sit at the surface
        end
        # a uniform atmosphere element is smooth and stays on DG (SpECTRE's TciOptions)
        far = findfirst(k -> abs(eng.elems[k].xlo) > 2*eng.R, 1:eng.K)
        @test eng.active[far] == 0
    end

    @testset "transfer between the two representations is exact and conservative" begin
        eng, st = setup_dgstarfd(eos, εc; grid=:I2)
        k = findfirst(k -> eng.active[k] == 0 && eng.elems[k].region == :interior, 1:eng.K)
        e = eng.elems[k]; w = eng.bases[e.p].w
        D0 = [st.D[e.off + a - 1] for a in 1:e.N]
        mass0 = sum(w[a]*e.J*D0[a] for a in 1:e.N)
        BDNKStar.DGStarFD._project_elem!(st, eng, k)
        massfd = sum(eng.shx[k]*st.sD[eng.soff[k] + c - 1] for c in 1:eng.sM[k])
        @test abs(massfd/mass0 - 1) < 1e-12                   # the projection conserves the integral
        BDNKStar.DGStarFD._reconstruct_elem!(st, eng, k)
        @test maximum(abs, [st.D[e.off + a - 1] - D0[a] for a in 1:e.N]) < 1e-10*maximum(abs, D0)
    end

    @testset "the hybrid conserves baryon mass across DG/FD interfaces" begin
        eng, st = setup_dgstarfd(eos, εc; grid=:I2)
        rec = evolve_dgstarfd!(st, eng; tmax=30.0, sample_dt=15.0)
        @test all(isfinite, rec.errD)
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-6            # measured 8e-12 … 4e-6
        @test 0.0 < rec.fd_fraction[end] < 0.5                 # only the surface is on subcells
    end

    @testset "pure finite volume (always_fd) is exactly conservative" begin
        eng, st = setup_dgstarfd(eos, εc; grid=:I2, always_fd=true)
        @test dgstarfd_fd_fraction(eng) == 1.0
        rec = evolve_dgstarfd!(st, eng; tmax=30.0, sample_dt=15.0)
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-11           # measured 2e-13: FV is conservative
        # the subcell centred on x = 0 must carry the central density: √γ = e^{λ/2}x² vanishes
        # there, so the recovery has to divide by the CELL AVERAGE of √γ, not its centre value
        @test abs(dgstarfd_central_density(st, eng)/1.28e-3 - 1) < 1e-2
    end

    @testset "the hybrid survives the seed that destroyed the pure-DG limiter" begin
        # VALIDATION.md §7.9: the I1 grid with the ΛΠ¹ minmod and a v = 1e-3 sin(πr/R) seed blows
        # a wind, the atmosphere reaching |v| → 1 and M_b growing by 7e-4. The hybrid holds.
        eng, st = setup_dgstarfd(eos, εc; grid=:I1)
        seed_dgstarfd_radial!(st, eng; A=1e-3)
        rec = evolve_dgstarfd!(st, eng; tmax=200.0, sample_dt=50.0)
        @test all(isfinite, rec.errD) && all(isfinite, rec.ρc)
        @test maximum(rec.vatm) < 0.5                          # no wind (pure-DG minmod reached 1.0)
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-4
    end

    @testset "seeded oscillation reproduces the radial spectrum of linear theory" begin
        # SpECTRE's configuration: uniform high-order elements, surface inside an element, no
        # limiter and no filter. The full scan (repro/dgstarfd_hybrid.jl) converges to
        # F = 2.685 (−0.04%), H1 = 4.533 (−0.35%) at K = 28 per side.
        F_lin = radial_cowling_spectrum(eos, εc; nmodes=2)[1] ./ BDNKStar.Units.Msun_to_km
        @test abs(F_lin[1] - 2.686) < 0.01
        R = 9.585578
        grid = [(xlo=0.0, xhi=3.0*R, K=14, p=5, region=:interior)]
        eng, st = setup_dgstarfd(eos, εc; grid=grid, p_center=5, center_fac=1)
        seed_dgstarfd_radial!(st, eng; A=1e-3)
        rec = evolve_dgstarfd!(st, eng; tmax=2000.0, sample_dt=0.5)
        pk, _, _ = dgstarfd_spectrum(rec.t, rec.ρc; npeaks=4)
        @test any(f -> abs(f/F_lin[1] - 1) < 0.02, pk)         # measured −0.4% at this resolution
        @test any(f -> abs(f/F_lin[2] - 1) < 0.04, pk)         # H1, measured −1.8%
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-8            # measured 4e-12
        @test maximum(rec.vatm) < 0.5
    end
end
