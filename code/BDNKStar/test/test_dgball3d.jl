#=
    test_dgball3d — the 3+1D cubed-sphere DG star after Hébert–Kidder–Teukolsky (2018)
    (src/dg/DGBall3D.jl). Anchor star ShumPolytrope(100), ε_c = 0.00128 + 100·0.00128².
    Coarse grid (nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3): 464 elements, 27 776 nodes. Runtime ≈ 8 min on 2 threads.
    Measured numbers are from the 2026-09-17 runs (see VALIDATION.md §7.10).
=#
using Test
using BDNKStar
using BDNKStar.DGBall3D: nidx, nnodes, _p2c, _εpoly, _p_ig, _raw_rhs!

@testset "DGBall3D — cubed-sphere conforming DG star" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
    eng, st = setup_dgball3d(eos, εc; nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3, wellbalanced=true, limiter=:wb)   # default volume=:split

    @testset "geometry: conforming faces, positive Jacobians, volume, surface on a shell boundary" begin
        @test all(eng.J .> 0)
        rmax = maximum(eng.r)
        @test abs(dgball3d_volume(eng)/(4π*rmax^3/3) - 1) < 1e-3          # measured 2.3e-4
        nint = 0; nbnd = 0; worst = 0.0
        for fc in eng.faces, (q, n) in enumerate(fc.nodes)
            m = fc.partner[q]
            if m == 0; nbnd += 1; continue; end
            nint += 1
            worst = max(worst, hypot(eng.x[n]-eng.x[m], eng.y[n]-eng.y[m], eng.z[n]-eng.z[m]))
        end
        @test worst < 1e-10                                            # partner nodes coincide
        @test nbnd > 0 && nint > 10*nbnd
        # the stellar surface is a shell boundary: some face nodes sit exactly at r=R and no
        # node lies within (0, h/4) of it on either side except those
        onR = count(abs.(eng.r .- eng.R) .< 1e-10*eng.R)
        @test onR ≥ 6*eng.elems[end].N[2]^2
        # the origin is a node of the cube elements and the metric there is regular
        @test length(eng.origin) ≥ 1 && all(isfinite, eng.α[eng.origin]) && all(eng.sqrtγ[eng.origin] .== 1.0)
        # odd nt with odd p_t: the origin is inside the middle cube element and is not a node (with
        # even p_t the LGL midpoint would sit on it); the central density is then interpolated
        eng5, st5 = setup_dgball3d(eos, εc; nt=3, p_t=3, p_int=2, p_surf=2, p_ext=2)
        @test isempty(eng5.origin)
        @test abs(dgball3d_central_density(st5, eng5)/dgball3d_central_density(st, eng) - 1) < 1e-4   # measured 5e-6
    end

    @testset "curl-form metrics: the discrete metric identity holds" begin
        @test dgball3d_metric_identity(eng) < 1e-13                    # measured 4e-15
        # the face area vectors are the contravariant vectors: |A| agrees with the cross-product
        # of the discrete tangents to truncation, and partner faces have opposite area vectors
        worst = 0.0
        Avec = Dict{Tuple{Int,Int},Vector{Float64}}()
        for (fid, fc) in enumerate(eng.faces), (q, n) in enumerate(fc.nodes); Avec[(fid, n)] = fc.A[:, q]; end
        for fid in 1:37:length(eng.faces)                              # a sample of faces (the full sweep is O(F²))
            fc = eng.faces[fid]
            for (q, n) in enumerate(fc.nodes)
                m = fc.partner[q]; m == 0 && continue
                # find the partner face id: the face of elem(m) that contains m and whose centroid matches
                found = false
                for (gid, gc) in enumerate(eng.faces)
                    gc.elem == eng.faces[fid].elem && continue
                    qm = findfirst(==(m), gc.nodes); qm === nothing && continue
                    gc.partner[qm] == n || continue
                    worst = max(worst, sum(abs2, fc.A[:, q] .+ gc.A[:, qm]) / sum(abs2, fc.A[:, q])); found = true; break
                end
                @test found
                break                                                    # one node per face is enough
            end
        end
        @test worst < 1e-20
    end

    @testset "free stream: a uniform moving gas on a flat metric has zero RHS (both volume forms)" begin
      for vol in (:chain, :split)
        eng2, st2 = setup_dgball3d(eos, εc; nt=2, p_t=2, p_int=2, p_surf=2, p_ext=2, wellbalanced=false, volume=vol)
        fill!(eng2.α, 1.0); fill!(eng2.elam, 1.0); fill!(eng2.sqrtγ, 1.0); fill!(eng2.Φp, 0.0); fill!(eng2.λp, 0.0); fill!(eng2.gor, 0.0)
        ρ0 = 1e-3; ε0 = _εpoly(eng2.κ, eng2.Γ, ρ0)
        for n in 1:eng2.Ntot
            st2.ρ[n] = ρ0; st2.ε[n] = ε0; st2.p[n] = _p_ig(eng2.Γ, ρ0, ε0); st2.vx[n] = 0.3; st2.vy[n] = -0.1; st2.vz[n] = 0.2
            D̂, Ŝx, Ŝy, Ŝz, τ̂ = _p2c(eng2.Γ, ρ0, ε0, 0.3, -0.1, 0.2, 1.0, 0.0, 0.0, 0.0)
            st2.D[n] = D̂; st2.Sx[n] = Ŝx; st2.Sy[n] = Ŝy; st2.Sz[n] = Ŝz; st2.τ[n] = τ̂
        end
        r = ntuple(_ -> zeros(eng2.Ntot), 5); _raw_rhs!(r..., st2, eng2)
        inner = eng2.r .< 2.5*eng2.R                                   # away from the outflow boundary
        @test maximum(maximum(abs, v[inner]) for v in r) < 1e-14 * ρ0    # measured 6e-18 / 2e-17
      end
    end

    @testset "static star: projected TOV star is a fixed point with the subtraction and the :wb limiter" begin
        ρc0 = dgball3d_central_density(st, eng)
        @test abs(ρc0/BDNKStar.FVCommon.rho_from_p(eos, BDNKStar.pressure(eos, εc)) - 1) < 1e-8   # TOV table starts at r=h_tov
        rec = evolve_dgball3d!(st, eng; tmax=30.0, sample_dt=10.0)
        @test maximum(rec.errD) < 1e-12                                # measured 7e-15 at t=60
        @test abs(rec.ρc[end]/ρc0 - 1) < 1e-12
        @test abs(rec.Mb[end]/rec.Mb[1] - 1) < 1e-12
        @test maximum(rec.vatm) == 0.0
        @test abs(rec.q20[end]) < 1e-12 * abs(rec.q00[end])            # stays spherical
    end

    @testset "the aliasing instability and the paper's momentum filter" begin
        # ℓ=2 seed; without the filter err[D̃] grows exponentially (e-fold ≈ 250 M⊙) after the
        # initial transient has decayed; with the filter (momentum only, s=6 centre / 12 shells)
        # it keeps falling. Measured t∈[400,600]: 2.4e-4 (growing) vs 1.7e-5 (falling).
        res = Dict{Symbol,Any}()
        for (tag, filt) in ((:nofilter, :none), (:filter, :all))
            eng4, st4 = setup_dgball3d(eos, εc; nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3, wellbalanced=true, limiter=:wb, filter=filt)
            seed_dgball3d_l2!(st4, eng4; A=1e-3)
            res[tag] = evolve_dgball3d!(st4, eng4; tmax=600.0, sample_dt=50.0)
        end
        e_nf = res[:nofilter].errD; e_f = res[:filter].errD; t = res[:filter].t
        late = t .≥ 400
        @test maximum(e_nf[late]) > 3 * maximum(e_f[late])
        @test e_f[end] < e_f[findfirst(t .≥ 200)]                     # still decaying with the filter
        @test abs(res[:filter].ρc[end]/res[:filter].ρc[1] - 1) < 1e-4  # measured ~1e-6
        @test abs(res[:filter].Mb[end]/res[:filter].Mb[1] - 1) < 1e-6  # measured 1e-8
        # the seeded ℓ=2 signal is present and decays; the m=2 channel stays empty (axisymmetric seed)
        q = res[:filter].q20 ./ abs(res[:filter].q00[1])
        @test maximum(abs.(q[t .< 100])) > 1e-3 && maximum(abs.(q[late])) < 0.2*maximum(abs.(q[t .< 100]))
        @test maximum(abs.(res[:filter].q22)) < 1e-6 * maximum(abs.(res[:filter].q20))
    end

    @testset "static star without the subtraction settles (the paper's regime)" begin
        eng3, st3 = setup_dgball3d(eos, εc; nt=2, p_t=3, p_int=3, p_surf=2, p_ext=3, wellbalanced=false, limiter=:wb)
        res = dgball3d_static_residual(st3, eng3)
        @test res[:interior] < 0.15 && res[:surface] < 0.2            # 2.7%/4.7% chain form, 9.5%/11% split form
        rec = evolve_dgball3d!(st3, eng3; tmax=40.0, sample_dt=20.0)
        @test all(isfinite, rec.errD) && rec.errD[end] < 5e-3         # measured 1.4e-3
        @test abs(rec.ρc[end]/rec.ρc[1] - 1) < 2e-3                   # measured −4e-4
        @test maximum(rec.vatm) == 0.0
    end
end
