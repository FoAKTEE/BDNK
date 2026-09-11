#=
    test_fvcartesian — Stage-1 (radial FV MVP) and Stage-2 (Cartesian FV) tests
    of the NONLINEAR flux-conservative Valencia GRHydro reproduction with an
    atmosphere (HLL + MinMod + SSP-RK2, well-balanced background).

    Validated behaviour (see HONEST ASSESSMENT in the module docs):
      Stage 1 (MVP, SOLID): the static TOV star is EXACTLY well-balanced (zero
        drift); a homologous radial (ℓ=0) perturbation produces a GRID-CONVERGED
        radial fundamental at ≈3.08 kHz (Hann-windowed FFT identical at N=200,400).
      Stage 2 (REPRODUCTION ATTEMPT — fv-insufficient): the Cartesian octant star
        is well-balanced STATICALLY (the atmosphere cures the staircase for the
        unperturbed star), but a seeded ℓ=2 perturbation excites a SECULAR
        surface (staircase) instability whose growth ACCELERATES with resolution,
        so the ℓ=2 f-mode is NOT cleanly reproduced at affordable resolution.
        The tests assert what is true: static stability + short-time boundedness.
=#
using Test
using BDNKStar

@testset "FV Cartesian/radial reproduction" begin
    κ   = 100.0
    εc  = 0.00128 + 100*0.00128^2
    eos = ShumPolytrope(κ)

    @testset "Stage 1: static well-balancing (the staircase cure, 1D)" begin
        eng, st = setup_fvradial(eos, εc; N=200, rmax_fac=1.3, atm_fac=1e-7, cfl=0.3)
        ts, probe, ρc, drift = evolve_fvradial!(st, eng; tmax=200.0, probe_frac=0.5)
        @test drift < 1e-10                       # well-balanced ⇒ machine-zero drift
        @test isfinite(ρc[end])
    end

    @testset "Stage 1: radial mode is GRID-CONVERGED (MVP validation)" begin
        fpk = Float64[]
        for N in (200, 400)
            eng, st = setup_fvradial(eos, εc; N=N, rmax_fac=1.3, atm_fac=1e-7, cfl=0.25)
            seed_radial_velocity!(st, eng; A=1e-4, profile=:linear)
            ts, probe, ρc, drift = evolve_fvradial!(st, eng; tmax=3000.0,
                                                    probe_frac=0.3, sample_dt=2.0)
            @test drift < 1e-2                    # perturbed star stays bounded
            # Hann-windowed central-density FFT → clean single radial fundamental
            fr, _, _ = radial_periodogram_freqs(ts, ρc; fmin_kHz=2.0, fmax_kHz=4.0,
                                                npeaks=1, window=true)
            push!(fpk, fr[1])
        end
        # the fundamental must be grid-converged (the FV-scheme-works proof)
        @test abs(fpk[1]-fpk[2])/fpk[1] < 0.03
        @test 2.5 < fpk[1] < 3.6                  # converged radial fundamental band
    end

    @testset "Stage 2: Cartesian static stability (staircase cure works statically)" begin
        eng, st = setup_fvcart(eos, εc; N=24, L_fac=1.3, atm_fac=1e-7, cfl=0.25)
        ρc0 = fvcart_central_density(st, eng)
        ts, q2, ρc = evolve_fvcart!(st, eng; tmax=40.0, sample_dt=4.0)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 1e-8   # well-balanced static
        @test all(isfinite, ρc)
    end

    @testset "Stage 2: Cartesian ℓ=2 oscillates and is short-time bounded" begin
        eng, st = setup_fvcart(eos, εc; N=24, L_fac=1.3, atm_fac=1e-7, cfl=0.25)
        seed_l2_velocity!(st, eng; A=1e-5)
        ρc0 = fvcart_central_density(st, eng)
        ts, q2, ρc = evolve_fvcart!(st, eng; tmax=120.0, sample_dt=4.0)
        @test all(isfinite, q2)
        @test any(q2 .> 0) && any(q2 .< 0)          # genuine oscillation
        # bounded over ~one dynamical time (the secular staircase drift grows later)
        @test maximum(abs.(ρc .- ρc0))/ρc0 < 0.1
    end
end
