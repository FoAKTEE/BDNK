#=
    test_dyngr — STAGE 2 dynamical-GR 1+1D engine (DynGR1D).

    The FIRST engine in the toolkit that DROPS the Cowling approximation and
    EVOLVES the spacetime coupled to the matter (constrained radial/areal gauge +
    polar slicing). These assertions cover the four decisive dynamical-GR tests,
    each of which a frozen-metric (Cowling) engine cannot do:

      1. static TOV is a genuine stationary solution of the FULL coupled system
         (no lake-at-rest subtraction beyond well-balancing) — tiny bounded drift
         over many dynamical times, no NaN, metric solve reproduces TOV m,α;
      2. the metric GENUINELY RESPONDS: the fundamental radial frequency from the
         dynamical evolution is distinct from (well below) the frozen-metric value,
         and the frozen-metric value tracks the validated relativistic Cowling
         eigenvalue;
      3. the corrected full-GR radial-pulsation eigensolver (KR 2001 convention)
         sits BELOW the Cowling eigensolver, hits the grid-converged fundamental
         ω²≈1.981e-3 km⁻² (F≈2.124 kHz), and lands in the published KR(2001)
         full-GR band — and the dynamical engine F_dyn matches it to ~5%;
      4. a supercritical / inward-kicked star COLLAPSES: collapse of the lapse
         (α_c → 0), 2m/r → 1 at the forming apparent horizon, central density grows.
=#
using Test
using BDNKStar
using BDNKStar: setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!, seed_dyngr_toroidal!,
                solve_tov, ShumPolytrope, dyngr_radial_freq,
                chandrasekhar_radial_omega2, cowling_radial_omega2,
                radial_cowling_spectrum
using BDNKStar.DynGR1D: _solve_metric!, prim2cons_toroidal, cons2prim_toroidal,
                        dyngr_metric!
using BDNKStar.FVCommon: AtmospherePars
using BDNKStar.Units: kHz_to_km

const _conv_kHz = 1/kHz_to_km          # km^-1 -> kHz (geometric, εc in km^-2)

_detrend(ts, y) = begin
    n = length(ts); x = (ts .- ts[1]) ./ (ts[end]-ts[1])
    A = hcat(ones(n), x, x.^2, x.^3); y .- A*(A\y)
end

@testset "DynGR1D — dynamical-GR 1+1D engine (Stage 2)" begin
    eos = ShumPolytrope(100.0)
    εc  = 0.0015
    star = solve_tov(eos, εc; h=2e-4)

    # --- the metric is a STATE that genuinely evolves (not frozen) -------------
    @testset "metric solve reproduces TOV" begin
        eng, st = setup_dyngr(eos, εc; N=400)
        _solve_metric!(st, eng)
        @test isapprox(eng.M, star.M; rtol=1e-3)        # ADM mass from constraint
        NG = eng.g.NG
        # interior α, m vs TOV at mid-radius
        ai = NG + clamp(round(Int, 0.5*star.R/eng.g.Δr + 0.5), 1, eng.g.N)
        r = eng.g.r[ai]
        j = clamp(searchsortedlast(star.r, r), 1, length(star.r)-1)
        t = (r-star.r[j])/(star.r[j+1]-star.r[j])
        mtov = star.m[j]+t*(star.m[j+1]-star.m[j])
        αtov = exp(0.5*(star.ν[j]+t*(star.ν[j+1]-star.ν[j])))
        @test isapprox(st.m[ai], mtov; rtol=2e-3)
        @test isapprox(st.α[ai], αtov; rtol=2e-3)
        @test st.X[ai] > 1.0                              # X = 1/√(1-2m/r) > 1
    end

    # --- 1. static TOV stationarity over many dynamical times -----------------
    @testset "static TOV stationary (full coupled system)" begin
        eng, st = setup_dyngr(eos, εc; N=400, cfl=0.3)
        res = evolve_dyngr!(st, eng; tmax=15*star.R)       # ~15 dynamical times
        @test res.ts[end] ≥ 14*star.R
        @test isfinite(res.ρc[end])
        @test res.drift < 1e-6                             # machine-small drift
        @test !res.collapsed
    end

    # --- 2. full-GR radial mode RESPONDS: dynamical ≠ frozen ------------------
    @testset "metric responds: dynamical radial mode below frozen" begin
        function fmode(fm, band)
            eng, st = setup_dyngr(eos, εc; N=500, cfl=0.25, freeze_metric=fm)
            seed_dyngr_velocity!(st, eng; A=2e-4, profile=:linear)
            res = evolve_dyngr!(st, eng; tmax=220*star.R, probe_frac=0.5,
                                sample_dt=0.5*eng.g.Δr)
            yv = _detrend(res.ts, copy(res.probe))
            fpk,_,_ = dyngr_radial_freq(res.ts, yv;
                        fmin=band[1]*kHz_to_km, fmax=band[2]*kHz_to_km, nf=8000)
            fpk*_conv_kHz
        end
        Fdyn    = fmode(false, (1.2, 4.0))      # dynamical metric
        Ffrozen = fmode(true,  (3.0, 6.5))      # frozen-metric control (≈ Cowling)
        Fcow = radial_cowling_spectrum(eos, εc; N=1500, h_tov=5e-5, nmodes=1)[1][1]
        # corrected SL eigenvalue (the GR ground truth the engine must match)
        ω2gr1 = chandrasekhar_radial_omega2(eos, εc; nmodes=1, N=2500)[1]
        Fgr   = sqrt(ω2gr1)/(2π)*_conv_kHz       # ≈ 2.124 kHz
        @test isfinite(Fdyn) && isfinite(Ffrozen)
        @test 1.4 < Fdyn < 2.6                  # dynamical fundamental ~2.13 kHz
        @test 3.6 < Ffrozen < 4.4               # frozen 4.00 kHz (= the Cowling eigenvalue)
        @test Fdyn < 0.7*Ffrozen                # metric response LOWERS the mode
        @test isapprox(Ffrozen, Fcow; rtol=0.03)# frozen tracks validated Cowling: −0.26%
        # GATE: with the well-balanced operator (VALIDATION §7.15) the engine reproduces BOTH
        # independent frequency-domain eigensolvers: F_dyn = 2.1304 vs the Chandrasekhar
        # full-GR 2.1236 kHz (+0.32%), F_frozen = 3.9994 vs the Cowling 4.0099 (−0.26%).
        # Before the fix the same runs gave 2.02 (−4.9%) and 4.6 (+15%).
        @test isapprox(Fdyn, Fgr; rtol=0.03)
    end

    # --- 3. full-GR and Cowling eigensolvers consistent -----------------------
    @testset "radial-pulsation eigensolvers" begin
        ω2gr  = chandrasekhar_radial_omega2(eos, εc; nmodes=2, N=2500)
        ω2cow = cowling_radial_omega2(eos, εc; nmodes=2)
        @test all(isfinite, ω2gr)
        @test all(isfinite, ω2cow)
        @test ω2gr[1]  > 0                       # stable star: fundamental ω²>0
        @test ω2cow[1] > 0
        @test ω2gr[2]  > ω2gr[1]                 # overtone above fundamental
        # corrected KR (2001) convention: full-GR fundamental is SOFTER than the
        # frozen-metric Cowling value (the spacetime response lowers ω²), matching
        # the dynamical engine (F_dyn ≈ 2.13 < F_frozen ≈ 4.00 kHz).
        @test ω2gr[1]  < ω2cow[1]
        # VALUE LOCK: the corrected fundamental (KR(2001) eqs.14-17, our convention)
        # is ω²≈1.981e-3 km⁻² → F≈2.124 kHz; the swapped-exponent bug gave ω²>ω²_Cow.
        @test isapprox(ω2gr[1], 1.981e-3; rtol=2e-2)
        Fgr = sqrt(ω2gr[1])/(2π)*_conv_kHz
        @test isapprox(Fgr, 2.124; atol=0.05)    # grid-converged fundamental
        # PUBLISHED ANCHOR: Kokkotas & Ruoff (2001) A&A 366,565 [gr-qc/0011093]
        # Table A.18, n=1 κ=100 km² polytrope (p=κρ², Γ=2) — OUR EOS family. Their
        # full-GR fundamentals for nearby stars span ν0≈2.15-2.32 kHz (M≈0.80-1.13
        # M⊙); our εc=0.0015 star (M≈0.96 M⊙) lands in that published band.
        @test 2.0 < Fgr < 2.45                   # within KR(2001) full-GR band
    end

    # --- 4. collapse to a black hole: lapse collapse + 2m/r → 1 ----------------
    # The collapsing star must be one that is ACTUALLY unstable. Until the operator was made
    # well-balanced (VALIDATION §7.15) this test used ρ_c = 2.4162×10⁻³ (ε_c = 0.003), which is
    # on the STABLE branch — the maximum mass of this EOS sits at ρ_c ≈ 3.16×10⁻³ — and it
    # "collapsed" only because the 5%-of-gravity momentum imbalance acted as a steady inward
    # force. With the corrected operator that star merely oscillates (ρ_c settles 8% above its
    # initial value, α_c = 0.51, max 2m/r = 0.44), which is the physical answer. The genuine
    # collapse test is the UNSTABLE-branch star of Font et al. 2002 (ρ_c = 7.993×10⁻³, past the
    # turning point, M = 1.448, R = 5.838) with an inward kick; it now collapses identically
    # with and without the equilibrium subtraction (ρ_c,end 1.614 vs 1.613×10⁻², max 2m/r 0.9595
    # vs 0.9598), where before the fix the two settings disagreed by a factor two in t_AH.
    @testset "collapse: lapse collapse + apparent horizon" begin
        ρc_c = 7.993e-3; εc_c = ρc_c + 100*ρc_c^2   # unstable branch (Font et al. 2002)
        starc = solve_tov(eos, εc_c; h=2e-4)
        eng, st = setup_dyngr(eos, εc_c; N=800, cfl=0.3, rmax_fac=2.0, atm_vbc=:outflow)
        seed_dyngr_velocity!(st, eng; A=-0.01, profile=:cubic)    # inward kick
        res = evolve_dyngr!(st, eng; tmax=120*starc.R, probe_frac=0.5,
                            sample_dt=eng.g.Δr)
        @test res.collapsed                                       # collapse flagged
        @test res.αc[end] < 1e-2                                  # lapse collapse
        @test res.max2mor[end] > 0.9                              # apparent horizon
        @test res.ρc[end] > 2*res.ρc[1]                           # central ρ grows

        # control 1: the STABLE-branch star takes the same kick and does NOT collapse
        eng2, st2 = setup_dyngr(eos, 0.003; N=800, cfl=0.3, rmax_fac=2.0, atm_vbc=:outflow)
        seed_dyngr_velocity!(st2, eng2; A=-0.03, profile=:cubic)
        res2 = evolve_dyngr!(st2, eng2; tmax=120*star.R, probe_frac=0.5)
        @test !res2.collapsed
        @test res2.max2mor[end] < 0.6
        # control 2: the SAME unstable star kicked OUTWARD migrates to the stable branch
        # instead of collapsing — the sign of the kick decides, which is the physics
        eng3, st3 = setup_dyngr(eos, εc_c; N=1600, cfl=0.3, rmax_fac=4.0, atm_vbc=:outflow)
        seed_dyngr_velocity!(st3, eng3; A=+0.05, profile=:cubic)
        res3 = evolve_dyngr!(st3, eng3; tmax=600.0, probe_frac=0.5, sample_dt=0.5)
        @test !res3.collapsed
        late = [res3.ρc[i] for i in eachindex(res3.ts) if res3.ts[i] > 300.0]
        @test 0.5e-3 < sum(late)/length(late) < 2.5e-3            # oscillates about the
        @test minimum(res3.ρc) > 1e-5                             # stable-branch star, ρ_c ≈ 1.3e-3
    end

    # ── TOROIDAL-FIELD GRMHD COUPLING (Stage-2 MHD) ───────────────────────────
    @testset "toroidal-field GRMHD: cons2prim + energy gravitates + induction" begin
        eos = ShumPolytrope(100.0); εc = 0.0015
        atm = AtmospherePars(1e-8, 1e-12, 1e-8, 5e-8, 0.999)
        # (1) magnetic cons2prim round-trip (b²=B²/W²), incl. strong field
        for (ρ,v,B) in ((0.01,0.3,0.05),(0.005,-0.5,0.1),(0.02,0.6,0.2))
            p = eos.κ*ρ^2
            D,S,τ = prim2cons_toroidal(eos, ρ, p, v, B)
            ρr,pr,εr,vr,Wr,isatm = cons2prim_toroidal(eos, D, S, τ, B, atm)
            @test isapprox(ρr, ρ; rtol=1e-8)
            @test isapprox(sign(S)*vr, v; rtol=1e-6, atol=1e-9)
        end
        # (2) STATIC field energy gravitates: ΔM = ∫4πr²·½B²dr exactly
        eng, st = setup_dyngr(eos, εc; N=300)
        M_unmag = eng.M
        M_mag = seed_dyngr_toroidal!(st, eng; B0=0.02, profile=:pressure)  # v=0
        g = eng.g; N=g.N; NG=g.NG
        Emag = sum(g.Δr*4π*g.r[NG+i]^2*0.5*st.B[NG+i]^2 for i in 1:N)
        @test isapprox(M_mag - M_unmag, Emag; rtol=1e-6)     # field energy gravitates
        @test M_mag > M_unmag                                # magnetised star heavier
        # (3) INDUCTION: evolve a kicked magnetised star — toroidal flux EXACTLY
        # conserved (pure advection), ADM mass conserved, field amplifies, stable.
        eng2, st2 = setup_dyngr(eos, εc; N=300)
        seed_dyngr_velocity!(st2, eng2; A=-0.01, profile=:linear)   # inward kick
        seed_dyngr_toroidal!(st2, eng2; B0=0.01, profile=:pressure)
        flux0 = sum(st2.B̃[NG+i] for i in 1:N)*g.Δr
        Mfield0 = eng2.M; Bc0 = st2.B[NG+1]
        d = evolve_dyngr!(st2, eng2; tmax=20.0)
        flux1 = sum(st2.B̃[NG+i] for i in 1:N)*g.Δr
        @test !d.collapsed && all(isfinite, st2.B)
        @test isapprox(flux0, flux1; rtol=1e-8)              # toroidal flux conserved
        @test isapprox(Mfield0, eng2.M; rtol=3e-3)           # ADM mass conserved (back-reaction)
        @test st2.B[NG+1] > 1.5*Bc0                          # field amplified by compression
        # (4) BDNK RESISTIVITY: the field decays MORE with larger r_b (resistive
        # diffusion of B), stable under the explicit parabolic dt cap.
        function magE_decay(r_b)
            e,s = setup_dyngr(eos, εc; N=100, freeze_metric=true)
            seed_dyngr_toroidal!(s, e; B0=0.02, profile=:sin, r_resist=r_b)
            mE(stt) = sum(e.g.Δr*4π*e.g.r[e.g.NG+i]^2*0.5*stt.B[e.g.NG+i]^2 for i in 1:e.g.N)
            E0 = mE(s); dd = evolve_dyngr!(s, e; tmax=8.0)
            (E0 - mE(s))/E0, all(isfinite, s.B) && !dd.collapsed
        end
        d0, ok0 = magE_decay(0.0); d4, ok4 = magE_decay(0.4)
        @test ok0 && ok4                                     # resistive run stable (dt cap)
        @test d4 > d0 + 0.005                                # resistivity ⇒ extra decay (above numerical floor)
    end

    @testset "viscous radial-mode damping (BDNK shear/bulk)" begin
        # Pluck a radial oscillation; viscous momentum diffusion (kinematic ν_visc=η̂)
        # damps it. ν=0 ⇒ ~undamped (numerical floor only); ν>0 ⇒ clear extra decay,
        # stable under the explicit parabolic dt cap dt≤0.4Δr²/(2ν). (The l=0 fundamental
        # is surface-boundary-layer dominated, γ∼√ν — see repro/viscous_damping.)
        eos = ShumPolytrope(100.0); εc = 0.0015
        star = solve_tov(eos, εc; h=2e-4); R = star.R
        function visc_decay(ν)
            e,s = setup_dyngr(eos, εc; N=200, cfl=0.25, ν_visc=ν, freeze_metric=true)
            seed_dyngr_velocity!(s, e; A=2e-4, profile=:linear)
            res = evolve_dyngr!(s, e; tmax=60*R, probe_frac=0.5)
            n=length(res.probe); rms(x)=sqrt(sum(abs2, x .- sum(x)/length(x))/length(x))
            (rms(res.probe[3n÷4:end])/rms(res.probe[1:n÷4]),
             all(isfinite,res.probe) && !res.collapsed)
        end
        r0, ok0 = visc_decay(0.0); rv, okv = visc_decay(0.06)
        @test ok0 && okv                                     # both runs stable (parabolic dt cap)
        @test r0 > 0.8                                       # inviscid: ~undamped (numerical floor only)
        @test rv < 0.4                                       # viscosity ⇒ oscillation clearly damped
        @test rv < 0.5*r0                                    # viscosity ⇒ substantially more decay than ν=0
    end
end
