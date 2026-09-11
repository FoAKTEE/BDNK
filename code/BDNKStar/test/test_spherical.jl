using Test
using BDNKStar
using BDNKStar: Units

# STAGE 3 — boundary-conforming (r,θ) spherical engine. The surface sits on the
# coordinate line r=R (no staircase), which is what makes the full-frame BDNK velocity
# recovery stable on the FULL star. Three layers: background metric, the ideal Cowling
# f-mode, and the full-frame BDNK viscous dissipation.
kHz2cyc(f) = f * Units.Msun_to_km * Units.kHz_to_km

@testset "Boundary-conforming (r,θ) spherical engine" begin
    eos = ShumPolytrope(100.0); εc = 0.00128 + 100 * 0.00128^2

    @testset "Gate 1 — (r,θ) TOV background" begin
        s = build_sphstar(eos, εc; Nr=64, Nθ=12)
        @test isapprox(s.ρ0[1], 0.00128; rtol=2e-2)              # central rest-mass density
        @test s.M > 0 && s.R > 0
        @test isapprox(s.α[end], sqrt(1 - 2*s.M/s.R); rtol=2e-2) # surface lapse = Schwarzschild
        @test all(>(0), s.cs2[1:end-1])                          # sound speed real in the interior
    end

    @testset "Gate 2-4 — ideal Cowling ℓ=2 f-mode (stable, near eigenmode)" begin
        s = build_sphstar(eos, εc; Nr=64, Nθ=12)
        e = setup_sphevo(s; σ_ko=0.02)
        st = SphState(s.grid.Nr, s.grid.Nθ); seed_sph_l2!(st, e; A=1e-3)
        dt = 0.25*s.grid.dr; nst = round(Int, 3000/dt)
        ts, q2, _ = evolve_sph!(st, e; dt=dt, nsteps=nst, sample=max(1,nst÷1200))
        @test all(isfinite, q2)
        @test maximum(abs, q2) < 5*abs(q2[1])                   # bounded full-star evolution
        fb = range(1.4, 2.4; length=1500); P = periodogram(ts, q2, fb .* kHz2cyc(1.0))
        f = fb[argmax(P)]
        @test 1.75 < f < 2.15                                   # ℓ=2 f-mode (eigensolver 1.883 kHz)
    end

    @testset "Gate 5 — full-frame BDNK: shear viscosity dissipates mode energy ∝ η̂" begin
        s = build_sphstar(eos, εc; Nr=48, Nθ=10)
        # energy ratio E_end/E0 for a quadrupole overtone seed; the faithful viscous
        # diagnostic is the quadratic ENERGY, not amplitude/quadrupole decay.
        function efin(η̂)
            ev = setup_sphbdnk(s; η̂=η̂, σ_ko=0.02)              # physical frame: all transport ∝ η̂
            st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, ev, 3; A=1e-3)
            _, _, en = evolve_sphbdnk!(st, ev; dt=0.005, nsteps=round(Int,300/0.005), sample=40)
            (en[end]/en[1], all(isfinite, en), maximum(en)/en[1])
        end
        r0, ok0, _   = efin(0.0)
        r4, ok4, m4  = efin(0.04)
        r8, ok8, m8  = efin(0.08)
        @test ok0 && ok4 && ok8                                 # all runs finite
        @test r8 < r4 < r0                                      # dissipation INCREASES with η̂
        @test r8 < 1.0                                          # net energy DECAY at η̂=0.08
        @test m8 < 2.0                                          # viscosity keeps it bounded
    end

    @testset "polar f+p spectrum vs eigensolver (ℓ=2,4)" begin
        # A broadband pluck on the stable (r,θ) engine reproduces the WHOLE polar spectrum
        # (f + p₁..p₄), for several ℓ. The benchmark is the frequency-domain eigensolver
        # itself (no hard-coded targets). p-modes have radial nodes ⇒ excite broadband and
        # read several point probes; the pole BCs are ℓ-independent (m=0), so only the
        # Legendre angular factor P_ℓ(cosθ) changes with ℓ.
        kHz2cyc(f) = f*Units.Msun_to_km*Units.kHz_to_km
        Pl(l,x) = l==2 ? (3x^2-1)/2 : (35x^4-30x^2+3)/8
        function td_match(l, Nθ; Nr=64, T=1800.0)
            # reference = eigensolver (no hard-coded targets); trimmed N here, full-accuracy
            # version is gated in test_nonradial.
            bench,_,_ = nonradial_cowling_spectrum(eos, εc; l=l, nmodes=5, N=3000, nscan=700)
            s = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ); ev = setup_sphevo(s; σ_ko=0.008)
            st = SphState(Nr, Nθ)
            for jj in 1:Nθ, ii in 1:Nr
                g = exp(-((s.grid.r[ii]-0.5*s.R)/(0.13*s.R))^2)
                st.δε[ii+1,jj+1] = 1e-3*s.ε0[ii]*g*Pl(l, s.grid.cosθ[jj])
            end
            dt = 0.25*s.grid.dr; nst = round(Int, T/dt)
            probes = [(round(Int,f*Nr)+1, 3) for f in (0.25,0.45,0.65,0.8)]
            k1=SphState(Nr,Nθ);k2=SphState(Nr,Nθ);k3=SphState(Nr,Nθ);k4=SphState(Nr,Nθ);tmp=SphState(Nr,Nθ)
            ts=Float64[]; prs=[Float64[] for _ in probes]; samp=max(1,round(Int,1.0/dt))
            for n in 0:nst
                if n%samp==0
                    push!(ts, n*dt)
                    for (pk,(ip,jp)) in enumerate(probes); push!(prs[pk], st.δε[ip,jp]); end
                end
                n==nst && break
                BDNKStar.SphEvolve._rk4!(st, ev, dt, k1,k2,k3,k4, tmp)
            end
            fb = collect(range(1.2, 1.15*maximum(bench); length=6000)); P = zeros(length(fb))
            for p in prs; P .+= periodogram(ts, p .- sum(p)/length(p), fb .* kHz2cyc(1.0)); end
            loc=Int[]; mx=maximum(P)
            for i in 2:length(P)-1; (P[i]>P[i-1] && P[i]>=P[i+1] && P[i]>0.01*mx) && push!(loc,i); end
            peaks = sort(fb[loc])
            (bench, [peaks[argmin(abs.(peaks .- b))] for b in bench], all(p->all(isfinite,p), prs))
        end
        for (l,Nθ) in ((2,10),(4,16))
            bench, matched, fin = td_match(l, Nθ)
            @test fin                                            # stable evolution
            @test issorted(matched)                             # f < p₁ < p₂ < p₃ < p₄ ordering
            for m in 1:5
                @test abs(matched[m]-bench[m])/bench[m] < 0.06   # each mode within 6%
            end
        end
    end
end
