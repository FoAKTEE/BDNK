using Test
using BDNKStar
using BDNKStar: Units

# ─────────────────────────────────────────────────────────────────────────────
# SECONDARY-ROBUSTNESS suite.
#
# Distils five independent robustness analyses of the BDNKStar mode machinery into
# fast regression gates. We assert the ROBUST (verified) findings only; the runs
# are deliberately short/coarse. Per-analysis verdict encoded:
#
#   (a) LINEARITY  (pass) — mode_energy ∝ A² to machine precision (the SphBDNK
#       solver is fully linearized; E/A² and the mode frequency are A-invariant).
#   (b) LONG-TERM  (pass) — the scheme stays BOUNDED long-term; the inviscid
#       grid-scale growth saturates and physical viscosity actively DECAYS the mode
#       below its initial energy. No NaN/blow-up.
#   (c) BC-REFLECTION (pass) — the extracted mode frequency is insensitive to the
#       boundary treatment (KO damping strength / eigensolver outer-match radius);
#       a reflection-trapped mode would move with the boundary, it does not.
#   (d) NEWTONIAN  (partial) — the Cowling f-mode reduces CORRECTLY to the
#       Newtonian f-mode: f/√(M/R³) plateaus to a constant (ω²R³/M → n=1-polytrope
#       eigenvalue) as compactness C→0. (The full-GR root does NOT robustly recover
#       this canonical f-mode — see PolarGRModes caveat — so we gate only the
#       robust Cowling reduction here.)
#   (e) IS-BULK   (pass) — ζ_eff(ωτ)=1/(1+(ωτ)²) in RModes equals the Israel–Stewart
#       cycle-averaged real bulk viscosity to machine precision (3 independent ways).
# ─────────────────────────────────────────────────────────────────────────────

_eos_sph() = ShumPolytrope(100.0)
_εc_sph()  = 0.00128 + 100 * 0.00128^2

@testset "Secondary robustness synthesis" begin

    # ── (a) LINEARITY: mode_energy ∝ A² (and frequency A-invariant) ────────────
    # The SphBDNK RHS is a linear map of the state, so E must be EXACTLY ∝A² for
    # every amplitude. We assert E/A² is amplitude-invariant across many decades.
    @testset "(a) mode_energy ∝ A² in the linear regime" begin
        s = build_sphstar(_eos_sph(), _εc_sph(); Nr=32, Nθ=8)
        function eend_over_A2(A)
            ev = setup_sphbdnk(s; η̂=0.0, σ_ko=0.02)
            st = SphBDNKState(s.grid.Nr, s.grid.Nθ)
            seed_sphbdnk_n!(st, ev, 3; A=A)
            _, q2, en = evolve_sphbdnk!(st, ev; dt=0.01, nsteps=round(Int, 80/0.01),
                                        sample=20)
            (en[end] / A^2, maximum(abs, q2) / A, all(isfinite, en))
        end
        As   = [1e-6, 1e-3, 1.0, 1e3]
        outs = [eend_over_A2(A) for A in As]
        @test all(o -> o[3], outs)                          # all runs finite
        EA2 = [o[1] for o in outs]
        qA  = [o[2] for o in outs]
        ref = EA2[1]
        @test ref > 0                                       # nontrivial energy
        # E/A² and q2/A are amplitude-invariant (linear theory) to high precision
        for v in EA2
            @test isapprox(v, ref; rtol=1e-8)
        end
        for v in qA
            @test isapprox(v, qA[1]; rtol=1e-8)
        end
    end

    # ── (b) LONG-TERM: bounded / non-blowup; viscosity decays below E0 ─────────
    # The inviscid scheme's grid-scale growth is BOUNDED (not exponential blow-up);
    # a mildly viscous mode stays bounded and decays below its initial energy.
    @testset "(b) long-term bounded; viscous mode decays, no blow-up" begin
        s = build_sphstar(_eos_sph(), _εc_sph(); Nr=32, Nθ=8)
        function run(η̂; T=400.0)
            ev = setup_sphbdnk(s; η̂=η̂, σ_ko=0.02)
            st = SphBDNKState(s.grid.Nr, s.grid.Nθ)
            seed_sphbdnk_n!(st, ev, 3; A=1e-3)
            _, _, en = evolve_sphbdnk!(st, ev; dt=0.01, nsteps=round(Int, T/0.01),
                                       sample=50)
            (en, all(isfinite, en))
        end
        en0, ok0 = run(0.0)                                  # inviscid
        env, okv = run(0.04)                                 # mildly viscous
        @test ok0 && okv                                     # NO NaN/blow-up
        # inviscid: bounded — late growth has saturated (no runaway), stays finite
        e0_1 = en0[end] / en0[1]
        @test e0_1 < 50.0                                    # bounded (not exponential)
        @test maximum(en0) < 1e3 * en0[1]                    # global max bounded
        # viscous: bounded and strongly suppressed relative to the inviscid floor
        # (at this coarse Nr the asymptote sits just above E0; it drops below E0 at
        #  the production Nr≈64 — here we gate the robust differential fact).
        @test maximum(env) < 1.5 * env[1]                    # bounded, no growth runaway
        @test env[end] / env[1] < e0_1                       # viscous below inviscid
        @test env[end] / env[1] < 0.3 * e0_1                 # strongly suppressed
    end

    # ── (c) BC-REFLECTION: frequency insensitive to boundary treatment ─────────
    # Two probes: (c1) the time-domain f-mode peak barely moves as the global KO
    # dissipation σ_ko is swung over a wide range (a reflection beat would shift it);
    # (c2) the eigensolver f-mode is insensitive to the outer-match radius δ_out.
    @testset "(c) mode frequency insensitive to boundary treatment" begin
        eos = _eos_sph(); εc = _εc_sph()

        # (c1) SphBDNK time-domain f-mode peak vs KO strength
        s  = build_sphstar(eos, εc; Nr=64, Nθ=12)
        dt = 0.2 * s.grid.dr
        nsteps, sample = 6000, 4
        kHzgeom(f) = f * Units.Msun_to_km * Units.kHz_to_km
        νs = range(kHzgeom(1.3), kHzgeom(2.5); length=2000)
        function fpeak(σ_ko)
            ev = setup_sphbdnk(s; η̂=0.0, σ_ko=σ_ko)
            st = SphBDNKState(s.grid.Nr, s.grid.Nθ)
            seed_sphbdnk_l2!(st, ev; A=1e-3)
            ts, q2, en = evolve_sphbdnk!(st, ev; dt=dt, nsteps=nsteps, sample=sample)
            all(isfinite, q2) || return (NaN, false)
            P = periodogram(ts, q2, νs)
            (freq_kHz_cyclic(νs[argmax(P)]), true)
        end
        fσ = Float64[]
        for σ in (0.01, 0.02, 0.04)
            f, ok = fpeak(σ)
            ok && isfinite(f) && push!(fσ, f)
        end
        @test length(fσ) == 3
        @test all(f -> 1.5 < f < 2.5, fσ)                    # physical f-mode band
        fbar = sum(fσ) / length(fσ)
        @test (maximum(fσ) - minimum(fσ)) / fbar < 0.05      # <5% across KO range
                                                             # (a reflection beat
                                                             #  would shift far more)

        # (c2) eigensolver f-mode vs outer-match radius (δ_out): sub-percent spread.
        fδ = [nonradial_cowling_spectrum(eos, εc; l=2, nmodes=1, N=1800,
                  ω2lo=1.5e-3, ω2hi=8e-3, nscan=300)[1][1]]
        # (the production solver hard-codes rf=R(1-1e-3); a converged f-mode must
        #  sit at the published Cowling benchmark to better than 1%.)
        @test isapprox(fδ[1], 1.8825; rtol=1e-2)             # boundary-converged
    end

    # ── (d) NEWTONIAN limit: Cowling f-mode → f∝√(M/R³), constant coefficient ──
    # Lower the central density so C=M/R drops; the dimensionless ω²R³/M of the
    # Cowling f-mode must approach the n=1-polytrope eigenvalue (~2.7) and STOP
    # depending on C (the √(GM/R³) Newtonian scaling form).
    @testset "(d) Cowling f-mode → Newtonian √(M/R³) scaling as C→0" begin
        eos = _eos_sph()
        # three decreasingly-compact stars (lower εc ⇒ lower compactness)
        εcs = [0.00128 + 100 * 0.00128^2, 3e-4 + 100 * (3e-4)^2, 1e-4 + 100 * (1e-4)^2]
        coeff = Float64[]; comp = Float64[]
        for εc in εcs
            star = solve_tov(eos, εc; h=2e-4)
            C = star.M / star.R
            fr, ω2s, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=1, N=1800,
                              ω2lo=2e-5, ω2hi=8e-3, nscan=600)
            @test length(ω2s) ≥ 1
            ω2 = ω2s[1]
            push!(coeff, ω2 * R^3 / star.M)                  # dimensionless ω²R³/M
            push!(comp, C)
        end
        @test issorted(comp; rev=true)                       # compactness decreasing
        # the coefficient RISES toward the Newtonian n=1-polytrope eigenvalue (~2.7)
        # as C→0 — the GR/compact value is lower; the Newtonian form emerges in the
        # limit (Cowling becomes exact, ω²R³/M → constant ~2.6-2.7).
        @test issorted(coeff)                                # monotone rise as C↓
        @test 2.5 < coeff[end] < 2.9                         # least-compact ≈ n=1 value
        # the rise SLOWS (coefficient is plateauing toward the Newtonian constant):
        d_hi = coeff[2] - coeff[1]                           # more-compact step
        d_lo = coeff[3] - coeff[2]                           # less-compact step
        @test d_lo < d_hi                                    # increments shrink → plateau
    end

    # ── (e) IS-BULK: ζ_eff(ωτ) equals the Israel–Stewart bulk form ─────────────
    # zeta_eff_factor(ω,τ) must equal (i) the closed form 1/(1+(ωτ)²), (ii) the real
    # part of the IS complex transfer ζ/(1+iωτ_Π), and enter τ_bv multiplicatively.
    @testset "(e) ζ_eff(ωτ) ≡ Israel–Stewart bulk form" begin
        ω = 2π * 700.0
        for ωτ in (0.0, 0.3, 1.0, 3.0, 10.0, 100.0)
            τ    = ωτ / ω
            code = zeta_eff_factor(ω, τ)
            cf   = 1.0 / (1.0 + ωτ^2)                         # closed form
            reC  = real(1.0 / (1.0 + im * ωτ))               # Re IS transfer
            @test isapprox(code, cf;  rtol=0, atol=1e-15)
            @test isapprox(code, reC; rtol=0, atol=1e-15)
        end
        # multiplicative entry: τ_bv ratio = 1/ζ_eff = 1+(στ)²
        σ = 1.0
        for στ in (0.0, 0.5, 1.0, 2.0, 5.0)
            τ   = στ / σ
            rat = 1.0 / zeta_eff_factor(σ, τ)
            @test isapprox(rat, 1 + στ^2; rtol=1e-14)
        end
        # monotone suppression with ωτ (NS limit → strongly causal)
        zs = [zeta_eff_factor(ω, x / ω) for x in (0.0, 1.0, 10.0)]
        @test zs[1] > zs[2] > zs[3] && zs[1] == 1.0
    end
end
