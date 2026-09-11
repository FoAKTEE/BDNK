using Test
using BDNKStar
using BDNKStar.PiecewisePolytrope: PUBLISHED

# ─────────────────────────────────────────────────────────────────────────────
# GENERAL / REALISTIC neutron-star EOS via the Read–Lackey–Owen–Friedman (2009)
# piecewise-polytrope parametrization [arXiv:0812.2163, Phys. Rev. D 79, 124032].
# A fixed SLy crust (Read Table II) joined to 3 high-density polytropes with the
# fixed dividing densities ρ1=10^14.7, ρ2=10^15.0 g/cm³; each EOS fixed by
# (log10 p1, Γ1, Γ2, Γ3) from Read Table III. Validation gate: TOV M_max & R_1.4
# vs the published Read 2009 values (cross-checked Lackey & Wade 2015 Table I).
# ─────────────────────────────────────────────────────────────────────────────

const _GCGS = BDNKStar.Units.gram_per_cm3_to_km_minus2   # g/cm³ → km⁻²
const _PCGS = BDNKStar.Units.dyne_per_cm2_to_km_minus2   # dyn/cm² → km⁻²

# Coarse TOV M–R scan over central ENERGY density [km⁻²]; returns the sequence
# and (M_max, R at M_max, R_1.4) extracted on the stable (rising) branch.
function _mr_scan(eos; εlo=4e14*_GCGS, εhi=4e15*_GCGS, n=36, h=0.012)
    εs = exp.(range(log(εlo), log(εhi); length=n))
    Rs = Float64[]; Ms = Float64[]; εok = Float64[]
    for εc in εs
        try
            st = solve_tov(eos, εc; h=h, rmax=60.0)
            push!(Rs, st.R); push!(Ms, mass_solar(st)); push!(εok, εc)
        catch
        end
    end
    imax = argmax(Ms)
    R14 = NaN
    for k in 1:imax-1
        if (Ms[k]-1.4)*(Ms[k+1]-1.4) ≤ 0
            t = (1.4 - Ms[k])/(Ms[k+1]-Ms[k]); R14 = Rs[k] + t*(Rs[k+1]-Rs[k]); break
        end
    end
    return εok, Rs, Ms, Ms[imax], Rs[imax], R14
end

@testset "General/realistic EOS — Read et al. 2009 piecewise polytropes" begin

    # ── 1. The presets build and reproduce p(ρ1)=10^logp1 dyn/cm² exactly ─────
    @testset "construction: p1 reproduced, named presets" begin
        ρ1g = 10.0^14.7 * _GCGS
        for (sym, logp1) in [(:SLy,34.384),(:APR4,34.269),(:H4,34.669),(:MS1,34.858)]
            eos = piecewise_polytrope(sym)
            ptarget = 10.0^logp1 * _PCGS
            # p at ρ1: round-trip through ε(p) then p(ε)
            pgot = pressure(eos, energy_from_pressure(eos, ptarget))
            @test isapprox(pgot, ptarget; rtol=1e-3)
        end
        @test_throws ErrorException piecewise_polytrope(:NOPE)
        # keyword builder agrees with the named preset
        e1 = piecewise_polytrope(:H4)
        e2 = piecewise_polytrope(; logp1=34.669, Γ1=2.909, Γ2=2.246, Γ3=2.144)
        ptest = 10.0^34.0 * _PCGS
        @test isapprox(energy_from_pressure(e1, ptest), energy_from_pressure(e2, ptest); rtol=1e-10)
    end

    # ── 2. Valid causal barotrope: p>0, monotone, 0<c_s²<1 up to ρ ≈ 10^15 ────
    # Honest caveat: the Read 2009 PP parametrization is acausal (c_s²>1) ABOVE
    # some density for the high-density-stiff fits (Read et al. note this; LAL
    # tracks `hMinAcausal`). It IS causal up to the central density of a 1.4 M⊙
    # star (≈10^15 g/cm³), the regime used for typical-star physics / the f-mode.
    @testset "causal barotrope up to 10^15 g/cm³ (p>0, monotone, 0<c_s²<1)" begin
        for sym in (:SLy, :APR4, :H4, :MS1)
            eos = piecewise_polytrope(sym)
            es  = exp.(range(log(1e13*_GCGS), log(1e15*_GCGS); length=200))
            ps  = [pressure(eos, e) for e in es]
            cs2 = [sound_speed2(eos, e) for e in es]
            @test all(ps .> 0)                         # positive pressure
            @test issorted(ps)                          # p monotone increasing in ε
            @test all(0 .< cs2)                          # genuinely a barotrope
            @test all(cs2 .< 1 + 1e-9)                   # subluminal in the typical-NS range
            @test all(is_thermodynamically_valid(eos, e) for e in es)
        end
    end

    # ── 2b. HONEST high-density acausality of the parametrization (documented) ─
    # The soft (low-radius) fits SLy/APR4 turn superluminal near the max-mass
    # central density; the stiff H4/MS1 stay causal everywhere. Report the onset.
    @testset "high-density acausality (Read 2009 parametrization caveat)" begin
        for sym in (:SLy, :APR4, :H4, :MS1)
            eos = piecewise_polytrope(sym)
            es  = exp.(range(log(1e14*_GCGS), log(5e15*_GCGS); length=400))
            cs2 = [sound_speed2(eos, e) for e in es]
            idx = findfirst(>(1.0), cs2)
            ρon = idx === nothing ? Inf : es[idx]/_GCGS
            @info "PP causality onset" EOS=sym max_cs2=maximum(cs2) acausal_above_g_cm3=ρon
            # H4 & MS1 are causal across the whole range; SLy & APR4 are not
            if sym in (:H4, :MS1)
                @test maximum(cs2) ≤ 1 + 1e-9
            else
                @test ρon > 1.0e15            # acausal ONLY above ~10^15 g/cm³
            end
        end
    end

    # ── 3. VALIDATION GATE: M_max & R_1.4 vs published Read et al. 2009 ───────
    @testset "M_max & R_1.4 vs Read 2009 (Lackey-Wade Table I)" begin
        results = Dict{Symbol,Any}()
        for sym in (:SLy, :APR4, :H4, :MS1)
            eos = piecewise_polytrope(sym)
            _, _, _, Mmax, RatMmax, R14 = _mr_scan(eos)
            pub = PUBLISHED[sym]
            dM = 100*(Mmax-pub.M_max)/pub.M_max
            dR = 100*(R14 -pub.R_14 )/pub.R_14
            results[sym] = (Mmax=Mmax, R14=R14, dM=dM, dR=dR)
            @info "Read 2009 PP validation" EOS=sym Mmax_got=Mmax Mmax_pub=pub.M_max pctM=dM R14_got=R14 R14_pub=pub.R_14 pctR=dR
        end
        # require ≥2 presets within 5% on BOTH M_max and R_1.4 (gate).
        # In practice all four agree to <2% with the coarse scan; assert all.
        for sym in (:SLy, :APR4, :H4, :MS1)
            r = results[sym]
            @test abs(r.dM) ≤ 5.0
            @test abs(r.dR) ≤ 5.0
        end
        # ordering sanity: MS1 (stiff) is the largest/most massive; APR4 (soft) compact
        @test results[:MS1].R14 > results[:H4].R14 > results[:SLy].R14
        @test results[:MS1].R14 > results[:APR4].R14
        @test results[:MS1].Mmax > results[:SLy].Mmax
    end

    # ── 4. Pipeline demo: SLy 1.4 M⊙ ℓ=2 Cowling f-mode (finite, in-band) ─────
    @testset "SLy 1.4 M⊙ ℓ=2 Cowling f-mode (realistic EOS pipeline)" begin
        eos = piecewise_polytrope(:SLy)
        εok, Rs, Ms, Mmax, _, _ = _mr_scan(eos; n=60)
        imax = argmax(Ms); ε14 = NaN
        for k in 1:imax-1
            if (Ms[k]-1.4)*(Ms[k+1]-1.4) ≤ 0
                t=(1.4-Ms[k])/(Ms[k+1]-Ms[k]); ε14 = εok[k]*(εok[k+1]/εok[k])^t; break
            end
        end
        @test isfinite(ε14)
        st = solve_tov(eos, ε14; h=0.005, rmax=60.0)
        @test isapprox(mass_solar(st), 1.4; atol=0.02)
        # km geometric units ⇒ Lunit_km = 1.0; scan ω² in km⁻² (f-mode ≈ 2-2.5 kHz)
        f, ω2, R = nonradial_cowling_spectrum(eos, ε14; l=2, nmodes=3, N=5000,
                                              nscan=700, Lunit_km=1.0,
                                              ω2lo=1e-4, ω2hi=0.06)
        @info "SLy 1.4 M⊙ Cowling f-mode" f_kHz=f[1] R_km=R M_Msun=mass_solar(st)
        @test length(f) ≥ 1
        @test all(isfinite, f)
        @test all(ω2 .> 0)                               # real ω² ⇒ stable
        @test issorted(f)                                # f < p1 < p2 …
        @test 1.5 ≤ f[1] ≤ 3.0                            # realistic-EOS Cowling f-mode in-band
    end
end
