#=
    bdnk_phase_diagram.jl  —  BDNK CAUSALITY / WELL-POSEDNESS PHASE DIAGRAM.

    The systematic map of where the BDNK frame is simultaneously CAUSAL and
    LINEARLY STABLE (well-posed).  Three independent frame parametrizations are
    mapped, then the M=1.4 production-frame star is checked radius-by-radius.

    GROUNDING (all formulas live in the package; this script only drives them):

    [A] SHUM hatted frame (ŝ,â,q̂), Shum et al. 2509.15303 eqs.67-71, implemented
        in src/transport/Transport.jl :: shum_frame_speeds / shum_frame_wellposed.
          V̂  = (4/3)η̂ + ζ̂
          c0 = cs √(q̂ η̂ /(â V̂))                              (diffusive)
          c± = cs √{ [â(1+ŝ)+q̂ ± √D]/(2â) },  D = q̂²+â²(4q̂+(ŝ-1)²)+2â q̂(1+ŝ)
        Well-posedness + linear stability (eq.71):     0 < q̂ < ŝ.
        Causality (subluminal modified sound):         c₊ ≤ 1  ⇔  cs²·argp ≤ 1.
        Production frame: (ŝ,â,q̂)=(1,1,0.999) ⇒ c₊=√3 cs, c₋≈0.0183 cs.

    [B] KOVTUN frame (Kovtun 1907.08191), the conditions enforced in
        repro/kovtun_sound.jl (at-rest linear stability eq:sound-c1 / eq:sound-c2)
        and the boosted-frame causality (large-k subluminal phase speed).  We map
        the (ε̄1, θ̄) plane of dimensionless first-order coefficients (ε̄1=cs²ε1/γs,
        θ̄=θ/γs, with the production point ε̄1=3, θ̄=4 of the kovtun_sound figure)
        and confirm the production point is stable+causal.

    [C] GENERAL BDNK biquadratic (Keeble & Redondo-Yuste), implemented in
        src/transport/Causality.jl :: causality_flag (Λ₂c⁴-2Λ₁c²+Λ₀=0), as an
        independent causal/real/subluminal monitor of the mapped Shum frame.

    [D] ALONG THE STAR.  For the production frame (ŝ,â,q̂)=(1,1,0.999) we solve the
        M=1.4 TOV star for (i) the Shum Γ=2 polytrope and (ii) a realistic Read
        et al. piecewise polytrope (SLy), and confirm at EVERY radius:
            c₊ = √3 cs(r) ≤ 1   AND   0 < q̂ < ŝ,
        reporting the radius-resolved causality margin 1 − c₊(r).

    Run: cd code/BDNKStar && JULIA_NUM_THREADS=6 julia --project=. \
            repro/bdnk_phase_diagram.jl
=#

include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.EquationOfState: ShumPolytrope, pressure, sound_speed2
using .BDNKStar.TOV: solve_tov
using .BDNKStar.Transport: shum_frame_speeds, shum_frame_wellposed, TransportCoefficients
using .BDNKStar.Causality: causality_flag
using Printf
using LinearAlgebra: eigvals

# Production frame + smallSB-F2 viscosity (Shum §V.A)
const ŝP, âP, q̂P = 1.0, 1.0, 0.999
const η̂P, ζ̂P     = 0.01, 0.01

# ===========================================================================
# Classification primitive (no throw on non-hyperbolic; that IS a class).
# Returns hyperbolic (real speeds), causal (c₊≤1 & c0≤1), well-posed (0<q̂<ŝ).
# ===========================================================================
function classify_shum(ŝ, â, q̂, η̂, ζ̂, cs)
    V̂   = (4/3)*η̂ + ζ̂
    D    = q̂^2 + â^2*(4q̂ + (ŝ-1)^2) + 2*â*q̂*(1+ŝ)
    base = â*(1+ŝ) + q̂
    arg0 = q̂*η̂/(â*V̂)
    sqD  = D ≥ 0 ? sqrt(D) : NaN
    argp = (base + sqD)/(2â)
    argm = (base - sqD)/(2â)
    hyp  = (D ≥ 0) && (argp ≥ 0) && (argm ≥ 0) && (arg0 ≥ 0)
    c0 = cs*sqrt(max(arg0,0.0)); cp = cs*sqrt(max(argp,0.0)); cm = cs*sqrt(max(argm,0.0))
    caus = hyp && (cp ≤ 1.0) && (c0 ≤ 1.0)
    wp   = shum_frame_wellposed(ŝ, â, q̂)
    return (D=D, c0=c0, cp=cp, cm=cm, hyp=hyp, caus=caus, wp=wp, cgood=(cp≤1 && c0≤1))
end

# ===========================================================================
# [A] SHUM frame-space phase diagram
# ===========================================================================
function section_A_shum(cs)
    println("="^78)
    println("[A] SHUM hatted-frame phase diagram  (cs = ", @sprintf("%.5f",cs),
            ",  η̂=ζ̂=", η̂P, ")")
    println("="^78)
    println("    Well-posed/stable  : 0 < q̂ < ŝ          (Shum eq.71)")
    println("    Causal             : c₊=cs√argp ≤ 1 and c0 ≤ 1")
    println()

    # ---- (A1) q̂ vs ŝ plane at â=1 : the well-posed wedge 0<q̂<ŝ ----
    println("  (A1) q̂–ŝ plane at â=1.  Symbol: '#'=WP+causal, '+'=WP only,")
    println("       'c'=causal not WP, '.'=neither.  Production (q̂,ŝ)=(0.999,1).")
    ŝgrid = collect(range(0.05, 2.0, length=21))   # rows (top = large ŝ)
    q̂grid = collect(range(0.05, 2.0, length=21))   # cols
    println("       ", " "^9, "q̂→  (", @sprintf("%.2f",q̂grid[1]), " … ",
            @sprintf("%.2f",q̂grid[end]), ")")
    for ŝ in reverse(ŝgrid)
        row = IOBuffer()
        print(row, @sprintf("    ŝ=%5.2f  ", ŝ))
        for q̂ in q̂grid
            r = classify_shum(ŝ, âP, q̂, η̂P, ζ̂P, cs)
            ch = r.wp ? (r.cgood ? '#' : '+') : (r.cgood ? 'c' : '.')
            print(row, ch)
        end
        println(String(take!(row)))
    end
    # boundary of well-posed: q̂<ŝ (line q̂=ŝ). boundary of causal: c₊=1.
    # find, per ŝ row, the largest q̂ that is BOTH WP and causal:
    println("\n       WP+causal q̂-boundary per ŝ (largest q̂ with 0<q̂<ŝ AND c₊≤1):")
    println("         ŝ      q̂_wp(=ŝ)   q̂_caus(c₊=1)   q̂_max(both)")
    for ŝ in [0.25, 0.5, 1.0, 1.5, 2.0]
        # well-posed limit is q̂ -> ŝ⁻.  causal limit: c₊(q̂)=1.
        q̂_caus = NaN
        qq = range(1e-4, 5.0, length=20000)
        for q̂ in qq
            r = classify_shum(ŝ, âP, q̂, η̂P, ζ̂P, cs)
            if r.cp > 1.0; q̂_caus = q̂; break; end
        end
        q̂both = min(ŝ, isnan(q̂_caus) ? Inf : q̂_caus)
        @printf("        %5.2f   %8.4f    %s   %s\n", ŝ, ŝ,
                isnan(q̂_caus) ? "  (c₊≤1 ∀q̂)" : @sprintf("%10.4f",q̂_caus),
                isfinite(q̂both) ? @sprintf("%8.4f",q̂both) : "   ∞")
    end

    # ---- (A2) â vs q̂ plane at ŝ=1 : WP is â-independent (0<q̂<1) ----
    println("\n  (A2) â–q̂ plane at ŝ=1.  WP (0<q̂<1) is â-independent; causality is")
    println("       NOT (c₊ depends on â).  Symbol legend as (A1).")
    âgrid = collect(range(0.1, 3.0, length=21))
    q̂grid2 = collect(range(0.05, 1.6, length=21))
    println("       ", " "^9, "q̂→  (", @sprintf("%.2f",q̂grid2[1]), " … ",
            @sprintf("%.2f",q̂grid2[end]), ")")
    for â in reverse(âgrid)
        row = IOBuffer()
        print(row, @sprintf("    â=%5.2f  ", â))
        for q̂ in q̂grid2
            r = classify_shum(ŝP, â, q̂, η̂P, ζ̂P, cs)
            ch = r.wp ? (r.cgood ? '#' : '+') : (r.cgood ? 'c' : '.')
            print(row, ch)
        end
        println(String(take!(row)))
    end

    # ---- production-frame margins ----
    rp = classify_shum(ŝP, âP, q̂P, η̂P, ζ̂P, cs)
    println("\n  PRODUCTION (ŝ,â,q̂)=(1,1,0.999):")
    @printf("     c0=%.5f  c₊=%.5f  c₋=%.5f   hyp=%s  causal=%s  well-posed=%s\n",
            rp.c0, rp.cp, rp.cm, rp.hyp, rp.caus, rp.wp)
    @printf("     WP margins:   q̂-0 = %.4f (>0) ;  ŝ-q̂ = %.4f (>0)\n", q̂P, ŝP-q̂P)
    @printf("     causal margin: 1 - c₊ = %.5f   (c₊=√3·cs=%.5f)\n", 1-rp.cp, rp.cp)
    return rp
end

# ===========================================================================
# [B] KOVTUN frame stable+causal region (conditions in repro/kovtun_sound.jl)
# ===========================================================================
# At-rest stability eq:sound-c1 / eq:sound-c2 (dimensionless ε̄1,ε̄2,θ̄,π̄1).
function kovtun_stability(; cs, ε̄1, ε̄2, θ̄, π̄1)
    cs2 = cs^2
    # eq:sound-c1 (e>0):  ε2+π1 > γs/cs² + cs² ε1  ⇒ in barred (÷γs):
    #   ε̄2 + π̄1 > 1/cs² + ε̄1            (since cs²ε1/γs = ε̄1, γs/cs²/γs=1/cs²)
    cond1 = (ε̄2 + π̄1) - (1/cs2 + ε̄1)
    # eq:sound-c2 Routh-Hurwitz (barred form, verbatim kovtun_sound.jl):
    cond2 = ε̄1^2/cs2 + cs2*(ε̄1-ε̄2)*(ε̄1+θ̄)^2*(ε̄1-π̄1) +
            (ε̄1+θ̄)*(2*ε̄1^2 - ε̄1*(ε̄2+π̄1) + (θ̄+ε̄2)*(θ̄+π̄1))
    return cond1, cond2
end

# Large-k (causal) phase speed of the SOUND quartic at rest:  the highest
# |Re ω/k| as k→∞.  From the rest quartic c4 ω⁴+…; the leading large-k balance
# c4 ω⁴ ~ (k² terms) gives the asymptotic group/phase speeds.  We evaluate the
# rest-frame quartic at very large k and take max |Re ω|/k (Kovtun eq:cs0largek).
function kovtun_largek_speed(; cs, ε1, ε2, θ, π1, γs=1.0, w0=1.0)
    cs2 = cs^2
    k = 1e6; k2 = k^2
    c4 = cs2*ε1*θ
    c3 = im*w0*(cs2*ε1 + θ)
    c2 = -( w0^2 + k2*cs2*( cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1 ) )
    c1 = -im*k2*w0*( γs + cs2^2*ε1 + cs2*θ )
    c0 = k2*cs2*( w0^2 + k2*θ*( cs2*(ε2+π1 - cs2*ε1) - γs ) )
    # companion roots of c4 x⁴+c3 x³+c2 x²+c1 x+c0
    a = ComplexF64[c0, c1, c2, c3, c4] ./ c4
    C = zeros(ComplexF64, 4, 4)
    C[2,1]=1; C[3,2]=1; C[4,3]=1
    for i in 1:4; C[i,4] = -a[i]; end
    rts = eigvals(C)
    return maximum(abs(real(z))/k for z in rts)
end

function section_B_kovtun()
    println("\n" * "="^78)
    println("[B] KOVTUN frame stable+causal region (conditions of repro/kovtun_sound.jl)")
    println("="^78)
    # Production point of the kovtun_sound figure: cs=0.5, ε̄1=3, θ̄=4, ε̄2=0, π̄1=3/cs²
    cs   = 0.5
    π̄1P  = 3/cs^2        # = 12
    PROD = (cs=cs, ε̄1=3.0, ε̄2=0.0, θ̄=4.0, π̄1=π̄1P)
    println("    Production point (Kovtun fig picresoundv09): cs=", cs,
            ", ε̄1=", PROD.ε̄1, ", θ̄=", PROD.θ̄, ", ε̄2=", PROD.ε̄2,
            ", π̄1=", @sprintf("%.1f",PROD.π̄1))
    println("    Conditions: c1 = (ε̄2+π̄1)-(1/cs²+ε̄1) > 0  [eq:sound-c1]")
    println("                c2 = Routh-Hurwitz > 0          [eq:sound-c2]")
    println("                causal: large-k max|Re ω/k| ≤ 1 [eq:cs0largek]")

    # ---- (B1) ε̄1 – θ̄ plane at fixed cs, ε̄2=0, π̄1=π̄1P : stable region ----
    println("\n  (B1) ε̄1–θ̄ plane (cs=0.5, ε̄2=0, π̄1=", @sprintf("%.1f",π̄1P),
            ").  '#'=stable+causal, 's'=stable only, '.'=unstable:")
    θ̄grid  = collect(range(0.5, 8.0, length=21))   # rows
    ε̄1grid = collect(range(0.5, 8.0, length=21))   # cols
    println("       ", " "^10, "ε̄1→  (", @sprintf("%.1f",ε̄1grid[1]), " … ",
            @sprintf("%.1f",ε̄1grid[end]), ")")
    for θ̄ in reverse(θ̄grid)
        row = IOBuffer(); print(row, @sprintf("    θ̄=%5.2f  ", θ̄))
        for ε̄1 in ε̄1grid
            c1, c2 = kovtun_stability(; cs=cs, ε̄1=ε̄1, ε̄2=PROD.ε̄2, θ̄=θ̄, π̄1=π̄1P)
            stable = (c1 > 0) && (c2 > 0)
            # causal check needs UNbarred ε1,θ,π1 (×γs, γs=1): ε1=ε̄1/cs²? -> ε̄1=cs²ε1/γs
            ε1 = ε̄1/cs^2; θ = θ̄; π1 = π̄1P
            v = kovtun_largek_speed(; cs=cs, ε1=ε1, ε2=PROD.ε̄2, θ=θ, π1=π1)
            causal = v ≤ 1 + 1e-6
            print(row, stable ? (causal ? '#' : 's') : '.')
        end
        println(String(take!(row)))
    end

    # ---- production point evaluation ----
    c1, c2 = kovtun_stability(; cs=cs, ε̄1=PROD.ε̄1, ε̄2=PROD.ε̄2, θ̄=PROD.θ̄, π̄1=π̄1P)
    ε1 = PROD.ε̄1/cs^2; v = kovtun_largek_speed(; cs=cs, ε1=ε1, ε2=PROD.ε̄2, θ=PROD.θ̄, π1=π̄1P)
    println("\n  PRODUCTION point evaluation:")
    @printf("     eq:sound-c1 = %.4f  (>0 ⇒ stable): %s\n", c1, c1>0)
    @printf("     eq:sound-c2 = %.4f  (>0 ⇒ stable): %s\n", c2, c2>0)
    @printf("     large-k max|Re ω/k| = %.5f  (≤1 ⇒ causal): %s\n", v, v≤1+1e-6)
    inside = (c1>0) && (c2>0) && (v≤1+1e-6)
    @printf("     production inside stable+causal region: %s\n", inside)
    # margins along ε̄1 at fixed θ̄=4: how far to the unstable boundary
    println("\n     stability margin scan along ε̄1 (θ̄=4 fixed): smallest c2:")
    minc2 = Inf; ε̄1_min = NaN
    for ε̄1 in range(0.1, 12.0, length=600)
        _, c2s = kovtun_stability(; cs=cs, ε̄1=ε̄1, ε̄2=PROD.ε̄2, θ̄=PROD.θ̄, π̄1=π̄1P)
        if c2s < minc2; minc2 = c2s; ε̄1_min = ε̄1; end
    end
    @printf("       min eq:sound-c2 over ε̄1∈[0.1,12] = %.4f at ε̄1=%.3f  (≥0 ⇒ whole line stable)\n",
            minc2, ε̄1_min)
    return (inside=inside, c1=c1, c2=c2, v=v)
end

# ===========================================================================
# [C] General BDNK biquadratic monitor of the mapped production frame
# ===========================================================================
function section_C_biquad(eos, εc)
    println("\n" * "="^78)
    println("[C] General BDNK biquadratic monitor (Causality.jl) of the production frame")
    println("="^78)
    # Map the production Shum frame to dimensionful (η,ζ,τε,τP,τQ) at the centre
    p   = pressure(eos, εc); cs2 = sound_speed2(eos, εc); cs = sqrt(cs2); ρ = εc + p
    V̂   = (4/3)*η̂P + ζ̂P; L = 1.0
    η   = q̂P*L*cs2*ρ*η̂P; ζ = q̂P*L*cs2*ρ*ζ̂P
    τε  = V̂*L; τp = ŝP*cs2*L*V̂; τQ = âP*L*V̂
    tc  = TransportCoefficients(η=η, ζ=ζ, κQ=0.0, τε=τε, τP=τp, τQ=τQ, L=L)
    f   = causality_flag(p, εc, cs2, tc)
    @printf("    centre cs²=%.5f  η=%.3e ζ=%.3e τε=%.3e τP=%.3e τQ=%.3e\n",
            cs2, η, ζ, τε, τp, τQ)
    @printf("    biquadratic: c²₋=%.4f c²₊=%.4f disc=%.3e\n", f.c2_minus, f.c2_plus, f.disc)
    @printf("    real_speeds=%s  nonneg=%s  subluminal(c²₊≤1)=%s  causal=%s\n",
            f.real_speeds, f.nonneg, f.subluminal, f.causal)
    println("    (NB: the biquadratic uses a DIFFERENT non-dimensionalisation than the")
    println("     Shum hatted frame — see shum_frame_analysis.jl §4 — so its c²₊ need")
    println("     not equal 3cs²; it is an independent real/hyperbolic monitor.)")
    return f
end

# ===========================================================================
# [D] ALONG THE STAR — radius-resolved causality margin, production frame
# ===========================================================================
function along_star(label, eos, εc; h, rmax, in_msun=false)
    star = solve_tov(eos, εc; h=h, ptol_rel=1e-12, rmax=rmax)
    Mrep = in_msun ? mass_solar(star.M) : star.M   # ShumPolytrope: M⊙=G=c=1 already
    N = length(star.r)
    allreal = true
    cp_max = 0.0; r_cpmax = 0.0; cs_cpmax = 0.0
    margin_min = Inf; r_margin = 0.0
    cs_max = 0.0; n_superlum = 0
    for i in 1:N
        ε = star.ε[i]; ε <= 0 && continue
        cs = sqrt(max(sound_speed2(eos, ε), 0.0))
        cs_max = max(cs_max, cs)
        c0, cp, cm = shum_frame_speeds(ŝP, âP, q̂P, η̂P, ζ̂P, cs)
        allreal &= isfinite(c0) && isfinite(cp) && isfinite(cm) && cp≥0 && cm≥0 && c0≥0
        if cp > cp_max; cp_max = cp; r_cpmax = star.r[i]; cs_cpmax = cs; end
        m = 1 - cp
        if m < margin_min; margin_min = m; r_margin = star.r[i]; end
        cp > 1 && (n_superlum += 1)
    end
    wp = shum_frame_wellposed(ŝP, âP, q̂P)
    println("\n  ", label, ":  M=", @sprintf("%.4f",Mrep), " M⊙",
            "  R=", @sprintf("%.4f",star.R), "  (", N, " radii)")
    @printf("     max cs = %.5f (1/√3=%.5f ⇒ c₊≤1 needs cs≤1/√3)\n", cs_max, 1/sqrt(3))
    @printf("     max c₊ = √3·cs = %.5f at r=%.4f (cs=%.5f)\n", cp_max, r_cpmax, cs_cpmax)
    @printf("     MIN causality margin (1-c₊) = %.5f at r=%.4f\n", margin_min, r_margin)
    @printf("     speeds real ∀r: %s   well-posed (0<q̂<ŝ) ∀r: %s   superluminal radii: %d\n",
            allreal, wp, n_superlum)
    # sampled profile
    println("       r/R      cs        c₊       1-c₊(margin)")
    idxs = unique(clamp.(round.(Int, range(1, N, length=7)), 1, N))
    for i in idxs
        ε = star.ε[i]; ε <= 0 && continue
        cs = sqrt(max(sound_speed2(eos, ε),0.0))
        _, cp, _ = shum_frame_speeds(ŝP, âP, q̂P, η̂P, ζ̂P, cs)
        @printf("      %5.3f   %7.5f  %7.5f   %8.5f\n", star.r[i]/star.R, cs, cp, 1-cp)
    end
    return (M=Mrep, R=star.R, cp_max=cp_max, margin_min=margin_min,
            allreal=allreal, wp=wp, n_superlum=n_superlum, cs_max=cs_max)
end

function section_D_star()
    println("\n" * "="^78)
    println("[D] ALONG THE STAR — production frame causality margin (M=1.4)")
    println("="^78)
    # (i) Shum Γ=2 polytrope, ρ0c=0.00128 ⇒ M_T=1.4
    κ = 100.0; eos1 = ShumPolytrope(κ); ρ0c = 0.00128; εc1 = ρ0c + κ*ρ0c^2
    s1 = along_star("ShumPolytrope (κ=100, M=1.4)", eos1, εc1; h=2e-4, rmax=50.0)
    # (ii) Realistic SLy piecewise polytrope, central ε tuned to M≈1.4 M⊙
    eos2 = piecewise_polytrope(:SLy)
    # SLy: M=1.4 near ε_c ≈ 8.9e-4 km^-2 (ρ_c ~ 9e14 g/cm³). Bracket-search M=1.4.
    function Mof(εc)
        st = solve_tov(eos2, εc; h=5e-3, ptol_rel=1e-10, rmax=40.0)
        return mass_solar(st.M)
    end
    lo, hi = 4e-4, 2.0e-3
    Mlo, Mhi = Mof(lo), Mof(hi)
    εc2 = lo
    for _ in 1:60
        mid = sqrt(lo*hi); Mm = Mof(mid)
        if Mm < 1.4; lo = mid; else hi = mid; end
        εc2 = mid
        abs(Mm-1.4) < 1e-3 && break
    end
    @printf("\n  SLy central-ε bracket for M=1.4 M⊙: εc=%.4e km⁻² (M=%.4f M⊙)\n",
            εc2, mass_solar(solve_tov(eos2, εc2; h=5e-3, rmax=40.0).M))
    s2 = along_star("PiecewisePolytrope SLy (M≈1.4 M⊙, realistic)", eos2, εc2;
                    h=5e-3, rmax=40.0, in_msun=true)
    return (shum=s1, sly=s2)
end

# ===========================================================================
# DRIVER
# ===========================================================================
function run()
    println("#"^78)
    println("# BDNK CAUSALITY / WELL-POSEDNESS PHASE DIAGRAM")
    println("#"^78)
    κ = 100.0; eos = ShumPolytrope(κ); ρ0c = 0.00128; εc = ρ0c + κ*ρ0c^2
    cs_centre = sqrt(sound_speed2(eos, εc))

    A = section_A_shum(cs_centre)
    B = section_B_kovtun()
    C = section_C_biquad(eos, εc)
    D = section_D_star()

    println("\n" * "="^78)
    println("SUMMARY — region boundaries, production-frame status, margins")
    println("="^78)
    println("  [A] SHUM frame:")
    println("      well-posed/stable region : 0 < q̂ < ŝ   (â-independent)")
    @printf( "      causal boundary at â=1   : c₊=cs√((â(1+ŝ)+q̂+√D)/2â)=1\n")
    @printf( "      production (1,1,0.999)   : WP=%s causal=%s  margins q̂-0=%.3f ŝ-q̂=%.3f 1-c₊=%.5f\n",
            A.wp, A.caus, q̂P, ŝP-q̂P, 1-A.cp)
    println("  [B] KOVTUN frame:")
    @printf( "      production point stable+causal : %s  (c1=%.3f c2=%.3f largek-speed=%.4f)\n",
            B.inside, B.c1, B.c2, B.v)
    println("  [C] BDNK biquadratic monitor (production frame, ShumPolytrope centre):")
    @printf( "      real=%s nonneg=%s subluminal=%s causal=%s\n",
            C.real_speeds, C.nonneg, C.subluminal, C.causal)
    println("  [D] ALONG THE STAR (production frame, c₊=√3 cs):")
    @printf( "      ShumPolytrope M=%.3f : real∀r=%s WP∀r=%s superlum=%d  min margin(1-c₊)=%.5f\n",
            D.shum.M, D.shum.allreal, D.shum.wp, D.shum.n_superlum, D.shum.margin_min)
    @printf( "      SLy realistic  M=%.3f : real∀r=%s WP∀r=%s superlum=%d  min margin(1-c₊)=%.5f\n",
            D.sly.M, D.sly.allreal, D.sly.wp, D.sly.n_superlum, D.sly.margin_min)
    safe = A.wp && A.caus && B.inside && D.shum.n_superlum==0 && D.sly.n_superlum==0 &&
           D.shum.wp && D.sly.wp
    @printf("\n  PRODUCTION FRAME + STAR SAFELY INSIDE CAUSAL+STABLE REGION : %s\n", safe)
    println("="^78)
    return (A=A, B=B, C=C, D=D, safe=safe)
end

run()
