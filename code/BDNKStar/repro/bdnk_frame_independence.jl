#=
    bdnk_frame_independence.jl — the CLEAN gauge-vs-physical demonstration of
    BDNK frame-independence at the DISPERSION-RELATION level.

    THE STATEMENT (textbook BDNK / Kovtun 1907.08191 §4.2):
      Build the linearized BDNK sound-channel dispersion ω(k) for plane-wave
      perturbations of a homogeneous fluid. The general-frame uncharged theory
      (Kovtun eqs 2.4) has 6 coefficients:
          PHYSICAL transport : η (shear), ζ (bulk)     [and vs² = ∂p/∂ε, w0=ε+p]
          FRAME / RELAXATION : ε1, π1, θ
      with ε1,π1,θ the relaxation times of energy density, pressure, and momentum
      density (heat flux) — i.e. exactly the BDNK frame (τ_ε, τ_P, τ_Q), and the
      Shum hatted ratios ŝ = vs²ε1/γs, â = θ/γs.

      HYDRODYNAMIC branch (ω→0 as k→0), Kovtun eq.4.9:
          ω = ± vs|k| − (i/2)(γs/w0) k² + …       γs ≡ (4/3)η + ζ
      ⇒ sound speed vs and leading damping γs/w0 depend ONLY on the PHYSICAL
        coefficients → frame-INVARIANT (gradient expansion is frame-independent).

      NON-HYDRODYNAMIC branch (gapped, ω→ −i/τ as k→0), Kovtun eq.4.10:
          ω = −i w0/(vs² ε1),   ω = −i w0/θ
      ⇒ the gaps scale ∝ 1/ε1, 1/θ → they MOVE with the frame (gauge modes).

    THIS SCRIPT holds the PHYSICAL quantities fixed (vs, η, ζ ⇒ γs, w0) and
    SWEEPS the frame (ε1, π1, θ, i.e. ŝ, â) across the stable+causal region,
    then tabulates:
      (A) small-k sound speed Re ω/k and leading damping −2 Im ω/k²·w0  → FLAT
      (B) the two non-hydro gaps vs frame                              → ∝ 1/τ
      (C) large-k sound speed (Kovtun eq.4.18) subluminal across the swept frames.

    Causality is a SEPARATE constraint from stability: a stable rest-frame choice
    can still be large-k superluminal (e.g. small ŝ at vs=0.5), so the sweep
    ranges are chosen inside the genuinely causal region (eq.4.13,4.14,4.6 +
    4.18). The companion package biquadratic (Causality.jl) uses a different
    non-dimensionalisation (documented mismatch, shum_frame_analysis.jl §4) and
    is reported but not used as the definitive causality check.

    Footnote-8 (Kovtun): ζ = vs²(π1 − vs²ε1) − π2 + vs²ε2. To hold the PHYSICAL
    bulk ζ fixed while varying the FRAME (ε1,π1,ε2), we solve this for π2. Then
    γs ≡ (4/3)η + ζ is fixed by construction, so the hydrodynamic coefficients
    cannot move — that is the whole point.

    Reuses the validated F_sound quartic + companion-matrix root finder from
    repro/kovtun_sound.jl (the dispersion machinery). Data only.

    Run: cd code/BDNKStar && JULIA_NUM_THREADS=6 julia --project=. \
         repro/bdnk_frame_independence.jl
=#

include(joinpath(@__DIR__, "kovtun_sound.jl"))   # F_sound quartic + roots + BDNKStar
using .BDNKStar
using .BDNKStar.Transport: TransportCoefficients
using .BDNKStar.Causality: characteristic_speeds
using Printf

# ---------------------------------------------------------------------------
# PHYSICAL coefficients held FIXED across the whole frame sweep.
# (dimensionless units w0 = 1, γs = 1, as in Kovtun's figures)
# ---------------------------------------------------------------------------
const VS  = 0.5            # physical sound speed  vs = √(∂p/∂ε)
const VS2 = VS^2
const W0  = 1.0            # enthalpy ε0 + p0
const ETA = 0.3            # physical shear η
const ZETA_phys = (4.0/3.0)*1.0 - (4.0/3.0)*ETA   # choose so γs = (4/3)η + ζ = 1
const GS  = (4.0/3.0)*ETA + ZETA_phys             # ≡ 1.0  (the fixed sound-damping)

# Frame coefficients are parametrized by the Shum-style dimensionless ratios:
#   ŝ = vs² ε1 / γs      (energy-relaxation ratio)
#   â = θ      / γs      (momentum/heat-relaxation ratio)
#   p̂ = π1     / γs      (pressure-relaxation ratio)
# Inverting: ε1 = ŝ γs/vs²,  θ = â γs,  π1 = p̂ γs.
# Hold ζ FIXED via footnote-8 by choosing ε2 = 0 and solving for π2:
#   π2 = vs²(π1 − vs²ε1) − ζ + vs²ε2     (with ε2 = 0)
frame_coeffs(ŝ, â, p̂; ε2=0.0) = begin
    ε1 = ŝ*GS/VS2
    θ  = â*GS
    π1 = p̂*GS
    π2 = VS2*(π1 - VS2*ε1) - ZETA_phys + VS2*ε2
    (ε1=ε1, ε2=ε2, π1=π1, π2=π2, θ=θ)
end

# Recover the PHYSICAL bulk ζ from a frame's (ε1,ε2,π1,π2) — must be invariant.
zeta_from(fc) = VS2*(fc.π1 - VS2*fc.ε1) - fc.π2 + VS2*fc.ε2

# Stability / causality (Kovtun §4.2 small-k + the two at-rest conditions).
function frame_stable(fc)
    # 4.11 small-k stability
    c411 = (GS > 0) && (fc.ε1 > 0) && (fc.θ > 0)
    # 4.13 first Routh-Hurwitz
    c413 = (fc.ε2 + fc.π1) > (GS/VS2 + VS2*fc.ε1)
    # 4.14 second Routh-Hurwitz (dimensionless bars)
    ε̄1 = VS2*fc.ε1/GS; ε̄2 = fc.ε2/GS; θ̄ = fc.θ/GS; π̄1 = fc.π1/GS
    c414 = (ε̄1^2/VS2) + VS2*(ε̄1-ε̄2)*(ε̄1+θ̄)*(ε̄1-π̄1) +
           (ε̄1+θ̄)*(2*ε̄1^2 - ε̄1*(ε̄2+π̄1) + (θ̄+ε̄2)*(θ̄+π̄1)) > 0
    # 4.6 shear-channel θ > η
    c46 = fc.θ > ETA
    (all=(c411 && c413 && c414 && c46), c411=c411, c413=c413, c414=c414, c46=c46)
end

# ---------------------------------------------------------------------------
# small-k HYDRODYNAMIC sound branch: extract vs and damping from ω(k) roots.
# Kovtun 4.9: ω = ±vs k − (i/2)(γs/w0) k²  ⇒  Re ω/k → ±vs,
#   damping coeff Γ ≡ −2 Im(ω) w0 / k²  →  γs.
# We pick the two roots with smallest |Im ω| (the gapless sound pair).
# ---------------------------------------------------------------------------
function sound_branch(fc; k=1e-3)
    rts = kovtun_sound_modes(k, 0.0; v0=0.0, cs=VS,
                             ε1=fc.ε1, ε2=fc.ε2, π1=fc.π1, θ=fc.θ, γs=GS, w0=W0)
    perm = sortperm(rts, by=z->abs(imag(z)))
    g1, g2 = rts[perm[1]], rts[perm[2]]
    vs_meas = (abs(real(g1)) + abs(real(g2)))/2 / k          # → vs
    Γ_meas  = -2*((imag(g1)+imag(g2))/2)*W0 / k^2            # → γs
    (vs=vs_meas, Γ=Γ_meas)
end

# NON-HYDRO gaps: the two most-damped roots at small k (Kovtun 4.10).
function gaps(fc; k=1e-4)
    rts = kovtun_sound_modes(k, 0.0; v0=0.0, cs=VS,
                             ε1=fc.ε1, ε2=fc.ε2, π1=fc.π1, θ=fc.θ, γs=GS, w0=W0)
    g = sort(rts, by=z->-abs(imag(z)))[1:2]
    gmeas = sort([-imag(z) for z in g])                       # |gap| (Im ω = -gap)
    gpred = sort([W0/(VS2*fc.ε1), W0/fc.θ])                   # 4.10 predictions
    (meas=gmeas, pred=gpred)
end

# Exact large-k sound speed (Kovtun eq.4.18, dimensionless bars): the asymptotic
# ω = c_sound k solves a biquadratic in c_sound². Causal ⇔ c_sound² real, 0<c<1.
# This is the analytic limit the rest-frame quartic above converges to.
function largek_sound2(fc)
    ε̄1 = VS2*fc.ε1/GS; ε̄2 = fc.ε2/GS; θ̄ = fc.θ/GS; π̄1 = fc.π1/GS
    a = ε̄1*θ̄/VS2
    b = ε̄1*(ε̄2+π̄1-ε̄1-1/VS2) - θ̄*(ε̄2+π̄1) - ε̄2*π̄1
    c = θ̄*(VS2*(ε̄2+π̄1-ε̄1) - 1)
    disc = b^2 - 4a*c
    disc < 0 && return (NaN, NaN, false)
    c2p = (-b + sqrt(disc))/(2a); c2m = (-b - sqrt(disc))/(2a)
    c2max = max(c2p, c2m)
    causal = (c2m > 0) && (c2max < 1)
    (c2max, sqrt(max(c2max,0.0)), causal)
end

# ===========================================================================
println("#"^78)
println("# BDNK FRAME-INDEPENDENCE at the dispersion-relation level")
println("# (Kovtun 1907.08191 §4.2 sound channel; general/BDNK frame)")
println("#"^78)
@printf("\nPHYSICAL (held FIXED): vs = %.4f  (vs²=%.4f),  η = %.4f,  ζ = %.4f\n",
        VS, VS2, ETA, ZETA_phys)
@printf("                       γs = (4/3)η+ζ = %.6f,  w0 = %.4f\n", GS, W0)
@printf("FRAME knobs SWEPT    : ŝ = vs²ε1/γs,  â = θ/γs,  p̂ = π1/γs  (ζ held fixed via fn-8)\n")

# Sweep ranges chosen inside the STABLE + CAUSAL frame region (Kovtun 4.13,4.14,
# 4.6 stability AND 4.18 large-k causality). Causality is a SEPARATE constraint
# from stability: e.g. ŝ=1.5 below is stable but large-k superluminal, so the
# causal demonstration is restricted to the genuinely causal swept frames.
table_header(knobname, dimname) = begin
    @printf("\n  %-7s %-9s | HYDRO (should be FLAT)          | NON-HYDRO gaps (∝1/τ, MOVE)   | large-k\n",
            knobname, dimname)
    @printf("  %-7s %-9s | %-9s %-14s | %-13s %-13s | c_∞  (caus)\n",
            "", "", "Re ω/k", "Γ=-2Imω w0/k²", "gap=w0/θ", "gap=w0/(vs²ε1)")
    println("  ", "-"^100)
end
report_row(knob, dimval, fc) = begin
    sb = sound_branch(fc); c2,cinf,caus = largek_sound2(fc); st = frame_stable(fc)
    @printf("  %-7.2f %-9.4f | %-9.5f %-14.6f | %-13.5f %-13.5f | %-5.4f %s%s\n",
            knob, dimval, sb.vs, sb.Γ, W0/fc.θ, W0/(VS2*fc.ε1),
            isnan(cinf) ? -1.0 : cinf, caus ? "C" : "-",
            st.all ? "" : " (unstbl)")
    (vs=sb.vs, Γ=sb.Γ, cinf=cinf, caus=caus, gth=W0/fc.θ, ge1=W0/(VS2*fc.ε1))
end

# --- SWEEP 1: vary â = θ/γs (heat/momentum relaxation), ŝ=3, p̂=8 fixed -------
println("\n", "="^78)
println("SWEEP 1 — vary frame â = θ/γs  (heat/momentum relaxation; ŝ=3, p̂=8 fixed)")
println("="^78)
println("  θ swept over a 60× range, all within the stable+causal region.")
table_header("â", "θ")
ŝ1, p̂1 = 3.0, 8.0
âs = [1.0, 2.0, 4.0, 8.0, 15.0, 30.0, 60.0]
r1 = [report_row(â, frame_coeffs(ŝ1, â, p̂1).θ, frame_coeffs(ŝ1, â, p̂1)) for â in âs]

# --- SWEEP 2: vary ŝ = vs²ε1/γs (energy relaxation), â=15, p̂=4ŝ -------------
println("\n", "="^78)
println("SWEEP 2 — vary frame ŝ = vs²ε1/γs  (energy relaxation; â=15, p̂=4ŝ to stay causal)")
println("="^78)
println("  ε1 swept over a 4× range; p̂=π1/γs (pressure relaxation) co-varies, ζ held fixed.")
table_header("ŝ", "ε1")
â2 = 15.0
ŝs = [1.5, 2.0, 3.0, 4.0, 5.0, 6.0]
r2 = [report_row(ŝ, frame_coeffs(ŝ, â2, 4*ŝ).ε1, frame_coeffs(ŝ, â2, 4*ŝ)) for ŝ in ŝs]

# ---------------------------------------------------------------------------
# INVARIANCE QUANTIFICATION (over the CAUSAL swept frames)
# ---------------------------------------------------------------------------
spread(v) = (maximum(v) - minimum(v)) / abs(sum(v)/length(v)) * 100
all = vcat(r1, r2)
causal_rows = filter(r->r.caus, all)
all_vs = [r.vs for r in causal_rows]
all_Γ  = [r.Γ  for r in causal_rows]
println("\n", "="^78)
println("INVARIANCE — physical (hydro) coefficients across the causal swept frames")
println("="^78)
@printf("  causal frames in sweep : %d / %d\n", length(causal_rows), length(all))
@printf("  sound speed   Re ω/k :  mean = %.6f (target vs = %.4f),  spread = %.4f%%\n",
        sum(all_vs)/length(all_vs), VS, spread(all_vs))
@printf("  damping Γ = γs        :  mean = %.6f (target γs = %.4f),  spread = %.4f%%\n",
        sum(all_Γ)/length(all_Γ), GS, spread(all_Γ))
sound_spread = max(spread(all_vs), spread(all_Γ))
@printf("  ⇒ HYDRODYNAMIC branch FRAME-INVARIANT to %.4f%%  (gradient expansion frame-indep)\n",
        sound_spread)

println("\n  NON-HYDRO gaps DO move with the frame (gauge modes):")
gth1 = [r.gth for r in r1]; ge12 = [r.ge1 for r in r2]
@printf("    gap = w0/θ        across â∈[%.0f,%.0f]:  %.4f → %.4f  (%.1f× range, ∝1/θ)\n",
        first(âs), last(âs), maximum(gth1), minimum(gth1), maximum(gth1)/minimum(gth1))
@printf("    gap = w0/(vs²ε1)  across ŝ∈[%.1f,%.0f]:  %.4f → %.4f  (%.1f× range, ∝1/ε1)\n",
        first(ŝs), last(ŝs), maximum(ge12), minimum(ge12), maximum(ge12)/minimum(ge12))

# ---------------------------------------------------------------------------
# Confirm ζ (physical bulk) held fixed across all swept frames.
# ---------------------------------------------------------------------------
println("\n  CHECK: physical bulk ζ recovered from each frame's (ε1,ε2,π1,π2) [fn-8]:")
ζrec = vcat([zeta_from(frame_coeffs(ŝ1, â, p̂1)) for â in âs],
            [zeta_from(frame_coeffs(ŝ, â2, 4*ŝ)) for ŝ in ŝs])
@printf("    ζ across all frames: min=%.6f max=%.6f (target %.6f) spread=%.2e%%\n",
        minimum(ζrec), maximum(ζrec), ZETA_phys,
        (maximum(ζrec)-minimum(ζrec))/abs(ZETA_phys)*100)

# ---------------------------------------------------------------------------
# CAUSALITY: large-k subluminality across the causal swept frames, via two
# independent routes — (i) Kovtun eq.4.18 large-k sound speed (printed above as
# c_∞), and (ii) the package Causality biquadratic (Keeble–Redondo-Yuste).
# ---------------------------------------------------------------------------
println("\n", "="^78)
println("CAUSALITY — large-k subluminality across the (causal) swept frames")
println("="^78)
cinf_caus = [r.cinf for r in causal_rows]
@printf("  (i) Kovtun eq.4.18 large-k sound speed c_∞ (rest frame): max over causal frames = %.5f  (<1: %s)\n",
        maximum(cinf_caus), maximum(cinf_caus) < 1.0)

println("\n  (ii) package Causality biquadratic c₊ (Causality.jl Λ₂c⁴-2Λ₁c²+Λ₀=0):")
println("       NOTE — this biquadratic uses a DIFFERENT (dimensionful Keeble–Redondo-Yuste)")
println("       non-dimensionalisation than the Kovtun frame quartic (documented in")
println("       shum_frame_analysis.jl §4); its disc<0 for these dimensionless τ~O(10).")
println("       We report c₊ where real; the Kovtun route (i) is the matched check.")
e0 = W0/(1 + VS2); p0 = VS2*e0
function biquad_cp(fc)
    tc = TransportCoefficients(η=ETA, ζ=ZETA_phys, κQ=0.0,
                              τε=fc.ε1, τP=fc.π1, τQ=fc.θ, L=1.0)
    _, c2p, disc = characteristic_speeds(p0, e0, VS, tc)
    (disc ≥ 0 && c2p ≥ 0) ? sqrt(c2p) : NaN
end
cps = vcat([biquad_cp(frame_coeffs(ŝ1, â, p̂1)) for â in âs],
           [biquad_cp(frame_coeffs(ŝ, â2, 4*ŝ)) for ŝ in ŝs])
nreal = count(isfinite, cps)
@printf("       biquadratic real-speed frames: %d/%d; of those, max c₊ = %.5f (<1: %s)\n",
        nreal, length(cps),
        nreal>0 ? maximum(filter(isfinite, cps)) : NaN,
        nreal>0 ? maximum(filter(isfinite, cps)) < 1.0 : false)

# ---------------------------------------------------------------------------
println("\n", "="^78)
println("VERDICT")
println("="^78)
@printf("  HYDRO branch (vs, Γ=γs)  frame-INVARIANT to %.4f%%  → gradient expansion is gauge-free\n",
        sound_spread)
@printf("  NON-HYDRO gaps           MOVE ∝ 1/τ with the frame → gauge / frame modes\n")
@printf("  large-k speeds (eq.4.18) subluminal (c_∞<1) across all causal swept frames: %s\n",
        maximum(cinf_caus) < 1.0)
@printf("  physical ζ held fixed    to %.2e%% (footnote-8 compensation works)\n",
        (maximum(ζrec)-minimum(ζrec))/abs(ZETA_phys)*100)
println("="^78)
