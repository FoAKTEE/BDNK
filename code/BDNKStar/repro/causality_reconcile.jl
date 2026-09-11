#=
    causality_reconcile.jl — RECONCILE the three characteristic-speed code paths.

    Paths:
      (K)  kovtun_sound.jl  : full F_sound quartic ω(k); LARGE-k slopes ω/k → c∞
           are the characteristic speeds (the object Causality.jl should equal).
      (K418) bdnk_frame_independence.jl largek_sound2: the ANALYTIC large-k
           biquadratic of that quartic (Kovtun eq.4.18, dimensionless bars).
      (C)  Causality.jl     : characteristic_speeds biquadratic Λ₂c⁴-2Λ₁c²+Λ₀=0
           (Keeble & Redondo-Yuste port), fed dimensionful (η,ζ,τε,τP,τQ).
      (S)  Transport.jl     : shum_frame_speeds (scale-free, c₊²+c₋²=3cs²).

    This script does the apples-to-apples large-k comparison the loose-end note
    asked for: for several frames, compute c∞ from (K) numerically, from (K418)
    analytically, and from (C) — and see where (C) diverges and WHY.  Data only.
=#

include(joinpath(@__DIR__, "kovtun_sound.jl"))   # F_sound quartic + roots + BDNKStar
using .BDNKStar
using .BDNKStar.Transport: TransportCoefficients
using .BDNKStar.Causality: characteristic_speeds
import .BDNKStar.Causality   # Λ0,Λ1,Λ2 are unexported; reach via the module
const Λ0 = Causality.Λ0
const Λ1 = Causality.Λ1
const Λ2 = Causality.Λ2
using Printf

# ---------------------------------------------------------------------------
# PHYSICAL block (dimensionless Kovtun units), same as bdnk_frame_independence.
# ---------------------------------------------------------------------------
const VS  = 0.5; const VS2 = VS^2; const W0 = 1.0
const ETA = 0.3
const ZETA = (4.0/3.0)*1.0 - (4.0/3.0)*ETA      # γs = (4/3)η+ζ = 1
const GS  = (4.0/3.0)*ETA + ZETA                # = 1

# Kovtun frame coeffs from Shum-style ratios (ŝ=vs²ε1/γs, â=θ/γs, p̂=π1/γs),
# ζ held fixed via footnote-8 (ε2=0).
frame_coeffs(ŝ, â, p̂; ε2=0.0) = begin
    ε1 = ŝ*GS/VS2; θ = â*GS; π1 = p̂*GS
    π2 = VS2*(π1 - VS2*ε1) - ZETA + VS2*ε2
    (ε1=ε1, ε2=ε2, π1=π1, π2=π2, θ=θ)
end

# (K) numeric large-k slopes from the FULL quartic: max |Re ω/k| as k→∞.
function largek_numeric(fc; k=1e7)
    rts = kovtun_sound_modes(k, 0.0; v0=0.0, cs=VS,
                             ε1=fc.ε1, ε2=fc.ε2, π1=fc.π1, θ=fc.θ, γs=GS, w0=W0)
    speeds = sort([abs(real(z))/k for z in rts])   # 4 slopes; sound pair = the 2 nonzero-ish
    (cmax = maximum(speeds), all = speeds)
end

# (K418) analytic large-k biquadratic of the quartic (Kovtun eq.4.18 bars).
function largek_analytic(fc)
    ε̄1 = VS2*fc.ε1/GS; ε̄2 = fc.ε2/GS; θ̄ = fc.θ/GS; π̄1 = fc.π1/GS
    a = ε̄1*θ̄/VS2
    b = ε̄1*(ε̄2+π̄1-ε̄1-1/VS2) - θ̄*(ε̄2+π̄1) - ε̄2*π̄1
    c = θ̄*(VS2*(ε̄2+π̄1-ε̄1) - 1)
    disc = b^2 - 4a*c
    disc < 0 && return (c2p=NaN, c2m=NaN, sum=NaN, ok=false)
    c2p = (-b + sqrt(disc))/(2a); c2m = (-b - sqrt(disc))/(2a)
    (c2p=max(c2p,c2m), c2m=min(c2p,c2m), sum=c2p+c2m, ok=true)
end

# (C) Causality.jl biquadratic, fed (η,ζ,τε,τP,τQ) = (ETA,ZETA, ε1,π1,θ).
# NOTE: This is the literal mapping used by bdnk_frame_independence.jl §(ii) and
# shum_frame_analysis §4: τε←ε1, τP←π1, τQ←θ.
function causality_biquad(fc)
    e0 = W0/(1 + VS2); p0 = VS2*e0
    tc = TransportCoefficients(η=ETA, ζ=ZETA, κQ=0.0,
                              τε=fc.ε1, τP=fc.π1, τQ=fc.θ, L=1.0)
    c2m, c2p, disc = characteristic_speeds(p0, VS, VS, tc)
    λ0 = Λ0(p0, e0, VS, ETA, ZETA, fc.ε1, fc.π1, fc.θ, 1.0)
    λ1 = Λ1(p0, e0, VS, ETA, ZETA, fc.ε1, fc.π1, fc.θ, 1.0)
    λ2 = Λ2(p0, e0, VS, ETA, ZETA, fc.ε1, fc.π1, fc.θ, 1.0)
    (c2m=c2m, c2p=c2p, disc=disc, sum=2λ1/λ2, prod=λ0/λ2, λ0=λ0, λ1=λ1, λ2=λ2)
end

# ---------------------------------------------------------------------------
# The KEY structural identities of the Causality.jl biquadratic.
# Λ2 c⁴ - 2Λ1 c² + Λ0 = 0  ⇒  c₊²+c₋² = 2Λ1/Λ2,  c₊²c₋² = Λ0/Λ2.
# With the reference coefficients (after the common prefactor cancels):
#   2Λ1/Λ2 = cs² (τε + τQ(τP+τε)) / (τQ τε)
#          = cs² ( 1/τQ + τP/τε + 1 )
#   Λ0/Λ2  = cs⁴ · 2 L² (3ζ+4η)² (τP-1) τQ / 9 ... let's just print it numerically.
# The "(-1+τP)" = (τP-1) factor and the bare "1" and "1/τQ" reveal that τε,τP,τQ
# here are DIMENSIONLESS ratios already (subtracting 1 from τP only makes sense if
# τP is dimensionless), NOT the dimensionful relaxation times.  We test that.
# ---------------------------------------------------------------------------

println("#"^80)
println("# RECONCILE characteristic-speed paths: Kovtun large-k  vs  Causality.jl")
println("#"^80)
@printf("\nPHYSICAL: vs=%.3f (vs²=%.3f) η=%.3f ζ=%.4f γs=%.4f w0=%.3f\n",
        VS, VS2, ETA, ZETA, GS, W0)

# Frames: production-like + several Kovtun frames (the 13-frame style sweep).
frames = [
    ("prod-like ŝ=â=1, p̂=8",   frame_coeffs(1.0, 1.0, 8.0)),
    ("ŝ=3 â=4 p̂=8",            frame_coeffs(3.0, 4.0, 8.0)),
    ("ŝ=3 â=15 p̂=8",           frame_coeffs(3.0, 15.0, 8.0)),
    ("ŝ=3 â=30 p̂=8",           frame_coeffs(3.0, 30.0, 8.0)),
    ("ŝ=4 â=15 p̂=16",          frame_coeffs(4.0, 15.0, 16.0)),
    ("ŝ=6 â=15 p̂=24",          frame_coeffs(6.0, 15.0, 24.0)),
]

println("\n", "="^80)
println("TABLE 1 — large-k sound speed²  c∞²  (the characteristic speed²)")
println("="^80)
@printf("  %-22s | %-10s | %-10s | %-22s\n",
        "frame", "K(numeric)", "K418(anal)", "Causality.jl biquad")
@printf("  %-22s | %-10s | %-10s | %-22s\n",
        "", "max c∞²", "c2p", "c2p  (disc sign)")
println("  ", "-"^78)
for (name, fc) in frames
    kn = largek_numeric(fc)
    ka = largek_analytic(fc)
    cb = causality_biquad(fc)
    @printf("  %-22s | %-10.5f | %-10.5f | %-8.4g  (disc %s)\n",
            name, kn.cmax^2, ka.c2p,
            isnan(cb.c2p) ? -1.0 : cb.c2p, cb.disc ≥ 0 ? "≥0" : "<0")
end

println("\n", "="^80)
println("TABLE 2 — WHY Causality.jl disagrees: the sum/prod of roots c₊²+c₋²")
println("="^80)
println("  Kovtun large-k (K418) sum c₊²+c₋² vs Causality.jl 2Λ1/Λ2.")
println("  Causality.jl 2Λ1/Λ2 = cs²(1/τQ + τP/τε + 1) with τ←(ε1,π1,θ).")
@printf("  %-22s | %-12s | %-26s | %-10s\n",
        "frame", "K418 sum", "Causality 2Λ1/Λ2", "Λ0/Λ2")
println("  ", "-"^78)
for (name, fc) in frames
    ka = largek_analytic(fc)
    cb = causality_biquad(fc)
    @printf("  %-22s | %-12.5f | %-26.5g | %-10.4g\n",
            name, ka.sum, cb.sum, cb.prod)
end

# ---------------------------------------------------------------------------
# TABLE 3 — the decisive test of the "τ are dimensionless ratios" reading.
# Reference docstring says Λ0,Λ1,Λ2 are functions of (…,τ_ε,τ_P,τ_Q,…) where the
# τ are DIMENSIONLESS frame numbers (only RATIOS map: ŝ=τP/(cs²τε), â=τQ/τε).
# Feed Causality.jl the DIMENSIONLESS Shum ratios directly:
#   τε_d = 1,  τP_d = ŝ cs² (so τP/(cs²τε)=ŝ),  τQ_d = â (so τQ/τε=â).
# Then 2Λ1/Λ2 = cs²(1/â + ŝ cs² + 1).  Compare to K418.
# ---------------------------------------------------------------------------
println("\n", "="^80)
println("TABLE 3 — feed Causality.jl the DIMENSIONLESS ratios (τε=1, τP=ŝcs², τQ=â)")
println("="^80)
println("  Does using dimensionless τ (the only reading where 'τP-1' is meaningful)")
println("  make Causality's c∞² match Kovtun's large-k?  (q̂ has no analogue in C.)")
@printf("  %-22s | %-12s | %-12s | %-12s | %-8s\n",
        "frame (ŝ,â)", "K418 c2p", "C-dimless c2p", "K418 sum", "C sum")
println("  ", "-"^78)
e0 = W0/(1 + VS2); p0 = VS2*e0
for (name, fc) in frames
    ŝ = VS2*fc.ε1/GS; â = fc.θ/GS
    τε_d = 1.0; τP_d = ŝ*VS2; τQ_d = â
    tc = TransportCoefficients(η=ETA, ζ=ZETA, κQ=0.0, τε=τε_d, τP=τP_d, τQ=τQ_d, L=1.0)
    c2m, c2p, disc = characteristic_speeds(p0, VS, VS, tc)
    ka = largek_analytic(fc)
    sumC = (c2m + c2p)
    @printf("  %-22s | %-12.5f | %-12s | %-12.5f | %-8.4g\n",
            @sprintf("(%.1f,%.1f)", ŝ, â), ka.c2p,
            disc ≥ 0 ? @sprintf("%.5f", c2p) : @sprintf("NaN(d<0)"),
            ka.sum, sumC)
end

# ---------------------------------------------------------------------------
# TABLE 4 — the analytic identity that settles bug-vs-convention.
# Build the LARGE-k limit of Kovtun's quartic SYMBOLICALLY in c=ω/k and read off
# its biquadratic coefficients; compare their RATIO structure to Causality.jl's.
# Kovtun F_sound highest-k powers (set ω=ck, drop lower-k): from fsound_rest_coeffs
#   c4-term: cs²ε1θ ω⁴                                   (k⁰ part — no k)
#   c2-term: -k² cs² A2 ω²   with A2 = cs⁴ε1²+γsε1+(ε2+π1)(θ-cs²ε1)+ε2π1
#   c0-term: k⁴ cs² θ (cs²(ε2+π1-cs²ε1) - γs)            (B1·k⁴)
# Leading large-k balance (ω=ck): k⁴[ cs²ε1θ c⁴ - cs²A2 c² + cs²θ(cs²(ε2+π1-cs²ε1)-γs) ]=0
#   ⇒  ε1θ c⁴ - A2 c² + θ(cs²(ε2+π1-cs²ε1) - γs) = 0.
# Divide by γs² and use bars (ε̄1=cs²ε1/γs etc.) to recover eq.4.18 (K418).  So:
#   Kovtun-largek :  (ε1θ) c⁴ - (A2) c² + θ(cs²(ε2+π1-cs²ε1)-γs) = 0
#   Causality.jl  :  Λ2 c⁴ - 2Λ1 c² + Λ0 = 0
# These match iff (2Λ1)/Λ2 = A2/(ε1θ) and Λ0/Λ2 = (cs²(ε2+π1-cs²ε1)-γs)/ε1.
# ---------------------------------------------------------------------------
println("\n", "="^80)
println("TABLE 4 — Kovtun-largek biquadratic (ε1θ)c⁴-(A2)c²+θ(…)=0  coefficient RATIOS")
println("           vs Causality.jl 2Λ1/Λ2 and Λ0/Λ2 (same frame, raw dimensionful τ)")
println("="^80)
@printf("  %-22s | %-16s %-16s | %-16s %-16s\n",
        "frame", "Kov sum A2/(ε1θ)", "C 2Λ1/Λ2", "Kov prod", "C Λ0/Λ2")
println("  ", "-"^78)
for (name, fc) in frames
    ε1=fc.ε1; ε2=fc.ε2; π1=fc.π1; θ=fc.θ
    A2 = VS2^2*ε1^2 + GS*ε1 + (ε2+π1)*(θ-VS2*ε1) + ε2*π1
    kov_sum  = A2/(ε1*θ)
    kov_prod = (VS2*(ε2+π1-VS2*ε1) - GS)/ε1
    cb = causality_biquad(fc)
    @printf("  %-22s | %-16.5g %-16.5g | %-16.5g %-16.5g\n",
            name, kov_sum, cb.sum, kov_prod, cb.prod)
end

println("\n", "#"^80)
println("# DONE — see tables above for the apples-to-apples large-k comparison.")
println("#"^80)
