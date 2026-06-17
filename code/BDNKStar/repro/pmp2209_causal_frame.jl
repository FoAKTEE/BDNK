# Reproduce PMP arXiv:2209.09265 ideal-gas causal-frame region
# (basis of picstabcaus / acaus_instab thresholds).
#
# GROUNDING (all from ref-paper/sources/arXiv-2209.09265/src/paper.tex):
#  - EOS (eq:EOS, line 395):      P = (Γ-1) m n e = n T
#  - specific int. energy (eq:e_defn, line 399): ε_tot = m n (1 + e)
#  - sound speed (eq:cs_sq, line 437):           c_s^2 = Γ P / ρ,  ρ ≡ ε_tot + P
#  - ω  (eq:omega, line 442):                    ω = m n P / (ε_tot ρ)
#  - simple frame constraints (eq:simple_constraints, line 506):
#        σ̂ ≤ 1/3   (linear stability)
#        τ̂ ≥ [ (Γ-1)(2 - c_s^2) + c_s^2 ] / (1 - c_s^2)   (causality)
#  - Γ→2 simpler bound (footnote, lines 510-512):  τ̂ ≥ 2 / (1 - c_s^2)
#  - second-law range  (eq:trans_coeff_ranges, line 1316): 0 < c_s^2 < 1,
#        0 < ω < 3 - 2√2 ≈ 0.2,  1 ≤ α
#
# TARGET (Table line 557 + caption line 1143, fig:acaus_instab):
#   Γ = 4/3, m = 0.1, V̂ = 4/3, σ̂ = 0, τ̂ ∈ {0.25, 0.4, 0.5, 1.5}.
#   Shockwave left state (eq:shockwave_params, line 1098): ε_L, v_L, n_L = {1, 0.6, 1}.
#   Classification: τ̂ = 0.25, 0.4, 0.5  ->  ACAUSAL (superluminal c_+ ~ 2,1.6,1.5)
#                   τ̂ = 1.5             ->  CAUSAL   (subluminal  c_+ ~ 0.9)

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar

# ---------------------------------------------------------------------------
# Predicate functions of (Γ, c_s^2, τ̂, σ̂)  [and ω for the 2nd law]
# ---------------------------------------------------------------------------

"Causality lower bound on τ̂  (eq:simple_constraints, line 506)."
tau_hat_causal_bound(Γ, cs2) = ((Γ - 1)*(2 - cs2) + cs2) / (1 - cs2)

"Simpler Γ→2 causality bound (footnote, line 511)."
tau_hat_causal_bound_simple(cs2) = 2 / (1 - cs2)

"Causality predicate: true if frame is CAUSAL (eq:simple_constraints)."
is_causal(Γ, cs2, τ̂) = τ̂ ≥ tau_hat_causal_bound(Γ, cs2)

"Linear-stability predicate (eq:simple_constraints): σ̂ ≤ 1/3."
is_stable(σ̂) = σ̂ ≤ 1/3

"Second-law predicate (eq:trans_coeff_ranges, line 1316): 0 < ω < 3-2√2."
second_law_ok(ω) = (0 < ω) && (ω < 3 - 2*sqrt(2))

# ---------------------------------------------------------------------------
# Build the figure state: Γ = 4/3, m = 0.1, left state ε_L=1, n_L=1.
# Package IdealGas uses (ρ = rest-mass density = m n, ϵ = specific int. energy e).
# ---------------------------------------------------------------------------
const Γ   = 4/3
const m   = 0.1
const σ̂  = 0.0            # Table line 557
const n_L = 1.0
const ε_L = 1.0           # total energy density of left state (eq:shockwave_params)

ρ_rest = m * n_L                       # = m n  (package ρ)
e_spec = ε_L / (m * n_L) - 1           # eq:e_defn  ε = m n (1+e)  ->  e = ε/(mn)-1

eos = IdealGas(Γ)
P   = pressure(eos, ρ_rest, e_spec)            # (Γ-1) m n e  (eq:EOS)
ε_tot = ρ_rest * (1 + e_spec)                  # = ε_L  (eq:e_defn)
ρ_tot = ε_tot + P                              # ρ ≡ ε + P  (line 371)
cs2 = sound_speed2(eos, ρ_rest, e_spec)        # package = Γ P / ρ_tot  (eq:cs_sq)
cs2_paper = Γ * P / ρ_tot                      # eq:cs_sq closed form (cross-check)
ω = m * n_L * P / (ε_tot * ρ_tot)              # eq:omega

println("=== fig:acaus_instab state (Γ=4/3, m=0.1, left state ε_L=1, n_L=1) ===")
println("P          = ", P)
println("ε_tot      = ", ε_tot, "   ρ_tot = ε+P = ", ρ_tot)
println("c_s^2 (pkg IdealGas)   = ", cs2)
println("c_s^2 (Γ P/ρ, eq:cs_sq) = ", cs2_paper, "   |Δ| = ", abs(cs2 - cs2_paper))
println("ω          = ", ω, "   (2nd law 0<ω<3-2√2≈", round(3-2*sqrt(2),digits=4), ": ",
        second_law_ok(ω), ")")
println("c_s        = ", sqrt(cs2))
println()

bound = tau_hat_causal_bound(Γ, cs2)
println("causality bound τ̂ ≥ ", bound,
        "   [eq:simple_constraints]")
println("simpler Γ→2 bound τ̂ ≥ ", tau_hat_causal_bound_simple(cs2),
        "   [footnote]")
println("σ̂ = ", σ̂, "  stable (σ̂≤1/3): ", is_stable(σ̂))
println()

# ---------------------------------------------------------------------------
# Classify the four figure τ̂ values
# ---------------------------------------------------------------------------
tauhats   = [0.25, 0.4, 0.5, 1.5]
expected  = ["ACAUSAL", "ACAUSAL", "ACAUSAL", "causal"]   # caption line 1143

println("=== classification (fig:acaus_instab) ===")
allok = true
for (τ̂, exp) in zip(tauhats, expected)
    causal = is_causal(Γ, cs2, τ̂)
    label  = causal ? "causal" : "ACAUSAL"
    ok = (label == exp)
    global allok &= ok
    println(rpad("τ̂ = $(τ̂)", 12),
            "-> ", rpad(label, 8),
            " (bound ", round(bound, digits=4), ")",
            "   expected ", rpad(exp, 8),
            ok ? "  ✓" : "  ✗  MISMATCH")
end
println()
println("ALL MATCH TARGET: ", allok)
