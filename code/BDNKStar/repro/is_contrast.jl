# =====================================================================
# repro/is_contrast.jl  —  node s2.is_contrast
#
# Reproduce the Chabanov–Rezzolla (arXiv:2311.13027) Israel–Stewart
# causality-fixed effective sound speed.  This is the IS side of the
# IS-vs-BDNK contrast.
#
# GROUNDING (all line numbers in
#   ref-paper/sources/arXiv-2311.13027/src/bulk_vis_in_bns_simulations.tex):
#
#  * Eq. (eq:limit_sound), line 510 — effective IS sound speed incl. bulk visc.:
#        c_s'^2 = (ζ/τ_Π)·(1/(ρ h')) + (∂p/∂e)_ρ + (1/h')·(∂p/∂ρ)_e
#    with subluminality requirement 0 ≤ c_s'^2 < 1  (eq:limit_sound, line 505).
#
#  * Eq. (eq:causal_limit), line 528 — the CAUSALITY FIX: when c_s'^2 > c_max^2,
#        τ_Π = (ζ/(ρ h')) · [ c_max^2 − (∂p/∂e)_ρ − (1/h')·(∂p/∂ρ)_e ]^{-1}
#    with free parameter 0 ≤ c_max^2 < 1 (line 533) "to explicitly ensure
#    causality".
#
#  * a', h' definitions (line 481):  a' := (p+Π)/(ρ(1+ε)),  h' := h + Π/ρ,
#    h := (e+p)/ρ  (line 476).
#
#  * IS conservative variables U (eq:bulk_implement1, line 337):
#        U = √γ ( ρW, (e+p+Π)W²v_j, (e+p+Π)W²−(p+Π)−ρW, ρWΠ )
#    — the 4th slot DΠ = ρWΠ is the extra IS variable (Maxwell–Cattaneo
#    constitutive relation eq:is_bulk, line 307: τ_Π u·∇Π = Π_NS − Π).
#
# PACKAGE REUSE:
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using Printf

# ---------------------------------------------------------------------
# Eq.30 (eq:limit_sound, line 510): effective IS sound speed.
# Inputs:
#   cs2_eq = (∂p/∂e)_ρ + (1/h')(∂p/∂ρ)_e   — the EQUILIBRIUM (adiabatic)
#            sound speed squared (the c_s'^2 with ζ/τ_Π → 0; cf. the
#            constant-Π lower bound B_Π = (∂p/∂e)_ρ + (1/h')(∂p/∂ρ)_e,
#            line 594-596).
#   ζ      — bulk viscosity
#   τΠ     — IS relaxation time
#   ρ      — rest-mass density
#   hp     — h' = h + Π/ρ  (specific enthalpy incl. bulk pressure)
# ---------------------------------------------------------------------
cs2_prime(cs2_eq, ζ, τΠ, ρ, hp) = (ζ / τΠ) * (1.0 / (ρ * hp)) + cs2_eq

# ---------------------------------------------------------------------
# Eq.31 (eq:causal_limit, line 528): causality-fixing relaxation time.
# Solve τ_Π so that the non-equilibrium part exactly fills the gap up to
# c_max^2:   (ζ/τ_Π)(1/(ρ h')) = c_max^2 − cs2_eq.
# Valid only when c_max^2 > cs2_eq (the equilibrium part must itself be
# subluminal, line 596: cs2_eq ≤ c_s'^2 < 1).
# ---------------------------------------------------------------------
function τΠ_causal(cmax2, cs2_eq, ζ, ρ, hp)
    return (ζ / (ρ * hp)) / (cmax2 - cs2_eq)
end

# ---------------------------------------------------------------------
# The limiting strategy item (iii), line 527: apply the fix only where
# the naive c_s'^2 would exceed c_max^2; otherwise keep τ_Π untouched.
# Returns (τΠ_used, cs2prime_after, fixed::Bool).
# ---------------------------------------------------------------------
function apply_causality_fix(cs2_eq, ζ, τΠ, ρ, hp; cmax2=0.9)
    cs2_naive = cs2_prime(cs2_eq, ζ, τΠ, ρ, hp)
    if cs2_naive > cmax2
        τfix = τΠ_causal(cmax2, cs2_eq, ζ, ρ, hp)
        cs2_after = cs2_prime(cs2_eq, ζ, τfix, ρ, hp)
        return (τfix, cs2_after, true)
    else
        return (τΠ, cs2_naive, false)
    end
end

# =====================================================================
# Demonstration on a realistic cold-NS barotrope (ShumPolytrope, the
# STAGE-1 EOS).  We evaluate at the central density of a TOV star so the
# equilibrium sound speed is a genuine NS value.
# =====================================================================
println("="^70)
println("Chabanov–Rezzolla 2311.13027 — IS causality-fixed sound speed")
println("="^70)

# Barotropic EOS (cold): for a barotrope (∂p/∂ρ)_e ≡ 0 and (∂p/∂e)_ρ = cs²,
# so cs2_eq = dp/de = sound_speed2(eos, e).   (Galeazzi limit, line 596.)
κ = 100.0
eos = ShumPolytrope(κ)

# central energy density of a representative star
εc  = 1.28e-3
star = solve_tov(eos, εc; h=1e-4, ptol_rel=1e-12, rmax=50.0)
ec   = star.ε[1]
ρc   = energy_from_pressure(eos, pressure(eos, ec)); # placeholder; recompute ρ below
# rest-mass density from e = ρ + p  (Γ=2 Shum):  ρ = √(p/κ)
pc   = pressure(eos, ec)
ρc   = sqrt(pc / κ)
hc   = (ec + pc) / ρc                 # h := (e+p)/ρ, line 476
cs2_eq = sound_speed2(eos, ec)        # (∂p/∂e)_ρ for a barotrope

println("Star:  M = $(round(mass_solar(star), digits=4)) Msun,",
        "  R = $(round(star.R, digits=3)) (geom. units)")
println("Central state:  e_c = $ec,  p_c = $(round(pc, sigdigits=5)),",
        "  ρ_c = $(round(ρc, sigdigits=5))")
println("  h    = $(round(hc, sigdigits=6))")
println("  cs2_eq = (∂p/∂e)_ρ = $(round(cs2_eq, sigdigits=6))   (equilibrium, subluminal)")
println()

# For the equilibrium part we take Π=0 (NS value at the limiter target),
# so h' = h + Π/ρ = h.
Π  = 0.0
hp = hc + Π/ρc

cmax2 = 0.9     # free causality parameter, 0 ≤ c_max² < 1 (line 533)

println("Causality parameter c_max² = $cmax2")
println("-"^70)
println("Scan ζ/τ_Π (bulk-viscosity / relaxation-time ratio):")
println()
@printf("%-14s %-16s %-16s %-8s\n", "ζ/τ_Π", "c_s'^2 naive", "c_s'^2 fixed", "fixed?")
println("-"^70)

# pick ζ fixed, vary τ_Π so that ζ/τ_Π spans below & above the luminal threshold
ζ = 1.0e-3
naive_super = false
fixed_all_sub = true
for ratio in (1e-4, 1e-3, 5e-3, (cmax2-cs2_eq)*ρc*hp, 1e-1, 1.0)
    τ = ζ / ratio
    cs2n = cs2_prime(cs2_eq, ζ, τ, ρc, hp)
    (τf, cs2f, fixed) = apply_causality_fix(cs2_eq, ζ, τ, ρc, hp; cmax2=cmax2)
    global naive_super  |= (cs2n ≥ 1.0)
    global fixed_all_sub &= (cs2f < 1.0)
    @printf("%-14.5g %-16.6f %-16.6f %-8s\n", ratio, cs2n, cs2f, fixed)
end
println("-"^70)

# ---------------------------------------------------------------------
# Algebraic self-consistency of the fix (the target reproduction):
# substituting Eq.31 back into Eq.30 must return EXACTLY c_max^2.
# ---------------------------------------------------------------------
ratio_super = 5.0                      # large ζ/τ_Π → naive grossly superluminal
τ_bad = ζ / ratio_super
cs2_bad = cs2_prime(cs2_eq, ζ, τ_bad, ρc, hp)
τ_fixed = τΠ_causal(cmax2, cs2_eq, ζ, ρc, hp)
cs2_recovered = cs2_prime(cs2_eq, ζ, τ_fixed, ρc, hp)
residual = abs(cs2_recovered - cmax2)

println()
println("CAUSALITY-FIX SELF-CONSISTENCY (Eq.31 → Eq.30 must give c_max²):")
@printf("  naive (ζ/τ_Π=%.1f):   c_s'^2 = %.6f   (superluminal: %s)\n",
        ratio_super, cs2_bad, cs2_bad > 1.0)
@printf("  fixed τ_Π          = %.6g\n", τ_fixed)
@printf("  recovered c_s'^2    = %.12f\n", cs2_recovered)
@printf("  target  c_max²      = %.12f\n", cmax2)
@printf("  |residual|          = %.3e\n", residual)
println()

# ---------------------------------------------------------------------
# IS conservative variable DΠ = ρWΠ (eq:bulk_implement1, line 346): the
# extra IS evolution variable that distinguishes IS from BDNK.  Show it
# is well-defined (W from z=Wv).
# ---------------------------------------------------------------------
v  = 0.1
W  = 1/sqrt(1 - v^2)
Πval = -ζ * 0.05            # an illustrative NS-like Π_NS = -ζΘ (Θ=0.05)
DΠ = ρc * W * Πval          # 4th IS conserved variable, line 346
@printf("IS extra conserved var  DΠ = ρWΠ = %.6g  (Π=%.4g, W=%.4f)\n", DΠ, Πval, W)
println()

# ---------------------------------------------------------------------
# Verdicts
# ---------------------------------------------------------------------
ok_residual = residual < 1e-12
ok_naive    = cs2_bad > 1.0          # naive branch IS superluminal
ok_fixed    = cs2_recovered < 1.0    # fixed branch IS subluminal
println("="^70)
@printf("PASS naive-superluminal-exists : %s\n", ok_naive)
@printf("PASS fixed-subluminal          : %s\n", ok_fixed)
@printf("PASS fix-recovers-cmax² (<1e-12): %s  (residual=%.2e)\n", ok_residual, residual)
allok = ok_residual && ok_naive && ok_fixed
@printf("ALL PASS                       : %s\n", allok)
println("="^70)
