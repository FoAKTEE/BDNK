# =============================================================================
# viscous_fmode_scaling.jl — MEASURED: the f-mode viscous damping exponent is
# SUB-LINEAR, γ ∝ η̂^p with p ≈ 0.47–0.50, cross-method.  *** STATUS: UNDER AUDIT —
# the MECHANISM is NOT established; do not cite this as physics yet. ***
#
# ⚠️ CORRECTION to an earlier version of this header, which claimed the result was
# "the classical Landau–Lifshitz free-surface Stokes layer". THAT IS WRONG:
#   * √(νω) Stokes-layer damping is the RIGID-WALL / NO-SLIP result;
#   * for a FREE surface (what a star has) the classical answer is BULK, γ = 2νk² ⇒ p = 1.
# So the measured p≈0.5 CONTRADICTS the standard free-surface expectation rather than
# confirming it, and is not explained by textbook theory.
#
# STRUCTURAL POINT (why p=0.5 is a strong claim): the dissipation integral
#   γ = ∫2η σ:σ dV / (2E),  η = η̂(ε+p)
# evaluated on the INVISCID eigenfunction is EXACTLY LINEAR in η̂. Getting p=0.5 REQUIRES
# the eigenfunction to deform with η̂ — a layer of width δ ∝ √η̂ gives ∫σ² ∝ 1/δ, hence
# γ ∝ η̂/√η̂ = √η̂. That layer is either genuine stellar-surface physics (η=η̂(ε+p)→0 while
# ν=η̂ stays finite, c_s→0, so the constant-density free-surface analysis breaks down) OR an
# under-resolution artifact (the layer thins as √η̂; Nr=140 and Nr=200 already DISAGREE for
# η̂ ≲ 0.04 — exactly the artifact signature).
#
# RESOLUTION (2026-07-19 audit): ❌ **THE √ν CLAIM IS RETRACTED — IT IS AN ARTIFACT.**
# A dedicated audit (convergence × dissipation integral × eigenfunction diagnostic ×
# BC sensitivity, plus adversarial verification) established:
#   1. The layer the claim requires DOES NOT EXIST. p=0.5 demands a viscous layer of
#      width δ∝√η̂; measured, the dissipation is GRID-LOCKED TO EXACTLY ONE CELL at
#      every η̂ — fixed CELL width, not fixed physical width, and not ∝√η̂.
#   2. The exponent is a function of UNPHYSICAL KNOBS: 0.42→0.99 when only the outer
#      ghost rule changes; 0.11→0.64 when only the floor `den_frac` changes.
#   3. This operator never implemented the invoked mechanism: `νm` is an r-INDEPENDENT
#      CONSTANT, never η(r)=η̂(ε+p)→0.
#   4. p NEVER CONVERGES. Fit over η̂∈[0.04,0.16] vs Nr=100,140,200,260,320,400 gives
#      0.560, 0.498, 0.436, −0.169, 0.031, 0.238 — it passes THROUGH 0.5 while drifting;
#      γ still moves 34–37% from Nr=320→400. The solver is documented valid only for
#      Nr≲200 (PolarViscousModes.jl:105-106), so the original points sat at that edge.
#   5. The 1D "cross-method confirmation" was NOT independent: DynGR1D hard-resets v=0
#      in the atmosphere = a numerical NO-SLIP RIGID WALL, and √(νω) is exactly the
#      rigid-wall Stokes law. See repro/viscous_damping.jl.
#   6. An η̂-INDEPENDENT OFFSET + linear damping fits 13–18× better than a power law —
#      that offset is the actual generator of the spurious p≈0.5.
# WHAT REPLACES IT: bulk p≈1, robust across four independent routes — dissipation
# integral (exactly linear), Dirichlet BC (0.985), extrapolation BC (0.966), 3/6/12-cell
# masking (1.004/1.004/1.003), axial shooting solver (0.959) — consistent with the
# classical FREE-surface result γ=2νk². Coefficient unsettled to ~factor 1.8.
# This script is retained as a REGRESSION/DEFECT probe, not as a physics result.
#
# COROLLARY for the 2+1D (r,θ) viscous engine (SphBDNK): the operator ν_mom∇²δS_i
# dissipates the ℓ=2 mode (qualitatively robust, test_spherical Gate 5), but the
# weakly-damped fundamental's RATE cannot be cleanly read from the linearized
# time-domain (numerical grid-scale growth dominates the low-k mode). The clean
# QUANTITATIVE f-mode damping lives in this converged eigenvalue solver.
# =============================================================================
using BDNKStar, Statistics
using BDNKStar: ShumPolytrope, polar_qnm, qnm_damping, qnm_freq_kHz

const EOS = ShumPolytrope(100.0)
const ΕC  = 0.00128 + 100*0.00128^2     # the (r,θ)-engine reference star

function fmode_scaling(Nr)
    ηs = (0.005, 0.01, 0.02, 0.04, 0.08, 0.16)
    g = Float64[]
    println("  Nr=$Nr:   η̂        γ           γ/η̂        γ/√η̂")
    for η̂ in ηs
        res = polar_qnm(EOS, ΕC; l=2, η̂=η̂, Nr=Nr, nmodes=1, warn=false)
        γ = abs(qnm_damping(res[1])); push!(g, γ)
        println("          ", rpad(η̂,9), rpad(round(γ,sigdigits=4),12),
                rpad(round(γ/η̂,sigdigits=4),11), round(γ/sqrt(η̂),sigdigits=4))
    end
    # power law over the RESOLVED tail (η̂≥0.04): γ ∝ η̂^p
    hi = collect(ηs)[4:end]; gh = g[4:end]
    lx=log.(hi); ly=log.(gh); m=length(lx)
    p=(m*sum(lx.*ly)-sum(lx)*sum(ly))/(m*sum(lx.^2)-sum(lx)^2)
    println("          → resolved-regime exponent p (γ∝η̂^p, η̂≥0.04) = ", round(p,digits=3),
            "   (0.5 ⇒ √ν boundary layer; 1.0 ⇒ bulk)")
    p
end

function main()
    println("ℓ=2 f-mode viscous damping γ(η̂): free-surface boundary layer (√ν) or bulk (ν)?\n")
    p140 = fmode_scaling(140); println()
    p200 = fmode_scaling(200)
    println("\n⇒ Large-η̂ exponent ≈ ", round((p140+p200)/2, digits=2),
            "  (bulk/dissipation-integral would be 1.0). MECHANISM NOT ESTABLISHED.")
    println("  The classical FREE-surface result is BULK γ=2νk² (p=1); √(νω) is the RIGID-WALL Stokes law,")
    println("  so p≈0.5 is NOT explained by textbook theory. p=0.5 requires an η̂-DEPENDENT layer (δ∝√η̂),")
    println("  which is either stellar-surface physics (η=η̂(ε+p)→0, c_s→0) or under-resolution.")
    println("  RED FLAG: Nr=140 vs Nr=200 disagree for η̂≲0.04. Treat as UNDER AUDIT, not as a result.")
end

main()
