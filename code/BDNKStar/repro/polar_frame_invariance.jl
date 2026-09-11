# Polar viscous QNM FRAME-INDEPENDENCE test.
#
# Question: at FIXED physical shear viscosity η̂, does sweeping the BDNK FRAME
# knobs (relaxation time τ̂/cτ, conduction coupling κ̂/cκ) move the PHYSICAL
# f and p1 modes? A genuine frame-gauge bug => the physical mode tracks the knob.
# A solver artifact => a fixed offset vs the shooting benchmark, present at ALL
# frames, independent of the knob.
#
# Frame knobs in polar_bdnk_operator:
#   τR = τ̂ + cτ η̂   (relaxation time; pure-gauge baseline τ̂)
#   κQ = κ̂ + cκ η̂   (reactive conduction coupling)
#   νm = ν̂ + cν η̂   (shear-viscous DAMPER; cν is the physical-transport channel)
# Physical transport: η̂ (shear). EOS/star held fixed.

using BDNKStar
using BDNKStar: PolarViscousModes, NonRadialModes
using Printf
using Statistics: mean, std

const eos = ShumPolytrope(100.0)
const εc  = 0.00128 + 100 * 0.00128^2          # n=1 K=100 polytrope, M≈1.4
const ETA = 0.03                                # FIXED physical shear viscosity
const NR  = 140                                 # FIXED, converged window
const L   = 2

freqkHz = PolarViscousModes.qnm_freq_kHz
damp    = PolarViscousModes.qnm_damping

# Shooting (NonRadialModes) inviscid benchmark — the "true" frame-independent f/p1.
fbench, _, _ = NonRadialModes.nonradial_cowling_spectrum(eos, εc; l=L, nmodes=2,
                                                          N=4000, nscan=900)
@printf("\n=== SHOOTING BENCHMARK (NonRadialModes, inviscid) ===\n")
@printf("  f-mode  = %.4f kHz\n", fbench[1])
@printf("  p1-mode = %.4f kHz\n", length(fbench)>=2 ? fbench[2] : NaN)

# helper: get (f_f, g_f, f_p1, g_p1) for a given frame-knob set at fixed η̂
function modes(; kw...)
    ev = PolarViscousModes.polar_qnm(eos, εc; l=L, η̂=ETA, Nr=NR, nmodes=2,
                                     warn=false, kw...)
    f = [(freqkHz(e), damp(e)) for e in ev]
    while length(f) < 2; push!(f, (NaN, NaN)); end
    (f[1]..., f[2]...)
end

# also inviscid (η̂=0) per-frame, to separate reactive-η̂ pull from the gauge.
function modes0(; kw...)
    ev = PolarViscousModes.polar_qnm(eos, εc; l=L, η̂=0.0, Nr=NR, nmodes=2,
                                     warn=false, kw...)
    f = [freqkHz(e) for e in ev]
    while length(f) < 2; push!(f, NaN); end
    (f[1], f[2])
end

println("\n", "="^78)
println("SWEEP 1: RELAXATION-TIME FRAME KNOB τ̂ (pure gauge), η̂=$ETA FIXED")
println("  causal region: larger τ̂ => smaller c_+ (more causal). Sweep across it.")
println("="^78)
@printf("  %-8s | %-10s %-10s | %-10s %-10s | %-12s %-12s\n",
        "τ̂", "f_f(kHz)", "γ_f", "f_p1(kHz)", "γ_p1", "f_f η̂=0", "f_p1 η̂=0")
tauhats = [0.06, 0.09, 0.12, 0.18, 0.25, 0.35, 0.50]
ff_tau=Float64[]; gf_tau=Float64[]; fp_tau=Float64[]; gp_tau=Float64[]
ff0_tau=Float64[]; fp0_tau=Float64[]
for τ in tauhats
    ff, gf, fp, gp = modes(τ̂=τ)
    f0f, f0p = modes0(τ̂=τ)
    push!(ff_tau,ff); push!(gf_tau,gf); push!(fp_tau,fp); push!(gp_tau,gp)
    push!(ff0_tau,f0f); push!(fp0_tau,f0p)
    @printf("  %-8.2f | %-10.4f %-10.4g | %-10.4f %-10.4g | %-12.4f %-12.4f\n",
            τ, ff, gf, fp, gp, f0f, f0p)
end

println("\n", "="^78)
println("SWEEP 2: RELAXATION-COUPLING cτ (frame knob on the η̂-dependent part), η̂=$ETA")
println("="^78)
@printf("  %-8s | %-10s %-10s | %-10s %-10s\n", "cτ", "f_f(kHz)", "γ_f", "f_p1(kHz)", "γ_p1")
ctaus = [0.0, 0.5, 1.0, 2.0, 4.0]
ff_ct=Float64[]; gf_ct=Float64[]; fp_ct=Float64[]; gp_ct=Float64[]
for c in ctaus
    ff, gf, fp, gp = modes(cτ=c)
    push!(ff_ct,ff); push!(gf_ct,gf); push!(fp_ct,fp); push!(gp_ct,gp)
    @printf("  %-8.2f | %-10.4f %-10.4g | %-10.4f %-10.4g\n", c, ff, gf, fp, gp)
end

println("\n", "="^78)
println("SWEEP 3: CONDUCTION-COUPLING FRAME KNOB (κ̂ baseline, cκ η̂-part), η̂=$ETA")
println("="^78)
@printf("  %-12s | %-10s %-10s | %-10s %-10s\n", "(κ̂,cκ)", "f_f(kHz)", "γ_f", "f_p1(kHz)", "γ_p1")
kappas = [(0.0,0.0), (0.0,0.5), (0.0,1.0), (0.05,0.0), (0.1,0.0), (0.05,1.0)]
ff_k=Float64[]; gf_k=Float64[]; fp_k=Float64[]; gp_k=Float64[]
for (kh,ck) in kappas
    ff, gf, fp, gp = modes(κ̂=kh, cκ=ck)
    push!(ff_k,ff); push!(gf_k,gf); push!(fp_k,fp); push!(gp_k,gp)
    @printf("  (%-4.2f,%-4.2f) | %-10.4f %-10.4g | %-10.4f %-10.4g\n", kh,ck, ff, gf, fp, gp)
end

# ----------------------------------------------------------------------------
# ANALYSIS
# ----------------------------------------------------------------------------
spread(x) = (m=mean(filter(isfinite,x)); 100*(maximum(filter(isfinite,x))-minimum(filter(isfinite,x)))/m)
offset(x,b) = (m=mean(filter(isfinite,x)); 100*(m-b)/b)

println("\n", "="^78)
println("ANALYSIS — is the ~5-8% systematic FRAME-GAUGE or SOLVER ARTIFACT?")
println("="^78)

@printf("\n[A] PHYSICAL-MODE SPREAD across the FRAME sweep (should be ~0 if frame-indep):\n")
@printf("  τ̂  sweep: f_f spread = %.2f%%   p1 spread = %.2f%%   (η̂=%.2f)\n",
        spread(ff_tau), spread(fp_tau), ETA)
@printf("  τ̂  sweep: γ_f spread = %.2f%%   γ_p1 spread = %.2f%%\n",
        spread(gf_tau), spread(gp_tau))
@printf("  cτ sweep: f_f spread = %.2f%%   p1 spread = %.2f%%\n",
        spread(ff_ct), spread(fp_ct))
@printf("  κ̂/cκ sweep: f_f spread = %.2f%%   p1 spread = %.2f%%\n",
        spread(ff_k), spread(fp_k))
@printf("  τ̂  sweep, η̂=0: f_f spread = %.2f%%   p1 spread = %.2f%%\n",
        spread(ff0_tau), spread(fp0_tau))

@printf("\n[B] OFFSET vs SHOOTING BENCHMARK (the fixed conservative-frame systematic):\n")
@printf("  inviscid (η̂=0) f_f offset across τ̂ frames: %.2f%% .. %.2f%%\n",
        offset([minimum(ff0_tau)],fbench[1]), offset([maximum(ff0_tau)],fbench[1]))
@printf("  inviscid (η̂=0) p1  offset across τ̂ frames: %.2f%% .. %.2f%%\n",
        offset([minimum(fp0_tau)],fbench[2]), offset([maximum(fp0_tau)],fbench[2]))
@printf("  viscous  (η̂=%.2f) f_f offset across τ̂ frames: %.2f%% .. %.2f%%\n",
        ETA, offset([minimum(ff_tau)],fbench[1]), offset([maximum(ff_tau)],fbench[1]))

@printf("\n[C] VERDICT LOGIC:\n")
@printf("  frame-spread (τ̂, η̂=0)   = %.2f%%  <-- pure gauge dependence of physical f\n", spread(ff0_tau))
@printf("  benchmark offset (mean)  = %.2f%%  <-- fixed conservative-frame systematic\n",
        offset(ff0_tau, fbench[1]))
@printf("  If spread << offset and offset≈const across frames => SOLVER ARTIFACT.\n")
@printf("  If spread tracks the knob (≳ offset)             => FRAME-GAUGE bug.\n")
println("="^78)
