#=
    frame_independence.jl — the VISUAL statement of BDNK frame-independence:
    PHYSICAL (hydrodynamic) observables are FLAT vs the frame knob, while the
    NON-HYDRODYNAMIC / frame ("gauge") modes are SLOPED (∝ 1/τ). The contrast
    between the flat and the sloped curves IS the result.

    Panels
      (A) DISPERSION (Kovtun 1907.08191 §4.2 sound channel; recomputed LIVE).
          At fixed physical (vs, η, ζ) — ζ pinned via fn-8 — sweep the frame
          â=θ/γs over a 60× range: small-k sound speed Re ω/k → vs and leading
          damping Γ → γs are FLAT; the non-hydro gap w0/θ is SLOPED ∝ 1/â.
      (B) AXIAL viscous QNM (Bussières 2604.13208). At fixed physical η_c=1e31
          cgs, sweep frame τ̂: the physical w-mode f & 1/τ are FLAT; the non-hydro
          η-mode second-sound speed c_η(0) ∝ 1/√τ̂ is SLOPED (gauge).
      (C) POLAR viscous QNM (Cowling BDNK). At fixed physical η̂=0.03 the f/p1
          frequencies are FLAT across the frame sweep, sitting at a FIXED solver
          offset vs the shooting benchmark (dashed) — a solver artifact, not a
          frame gauge.

    Panel A is computed live; B & C plot the converged values reported by
    repro/axial_frame_independence.jl and repro/polar_frame_invariance.jl.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/frame_independence.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
using LinearAlgebra: eigvals
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# ----------------------------------------------------------------------
# (A) DISPERSION — live quartic (mirrors repro/kovtun_sound.jl, rest frame)
# ----------------------------------------------------------------------
function polyroots(c::Vector{ComplexF64})
    while length(c) > 1 && abs(c[end]) < 1e-300; c = c[1:end-1]; end
    n = length(c) - 1; n == 0 && return ComplexF64[]
    a = c ./ c[end]; C = zeros(ComplexF64, n, n)
    for i in 1:n-1; C[i+1, i] = 1.0; end
    for i in 1:n;   C[i, n] = -a[i]; end
    return eigvals(C)
end
fsound_rest(k2; cs2, ε1, ε2, π1, θ, γs, w0) = ComplexF64[
    k2*cs2*( w0^2 + k2*θ*( cs2*(ε2+π1 - cs2*ε1) - γs ) ),
    -im*k2*w0*( γs + cs2^2*ε1 + cs2*θ ),
    -( w0^2 + k2*cs2*( cs2^2*ε1^2 + γs*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1 ) ),
    im*w0*(cs2*ε1 + θ),
    cs2*ε1*θ ]

const VS=0.5; const VS2=VS^2; const W0=1.0; const ETA=0.3
const ZETA=(4.0/3.0)*1.0 - (4.0/3.0)*ETA; const GS=(4.0/3.0)*ETA + ZETA   # γs=1
frame(ŝ, â, p̂; ε2=0.0) = (ε1=ŝ*GS/VS2, ε2=ε2, π1=p̂*GS,
                          π2=VS2*(p̂*GS - VS2*ŝ*GS/VS2) - ZETA + VS2*ε2, θ=â*GS)
function sound_branch(fc; k=1e-3)
    rts = polyroots(fsound_rest(k^2; cs2=VS2, ε1=fc.ε1, ε2=fc.ε2, π1=fc.π1,
                                θ=fc.θ, γs=GS, w0=W0))
    p = sortperm(rts, by=z->abs(imag(z))); g1, g2 = rts[p[1]], rts[p[2]]
    ((abs(real(g1))+abs(real(g2)))/2/k, -2*((imag(g1)+imag(g2))/2)*W0/k^2)
end

âs = [1.0, 2.0, 4.0, 8.0, 15.0, 30.0, 60.0]
fcsA = [frame(3.0, â, 8.0) for â in âs]
vs_A = Float64[]; Γ_A = Float64[]
for fc in fcsA; v, g = sound_branch(fc); push!(vs_A, v); push!(Γ_A, g); end
gap_A = [W0/fc.θ for fc in fcsA]        # non-hydro gap, ∝ 1/â (gauge)

# ----------------------------------------------------------------------
# (B) AXIAL — converged values (repro/axial_frame_independence.jl)
# ----------------------------------------------------------------------
ax_tau = [5.0, 10.0, 20.0, 40.0, 80.0]
ax_wf  = [10.08276, 10.08982, 10.09286, 10.09425, 10.09493]   # w-mode f [kHz], A
ax_wg  = [0.032395, 0.032389, 0.032382, 0.032378, 0.032376]   # w-mode 1/τ [1/μs], A
ax_finv = 10.50135                                            # inviscid ref
ce_tau = [10.0, 20.0, 40.0, 80.0]
ce     = [0.2146, 0.1518, 0.1073, 0.0759]                     # η-mode c_η(0) [1/km], B ∝1/√τ̂

# ----------------------------------------------------------------------
# (C) POLAR — converged values (repro/polar_frame_invariance.jl)
# ----------------------------------------------------------------------
po_tau = [0.06, 0.12, 0.50]
po_ff0 = [2.0125, 2.0118, 1.9996]     # f-mode (η̂=0) — flat
po_tp1 = [0.06, 0.50]
po_fp0 = [3.9798, 4.0109]             # p1-mode (η̂=0) — flat
f_bench, p1_bench = 1.8829, 4.1067    # shooting benchmark (frame-indep reference)

# ----------------------------------------------------------------------
# FIGURE
# ----------------------------------------------------------------------
fig = Figure(size=(1500, 540))
cFLAT = :seagreen; cFLAT2 = :dodgerblue; cGAUGE = :crimson; cBENCH = :black

# Panel A: dispersion
axA = Axis(fig[1,1], xlabel="frame knob  â = θ/γs  (log)", xscale=log10,
           title="(A) Dispersion (Kovtun §4.2)\nhydro FLAT · non-hydro gap ∝ 1/τ")
lines!(axA, âs, vs_A; color=cFLAT, linewidth=2.5)
scatter!(axA, âs, vs_A; color=cFLAT, markersize=9, label="sound speed Re ω/k → vs (FLAT)")
lines!(axA, âs, Γ_A; color=cFLAT2, linewidth=2.5)
scatter!(axA, âs, Γ_A; color=cFLAT2, markersize=9, label="damping Γ → γs (FLAT)")
lines!(axA, âs, gap_A; color=cGAUGE, linewidth=2.5, linestyle=:dash)
scatter!(axA, âs, gap_A; color=cGAUGE, marker=:diamond, markersize=10,
         label="non-hydro gap w0/θ ∝ 1/τ (GAUGE)")
axislegend(axA; position=:rc, labelsize=10, framevisible=true)
ylims!(axA, -0.05, 1.4)

# Panel B: axial — w-mode flat, η-mode sloped (twin axis for differing scales)
axB = Axis(fig[1,2], xlabel="frame knob  τ̂", xscale=log10,
           ylabel="w-mode f [kHz]",
           title="(B) Axial viscous QNM (Bussières)\nw-mode FLAT · η-mode c_η ∝ 1/√τ̂")
lines!(axB, ax_tau, ax_wf; color=cFLAT, linewidth=2.5)
scatter!(axB, ax_tau, ax_wf; color=cFLAT, markersize=9, label="w-mode f (FLAT)")
hlines!(axB, [ax_finv]; color=cBENCH, linestyle=:dot, linewidth=1.5,
        label="inviscid w-mode ref")
ylims!(axB, 10.05, 10.55)
axB2 = Axis(fig[1,2], xscale=log10, yaxisposition=:right,
            ylabel="η-mode c_η(0) [1/km]  (gauge)", ygridvisible=false)
hidespines!(axB2); hidexdecorations!(axB2)
lines!(axB2, ce_tau, ce; color=cGAUGE, linewidth=2.5, linestyle=:dash)
scatter!(axB2, ce_tau, ce; color=cGAUGE, marker=:diamond, markersize=10,
         label="η-mode c_η(0) ∝ 1/√τ̂ (GAUGE)")
ylims!(axB2, 0.0, 0.26)
linkxaxes!(axB, axB2)
axislegend(axB; position=:lc, labelsize=9, framevisible=true)
axislegend(axB2; position=:rt, labelsize=9, framevisible=true)

# Panel C: polar — f/p1 flat vs frame, fixed offset vs benchmark
axC = Axis(fig[1,3], xlabel="frame knob  τ̂",
           ylabel="frequency [kHz]",
           title="(C) Polar viscous QNM (Cowling BDNK)\nf/p1 FLAT · fixed solver offset (not gauge)")
lines!(axC, po_tau, po_ff0; color=cFLAT, linewidth=2.5)
scatter!(axC, po_tau, po_ff0; color=cFLAT, markersize=9, label="f-mode (FLAT)")
hlines!(axC, [f_bench]; color=cFLAT, linestyle=:dash, linewidth=1.8,
        label="f shooting benchmark")
lines!(axC, po_tp1, po_fp0; color=cFLAT2, linewidth=2.5)
scatter!(axC, po_tp1, po_fp0; color=cFLAT2, markersize=9, label="p1-mode (FLAT)")
hlines!(axC, [p1_bench]; color=cFLAT2, linestyle=:dash, linewidth=1.8,
        label="p1 shooting benchmark")
text!(axC, 0.28, 1.94; text="fixed +6.7% offset\n(solver, not frame)",
      color=cGAUGE, fontsize=10, align=(:center,:top))
axislegend(axC; position=:rc, labelsize=9, framevisible=true)
ylims!(axC, 1.8, 4.2)

outfile = joinpath(outdir, "frame_independence.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# console summary
spread(v) = (maximum(v)-minimum(v))/abs(sum(v)/length(v))*100
println("\n=== FRAME-INDEPENDENCE (visual) ===")
@printf("(A) dispersion : Re ω/k spread = %.4f%%, Γ spread = %.4f%% (FLAT);  gap 60× ∝1/τ\n",
        spread(vs_A), spread(Γ_A))
@printf("(B) axial      : w-mode f spread = %.3f%% (FLAT);  c_η drops %.4f→%.4f ∝1/√τ̂ (GAUGE)\n",
        spread(ax_wf), maximum(ce), minimum(ce))
@printf("(C) polar      : f-mode frame spread = %.3f%% (FLAT);  fixed +%.1f%% offset vs benchmark\n",
        spread(po_ff0), 100*(sum(po_ff0)/length(po_ff0)-f_bench)/f_bench)
