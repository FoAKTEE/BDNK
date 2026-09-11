#=
    CONVERGENCE + NUMERICAL-vs-PHYSICAL VISCOSITY synthesis figure.

    Four panels distilling the convergence studies of the BDNKStar solver tracks
    and the time-domain SphBDNK physical-vs-numerical discrimination:

    (A) SphBDNK Γ vs σ_ko at η̂=0 vs η̂>0.
        PHYSICAL  = flat line (rate independent of the numerical KO knob).
        NUMERICAL = sloped line.  The η̂=0 baseline is a steeply σ_ko-dependent
        grid-scale GROWTH (Γ<0); adding shear viscosity FLATTENS the slope —
        the rate becomes more physical (but not perfectly flat ⇒ "mixed").

    (B) SphBDNK Γ vs spatial resolution Nr at fixed η̂, with the η̂=0 NUMERICAL
        FLOOR shrinking with resolution.  Honest: the floor weakens but does not
        cleanly vanish at scheme order.

    (C) Γ(η̂) scaling at fixed (Nr,σ_ko): the energy ratio / damping grows with
        η̂; we overlay the inviscid floor Γ₀ so (Γ−Γ₀)>0 is visible.

    (D) Solver convergence panel (log-log): the RMode LOM98 Lane–Emden structure
        integral J̃ converges 2nd-order in N (Richardson p≈2), and the polar-
        viscous f-mode frequency is stable across its Nr≲200 window.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/convergence_viscous.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const EOS = ShumPolytrope(100.0)
const ΕC  = 0.00128 + 100 * 0.00128^2

# ── helpers ──────────────────────────────────────────────────────────────────
"Least-squares exponential rate Γ of an energy series: E≈E0·exp(−Γt). Γ>0 decay."
function energy_rate(ts, en)
    m = (en .> 0) .& isfinite.(en)
    t = collect(ts[m]); y = log.(en[m]); n = length(t); n < 3 && return NaN
    t̄=sum(t)/n; ȳ=sum(y)/n
    return -(sum((t.-t̄).*(y.-ȳ)) / sum((t.-t̄).^2))
end

function sph_gamma(s; η̂, σ_ko, T=120.0, dt=0.01)
    ev = setup_sphbdnk(s; η̂=η̂, σ_ko=σ_ko)
    st = SphBDNKState(s.grid.Nr, s.grid.Nθ); seed_sphbdnk_n!(st, ev, 3; A=1e-3)
    nst = round(Int, T/dt)
    ts, _, en = evolve_sphbdnk!(st, ev; dt=dt, nsteps=nst, sample=max(1, nst÷200))
    return energy_rate(ts, en)
end

println("building SphBDNK Nr=48 / Nr=64 backgrounds …")
s48 = build_sphstar(EOS, ΕC; Nr=48, Nθ=10)
s64 = build_sphstar(EOS, ΕC; Nr=64, Nθ=12)

fig = Figure(size=(1240, 980))

# ── Panel A: Γ vs σ_ko at η̂=0 vs η̂>0 ────────────────────────────────────────
println("Panel A: Γ vs σ_ko …")
σgrid = [0.005, 0.01, 0.02, 0.04]
Γ_eta0 = [sph_gamma(s48; η̂=0.0,  σ_ko=σ) for σ in σgrid]
Γ_etav = [sph_gamma(s48; η̂=0.08, σ_ko=σ) for σ in σgrid]

axA = Axis(fig[1,1], xlabel="Kreiss–Oliger σ_ko", ylabel="energy rate  Γ   (Γ>0 = decay)",
           title="(A) SphBDNK Γ vs σ_ko (Nr=48)\nphysical = flat, numerical = sloped")
hlines!(axA, [0.0]; color=:gray, linestyle=:dash, linewidth=1)
scatterlines!(axA, σgrid, Γ_eta0; color=:crimson,    markersize=10, linewidth=2,
              label="η̂=0 (numerical floor: steep, GROWTH)")
scatterlines!(axA, σgrid, Γ_etav; color=:dodgerblue, markersize=10, linewidth=2,
              label="η̂=0.08 (viscous: flatter ⇒ physical)")
axislegend(axA; position=:rb, labelsize=10, framevisible=true)

# ── Panel B: Γ vs resolution at fixed η̂, plus η̂=0 floor ──────────────────────
println("Panel B: Γ vs resolution …")
Γ0_48 = sph_gamma(s48; η̂=0.0,  σ_ko=0.005)
Γ0_64 = sph_gamma(s64; η̂=0.0,  σ_ko=0.005)
Γv_48 = sph_gamma(s48; η̂=0.08, σ_ko=0.02)
Γv_64 = sph_gamma(s64; η̂=0.08, σ_ko=0.02)
Nrs = [48, 64]

axB = Axis(fig[1,2], xlabel="radial resolution  Nr", ylabel="energy rate  Γ",
           title="(B) SphBDNK Γ vs resolution\nη̂=0 numerical floor SHRINKS with Nr",
           xticks=Nrs)
hlines!(axB, [0.0]; color=:gray, linestyle=:dash, linewidth=1)
scatterlines!(axB, Nrs, [Γ0_48, Γ0_64]; color=:crimson, markersize=11, linewidth=2,
              label="η̂=0 floor (|Γ| weakens)")
scatterlines!(axB, Nrs, [Γv_48, Γv_64]; color=:dodgerblue, markersize=11, linewidth=2,
              label="η̂=0.08 (net decay)")
axislegend(axB; position=:rc, labelsize=10, framevisible=true)

# ── Panel C: Γ(η̂) scaling vs the inviscid floor ──────────────────────────────
println("Panel C: Γ(η̂) scaling …")
ηgrid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10]
Γ_of_η = [sph_gamma(s48; η̂=η, σ_ko=0.02) for η in ηgrid]
Γ0 = Γ_of_η[1]

axC = Axis(fig[2,1], xlabel="shear viscosity  η̂", ylabel="energy rate  Γ",
           title="(C) SphBDNK Γ(η̂) at Nr=48, σ_ko=0.02\n(Γ−Γ₀) > 0 and grows with η̂")
hlines!(axC, [0.0]; color=:gray, linestyle=:dash, linewidth=1)
hlines!(axC, [Γ0];  color=:crimson, linestyle=:dot, linewidth=1.5,
        label=@sprintf("Γ₀ = %.2e (inviscid floor)", Γ0))
scatterlines!(axC, ηgrid, Γ_of_η; color=:purple, markersize=10, linewidth=2,
              label="Γ(η̂)")
axislegend(axC; position=:lt, labelsize=10, framevisible=true)

# ── Panel D: solver convergence (log-log) ─────────────────────────────────────
println("Panel D: solver convergence …")
Ns = [500, 1000, 2000, 4000, 8000, 16000]
Jref = rmode_validate_lom98(N=32000).J̃
Jerr = [abs(rmode_validate_lom98(N=N).J̃ - Jref) for N in Ns]

axD = Axis(fig[2,2], xlabel="resolution N", ylabel="error",
           xscale=log10, yscale=log10,
           title="(D) Solver convergence\nRMode LOM98 J̃: 2nd-order (Richardson p≈2)")
scatterlines!(axD, Ns, Jerr; color=:seagreen, markersize=10, linewidth=2,
              label="|J̃(N) − J̃(32000)|")
# reference 2nd-order slope guide through the first point
g2 = Jerr[1] .* (Float64.(Ns) ./ Ns[1]).^(-2)
lines!(axD, Ns, g2; color=:black, linestyle=:dash, linewidth=1.5, label="∝ N⁻² (theory)")
axislegend(axD; position=:lb, labelsize=10, framevisible=true)

# polar-viscous f-mode stability annotation (text inset — cheap, robust)
f96  = qnm_freq_kHz(polar_qnm(EOS, ΕC; l=2, η̂=0.0, Nr=96,  nmodes=1)[1])
f160 = qnm_freq_kHz(polar_qnm(EOS, ΕC; l=2, η̂=0.0, Nr=160, nmodes=1)[1])
text!(axD, 600, minimum(Jerr)*3;
      text=@sprintf("Polar f-mode (Nr≲200 window):\nf(96)=%.3f kHz, f(160)=%.3f kHz\nΔ=%.1f%% (stable)",
                    f96, f160, 100*abs(f160-f96)/f96),
      align=(:left,:bottom), fontsize=9, color=:gray25)

outfile = joinpath(outdir, "convergence_viscous.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# console summary
println("\n── SphBDNK discrimination (Nr=48) ──")
@printf("σ_ko-slope of Γ:  η̂=0 → %.2e ;  η̂=0.08 → %.2e   (viscosity flattens)\n",
        abs(Γ_eta0[end]-Γ_eta0[1]), abs(Γ_etav[end]-Γ_etav[1]))
@printf("η̂=0 floor:  Nr48 Γ=%.2e ,  Nr64 Γ=%.2e   (|Γ| shrinks with Nr)\n", Γ0_48, Γ0_64)
@printf("Γ(η̂)−Γ₀:  η̂=0.06 → %.2e ,  η̂=0.10 → %.2e   (grows with η̂)\n",
        Γ_of_η[4]-Γ0, Γ_of_η[6]-Γ0)
println("RMode J̃ error N=500→16000: ", join((@sprintf("%.1e",e) for e in Jerr), ", "))
