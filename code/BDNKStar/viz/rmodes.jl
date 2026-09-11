#=
    r-MODE (l=m=2) CFS GRAVITATIONAL-WAVE INSTABILITY WINDOW + BDNK causal ζ_eff.

    PRIMARY  Lindblom, Owen & Morsink 1998, PRL 80, 4843 [gr-qc/9803053].
    CROSS    Owen, Lindblom, Cutler, Schutz, Vecchio & Andersson 1998,
             PRD 58, 084020 [gr-qc/9804044], Table I.
    CROSS    Andersson & Kokkotas 2001 review, IJMPD 10, 381, Fig. 5.

    (A) The three r-mode timescales τ_GW, τ_sv, τ_bv vs core temperature T at a
        fixed spin (near Kepler) for the n=1 polytrope, showing the GW–shear and
        GW–bulk crossings that bound the instability window.  |τ| on a log axis;
        GW driving (τ<0) drawn dashed.
    (B) The instability window: critical spin ν_crit(T) [Hz] (and Ω_c/Ω_K) vs
        core temperature for the n=1 polytrope + realistic SLy/APR4 EOS, with the
        Navier–Stokes (τ=0) vs causal-ζ_eff (τ>0) high-T edge contrast; the
        unstable region (Ω>Ω_crit) is shaded.  The classic minimum sits near
        T~10⁹·³ K at Ω_c/Ω_K~0.06 (LOM98).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/rmodes.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
const GC = BDNKStar.Units.gram_per_cm3_to_km_minus2

"n=1 (Γ=2) relativistic TOV star near M=1.4 M⊙."
function n1_star(; K=190.0)
    eos = isentropic_idealgas(Γ=2.0, K=K, ρ_lo=1e-10, ρ_hi=1e-2, N=3000)
    best=nothing; bd=Inf
    for ρc in exp.(range(log(3e-4), log(2e-3); length=120))
        pc=K*ρc^2; s=solve_tov(eos, ρc+pc; h=0.0025)
        d=abs(mass_solar(s)-1.4); d<bd && (bd=d; best=s)
    end
    return best, eos
end

"Realistic-EOS star near 1.4 M⊙."
function eos_star(sym)
    eos = piecewise_polytrope(sym)
    best=nothing; bd=Inf
    for εc in exp.(range(log(6e14*GC), log(2.2e15*GC); length=60))
        s=solve_tov(eos, εc; h=0.004)
        d=abs(mass_solar(s)-1.4); d<bd && (bd=d; best=s)
    end
    return best, eos
end

star1, eos1 = n1_star()
ΩK1 = kepler_frequency(star1, eos1)

fig = Figure(size=(1240, 540))

# ── Panel A: timescales vs T at fixed near-Kepler spin ───────────────────────
axA = Axis(fig[1,1], xscale=log10, yscale=log10,
           xlabel="core temperature  T  [K]", ylabel="|τ|  [s]",
           title="r-mode timescales (n=1, M=1.4 M⊙, Ω=0.8 Ω_K)\nGW driving (dashed) vs shear/bulk damping (solid)")
Tg = 10 .^ range(8.0, 10.6; length=120)
Ωfix = 0.8*ΩK1
τGWv=Float64[]; τsvv=Float64[]; τbvv=Float64[]
for T in Tg
    g,s,b,_ = rmode_timescales(star1, eos1; Ω=Ωfix, T=T)
    push!(τGWv, abs(g)); push!(τsvv, s); push!(τbvv, b)
end
lines!(axA, Tg, τGWv; color=:crimson,   linewidth=2.5, linestyle=:dash, label="|τ_GW| (driving)")
lines!(axA, Tg, τsvv; color=:dodgerblue, linewidth=2.5, label="τ_sv (shear, ∝T⁻²)")
lines!(axA, Tg, τbvv; color=:seagreen,   linewidth=2.5, label="τ_bv (bulk, ∝T⁻⁶)")
axislegend(axA; position=:lb, framevisible=true, labelsize=11)

# ── Panel B: instability window ν_crit(T) — NS vs causal, EOS overlay ─────────
axB = Axis(fig[1,2], xscale=log10, yscale=log10,
           xlabel="core temperature  T  [K]", ylabel="critical spin  ν_crit  [Hz]",
           title="r-mode instability window (l=m=2)\nunstable above the curve; NS vs BDNK causal ζ_eff")
Tw = 10 .^ range(8.5, 10.7; length=90)

# n=1 polytrope, Navier–Stokes
wNS = rmode_instability_window(star1, eos1; Tgrid=Tw, τ_bulk_relax=0.0)
# shade the unstable region (between ν_crit and ν_Kepler) for the NS n=1 case
band!(axB, wNS.T, wNS.ν_crit_Hz, fill(wNS.ν_K_Hz, length(Tw));
      color=(:crimson, 0.12))
lines!(axB, wNS.T, wNS.ν_crit_Hz; color=:black, linewidth=3,
       label="n=1 (Navier–Stokes, τ=0)")

# n=1 polytrope, BDNK causal ζ_eff (τ_bulk = 1e-3 s)
wC = rmode_instability_window(star1, eos1; Tgrid=Tw, τ_bulk_relax=1e-3)
lines!(axB, wC.T, wC.ν_crit_Hz; color=:purple, linewidth=3, linestyle=:dash,
       label="n=1 (causal ζ_eff, τ=10⁻³ s)")

# realistic EOS
for (sym, col) in ((:SLy,:dodgerblue), (:APR4,:seagreen))
    s, e = eos_star(sym)
    w = rmode_instability_window(s, e; Tgrid=Tw)
    lines!(axB, w.T, w.ν_crit_Hz; color=col, linewidth=2.2,
           label=@sprintf("%s (M=%.2f, ν_K=%.0f Hz)", String(sym), mass_solar(s), w.ν_K_Hz))
end

# Kepler line + minimum marker
hlines!(axB, [wNS.ν_K_Hz]; color=:gray, linestyle=:dot, linewidth=1.5)
text!(axB, 1.2e9, wNS.ν_K_Hz*0.9; text="ν_Kepler (mass-shedding)", color=:gray, fontsize=10)
imin = argmin(wNS.ν_crit_Hz)
scatter!(axB, [wNS.T[imin]], [wNS.ν_crit_Hz[imin]]; color=:crimson, markersize=12)
text!(axB, wNS.T[imin]*1.1, wNS.ν_crit_Hz[imin]*1.15;
      text=@sprintf("min: T=%.1e K\nΩ_c/Ω_K=%.3f", wNS.T[imin], wNS.Ω_crit_over_ΩK[imin]),
      fontsize=10)
axislegend(axB; position=:rt, framevisible=true, labelsize=10)
ylims!(axB, 8, 2000)

outfile = joinpath(outdir, "rmodes_window.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# ── console: LOM98 validation table + window summary ─────────────────────────
v = rmode_validate_lom98()
pct(c,p)=100*(c-p)/abs(p)
println("\nLOM98 n=1 r-mode validation (Lane–Emden, Ω=√(πGρ̄), T=1e9 K):")
@printf("  %-6s %12s %12s %8s\n", "qty", "computed", "published", "%diff")
@printf("  %-6s %12.4e %12.4e %+7.2f%%\n", "J̃",   v.J̃,  JTILDE_PUB, pct(v.J̃,JTILDE_PUB))
@printf("  %-6s %12.4e %12.4e %+7.2f%%\n", "Ĩ",   v.Ĩ,  ITILDE_PUB, pct(v.Ĩ,ITILDE_PUB))
@printf("  %-6s %12.4e %12.4e %+7.2f%%\n", "τ_GW", v.τ_GW, TAU_GW_PUB, pct(v.τ_GW,TAU_GW_PUB))
@printf("  %-6s %12.4e %12.4e %+7.2f%%\n", "τ_sv", v.τ_sv, TAU_SV_PUB, pct(v.τ_sv,TAU_SV_PUB))
@printf("  %-6s %12.4e %12.4e %+7.2f%%\n", "τ_bv", v.τ_bv, TAU_BV_PUB, pct(v.τ_bv,TAU_BV_PUB))

σK = (2/3)*ΩK1
@printf("\nωτ regime at near-Kepler: σ=(2/3)Ω_K=%.0f s⁻¹ (%.0f Hz); ζ_eff/ζ_NS(τ=1e-3 s)=%.3f\n",
        σK, σK/(2π), zeta_eff_factor(σK, 1e-3))
println("Window minimum (n=1, NS): T=", @sprintf("%.2e K", wNS.T[imin]),
        "  Ω_c/Ω_K=", @sprintf("%.3f", wNS.Ω_crit_over_ΩK[imin]),
        "  ν_c=", @sprintf("%.1f Hz", wNS.ν_crit_Hz[imin]))
