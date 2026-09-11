#=
    viz/dyngr_1d.jl — STAGE 2 dynamical-GR figure.

    Three panels:
      (a) static-TOV stationarity: central density drift over many dynamical
          times for the FULL coupled system (unperturbed);
      (b) radial mode: the central-density spectrum for the DYNAMICAL metric vs
          the FROZEN-metric control, with the validated relativistic Cowling
          eigenvalue marked — the dynamical mode sits well below frozen/Cowling
          (the spacetime-response signature);
      (c) collapse to a black hole: α_c(t) → 0 (collapse of the lapse),
          max(2m/r) → 1 (apparent horizon), ρ_c(t) ↑ (central density grows).
=#
#   Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/dyngr_1d.jl
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar: setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!, solve_tov,
                 ShumPolytrope, dyngr_radial_freq, radial_cowling_spectrum
using .BDNKStar.Units: kHz_to_km
using CairoMakie
using Printf

const CONV = 1/kHz_to_km     # km^-1 -> kHz (geometric)

_detrend(ts, y) = begin
    n=length(ts); x=(ts.-ts[1])./(ts[end]-ts[1]); A=hcat(ones(n),x,x.^2,x.^3); y.-A*(A\y)
end

eos = ShumPolytrope(100.0)
εc  = 0.0015
star = solve_tov(eos, εc; h=2e-4)

# (a) static stationarity ----------------------------------------------------
eng, st = setup_dyngr(eos, εc; N=400, cfl=0.3)
resS = evolve_dyngr!(st, eng; tmax=20*star.R, sample_dt=eng.g.Δr)
ρc0 = resS.ρc[1]
relρ = abs.(resS.ρc .- ρc0)./ρc0

# (b) radial mode spectra: dynamical vs frozen -------------------------------
function spectrum(fm, band)
    eng, st = setup_dyngr(eos, εc; N=600, cfl=0.25, freeze_metric=fm)
    seed_dyngr_velocity!(st, eng; A=2e-4, profile=:linear)
    res = evolve_dyngr!(st, eng; tmax=250*star.R, probe_frac=0.5, sample_dt=0.5*eng.g.Δr)
    yv = _detrend(res.ts, copy(res.probe))
    fpk,P,fg = dyngr_radial_freq(res.ts, yv; fmin=0.3*kHz_to_km, fmax=8*kHz_to_km, nf=6000)
    # restrict reported peak to the band
    fpk2,_,_ = dyngr_radial_freq(res.ts, yv; fmin=band[1]*kHz_to_km, fmax=band[2]*kHz_to_km, nf=6000)
    (fg.*CONV, P./maximum(P), fpk2*CONV)
end
fgD, PD, FD = spectrum(false, (1.2,4.0))
fgF, PF, FF = spectrum(true,  (3.5,6.5))
Fcow = radial_cowling_spectrum(eos, εc; N=1500, h_tov=5e-5, nmodes=1)[1][1]

# (c) collapse ----------------------------------------------------------------
εc_c = 0.003
starc = solve_tov(eos, εc_c; h=2e-4)
engc, stc = setup_dyngr(eos, εc_c; N=500, cfl=0.3, rmax_fac=1.5)
seed_dyngr_velocity!(stc, engc; A=-0.03, profile=:linear)
resC = evolve_dyngr!(stc, engc; tmax=120*starc.R, probe_frac=0.5, sample_dt=engc.g.Δr)
tC = resC.ts ./ starc.R

# ---------------------------------------------------------------------------
fig = Figure(size=(1500,460))

ax1 = Axis(fig[1,1], xlabel="t / R", ylabel="|Δρ_c|/ρ_c",
           title="(a) static TOV stationary (full coupled GR)", yscale=log10)
lines!(ax1, resS.ts./star.R, max.(relρ,1e-16); color=:navy)
text!(ax1, 0.05, 0.92; text=@sprintf("max drift = %.1e", maximum(relρ)),
      space=:relative, align=(:left,:top))

ax2 = Axis(fig[1,2], xlabel="f [kHz]", ylabel="normalized power",
           title="(b) radial mode: metric responds")
lines!(ax2, fgD, PD; color=:crimson, label=@sprintf("dynamical  F=%.2f kHz", FD))
lines!(ax2, fgF, PF; color=:teal,    label=@sprintf("frozen     F=%.2f kHz", FF))
vlines!(ax2, [Fcow]; color=:black, linestyle=:dash, label=@sprintf("Cowling eigen %.2f", Fcow))
xlims!(ax2, 0, 8); axislegend(ax2; position=:rt)

ax3 = Axis(fig[1,3], xlabel="t / R", title="(c) collapse to a black hole")
lines!(ax3, tC, resC.αc; color=:purple, label="α_c (lapse)")
lines!(ax3, tC, resC.max2mor; color=:orange, label="max(2m/r)")
lines!(ax3, tC, resC.ρc./resC.ρc[1]./maximum(resC.ρc./resC.ρc[1]); color=:green,
       label="ρ_c (norm)")
hlines!(ax3, [1.0]; color=:gray, linestyle=:dot)
axislegend(ax3; position=:lc)

outdir = joinpath(@__DIR__, "..", "figures")
mkpath(outdir)
save(joinpath(outdir, "dyngr_1d.png"), fig)
println("wrote figures/dyngr_1d.png")
@printf("static drift=%.2e | F_dyn=%.3f F_frozen=%.3f Cowling=%.3f kHz | collapse: αc %.4f→%.4f, 2m/r %.3f→%.3f, ρc×%.2f\n",
        maximum(relρ), FD, FF, Fcow, resC.αc[1], resC.αc[end],
        resC.max2mor[1], resC.max2mor[end], resC.ρc[end]/resC.ρc[1])
