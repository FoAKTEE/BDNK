#=
    UNIVERSAL RELATIONS ACROSS CENTRAL-DENSITY SEQUENCES AND ACROSS EOS.

    Three published quasi-EOS-independent neutron-star relations, each computed on
    realistic Read et al. (2009) piecewise-polytrope sequences (SLy/APR4/H4[/MS1])
    and overlaid on the published fit:

      (A) Andersson-Kokkotas (1998) f-mode frequency relation
          [MNRAS 299, 1059; coeffs from Kokkotas, Apostolatos & Andersson 1999,
           gr-qc/9901072 App.A Eq.(A1)]:
              f_f[kHz] = 0.78 + 1.635 * sqrt(Mbar/Rbar^3),  Mbar=M/1.4Msun, Rbar=R/10km.
          Plotted: our relativistic-Cowling f-mode (NonRadialModes) dense
          sequence + full-GR (GW-damped) anchors (PolarGRModes.polar_gr_qnm).

      (B) Breu & Rezzolla (2016) moment-of-inertia relation
          [MNRAS 459, 646, arXiv:1601.06083, Eq.(20)/Table 2 slow-rot. row]:
              Ibar = I/M^3 = 0.8134/C + 0.2101/C^2 + 3.175e-3/C^3 - 2.717e-4/C^4.

      (C) Lindblom, Owen & Morsink (1998) r-mode (l=m=2) CFS instability window
          [PRL 80, 4843]: window-minimum critical spin Omega_c/Omega_K, quasi-EOS-
          independent band ~0.05-0.06.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/universal_relations.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.PolarGRModes
const PG = BDNKStar.PolarGRModes
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
const GCGS = BDNKStar.Units.gram_per_cm3_to_km_minus2

colors = Dict(:SLy=>:dodgerblue, :APR4=>:seagreen, :H4=>:darkorange, :MS1=>:crimson)

ak_f(M,R)     = (Mb=M/1.4; Rb=R/10.0; 0.78 + 1.635*sqrt(Mb/Rb^3))
ak_invtau(M,R)= (Mb=M/1.4; Rb=R/10.0; (Mb^3/Rb^4)*(22.85 - 14.65*(Mb/Rb)))

fig = Figure(size=(1500, 480))

# ============================================================================
# Panel A: AK1998 f-mode  f vs sqrt(Mbar/Rbar^3)
# ============================================================================
axA = Axis(fig[1,1],
    xlabel="√(M̄ / R̄³)   (M̄=M/1.4M⊙, R̄=R/10km)", ylabel="f-mode frequency  f [kHz]",
    title="(A) f-mode AK1998 universal relation\n(Cowling sequence + full-GR anchors)")

windowsA = Dict(:SLy=>(9e14*GCGS, 2.3e15*GCGS),
                :APR4=>(9e14*GCGS, 2.6e15*GCGS),
                :H4 =>(6e14*GCGS, 1.7e15*GCGS))

# AK1998 published fit line over the x-range we cover
xline = range(0.45, 1.30; length=200)
lines!(axA, collect(xline), 0.78 .+ 1.635 .* xline;
       color=:black, linewidth=2.5, label="AK1998 fit (gr-qc/9901072)")

for sym in (:SLy, :APR4, :H4)
    eos = piecewise_polytrope(sym); εlo, εhi = windowsA[sym]
    xs = Float64[]; fs = Float64[]
    for εc in exp.(range(log(εlo), log(εhi); length=8))
        s = solve_tov(eos, εc; h=0.01); M = mass_solar(s); R = s.R
        M < 1.0 && continue
        f,_,_ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=1, N=2800,
                    nscan=450, ω2lo=1e-4, ω2hi=0.2)
        push!(xs, sqrt((M/1.4)/(R/10)^3)); push!(fs, f[1])
    end
    scatter!(axA, xs, fs; color=colors[sym], markersize=9,
             label="$(String(sym)) Cowling-f")
end

# full-GR (GW-damped) anchors: SLy + APR4
for (sym, εc) in [(:SLy, 1.5e15*GCGS), (:APR4, 1.6e15*GCGS)]
    eos = piecewise_polytrope(sym); s = solve_tov(eos, εc; h=0.008)
    M, R = mass_solar(s), s.R
    res = PG.polar_gr_qnm(eos, εc; l=2, ω0=PG.polar_ftau_to_omega(2.2, 0.18), star=s)
    x = sqrt((M/1.4)/(R/10)^3)
    scatter!(axA, [x], [res.f_kHz]; color=colors[sym], marker=:star5,
             markersize=22, strokecolor=:black, strokewidth=1.0,
             label="$(String(sym)) full-GR-f")
    @printf("full-GR-f %-4s  M=%.3f R=%.2f  f=%.4f (AK %.4f)  1/τ=%.3f (AK %.3f)\n",
            String(sym), M, R, res.f_kHz, ak_f(M,R), 1/res.tau_s, ak_invtau(M,R))
end
axislegend(axA; position=:lt, framevisible=true, labelsize=10, nbanks=2)

# ============================================================================
# Panel B: Breu-Rezzolla  Ibar vs C
# ============================================================================
axB = Axis(fig[1,2], xlabel="compactness  C = M/R", ylabel="Ī = I / M³",
    title="(B) Moment of inertia Ī–C\n(slow rotation; Breu & Rezzolla 2016)")

windowsB = Dict(:SLy=>(6e14*GCGS, 2.4e15*GCGS),
                :APR4=>(6e14*GCGS, 2.6e15*GCGS),
                :H4 =>(4e14*GCGS, 1.8e15*GCGS),
                :MS1=>(2.5e14*GCGS, 1.4e15*GCGS))

Cgrid = range(0.08, 0.32; length=200)
lines!(axB, collect(Cgrid), IBAR_FROM_C_BREU.(Cgrid);
       color=:black, linewidth=2.5, label="Breu & Rezzolla (2016)")

for sym in (:SLy, :APR4, :H4, :MS1)
    eos = piecewise_polytrope(sym); εlo, εhi = windowsB[sym]
    Cs = Float64[]; Ib = Float64[]
    for εc in exp.(range(log(εlo), log(εhi); length=14))
        r = moment_of_inertia(eos, εc; h=0.01)
        r.M_Msun < 1.0 && continue
        push!(Cs, r.C); push!(Ib, r.Ibar)
    end
    scatter!(axB, Cs, Ib; color=colors[sym], markersize=8, label=String(sym))
end
axislegend(axB; position=:rt, framevisible=true, labelsize=10)
xlims!(axB, 0.08, 0.32)

# ============================================================================
# Panel C: r-mode window minimum  Omega_c/Omega_K vs M
# ============================================================================
axC = Axis(fig[1,3], xlabel="gravitational mass  M [M⊙]",
    ylabel="window minimum  Ω_c / Ω_K",
    title="(C) r-mode CFS window minimum\n(LOM98 near-universal band)")

# LOM98 universal band 0.05-0.06
band!(axC, [1.05, 1.65], fill(0.05, 2), fill(0.06, 2);
      color=(:gray, 0.25))
lines!(axC, [1.05, 1.65], [0.055, 0.055]; color=:black, linestyle=:dash,
       linewidth=1.5, label="LOM98 band 0.05–0.06")

Tg = 10 .^ range(8.5, 10.5; length=120)
for sym in (:SLy, :APR4, :H4)
    eos = piecewise_polytrope(sym)
    Ms = Float64[]; mins = Float64[]
    for Mt in (1.2, 1.4, 1.6)
        best=nothing; bd=Inf
        for εc in exp.(range(log(6e14*GCGS), log(2.4e15*GCGS); length=36))
            s = solve_tov(eos, εc; h=0.008); d = abs(mass_solar(s)-Mt)
            if d<bd; bd=d; best=s; end
        end
        star = best
        w = rmode_instability_window(star, eos; Tgrid=Tg)
        imin = argmin(w.Ω_crit_over_ΩK)
        push!(Ms, mass_solar(star)); push!(mins, w.Ω_crit_over_ΩK[imin])
    end
    scatterlines!(axC, Ms, mins; color=colors[sym], markersize=11,
                  linewidth=2, label=String(sym))
end
axislegend(axC; position=:rt, framevisible=true, labelsize=10)
ylims!(axC, 0.04, 0.075)

outfile = joinpath(outdir, "universal_relations.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)
