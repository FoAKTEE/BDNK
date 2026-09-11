#=
    GENERAL / REALISTIC neutron-star EOS — Read–Lackey–Owen–Friedman (2009)
    piecewise-polytrope M–R sequences [arXiv:0812.2163, PRD 79, 124032].

    A fixed SLy crust joined to 3 high-density polytropes; each EOS fixed by
    (log10 p1, Γ1, Γ2, Γ3) from Read Table III. The TOV M–R curves for the four
    presets (SLy, APR4, H4, MS1) are overlaid with their 1.4 M⊙ points and the
    observed-NS maximum-mass band (PSR J0740+6620, Fonseca et al. 2021:
    M = 2.08 ± 0.07 M⊙) marked. Published Read 2009 / Lackey–Wade values are
    plotted as ✕ for direct visual validation.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/general_eos_mr.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.PiecewisePolytrope: PUBLISHED
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const GCGS = BDNKStar.Units.gram_per_cm3_to_km_minus2

"""TOV M–R sequence over central energy density; returns (R[km], M[M⊙], M_max, R14)."""
function mr_sequence(eos; εlo=3.5e14*GCGS, εhi=4.5e15*GCGS, n=70, h=0.008)
    εs = exp.(range(log(εlo), log(εhi); length=n))
    Rs = Float64[]; Ms = Float64[]
    for εc in εs
        st = solve_tov(eos, εc; h=h, rmax=60.0)
        push!(Rs, st.R); push!(Ms, mass_solar(st))
    end
    imax = argmax(Ms)
    R14 = NaN
    for k in 1:imax-1
        if (Ms[k]-1.4)*(Ms[k+1]-1.4) ≤ 0
            t=(1.4-Ms[k])/(Ms[k+1]-Ms[k]); R14 = Rs[k]+t*(Rs[k+1]-Rs[k]); break
        end
    end
    return Rs, Ms, Ms[imax], R14
end

presets = (:SLy, :APR4, :H4, :MS1)
colors  = Dict(:SLy=>:dodgerblue, :APR4=>:seagreen, :H4=>:darkorange, :MS1=>:crimson)

fig = Figure(size=(820, 660))
ax  = Axis(fig[1,1], xlabel="R  [km]", ylabel="M  [M⊙]",
           title="Read et al. (2009) piecewise-polytrope M–R sequences\n(SLy crust + 3 high-density polytropes; TOV)")

# observed max-mass band: PSR J0740+6620, M = 2.08 ± 0.07 M⊙ (Fonseca et al. 2021)
band!(ax, [8.5, 16.0], fill(2.08-0.07, 2), fill(2.08+0.07, 2);
      color=(:gray, 0.25))
text!(ax, 8.7, 2.09; text="PSR J0740+6620  2.08 ± 0.07 M⊙", align=(:left,:bottom),
      fontsize=11, color=:gray30)

for sym in presets
    eos = piecewise_polytrope(sym)
    Rs, Ms, Mmax, R14 = mr_sequence(eos)
    c = colors[sym]
    lines!(ax, Rs, Ms; color=c, linewidth=2.5,
           label=@sprintf("%s  (Mₘₐₓ=%.2f, R₁.₄=%.1f km)", String(sym), Mmax, R14))
    # 1.4 M⊙ point (computed)
    scatter!(ax, [R14], [1.4]; color=c, marker=:circle, markersize=12, strokewidth=1, strokecolor=:black)
    # published Read 2009 / Lackey–Wade reference (✕)
    pub = PUBLISHED[sym]
    scatter!(ax, [pub.R_14], [1.4]; color=c, marker=:xcross, markersize=14)
    scatter!(ax, [Rs[argmax(Ms)]], [Mmax]; color=c, marker=:utriangle, markersize=11)
end

hlines!(ax, [1.4]; color=:gray70, linestyle=:dash, linewidth=1)
text!(ax, 15.6, 1.40; text="1.4 M⊙", align=(:right,:bottom), fontsize=11, color=:gray50)
xlims!(ax, 8.5, 16.0); ylims!(ax, 0.0, 3.0)
axislegend(ax; position=:rb, framevisible=true, labelsize=11)

# annotate the marker legend
text!(ax, 8.7, 0.18;
      text="● computed R₁.₄    ✕ published R₁.₄ (Read 2009 / Lackey–Wade)    ▲ Mₘₐₓ",
      align=(:left,:bottom), fontsize=10, color=:gray20)

outfile = joinpath(outdir, "general_eos_mr.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# console summary (validation table)
println("\nEOS    M_max(got)  M_max(pub)   R_1.4(got)  R_1.4(pub)")
for sym in presets
    eos = piecewise_polytrope(sym)
    _, _, Mmax, R14 = mr_sequence(eos)
    pub = PUBLISHED[sym]
    @printf("%-5s   %8.3f    %8.3f    %8.3f    %8.3f\n",
            sym, Mmax, pub.M_max, R14, pub.R_14)
end
