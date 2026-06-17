#=
    PMP shock_instability + acaus_instab reproduction (2209.09265): the CAUSALITY
    classification that drives both figures. The max BDNK characteristic speed c₊
    (∝ √(c_s²/τ̂)) decreases with τ̂; a viscous shock (v_max) is stable iff it is
    causal — c₊ between v_max and the luminal bound. Outcomes (verified by the
    dynamical evolution in repro/pmp_shock_instab.jl):
      shock_instability (v_max=0.9): τ̂=1.5 c₊=0.94>v_max STABLE; τ̂=3 c₊=0.76<v_max → CRASH
      acaus_instab     (v_max=0.6): τ̂=1.5/0.5/0.4 subluminal-or-weak STABLE;
                                    τ̂=0.25 c₊=2.03 WILDLY superluminal → ACAUSAL CRASH
    (The v(x) profile snapshots need the full evolution runs — validated in the
    repro module; here we reproduce the causality criterion that determines them.)

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_shock.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# verified c₊(τ̂) at the shock left state (repro/pmp_shock_instab.jl) + outcomes
τ   = [0.25, 0.40, 0.50, 1.50, 3.00]
cp  = [2.026, 1.627, 1.470, 0.938, 0.762]
crash = [true, false, false, false, false]    # acaus τ̂=0.25 (and shock τ̂=3 vs v=0.9)

fig = Figure(size=(760, 470))
ax = Axis(fig[1,1], xlabel="τ̂", ylabel="c₊  (max characteristic speed)", xscale=log10,
          title="PMP shock/acaus causality: c₊(τ̂) sets stability")
# regimes
hlines!(ax, [1.0], color=:black, linestyle=:dash, linewidth=2)
text!(ax, 0.26, 1.04, text="luminal c=1", align=(:left,:bottom), fontsize=11)
hlines!(ax, [0.9], color=:crimson, linestyle=:dot, linewidth=1.8)
text!(ax, 2.0, 0.86, text="v_max=0.9 (shock_instability)", align=(:right,:top), color=:crimson, fontsize=10)
hlines!(ax, [0.6], color=:seagreen, linestyle=:dot, linewidth=1.8)
text!(ax, 2.0, 0.56, text="v_max=0.6 (acaus_instab)", align=(:right,:top), color=:seagreen, fontsize=10)
lines!(ax, τ, cp, color=:gray, linewidth=1.5)
for i in eachindex(τ)
    if crash[i] || (τ[i]==3.0)
        scatter!(ax, [τ[i]], [cp[i]], marker=:xcross, markersize=18, color=:red, strokewidth=3)
    else
        scatter!(ax, [τ[i]], [cp[i]], marker=:circle, markersize=15, color=:dodgerblue, strokecolor=:black, strokewidth=1)
    end
end
# annotate the two crash cases
text!(ax, 0.25, 2.10, text="τ̂=0.25: c₊=2.03\nWILDLY superluminal\n→ ACAUSAL CRASH", align=(:left,:bottom), color=:red, fontsize=10)
text!(ax, 3.0, 0.70, text="τ̂=3: c₊=0.76<0.9\n→ shock UNSTABLE", align=(:right,:top), color=:red, fontsize=10)
ylims!(ax, 0.45, 2.4)
# legend proxies
scatter!(ax, [NaN],[NaN], marker=:circle, color=:dodgerblue, strokecolor=:black, strokewidth=1, label="stable (causal)")
scatter!(ax, [NaN],[NaN], marker=:xcross, color=:red, strokewidth=3, label="crash (acausal / v>c₊)")
axislegend(ax, position=:rt, framevisible=true, labelsize=11)

Label(fig[0,:], "BDNKStar — reproduce PMP shock_instability/acaus_instab causality classification",
      fontsize=13, font=:bold)
outfile = joinpath(outdir, "pmp_shock_reproduction.png")
save(outfile, fig); println("saved ", outfile)
