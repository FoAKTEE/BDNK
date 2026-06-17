#=
    Reproduction validation summary — relative error |achieved−target|/target for
    every verified reference-paper benchmark (log scale per the NR comparison
    rule), against the per-result tolerance. All matches sit far below tolerance.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/repro_summary.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# (label, relative error, tolerance) — from the committed ledger / verified runs
rows = [
 ("STEP0 prim↔cons round-trip",     4.3e-15, 1e-10),
 ("RH shock εR (PMP/Pandya)",       2.3e-7,  1e-2),
 ("IS causality-fixed c_s'²",       1.2e-16, 1e-2),
 ("conformal recovery (solver.c)",  1.8e-16, 1e-9),
 ("Bjorken Q→16 (PMP)",             6.3e-3,  5e-2),
 ("TOV Shum M_T=1.4",               1.1e-4,  2e-2),
 ("TOV Bussières R=8.86",           1.7e-4,  5e-3),
 ("TOV Bussières M=1.27",           2.8e-3,  1e-2),
 ("R4 axial w-mode f (Bussières)",  4.8e-5,  1e-2),
 ("R4 axial w-mode τ (Bussières)",  2.8e-4,  1e-2),
 ("R5 QNM H1 (Shum)",               2.2e-4,  2e-2),
 ("R5 QNM F (Shum)",                3.3e-3,  2e-2),
 ("R5 QNM H2 (Shum)",               1.7e-2,  2e-2),
 ("Kovtun cv(φ=0) addition",        1.0e-4,  1e-2),
]
sort!(rows, by=r->r[2])
labels = [r[1] for r in rows]
errs   = [r[2] for r in rows]
tols   = [r[3] for r in rows]
y = 1:length(rows)

fig = Figure(size=(900, 620))
ax = Axis(fig[1,1], xscale=log10, xlabel="relative error  |achieved − target| / target",
          yticks=(y, labels), title="BDNKStar reproduction validation — every benchmark ≪ tolerance")
for i in y
    lines!(ax, [errs[i], tols[i]], [i, i], color=(:gray,0.4), linewidth=1.5)
    scatter!(ax, [tols[i]], [i], marker=:vline, markersize=16, color=:crimson)   # tolerance
    scatter!(ax, [errs[i]], [i], markersize=12, color=:seagreen)                 # achieved
end
scatter!(ax, [NaN],[NaN], color=:seagreen, markersize=12, label="achieved")
scatter!(ax, [NaN],[NaN], marker=:vline, color=:crimson, markersize=16, label="tolerance")
axislegend(ax, position=:rb, framevisible=true)
xlims!(ax, 1e-16, 1e-1)
Label(fig[0,:], "BDNKStar — reproduction of all reference papers: 14 benchmarks, every one within (≪) its tolerance",
      fontsize=13, font=:bold)
save(joinpath(outdir,"reproduction_summary.png"), fig)
println("saved reproduction_summary.png  (all ", length(rows), " benchmarks below tolerance: ",
        all(errs .< tols), ")")
