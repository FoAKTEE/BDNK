#=
    Shum fine-Δr instability diagnostic — documents why error_fit / decay-continuum
    / stable_evol_resolutions(inset) are out of reach.  The nonlinear Cowling BDNK
    engine is STABLE at coarse Δr (0.04, used for the matched QNM + per-Δr decay)
    but blows up early at the paper's fine Δr (≲0.0032), at t≈37/27/4 for
    Δr=0.0032/0.0024/0.0020 — and elevated Kreiss–Oliger dissipation (σKO=2–6)
    only delays it.  The fine-Δr ladder needs t≳1500 for a decay fit, so those
    three figures require a deeper stability fix (not KO) — documented HOLE.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_fine_instability.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

load(f)=( rows=[split(l) for l in readlines(f) if !startswith(l,"#") && !isempty(strip(l))];
          ([parse(Float64,r[1]) for r in rows], [parse(Float64,r[2]) for r in rows]) )
tc,ec = load("/tmp/r5_coarse_0.04_head.txt")
t32,e32 = load("/tmp/r5_fine_0.0032.txt")
t24,e24 = load("/tmp/r5_fine_0.0024.txt")
t20,e20 = load("/tmp/r5_fine_0.002.txt")

fig=Figure(size=(820,500))
ax=Axis(fig[1,1], xlabel="t  [M_⊙]", ylabel="ε_c(t)",
        title="Shum fine-Δr instability — engine blows up at Δr≲0.0032 (σKO=2; KO only delays it)")
lines!(ax, tc,  ec,  color=:black,     linewidth=2.2, label="Δr=0.04 (stable — used for QNM+decay)")
lines!(ax, t32, e32, color=:dodgerblue,linewidth=1.8, label="Δr=0.0032 → blow-up t≈37")
lines!(ax, t24, e24, color=:seagreen,  linewidth=1.8, label="Δr=0.0024 → blow-up t≈27")
scatter!(ax, t20, e20, color=:crimson, markersize=9, label="Δr=0.0020 → blow-up t≈4")
xlims!(ax, 0, 50); ylims!(ax, 0.0, 0.0027)
axislegend(ax, position=:rt, framevisible=true)
Label(fig[0,:], "BDNKStar — Shum fine-Δr instability: blocks error_fit / decay-continuum / stable_evol_resolutions(inset) [documented HOLE; needs deeper stability fix]",
      fontsize=10.5, font=:bold)
save(joinpath(outdir,"shum_fine_instability.png"), fig)
println("saved shum_fine_instability.png | blow-up rows: 0.0032=",length(t32)," 0.0024=",length(t24)," 0.002=",length(t20))
