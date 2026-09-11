#=
    BDNK viscoresistive relativistic MHD — 1D SHOCK-TUBE figure (Stage 1b).
    Reproduces the perpendicular-shock scan of Lier, Armas, Porth (2026,
    arXiv:2606.22691), Eq.22 / Fig.3, Table I (ST-a..j).

    Panels:
      (A) the out-of-plane magnetic component J^{ty}(x, t=0.4) across the scan,
      (B) the pressure p=ε/3 (x, t=0.4),
      (C) the r_b trend  : at fixed D_u=D_ε, raising r_b smooths J^{ty},
      (D) the D_u trend  : at fixed r_b, raising D_u,D_ε smooths the pressure
          fast/rarefaction waves;
      and a causality strip: the per-cell front-velocity v_max field.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bdnk_mhd_shocktube.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# evolve a case to t=0.4 at the paper resolution and return the profiles
function run_case(tag; N=1024, tmax=0.4)
    eng, st = setup_mhd1d_shocktube(tag; N=N)
    d = evolve_mhd1d!(st, eng; tmax=tmax, monitor_causality=false)
    pr = mhd1d_primitives(st, eng)
    vf = mhd1d_vmax_field(st, eng; stride=24, nθ=61)
    return pr, vf, d
end

fig = Figure(size=(1500, 980))

# ── Panel A: J^{ty} across the full scan ─────────────────────────────────────
axA = Axis(fig[1,1], xlabel="x", ylabel="J^{ty}",
           title="(A) out-of-plane B  J^{ty}(x, t=0.4) — scan ST-a..j (Eq.22)")
cols = cgrad(:viridis, 10; categorical=true)
results = Dict{Symbol,Any}()
for (k, tag) in enumerate((:a,:b,:c,:d,:e,:f,:g,:h,:i,:j))
    pr, vf, d = run_case(tag)
    results[tag] = (pr=pr, vf=vf, d=d)
    lines!(axA, pr.x, pr.Jty; color=cols[k], linewidth=1.6, label="ST-$tag")
    @printf("ST-%s: nan=%s εmax=%.3f |Jty|max=%.3f TV(Jty)=%.3f\n",
            tag, d.nan, maximum(pr.ε), maximum(abs.(pr.Jty)),
            sum(abs.(diff(pr.Jty))))
end
axislegend(axA; position=:lb, nbanks=2, labelsize=8, framevisible=true)

# ── Panel B: pressure across the full scan ───────────────────────────────────
axB = Axis(fig[1,2], xlabel="x", ylabel="p = ε/3",
           title="(B) pressure p(x, t=0.4) — scan ST-a..j")
for (k, tag) in enumerate((:a,:b,:c,:d,:e,:f,:g,:h,:i,:j))
    lines!(axB, results[tag].pr.x, results[tag].pr.p; color=cols[k], linewidth=1.6)
end
hlines!(axB, [1.0, 0.1]; color=:gray, linestyle=:dot, linewidth=1.0)

# ── Panel C: the r_b trend (smooths J^{ty}) ──────────────────────────────────
axC = Axis(fig[2,1], xlabel="x", ylabel="J^{ty}",
           title="(C) r_b TREND at D_u=D_ε=1e-3 :\nr_b=1e-4 (d) → 1e-3 (e) → 1e-2 (f) smooths the magnetic profile")
for (tag, col, lab) in ((:d, :navy, "r_b=1e-4"), (:e, :seagreen, "r_b=1e-3"),
                        (:f, :crimson, "r_b=1e-2"))
    lines!(axC, results[tag].pr.x, results[tag].pr.Jty; color=col, linewidth=2.2, label=lab)
end
axislegend(axC; position=:lb, framevisible=true, labelsize=9)

# ── Panel D: the D_u,D_ε trend (smooths pressure waves) ──────────────────────
axD = Axis(fig[2,2], xlabel="x", ylabel="p = ε/3",
           title="(D) D_u=D_ε TREND at r_b=1e-2 :\nD=1e-4 (c) → 1e-3 (f) → 1e-2 (i) smooths the fast/rarefaction waves")
for (tag, col, lab) in ((:c, :navy, "D=1e-4"), (:f, :seagreen, "D=1e-3"),
                        (:i, :crimson, "D=1e-2"))
    lines!(axD, results[tag].pr.x, results[tag].pr.p; color=col, linewidth=2.2, label=lab)
end
axislegend(axD; position=:rt, framevisible=true, labelsize=9)

# ── Panel E: ST-i vs ST-j (the extra-τ_X case) overlay ───────────────────────
axE = Axis(fig[3,1], xlabel="x", ylabel="J^{ty}",
           title="(E) ST-i (τ_X=2e-2) vs ST-j (τ_X=4e-2, extra-τ_X-for-stability)")
lines!(axE, results[:i].pr.x, results[:i].pr.Jty; color=:navy, linewidth=2.2, label="ST-i τ_X=2e-2")
lines!(axE, results[:j].pr.x, results[:j].pr.Jty; color=:darkorange, linewidth=2.2,
       linestyle=:dash, label="ST-j τ_X=4e-2")
axislegend(axE; position=:lb, framevisible=true, labelsize=9)

# ── Panel F: causality strip — front-velocity v_max field ────────────────────
axF = Axis(fig[3,2], xlabel="x", ylabel="front velocity v_max(x)",
           title="(F) causality monitor (foundation v_max field, ST-e)")
vf = results[:e].vf
lines!(axF, vf.x, vf.vmax; color=:purple, linewidth=2.2, label="v_max")
hlines!(axF, [1.0]; color=:black, linestyle=:dash, linewidth=1.2, label="luminal")
band!(axF, vf.x, fill(0.0, length(vf.x)), vf.imW; color=(:red,0.2))
ylims!(axF, 0.0, max(1.1, maximum(vf.vmax)*1.05))
axislegend(axF; position=:rt, framevisible=true, labelsize=9)

Label(fig[0, :],
      "BDNKStar — reproduce Lier–Armas–Porth (2026): BDNK viscoresistive MHD 1D shock-tube scan (Eq.22 / Fig.3 / Table I)",
      fontsize=14, font=:bold)

outfile = joinpath(outdir, "bdnk_mhd_shocktube.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)
