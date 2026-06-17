#=
    Pandya 2201.12317 Tab_cons — discrete conservation of ∫T^tt.
    FV (conserved variable, ×) holds at machine precision; FD (reconstructed from
    the non-conservatively evolved primitives, line) drifts; both jump when the
    pulse reaches the outflow boundary.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_tab_cons.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

rows=[split(l) for l in readlines(joinpath(repo,"conf_tab_cons.txt")) if !startswith(l,"#") && !isempty(strip(l))]
t  =[parse(Float64,r[1]) for r in rows]
fv =[max(parse(Float64,r[2]),1e-16) for r in rows]
fd =[max(parse(Float64,r[3]),1e-16) for r in rows]

fig=Figure(size=(820,520))
ax=Axis(fig[1,1], yscale=log10, xlabel="t", ylabel="|∫ [T^tt(t,x) − T^tt(0,x)] dx|",
        title="Pandya Tab_cons — FV conserves ∫T^tt to machine precision, FD drifts")
lines!(ax, t, fd, color=:gray55, linewidth=3.0, label="FD (primitive-reconstructed)")
scatter!(ax, t, fv, color=:black, marker=:xcross, markersize=7, label="FV (conserved variable)")
# boundary marker where both jump
ib = findfirst(i->fv[i]>1e-6, eachindex(fv))
ib !== nothing && vlines!(ax, [t[ib]], color=:red, linestyle=:dot, linewidth=1.2)
ylims!(ax, 1e-16, 5e1)
axislegend(ax, position=:rc, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Pandya Tab_cons: FV ∫T^tt machine-precision conserved (~1e-14) vs FD drift (~1e-2); both jump at outflow boundary",
      fontsize=11, font=:bold)
save(joinpath(outdir,"conf_tab_cons.png"), fig)
@printf("saved conf_tab_cons.png | FV plateau=%.2e  FD plateau=%.2e  boundary-jump t≈%.0f\n",
        fv[20], fd[20], ib===nothing ? -1 : t[ib])
