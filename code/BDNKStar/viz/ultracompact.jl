#=
    Bussières 2604.13208 Sec.V-D ultracompact trapped ℓ=2 w-modes (from the
    verified solver repro/axial_ultracompact.jl). Complex-ω plane: as the
    compactness 𝒞=M/R rises (0.40→0.44, constant-density frame-B stars), the
    trapped-mode ladder deepens (2→6 long-lived modes) and the fundamental
    becomes longer-lived (|ω_I| drops); shear viscosity increases the damping.
    [Complementary to the paper's ω(ℓ) view in plot_ultracompact.]

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/ultracompact.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# verified solver output (repro/axial_ultracompact.jl)
ideal = Dict(
 0.40 => [(0.103087,0.005498),(0.136518,0.022529)],
 0.42 => [(0.087979,0.013026),(0.128226,0.021146)],
 0.44 => [(0.036606,0.000670),(0.064275,0.003981),(0.078842,0.005427),
          (0.093567,0.006610),(0.108376,0.007580),(0.123246,0.008383)])
visc  = Dict(0.40=>(0.101072,0.009176), 0.42=>(0.088796,0.016429), 0.44=>(0.036659,0.001844))
cols  = Dict(0.40=>:crimson, 0.42=>:seagreen, 0.44=>:dodgerblue)

fig = Figure(size=(760, 470))
ax = Axis(fig[1,1], xlabel="ω_R M  (geom.)", ylabel="−ω_I M  (damping)", yscale=log10,
          title="Ultracompact ℓ=2 trapped w-mode ladder vs compactness 𝒞 (Bussières Sec.V-D)")
for C in (0.40, 0.42, 0.44)
    pts = ideal[C]
    scatter!(ax, [p[1] for p in pts], [p[2] for p in pts], marker=:rect, markersize=13,
             color=cols[C], label="𝒞=$C inviscid ($(length(pts)) trapped)")
    v = visc[C]
    scatter!(ax, [v[1]], [v[2]], marker=:circle, markersize=15, color=cols[C],
             strokecolor=:black, strokewidth=1.5)
end
# guide the eye along each ladder
for C in (0.44,)
    pts = ideal[C]
    lines!(ax, [p[1] for p in pts], [p[2] for p in pts], color=(cols[C],0.4), linewidth=1.2)
end
axislegend(ax, position=:rb, framevisible=true, labelsize=11)
Label(fig[0,:], "BDNKStar — ultracompact trapped w-modes: ladder deepens (2→6) + longer-lived as 𝒞 rises; ○ = viscous (more damped)",
      fontsize=12.5, font=:bold)
outfile = joinpath(outdir, "ultracompact_reproduction.png")
save(outfile, fig); println("saved ", outfile)
