#=
    Bussières 2604.13208 — axial ℓ=2 w-mode (f,τ) VISCOUS TRAJECTORY, frame A1, as
    η_c sweeps 3e29→1e31 g/cm/s (same kind of trace as complex_plane_2, for the
    frame validated against Table II).  Points coloured by η_c; the curve passes
    through the two published Table-II anchors (×) and starts at the inviscid star.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/axial_qnm_trajectory.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

rows=[split(l) for l in readlines(joinpath(@__DIR__,"..","repro","axial_qnm_trajectory.txt")) if !startswith(l,"#") && !isempty(strip(l))]
ηc=[parse(Float64,r[1]) for r in rows]; f=[parse(Float64,r[2]) for r in rows]; τ=[parse(Float64,r[3]) for r in rows]
visc = ηc .> 0
# Table-II anchors
anchf=[10.4884,10.0898]; ancht=[29.5870,30.8857]; anchη=[3e29,1e31]

fig=Figure(size=(760,540))
ax=Axis(fig[1,1], xlabel="f  [kHz]", ylabel="τ  [µs]",
        title="Bussières axial ℓ=2 w-mode (f,τ) viscous trajectory — frame A1, η_c sweep")
lines!(ax, f[visc], τ[visc], color=:gray60, linewidth=1.5)
sc=scatter!(ax, f[visc], τ[visc], color=log10.(ηc[visc]), colormap=:viridis, markersize=13)
scatter!(ax, [f[1]], [τ[1]], color=:black, marker=:star5, markersize=20, label="inviscid (η_c=0)")
scatter!(ax, anchf, ancht, color=:red, marker=:xcross, markersize=18, label="Bussières Table II")
Colorbar(fig[1,2], sc, label="log₁₀ η_c [g cm⁻¹ s⁻¹]")
axislegend(ax, position=:lt, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce axial w-mode viscous (f,τ) trajectory (frame A1): curve passes through Table-II anchors (<0.04%) [complex_plane_2 physics, validated frame]",
      fontsize=10, font=:bold)
save(joinpath(outdir,"axial_qnm_trajectory.png"), fig)
@printf("saved axial_qnm_trajectory.png | %d viscous pts; anchors: 3e29→(%.4f,%.4f) vs (10.4884,29.587); 1e31→(%.4f,%.4f) vs (10.0898,30.886)\n",
        count(visc), f[findmin(abs.(ηc.-3e29))[2]], τ[findmin(abs.(ηc.-3e29))[2]],
        f[end], τ[end])
