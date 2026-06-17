#=
    Bussières 2604.13208 Sec.V-C — axial ℓ=2 avoided crossing (frame B1, EOS1):
    Re ω vs η_c for the w-mode (continued from the inviscid mode, blue ●) and the
    η-mode (viscosity-driven second-sound resonance with NO inviscid counterpart,
    red ▲).  The η-mode is located by a complex scan at large η_c and down-tracked
    until the AVOIDED CROSSING terminates the track (~3.96e31) — exactly the
    obstruction the paper reports [main.tex 588,591].  The w-mode frequency bends
    (repulsion) near the same η_c — the avoided-crossing signature.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/axial_avoided_crossing.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

W=Tuple{Float64,Float64}[]; E=Tuple{Float64,Float64}[]
for l in readlines(joinpath(@__DIR__,"..","repro","axial_avoided_crossing.txt"))
    startswith(l,"#") && continue; isempty(strip(l)) && continue
    p=split(l)
    p[1]=="W" && push!(W,(parse(Float64,p[2]),parse(Float64,p[3])))
    p[1]=="E" && push!(E,(parse(Float64,p[2]),parse(Float64,p[3])))
end
# w-mode: drop the η_c=0 inviscid point for the log-x plot (mark separately)
w0 = W[1][2]; Wv=[x for x in W if x[1]>0]

fig=Figure(size=(780,520))
ax=Axis(fig[1,1], xscale=log10, xlabel="η_c  [g cm⁻¹ s⁻¹]", ylabel="Re ω  [1/km]",
        title="Bussières axial ℓ=2 avoided crossing (frame B1): w-mode (●) vs η-mode (▲)")
scatterlines!(ax, [x[1] for x in Wv], [x[2] for x in Wv], color=:dodgerblue, markersize=11, linewidth=2.0, label="w-mode (← inviscid)")
hlines!(ax, [w0], color=:dodgerblue, linestyle=:dot, linewidth=1.0)
scatterlines!(ax, [x[1] for x in E], [x[2] for x in E], color=:crimson, marker=:utriangle, markersize=12, linewidth=2.0, label="η-mode (no inviscid counterpart)")
xc = E[end][1]
vlines!(ax, [xc], color=:gray, linestyle=:dash, linewidth=1.2)
text!(ax, xc*0.7, 0.16, text="avoided crossing\n(η-track lost ≈$(round(xc/1e31,digits=2))e31)\n[paper too, 588]", fontsize=10, color=:gray25)
axislegend(ax, position=:lc, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Bussières η-mode + avoided crossing: η-mode located (complex scan), down-tracked until the avoided crossing terminates it [main.tex 588,591]",
      fontsize=9.5, font=:bold)
save(joinpath(outdir,"axial_avoided_crossing.png"), fig)
@printf("saved axial_avoided_crossing.png | w-mode %d pts (Re ω %.3f→bend→%.3f); η-mode %d pts (%.3f→%.3f), avoided crossing ≈%.2e\n",
        length(Wv), Wv[1][2], Wv[end][2], length(E), E[1][2], E[end][2], xc)
