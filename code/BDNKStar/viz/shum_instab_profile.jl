#=
    Shum fine-Δr instability — SPATIAL localization.  |δε(r)|=|ε−ε_bg| just before
    blow-up: the fine Δr=0.0024 run's perturbation is a huge spike at the ORIGIN
    (r≈0), while the stable coarse Δr=0.04 run carries only the physical ~1e-6
    surface QNM.  ⇒ the instability is a CENTRAL r=0 regularity problem (the 1/r
    terms / origin boundary), NOT the surface — refines the fix hypothesis.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_instab_profile.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

function parse_blocks(path)
    rF=Float64[]; dF=Float64[]; rC=Float64[]; dC=Float64[]; sect=0
    for l in readlines(path)
        if startswith(l,"# FINE"); sect=1; continue; end
        if startswith(l,"# COARSE"); sect=2; continue; end
        (isempty(strip(l)) || startswith(l,"#")) && continue
        p=split(l); r=parse(Float64,p[1]); de=parse(Float64,p[2])
        if sect==1; push!(rF,r); push!(dF,de); elseif sect==2; push!(rC,r); push!(dC,de); end
    end
    return rF,dF,rC,dC
end
rF,dF,rC,dC = parse_blocks(joinpath(@__DIR__, "..", "repro", "shum_instab_profile.txt"))
af(v)=max.(abs.(v),1e-18)

fig=Figure(size=(820,500))
ax=Axis(fig[1,1], yscale=log10, xlabel="r / M_⊙", ylabel="|δε(r)| = |ε − ε_bg|",
        title="Shum fine-Δr instability is at the ORIGIN (r≈0), not the surface")
lines!(ax, rF, af(dF), color=:crimson,    linewidth=2.0, label="Δr=0.0024 (fine, near blow-up) — spike at r≈0")
lines!(ax, rC, af(dC), color=:black,      linewidth=2.0, label="Δr=0.04 (coarse, stable) — physical ~1e-6 surface QNM")
vlines!(ax, [0.0], color=:gray, linestyle=:dot)
xlims!(ax, -0.3, 10); ylims!(ax, 1e-9, 1e4)
axislegend(ax, position=:rt, framevisible=true)
Label(fig[0,:], "BDNKStar — Shum fine-Δr instability LOCALIZED to the origin (r=0 regularity / 1/r terms) [refines fix hypothesis]",
      fontsize=10.5, font=:bold)
save(joinpath(outdir,"shum_instab_profile.png"), fig)
println("saved shum_instab_profile.png | fine max|δε|=", round(maximum(abs.(dF)),sigdigits=4),
        " at r=", round(rF[argmax(abs.(dF))],digits=3), " ; coarse max|δε|=", round(maximum(abs.(dC)),sigdigits=3),
        " at r=", round(rC[argmax(abs.(dC))],digits=2))
