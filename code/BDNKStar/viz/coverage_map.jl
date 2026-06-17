#=
    BDNKStar reproduction COVERAGE MAP — per reference paper, the count of original
    figures reproduced vs blocked (by reason).  Summarizes the comprehensive scope:
    5 of 6 papers fully covered for their feasible figures; the blocked ones are
    2D (no 2D reference code), fine-Δr (origin instability), config-halt, or the
    3D GRMHD merger (GPU).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/coverage_map.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# (paper, reproduced, 2D-blocked, fineΔr-blocked, config-halt, GPU-merger)
papers = [
 ("Kovtun 1907.08191",        7, 0, 0, 0, 0),
 ("Pandya 2201.12317",        7, 4, 0, 0, 0),
 ("PMP 2209.09265",           7, 0, 0, 0, 0),
 ("Shum 2509.15303",          4, 0, 2, 0, 0),
 ("Bussières 2604.13208",     4, 0, 0, 1, 0),
 ("Chabanov-Rezzolla 2311…",  0, 0, 0, 0, 26),
]
labels=[p[1] for p in papers]; y=1:length(papers)
rep =[p[2] for p in papers]; b2d=[p[3] for p in papers]
bfd =[p[4] for p in papers]; bcf=[p[5] for p in papers]; bgpu=[p[6] for p in papers]

fig=Figure(size=(960,520))
ax=Axis(fig[1,1], xlabel="# original figures", yticks=(y,labels),
        title="BDNKStar reproduction coverage — reproduced vs blocked (by reason)")
for i in y
    x0=0.0
    for (val,col,lab) in [(rep[i],:seagreen,"reproduced"),(b2d[i],:orange,"2D: no ref code"),
                          (bfd[i],:crimson,"fine-Δr instability"),(bcf[i],:purple,"config halt"),
                          (bgpu[i],:gray60,"3D merger: GPU")]
        val>0 && (barplot!(ax,[i],[val],direction=:x,offset=[x0],color=col); x0+=val)
    end
end
# legend proxies
for (col,lab) in [(:seagreen,"reproduced"),(:orange,"2D: no 2D ref code"),(:crimson,"fine-Δr instability"),(:purple,"config halt"),(:gray60,"3D merger: GPU")]
    barplot!(ax,[NaN],[1],direction=:x,color=col,label=lab)
end
axislegend(ax, position=:rb, framevisible=true)
xlims!(ax,0,28)
tot_rep=sum(rep); tot_blk=sum(b2d)+sum(bfd)+sum(bcf)+sum(bgpu)
Label(fig[0,:], "BDNKStar — $(tot_rep) original figures reproduced + verified; $(tot_blk) blocked (26 = infeasible 3D merger; rest = 2D-no-ref / fine-Δr / config)",
      fontsize=11, font=:bold)
save(joinpath(outdir,"coverage_map.png"), fig)
println("saved coverage_map.png | reproduced=$tot_rep blocked=$tot_blk (merger=", sum(bgpu), ")")
