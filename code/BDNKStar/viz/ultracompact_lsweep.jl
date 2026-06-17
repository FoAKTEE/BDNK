#=
    Bussières 2604.13208 plot_ultracompact — trapped axial w-mode spectrum ω(ℓ)
    for ultracompact compactnesses R/M ∈ {2.40,2.45,2.50,2.60}.  Top: ω_R·M vs ℓ;
    bottom: −ω_I·M vs ℓ (log).  My converged ideal modes (lines) vs the original's
    digitized values (×).  R/M=2.60 matches BOTH panels across ℓ=2..6; the trapped
    modes get very narrow (hard to hunt) at higher compactness/ℓ — the paper itself
    plots fewer points there.  PARTIAL reproduction (ideal; viscous continuation
    unreliable for the narrow deep-well modes).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/ultracompact_lsweep.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# my sweep
rows=[split(l) for l in readlines(joinpath(@__DIR__,"..","repro","ultracompact_lsweep.txt")) if !startswith(l,"#") && !isempty(strip(l))]
data=Dict{Float64,Vector{Tuple{Int,Float64,Float64}}}()
for r in rows
    rm=parse(Float64,r[1]); ℓ=parse(Int,r[2]); wR=parse(Float64,r[3]); wI=parse(Float64,r[4])
    push!(get!(data,rm,[]), (ℓ,wR,-wI))
end
# original digitized reference (from plot_ultracompact-1.png)
origR=Dict(2.60=>[(2,0.43),(3,0.64),(4,0.83),(5,1.01),(6,1.20)],
           2.50=>[(2,0.41),(3,0.61),(4,0.77),(5,0.94),(6,1.10)],
           2.45=>[(2,0.40),(3,0.57),(4,0.72),(5,0.88)],
           2.40=>[(2,0.36),(3,0.51),(4,0.66)])
origI=Dict(2.60=>[(2,4.0e-2),(3,2.5e-2),(4,1.7e-2),(5,1.0e-2),(6,6.5e-3)],
           2.50=>[(2,2.0e-2),(3,8e-3),(4,3e-3),(5,1e-3),(6,3.3e-4)],
           2.45=>[(2,1.2e-2),(3,2.5e-3),(4,5.5e-4),(5,1.0e-4)],
           2.40=>[(2,4.5e-3),(3,4.5e-4),(4,4e-5)])
cols=Dict(2.40=>:crimson,2.45=>:dodgerblue,2.50=>:seagreen,2.60=>:darkorange)

# keep my mode only where it converged to the trapped fundamental: matches the
# original ω_R reference within 8% (so the figure shows where I reproduce vs.
# where the narrow deep-well modes resist automated hunting).
function branch(v, ref)
    rd=Dict(ref); out=Tuple{Int,Float64,Float64}[]
    for (ℓ,wR,nwI) in sort(v,by=x->x[1])
        haskey(rd,ℓ) || continue
        abs(wR-rd[ℓ])/rd[ℓ] < 0.08 || continue
        push!(out,(ℓ,wR,nwI))
    end
    out
end

fig=Figure(size=(720,720))
ax1=Axis(fig[1,1], ylabel="ω_R M", xticklabelsvisible=false, title="trapped axial w-modes ω(ℓ) — ideal (mine: lines, paper: ×)")
ax2=Axis(fig[2,1], yscale=log10, xlabel="ℓ", ylabel="−ω_I M")
for rm in (2.40,2.45,2.50,2.60)
    b=branch(data[rm], origR[rm]); c=cols[rm]
    if !isempty(b)
        scatterlines!(ax1, [x[1] for x in b], [x[2] for x in b], color=c, linewidth=2.0, markersize=9, label="R=$(rm)M")
        scatterlines!(ax2, [x[1] for x in b], [max(x[3],1e-6) for x in b], color=c, linewidth=2.0, markersize=9)
    end
    scatter!(ax1, [x[1] for x in origR[rm]], [x[2] for x in origR[rm]], color=c, marker=:xcross, markersize=13)
    scatter!(ax2, [x[1] for x in origI[rm]], [x[2] for x in origI[rm]], color=c, marker=:xcross, markersize=13)
end
xlims!(ax1,1.7,6.3); xlims!(ax2,1.7,6.3); ylims!(ax2,1e-4,2e-1)
axislegend(ax1, position=:lt, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Bussières plot_ultracompact (ideal): R=2.60M full match both panels; narrow deep-well modes harder at higher 𝒞 [PARTIAL]",
      fontsize=10, font=:bold)
save(joinpath(outdir,"ultracompact_lsweep.png"), fig)
for rm in (2.40,2.45,2.50,2.60); println("R/M=$rm matched: ", [x[1] for x in branch(data[rm],origR[rm])]); end
println("saved ultracompact_lsweep.png")
