#=
    Pandya 2201.12317 Conv_plot — self-convergence of the conformal-BDNK Gaussian
    clump.  Using the finest grid (N=513) as the reference, the dx-weighted L1
    error of N=129 and N=257 vs the reference gives the convergence order for BOTH
    the Julia engine and the reference C code; they should match.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conf_convergence.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

# load Julia ladder (x,eps columns)
function loadjl(N)
    rows=[split(l) for l in readlines(joinpath(repo,"conf_conv_N$(N).txt")) if !startswith(l,"#") && !isempty(strip(l))]
    [parse(Float64,r[2]) for r in rows]
end
# load C last-rows (single row of eps)
loadc(R) = [parse(Float64,t) for t in split(strip(read(joinpath(repo,"conf_cref_conv_RES$(R).txt"),String)))]

u = Dict(:jl=>Dict(129=>loadjl(129),257=>loadjl(257),513=>loadjl(513)),
         :c =>Dict(129=>loadc(1),    257=>loadc(2),    513=>loadc(4)))

# dx-weighted L1 error vs the N=513 reference (restricted to the coarse grid)
function errN(d, N)
    uN = d[N]; uR = d[513]; step = 512 ÷ (N-1)
    ref = uR[1:step:end]                       # restrict 513-grid to N-grid
    dx = 400.0/(N-1)
    sum(abs.(uN .- ref)) * dx
end
res = Dict()
for code in (:jl,:c)
    e1 = errN(u[code],129); e2 = errN(u[code],257)
    p  = log2(e1/e2)
    res[code] = (e1=e1, e2=e2, p=p)
    @printf("  %-5s: err(N=129)=%.4e  err(N=257)=%.4e  order p=%.2f\n", code, e1, e2, p)
end

Ns = [129, 257]
fig = Figure(size=(820, 540))
ax = Axis(fig[1,1], xscale=log10, yscale=log10, xlabel="N (cells)",
          ylabel="dx-weighted L1 error vs N=513 reference",
          title="Pandya Conv_plot — self-convergence (Julia engine vs reference C code)")
scatterlines!(ax, Ns, [res[:jl].e1, res[:jl].e2], color=:dodgerblue, markersize=12,
              linewidth=2.4, label="Julia engine  (p=$(round(res[:jl].p,digits=2)))")
scatterlines!(ax, Ns, [res[:c].e1, res[:c].e2], color=:crimson, marker=:rect, markersize=11,
              linewidth=2.0, linestyle=:dash, label="reference C code  (p=$(round(res[:c].p,digits=2)))")
# reference slopes p=1 and p=2 anchored at the N=129 Julia error
e0 = res[:jl].e1
lines!(ax, Ns, [e0, e0*(129/257)^1], color=(:gray,0.6), linewidth=1.2, linestyle=:dot, label="slope p=1")
lines!(ax, Ns, [e0, e0*(129/257)^2], color=(:gray,0.35), linewidth=1.2, linestyle=:dashdot, label="slope p=2")
axislegend(ax, position=:rt, framevisible=true)
Label(fig[0,:], "BDNKStar — reproduce Pandya Conv_plot: both codes converge at order p≈$(round(res[:jl].p,digits=1)) (Julia=$(round(res[:jl].p,digits=2)), C=$(round(res[:c].p,digits=2)))",
      fontsize=12, font=:bold)
save(joinpath(outdir,"conf_convergence.png"), fig)
println("saved conf_convergence.png | Julia p=", round(res[:jl].p,digits=2), " C p=", round(res[:c].p,digits=2))
