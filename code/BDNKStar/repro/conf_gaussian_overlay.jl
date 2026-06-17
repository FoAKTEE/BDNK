#=
    Pandya 2201.12317 — EXACT overlay of the Julia conformal-BDNK engine against
    the published reference C code (ref-code/1D_conformal_bdnk, A. Pandya).

    Same setup as the reference parameters.h:
      conformal PMP-luminal frame  eta0 = 10^{1/4}/(3π),  (λ0,χ0)=(25/7,25/4)·eta0
      Gaussian ID  ε = exp(−x²/25) + 0.1,  u^x=0,  x∈[−200,200],  N=129,  CFL=0.1
      outflow (ghost) BC;  dt = CFL·dx = 0.3125;  save every 10 steps.

    We evolve with the Julia engine (WENO5 + Kurganov–Tadmor + Heun-RK2 + BDNK
    primitive solve) and dump ε(x) at the reference snapshot rows {1,16,31,51,102}
    → t = {0, 46.875, 93.75, 156.25, 315.625}.  The viz overlays these on the C
    eps.txt to show the two independent codes agree.

    Run: julia code/BDNKStar/repro/conf_gaussian_overlay.jl
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using .BDNKStar.ConformalBDNK: ConformalFrame

eta0 = 10.0^0.25/(3π)                       # = 0.188681 (matches C exactly)
fr   = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)
s    = init_gaussian(fr; N=129, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)

# reference rows (1-based) -> cumulative C timesteps = 10·(row-1)
rows  = [1, 16, 31, 51, 102]
nstep = [10*(r-1) for r in rows]            # 0,150,300,500,1010
labels= ["t=$(round(n*s.dt,digits=3))" for n in nstep]

x = copy(s.x)
cols = Vector{Vector{Float64}}()
done = 0
push!(cols, copy(energy_density(s)))        # t=0
for k in 2:length(nstep)
    evolve!(s, nstep[k]-done); global done = nstep[k]
    push!(cols, copy(energy_density(s)))
    println("  evolved to step $done (t=$(round(done*s.dt,digits=3)))  εmax=",
            round(maximum(cols[end]),sigdigits=5), " nan=", any(!isfinite, cols[end]))
end

open(joinpath(@__DIR__,"conf_gaussian_overlay.txt"),"w") do io
    println(io, "# x  ", join(labels, "  "), "   (Julia conformal-BDNK engine vs C reference rows ", rows, ")")
    for i in eachindex(x)
        print(io, x[i]); for c in cols; print(io, "  ", c[i]); end; println(io)
    end
end
println("SAVED conf_gaussian_overlay.txt  (",length(x)," x-points, ",length(cols)," times)")
