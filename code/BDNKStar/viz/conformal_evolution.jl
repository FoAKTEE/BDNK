#=
    Conformal BDNK flat-space evolution figure (reproduces the Pandya 2201.12317
    flat-space tests with the solver.c engine):
      A  Gaussian energy clump — initial vs evolved (viscous spreading / two
         outgoing pulses)
      B  steady planar shock — erf initial data relaxing to the BDNK steady
         profile; asymptotic states = Rankine–Hugoniot (εL=1, εR=4.4074)

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/conformal_evolution.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
fr = pmp_luminal_frame(10.0)

# Gaussian
sg = init_gaussian(fr; N=257, A=1.0, x0=0.0, w=25.0, c=0.1, cfl=0.1)
xg = copy(sg.x); εg0 = copy(energy_density(sg))
evolve!(sg, 250); εgT = energy_density(sg)

# Steady shock
ss = init_smooth_shock(fr; N=257, εL=1.0, vL=0.8, cfl=0.1)
xs = copy(ss.x); εs0 = copy(energy_density(ss))
evolve!(ss, 300); εsT = energy_density(ss)

fig = Figure(size=(1150, 460))
axA = Axis(fig[1,1], title="A. Gaussian clump (conformal BDNK)", xlabel="x", ylabel="ε")
lines!(axA, xg, εg0, color=:gray, linewidth=2, linestyle=:dash, label="t = 0")
lines!(axA, xg, εgT, color=:dodgerblue, linewidth=2.5, label="evolved")
xlims!(axA, -200, 200); axislegend(axA, position=:rt, framevisible=true, labelsize=11)

axB = Axis(fig[1,2], title="B. Steady planar shock (εL=1, vL=0.8)", xlabel="x", ylabel="ε")
lines!(axB, xs, εs0, color=:gray, linewidth=2, linestyle=:dash, label="t=0 (erf guess)")
lines!(axB, xs, εsT, color=:crimson, linewidth=2.5, label="relaxed BDNK shock")
hlines!(axB, [1.0, 4.40741], color=:black, linestyle=:dot, label="RH: εL=1, εR=4.4074")
xlims!(axB, -100, 100); axislegend(axB, position=:rt, framevisible=true, labelsize=10)

Label(fig[0, :], "BDNKStar conformal flat-space evolution (WENO5+KT+Heun, solver.c port) — stable; RH-consistent shock",
      fontsize=14, font=:bold)

outfile = joinpath(outdir, "conformal_evolution.png")
save(outfile, fig); println("saved ", outfile)
