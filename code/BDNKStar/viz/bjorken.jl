#=
    Bjorken flow figure (PMP 2209.09265 Sec. Bjorken_flow, inviscid baseline):
      A  ε(τ) RK4 vs analytic m n0 τ^{-1}[1+e0 τ^{-(Γ-1)}] (log–log)
      B  RK4 self-convergence: max error vs N (∝ N⁻⁴, Q→16)

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bjorken.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
Γ, m, n0, e0 = 4/3, 1.0, 1.0, 1.0
ε0 = bjorken_inviscid_analytic(1.0; Γ=Γ, m=m, n0=n0, e0=e0)

τs, εs = bjorken_evolve_rk4(1.0, 20.0, ε0; N=400, Γ=Γ, m=m, n0=n0)
εa = [bjorken_inviscid_analytic(τ; Γ=Γ, m=m, n0=n0, e0=e0) for τ in τs]

Ns = [125, 250, 500, 1000, 2000, 4000]
errs = Float64[]
for N in Ns
    t, e = bjorken_evolve_rk4(1.0, 20.0, ε0; N=N, Γ=Γ, m=m, n0=n0)
    ea = [bjorken_inviscid_analytic(τ; Γ=Γ, m=m, n0=n0, e0=e0) for τ in t]
    push!(errs, maximum(abs.(e .- ea)))
end

fig = Figure(size=(1150, 460))
axA = Axis(fig[1,1], title="A. Bjorken ε(τ): RK4 vs analytic (Γ=4/3)",
           xlabel="τ", ylabel="ε", xscale=log10, yscale=log10)
lines!(axA, τs, εa, color=:black, linewidth=3, label="analytic")
scatter!(axA, τs[1:20:end], εs[1:20:end], color=:crimson, markersize=8, label="RK4")
axislegend(axA, position=:rt, framevisible=true, labelsize=11)

axB = Axis(fig[1,2], title="B. RK4 self-convergence (4th order, Q→16)",
           xlabel="N steps", ylabel="max |ε_RK4 − ε_analytic|", xscale=log10, yscale=log10)
scatterlines!(axB, Float64.(Ns), errs, color=:purple, linewidth=2.5, markersize=11)
ref = errs[1] .* (Float64.(Ns) ./ Ns[1]).^(-4)
lines!(axB, Float64.(Ns), ref, color=:gray, linestyle=:dash, label="∝ N⁻⁴ (Q=16)")
axislegend(axB, position=:rt, framevisible=true, labelsize=11)

Label(fig[0, :], "BDNKStar Bjorken flow (PMP 2209.09265) — RK4 reproduces analytic ε(τ); 4th-order Q→16",
      fontsize=14, font=:bold)
outfile = joinpath(outdir, "bjorken.png")
save(outfile, fig); println("saved ", outfile)
