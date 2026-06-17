#=
    PMP shock_instability v(x) profiles (the actual evolution snapshots) — from
    the detached viscous-BDNK runs (repro/shock_profiles_run.jl). Top: τ̂=3
    (c₊<v_max ⇒ high-freq instability where v>c₊) at t=27. Bottom: τ̂=1.5
    (c₊>v_max ⇒ stable steady shock) at t=100. Dotted line = c₊.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pmp_shock_profiles.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

function load(f)
    rows = [split(l) for l in readlines(f) if !startswith(l,"#") && !isempty(strip(l))]
    x  = [parse(Float64,r[1]) for r in rows]
    v  = [parse(Float64,r[2]) for r in rows]
    cp = [parse(Float64,r[3]) for r in rows]
    return x, v, cp
end
x3,v3,cp3 = load(joinpath(repo,"shock_prof_t3.txt"))
x1,v1,cp1 = load(joinpath(repo,"shock_prof_t1p5.txt"))

fig = Figure(size=(720, 560))
ax1 = Axis(fig[1,1], ylabel="v", title="τ̂ = 3  (c₊ < v_max ⇒ unstable)  t=27")
lines!(ax1, x3, v3, color=:black, linewidth=1.6)
lines!(ax1, x3, cp3, color=:black, linestyle=:dot, linewidth=1.6)
xlims!(ax1,-50,50); ylims!(ax1,0.3,0.95); text!(ax1, 26, 0.80, text="⋯ c₊", fontsize=12)
ax2 = Axis(fig[2,1], xlabel="x", ylabel="v", title="τ̂ = 1.5  (c₊ > v_max ⇒ stable)  t=100")
lines!(ax2, x1, v1, color=:black, linewidth=1.6)
lines!(ax2, x1, cp1, color=:black, linestyle=:dot, linewidth=1.6)
xlims!(ax2,-50,50); ylims!(ax2,0.3,1.0); text!(ax2, 26, 0.92, text="⋯ c₊", fontsize=12)
Label(fig[0,:], "BDNKStar — reproduce PMP shock_instability v(x) profiles (viscous BDNK evolution)",
      fontsize=13, font=:bold)
save(joinpath(outdir,"pmp_shock_profiles.png"), fig)
println("saved  τ̂=3 vmax=",round(maximum(v3),digits=3)," c₊min=",round(minimum(cp3),digits=3),
        "  | τ̂=1.5 vmax=",round(maximum(v1),digits=3)," c₊min=",round(minimum(cp1),digits=3))
