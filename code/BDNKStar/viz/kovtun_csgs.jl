#=
    Kovtun 1907.08191 piccsGs01 + piccsGs09 (Fig. soundspeeds): the two sound
    branches' velocity c_s(φ) (solid) and damping Γ_s(φ) (dashed) vs the angle φ
    between k and v₀, for v_s=1/2 at v₀=0.1 and v₀=0.9. c_s(φ) is the cvphi
    formula (kovtun_cv, c₀=v_s); Γ_s from eq below eq:csv-eqn (in units γs/w₀).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/kovtun_csgs.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# Γ_s(c_s, φ; v0, vs)  (paper eq after eq:csv-eqn), in units γs/w₀
function Γs(cs, φ, v0, vs)
    num = cs - v0*cos(φ)
    f2  = 1 + cs^2*v0^2 - 2cs*v0*cos(φ) - v0^2*sin(φ)^2
    den = cs*(1 - v0^2*vs^2) - v0*(1 - vs^2)*cos(φ)
    return (num/sqrt(1-v0^2)) * (f2/den)
end

vs = 0.5
φs = range(0, π; length=400)
fig = Figure(size=(1050, 430))
for (k, v0) in enumerate((0.1, 0.9))
    ax = Axis(fig[1,k], xlabel="φ", ylabel = k==1 ? "c_s,  (w₀/γs)Γ_s" : "",
              title="piccsGs$(v0==0.1 ? "01" : "09")  (v_s=0.5, v₀=$v0)")
    hlines!(ax, [0.0], color=:black, linewidth=0.8)
    cp = [kovtun_cv(v0, φ; c0=vs)[1] for φ in φs]   # branch +
    cm = [kovtun_cv(v0, φ; c0=vs)[2] for φ in φs]   # branch −
    lines!(ax, φs, cp, color=:dodgerblue, linewidth=2.2)               # c_s solid
    lines!(ax, φs, cm, color=:orange,     linewidth=2.2)
    lines!(ax, φs, [Γs(cp[i], φs[i], v0, vs) for i in eachindex(φs)], color=:dodgerblue, linestyle=:dash, linewidth=2)   # Γ_s dashed
    lines!(ax, φs, [Γs(cm[i], φs[i], v0, vs) for i in eachindex(φs)], color=:orange,     linestyle=:dash, linewidth=2)
    ylims!(ax, -1.1, 1.55)
end
Label(fig[0,:], "BDNKStar — reproduce Kovtun piccsGs: sound velocity c_s (solid) + damping Γ_s (dashed) vs φ",
      fontsize=13, font=:bold)
save(joinpath(outdir, "kovtun_csgs.png"), fig)
println("saved kovtun_csgs.png  | c±(φ=0,v0=0.9)=", round.(kovtun_cv(0.9, 0.0; c0=0.5), digits=4))
