#=
    R4 reproduction figure — Bussières 2604.13208 axial ℓ=2 fundamental w-mode:
    frequency f and damping time τ vs heat-conduction/viscosity η_c, comparing
    the BDNKStar shooting+Leaver solver (repro/axial_qnm.jl, independently
    re-verified) against Bussières Table II. (Values are the verified solver
    output; the solver itself is repro/axial_qnm.jl.)

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/axial_qnm.jl
=#
using Pkg
Pkg.activate(@__DIR__)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# verified solver output (repro/axial_qnm.jl) vs Bussières Table II targets
ηc      = [3e29, 1e31]
f_ach   = [10.4879, 10.0898];  f_tgt = [10.4884, 10.0898]
τ_ach   = [29.5786, 30.8748];  τ_tgt = [29.5870, 30.8857]
f_invisc_ach, f_invisc_tgt = 10.5014, 10.50
τ_invisc_ach, τ_invisc_tgt = 29.5338, 29.54

fig = Figure(size=(1150, 460))
axf = Axis(fig[1,1], title="A. ℓ=2 w-mode frequency vs η_c (frame A1)",
           xlabel="η_c", ylabel="f [kHz]", xscale=log10)
scatter!(axf, ηc, f_tgt, marker=:star5, markersize=20, color=:crimson, label="Bussières Table II")
scatter!(axf, ηc, f_ach, marker=:circle, markersize=11, color=:black, label="BDNKStar (shooting+Leaver)")
hlines!(axf, [f_invisc_tgt], color=:gray, linestyle=:dash, label="inviscid (10.50)")
axislegend(axf, position=:lb, framevisible=true, labelsize=11)

axt = Axis(fig[1,2], title="B. w-mode damping time vs η_c",
           xlabel="η_c", ylabel="τ [μs]", xscale=log10)
scatter!(axt, ηc, τ_tgt, marker=:star5, markersize=20, color=:crimson, label="Bussières Table II")
scatter!(axt, ηc, τ_ach, marker=:circle, markersize=11, color=:black, label="BDNKStar")
hlines!(axt, [τ_invisc_tgt], color=:gray, linestyle=:dash, label="inviscid (29.54)")
axislegend(axt, position=:lt, framevisible=true, labelsize=11)

Label(fig[0, :], "BDNKStar R4 — reproduce Bussières 2604.13208 axial w-mode (f,τ) to <0.04% (shooting + Leaver continued fraction)",
      fontsize=14, font=:bold)
outfile = joinpath(outdir, "axial_qnm_reproduction.png")
save(outfile, fig); println("saved ", outfile)
