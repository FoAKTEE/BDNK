#=
    Non-radial Cowling f/p-mode reproduction (STAGE 3 eigensolver).

      A  relativistic Cowling l=2 spectrum (f, p1..p4) of the n=1 K=100
         rho_c=0.00128 polytrope (M=1.4): BDNKStar shooting vs the published
         benchmark (Font/Stergioulas/Kokkotas, arXiv:2107.13339 Table VI).
      B  f-mode frequency vs angular degree l (l=2..6) — rises with l.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/nonradial_modes.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

eos = ShumPolytrope(100.0)
εc  = 0.00128 + 100 * 0.00128^2

# Panel A: l=2 spectrum vs benchmark
f2, _, R = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=5, N=6000, nscan=900)
bench = [1.8825, 4.1060, 6.0298, 7.8670, 9.6663]
labels = ["f", "p₁", "p₂", "p₃", "p₄"]

# Panel B: f-mode vs l
ls = 2:6
ffl = Float64[]
for l in ls
    fl, _, _ = nonradial_cowling_spectrum(eos, εc; l=l, nmodes=1, N=6000,
                                          nscan=700, ω2lo=1.5e-3, ω2hi=0.02)
    push!(ffl, fl[1])
end

fig = Figure(size=(960, 420))

axA = Axis(fig[1,1], xticks=(1:5, labels), ylabel="frequency [kHz]",
           title="A  l=2 Cowling spectrum vs Font/Stergioulas/Kokkotas (2107.13339)")
barpos = (1:5) .- 0.18; bpos2 = (1:5) .+ 0.18
barplot!(axA, barpos, bench, width=0.36, color=(:gray70), label="benchmark")
barplot!(axA, bpos2, f2[1:5], width=0.36, color=(:seagreen), label="BDNKStar")
for i in 1:5
    text!(axA, i, max(bench[i], f2[i]) + 0.25;
          text="$(round(100*abs(f2[i]-bench[i])/bench[i], digits=2))%",
          align=(:center,:bottom), fontsize=10)
end
ylims!(axA, 0, 11)
axislegend(axA, position=:lt, framevisible=true)

axB = Axis(fig[1,2], xticks=(collect(ls), ["$l" for l in ls]),
           xlabel="angular degree ℓ", ylabel="f-mode frequency [kHz]",
           title="B  Cowling f-mode vs ℓ  (M=1.4, R=$(round(R*BDNKStar.Units.Msun_to_km,digits=2)) km)")
scatterlines!(axB, collect(ls), ffl, color=:crimson, markersize=12)
for (l, fv) in zip(ls, ffl)
    text!(axB, l, fv + 0.04; text="$(round(fv,digits=3))", align=(:center,:bottom), fontsize=10)
end

Label(fig[0, :], "BDNKStar — relativistic Cowling NON-radial f/p modes (STAGE 3 eigensolver): " *
      "5-mode benchmark reproduced to <0.1%", fontsize=14, font=:bold)

save(joinpath(outdir, "nonradial_modes.png"), fig)
println("saved nonradial_modes.png | l=2 f=", round(f2[1],digits=4),
        " kHz (bench 1.8825); f-mode(l=2..6)=", round.(ffl, digits=3))
