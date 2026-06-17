#=
    TOV reproduction figure:
      A  mass–radius sequence (energy polytrope κ=100, n=1) with the Bussières
         et al. (2604.13208) EOS1 benchmark point M=1.27 M☉, R=8.86 km marked
      B  interior profiles p/p_c, m/M, ε/ε_c vs r/R for the benchmark star

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/tov_reproduction.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
const G2KM = BDNKStar.Units.gram_per_cm3_to_km_minus2
eos = PolytropeEnergy(100.0, 1.0)

# --- M-R sequence ---
εcs = [k * G2KM for k in 10 .^ range(log10(5e14), log10(5e16); length=60)]
Ms = Float64[]; Rs = Float64[]
for εc in εcs
    st = solve_tov(eos, εc; h=1e-3)
    push!(Ms, mass_solar(st)); push!(Rs, st.R)
end

# --- benchmark star profiles ---
stB = solve_tov(eos, 3e15 * G2KM; h=5e-4)
MB = mass_solar(stB); RB = stB.R

fig = Figure(size=(1150, 460))
axA = Axis(fig[1,1], title="A. Mass–radius (energy polytrope κ=100, n=1)",
           xlabel="R [km]", ylabel="M [M☉]")
lines!(axA, Rs, Ms, color=:dodgerblue, linewidth=2.5, label="TOV sequence")
scatter!(axA, [8.86], [1.27], color=:crimson, markersize=16, marker=:star5,
         label="Bussières EOS1: (8.86, 1.27)")
scatter!(axA, [RB], [MB], color=:black, markersize=9,
         label="this work: ($(round(RB,digits=2)), $(round(MB,digits=3)))")
axislegend(axA, position=:lt, framevisible=true, labelsize=11)

axB = Axis(fig[1,2], title="B. Interior profiles (Bussières EOS1 star)",
           xlabel="r / R", ylabel="normalized")
rr = stB.r ./ RB
lines!(axB, rr, stB.p ./ stB.p[1], color=:crimson, linewidth=2.5, label="p/p_c")
lines!(axB, rr, stB.m ./ stB.M, color=:seagreen, linewidth=2.5, label="m/M")
lines!(axB, rr, stB.ε ./ stB.ε[1], color=:dodgerblue, linewidth=2.5, label="ε/ε_c")
axislegend(axB, position=:rt, framevisible=true, labelsize=11)

Label(fig[0, :], "BDNKStar TOV — reproduce Bussières EOS1: M=$(round(MB,digits=3)) M☉ (target 1.27), R=$(round(RB,digits=2)) km (target 8.86)",
      fontsize=15, font=:bold)

outfile = joinpath(outdir, "tov_reproduction.png")
save(outfile, fig)
println("saved ", outfile, "  M=", MB, " R=", RB)
