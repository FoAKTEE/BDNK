#=
    Radial Cowling mode figure:
      A  radial pulsation spectrum f_n (n=0..3) for the κ=100,n=1 polytrope star
      B  fundamental f-mode convergence: |f0(N) - f0(2N)| vs N (log–log)

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/radial_modes.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
const G2KM = BDNKStar.Units.gram_per_cm3_to_km_minus2
eos = PolytropeEnergy(100.0, 1.0); εc = 5.5e15 * G2KM

f, ω2, R = radial_cowling_spectrum(eos, εc; N=1600, h_tov=5e-5, nmodes=4)

Ns = [100, 200, 400, 800, 1600]
f0s = [radial_cowling_spectrum(eos, εc; N=N, h_tov=5e-5, nmodes=1)[1][1] for N in Ns]
errs = abs.(diff(f0s))                      # |f0(N)-f0(2N)|

fig = Figure(size=(1150, 460))
axA = Axis(fig[1,1], title="A. Radial Cowling spectrum (κ=100, n=1; R=$(round(R,digits=2)) km)",
           xlabel="mode number n", ylabel="f_n [kHz]", xticks=0:3)
stem!(axA, 0:3, f, color=:dodgerblue, markersize=14, stemwidth=2)
for (i, fi) in enumerate(f)
    text!(axA, i-1, fi+0.4, text="$(round(fi,digits=3))", align=(:center,:bottom), fontsize=12)
end
ylims!(axA, 0, maximum(f)*1.15)
text!(axA, 0.0, maximum(f)*1.05, text="fundamental f₀=$(round(f[1],digits=3)) kHz",
      align=(:left,:top), fontsize=12, color=:crimson)

axB = Axis(fig[1,2], title="B. Fundamental-mode convergence",
           xlabel="grid N", ylabel="|f₀(N) − f₀(2N)| [kHz]", xscale=log10, yscale=log10)
scatterlines!(axB, Float64.(Ns[1:end-1]), errs, color=:purple, linewidth=2.5, markersize=11)
ref = errs[1] .* (Float64.(Ns[1:end-1]) ./ Ns[1]).^(-2)
lines!(axB, Float64.(Ns[1:end-1]), ref, color=:gray, linestyle=:dash, label="∝ N⁻² (2nd order FD)")
axislegend(axB, position=:rt, framevisible=true, labelsize=11)

Label(fig[0, :], "BDNKStar radial Cowling eigensolver (Caballero–Yunes 2506.09149 / NSO sector) — converged, stable spectrum",
      fontsize=14, font=:bold)

outfile = joinpath(outdir, "radial_modes.png")
save(outfile, fig)
println("saved ", outfile, "  f[kHz]=", round.(f, digits=4))
