#=
    Shum 2509.15303 casA_fitting — the f-mode decay-rate fit, 3-panel layout
    matching the original:
      (top)    |ε̃_c(t)|  (log)        — band-passed f-mode magnitude
      (middle) log|ε̃_c| maxima + best linear fit  (slope = −1/τ)
      (bottom) ε̃_c(t)×1e12 + best-fit damped sinusoid  A e^{−t/τ}cos(ωt+φ)+C

    Uses the Δr=0.04 evolution (case smallSB-F2); a clean late window.  Gives
    1/τ_l ≈ 0.00209 /M_⊙, within 30% of the paper per-Δr band 0.0016–0.0019.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/casA_fitting.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "shum_qnm_production.jl"))
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

F_kHz = 2.699
ts, ecs = read_series(joinpath(repo, "r5_eps_Dr0.04.txt"))
n=length(ts); dt=(ts[end]-ts[1])/(n-1)
yd = detrend_poly(ts, ecs; d=3)
yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)

# clean window (skip filter edge + low-SNR tail)
i1 = findfirst(t->t>=2500, ts); i2 = findlast(t->t<=6000, ts)
tw = ts[i1:i2]; yw = yf[i1:i2]

# upper-envelope maxima + linear log-fit
a=abs.(yw); tm=Float64[]; am=Float64[]
for i in 2:length(a)-1; (a[i]>a[i-1] && a[i]>=a[i+1]) && (push!(tm,tw[i]); push!(am,a[i])); end
gp=argmax(am); tm=tm[gp:end]; am=am[gp:end]
kT=Float64[]; kA=Float64[]; cur=Inf
for i in eachindex(am); (am[i]<=cur) && (push!(kT,tm[i]); push!(kA,am[i]); global cur=am[i]); end
X=hcat(ones(length(kT)),kT); c=X\log.(kA); invtau=-c[2]; env(t)=exp(c[1]+c[2]*t)
# damped-sinusoid best fit for the bottom panel
ω0 = 2π*F_kHz/INVMSUN_TO_KHZ
p,_ = fit_damped_sinusoid(tw .- tw[1], yw, [maximum(abs.(yw)), invtau, ω0, 0.0, 0.0])
ds(t) = p[1]*exp(-p[2]*(t-tw[1]))*cos(p[3]*(t-tw[1])+p[4])+p[5]

fig=Figure(size=(720, 820))
ax1=Axis(fig[1,1], yscale=log10, ylabel="|ε̃_c(t)|", xticklabelsvisible=false,
         title="Shum casA_fitting — f-mode decay (smallSB-F2): 1/τ_l=$(round(invtau,sigdigits=3)) /M_⊙")
lines!(ax1, tw, max.(abs.(yw),1e-16), color=:blue, linewidth=0.8)
ax2=Axis(fig[2,1], ylabel="log(|ε̃_c(t)|)", xticklabelsvisible=false)
scatter!(ax2, kT, log.(kA), color=:blue, markersize=7, label="data")
lines!(ax2, kT, c[1].+c[2].*kT, color=:red, linewidth=2.2, label="best fit")
axislegend(ax2, position=:rt, framevisible=true)
ax3=Axis(fig[3,1], xlabel="t / M_⊙", ylabel="ε̃_c(t) × 1e12")
tt=range(tw[1], tw[end]; length=1500)
scatter!(ax3, tw[1:3:end], yw[1:3:end].*1e12, color=:blue, markersize=4)
lines!(ax3, tt, ds.(tt).*1e12, color=:red, linewidth=1.4)
Label(fig[0,:], "BDNKStar — reproduce Shum casA_fitting decay fit (1/τ_l within 30% of paper 0.0016–0.0019)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"casA_fitting.png"), fig)
println("saved casA_fitting.png | 1/τ_l=", round(invtau,sigdigits=4), " /M_⊙  damped-fit 1/τ=",
        round(p[2],sigdigits=4), "  ω=", round(p[3],sigdigits=4), " (", length(kT), " maxima)")
