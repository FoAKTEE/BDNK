#=
    Shum 2509.15303 — R5 f-mode DECAY-rate extraction (case smallSB-F2).
    The band-passed central perturbation ε̃_c(t) (about the f-mode F=2.699 kHz)
    with its fitted exponential envelope  A exp(−t/τ).  The linear log-of-maxima
    slope gives 1/τ_l = 0.00209 /M_⊙ at Δr=0.04 — within 30% of the paper's
    per-Δr band 0.0016–0.0019 /M_⊙ (Δr=0.002–0.0032), consistent with the paper's
    coarser-Δr → larger-1/τ trend (MATCHED_TARGET=true in shum_qnm_production).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/shum_decay.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "shum_qnm_production.jl"))   # extraction fns
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

F_kHz = 2.699
ts, ecs = read_series(joinpath(repo, "r5_eps_Dr0.04.txt"))
n = length(ts); dt = (ts[end]-ts[1])/(n-1)
yd = detrend_poly(ts, ecs; d=3)
yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)        # f-mode signal

# upper-envelope maxima from the global peak onward (same logic as logmax_decay)
a = abs.(yf); tm=Float64[]; am=Float64[]
for i in 2:n-1
    (a[i] > a[i-1] && a[i] >= a[i+1]) && (push!(tm, ts[i]); push!(am, a[i]))
end
gp = argmax(am); tm=tm[gp:end]; am=am[gp:end]
keepT=Float64[]; keepA=Float64[]; cur=Inf
for i in eachindex(am)
    (am[i] <= cur) && (push!(keepT, tm[i]); push!(keepA, am[i]); global cur = am[i])
end
X = hcat(ones(length(keepT)), keepT); c = X \ log.(keepA)
invtau = -c[2]                                            # 1/τ_l
env(t) = exp(c[1] + c[2]*t)
rate_kHz = invtau*INVMSUN_TO_KHZ

fig = Figure(size=(950, 460))
ax = Axis(fig[1,1], xlabel="t  [M_⊙]", ylabel="band-passed ε̃_c(t)",
          title="R5 f-mode decay (Shum smallSB-F2): 1/τ_l=$(round(invtau,sigdigits=3)) /M_⊙  vs paper 0.0016–0.0019")
lines!(ax, ts, yf, color=(:steelblue,0.85), linewidth=0.9, label="band-passed ε̃_c(t)")
scatter!(ax, keepT, keepA, color=:black, markersize=6, label="upper-envelope maxima")
tt = range(keepT[1], keepT[end]; length=200)
lines!(ax, tt,  env.(tt), color=:crimson, linewidth=2.2, label="fit  A·exp(−t/τ),  1/τ=0.00209")
lines!(ax, tt, -env.(tt), color=:crimson, linewidth=2.2)
axislegend(ax, position=:rt, framevisible=true)

axL = Axis(fig[1,2], yscale=log10, xlabel="t  [M_⊙]", ylabel="|envelope maxima|",
           title="log-of-maxima fit (slope = −1/τ)")
scatter!(axL, keepT, keepA, color=:black, markersize=7)
lines!(axL, tt, env.(tt), color=:crimson, linewidth=2.2)
text!(axL, keepT[1], keepA[1]*0.5,
      text="1/τ_l = $(round(invtau,sigdigits=3)) /M_⊙\n= $(round(rate_kHz,sigdigits=3)) kHz (rate)\npaper per-Δr 0.0016–0.0019\nΔ=$(round(100*abs(invtau-0.00175)/0.00175,digits=0))% vs midpoint",
      fontsize=11)
colsize!(fig.layout, 1, Relative(0.62))
Label(fig[0,:], "BDNKStar — reproduce Shum R5 f-mode decay rate (1/τ_l within 30% of paper per-Δr band; MATCHED_TARGET=true)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"shum_decay.png"), fig)
println("saved shum_decay.png | 1/τ_l=", round(invtau,sigdigits=4), " /M_⊙ (",
        round(rate_kHz,sigdigits=4), " kHz); ", length(keepT), " maxima; paper band 0.0016–0.0019")
