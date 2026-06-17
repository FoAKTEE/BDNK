#=
    Pandya 2201.12317 CC_plot — WENO5 mixed-derivative commutator diagnostic.
    (A) ∫∫|∂x∂y ξ − ∂y∂x ξ| vs t for ε_W ∈ {1e-3, 1, 1e15} (original's format).
    (B) the same commutator vs ε_W at t=0: nonlinear plateau → peak → collapse to
        the machine-precision floor as ε_W→∞ (the linear-weight limit commutes).

    Mechanism reproduction (PRELIMINARY): the provided reference code is 1D-slab
    only, so there is NO 2D BDNK solver to extract Pandya's exact curves — this
    reproduces the diagnostic's mechanism on a representative 2D smoothing flow.
    Faithfully matched: ε_W=1e15 machine floor + flat; exponential decay in t.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/pandya_cc.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "cc_commutator.jl"))   # engine
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
repo = joinpath(@__DIR__, "..", "repro")

# ---- panel A data: vs t ----
rows = [split(l) for l in readlines(joinpath(repo,"cc_commutator.txt")) if !startswith(l,"#") && !isempty(strip(l))]
t   = [parse(Float64,r[1]) for r in rows]
c3  = [parse(Float64,r[2]) for r in rows]   # ε_W=1e-3
c1  = [parse(Float64,r[3]) for r in rows]   # ε_W=1
c15 = [parse(Float64,r[4]) for r in rows]   # ε_W=1e15
flr(v) = max.(v, 1e-16)

# ---- panel B data: vs ε_W at t=0 ----
eWs = [10.0^p for p in -6:0.5:15]
cE  = [commutator_series(eW; T=0.0, nout=1)[2][1] for eW in eWs]

fig = Figure(size=(1040, 430))
axA = Axis(fig[1,1], yscale=log10, xlabel="t",
           ylabel="∫∫ dy dx |∂x∂y ξ − ∂y∂x ξ|", title="(A) commutator vs t")
lines!(axA, t, flr(c3),  color=:black, linestyle=:dot,   linewidth=2.2, label="ε_W = 1e-3")
lines!(axA, t, flr(c1),  color=:black, linestyle=:solid, linewidth=1.8, label="ε_W = 1")
lines!(axA, t, flr(c15), color=:gray35, linestyle=:dash, linewidth=2.0, label="ε_W = 1e15")
ylims!(axA, 1e-16, 1e1); axislegend(axA, position=:lb, framevisible=true)

axB = Axis(fig[1,2], xscale=log10, yscale=log10, xlabel="ε_W",
           ylabel="commutator (t=0)", title="(B) collapse to machine floor as ε_W→∞")
lines!(axB, eWs, flr(cE), color=:dodgerblue, linewidth=2.4)
scatter!(axB, [1e-3,1.0,1e15], flr([commutator_series(e; T=0.0, nout=1)[2][1] for e in (1e-3,1.0,1e15)]),
         color=[:black,:black,:gray35], markersize=11)
hlines!(axB, [2e-12], color=:crimson, linestyle=:dash, linewidth=1.2)
text!(axB, 1e-5, 4e-12, text="linear-weight floor (∂x,∂y commute)", color=:crimson, fontsize=10)
ylims!(axB, 1e-13, 5e0)

Label(fig[0,:], "BDNKStar — reproduce Pandya CC_plot mechanism: WENO5 mixed-derivative commutator (PRELIMINARY: 2D-surrogate, no 2D ref code)",
      fontsize=12, font=:bold)
save(joinpath(outdir,"pandya_cc_plot.png"), fig)
println("saved pandya_cc_plot.png | floor(ε_W=1e15)=", round(c15[1],sigdigits=2),
        " flat; nonlinear decay c3 ", round(c3[1],sigdigits=2), "→", round(c3[end],sigdigits=2),
        "; vs-εW peak=", round(maximum(cE),sigdigits=3), " @ε_W=", eWs[argmax(cE)])
