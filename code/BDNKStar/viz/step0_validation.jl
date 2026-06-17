#=
    STEP 0 validation figure (3 panels, log scale):
      A  prim<->cons round-trip max rel error vs |v|, per EOS, with the 1e-10 gate
      B  tabulated EOS interpolation error vs table size N (convergence)
      C  BDNK characteristic speeds c²∓ vs τ_P (causality boundary τ_P ≥ 1)

    Run:  julia --project=code/BDNKStar/viz code/BDNKStar/viz/step0_validation.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures")
isdir(outdir) || mkpath(outdir)

# ---- Panel A data: round-trip error vs |v| -------------------------------
vs = range(0.0, 0.85; length=60)
eos_list = [("polytrope p=κε^{1+1/n}", PolytropeEnergy(100.0, 1.0), :barotropic),
            ("Shum Γ=2 (1C target)",  ShumPolytrope(100.0),       :barotropic),
            ("ideal gas Γ=2",          IdealGas(2.0),              :general)]
roundtrip = Dict{String,Vector{Float64}}()
for (name, eos, kind) in eos_list
    errs = Float64[]
    for v in vs
        if kind == :barotropic
            e = 1e-3
            E, S, _ = prim2cons_barotropic(eos, e, v)
            e2, v2, _, _ = cons2prim_barotropic(eos, E, S)
            push!(errs, max(abs(e2-e)/e, abs(v2-v)) + 1e-18)
        else
            ρ, ϵ = 1e-3, 0.5
            D, S, τ, _ = prim2cons_general(eos, ρ, v, ϵ)
            ρ2, v2, ϵ2, _, _ = cons2prim_general(eos, D, S, τ)
            push!(errs, max(abs(ρ2-ρ)/ρ, abs(v2-v), abs(ϵ2-ϵ)/ϵ) + 1e-18)
        end
    end
    roundtrip[name] = errs
end

# ---- Panel B data: tabulated convergence ---------------------------------
base = PolytropeEnergy(100.0, 1.0)
es_test = range(2e-4, 1.8e-3; length=41)
Ns = [25, 50, 100, 200, 400, 800]
tab_err = Float64[]
for N in Ns
    tb = tabulate(base, 1e-4, 2e-3, N)
    push!(tab_err, maximum(abs(pressure(tb, e) - pressure(base, e))/pressure(base, e) for e in es_test))
end

# ---- Panel C data: characteristic speeds vs τ_P --------------------------
e = 1e-3; p = pressure(base, e); cs2 = sound_speed2(base, e); cs = sqrt(cs2)
τPs = range(0.4, 3.0; length=140)
c2p = Float64[]; c2m = Float64[]; disc = Float64[]
for τP in τPs
    tc = TransportCoefficients(η=1e-2, ζ=1e-2, κQ=0.0, τε=1.0, τP=τP, τQ=1.0, L=1.0)
    cm, cp, d = characteristic_speeds(p, e, cs, tc)
    push!(c2m, cm); push!(c2p, cp); push!(disc, d)
end

# ---- plot ----------------------------------------------------------------
fig = Figure(size=(1500, 470))

axA = Axis(fig[1,1], title="A. prim↔cons round-trip (STEP 0 gate)",
           xlabel="|v|", ylabel="max relative error", yscale=log10)
cols = [:dodgerblue, :crimson, :seagreen]
for (i, (name, _, _)) in enumerate(eos_list)
    lines!(axA, vs, roundtrip[name], color=cols[i], linewidth=2.5, label=name)
end
hlines!(axA, [1e-10], color=:black, linestyle=:dash, linewidth=2)
text!(axA, 0.02, 1.4e-10, text="1e-10 gate", align=(:left,:bottom), fontsize=12)
ylims!(axA, 1e-17, 1e-8)
axislegend(axA, position=:lt, framevisible=true, labelsize=11)

axB = Axis(fig[1,2], title="B. tabulated EOS convergence",
           xlabel="table size N", ylabel="max rel error in p(ε)",
           xscale=log10, yscale=log10)
scatterlines!(axB, Float64.(Ns), tab_err, color=:purple, linewidth=2.5, markersize=11)
# 4th-order reference slope
ref = tab_err[1] .* (Float64.(Ns) ./ Ns[1]).^(-4)
lines!(axB, Float64.(Ns), ref, color=:gray, linestyle=:dash, label="∝ N⁻⁴ (Hermite)")
axislegend(axB, position=:rt, framevisible=true, labelsize=11)

axC = Axis(fig[1,3], title="C. BDNK characteristic speeds vs τ_P",
           xlabel="τ_P", ylabel="c²")
band!(axC, [0.4, 1.0], -0.2, 1.3, color=(:red, 0.08))
text!(axC, 0.45, 1.15, text="τ_P<1: acausal", align=(:left,:top), color=:red, fontsize=12)
lines!(axC, τPs, c2p, color=:crimson, linewidth=2.5, label="c²₊")
lines!(axC, τPs, c2m, color=:dodgerblue, linewidth=2.5, label="c²₋")
hlines!(axC, [1.0], color=:black, linestyle=:dash, label="luminal c²=1")
hlines!(axC, [cs2], color=:seagreen, linestyle=:dot, label="c_s²")
vlines!(axC, [1.0], color=:gray, linestyle=:dashdot)
ylims!(axC, -0.05, 1.25)
axislegend(axC, position=:rb, framevisible=true, labelsize=10)

Label(fig[0, :], "BDNKStar STEP 0 — EOS + primitive recovery + causality monitor (commit-tracked validation)",
      fontsize=16, font=:bold)

outfile = joinpath(outdir, "step0_validation.png")
save(outfile, fig)
println("saved ", outfile)
println("round-trip max over all EOS/|v|: ", maximum(maximum(v) for v in values(roundtrip)))
println("tabulated err N=25 -> N=800: ", tab_err[1], " -> ", tab_err[end])
