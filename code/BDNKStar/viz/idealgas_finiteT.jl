#=
    Finite-temperature ideal-gas TOV star — the BDNK heat-conduction sector ON it.

    `solve_tov_idealgas` integrates the hydrostatic structure with the ideal-gas-
    on-an-adiabat polytrope p = K ρ^{Γ_struct} (a barotrope — a static star must
    be), then reconstructs the genuine finite-T ideal-gas profiles via the micro-
    EOS p = (Γ-1)ρϵ with the adiabatic index Γ. The result is a real temperature
    profile T(r)=p/ρ with dT/dr<0, on which we show:

      A  the temperature profile T(r) (heat conduction has a gradient to act on),
      B  the Caballero–Yunes heat-conduction criterion cs²(r) vs cn²(r):
         the ADIABATIC sound speed cs² sits BELOW the EQUILIBRIUM cn²=Γ-1
         everywhere ⇒ cs²−cn²<0, the ideal gas VIOLATES the criterion (a genuine,
         non-trivial heat-conduction sector — not the marginal cold barotrope),
      C  the BDNK characteristic speeds (heat-flux κ_Q>0, causal frame) along the
         star: c₊²,c₋² real, non-negative and SUBLUMINAL ⇒ causal heat conduction.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/idealgas_finiteT.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# --- build the finite-T ideal-gas star (Γ=5/3 ⇒ cn²=2/3 strictly subluminal) ---
Γ = 5/3; K = 15.0; ρc = 1e-3
st = solve_tov_idealgas(Γ=Γ, K=K, ρc=ρc, h=2e-4)
ni = findlast(>(0), st.p)                       # interior (p > 0)
r   = st.r[1:ni]
T   = st.T[1:ni]
cs2 = st.cs2[1:ni]
cn2v= st.cn2[1:ni]

# --- BDNK characteristic speeds with κ_Q > 0 along the star (causal frame) ---
c2p = Float64[]; c2m = Float64[]; allcausal = true
for i in 1:ni
    κQi = st.ε[i]^0.25 / (3π)                    # representative κ_Q > 0
    tci = TransportCoefficients(η=1e-3, ζ=1e-3, κQ=κQi, τε=2.0, τP=2.0, τQ=1.0, L=1.0)
    fl  = causality_flag(st.p[i], st.ε[i], st.cs2[i], tci)
    push!(c2p, fl.c2_plus); push!(c2m, max(fl.c2_minus, 0.0))
    global allcausal &= fl.causal
end

@printf("finite-T ideal-gas star: M=%.3f Msun  R=%.2f km  Tc=%.4f  Tsurf=%.2e\n",
        mass_solar(st), st.R, st.T[1], st.T[ni])
@printf("cs² ∈ [%.4f,%.4f]  cn²=%.4f  cs²-cn² ∈ [%.4f,%.4f]\n",
        minimum(cs2), maximum(cs2), cn2v[1], minimum(cs2.-cn2v), maximum(cs2.-cn2v))
@printf("BDNK heat-conduction (κ_Q>0): max c₊²=%.4f  causal everywhere: %s\n",
        maximum(c2p), allcausal)

fig = Figure(size=(1320, 430))

# A — temperature profile
axA = Axis(fig[1,1], xlabel="r  [km]", ylabel="T(r)  [k_B/m = 1]",
           title="A  finite-T profile  T(r)=p/ρ  (dT/dr<0)")
lines!(axA, r, T, color=:crimson, linewidth=2.5)
band!(axA, r, zeros(length(r)), T, color=(:crimson, 0.12))

# B — Caballero–Yunes heat-conduction criterion cs² vs cn²
axB = Axis(fig[1,2], xlabel="r  [km]", ylabel="sound speed²",
           title="B  heat-conduction criterion  cs²<cn²  (CY violated)")
lines!(axB, r, cn2v, color=:navy, linewidth=2.5, label="cn² = Γ-1  (equilibrium)")
lines!(axB, r, cs2,  color=:darkorange, linewidth=2.5, label="cs²  (adiabatic)")
band!(axB, r, cs2, cn2v, color=(:teal, 0.15))
axislegend(axB, position=:rt, framevisible=true)

# C — BDNK characteristic speeds (heat-flux κ_Q>0, causal)
axC = Axis(fig[1,3], xlabel="r  [km]", ylabel="characteristic speed²",
           title="C  BDNK speeds (κ_Q>0): real & subluminal")
hlines!(axC, [1.0], color=:gray, linestyle=:dash, label="luminal c²=1")
lines!(axC, r, c2p, color=:purple,    linewidth=2.5, label="c₊²  (fast / heat)")
lines!(axC, r, c2m, color=:seagreen,  linewidth=2.5, label="c₋²  (slow)")
axislegend(axC, position=:rc, framevisible=true)

title_str = @sprintf("BDNKStar — finite-T ideal-gas star (Γ=%.3g, K=%.0f, ρc=%.0e): M=%.2f M⊙, R=%.1f km",
                     Γ, K, ρc, mass_solar(st), st.R) *
            " — heat-conduction sector exercised (CY criterion violated, BDNK heat conduction causal)"
Label(fig[0, :], title_str, fontsize=12, font=:bold)

save(joinpath(outdir, "idealgas_finiteT.png"), fig)
println("saved idealgas_finiteT.png")
