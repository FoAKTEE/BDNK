#=
    ℓ=2 TIDAL DEFORMABILITY figure (the I-Love-Q completion + dynamical tide).

      A  Λ vs M along the mass sequence for the realistic Read et al. (2009)
         piecewise-polytrope EOS (SLy/APR4/H4/MS1), with the GW170817 (Abbott
         et al. 2018, PRL 121,161101) Λ_1.4 = 190 +390/−120 band and the ≲800
         90%-upper bound: soft SLy/APR4 inside, stiff H4/MS1 disfavored.
      B  I-Love universal relation — our (Ī, Λ) on the Yagi–Yunes (2013,
         Science 341,365; arXiv:1302.4499) fit; all EOS collapse onto one curve.
      C  DYNAMICAL TIDE — Λ_eff(ω)/Λ_static vs GW frequency, the f-mode
         resonance enhancement (Hinderer et al. 2016 / Steinhoff et al. 2016).

    STATIC method: Postnikov–Prakash–Lattimer 2010 y(r) Riccati ODE + Hinderer
    2008 k₂ closed form (TidalDeformability.tidal_love_number); Ī from
    SlowRotation.moment_of_inertia.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/tidal_deformability.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# central energy density (km⁻²) tuned to a target gravitational mass (M⊙)
function εc_for_mass(eos, Mt; lo=3e-4, hi=4e-3)
    f(εc) = mass_solar(solve_tov(eos, εc; h=2.5e-3)) - Mt
    a, b = lo, hi; fa, fb = f(a), f(b); it = 0
    while fa*fb > 0 && it < 40; b *= 1.15; fb = f(b); it += 1; end
    for _ in 1:70
        c = 0.5*(a+b); fc = f(c)
        fa*fc ≤ 0 ? (b = c; fb = fc) : (a = c; fa = fc)
        abs(b-a) < 1e-7 && break
    end
    return 0.5*(a+b)
end

const EOSES = (:SLy, :APR4, :H4, :MS1)
const COL = Dict(:SLy=>:steelblue, :APR4=>:seagreen, :H4=>:darkorange, :MS1=>:crimson)
const MASSES = (1.0, 1.2, 1.4, 1.6, 1.8, 2.0)

# ── gather Λ(M), Ī(M), Λ(1.4) per EOS ───────────────────────────────────────
seqΛ = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()   # (M, Λ)
seqI = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()   # (Λ, Ībar)
Λ14  = Dict{Symbol,Float64}()
for sym in EOSES
    eos = piecewise_polytrope(sym; N=4000)
    Ms = Float64[]; Λs = Float64[]; Ibars = Float64[]
    for Mt in MASSES
        εc = try; εc_for_mass(eos, Mt); catch; continue; end
        star = solve_tov(eos, εc; h=1.5e-3)
        (mass_solar(star) < 0.9) && continue
        _, Λ, _ = tidal_love_number(eos, star)
        Ibar = moment_of_inertia(star).Ibar
        push!(Ms, mass_solar(star)); push!(Λs, Λ); push!(Ibars, Ibar)
    end
    seqΛ[sym] = (Ms, Λs); seqI[sym] = (Λs, Ibars)
    εc14 = εc_for_mass(eos, 1.4); s14 = solve_tov(eos, εc14; h=1.5e-3)
    _, l14, _ = tidal_love_number(eos, s14); Λ14[sym] = l14
end

fig = Figure(size=(1500, 460))

# ── Panel A: Λ vs M with the GW170817 band ──────────────────────────────────
axA = Axis(fig[1,1], xlabel="M [M⊙]", ylabel="dimensionless Λ", yscale=log10,
           title="A  Λ(M) — Read 2009 piecewise polytropes + GW170817")
# GW170817 Λ_1.4 = 190 +390/-120 ⇒ band [70,580]; 90% upper ≲800
hspan!(axA, 70, 580, color=(:gray70, 0.35))
hlines!(axA, [800], color=:black, linestyle=:dash, linewidth=1.5)
vlines!(axA, [1.4], color=:gray60, linestyle=:dot)
text!(axA, 1.02, 800; text="GW170817 ≲800 (90%)", align=(:left,:bottom), fontsize=10)
text!(axA, 1.02, 190; text="Λ₁.₄ 90% band", align=(:left,:center), fontsize=10, color=:gray30)
for sym in EOSES
    Ms, Λs = seqΛ[sym]
    scatterlines!(axA, Ms, Λs, color=COL[sym], markersize=8, label=string(sym))
end
ylims!(axA, 30, 4000); xlims!(axA, 0.95, 2.05)
axislegend(axA, position=:rt, framevisible=true)

# ── Panel B: I-Love universal relation (Ī vs Λ) ─────────────────────────────
axB = Axis(fig[1,2], xlabel="Λ", ylabel="Ī = I/M³", xscale=log10, yscale=log10,
           title="B  I-Love (Yagi-Yunes 2013) — EOS-independent")
# the YY fit curve
Λgrid = 10 .^ range(log10(20), log10(3000); length=200)
lines!(axB, Λgrid, ibar_yagi_yunes.(Λgrid), color=:black, linewidth=2.5,
       label="Yagi-Yunes fit")
maxdev = 0.0
for sym in EOSES
    Λs, Ibars = seqI[sym]
    scatter!(axB, Λs, Ibars, color=COL[sym], markersize=11, label=string(sym))
    for (Λ, Ib) in zip(Λs, Ibars)
        global maxdev = max(maxdev, abs(Ib - ibar_yagi_yunes(Λ))/ibar_yagi_yunes(Λ))
    end
end
text!(axB, 30, 6.5; text="max |Δ| = $(round(100*maxdev,digits=2))%",
      align=(:left,:bottom), fontsize=11)
axislegend(axB, position=:rt, framevisible=true)

# ── Panel C: dynamical Λ_eff(ω)/Λ_static — f-mode resonance ──────────────────
axC = Axis(fig[1,3], xlabel="GW frequency f_GW [kHz]", ylabel="Λ_eff / Λ_static",
           title="C  Dynamical tide — f-mode resonance enhancement")
# reference SLy 1.4 star: full-GR f-mode f_f≈2.0 kHz region; resonance f_GW=2 f_f
Λs_ref = Λ14[:SLy]
for (f_f, lab, col) in ((2.86, "full-GR f-mode (2.86 kHz)", :crimson),
                        (2.0,  "f-mode 2.0 kHz", :steelblue))
    m = lambda_eff_modes(Λs_ref, f_f; g_freqs_kHz=Float64[], g_frac=1e-3)
    f_res = 2*f_f
    fGW = range(0.05, 0.985*f_res; length=400)
    ratio = [lambda_eff(m; f_GW_kHz=fg)/Λs_ref for fg in fGW]
    lines!(axC, collect(fGW), ratio, color=col, linewidth=2.0, label=lab)
    vlines!(axC, [f_res], color=col, linestyle=:dash, linewidth=1.0)
end
hlines!(axC, [1.0], color=:gray60, linestyle=:dot)
text!(axC, 0.1, 1.05; text="static limit Λ_eff(ω→0)=Λ_static",
      align=(:left,:bottom), fontsize=10, color=:gray30)
ylims!(axC, 0.0, 6.0); xlims!(axC, 0, 6.2)
axislegend(axC, position=:lt, framevisible=true)

Label(fig[0, :], "BDNKStar — ℓ=2 TIDAL DEFORMABILITY: static k₂/Λ (PPL2010 y-ODE), " *
      "I-Love-Q completion, GW170817 confrontation, and the dynamical Λ_eff(ω) f-mode resonance",
      fontsize=15, font=:bold)

save(joinpath(outdir, "tidal_deformability.png"), fig)
println("saved tidal_deformability.png")
println("Λ(1.4): ", join(["$s=$(round(Λ14[s],digits=1))" for s in EOSES], "  "))
println("I-Love max |Δ| = ", round(100*maxdev, digits=2), "%")
