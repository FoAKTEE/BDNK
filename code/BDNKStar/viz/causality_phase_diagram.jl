#=
    causality_phase_diagram.jl  —  BDNK CAUSALITY / WELL-POSEDNESS PHASE DIAGRAM
    figure (companion to repro/bdnk_phase_diagram.jl + repro/causality_reconcile.jl).

      A  SHUM hatted-frame (q̂,ŝ) phase diagram at â=1, η̂=ζ̂=0.01, cs at the
         ShumPolytrope M=1.4 centre.  Shaded by class: causal+well-posed,
         well-posed only, causal only, neither.  Boundaries q̂=ŝ (well-posed,
         Shum eq.71) and c₊=1 (causal) overdrawn; production frame (0.999,1)
         marked at the q̂=ŝ corner.
      B  RADIUS-RESOLVED causality margin 1-c₊(r) along the production-frame
         M=1.4 star, for the soft Shum Γ=2 polytrope (causal everywhere) and a
         realistic stiff SLy star (superluminal core).  c₊=√3·cs(r); the
         well-posed flag 0<q̂<ŝ holds at every radius for both (EOS-independent).
      C  RECONCILIATION: the Kovtun large-k limit IS the characteristic speed —
         numeric large-k slopes of the full sound quartic vs the analytic eq.4.18
         biquadratic agree across frames (the genuine identity), whereas the
         Causality.jl Keeble–Redondo-Yuste biquadratic is a documented split.

    Run: cd code/BDNKStar && julia --project=viz viz/causality_phase_diagram.jl
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "viz"))
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar.EquationOfState: ShumPolytrope, sound_speed2
using .BDNKStar.TOV: solve_tov, mass_solar
using .BDNKStar.Transport: shum_frame_speeds, shum_frame_wellposed
using CairoMakie
using LinearAlgebra: eigvals
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

const ŝP, âP, q̂P = 1.0, 1.0, 0.999
const η̂P, ζ̂P     = 0.01, 0.01

# ---------------------------------------------------------------------------
# Stellar-centre sound speed of the Shum benchmark star.
# ---------------------------------------------------------------------------
eos_shum = ShumPolytrope(100.0); ρ0c = 0.00128; εc = ρ0c + 100*ρ0c^2
cs_centre = sqrt(sound_speed2(eos_shum, εc))

# Shum-frame classifier at â=1 (cp,cm,c0 and well-posed).  Guarded sqrt: in the
# grid scan non-hyperbolic frames (disc<0 or base-√D<0) ARE a class, not an
# error — mirrors classify_shum in repro/bdnk_phase_diagram.jl.
function classify(ŝ, â, q̂, cs)
    V̂   = (4/3)*η̂P + ζ̂P
    D    = q̂^2 + â^2*(4q̂ + (ŝ-1)^2) + 2*â*q̂*(1+ŝ)
    base = â*(1+ŝ) + q̂
    sqD  = D ≥ 0 ? sqrt(D) : NaN
    argp = (base + sqD)/(2â); argm = (base - sqD)/(2â); arg0 = q̂*η̂P/(â*V̂)
    hyp  = (D ≥ 0) && (argp ≥ 0) && (argm ≥ 0) && (arg0 ≥ 0)
    cp   = cs*sqrt(max(argp,0.0)); c0 = cs*sqrt(max(arg0,0.0))
    caus = hyp && cp ≤ 1.0 && c0 ≤ 1.0
    wp   = shum_frame_wellposed(ŝ, â, q̂)
    (cp=cp, caus=caus, wp=wp)
end

# ===========================================================================
# PANEL A — Shum (q̂,ŝ) phase diagram at â=1
# ===========================================================================
nq, ns = 240, 240
q̂s = range(0.02, 2.0, length=nq)
ŝs = range(0.02, 2.0, length=ns)
# class code: 0 neither, 1 causal-only, 2 wp-only, 3 causal+wp
Z = Array{Float64}(undef, nq, ns)
for (i,q̂) in enumerate(q̂s), (j,ŝ) in enumerate(ŝs)
    r = classify(ŝ, âP, q̂, cs_centre)
    Z[i,j] = r.wp ? (r.caus ? 3.0 : 2.0) : (r.caus ? 1.0 : 0.0)
end

# discrete 4-colour map
cmapA = [RGBAf(0.85,0.85,0.88,1), RGBAf(0.99,0.80,0.55,1),
         RGBAf(0.62,0.78,0.95,1), RGBAf(0.40,0.72,0.45,1)]

# ===========================================================================
# PANEL B — radius-resolved causality margin along the M=1.4 star
# ===========================================================================
function margin_profile(eos, εc; h, rmax, in_msun=false)
    star = solve_tov(eos, εc; h=h, ptol_rel=1e-11, rmax=rmax)
    rR = Float64[]; marg = Float64[]
    for i in eachindex(star.r)
        ε = star.ε[i]; ε > 0 || continue
        cs = sqrt(max(sound_speed2(eos, ε), 0.0))
        _, cp, _ = shum_frame_speeds(ŝP, âP, q̂P, η̂P, ζ̂P, cs)
        push!(rR, star.r[i]/star.R); push!(marg, 1 - cp)
    end
    M = in_msun ? mass_solar(star.M) : star.M
    (rR=rR, marg=marg, M=M)
end

Bshum = margin_profile(eos_shum, εc; h=5e-4, rmax=50.0)
# realistic SLy, central ε tuned to M≈1.4 M⊙ (same bracket as bdnk_phase_diagram.jl)
eos_sly = piecewise_polytrope(:SLy)
Mof(e) = mass_solar(solve_tov(eos_sly, e; h=5e-3, ptol_rel=1e-10, rmax=40.0).M)
lo, hi = 4e-4, 2.0e-3; εsly = lo
for _ in 1:60
    mid = sqrt(lo*hi); Mm = Mof(mid)
    Mm < 1.4 ? (global lo = mid) : (global hi = mid)
    global εsly = mid; abs(Mm-1.4) < 1e-3 && break
end
Bsly = margin_profile(eos_sly, εsly; h=5e-3, rmax=40.0, in_msun=true)

# ===========================================================================
# PANEL C — reconciliation: Kovtun large-k numeric vs analytic eq.4.18
# ===========================================================================
const CS_K = 0.5; const GS_K = 1.0; const W0_K = 1.0
function k_largek_num(ε1, ε2, θ, π1; k=1e7)
    cs2 = CS_K^2; k2 = k^2
    c4 = cs2*ε1*θ
    c3 = im*W0_K*(cs2*ε1 + θ)
    c2 = -(W0_K^2 + k2*cs2*(cs2^2*ε1^2 + GS_K*ε1 + (ε2+π1)*(θ - cs2*ε1) + ε2*π1))
    c1 = -im*k2*W0_K*(GS_K + cs2^2*ε1 + cs2*θ)
    c0 = k2*cs2*(W0_K^2 + k2*θ*(cs2*(ε2+π1 - cs2*ε1) - GS_K))
    a = ComplexF64[c0, c1, c2, c3, c4] ./ c4
    C = zeros(ComplexF64,4,4); C[2,1]=1; C[3,2]=1; C[4,3]=1
    for i in 1:4; C[i,4] = -a[i]; end
    maximum((real(z)/k)^2 for z in eigvals(C))
end
function k_largek_an(ε1, ε2, θ, π1)
    cs2 = CS_K^2
    ε̄1 = cs2*ε1/GS_K; ε̄2 = ε2/GS_K; θ̄ = θ/GS_K; π̄1 = π1/GS_K
    a = ε̄1*θ̄/cs2
    b = ε̄1*(ε̄2+π̄1-ε̄1-1/cs2) - θ̄*(ε̄2+π̄1) - ε̄2*π̄1
    c = θ̄*(cs2*(ε̄2+π̄1-ε̄1) - 1)
    (-b + sqrt(b^2 - 4a*c))/(2a)
end
# 6-frame sweep from Shum ratios (ŝ=cs²ε1/γs, â=θ/γs, p̂=π1/γs)
fr = [(1.,1.,8.), (3.,4.,8.), (3.,15.,8.), (3.,30.,8.), (4.,15.,16.), (6.,15.,24.)]
labels = ["(1,1,8)","(3,4,8)","(3,15,8)","(3,30,8)","(4,15,16)","(6,15,24)"]
c2num = Float64[]; c2an = Float64[]
for (ŝ,â,p̂) in fr
    ε1 = ŝ*GS_K/CS_K^2; θ = â*GS_K; π1 = p̂*GS_K
    push!(c2num, k_largek_num(ε1, 0.0, θ, π1))
    push!(c2an,  k_largek_an(ε1, 0.0, θ, π1))
end

# ===========================================================================
# FIGURE
# ===========================================================================
fig = Figure(size=(1500, 470))

# --- Panel A ---
axA = Axis(fig[1,1], xlabel="q̂", ylabel="ŝ",
           title="A  Shum frame phase diagram (â=1, cs=$(round(cs_centre,digits=3)))")
heatmap!(axA, q̂s, ŝs, Z, colormap=cmapA, colorrange=(-0.5,3.5))
lines!(axA, [0,2], [0,2], color=:black, linestyle=:dash, linewidth=2)  # q̂=ŝ WP bdry
# causal boundary c₊=1 contour at â=1 (per ŝ, the q̂ where cp crosses 1)
qcaus = Float64[]; scaus = Float64[]
for ŝ in ŝs
    prev = classify(ŝ, âP, q̂s[1], cs_centre).cp
    for q̂ in q̂s[2:end]
        cp = classify(ŝ, âP, q̂, cs_centre).cp
        if (prev-1)*(cp-1) ≤ 0
            push!(scaus, ŝ); push!(qcaus, q̂); break
        end
        prev = cp
    end
end
!isempty(qcaus) && lines!(axA, qcaus, scaus, color=:firebrick, linewidth=2.5)
scatter!(axA, [q̂P], [ŝP], color=:black, marker=:star5, markersize=22,
         strokecolor=:white, strokewidth=1)
text!(axA, q̂P-0.05, ŝP+0.06; text="production\n(0.999, 1)", align=(:right,:bottom),
      fontsize=11)
text!(axA, 0.30, 1.65; text="causal+well-posed", color=RGBAf(0.15,0.45,0.2,1),
      fontsize=12, font=:bold)
text!(axA, 1.45, 0.35; text="well-posed only", color=RGBAf(0.25,0.4,0.6,1),
      fontsize=11, rotation=0)
xlims!(axA, 0, 2); ylims!(axA, 0, 2)
# legend proxies
elems = [PolyElement(color=cmapA[4]), PolyElement(color=cmapA[3]),
         PolyElement(color=cmapA[2]), PolyElement(color=cmapA[1]),
         LineElement(color=:black, linestyle=:dash),
         LineElement(color=:firebrick)]
Legend(fig[2,1], elems,
       ["causal+well-posed","well-posed only","causal only","neither",
        "q̂=ŝ (eq.71)","c₊=1 (causal)"],
       orientation=:horizontal, nbanks=2, framevisible=true, labelsize=9)

# --- Panel B ---
axB = Axis(fig[1,2], xlabel="r / R", ylabel="causality margin  1 − c₊(r)",
           title="B  radius-resolved margin, production frame (M=1.4)")
hlines!(axB, [0.0], color=:black, linestyle=:dot)
lines!(axB, Bshum.rR, Bshum.marg, color=:seagreen, linewidth=2.5,
       label="Shum Γ=2  (M=$(round(Bshum.M,digits=2)))")
lines!(axB, Bsly.rR, Bsly.marg, color=:crimson, linewidth=2.5,
       label="SLy realistic  (M=$(round(Bsly.M,digits=2)) M⊙)")
# shade the acausal (margin<0) region
band!(axB, [0,1], [minimum(Bsly.marg),minimum(Bsly.marg)], [0,0],
      color=(:crimson,0.08))
text!(axB, 0.04, -0.06; text="superluminal core (c₊>1)", color=:crimson, fontsize=10)
text!(axB, 0.55, 0.30; text="c₊ = √3·cs(r);  0<q̂<ŝ holds ∀r (both)",
      fontsize=10, color=:gray30)
axislegend(axB, position=:rt, framevisible=true, labelsize=10)

# --- Panel C ---
axC = Axis(fig[1,3], xlabel="frame (ŝ,â,p̂)", ylabel="large-k  c∞²",
           title="C  reconciliation: char-speed = Kovtun large-k limit",
           xticks=(1:length(fr), labels), xticklabelrotation=π/6)
scatter!(axC, 1:length(fr), c2num, color=:navy, markersize=15,
         label="numeric  ω/k|_{k→∞}  (full quartic)")
scatter!(axC, 1:length(fr), c2an, color=:orange, marker=:xcross, markersize=15,
         label="analytic eq.4.18 biquadratic")
hlines!(axC, [1.0], color=:black, linestyle=:dash, label="luminal c²=1")
axislegend(axC, position=:rt, framevisible=true, labelsize=9)
text!(axC, 1.0, maximum(c2num)*0.55;
      text="numeric ≡ analytic\n(max rel.err $(round(maximum(abs.(c2num.-c2an)./c2an),sigdigits=2)))",
      fontsize=10, color=:gray25)

Label(fig[0, :], "BDNKStar — Causality phase diagram + reconciliation " *
      "(Shum 2509.15303 eqs.67-71;  Kovtun 1907.08191 eq.4.18)",
      fontsize=15, font=:bold)

save(joinpath(outdir, "causality_phase_diagram.png"), fig)
println("saved causality_phase_diagram.png")
println("  PanelA: production (q̂,ŝ)=(0.999,1) class=", Int(Z[argmin(abs.(q̂s.-q̂P)), argmin(abs.(ŝs.-ŝP))]),
        " (3=causal+wp)")
println("  PanelB: Shum min margin=", round(minimum(Bshum.marg),digits=4),
        "  SLy min margin=", round(minimum(Bsly.marg),digits=4))
println("  PanelC: max |num-an| rel.err=", maximum(abs.(c2num.-c2an)./c2an))
