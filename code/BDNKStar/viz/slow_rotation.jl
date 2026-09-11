#=
    SLOW-ROTATION (Hartle 1967) frame dragging + moment of inertia.

    (A) The Ī–C universal relation: our computed (C, Ī=I/M³) points for the four
        Read et al. (2009) realistic EOS (SLy/APR4/H4/MS1) across a mass range,
        overlaid on the PUBLISHED Breu & Rezzolla (2016) slow-rotation fit
        [MNRAS 459, 646, arXiv:1601.06083, Eq. (20) + Table 2 slow-rot. row]:
            Ī = 0.8134 C⁻¹ + 0.2101 C⁻² + 3.175e-3 C⁻³ − 2.717e-4 C⁻⁴.
    (B) The Lense-Thirring frame-dragging profile ω(r)/Ω vs r/R for a 1.4 M⊙ SLy
        star: dragging is maximal at the center (ω(0)/Ω = central drag fraction)
        and falls off ∝ 1/r³ outside the surface.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/slow_rotation.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
const GCGS = BDNKStar.Units.gram_per_cm3_to_km_minus2

presets = (:SLy, :APR4, :H4, :MS1)
colors  = Dict(:SLy=>:dodgerblue, :APR4=>:seagreen, :H4=>:darkorange, :MS1=>:crimson)
# per-EOS central-density window (km⁻²) spanning ~1 M⊙ up to near max mass
εwin = Dict(
    :SLy  => (6e14*GCGS, 2.4e15*GCGS),
    :APR4 => (6e14*GCGS, 2.6e15*GCGS),
    :H4   => (4e14*GCGS, 1.8e15*GCGS),
    :MS1  => (2.5e14*GCGS, 1.4e15*GCGS),
)

"(C, Ī) sequence over a central-density range for one EOS."
function ic_sequence(eos, εlo, εhi; n=24, h=0.01)
    Cs = Float64[]; Ibars = Float64[]; Ms = Float64[]
    for εc in exp.(range(log(εlo), log(εhi); length=n))
        r = moment_of_inertia(eos, εc; h=h)
        push!(Cs, r.C); push!(Ibars, r.Ibar); push!(Ms, r.M_Msun)
    end
    return Cs, Ibars, Ms
end

"Model nearest to a target mass (for the 1.4 M⊙ profile + markers)."
function near_mass(eos, εlo, εhi, Mt; n=40, h=0.006)
    best=nothing; bd=Inf
    for εc in exp.(range(log(εlo), log(εhi); length=n))
        r = moment_of_inertia(eos, εc; h=h); d=abs(r.M_Msun-Mt)
        if d<bd; bd=d; best=r; end
    end
    return best
end

fig = Figure(size=(1180, 520))

# ── Panel A: Ī–C universal relation ──────────────────────────────────────────
axA = Axis(fig[1,1], xlabel="compactness  C = M/R", ylabel="Ī = I / M³",
           title="Moment of inertia: Ī–C universal relation\n(slow rotation, Hartle 1967; fit: Breu & Rezzolla 2016)")

Cgrid = range(0.08, 0.30; length=200)
lines!(axA, collect(Cgrid), IBAR_FROM_C_BREU.(Cgrid);
       color=:black, linewidth=2.5, label="Breu & Rezzolla (2016) slow-rot. fit")

for sym in presets
    eos = piecewise_polytrope(sym)
    εlo, εhi = εwin[sym]
    Cs, Ibars, _ = ic_sequence(eos, εlo, εhi)
    scatter!(axA, Cs, Ibars; color=colors[sym], markersize=8, label=String(sym))
end
axislegend(axA; position=:rt, framevisible=true, labelsize=11)
xlims!(axA, 0.08, 0.30)

# ── Panel B: frame-dragging profile for a 1.4 M⊙ SLy star ────────────────────
axB = Axis(fig[1,2], xlabel="r / R", ylabel="ω(r) / Ω   (Lense–Thirring dragging)",
           title="Frame dragging in a 1.4 M⊙ star\nω(0)/Ω = central drag fraction")

for sym in presets
    eos = piecewise_polytrope(sym)
    εlo, εhi = εwin[sym]
    r = near_mass(eos, εlo, εhi, 1.4)
    rR = r.r ./ r.R_km
    lines!(axB, rR, r.ω_over_Ω; color=colors[sym], linewidth=2.2,
           label=@sprintf("%s  (M=%.2f, ω₀/Ω=%.2f)", String(sym), r.M_Msun, r.drag_ratio))
end
axislegend(axB; position=:rt, framevisible=true, labelsize=10)
xlims!(axB, 0.0, 1.0); ylims!(axB, 0.0, nothing)

outfile = joinpath(outdir, "slow_rotation.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# ── console validation table ─────────────────────────────────────────────────
println("\nEOS    M[M⊙]   R[km]   C      I[1e45 gcm2]  Ibar   Breu    %Breu   drag ω0/Ω")
for sym in presets
    eos = piecewise_polytrope(sym)
    εlo, εhi = εwin[sym]
    r = near_mass(eos, εlo, εhi, 1.4)
    pB = IBAR_FROM_C_BREU(r.C)
    @printf("%-5s  %5.3f  %6.2f  %.4f  %10.3f   %6.2f  %6.2f  %+5.1f%%  %.3f\n",
            sym, r.M_Msun, r.R_km, r.C, r.I_cgs/1e45, r.Ibar, pB,
            100*(r.Ibar-pB)/pB, r.drag_ratio)
end
println("\nAnchor: PSR J0737-3039A I_1.338 = 1.425e45 g cm² (Reed et al. 2022, arXiv:2204.09000)")
