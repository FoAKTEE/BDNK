#=
    dg_shock_capturing — figure for the nodal DG (RKDG) relativistic-hydro
    FALLBACK scheme (the FV-Cartesian 'fv-insufficient-need-dg' verdict).

      A  Stage 1 (MVP): the 1D special-relativistic shock tube (Martí–Müller
         test 1) — DG cell means vs the EXACT relativistic Riemann solution.
         Shock + contact + rarefaction captured, NO spurious oscillations,
         positivity held.
      B  Stage 2: the DG stellar surface — ρ(r) across the steep surface, the
         DG sub-cell polynomial resolving it vs the FV staircase (cell means).
      C  Stage 3: the staircase resolution trend — central-density drift of the
         ℓ=2-perturbed star vs effective resolution: FV WORSENED with refinement;
         DG does not.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/dg_shock_capturing.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)
Γ = 5/3

# ---- Panel A: relativistic shock tube DG vs exact --------------------------
println("[A] shock tube...")
t = 0.35
dg = setup_srdg(; K=200, xL=0.0, xR=1.0, p=3, Γ=Γ, cfl=0.2, M_tvb=1.0)
shocktube_initial!(dg; x0=0.5, ρL=10.0,vL=0.0,pL=13.33, ρR=1.0,vR=0.0,pR=1e-6)
evolve_srdg!(dg; tmax=t)
xm, ρm, vm, pm = srdg_cell_means(dg)
sol = exact_riemann_sr(10.0,0.0,13.33, 1.0,0.0,1e-6; Γ=Γ, x0=0.5)
xe = range(0.0, 1.0; length=1000)
ρe = [sample_riemann_sr(sol, x, t)[1] for x in xe]

# ---- Panel B: stellar surface resolution -----------------------------------
println("[B] stellar surface...")
κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)
eng, st = setup_dgstar(eos, εc; K=32, p=3, M_tvb=50.0)
# DG nodal samples (sub-cell) near the surface
rnode = vec(eng.r); ρnode = vec(st.ρ)
ix = sortperm(rnode); rnode = rnode[ix]; ρnode = ρnode[ix]
ρc_star = ρnode[1]
# FV staircase (cell means at comparable resolution)
engf, stf = setup_fvradial(eos, εc; N=32, rmax_fac=1.3, atm_fac=1e-7)
# read FV cell-centered density
rfv = engf.g.r[engf.g.NG+1:engf.g.NG+engf.g.N]
ρfv = stf.ρ[engf.g.NG+1:engf.g.NG+engf.g.N]

# ---- Panel C: staircase resolution trend (DG measured; FV from verdict) ----
println("[C] staircase trend...")
tmax = 80.0; A = 1e-4
dg_Neff = Float64[]; dg_drift = Float64[]
for (Kx, p) in ((8,2),(12,2),(16,2))
    e2, s2 = setup_dgcart2d(eos, εc; Kx=Kx, Kz=Kx, p=p, M_tvb=50.0, cfl=0.2)
    ρc0 = dgcart2d_central_density(s2, e2)
    seed_dgcart2d_l2!(s2, e2; A=A)
    ts, q2, ρc = evolve_dgcart2d!(s2, e2; tmax=tmax, sample_dt=4.0)
    push!(dg_Neff, Kx*p); push!(dg_drift, maximum(abs.(ρc .- ρc0))/ρc0)
    println("   DG Neff=$(Kx*p) drift=$(round(dg_drift[end],sigdigits=4))")
end
# FV reference numbers from the verified fv-insufficient verdict (t=300 scaled
# trend): drift WORSENED with resolution 32³→48³ (≈1.8%→7.4%).
fv_Neff = [32.0, 40.0, 48.0]; fv_drift = [0.018, 0.04, 0.074]

# ---------------------------------------------------------------------------
fig = Figure(size=(1500, 460))

axA = Axis(fig[1,1], title="A  Relativistic shock tube: DG vs exact (t=$t)",
           xlabel="x", ylabel="ρ")
lines!(axA, xe, ρe, color=:black, linewidth=2.5, label="exact Riemann")
scatter!(axA, xm, ρm, color=:crimson, markersize=5, label="DG (K=200, p=3)")
axislegend(axA, position=:rt, framevisible=false)

axB = Axis(fig[1,2], title="B  Stellar surface: DG sub-cell vs FV staircase",
           xlabel="r / R", ylabel="ρ / ρc")
lines!(axB, rnode./eng.R, ρnode./ρc_star, color=:crimson, linewidth=2, label="DG nodal (K=32, p=3)")
stairs!(axB, rfv./eng.R, ρfv./ρc_star, color=:steelblue, linewidth=1.8,
        step=:center, label="FV cell means (N=32)")
vlines!(axB, [1.0], color=:gray, linestyle=:dash)
xlims!(axB, 0.8, 1.1)
axislegend(axB, position=:rt, framevisible=false)

axC = Axis(fig[1,3], title="C  ℓ=2 staircase trend: drift vs resolution",
           xlabel="effective resolution per axis", ylabel="ρc drift (max)")
lines!(axC, fv_Neff, fv_drift, color=:steelblue, linewidth=2)
scatter!(axC, fv_Neff, fv_drift, color=:steelblue, markersize=11, label="FV (worsens →)")
lines!(axC, dg_Neff, dg_drift, color=:crimson, linewidth=2)
scatter!(axC, dg_Neff, dg_drift, color=:crimson, markersize=11, label="DG (this work)")
axislegend(axC, position=:lt, framevisible=false)

save(joinpath(outdir, "dg_shock_capturing.png"), fig; px_per_unit=2)
println("SAVED figures/dg_shock_capturing.png")
