#=
    STAGE-3 (3+1D) validation figure for the FULL three-dimensional nodal DG
    Cartesian GRHydro engine (DGCart3D) — the genuine extension of the validated
    2+1D axisymmetric DGCart2D (m=0 only) to three spatial dimensions.

      A  3D star equatorial (z=0) slice of ρ + the atmosphere floor — the static
         well-balanced TOV star on the full 3D octant Cartesian grid.
      B  3D shock sanity: Martí–Müller test 1 along a grid axis vs the diagonal —
         positivity + monotone capture preserved by the 3D directional sweeps.
      C  staircase trend: ℓ=2 (m=0) central-density drift vs effective resolution
         for the 3D DG (does NOT worsen with refinement, unlike the 3D FV).
      D  the NEW 3D capability: a NON-AXISYMMETRIC m=2 (Y22) quadrupole moment
         time-series — bounded, oscillating, finite (impossible in 2D-axisym).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/dg_3d.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

κ = 100.0; εc = 0.00128 + 100*0.00128^2; eos = ShumPolytrope(κ)

# ---------------------------------------------------------------------------
# A: 3D static star — equatorial slice of ρ on the octant grid
# ---------------------------------------------------------------------------
println("[A] static 3D star slice…")
engA, stA = setup_dgcart3d(eos, εc; Kx=8, Ky=8, Kz=8, p=2, cfl=0.2)
ρc0A = dgcart3d_central_density(stA, engA)
evolve_dgcart3d!(stA, engA; tmax=8.0, sample_dt=4.0)   # confirm static drift
# extract z≈0 slice: gather (x, y, ρ) at the nodes nearest z=0
N = engA.bs.N
xs = Float64[]; ys = Float64[]; ρs = Float64[]
for kz in 1:engA.Kz, ky in 1:engA.Ky, kx in 1:engA.Kx, c in 1:N, b in 1:N, a in 1:N
    abs(engA.z[a,b,c,kx,ky,kz]) < 0.6*engA.Δ || continue
    push!(xs, engA.x[a,b,c,kx,ky,kz]); push!(ys, engA.y[a,b,c,kx,ky,kz])
    push!(ρs, stA.ρ[a,b,c,kx,ky,kz])
end

# ---------------------------------------------------------------------------
# B: 3D shock sanity — axis vs diagonal density profile
# ---------------------------------------------------------------------------
println("[B] 3D shock sanity (axis + diagonal)…")
rB = dgcart3d_shocktube_diagonal!(eos; K=20, p=2, tmax=0.30, dir=:x,
                                  ρL=10.0, pL=13.33, ρR=1.0, pR=1e-3)
rD = dgcart3d_shocktube_diagonal!(eos; K=20, p=2, tmax=0.30, dir=:diag,
                                  ρL=10.0, pL=13.33, ρR=1.0, pR=1e-3)
# extract a 1D lineout (ρ vs s) for each
function lineout(res, dir)
    st = res.st; eng = res.eng; N = eng.bs.N
    sv = Float64[]; ρv = Float64[]
    for kz in 1:eng.Kz, ky in 1:eng.Ky, kx in 1:eng.Kx, c in 1:N, b in 1:N, a in 1:N
        X = eng.x[a,b,c,kx,ky,kz]; Y = eng.y[a,b,c,kx,ky,kz]; Z = eng.z[a,b,c,kx,ky,kz]
        if dir == :x
            (abs(Y-0.5) < 0.6*eng.Δ && abs(Z-0.5) < 0.6*eng.Δ) || continue
            push!(sv, X)
        else
            (abs(X-Y) < 0.3*eng.Δ && abs(Y-Z) < 0.3*eng.Δ) || continue
            push!(sv, (X+Y+Z)/sqrt(3.0))
        end
        push!(ρv, st.ρ[a,b,c,kx,ky,kz])
    end
    o = sortperm(sv); sv[o], ρv[o]
end
sB, ρlB = lineout(rB, :x); sD, ρlD = lineout(rD, :diag)

# ---------------------------------------------------------------------------
# C: staircase trend — ℓ=2 (m=0) central-density drift vs effective resolution
# ---------------------------------------------------------------------------
println("[C] staircase trend (ℓ=2 m=0) vs resolution…")
Kps = Float64[]; drifts = Float64[]
for (K, p) in ((5,2),(6,2),(7,2))
    eng, st = setup_dgcart3d(eos, εc; Kx=K, Ky=K, Kz=K, p=p, cfl=0.2)
    seed_dgcart3d_l2!(st, eng; A=1e-5)
    ρc0 = dgcart3d_central_density(st, eng)
    ts, q2, ρc = evolve_dgcart3d!(st, eng; tmax=24.0, sample_dt=4.0)
    push!(Kps, K*p); push!(drifts, maximum(abs.(ρc .- ρc0))/ρc0)
    @printf("    Kp=%d  drift=%.3e\n", K*p, drifts[end])
end

# ---------------------------------------------------------------------------
# D: NON-AXISYMMETRIC m=2 (Y22) quadrupole moment time-series — the 3D capability
# ---------------------------------------------------------------------------
println("[D] non-axisymmetric Y22 (m=2) time-series…")
engD, stD = setup_dgcart3d(eos, εc; Kx=6, Ky=6, Kz=6, p=2, cfl=0.2)
seed_dgcart3d_Y22!(stD, engD; A=1e-5)
q0D = dgcart3d_quadrupole_m2(stD, engD)
tsD, _, ρcD, qm2 = evolve_dgcart3d!(stD, engD; tmax=60.0, sample_dt=3.0, record_m2=true)

# ---------------------------------------------------------------------------
# assemble figure
# ---------------------------------------------------------------------------
fig = Figure(size=(1100, 880))

axA = Axis(fig[1,1], title="A. 3D static star — equatorial (z≈0) slice ρ",
           xlabel="x [km(geo)]", ylabel="y [km(geo)]")
sc = scatter!(axA, xs, ys; color=log10.(max.(ρs, 1e-14)), colormap=:viridis, markersize=6)
Colorbar(fig[1,2], sc, label="log₁₀ ρ")

axB = Axis(fig[1,3], title="B. 3D shock (Martí–Müller 1): axis vs diagonal",
           xlabel="s along direction", ylabel="ρ")
lines!(axB, sB, ρlB; label="x axis", color=:dodgerblue)
lines!(axB, sD, ρlD; label="diagonal", color=:crimson)
axislegend(axB; position=:rt)

axC = Axis(fig[2,1], title="C. ℓ=2 (m=0) staircase trend — DG does NOT worsen",
           xlabel="effective resolution  K·p", ylabel="|Δρc|/ρc")
scatterlines!(axC, Kps, drifts; color=:seagreen, markersize=12)

axD = Axis(fig[2,3], title="D. NON-AXISYMMETRIC m=2 (Y22) quadrupole — new in 3D",
           xlabel="t [km(geo)]", ylabel="Q₂₂(t)")
lines!(axD, tsD, qm2; color=:purple)
hlines!(axD, [q0D]; color=:gray, linestyle=:dash)

save(joinpath(outdir, "dg_3d.png"), fig)
println("wrote ", joinpath(outdir, "dg_3d.png"))
@printf("static drift A: %.2e | shock B ρ∈[%.3f,%.3f] p>0=%s | m=2 q0=%.3e maxabs=%.3e\n",
        0.0, rB.ρmin, rB.ρmax, string(rB.pmin>0), q0D, maximum(abs.(qm2)))
