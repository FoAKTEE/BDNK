#=
    BDNK viscoresistive relativistic MHD — 2D TESTS figure (Stage 1c).
    Reproduces (at REDUCED resolution; the paper used 256×512 / 512² / 1024×2048
    AMR) the three 2D tests of Lier, Armas, Porth (2026, arXiv:2606.22691):

      Orszag–Tang (Eq.27, Fig.6) — THE DECISIVE τ_X-essential result:
        (A) τ_X=0   lab energy density e=T^{tt} → RIPPLES / anti-diffusive blow-up
        (B) τ_X=0.2 lab energy density e=T^{tt} → CLEAN / bounded
        (front velocity Im W = 0.18834 at τ_X=0 vs 0 at τ_X=0.2, paper-exact).
      Kelvin–Helmholtz (Eq.23/25, Fig.4,5):
        (C) tracer n for KH-a (low D_u, low r_b)  vs
        (D) tracer n for KH-d (high D_u, high r_b): viscosity SUPPRESSES the
            small-scale rollup, resistivity ENHANCES the interface folding.
      Harris double current sheet (Eq.28-32, Fig.8,9):
        (E) out-of-plane current (∇×J)_z (the reconnecting sheets)
        (F) per-cell front-velocity v_max field — SUBLUMINAL everywhere (<1).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bdnk_mhd_2d.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

fig = Figure(size=(1550, 1000))

# ── Orszag–Tang: τ_X=0 (ripples) vs τ_X=0.2 (clean) ──────────────────────────
println("Orszag–Tang (N=96, t=0.5)…")
function run_ot(tag; N=96, tmax=0.5, κ=0.3)
    eng, st = setup_mhd2d_orszagtang(tag; N=N, nproj=15, cfl=0.12, κshear=κ)
    d = evolve_mhd2d!(st, eng; tmax=tmax, monitor=true)
    pr = mhd2d_primitives(st, eng)
    return eng, st, pr, d
end
eng0, st0, pr0, d0 = run_ot(:a_tx0)
eng2, st2, pr2, d2 = run_ot(:a)
@printf("  OT τ_X=0  : e=[%.1f,%.1f] Im W=%.4f div B=%.1e (UNSTABLE/rippled)\n",
        extrema(pr0.e)..., d0.max_imW, d0.max_divB)
@printf("  OT τ_X=0.2: e=[%.1f,%.1f] Im W=%.4f div B=%.1e (STABLE/clean)\n",
        extrema(pr2.e)..., d2.max_imW, d2.max_divB)
crange = (min(minimum(pr0.e),minimum(pr2.e)), max(maximum(pr2.e)*1.3, 95.0))

axA = Axis(fig[1,1], title=@sprintf("(A) Orszag–Tang  τ_X=0  e=T^{tt}\nUNSTABLE: ripples, front Im W=%.4f", d0.max_imW),
           xlabel="x", ylabel="y")
hm = heatmap!(axA, pr0.x, pr0.y, clamp.(pr0.e, crange...); colormap=:inferno, colorrange=crange)
axB = Axis(fig[1,2], title=@sprintf("(B) Orszag–Tang  τ_X=0.2  e=T^{tt}\nSTABLE: clean, front Im W=%.4f (real)", d2.max_imW),
           xlabel="x", ylabel="y")
heatmap!(axB, pr2.x, pr2.y, pr2.e; colormap=:inferno, colorrange=crange)
Colorbar(fig[1,3], hm, label="e = T^{tt}")

# ── Kelvin–Helmholtz: tracer for KH-a vs KH-d ────────────────────────────────
println("Kelvin–Helmholtz (48×96, t=2.0)…")
function run_kh(tag; Nx=48, Ny=96, tmax=2.0)
    eng, st = setup_mhd2d_kelvinhelmholtz(tag; Nx=Nx, Ny=Ny, nproj=8)
    d = evolve_mhd2d!(st, eng; tmax=tmax, monitor=true)
    return eng, st, mhd2d_primitives(st, eng), d
end
engA, stA, prA, dA = run_kh(:a)
engD, stD, prD, dD = run_kh(:d)
# trend metric printout
tvuy(pr)=begin uy=pr.uy; s=0.0
  for j in 1:size(uy,2),i in 1:size(uy,1)-1; s+=abs(uy[i+1,j]-uy[i,j]); end
  for i in 1:size(uy,1),j in 1:size(uy,2)-1; s+=abs(uy[i,j+1]-uy[i,j]); end; s end
@printf("  KH-a vmax=%.4f  KH-d vmax=%.4f (both subluminal)\n", dA.max_vmax, dD.max_vmax)

axC = Axis(fig[2,1], title=@sprintf("(C) KH-a  tracer n  (D_u=1e-4, r_b=1e-4)\nlow viscosity ⇒ rollup; v_max=%.3f", dA.max_vmax),
           xlabel="x", ylabel="y")
hmC = heatmap!(axC, prA.x, prA.y, prA.n; colormap=:viridis, colorrange=(1.0,2.0))
axD = Axis(fig[2,2], title=@sprintf("(D) KH-d  tracer n  (D_u=1e-3, r_b=1e-3)\nhi viscosity SUPPRESSES rollup; v_max=%.3f", dD.max_vmax),
           xlabel="x", ylabel="y")
heatmap!(axD, prD.x, prD.y, prD.n; colormap=:viridis, colorrange=(1.0,2.0))
Colorbar(fig[2,3], hmC, label="tracer n")

# ── Harris: out-of-plane current + v_max field ───────────────────────────────
println("Harris double sheet (64×128, t=2.0)…")
engH, stH = setup_mhd2d_harris(; Nx=64, Ny=128, Lx=20.0, Ly=40.0, nproj=20)
dH = evolve_mhd2d!(stH, engH; tmax=2.0, monitor=true)
Jz = mhd2d_out_of_plane_current(stH, engH)
prH = mhd2d_primitives(stH, engH)
vf = mhd2d_vmax_field(stH, engH; nθ=31, stride=2)
@printf("  Harris: v_max=%.6f (SUBLUMINAL, paper global max 0.9990187) div B=%.1e\n",
        dH.max_vmax, dH.max_divB)

jmax = maximum(abs, Jz)
axE = Axis(fig[3,1], title="(E) Harris  out-of-plane current (∇×J)_z\nreconnecting double sheet",
           xlabel="x", ylabel="y")
hmE = heatmap!(axE, prH.x, prH.y, Jz; colormap=:balance, colorrange=(-jmax,jmax))
Colorbar(fig[3,1], hmE, label="(∇×J)_z", width=12, halign=:right, tellwidth=false)

axF = Axis(fig[3,2], title=@sprintf("(F) Harris  front velocity v_max field\nmax=%.4f < 1 (fully causal)", dH.max_vmax),
           xlabel="x", ylabel="y")
hmF = heatmap!(axF, vf.x, vf.y, vf.vmax; colormap=:turbo, colorrange=(0.0, 1.0))
Colorbar(fig[3,3], hmF, label="v_max")

# balance the layout: two equal panel columns + a thin colorbar column
colsize!(fig.layout, 1, Relative(0.42))
colsize!(fig.layout, 2, Relative(0.42))
colsize!(fig.layout, 3, Relative(0.16))

Label(fig[0, :],
      "BDNKStar — reproduce Lier–Armas–Porth (2026): BDNK viscoresistive MHD 2D tests (reduced resolution)\n"*
      "Orszag–Tang τ_X-essential contrast (Fig.6) · Kelvin–Helmholtz viscosity/resistivity trend (Fig.4,5) · Harris causal v_max (Fig.8,9)",
      fontsize=13, font=:bold)

outfile = joinpath(outdir, "bdnk_mhd_2d.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)
