#=
    Kovtun (1907.08191) dispersion-relation figures — reproductions of
    piccvphi.pdf, picreshearv09.pdf, picimshearv09.pdf.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/kovtun.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

# ---- piccvphi: c_v(φ) vs v0, c0=1/2 ----
fig1 = Figure(size=(560, 380))
ax1 = Axis(fig1[1,1], xlabel="v₀", ylabel="c_v(φ)", title="Kovtun fig. cvphi (c₀=1/2)")
v0s = range(-1, 1; length=400)
for φ in range(0, π; length=13)
    cp = [kovtun_cv(v0, φ; c0=0.5)[1] for v0 in v0s]
    cm = [kovtun_cv(v0, φ; c0=0.5)[2] for v0 in v0s]
    lines!(ax1, v0s, cp, linewidth=1.6); lines!(ax1, v0s, cm, linewidth=1.6)
end
hlines!(ax1, [0.0], color=:black, linewidth=0.8); vlines!(ax1, [0.0], color=:black, linewidth=0.8)
ylims!(ax1, -1.05, 1.05)
save(joinpath(outdir, "kovtun_cvphi.png"), fig1)

# ---- shear channel Re/Im ω(k), v0=0.9, θ/η=2 ----
ks = range(0, 1.5; length=300)
φs = range(0, π/2; length=9)
figR = Figure(size=(560, 380))
axR = Axis(figR[1,1], xlabel="(η/w₀) k", ylabel="(η/w₀) Re ω", title="Kovtun fig. reshearv09 (v₀=0.9, θ/η=2)")
figI = Figure(size=(560, 380))
axI = Axis(figI[1,1], xlabel="(η/w₀) k", ylabel="(η/w₀) Im ω", title="Kovtun fig. imshearv09 (v₀=0.9, θ/η=2)")
for φ in φs
    w1 = [kovtun_shear_modes(k, φ; v0=0.9, θη=2.0)[1] for k in ks]
    w2 = [kovtun_shear_modes(k, φ; v0=0.9, θη=2.0)[2] for k in ks]
    lines!(axR, ks, real.(w1), linewidth=1.8); lines!(axR, ks, real.(w2), linewidth=1.8)
    lines!(axI, ks, imag.(w1), linewidth=1.8); lines!(axI, ks, imag.(w2), linewidth=1.8)
end
lines!(axR, ks, collect(ks), color=:black, linestyle=:dash); lines!(axR, ks, -collect(ks), color=:black, linestyle=:dash)
ylims!(axR, -1.5, 1.6); ylims!(axI, -0.4, 0.02)
save(joinpath(outdir, "kovtun_reshearv09.png"), figR)
save(joinpath(outdir, "kovtun_imshearv09.png"), figI)
println("saved kovtun_cvphi / reshearv09 / imshearv09")
