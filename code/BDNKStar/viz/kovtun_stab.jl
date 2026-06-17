#=
    Kovtun 1907.08191 picstab — the sound-channel STABILITY region in the
    (v_s²ε₁/γs, θ/γs) plane for v_s ∈ {0.1,0.3,0.5,0.7}, with ε₂=0, π₁/γs=3/v_s²
    (Fig. caption). Reuses the verified at-rest Routh–Hurwitz stability
    conditions `stability_at_rest` from repro/kovtun_sound.jl. The stable region
    is larger for smaller v_s (paper).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/kovtun_stab.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "kovtun_sound.jl"))   # stability_at_rest (+ validation on load)
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

ε̄1 = range(0.0, 11.0; length=320)        # x = v_s²ε₁/γs
θ̄  = range(0.0, 21.0; length=320)        # y = θ/γs
vss = [(0.7, RGBAf(0.85,0.45,0.10,0.55)), (0.5, RGBAf(0.45,0.65,0.10,0.50)),
       (0.3, RGBAf(0.70,0.55,0.30,0.50)), (0.1, RGBAf(0.25,0.45,0.80,0.50))]

fig = Figure(size=(720, 660))
ax = Axis(fig[1,1], xlabel="v_s² ε₁ / γs", ylabel="θ / γs",
          title="Kovtun picstab — sound-channel stability region (ε₂=0, π₁/γs=3/v_s²)")
for (vs, col) in vss
    γs = 1.0; cs = vs
    stable = fill(NaN, length(ε̄1), length(θ̄))
    for (i, eb) in enumerate(ε̄1), (j, tb) in enumerate(θ̄)
        ε1 = eb/vs^2; π1 = 3/vs^2; θ = tb              # ε̄1=v_s²ε1, π̄1=π1=3/v_s²
        c1, c2 = stability_at_rest(; cs=cs, ε1=ε1, ε2=0.0, π1=π1, θ=θ, γs=γs)
        (c1 > 0 && c2 > 0) && (stable[i,j] = 1.0)        # stable ⇔ both conditions > 0
    end
    heatmap!(ax, ε̄1, θ̄, stable, colormap=[col, col], colorrange=(0.9,1.1))
    # label near the top of each band
    text!(ax, vs==0.7 ? 1.5 : vs==0.5 ? 5.0 : vs==0.3 ? 8.5 : 9.5, vs==0.1 ? 9.0 : 15.0,
          text="v_s=$vs", color=:black, fontsize=13)
end
xlims!(ax, 0, 11); ylims!(ax, 0, 21)
Label(fig[0,:], "BDNKStar — reproduce Kovtun picstab (stable region grows as v_s decreases)",
      fontsize=13, font=:bold)
save(joinpath(outdir, "kovtun_picstab.png"), fig)
println("saved kovtun_picstab.png")
