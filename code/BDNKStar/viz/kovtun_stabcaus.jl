#=
    Kovtun 1907.08191 picstabcaus — sound-channel region that is BOTH stable
    (Routh–Hurwitz) AND short-wavelength causal (lim_{k→∞}|ω/k|<1), in the
    (v_s²ε₁/γs, θ/γs) plane for v_s∈{0.1,0.3,0.5,0.7} (ε₂=0, π₁/γs=3/v_s²). The
    causality cut removes the small-ε̄1 part of picstab (origin ε₁=θ=0 excluded).
    Reuses the verified stability_at_rest + kovtun_sound_modes (large-k speed).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/kovtun_stabcaus.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "repro", "kovtun_sound.jl"))
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

ε̄1 = range(0.0, 11.0; length=240)
θ̄  = range(0.0, 21.0; length=240)
vss = [(0.7, RGBAf(0.85,0.45,0.10,0.55)), (0.5, RGBAf(0.45,0.65,0.10,0.50)),
       (0.3, RGBAf(0.70,0.55,0.30,0.50)), (0.1, RGBAf(0.25,0.45,0.80,0.50))]
const KBIG = 1.0e4

fig = Figure(size=(720, 660))
ax = Axis(fig[1,1], xlabel="v_s² ε₁ / γs", ylabel="θ / γs",
          title="Kovtun picstabcaus — stable AND causal region (ε₂=0, π₁/γs=3/v_s²)")
for (vs, col) in vss
    γs=1.0; cs=vs
    reg = fill(NaN, length(ε̄1), length(θ̄))
    for (i, eb) in enumerate(ε̄1), (j, tb) in enumerate(θ̄)
        ε1 = eb/vs^2; π1 = 3/vs^2; θ = tb
        c1, c2 = stability_at_rest(; cs=cs, ε1=ε1, ε2=0.0, π1=π1, θ=θ, γs=γs)
        stable = (c1 > 0 && c2 > 0)
        if stable
            z = kovtun_sound_modes(KBIG, 0.0; v0=0.0, cs=cs, ε1=ε1, ε2=0.0, π1=π1, θ=θ, γs=γs)
            causal = maximum(abs(real(w))/KBIG for w in z) < 1.0
            (causal) && (reg[i,j] = 1.0)
        end
    end
    heatmap!(ax, ε̄1, θ̄, reg, colormap=[col, col], colorrange=(0.9,1.1))
    text!(ax, vs==0.7 ? 2.7 : vs==0.5 ? 5.5 : vs==0.3 ? 8.7 : 9.5, vs==0.1 ? 9.0 : 15.0,
          text="v_s=$vs", color=:black, fontsize=13)
end
xlims!(ax, 0, 11); ylims!(ax, 0, 21)
Label(fig[0,:], "BDNKStar — reproduce Kovtun picstabcaus (causality removes small-ε₁; origin excluded)",
      fontsize=13, font=:bold)
save(joinpath(outdir, "kovtun_picstabcaus.png"), fig)
println("saved kovtun_picstabcaus.png")
