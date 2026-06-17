#=
    Kovtun 1907.08191 sound-channel figures — reproductions of picresoundv09.pdf
    and picimsoundv09.pdf. Reuses the verified quartic F_sound solver
    `kovtun_sound_modes` + `PARAMS` from repro/kovtun_sound.jl (v0=0.9, c_s=0.5,
    ε1=π1=3/c_s², ε2=0, θ=4, γs=w0=1; Im ω ≤ 0 ⇒ stable frame).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/kovtun_sound.jl
=#
using Pkg
Pkg.activate(@__DIR__)
# the repro script defines kovtun_sound_modes + PARAMS (and prints a validation
# block on load — harmless); it include()s BDNKStar + LinearAlgebra itself.
include(joinpath(@__DIR__, "..", "repro", "kovtun_sound.jl"))
using CairoMakie
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

ks = collect(range(1e-3, 3.0; length=240))
φs = range(0, π/2; length=9)

figR = Figure(size=(580, 400))
axR = Axis(figR[1,1], xlabel="(γₛ/w₀) k", ylabel="(γₛ/w₀) Re ω",
           title="Kovtun picresoundv09 (v₀=0.9 sound channel)")
figI = Figure(size=(580, 400))
axI = Axis(figI[1,1], xlabel="(γₛ/w₀) k", ylabel="(γₛ/w₀) Im ω",
           title="Kovtun picimsoundv09 (v₀=0.9 sound channel)")
for φ in φs
    # continuous mode-tracking: assign each k's 4 roots to 4 tracks by nearest
    # complex distance to the previous k (avoids sort-by-Re jumps at crossings)
    tracks = [ComplexF64[] for _ in 1:4]
    prev = sort(kovtun_sound_modes(ks[1], φ; PARAMS...), by=real)
    for m in 1:4; push!(tracks[m], prev[m]); end
    for ki in 2:length(ks)
        z = kovtun_sound_modes(ks[ki], φ; PARAMS...)
        used = falses(4)
        for m in 1:4
            best = 0; bd = Inf
            for j in 1:4
                used[j] && continue
                d = abs(z[j] - prev[m]); if d < bd; bd = d; best = j; end
            end
            used[best] = true; push!(tracks[m], z[best]); prev[m] = z[best]
        end
    end
    for m in 1:4
        lines!(axR, ks, real.(tracks[m]), linewidth=1.4)
        lines!(axI, ks, imag.(tracks[m]), linewidth=1.4)
    end
end
lines!(axR, ks, collect(ks), color=:black, linestyle=:dash)   # light cone ω=±k
lines!(axR, ks, -collect(ks), color=:black, linestyle=:dash)
ylims!(axR, -3.1, 3.1); ylims!(axI, -0.25, 0.02)

save(joinpath(outdir, "kovtun_resoundv09.png"), figR)
save(joinpath(outdir, "kovtun_imsoundv09.png"), figI)
println("saved kovtun_resoundv09 / kovtun_imsoundv09")
