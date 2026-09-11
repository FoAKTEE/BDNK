#=
    STAGE-1B/3 result: AXIAL (odd-parity) VISCOUS quasi-normal modes — the viscosity-driven
    η-mode family, via the module `AxialViscousModes` (full GR: interior shooting + exterior
    vacuum Leaver continued fraction + Wronskian matching + complex-ω root finding, Bussières
    et al. 2604.13208). Contrast with the polar sector: the axial fluid sector is otherwise
    trivial, so viscosity DRIVES new dynamics (no perfect-fluid counterpart) and shifts the
    spacetime w-mode.

      A  fundamental w-mode frequency f vs central shear viscosity η_c (frame A): the
         viscous frequency shift.
      B  damping time τ vs η_c: how dissipation reshapes the mode lifetime.
    Markers tracked by complex continuation in η_c from the inviscid w-mode.

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/axial_viscous_qnm.jl
=#
using Pkg; Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie, Printf
CairoMakie.activate!(type="png")
outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

G = 6.6743015e-11; c = 299_792_458.0
εc = (3e15 * 1e3) * G / c^2 * 1e6                 # EOS1 reference star ρc=3e15 g/cc
eos = PolytropeEnergy(100.0, 1.0)

inv = axial_qnm(eos, εc; l=2, ηc_cgs=0.0)
@printf("inviscid w-mode: f=%.4f kHz, τ=%.3f μs\n", inv.f_kHz, inv.tau_us)

ηcs_all = [3e29, 1e30, 3e30, 1e31, 3e31]         # central shear viscosity [cgs], frame A
ηcs, fs, τs = let ω = inv.omega, xs = Float64[], fs = Float64[], τs = Float64[]
    for ηc in ηcs_all
        r = axial_qnm(eos, εc; l=2, ηc_cgs=ηc, τ̂=10.0, ω0=ω)   # continuation in η_c
        @printf("  η_c=%.1e : f=%.4f kHz, τ=%.3f μs (converged=%s)\n", ηc, r.f_kHz, r.tau_us, r.converged)
        r.converged || continue                                  # keep only converged roots
        push!(xs, ηc); push!(fs, r.f_kHz); push!(τs, r.tau_us); ω = r.omega
    end
    xs, fs, τs
end

fig = Figure(size=(960, 420))
axA = Axis(fig[1,1], xscale=log10, xlabel="central shear viscosity η_c [g cm⁻¹ s⁻¹]",
           ylabel="w-mode frequency f [kHz]", title="A  axial w-mode frequency vs viscosity")
hlines!(axA, [inv.f_kHz], color=:gray, linestyle=:dash, label="inviscid")
scatter!(axA, ηcs, fs, color=:crimson, markersize=11); lines!(axA, ηcs, fs, color=:crimson)
axislegend(axA, position=:lb, framevisible=true)

axB = Axis(fig[1,2], xscale=log10, xlabel="central shear viscosity η_c [g cm⁻¹ s⁻¹]",
           ylabel="damping time τ [μs]", title="B  w-mode damping time vs viscosity")
hlines!(axB, [inv.tau_us], color=:gray, linestyle=:dash, label="inviscid")
scatter!(axB, ηcs, τs, color=:navy, markersize=11); lines!(axB, ηcs, τs, color=:navy)
axislegend(axB, position=:lb, framevisible=true)

Label(fig[0,:], "BDNKStar — STAGE 1B: axial (odd-parity) viscous QNMs of the M≈1.27 M☉ star — " *
      "the η-mode family (full-GR shooting; Bussières et al. 2604.13208)", fontsize=12, font=:bold)
save(joinpath(outdir, "axial_viscous_qnm.png"), fig)
println("saved axial_viscous_qnm.png")
