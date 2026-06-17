#=
    Bussières 2604.13208 — axial ℓ=2 fundamental w-mode (f,τ) VISCOUS TRAJECTORY
    in frame A1 (τ̂=10) as the central viscosity η_c sweeps 1e29→1e31 g/cm/s.
    Same physics as complex_plane_2 (a continuous (f,τ) trace vs η_c); here for the
    frame I validated against Bussières Table II.  The trajectory is ANCHORED by
    the two published Table-II points (η_c=3e29→(10.4884,29.587); 1e31→(10.0898,
    30.886)) — the curve must pass through them.  Continuation tracking from the
    inviscid w-mode (10.50 kHz, 29.54 µs).

    Run: julia code/BDNKStar/repro/axial_qnm_trajectory.jl
=#
include(joinpath(@__DIR__, "axial_qnm.jl"))   # main-guarded → no driver run
using Printf

eos, star = build_star(); R = star.R; ℓ = 2
a = 1.6*R; rmin = 1e-3; nint = 8000; next = 4000; Ncf = 800; surf_cut = 1e-3

# inviscid anchor
g0 = ω -> matching_residual(star, eos, ω, ℓ, inviscid(); a=a, rmin=rmin, nint=nint, next=next, Ncf=Ncf)
ω0, res0, ok0 = find_qnm(g0, ftau_to_omega(10.5, 29.5); tol=1e-10, maxit=100)
f0, τ0 = omega_to_ftau(ω0)
@printf("INVISCID: (f,τ)=(%.4f, %.4f)  ok=%s\n", f0, τ0, ok0)

function sweep(ω0, f0, τ0)
    # the paper's strategy: guess = the INVISCID w-mode for every viscous solve
    # (continuation drifts onto spurious roots).  Validated range η_c∈[3e29,1e31].
    etas = 10 .^ range(log10(3e29), 31.0; length=12)
    open(joinpath(@__DIR__,"axial_qnm_trajectory.txt"),"w") do io
        println(io, "# eta_c  f_kHz  tau_us   (frame A1, ℓ=2 w-mode; guess=inviscid)")
        println(io, "0.0  ", f0, "  ", τ0)
        for ηc in etas
            visc,_,_ = frameA_viscosity(star, eos, ηc, 10.0)
            g = ω -> matching_residual(star, eos, ω, ℓ, visc; a=a, rmin=rmin, nint=nint,
                                       next=next, Ncf=Ncf, surf_cut=surf_cut)
            ωn, res, ok = find_qnm(g, ω0; tol=1e-9, maxit=130)
            if ok && isfinite(real(ωn)) && real(ωn)>0 && 9.5<real(ωn)/(2π)/KHZ_TO_KM<10.8
                f, τ = omega_to_ftau(ωn)
                println(io, ηc, "  ", f, "  ", τ); flush(io)
                @printf("η_c=%.3e : (f,τ)=(%.4f, %.4f)\n", ηc, f, τ)
            else
                @printf("η_c=%.3e : off-branch (f=%.3f) — skipped\n", ηc,
                        real(ωn)/(2π)/KHZ_TO_KM)
            end
        end
    end
end
sweep(ω0, f0, τ0)
println("AXIAL_TRAJECTORY_DONE")
