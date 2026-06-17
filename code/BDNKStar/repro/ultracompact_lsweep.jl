#=
    Bussières 2604.13208 plot_ultracompact — trapped axial w-modes ω(ℓ) for four
    ultracompact compactnesses R/M ∈ {2.40,2.45,2.50,2.60} (𝒞=M/R), ℓ=2..6,
    inviscid + viscous (frame B, η_c=1e31, τ̂=10).  Reuses the verified axial QNM
    solver (axial_ultracompact.jl); finds the fundamental (longest-lived) trapped
    mode per (𝒞,ℓ) and records ω_R·M, −ω_I·M.

    Run: julia code/BDNKStar/repro/ultracompact_lsweep.jl
=#
include(joinpath(@__DIR__, "axial_ultracompact.jl"))   # main-guarded → no driver run
using Printf

const R = 10.0
const rmin=1e-3; const nint=9000; const next=4000; const Ncf=1200
const surf_cut=5e-3; const ηc_cgs=1e31
eos = ConstDensityEOS()
RoverM = [2.40, 2.45, 2.50, 2.60]
ells   = [2,3,4,5,6]

# target_wR = expected Re(ω) for the fundamental branch (branch-continuation):
# eikonal seed for ℓ=2, then previous-ℓ value + typical increment.
function fundamental_trapped(𝒞, ℓ, target_wR)
    star, ρ = const_density_star(𝒞, R); M = star.M; a = 1.6*R
    ΩLR = schwarzschild_OmegaLR(M)
    visc0 = inviscid()
    ωlo = 0.10*(ℓ+0.5)*ΩLR; ωhi = 1.25*(ℓ+0.5)*ΩLR
    seeds,_,_ = scan_trapped_real(star, eos, ℓ, visc0; a=a, rmin=rmin, nint=nint,
        next=next, Ncf=Ncf, surf_cut=0.0, ωre_lo=ωlo, ωre_hi=ωhi, nω=400, imω=-1e-3)
    roots = ComplexF64[]
    for (s,_) in seeds
        ωr,_,ok = find_trapped_mode(star, eos, ℓ, visc0; a=a, rmin=rmin, nint=nint,
            next=next, Ncf=Ncf, surf_cut=0.0, ωguess=s, tol=1e-10, maxit=150)
        if ok && imag(ωr)<0 && imag(ωr)>-0.1 && real(ωr)>0.05 &&
           all(abs(ωr-z)>1e-3 for z in roots); push!(roots, ωr); end
    end
    isempty(roots) && return (M, NaN+0im, NaN+0im)
    # BRANCH SELECTION: root whose ω_R·M is closest to the expected branch value
    sort!(roots, by=z->abs(real(z)*M - target_wR)); ω_id = roots[1]
    visc,_,_ = frameB_viscosity(star, eos, ηc_cgs, 10.0)
    bestv = NaN+0im
    for sc2 in (surf_cut, 1e-2, 3e-3, 7e-3, 2e-3)
        ω_v,_,ok_v = find_trapped_mode(star, eos, ℓ, visc; a=a, rmin=rmin, nint=nint,
            next=next, Ncf=Ncf, surf_cut=sc2, ωguess=ω_id, tol=1e-9, maxit=150)
        # accept only if it stays on the branch (Re near ideal) and physically damped
        if ok_v && !isnan(real(ω_v)) && abs(real(ω_v)-real(ω_id)) < 0.15/M &&
           imag(ω_v) < 0 && imag(ω_v) > -0.2/M
            bestv = ω_v; break
        end
    end
    return (M, ω_id, bestv)
end

const EIK = 1.0/(3*sqrt(3))            # Ω_LR·M = 1/(3√3); eikonal ω·M=(ℓ+0.5)·EIK
open(joinpath(@__DIR__,"ultracompact_lsweep.txt"),"w") do io
    println(io, "# RoverM  ell  wR_id_M  wI_id_M  wR_v_M  wI_v_M")
    for rm in RoverM
        𝒞 = 1.0/rm
        target = 0.80*2.5*EIK            # branch start for ℓ=2 (fundamental ≈0.8 eikonal)
        for ℓ in ells
            M, ωid, ωv = fundamental_trapped(𝒞, ℓ, target)
            isfinite(real(ωid)) && (target = real(ωid)*M + 0.19)   # branch-continue
            @printf(io, "%.2f  %d  %.6f  %.6e  %.6f  %.6e\n", rm, ℓ,
                    real(ωid)*M, imag(ωid)*M, real(ωv)*M, imag(ωv)*M)
            @printf("R/M=%.2f ℓ=%d : ωR·M=%.4f −ωI·M(id)=%.3e (visc)=%.3e\n",
                    rm, ℓ, real(ωid)*M, -imag(ωid)*M, -imag(ωv)*M); flush(io); flush(stdout)
        end
    end
end
println("ULTRACOMPACT_LSWEEP_DONE")
