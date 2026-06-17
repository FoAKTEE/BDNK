#=
    Bussières 2604.13208 Sec.V-C — axial ℓ=2 AVOIDED CROSSING between the w-mode
    and the (viscosity-driven) η-mode, frame B1, EOS1 star, as η_c varies.  The
    w-mode (continued from the inviscid mode) and the η-mode (a second-sound
    resonance, located by a complex scan at large η_c and down-tracked) approach
    in frequency but DO NOT cross [main.tex 591].  Reuses the verified η-mode
    machinery from axial_ultracompact.jl.

    Run: julia code/BDNKStar/repro/axial_avoided_crossing.jl
=#
include(joinpath(@__DIR__, "axial_ultracompact.jl"))   # main-guarded → no driver run
using Printf

eos, star = eta_mode_setup(); R = star.R; M = star.M; ℓ = 2
a = 1.6*R; rmin=1e-3; nint=7000; next=3500; Ncf=800; surf_cut=1e-3; τ̂=10.0
gof(v) = (ω -> matching_residual(star, eos, ω, ℓ, v; a=a, rmin=rmin, nint=nint,
                                 next=next, Ncf=Ncf, surf_cut=surf_cut))

# ---- w-mode (frame B1): continue from the inviscid mode, sweep η_c up ----------
ω0,_,_ = find_qnm(gof(inviscid()), ftau_to_omega(10.5,29.5); tol=1e-10, maxit=100)
fW = omega_to_ftau(ω0); @printf("inviscid w-mode (f,τ)=(%.4f,%.4f)\n", fW[1], fW[2])
ηc_w = 10 .^ range(log10(3e29), log10(5e31); length=10)
wtraj = Tuple{Float64,Float64}[(0.0, real(ω0))]    # (η_c, Re ω)
for ηc in ηc_w
    v,_,_ = frameB_viscosity(star, eos, ηc, τ̂)
    ωr,_,ok = find_qnm(gof(v), ω0; tol=1e-9, maxit=130)   # inviscid-mode guess each
    (ok && real(ωr)>0 && 0.2<real(ωr)<0.45) && push!(wtraj,(ηc, real(ωr)))
end

# ---- η-mode: locate at large η_c, down-track (second-sound, terminates at xing)-
function track_eta()
    v0,η̂0,_ = frameB_viscosity(star, eos, 5e31, τ̂)
    seeds = scan_complex_for_root(gof(v0); re_lo=0.02, re_hi=0.55, im_lo=-0.08, im_hi=-0.0008, nre=80, nim=40)
    ηcands = ComplexF64[]
    for (s,_) in seeds
        ωr,r,ok = find_qnm(gof(v0), s; tol=1e-8, maxit=120)
        (ok && real(ωr)>0.01 && -0.12<imag(ωr)<-1e-6 && real(ωr)<0.6 &&
         all(abs(ωr-z)>3e-3 for z in ηcands)) && push!(ηcands,ωr)
    end
    @printf("    η-mode candidates at 5e31: %d\n", length(ηcands))
    etraj = Tuple{Float64,Float64}[]
    isempty(ηcands) && return etraj
    sort!(ηcands, by=z->abs(imag(z))); ωη = ηcands[1]
    push!(etraj, (5e31, real(ωη)))
    for ηc in 10 .^ range(log10(4.7e31), log10(2.5e31); length=12)
        v,_,_ = frameB_viscosity(star, eos, ηc, τ̂)
        ωr,r,ok = _local_root(gof(v), ωη; w=0.012, n=11)
        if ok && imag(ωr)<-1e-4 && imag(ωr)>-0.25 && real(ωr)>0.02 && abs(ωr-ωη)<0.05
            ωη = ωr; push!(etraj,(ηc, real(ωr)))
        else
            @printf("    η-branch lost at avoided crossing near η_c=%.2e [paper too, 588]\n", ηc); break
        end
    end
    return etraj
end
etraj = track_eta()

open(joinpath(@__DIR__,"axial_avoided_crossing.txt"),"w") do io
    println(io, "# WMODE  eta_c  Re_omega")
    for (e,w) in wtraj; println(io, "W  ", e, "  ", w); end
    println(io, "# ETAMODE  eta_c  Re_omega")
    for (e,w) in etraj; println(io, "E  ", e, "  ", w); end
end
@printf("SAVED axial_avoided_crossing.txt | w-mode pts=%d, η-mode pts=%d\n", length(wtraj), length(etraj))
