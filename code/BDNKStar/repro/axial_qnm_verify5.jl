# Final check: viscosity sweep with the file's (-C1) convention should reproduce
# the MONOTONIC trend of Bussieres Table II A1 column (f decreases, tau increases
# with eta_c). Compares against ALL eight tabulated A1 rows (tex l.495-534).
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/axial_qnm.jl")
using Printf

function main()
    eos, star = build_star()
    ℓ=2; R=star.R; a=1.6*R
    # full A1 column of Table II (tex lines 495-534)
    tableA1 = [(3e29,10.4884,29.5870),(5e29,10.4795,29.6169),(8e29,10.4661,29.6619),
               (1e30,10.4571,29.6917),(3e30,10.3692,29.9752),(5e30,10.2854,30.2463),
               (8e30,10.1659,30.6362),(1e31,10.0898,30.8857)]
    # inviscid seed
    g0 = ω -> matching_residual(star,eos,ω,ℓ,inviscid(); a=a,rmin=1e-3,nint=8000,next=4000,Ncf=800)
    ω0,_,_ = find_qnm(g0, ftau_to_omega(10.5,29.5); tol=1e-10,maxit=100)
    f0,τ0 = omega_to_ftau(ω0)
    @printf("inviscid: (f,τ)=(%.4f, %.4f)\n", f0,τ0)
    @printf("%-9s %-20s %-20s %-16s\n","eta_c","achieved (f,τ)","Table II A1","Δf%% / Δτ%%")
    ωprev = ω0
    maxdf=0.0; maxdτ=0.0
    for (ηc,ft,τt) in tableA1
        visc,_,_ = frameA_viscosity(star, eos, ηc, 10.0)
        gv = ω -> matching_residual(star,eos,ω,ℓ,visc; a=a,rmin=1e-3,nint=8000,next=4000,Ncf=800,surf_cut=1e-3)
        ωv,res,_ = find_qnm(gv, ωprev; tol=1e-9,maxit=150)
        fv,τv = omega_to_ftau(ωv)
        df=100*(fv-ft)/ft; dτ=100*(τv-τt)/τt
        maxdf=max(maxdf,abs(df)); maxdτ=max(maxdτ,abs(dτ))
        @printf("%-9.0e (%.4f, %.4f)   (%.4f, %.4f)   %+.3f / %+.3f\n",
                ηc, fv,τv, ft,τt, df,dτ)
        ωprev = ωv
    end
    @printf("\nMAX |Δf|=%.3f%%  MAX |Δτ|=%.3f%%  over ALL 8 A1 rows\n", maxdf, maxdτ)
    println(maxdf<1.0 && maxdτ<1.0 ? "=> ALL 8 A1 rows reproduced within 1%." :
                                     "=> some rows exceed 1%.")
end
main()
