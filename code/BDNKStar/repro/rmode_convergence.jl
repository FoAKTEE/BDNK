#=
    Convergence analysis of the r-mode / CFS results in RModes.jl.

    (a) Lane–Emden profile resolution N for rmode_validate_lom98 — J̃,Ĩ,τ_GW,τ_sv
    (b) TOV grid resolution h for a realistic-EOS star — window min Ω_c/Ω_K, ν_c
    (c) window temperature-grid density / bisection tolerance
    (d) independent confirmation of the trapezoid integration order
=#
using BDNKStar
using BDNKStar.RModes: trapz, _BULK_CAL
using Printf

richardson_rate(e1, e2, ref_ratio) = log(abs(e1)/abs(e2)) / log(ref_ratio)

# generic "rate from successive Richardson with reference value" using triples
# p = log2( |f_N - f_2N| / |f_2N - f_4N| ) when each step doubles resolution
function step_rates(vals)
    rs = Float64[]
    for i in 1:length(vals)-2
        d1 = vals[i+1]-vals[i]; d2 = vals[i+2]-vals[i+1]
        push!(rs, log2(abs(d1)/abs(d2)))
    end
    return rs
end

println("="^78)
println("(a) LANE–EMDEN PROFILE RESOLUTION N  (rmode_validate_lom98)")
println("="^78)
Ns = [1000, 2000, 4000, 8000, 16000]
Js = Float64[]; Is = Float64[]; tGW = Float64[]; tSV = Float64[]; tBV = Float64[]
for N in Ns
    v = rmode_validate_lom98(; N=N)
    push!(Js, v.J̃); push!(Is, v.Ĩ); push!(tGW, v.τ_GW); push!(tSV, v.τ_sv); push!(tBV, v.τ_bv)
    @printf("  N=%6d  J̃=%.10e  Ĩ=%.10e  τ_GW=%.8e  τ_sv=%.8e  τ_bv=%.8e\n",
            N, v.J̃, v.Ĩ, v.τ_GW, v.τ_sv, v.τ_bv)
end
println("\n  Richardson rates (p = log2(|Δ_N|/|Δ_2N|), resolution doubling):")
for (nm, vv) in (("J̃",Js),("Ĩ",Is),("τ_GW",tGW),("τ_sv",tSV),("τ_bv",tBV))
    rs = step_rates(vv)
    @printf("    %-5s p = %s   (Δ tail: %.3e)\n", nm, join((@sprintf("%.3f",r) for r in rs), ", "),
            vv[end]-vv[end-1])
end

println("\n" * "="^78)
println("(b) TOV GRID RESOLUTION h  (realistic EOS, window min Ω_c/Ω_K and ν_c)")
println("="^78)
const _GC = BDNKStar.Units.gram_per_cm3_to_km_minus2
function eos_star(eos, εc; h)
    solve_tov(eos, εc; h=h)
end
# fixed central density tuned to ~1.4 Msun (find once at fine h, reuse εc so the
# ONLY thing changing is the integration step — clean grid-convergence study)
function tune_center(eos; h=0.004)
    best=nothing; bd=Inf; bestε=0.0
    for εc in exp.(range(log(6e14*_GC), log(2.2e15*_GC); length=60))
        s = solve_tov(eos, εc; h=h)
        d = abs(mass_solar(s)-1.4); if d<bd; bd=d; best=s; bestε=εc; end
    end
    return bestε
end
for sym in (:SLy, :APR4)
    eos = piecewise_polytrope(sym)
    εc = tune_center(eos)
    println("\n  EOS=$sym  (fixed εc=$(εc), tuned to ~1.4 Msun)")
    hs = [0.02, 0.01, 0.005, 0.0025, 0.00125]
    Tg = 10 .^ range(8.5, 10.4; length=120)
    Mvals=Float64[]; Rvals=Float64[]; ratios=Float64[]; nus=Float64[]; Tmins=Float64[]
    for h in hs
        s = solve_tov(eos, εc; h=h)
        w = rmode_instability_window(s, eos; Tgrid=Tg)
        imin = argmin(w.Ω_crit_over_ΩK)
        push!(Mvals, mass_solar(s)); push!(Rvals, s.R)
        push!(ratios, w.Ω_crit_over_ΩK[imin]); push!(nus, w.ν_crit_Hz[imin]); push!(Tmins, w.T[imin])
        @printf("    h=%.5f  Npts=%5d  M=%.5f  R=%.4f  Ω_c/Ω_K=%.8e  ν_c=%.6f Hz  T_min=%.4e\n",
                h, length(s.r), mass_solar(s), s.R, w.Ω_crit_over_ΩK[imin], w.ν_crit_Hz[imin], w.T[imin])
    end
    println("    Richardson rates (p=log2(|Δ_h|/|Δ_{h/2}|), step halving):")
    for (nm,vv) in (("M",Mvals),("R",Rvals),("Ω_c/Ω_K",ratios),("ν_c",nus))
        rs = step_rates(vv)
        @printf("      %-8s p = %s   (Δ tail: %.3e)\n", nm,
                join((@sprintf("%.3f",r) for r in rs), ", "), vv[end]-vv[end-1])
    end
end

println("\n" * "="^78)
println("(c) WINDOW TEMPERATURE-GRID DENSITY + BISECTION TOLERANCE")
println("="^78)
# Use a fixed fine star so only the window discretization varies.
let eos = piecewise_polytrope(:SLy)
    εc = tune_center(eos)
    s = solve_tov(eos, εc; h=0.0025)
    println("\n  (c1) Temperature-grid density NT (fixed [10^8.5,10^10.4] K):")
    NTs = [20, 40, 80, 160, 320, 640]
    rmin=Float64[]; numin=Float64[]; Tmin=Float64[]
    for NT in NTs
        Tg = 10 .^ range(8.5, 10.4; length=NT)
        w = rmode_instability_window(s, eos; Tgrid=Tg)
        imin = argmin(w.Ω_crit_over_ΩK)
        push!(rmin, w.Ω_crit_over_ΩK[imin]); push!(numin, w.ν_crit_Hz[imin]); push!(Tmin, w.T[imin])
        @printf("    NT=%4d  Ω_c/Ω_K(min)=%.8e  ν_c=%.6f Hz  T_min=%.5e K\n",
                NT, w.Ω_crit_over_ΩK[imin], w.ν_crit_Hz[imin], w.T[imin])
    end
    println("    rates (T-grid doubling):")
    for (nm,vv) in (("Ω_c/Ω_K",rmin),("ν_c",numin))
        rs = step_rates(vv)
        @printf("      %-8s p = %s   (Δ tail: %.3e)\n", nm,
                join((@sprintf("%.3f",r) for r in rs), ", "), vv[end]-vv[end-1])
    end

    println("\n  (c2) Bisection iteration count (fixed NT=160) at fixed T near the min:")
    # The module hard-codes 80 bisection iters. Independently bisect inv_tau here
    # to measure how Ω_crit converges with iteration count (tolerance proxy).
    # Reconstruct inv_tau via rmode_timescales at fixed T, scanning Ω.
    Tstar = 2e9
    ΩK = kepler_frequency(s, eos)
    f(Ω) = (1.0 ./ rmode_timescales(s, eos; Ω=Ω, T=Tstar)[4])  # 1/τ_tot
    a0 = 1e-3*ΩK; b0 = ΩK
    fa = f(a0); fb = f(b0)
    @printf("    f(Ωlo)=%.4e  f(Ωhi)=%.4e  (sign change: %s)\n", fa, fb, fa*fb<0)
    if fa*fb < 0
        prev = NaN
        for niter in (10,20,30,40,50,60,70,80)
            a,b,la,lb = a0,b0,fa,fb
            for _ in 1:niter
                m=0.5*(a+b); fm=f(m)
                (la*fm ≤ 0) ? (b=m; lb=fm) : (a=m; la=fm)
            end
            Ωc=0.5*(a+b)
            d = isnan(prev) ? NaN : Ωc-prev
            @printf("    iters=%3d  Ω_c=%.12e  ΔΩ_c=%.3e  bracket=%.3e\n",
                    niter, Ωc, d, b-a)
            prev = Ωc
        end
    end
end

println("\n" * "="^78)
println("(d) INDEPENDENT TRAPEZOID INTEGRATION ORDER CHECK")
println("="^78)
# Integrate a smooth analytic function with known exact value on a uniform grid,
# refine, and confirm error ∝ N^{-2} (second order). Use ∫_0^1 e^x dx = e-1.
let
    exact = exp(1)-1
    Ns = [50,100,200,400,800,1600]
    errs=Float64[]
    for N in Ns
        x = collect(range(0.0,1.0;length=N)); y = exp.(x)
        I = trapz(x,y)
        push!(errs, abs(I-exact))
        @printf("    N=%5d  I=%.12f  err=%.6e\n", N, I, abs(I-exact))
    end
    println("    observed order p = log2(err_N/err_2N):")
    for i in 1:length(errs)-1
        @printf("      N=%5d→%5d  p=%.5f\n", Ns[i], Ns[i+1], log2(errs[i]/errs[i+1]))
    end
    # second smooth test: ∫_0^π sin x dx = 2
    exact2=2.0; errs2=Float64[]
    for N in Ns
        x=collect(range(0.0,π;length=N)); y=sin.(x)
        push!(errs2, abs(trapz(x,y)-exact2))
    end
    println("    sin check order:")
    for i in 1:length(errs2)-1
        @printf("      N=%5d→%5d  p=%.5f\n", Ns[i], Ns[i+1], log2(errs2[i]/errs2[i+1]))
    end
end

println("\nDONE.")
