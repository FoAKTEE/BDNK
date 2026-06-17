#=
    shum_qnm_fmode.jl — focused f-mode decay extraction from a dumped ε_c(t)
    series (paper eq. l.725).  Steps, mirroring paper l.713–728:
      1. remove the slow secular drift (low-frequency content) — the paper uses a
         4th-order Butterworth band-pass; here we (a) subtract a smooth cubic
         drift, then (b) band-pass tightly around the f-mode F-band.
      2. fit ε̃_c(t) = A e^{-t/τ} cos(ω t + φ0) + C over an automatically chosen
         window that brackets the decaying portion of the f-mode envelope.
    We FORCE the analysis to the F-band [2.2,3.1] kHz so we extract the
    FUNDAMENTAL (not the dominant overtone, which at fine Δr can outweigh F).
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_qnm.jl")
using Printf
using LinearAlgebra: qr

function read_series(path)
    dr=NaN; ec0=NaN; ts=Float64[]; ys=Float64[]
    for ln in eachline(path)
        if startswith(ln,"#")
            m = match(r"dr=([0-9.]+).*ec0=([0-9.eE+-]+)", ln)
            m !== nothing && (dr=parse(Float64,m.captures[1]); ec0=parse(Float64,m.captures[2]))
            continue
        end
        sp=split(strip(ln)); length(sp)<2 && continue
        push!(ts,parse(Float64,sp[1])); push!(ys,parse(Float64,sp[2]))
    end
    return dr,ec0,ts,ys
end

# subtract a degree-d polynomial least-squares drift (removes secular trend)
function detrend_poly(ts, ys; d=3)
    t0=ts[1]; tn=(ts.-t0)./(ts[end]-t0)
    V = hcat([tn.^k for k in 0:d]...)
    c = V \ ys
    return ys .- V*c
end

# strongest PSD peak in [lo,hi] kHz (parabolic refine)
function band_peak(fk,psd,lo,hi)
    best=-1.0; fb=NaN
    for i in 2:length(psd)-1
        fk[i]<lo && continue; fk[i]>hi && break
        if psd[i]>psd[i-1] && psd[i]>psd[i+1] && psd[i]>best
            y0,y1,y2=psd[i-1],psd[i],psd[i+1]; den=(y0-2y1+y2)
            δ=den!=0 ? 0.5*(y0-y2)/den : 0.0; fb=fk[i]+δ*(fk[i+1]-fk[i]); best=psd[i]
        end
    end
    return fb
end

function fmode_extract(path; Fband=(2.2,3.1), Fguess=2.69)
    dr,ec0,ts,ys = read_series(path)
    dt=(ts[end]-ts[1])/(length(ts)-1)
    # PSD (Blackman) of the detrended signal -> locate the F-band peak
    yd = detrend_poly(ts, ys; d=3)
    fk,psd = psd_blackman(ts, yd; fmin_kHz=0.5, fmax_kHz=8.0, nf=8000)
    Fpk = band_peak(fk,psd,Fband[1],Fband[2]); isnan(Fpk) && (Fpk=Fguess)
    H1 = band_peak(fk,psd,3.8,5.2); H2 = band_peak(fk,psd,5.6,6.9)
    # tight band-pass around the F-band peak, then damped-sinusoid fit
    yf = butter4_bandpass(yd, Fpk, dt; bw_frac=0.35)
    # window: skip 2 periods, fit until envelope falls below 5% of its early max
    P = INVMSUN_TO_KHZ/Fpk/dt
    aenv = abs.(yf); n=length(yf)
    i1 = max(2, Int(round(2*P)))
    pk = maximum(aenv[i1:min(n,i1+Int(round(8*P)))])
    i2 = n - 3
    for k in i1:n-3
        if aenv[k] < 0.05*pk && k > i1+Int(round(3*P)); i2=k; break; end
    end
    i2 = max(i2, i1+Int(round(4*P))); i2=min(i2,n-3)
    tt = ts[i1:i2].-ts[i1]; yy=yf[i1:i2]
    ω0 = 2π*Fpk/INVMSUN_TO_KHZ
    p0 = [maximum(abs.(yy)), 1e-3, ω0, 0.0, 0.0]
    p,_ = fit_damped_sinusoid(tt, yy, p0)
    A,λ,ω,φ0,C = p
    fnl = ω*INVMSUN_TO_KHZ/(2π)
    @printf("[%s Δr=%.3f] F-band peak=%.3f kHz  H1=%.3f H2=%.3f  | fit window t∈[%.0f,%.0f] (%d pts)\n",
            basename(path), dr, Fpk, H1, H2, ts[i1], ts[i2], length(tt))
    @printf("    -> ω_nl=%.5f /M_⊙  f_nl=%.4f kHz   1/τ=%.6g /M_⊙ (=%.4f kHz)\n",
            ω, fnl, λ, λ*INVMSUN_TO_KHZ)
    return (dr=dr, F=Fpk, H1=H1, H2=H2, omega_nl=ω, f_nl=fnl, invtau=λ)
end

function main_f()
    res = [fmode_extract(f) for f in ARGS]
    println("\n"*"="^72)
    println("F-MODE SUMMARY (smallSB-F2)   vs paper Δr=0.002: F=2.69 ω_nl=0.0834 1/τ=0.00157")
    @printf("%-8s %-8s %-8s %-8s %-10s %-12s\n","Δr","F_kHz","H1","H2","ω_nl","1/τ")
    for r in res
        @printf("%-8.4f %-8.3f %-8.3f %-8.3f %-10.5f %-12.6g\n",
                r.dr, r.F, r.H1, r.H2, r.omega_nl, r.invtau)
    end
    # Richardson extrapolation of F and ω_nl to Δr->0 (assume leading error ~Δr^q)
    if length(res) == 2
        d1,d2 = res[1].dr, res[2].dr
        for (nm,a1,a2) in (("F",res[1].F,res[2].F),("omega_nl",res[1].omega_nl,res[2].omega_nl))
            # 1st-order Richardson with the two Δr (assume p=1): a0 = a1 + (a1-a2)*d1/(d2-d1)
            a0 = a1 + (a1-a2)*d1/(d2-d1)
            @printf("   Richardson(Δr->0, p=1) %s ≈ %.4f\n", nm, a0)
        end
    end
    println("="^72)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_f()
end
