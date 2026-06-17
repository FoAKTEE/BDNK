#=
    shum_qnm_reanalyze.jl — offline post-processing of the ε_c(t) series dumped
    by repro/shum_qnm.jl (env QNM_DUMP).  Recomputes the Blackman-window PSD and
    the damped-sinusoid decay fit (paper eq. l.725) with finer control over the
    analysis windows, WITHOUT re-running the expensive evolution.  Also performs
    the continuum extrapolation of the f-mode decay rate over Δr using the paper's
    eq. l.757   1/τ_Δr = 1/τ_0 + m (Δr)^p   (Chabanov et al.), to show the trend
    toward the paper's Δr=0.002 value.

    USAGE:  julia --project=. repro/shum_qnm_reanalyze.jl  series1.txt series2.txt ...
    (each file: "# dr=.. t_f=.. ec0=.."  then rows  "t  eps_c").
=#
include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_qnm.jl")
using Printf

function read_series(path)
    dr=NaN; ec0=NaN; ts=Float64[]; ys=Float64[]
    for ln in eachline(path)
        if startswith(ln,"#")
            m = match(r"dr=([0-9.]+).*ec0=([0-9.eE+-]+)", ln)
            if m !== nothing; dr=parse(Float64,m.captures[1]); ec0=parse(Float64,m.captures[2]); end
            continue
        end
        sp = split(strip(ln))
        length(sp) < 2 && continue
        push!(ts, parse(Float64, sp[1])); push!(ys, parse(Float64, sp[2]))
    end
    return dr, ec0, ts, ys
end

# strongest PSD peak inside a kHz band [lo,hi] (parabolically refined)
function band_peak(fk, psd, lo, hi)
    best=-1.0; fbest=NaN
    for i in 2:length(psd)-1
        fk[i] < lo && continue; fk[i] > hi && break
        if psd[i] > psd[i-1] && psd[i] > psd[i+1] && psd[i] > best
            y0,y1,y2 = psd[i-1],psd[i],psd[i+1]; den=(y0-2y1+y2)
            δ = den != 0 ? 0.5*(y0-y2)/den : 0.0
            fbest = fk[i] + δ*(fk[i+1]-fk[i]); best=psd[i]
        end
    end
    return fbest
end

function analyze(path)
    dr, ec0, ts, ys = read_series(path)
    @printf("\n===== %s   (Δr=%.4f, n=%d, t_f=%.0f) =====\n", path, dr, length(ts), ts[end])
    fk, psd = psd_blackman(ts, ys; fmin_kHz=0.5, fmax_kHz=8.0, nf=8000)
    pks = psd_peaks(fk, psd; npk=10)
    @printf("Blackman-PSD global peaks (kHz): ")
    for (f,p) in pks[1:min(6,end)]; @printf("%.3f(%.2f) ", f, p/pks[1][2]); end
    println()
    # band-restricted mode identification (paper Table l.694: F≈2.7, H1≈4.6, H2≈6.36)
    F  = band_peak(fk, psd, 2.0, 3.2)
    H1 = band_peak(fk, psd, 3.8, 5.2)
    H2 = band_peak(fk, psd, 5.6, 6.9)
    @printf("band-restricted: F=%.3f H1=%.3f H2=%.3f kHz  (paper 2.69/4.60/6.36)\n", F,H1,H2)
    λ, ω, fnl, _ = fit_fmode_decay(ts, ys, F; tstart_frac=0.10, label="reanalyze")
    return (dr=dr, F=F, H1=H1, H2=H2, invtau=λ, omega_nl=ω, f_nl=fnl)
end

# continuum extrapolation 1/τ_Δr = 1/τ_0 + m Δr^p   (paper eq. l.757)
# fit (τ0, m, p) by nonlinear least squares over the (Δr, 1/τ) data points.
function extrapolate_decay(drs, rates)
    n = length(drs)
    n < 3 && return (NaN, NaN, NaN)
    # parametrize, fit via simple LM on (τ0, m, p)
    model(q) = q[1] .+ q[2].*(drs.^q[3])
    resid(q) = model(q) .- rates
    q = [minimum(rates)*0.5, 1.0, 1.0]
    λlm=1e-2
    for _ in 1:2000
        r = resid(q); J=zeros(n,3)
        for i in 1:n
            J[i,1]=1.0; J[i,2]=drs[i]^q[3]; J[i,3]=q[2]*drs[i]^q[3]*log(drs[i])
        end
        JtJ=J'J; g=J'r
        for _ in 1:30
            A = JtJ + λlm*Diagonal(diag(JtJ).+1e-30)
            δ = -(A\g); qn=q.+δ
            if sum(abs2,resid(qn)) < sum(abs2,r); q=qn; λlm=max(λlm/3,1e-12); break
            else; λlm*=4; end
        end
    end
    return (q[1], q[2], q[3])
end

function main_re()
    files = ARGS
    isempty(files) && (println("no series files given"); return)
    results = [analyze(f) for f in files]
    println("\n"*"="^70)
    println("CONVERGENCE SUMMARY (case smallSB-F2):")
    @printf("%-8s %-8s %-8s %-8s %-12s %-10s\n","Δr","F","H1","H2","1/τ","ω_nl")
    for r in results
        @printf("%-8.4f %-8.3f %-8.3f %-8.3f %-12.6g %-10.5f\n",
                r.dr, r.F, r.H1, r.H2, r.invtau, r.omega_nl)
    end
    println("paper Δr=0.002:  F=2.69 H1=4.60 H2=6.36  1/τ=0.00157  ω_nl=0.0834")
    drs = [r.dr for r in results]; rates=[r.invtau for r in results]
    if length(results) >= 3
        τ0,m,p = extrapolate_decay(drs, rates)
        @printf("\nEXTRAP (eq. l.757) 1/τ_Δr=1/τ_0+m Δr^p:  1/τ_0=%.6g  m=%.4g  p=%.3f\n", τ0,m,p)
        @printf("   paper continuum 1/τ_0(smallSB-F2)=0.0011 (Table l.777)\n")
    end
    println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_re()
end
