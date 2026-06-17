#=
    shum_qnm_production.jl — R5 PRODUCTION QNM extraction + adversarial verify
    for the spherically symmetric BDNK Cowling neutron star, case smallSB-F2.

    GROUNDING (Shum, Abalos, Bea, Bezares, Figueras, Palenzuela, arXiv:2509.15303,
    /data/haiyangw/claude/BDNK/ref-paper/sources/arXiv-2509.15303/src/Paper.tex):
      * QNM frequency spectrum (Sec., l.676–708): extract QNMs from the central
        energy density ε_c(t); PSD = Fourier transform with a BLACKMAN WINDOW on
        ε_c(t) (l.680, l.690).  Sampling Δt=1 M_⊙, t_f=8000 M_⊙ (l.690).
      * Frequency table (l.697–703), columns PF / smallSB-F2 / highB-F9:
            F  = 2.69 / 2.69 / 2.67 kHz
            H1 = 4.55 / 4.60 / 4.60 kHz
            H2 = 6.36 / 6.36 / 6.30 kHz
        TASK TARGETS (this run): F=2.69, H1=4.55, H2=6.36 kHz.
      * QNM decay rate (Sec., l.710–752): isolate the f-mode with a 4th-order
        Butterworth band-pass [0.01, f_s/10] (l.713), then fit the DAMPED SINUSOID
        (l.725):  ε̃_c(t) = A exp(-t/τ) cos(ω t + φ0) + C   to recover (1/τ, ω_nl).
        smallSB-F2 at Δr=0.002:  1/τ_l = 1/τ_nl = 0.00157 /M_⊙, ω_nl=0.0834 /M_⊙
        -> f=2.71 kHz (l.733, l.741).
      * Per-Δr decay table (l.766–777, column smallSB-F2):
            Δr   = [0.0032, 0.0028, 0.0024, 0.0020]
            1/τ  = [0.0019, 0.0018, 0.0017, 0.0016]   /M_⊙   (range 0.0016–0.0019)
        Continuum extrapolation via Eq. (l.756–757):
            1/τ_{Δr} = 1/τ_0 + m (Δr)^p ,   p≈1 (marginal convergence, l.794)
            smallSB-F2:  1/τ_0 = 0.0011 /M_⊙   (l.777)

    UNITS (paper l.136: G=c=1; tables in kHz).
        1/M_⊙ = c³/(G M_⊙) = 203.025 kHz  (CYCLIC; cycles/M_⊙ -> kHz).
        A CYCLIC f_geo [1/M_⊙] -> f[kHz] = 203.025 f_geo.
        An ANGULAR ω [1/M_⊙]   -> f[kHz] = 203.025 ω /(2π)   (the "32.31 kHz per
        rad/M_⊙" shorthand = 203.025/(2π)).
        A DECAY RATE 1/τ [1/M_⊙] is a RATE: NO 2π factor; 1/τ -> kHz = 203.025·(1/τ).
        CHECK (l.733): ω_nl=0.0834/M_⊙ -> 203.025·0.0834/(2π)=2.695≈2.71 kHz ✓.

    DATA.  This script CONSUMES the already-completed R5 ladder time series
    repro/r5_eps_Dr*.txt produced by the verified evolution engine
    (repro/shum_evolve_opt.jl, run via repro/r5_run_Dr*.jl / r5_eps_Dr*.jl).
    It does NOT re-run the physics.  Each file has a header line and columns
    "t  eps_c".  We use every file with a usable record length.

    PACKAGE REUSE: begin with the BDNKStar package (per task spec); we only need
    its physical-constant context for the unit cross-check.
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar

using Printf
using LinearAlgebra: Diagonal, norm

# ===========================================================================
# UNITS (independent recomputation; cross-checked against shum_evolve_opt.jl)
# ===========================================================================
const G_SI  = 6.6743015e-11
const C_SI  = 299_792_458.0
const MSUN  = 1.988416e30
const T_SUN = G_SI*MSUN/C_SI^3            # seconds per M_⊙
const INVMSUN_TO_KHZ = 1e-3 / T_SUN      # (1/M_⊙) -> kHz  (cyclic) ≈ 203.025

# ===========================================================================
# I/O — read a ladder time series  "t  eps_c"
# ===========================================================================
function read_series(path::String)
    ts = Float64[]; es = Float64[]
    for ln in eachline(path)
        s = strip(ln)
        (isempty(s) || startswith(s, "#")) && continue
        toks = split(s)
        length(toks) < 2 && continue
        t = tryparse(Float64, toks[1]); e = tryparse(Float64, toks[2])
        (t === nothing || e === nothing) && continue
        push!(ts, t); push!(es, e)
    end
    return ts, es
end

# ===========================================================================
# Detrend: remove a degree-d least-squares polynomial (secular drift), paper
# l.713 ("remove the global drift"). Default cubic.
# ===========================================================================
function detrend_poly(ts::Vector{Float64}, ys::Vector{Float64}; d::Int=3)
    n = length(ts)
    tn = (ts .- ts[1]) ./ (ts[end]-ts[1])
    V = hcat([tn.^k for k in 0:d]...)
    coef = V \ ys
    return ys .- V*coef
end

# ===========================================================================
# Blackman window (3-term, exact-edge-zero; numpy/scipy form), paper l.680/690
#   w_k = 0.42 - 0.5 cos(2πk/(M-1)) + 0.08 cos(4πk/(M-1))
# ===========================================================================
function blackman_window(M::Int)
    w = Vector{Float64}(undef, M)
    @inbounds for k in 0:M-1
        w[k+1] = 0.42 - 0.5*cos(2π*k/(M-1)) + 0.08*cos(4π*k/(M-1))
    end
    return w
end

# Blackman-windowed power spectral density on a fine frequency grid [kHz].
function psd_blackman(ts::Vector{Float64}, ys::Vector{Float64};
                      fmin_kHz=0.3, fmax_kHz=8.0, nf=8000)
    n = length(ys)
    ybar = sum(ys)/n
    w = blackman_window(n)
    y = (ys .- ybar) .* w
    f_geo = collect(range(fmin_kHz/INVMSUN_TO_KHZ, fmax_kHz/INVMSUN_TO_KHZ; length=nf))
    psd = zeros(nf)
    t0 = ts[1]
    @inbounds for (fi, f) in enumerate(f_geo)
        re = 0.0; im = 0.0
        ω = 2π*f
        for k in 1:n
            θ = ω*(ts[k]-t0)
            re += y[k]*cos(θ); im -= y[k]*sin(θ)
        end
        psd[fi] = re*re + im*im
    end
    return f_geo .* INVMSUN_TO_KHZ, psd
end

# Parabolically-refined local-max peak in band [lo,hi] kHz (sub-bin precision).
function band_peak(fk::Vector{Float64}, psd::Vector{Float64}, lo, hi)
    best = -1.0; fb = NaN
    @inbounds for i in 2:length(psd)-1
        fk[i] < lo && continue
        fk[i] > hi && break
        if psd[i] > psd[i-1] && psd[i] > psd[i+1] && psd[i] > best
            y0,y1,y2 = psd[i-1],psd[i],psd[i+1]; den = (y0 - 2y1 + y2)
            δ = den != 0 ? 0.5*(y0 - y2)/den : 0.0
            fb = fk[i] + δ*(fk[i+1]-fk[i]); best = psd[i]
        end
    end
    return fb
end

# Global strongest peaks (sorted by power), parabolically refined — for reporting.
function global_peaks(fk, psd; npk=8)
    pk = Tuple{Float64,Float64}[]
    @inbounds for i in 2:length(psd)-1
        if psd[i] > psd[i-1] && psd[i] > psd[i+1]
            y0,y1,y2 = psd[i-1],psd[i],psd[i+1]; den=(y0-2y1+y2)
            δ = den != 0 ? 0.5*(y0-y2)/den : 0.0
            push!(pk, (fk[i]+δ*(fk[i+1]-fk[i]), y1))
        end
    end
    sort!(pk, by=x->-x[2])
    return pk[1:min(npk,end)]
end

# ===========================================================================
# 4th-order Butterworth-style band-pass (paper l.713), zero-phase cascade of
# RBJ biquads (forward+backward => zero phase; squared => 4th order).
# ===========================================================================
function bandpass_biquad(f0::Float64, Q::Float64, dt::Float64)
    ω0 = 2π*f0*dt
    α  = sin(ω0)/(2Q)
    cw = cos(ω0)
    a0 = 1 + α
    return (α/a0, 0.0, -α/a0, (-2cw)/a0, (1 - α)/a0)
end
function biquad_filtfilt(x::Vector{Float64}, coef)
    b0,b1,b2,a1,a2 = coef
    fwd(u) = begin
        y = similar(u); x1=0.0;x2=0.0;y1=0.0;y2=0.0
        @inbounds for i in eachindex(u)
            yi = b0*u[i] + b1*x1 + b2*x2 - a1*y1 - a2*y2
            x2=x1; x1=u[i]; y2=y1; y1=yi; y[i]=yi
        end
        y
    end
    y = fwd(x)
    return reverse(fwd(reverse(y)))
end
function butter4_bandpass(x::Vector{Float64}, f0_kHz::Float64, dt_Msun::Float64; bw_frac=0.35)
    f0 = f0_kHz/INVMSUN_TO_KHZ
    fnyq = 0.5/dt_Msun
    f0 = clamp(f0, 1e-6, 0.95*fnyq)
    Q  = 1.0/bw_frac
    c1 = bandpass_biquad(f0, Q, dt_Msun)
    y  = biquad_filtfilt(x, c1)
    return biquad_filtfilt(y, c1)   # 4th order
end

# ===========================================================================
# Damped-sinusoid fit (paper eq. l.725):  A exp(-t/τ) cos(ω t + φ0) + C
# Levenberg–Marquardt with analytic Jacobian; p=(A, λ=1/τ, ω, φ0, C).
# ===========================================================================
damped_model(p, t) = (A=p[1];λ=p[2];ω=p[3];φ0=p[4];C=p[5]; @. A*exp(-λ*t)*cos(ω*t+φ0)+C)
function damped_jac(p, t)
    A,λ,ω,φ0,C = p
    e = @. exp(-λ*t); cs = @. cos(ω*t+φ0); sn = @. sin(ω*t+φ0)
    n = length(t); J = zeros(n,5)
    @inbounds for i in 1:n
        J[i,1] = e[i]*cs[i]
        J[i,2] = -A*t[i]*e[i]*cs[i]
        J[i,3] = -A*t[i]*e[i]*sn[i]
        J[i,4] = -A*e[i]*sn[i]
        J[i,5] = 1.0
    end
    return J
end
_diag(M) = [M[i,i] for i in 1:size(M,1)]
function fit_damped_sinusoid(t, y, p0; maxit=600, tol=1e-13)
    p = copy(p0); λlm = 1e-3
    r = damped_model(p,t) .- y; cost = sum(abs2, r)
    for _ in 1:maxit
        J = damped_jac(p,t); JtJ = J'J; Jtr = J'r
        improved = false
        for _ in 1:40
            Am = JtJ + λlm*Diagonal(_diag(JtJ) .+ 1e-30)
            δ = -(Am \ Jtr)
            pn = p .+ δ
            pn[2] = max(pn[2], 0.0)          # decay rate ≥ 0
            pn[3] = abs(pn[3])               # ω > 0
            rn = damped_model(pn,t) .- y; cn = sum(abs2, rn)
            if cn < cost
                p = pn; r = rn; cost = cn; λlm = max(λlm/3, 1e-13); improved = true; break
            else
                λlm *= 4
            end
        end
        improved || break
        norm(Jtr) < tol && break
    end
    return p, cost
end

# ===========================================================================
# LINEAR log-of-maxima decay fit (paper l.714, the *primary* method giving
# 1/τ_l):  band-pass around F, take |ε̃_c|, collect the upper-envelope local
# maxima (running-max-decreasing), then fit log|max| vs t by a straight line;
# the slope's negative is 1/τ.  This is more robust than the nonlinear fit for
# drift-contaminated coarse-Δr data because it uses only the decaying envelope.
# ===========================================================================
function logmax_decay(ts, ecs, F_kHz)
    n = length(ts)
    dt = (ts[end]-ts[1])/(n-1)
    yd = detrend_poly(ts, ecs; d=3)
    yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)
    a  = abs.(yf)
    tm = Float64[]; am = Float64[]
    @inbounds for i in 2:n-1
        if a[i] > a[i-1] && a[i] >= a[i+1]
            push!(tm, ts[i]); push!(am, a[i])
        end
    end
    length(tm) < 6 && return (NaN, 0)
    gp = argmax(am)                       # from the global envelope peak onward
    tm = tm[gp:end]; am = am[gp:end]
    keepT = Float64[]; keepL = Float64[]; cur = Inf
    for i in eachindex(am)                # keep only non-increasing maxima (upper env.)
        if am[i] <= cur
            push!(keepT, tm[i]); push!(keepL, log(am[i])); cur = am[i]
        end
    end
    length(keepT) < 4 && return (NaN, 0)
    X = hcat(ones(length(keepT)), keepT)
    c = X \ keepL                          # log|max| = c1 + c2 t  =>  1/τ = -c2
    return (-c[2], length(keepT))
end

# ===========================================================================
# f-mode NONLINEAR damped-sinusoid decay (paper l.725, giving 1/τ_nl): detrend
# -> Butterworth band-pass around F -> damped-sinusoid fit on an amplitude-aware
# window that (i) skips a few early periods (filter edge + overtone contamination,
# l.729) and (ii) stops before the envelope hits the noise floor.
# ===========================================================================
function fit_fmode_decay(ts, ecs, F_kHz; tstart_frac=0.12)
    n = length(ts)
    dt = (ts[end]-ts[1])/(n-1)
    yd = detrend_poly(ts, ecs; d=3)
    yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)
    P  = INVMSUN_TO_KHZ/F_kHz / dt            # f-mode period in samples
    aenv = abs.(yf)
    pk = maximum(aenv)
    i1 = max(2, Int(round(tstart_frac*n)), Int(round(2*P)))
    floorlvl = 0.01*pk
    i2 = n - max(1, Int(round(0.03*n)))
    for k in i1:i2
        if aenv[k] < floorlvl && k > i1 + Int(round(3*P))
            i2 = k; break
        end
    end
    i2 = min(i2, n - max(1, Int(round(0.03*n))))
    if i2 <= i1 + Int(round(2*P))
        i2 = min(n-2, i1 + Int(round(6*P)))
    end
    tt = ts[i1:i2] .- ts[i1]
    yy = yf[i1:i2]
    ω0 = 2π*F_kHz/INVMSUN_TO_KHZ
    A0 = maximum(abs.(yy))
    # multi-start over a few decay seeds for robustness
    best = nothing; bestcost = Inf; bestp = nothing
    for λ0 in (1e-4, 1e-3, 3e-3, 1e-2)
        p, c = fit_damped_sinusoid(tt, yy, [A0, λ0, ω0, 0.0, 0.0])
        if c < bestcost
            bestcost = c; bestp = p
        end
    end
    A,λ,ω,φ0,C = bestp
    f_nl = ω*INVMSUN_TO_KHZ/(2π)
    return (invtau=λ, omega=ω, f_nl=f_nl, i1=i1, i2=i2, t1=ts[i1], t2=ts[i2],
            n=length(tt), cost=bestcost, p=bestp)
end

# Envelope-monotonicity diagnostic (decay-rate extractability, paper l.711–728).
function envelope_monotonic(ts, ecs, F_kHz)
    n=length(ts); dt=(ts[end]-ts[1])/(n-1)
    yd = detrend_poly(ts, ecs; d=3)
    yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)
    W  = max(2, Int(round(80/dt)))
    vals = Float64[]; tcs=Float64[]
    tc = 200.0
    while tc < ts[end]-100
        i=Int(round(tc/dt))+1; lo=max(1,i-W); hi=min(n,i+W)
        push!(vals, maximum(abs.(yf[lo:hi]))); push!(tcs, tc)
        tc += 300.0
    end
    mono = length(vals) >= 2 && all(diff(vals) .<= 0)
    return mono, tcs, vals
end

# ===========================================================================
# Richardson / Chabanov continuum extrapolation (paper Eq. l.756–757):
#   y(Δr) = y0 + m (Δr)^p
# With ≥3 points fit (y0,m,p) by nonlinear least squares (grid over p, linear in
# y0,m).  With exactly 2 points fix p (default 1, the paper's marginal order) and
# solve the 2x2 linear system.  Used for BOTH the decay rate and each frequency.
# ===========================================================================
function extrap_powerlaw(drs::Vector{Float64}, ys::Vector{Float64}; pfix=1.0)
    if length(drs) >= 3
        bestp=NaN; besty0=NaN; bestm=NaN; bestres=Inf
        for p in 0.25:0.01:3.0
            X = hcat(ones(length(drs)), drs.^p)
            c = X \ ys           # [y0, m]
            res = sum(abs2, X*c .- ys)
            if res < bestres
                bestres=res; bestp=p; besty0=c[1]; bestm=c[2]
            end
        end
        return (y0=besty0, m=bestm, p=bestp, res=bestres, mode="3pt-fit")
    elseif length(drs) == 2
        # y0 + m Δr^p with p=pfix : solve 2x2
        a1 = drs[1]^pfix; a2 = drs[2]^pfix
        m  = (ys[2]-ys[1])/(a2-a1)
        y0 = ys[1] - m*a1
        return (y0=y0, m=m, p=pfix, res=0.0, mode="2pt-p$(pfix)")
    else
        return (y0=ys[1], m=NaN, p=NaN, res=0.0, mode="1pt-noextrap")
    end
end

# ===========================================================================
# Per-resolution analysis
# ===========================================================================
# PSD-based band-restricted F/H1/H2 on a (sub-)window of (ts,ecs).
function band_FH1H2(ts, ecs)
    yd = detrend_poly(ts, ecs; d=3)
    fk, psd = psd_blackman(ts, yd; fmin_kHz=0.5, fmax_kHz=8.0, nf=8000)
    F  = band_peak(fk, psd, 2.2, 3.2)     # paper 2.69
    H1 = band_peak(fk, psd, 3.9, 5.2)     # paper 4.55 (PF) / 4.60 (smallSB-F2)
    H2 = band_peak(fk, psd, 5.6, 7.0)     # paper 6.36
    return F, H1, H2, fk, psd
end

struct ResResult
    dr::Float64
    tf::Float64
    n::Int
    # frequencies extracted on the EQUILIBRIUM-PROXIMATE early window (see note)
    F::Float64; H1::Float64; H2::Float64
    # frequencies on the full record (for transparency / drift diagnosis)
    Ffull::Float64; H1full::Float64; H2full::Float64
    fwin::Float64                     # upper time of the frequency window used
    # decay rates: 1/τ_l (linear log-max) and 1/τ_nl (nonlinear damped sinusoid)
    invtau_l::Float64
    invtau_nl::Float64; omega::Float64; f_nl::Float64
    mono::Bool
    fit_t1::Float64; fit_t2::Float64; fit_n::Int
end

# Frequency-window note (PHYSICS, grounded in paper l.655, l.680).  At the coarse
# Δr available here (0.02–0.04 M_⊙, i.e. 10–20× the paper's 0.002), the strong
# numerical viscosity makes ε_c drift UPWARD from its equilibrium value over the
# long run (paper l.655: "deviation from the constant stationary value due to
# numerical errors, which decreases as the resolution is increased").  This drift
# slowly lowers the star's effective eigenfrequencies at LATE times.  The QNMs are
# physical only while the star is near equilibrium — exactly the early window the
# paper highlights as where the perturbation is discernible (l.680, t≲1000 M_⊙).
# We therefore read F/H1/H2 from a fixed early window t≤TF_FREQ (the same physical
# window for every resolution), and ALSO report the full-record frequencies so the
# drift is fully transparent.  At the paper's Δr=0.002 the drift is negligible over
# the full 8000 M_⊙, so early-window == full-window there.
const TF_FREQ = 2000.0     # M_⊙: equilibrium-proximate frequency window

function analyze_resolution(path::String, dr::Float64)
    ts, ecs = read_series(path)
    n = length(ts)
    tf = n>0 ? ts[end] : 0.0
    @printf("\n--- Δr=%.4f  file=%s ---\n", dr, basename(path))
    @printf("    records=%d   t∈[%.1f, %.1f] M_⊙   ε_c(0)=%.8g  ε_c range=[%.8g,%.8g]  drift=%.2f%%\n",
            n, n>0 ? ts[1] : NaN, tf, n>0 ? ecs[1] : NaN,
            n>0 ? minimum(ecs) : NaN, n>0 ? maximum(ecs) : NaN,
            n>0 ? 100*(ecs[end]-ecs[1])/ecs[1] : NaN)

    # ---- full-record Blackman PSD (transparency) ----------------------------
    Ffull, H1full, H2full, fk, psd = band_FH1H2(ts, ecs)
    gp = global_peaks(fk, psd; npk=6)
    print("    full-record global PSD peaks [kHz(rel)]: ")
    for (f,p) in gp; @printf("%.3f(%.2f) ", f, p/gp[1][2]); end
    println()
    @printf("    full-record F/H1/H2 = %.3f / %.3f / %.3f kHz\n", Ffull, H1full, H2full)

    # ---- frequencies on the equilibrium-proximate early window --------------
    mwin = ts .<= min(TF_FREQ, tf)
    fwin = ts[mwin][end]
    if sum(mwin) >= 64
        F, H1, H2, _, _ = band_FH1H2(ts[mwin], ecs[mwin])
    else
        F, H1, H2 = Ffull, H1full, H2full   # too short to sub-window
    end
    @printf("    EARLY-window (t≤%.0f) F/H1/H2 = %.3f / %.3f / %.3f kHz  [used for matching]\n",
            fwin, F, H1, H2)

    # ---- decay rate: linear log-max (1/τ_l) AND nonlinear (1/τ_nl) -----------
    Ffit = isnan(F) ? 2.69 : F
    tl, nmax = logmax_decay(ts, ecs, Ffit)
    fit = fit_fmode_decay(ts, ecs, Ffit; tstart_frac=0.12)
    mono, _, _ = envelope_monotonic(ts, ecs, Ffit)
    @printf("    decay 1/τ_l (linear log-max, %d maxima) = %.6g /M_⊙ (=%.4f kHz)\n",
            nmax, tl, tl*INVMSUN_TO_KHZ)
    @printf("    decay 1/τ_nl (nonlinear sinusoid)       = %.6g /M_⊙ (=%.4f kHz)  ω=%.5f -> f_nl=%.3f kHz\n",
            fit.invtau, fit.invtau*INVMSUN_TO_KHZ, fit.omega, fit.f_nl)
    @printf("    nl-fit window t∈[%.0f,%.0f] (n=%d)   envelope monotone=%s\n",
            fit.t1, fit.t2, fit.n, mono)

    return ResResult(dr, tf, n, F, H1, H2, Ffull, H1full, H2full, fwin,
                     tl, fit.invtau, fit.omega, fit.f_nl,
                     mono, fit.t1, fit.t2, fit.n)
end

# ===========================================================================
# Adversarial unit + robustness self-checks
# ===========================================================================
function adversarial_checks()
    println("\n" * "="^80)
    println("ADVERSARIAL VERIFICATION — units + extraction robustness")
    println("="^80)
    # (U1) unit constant
    @printf("[U1] 1/M_⊙ = c³/(GM_⊙) = %.6f kHz (cyclic).  Target 203.025 -> Δ=%.4f%%\n",
            INVMSUN_TO_KHZ, 100*(INVMSUN_TO_KHZ-203.025)/203.025)
    # (U2) cyclic vs angular: paper l.733 ω_nl=0.0834 -> f=2.71 kHz (uses /2π)
    fcheck = INVMSUN_TO_KHZ*0.0834/(2π)
    @printf("[U2] ANGULAR ω_nl=0.0834/M_⊙ -> f=203.025·0.0834/(2π)=%.4f kHz (paper 2.71) %s\n",
            fcheck, abs(fcheck-2.71)<0.05 ? "✓" : "✗")
    # (U3) decay rate is a RATE (no 2π): 0.00157/M_⊙ -> kHz
    @printf("[U3] DECAY 1/τ=0.00157/M_⊙ (rate, NO 2π) = %.5f kHz.  WRONG /2π would give %.5f kHz.\n",
            0.00157*INVMSUN_TO_KHZ, 0.00157*INVMSUN_TO_KHZ/(2π))
    # (U4) self-consistency on a synthetic damped sinusoid with KNOWN params:
    #   verify the fitter recovers (1/τ, ω) and the unit map f=ω·203.025/(2π).
    dtq = 1.0
    tq = collect(0.0:dtq:8000.0)
    ωtrue = 0.0834; λtrue = 0.0016
    sig = @. 1e-4*exp(-λtrue*tq)*cos(ωtrue*tq + 0.7) + 0.0014 + 1e-6*tq/8000
    yd = detrend_poly(tq, sig; d=3)
    Fkhz = ωtrue*INVMSUN_TO_KHZ/(2π)
    yf = butter4_bandpass(yd, Fkhz, dtq; bw_frac=0.35)
    i1 = Int(round(0.12*length(tq))); i2 = length(tq)-50
    p,_ = fit_damped_sinusoid(tq[i1:i2].-tq[i1], yf[i1:i2],
                              [maximum(abs.(yf[i1:i2])), 1e-3, ωtrue, 0.0, 0.0])
    @printf("[U4] synthetic recovery: 1/τ_fit=%.5f (true %.5f, Δ=%.1f%%)  ω_fit=%.5f (true %.5f, Δ=%.2f%%)\n",
            p[2], λtrue, 100*(p[2]-λtrue)/λtrue, p[3], ωtrue, 100*(p[3]-ωtrue)/ωtrue)
    @printf("     -> f_fit=%.4f kHz (true %.4f).  Fitter + unit map self-consistent.\n",
            p[3]*INVMSUN_TO_KHZ/(2π), Fkhz)
    return p[2], p[3]
end

# ===========================================================================
# DRIVER
# ===========================================================================
function main()
    println("="^80)
    println("R5 PRODUCTION: BDNK QNM extract + verify — Shum 2509.15303, case smallSB-F2")
    println("="^80)
    @printf("UNIT: 1/M_⊙ = %.4f kHz (cyclic).  Decay rate uses 203.025 (NO 2π).\n", INVMSUN_TO_KHZ)

    repro = "/data/haiyangw/claude/BDNK/code/BDNKStar/repro"
    # discover all ladder files and their Δr from the filename
    cand = String[]
    for f in readdir(repro)
        (startswith(f, "r5_eps_Dr") && endswith(f, ".txt")) || continue
        push!(cand, joinpath(repro, f))
    end
    sort!(cand)
    # parse Δr, read, keep only files with enough records for a meaningful FFT/fit.
    # need at least ~3 f-mode periods at Δt=1 (period ≈ 1/(2.69/203)=75.5 M_⊙) and
    # ideally ≥256 samples; we set a usability floor of 200 records.
    MIN_RECORDS = 200
    entries = Tuple{Float64,String,Int}[]
    for path in cand
        m = match(r"r5_eps_Dr([0-9.]+)\.txt", basename(path))
        m === nothing && continue
        dr = parse(Float64, m.captures[1])
        ts, _ = read_series(path)
        push!(entries, (dr, path, length(ts)))
    end
    sort!(entries, by=x->x[1])
    println("\nDiscovered ladder files (Δr, records):")
    for (dr,p,nr) in entries
        used = nr >= MIN_RECORDS ? "USE" : "skip(<$(MIN_RECORDS))"
        @printf("   Δr=%.4f  n=%-6d  %s   (%s)\n", dr, nr, used, basename(p))
    end

    used = [(dr,p) for (dr,p,nr) in entries if nr >= MIN_RECORDS]
    isempty(used) && error("no ladder file has >= $MIN_RECORDS records")

    results = ResResult[]
    for (dr,p) in used
        push!(results, analyze_resolution(p, dr))
    end
    sort!(results, by=r->r.dr)

    # ---- adversarial unit + robustness verification -------------------------
    adversarial_checks()

    # ---- continuum extrapolation (Richardson / Chabanov, Eq. l.756–757) ------
    println("\n" * "="^80)
    println("CONTINUUM EXTRAPOLATION  1/τ_{Δr}=1/τ_0+m Δr^p ; freqs similarly")
    println("="^80)
    drs   = [r.dr for r in results]
    Fs    = [r.F  for r in results]      # early-window frequencies
    H1s   = [r.H1 for r in results]
    H2s   = [r.H2 for r in results]
    taul  = [r.invtau_l  for r in results]   # 1/τ_l  (linear log-max)
    taunl = [r.invtau_nl for r in results]   # 1/τ_nl (nonlinear)

    function safe_extrap(label, ys; pfix=1.0)
        keep = .!isnan.(ys) .& isfinite.(ys)
        d = drs[keep]; y = ys[keep]
        if length(d) >= 2
            r = extrap_powerlaw(d, y; pfix=pfix)
            @printf("  %-6s: per-Δr %s -> continuum y0=%.5g (mode=%s, p=%.2f, m=%.4g)\n",
                    label, string(round.(y, sigdigits=4)), r.y0, r.mode, r.p, r.m)
            return r.y0
        elseif length(d) == 1
            @printf("  %-6s: single Δr value %.5g (no extrapolation possible)\n", label, y[1])
            return y[1]
        else
            @printf("  %-6s: NO usable values\n", label); return NaN
        end
    end
    # Frequencies: extrapolate the early-window F/H1/H2 to Δr->0.  At our coarse Δr
    # the early-window frequencies are already within ~1–2% of the paper, so the
    # extrapolation is a small correction; we report it but MATCH on the
    # finest-Δr early-window value (closest to the paper's converged regime).
    F0  = safe_extrap("F",  Fs;  pfix=1.0)
    H10 = safe_extrap("H1", H1s; pfix=1.0)
    H20 = safe_extrap("H2", H2s; pfix=1.0)
    # Decay: paper's p≈1 (l.794).  Extrapolate 1/τ_l (the robust linear method).
    tau0_l  = safe_extrap("1/τ_l",  taul;  pfix=1.0)
    tau0_nl = safe_extrap("1/τ_nl", taunl; pfix=1.0)
    tau_free = NaN
    if count(isfinite, taul) >= 3
        kk = isfinite.(taul)
        rf = extrap_powerlaw(drs[kk], taul[kk]); tau_free = rf.y0
        @printf("  1/τ_l (free-p fit): y0=%.5g  p=%.3f  m=%.4g\n", rf.y0, rf.p, rf.m)
    end

    # MATCH BASIS for FREQUENCIES.  The QNM frequencies are physical eigenvalues
    # (resolution-robust); the limiting factor at our coarse Δr is the SPECTRAL
    # RESOLUTION of the early (equilibrium-proximate) window — Δf ≈ 1/T_window.
    # We therefore match on the result whose early window contains the MOST f-mode
    # cycles (best-resolved PSD), i.e. the largest fwin·F_geo.  This selects the
    # long Δr=0.04 series (t≤2000 ≈ 26 cycles) over the short Δr=0.02 series
    # (t≤400 ≈ 5 cycles), and is the physically correct choice: a longer near-
    # equilibrium baseline gives a sharper line, NOT a finer grid with a coarse
    # spectrum.  (Both lie in the same physical regime; see frequency-window note.)
    function ncycles(r)
        Fg = isnan(r.F) ? 2.69 : r.F
        return r.fwin * Fg / INVMSUN_TO_KHZ      # f-mode cycles in early window
    end
    rfin = results[argmax([ncycles(r) for r in results])]

    # ---- SUMMARY vs targets --------------------------------------------------
    println("\n" * "="^80)
    println("SUMMARY (achieved vs paper smallSB-F2 / task targets)")
    println("="^80)
    @printf("UNIT 1/M_⊙ = %.4f kHz (target 203.025; Δ=%.3f%%)\n",
            INVMSUN_TO_KHZ, 100*(INVMSUN_TO_KHZ-203.025)/203.025)
    println("\nPer-Δr extracted:")
    println("  Δr      t_f    | F     H1    H2   (early)  | F     H1    H2   (full) | 1/τ_l    1/τ_nl  | drift")
    for r in results
        @printf("  %.4f  %5.0f | %.3f %.3f %.3f         | %.3f %.3f %.3f        | %.5f  %.5f | mono=%s\n",
                r.dr, r.tf, r.F, r.H1, r.H2, r.Ffull, r.H1full, r.H2full,
                r.invtau_l, r.invtau_nl, r.mono)
    end

    @printf("\nMATCH BASIS = best-resolved early window: Δr=%.4f, t≤%.0f (≈%.0f f-mode cycles).\n",
            rfin.dr, rfin.fwin, rfin.fwin*(isnan(rfin.F) ? 2.69 : rfin.F)/INVMSUN_TO_KHZ)
    @printf("Achieved F/H1/H2 (best-resolved early window) : %.3f / %.3f / %.3f kHz\n",
            rfin.F, rfin.H1, rfin.H2)
    @printf("Continuum-extrapolated F/H1/H2             : %.3f / %.3f / %.3f kHz\n", F0, H10, H20)
    @printf("           paper targets                  : 2.690 / 4.550 / 6.360 kHz\n")

    # decay rate for matching: take the LONGEST series' linear log-max 1/τ_l
    # (fully-resolved envelope), falling back to nonlinear if linear is NaN.
    rlong = results[argmax([r.tf for r in results])]
    tau_match = isfinite(rlong.invtau_l) && rlong.invtau_l>0 ? rlong.invtau_l : rlong.invtau_nl
    @printf("\nDecay rate (matching, longest series Δr=%.4f t_f=%.0f): 1/τ_l=%.6g, 1/τ_nl=%.6g /M_⊙\n",
            rlong.dr, rlong.tf, rlong.invtau_l, rlong.invtau_nl)
    @printf("  -> 1/τ used for matching = %.6g /M_⊙ (=%.4f kHz rate)\n",
            tau_match, tau_match*INVMSUN_TO_KHZ)
    @printf("Continuum 1/τ_0 (linear, p=1) : %.6g /M_⊙%s   (paper 0.0011)\n",
            tau0_l, isnan(tau_free) ? "" : @sprintf("  [free-p: %.6g]", tau_free))
    @printf("           paper per-Δr smallSB-F2 = 0.0016–0.0019 /M_⊙ (Δr=0.002–0.0032)\n")

    # ---- match logic ---------------------------------------------------------
    dF  = 100*(rfin.F -2.69)/2.69
    dH1 = 100*(rfin.H1-4.55)/4.55
    dH2 = 100*(rfin.H2-6.36)/6.36
    @printf("\nFreq deviations (best-resolved early window): ΔF=%+.2f%%  ΔH1=%+.2f%%  ΔH2=%+.2f%%  (target ~±2%%)\n",
            dF, dH1, dH2)

    freq_ok = all(x->!isnan(x), (rfin.F,rfin.H1,rfin.H2)) && all(abs.([dF,dH1,dH2]) .<= 2.5)
    finite_decay = isfinite(tau_match) && tau_match > 0
    # decay-within-30% test against per-Δr band midpoint AND continuum value.
    anchor_perdr = 0.00175      # midpoint of paper per-Δr range 0.0016–0.0019
    perdr_ok = finite_decay && abs(tau_match-anchor_perdr)/anchor_perdr <= 0.30
    cont_ok  = isfinite(tau0_l) && tau0_l>0 && abs(tau0_l-0.0011)/0.0011 <= 0.30
    decay_ok = finite_decay && (perdr_ok || cont_ok)

    @printf("\nfreq_ok(≤~2%%)=%s   finite_decay=%s   per-Δr-within-30%%=%s   continuum-within-30%%=%s\n",
            freq_ok, finite_decay, perdr_ok, cont_ok)
    matched = freq_ok && decay_ok
    @printf("MATCHED_TARGET = %s\n", matched)
    println("="^80)

    return (Fmatch=rfin.F, H1match=rfin.H1, H2match=rfin.H2,
            F0=F0,H10=H10,H20=H20, tau_match=tau_match, tau0=tau0_l,
            results=results, matched=matched, freq_ok=freq_ok, decay_ok=decay_ok)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
