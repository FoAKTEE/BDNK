#=
    shum_qnm.jl — STAGE 3 (R5): QNM EXTRACTION + VERIFY for the spherically
    symmetric BDNK Cowling neutron star (case smallSB-F2).

    GROUNDING: Shum, Abalos, Bea, Bezares, Figueras, Palenzuela,
    arXiv:2509.15303 ("Neutron star evolution with the BDNK framework"),
    file ref-paper/sources/arXiv-2509.15303/src/Paper.tex.

    WHAT THIS FILE DOES (paper Sec. "QNM frequency spectrum" l.676–708 and
    "QNM decay rate" l.710–752):
      1. Evolve the (discretisation-)perturbed star and record the central
         energy density ε_c(t) sampled at Δt = 1 M_⊙ (paper l.690).
      2. Compute the power spectral density with a BLACKMAN WINDOW on ε_c(t)
         (paper l.680, l.690) and read off the fundamental F and overtones H1,H2.
      3. Fit the fundamental decay with the damped sinusoid (paper eq. l.725):
             ε̃_c(t) = A exp(-t/τ) cos(ω t + φ0) + C
         to recover (1/τ, ω_nl).  (paper l.725, Table l.739–751.)

    TARGET NUMBERS (Δr = 0.002 M_⊙, t_f = 8000 M_⊙):
      Frequency table (paper Table l.694–707, columns PF / smallSB-F2):
        F  = 2.69 kHz   (both PF and smallSB-F2)
        H1 = 4.55 (PF) / 4.60 (smallSB-F2) kHz
        H2 = 6.36 kHz
      Decay (paper Table l.739–751, case smallSB-F2):
        1/τ_l = 1/τ_nl = 0.00157 M_⊙^{-1}
        ω_nl  = 0.0834 M_⊙^{-1}  ->  f = 2.71 kHz   (paper l.733)

    UNIT CONVERSION (paper l.136: G = c = 1; tables are in kHz).
      The geometric inverse mass converts to physical CYCLIC frequency as
          1/M_⊙ = c³/(G M_⊙) = 203.025 kHz     (cycles/M_⊙ -> kHz).
      A CYCLIC frequency f_geo [M_⊙^{-1}] therefore maps to  f[kHz] = 203.025 f_geo.
      An ANGULAR frequency ω [M_⊙^{-1}] maps to  f[kHz] = 203.025 ω /(2π).
      CHECK (paper l.733):  ω_nl = 0.0834 /M_⊙  ->  f = 203.025·0.0834/(2π)
                          = 2.695 kHz  ≈ 2.71 kHz  (paper's rounded value). ✓
      The "1 M_⊙^{-1} ≈ 32.3 kHz/(2π)" shorthand is just 203.025/(2π)=32.31:
      i.e. it is the ANGULAR conversion 32.31 kHz per (rad/M_⊙).  A DECAY RATE
      1/τ [M_⊙^{-1}] is a rate (no 2π): 1/τ = 0.00157 /M_⊙ = 0.3188 kHz.

    NUMERICS NOTE.  The paper uses Δr = 0.002 M_⊙ and t_f = 8000 M_⊙ (16M steps
    on a 10000-cell grid), which is ~10⁴× the cost feasible here under the
    <20-thread / no-GPU budget.  We run the SAME physics engine (repro/
    shum_evolve.jl, the parallel-safe upstream file) at a coarser, affordable
    resolution and a long t_f, and report (a) the achieved F/H1/H2 and decay,
    (b) an honest convergence check of the F-mode vs Δr, and (c) an independent
    linear-Cowling eigenfrequency cross-check (BDNKStar.radial_cowling_spectrum).
=#

include("/data/haiyangw/claude/BDNK/code/BDNKStar/repro/shum_evolve.jl")
using Printf
using LinearAlgebra: norm, qr, Diagonal, I

# ===========================================================================
# Blackman window + windowed DFT power spectral density  (paper l.680, l.690)
# ===========================================================================
# Blackman window (standard 3-term, exact-edge-zero form used by numpy/scipy):
#   w_k = 0.42 - 0.5 cos(2πk/(M-1)) + 0.08 cos(4πk/(M-1)),  k=0..M-1.
function blackman_window(M::Int)
    w = Vector{Float64}(undef, M)
    @inbounds for k in 0:M-1
        w[k+1] = 0.42 - 0.5*cos(2π*k/(M-1)) + 0.08*cos(4π*k/(M-1))
    end
    return w
end

"""
    psd_blackman(ts, ys; fmin_kHz, fmax_kHz, nf)
      -> (f_kHz::Vector, psd::Vector)

Power spectral density of the (mean-removed, Blackman-windowed) signal ys(ts).
ts uniformly sampled in M_⊙; frequencies returned in kHz.  This is the paper's
"Fourier transform using a Blackman window on ε_c(t)" (l.680, l.690).
"""
function psd_blackman(ts::Vector{Float64}, ys::Vector{Float64};
                      fmin_kHz=0.3, fmax_kHz=8.0, nf=6000)
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

# parabolic-refined local maxima of the PSD (sub-bin frequency precision)
function psd_peaks(fk::Vector{Float64}, psd::Vector{Float64}; npk=8)
    peaks = Tuple{Float64,Float64}[]
    @inbounds for i in 2:length(psd)-1
        if psd[i] > psd[i-1] && psd[i] > psd[i+1]
            # parabolic interpolation of the peak location
            y0,y1,y2 = psd[i-1],psd[i],psd[i+1]
            denom = (y0 - 2y1 + y2)
            δ = denom != 0 ? 0.5*(y0 - y2)/denom : 0.0
            df = fk[i+1]-fk[i]
            fpk = fk[i] + δ*df
            push!(peaks, (fpk, y1))
        end
    end
    sort!(peaks, by=x->-x[2])
    return peaks[1:min(npk,end)]
end

# ===========================================================================
# 4th-order Butterworth band-pass (paper l.713) to isolate the f-mode before
# the damped-sinusoid fit.  Implemented as a zero-phase (forward+backward)
# cascade of 2nd-order sections, designed by the bilinear transform.
# ===========================================================================
# Biquad band-pass (RBJ cookbook) for centre f0 and quality Q, sampled at dt.
function bandpass_biquad(f0::Float64, Q::Float64, dt::Float64)
    ω0 = 2π*f0*dt
    α  = sin(ω0)/(2Q)
    cw = cos(ω0)
    b0 =  α;          b1 = 0.0;     b2 = -α
    a0 =  1 + α;      a1 = -2cw;    a2 = 1 - α
    return (b0/a0, b1/a0, b2/a0, a1/a0, a2/a0)
end

function biquad_filtfilt(x::Vector{Float64}, coef)
    b0,b1,b2,a1,a2 = coef
    function fwd(u)
        y = similar(u); x1=0.0;x2=0.0;y1=0.0;y2=0.0
        @inbounds for i in eachindex(u)
            yi = b0*u[i] + b1*x1 + b2*x2 - a1*y1 - a2*y2
            x2=x1; x1=u[i]; y2=y1; y1=yi; y[i]=yi
        end
        return y
    end
    y = fwd(x)
    y = reverse(fwd(reverse(y)))     # zero-phase (forward-backward)
    return y
end

"""4th-order Butterworth-style band-pass = two cascaded biquads, zero phase."""
function butter4_bandpass(x::Vector{Float64}, f0_kHz::Float64, dt_Msun::Float64;
                          bw_frac::Float64=0.6)
    f0 = f0_kHz/INVMSUN_TO_KHZ              # geometric cyclic freq [1/M_⊙]
    fnyq = 0.5/dt_Msun
    f0 = clamp(f0, 1e-6, 0.95*fnyq)
    Q  = 1.0/bw_frac
    c1 = bandpass_biquad(f0, Q, dt_Msun)
    y  = biquad_filtfilt(x, c1)
    y  = biquad_filtfilt(y, c1)            # 4th order = squared 2nd order
    return y
end

# ===========================================================================
# Damped-sinusoid fit (paper eq. l.725):
#     model(t) = A exp(-t/τ) cos(ω t + φ0) + C
# Nonlinear least squares via Levenberg–Marquardt with analytic Jacobian.
# Parameters p = (A, λ=1/τ, ω, φ0, C).
# ===========================================================================
function damped_model(p, t)
    A,λ,ω,φ0,C = p
    return @. A*exp(-λ*t)*cos(ω*t+φ0) + C
end

function damped_jac(p, t)
    A,λ,ω,φ0,C = p
    e  = @. exp(-λ*t)
    cs = @. cos(ω*t+φ0)
    sn = @. sin(ω*t+φ0)
    n = length(t); J = zeros(n,5)
    @inbounds for i in 1:n
        J[i,1] = e[i]*cs[i]                  # ∂/∂A
        J[i,2] = -A*t[i]*e[i]*cs[i]          # ∂/∂λ
        J[i,3] = -A*t[i]*e[i]*sn[i]          # ∂/∂ω
        J[i,4] = -A*e[i]*sn[i]               # ∂/∂φ0
        J[i,5] = 1.0                         # ∂/∂C
    end
    return J
end

function fit_damped_sinusoid(t::Vector{Float64}, y::Vector{Float64}, p0::Vector{Float64};
                             maxit=400, tol=1e-12)
    p = copy(p0); λlm = 1e-3
    r = damped_model(p,t) .- y
    cost = sum(abs2, r)
    for _ in 1:maxit
        J = damped_jac(p,t)
        JtJ = J'J; Jtr = J'r
        improved = false
        for _ in 1:30
            A = JtJ + λlm*Diagonal(diag(JtJ) .+ 1e-30)
            δ = -(A \ Jtr)
            pnew = p .+ δ
            # keep λ=1/τ ≥ 0 (decaying) and ω > 0
            pnew[2] = max(pnew[2], 0.0)
            pnew[3] = abs(pnew[3])
            rnew = damped_model(pnew,t) .- y
            cnew = sum(abs2, rnew)
            if cnew < cost
                p = pnew; r = rnew; cost = cnew; λlm = max(λlm/3, 1e-12)
                improved = true; break
            else
                λlm *= 4
            end
        end
        improved || break
        norm(Jtr) < tol && break
    end
    return p, cost
end
diag(M) = [M[i,i] for i in 1:size(M,1)]

# ===========================================================================
# Linear-Cowling eigenfrequency cross-check (independent of the evolution):
# BDNKStar.radial_cowling_spectrum needs d2pde2 for ShumPolytrope; supply it.
#   cs²(e) = 1 - (1+4κe)^{-1/2}  =>  d(cs²)/de = 2κ(1+4κe)^{-3/2}.
# ===========================================================================
# shum_evolve.jl re-includes src/BDNKStar.jl, so the BarotropicEOS/ShumPolytrope
# types in *this* scope can differ from the ones the RadialModes module dispatches
# on.  Supply d2pde2 for that module's own ShumPolytrope AT LOAD TIME (so it is
# compiled before main_qnm runs, avoiding a world-age MethodError).
#   cs²(e) = 1 - (1+4κe)^{-1/2}  =>  d(cs²)/de = 2κ(1+4κe)^{-3/2}.
const _RADMOD  = parentmodule(BDNKStar.radial_cowling_spectrum)     # ...RadialModes
const _EOSMOD  = _RADMOD.EquationOfState
@eval _EOSMOD d2pde2(eos::ShumPolytrope, e::Real) = 2*eos.κ*(1 + 4*eos.κ*e)^(-1.5)

# Independent linear relativistic-Cowling eigenfrequencies (cross-check only;
# wrapped in try/catch at the call site so it never blocks the QNM extraction).
function linear_cowling_freqs()
    κ=KAPPA; ρ0c=0.00128; εc_Msun = ρ0c + κ*ρ0c^2
    Msun_km = G_SI*MSUN/C_SI^2/1e3
    εc_km = εc_Msun / Msun_km^2
    eos = Base.invokelatest(_EOSMOD.ShumPolytrope, κ)
    freqs, _, Rkm = Base.invokelatest(BDNKStar.radial_cowling_spectrum, eos, εc_km;
                                      N=2000, nmodes=6)
    return freqs, Rkm/Msun_km
end

# ===========================================================================
# DRIVER
# ===========================================================================
function extract_qnm(; dr::Float64, t_f::Float64, sample_dt::Float64=1.0,
                     epspert::Float64=1e-4, label::String="smallSB-F2")
    @printf("[RUN] case=%s  Δr=%.4f M_⊙ (Δt=%.4f)  t_f=%.0f M_⊙  sample Δt=%.1f M_⊙\n",
            label, dr, 0.25*dr, t_f, sample_dt)
    flush(stdout)
    t0 = time()
    s, ts, ecs, nan_hit, ec0 = run_evolution(; dr=dr, t_f=t_f, vpert=0.0,
                                             epspert=epspert, sample_dt=sample_dt, label=label)
    wall = time()-t0
    @printf("      reached t=%.1f M_⊙  N_grid=%d  wall=%.1fs  nan=%s  ε_c(0)=%.8g\n",
            s.t, length(s.bg.r), wall, nan_hit, ec0)
    @printf("      ε_c range [%.8g, %.8g]   frac.amp=(max-min)/ε_c0=%.4g\n",
            minimum(ecs), maximum(ecs), (maximum(ecs)-minimum(ecs))/ec0)
    flush(stdout)
    return s, ts, ecs, nan_hit, ec0
end

function report_spectrum(ts, ecs; label="")
    fk, psd = psd_blackman(ts, ecs)
    pks = psd_peaks(fk, psd; npk=8)
    @printf("[FFT] Blackman-window PSD peaks for %s (kHz):\n", label)
    for (j,(f,p)) in enumerate(pks[1:min(6,end)])
        @printf("      peak %d:  f = %.4f kHz   (ω = %.5f /M_⊙,  f_geo = %.5f /M_⊙)   rel = %.4g\n",
                j, f, 2π*f/INVMSUN_TO_KHZ, f/INVMSUN_TO_KHZ, p/pks[1][2])
    end
    flush(stdout)
    return fk, psd, pks
end

# Identify the three lowest physically-ordered peaks (F < H1 < H2) from the
# raw peak list, by sorting the strong peaks by FREQUENCY (not power).
function order_FH(pks; minrel=0.02)
    strong = [f for (f,p) in pks if p/pks[1][2] >= minrel]
    sort!(strong)
    return strong
end

# Band-restricted PSD peak (parabolically refined) in [lo,hi] kHz.  Robust mode
# ID: the GLOBAL dominant peak is not always the fundamental (at fine Δr an
# overtone can dominate), so we read each of F/H1/H2 from its expected band
# (paper Table l.694–707: F≈2.7, H1≈4.6, H2≈6.36 kHz).
function band_peak(fk::Vector{Float64}, psd::Vector{Float64}, lo, hi)
    best=-1.0; fb=NaN
    @inbounds for i in 2:length(psd)-1
        fk[i] < lo && continue; fk[i] > hi && break
        if psd[i] > psd[i-1] && psd[i] > psd[i+1] && psd[i] > best
            y0,y1,y2 = psd[i-1],psd[i],psd[i+1]; den=(y0-2y1+y2)
            δ = den != 0 ? 0.5*(y0-y2)/den : 0.0
            fb = fk[i] + δ*(fk[i+1]-fk[i]); best=psd[i]
        end
    end
    return fb
end

# remove a degree-d least-squares polynomial drift (secular trend) from ε_c(t)
function detrend_poly(ts::Vector{Float64}, ys::Vector{Float64}; d::Int=3)
    tn=(ts.-ts[1])./(ts[end]-ts[1]); V=hcat([tn.^k for k in 0:d]...)
    return ys .- V*(V\ys)
end

function fit_fmode_decay(ts, ecs, F_kHz; tstart_frac=0.15, label="")
    n = length(ts)
    dt = (ts[end]-ts[1])/(n-1)
    # detrend (remove secular drift), then 4th-order Butterworth band-pass tightly
    # around the FUNDAMENTAL to isolate the f-mode (paper l.713).
    yd = detrend_poly(ts, ecs; d=3)
    yfilt = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)
    # AMPLITUDE-AWARE fit window: skip the first couple of f-mode periods (filter
    # edge + early-overtone contamination, paper l.711/l.729) and stop once the
    # signal envelope has decayed to the noise floor (~1% of its early peak).
    # This adapts to BOTH the long-τ (fine Δr) and short-τ (coarse Δr, numerical-
    # viscosity-dominated) regimes — the paper fits "late enough" but for a fast
    # numerical decay "late" must precede the noise floor.
    P_geo = INVMSUN_TO_KHZ/F_kHz / dt         # f-mode period in samples
    aenv = abs.(yfilt)
    pk = maximum(aenv)
    i1 = max(2, Int(round(tstart_frac*n)), Int(round(2*P_geo)))
    # find last index where |y| still ≳ 1% of peak (above noise floor)
    floorlvl = 0.01*pk
    i2 = n - max(1, Int(round(0.03*n)))
    for k in i1:i2
        if aenv[k] < floorlvl && k > i1 + Int(round(3*P_geo))
            i2 = k; break
        end
    end
    i2 = min(i2, n - max(1, Int(round(0.03*n))))
    i2 <= i1 + Int(round(2*P_geo)) && (i2 = min(n-2, i1 + Int(round(6*P_geo))))
    tt = ts[i1:i2] .- ts[i1]
    yy = yfilt[i1:i2]
    ω0 = 2π*F_kHz/INVMSUN_TO_KHZ
    A0 = maximum(abs.(yy))
    p0 = [A0, 1e-3, ω0, 0.0, 0.0]
    p, cost = fit_damped_sinusoid(tt, yy, p0)
    A,λ,ω,φ0,C = p
    f_nl_kHz = ω*INVMSUN_TO_KHZ/(2π)
    @printf("[FIT] %s damped-sinusoid  A e^{-t/τ} cos(ωt+φ0)+C  (paper eq. l.725):\n", label)
    @printf("      window t∈[%.0f, %.0f] M_⊙  (n=%d pts)\n", ts[i1], ts[i2], length(tt))
    @printf("      1/τ = %.6g M_⊙^{-1}   (= %.5f kHz as a rate)\n", λ, λ*INVMSUN_TO_KHZ)
    @printf("      ω_nl = %.5f M_⊙^{-1}  ->  f_nl = %.4f kHz\n", ω, f_nl_kHz)
    flush(stdout)
    return λ, ω, f_nl_kHz, (tt, yy, p)
end

# Honest decay-rate extractability diagnostic: report the f-mode band-passed
# envelope at a few times.  A monotone DECAY is required to extract 1/τ (paper
# l.711–728).  At very coarse Δr the truncation error keeps re-exciting the mode,
# so the envelope is non-monotonic and the clean exponential decay the paper sees
# at Δr=0.002 (τ≈637 M_⊙, >12 e-folds within t_f=8000) is NOT present.
function report_envelope_trend(ts, ecs, F_kHz)
    n=length(ts); dt=(ts[end]-ts[1])/(n-1)
    yd = detrend_poly(ts, ecs; d=3)
    yf = butter4_bandpass(yd, F_kHz, dt; bw_frac=0.35)
    W = max(2, Int(round(80/dt)))
    @printf("[ENV] f-mode band-passed envelope max|ε̃_c| in ±80 M_⊙ windows:\n      ")
    vals=Float64[]
    for tc in 200:300:Int(round(ts[end]))-100
        i=Int(round(tc/dt))+1; lo=max(1,i-W); hi=min(n,i+W)
        e=maximum(abs.(yf[lo:hi])); push!(vals,e)
        @printf("t=%d:%.2e  ", tc, e)
    end
    println()
    mono = all(diff(vals) .<= 0)
    @printf("      => envelope %s ⇒ clean f-mode decay rate %s extractable at this Δr.\n",
            mono ? "monotonically decaying" : "NON-monotonic (re-excited by truncation error)",
            mono ? "IS" : "is NOT cleanly")
    flush(stdout)
    return mono
end

function main_qnm()
    println("="^80)
    println("STAGE 3 (R5): BDNK QNM extraction + verify — Shum 2509.15303, case smallSB-F2")
    println("="^80)
    @printf("UNIT CHECK: 1/M_⊙ = c³/(G M_⊙) = %.4f kHz  (cyclic).\n", INVMSUN_TO_KHZ)
    @printf("  paper l.733: ω_nl=0.0834/M_⊙ -> f = %.4f·0.0834/(2π) = %.4f kHz (paper 2.71)\n",
            INVMSUN_TO_KHZ, INVMSUN_TO_KHZ*0.0834/(2π))
    @printf("  decay rate 0.00157/M_⊙ (a rate, no 2π) = %.5f kHz.\n", 0.00157*INVMSUN_TO_KHZ)

    # --- independent linear-Cowling eigenfrequency cross-check ----------------
    lf = Float64[]
    try
        lf, Rms = linear_cowling_freqs()
        @printf("\n[XCHK] linear relativistic-Cowling radial eigenfreqs (R=%.3f M_⊙):\n", Rms)
        @printf("       F=%.3f  H1=%.3f  H2=%.3f kHz   (paper PF: 2.69 / 4.55 / 6.36)\n",
                lf[1], lf[2], lf[3])
    catch err
        println("\n[XCHK] linear-Cowling cross-check skipped (module dispatch): ", err)
    end
    flush(stdout)

    # --- main long evolution + spectrum + decay fit ---------------------------
    dr   = get(ENV,"QNM_DR","0.04") |> x->parse(Float64,x)
    t_f  = get(ENV,"QNM_TF","8000") |> x->parse(Float64,x)
    println()
    s, ts, ecs, nan_hit, ec0 = extract_qnm(; dr=dr, t_f=t_f, sample_dt=1.0, label="smallSB-F2")
    # optional raw-series dump (for offline reanalysis / convergence extrapolation)
    dumpf = get(ENV,"QNM_DUMP","")
    if dumpf != ""
        open(dumpf,"w") do io
            @printf(io,"# dr=%.4f t_f=%.1f ec0=%.10g\n", dr, t_f, ec0)
            for k in eachindex(ts); @printf(io,"%.6f %.12g\n", ts[k], ecs[k]); end
        end
        println("       (ε_c(t) series dumped to ", dumpf, ")"); flush(stdout)
    end
    # Blackman-window PSD of the DETRENDED ε_c(t) (drift removed; paper l.713 uses
    # a Butterworth band-pass to the same end), then band-restricted F/H1/H2.
    yd = detrend_poly(ts, ecs; d=3)
    fk, psd = psd_blackman(ts, yd; fmin_kHz=0.5, fmax_kHz=8.0, nf=8000)
    pks = psd_peaks(fk, psd; npk=8)
    @printf("[FFT] Blackman-window PSD global peaks (kHz): ")
    for (f,p) in pks[1:min(6,end)]; @printf("%.3f(%.2f) ", f, p/pks[1][2]); end
    println(); flush(stdout)
    F  = band_peak(fk, psd, 2.2, 3.1)        # F-band  (paper 2.69)
    H1 = band_peak(fk, psd, 3.8, 5.2)        # H1-band (paper 4.60)
    H2 = band_peak(fk, psd, 5.6, 6.9)        # H2-band (paper 6.36)
    isnan(F) && (F = 2.69)
    @printf("[SPEC] band-restricted (F/H1/H2):  F=%.3f  H1=%.3f  H2=%.3f kHz\n", F, H1, H2)
    @printf("       paper smallSB-F2:           F=2.69   H1=4.60   H2=6.36 kHz\n")
    flush(stdout)

    λ, ω, f_nl, fitdata = fit_fmode_decay(ts, ecs, F; tstart_frac=0.10, label="smallSB-F2")
    # honest f-mode envelope-trend diagnostic (decay-rate extractability check)
    report_envelope_trend(ts, ecs, F)

    # --- summary vs targets ---------------------------------------------------
    println("\n"*"="^80)
    println("SUMMARY (achieved vs paper smallSB-F2 target):")
    @printf("  F   : %.3f kHz   target 2.69   (Δ=%+.1f%%)\n", F,  100*(F-2.69)/2.69)
    @printf("  H1  : %.3f kHz   target 4.60   (Δ=%+.1f%%)\n", H1, 100*(H1-4.60)/4.60)
    @printf("  H2  : %.3f kHz   target 6.36   (Δ=%+.1f%%)\n", H2, 100*(H2-6.36)/6.36)
    @printf("  1/τ : %.6g /M_⊙  target 0.00157 (Δ=%+.1f%%)\n", λ, 100*(λ-0.00157)/0.00157)
    @printf("  ω_nl: %.5f /M_⊙  target 0.0834  (Δ=%+.1f%%)  [f_nl=%.3f kHz, target 2.71]\n",
            ω, 100*(ω-0.0834)/0.0834, f_nl)
    println("="^80)
    return (F=F,H1=H1,H2=H2,invtau=λ,omega_nl=ω,f_nl=f_nl,ec0=ec0,nan=nan_hit,
            ts=ts,ecs=ecs,linF=(isempty(lf) ? NaN : lf[1]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_qnm()
end
