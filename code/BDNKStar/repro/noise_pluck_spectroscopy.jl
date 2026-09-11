#=
    noise_pluck_spectroscopy.jl  —  THE EXPERIMENT.

    Initialize a VALID + CAUSAL band-limited Gaussian-noise perturbation on the
    causal Shum benchmark star and spectroscope WHICH stellar modes a generic
    random pluck excites, using the full-frame (causal) SphBDNK engine.

    Cross-check: time-domain ℓ-projected power spectra  vs  the NonRadialModes
    Cowling eigensolver (ℓ=2,4 f/p) and the ℓ=0 time-domain peak.

    Pipeline:
      1. CAUSAL star: ShumPolytrope(100), εc=ρ0c+100ρ0c² (M=1.4); production frame
         c₊=√3·cs subluminal everywhere on this soft star.
      2. NOISE ID: band-limited GRF on radial Fourier × Legendre P_ℓ(cosθ),
         ℓ=0..6, Gaussian spectral cutoff → power on RESOLVED scales. Each field
         component projected onto its correct centre/pole PARITY.
      3. VALID: ε0+δρ>0, velocities subluminal, recoverable (rescale A).
      4. CAUSALITY (Shum route): shum_frame_wellposed(1,1,0.999)=true and
         c₊=√3·cs(r)≤1 at every radius; report margin.
      5. EVOLVE ideal (η̂=0) long (T~1500), sample ℓ=0,2,4 projected moments + probes.
      6. FFT → power spectra; overlay eigenfreqs; report peaks + match%.
      7. EVOLVE with causal viscosity (η̂=0.05); report spectral change + energy.
=#
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar: SphBDNK, SphBackground
using Printf, Statistics, Random

const EOS  = ShumPolytrope(100.0)
const ρ0c  = 0.00128
const εc   = ρ0c + 100*ρ0c^2
kHz2cyc(f) = f*BDNKStar.Units.Msun_to_km*BDNKStar.Units.kHz_to_km   # kHz -> cyclic geom freq

# Legendre P_ℓ(x), ℓ=0..6
function Pl(l, x)
    l==0 && return 1.0
    l==1 && return x
    l==2 && return (3x^2-1)/2
    l==3 && return (5x^3-3x)/2
    l==4 && return (35x^4-30x^2+3)/8
    l==5 && return (63x^5-70x^3+15x)/8
    l==6 && return (231x^6-315x^4+105x^2-5)/16
    return 0.0
end

# ============================================================================
# 1. CAUSALITY CHECK (Shum route) — authoritative, BEFORE evolving
# ============================================================================
function causality_check(s)
    wp = shum_frame_wellposed(1.0, 1.0, 0.999)        # 0<q̂<ŝ
    cpmax = 0.0; rmax = 0.0; csmax = 0.0; nbad = 0; margin = Inf
    for i in 1:s.grid.Nr
        ε = s.ε0[i]; ε ≤ 0 && continue
        cs = sqrt(max(sound_speed2(EOS, ε), 0.0))
        _, cp, _ = shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, cs)   # c₊=√3 cs
        cp > cpmax && (cpmax = cp; rmax = s.grid.r[i]; csmax = cs)
        cp > 1.0 && (nbad += 1)
        margin = min(margin, 1.0 - cp)
    end
    @printf("CAUSALITY (Shum route, production frame ŝ,â,q̂=1,1,0.999):\n")
    @printf("   well-posed (0<q̂<ŝ)              : %s\n", wp)
    @printf("   max c₊=√3·cs = %.5f  at r/R=%.3f (cs=%.4f)\n", cpmax, rmax/s.R, csmax)
    @printf("   superluminal radii (c₊>1)       : %d / %d\n", nbad, s.grid.Nr)
    @printf("   causality margin  min(1−c₊)     : %.5f\n", margin)
    causal = wp && nbad == 0
    @printf("   => CAUSAL + WELL-POSED          : %s\n\n", causal)
    return (causal=causal, wp=wp, margin=margin, cpmax=cpmax)
end

# ============================================================================
# 2. BAND-LIMITED GAUSSIAN-NOISE SEED, parity-projected per component
# ============================================================================
# Build each component as Σ_{k,l} g_kl * c_kl * radial_k(r) * Pl(l,cosθ),
# with c_kl ~ N(0,1), g_kl a smooth Gaussian spectral cutoff on (k,l), and
# radial_k chosen to enforce the correct CENTRE parity per component:
#   δρ  : even at centre (kept finite), even in θ  -> cos-like radial, ℓ even or all
#   δvr : ODD  at centre (∝ r near 0), even in θ   -> sin-like radial
#   δvθ : even at centre, ODD in θ                 -> use ℓ odd Legendre? No:
#         δvθ pole-parity is ODD (antisym about pole). We realise the θ-structure
#         via dPl/dθ-like factor sinθ*P'_ℓ which is naturally odd about poles.
# We respect pole parity by component:
#   δρ even-θ, δvr even-θ : use P_ℓ(cosθ)            (even about θ=0,π for ℓ even/all)
#   δvθ odd-θ             : use sinθ * dPl/dx (cosθ)  (vanishes & flips at poles)
function dPl(l, x)   # d/dx P_l
    l==0 && return 0.0
    l==1 && return 1.0
    l==2 && return 3x
    l==3 && return (15x^2-3)/2
    l==4 && return (35x^3 - 15x)/2      # d/dx (35x^4-30x^2+3)/8 = (140x^3-60x)/8
    l==5 && return (315x^4-210x^2+15)/8
    l==6 && return (1386x^5-1260x^3+210x)/16
    return 0.0
end

const RNGSEED = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 20260620

function seed_noise!(st, s; A=1e-4, kmax=6, lmax=6, kcut=4.0, lcut=4.0, rng=MersenneTwister(RNGSEED))
    Nr=s.grid.Nr; Nθ=s.grid.Nθ; R=s.R
    # random coefficients per (component,k,l), with Gaussian spectral envelope
    g(k,l) = exp(-0.5*((k/kcut)^2 + (l/lcut)^2))
    # draw coefficient tensors
    cρ  = randn(rng, kmax, lmax+1)
    cvr = randn(rng, kmax, lmax+1)
    cvθ = randn(rng, kmax, lmax+1)
    # zero the fields' interior
    fill!(st.δρ,0.0); fill!(st.δSr,0.0); fill!(st.δSθ,0.0); fill!(st.δvr,0.0); fill!(st.δvθ,0.0)
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        r = s.grid.r[ii]; x = s.grid.cosθ[jj]; ξ = r/R
        dρ = 0.0; dvr = 0.0; dvθ = 0.0
        for k in 1:kmax
            # CENTRE-parity radial basis:
            #   even-at-centre: cos(kπξ/?)  -> use (1-ξ) tapered even cos: cos((k-1/2)πξ)?
            #   simplest even basis vanishing-derivative at 0: cos(kπξ)  (even in r)
            #   odd-at-centre  : sin(kπξ)   (∝ r near 0)
            re_even = cos((k-1)*π*ξ)        # k=1 -> 1 (DC), even about r=0
            re_odd  = sin(k*π*ξ)            # ∝ r near 0, odd about r=0; vanishes at surface
            for l in 0:lmax
                w = g(k,l)
                Plx = Pl(l, x)
                # δρ : even-centre, even-θ
                dρ  += w*cρ[k,l+1]  * re_even * Plx
                # δvr: odd-centre,  even-θ
                dvr += w*cvr[k,l+1] * re_odd  * Plx
                # δvθ: even-centre, odd-θ (sinθ * P'_l(cosθ))
                dvθ += w*cvθ[k,l+1] * re_even * (s.grid.sinθ[jj]*dPl(l,x))
            end
        end
        st.δρ[ii+1,jj+1]  = dρ
        st.δvr[ii+1,jj+1] = dvr
        st.δvθ[ii+1,jj+1] = dvθ
    end
    # ---- normalise component RMS to comparable scale, then enforce VALID ----
    # scale δρ to be a small fraction of ε0; scale velocities to a small value.
    ρrms = sqrt(mean(st.δρ[2:Nr+1,2:Nθ+1].^2))
    vrrms = sqrt(mean(st.δvr[2:Nr+1,2:Nθ+1].^2))
    vθrms = sqrt(mean(st.δvθ[2:Nr+1,2:Nθ+1].^2))
    ε0rms = sqrt(mean(s.ε0.^2))
    sρ  = A*ε0rms/max(ρrms,1e-300)
    sv  = A/max(max(vrrms,vθrms),1e-300)
    st.δρ  .*= sρ
    st.δvr .*= sv
    st.δvθ .*= sv
    # initial momentum consistent with velocity: δS_i ≈ w₀ δv_i (lowered)
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        w0 = s.ε0[ii]+s.p0[ii]
        st.δSr[ii+1,jj+1] = w0*s.eΛ[ii]*st.δvr[ii+1,jj+1]          # δS_r = w₀ g_rr δv^r
        st.δSθ[ii+1,jj+1] = w0*s.grid.r[ii]^2*st.δvθ[ii+1,jj+1]    # δS_θ = w₀ g_θθ δv^θ
    end
    # ---- enforce ε0+δρ>0 and subluminal: rescale globally if violated ----
    minfac = 1.0
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        # δε = h0 δρ ; require ε0 + δε > 0 -> here δρ already energy-density-scaled
        h0 = (s.ε0[ii]+s.p0[ii])/s.ρ0[ii]
        δε = h0*st.δρ[ii+1,jj+1]
        if s.ε0[ii] + δε ≤ 0
            f = 0.5*s.ε0[ii]/abs(δε)
            minfac = min(minfac, f)
        end
        vmag = sqrt(s.eΛ[ii]*st.δvr[ii+1,jj+1]^2 + s.grid.r[ii]^2*st.δvθ[ii+1,jj+1]^2)
        if vmag ≥ 0.5
            minfac = min(minfac, 0.4/vmag)
        end
    end
    if minfac < 1.0
        for F in (st.δρ,st.δSr,st.δSθ,st.δvr,st.δvθ); F .*= minfac; end
    end
    return (sρ=sρ, sv=sv, minfac=minfac)
end

# validity report
function validity_report(st, s)
    Nr=s.grid.Nr; Nθ=s.grid.Nθ
    minε = Inf; maxv = 0.0
    @inbounds for jj in 1:Nθ, ii in 1:Nr
        h0 = (s.ε0[ii]+s.p0[ii])/s.ρ0[ii]
        δε = h0*st.δρ[ii+1,jj+1]
        minε = min(minε, (s.ε0[ii]+δε)/s.ε0[ii])
        vmag = sqrt(s.eΛ[ii]*st.δvr[ii+1,jj+1]^2 + s.grid.r[ii]^2*st.δvθ[ii+1,jj+1]^2)
        maxv = max(maxv, vmag)
    end
    @printf("VALID constraints:\n")
    @printf("   min (ε0+δε)/ε0  (>0 ⇒ positivity) : %.6f\n", minε)
    @printf("   max |δv|        (<1 ⇒ subluminal)  : %.3e\n", maxv)
    @printf("   positivity OK   : %s ;  subluminal OK : %s\n\n", minε>0, maxv<1)
    return (minε=minε, maxv=maxv)
end

# ============================================================================
# 3. ℓ-projected moments
# ============================================================================
# M_ℓ(t) = ∫ δρ * P_ℓ(cosθ) * √γ sinθ dr dθ  (over the star)
function lproj(st, s, l)
    acc = 0.0
    @inbounds for jj in 1:s.grid.Nθ, ii in 1:s.grid.Nr
        acc += st.δρ[ii+1,jj+1]*Pl(l, s.grid.cosθ[jj])*s.sqγr[ii]*s.grid.sinθ[jj]
    end
    acc*s.grid.dr*s.grid.dθ
end

# evolve and record ℓ=0,2,4 moments + 4 point probes + energy
function evolve_record(st, e; dt, nsteps, sample, probes)
    s = e.s; Nr=s.grid.Nr; Nθ=s.grid.Nθ
    k1=SphBDNKState(Nr,Nθ);k2=SphBDNKState(Nr,Nθ);k3=SphBDNKState(Nr,Nθ)
    k4=SphBDNKState(Nr,Nθ);tmp=SphBDNKState(Nr,Nθ)
    ts=Float64[]; m0=Float64[]; m2=Float64[]; m4=Float64[]
    en=Float64[]; prs=[Float64[] for _ in probes]
    for n in 0:nsteps
        if n % sample == 0
            push!(ts, n*dt)
            push!(m0, lproj(st,s,0)); push!(m2, lproj(st,s,2)); push!(m4, lproj(st,s,4))
            push!(en, mode_energy(st,e))
            for (pk,(ip,jp)) in enumerate(probes); push!(prs[pk], st.δρ[ip,jp]); end
        end
        n == nsteps && break
        SphBDNK._rk4!(st, e, dt, k1,k2,k3,k4, tmp)
    end
    return (ts=ts, m0=m0, m2=m2, m4=m4, en=en, prs=prs)
end

# linear detrend (kills the slow drift / DC-leakage forest)
function detrend(ts, q)
    n=length(ts); sx=sum(ts); sy=sum(q); sxx=sum(abs2,ts); sxy=sum(ts.*q)
    b=(n*sxy-sx*sy)/(n*sxx-sx^2); a=(sy-b*sx)/n
    q .- (a .+ b.*ts)
end

# find spectral peaks of series q(t) over kHz grid fb; return (fb,Pnorm,peakfreqs)
function spectrum(ts, q, fb; thr=0.02)
    P = periodogram(ts, detrend(ts,q), fb .* kHz2cyc(1.0))
    mx = maximum(P); mx>0 && (P ./= mx)
    loc = Int[]
    for i in 2:length(P)-1
        (P[i]>P[i-1] && P[i] ≥ P[i+1] && P[i] > thr) && push!(loc, i)
    end
    peaks = fb[loc]; pw = P[loc]
    order = sortperm(pw, rev=true)
    (fb=fb, P=P, peaks=peaks[order], pw=pw[order])
end

# ============================================================================
# MAIN
# ============================================================================
function main()
    Nr=96; Nθ=24
    println("#"^72)
    println("# GAUSSIAN-NOISE PLUCK SPECTROSCOPY — causal Shum star (M=1.4 M⊙)")
    println("#"^72, "\n")
    s = build_sphstar(EOS, εc; Nr=Nr, Nθ=Nθ)
    @printf("STAR: M=%.4f M⊙  R=%.4f M⊙  (Nr=%d, Nθ=%d)\n\n", s.M, s.R, Nr, Nθ)

    cz = causality_check(s)

    # eigensolver overlay
    bench2,_,_ = nonradial_cowling_spectrum(EOS, εc; l=2, nmodes=4)
    bench4,_,_ = nonradial_cowling_spectrum(EOS, εc; l=4, nmodes=4)
    @printf("EIGENSOLVER (NonRadialModes Cowling f/p):\n")
    @printf("   ℓ=2 : f=%.4f  p1=%.4f  p2=%.4f  p3=%.4f kHz\n", bench2[1],bench2[2],bench2[3],bench2[4])
    @printf("   ℓ=4 : f=%.4f  p1=%.4f  p2=%.4f  p3=%.4f kHz\n\n", bench4[1],bench4[2],bench4[3],bench4[4])

    # seed the noise
    st = SphBDNKState(Nr,Nθ)
    sc = seed_noise!(st, s; A=1e-4)
    @printf("SEED: band-limited GRF (k≤6 Fourier-r × P_ℓ ℓ≤6, Gaussian cutoff kcut=lcut=4)\n")
    @printf("   scale δρ=%.3e  δv=%.3e  rescale-for-validity=%.4f\n", sc.sρ, sc.sv, sc.minfac)
    val = validity_report(st, s)

    # keep a pristine copy of the seed for the viscous run
    st0 = SphBDNKState(Nr,Nθ)
    for (a,b) in zip((st0.δρ,st0.δSr,st0.δSθ,st0.δvr,st0.δvθ),
                     (st.δρ,st.δSr,st.δSθ,st.δvr,st.δvθ)); copyto!(a,b); end

    # frequency grid up to ~1.15× the highest eigenfreq we overlay
    fhi = 1.15*max(bench2[end], bench4[end])
    fb = collect(range(0.3, fhi; length=6000))

    # ---- IDEAL run (η̂=0) ----
    T = 1500.0
    e0 = setup_sphbdnk(s; η̂=0.0, ν̂=0.0, σ_ko=0.01)
    dt = 0.25*s.grid.dr; nst = round(Int, T/dt); samp = max(1, round(Int, 1.0/dt))
    probes = [(round(Int,f*Nr)+1, jp) for (f,jp) in ((0.3,4),(0.5,8),(0.7,12),(0.85,18))]
    @printf("\nIDEAL EVOLUTION: η̂=0, T=%.0f, dt=%.4f, nsteps=%d, Δf≈%.4f kHz\n",
            T, dt, nst, 1/(T*kHz2cyc(1.0)))
    rec = evolve_record(st, e0; dt=dt, nsteps=nst, sample=samp, probes=probes)
    @printf("   energy E_end/E0 = %.3f\n", rec.en[end]/rec.en[1])

    # spectra of ℓ-moments and (probes summed)
    sp0 = spectrum(rec.ts, rec.m0, fb; thr=0.01)
    sp2 = spectrum(rec.ts, rec.m2, fb; thr=0.01)
    sp4 = spectrum(rec.ts, rec.m4, fb; thr=0.01)
    # also a combined probe spectrum (captures all ℓ at a point)
    Pp = zeros(length(fb)); for p in rec.prs; Pp .+= periodogram(rec.ts, p.-mean(p), fb.*kHz2cyc(1.0)); end
    Pp ./= maximum(Pp)

    # ---- match peaks to eigenfreqs ----
    function match_table(label, bench, sp, modes)
        println("\n  ", label, " peaks vs eigensolver (distinct peak per mode):")
        println("    mode   eigen[kHz]   peak[kHz]   match%   (all peaks: ",
                join((@sprintf("%.3f", x) for x in sort(sp.peaks)[1:min(10,end)]), ", "), ")")
        out = NamedTuple[]
        avail = collect(sp.peaks)                       # assign WITHOUT reuse
        for (mi,b) in enumerate(bench)
            if isempty(avail)
                @printf("    %-5s  %9.4f   %9s   %6s\n", modes[mi], b, "--", "--")
                continue
            end
            j = argmin(abs.(avail .- b)); pk = avail[j]; deleteat!(avail, j)
            mpct = 100*(1 - abs(pk-b)/b)
            @printf("    %-5s  %9.4f   %9.4f   %6.1f\n", modes[mi], b, pk, mpct)
            push!(out, (mode=modes[mi], eig=b, peak=pk, mpct=mpct))
        end
        out
    end
    m2tab = match_table("ℓ=2", bench2, sp2, ["f","p1","p2","p3"])
    m4tab = match_table("ℓ=4", bench4, sp4, ["f","p1","p2","p3"])
    println("\n  ℓ=0 (radial) top peaks [kHz]: ",
            join((@sprintf("%.4f", x) for x in sp0.peaks[1:min(5,end)]), ", "))

    # ---- VISCOUS run (η̂=0.05) from the same seed ----
    for (a,b) in zip((st.δρ,st.δSr,st.δSθ,st.δvr,st.δvθ),
                     (st0.δρ,st0.δSr,st0.δSθ,st0.δvr,st0.δvθ)); copyto!(a,b); end
    eV = setup_sphbdnk(s; η̂=0.05, ν̂=0.0, cν=1.0, σ_ko=0.01)
    @printf("\nVISCOUS EVOLUTION: η̂=0.05 (τ_R,κ_Q,ν_mom all ∝η̂), T=%.0f\n", T)
    recV = evolve_record(st, eV; dt=dt, nsteps=nst, sample=samp, probes=probes)
    @printf("   energy E_end/E0 = %.3f  (ideal was %.3f)\n", recV.en[end]/recV.en[1], rec.en[end]/rec.en[1])
    sp2V = spectrum(recV.ts, recV.m2, fb)
    sp0V = spectrum(recV.ts, recV.m0, fb)

    println("\n  ℓ=2 VISCOUS peaks [kHz]: ",
            join((@sprintf("%.3f", x) for x in sp2V.peaks[1:min(6,end)]), ", "))
    # peak-power comparison at the f-mode location, ideal vs viscous
    function pkpow(sp, f)
        i = argmin(abs.(sp.fb .- f)); sp.P[i]
    end
    @printf("\n  f-mode (ℓ=2, %.3f kHz) relative peak power: ideal=1.000  viscous=%.3f\n",
            bench2[1], (pkpow(sp2V,bench2[1])*maximum(periodogram(recV.ts,recV.m2.-mean(recV.m2),fb.*kHz2cyc(1.0)))) /
                       (pkpow(sp2 ,bench2[1])*maximum(periodogram(rec.ts ,rec.m2 .-mean(rec.m2 ),fb.*kHz2cyc(1.0)))) )

    # spectral-width (FWHM-ish) of the f peak ideal vs viscous: fraction of grid above half-max near f
    function peakwidth(sp, f, half=0.5)
        win = abs.(sp.fb .- f) .< 0.4
        any(win) || return NaN
        Pw = sp.P[win]; fw = sp.fb[win]
        loc = argmax(Pw); pk = Pw[loc]
        above = fw[Pw .> half*pk]
        isempty(above) ? 0.0 : maximum(above)-minimum(above)
    end
    @printf("  f-mode spectral width (kHz, ~FWHM): ideal=%.4f  viscous=%.4f\n",
            peakwidth(sp2,bench2[1]), peakwidth(sp2V,bench2[1]))

    println("\nDONE.")
    return (cz=cz, val=val, bench2=bench2, bench4=bench4,
            m2tab=m2tab, m4tab=m4tab,
            sp0=sp0, sp2=sp2, sp4=sp4, sp2V=sp2V,
            en_ideal=rec.en[end]/rec.en[1], en_visc=recV.en[end]/recV.en[1])
end

main()
