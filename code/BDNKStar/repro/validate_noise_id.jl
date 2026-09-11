#=
    INDEPENDENT VALIDATION of the band-limited Gaussian-noise initial data +
    spectroscopy robustness for the SphBDNK mode-excitation experiment.
    (Stress-test, not part of the build.)
=#
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "viz"))
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using .BDNKStar: SphBackground, Units
using Random, Printf, Statistics

const Msun_to_km = BDNKStar.Units.Msun_to_km
const kHz_to_km  = BDNKStar.Units.kHz_to_km

# ---------------------------------------------------------------- Shum star
eos = ShumPolytrope(100.0)
εc  = 0.00128 + 100*0.00128^2
Nr, Nθ = 96, 24
s = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ)
println(@sprintf("Shum star: R=%.4f (%.2f km)  M=%.4f Msun  Nr=%d Nθ=%d  εc=%.6g",
        s.R, s.R*Msun_to_km, s.M, Nr, Nθ, εc))

# ---------------------------------------------------- band-limited Gaussian ID
# Build noise on the interior (Nr×Nθ), keep only low-k radial & angular Fourier
# content (band-limited), respect the centre/axis parity that _fill_ghosts!
# expects: δρ EVEN at centre & poles; δSr ODD at centre, δvr ODD at centre;
# δSθ ODD at poles, δvθ ODD at poles. We seed δρ (the density mode) — a scalar
# field that is even everywhere — and recover δv via the evolution.
# Band limit: keep radial modes n=1..Kr, angular Legendre-like cos(mθ) m=1..Kθ.
function bandlimited_noise!(st, s; A=1e-4, Kr=6, Kθ=3, seed=1, parity_even=true)
    rng = MersenneTwister(seed)
    Nr=s.grid.Nr; Nθ=s.grid.Nθ; R=s.R
    # random amplitudes/phases for a low-k cosine/sine basis
    ar = randn(rng, Kr); aθ = randn(rng, Kθ)
    # δρ even at centre (sin(nπr/R) → 0 at r=0, even reflection ok since field→0)
    # and even at poles via cos(mθ) basis with m chosen so ∂θ=0 at poles? We use
    # Legendre P_l(cosθ) for l=2..(Kθ+1): even in cosθ → even at both poles.
    for jj in 1:Nθ, ii in 1:Nr
        r=s.grid.r[ii]; c=s.grid.cosθ[jj]
        radial = 0.0
        for n in 1:Kr
            radial += ar[n]*sin(n*π*r/R)        # vanishes at r=0 and r=R: regular + Δp~0
        end
        ang = 0.0
        for l in 2:(Kθ+1)
            # Legendre P_l(cosθ), even-parity (even l) and odd mixed; build a few
            Pl = l==2 ? (3c^2-1)/2 :
                 l==3 ? (5c^3-3c)/2 :
                 l==4 ? (35c^4-30c^2+3)/8 :
                        (63c^5-70c^3+15c)/8
            ang += aθ[l-1]*Pl
        end
        st.δρ[ii+1,jj+1] = A*s.ε0[ii]*radial*ang
    end
    st
end

# grid-scale (un-band-limited) white noise: independent per cell, full spectrum
function gridscale_noise!(st, s; A=1e-4, seed=1)
    rng = MersenneTwister(seed)
    Nr=s.grid.Nr; Nθ=s.grid.Nθ
    for jj in 1:Nθ, ii in 1:Nr
        st.δρ[ii+1,jj+1] = A*s.ε0[ii]*randn(rng)
    end
    st
end

# ===================================================================
# 1. CONSTRAINT VALIDITY
# ===================================================================
println("\n========== 1. CONSTRAINT VALIDITY ==========")
e = setup_sphbdnk(s; η̂=0.0)             # inviscid baseline for ID checks
st = SphBDNKState(Nr, Nθ)
bandlimited_noise!(st, s; A=1e-4, seed=1)

# ε0 + δε > 0 everywhere (δε = h0 δρ)
minε, minε0 = let mε=Inf, mε0=Inf
    for jj in 1:Nθ, ii in 1:Nr
        h0 = (s.ε0[ii]+s.p0[ii])/s.ρ0[ii]
        δε = h0*st.δρ[ii+1,jj+1]
        mε  = min(mε, s.ε0[ii]+δε)
        mε0 = min(mε0, s.ε0[ii])
    end
    (mε, mε0)
end
maxδρ = maximum(abs, st.δρ)
println(@sprintf("min ε0           = %.6e", minε0))
println(@sprintf("min(ε0+δε)       = %.6e   (positivity %s)", minε, minε>0 ? "OK" : "VIOLATED"))
println(@sprintf("max|δρ|/ε        = %.3e   (relative perturbation amplitude)", maxδρ/minε0))

# parity / regularity: fill ghosts and check no NaN/Inf and bounded reflection
BDNKStar.SphBDNK._fill_ghosts!(st, Nr, Nθ)
finite_ghost = all(isfinite, st.δρ) && all(isfinite, st.δSr) && all(isfinite, st.δvr)
# centre regularity: δρ even => ghost(1,j)=+δρ(2,j); confirm
ctr_ok = all(st.δρ[1,j] ≈ st.δρ[2,j] for j in 2:Nθ+1)
axis_ok = all(st.δρ[i,1] ≈ st.δρ[i,2] for i in 2:Nr+1)
println(@sprintf("ghost fill finite = %s   centre δρ even = %s   axis δρ even = %s",
        finite_ghost, ctr_ok, axis_ok))

# subluminal velocities: ID seeds δρ only (δv=0 initially → trivially subluminal).
# After a few steps, recover δv and check |v|<1. Run a short evolution.
dt = 0.2*s.grid.dr
evolve_sphbdnk!(st, e; dt=dt, nsteps=200, sample=200)
vmax = let vm=0.0
    for jj in 1:Nθ, ii in 1:Nr
        vr=st.δvr[ii+1,jj+1]; vθ=st.δvθ[ii+1,jj+1]
        vphys = sqrt(s.eΛ[ii]*vr^2 + s.grid.r[ii]^2*vθ^2)
        vm = max(vm, vphys)
    end
    vm
end
allfin = all(isfinite, st.δρ)
println(@sprintf("after 200 steps: max|δv|_phys = %.3e (<1 subluminal=%s)  fields finite=%s",
        vmax, vmax<1, allfin))

# ===================================================================
# 2. CAUSALITY (Shum/Kovtun route)
# ===================================================================
println("\n========== 2. CAUSALITY (Shum/Kovtun route) ==========")
wp = shum_frame_wellposed(1.0, 1.0, 0.999)
println(@sprintf("shum_frame_wellposed(1,1,0.999) = %s  (0<q̂<ŝ ⇒ well-posed+stable)", wp))
# c₊ = √3 cs(r) at every radius; require ≤ 1
maxcp, rcp, cs_at = let mc=0.0, rc=0.0, ca=0.0
    for ii in 1:Nr
        cs = sqrt(max(s.cs2[ii], 0.0))
        _,cp,cm = shum_frame_speeds(1.0, 1.0, 0.999, 0.01, 0.01, cs)  # c₊=√3·cs
        if cp > mc; mc=cp; rc=s.grid.r[ii]; ca=cs; end
    end
    (mc, rc, ca)
end
csmax = sqrt(maximum(s.cs2))
println(@sprintf("max cs(r)        = %.6f", csmax))
println(@sprintf("max c₊ = √3·cs   = %.6f  at r=%.4f (cs=%.6f)", maxcp, rcp, cs_at))
println(@sprintf("subluminal margin 1-max c₊ = %.6f  (%s)", 1-maxcp, maxcp<=1 ? "CAUSAL" : "ACAUSAL"))
# cross-check c₊/cs ratio: at q̂=0.999 it is ≈√3 (exactly √3 only as q̂→1)
ratio = maxcp/cs_at
println(@sprintf("c₊/cs(at max) = %.6f   √3 = %.6f   (q̂=0.999 ⇒ %.3f%% below √3, expected)",
        ratio, sqrt(3), 100*(sqrt(3)-ratio)/sqrt(3)))

# ===================================================================
# 3. SEED-INDEPENDENCE of peak frequencies
# ===================================================================
println("\n========== 3. SEED-INDEPENDENCE (≥3 seeds) ==========")
# Cowling benchmark f-mode for this star (independent eigensolver) for context:
fbench,_,_ = nonradial_cowling_spectrum(eos, εc; l=2, nmodes=5, N=6000, nscan=900)
println(@sprintf("Cowling eigensolver l=2 (f,p1..): %s kHz  (independent reference)",
        string(round.(fbench, digits=3))))

# Long evolution so the periodogram resolves modes: with dt=0.2dr the f-mode
# (~1.9 kHz) period is ~108 Msun; nsteps=80000 → T~1600 Msun (~15 cycles),
# Δν~0.13 kHz frequency resolution.
const NSTEPS = 80000
const SAMPLE = 4
νs_kHz = collect(range(0.3, 9.0; length=2400))
νs_geo = νs_kHz .* (Msun_to_km*kHz_to_km)

# Hann-windowed periodogram to suppress spectral leakage (cleaner peaks)
function wperiodogram(ts, q, νgeo)
    n=length(q); q̄=q .- sum(q)/n
    w=[0.5-0.5*cos(2π*(k-1)/(n-1)) for k in 1:n]; q̄ .*= w
    P=similar(νgeo)
    @inbounds for (m,ν) in enumerate(νgeo)
        ω=2π*ν; re=0.0; im=0.0
        for k in 1:n; sc=sincos(ω*ts[k]); re+=q̄[k]*sc[2]; im-=q̄[k]*sc[1]; end
        P[m]=re^2+im^2
    end
    P
end

function run_spectroscopy(seed; A=1e-4, idfun=bandlimited_noise!, nsteps=NSTEPS, sample=SAMPLE)
    e = setup_sphbdnk(s; η̂=0.0)
    st = SphBDNKState(Nr, Nθ)
    idfun(st, s; A=A, seed=seed)
    dt = 0.2*s.grid.dr
    ts, q2, en = evolve_sphbdnk!(st, e; dt=dt, nsteps=nsteps, sample=sample)
    finite = all(isfinite, q2)
    P = finite ? wperiodogram(ts, q2, νs_geo) : fill(NaN, length(νs_geo))
    return ts, q2, en, P, finite
end

# locate dominant peaks (local maxima), strongest first
function peaks(νkHz, P; npk=5, minsep=0.3, relthr=0.02)
    isfinite(sum(P)) || return Float64[]
    Pm=maximum(P)
    idx = Int[]
    for n in 2:length(P)-1
        if P[n]>P[n-1] && P[n]≥P[n+1] && P[n]>relthr*Pm; push!(idx,n); end
    end
    sort!(idx, by=i->-P[i])
    out = Float64[]
    for i in idx
        f = νkHz[i]
        all(abs(f-g)>minsep for g in out) || continue
        push!(out, f)
        length(out)≥npk && break
    end
    sort(out)
end

Δν_kHz = (νs_kHz[end]-νs_kHz[1])/(length(νs_kHz)-1)
seeds = [1, 7, 42, 123]
peaktab = Vector{Vector{Float64}}()
for sd in seeds
    ts,q2,en,P,fin = run_spectroscopy(sd)
    pk = peaks(νs_kHz, P)
    push!(peaktab, pk)
    Pmax = fin ? maximum(P) : NaN
    @printf("seed=%-4d finite=%s  peaks[kHz]=%s   Pmax=%.3e  q2max=%.3e\n",
            sd, fin, isempty(pk) ? "[]" : string(round.(pk,digits=3)),
            Pmax, fin ? maximum(abs,q2) : NaN)
end

# A physical mode = a peak present in EVERY seed within the freq resolution.
# Cluster all seeds' peaks; report clusters that appear in all 4 seeds.
println(@sprintf("\nfreq resolution Δν ≈ %.3f kHz; cluster tolerance = 2Δν = %.3f kHz", Δν_kHz, 2Δν_kHz))
allpk = sort(vcat(peaktab...))
tol = max(5Δν_kHz, 0.05)
clusters = Vector{Vector{Float64}}()
for f in allpk
    placed=false
    for c in clusters
        if abs(f - mean(c)) ≤ tol; push!(c,f); placed=true; break; end
    end
    placed || push!(clusters, [f])
end
println("common (seed-independent) mode peaks [kHz | #seeds | scatter]:")
for c in sort(clusters, by=mean)
    # count distinct seeds contributing
    nseed = count(any(abs(p-mean(c))≤tol for p in pk) for pk in peaktab)
    flag = nseed==length(seeds) ? "★ ALL SEEDS" : "($nseed/$(length(seeds)) seeds)"
    @printf("  f=%.3f  scatter=%.3f kHz  %s\n", mean(c), maximum(c)-minimum(c), flag)
end

# ===================================================================
# 4. BAND-LIMIT SENSITIVITY: grid-scale noise excites numerical growth
# ===================================================================
println("\n========== 4. BAND-LIMIT SENSITIVITY ==========")
# Same length/window for both; energy series is the faithful growth diagnostic.
tsb,q2b,enb,Pb,finb = run_spectroscopy(1; idfun=bandlimited_noise!)
tsg,q2g,eng,Pg,fing = run_spectroscopy(1; idfun=gridscale_noise!)

# energy end/start; a clean mode bath stays O(1), grid-scale noise feeds the
# scheme's slow grid-scale numerical GROWTH (documented in SphBDNK header).
egrow(en) = (isempty(en)||en[1]==0) ? NaN : en[end]/en[1]
pkb = peaks(νs_kHz, Pb); pkg = peaks(νs_kHz, Pg)
@printf("band-limited : finite=%s  E_end/E0=%.3e  q2 end/start=%.3e  peaks=%s\n",
        finb, egrow(enb), finb ? abs(q2b[end])/max(abs(q2b[1]),1e-300) : NaN,
        isempty(pkb) ? "[]" : string(round.(pkb,digits=3)))
@printf("grid-scale   : finite=%s  E_end/E0=%.3e  q2 end/start=%.3e  peaks=%s\n",
        fing, egrow(eng), fing ? abs(q2g[end])/max(abs(q2g[1]),1e-300) : NaN,
        isempty(pkg) ? "[]" : string(round.(pkg,digits=3)))
# spectral sharpness: peak power / median power. Clean modes ≫ broadband floor.
sharp(P) = isfinite(sum(P)) ? maximum(P)/median(P) : NaN
@printf("spectral sharpness (max/median P): band-limited=%.2e   grid-scale=%.2e\n",
        sharp(Pb), sharp(Pg))
# high-freq (grid-scale, ν>5 kHz) power fraction: large ⇒ noise excites the grid,
# not clean low-order modes.
hf(P,νk) = (m=νk.>5.0; isfinite(sum(P)) ? sum(P[m])/sum(P) : NaN)
@printf("high-freq (>5 kHz) power fraction: band-limited=%.3f   grid-scale=%.3f\n",
        hf(Pb,νs_kHz), hf(Pg,νs_kHz))
println("\n========== DONE ==========")
