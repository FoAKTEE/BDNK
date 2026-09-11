#=
    radial_overtone_crosscheck.jl

    INDEPENDENT cross-check of the full-GR radial OVERTONE spectrum (F,H1,H2,H3,H4)
    of the BDNKStar benchmark star, by EXTENDING the GHZ(1997) shooting cross-check
    (repro/radial_shoot_crosscheck.jl) from the fundamental to the overtones — a
    DIFFERENT method than the Sturm–Liouville matrix eigensolver
    `chandrasekhar_radial_omega2` / `sl_modes` and the time-domain DynGR1D engine.

    PART 1 — SHOOTING OVERTONES
      Each overtone H_n is the (n+1)-th ω² root of the surface condition Δp(R)=0,
      with the displacement eigenfunction ξ(r) having exactly n interior nodes.
      We scan ω², bracket every sign change of Δp(R), bisect each, and node-count
      ξ at each eigenvalue. Cross-check each H_n against the SL eigensolver to ≳4
      digits, and confirm node count == mode number.

    Two reference targets:
      (a) the engine benchmark star  εc = 0.0015 km^-2 (NOT an exact KR table row;
          published ν only by interpolation). Primary check here is shoot-vs-SL.
      (b) an EXACT KR(2001) Table A.18 row (ρ_c = 2.0e15 g/cm^3, n=1 κ=100 polytrope)
          for a DIRECT published-frequency comparison of F, H1, H2.

    PART 2 — DYNAMICAL ENGINE OVERTONES
      Excite the engine (DynGR1D) with a node-bearing radial velocity profile
      (sin(2πr/R), one interior node — projects onto H1) and FFT ρ_c(t) and a
      probe near a velocity antinode. Report which radial modes (F, H1, ...) the
      time-domain dynamical-GR engine resolves and whether the peak sits at the
      eigenvalue H1.
=#

using BDNKStar
using BDNKStar: solve_tov, ShumPolytrope, pressure, sound_speed2,
                energy_from_pressure, rho_from_p,
                setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!,
                dyngr_radial_freq
using BDNKStar.Units: kHz_to_km, gram_per_cm3_to_km_minus2
using LinearAlgebra: eigen, Symmetric
using Printf

const CONV_kHz = 1.0 / kHz_to_km
fkHz(ω2) = sqrt(max(ω2, 0.0)) / (2π) * CONV_kHz

const EOS  = ShumPolytrope(100.0)
const KAPPA = 100.0
const HTOV = 2e-5

# ρ_c[1e15 g/cm3] -> εc[km^-2] via ε = ρ + κ ρ²  (exact for KR n=1 κ=100 polytrope)
epsc_from_rhoc(rhoc_1e15) = let ρ = rhoc_1e15*1e15*gram_per_cm3_to_km_minus2
    ρ + KAPPA*ρ^2
end

# ===========================================================================
# SHOOTING (co-integrated GHZ (ξ,Δp)) — generalized to a chosen star (εc).
# Returns Δp(R), ξ(R), interior node count of ξ, ν(R), R_int.
# ===========================================================================
@inline function deriv(r, y, ω2, εc)
    m, p, ν, ξ, Δp = y
    p ≤ 0 && return (0.0, 0.0, 0.0, 0.0, 0.0)
    ε   = energy_from_pressure(EOS, p)
    den = r * (r - 2m); fac = m + 4π * r^3 * p
    dm  = 4π * r^2 * ε
    dν  = 2 * fac / den
    pp  = -(ε + p) * fac / den
    eλ  = 1.0 / max(1 - 2m / r, 1e-12)
    eνm = exp(-ν)
    cs2 = clamp(sound_speed2(EOS, ε), 1e-12, 1.0)
    Γ1  = cs2 * (ε + p) / p
    V = -3.0 / r - pp / (ε + p)
    W = -(1.0 / r) / (Γ1 * p)
    X = ω2 * eλ * eνm * (ε + p) * r - 4.0 * pp + pp^2 * r / (ε + p) -
        8π * eλ * (ε + p) * p * r
    Y = pp / (ε + p) - 4π * (ε + p) * r * eλ
    return (dm, pp, dν, V*ξ + W*Δp, X*ξ + Y*Δp)
end

# in-place 5-tuple update a + s*b (avoids the closure ALLOCATIONS that ntuple(i->..)
# incurs — those boxed ~61 MB/shoot and made each shoot ~65x slower).
@inline _add5(a, b, s) =
    (a[1]+s*b[1], a[2]+s*b[2], a[3]+s*b[3], a[4]+s*b[4], a[5]+s*b[5])

function shoot(ω2, εc; h::Float64=3e-4, ptol_rel::Float64=1e-8)
    pc = pressure(EOS, εc); ptol = ptol_rel * pc
    r0 = h
    m = (4π/3) * εc * r0^3
    p = pc - 2π*(εc + pc)*(εc/3 + pc) * r0^2
    ν = 0.0
    ε0 = energy_from_pressure(EOS, p)
    Γ10 = clamp(sound_speed2(EOS, ε0), 1e-12, 1.0) * (ε0 + p) / p
    ξ = 1.0; Δp = -3.0 * Γ10 * p * ξ
    y = (m, p, ν, ξ, Δp); r = r0; nodes = 0; prevξ = ξ
    Δp_s = Δp; ξ_s = ξ; ν_s = ν; R_int = r0
    while p > ptol && r < 100.0
        k1 = deriv(r,     y, ω2, εc)
        k2 = deriv(r+h/2, _add5(y, k1, h/2), ω2, εc)
        k3 = deriv(r+h/2, _add5(y, k2, h/2), ω2, εc)
        k4 = deriv(r+h,   _add5(y, k3, h),   ω2, εc)
        yn = (y[1] + h/6*(k1[1] + 2k2[1] + 2k3[1] + k4[1]),
              y[2] + h/6*(k1[2] + 2k2[2] + 2k3[2] + k4[2]),
              y[3] + h/6*(k1[3] + 2k2[3] + 2k3[3] + k4[3]),
              y[4] + h/6*(k1[4] + 2k2[4] + 2k3[4] + k4[4]),
              y[5] + h/6*(k1[5] + 2k2[5] + 2k3[5] + k4[5]))
        r += h; pn = yn[2]
        if pn ≤ ptol
            frac = y[2] / (y[2] - pn)
            ξ_s  = y[4] + frac*(yn[4] - y[4])
            Δp_s = y[5] + frac*(yn[5] - y[5])
            ν_s  = y[3] + frac*(yn[3] - y[3])
            R_int = (r - h) + frac*h
            break
        end
        if pn > 1e-4*pc && sign(yn[4]) != sign(prevξ)
            nodes += 1
        end
        prevξ = yn[4]
        p = pn; y = yn
    end
    return Δp_s, ξ_s, nodes, ν_s, R_int
end

# Schwarzschild-ν shift: ω²_phys = ω²_raw · e^{Δν}, Δν = ln(1-2M/R_int) - ν_raw(R)
function nu_shift(εc, M)
    (_, _, _, νR, Rint) = shoot(1e-3, εc)
    return log(1 - 2M/Rint) - νR
end

# Scan ω²_raw, bracket EVERY Δp(R) sign change, bisect each → ordered ω² roots
# with node counts. Returns Vector of (ω2_raw, nodes).
function shoot_spectrum(εc; ω2lo=1e-4, ω2hi=0.10, nscan=1200)
    surf(ω2) = shoot(ω2, εc)[1]
    grid = collect(range(ω2lo, ω2hi; length=nscan))
    prevv = surf(grid[1]); prevω = grid[1]
    roots = Tuple{Float64,Int}[]
    for ω2 in grid[2:end]
        cur = surf(ω2)
        if sign(cur) != sign(prevv) && isfinite(cur) && isfinite(prevv)
            a, c = prevω, ω2; fa = surf(a)
            for _ in 1:200
                mid = 0.5*(a+c); fm = surf(mid)
                sign(fm) == sign(fa) ? (a = mid; fa = fm) : (c = mid)
                (c - a) < 1e-12 && break
            end
            ω2r = 0.5*(a+c); nd = shoot(ω2r, εc)[3]
            push!(roots, (ω2r, nd))
        end
        prevv = cur; prevω = ω2
    end
    return roots
end

# ===========================================================================
# SL eigensolver returning ω² AND node counts of the ζ eigenvector.
# (verbatim from chandrasekhar_radial_omega2 / kr2001_radial_spectrum.sl_modes)
# ===========================================================================
@inline _lin(xs, ys, x) = begin
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end
function sl_modes(εc; N=3000, h_tov=HTOV, nmodes=6)
    star = solve_tov(EOS, εc; h=h_tov); R=star.R
    rt=star.r; mt=star.m; νt=star.ν; εt=star.ε; pt=star.p
    bg(r)=(max(_lin(rt,mt,r),0.0), _lin(rt,νt,r),
           max(_lin(rt,εt,r),1e-20), max(_lin(rt,pt,r),1e-30))
    n=N; r=collect(range(R/n, R*(1-1e-6); length=n)); dr=r[2]-r[1]
    P=zeros(n); Wt=zeros(n); Q=zeros(n)
    for i in 1:n
        m,ν,ε,p=bg(r[i]); eλ=1/max(1-2m/r[i],1e-12); λ=log(eλ)
        cs2=clamp(sound_speed2(EOS,ε),1e-12,1.0); Γ1=cs2*(ε+p)/p
        νp=2*(m+4π*r[i]^3*p)/(r[i]*(r[i]-2m))
        ePQ=exp((λ+3ν)/2); eW=exp((3λ+ν)/2)
        P[i]=Γ1*p*ePQ/r[i]^2; Wt[i]=(ε+p)*eW/r[i]^2
        Q[i]=ePQ*(ε+p)/r[i]^2*(νp^2/4 + 2νp/r[i] - 8π*eλ*p)
    end
    A=zeros(n,n); B=zeros(n,n); Pf=zeros(n+1)
    for i in 1:n-1; Pf[i+1]=0.5*(P[i]+P[i+1]); end
    Pf[1]=P[1]; Pf[n+1]=0.0
    for i in 1:n
        aw=Pf[i]/dr^2; ae=Pf[i+1]/dr^2
        A[i,i]=aw+ae-Q[i]; i>1 && (A[i,i-1]=-aw); i<n && (A[i,i+1]=-ae); B[i,i]=Wt[i]
    end
    F=eigen(Symmetric(A),Symmetric(B))
    idx=sortperm(real.(F.values)); ω2=real.(F.values)[idx]
    k=min(nmodes,n); ω2k=ω2[1:k]
    nodes=Int[]
    for j in 1:k
        v=real.(F.vectors[:,idx[j]]); nc=0
        for i in 2:n-1; (v[i]*v[i-1] < 0) && (nc+=1); end
        push!(nodes,nc)
    end
    return ω2k, nodes, star
end

modelabel(n) = n==0 ? "F" : "H$n"

# ===========================================================================
# PART 1A — benchmark star εc=0.0015 : shooting overtones vs SL
# ===========================================================================
const EPSC_BENCH = 0.0015
star_b = solve_tov(EOS, EPSC_BENCH; h=HTOV)
const Rb = star_b.R; const Mb = star_b.M

println("# ==========================================================================")
println("# PART 1 — SHOOTING OVERTONE SPECTRUM vs SL EIGENSOLVER")
println("# ==========================================================================")
@printf("# benchmark star: ShumPolytrope(kappa=100), eps_c=%.5g km^-2\n", EPSC_BENCH)
@printf("# R=%.5f km  M=%.5f km  2M/R=%.5f\n\n", Rb, Mb, 2Mb/Rb)

Δν_b = nu_shift(EPSC_BENCH, Mb)
# NOTE: shoot_spectrum scans in RAW ω² (pre Schwarzschild-ν shift, factor e^Δν≈0.44);
# the 5 physical modes F,H1..H4 are the first 5 raw roots, spanning raw ω²∈[0.004,0.21].
roots_b = shoot_spectrum(EPSC_BENCH; ω2lo=1e-4, ω2hi=0.25, nscan=5000)
print("# PART 1A shooting scan done\n"); flush(stdout)
ω2sl_b, nodes_sl_b, _ = sl_modes(EPSC_BENCH; N=2000, nmodes=6)
print("# PART 1A SL solve done\n"); flush(stdout)

println("=== (1A) benchmark eps_c=0.0015 : SHOOTING vs SL ===")
@printf("%5s | %14s %8s %6s | %14s %8s %6s | %8s %8s\n",
        "mode","w2_shoot","F_shoot","nodes","w2_SL","F_SL","nodes","dF%","match")
nshow = min(5, length(roots_b))
bench_rows = NamedTuple[]
for n in 0:nshow-1
    ω2raw, nd = roots_b[n+1]
    ω2sh = ω2raw * exp(Δν_b)
    Fsh = fkHz(ω2sh)
    ω2s = ω2sl_b[n+1]; Fs = fkHz(ω2s); nds = nodes_sl_b[n+1]
    pct = 100*(Fsh-Fs)/Fs
    ok = abs(pct) < 0.1 ? "OK" : (abs(pct)<1 ? "~" : "X")
    @printf("%5s | %14.6e %8.4f %6d | %14.6e %8.4f %6d | %+7.3f %8s\n",
            modelabel(n), ω2sh, Fsh, nd, ω2s, Fs, nds, pct, ok)
    push!(bench_rows, (mode=modelabel(n), Fsh=Fsh, Fsl=Fs, nodes=nd, nodes_sl=nds, pct=pct))
end

# ===========================================================================
# PART 1B — EXACT KR(2001) Table A.18 row ρ_c=2.0e15 g/cm^3 : DIRECT published
#           comparison of F=ν0, H1=ν1, H2=ν2.
# ===========================================================================
const RHOC_KR = 2.000           # 1e15 g/cm^3 (KR Table A.18 row)
const KR_PUB  = (R=9.673, M=1.126, ν0=2.323, ν1=6.237, ν2=9.295)  # KR A.18
εc_kr = epsc_from_rhoc(RHOC_KR)
star_kr = solve_tov(EOS, εc_kr; h=HTOV)
Δν_kr = nu_shift(εc_kr, star_kr.M)
roots_kr = shoot_spectrum(εc_kr; ω2lo=1e-4, ω2hi=0.25, nscan=5000)
ω2sl_kr, nodes_sl_kr, _ = sl_modes(εc_kr; N=2000, nmodes=6)
print("# PART 1B done\n"); flush(stdout)

println("\n=== (1B) EXACT KR(2001) Table A.18 row rho_c=2.0e15 g/cm^3 ===")
@printf("# KR published: R=%.3f km  M=%.3f Msun  nu0=%.3f nu1=%.3f nu2=%.3f kHz\n",
        KR_PUB.R, KR_PUB.M, KR_PUB.ν0, KR_PUB.ν1, KR_PUB.ν2)
@printf("# our star    : R=%.3f km  M=%.4f km\n", star_kr.R, star_kr.M)
@printf("%5s | %9s %9s %9s | %8s %8s\n",
        "mode","F_shoot","F_SL","nu_KR","sh-SL%","sh-KR%")
kr_pub_vec = [KR_PUB.ν0, KR_PUB.ν1, KR_PUB.ν2]
kr_rows = NamedTuple[]
for n in 0:min(4, length(roots_kr)-1)
    ω2raw, nd = roots_kr[n+1]
    Fsh = fkHz(ω2raw * exp(Δν_kr))
    Fs  = fkHz(ω2sl_kr[n+1])
    νpub = n < 3 ? kr_pub_vec[n+1] : NaN
    pSL  = 100*(Fsh-Fs)/Fs
    pKR  = isnan(νpub) ? NaN : 100*(Fsh-νpub)/νpub
    @printf("%5s | %9.4f %9.4f %9s | %+7.3f %s\n",
            modelabel(n), Fsh, Fs, isnan(νpub) ? "-" : @sprintf("%.4f",νpub),
            pSL, isnan(pKR) ? "    -" : @sprintf("%+7.3f", pKR))
    push!(kr_rows, (mode=modelabel(n), Fsh=Fsh, Fsl=Fs, nodes=nd,
                    nodes_sl=nodes_sl_kr[n+1], pub=isnan(νpub) ? nothing : νpub,
                    pctKR=pKR))
end

# ordering / node sanity
ωsh_b = [roots_b[i+1][1]*exp(Δν_b) for i in 0:nshow-1]
ord_ok = all(diff(ωsh_b) .> 0)
nodes_ok = all(roots_b[i+1][2] == i for i in 0:nshow-1)
@printf("\n# ordering F<H1<H2<...: %s   node count == mode number: %s\n",
        ord_ok ? "OK" : "FAIL", nodes_ok ? "OK" : "FAIL")
flush(stdout)   # ensure PART 1 results are emitted before the slow engine runs

# ===========================================================================
# PART 2 — DYNAMICAL ENGINE OVERTONES
#   Excite a node-bearing velocity profile and FFT. Can the engine resolve H1?
# ===========================================================================
println("\n# ==========================================================================")
println("# PART 2 — DYNAMICAL ENGINE (DynGR1D) OVERTONE EXCITATION")
println("# ==========================================================================")

F_eig  = fkHz(ω2sl_b[1])
H1_eig = fkHz(ω2sl_b[2])
H2_eig = fkHz(ω2sl_b[3])
@printf("# eigenvalue targets (SL): F=%.4f  H1=%.4f  H2=%.4f kHz\n\n", F_eig, H1_eig, H2_eig)

# Custom node-bearing seed: v = A * sin(2π r/R) (ONE interior node @ r=R/2),
# orthogonal to the nodeless fundamental → should project strongly onto H1.
function seed_nodebearing!(st, eng; A=2e-3, k=2)
    g=eng.g; eos=eng.eos; N=g.N; NG=g.NG
    @inbounds for i in 1:N
        ai=NG+i; r=g.r[ai]
        if r < eng.R
            v = A*sin(k*π*r/eng.R)
            st.v[ai]=v
            sg=st.sqrtg[ai]
            D̂,Ŝ,τ̂ = BDNKStar.DynGR1D.prim2cons_barotrope(eos, st.ρ[ai], st.p[ai], abs(v))
            st.D[ai]=sg*D̂; st.S[ai]=sg*sign(v)*Ŝ; st.τ[ai]=sg*τ̂
        end
    end
    BDNKStar.DynGR1D._fill_ghosts!(st, eng)
    BDNKStar.DynGR1D.dyngr_metric!(st, eng)
end

# Period of fundamental in km (geometric): T ~ 1/F[km^-1]
F_geom = sqrt(ω2sl_b[1])/(2π)
H1_geom = sqrt(ω2sl_b[2])/(2π)
# The dynamical engine is only LONG-TERM STABLE for ~15-20 fundamental periods
# (the validated benchmark used tmax≈220 R ≈ 15 periods; beyond that the
# constrained-evolution central density drifts and collapses). Use 18 periods.
const TPER = 1.0/F_geom                       # one fundamental period in km
                                              # (F_geom = sqrt(ω²)/(2π) is cycles/km)
tmax = 18.0 * TPER

# Drop any post-collapse / non-finite tail BEFORE the periodogram (a diverging
# ρ_c tail otherwise swamps the windowed FFT and erases the mode peaks).
function _clean_series(ts, q, ρc0)
    n = length(ts); keep = n
    @inbounds for i in 1:n
        if !isfinite(q[i]) || !isfinite(ts[i]) || abs(q[i]) > 50*abs(ρc0)+1e-30
            keep = i-1; break
        end
    end
    keep = max(keep, 0)
    return ts[1:keep], q[1:keep]
end

# Probe near a velocity ANTINODE of H1 (r≈R/4, where sin(2πr/R) is near max)
function run_engine(seed_fn, label; A, k=2, Ngrid=600, probe_frac=0.25)
    eng, st = setup_dyngr(EOS, EPSC_BENCH; N=Ngrid, rmax_fac=1.3, cfl=0.3,
                          h_tov=HTOV)
    seed_fn(st, eng; A=A, k=k)
    println("  [running engine: ", label, ", N=", Ngrid, ", tmax=",
            round(tmax,digits=2), " km, ", round(tmax/TPER,digits=1),
            " periods] ..."); flush(stdout)
    res = evolve_dyngr!(st, eng; tmax=tmax, probe_frac=probe_frac,
                        sample_dt=0.5*eng.g.Δr)
    println("    -> steps done, ", length(res.ts), " samples, drift=",
            round(res.drift,sigdigits=3), " collapsed=", res.collapsed); flush(stdout)
    fmaxkm = 1.5*H1_geom*1.6    # span a bit above H2 too
    ρc0 = length(res.ρc) > 0 ? res.ρc[1] : 1.0
    ts_rc, q_rc = _clean_series(res.ts, res.ρc, ρc0)
    ts_pr, q_pr = _clean_series(res.ts, res.probe, ρc0)
    fpk_rc, Prc, fg = dyngr_radial_freq(ts_rc, q_rc; nf=4000, fmax=fmaxkm)
    fpk_pr, Ppr, _  = dyngr_radial_freq(ts_pr, q_pr; nf=4000, fmax=fmaxkm)
    println("    -> usable samples: rho_c=", length(ts_rc),
            " probe=", length(ts_pr)); flush(stdout)
    return res, fg, Prc, Ppr
end

# find all significant peaks (local maxima above frac*global max) in kHz
function peaks_kHz(fg, P; relthresh=0.12)
    Pmx = maximum(P); out=Tuple{Float64,Float64}[]
    for k in 2:length(P)-1
        if P[k]>P[k-1] && P[k]≥P[k+1] && P[k] > relthresh*Pmx
            push!(out, (fg[k]*CONV_kHz, P[k]/Pmx))
        end
    end
    sort!(out, by=x->-x[2])
    return out
end

println("=== (2A) HOMOLOGOUS seed v ~ sin(π r/R) (nodeless → fundamental F) ===")
res1, fg1, Prc1, Ppr1 = run_engine(
    (st,eng;A,k)->seed_dyngr_velocity!(st,eng;A=A,profile=:sin), "homol"; A=2e-4, k=1)
pk_rc1 = peaks_kHz(fg1, Prc1); pk_pr1 = peaks_kHz(fg1, Ppr1)
@printf("  rho_c peaks (kHz, rel.power): %s\n",
        join([@sprintf("%.4f(%.2f)",f,p) for (f,p) in pk_rc1[1:min(4,end)]], "  "))
@printf("  probe peaks (kHz, rel.power): %s\n",
        join([@sprintf("%.4f(%.2f)",f,p) for (f,p) in pk_pr1[1:min(4,end)]], "  "))
flush(stdout)

println("\n=== (2B) NODE-BEARING seed v ~ sin(2π r/R) (1 node → projects onto H1) ===")
res2, fg2, Prc2, Ppr2 = run_engine(seed_nodebearing!, "node"; A=3e-4, k=2, probe_frac=0.25)
pk_rc2 = peaks_kHz(fg2, Prc2); pk_pr2 = peaks_kHz(fg2, Ppr2)
@printf("  rho_c peaks (kHz, rel.power): %s\n",
        join([@sprintf("%.4f(%.2f)",f,p) for (f,p) in pk_rc2[1:min(5,end)]], "  "))
@printf("  probe peaks (kHz, rel.power): %s\n",
        join([@sprintf("%.4f(%.2f)",f,p) for (f,p) in pk_pr2[1:min(5,end)]], "  "))

# Identify which eigenmode each engine peak corresponds to
function classify(pk, targets)
    out=String[]
    for (f,p) in pk
        best=argmin([abs(f-t) for (_,t) in targets])
        lbl,tf=targets[best]; d=100*(f-tf)/tf
        push!(out, @sprintf("%s@%.4f(eng=%.4f, %+.1f%%, pow=%.2f)", lbl, tf, f, d, p))
    end
    return out
end
targets = [("F",F_eig),("H1",H1_eig),("H2",H2_eig)]

println("\n=== (2C) ENGINE PEAKS classified vs eigenvalues (F,H1,H2) ===")
println("  homol  rho_c : ", join(classify(pk_rc1[1:min(3,end)], targets), " | "))
println("  node   rho_c : ", join(classify(pk_rc2[1:min(3,end)], targets), " | "))
println("  node   probe : ", join(classify(pk_pr2[1:min(3,end)], targets), " | "))

# Does the node-bearing run show a peak near H1?
function nearest(pk, tf; tol=0.06)
    best=nothing; bd=Inf
    for (f,p) in pk
        d=abs(f-tf)/tf
        if d<bd; bd=d; best=(f,p,d); end
    end
    return (best !== nothing && best[3] < tol) ? best : nothing
end
h1_rc = nearest(pk_rc2, H1_eig); h1_pr = nearest(pk_pr2, H1_eig)
f_rc2 = nearest(pk_rc2, F_eig);  f_pr2 = nearest(pk_pr2, F_eig)

println("\n=== (2D) DOES THE ENGINE RESOLVE H1? ===")
@printf("  H1 eigenvalue = %.4f kHz\n", H1_eig)
if h1_rc !== nothing
    @printf("  rho_c: H1 peak at %.4f kHz (%+.2f%%, rel.power %.2f) -> RESOLVED\n",
            h1_rc[1], 100*(h1_rc[1]-H1_eig)/H1_eig, h1_rc[2])
else
    println("  rho_c: NO peak within 6% of H1")
end
if h1_pr !== nothing
    @printf("  probe: H1 peak at %.4f kHz (%+.2f%%, rel.power %.2f) -> RESOLVED\n",
            h1_pr[1], 100*(h1_pr[1]-H1_eig)/H1_eig, h1_pr[2])
else
    println("  probe: NO peak within 6% of H1")
end

println("\n# DONE")
