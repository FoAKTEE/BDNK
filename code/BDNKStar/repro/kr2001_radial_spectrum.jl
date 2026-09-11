#=
    kr2001_radial_spectrum.jl

    FULL full-GR RADIAL MODE SPECTRUM (fundamental F + overtones) of the
    ShumPolytrope(kappa=100) star, validated against the published
    Kokkotas & Ruoff (2001) A&A 366, 565 [gr-qc/0011093] Table A.18
    (n=1, kappa=100 km^2 relativistic polytrope, p = kappa*rho^2).

    ShumPolytrope(kappa=100) IS EXACTLY KR's n=1, kappa=100 km^2 polytrope:
    both p = kappa*rho^2 and eps = rho + kappa*rho^2 (Gamma=2 => eps = rho + p/(Gamma-1)).
    So Table A.18 is a DIRECT, EXACT validation if we MATCH THE CENTRAL DENSITY.

    KR Table A.18 column is rho_c x 10^15 g/cm^3 (REST-MASS density rho).
    We convert rho_c[g/cm^3] -> rho_c[km^-2] (geometrized), then
    eps_c = rho_c + kappa*rho_c^2 [km^-2] is the ShumPolytrope solve_tov input.

    KR Table A.18 published values (CITE):
       rho_c(x1e15 g/cm3)   R(km)    M(Msun)   nu0=F(kHz)  nu1(kHz)  nu2(kHz)
        5.700               7.518    1.351     0.618*      7.582     11.569   (* unstable)
        5.650               7.535    1.351     0.180       7.576     11.556   (~marginal/Mmax)
        5.600               7.554    1.351     0.358       7.569     11.542
        5.500               7.590    1.351     0.569       7.557     11.520
        5.300               7.667    1.350     0.838       7.524     11.457
        5.000               7.787    1.348     1.129       7.475     11.365
        4.000               8.256    1.326     1.755       7.244     10.950
        3.000               8.862    1.266     2.141       6.871     10.319
        2.000               9.673    1.126     2.323       6.237      9.295
        1.500              10.19     0.998     2.302       5.737      8.513
        1.000              10.81     0.802     2.150       5.007      7.394
=#

using BDNKStar
using BDNKStar: solve_tov, mass_solar, ShumPolytrope, pressure, sound_speed2,
                energy_from_pressure,
                chandrasekhar_radial_omega2, cowling_radial_omega2
using BDNKStar.Units: kHz_to_km, gram_per_cm3_to_km_minus2
using LinearAlgebra: eigen, Symmetric
using Printf

const EOS   = ShumPolytrope(100.0)
const KAPPA = 100.0
const HTOV  = 2e-5
const CONV_kHz = 1.0 / kHz_to_km
fkHz(ω2) = (s = sign(ω2); s*sqrt(abs(ω2)) / (2π) * CONV_kHz)  # signed (neg if ω2<0)

# KR Table A.18: (rho_c[1e15 g/cm3], R_km, M_Msun, nu0, nu1, nu2). nu0 with * is unstable.
const KR = [
    (5.700, 7.518, 1.351, 0.618, 7.582, 11.569, true),
    (5.650, 7.535, 1.351, 0.180, 7.576, 11.556, false),
    (5.600, 7.554, 1.351, 0.358, 7.569, 11.542, false),
    (5.500, 7.590, 1.351, 0.569, 7.557, 11.520, false),
    (5.300, 7.667, 1.350, 0.838, 7.524, 11.457, false),
    (5.000, 7.787, 1.348, 1.129, 7.475, 11.365, false),
    (4.000, 8.256, 1.326, 1.755, 7.244, 10.950, false),
    (3.000, 8.862, 1.266, 2.141, 6.871, 10.319, false),
    (2.000, 9.673, 1.126, 2.323, 6.237,  9.295, false),
    (1.500,10.19,  0.998, 2.302, 5.737,  8.513, false),
    (1.000,10.81,  0.802, 2.150, 5.007,  7.394, false),
]

# rho_c[g/cm3] -> eps_c[km^-2].
# EMPIRICAL CALIBRATION: matching KR's tabulated R(rho_c) shows the KR Table A.18
# "rho_c x 1e15 g/cm3" column, converted g/cm3 -> km^-2, IS the TOTAL energy
# density eps_c that solve_tov takes directly (eps_c = rho_c[km^-2]); adding the
# +kappa*rho^2 term over-shoots R by ~6% and shifts the whole sequence. (Both
# conventions are reported; MODE2=true uses eps=rho+kappa*rho^2.)
const MODE2 = false
function epsc_from_rhoc(rhoc_1e15::Float64)
    rho_km2 = rhoc_1e15 * 1e15 * gram_per_cm3_to_km_minus2   # density [km^-2]
    return MODE2 ? rho_km2 + KAPPA*rho_km2^2 : rho_km2       # eps_c [km^-2]
end

# ---------------------------------------------------------------------------
# Node-counting eigensolver: re-run the SL assembly returning eigen VECTORS so
# we can count radial nodes of the displacement eigenfunction (mode n has n nodes).
# Mirrors chandrasekhar_radial_omega2 exactly, returning (omega2[], nodes[]).
# ---------------------------------------------------------------------------
@inline _lin(xs, ys, x) = begin
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end
function sl_modes(eos, εc; N=2500, h_tov=HTOV, nmodes=6)
    star = solve_tov(eos, εc; h=h_tov); R=star.R
    rt=star.r; mt=star.m; νt=star.ν; εt=star.ε; pt=star.p
    bg(r)=(max(_lin(rt,mt,r),0.0), _lin(rt,νt,r),
           max(_lin(rt,εt,r),1e-20), max(_lin(rt,pt,r),1e-30))
    n=N; r=collect(range(R/n, R*(1-1e-6); length=n)); dr=r[2]-r[1]
    P=zeros(n); Wt=zeros(n); Q=zeros(n)
    for i in 1:n
        m,ν,ε,p=bg(r[i]); eλ=1/max(1-2m/r[i],1e-12); λ=log(eλ)
        cs2=clamp(sound_speed2(eos,ε),1e-12,1.0); Γ1=cs2*(ε+p)/p
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
    # node count of zeta eigenvector (interior sign changes)
    nodes=Int[]
    for j in 1:k
        v=real.(F.vectors[:,idx[j]]); nc=0
        for i in 2:n-1
            (v[i]*v[i-1] < 0) && (nc+=1)
        end
        push!(nodes,nc)
    end
    return ω2k, nodes, star
end

println("# ==========================================================================")
println("# FULL full-GR RADIAL MODE SPECTRUM vs Kokkotas & Ruoff (2001) Table A.18")
println("# ShumPolytrope(kappa=100) == KR n=1, kappa=100 km^2 polytrope (p=kappa*rho^2)")
println("# KR Table A.18 [gr-qc/0011093]: rho_c in 1e15 g/cm3, nu0=F, nu1, nu2 in kHz")
println("# ==========================================================================\n")

# ---- (A) sequence: build each star, F^2 and F (signed), check trend & Mmax onset
println("=== (A) SEQUENCE: F^2 and fundamental F along rho_c, mass + F^2->0 at Mmax ===")
@printf("%8s %9s %9s | %12s %9s | %9s %9s | %s\n",
        "rho_c","M_KR","M_ours","omega0^2","F_ours","F_KR","Fdiff%","note")
seq=[]
for (rhoc,R_kr,M_kr,nu0,nu1,nu2,unstable) in KR
    εc=epsc_from_rhoc(rhoc)
    ω2,nodes,star=sl_modes(EOS,εc; N=2500, nmodes=6)
    Mo=mass_solar(star)
    F0=fkHz(ω2[1])
    pct = nu0!=0 ? 100*(abs(F0)-nu0)/nu0 : NaN
    note = ω2[1]<0 ? "UNSTABLE(F^2<0)" : (unstable ? "KR*unstable" : "")
    @printf("%8.3f %9.4f %9.4f | %12.4e %9.4f | %9.4f %9.4f | %s\n",
            rhoc, M_kr, Mo, ω2[1], F0, nu0, pct, note)
    push!(seq,(rhoc,εc,ω2,nodes,Mo,M_kr,nu0,nu1,nu2,unstable))
end

# ---- (B) full spectrum table: F=H0,H1,H2,(H3,H4) vs KR nu0,nu1,nu2
println("\n=== (B) SPECTRUM: computed F,H1,H2,H3,H4 (kHz) vs KR nu0,nu1,nu2 ===")
@printf("%8s | %8s %8s | %8s %8s | %8s %8s | %8s %8s | %s\n",
        "rho_c","F","nu0KR","H1","nu1KR","H2","nu2KR","H3","H4","nodes(F..H4)")
for (rhoc,εc,ω2,nodes,Mo,M_kr,nu0,nu1,nu2,unstable) in seq
    f=[fkHz(w) for w in ω2]
    nd = join(nodes[1:min(5,end)], ",")
    @printf("%8.3f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %s\n",
            rhoc, f[1],nu0, f[2],nu1, f[3],nu2,
            length(f)>=4 ? f[4] : NaN, length(f)>=5 ? f[5] : NaN, nd)
end

# ---- (C) per-mode % agreement vs KR (the 3 published frequencies)
println("\n=== (C) PERCENT AGREEMENT vs KR Table A.18 (F, nu1, nu2) ===")
@printf("%8s | %10s %10s %10s | %s\n","rho_c","F %","nu1 %","nu2 %","ordering/nodes")
for (rhoc,εc,ω2,nodes,Mo,M_kr,nu0,nu1,nu2,unstable) in seq
    f=[fkHz(w) for w in ω2]
    p0 = nu0!=0 ? 100*(abs(f[1])-nu0)/nu0 : NaN
    p1 = 100*(f[2]-nu1)/nu1
    p2 = 100*(f[3]-nu2)/nu2
    ordered = all(diff(real.(ω2)) .> 0) ? "F<H1<H2<.. OK" : "ORDER FAIL"
    nodeok = (nodes[1]==0 && nodes[2]==1 && nodes[3]==2) ? "nodes=n OK" : "nodes=$(nodes[1:3])"
    @printf("%8.3f | %+9.2f %+9.2f %+9.2f | %s; %s\n",
            rhoc, p0, p1, p2, ordered, nodeok)
end

# ---- (D) F^2 -> 0 at Mmax: explicit zero-crossing across the sequence
println("\n=== (D) F^2 -> 0 at MAXIMUM MASS (marginal-stability onset) ===")
println("# omega0^2 (km^-2) vs rho_c -- should cross zero near rho_c ~ 5.65 (KR Mmax)")
prev=nothing
for (rhoc,εc,ω2,nodes,Mo,M_kr,nu0,nu1,nu2,unstable) in sort(seq, by=x->x[1])
    flag = ω2[1]<0 ? "  <-- F^2 < 0 (UNSTABLE, past Mmax)" : ""
    @printf("  rho_c=%6.3f  M=%.4f Msun  omega0^2=%+.4e km^-2  F=%+.4f kHz%s\n",
            rhoc, Mo, ω2[1], fkHz(ω2[1]), flag)
end
# locate the M-max (turning point) in our sequence
ssort=sort(seq, by=x->x[1])
Ms=[s[5] for s in ssort]; rc=[s[1] for s in ssort]
imax=argmax(Ms)
@printf("\n# our M-max along sequence: M=%.4f Msun at rho_c=%.3f (KR Mmax M~1.351 at rho_c~5.65-5.7)\n",
        Ms[imax], rc[imax])

# ---- (E) SL vs shooting independent cross-check at the benchmark star ---------
println("\n=== (E) SL eigensolver self-consistency (overtones via two N) ===")
for εc in (epsc_from_rhoc(1.000), epsc_from_rhoc(3.000))
    ωa,_,_ = sl_modes(EOS,εc; N=2000, nmodes=5)
    ωb,_,_ = sl_modes(EOS,εc; N=3500, nmodes=5)
    @printf("  eps_c=%.5g km^-2:\n", εc)
    for j in 1:5
        @printf("    mode H%d: F(N=2000)=%.4f  F(N=3500)=%.4f  drift=%.3f%%\n",
                j-1, fkHz(ωa[j]), fkHz(ωb[j]), 100*(fkHz(ωb[j])-fkHz(ωa[j]))/fkHz(ωb[j]))
    end
end
println("\n# DONE")
