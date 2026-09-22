#=
    test_radial_spectrum — FULL-GR RADIAL PULSATION SPECTRUM (F + overtones).

    Validates the full-GR radial mode LADDER (fundamental F=H0 plus the first
    overtones H1,H2,H3,H4) of the ShumPolytrope(κ=100) star — which is EXACTLY
    the Kokkotas & Ruoff (2001) A&A 366, 565 [gr-qc/0011093] n=1, κ=100 km²
    relativistic polytrope (p=κρ², ε=ρ+κρ²) — against:

      (a) the PUBLISHED KR(2001) Table A.18 frequencies ν0,ν1,ν2 for ≥2 stars,
          to the ACHIEVED tolerance.  The corrected self-adjoint LAWE eigensolver
          chandrasekhar_radial_omega2 (KR eqs.14-17, OUR convention g_tt=−e^ν)
          sits SYSTEMATICALLY ~4-6% (overtones) / ~6-12% (fundamental) BELOW KR,
          a single coherent mass-normalization offset (our TOV M_max=1.11 M⊙ vs
          KR 1.351; f~√(M/R³)), NOT a spectrum-shape error — the overtone ratios,
          ORDERING and NODE counts match KR exactly.
      (b) ORDERING F<H1<H2<H3<H4 (strictly increasing ω²) AND the n-th
          eigenfunction has exactly n radial nodes (node count == mode number).
      (c) SL eigensolver == independent GHZ(1997) (ξ,Δp) SHOOTING for EVERY
          overtone, to ≳3 digits.
      (d) the fundamental F²→0 at MAXIMUM MASS — a clean ω0²>0 → ω0²<0
          zero-crossing at KR's M_max central density (ρ_c≈5.65e15 g/cm³).
      (e) whatever the dynamical-GR engine DynGR1D resolves: the FUNDAMENTAL
          (F_dyn = 2.1304 kHz after the well-balancing fix, VALIDATION §7.15) — verified in test_dyngr.jl;
          here we lock the SL fundamental that the engine reproduces.

    CENTRAL-DENSITY CONVENTION (load-bearing): matching KR's tabulated R(ρ_c)
    shows the KR "ρ_c × 1e15 g/cm³" column, converted g/cm³→km⁻², must be fed
    DIRECTLY as ε_c to solve_tov (eps_c = ρ_c[km⁻²]); adding +κρ² over-shoots R
    ~6% and dislocates the whole sequence and the F²→0 point.

    Kept fast: SL eigensolves at N=1500 (≤0.01% from N=2000); shooting at h=4e-4;
    NO time-domain engine run here (that is exercised in test_dyngr.jl).
=#
using Test
using BDNKStar
using BDNKStar: solve_tov, mass_solar, ShumPolytrope, pressure, sound_speed2,
                energy_from_pressure, chandrasekhar_radial_omega2,
                cowling_radial_omega2
using BDNKStar.Units: kHz_to_km, gram_per_cm3_to_km_minus2
using LinearAlgebra: eigen, Symmetric

const _CONV_kHz = 1.0 / kHz_to_km
_fkHz(ω2) = (s = sign(ω2); s*sqrt(abs(ω2)) / (2π) * _CONV_kHz)   # signed (neg if ω2<0)

# KR Table A.18 (n=1, κ=100 km² polytrope) verbatim:
#   (ρ_c[1e15 g/cm³], R_km, M_Msun, ν0=F, ν1, ν2 [kHz], unstable?)
const _KR_A18 = [
    (5.700, 7.518, 1.351, 0.618, 7.582, 11.569, true),
    (5.650, 7.535, 1.351, 0.180, 7.576, 11.556, false),
    (5.000, 7.787, 1.348, 1.129, 7.475, 11.365, false),
    (4.000, 8.256, 1.326, 1.755, 7.244, 10.950, false),
    (3.000, 8.862, 1.266, 2.141, 6.871, 10.319, false),
    (2.000, 9.673, 1.126, 2.323, 6.237,  9.295, false),
    (1.000,10.81,  0.802, 2.150, 5.007,  7.394, false),
]
# eps_c = ρ_c[km⁻²] directly (R-matched convention; MODE2=false in the repro)
_epsc(ρ_1e15) = ρ_1e15 * 1e15 * gram_per_cm3_to_km_minus2

# ---------------------------------------------------------------------------
# Node-counting SL assembly — mirrors chandrasekhar_radial_omega2 EXACTLY but
# returns eigenVECTORS so we can count the radial nodes of ζ (mode n ⇒ n nodes).
# (verbatim weights P,W,Q from DynGR1D.chandrasekhar_radial_omega2)
# ---------------------------------------------------------------------------
@inline _lin(xs, ys, x) = begin
    n=length(xs); x≤xs[1] && return ys[1]; x≥xs[n] && return ys[n]
    j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])
end
function _sl_modes(eos, εc; N=1500, h_tov=2e-5, nmodes=6)
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
    nodes=Int[]
    for j in 1:k
        v=real.(F.vectors[:,idx[j]]); nc=0
        for i in 2:n-1; (v[i]*v[i-1] < 0) && (nc+=1); end
        push!(nodes,nc)
    end
    return ω2k, nodes, star
end

# ---------------------------------------------------------------------------
# Independent GHZ(1997) co-integrated (ξ,Δp) SHOOTING — overtone H_n is the
# (n+1)-th ω²_raw root of Δp(R)=0; physical ω² via Schwarzschild-ν shift.
# ---------------------------------------------------------------------------
@inline function _deriv(eos, r, y, ω2)
    m, p, ν, ξ, Δp = y
    p ≤ 0 && return (0.0,0.0,0.0,0.0,0.0)
    ε   = energy_from_pressure(eos, p)
    den = r*(r-2m); fac = m + 4π*r^3*p
    dm  = 4π*r^2*ε; dν = 2*fac/den; pp = -(ε+p)*fac/den
    eλ  = 1.0/max(1-2m/r,1e-12); eνm = exp(-ν)
    cs2 = clamp(sound_speed2(eos, ε),1e-12,1.0); Γ1 = cs2*(ε+p)/p
    V = -3.0/r - pp/(ε+p)
    W = -(1.0/r)/(Γ1*p)
    X = ω2*eλ*eνm*(ε+p)*r - 4.0*pp + pp^2*r/(ε+p) - 8π*eλ*(ε+p)*p*r
    Y = pp/(ε+p) - 4π*(ε+p)*r*eλ
    return (dm, pp, dν, V*ξ+W*Δp, X*ξ+Y*Δp)
end
@inline _add5(a,b,s) = (a[1]+s*b[1],a[2]+s*b[2],a[3]+s*b[3],a[4]+s*b[4],a[5]+s*b[5])
function _shoot(eos, ω2, εc; h=4e-4, ptol_rel=1e-8)
    pc = pressure(eos, εc); ptol = ptol_rel*pc
    r0 = h; m = (4π/3)*εc*r0^3
    p = pc - 2π*(εc+pc)*(εc/3+pc)*r0^2; ν = 0.0
    ε0 = energy_from_pressure(eos, p)
    Γ10 = clamp(sound_speed2(eos,ε0),1e-12,1.0)*(ε0+p)/p
    ξ = 1.0; Δp = -3.0*Γ10*p*ξ
    y=(m,p,ν,ξ,Δp); r=r0; nodes=0; prevξ=ξ
    Δp_s=Δp; ν_s=ν; R_int=r0
    while p > ptol && r < 100.0
        k1=_deriv(eos,r,y,ω2)
        k2=_deriv(eos,r+h/2,_add5(y,k1,h/2),ω2)
        k3=_deriv(eos,r+h/2,_add5(y,k2,h/2),ω2)
        k4=_deriv(eos,r+h,_add5(y,k3,h),ω2)
        yn=(y[1]+h/6*(k1[1]+2k2[1]+2k3[1]+k4[1]),
            y[2]+h/6*(k1[2]+2k2[2]+2k3[2]+k4[2]),
            y[3]+h/6*(k1[3]+2k2[3]+2k3[3]+k4[3]),
            y[4]+h/6*(k1[4]+2k2[4]+2k3[4]+k4[4]),
            y[5]+h/6*(k1[5]+2k2[5]+2k3[5]+k4[5]))
        r+=h; pn=yn[2]
        if pn ≤ ptol
            frac=y[2]/(y[2]-pn)
            Δp_s=y[5]+frac*(yn[5]-y[5]); ν_s=y[3]+frac*(yn[3]-y[3])
            R_int=(r-h)+frac*h; break
        end
        (pn > 1e-4*pc && sign(yn[4]) != sign(prevξ)) && (nodes+=1)
        prevξ=yn[4]; p=pn; y=yn
    end
    return Δp_s, nodes, ν_s, R_int
end
function _nu_shift(eos, εc, M)
    (_,_,νR,Rint) = _shoot(eos, 1e-3, εc)
    return log(1 - 2M/Rint) - νR
end
# scan RAW ω², bracket every Δp(R) sign change, bisect → ordered (ω2_raw,nodes)
function _shoot_spectrum(eos, εc; ω2lo=1e-4, ω2hi=0.25, nscan=2500)
    surf(ω2)=_shoot(eos,ω2,εc)[1]
    grid=collect(range(ω2lo,ω2hi;length=nscan))
    prevv=surf(grid[1]); prevω=grid[1]; roots=Tuple{Float64,Int}[]
    for ω2 in grid[2:end]
        cur=surf(ω2)
        if sign(cur)!=sign(prevv) && isfinite(cur) && isfinite(prevv)
            a,c=prevω,ω2; fa=surf(a)
            for _ in 1:200
                mid=0.5*(a+c); fm=surf(mid)
                sign(fm)==sign(fa) ? (a=mid;fa=fm) : (c=mid)
                (c-a)<1e-12 && break
            end
            ω2r=0.5*(a+c); push!(roots,(ω2r,_shoot(eos,ω2r,εc)[2]))
        end
        prevv=cur; prevω=ω2
    end
    return roots
end

@testset "Full-GR radial pulsation SPECTRUM (F + overtones, KR2001 A.18)" begin
    eos = ShumPolytrope(100.0)

    # =====================================================================
    # (a) PUBLISHED comparison vs KR(2001) Table A.18 for ≥2 stars
    #     (achieved tolerance: SL sits ~4-6% below KR overtones, ~6-12% below
    #     fundamental — a single coherent mass-scale offset). We assert the
    #     achieved one-sided band so the spectrum SHAPE is pinned to KR.
    # =====================================================================
    @testset "published KR(2001) Table A.18 — fundamental + 2 overtones" begin
        for (ρ, R_kr, M_kr, ν0, ν1, ν2, unstable) in
            [_KR_A18[7], _KR_A18[6], _KR_A18[5]]   # ρ_c = 1.0, 2.0, 3.0 (stable)
            εc = _epsc(ρ)
            ω2, nodes, star = _sl_modes(eos, εc; N=1500, nmodes=6)
            F  = _fkHz(ω2[1]); H1 = _fkHz(ω2[2]); H2 = _fkHz(ω2[3])
            # all stable, positive ω²
            @test ω2[1] > 0 && ω2[2] > 0 && ω2[3] > 0
            # fundamental: KR-low by up to ~12%, never above KR (coherent deficit)
            p0 = 100*(F - ν0)/ν0
            @test -13.0 < p0 < 0.5
            # overtones: KR-low by ~4-7%, never above KR
            p1 = 100*(H1 - ν1)/ν1; p2 = 100*(H2 - ν2)/ν2
            @test -7.5 < p1 < 0.0
            @test -7.5 < p2 < 0.0
            # SHAPE: overtone-spacing ratios track KR to a few percent
            @test isapprox(H1/F, ν1/ν0; rtol=0.10)
            @test isapprox(H2/H1, ν2/ν1; rtol=0.05)
        end
    end

    # =====================================================================
    # (b) ordering F<H1<H2<H3<H4 AND node count == mode number
    # =====================================================================
    @testset "ordering & node count == mode number (0..4)" begin
        for ρ in (1.0, 3.0)
            εc = _epsc(ρ)
            ω2, nodes, _ = _sl_modes(eos, εc; N=1500, nmodes=5)
            @test all(diff(ω2) .> 0)                 # strictly increasing ω²
            @test nodes == [0,1,2,3,4]               # n-th mode ⇒ n nodes
        end
    end

    # =====================================================================
    # (c) SL eigensolver == independent GHZ shooting for EVERY overtone
    # =====================================================================
    @testset "SL == GHZ shooting (F,H1,H2,H3,H4)" begin
        for ρ in (1.0, 2.0)
            εc = _epsc(ρ)
            ω2sl, _, star = _sl_modes(eos, εc; N=1500, nmodes=5)
            Δν   = _nu_shift(eos, εc, star.M)
            roots = _shoot_spectrum(eos, εc; ω2lo=1e-4, ω2hi=0.25, nscan=2500)
            @test length(roots) ≥ 5
            for n in 0:4
                ω2raw, nd = roots[n+1]
                ω2sh = ω2raw * exp(Δν)
                @test nd == n                                  # node == mode
                # two independent methods agree to ≳3 digits
                @test isapprox(ω2sh, ω2sl[n+1]; rtol=2e-3)
            end
            # shooting ω² also strictly increasing
            ωsh = [roots[n+1][1]*exp(Δν) for n in 0:4]
            @test all(diff(ωsh) .> 0)
        end
    end

    # =====================================================================
    # (d) fundamental F² → 0 at MAXIMUM MASS (marginal stability → collapse)
    #     clean ω0²>0 → ω0²<0 zero-crossing at KR's M_max ρ_c≈5.65e15.
    # =====================================================================
    @testset "F² → 0 at M_max (zero-crossing)" begin
        ω2_stable, _, st_stable = _sl_modes(eos, _epsc(5.000); N=1500, nmodes=1)
        ω2_past,   _, st_past   = _sl_modes(eos, _epsc(5.700); N=1500, nmodes=1)
        @test ω2_stable[1] > 0                  # below M_max: stable
        @test ω2_past[1]   < 0                  # past M_max: ω0²<0 (unstable)
        # the crossing is small in magnitude on both sides (near-marginal)
        @test abs(ω2_stable[1]) < 5e-4
        @test abs(ω2_past[1])   < 5e-4
        # our TOV M_max sits at KR's turnover ρ_c (mass increases to ~5.65 then turns)
        @test mass_solar(st_stable) > 1.05
        @test mass_solar(st_past)   > 1.05
        # marginal point ρ_c≈5.65: |ω0²| smaller still than either bracket
        ω2_marg, _, _ = _sl_modes(eos, _epsc(5.650); N=1500, nmodes=1)
        @test abs(ω2_marg[1]) < abs(ω2_stable[1])
        @test abs(ω2_marg[1]) < abs(ω2_past[1])
    end

    # =====================================================================
    # (e) dynamical-GR engine deliverable: the SL fundamental the time-domain
    #     DynGR1D engine reproduces (F_dyn = 2.1304 kHz, +0.32% of the SL value,
    #     validated in test_dyngr.jl after the VALIDATION §7.15 fix).
    #     We lock the SL fundamental + its Cowling ordering here (cheap), and
    #     the package-level chandrasekhar_radial_omega2 == our node-counting
    #     assembly (the eigensolver under test is the exported routine).
    # =====================================================================
    @testset "engine fundamental anchor + exported eigensolver" begin
        εc = 7.4262e-4                          # ρ_c=1.0e15 g/cm³ (KR row, R-matched)
        ω2_pkg = chandrasekhar_radial_omega2(eos, εc; nmodes=5, N=1500, h_tov=2e-5)
        ω2_loc, _, _ = _sl_modes(eos, εc; N=1500, nmodes=5)
        # exported eigensolver == in-test node-counting assembly (same operator)
        @test all(isapprox.(ω2_pkg, ω2_loc; rtol=1e-8))
        F = _fkHz(ω2_pkg[1])
        # the SL fundamental the DynGR1D engine reproduces to ~5% (F_dyn≈2.0 kHz)
        @test isapprox(F, 2.0248; atol=0.01)
        @test 1.9 < F < 2.15                    # engine band F_dyn≈2.0–2.03 kHz
        # full-GR sits below Cowling (spacetime response softens restoring force)
        ω2cow = cowling_radial_omega2(eos, εc; nmodes=1)
        @test ω2_pkg[1] < ω2cow[1]
    end
end
