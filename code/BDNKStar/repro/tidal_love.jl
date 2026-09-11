#=
    tidal_love.jl — STATIC ℓ=2 TIDAL DEFORMABILITY (the I-Love-Q completion).

    Computes the ℓ=2 tidal Love number k₂ and dimensionless tidal deformability
    Λ = (2/3) k₂ (R/M)^5 = (2/3) k₂ C^{-5} for the realistic Read et al. (2009)
    piecewise-polytrope EOS (SLy/APR4/H4/MS1), and combines with the slow-rotation
    moment of inertia already in the code to test the I-Love universal relation.

    METHOD (source-verified — Postnikov, Prakash & Lattimer 2010, PRD 82, 024016
    [arXiv:1004.5098], Eqs. (eqH),(eq:Q),(eq:met),(eq:y); k₂ closed form Hinderer
    2008 ApJ 677,1216 + erratum):

      Integrate the first-order Riccati ODE for y(r)=rH'/H from the centre y(0)=2:
        r y' + y² + y·F(r) + r²·Q(r) = 0,
      with   e^λ = 1/(1−2m/r),
        F(r) = e^λ [1 + 4π r²(p−ε)],
        Q(r) = 4π e^λ (5ε + 9p + (ε+p)/c_s²) − 6 e^λ/r² − (ν')²,
        ν'  = 2 e^λ (m + 4π r³ p)/r²,     c_s² = dp/dε.
      Surface y_R = y(R). With C = M/R,

        k₂ = (8C⁵/5)(1−2C)²[2 + 2C(y_R−1) − y_R]
             / { 2C[6 − 3y_R + 3C(5y_R−8)]
               + 4C³[13 − 11y_R + C(3y_R−2) + 2C²(1+y_R)]
               + 3(1−2C)²[2 − y_R + 2C(y_R−1)] ln(1−2C) },
        Λ = (2/3) k₂ C^{-5}.

    The c_s² for the piecewise polytrope is the EOS sound_speed2(ε)=dp/dε; we
    re-derive ε,p,c_s² from the EOS at each grid point (the TOV grid is reused).
=#
using BDNKStar
using BDNKStar: piecewise_polytrope, solve_tov, sound_speed2, pressure,
                energy_from_pressure, mass_solar
using BDNKStar.SlowRotation: moment_of_inertia, SlowRotResult

# ── k₂ closed form (Hinderer 2008 + erratum; matches PPL2010 eq:k2C) ──────────
function k2_of(C::Float64, y::Float64)
    num = (8C^5/5)*(1-2C)^2*(2 + 2C*(y-1) - y)
    den = 2C*(6 - 3y + 3C*(5y-8)) +
          4C^3*(13 - 11y + C*(3y-2) + 2C^2*(1+y)) +
          3*(1-2C)^2*(2 - y + 2C*(y-1))*log(1-2C)
    return num/den
end

# ── y(r) integration on the TOV grid (RK4, midpoint EOS lookups) ─────────────
# RHS:  y' = -[ y² + y F + r² Q ] / r
function tidal_yR(eos, star)
    r = star.r; m = star.m
    # helper: F, Q at a radius given (r,m,p,ε,cs2)
    @inline function FQ(rr, mm, p, ε, cs2)
        eλ = 1/(1 - 2mm/rr)
        F  = eλ*(1 + 4π*rr^2*(p - ε))
        νp = 2*eλ*(mm + 4π*rr^3*p)/rr^2
        Q  = 4π*eλ*(5ε + 9p + (ε+p)/cs2) - 6*eλ/rr^2 - νp^2
        return F, Q
    end
    @inline function dydr(rr, yy, mm, p, ε, cs2)
        F, Q = FQ(rr, mm, p, ε, cs2)
        return -(yy^2 + yy*F + rr^2*Q)/rr
    end
    # state at a grid point i: (p,ε,cs2) from the EOS at ε[i]
    @inline function state(i)
        εi = star.ε[i]; pi = star.p[i]
        cs2 = εi > 0 ? sound_speed2(eos, εi) : 1.0
        return pi, εi, cs2
    end
    # interpolated midpoint state between i and i+1 (linear in r on m,p,ε)
    @inline function state_mid(i)
        εi = 0.5*(star.ε[i] + star.ε[i+1])
        pmid = εi > 0 ? pressure(eos, εi) : 0.0
        cs2  = εi > 0 ? sound_speed2(eos, εi) : 1.0
        mmid = 0.5*(m[i] + m[i+1])
        return mmid, pmid, εi, cs2
    end

    y = 2.0
    n = length(r)
    # stop a hair inside the surface to avoid p→0,ε→0 (c_s² ill-defined at edge);
    # integrate up to the last node with ε>0, then carry y to R.
    ilast = n
    while ilast > 1 && star.ε[ilast] ≤ 0
        ilast -= 1
    end
    @inbounds for i in 1:ilast-1
        hh = r[i+1] - r[i]
        p0,ε0,c0 = state(i)
        k1 = dydr(r[i], y, m[i], p0, ε0, c0)
        mm,pm,εm,cm = state_mid(i)
        rm = r[i] + hh/2
        k2 = dydr(rm, y + hh/2*k1, mm, pm, εm, cm)
        k3 = dydr(rm, y + hh/2*k2, mm, pm, εm, cm)
        p1,ε1,c1 = state(i+1)
        k4 = dydr(r[i+1], y + hh*k3, m[i+1], p1, ε1, c1)
        y += hh/6*(k1 + 2k2 + 2k3 + k4)
    end
    # carry y across the (tiny) crust edge to the surface with the vacuum eq.
    # (p=ε=0 ⇒ F=1, Q=-(ν')²): y' = -(y²+y)/r + r(ν')² ; over a thin shell this is
    # a negligible correction, but we include it for completeness.
    @inbounds for i in ilast:n-1
        hh = r[i+1] - r[i]
        eλ = 1/(1 - 2m[i]/r[i])
        νp = 2*eλ*m[i]/r[i]^2
        y += hh*(-(y^2 + y*eλ)/r[i] + r[i]*νp^2)
    end
    return y
end

# central energy density (km⁻²) tuned to a target gravitational mass (M⊙)
function εc_for_mass(eos, Mtarget_Msun; εlo=3e-4, εhi=4e-3)
    f(εc) = mass_solar(solve_tov(eos, εc; h=2e-3)) - Mtarget_Msun
    a, b = εlo, εhi
    fa, fb = f(a), f(b)
    # expand bracket if needed
    it = 0
    while fa*fb > 0 && it < 40
        b *= 1.15; fb = f(b); it += 1
    end
    for _ in 1:80
        c = 0.5*(a+b); fc = f(c)
        if fa*fc ≤ 0; b = c; fb = fc; else; a = c; fa = fc; end
        abs(b-a) < 1e-8 && break
    end
    return 0.5*(a+b)
end

# ── Yagi–Yunes (2013) I–Love fit:  ln Ī = Σ coeff (ln Λ)^k  (Table I, Ī–Λ̄) ──
const YY_ILOVE = (1.47, 0.0817, 0.0149, 2.87e-4, -3.64e-5)
function ibar_yy(Λ)
    x = log(Λ)
    a,b,c,d,e = YY_ILOVE
    return exp(a + b*x + c*x^2 + d*x^3 + e*x^4)
end

# ── Published Λ(1.4) references (CITED in report) ────────────────────────────
const PUB_L14 = Dict(:SLy=>300.0, :APR4=>245.0, :H4=>900.0, :MS1=>1250.0)

const EOSES = (:SLy, :APR4, :H4, :MS1)
const h_fine = 1e-3

println("# ℓ=2 tidal deformability — Read et al. (2009) piecewise polytropes")
println("# method: Postnikov-Prakash-Lattimer 2010 y-ODE + Hinderer 2008 k₂ closed form\n")

# ---- (1) Λ(1.4) table vs published ----
println("## Λ(1.4 M⊙) table")
println(rpad("EOS",6), rpad("M[M⊙]",9), rpad("R[km]",9), rpad("C",9),
        rpad("k2",9), rpad("Λ",11), rpad("Λ_pub",9), "Δ%")
results_14 = Dict{Symbol,Any}()
for eos_sym in EOSES
    eos = piecewise_polytrope(eos_sym; N=4000)
    εc  = εc_for_mass(eos, 1.4)
    star = solve_tov(eos, εc; h=h_fine)
    C = star.M/star.R
    yR = tidal_yR(eos, star)
    k2 = k2_of(C, yR)
    Λ  = (2/3)*k2*C^(-5)
    pub = PUB_L14[eos_sym]
    dpc = 100*(Λ-pub)/pub
    results_14[eos_sym] = (M=mass_solar(star), R=star.R, C=C, k2=k2, Λ=Λ, yR=yR)
    println(rpad(string(eos_sym),6),
            rpad(round(mass_solar(star),digits=3),9),
            rpad(round(star.R,digits=3),9),
            rpad(round(C,digits=4),9),
            rpad(round(k2,digits=4),9),
            rpad(round(Λ,digits=2),11),
            rpad(round(pub,digits=1),9),
            round(dpc,digits=1))
end

# ---- (2) I-Love universal relation across the mass sequence ----
println("\n## I-Love (Yagi-Yunes 2013) — Ī vs Λ along the mass sequence")
println(rpad("EOS",6), rpad("M[M⊙]",9), rpad("Λ",11), rpad("Ibar",9),
        rpad("Ibar_YY",10), "Δ%")
iloveq = Tuple[]
for eos_sym in EOSES
    eos = piecewise_polytrope(eos_sym; N=4000)
    for Mt in (1.2, 1.4, 1.6, 1.8)
        εc = εc_for_mass(eos, Mt)
        star = solve_tov(eos, εc; h=h_fine)
        C = star.M/star.R
        yR = tidal_yR(eos, star)
        k2 = k2_of(C, yR)
        Λ  = (2/3)*k2*C^(-5)
        rot = moment_of_inertia(star)
        Ibar = rot.Ibar
        Iyy  = ibar_yy(Λ)
        dpc  = 100*(Ibar-Iyy)/Iyy
        push!(iloveq, (eos_sym, mass_solar(star), Λ, Ibar, Iyy, dpc))
        println(rpad(string(eos_sym),6),
                rpad(round(mass_solar(star),digits=3),9),
                rpad(round(Λ,digits=2),11),
                rpad(round(Ibar,digits=4),9),
                rpad(round(Iyy,digits=4),10),
                round(dpc,digits=2))
    end
end
maxdev = maximum(abs(t[6]) for t in iloveq)
println("# max |Δ%| vs Yagi-Yunes I-Love fit = ", round(maxdev,digits=2), "%")

# ---- (3) GW170817 confrontation ----
println("\n## GW170817 confrontation (Λ_1.4 = 190 +390/−120; 90% upper ≲ 800)")
for eos_sym in EOSES
    Λ = results_14[eos_sym].Λ
    status = Λ ≤ 800 ? "CONSISTENT (≤800)" : "DISFAVORED (>800)"
    inband = (70 ≤ Λ ≤ 580) ? "within 90% band [70,580]" : "outside 90% band"
    println(rpad(string(eos_sym),6), "Λ=", rpad(round(Λ,digits=1),9),
            status, " | ", inband)
end
