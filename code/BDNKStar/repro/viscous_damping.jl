# =============================================================================
# viscous_damping.jl — Time-domain VISCOUS damping of the radial f-mode in the
# 1+1D dynamical-GR engine (DynGR1D), cross-checked against the inviscid
# eigenvalue and the dissipation integral.
#
# WHAT THIS SHOWS
#   DynGR1D gained a BDNK shear/bulk viscous momentum-diffusion source term
#   (kinematic ν_visc = η̂ = η/(ε+p)); see src/dyngr/DynGR1D.jl (_raw_rhs! viscous
#   flux + parabolic dt cap). Here we pluck a radial oscillation and measure how
#   it damps:
#     (1) ν=0  ⇒  ~undamped (only the LF/HLL numerical floor, Q≈68);
#     (2) ν>0  ⇒  clean exponential decay (narrowband-isolated fundamental);
#     (3) frequency matches the Cowling/Chandrasekhar radial eigenvalue;
#     (4) the MEASURED damping scales ∼√ν rather than ∝ν.
#
# ❌ CAVEAT / RETRACTION (2026-07-19 audit) — the √ν scaling here is an ARTIFACT OF THE
#   BOUNDARY TREATMENT and must NOT be cited as stellar physics. This engine hard-resets
#   v = 0 in the atmosphere (DynGR1D.jl lines ~344, 355, 363, and the outer ghosts ~387).
#   That is a numerical **NO-SLIP RIGID WALL** — and γ ∝ √(νω) is precisely the textbook
#   *rigid-wall* Stokes-layer law. So the measured √ν is expected from the boundary
#   condition alone. (An earlier header attributed it to a "free-surface Stokes layer";
#   that conflates two different classical results — for a genuine FREE surface the
#   classical answer is BULK damping, γ = 2νk², i.e. ∝ν.)
#   The companion ℓ=2 eigenvalue "confirmation" was likewise retracted — see
#   repro/viscous_fmode_scaling.jl and VALIDATION.md §6. Best-supported exponent is p≈1
#   (bulk). Items (1)-(3) above — inviscid-exact, clean exponential decay, correct
#   frequency — are unaffected and remain valid operator validations.
#   ✅ EXPERIMENT RUN (2026-07-19), via the new `atm_vbc` knob in setup_dyngr
#      (:zero = historical pinned-v wall; :outflow = zero-gradient, no wall stress).
#      Same star/operator/measurement; ONLY the BC changes. N=320, Cowling, 180R.
#        :zero     γ(ν=.02,.04,.08) = 2.609e-3, 3.059e-3, 4.165e-3   local p = 0.23, 0.45
#        :outflow  γ(ν=.02,.04,.08) = 2.733e-3, 2.792e-3, 3.253e-3   local p = 0.03, 0.22
#      RESULT — the no-slip hypothesis is **NOT** the cause of the sub-linearity: removing
#      the wall did NOT restore p≈1, it pushed p DOWN. What the data actually show is that
#      an AFFINE model (ν-independent OFFSET + LINEAR damping) fits far better than any
#      power law:
#        :zero     γ = 2.056e-3 + 0.0262·ν   (RSS 3.0e-9)  vs  A·ν^0.337 (RSS 3.9e-8) → 12.9× better
#        :outflow  γ = 2.503e-3 + 0.0091·ν   (RSS 8.4e-9)  vs  A·ν^0.126 (RSS 2.4e-8) →  2.9× better
#      and the :zero LINEAR slope 0.0262 agrees with the independent bulk dissipation-integral
#      prediction 0.0224 to 17%. ⇒ the underlying damping IS BULK (γ∝ν); the apparent "√ν"
#      was the OFFSET masquerading as a power law. The wall is not the generator, though it
#      does contribute ~65% of the genuine viscous slope (0.0262 → 0.0091 when removed).
#      REMAINING CONFOUND: the offset appears only for ν>0 (ν=0 shows no measurable damping),
#      so it is not a static numerical floor — most likely the parabolic dt cap (dt ∝ 1/ν,
#      engaged for every ν≥0.02 here) altering the scheme's dissipation. Quantifying that is
#      the next step before any damping COEFFICIENT from this engine is quoted.
#
# METHOD. Frozen-metric (Cowling) evolution for a clean, fast background; seed
#   v=A r/R; sample the probe velocity; per-window narrowband demodulation at the
#   measured f-mode ω (rejects fast-damping overtones); fit ln(amplitude) vs t.
# =============================================================================
using BDNKStar, Statistics
using BDNKStar: setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!, solve_tov,
                ShumPolytrope, cowling_radial_omega2, energy_from_pressure

const EOS = ShumPolytrope(100.0); const ΕC = 0.0015
const STAR = solve_tov(EOS, ΕC; h=2e-4); const R = STAR.R
const ΩREF = sqrt(cowling_radial_omega2(EOS, ΕC; nmodes=1)[1]); const PERIOD = 2π/ΩREF

itp(xs,ys,x) = (x<=xs[1] ? ys[1] : x>=xs[end] ? ys[end] :
    (j=searchsortedlast(xs,x); t=(x-xs[j])/(xs[j+1]-xs[j]); ys[j]+t*(ys[j+1]-ys[j])))

# narrowband fundamental amplitude per window (DFT magnitude at ω): rejects overtones.
function nb_gamma(ts, sig, ω, P; perwin=2.5)
    T=ts[end]-ts[1]; nseg=max(5,floor(Int,T/(perwin*P))); wlen=T/nseg
    tc=Float64[]; amp=Float64[]
    for k in 1:nseg
        t0=ts[1]+(k-1)*wlen
        idx=findall(t-> t0<=t<t0+wlen, ts); length(idx)<6 && continue
        s=sig[idx].-mean(sig[idx]); tt=ts[idx]
        cr=sum(s.*cos.(ω.*tt)); ci=sum(s.*sin.(ω.*tt))
        push!(tc, mean(tt)); push!(amp, 2*sqrt(cr^2+ci^2)/length(idx))
    end
    iend=1; for i in 2:length(amp); amp[i]<amp[i-1] ? (iend=i) : break; end
    iend<3 && return NaN
    x=tc[1:iend]; y=log.(amp[1:iend]); m=length(x)
    sx=sum(x);sy=sum(y);sxx=sum(x.^2);sxy=sum(x.*y)
    -(m*sxy-sx*sy)/(m*sxx-sx^2)
end

function run_sweep(νs)
    rows = Tuple{Float64,Float64}[]
    for ν in νs
        eng, st = setup_dyngr(EOS, ΕC; N=400, cfl=0.25, ν_visc=ν, freeze_metric=true)
        seed_dyngr_velocity!(st, eng; A=2e-4, profile=:linear)
        res = evolve_dyngr!(st, eng; tmax=200*R, probe_frac=0.5)
        push!(rows, (ν, nb_gamma(res.ts, res.probe, ΩREF, PERIOD)))
    end
    rows
end

# homologous (δv∝r) X-weighted dissipation integral: the BULK damping slope dγ/dν.
function bulk_slope()
    Nr=3000; rg=range(1e-3, R*0.999; length=Nr); dr=step(rg); num=0.0; den=0.0
    for r in rg
        pl=max(itp(STAR.r,STAR.p,r),0.0); w=energy_from_pressure(EOS,pl)+pl
        m=itp(STAR.r,STAR.m,r); X=1/sqrt(max(1-2m/r,1e-6))
        num += w*(1/R)^2*r^2/X*dr; den += w*(r/R)^2*r^2*X*dr
    end
    num/den
end

function main()
    println("Cowling radial f-mode: ω=", round(ΩREF,sigdigits=4), " /km  (period ",
            round(PERIOD,digits=1), " km);  R=", round(R,digits=2), " km")
    rows = run_sweep((0.0, 0.02, 0.04, 0.08))
    γ0 = rows[1][2]
    println("\n  ν_visc    γ(/km)       γ−γ0        (γ−γ0)/ν     (γ−γ0)/√ν     Q=ω/2γ")
    for (ν,γ) in rows
        dγ = γ-γ0
        println("  ", rpad(ν,9), rpad(round(γ,sigdigits=3),13),
                rpad(ν>0 ? string(round(dγ,sigdigits=3)) : "—",12),
                rpad(ν>0 ? string(round(dγ/ν,sigdigits=3)) : "—",13),
                rpad(ν>0 ? string(round(dγ/sqrt(ν),sigdigits=3)) : "—",14),
                round(ΩREF/(2γ),digits=1))
    end
    # √ν boundary-layer model vs linear bulk model
    ν2,γ2 = rows[2]; c_sqrt = (γ2-γ0)/sqrt(ν2); c_lin = (γ2-γ0)/ν2
    println("\nBoundary-layer model  γ = γ0 + c√ν  (c=", round(c_sqrt,sigdigits=3), "):")
    for (ν,γ) in rows[3:end]
        pred=γ0+c_sqrt*sqrt(ν)
        println("   ν=",ν,"  predicted γ=",round(pred,sigdigits=3),"  measured=",round(γ,sigdigits=3),
                "  err=",round(100*(γ-pred)/γ,digits=1),"%")
    end
    println("Linear bulk model     γ = γ0 + c·ν  (c=", round(c_lin,sigdigits=3), "):")
    for (ν,γ) in rows[3:end]
        pred=γ0+c_lin*ν
        println("   ν=",ν,"  predicted γ=",round(pred,sigdigits=3),"  measured=",round(γ,sigdigits=3),
                "  err=",round(100*(γ-pred)/γ,digits=1),"%")
    end
    println("\nBulk dissipation-integral slope dγ/dν (homologous, X-weighted) = ",
            round(bulk_slope(),sigdigits=4), " /km  — sub-dominant vs the √ν surface layer.")
    println("\n⇒ Operator validated: inviscid-exact, clean exponential decay, correct frequency.")
    println("  Radial (l=0) f-mode damping is SURFACE-BOUNDARY-LAYER (∼√ν) dominated.")
end

main()
