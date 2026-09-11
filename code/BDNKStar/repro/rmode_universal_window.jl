using BDNKStar
using BDNKStar: piecewise_polytrope, solve_tov, mass_solar
using BDNKStar.RModes: rmode_instability_window, kepler_frequency

# central-density bisection to hit a target gravitational mass [M⊙]
function star_for_mass(eos, Mtarget; εlo=2e-4, εhi=4e-3, tol=1e-4)
    flo = mass_solar(solve_tov(eos, εlo)) - Mtarget
    fhi = mass_solar(solve_tov(eos, εhi)) - Mtarget
    if flo*fhi > 0
        return nothing
    end
    a, b = εlo, εhi
    local star
    for _ in 1:60
        εm = 0.5*(a+b)
        star = solve_tov(eos, εm)
        fm = mass_solar(star) - Mtarget
        (flo*fm ≤ 0) ? (b=εm; fhi=fm) : (a=εm; flo=fm)
        abs(fm) < tol && break
    end
    return star
end

eoses = [(:SLy, piecewise_polytrope(:SLy)),
         (:APR4, piecewise_polytrope(:APR4)),
         (:H4,   piecewise_polytrope(:H4))]

masses = [1.2, 1.4, 1.6]
# core-temperature grid spanning the LOM98 instability window (shear-low to bulk-high edge)
Tgrid = 10 .^ range(6.0, 11.0; length=400)

println("eos,M_Msun,R_km,compactness,OmegaK_s,nuK_Hz,minRatio,nu_c_Hz,T_at_min_K")
for (name, eos) in eoses
    for M in masses
        star = star_for_mass(eos, M)
        star === nothing && (println("$name,$M,NA,NA,NA,NA,NA,NA,NA"); continue)
        R = star.R
        C = star.M / R                      # GM/(Rc^2), geometric
        ΩK = kepler_frequency(star, eos)
        w  = rmode_instability_window(star, eos; Tgrid=Tgrid)
        imin = argmin(w.Ω_crit_over_ΩK)
        minratio = w.Ω_crit_over_ΩK[imin]
        νc  = w.ν_crit_Hz[imin]
        Tmin = w.T[imin]
        Mreal = mass_solar(star)
        println("$name,$(round(Mreal,digits=4)),$(round(R,digits=4)),$(round(C,digits=4)),",
                "$(round(ΩK,digits=2)),$(round(w.ν_K_Hz,digits=2)),",
                "$(round(minratio,digits=5)),$(round(νc,digits=2)),$(round(Tmin,sigdigits=4))")
    end
end

# LOM98 n=1 benchmark validation (sanity)
v = BDNKStar.RModes.rmode_validate_lom98()
println("# LOM98 n=1 check: Jt=$(round(v.J̃,sigdigits=5)) It=$(round(v.Ĩ,sigdigits=5)) ",
        "tauGW=$(round(v.τ_GW,digits=3)) tauSV=$(round(v.τ_sv,sigdigits=4)) tauBV=$(round(v.τ_bv,sigdigits=4))")
