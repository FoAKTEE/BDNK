# ======================================================================
# Phase-1 counterpart runs in the genuine 1+1D engine (code/BDNKStar/src/dyngr/DynGR1D.jl):
# the same two stars and the same cubic kick as the AthenaK inputs, so the collapse/migration
# verdict, the central proper time and the horizon mass can be compared across the two codes
# and the two gauges (areal/polar slicing here, puncture gauge there).
# Usage (from code/BDNKStar):  julia --project=. ../athenak-bdnk/analysis/p1_dyngr1d.jl
# Writes ../athenak-bdnk/runs/dyngr1d_p1_<tag>.csv and prints the summary table.
# ======================================================================
using BDNKStar, Printf
using BDNKStar: setup_dyngr, evolve_dyngr!, seed_dyngr_velocity!, solve_tov, ShumPolytrope
eos = ShumPolytrope(100.0)
out = normpath(joinpath(@__DIR__, "..", "runs"))
ρ_of_ε(ε) = (-1 + sqrt(1 + 4*100*ε))/200
ε_of_ρ(ρ) = ρ + 100*ρ^2
# The same signed kicks as the AthenaK inputs. An UNPERTURBED unstable star is decided by
# whatever the scheme's own error happens to be — with the §7.15 operator DynGR1D collapsed it at
# t ≈ 194, with the hydrostatic reconstruction of §7.16 it expands instead (ρ_c down to 8% of its
# initial value by t = 300) — so both tests carry an EXPLICIT kick and neither relies on that.
# The migration run needs a wider box because the star expands well past the stable-branch radius
# on its first overshoot (rmax = 4R, N doubled).
runs = [ ("p1a_collapse_eps0003", 2.4162e-3, -0.03, 400.0, 2.0,  800),
         ("p1b_collapse_font",    7.993e-3,  -0.01, 400.0, 2.0,  800),
         ("p1b_migration_font",   7.993e-3,  +0.05, 600.0, 4.0, 1600) ]
@printf("%-22s %8s %8s %8s | %18s %8s %10s %8s %8s %10s %8s
", "run", "rho_c", "M", "R",
        "verdict", "t_AH", "tau_c(AH)", "M_AH", "max2m/r", "rhoc_late", "alpha_c")
for (tag, ρc, A, T, rfac, Ncell) in runs
    εc = ε_of_ρ(ρc); star = solve_tov(eos, εc; h=2e-4)
    eng, st = setup_dyngr(eos, εc; N=Ncell, cfl=0.3, rmax_fac=rfac, atm_vbc=:outflow)
    A != 0.0 && seed_dyngr_velocity!(st, eng; A=A, profile=:cubic)
    res = evolve_dyngr!(st, eng; tmax=T, probe_frac=0.5, sample_dt=0.25)
    # central proper time tau_c = ∫ alpha_c dt, and the horizon: first 2m/r > 0.98
    τ = zeros(length(res.ts))
    for i in 2:length(res.ts); τ[i] = τ[i-1] + 0.5*(res.αc[i]+res.αc[i-1])*(res.ts[i]-res.ts[i-1]); end
    # In polar slicing the lapse collapses BEFORE 2m/r reaches 1 (that is what makes the slicing
    # singularity-avoiding), so the horizon is approached asymptotically: take the trapped-surface
    # threshold at 2m/r > 0.9 and read the horizon mass off the mass function there, M_AH = r/2.
    iAH = findfirst(>(0.9), res.max2mor)
    tAH = iAH === nothing ? NaN : res.ts[iAH]; τAH = iAH === nothing ? NaN : τ[iAH]
    MAH = NaN
    if iAH !== nothing
        g = eng.g
        for i in 1:g.N
            ai = g.NG + i
            if 2*st.m[ai]/g.r[ai] > 0.9; MAH = 0.5*g.r[ai]; break; end
        end
    end
    m2max = maximum(res.max2mor); ρmax_over = maximum(res.ρc)/res.ρc[1]
    Mb = 4π*sum(st.D[eng.g.NG+i] for i in 1:eng.g.N)*eng.g.Δr
    late = [res.ρc[i] for i in eachindex(res.ts) if res.ts[i] > 0.5*res.ts[end]]
    open(joinpath(out, "dyngr1d_p1_$(tag).csv"), "w") do io
        println(io, "t,tau_c,rhoc,alphac,max2mor")
        for i in eachindex(res.ts)
            @printf(io, "%.4f,%.4f,%.10e,%.8f,%.6f\n", res.ts[i], τ[i], res.ρc[i], res.αc[i], res.max2mor[i])
        end
    end
    verdict = res.collapsed ? "COLLAPSE" : (minimum(res.ρc) < 1e-6 ? "disperse" : "migrate/oscillate")
    @printf("%-22s %8.3e %8.5f %8.4f | %18s %8.2f %10.2f %8.4f %8.4f %10.3e %8.2e
", tag, ρc, star.M, star.R,
            verdict, tAH, τAH, MAH, m2max, isempty(late) ? NaN : sum(late)/length(late), res.αc[end])

end
println("DYNGR1D P1 DONE")
