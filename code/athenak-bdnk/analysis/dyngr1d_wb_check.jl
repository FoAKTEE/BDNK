# ======================================================================
# The well-balancing benchmark for DynGR1D (VALIDATION.md §7.15, §7.16).
#
# (1) RESIDUAL. The raw static momentum residual |∂_t S̃| over the stellar interior, divided by
#     the local gravity scale √γ̃(ρhW²−p)α ∂_r lnα, for each scheme. It is reported BOTH over the
#     whole star and over r < 0.95R, because the schemes fail in different PLACES and one number
#     over the whole star cannot tell them apart:
#       • plain operator — the residual is spread through the bulk and falls as Δr²;
#       • hydrostatic reconstruction — the bulk residual is ZERO (round-off, and therefore
#         Δr-INDEPENDENT), and what is left sits in the one cell where the star meets the
#         constant-density artificial atmosphere. That cell is not an equilibrium at all, so no
#         flux/source pairing can balance it; reporting it mixed in with the bulk is what hid
#         the difference between an exactly balanced operator and a second-order one.
#     The `before_fix` column is the pre-2026-09-22 operator, whose residual did not fall with
#     Δr at all (5.346e-2 and 1.578e-1 at every resolution — a consistency error, not truncation).
#     `plain+dic` isolates the initial data: the plain operator on the SAME discrete-equilibrium
#     initial data the hydrostatic scheme uses, which shows the initial data alone fixes nothing.
#
# (2) STATIC EVOLUTION. The unsubtracted star over 200 M⊙, with the central density as the
#     probe. ρ_c does not move at all (round-off) until a signal has had time to travel in from
#     the surface; the arrival time is reported, and it is what dates the remaining drift to the
#     star/atmosphere boundary rather than to the operator.
#
# Usage (from code/BDNKStar): julia --project=. ../athenak-bdnk/analysis/dyngr1d_wb_check.jl
# Writes ../athenak-bdnk/runs/dyngr1d_wb_residual.csv and dyngr1d_wb_static.csv.
# ======================================================================
using BDNKStar, Printf
using BDNKStar: setup_dyngr, evolve_dyngr!, ShumPolytrope
using BDNKStar.DynGR1D: _raw_rhs!

eos = ShumPolytrope(100.0)
out = normpath(joinpath(@__DIR__, "..", "runs"))

# the three schemes, as (tag, setup kwargs). Nothing is ever subtracted here: the point is
# what the OPERATOR does on its own.
const SCHEMES = (("plain",       (wellbalanced = :none,        discrete_ic = false)),
                 ("plain+dic",   (wellbalanced = :none,        discrete_ic = true)),
                 ("hydrostatic", (wellbalanced = :hydrostatic, discrete_ic = true)))

"""Momentum residual over r < `frac`·R, normalised by the local gravity scale."""
function residual(eng, st, frac)
    n = length(st.D); r = ntuple(_ -> zeros(n), 4)
    _raw_rhs!(r..., st, eng)
    g = eng.g; num = 0.0; den = 0.0
    for i in 1:g.N
        ai = g.NG + i; g.r[ai] < frac*eng.R || continue
        num += abs(r[2][ai])
        den += max(abs(st.sqrtg[ai]*st.ε[ai]*st.α[ai]*st.dlnα[ai]), 1e-300)
    end
    return num/den
end

open(joinpath(out, "dyngr1d_wb_residual.csv"), "w") do io
    println(io, "star,rho_c,scheme,N,dr,residual,residual_bulk,before_fix")
    for (tag, ρc, before) in (("stable", 1.28e-3, 5.346e-2), ("unstable", 7.993e-3, 1.578e-1))
        for (scheme, kw) in SCHEMES, N in (400, 800, 1600, 3200)
            eng, st = setup_dyngr(eos, ρc + 100*ρc^2; N=N, cfl=0.3, rmax_fac=2.0,
                                  atm_vbc=:outflow, kw...)
            full = residual(eng, st, 1.0); bulk = residual(eng, st, 0.95)
            @printf(io, "%s,%.5e,%s,%d,%.6f,%.6e,%.6e,%.6e\n",
                    tag, ρc, scheme, N, eng.g.Δr, full, bulk, before)
            @printf("%-8s %-12s N=%4d Δr=%.4f: whole star %.3e   r<0.95R %.3e\n",
                    tag, scheme, N, eng.g.Δr, full, bulk)
        end
        @printf("%-8s before the 2026-09-22 flux/source fix: %.3e at EVERY Δr\n", tag, before)
    end
end

open(joinpath(out, "dyngr1d_wb_static.csv"), "w") do io
    println(io, "scheme,N,t,rhoc_over_rhoc0")
    for (scheme, kw) in SCHEMES, N in (400, 800, 1600)
        eng, st = setup_dyngr(eos, 1.28e-3 + 100*1.28e-3^2; N=N, cfl=0.3, rmax_fac=2.0,
                              atm_vbc=:outflow, kw...)
        res = evolve_dyngr!(st, eng; tmax=200.0, sample_dt=1.0)
        d = abs.(res.ρc ./ res.ρc[1] .- 1)
        for i in eachindex(res.ts)
            @printf(io, "%s,%d,%.4f,%.12f\n", scheme, N, res.ts[i], res.ρc[i]/res.ρc[1])
        end
        # when does the centre first learn that the surface is not balanced?
        k = findfirst(>(1e-10), d)
        tsig = k === nothing ? Inf : res.ts[k]
        @printf("static, NOTHING subtracted, %-12s N=%4d: max |ρc/ρc0−1| = %.3e; centre at round-off until t = %.1f (before the fix: 8.7e-1)\n",
                scheme, N, maximum(d), tsig)
    end
end
println("DYNGR1D WB CHECK DONE")
