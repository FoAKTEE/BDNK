# ======================================================================
# repro/dgstarfd_hybrid.jl — THE DG / FINITE-DIFFERENCE HYBRID ON THE RADIAL COWLING STAR
# (src/dg/DGStarFD.jl, machinery in src/dg/DGSubcell.jl), the scheme Deppe et al. (PRD 105,
# 123031; arXiv:2109.12033) endorse after every classical DG limiter fails on a neutron star.
#
# Run on the SAME grids as the pure-DG study (VALIDATION.md §7.9) so the comparison is
# limiter-for-limiter, plus SpECTRE-style uniform high-order grids:
#   I1  10 linear surface elements   — the configuration whose pure-DG minmod blew a wind
#   I2   5 quadratic surface elements — the configuration the pure-DG study favoured
#   uniform p=5, K per side          — SpECTRE's setup (surface inside an element, no grading)
#   always_fd                        — every element on the subgrid: pure second-order finite
#                                      volume, the reference for what the DG elements buy
#
# Output: repro/data/dgstarfd_hybrid_summary.csv (one row per run) and
# repro/data/dgstarfd_hybrid_series_<tag>.csv. Linear radial Cowling modes of this star:
# F = 2.686, H1 = 4.550, H2 = 6.341, H3 = 8.108 kHz. Runtime ≈ 25 min single-threaded.
# ======================================================================
using BDNKStar, Printf
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
datadir = joinpath(@__DIR__, "data")
R = 9.585578
uni(K, p; xmax=3.0) = [(xlo=0.0, xhi=xmax*R, K=K, p=p, region=:interior)]
Flin = radial_cowling_spectrum(eos, εc; nmodes=4)[1] ./ BDNKStar.Units.Msun_to_km
@printf("linear Cowling radial modes [kHz]: %s\n", join([@sprintf("%.4f", f) for f in Flin], "  "))

runs = [  # tag, grid kwargs, seed amplitude, T
    # SpECTRE's configuration: uniform high-order elements, the surface inside an element, no
    # grading and no limiter — robustness comes from the subcells alone. Resolution scan.
    ("uni_p5_K10_unseeded", (grid=uni(10,5), p_center=5, center_fac=1), 0.0, 10000.0),
    ("uni_p5_K14_unseeded", (grid=uni(14,5), p_center=5, center_fac=1), 0.0, 10000.0),
    ("uni_p5_K20_unseeded", (grid=uni(20,5), p_center=5, center_fac=1), 0.0, 10000.0),
    ("uni_p5_K28_unseeded", (grid=uni(28,5), p_center=5, center_fac=1), 0.0, 10000.0),
    ("uni_p5_K10_seeded", (grid=uni(10,5), p_center=5, center_fac=1), 1e-3, 4000.0),
    ("uni_p5_K14_seeded", (grid=uni(14,5), p_center=5, center_fac=1), 1e-3, 4000.0),
    ("uni_p5_K20_seeded", (grid=uni(20,5), p_center=5, center_fac=1), 1e-3, 4000.0),
    ("uni_p5_K28_seeded", (grid=uni(28,5), p_center=5, center_fac=1), 1e-3, 4000.0),
    # pure second-order finite volume on the same grid: what the DG elements buy
    ("uni_p5_K20_purefv", (grid=uni(20,5), p_center=5, center_fac=1, always_fd=true), 1e-3, 4000.0),
    # the hp grids of the pure-DG study (VALIDATION §7.9), for the limiter-for-limiter comparison
    ("I2_seeded", (grid=:I2,), 1e-3, 4000.0),
    ("I1_seeded", (grid=:I1,), 1e-3, 4000.0),
]

open(joinpath(datadir, "dgstarfd_hybrid_summary.csv"), "w") do io
    println(io, "tag,elements,dg_nodes,subcells,dt,T,seed,errD_end,rhoc_drift_end,Mb_drift_end,vatm_max,fd_fraction_end,ndrop,nback,peak1_kHz,peak2_kHz,peak3_kHz,peak4_kHz,seconds")
    for (tag, kw, A, T) in runs
        eng, st = setup_dgstarfd(eos, εc; kw...)
        A > 0 && seed_dgstarfd_radial!(st, eng; A=A)
        rec = evolve_dgstarfd!(st, eng; tmax=T, sample_dt=0.5)
        open(joinpath(datadir, "dgstarfd_hybrid_series_$(tag).csv"), "w") do ios
            println(ios, "t,rhoc,errD,Mb,fd_fraction,vatm")
            for i in eachindex(rec.t)
                @printf(ios, "%.3f,%.10e,%.6e,%.10e,%.4f,%.4e\n", rec.t[i], rec.ρc[i], rec.errD[i], rec.Mb[i], rec.fd_fraction[i], rec.vatm[i])
            end
        end
        pk, _, _ = dgstarfd_spectrum(rec.t, rec.ρc; npeaks=4, tmax=min(T, 4000.0))
        pk4 = vcat(pk, fill(NaN, 4 - length(pk)))
        am = dgstarfd_active_map(eng)
        @printf(io, "%s,%d,%d,%d,%.5f,%.0f,%.0e,%.3e,%.3e,%.3e,%.3e,%.3f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.0f\n",
                tag, eng.K, eng.Ntot, eng.Mtot, rec.dt, T, A, rec.errD[end], rec.ρc[end]/rec.ρc[1]-1,
                rec.Mb[end]/rec.Mb[1]-1, maximum(rec.vatm), rec.fd_fraction[end], rec.ndrop, rec.nback,
                pk4..., rec.seconds); flush(io)
        @printf("%-22s K=%3d Ndg=%3d Nfd=%3d [%4.0f s]: errD=%.2e ρc=%+.2e Mb=%+.1e vatm=%.1e FD=%.0f%% (r/R %s) peaks=%s\n",
                tag, eng.K, eng.Ntot, eng.Mtot, rec.seconds, rec.errD[end], rec.ρc[end]/rec.ρc[1]-1,
                rec.Mb[end]/rec.Mb[1]-1, maximum(rec.vatm), 100*rec.fd_fraction[end],
                am.n == 0 ? "-" : @sprintf("%.2f..%.2f", minimum(am.r_over_R), maximum(am.r_over_R)),
                join([@sprintf("%.3f", x) for x in pk], ",")); flush(stdout)
    end
end
println("DGSTARFD DONE")
