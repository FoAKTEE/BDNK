# ======================================================================
# repro/dgstarhp_hkt.jl — THE 1D hp-ADAPTED DG STAR OF HÉBERT–KIDDER–TEUKOLSKY (2018), Sec. VI.A,
# rebuilt in BDNKStar (src/dg/DGStarHP.jl) and compared limiter-by-limiter.
#
# Produces repro/data/dgstarhp_hkt_summary.csv (one row per run) and
# repro/data/dgstarhp_hkt_series_<tag>.csv (time series for the plotted runs; the figure is
# repro/dgstarhp_hkt_figure.py → paper/figs/dgstarhp_hkt.png). Write-up: VALIDATION.md §7.9.
#
# Grids (their Table I rescaled to the areal radius R=9.586 M⊙ of the κ=100, Γ=2, ρ_c=1.28e-3
# star — the same star they evolve): I1 = 12½ interior p=3 + 10 linear surface + 7 exterior;
# I2 = 5 quadratic surface elements; I1R/I2R refined. The surface sits on an element boundary.
# Central element: p=7 over 3 nominal widths (the p=3 centre has an O(1) static force error).
# Limiters: :minmod = the paper's ΛΠ¹ on the surface elements; :wb = positivity scaling about
# the equilibrium-preserving reference of DGCart3D; :mean = mean-based scaling.
# wellbalanced=false is the paper's scheme (the star settles to its numerical equilibrium);
# wellbalanced=true subtracts the static residual (exact fixed point).
# Runtime ≈ 15 min single-threaded.
# ======================================================================
using BDNKStar, Printf
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
datadir = joinpath(@__DIR__, "data")
runs = [ # (tag, grid, limiter, wellbalanced, A, T, save_series); seeds are v̂ = A sin(πr/R)
    ("I1_minmod_noWB_static", :I1, :minmod, false, 0.0, 10000.0, true),     # the paper's I1
    ("I1_minmod_noWB_seeded", :I1, :minmod, false, 1e-3, 4000.0, true),
    ("I1R_minmod_noWB_static", :I1R, :minmod, false, 0.0, 10000.0, true),   # the paper's I1R
    ("I1R_minmod_noWB_seeded", :I1R, :minmod, false, 1e-3, 4000.0, true),
    ("I2_minmod_noWB_static", :I2, :minmod, false, 0.0, 4000.0, true),      # the paper's I2 (quadratic surface)
    ("I2_wb_noWB_static", :I2,  :wb,     false, 0.0,  10000.0, true),       # our limiter, paper's scheme
    ("I2_wb_noWB_seeded", :I2,  :wb,     false, 1e-3, 4000.0, true),
    ("I2_wb_WB_static",   :I2,  :wb,     true,  0.0,  4000.0, true),        # our limiter + subtraction
    ("I2_wb_WB_seeded",   :I2,  :wb,     true,  1e-3, 4000.0, true),
    ("I2R_wb_noWB_static", :I2R, :wb,    false, 0.0,  4000.0, true),
    ("I2R_wb_WB_seeded",  :I2R, :wb,     true,  1e-3, 4000.0, true),
    ("I1_wb_noWB_static", :I1,  :wb,     false, 0.0,  4000.0, false),       # our limiter on linear elements
    ("I1_wb_WB_seeded",   :I1,  :wb,     true,  1e-3, 4000.0, false),
    ("I2_mean_WB_seeded", :I2,  :mean,   true,  1e-3, 4000.0, false),       # mean-based scaling
    ("I2_none_WB_static", :I2,  :none,   true,  0.0,  4000.0, false),       # no limiter at all
    ("I1_minmod_noWB_seededlin", :I1, :minmod, false, 1e-3, 4000.0, false), # homologous seed hits the surface
]
Fl = radial_cowling_spectrum(eos, εc; nmodes=4)[1] ./ BDNKStar.Units.Msun_to_km   # kHz (εc given in M⊙ units)
@printf("linear Cowling radial modes [kHz]: %s\n", join([@sprintf("%.4f",f) for f in Fl], "  "))
open(joinpath(datadir, "dgstarhp_hkt_summary.csv"), "w") do io
    println(io, "tag,grid,limiter,wellbalanced,A,T,nodes,dt,errD_end,rhoc_drift_end,Mb_drift_end,Snorm_end,vatm_max,rstar_max_R,amp_first300,amp_last300,peak1_kHz,peak2_kHz,peak3_kHz,peak4_kHz,nfix,seconds")
    for (tag, grid, lim, wb, A, T, save) in runs
        eng, st = setup_dgstarhp(eos, εc; grid=grid, limiter=lim, wellbalanced=wb)
        A > 0 && seed_dgstarhp_radial!(st, eng; A=A, profile=(endswith(tag,"lin") ? :linear : :sine))
        t0 = time(); rec = evolve_dgstarhp!(st, eng; tmax=T, sample_dt=0.5); secs = time()-t0
        pk, _, _ = dgstarhp_spectrum(rec.t, rec.ρc; npeaks=4, tmax=4000.0)
        pk4 = vcat(pk, fill(NaN, 4-length(pk)))
        q = rec.ρc ./ rec.ρc[1] .- 1; t = rec.t
        m0 = t .< 300; m1 = t .> T-300
        a0 = maximum(abs.(q[m0] .- sum(q[m0])/count(m0))); a1 = maximum(abs.(q[m1] .- sum(q[m1])/count(m1)))
        @printf(io, "%s,%s,%s,%s,%.0e,%.0f,%d,%.5f,%.3e,%.3e,%.3e,%.3e,%.3e,%.4f,%.3e,%.3e,%.4f,%.4f,%.4f,%.4f,%d,%.0f\n", tag, grid, lim, wb, A, T, eng.Ntot, rec.dt,
                rec.errD[end], q[end], rec.Mb[end]/rec.Mb[1]-1, rec.Snorm[end], maximum(rec.vatm), maximum(rec.rstar)/eng.R, a0, a1, pk4..., rec.nfix, secs); flush(io)
        @printf("%-26s errD=%.2e ρc drift=%+.2e Mb drift=%+.1e vatm_max=%.2f rstar_max=%.2fR amp %.1e→%.1e peaks=%s [%.0fs]\n", tag,
                rec.errD[end], q[end], rec.Mb[end]/rec.Mb[1]-1, maximum(rec.vatm), maximum(rec.rstar)/eng.R, a0, a1, join([@sprintf("%.3f",x) for x in pk],","), secs)
        if save
            open(joinpath(datadir, "dgstarhp_hkt_series_$(tag).csv"), "w") do ios
                println(ios, "t,rhoc_rel,errD,Snorm,Mb_rel,vatm,rstar_R")
                for i in eachindex(rec.t)
                    @printf(ios, "%.2f,%.8e,%.6e,%.6e,%.10e,%.4e,%.4f\n", rec.t[i], q[i], rec.errD[i], rec.Snorm[i], rec.Mb[i]/rec.Mb[1]-1, rec.vatm[i], rec.rstar[i]/eng.R)
                end
            end
        end
    end
end
println("HKT DONE")
