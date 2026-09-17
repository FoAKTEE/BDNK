# Live runs for paper/figs: f-mode convergence (cut vs mask), m-seeds on the full cube, and the
# time-series figure. Writes repro/data/{cowling3d_anchor_convergence,cowling3d_timeseries_cut_N24,
# bdnk3d_timeseries_F1_N24}.csv. Cost ~1 CPU-hour (N=48 dominates). Run: julia --project=. repro/paper_figure_runs.jl
using BDNKStar, Printf
using BDNKStar: Msun_to_km, kHz_to_km
const D = joinpath(@__DIR__, "data")
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
const F1D = 1.88291; ν_ref = F1D*kHz_to_km*Msun_to_km; P = 1/ν_ref; NPER = 9
conv = open(joinpath(D, "cowling3d_anchor_convergence.csv"), "w")
println(conv, "# CowlingEvolve3D, sigma_KO=0.01, kappa_min=0.5, dt=0.30dx, $(NPER) f-mode periods, seed_ylm! (l=2,m), anchor star. surface=cut|mask.")
println(conv, "surface,N,m,f_pgram,f_pencil,gamma,stable,ratio,agree_pct,df_kHz"); flush(conv)
function linrun(N, cut, m; save_ts=nothing)
    s = build_star3d(eos, εc; N=N, Lfac=1.2); e = setup_evo3d(s; σ_ko=0.01, cut=cut)
    dt = 0.30*s.grid.dx; nsteps = ceil(Int, NPER*P/dt); sample = max(1, floor(Int, P/(50dt)))
    st = EvolState(N); seed_ylm!(st, e; l=2, m=m, A=1e-3)
    ts, Q = evolve3d_moments!(st, e; dt=dt, nsteps=nsteps, sample=sample, moments=[(2, m)])
    a = analyze_qnm(ts, Q[:,1]; ν_ref=ν_ref)
    @printf(conv, "%s,%d,%d,%.7f,%.7f,%.6e,%s,%.4g,%.4f,%.4f\n", cut ? "cut" : "mask", N, m, a.f_kHz, a.f_pencil_kHz, a.γ, a.stable, a.envelope_ratio, a.agree_pct, a.df_kHz); flush(conv)
    @printf("%-4s N=%d m=%d  f=%.6f/%.6f  stable=%s\n", cut ? "cut" : "mask", N, m, a.f_kHz, a.f_pencil_kHz, a.stable)
    if save_ts !== nothing
        open(joinpath(D, save_ts), "w") do io
            println(io, "# q20(t) of the ideal cut-cell run, N=24, sigma_KO=0.01, seed (2,0). t in f-mode periods (P from the 1D f=1.88291 kHz).")
            println(io, "t_over_P,q"); for (t, q) in zip(ts, Q[:,1]); @printf(io, "%.6f,%.8e\n", t/P, q); end
        end
    end
end
for N in (24, 32, 40, 48); linrun(N, true, 0; save_ts = N == 24 ? "cowling3d_timeseries_cut_N24.csv" : nothing); end
for N in (24, 32, 40);     linrun(N, false, 0); end
for N in (24, 32, 40), m in (2, 1); linrun(N, true, m); end
close(conv)
# BDNK F1 time series at N=24 (viscous run of the driver), sigma_KO=0.005
r = bdnk3d_qnm(eos, εc; N=24, η̂=0.03, frame=(25/4, 25/7), η̂_frame=0.03, control=false, nperiods=NPER, f_ref_kHz=F1D, σ_ko=0.005, dt_fac=0.20)
open(joinpath(D, "bdnk3d_timeseries_F1_N24.csv"), "w") do io
    println(io, "# q20(t) of the BDNK run, frame F1, eta_hat=0.03, N=24, sigma_KO=0.005, dt=0.20dx. gamma_per_period = envelope damping rate x P.")
    println(io, "t_over_P,q,gamma_per_period"); for (t, q) in zip(r.ts, r.q); @printf(io, "%.6f,%.8e,%.6f\n", t/P, q, r.γ*P); end
end
println("RUNS DONE")
