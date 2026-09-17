# Masked (cut=false, rigid-wall) runs at the low dissipation sigma_KO=0.004 where the mask's answer is
# estimator-unambiguous (-6.9%), 23 periods as in the production cut-cell runs. Appends to
# repro/data/cowling3d_mask_sigko.csv. Run: julia --project=. repro/paper_figure_runs_mask.jl
using BDNKStar, Printf
using BDNKStar: Msun_to_km, kHz_to_km
const D = joinpath(@__DIR__, "data")
eos = ShumPolytrope(100.0); εc = 0.00128 + 100*0.00128^2
const F1D = 1.88291; ν_ref = F1D*kHz_to_km*Msun_to_km; P = 1/ν_ref; NPER = 23
out = open(joinpath(D, "cowling3d_mask_sigko.csv"), "w")
println(out, "# CowlingEvolve3D, cut=false (binary interior mask), sigma_KO=0.004, dt=0.30dx, 23 f-mode periods, seed (2,0), anchor star.")
println(out, "surface,sigko,nper,N,m,f_pgram,f_pencil,gamma,stable,ratio,agree_pct,df_kHz"); flush(out)
for N in (24, 32, 40)
    s = build_star3d(eos, εc; N=N, Lfac=1.2); e = setup_evo3d(s; σ_ko=0.004, cut=false)
    dt = 0.30*s.grid.dx; nsteps = ceil(Int, NPER*P/dt); sample = max(1, floor(Int, P/(50dt)))
    st = EvolState(N); seed_ylm!(st, e; l=2, m=0, A=1e-3)
    ts, Q = evolve3d_moments!(st, e; dt=dt, nsteps=nsteps, sample=sample, moments=[(2, 0)])
    a = analyze_qnm(ts, Q[:,1]; ν_ref=ν_ref)
    @printf(out, "mask,0.004,%d,%d,0,%.7f,%.7f,%.6e,%s,%.4g,%.4f,%.4f\n", NPER, N, a.f_kHz, a.f_pencil_kHz, a.γ, a.stable, a.envelope_ratio, a.agree_pct, a.df_kHz); flush(out)
    @printf("mask σ=0.004 N=%d  f=%.6f/%.6f  stable=%s ratio=%.3g\n", N, a.f_kHz, a.f_pencil_kHz, a.stable, a.envelope_ratio)
end
close(out); println("MASK RUNS DONE")
