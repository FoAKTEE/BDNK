#!/usr/bin/env julia
# ======================================================================
# hpc/bench_dyngr.jl — measure what a 1+1D dynamical-GR run actually costs, so that the
# "should this go on a GPU?" question is answered by numbers rather than by intuition.
#
#     julia --project=.. -t 16 bench_dyngr.jl
#
# Reports, for a range of radial resolutions: wall time, cell-updates per second, the size of
# the whole evolved state, and the parallel work available inside ONE run. Then measures the
# ensemble throughput across threads, which is the parallelism that actually exists in this
# problem. Re-run it on the cluster before trusting the numbers in README.md.
# ======================================================================
using BDNKStar, Printf
using Base.Threads: @threads, nthreads

const eos = ShumPolytrope(100.0)
const εc  = 0.00128 + 100*0.00128^2
const TB  = 20.0                        # short, only to time the operator

function one_run(N, tmax)
    eng, st = setup_dyngr(eos, εc; N=N)
    t0 = time(); res = evolve_dyngr!(st, eng; tmax=tmax); (time()-t0, length(res[1]))
end

println("=== single run: cost and the parallel work available inside it ===")
@printf("%-7s %-10s %-9s %-14s %-14s %s\n", "N", "wall[s]", "samples", "cell-upd/s", "state bytes", "concurrent cells")
one_run(200, 1.0)                        # warm up the JIT
for N in (200, 400, 800, 1600, 3200)
    wall, ns = one_run(N, TB)
    @printf("%-7d %-10.3f %-9d %-14.2e %-14d %d\n", N, wall, ns, N*ns/wall, 8*8*N, N)
end

println("""
A GPU needs O(10^4–10^5) concurrent work items to beat a single CPU core, and each kernel
launch costs ~5–10 microseconds. The rightmost column is the entire parallelism inside one
run. Two further obstacles are structural, not tuning problems:
  * the constraint solve integrates the Hamiltonian constraint and the polar-slicing lapse
    ODE OUTWARD in r by trapezoid, twice per Runge-Kutta substage. That is a sequential scan
    along the only spatial dimension — the access pattern GPUs are worst at;
  * the primitive recovery is a per-cell root find, so threads in a warp diverge.
The parallelism that does exist is ACROSS runs, measured next.""")

println("\n=== ensemble: many independent stars, which is the real parallelism ===")
nmem = 4*nthreads()
εs = collect(range(0.0006, 0.0028; length=nmem))
t0 = time()
@threads for i in 1:nmem
    eng, st = setup_dyngr(eos, εs[i] + 100*εs[i]^2; N=800)
    evolve_dyngr!(st, eng; tmax=TB)
end
wall = time() - t0
@printf("%d members, N=800, t=%.0f: %.2f s wall on %d threads → %.2f s/member/core\n",
        nmem, TB, wall, nthreads(), wall*nthreads()/nmem)
@printf("extrapolated: a 10,000-member sweep at t=2000 M⊙ costs %.0f core-hours at this efficiency\n",
        1e4 * (wall*nthreads()/nmem) * (2000/TB) / 3600)
# single-thread reference, to expose threading efficiency
t1 = time()
for i in 1:max(nthreads(), 2)
    eng, st = setup_dyngr(eos, εs[i] + 100*εs[i]^2; N=800); evolve_dyngr!(st, eng; tmax=TB)
end
ser = (time() - t1)/max(nthreads(), 2)
@printf("serial reference: %.2f s/member → threading efficiency %.0f%% on %d threads\n",
        ser, 100*ser/(wall*nthreads()/nmem), nthreads())
println("""
Read this as: the ensemble is embarrassingly parallel, but Julia's stop-the-world garbage
collector serialises the threads because the right-hand side allocates work arrays on every
call. Preallocating them in the engine is worth roughly a factor two on a node and is a far
better use of effort than a GPU port. See README.md, "Where the speed actually is".""")
