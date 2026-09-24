# Running the 1+1D dynamical-GR star (`DynGR1D`) on Caltech HPC

Everything here targets `src/dyngr/DynGR1D.jl`, the only engine in the toolkit that drops the
Cowling approximation and evolves the spacetime with the matter.

## The short answer about GPUs

**A single 1+1D run must not go on a GPU — it would be slower than one laptop core.** This is
not a tuning problem, it is the shape of the problem. Measured on this machine with
`hpc/bench_dyngr.jl`:

| radial cells | wall for t = 20 M⊙ | cell-updates/s | whole evolved state | concurrent work inside one run |
|---|---|---|---|---|
| 200 | 0.049 s | 2.2×10⁵ | 12.8 kB | 200 cells |
| 800 | 0.742 s | 2.3×10⁵ | 51 kB | 800 cells |
| 3200 | 10.9 s | 2.5×10⁵ | 205 kB | 3200 cells |

A GPU needs 10⁴–10⁵ concurrent work items to beat a single core, and a kernel launch costs
5–10 µs. One run offers a few hundred. Two further obstacles are structural:

* the constraint solve integrates the Hamiltonian constraint and the polar-slicing lapse ODE
  **outward in r by trapezoid, twice per Runge–Kutta substage**. That is a sequential scan along
  the only spatial dimension — precisely the access pattern GPUs handle worst;
* the primitive recovery is a per-cell root find, so threads in a warp diverge.

The parallelism in this problem is **across runs**, not inside one.

## Where the speed actually is

Measured, 16 members at N = 800:

| | per member | efficiency |
|---|---|---|
| 1 thread | 0.61–0.73 s | — |
| 8 threads | 1.75–1.83 s of core time | **35–40%** |

Garbage collection is 19% of wall time at 8 threads against 3% at one thread. The right-hand
side allocates its work arrays (`Ec`, `qcenter`, `mface`, …) on **every call**, and Julia's
stop-the-world collector then serialises every thread. Preallocating those in `DynGREngine` is
worth roughly a factor two on a node — a much better use of effort than any GPU port, and the
first thing to do if these sweeps become a bottleneck.

Cost of a 10,000-member sweep at t = 2000 M⊙, N = 800:

* 203 core-hours at single-thread efficiency;
* ≈ 32 h wall on one 16-core node at the current 40% threading efficiency;
* ≈ 4 h wall on an 8-task array of 16 cores (the `dyngr_array.sbatch` default);
* ≈ 1.8 h on the same array once the allocations are fixed.

## When a GPU *would* be the right tool

Only for an ensemble large enough that the numbers above become inconvenient — roughly
10⁵–10⁶ members, e.g. a Bayesian inference over equation-of-state and viscosity parameters. The
design would then be **one star per thread block**, not one cell per thread: a block of ~128
threads holds the radial grid in shared memory, the constraint solve becomes a block-level
parallel scan, and thousands of blocks run independent stars. That is a genuine rewrite
(`KernelAbstractions.jl` or `CUDA.jl`, plus a branch-free primitive recovery), not a port.
Ask before starting it; at 10⁴ members the CPU array wins on every axis including your time.

Note also that this repository's own plan (`progress/plan_athenak_bdnk_collapse.md`) already
designates a *different* code for the GPU allocation: AthenaK is Kokkos-based, and the 128³ Z4c
collapse runs are described there as "a GPU job (hours per run)". If the goal is to use a GPU
allocation productively, that is the workload that needs it — not this one.

## Deploying

`BDNKStar` is **stdlib-only** (`LinearAlgebra`, `Printf`), so there is nothing to download on a
compute node and no network access is required. You need a Julia ≥ 1.10 binary and the repo.

1. **Discover the cluster's names** — run once on a login node, submits nothing:

   ```bash
   bash hpc/discover.sh | tee hpc/cluster_facts.txt
   ```

   It prints your accounts, the partitions and their time limits, which partitions actually have
   GPUs, cores and memory per node, whether Julia is a module, and the scratch filesystems.

2. **Fill the placeholders** in `hpc/dyngr_array.sbatch`: `{{ACCOUNT}}`, `{{PARTITION}}`,
   `{{JULIA_LOAD}}`, `{{SCRATCH}}`.

3. **Sanity-check one member interactively** before burning an allocation:

   ```bash
   SWEEP_NEPS=1 SWEEP_EPS_LO=0.00128 SWEEP_EPS_HI=0.00128 SWEEP_N=800 SWEEP_TMAX=2000 \
   SWEEP_OUT=/tmp/chk julia --project=. hpc/sweep_dyngr.jl
   ```

   The anchor star must come back with M = 1.4001, R = 9.5856 and F = 2.17 kHz at a spectral
   resolution df = 0.15 kHz — consistent with the Chandrasekhar eigenvalue 2.1236 kHz that
   `VALIDATION.md` row 3 records. If it does not, stop.

4. **Submit**: `sbatch hpc/dyngr_array.sbatch`, then concatenate the per-task CSVs with the
   command in the footer of that file.

## Files

| file | what it is |
|---|---|
| `discover.sh` | read-only probe of the cluster; prints what to fill in |
| `sweep_dyngr.jl` | the sweep driver: SLURM-array aware, threads members across the cores it is given, one CSV per task |
| `dyngr_array.sbatch` | CPU array template |
| `bench_dyngr.jl` | reproduces every number above; re-run it on the cluster before trusting them |

## Two things to know about the physics before sweeping

* **The record must be long enough to resolve the mode.** The fundamental is ≈ 2.1 kHz, a
  period of ≈ 97 M⊙; the periodogram resolution is df = 1/T. At the default t = 2000 M⊙ that is
  0.15 kHz. A short run returns a meaningless peak — `sweep_dyngr.jl` records `df_kHz` in every
  row so this is visible rather than silent.
* **`DynGR1D` has a known caveat** recorded in `VALIDATION.md`: its atmosphere velocity boundary
  condition was implicated in a spurious damping exponent, and `setup_dyngr` carries an
  `atm_vbc` switch (`:zero` = pinned-velocity wall, `:outflow` = zero-gradient) for exactly that
  reason. Sweeps that measure *damping* rather than frequency should run both.
