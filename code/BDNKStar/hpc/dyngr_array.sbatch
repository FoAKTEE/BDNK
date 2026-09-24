#!/usr/bin/env bash
# ======================================================================
# hpc/dyngr_array.sbatch — run the 1+1D dynamical-GR parameter sweep as a SLURM array.
#
# Fill the {{PLACEHOLDERS}} from `bash hpc/discover.sh`, then:
#     sbatch hpc/dyngr_array.sbatch
#
# Sizing: one member (N = 800 cells, t = 2000 M⊙) costs ≈ 77 core-seconds. With
# --cpus-per-task=16 a task finishes 16 members in ~80 s, so the default 24-member sweep is
# done in a couple of minutes on a single task, and a 10,000-member sweep takes about
# 13 core-hours — i.e. ~10 minutes wall on 8 tasks of 16 cores. There is no GPU in this path
# and none is wanted: see README.md ("Why not a GPU").
# ======================================================================
#SBATCH --job-name=dyngr_sweep
#SBATCH --account={{ACCOUNT}}
#SBATCH --partition={{PARTITION}}
#SBATCH --array=0-7                     # 8 tasks; each takes every 8th member
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G                        # the whole state is ~100 kB per member; this is ample
#SBATCH --time=02:00:00
#SBATCH --output=logs/dyngr_%A_%a.out
#SBATCH --error=logs/dyngr_%A_%a.err

set -euo pipefail
mkdir -p logs

{{JULIA_LOAD}}                          # e.g. module load julia/1.10   (delete if julia is on PATH)

# Keep the Julia depot OFF the home filesystem: precompilation is metadata-heavy and home
# directories on shared clusters are usually the slowest mount. BDNKStar is stdlib-only, so
# nothing is downloaded — the depot only holds precompiled caches.
export JULIA_DEPOT_PATH={{SCRATCH}}/julia_depot
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# one BLAS thread per Julia thread would oversubscribe; the sweep is not BLAS-bound
export OPENBLAS_NUM_THREADS=1

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SWEEP_OUT={{SCRATCH}}/dyngr_sweep
export SWEEP_N=800
export SWEEP_TMAX=2000
export SWEEP_NEPS=24
export SWEEP_AMPS=1e-3

echo "task ${SLURM_ARRAY_TASK_ID:-0} of ${SLURM_ARRAY_TASK_COUNT:-1} on $(hostname), ${JULIA_NUM_THREADS} threads"
srun julia --project="$REPO" "$REPO/hpc/sweep_dyngr.jl"

# Combine afterwards (from the login node):
#     cd {{SCRATCH}}/dyngr_sweep
#     head -1 "$(ls sweep_*.csv | head -1)" > sweep_all.csv
#     tail -q -n +2 sweep_*.csv | sort -t, -k1,1n >> sweep_all.csv
