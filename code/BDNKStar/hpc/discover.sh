#!/usr/bin/env bash
# ======================================================================
# hpc/discover.sh — run this ONCE on a Caltech HPC login node. It prints the facts needed to
# fill the {{PLACEHOLDERS}} in dyngr_array.sbatch and dyngr_gpu.sbatch. It submits nothing and
# changes nothing.
#
#     bash hpc/discover.sh | tee hpc/cluster_facts.txt
# ======================================================================
set -u
line() { printf '\n===== %s =====\n' "$1"; }
# print a command's output, or a fallback note when it produces nothing (an empty pipeline
# still exits 0, so "cmd || echo" would never fire)
try() { local out; out="$(eval "$1" 2>/dev/null)"; if [ -n "$out" ]; then printf '%s\n' "$out"; else printf '(nothing returned) %s\n' "$2"; fi; }

line "host / scheduler"
hostname -f 2>/dev/null || hostname
sinfo --version 2>/dev/null || echo "sinfo not found — is this a SLURM cluster?"

line "accounts you can charge to (use one as {{ACCOUNT}})"
try 'sacctmgr -nP show assoc user="$USER" format=Account,Partition | sort -u' \
    "— sacctmgr unavailable; try 'sshare -U' or ask the help desk"

line "partitions: name, time limit, nodes, state (pick {{PARTITION}} / {{GPU_PARTITION}})"
try "sinfo -o '%20P %12l %8D %10t %N' | head -30" "— no SLURM here"

line "GPU resources actually configured (empty output = no GPUs in these partitions)"
try "sinfo -o '%20P %10G %8D %N' | grep -v 'null' | head -20" "— no GPU GRES visible"
try "scontrol show config | grep -i gres | head" "— no GRES in the SLURM config"

line "cores and memory per node, by partition"
try "sinfo -o '%20P %6c %10m %N' | head -20" "— no SLURM here"

line "Julia: module or binary?"
try 'module avail julia 2>&1 | head -20' "— no module system, or no julia module"
try 'command -v julia && julia --version' "— no julia on PATH; load a module or install one (see README)"

line "scratch / project space (put JULIA_DEPOT_PATH and the sweep output on the fastest one)"
for d in "$HOME" /central/scratch/"$USER" /central/groups /scratch/"$USER" /resnick/scratch/"$USER"; do
  [ -d "$d" ] && printf '%-32s %s\n' "$d" "$(df -h "$d" 2>/dev/null | awk 'NR==2{print $4" free of "$2}')"
done

line "what to fill in"
cat <<'EOF'
  {{ACCOUNT}}        an Account from the first block
  {{PARTITION}}      a CPU partition with enough time limit
  {{GPU_PARTITION}}  a partition whose GRES column shows gpu:<type>:<n>   (only if you go the GPU route)
  {{JULIA_LOAD}}     e.g. "module load julia/1.10" or "" if julia is already on PATH
  {{SCRATCH}}        a large, fast directory you can write to
EOF
