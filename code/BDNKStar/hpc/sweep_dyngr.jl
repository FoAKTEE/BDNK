#!/usr/bin/env julia
# ======================================================================
# hpc/sweep_dyngr.jl — parameter sweep over the 1+1D dynamical-GR star (DynGR1D) sized for a
# cluster: a SLURM array task takes a stride of the member list and runs its members
# concurrently across the cores it was given.
#
# DynGR1D holds no module-level mutable state and BDNKStar is stdlib-only, so members are
# independent and nothing needs to be downloaded on the compute node.
#
# Usage (one node, interactive):
#     julia --project=.. -t 16 sweep_dyngr.jl
# Usage (SLURM array, see dyngr_array.sbatch):
#     SLURM_ARRAY_TASK_ID / SLURM_ARRAY_TASK_COUNT select this task's stride.
#
# Environment overrides (all optional):
#     SWEEP_OUT      output directory                          (default ./sweep_out)
#     SWEEP_N        radial cells per member                   (default 800)
#     SWEEP_TMAX     evolution time per member, M⊙             (default 2000)
#     SWEEP_NEPS     number of central densities               (default 24)
#     SWEEP_EPS_LO   lowest central energy density             (default 4.0e-4)
#     SWEEP_EPS_HI   highest central energy density            (default 3.2e-3)
#     SWEEP_AMPS     comma-separated seed amplitudes           (default 1e-3)
#     SWEEP_KAPPA    polytropic constant                       (default 100)
#
# Each task writes ONE csv, sweep_<taskid>.csv, so array tasks never collide. Concatenate with
#     cat sweep_out/sweep_*.csv | sort -u -t, -k1,1n > sweep_all.csv
# ======================================================================
using BDNKStar
using Printf
using Base.Threads: @threads, nthreads, threadid

# ---------------------------------------------------------------- configuration
getenvf(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d
getenvi(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d

const OUTDIR = get(ENV, "SWEEP_OUT", joinpath(@__DIR__, "sweep_out"))
const NCELL  = getenvi("SWEEP_N", 800)
const TMAX   = getenvf("SWEEP_TMAX", 2000.0)
const NEPS   = getenvi("SWEEP_NEPS", 24)
const EPSLO  = getenvf("SWEEP_EPS_LO", 4.0e-4)
const EPSHI  = getenvf("SWEEP_EPS_HI", 3.2e-3)
const KAPPA  = getenvf("SWEEP_KAPPA", 100.0)
const AMPS   = [parse(Float64, s) for s in split(get(ENV, "SWEEP_AMPS", "1e-3"), ",")]

const TASK_ID = getenvi("SLURM_ARRAY_TASK_ID", 0)
const TASK_N  = getenvi("SLURM_ARRAY_TASK_COUNT", 1)

# ---------------------------------------------------------------- member list
struct Member
    idx::Int
    εc::Float64
    amp::Float64
end
function build_members()
    ms = Member[]
    i = 0
    for εc0 in range(EPSLO, EPSHI; length=NEPS), amp in AMPS
        i += 1
        # the engine takes the central ENERGY density; for the p = κρ² barotrope the TOV
        # convention used throughout this repo is ε_c = ρ_c + κρ_c²
        push!(ms, Member(i, εc0 + KAPPA*εc0^2, amp))
    end
    ms
end

# members this task owns: a stride, so a short task list still covers the whole range
mine(ms) = TASK_N <= 1 ? ms : [m for m in ms if (m.idx - 1) % TASK_N == TASK_ID % TASK_N]

# ---------------------------------------------------------------- one member
function run_member(m::Member, eos)
    t0 = time()
    row = Dict{String,Any}("idx" => m.idx, "eps_c" => m.εc, "amp" => m.amp, "N" => NCELL, "tmax" => TMAX)
    try
        eng, st = setup_dyngr(eos, m.εc; N=NCELL)
        row["M"] = eng.M; row["R"] = eng.R
        m.amp > 0 && seed_dyngr_velocity!(st, eng; A=m.amp)
        ts, probe, ρc, αc, max2mor, drift, collapsed = evolve_dyngr!(st, eng; tmax=TMAX)
        row["drift"] = drift
        row["collapsed"] = collapsed ? 1 : 0
        row["alpha_c_min"] = minimum(αc)
        row["max_2mor"] = maximum(max2mor)
        row["rhoc_end_over_start"] = ρc[end]/ρc[1]
        # dyngr_radial_freq returns (f_geom, P, fgrid) with f_geom in km⁻¹; the conversion to
        # kHz is f[kHz] = f[km⁻¹]/kHz_to_km with NO M⊙ rescaling (see its docstring). The
        # spectral resolution of a record of length T is df = 1/T, recorded alongside.
        fpk, _, _ = dyngr_radial_freq(ts, ρc)
        row["F_kHz"] = fpk/BDNKStar.Units.kHz_to_km
        row["df_kHz"] = (1/max(ts[end]-ts[1], eps()))/BDNKStar.Units.kHz_to_km
        row["status"] = "ok"
    catch err
        row["status"] = "FAILED: " * first(sprint(showerror, err), 120)
        for k in ("M","R","drift","collapsed","alpha_c_min","max_2mor","rhoc_end_over_start","F_kHz","df_kHz")
            get!(row, k, NaN)
        end
    end
    row["wall_s"] = time() - t0
    row
end

const COLS = ["idx","eps_c","amp","N","tmax","M","R","F_kHz","df_kHz","drift","collapsed",
              "alpha_c_min","max_2mor","rhoc_end_over_start","wall_s","status"]
fmt(v) = v isa AbstractString ? v : (v isa Integer ? string(v) : @sprintf("%.10g", v))

function main()
    mkpath(OUTDIR)
    eos = ShumPolytrope(KAPPA)
    all = build_members(); ms = mine(all)
    @info "sweep: $(length(all)) members total, $(length(ms)) on this task " *
          "(task $TASK_ID of $TASK_N), $(nthreads()) threads, N=$NCELL, tmax=$TMAX"
    rows = Vector{Any}(undef, length(ms))
    @threads for i in eachindex(ms)
        rows[i] = run_member(ms[i], eos)
    end
    out = joinpath(OUTDIR, @sprintf("sweep_%04d.csv", TASK_ID))
    open(out, "w") do io
        println(io, join(COLS, ","))
        for r in rows
            println(io, join([fmt(get(r, c, "")) for c in COLS], ","))
        end
    end
    nfail = count(r -> !startswith(String(r["status"]), "ok"), rows)
    tot = sum(r -> r["wall_s"], rows; init=0.0)
    @info "wrote $out — $(length(rows)) members, $nfail failed, $(round(tot, digits=1)) core-seconds"
end

main()
