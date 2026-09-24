# athenak-bdnk — AthenaK runs for the GR-coupled BDNK collapse programme

Plan: `progress/plan_athenak_bdnk_collapse.md`. AthenaK itself lives in `ref-code/athenak`
(gitignored, commit c5a0d7f "Cyclic zoom (#713)", Kokkos 4.7.2); this directory holds what is
ours — the Windows build recipe, inputs, run driver, analysis scripts and the record.

## Build on this laptop (Windows 11, no WSL) — 2026-09-21

Toolchain already present: CMake 4.4.3, MSYS2 UCRT64 GCC 16.2, Ninja 1.13. Kokkos OpenMP
backend (8 cores, no NVIDIA GPU). **No AthenaK source is modified**: three POSIX-isms are
shimmed by a force-included header, `win/win_compat.hpp` (`M_PI` under strict `-std=c++17`,
two-argument `mkdir`, `asctime_r`/`sleep` in the watchdog; `<windows.h>` must NOT be
included — it defines `r2`, `min`, `max`, … and breaks `gr_monopole.cpp`).

```bash
cd ref-code/athenak && mkdir -p build-win && cd build-win
export PATH="/c/msys64/ucrt64/bin:$PATH"
cmake -G Ninja -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_C_COMPILER=/c/msys64/ucrt64/bin/gcc.exe -D CMAKE_CXX_COMPILER=/c/msys64/ucrt64/bin/g++.exe \
  -D CMAKE_MAKE_PROGRAM="$(cygpath -m "$(which ninja)")" \
  -D PROBLEM=dyn_grmhd/dyngr_tov -D Athena_ENABLE_OPENMP=ON -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ENABLE_SERIAL=ON -D Kokkos_ARCH_NATIVE=ON \
  -D CMAKE_CXX_FLAGS="-include <abs path>/code/athenak-bdnk/win/win_compat.hpp" ..
ninja -j8 athena          # 3.5 min; src/athena.exe (18 MB)
```

Running: `OMP_NUM_THREADS=8 athena.exe -i <input>` from the run directory. Do not set
`OMP_PROC_BIND`/`OMP_PLACES` (libgomp affinity is unsupported on Windows). Input files must
have LF line endings, since the parser does not strip `\r`; this checkout has
`core.autocrlf = true`, so `.gitattributes` here pins `*.athinput` and `*.sh` to LF. Parameter
lines must also carry at most one `=`, inline comments included: AthenaK copies each line into
the header of every binary dump, and `vis/python/bin_convert.py` splits it on every `=` (the
shipped whisky input trips this itself with `kappa = 100.0 # P = kappa*rho^gamma`).
The shipped `<coord>` keys `general_rel`/`m` are reported unused by `dyn_grmhd` — harmless.

## Phase 0 runs (`inputs/`, driver `runs/run_p0.sh`, analysis `analysis/p0_history.py`)

The star of the whole project: K=100, Γ=2, ρ_c = 1.28×10⁻³ (M = 1.40016, R = 9.586, R_iso = 8.125),
on the octant [0, 20]³ of Deppe et al. 2021's box, isotropic initial slice, reflect/diode,
PPMX + HLLE + FOFC, RK3, the shipped floors (dfloor 10⁻¹⁰, tfloor 10⁻⁸, dthreshold 1.02).

| tag | spacetime | grid | kick | tlim | purpose |
|---|---|---|---|---|---|
| `p0a_tov_cowling_K32_timing` | `<adm>` (Cowling) | 32³ | — | 50 | cost: 3.8×10⁵ zone-cycles/s, 17 s |
| `p0b_tov_z4c_K32_timing` | Z4c, puncture gauge | 32³ | — | 50 | cost: 1.15×10⁵ zone-cycles/s, 91 s |
| `p0a_tov_cowling_K64_static` | Cowling | 64³ (Δx = 0.3125, 26 cells/R) | — | 2000 | settling, spectrum vs F/H₁/H₂ and vs the DG-FD hybrid of VALIDATION §7.13 |
| `p0a_tov_cowling_K64_kick` | Cowling | 64³ | v_pert = −10⁻³ | 2000 | clean radial spectrum |
| `p0b_tov_z4c_K64_static` | Z4c | 64³ | — | 1000 | gate: static star holds, full-GR F below Cowling |

Projected cost (from the 32³ timings, 8×cells × 2×steps): Cowling 64³ ≈ 3 h per 2000 M⊙,
Z4c 64³ ≈ 8 h per 1000 M⊙. Results are appended below when the runs finish.

### Results

**`p0a_tov_cowling_K64_static`** (2026-09-21, 3.34 h wall, 3.49×10⁵ zone-cycles/s, 16 000 cycles):
the projected star settles and stays: max ρ within [−5.9×10⁻⁴, +1.5×10⁻³] of its initial value over
2000 M⊙ (1.5×10⁻⁴ at the end), settling ringing of 3.4×10⁻⁴ peak-to-peak after t = 200, baryon
mass drift −3.9×10⁻⁶, min α constant (Cowling). Spectrum of max ρ(t) against the linear radial
Cowling modes of this star: **F 2.676 kHz (−0.4%)**, **H₁ 4.503 (−1.0%)**, H₂ 6.51 (+2.7%), a weak
8.33 (H₃ +2.7%), plus a 4.70/4.83 kHz doublet next to H₁ that is not a radial mode (the
Cartesian grid's non-radial settling; the octant box has no ℓ=2 seed but the staircase does).
Compare the DG–FD hybrid on the same star (VALIDATION §7.13, K=3: F −1.2%; K=6 pending) and
Deppe et al. Fig. 9 (their 12-element run resolves F, H₁ and a third mode). Gate for the Cowling
baseline: passed.

**`p0a_tov_cowling_K64_kick`** (2.75 h, 4.23×10⁵ zone-cycles/s): v_pert = −10⁻³ excites the
fundamental cleanly — F 2.683 kHz in [0, 700] and 2.677 in [700, 2000] (−0.1 / −0.3% of 2.6861),
H₁ 4.529 (−0.5%); ρ_max stays within ±4×10⁻³, M_b drift −3.9×10⁻⁶.

**`p0b_tov_z4c_K64_static`** (5.6 h, 1.68×10⁵ zone-cycles/s, 1000 M⊙ = 4.9 ms): **the Phase-0
gate**. The dynamical run is stable (ρ_max within ±0.8%, M_b drift −3.6×10⁻⁶); min α rises
0.6703 → 0.7065 and settles, which is the 1+log slicing relaxing off the TOV lapse, not a
collapse; the settling transient rings at **1.49 kHz** against our own full-GR fundamental for
this star, **F_GR = 1.4425 kHz** (`chandrasekhar_radial_omega2`, the corrected Chandrasekhar
LAWE of `test_dyngr`; the published Font et al. 2002 value for this benchmark star is ≈1.44 kHz):
**+3.3%**, and 45% below the Cowling F the same code reproduces to −0.4%. The metric responds —
this is what no Cowling engine in the toolkit can do, and it is the ideal-fluid baseline the
BDNK module must reproduce at η̂ = ζ̂ = 0.

| run | F measured | reference | Δ |
|---|---|---|---|
| Cowling, static (settling transient) | 2.676 kHz | 2.6861 (Cowling LAWE) | −0.4% |
| Cowling, kicked | 2.683 / 2.677 kHz | 2.6861 | −0.1 / −0.3% |
| Z4c, static (settling transient) | 1.49 kHz | 1.4425 (full-GR LAWE) | +3.3% |

Cost on this laptop (8 OpenMP threads, no usable GPU — see below): 64³ Cowling 3.3 h / 2000 M⊙,
64³ Z4c 5.6 h / 1000 M⊙. Phase 0 is complete.

## Can this laptop's GPU run AthenaK? Measured, 2026-09-22

The machine has an **Intel Arc 140V** — the Lunar Lake *integrated* GPU, which shares the one
LPDDR5X memory controller with the CPU (the "16 GB" is system RAM, not VRAM). Two questions were
settled with measurements rather than assumptions (`win/fp64_probe.c`, OpenCL, built with the same
MinGW GCC as AthenaK: `gcc -O2 -I<oneAPI>/compiler/2025.1/include fp64_probe.c C:/Windows/System32/OpenCL.dll`):

| device | FP64 native? | STREAM triad (fp64) | FMA chain (fp64) |
|---|---|---|---|
| Arc 140V iGPU (64 EU @ 2 GHz) | yes — `double_fp_config = 0x3f`, preferred vector width 1 | **96.1 GB/s** | **201.8 GFLOP/s** |
| Core Ultra 7 268V (8 cores) | yes | **74.9 GB/s** | **112.3 GFLOP/s** |

So FP64 is *native*, not emulated (an earlier guess to the contrary was wrong), but the iGPU's
ceiling over the CPU is only **1.28× on memory bandwidth and 1.8× on FP64 arithmetic** — both
devices are drinking from the same ~136 GB/s memory bus at 55–70% of its peak. A finite-volume
GRHD + Z4c code sits between those two limits, so **no GPU port of AthenaK can beat ~1.3–1.8× on
this machine**, before paying for the port.

And the port is not cheap. Kokkos' only Intel-GPU backend is SYCL, which Kokkos supports on Linux
(no Windows/SYCL path in the 4.7.2 CMake); oneAPI's `icpx` on Windows targets the MSVC ABI and has
no C++ standard library of its own (`fatal error: 'climits' file not found` — it needs Visual
Studio, not installed here, ~5 GB); and `src/utils/watchdog.cpp` uses POSIX `pthread`, which the
MinGW build satisfies with winpthreads and an MSVC build would not.

**Verdict: run on the CPU here, and take the GPU on the cluster.** AthenaK's well-trodden GPU path
is CUDA (`-D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_<arch>=ON`, same inputs, no source changes),
which is where the Phase-3 matrix belongs anyway. The oneAPI Base Toolkit installed for this test
(≈4 GB under `C:\Program Files (x86)\Intel\oneAPI`) is not needed by anything else here and can be
removed with `winget uninstall Intel.OneAPI.BaseToolkit`.

## Phase 1 runs (`inputs/p1*`, drivers `runs/run_p1.sh` and `runs/run_tags.sh`, 1+1D counterpart `analysis/p1_dyngr1d.jl`)

Ideal-fluid collapse and migration, Z4c, isotropic initial slice, octant [0,16]³ at 64³ (Δx = 0.25)
unless noted. `runs/run_tags.sh <tag>...` runs any inputs in sequence and honours
`OMP_NUM_THREADS`. The kick is the pgen's v_r = ½·v_pert·(3x − x³), x = r/R — but **the same
v_pert is a stronger kick than `DynGR1D`'s `profile=:cubic` gives for the same number**: the pgen
writes the Valencia primitive on the isotropic grid, so the orthonormal speed is larger by
ψ² = r_schw/r_iso (1.69 at the centre, 1.26 at the surface, least-squares 1.41×). Compare the codes
on the orthonormal velocity, not the input parameter.

| tag | setup | result |
|---|---|---|
| `p1a_collapse_eps0003_K64`, `_L1`, `_L2` | stable-branch star, v_pert = −0.03 | collapses (ρ ×25.6) where `DynGR1D` bounces; refinement changes nothing. **Dropped as a gate**: the star sits on a bifurcation that the 1.41× kick straddles |
| `p1b_collapse_font_K64`, `_L1` | Font et al. 2002 unstable star, −0.01 | collapses, as in `DynGR1D`. τ_c at lapse collapse 10.43 (Δx 0.25) and 8.39 (0.125) against `DynGR1D`'s 10.01 at the horizon: **not converged**; both lose the 3-metric at the centre at t ≈ 65 |
| `p1b_migration_font_K64` | same star, +0.05 | migrates, as in `DynGR1D`; clean to t = 800, but loses 13.3% of M_b through the diode boundary |
| `p1b_migration_font_box32` | [0,32], one level back on [0,16] | outflow cured (M_b −5.5×10⁻⁶ at t = 200 against −8.0×10⁻³), then a Z4c blow-up at its refinement interface at t = 221.5 |
| `p1b_migration_font_box32_uni64` | [0,32] uniform 64³, no interface | clean to t = 300, so the interface killed box32; the outflow returns once the second expansion reaches r = 32 |
| `p1b_collapse_font_ahf_cluster` | full box, four centre levels, `<fastflow>` | **deferred to the cluster** — smoke-tested; 1.2×10¹⁰ zone-cycles, about two days here |

The horizon run needs a full box because `<fastflow>` interpolates over a whole sphere around the
centre, which an octant cannot supply beyond its ghost cells, and four centre levels because the
horizon (isotropic radius M/2 = 0.64) is only 2.5 cells across at Δx = 0.25 — which is also why the
collapse runs break down as the lapse collapses.

**Status:** the verdict half of the Phase-1 gate passes on both Font setups; the quantitative half
does not. Every number, and how each was established: `progress/plan_athenak_bdnk_collapse.md`,
Phase 1.

The 1+1D counterpart, after the well-balancing fixes of VALIDATION.md §7.15–7.16
(`analysis/p1_dyngr1d.jl`, `runs/dyngr1d_p1_*.csv`):

| setup | verdict | t_AH | τ_c(AH) | M_AH | ρ_c late |
|---|---|---|---|---|---|
| stable star, kick −0.03 | oscillates (no collapse) | — | — | — | 2.56×10⁻³ |
| unstable star, kick −0.01 | **collapse** | 53.71 | 10.01 | 1.2734 | 1.38×10⁻² |
| unstable star, kick +0.05 | **migration** | — | — | — | 1.47×10⁻³ (stable-branch value) |

Compare in central proper time τ = ∫α dt, not coordinate time: the two codes use different gauges.
