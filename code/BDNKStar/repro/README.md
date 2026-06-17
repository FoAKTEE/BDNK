# BDNKStar reproduction suite

Verified reference-paper reproductions, each a self-contained Julia script that
`include`s the `BDNKStar` package (`../src/BDNKStar.jl`), implements the paper's
method on top of the package API, runs, and prints a numeric comparison against
the published target. **No fabrication** — every number is computed.

**Run convention** (from the package root, so the relative `include`s resolve):

```bash
cd code/BDNKStar
JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 julia --project=. repro/<module>.jl
```

`<20` threads, no GPU. Dependency `include`s between modules are noted below.

## Stage 1B — axial QNM (Bussières 2604.13208 / Redondo-Yuste 2411.16841)
| module | reproduces | achieved | deps |
|---|---|---|---|
| `axial_waveeqs.jl` | coupled axial wave eqs (14a/14b) on TOV bg | inviscid limit → RW w-mode | — |
| `axial_qnm.jl` | ℓ=2 w-mode (f,τ) via shooting + Leaver | **<0.04%** vs Table II (10.4879,29.5786) | `axial_waveeqs.jl` |
| `axial_ultracompact.jl` | ultracompact trapped w-modes + η-mode | 6-mode ladder at C=0.44 | `axial_qnm.jl` |

## Stage 1C — nonlinear Cowling BDNK (Shum 2509.15303)
| module | reproduces | achieved | deps |
|---|---|---|---|
| `shum_core.jl` | general-EOS BDNK recovery + areal→isotropic | round-trip closes | — |
| `shum_evolve_opt.jl` | SSP-RK3 stellar evolution (optimized) | Δr=0.04 t_f=8000 stable, 3064 steps/s | `shum_core.jl` |
| `shum_qnm_production.jl` | FFT QNM + decay extraction | **F=2.699/H1=4.551/H2=6.468 kHz (<2%)** | `shum_evolve_opt.jl` |

## PMP viscous 1D (Pandya–Most–Pretorius 2209.09265 / 2201.12317)
| module | reproduces | achieved | deps |
|---|---|---|---|
| `pmp_viscous_core.jl` | general-EOS ideal-gas viscous BDNK 1D engine | const exact; Bjorken limit | — |
| `pmp_telegrapher.jl` | telegrapher dispersion + d'Alembert split | match 8.9e-16 | `pmp_viscous_core.jl` |
| `pmp_shock_instab.jl` | shock/acausal causality-crash classification | c₊ to 2 sig figs | `pmp_viscous_core.jl` |
| `pmp_heat_stationary.jl` | stationary heat-conduction profile | σ̂=0 → zero flux | `pmp_viscous_core.jl` |

## Analytic fleet (dispersion / causality / IS)
| module | reproduces | achieved | deps |
|---|---|---|---|
| `kovtun_sound.jl` | Kovtun sound channel (picresound/picstab) | Im ω≤0 stable, ±c_s pair | — |
| `is_contrast.jl` | Chabanov–Rezzolla IS causality-fixed c_s′² | =0.9 exact (1e-16) | — |
| `pmp2209_causal_frame.jl` | PMP causal/acausal frame classification | 4/4 (τ̂=0.25–1.5) | — |
| `conf_cc_cons.jl` | Pandya conformal char speeds + conservation | luminal=1, drift 3e-16 | — |

## Data
- `r5_eps_Dr0.04.txt` — canonical R5 central-ε(t) series (Δr=0.04, t_f=8000),
  the data behind `figures/shum_qnm_reproduction.png`. Regenerate via
  `shum_evolve_opt.jl run_shum(0.04, 8000.0)`. Other `r5_eps_Dr*.txt`
  (scratch resolutions) are gitignored.

## Promotion roadmap
These are verified, reusable reproduction modules. The next code-clarity step
(per `progress/CLAUDE.md`) is to promote the engines into the package `src/`
as proper submodules — `src/perturbations/AxialQNM.jl` (axial_*),
`src/cowling/CowlingBDNK.jl` (shum_*), `src/viscous/ViscousBDNK1D.jl`
(pmp_viscous_core) — moving each validation block into `test/`. Deferred here as
a deliberate refactor to avoid destabilising the verified state in a cleanup pass.
