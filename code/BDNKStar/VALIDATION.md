# BDNKStar — Validation & Scope

Companion to `README.md`. This document states (1) the unit conventions, (2) what has
been validated against published results and to what tolerance, (3) reproduction
coverage, (4) known limitations and honest negatives, and (5) which results are novel
rather than reproductions, with their maturity. Every number below is either measured
by the test suite / a script in this repo, or cited to its source.

---

## 1. Test-suite status

**1472 / 1472 passing** (full suite, `julia --project=. test/runtests.jl`, 29 m 22 s at
6 threads). 45 test files.

**Reproduction scripts: 58 / 58 execute with no errors** (audited 2026-07-19, wall-clock at
3 threads/script). 50 finish < 330 s; 6 are slow-but-complete (330–1335 s); 2 are finite
τ̂-sweeps exceeding 1500 s (`axial_frame_independence`, `pmp_shock_instab`) — very long, not
hung. A driver, **`repro/RUNALL.jl`**, runs them in tiers: fast lane per-commit
(`julia --project=. repro/RUNALL.jl`), `--all` for the heavy set (nightly), `--xheavy` for the
two long sweeps (manual). It exits non-zero only on a genuine error, so it is CI-usable to keep
the "figures reproduced" claim continuously verified. **Note: the reproduction scripts are NOT
part of `runtests.jl`** — the 1472 tests validate the *library*; `RUNALL.jl` validates the
*reproductions*. Wire both into CI before release.

---

## 2. Unit conventions — read this first

`solve_tov` is **scale-agnostic**: it integrates in whatever geometric units the EOS and
central density are expressed in. Two conventions are in use, and confusing them changes
every mass and frequency by exactly `Msun_to_km = 1.4766`:

| convention | set by | `star.M` | `star.R` | mass → M⊙ | freq → kHz |
|---|---|---|---|---|---|
| **km-geometric** | εc in km⁻² (e.g. `piecewise_polytrope` via ρ[g/cm³]·7.4237e-19) | km | km | `mass_solar(star)` | `f[km⁻¹]/kHz_to_km` |
| **M⊙-geometric** | G=c=M⊙=1 (NS-oscillation literature standard) | already M⊙ | M⊙ (×1.4766 → km) | **identity** — do *not* call `mass_solar` | `f·Msun_to_km·kHz_to_km` |

Because the Γ=2 (n=1) polytrope is **scale-invariant** (M ∝ √κ), the *same* raw
`ShumPolytrope(100.0)` TOV solution legitimately reproduces **two different published
benchmark stars** depending on the units assigned to κ:

| reading | κ | M | R | anchor |
|---|---|---|---|---|
| M⊙-geometric | 100 M⊙² | **1.400 M⊙** | **14.154 km** | Font/Stergioulas/Kokkotas ℓ=2 benchmark (`test_nonradial`, `test_spherical`) |
| km-geometric | 100 km² | **0.948 M⊙** | **9.586 km** | Kokkotas & Ruoff 2001, A&A 366, 565, Tab. A.18 radial band (`test_dyngr`) |

Both are physically valid and each module is internally consistent with the reference it
targets. **Dimensionless results (compactness M/R, Mω, Λ, Q, damping exponents) are
identical under both readings** and are unaffected. Conventions are documented on
`ShumPolytrope` (`src/eos/EquationOfState.jl`) and `mass_solar` (`src/tov/TOV.jl`).

---

## 3. Quantitative validation against published results

Tolerances below are the ones **enforced by the test suite** (all currently pass).

| # | Benchmark | Published reference | Enforced tolerance / result | Where |
|---|---|---|---|---|
| 1 | n=1, K=100 TOV star | Font/Stergioulas/Kokkotas | M = 1.400 M⊙, R = 14.154 km vs published 1.4 / 14.15 km (atol 0.1 km) | `test_nonradial.jl:23` |
| 2 | ℓ=2 Cowling f + p₁…p₄ spectrum | Font/Stergioulas/Kokkotas table | **every mode within 0.3 %** (`rtol=3e-3`) | `test_nonradial.jl:19` |
| 3 | Radial fundamental (full-GR, Chandrasekhar LAWE) | Kokkotas & Ruoff 2001 Tab. A.18 | F = 2.124 kHz, inside the published ν₀ band (asserted 2.0–2.45 kHz) | `test_dyngr.jl` |
| 4 | Radial spectrum F | in-repo cross-solver anchor | F = 2.0248 ± 0.01 | `test_radial_spectrum.jl:277` |
| 5 | Λ(1.4 M⊙), SLy/APR4/H4/MS1 | Read et al. piecewise polytropes | within **5 %** of published Λ(1.4) | `test_tidal.jl:65` |
| 6 | I–Love universal relation | Yagi & Yunes | on the universal curve to **< 3 %**; EOS-independent to < 3 % | `test_tidal.jl:86,91` |
| 7 | Ī–C sequence | Yagi & Yunes | dev < 8 %; measured worst 0.021 (SLy), 0.052 (APR4), 0.049 (H4) | `test_universal_relations.jl:153` |
| 8 | Ī EOS-collapse at fixed C | Yagi & Yunes | spread < 8 %; measured 0.071 (C=0.16), 0.070 (0.18), 0.035 (0.20) | `test_universal_relations.jl:175` |
| 9 | f-mode vs AK1998 relation | Andersson & Kokkotas 1998 | full-GR anchor, fit residual ~1e-7 (see §5 for the known calibration offset) | `test_universal_relations.jl` |
| 10 | r-mode GW coefficient | Lindblom–Owen–Morsink 1998 | `RMODE_GW_COEFF` = 2.5104 ± 1e-3 | `test_rmodes.jl:61` |
| 11 | r-mode τ_GW, τ_sv, τ_bv (n=1) | LOM98 Table I | reproduced with calibrated bulk eigenfunction | `src/rotation/RModes.jl` |
| 12 | r-mode instability window | LOM98 / realistic EOS | M = 1.4 ± 0.06 M⊙; window minima spread < 0.03 | `test_rmodes.jl:148,158` |
| 13a | TOV star, EOS1 (κ=100 km², n=1, ρc=3e15 g/cm³) | Bussières et al. 2604.13208 | M = 1.27 M⊙ (atol 0.01), R = 8.86 km (atol 0.02) — direct published values | `test_tov.jl:6` |
| 13b | Axial w-mode QNM, inviscid | Bussières Tab. II | f = 10.50 kHz, τ = 29.54 µs, both **rtol 1 %** | `test_axial_viscous.jl:18` |
| 13c | Axial w-mode QNM, viscous (BDNK frame A, η_c=1e31 cgs) | Bussières Tab. II | f = 10.0898 kHz, τ = 30.8857 µs, both **rtol 1 %**; viscosity lowers f | `test_axial_viscous.jl:24` |
| 14 | BDNK causality / frame well-posedness | Shum et al. 2509.15303 | production-frame speeds real, ordered, subluminal at every radius | `test_causality_map.jl`, `test_shum_frame.jl` |
| 15 | BDNK viscoresistive MHD constitutive | Lier, Armas & Porth 2606.22691 (Eq. 5,6,8,20) | **verbatim** constitutive; analytic Eq.20 Jacobian byte-identical to the finite-basis form to ~1e-15 | `src/mhd/BDNKMHDConstitutive.jl` |
| 16 | PMP telegrapher / heat / shock | PMP 2209.09265 | 7 figures reproduced | `repro/pmp_*.jl` |
| 17 | **Viscous damping, cross-method (eigenvalue vs dissipation integral)** — internal, and the strongest self-consistency check in the package | in-repo, two independent methods | p₁: eigenvalue **slope** dγ/dη̂ vs the linear dissipation integral agree to **85.4 %** (asserted > 75 %, `rtol 0.25`); dissipation-integral values γd_f = 0.0875, γd_p = 0.585 (rtol 5 %); ordering γd_p/γd_f > 4 | `test_cross_method.jl:125-153` |

> **Note on row 17 — this is the anchor that governs viscous-damping claims.** It fits a *slope*
> (immune to an additive offset) over η̂ ∈ [0.01, 0.02, 0.04] at Nr=110, i.e. inside **both** the
> documented Nr≈96–160 window **and** the η̂ < τ̂ ≈ 0.12 recovery limit. Within that regime the
> damping is **linear in η̂** and matches the dissipation integral for p₁ to ~15 %. The suite also
> explicitly asserts `mf_win > γd_f` — i.e. **the ℓ=2 f-mode eigenvalue is known to be
> over-estimated/contaminated**, and only its qualitative behaviour is claimed. Any future
> viscous-damping claim should be checked against this test *first*; the retracted √ν result
> (§6) was built on the contaminated f-mode over an η̂ range extending past the validity limit,
> and this pre-existing test already contradicted it.

---

## 4. Reproduction coverage

Per-paper figure counts (`viz/coverage_map.jl` → `figures/coverage_map.png`):
**29 original figures reproduced** across 6 reference papers.

| paper | reproduced | blocked | blocking reason |
|---|---|---|---|
| Kovtun 1907.08191 | 7 | 0 | — |
| Pandya 2201.12317 | 7 | 4 | 2D: no reference code |
| PMP 2209.09265 | 7 | 0 | — |
| Shum 2509.15303 | 4 | 2 | fine-Δr origin instability |
| Bussières 2604.13208 | 4 | 1 | configuration halt |
| Chabanov–Rezzolla 2311.xxxxx | 0 | 26 | 3+1D GRMHD merger (GPU) |

---

## 5. Known limitations and honest negatives

These bound what may be claimed. None is hidden by the test suite; several are asserted
*as* limitations.

**Numerics / stability**
- **Stiff BDNK-MHD shock tubes** (high `D_u`, cases ST-g…j) blow up under explicit
  time-stepping; this is a genuine explicit-stiffness limit, not a CFL or dissipation
  tuning issue. **An IMEX integrator is not implemented.** Tests are split so a–f run
  nonlinearly and g–j are asserted only as linearly well-posed.
- **2D Kelvin–Helmholtz evolution** is superluminally unstable; front velocities are
  therefore validated **algebraically**, not by evolution.
- **Time-domain magnetized mode extraction** is unstable for B₀ ≳ 0.03 (local
  p_mag > p_gas); magnetized mode frequencies come from the **eigenvalue** solvers.
- **Spherical (r,θ) full-star evolution** requires graded excision; the linearized
  engine exhibits slow grid-scale growth (see below).

**Quantitative accuracy of dissipation**
- **The polar viscous QNM operator has no physical outer boundary condition, and this is a
  real defect — not a tolerance.** `PolarViscousModes.jl:63-64` imposes a bare zero-gradient
  (Neumann) ghost on all five fields, at a surface that is **not even a grid point**
  (`SphBackground.jl:71` puts the outermost cell centre at R−dr/2, while the docstring at
  `:60-62` says the surface lies *on* the last point). The result is a **one-cell dissipation
  artifact** contributing an η̂-independent offset to γ. Consequences, all measured:
  a zero-gradient ghost on a truncated free surface **mimics a no-slip wall**; the fitted
  damping exponent moves 0.42→0.99 with the ghost rule alone and 0.11→0.64 with the numerical
  floor `den_frac` alone; and the module's self-reported "~5–8 % free-surface systematic" and
  "+20–30 % reactive viscous pull" are consequences of this BC, not frame effects. **Damping
  sign/trend is trustworthy; the shipped η̂-scaling and the absolute rate are not.**
- **The polar viscous operator uses a constant scalar viscosity**, not the variable-coefficient
  ∇_μ(2ησ^{μν}) with η(r) = η̂(ε₀+p₀) → 0 at the surface. Any claim resting on the surface
  behaviour of η(r) is therefore untested by this solver as shipped.
- **Validity window.** `polar_qnm` is documented for **Nr ≲ 200** (`PolarViscousModes.jl:105-106`,
  recommending Nr≈96–160) with γ drifting ~tens-of-% over Nr=96–200; beyond Nr≈200 the mode
  continuation jumps eigen-branches. Separately, η̂ must satisfy the BDNK recovery condition —
  at η̂=0.16 > τ̂=0.12 the recovery denominator goes negative and is floored over the **entire**
  star, so η̂ ≳ 0.12 sweeps are outside the operator's valid regime.
- **Polar viscous QNM damping rates are semi-quantitative** (tens-of-% level). Frequencies are
  better converged than damping rates, but the inviscid ℓ=2 f-mode is itself only ~4.6 % from its
  Richardson limit at Nr=400.
- **2+1D time-domain viscous damping is qualitative only.** In `SphBDNK` the shear
  friction `ν_mom∇²δS` demonstrably removes mode energy monotonically in η̂, but the
  extracted *rate* is contaminated by Kreiss–Oliger dissipation, grid-scale growth and an
  unconverged eigenfunction (~100× below the analytic ν_mom·k²/2; linearity R²=0.21).
  The weakly damped ℓ=2 **fundamental** cannot be measured this way at feasible
  resolution — the ℓ=2 quadrupole moment is dephasing-dominated and the *energy* is the
  only faithful diagnostic. Quantitative viscous damping must come from the
  frequency-domain solvers.
- **1+1D radial viscous damping** has a numerical floor γ₀ ≈ 6.2e-4 (Q ≈ 68) from the
  LF/HLL scheme, which must be subtracted.

**Physics scope**
- **No temperature-dependent transport.** Transport is parametrized by a constant η̂;
  there is **no microphysical ζ(ρ,T)** outside the r-mode module (`RModes.jl` carries the
  Sawyer modified-Urca ζ ∝ ρ²T⁶ and Cutler–Lindblom η ∝ ρ^{9/4}T⁻², not exposed to the
  transport layer or the evolvers). All merger-temperature transport statements are
  therefore **literature-cited or order-of-magnitude estimates, not in-code results**.
- **No 3+1D GRMHD.** The Chabanov–Rezzolla merger figures are entirely out of scope.
- **Cowling-vs-GR f-mode calibration offset**: Cowling f-mode frequencies carry a known
  ~10–20 % offset vs full GR; the f-mode/AK1998 comparison offset was diagnosed as this
  calibration effect (Cowling-vs-GR plus the AK fit), **not** a nonequilibrium effect.
- **g-modes require stratification.** A strictly barotropic EOS has none; all g-mode
  frequencies, overlaps and damping used anywhere are model inputs.
- **`mass_solar` is convention-dependent** — see §2.

---

## 6. Novel results (not reproductions) and their maturity

| result | status | what it still needs |
|---|---|---|
| **f-mode viscous damping is LINEAR in η̂ (bulk, p ≈ 1)** — an earlier sub-linear claim (p ≈ 0.47–0.50, "γ∝√ν") is **RETRACTED** as a boundary artifact | **RETRACTED — see the retraction box** | Coefficient still unsettled to ~factor 1.8; needs a derived free-surface BC and the true variable-coefficient operator (below) |

> ❌ **RETRACTION (2026-07-19) — the "f-mode damps as γ ∝ √ν" claim is withdrawn.**
> A dedicated audit (convergence × dissipation integral × eigenfunction diagnostic × BC sensitivity, with
> adversarial verification) established that the sub-linear exponent is a **numerical artifact of the polar
> solver's outer boundary treatment**, not physics. The evidence:
>
> 1. **The layer the claim requires does not exist.** p = 0.5 demands an η̂-dependent viscous layer of width
>    δ ∝ √η̂. Measured directly, the dissipation is **grid-locked to exactly one cell at every η̂** — fixed
>    *cell* width, not fixed physical width, and it does not scale as √η̂.
> 2. **The exponent is a function of unphysical knobs.** It swings **0.42 → 0.99** when only the outer ghost
>    rule is changed, and **0.11 → 0.64** when only the numerical floor `den_frac` is changed.
> 3. **The operator does not implement the invoked mechanism.** The claim appealed to η = η̂(ε+p) → 0 at the
>    surface; the shipped operator uses an **r-independent constant** `νm`, never η(r).
> 4. **The exponent never converges.** Fitted p over η̂∈[0.04,0.16] vs Nr = 100,140,200,260,320,400 is
>    0.560, 0.498, 0.436, **−0.169, 0.031, 0.238** — it passes *through* 0.5 while drifting, and γ itself still
>    moves 34–37 % from Nr=320→400. The solver's own docstring (`PolarViscousModes.jl:105-106`) limits it to
>    **Nr ≲ 200**, so the original Nr=140/200 points sat at the documented edge, and the continuation jumps
>    eigen-branches beyond it (the viscous frequency shift flips sign between Nr=200 and 260).
> 5. **The cross-method "confirmation" was not independent.** `DynGR1D` hard-resets v = 0 in the atmosphere —
>    a numerical **no-slip rigid wall**, and √(νω) is precisely the *rigid-wall* Stokes law. The 1D l=0 √ν was
>    produced by the boundary treatment, not by stellar physics.
> 6. **An η̂-independent offset plus linear damping fits 13–18× better** than a power law — that offset is the
>    actual generator of the spurious p ≈ 0.5.
>
> **What replaces it:** the bulk exponent **p ≈ 1** is robust across four independent routes — the dissipation
> integral (exactly linear by construction), a Dirichlet outer BC (0.985), an extrapolation BC (0.966), 3/6/12-cell
> surface masking (1.004/1.004/1.003), and the independent axial shooting solver (0.959). This agrees with the
> classical free-surface result γ = 2νk². The **coefficient** is *not* settled better than a factor ≈1.8.
>
> **Follow-up experiment (2026-07-19), and it refined the diagnosis.** A dedicated `atm_vbc` knob was added to
> `setup_dyngr` (`:zero` = historical pinned-v atmosphere wall; `:outflow` = zero-gradient, no wall stress) and the
> 1D sweep re-run changing *only* the BC. The no-slip hypothesis **was not confirmed**: removing the wall did not
> restore p≈1, it pushed the apparent exponent *down* (local p 0.23,0.45 → 0.03,0.22). The data instead show that
> an **affine model — a ν-independent offset plus strictly linear damping — fits far better than any power law**:
> `:zero` γ = 2.056e-3 + 0.0262·ν (RSS 3.0e-9) vs A·ν^0.337 (RSS 3.9e-8), i.e. **12.9× better**. Crucially the
> fitted **linear slope 0.0262 agrees with the independent bulk dissipation-integral prediction 0.0224 to 17 %**.
> So the underlying damping is genuinely **bulk (γ∝ν)** and the apparent "√ν" was an **offset masquerading as a
> power law** — the mechanism, now measured directly in a second code. The wall is not the generator, though it
> does supply ~65 % of the real viscous slope (0.0262 → 0.0091 when removed).
> **Open confound:** the offset appears only for ν>0 (ν=0 shows no measurable damping), so it is not a static
> numerical floor; the leading suspect is the parabolic dt cap (dt ∝ 1/ν, active for every ν ≥ 0.02 here) changing
> the scheme's dissipation. No damping *coefficient* should be quoted from this engine until that is quantified.
>
> *Terminology note, also corrected:* √(νω) Stokes damping is the **rigid-wall/no-slip** law; for a **free
> surface** the classical answer is bulk, γ = 2νk². An earlier draft of this file attributed our result to a
> "classical Landau–Lifshitz free-surface Stokes layer," which conflated the two.
| **BDNK causal resistive QPO damping**, r_eff(ω) = r_b/(1+(ωτ_b)²) | prediction | a validated resistivity/microphysics model |
| **Dynamical-tide reversible/irreversible split** — dissipation enters the reversible (universal) channel only at O(Kn²), the irreversible channel at O(Kn¹) | analytically clean | absolute numbers rest on semi-quantitative QNM damping |
| **Merger: irreversible tide becomes important via the Urca ζ(T) resonance** | order-of-magnitude | an in-code ζ(ρ,T) module — currently literature-only |

---

## 7. 3+1D engine mode census

Everything in §§3–6 above that carries a frequency was produced by the **1D radial solvers**
(`NonRadialModes`, `PolarGRModes`, `AxialViscousModes`, `TidalDeformability`, `GravityModes`,
`RadialModes`) on spherically symmetric TOV backgrounds. This section records, separately and
completely, every mode that was actually extracted from a **3+1D engine** (2026-09-10/11;
4 compute blocks, each re-run bit-for-bit by an adversarial verifier). Voided rows are kept:
a table that omits what failed misrepresents the engine.

**Anchor star throughout:** `ShumPolytrope(100.0)`, ε_c = 0.00128 + 100·0.00128², M = 1.40016 M⊙,
R = 14.1544 km. 1D reference recomputed in-session: f = 1.88291, p₁ = 4.10672, p₂ = 6.03113 kHz.

### 7.1 What the 3+1D engines cannot compute — do not look for it below

All three engines (`CowlingEvolve3D`, `CowlingBDNK3D`, `DGCart3D`) are **Cowling: the metric is
frozen.** There is therefore no gravitational radiation and no metric perturbation, so the
following are **structurally absent**, not merely unconverged: axial w-modes, polar GR damping
times τ_f, tidal deformability Λ (hence f-Love and w-Love), r-modes (no rotation), g-modes (no
stratification), i-modes (no density discontinuity). Every figure and table in §§3–6 that uses
these quantities is 1D-sourced permanently.

### 7.2 The single number that governs every row

Time-domain frequency extraction has spectral resolution df = 1/T. At affordable 3D resolution
the usable record is **23 f-mode periods** (linear/BDNK engines, df/f = 4.4%) or **4.6 periods**
(`DGCart3D`, df/f = 21.6% — a cap set at the time not by cost but by the engine's background
eroding 10% in that time; the erosion is diagnosed and fixed in §7.8, and the §7.5–7.6 rows below
are from the pre-fix engine). **Every sub-percent figure below is a parametric peak estimate** (parabolic vertex
+ matrix pencil), not a resolved spectral line, and is quoted only where the two estimators
agree. No third decimal is meaningful on any 3D frequency.

### 7.3 Linear Cowling engine — `CowlingEvolve3D`, ℓ=2, m=0, cut-cell surface

Record 23 periods, df/f ≈ 4.4%, dt = 0.30 dx, κ_min = 0.5. Deviation is against the 1D shooting
value for the identical star.

| EOS (M/M⊙, R km) | N | σ_ko | mode | f₃D [kHz] | f₁D [kHz] | dev | status |
|---|---|---|---|---|---|---|---|
| Γ=2 anchor (1.400, 14.154) | 24 | 0.01 | f | 1.88301 | 1.88291 | +0.005% | ok |
| | 32 | 0.01 | f | 1.87560 | | −0.388% | ok |
| | 40 | 0.01 | f | 1.87701 | | −0.313% | ok |
| | **48** | 0.01 | **f** | **1.88135** | | **−0.083%** | ok — best anchor |
| | 48 | 0.01 | p₁ | 4.0665 | 4.10672 | −0.98% | ok |
| | 48 | 0.02 | f | 1.89076 | 1.88291 | +0.417% | σ_ko systematic, brackets from above |
| | 32, 40 | 0.004 | f | 1.75235, 1.75205 | | **−6.93%, −6.95%** | **cut=false — the `setup_evo3d` DEFAULT — wrong and non-converging** |
| Γ=2.5 polytrope | 48 | 0.01 | f | 2.23448 | 2.23404 | **+0.020%** | ok — best of study |
| | 48 | 0.01 | p₁ | 5.8517 | 5.9007 | −0.83% | ok |
| Γ=3 polytrope | 48 | 0.01 | f | 2.40545 | 2.40739 | −0.080% | ok |
| | 48 | 0.01 | p₁ | 7.0718 | 7.1457 | −1.03% | ok |
| SLy (1.4205, 11.675) | 32 | 0.01 | f | 2.37154 | 2.41928 | −1.973% | ok — highest stable N |
| | 24→32 Rich. | 0.01 | f | 2.39450 | | −1.024% | monotone pair, Richardson valid |
| | 40, 48 | 0.01 | f | **VOID** | | — | surface mode +0.71, +1.23 e-folds/period; end/start 4×10⁴, 4×10¹¹ |
| | 48 | 0.04 | f | 2.31598 | | −4.270% | stable but KO-biased; bias grows with N |
| APR4 (1.3595, 11.309) | 32 | 0.01 | f | 2.45373 | 2.49463 | −1.640% | windowed to 12.2 periods |
| | 24→32 Rich. | 0.01 | f | 2.48193 | | **−0.509%** | best crusted-EOS result |
| | 40, 48 | 0.01 | f | **VOID** | | — | +0.89, +1.59 e-folds/period; end/start 8×10¹⁴ (worst in study) |
| | 48 | 0.04 | f | 2.38693 | | −4.317% | KO-biased |
| H4 (1.4059, 13.926) | 32 | 0.01 | f | 1.89169 | 1.94832 | −2.907% | marginal (−0.02 e-folds/period) |
| | 24→32 Rich. | 0.01 | f | 1.89821 | | −2.572% | Richardson barely helps |
| | 48 | 0.04 | f | 1.83965 | | −5.578% | KO-biased |
| MS1 (1.3588, 14.895) | 24 | 0.01 | f | 1.69672 | 1.76756 | −4.008% | **only stable resolution** at σ_ko=0.01 |
| | 40, 48 | 0.04 | f | 1.64607, 1.64294 | | −6.873%, −7.050% | KO-biased; largest bias in study |
| all 4 crusted EOS | 24, 32 | 0.01 | p₁ | 4.67–6.43 | 5.28–7.26 | −10% to −18% | **not resolved** |

**Reading this table.** Smooth analytic barotropes converge to <0.1% at N=48, but the N=24→48
sequence is non-monotone (+0.005/−0.39/−0.31/−0.08%), so the honest error bar is the N-spread,
±0.4%, and Richardson on (40,48) makes every analytic EOS *worse*. The four realistic crusted
EOS excite a growing surface mode for N ≥ 40 at σ_ko = 0.01. The discriminant is the surface
index n in ε ∝ (R−r)ⁿ — 0.50/0.67/1.005 (stable) versus 1.84/1.87/1.95/1.95 (unstable),
monotone over 19 decades of end/start ratio. Mechanism **open**: the atmosphere-sound-speed
hypothesis was tested (all 77 040 atmosphere cells zeroed, SLy N=48) and refuted (+0.596 vs
+0.595 e-folds/period). No affordable configuration gives a converging crusted-EOS f-mode.

### 7.4 BDNK viscous engine — `CowlingBDNK3D`, ℓ=2, m=0, three Clarisse frames

η̂ = 0.03, ζ = 0, cut-cell, σ_ko = 0.01, dt = 0.20 dx. Frames per Clarisse et al.
(arXiv:2510.16603 Tab. 1) mapped cell-by-cell as τ_ε = a₁η/w₀, τ_Q = a₂η/w₀, τ_P = c_s²τ_ε (the
last is pinned by the module's Π = c_s²𝒜 closure, not free). Damping by the difference protocol
γ_visc = γ(η̂=0.03) − γ(η̂=0) at identical N and frame, which cancels the numerical floor.

| frame (a₁, a₂) | N | T [M⊙] / periods | f₃D [kHz] (pencil / pgram) | f₁D polar | γ_visc [10⁻³/M⊙] | τ_damp, Q | BDN bound violated | stable? |
|---|---|---|---|---|---|---|---|---|
| F₁ (25/4, 25/7) | 40 | 2500 / 22.8 | 1.87954 / 1.88027 | 1.9038 | 2.48251 | 1.984 ms, 11.72 | **100.0%** of cells | **yes** (end/start 6.4×10⁻⁵) |
| F₂ (25/2, 25/3) | 40 | 2500 / 22.8 | 1.87973 / 1.88036 | 1.9105 | 2.48229 | 1.984 ms, 11.72 | 94.4% cells (94.9% vol) | yes |
| F₃ (25, 25) | 40 | 2500 / 22.8 | 1.88040 / 1.88082 | 1.9272 | 2.48139 | 1.985 ms, 11.73 | 51.3% cells (51.3% vol) | yes |
| η̂=0 control, all frames | 40 | 2500 | 1.87681 (identical to 6 digits) | 1.887 | — | — | — | yes |
| F₁, F₃ long record | 32 | 5000 / 46.4 | — | — | decay −4.4×10⁻³/M⊙ constant to 1% | — | — | **stable through 12 decades** |

| derived quantity | N=24 | N=32 | N=40 | note |
|---|---|---|---|---|
| frame spread of f | 0.052% | 0.048% | **0.046%** | vs 1.22% in the 1D polar operator; upper bound — df/f = 4.4% |
| frame spread of γ_visc | 0.150% | 0.130% | **0.045%** | shrinks with N; 1D gave 14.28% (withdrawn — missing lapse) |
| γ_visc, F₁ | 2.35482 | 2.44941 | 2.48251 | order ≈ 3 at σ_KO=0.01; at σ_KO=0.005 the pencil rate is 2.597 (N=40) / 2.607 (N=48) = **99.2% of the NS integral 2.6238** — the 96–97% was KO-biased |
| numerical floor γ(η̂=0) | 3.2697 | 1.9551 | 1.2993 | ∝ dx^1.8; still 52% of γ_visc at N=40 — difference protocol mandatory |
| viscous frequency shift | +0.40% | +0.28% | +0.15% | **KO artifact** — at σ_KO=0.005 it converges to the damped-oscillator value −(γ²−γ₀²)/2ω²: −0.125% measured vs −0.133% predicted at N=48 (`repro/bdnk3d_viscous_shift.jl`) |

**Two physics results from this block.** (i) The stability bound τ_Q w₀ c_s² > (4/3)η reduces
under the Clarisse mapping to the η-independent a₂c_s² > 4/3, and **its violation is sufficient-
not-necessary**: F₁ violates it in every cell and is as stable as F₃. The local 4×4 sound-sector
analysis shows why — the unstable band lies above k_c = C/η̂ (15.3 / 10.1 / 5.8 per M⊙ for
F₁/F₂/F₃), which the 2nd-order stencil only resolves for N > 353 / 232 / 135. Moved to the
cheap η̂ axis at N=24 the criterion predicts 6/6 blow-ups correctly in sign and to 24–28% in
growth rate, and **F₃ destabilises first** (k_c 2.6× smaller) despite violating in fewer cells.
(ii) a₁ is inert: a₁ = 6.25/25/100 at fixed a₂ gives f = 1.89112/1.89113/1.89123 kHz. The three
frames differ through a₂ alone, in 3D as in the axial sector.

### 7.5 Nonlinear DG engine — `DGCart3D`, ℓ=2, m = 0 and 2, octant grid, p=2 (pre-fix engine)

**These rows were produced with the mean-based positivity limiter that eroded the background
(§7.7 item 5, fixed in §7.8).** The m-degeneracy contrasts stand as recorded; the absolute
frequencies do not (see §7.8 for the fixed engine's f). The only engine that can carry a
non-axisymmetric perturbation. Analysis window [0, 600] M⊙ =
4.6 periods at every K, df/f = 21.6%. Seed amplitude A = 10⁻² (v^r = A(r/R)Y).

| K (nodes/axis) | f(m=0) [kHz] | f(m=2) [kHz] | dev vs 1D | ∣Δf∣/f (m=2 vs 0) | ρ_c drift by t=600 | status |
|---|---|---|---|---|---|---|
| 4 (12) | 1.85758 | 1.49994 | −1.3% / −20.3% | 19.25% | +6.1%, then collapses | **void** — background collapses |
| 5 (15) | 1.47898 | 1.57348 | −21.5% / −16.4% | 6.39% | | coarse |
| 6 (18) | 1.61488 | 1.55356 | −14.2% / −17.5% | 3.80% | | usable |
| 8 (24) | 1.58933 | 1.55637 | −15.6% / −17.3% | 2.07% | −13% | usable |
| 10 (30) | 1.57218 | 1.58496 | −16.5% / −15.8% | **0.81%** | −8.3% / −12.8% | best |
| 6, Y₂₁ seed, one run | 1.56352 (q₂₀) | 1.56366 (q₂₂) | −17.0% | **0.0090%** | | same-run channels: degeneracy ≤ 10⁻⁴ |

**Reading this table.** On a non-rotating star the m-degeneracy is exact, so every splitting
is a grid artifact; it falls as h^3.2 (K=4–10), consistent with the scheme's O(h³). The
within-run Y₂₁ bound (0.009%) is 90× tighter than the best cross-run value and shows the grid's
genuine degeneracy violation is ≤10⁻⁴ — the cross-run "splitting" is run-to-run systematics
(background drift, 4.6-period record, nonlinear back-reaction). **Two caveats are load-bearing:**
the absolute f is −14 to −18% at every K (common-mode between m=0 and m=2; diagnosed in §7.8 as
the eroding background — the fixed engine gives −2 to −10%, its residual staircase-surface
systematic); and **m=1 is structurally absent** — Re Y₂₁ ∝ n_x n_z has parity
(−,+,−) under the octant's three separate reflections, so the whole T₂g triplet cannot live on
this grid. The Y₂₁ seed's mirror extension projects onto the E_g doublet {Y₂₀, Re Y₂₂}
(q₂₂/q₂₀ = 1.0000 measured, confirmed to 1.2×10⁻⁴), which is why it is a degeneracy diagnostic
and not an m=1 measurement. `dgcart3d_quadrupole_m2` discriminates m=2 from m=0 by 3×10⁵–5×10⁶
at every K.

### 7.6 Nonlinear amplitude dependence — `DGCart3D`, K=8, p=2, window [0, 600]

| seed A | ξ_max/R | f(A)/f(A→0) − 1 | vs 0.82% floor | status |
|---|---|---|---|---|
| 0, 10⁻⁴, 3×10⁻⁴ | — | — | SNR 0.33, 0.62, 2.06 | unusable (below engine's own q₂ noise) |
| 10⁻³ | 0.0041 | +5.04 ± 8.99% | — | not resolved (SNR 4.3) |
| 3×10⁻³ | 0.0123 | −0.85 ± 0.74% | 1.0× | not resolved |
| 10⁻² | 0.0412 | +0.54 ± 0.33% | 0.7× | not resolved |
| 3×10⁻² | 0.124 | +0.31 ± 0.57% | 0.4× | not resolved (predicted −0.16%) |
| 10⁻¹ | 0.412 | −1.38 ± 1.02% | 1.7× | **suggestive, not resolved**; K=6 gives −7.48 ± 0.90% — **not converged** (5.4× change for 25% refinement) |
| 3×10⁻¹ | 1.23 | −10.1 ± 7.2% | — | invalid: limiter-clipped (A_eff/A = 0.28) |

Second harmonic at 2f: **not detected** in either the ℓ=2 or ℓ=0 channel, A(2f)/A(f) < 0.09–0.19
at ξ/R = 0.41. Combination tone with p₁: not measurable (p₁ unresolved in this engine). Linear-
engine control: `CowlingEvolve3D` reproduces exact amplitude proportionality to 10⁻¹⁴, so the
pipeline floor is ~10⁻¹⁶ and nothing above is a diagnostic artifact. **Verdict: no nonlinear
frequency shift is resolved at affordable cost**; the one suggestive point collapses with
resolution, the signature of limiter/Rusanov dissipation rather than fluid nonlinearity. K=12
(first resolution with a quiet background, drift −6×10⁻³ at t=500) is the experiment that would
settle it and was not affordable. With the fixed limiter (§7.8) the background is quiet at every
K (drift −5×10⁻⁴ at t=600 for A=10⁻² at K=6), so the experiment is now affordable; it has not
been rerun.

### 7.7 Engine defects found during this census (status noted per item)

1. **`setup_evo3d` default `cut=false` gives −6.95%** on f from a run that looks healthy (stable,
   55-period record, window-independent peak) and does not converge away between N=32 and 40.
2. **σ_ko is an uncontrolled systematic**: ±0.25% on f with cut=true, >7 pp with cut=false. Any 3D
   frequency must state σ_ko.
3. Growing late-time surface mode on the cut-cell path at σ_ko = 0.004 for N ≥ 32 (growth
   0 → 3.3×10⁻⁴ → 2.5×10⁻³ → 8.8×10⁻³ /M⊙ at N=24/32/40/56).
4. `Background3D.build_star3d:116` sets ρ₀ = ε − p, exact only for Γ=2. Harmless for the linear
   engines (they never read ρ₀); would bite any nonlinear user.
5. **`DGCart3D`'s background is not static** — **FIXED (§7.8)**. As found: ρ_c eroded −10% by
   t=600 M⊙, −31 to −43% by t=3000, absolute f −16% vs 1D at every K, which capped the record at
   4.6 periods and voided the m-degeneracy and nonlinear tests as precision measurements. Cause:
   the well-balanced RHS subtraction is exact, but the mean-based Zhang–Shu positivity limiter
   flattens the equilibrium's own intra-element profile whenever it engages, and the Rusanov
   dissipation acts on the equilibrium jump with a state-dependent wave speed. Both now act on
   the deviation from equilibrium; the equilibrium is a bitwise fixed point of limiter + RHS.
6. Numerical damping floor γ_num ≈ 1.2×10⁻³ /M⊙ at N=48 — a hard floor under any viscous
   measurement without the η=0 difference protocol.
7. Julia: `2f1` parses as the `Float32` literal `20.0`, not `2*f1`. A harmonic search silently
   sampled 15–30 kHz instead of 2f and returned a fake null until caught by synthetic injection.

### 7.8 `DGCart3D` limiter fix — equilibrium-preserving Zhang–Shu, deviation-form Rusanov (2026-09-17)

**Attribution.** With the limiter switched off the unseeded star is static to 10⁻¹² (the RHS
subtraction is exact); with the mean-based limiter on, the same unseeded star is thrown into a
±3% radial oscillation with 0.6c surface velocities within 50 M⊙. The limiter was the erosion
source. Mechanism: Zhang–Shu scales every field toward its cell mean, U ← Ū + θ(U−Ū); in a
surface element D_eq spans orders of magnitude, so any θ<1 flattens the equilibrium's own
profile and pushes mass from the stellar edge into the atmosphere nodes, every stage of every
step.

**Fix** (`_limit_wb!`, `_rus5wb`; `wb=true` whenever a stellar background is stored). The
limiter scales about a reference R that (i) equals U_eq when the deviation δU = U − U_eq
vanishes and (ii) has the cell mean of U: R_D = D_eq + δD̄·w_D, R_S = S_eq + δS̄·w_S,
R_τ = τ_eq + ΔK + (δτ̄ − ΔK̄)·w_τ, with mean-one weights w ∝ the background (w_D = D_eq/D̄_eq,
w_τ = τ_eq/τ̄_eq, w_S ∝ D_eq − D_floor) and ΔK the cold kinetic energy √(R_D²+|R_S|²) − R_D of
the reference momentum. θ enforces D ≥ ½√γρ_atm and the pressure proxy q ≥ 10⁻¹² q(R) at every
node (chord root; q is concave so the bound is guaranteed). Only an element whose reference is
itself infeasible (lost more than half its mass or all of its thermal energy) falls back to the
mean-based limiter. The Rusanov flux dissipates a_max on δU only and a frozen equilibrium sound
speed on the equilibrium jump. Two false starts are recorded in the source comments: a uniform
shift R = U_eq + δŪ swamps the atmosphere nodes and fell back in every surface element
(ρ_c −3.6% at t=600); a reference without ΔK has negative pressure at the outermost stellar nodes
for a 1% seed (τ_eq ∝ ρ² vanishes faster than ½Dv̄²), so the three surface elements of the K=6
grid fell back on every stage (ρ_c −0.7% by t=500, accelerating).

**Fixed-point checks** (`test/test_dg3d.jl`): `_limit!(U_eq)` returns U_eq bitwise;
`rhs(U_eq) ≡ 0` exactly; the unseeded star with the limiter ON stays at |ρ_c/ρ_c0 − 1| ≤ 2×10⁻¹³
and |v| ≤ 4×10⁻¹¹ over 250 M⊙ (was −10% eroding); cell means are conserved to 8×10⁻¹⁶ while
the limiter engages on 67/216 elements under a violent perturbation.

**Seeded star, K=6, p=2, v^r = A(r/R)Y₂₀** (`dgcart3d_limiter_census` classifies what the limiter
does to each element; "star fallback" = elements with stellar nodes handed to the flattening
limiter):

| A | t [M⊙] | ρ_c/ρ_c0 − 1 | Σ D/Σ D₀ − 1 | max∣v∣ | elements free / scaled / fallback | star fallback |
|---|---|---|---|---|---|---|
| 10⁻² | 100 | −1.4×10⁻⁴ | +1.6×10⁻⁵ | 0.075 | 209 / 0 / 7 (all r/R = 1.29–1.36) | 0 |
| 10⁻² | 600 | −5.8×10⁻⁴ | −1.6×10⁻⁴ | 0.003 | 211 / 0 / 5 | 0 |
| 10⁻² | 1000 | −7.3×10⁻⁴ | −3.6×10⁻⁴ | 0.003 | 211 / 0 / 5 | 0 |
| 10⁻³ | 600 | −1.3×10⁻⁵ | +1.5×10⁻⁶ | 0.001 | 216 / 0 / 0 | 0 |

The A=10⁻² drift decelerates (−1.2, −1.3, −1.0, −0.9, −0.7, −0.4, −0.2, −0.07 ×10⁻⁴ per 100 M⊙
from t=200) toward ≈ −7.5×10⁻⁴, the order of the mode's kinetic energy fraction thermalised by
the scheme's dissipation as the oscillation decays; it scales as A² (−1.3×10⁻⁵ at A=10⁻³).
Pre-fix values at the same (K, A, t=600): −1.0×10⁻¹; first attempt −3.6×10⁻²; second −8.5×10⁻³.
The fallback elements are the pure-atmosphere corner elements of the octant box, where ejecta
arrive with a negative τ error; no stellar element is ever flattened after the initial half
period (a transient of ≤2 surface elements at t≈50 is recorded in the test).

**f-mode from the fixed engine** (`analyze_qnm`, window [0, T], periodogram / matrix pencil;
1D reference 1.88291 kHz):

| K (nodes/axis) | A | T [M⊙] | periods | df/f | f_pgram [kHz] | f_pencil [kHz] | vs 1D | envelope | ρ_c drift |
|---|---|---|---|---|---|---|---|---|---|
| 6 (18) | 10⁻² | 600 | 5.6 | 18% | 1.8142 | 1.7568 | −3.6% / −6.7% | decaying (0.55) | −5.8×10⁻⁴ |
| 6 (18) | 10⁻³ | 600 | 5.6 | 18% | 1.8228 | 1.7581 | −3.2% / −6.6% | flat at noise floor | −1.3×10⁻⁵ |
| 8 (24) | 10⁻² | 600 | 5.6 | 18% | 1.8104 | 1.8512 | −3.9% / −1.7% | decaying (0.58) | −5.2×10⁻⁴ |
| 10 (30) | 10⁻² | 600 | 5.6 | 18% | 1.6935 | 1.6887 | −10.1% / −10.3% | decaying (0.57) | −3.7×10⁻⁴ |
| 8 (24) | 10⁻² | 1300 | 12.1 | 8.7% | 1.7985 | 1.7708 | −4.5% / −6.0% | decaying (0.58) | +1.3×10⁻⁴ (min −5.3×10⁻⁴ at t≈700, recovering) |

A third estimator, the median half-period between zero crossings of the cubic-detrended q₂
record after t=100, gives 1.70/1.81/1.75/1.71/1.80 kHz for the five rows (−10 to −4%) with a
6–12% half-period scatter; the detrended-away slow quadrupole is comparable to the oscillation
amplitude in every run. Pre-fix: 1.61/1.59/1.57 kHz (−14 to −16%) at K=6/8/10 with a growing
envelope. **Reading.** The fixed engine's f is low by 2–10% with estimator scatter of the same
size and no monotone trend in K (K=10 is the lowest). That residual is the staircase/atmosphere
surface of a nodal DG star at 14–23 nodes across R — the same systematic the linear engine shows
with a masked (rigid-wall) surface, −7% (§7.7 item 1) — and not the limiter, whose background is
now static. It is not a precision measurement of anything; the linear cut-cell engine is.



### 7.9 1D hp-adapted DG star — the Hébert–Kidder–Teukolsky (2018) method in `DGStarHP` (2026-09-17)

`src/dg/DGStarHP.jl` rebuilds Sec. VI.A of Hébert, Kidder & Teukolsky, PRD 98, 044041
(arXiv:1804.02003): symmetric staggered domain (no node at r=0), regions of different element
size and polynomial order (p=3 inside and outside, thin p=1 or p=2 elements across the surface,
which here sits on an element boundary), Valencia GRHD on the frozen TOV metric with the
covariant momentum and a Γ-law EOS, Galeazzi-type atmosphere fixing (ρ_atm = 10⁻¹³ρ_c, entropy
bounds κρ ≤ ε ≤ 100κρ), HLL flux, SSP-RK3, the minmod ΛΠ¹ limiter on the surface elements only.
Their star is the anchor star of this file. Their grids I1/I2/I1R (Table I) are rescaled to
R = 9.586 M⊙; I2R (quadratic surface, refined) is ours. Data: `repro/data/dgstarhp_hkt_*.csv`
from `repro/dgstarhp_hkt.jl`; figure `paper/figs/dgstarhp_hkt.png` (`repro/dgstarhp_hkt_figure.py`);
tests `test/test_dgstarhp.jl` (46).

**Four things that had to be understood before the method worked** (each is a switch in
`setup_dgstarhp` and a test):

1. **The doubled domain has an unphysical antisymmetric mode** (odd δρ, even v — a translation
   of the star) that the spherical problem does not have. It grows from round-off at 0.049/M⊙
   (e-fold 20 M⊙, the dynamical time) on every grid and for every flux, CFL and limiter, is
   first visible at the central element, and destroys the star by t ≈ 600 M⊙; its measured
   antisymmetric fraction is > 1 (pure antisymmetry gives √2). Projecting the state onto the
   physical parity after every stage (`symmetrize=true`) removes it exactly. The paper's
   exponential filter (α = 36, s = 6 on all high-order elements — their central-cube setting)
   also suppresses it, but only by damping the linear mode of every element 5% per step; s ≥ 12
   does not. This is presumably the "numerical instability in S̃_i on O(100 M⊙) timescales" they
   cure with the filter.
2. **A p=3 central element has an O(1) static force error**: the flux x²p(x) of the static star
   has an x⁴ term whose aliased derivative is as large as the (∝x³) gravitational force at the
   innermost nodes — residual 3× gravity, scale-free (same on I1 and I1R); the adjacent element
   has 0.6. With a p=7 central element spanning three nominal widths: 1.7×10⁻⁴ and 1.8×10⁻².
3. **Linear surface elements carry a static residual of 0.5× gravity** (p ∝ (R−r)² is not linear);
   quadratic ones hold it to 2.6×10⁻³. With linear elements and no slope limiter the star launches
   a wind; the paper's minmod flattens S̃ every stage and turns those elements into first-order
   finite volumes whose mean force balance is O(h²) — that is why their scheme holds.
4. **The minmod slope must be the exact linear modal coefficient.** The LGL-quadrature formula
   1.5 Σ w ξ u is exact only for p ≥ 2 and is 3× too large on p=1 elements; with it the limiter
   rewrote every surface slope every stage and the star was destroyed (ρ_c −97%, M_b +25% from
   atmosphere resets of negative nodes). Fixed, the paper's I1 reproduces their Figs. 9 and 12.

**Results** (unseeded = the paper's protocol: the truncation-level settling transient is the
only excitation; seeded = v̂ = 10⁻³ sin(πr/R); linear radial Cowling modes F = 2.686, H1 = 4.550,
H2 = 6.341, H3 = 8.108 kHz):

| grid, limiter, subtraction | seed | T [M⊙] | err[D̃] at T | ρ_c/ρ_c0 − 1 | M_b/M_b0 − 1 | atmosphere | ρ_c spectrum peaks [kHz] |
|---|---|---|---|---|---|---|---|
| I1, minmod, none (the paper's) | — | 10⁴ | 1.1×10⁻³ | −3.3×10⁻⁴ | −3×10⁻¹¹ | transient v≤0.32, then quiet | 2.655 (F −1.2%), 8.19 (H3) |
| I1, minmod, none | 10⁻³ | 4000 | 4.4×10⁻² | −3.2×10⁻³ | +6.8×10⁻⁴ | wind, v→1 | 2.652 |
| I1R, minmod, none | — | 10⁴ | 5.1×10⁻³ (growing) | +1.6×10⁻³ | −9×10⁻¹⁰ | v≤0.19 | 4.585, 6.319, 8.069 (H1–H3 ≤0.8%) |
| I1R, minmod, none | 10⁻³ | 4000 | 7.1×10⁻³ | +1.2×10⁻⁴ | −2×10⁻¹⁰ | v≤0.19 | **2.685, 4.543, 6.298** (F, H1, H2 ≤0.7%) |
| I2, minmod, none | — | 4000 | 3.5×10⁻² | +1.5×10⁻² | +6.4×10⁻⁴ | wind | — |
| I2, wb, none | — | 10⁴ | 1.0×10⁻³ | −3.3×10⁻⁴ | −2×10⁻¹¹ | quiet | — |
| I2R, wb, none | — | 4000 | **1.6×10⁻⁵** | −2.9×10⁻⁶ | −8×10⁻¹² | quiet | 2.692, 4.573 |
| I2, wb, subtraction | — | 4000 | **1.3×10⁻¹²** | −1.4×10⁻¹² | −1×10⁻¹² | quiet | — |
| I2, wb, subtraction | 10⁻³ | 4000 | 9.3×10⁻⁴ | +1.8×10⁻³ | −6×10⁻¹⁰ | v≤0.47 | 2.690 (F +0.1%) |
| I2R, wb, subtraction | 10⁻³ | 4000 | 2.7×10⁻³ | +4.1×10⁻³ | −2×10⁻¹⁰ | v≤0.25 | **2.685, 4.547, 6.326, 8.050** (F–H3 ≤0.7%) |
| I1, wb, none | — | 4000 | 1.2×10⁻¹ | −2.9×10⁻⁴ | +1.9×10⁻² | wind | — |
| I2, mean scaling, subtraction | 10⁻³ | 4000 | 2.7×10⁻² | +1.4×10⁻³ | +2×10⁻⁹ | v≤0.43 | 2.666 |
| I2, no limiter, subtraction | — | 4000 | 2.3×10⁻³ | +3.8×10⁻³ | +8×10⁻⁴ | v≤0.95 | — |

**Reading.** (i) The paper's I1 scheme reproduces their numbers: err[D̃] settles at 1×10⁻³
(theirs 7×10⁻⁴ at 10⁴), ρ_c −3×10⁻⁴ (theirs −5×10⁻⁴), F and the first overtones appear in the
spectrum of the settling transient, and M_b is conserved to 10⁻¹¹ because our minmod acts on the
densitized D̃ (theirs loses M_b at 10⁻⁴). (ii) Our refined I1R does **not** reproduce their
order-of-magnitude improvement (5×10⁻³ and slowly growing, e-fold ≈ 3000 M⊙, against their
3×10⁻⁵): thin linear surface elements plus minmod are marginal in this implementation. (iii) The
equilibrium-preserving scaling limiter of `DGCart3D` on quadratic surface elements reaches that
quality instead: I2R settles to 1.6×10⁻⁵ with a quiet atmosphere, and with the residual
subtraction the I2 star is an exact fixed point (10⁻¹²). (iv) Seeded finite-amplitude
oscillations (v = 10⁻³c) are where the schemes differ most: linear elements with minmod blow a
wind whatever the seed profile, while the quadratic-element wb configuration returns F, H1, H2,
H3 to 0.7% or better on I2R. (v) The mean-based scaling (the paper's MRS pathology, and the
old `DGCart3D` limiter) erodes the surface; with no limiter at all the surface goes slowly
unstable. (vi) The subtraction is exact at the projected equilibrium and is a small fake force
elsewhere (of the size of the raw residual), so it should only be used with a grid whose raw
residual is small (quadratic surface, p=7 centre) — under a finite-amplitude oscillation it
still drifts ρ_c by +2–4×10⁻³ in 4000 M⊙, where the same runs without it drift by −2×10⁻⁴ but
damp the mode 100× faster.

---

*Generated as part of the pre-publication audit. Suite state and all tabulated numbers
correspond to the commit in which this file was added.*
