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
8. **`DGCart3D`'s momentum and energy sources were incomplete — FIXED 2026-09-17** (found while
   deriving `DGBall3D`): it applied only −α(ε+p)W²Φ′ n_i, without the √γ factor, without the
   metric-derivative term (α/2)√γ T^{ab}∂_iγ_ab and without the +√γ p αΦ′ piece of −√γ E ∂_iα, so
   its static residual was not TOV (hidden by the well-balanced subtraction). With the complete
   sources the interior static residual is 1% of gravity (the p=2 derivative error; the analytic
   balance of the coefficient arrays closes to 10⁻¹⁶); the old sources were off by 21% at the
   centre. Effect on the ℓ=2 f-mode (same runs as §7.8, `repro/data/dgcart3d_fmode_scan_srcfix.csv`):

   | K | A | T | periodogram / pencil, old sources | periodogram / pencil, fixed sources |
   |---|---|---|---|---|
   | 6 | 10⁻² | 600 | −3.6% / −6.7% | −0.55% / −0.04% |
   | 8 | 10⁻² | 600 | −3.9% / −1.7% | −1.1% / −1.2% |
   | 10 | 10⁻² | 600 | −10.1% / −10.3% | −3.0% / −3.5% |
   | 8 | 10⁻² | 1300 | −4.5% / −6.0% | −2.9% / −1.3% |

   Two thirds of the "staircase systematic" of §7.8 was the wrong source term. What remains
   (−0.5 to −3.5%, still non-monotone in K) is the staircase surface. test_dg3d.jl: 25/25.
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
**Addendum (same day):** two thirds of that offset was the incomplete source term of §7.7 item 8;
with the complete sources the same scan gives −0.5 / −1.1 / −3.0% at K=6/8/10 (table there).



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


### 7.10 3+1D cubed-sphere DG star — the Hébert–Kidder–Teukolsky method in `DGBall3D` (2026-09-17)

`src/dg/DGBall3D.jl` implements their 3D method (Secs. II–IV, VI.B, App. A–B): a rounded central
cube surrounded by six-wedge cubed-sphere shells that conform to the star (their Fig. 14 /
Table II rescaled to R = 9.586 M⊙), strong-form nodal DG on mapped hexahedra with the chain-rule
flux divergence and the discrete Jacobian, geometric face matching, the complete Cartesian
Cowling sources (which reduce to TOV; `DGCart3D`'s do not — §7.7 item 8 below), Γ-law EOS with
Galeazzi fixing, HLL, SSP-RK3, the `DGCart3D` equilibrium-preserving scaling limiter on the surface
shells, their exponential momentum filter, and an optional static-residual subtraction. The
surface sits on a shell boundary and the surface shells are quadratic (the I2 lesson of §7.9).
Coarse grid `nt=2, p=3`: 464 elements, 27 776 nodes, Δt = 0.10 M⊙, 0.03 s/step on 8 threads.
Tests `test/test_dgball3d.jl` (23). Data `repro/data/dgball3d_hkt_*` from `repro/dgball3d_hkt.jl`;
figure `paper/figs/dgball3d_hkt.png`.

**Geometry and consistency.** All Jacobians positive; the ball's volume to 2.3×10⁻⁴; all 42 240
interior face nodes matched to their partner at identical positions (the wedge rotations and the
cube–wedge interfaces need no bookkeeping); a uniform moving gas on a flat metric has RHS
6×10⁻¹⁸ (exact free stream); the origin is a regular node. Static momentum residual / gravity:
1.9–2.7% in the cube, centre and interior shells, 4.7% in the quadratic surface shells (curved
p=3 elements; the 1D I2 grid had 10⁻³).

**The aliasing instability and the filter (load-bearing).** Every seeded run without the filter
develops an exponential instability with e-folding time ≈ 250 M⊙ once the seeded transient has
decayed (ℓ=0 and ℓ=2 alike; HLL or LLF; with or without the entropy floor or the subtraction):
err[D̃] and the central density drift at an accelerating rate, mass is created (the chain-rule
strong form is not exactly conservative on curved elements: ΔM_b/M_b ∝ seed amplitude), and the
star blows up at t ≈ 600 (no subtraction) to 1900 M⊙ (with). This is the paper's "numerical
instability in S̃_i on O(100 M⊙) timescales" and their cure works verbatim: the exponential
modal filter exp[−α(i/p)^s] applied after every step to the **momentum only**, α=36, s=6 in the
cube and centre shells and s=12 in the cubed-sphere shells. With it, err[D̃] keeps falling (ℓ=2
seed, t∈[600,800]: 6.0×10⁻⁴ growing → 1.1×10⁻⁵ falling), ρ_c is flat to 10⁻⁶ and M_b to 10⁻⁸.
Filtering all five variables destroys the star. The filter is on by default and is not optional.

**Results** (nt=2; seeds v = 10⁻³ sin(πr/R) n̂ (ℓ=0) and 10⁻³ sin(πr/R) Y₂₀ n̂ (ℓ=2); 1D references
F = 2.686 kHz, f = 1.88291 kHz):

| run | T [M⊙] | err[D̃] at T | ρ_c/ρ_c0 − 1 | M_b/M_b0 − 1 | frequency |
|---|---|---|---|---|---|
| static, subtraction | 200 | 5×10⁻¹⁶ | 0 | 0 | — |
| static, no subtraction (paper's scheme) | 2000 | 5.1×10⁻⁴ | −1.3×10⁻⁴ | +1.6×10⁻⁵ | — (paper B1: 6×10⁻⁴, +2.5×10⁻⁴ at 4000) |
| ℓ=0 seed, subtraction | 2000 | 7.8×10⁻⁶ | −6.9×10⁻⁶ | −4×10⁻⁸ | F = 2.674 kHz (−0.45%); 2F harmonic at 5.3 |
| ℓ=2 seed, subtraction | 2000 | 1.1×10⁻⁵ | −1.7×10⁻⁶ | −2×10⁻⁸ | f = 1.845 / 1.888 / 1.912 kHz (−2.0 / +0.2 / +1.5%) |
| ℓ=2 seed, no subtraction | 2000 | 5.1×10⁻⁴ | −1.3×10⁻⁴ | +1.6×10⁻⁵ | f = 1.840 / 1.884 / 1.916 kHz (−2.3 / +0.1 / +1.8%) |

The three ℓ=2 numbers are the periodogram peak, the matrix-pencil pole and a damped-sinusoid fit
over [30, 600] M⊙; they scatter by ±2% because the mode is heavily damped on this grid (fit
Q ≈ 3.3, e-folding 115 M⊙ ≈ one period): the record carries only 5–8 usable periods and the
s=12 filter removes 24% of every element's quadratic mode per step, which on a two-element-per-
quadrant grid is a large part of an ℓ=2 pattern. The result is therefore "f within ±2% of the
1D Cowling value" — the same conclusion the linear cut-cell engine reaches at 0.1%, and a
quantitative improvement over `DGCart3D`'s −2 to −10% staircase systematic (§7.8) on a grid with
half the nodes. **Filter strength is not a free knob.** Weakening the shell filter to s = 16 (5% per step on
the quadratic momentum mode instead of 24%) does not merely reduce the dissipation: the ℓ=2 run
then develops a *faster* instability (err[D̃] 7×10⁻⁵ at t=180 → 1.2×10⁻² at t=500, e-fold ≈ 50 M⊙,
blow-up at 540 M⊙) — faster than with no filter at all (e-fold 250 M⊙). Removing the top
momentum mode while leaving the quadratic mode nearly undamped is destabilizing; the paper's
s = 12 works because its 24% quadratic-mode damping outweighs the growth. Their parameters are
kept as the defaults and documented as such.
**Resolution (nt=3, 1053 elements, 63 072 nodes, same radial shells).** With the paper's filter
and the subtraction, the ℓ=2 seed decays as at nt=2 until t ≈ 450 M⊙, and then a *spurious
quadrupolar deformation* grows (q₂₀/q₀₀ 4.8×10⁻³ → 5.3×10⁻² between 450 and 1200, e-fold ≈ 310
M⊙) and saturates at q₂₀/q₀₀ ≈ 7.5% with err[D̃] 1.0×10⁻², ρ_c −1.2×10⁻³, M_b −3×10⁻⁵; the atmosphere
stays at rest. At nt=2 the unfiltered instability was not quadrupolar (q₂₀ stayed at 3×10⁻⁴
while ρ_c collapsed), so these are grid modes of the cubic-symmetric discretization, not one
physical mode. The static momentum residual of 2–5% of gravity identifies their driver: it is
the geometric error of representing a 45°-wide spherical patch by a cubic (the discrete
Jacobian; ~(Δθ)⁴/4! ≈ 1.6%), far above the 10⁻³ of the 1D grid, and the same error acts on every
perturbation as a spurious force of that size with the grid's angular pattern. The paper's B1 grid
has 15°-wide elements (6×6 per wedge, (Δθ)⁴/4! ≈ 2×10⁻⁴) and shows no such growth over 10⁴ M⊙.
**What refinement does and does not do** (ℓ=2 seed, subtraction unless noted; three
frequency estimators: periodogram / matrix pencil; growth windows of err[D̃] and q₂₀/q₀₀):

| grid | nodes | static residual int./surf. | err[D̃], q₂₀/q₀₀ at t∈[600,800] | at the end | f [kHz] |
|---|---|---|---|---|---|
| nt=2, p_t=3 (filter) | 27 776 | 2.7% / 4.7% | 1.1×10⁻⁵, 1.8×10⁻⁴ (falling) | 1.1×10⁻⁵ at 2000 | 1.845 / 1.888 |
| nt=2, p_t=3, no filter | 27 776 | 2.7% / 4.7% | 6.0×10⁻⁴, 3×10⁻⁴ (growing) | blow-up at 1900 | — |
| nt=3, p_t=3 (filter) | 63 072 | 2.7% / 4.7% | 3.7×10⁻⁴, 7.9×10⁻³ (growing) | 1.0×10⁻², 7.5% at 2000 (saturated) | — |
| nt=3, p_t=3, no subtraction | 63 072 | 2.7% / 4.7% | 5.0×10⁻⁴, 8.4×10⁻³ (growing) | 1.3×10⁻³, 2.0% at 1000 | 1.825 / 1.855 |
| nt=2, p_t=5 (filter) | 63 072 | 0.43% / 3.8% | 4.6×10⁻⁵, 6.9×10⁻⁴ (falling) | 2.1×10⁻³, 2.7% at 1500 (growing from ~800) | 1.891 / 1.917 |
| nt=3, p_t=3, shell filter s=6 | 63 072 | 2.7% / 4.7% | 2.7×10⁻⁶, 4.9×10⁻⁵ (falling) | 2.4×10⁻⁶, 1.0×10⁻⁵ at 1200 (stable) | 1.707 / 1.809 (over-damped) |

Raising the tangential order to 5 cuts the interior geometric residual six-fold (2.7% → 0.43%,
as (Δθ)^{p+1} predicts) and brings the periodogram f to +0.4%, but the quadrupolar grid mode
still appears, later (t ≈ 800) and more slowly (e-fold ≈ 190 M⊙); removing the subtraction
changes nothing. The paper's central-cube strength (s = 6) applied to every shell does hold
nt=3 stable for 1200 M⊙ (err[D̃] 2×10⁻⁶ and falling, M_b to 3×10⁻⁹), but at that strength the
mode is over-damped and its frequency estimators scatter to −4…−9%: the filter that controls
the grid mode is the same dissipation that destroys the measurement. **Verdict on the 3D method
as implemented:** on the coarsest grid, with the
paper's filter, it is stable for 2000 M⊙ and reproduces their static results and the 1D Cowling
F and f to 0.5% and ±2%; refinement in either h or p exposes a slowly growing, saturating
quadrupolar grid mode that the exponential filter does not control. The likely cure is the
one the DG literature uses for exactly this symptom on curved elements — a split-form
(entropy-stable) or over-integrated volume term in place of the chain-rule strong form — and
that is the next engineering step, not a parameter choice.

**Split form and localization (2026-09-18).** Two follow-ups were done. (i) The volume term was
re-implemented in the flux-differencing split form of Gassner–Winters–Kopriva with Kopriva's
curl-form metric vectors and a Kennedy–Gruber-type two-point flux for the Valencia system
(`volume=:split`): the discrete metric identity closes to 4×10⁻¹⁵, free stream is exact
(2×10⁻¹⁷), the static star is a fixed point with the subtraction and settles to 1.5×10⁻³ without
it. It did **not** cure the nt=3 grid mode: its static residual is larger (7–11% of gravity
against 2–5%) and the nt=3 seeded run blows up before 1500 M⊙ where the chain form saturates,
at twice the cost. The chain form stays the default. (ii) The Y₂₀ moment and the density error
were tracked *per region* on the chain-form nt=3 run (`repro/data/dgball3d_gridmode_localization_nt3.txt`):
the growth starts in the **surface shells** (their q₂₀/q₀₀ −1.8×10⁻³ → −5.7×10⁻³ and their
relative density error 3.0×10⁻³ → 1.1×10⁻² between t=500 and 800) and spreads inward with a
lag (interior error 2.3×10⁻⁴ → 1.0×10⁻³), while the cube and centre shells stay at 10⁻⁵. The
grid mode is therefore a surface-shell instability of the star–atmosphere interface on the
curved shells, not an interior volume-discretization effect — which is also why the strong s=6
filter in every shell (it damps the surface shells' momentum modes) is what stabilizes nt=3,
and why the split form, which changes only the volume term, cannot help. The paper controls
exactly this region with the minmod slope limiter on *linear* surface shells; our 3D surface
shells are quadratic with the positivity-only scaling limiter, the configuration the 1D study
favoured. The corrective step is therefore the paper's own: linear surface shells with the
per-direction ΛΠ¹ minmod and the physical-state check, in 3D.

**The paper's surface recipe in 3D (2026-09-18).** `limiter=:minmod` with `surface=:linear`
implements it: 10 linear shells of h=R/32 across [R−2h, R+8h], the ΛΠ¹ slope limiter per
reference direction with neighbour means across the six sides (element adjacency from the
geometric face matching), higher modes dropped when any direction is limited, then the
physical-state check with slope halving; means and slopes in the reference coordinates, as
they do. Static residual of the linear shells 0.59 of gravity (the 1D value). At nt=2 the
unseeded star settles to err[D̃] 1.3×10⁻³ with M_b drifting at 10⁻⁴ (their reported level for
this limiter on deformed elements). At nt=3 with the ℓ=2 seed and the momentum filter, no
subtraction — their B1 recipe except for the angular resolution — err[D̃] is already 4.9×10⁻³
at t∈[200,400], 2.9×10⁻² at [600,800], and the run blows up before 1500 M⊙: worse than the
quadratic shells with the scaling limiter, which saturate. The one remaining difference from
their B1 grid was the element width, 30° here against their 15° (6×6 per wedge), so the same
recipe was run on their grid: nt=6, 5400 elements, 276 480 nodes (their B1 has 5184), Δt = 0.08.
The interior static residual falls to 0.4% as the geometry predicts, the surface shells stay at
0.59, and the seeded run fails the same way — err[D̃] 4.7×10⁻³ at [200,400], 3.5×10⁻² at [600,800],
blow-up before 1200 M⊙. **Angular resolution is not the missing ingredient.** Every 3D surface
treatment tried here is unstable to a seeded non-radial perturbation beyond a few hundred
dynamical times except two: quadratic shells with the equilibrium-preserving scaling limiter and
the paper's filter at nt=2 (stable, f to ±2%), and the s=6 filter in every shell at nt=3 (stable,
over-damped). What differs from the paper's stable B1 evolution that we have not reproduced:
their isotropic coordinates (a conformally flat spatial metric; ours is the anisotropic areal
Cartesian one), their exact atmosphere and inversion-fixing details (Galeazzi et al. App. C),
their surface inside a shell rather than on a boundary, a smaller Courant number
(Δt·a_max/Δx_min ≈ 0.18 against our 0.25), and whatever their ΛΠ^N implementation does beyond
the description — their code (SpECTRE) is open and its `Minmod` limiter and TOV test are the
right reference to check against next. Data: `repro/data/dgball3d_hkt_series_mm_*.csv`.

**Step-2 verdict.** The 3D method is implemented and verified in its static and coarse-grid
behaviour; at coarse resolution it reproduces the paper's static settling and gives F to 0.5%
and f to ±2%; it is not yet a refinable precision tool, and the obstacle has been localized to
the star–atmosphere interface on the curved surface shells, not to the volume discretization.

### 7.11 The 3+1D Cowling validation against SpECTRE (2026-09-18)

SpECTRE (sxs-collaboration/spectre, `develop` at 00f341c) is the open successor of the code
used in the paper; it could not be run here (no container runtime or CMake on this machine, and
3.7 GB of free disk against a build that needs tens), so the comparison is against its
regression setup, its published results and its source.

**What SpECTRE's own Cowling TOV test is.** `tests/InputFiles/GrMhd/ValenciaDivClean/TovStar.yaml`
evolves our star (K=100, Γ=2, ρ_c=1.28×10⁻³, isotropic coordinates) on a filled sphere of
inner radius 3 and **outer radius 7 — inside the star** (isotropic R = 8.125), with the analytic
TOV solution as Dirichlet boundary data: "these domains are chosen so that we only need to
simulate the interior of the star, avoiding the surface discontinuity". Strong-inertial
(chain-rule) DG formulation, Gauss–Lobatto, 5 points per element, RK3 SSP with a CFL safety
factor of 0.5, no filter, DG–FD subcell fallback armed but idle on a smooth interior, an
atmosphere at ρ_atm = 10⁻¹⁵ with cutoff 1.1×10⁻¹⁵, Kastaun inversion, and velocity and entropy
("kappa") limiting near the atmosphere; observed: L2 errors of ρ, ε, v, p against the analytic
solution. The CI run lasts three slabs. The full-star Cowling evolutions of Deppe et al. 2021
(arXiv:2109.12033, Sec. IV.7: [−20,20]³, 6/12/24 P5 elements, ρ_atm 10⁻¹⁵, cutoff 1.01×10⁻¹⁵)
use the DG–finite-difference hybrid; of the classical DG limiters they tried on that star,
ΛΠ^N "falls back to the linear approximation", Krivodonova "succeeded at some resolutions (3 of
16 attempted runs)", HWENO and simple WENO were stable only at P2, and "the only limiting
strategy we can endorse is a discontinuous Galerkin–finite-difference hybrid method". Their
GH+MHD full-star inputs run with `AlwaysUseSubcells: true` — finite differences everywhere —
and a `Hypercube` exponential filter with half-power 64, i.e. c_i → c_i exp[−36 (i/N)^{128}]:
the top mode only. The 2018 paper's s = 6/12 filters are far stronger. Their `Minmod` acts in
element-logical coordinates, is non-conservative on deformed elements (documented), and divides
the neighbour-mean differences by ½(1 + h_neighbour/h_self); that factor is now in ours too.

**Their interior-only setup in our engine** (`ball_grid(…; interior_only=true, rmax_int_fac=0.85)`,
`boundary=:tov`, no subtraction, no limiter; nt=2, 152 elements, degree p in every direction;
errors are the point-wise L2 norms SpECTRE observes; `repro/data/dgball3d_spectre_interior*.{csv,txt}`):

| p | filter | err(ρ) t=100 | t=300 | t=500 | ρ_c/ρ_c0 − 1 at 500 |
|---|---|---|---|---|---|
| 2 | none | 4.6×10⁻² | blow-up | | |
| 3 | none | 1.5×10⁻³ | 4.8×10⁻³ | 1.6×10⁻² | −1.2×10⁻² (e-fold ≈ 170 M⊙) |
| 4 | none | 8.8×10⁻⁴ | 4.0×10⁻³ | blow-up | |
| 3 | none, CFL halved | identical to the unfiltered p=3 row | | | |
| 3 | none, LLF instead of HLL | identical to the unfiltered p=3 row | | | |
| 3 | none, split form | 5.3×10⁻² | 1.0×10⁻¹ | 1.1×10⁻¹ | −0.26 |
| 3 | none, residual subtraction | 6×10⁻¹⁵ | 2×10⁻¹⁴ | 4×10⁻¹⁴ | fixed point (nothing seeds it) |
| 2 | momentum, s=6/12 | 1.2×10⁻² | 2.3×10⁻² | 3.6×10⁻² | −2.2×10⁻² (still growing) |
| 3 | momentum, s=6/12 | 2.3×10⁻⁴ | 2.7×10⁻⁴ | 3.2×10⁻⁴ | −1.5×10⁻⁴ |
| 4 | momentum, s=6/12 | 3.8×10⁻⁵ | 4.2×10⁻⁵ | 5.1×10⁻⁵ | −8.4×10⁻⁶ |
| 3 | all variables, s=6/12 | 1.7×10⁻⁴ | 1.9×10⁻⁴ | 2.0×10⁻⁴ | −1.3×10⁻⁴ |
| 3 | momentum, with subtraction | 6×10⁻¹⁶ | 6×10⁻¹⁶ | 6×10⁻¹⁶ | 0 |

**Reading.** With no surface, no atmosphere and no limiter in play, the unfiltered chain-rule
scheme is still unstable on the smooth star interior, with the same e-folding time as the
full-star runs (≈ 170–250 M⊙), independent of the Courant number and of the Riemann solver. The
instability of §7.10 is therefore a property of the volume scheme on the curved elements — the
paper's "numerical instability in S̃_i caused by aliasing", now isolated — and the surface
shells were where it showed first, not where it comes from. The paper's momentum filter is the
cure here too: with it the interior error is 3.2×10⁻⁴ at p=3 and 5.1×10⁻⁵ at p=4 after 500 M⊙,
a factor six per degree, which is the spectral convergence SpECTRE's interior test is built to
show; at p=2 the same filter is too weak to hold the scheme. The split-form volume term makes
the unfiltered instability worse, not better. SpECTRE's three-slab regression run cannot see a
170 M⊙ instability; its production configuration removes the top mode every step and hands
troubled elements to finite differences, which is the direction any further work on the
nonlinear 3D star should take.

**The paper's Sec. VI.B protocol, literally** (`repro/data/dgball3d_hkt_series_B1_unseeded_nt6_p3.csv`,
figure `paper/figs/dgball3d_b1.png`, analysis `repro/dgball3d_b1_analysis.py`): their B1 grid
structure (nt=6, 5400 elements, 276 480 nodes, Δt = 0.08), linear surface shells with the
per-direction minmod, the momentum filter with s=6 in the centre and s=12 in the shells, no
subtraction, **unseeded** — the settling transient is the only excitation, as in their Figs. 15–16.
Ours does not settle: err[D̃] rises from 9×10⁻⁴ (t∈[0,200]) through 2.1×10⁻³, 5.7×10⁻³ and
1.5×10⁻² to 2.2×10⁻² at [800,1000] while the ρ_c oscillation amplitude grows 2.5×10⁻³ → 6×10⁻³ →
1.8×10⁻² → 4.4×10⁻² → 5.7×10⁻² (e-fold ≈ 200 M⊙), the atmosphere reaches 0.9c after t ≈ 800, M_b
drifts +2.7×10⁻³, and the run blows up at t ≈ 950; the pre-blow-up spectrum is dominated by a
4.38 kHz line (H₁ −3.7%) with F at 2.694 (+0.3%) weak. Their B1: err[D̃] ≈ 6×10⁻⁴, ρ_c settling at
+2.5×10⁻⁴, M_b error 10⁻⁴, stable to 10⁴ M⊙, F and H₁ sharp. **The validation of Sec. VI.B is not
reproduced.** Together with §7.11's interior result (the filtered interior converges) and
§7.10's localization, the failure sits at the star–atmosphere interface: linear radial shells on
which the s=6/12 momentum filter annihilates the radial slope of the momentum deviation every
step and the minmod acts on reference-coordinate means, driven from a surface whose density
cut sits nine orders of magnitude below the last stellar node. The ingredients of their setup
not reproduced here are the isotropic coordinates, the surface inside a shell rather than on a
boundary, the exact inversion-fixing recipe, and a Courant number of 0.18 against our 0.25;
none of these changed the instability in the variants tried (§7.10). SpECTRE, the successor
code, no longer evolves a DG star surface at all: its Cowling test excludes the surface and its
full-star runs use finite-difference subcells. That is the state of the art and it is where a
nonlinear 3D star for this project would have to go.

### 7.12 The DG / finite-difference hybrid — `DGSubcell` + `DGStarFD`, radial star (2026-09-19)

The scheme Deppe et al. (PRD 105, 123031; arXiv:2109.12033) endorse after finding that every
classical DG limiter fails on a neutron star: each element carries both a DG polynomial on N
Legendre–Gauss–Lobatto nodes and M = 2N−1 finite-volume subcells, a troubled-cell indicator
chooses which one is evolved, and the two are exchanged conservatively. `src/dg/DGSubcell.jl`
holds the grid-agnostic machinery (projection, reconstruction, indicators, MC slope),
`src/dg/DGStarFD.jl` the radial star. Tests `test/test_dgstarfd.jl` (62). Data
`repro/data/dgstarfd_hybrid_*`; figure `paper/figs/dgstarfd_hybrid.png`.

**The machinery, verified to machine precision.** The projection P (DG → subcell averages) is
built from exact Gauss–Legendre integrals of the Lagrange basis over each subcell; the
reconstruction R (subcell → DG) is the constrained least squares of SpECTRE's
`reconstruction_matrix`, minimising ‖Pu − v‖² subject to the element integral being preserved.
At every degree tested (p = 1…7): ‖R P − I‖ ≤ 9×10⁻¹⁶, a constant projects to itself to 10⁻¹⁶,
and P carries the element integral exactly. The Persson indicator gives a top-mode fraction of
1.3×10⁻³ for smooth data against 1.5×10⁻¹ for a jump, with threshold p⁻⁴ = 1.2×10⁻².

**Three implementation points that had to be right.**
1. *A uniform atmosphere element is smooth.* The first indicator forced every element below the
   density cutoff onto the subgrid — 56% of the grid — because they are "unphysical". SpECTRE's
   `TciOptions` does the opposite: an element lying entirely in the atmosphere is uniform and
   stays on DG; only elements that *straddle* the cutoff are troubled. With that correction the
   subgrid holds 4–10% of the elements, exactly those containing the stellar surface.
2. *The subcell recovery must divide by the cell average of √γ.* M = 2N−1 is odd, so the central
   element's middle subcell is centred on x = 0, where √γ = e^{λ/2}x² vanishes identically: the
   densitized state carries no information there and the recovery returned atmosphere, emptying
   the star from the centre (ρ_c → ρ_atm in the pure-FV run). The cell average is both regular
   and the consistent finite-volume reading of q̄ = (1/h)∫√γρW dx.
3. *Closure capture in the setup.* A nested geometry closure assigning to `α` wrote into the
   enclosing array of the same name — Julia binds an inner assignment to an existing enclosing
   local. Renaming the closure's locals fixed it.

**Results on SpECTRE's configuration** (uniform degree-5 elements, the surface inside an element
rather than on a boundary, no limiter and no filter; K elements per side over [0, 3R]; unseeded
runs to 10⁴ M⊙, seeded runs v = 10⁻³ sin(πr/R) to 4×10³ M⊙; linear Cowling modes F = 2.6861,
H₁ = 4.5495, H₂ = 6.3414, H₃ = 8.1082 kHz):

| K | DG nodes | subcells | subgrid | err[D̃] at 10⁴ | F | H₁ | H₂ | H₃ |
|---|---|---|---|---|---|---|---|---|
| 10 | 126 | 231 | 10% | 3.14×10⁻² | −0.75% | −3.53% | — | — |
| 14 | 174 | 319 | 7% | 1.79×10⁻² | −0.38% | −1.84% | −4.77% | −5.06% |
| 20 | 246 | 451 | 5% | 1.08×10⁻² | −0.12% | −0.91% | −2.94% | −5.35% |
| 28 | 342 | 627 | 4% | 6.86×10⁻³ | **−0.04%** | **−0.36%** | **−1.17%** | −2.62% |
| 20, pure finite volume | — | 451 | 100% | 2.10×10⁻² | −0.56% | −2.45% | −4.96% | — |

err[D̃] settles (it does not grow at any resolution — the left panel of the figure is flat from
t ≈ 2000 to 10⁴) and converges as h^1.4–1.7, the rate being set by the second-order subcells at
the surface rather than by the degree-5 interior. The mode frequencies converge much faster and
reach the linear values to 0.04% (F) and 0.36% (H₁) at K = 28. Baryon mass is conserved to
10⁻¹² — the common numerical flux is shared across every DG/FD interface, so the hybrid is
conservative by construction, and the pure-FV limit is conservative to 2×10⁻¹³.

**What the DG elements buy.** At K = 20 the hybrid and the pure finite-volume run use the same
451 subcells and cost within a factor 1.6; the hybrid's err[D̃] is 1.9× smaller and its mode
frequencies are 3–5× closer to linear theory (F −0.12% against −0.56%, H₁ −0.91% against
−2.45%). Evolving 95% of the star spectrally is worth that.

**Robustness — the point of the method.** Nothing tried blew up. In particular the I1 grid with
the v = 10⁻³ seed, which under the paper's own ΛΠ¹ minmod blows a wind with the atmosphere
reaching |v| → 1 and M_b growing by 7×10⁻⁴ (§7.9), runs with the atmosphere bounded and M_b to
10⁻⁴ under the hybrid. That is the claim of arXiv:2109.12033 reproduced.

**Where the hybrid is NOT the best choice.** On the hp grids designed for limiters it is worse
than the pure-DG scheme of §7.9. The thin surface elements of I1/I2 sit next to interior
elements several times their size, and the subcell ghost exchange across that jump rings the
star: I2 seeded settles at err[D̃] = 2.5×10⁻³ but keeps a 1.3% central-density oscillation whose
spectrum is not the radial tower, and I1 still reaches |v| ≈ 0.9 in the atmosphere. The
equilibrium-preserving scaling limiter of §7.8–7.9 on the I2 grid remains the most accurate 1D
configuration (an exact fixed point with the residual subtraction, err[D̃] = 10⁻³ without it),
because it is built to preserve the TOV equilibrium, whereas the hybrid deliberately lets the
star settle. The two answer different questions: the limiter is for precision on a star that
stays near equilibrium, the hybrid is for robustness when it does not.

### 7.15 `DynGR1D` was not well-balanced — diagnosis and fix (2026-09-22)

Found while setting up the Phase-1 cross-code collapse comparison against AthenaK
(`progress/plan_athenak_bdnk_collapse.md`, `code/athenak-bdnk/`). `setup_dyngr` stores the raw
RHS of the initial TOV state in `Seq_*` and subtracts it at every step (`wellbalanced=true`,
the default). That made the projected equilibrium an exact fixed point — and hid the fact that
the operator underneath it was not balanced at all.

**The diagnosis.** The raw static momentum residual (|∂_t S̃| over the stellar interior divided
by the local gravity scale √γ̃(ρhW²−p)α ∂_r lnα) was **resolution-independent**: 5.347×10⁻² at
Δr = 0.048, 5.346×10⁻² at 0.024 and 5.346×10⁻² at 0.012 for the ρ_c = 1.28×10⁻³ star, and
1.578×10⁻¹ at every resolution for ρ_c = 7.993×10⁻³ (worst cells 22% and 41%). A truncation
error falls as Δr²; this did not fall at all, so it was a consistency error. Subtracting the two
terms the source was missing accounted for **99.99%** of it:

* the momentum source read `−αXr²(ε+p)W²Φ' + 2αXr p`, i.e. it used the energy density
  E = (ε+p)W² where the Valencia source calls for ρhW² − p, dropping **αXr²pΦ'**;
* and it omitted the radial-stress term of (α/2)√γT^{jk}∂_rγ_jk, **αr²X′(ρhW²v²+p)**, entirely.

Both are of order p/ε — which is exactly the size measured (p/ε = 0.44 at the centre of the
unstable star against a 41% worst-cell residual). Separately, the D and τ rows carried the
momentum row's area weight αXr² instead of αr²: the coordinate flux is √γ α F(v^r) with
v^r = v/X, so only the momentum row keeps the X.

**The fix** (`_raw_rhs!`). Every pressure piece of the Valencia source is exactly
p ∂_r(αXr²) — since αr²X′ + 2αXr + αXr²Φ' ≡ ∂_r(αXr²) — so the source is now grouped as
`p·A′ + αr²X′ρhW²v² − αXr²ρhW²Φ'` with **A′ discretised as (A_{k+1}−A_k)/Δr from the same face
areas the flux difference uses**. The pressure terms then cancel analytically against the flux,
−∂_r(Ap) + pA′ = −Ap′, leaving −A[p′ + (ε+p)Φ'] — the TOV equation — so hydrostatic balance
holds to the order of the scheme by construction. The D and τ rows now carry αr²; the evolved
orthonormal momentum S̃ = √γS_r/X carries the 1/X and the frame term −S̃∂_t lnX, with
∂_t lnX = −αK^r_r fixed algebraically by the momentum constraint K^r_r = 4πrXρhW²v; X′ comes
from the Hamiltonian constraint, X′ = X³(4πrE − m/r²). The energy row gains the
α√γT^{ij}K_ij term. `wellbalanced` is kept as an option.

**Verification.**

| | before | after |
|---|---|---|
| raw static residual, Δr = 0.048 / 0.024 / 0.012 / 0.006 | 5.35×10⁻² at every Δr | 6.15×10⁻⁵ → 1.54×10⁻⁵ → 3.87×10⁻⁶ → 9.70×10⁻⁷, **order 2.00** |
| same, unstable star | 1.578×10⁻¹ at every Δr | 8.81×10⁻⁵ → … → 1.51×10⁻⁶, order 1.9–2.0 |
| unkicked star, no subtraction, 200 M⊙ (N = 400/800/1600) | ρ_c ×1.867 | max \|ρ_c/ρ_c0−1\| = 4.3×10⁻³ / 2.6×10⁻³ / 1.5×10⁻³ (sampled every 1 M⊙) |
| dynamical radial mode vs the Chandrasekhar eigenvalue 2.1236 kHz | 2.02 kHz, **−4.9%** | **2.1304 kHz, +0.32%** |
| frozen-metric mode vs the Cowling eigenvalue 4.0099 kHz | 4.6 kHz, +15% | **3.9994 kHz, −0.26%** |
| unstable star −1% kick, `wellbalanced` true vs false | t_AH 26.5 vs 13.5 (factor 2) | ρ_c,end 1.614 vs 1.613×10⁻², max 2m/r 0.9595 vs 0.9598 |
| unstable star +5% outward kick | dispersed to the floor (M_b −24%) | **migrates**, ρ_c oscillating about 1.40×10⁻³ — the stable-branch star of this baryon mass |

The engine now reproduces **both** independent frequency-domain eigensolvers to better than
0.35%, which is the real check on the corrected flux weights. `test_dyngr.jl`: 50/50.

**A recorded result that the fix overturns.** The collapse test used ρ_c = 2.4162×10⁻³
(ε_c = 0.003), described as "moderately compact"; that star is on the **stable** branch (the
maximum mass of this EOS is at ρ_c ≈ 3.16×10⁻³), and a −3% kick cannot unbind it. It
"collapsed" only because the 5%-of-gravity imbalance acted as a steady inward force. With the
corrected operator it oscillates and settles 9% above its initial central density (α_c = 0.55,
max 2m/r = 0.47) — the physical answer. **The former "collapse-to-BH (α_c → 0.001, 2m/r → 0.95,
ρ_c × 4.8)" was an artifact and is withdrawn.** The genuine collapse test is now the
unstable-branch star of Font et al. 2002 (ρ_c = 7.993×10⁻³, M = 1.448, R = 5.838) with a −1%
kick: it collapses with t_AH = 53.7, central proper time τ_c = 10.0, horizon mass M_AH = 1.273
(88% of the gravitational mass), identically with and without the subtraction. Two controls are
now asserted alongside it: the stable star takes the same kick without collapsing, and the same
unstable star kicked **outward** migrates instead — the sign of the kick decides, which is the
physics.

**What still needs the subtraction.** *(Superseded by §7.16 — hydrostatic reconstruction
makes the raw operator exact in the stellar interior, and the subtraction is off by default
from 2026-09-22.)* Balance at this point was to truncation order, not machine
precision: the raw operator leaves ρ_c ringing at a few 10⁻³ over 200 M⊙ (and that amplitude
converges only slowly, because it is set by the mismatch between the interpolated TOV data and
the discrete equilibrium rather than by the residual alone) where `wellbalanced=true` is exact.
The task's 10⁻⁶ target for the unsubtracted operator is therefore **not** met. Exact balance would need hydrostatic reconstruction (reconstructing
the deviation from the local hydrostatic profile). The subtraction is now safe to use, because
what it stores is a truncation-size residual rather than a 5–16% fake force; both settings agree
on every Phase-1 verdict and on the collapse numbers above. The magnetised (toroidal) sector was
transformed consistently with the fluid rows but **not** re-derived from the GRMHD source terms;
its tests (flux conservation, field amplification, ADM-mass back-reaction) still pass and it
remains the least-validated part of this engine.

**Follow-up owed.** Four reproduction scripts and two figures consume this engine and were not
re-run here: `repro/radial_overtone_crosscheck.jl`, `repro/radial_shoot_crosscheck.jl`,
`repro/viscous_damping.jl`, `repro/viscous_fmode_scaling.jl`, `viz/dyngr_1d.jl` and
`viz/dyngr_radial_validation.jl`. Their recorded numbers predate the fix; the frequencies in
particular will move by the same ~5% the fundamental did. `test_radial_spectrum.jl` is
unaffected (it runs no time-domain engine).

### 7.16 `DynGR1D` hydrostatic reconstruction — the raw operator now holds the star exactly (2026-09-22)

§7.15 fixed the flux/source pair and left one thing open: balance held only to truncation
order, so the unsubtracted operator still let ρ_c drift at a few 10⁻³ and the 10⁻⁶ target was
**not** met. The remedy named there was hydrostatic reconstruction. This is it.

**The construction.** A static barotropic star obeys an exact first integral. With
p′ = −(ε+p)(ln α)′ and dp/(ε+p) = dH for the pseudo-enthalpy H(p) = ∫₀^p dp′/(ε+p′) = ln h,

>  **q ≡ H(p) + ln α = const**    (the relativistic Bernoulli integral)

so *q*, not (ρ,p), is the variable the scheme should limit.

1. **Faces.** Build q at cell centres, limit it with the same minmod slope, and map each face
   state back onto the local hydrostatic profile, p_face = H⁻¹(q_face − ln α_face); ρ and ε then
   follow from p through the EOS, so the face state is barotropically **consistent** (limiting ρ
   and p independently is not). A star in equilibrium has q ≡ const, so the slope is exactly
   zero, the two face states are the same number, and the HLL dissipation vanishes identically
   instead of being O(Δr²).
2. **Source.** Gravity is discretised as the cell's OWN hydrostatic profile differenced across
   its two faces, `[A_R(p_eq,R − p) − A_L(p_eq,L − p)]/Δr`, in place of the algebraic
   −A(ε+p)Φ′. It is second-order consistent (the two ±Δr/2 extrapolations are centred) and,
   because those are the very p_eq the flux reconstructed, it cancels the flux difference
   **identically**, not to truncation order. What is left of the Φ′ term is the kinetic piece
   (ε+p)(W²−1), written so that it is identically zero at v = 0. The p·∂_r(αXr²) grouping of
   §7.15 is unchanged, and the D and τ rows need nothing: with v = 0 and equal face states their
   fluxes and sources are already exactly zero.
3. **Initial data.** A well-balanced operator can only hold an equilibrium *of the
   discretisation*. Sampling the continuum TOV solution at cell centres is not one — the
   discrete constraint solve returns a lapse differing from the continuum one by its own
   truncation error, so q carries an O(Δr²) ripple that the scheme then faithfully accelerates.
   `_discrete_equilibrium!` iterates the first integral against that same constraint solve,
   p_i ← H⁻¹(C − ln α_i) with C anchored on the central cell and α ← metric(ρ,p), to a fixed
   point. It converges in 12–17 iterations. It is also not optional: for the ρ_c = 7.993×10⁻³
   star at N = 800, **without** it the spread of q across the star is 1.05×10⁻⁴ and the interior
   residual 6.8×10⁻⁶; **with** it, 3.6×10⁻¹³ and 5.8×10⁻¹⁴.

`ShumPolytrope` (p = κρ², ε = ρ+p ⇒ h = 1+2κρ) has H = `log1p(2√(κp))` and
H⁻¹ = `expm1(H)²/4κ` in closed form, which hold full precision down to the floor; any other
barotrope falls back to H = ln((ε+p)/ρ) with a bisection inverse.

**Verification.** Everything below with **nothing subtracted**
(`../athenak-bdnk/analysis/dyngr1d_wb_check.jl`, figure `dyngr1d_wb.png`).

The residual must be read separately in the bulk and at the surface, because the two schemes
fail in different *places* and one number over the whole star cannot tell an exactly balanced
operator from a second-order one — which is why §7.15 could not see the difference. Momentum
residual |∂_t S̃| over the stellar **interior**, r < 0.95R, divided by the local gravity scale,
for the ρ_c = 1.28×10⁻³ star:

| Δr | 0.0479 | 0.0240 | 0.0120 | 0.0060 | order |
|---|---|---|---|---|---|
| plain operator (§7.15) | 3.50×10⁻⁵ | 8.72×10⁻⁶ | 2.18×10⁻⁶ | 5.45×10⁻⁷ | 2.00 |
| plain + discrete-equilibrium data | 4.84×10⁻⁵ | 1.20×10⁻⁵ | 3.00×10⁻⁶ | 7.49×10⁻⁷ | 2.00 |
| **hydrostatic reconstruction** | **3.85×10⁻¹³** | **8.54×10⁻¹³** | **1.99×10⁻¹²** | **3.79×10⁻¹²** | round-off |

Eight decades, and the last row does not converge — it *rises* linearly with the cell count,
which is the signature of accumulated round-off rather than of a truncation error. The
ρ_c = 7.993×10⁻³ star behaves the same way (7.84×10⁻⁵ → 1.35×10⁻⁶ against 3.08×10⁻¹⁴ →
2.31×10⁻¹³). The middle row is the control that matters: the new initial data on the old
operator is no better than before, so the gain is the scheme, not the data.

Evolved, over 200 M⊙, unkicked, nothing subtracted (ρ_c = 1.28×10⁻³ star):

| | N = 400 | N = 800 | N = 1600 | centre first moves at |
|---|---|---|---|---|
| plain operator | 4.31×10⁻³ | 2.55×10⁻³ | 1.48×10⁻³ | t = 1.0–1.1 |
| plain + discrete-equilibrium data | 4.19×10⁻³ | 2.52×10⁻³ | 1.47×10⁻³ | t = 1.0–1.1 |
| **hydrostatic reconstruction** | **9.51×10⁻⁴** | **5.52×10⁻⁴** | **2.65×10⁻⁴** | **t = 41.9 / 47.0 / 50.5** |

The last column is the sharpest statement available: with the plain operator the centre starts
moving within one M⊙, because the imbalance is spread through the bulk; with hydrostatic
reconstruction ρ_c sits at round-off (∼3×10⁻¹⁵) for ~45 M⊙, which is about the time a sound
signal needs to cross this star (R = 9.6, coordinate sound speed ≈ 0.3), and the delay **grows**
with resolution — consistent with the surface injecting less as it is better resolved, which is
also what the drift column does. Nothing in the interior moves the star; the drift is injected
at the surface and has to travel in.

Against the two independent frequency-domain eigensolvers (ε_c = 0.0015, N = 500, 220 R):

| operator | F_dyn vs Chandrasekhar 2.1236 kHz | F_frozen vs Cowling 4.0099 kHz |
|---|---|---|
| before §7.15 | 2.02, −4.9% | 4.6, +15% |
| §7.15, subtraction | 2.1304, +0.32% | 3.9994, −0.26% |
| **hydrostatic, nothing subtracted** | **2.1234, −0.01%** | **4.0020, −0.20%** |
| hydrostatic + subtraction | 2.1189, −0.22% | 3.9950, −0.37% |

The pure hydrostatic operator is the most accurate of the four, and by a factor ~30 on the
dynamical fundamental. That is the decisive evidence, because it is a comparison against an
independent method rather than against the engine's own initial data. It is also why the
subtraction is now **off** by default: it is the least accurate of the three post-fix rows,
since what it freezes into place is a state-specific residual that the star then oscillates
around. `test_dyngr.jl`: **56/56**.

**The overtones too.** `repro/radial_overtone_crosscheck.jl`, one of the scripts §7.15 left
owed, was re-run. It compares three independent routes — GHZ(1997) shooting, the
Sturm–Liouville matrix eigensolver and the time-domain engine — and the engine now lands on
every mode the two seeds excite:

| | F | H1 | H2 | H3 |
|---|---|---|---|---|
| SL eigenvalue (kHz) | 2.1236 | 5.8946 | 8.8243 | 11.6005 |
| DynGR1D periodogram peak | +0.05% | −0.04% | −0.03% | −0.03% |

The fundamental turns up in both seeds, at +0.05% from the homologous one and −0.12% from the
node-bearing one; H4 (14.3091 kHz) is not excited by either seed and so is not measured.

`repro/radial_shoot_crosscheck.jl` was re-run too, and its three routes close: shooting 2.1236,
SL 2.1236, engine 2.1234 — shoot-vs-engine **+0.01%**, where the value recorded in that script
had been 2.03 and would have printed a 4.6% disagreement.

**That benchmark had to be repaired before it could say so**, and both of its defects turned a
correct engine result into an apparent failure. Its classifier held only (F, H1, H2) while the
periodogram window reaches past H3, so the genuine H3 peak at 11.5967 kHz was assigned to H2 and
printed as a **31.4% error**. And its closing verdict asked whether H1 was resolved while looking
only at the node-bearing seed — which is orthogonal to H1's eigenfunction by construction — so it
printed "NO peak within 6% of H1" in the same report whose homologous probe carried H1 at −0.04%
with 0.95 of the peak power. The target list now runs to H4 and the verdict reports every mode
across both seeds and both probes.

**Nothing physical moved.** The Phase-1 runs (`../athenak-bdnk/analysis/p1_dyngr1d.jl`) were
regenerated with the new operator and reproduce §7.15 verdict for verdict and digit for digit:
the Font et al. unstable star given a −1% kick still collapses at t_AH = 53.71 with central
proper time τ_c = 10.01, horizon mass M_AH = 1.2734 and max 2m/r = 0.9598 (§7.15: 53.7, 10.0,
1.273, 0.9595–0.9598); kicked +5% it still migrates, settling at ρ_c = 1.47×10⁻³ against the
1.3×10⁻³ stable-branch star of that baryon mass; and the ε_c = 0.003 stable-branch star still
merely oscillates, 6% above its initial central density with max 2m/r = 0.47. A scheme change
that leaves an eight-decade mark on the residual and none at all on the physics is the outcome
to want.

**API.** `setup_dyngr(...; wellbalanced=:hydrostatic)` is the default; `:plain` (also `:none`,
`false`) is the §7.15 operator. `subtract` is now a separate, orthogonal keyword defaulting to
`false`; `wellbalanced=:subtract` reproduces the pre-2026-09-22 default (plain + subtraction).
`discrete_ic` controls the initial-data projection and follows the scheme.

**The honest negative: the surface.** The residual is exactly zero in the bulk and entirely
concentrated in the **single cell** where the star meets the constant-density artificial
atmosphere. That cell is not an equilibrium at all — ρ falls from ∼10⁻³ρ_c to the floor across
it — so no flux/source pairing can balance it, and summed over the whole star the residual comes
back up into the 10⁻⁶–10⁻⁸ range (unstable star 1.69×10⁻⁶ → 2.68×10⁻⁸ across the four
resolutions, second order; stable star 4.78×10⁻⁶ → 5.05×10⁻⁷ but **not** monotone — the N = 800
point is the worst of the four, because a one-cell residual depends on exactly where the surface
happens to fall between cell faces). That one cell is what launches the drift in the table above. So: the 10⁻⁶ target is met **for the operator** (interior residual
10⁻¹³, eight decades below it) and **not** for the time-integrated central density, which stays
at a few 10⁻⁴ and converges only at ≈Δr¹ because it is set by how well the surface is resolved.
Pass `subtract=true` to remove it if a strictly frozen background is wanted; it costs the 0.2%
accuracy shown in the frequency table.

Two further limits. The hydrostatic map is applied only where there is hydrostatic support to
preserve: it is switched off in the atmosphere (ρ ≤ ρ_cut), because H(p_atm) ∼ 10⁻⁷ against half
a cell of ln α ∼ 10⁻³ would return a face pressure eight orders of magnitude above the floor and
pump mass into the surface — that was a real failure mode during development, worth 10⁻³ of
spurious baryon-mass gain per 200 M⊙. A second guard (`hyd_fac`, default 10⁴) distrusts the map
where the profile steepens by more than four decades across half a cell; across sixteen
(star, resolution, box) combinations it excluded **at most one** cell — always the outermost
fluid cell, the very one the star/atmosphere face already leaves unbalanced — and none at all in
half of them. It is a safety valve against an unresolved profile, not a tuning knob. And the
magnetised (toroidal) sector is unchanged from §7.15: B ≠ 0 is not part of the first integral,
so a magnetised star is not well-balanced by this scheme and still needs `rebalance=true`.

**Follow-up owed.** Two of the six scripts §7.15 listed are now discharged:
`repro/radial_overtone_crosscheck.jl` and `repro/radial_shoot_crosscheck.jl`, above. A third,
`repro/viscous_fmode_scaling.jl`, was on that list in error — it names DynGR1D only in a comment
about a past audit and never runs it. What is still owed is `repro/viscous_damping.jl`, and the
two figures `viz/dyngr_1d.png` and `viz/dyngr_radial_validation.png`, which still show the
pre-§7.15 frequencies. The figures were **attempted and could not be produced here**: their
scripts need `viz/Project.toml`'s CairoMakie, which is declared but not installed in this
checkout, and instantiating it was left to the owner of the environment rather than done in
passing. Neither script hard-codes a frequency — both call the engine — so re-running them is
all that is required.

---

*Generated as part of the pre-publication audit. Suite state and all tabulated numbers
correspond to the commit in which this file was added.*
