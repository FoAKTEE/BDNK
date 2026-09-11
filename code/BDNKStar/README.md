# BDNKStar

A causal first-order (BDNK) viscous neutron-star and relativistic-hydrodynamics
toolkit in Julia. It **reproduces** published BDNK / first-order-hydro results
(then extends them), with every reproduction computed from the package API — no
hard-coded target numbers.

BDNK = Bemfica–Disconzi–Noronha–Kovtun: first-order viscous relativistic
hydrodynamics made causal/stable by a "frame" choice (relaxation times on the
hydrodynamic fields) rather than by adding new dynamical fields (Israel–Stewart).

Conventions: geometrized units `G = c = 1`. **Two length conventions are in use** and
must not be conflated — `solve_tov` is scale-agnostic, so the convention is fixed by the
EOS + central density you supply:

* **km-geometric** — εc in km⁻² (e.g. `piecewise_polytrope` via ρ[g/cm³]·7.4237e-19, or
  Bussières EOS1 with κ=100 km²). `star.M`/`star.R` are km; use `mass_solar(star)` and
  `f[km⁻¹]/kHz_to_km`.
* **M⊙-geometric** — the G=c=M⊙=1 convention standard in the NS-oscillation literature
  (the n=1, K=100 Font/Stergioulas/Kokkotas benchmark). `star.M` is **already in M⊙**
  (do *not* call `mass_solar`); `star.R × 1.4766` → km; `f·Msun_to_km·kHz_to_km` → kHz.

Because the Γ=2 polytrope is scale-invariant (M ∝ √κ), the same `ShumPolytrope(100.0)`
solution reproduces two different published stars (1.400 M⊙/14.15 km vs 0.948 M⊙/9.59 km)
depending on the reading. Dimensionless results are identical under both.

📋 **See [`VALIDATION.md`](VALIDATION.md)** for the unit conventions in full, the
quantitative validation table against published benchmarks (with enforced tolerances),
reproduction coverage, and the explicit list of known limitations / honest negatives.

## Layout at a glance

| Path | Role |
|---|---|
| `Project.toml` / `Manifest.toml` | Package manifest. **Stdlib-only** core (LinearAlgebra, Printf) so the test gate instantiates and runs fast. |
| `src/` | The physics package — one concern per file (see below). |
| `test/` | `Pkg.test` suite: the correctness gates (round-trip, causality, convergence, reproduction targets). |
| `repro/` | Self-contained reference-paper reproduction scripts (+ their `.txt` data dumps). Has its own [README](repro/README.md). |
| `viz/` | Figure generators (separate `CairoMakie` environment) → `figures/`. |
| `figures/` | 42 rendered reproduction figures (`.png`). |

The core package is deliberately dependency-light. Plotting (`CairoMakie`) lives
in a separate `viz/` environment so the test gate stays lean.

---

## `src/` — the physics package

Top-level module [src/BDNKStar.jl](src/BDNKStar.jl) wires the submodules below
and re-exports the public surface. Modular by design so the EOS/recovery trunk
is reusable by every later stage (radial → axial → nonlinear Cowling → 1+1D GR).

### Foundations
- **[Numerics.jl](src/Numerics.jl)** — dependency-free root finding. A guarded
  **Brent solver** (bracketing, derivative-free, superlinear), the workhorse for
  primitive recovery (a 1-D root find in the pressure).
- **[Units.jl](src/Units.jl)** — geometrized-unit constants ported from
  `NeutronStarOscillations.jl`: `Msun_to_km`, `kHz_to_km`, g/cm³ and dyne/cm² →
  km⁻² conversions, so backgrounds and transport coefficients are comparable to
  the benchmark codes.

### Equation of state & background star
- **[eos/EquationOfState.jl](src/eos/EquationOfState.jl)** — the shared EOS
  driver (STEP-0 trunk, the project's highest-leverage module). Two families
  behind one interface: `BarotropicEOS` (`p = p(e)`, cold stars, reference
  polytrope) and `GeneralEOS` (`p = p(ρ,ϵ)`, finite-T ideal gas, the BDNK
  heat-conduction sector). Exposes thermodynamically consistent sound speed
  `cs²`, partials `∂p/∂ρ|ϵ`, `∂p/∂ϵ|ρ`, enthalpy, temperature; enforces
  monotonicity + density floor; supports monotone-cubic **tabulated** EOS.
  `isentropic_idealgas(; Γ, K)` tabulates the ideal gas p=(Γ-1)ρϵ on an adiabat
  (p=Kρ^Γ, ε=ρ+p/(Γ-1)) as a barotrope, so ideal-gas microphysics flows through the
  whole TOV→f/p-eigensolver→viscous-QNM pipeline.
- **finite-T ideal-gas star** — `solve_tov_idealgas(; Γ, Γ_struct, K, ρc)` →
  `IdealGasStar` builds the hydrostatic structure (barotropic, as any static star must
  be) and reconstructs the genuine finite-temperature profiles T(r)=p/ρ, adiabatic
  c_s²(r) and equilibrium c_n²(r). It **exercises the BDNK heat-conduction sector**:
  a real temperature gradient dT/dr<0, the Caballero–Yunes criterion c_s²−c_n²(r), and
  causal BDNK heat-flux characteristic speeds (κ_Q>0, real & subluminal in a τ_P>1 frame).
  Figure: `idealgas_finiteT.png`. (The full thermal-perturbation δT/heat-flux QNM is a
  further extension.)
- **[eos/PiecewisePolytrope.jl](src/eos/PiecewisePolytrope.jl)** — GENERAL /
  REALISTIC cold-NS equations of state via the standard **Read–Lackey–Owen–Friedman
  (2009)** piecewise-polytrope parametrization (arXiv:0812.2163): a fixed SLy crust
  joined to 3 high-density polytropes `p=K_iρ^{Γ_i}` (fixed dividing densities
  ρ1=10^14.7, ρ2=10^15.0 g/cm³), each EOS fixed by four numbers (log₁₀p1, Γ1, Γ2, Γ3)
  with first-law ε-continuity. `piecewise_polytrope(:SLy|:APR4|:H4|:MS1)` (or the
  keyword builder) tabulates ρ→(p,ε,c_s²) into a `TabulatedBarotrope`, so any realistic
  EOS plugs straight into TOV → M–R and every QNM solver. **Validated:** TOV M_max and
  R_1.4 reproduce the published Read 2009 / Lackey–Wade values to **≤1.1%** for all four
  EOS (SLy 2.05 M☉/11.7 km, APR4 2.19/11.3, H4 2.01/13.9, MS1 2.75/14.9); the SLy
  1.4 M☉ ℓ=2 Cowling f-mode = **2.41 kHz**. Honest caveat (documented + tested): the
  soft fits SLy/APR4 turn superluminal only above ~10^15 g/cm³ (a known Read-parametrization
  feature; H4/MS1 stay causal everywhere). Figure: `general_eos_mr.png`.
- **[tov/TOV.jl](src/tov/TOV.jl)** — Tolman–Oppenheimer–Volkoff static
  background star. RK4 outward integration from a Taylor-seeded regular center,
  surface at `p→0`, Schwarzschild exterior match. The shared background for the
  linear and nonlinear analyses.
- **[rotation/SlowRotation.jl](src/rotation/SlowRotation.jl)** — **ROTATION** to
  O(Ω): the Hartle (1967) slow-rotation **moment of inertia** + Lense–Thirring
  frame dragging. Integrates the frame-drag ODE `d/dr[r⁴j ϖ′] + 4r³ j′ ϖ = 0`
  (ϖ≡Ω−ω, j=e^{−ν/2}√(1−2m/r)) on the TOV background and matches to the vacuum
  exterior → `J=R⁴ϖ′(R)/6`, `Ω=ϖ(R)+Rϖ′(R)/3`, `I=J/Ω` (central normalization
  cancels). `moment_of_inertia(star) -> SlowRotResult` (I in km³ and g cm², the
  dimensionless Ī=I/M³, the dragging profile ω(r)/Ω). **Validated:** Ī of the
  realistic EOS (SLy, APR4, H4, MS1) lies on the **Breu–Rezzolla (2016)** universal
  Ī–C relation to **≤1.6%** (SLy −1.3%, APR4 −0.2%, H4 −1.6%, MS1 +0.9%), with
  Lattimer–Schutz (2005) as cross-check and the constrained PSR J0737−3039A
  I=1.425×10⁴⁵ g cm² as the absolute-scale anchor (SLy 1.4 M⊙ ⇒ I=1.36×10⁴⁵,
  Ī=11.5; central drag ω(0)/Ω≈0.45). Verified two independent ways (volume-integral
  form and a different ODE formulation agree to ≤1.4×10⁻⁶). Figure: `slow_rotation.png`.
- **[rotation/RModes.jl](src/rotation/RModes.jl)** — **r-modes + viscous CFS
  instability window**: the l=m=2 Rossby mode (corotating σ=⅔Ω, inertial ω=−4⁄3Ω)
  is CFS-unstable to gravitational radiation at *all* spins and damped by viscosity,
  so the instability window is the cleanest test of the viscous sector. Implements
  the Lindblom–Owen–Morsink (1998) current-quadrupole GW growth `τ_GW`, shear `τ_sv`
  (η∝ρ^{9/4}T⁻²) and bulk `τ_bv` (ζ∝ρ²T⁶) damping integrals, and the critical
  spin `Ω_c(T)` window. **Validated** on the exact n=1 Lane–Emden benchmark (LOM98's
  Newtonian setup) to ≤0.2%: J̃ 0.02%, Ĩ 0.15%, τ_GW 0.09%, τ_sv 0.16%; window
  minimum Ω_c/Ω_K≈0.057 near T~2×10⁹ K (LOM98 ~0.06). The empirical gate **caught a
  real coefficient error** both literature-derivations shared (shear prefactor is
  (l−1)(2l+1)=5, not 15). Realistic-EOS windows (SLy/APR4/H4) confirm the minimum
  Ω_c/Ω_K is near-EOS-independent (0.055–0.058). **BDNK headline:** the causal
  frequency-dependent bulk viscosity `ζ_eff(ω)=ζ_NS/(1+(ωτ)²)` (relaxation-theory
  class, of which BDNK is an instance) suppresses high-T damping and *widens* the
  high-T edge of the window vs Navier–Stokes; for the r-mode στ~1 at τ~3×10⁻⁴ s.
  Honest caveat: τ_bv carries one O(1) eigenfunction normalization calibrated once to
  the LOM98 benchmark (the dropped O(Ω²) δΨ piece; LOM98 call this channel factor-2
  approximate) — so the τ_bv=6.99e8 s test is a self-consistency check, not an
  independent validation; the EOS-dependence and ζ_eff replacement act on the
  genuinely-computed integral. Figure: `rmodes_window.png`.

### BDNK transport, causality, recovery
- **[transport/Transport.jl](src/transport/Transport.jl)** — the BDNK transport
  coefficient container (shear η, bulk ζ, heat conductivity κ_Q, relaxation
  times τ_ε, τ_P, τ_Q = the "frame") plus named frames (conformal PMP, Shum).
- **[transport/Causality.jl](src/transport/Causality.jl)** — pointwise causality
  / stability monitor. Solves the BDNK characteristic-speed **biquadratic**
  `Λ₂c⁴ − 2Λ₁c² + Λ₀ = 0` (coefficients ported verbatim from the reference
  `BDNKCharacteristicSpeeds.jl`); flags real, non-negative, subluminal speeds.
- **[recovery/Recovery.jl](src/recovery/Recovery.jl)** — the conservative ↔
  primitive map (the recurring "wall"). Flat-space barotropic and general
  inversions, each reduced to one 1-D pressure root find. BDNK twist: conserved
  densities carry the first-derivative dissipative corrections, so recovery is
  the ideal inversion shifted by frozen sources; in the ideal limit it closes to
  machine precision (the STEP-0 gate, ≤ 1e-10).

### Conformal flat-space BDNK (the nonlinear engine)
- **[conformal/ConformalBDNK.jl](src/conformal/ConformalBDNK.jl)** — conformal
  (traceless, ζ=0, cs²=1/3) BDNK in flat space, slab symmetry. Verbatim port of
  Pandya's reference C solver (arXiv:2201.12317): BDNK stress tensor with
  dissipative corrections, the time-derivative primitive recovery (ξ̇, u̇ linear
  solve), the PMP luminal frame, and the Rankine–Hugoniot steady-shock state.
- **[conformal/ConformalEvolution.jl](src/conformal/ConformalEvolution.jl)** —
  full conformal-BDNK time evolution: 5th-order WENO reconstruction,
  Kurganov–Tadmor central flux, Heun (SSP-RK2) stepping, BDNK primitive solve
  each substage. Initial data: Gaussian, smooth shock, step. Reproduces the
  steady shock and smooth Gaussian evolution.

### Specialized physics
- **[perturbations/RadialModes.jl](src/perturbations/RadialModes.jl)** — radial
  (ℓ=0) stellar pulsations in the relativistic **Cowling** approximation
  (frozen metric). Frequency-domain matrix eigensolver of the pulsation ODE on
  the TOV background; eigenvalues ω² give mode frequencies (STAGE 1A,
  Caballero–Yunes).
- **[flows/Bjorken.jl](src/flows/Bjorken.jl)** — boost-invariant (0+1D) Bjorken
  flow in Milne coordinates, the PMP (2209.09265) test problem. Inviscid
  analytic solution + RK4 integrator (4th-order self-convergence, Q→16).
- **[dispersion/Kovtun.jl](src/dispersion/Kovtun.jl)** — linearized dispersion
  relations of first-order (general-frame) hydro (Kovtun 1907.08191): phase
  velocity `c_v(φ)` and shear-channel eigenfrequencies; reproduces the paper's
  stability claim (Im ω ≤ 0).
- **[viscous/IsraelStewart.jl](src/viscous/IsraelStewart.jl)** — Israel–Stewart
  bulk-viscous closure + causal limiter (Chabanov–Rezzolla 2311.13027): the
  contrast against which the first-order BDNK frame is measured (`is_contrast`).

### STAGE 3 — 3+1D Cowling non-radial modes
The polar (even-parity) f/p spectrum on the frozen-metric background, built four ways.
- **[perturbations/NonRadialModes.jl](src/perturbations/NonRadialModes.jl)** —
  frequency-domain eigensolver for polar (ℓ≥2) Cowling f/p modes. The reference
  spectrum: ℓ=2 f-mode **1.883 kHz** for the M=1.4 M☉ Shum star (matches the
  benchmark to 0.02%); the time-domain engines are validated against it.
- **[cowling3d/](src/cowling3d/)** — Cartesian 3+1D evolution. `Background3D` +
  `CowlingEvolve3D` (ideal + leading/Navier–Stokes-limit viscosity, `1/τ ∝ η̂`)
  and `CowlingBDNK3D` (full-frame causal BDNK). Exposes the surface problem: the
  staircased Cartesian boundary makes the full-frame velocity recovery blow up.
- **[spherical/](src/spherical/)** — the fix: a **boundary-conforming (r,θ)**
  engine with the surface on the coordinate line `r=R` (no staircase).
  `SphBackground` (TOV metric on the (r,θ) grid), `SphEvolve` (ideal Cowling —
  stable on the FULL star with no excision; a broadband Pₗ(cosθ) pluck reproduces
  the **whole polar spectrum** — f and the overtones p₁–p₄ — for **ℓ=2–6** to
  ≤4% of the eigensolver (ℓ≥3 to <1%), fig `sph_fp_spectrum.png`; the pole BCs are
  ℓ-independent for m=0, so only the angular factor changes), and **`SphBDNK`**
  (full-frame causal BDNK: continuity + momentum + the velocity **recovery**,
  the part the Cartesian surface destabilized — now stable). The viscous result
  lives here: of the three frame channels, only the shear-viscous **momentum**
  force `ν_mom∇²δS` dissipates the mode (it makes the (δρ,δS) oscillator damped);
  the conduction/heat term `κ_Q c_s²∂ε` is **reactive** (frequency/stability,
  not damping) and `τ_R` is the causal-frame stiffness. The faithful diagnostic
  is the quadratic **mode energy** (`mode_energy`), not amplitude decay — the
  latter is masked by Kreiss–Oliger dissipation and by dephasing of a
  non-eigenmode seed. Energy dissipation increases monotonically with η̂ while
  the evolution stays bounded on the full star. Figure: `sphbdnk_viscous.png`.
- **[perturbations/PolarViscousModes.jl](src/perturbations/PolarViscousModes.jl)**
  — the polar **viscous QNM solver** (even-parity counterpart of the axial track).
  Frequency-domain matrix eigenproblem: the 1D ℓ-reduced linearised polar BDNK–Cowling
  operator `L_ℓ` is assembled and `eigvals(L_ℓ)` gives the complex QNMs (`polar_qnm`,
  `polar_bdnk_operator`; freq=|Im λ|, damping γ=−Re λ), modes tracked by complex
  continuation from the η̂=0 ideal modes. η̂→0 reproduces the f/p modes with γ→0; η̂>0
  gives **γ>0 rising ∝η̂, p₁ damping more than f** (higher k ⇒ more shear), consistent
  with the `SphBDNK` time-domain dissipation. The full spectrum also shows the
  heavily-damped frame/relaxation (non-hydrodynamic) modes — the only "viscosity-driven"
  addition in the polar sector (no long-lived η-mode analog). **SEMI-QUANTITATIVE**
  (adversarially reviewed): the damping and its η̂-scaling are robust, but it reliably
  tracks only **f and p₁** (higher overtones mix with frame modes — use NonRadialModes),
  is valid for **Nr≲200** (marginally-stable operator), and the frequency carries a
  ~5–8% frame systematic at η̂→0 plus a large reactive viscous pull (~20–30%) at finite
  η̂; for precision QNMs use the axial full-GR shooting solver. Figure: `polar_viscous_qnm.png`.
- **[perturbations/AxialViscousModes.jl](src/perturbations/AxialViscousModes.jl)** —
  the **axial (odd-parity) viscous QNM solver**, full GR (no Cowling), ported from the
  validated `repro/axial_*` track (Bussières et al. 2604.13208): interior shooting +
  exterior vacuum (Leaver continued fraction) + Wronskian/log-derivative matching +
  complex-ω root finding (`axial_qnm`; convention ω=2πf−i/τ). Reproduces Bussières
  Table II to **<0.05%** (inviscid w-mode 10.50 kHz/29.54 μs; the viscosity-driven
  **η-modes**). Here viscosity DRIVES new dynamics (the axial fluid sector is otherwise
  trivial) — the contrast with the polar sector where it only damps. Figure
  `axial_viscous_qnm.png` (w-mode f, τ vs central η_c: f falls, τ rises with viscosity).
- **[perturbations/GravityModes.jl](src/perturbations/GravityModes.jl)** — the
  relativistic Cowling polar eigensolver WITH buoyancy (the Schwarzschild-discriminant
  term the f/p solver drops), so it yields f, p AND **g-modes** (`gmode_spectrum`,
  `brunt_vaisala`). g-modes are a **genuinely new mode family** here: they appear ONLY
  when the finite-T star is stratified (adiabatic Γ ≠ structure Γ_struct ⇒ c_s²≠c_e² ⇒
  N²≠0). Validated: barotropic (Γ=Γ_struct) ⇒ N²≈0 and NO g-modes, f/p matching
  NonRadialModes to <3%; stratified, convectively stable (Γ>Γ_struct ⇒ N²>0) ⇒ a real
  g-mode tower g₁>g₂>… BELOW the f-mode (p-modes above), converged and rising with the
  stratification. (U,V) system of Jaikumar/Shirke (2506.18892) with the Gaertig–Kokkotas
  discriminant sign. Cowling ⇒ ~15% absolute accuracy (qualitative-to-semiquantitative).
  Figure: `gravity_modes.png`.
- **[perturbations/PolarGRModes.jl](src/perturbations/PolarGRModes.jl)** — the
  **FULL general-relativistic** polar QNM solver (the Cowling approximation DROPPED):
  the Lindblom–Detweiler interior {H₁,K,W,X} metric+fluid system (DL85 1985) matched
  to the exterior **Zerilli** equation with the outgoing-wave BC (Chandrasekhar–Detweiler
  transform of the Leaver continued fraction), complex-ω root find (`polar_gr_qnm`). The
  even-parity analog of the axial full-GR solver — now the f-mode acquires a **gravitational-
  wave damping time** (Im ω ≠ 0, impossible in Cowling). Validated: the f-mode is
  **f=2.86 kHz, τ_GW=0.115 s** — *below* the Cowling f-mode (3.30; ratio 0.87, the GR
  overestimate-correction) and agreeing with the Andersson–Kokkotas (1998) universal
  relations to ≈8% (frequency) / ≈10% (damping); a polar w-mode (9.2 kHz, 64 μs) is also
  recovered. Figure: `polar_gr_qnm.png`. (Getting the f-mode right required correcting the
  LD X′ eq. against the DL85 paper: the H₀-term sign, the eq-(5) X-term e^{+ν/2}, and the
  ½ on the V coefficient — the f-mode's GW damping, Im ω/Re ω ~ 10⁻⁴, is exquisitely
  sensitive to these.)

## `test/` — correctness gates

Run with `julia --project=. -e 'using Pkg; Pkg.test()'`.
[runtests.jl](test/runtests.jl) drives:

| Test | Gate |
|---|---|
| `test_eos` | EOS thermodynamic consistency |
| `test_recovery` | Barotropic polytrope round-trip ≤ 1e-10 |
| `test_causality` | Characteristic speeds solve the biquadratic (Vieta residual) |
| `test_tov` | Reproduce Bussières EOS1 star (M = 1.27 M☉, R = 8.86 km) |
| `test_conformal` | Rankine–Hugoniot steady shock (PMP/Pandya) |
| `test_conformal_evolution` | Uniform state exactly preserved |
| `test_conformal_convergence` | Ordered self-convergence (2nd–5th order) |
| `test_radial` | Radial Cowling modes: convergence + stability |
| `test_heat_criterion` | Caballero–Yunes heat-conduction criterion `cs² − cn²` |
| `test_bjorken` | Bjorken RK4 Q→16 + analytic + diagnostic |
| `test_kovtun` | Kovtun dispersion stability + analytic limits |
| `test_is_contrast` | Israel–Stewart vs BDNK Reynolds-limit contrast |
| `test_nonradial` | ℓ≥2 Cowling f/p eigensolver (5-mode benchmark) |
| `test_cowling3d` | Cartesian 3+1D Cowling f-mode + leading viscosity |
| `test_bdnk_viscous` / `test_bdnk_fullframe` | leading-/full-frame BDNK: hyperbolic stability + `1/τ ∝ η̂` |
| `test_spherical` | boundary-conforming (r,θ): metric, ideal f-mode, **f+p₁–p₄ polar spectrum (ℓ=2,4) ≤6%**, full-frame BDNK energy dissipation ↑ with η̂ |
| `test_polar_viscous` | polar viscous QNM eigenproblem: η̂→0 non-dissipative f/p; **γ>0 ∝η̂**, p₁ damps more than f |
| `test_axial_viscous` | axial full-GR viscous QNM: Bussières Table II w-mode (10.50 kHz/29.54 μs) + viscous η-mode shift, <1% |
| `test_idealgas` | ideal-gas EOS through TOV + f/p eigensolver + polar & axial viscous QNM (EOS-generality) |
| `test_idealgas_finiteT` | finite-T ideal-gas star: T(r) gradient, Caballero–Yunes criterion, causal BDNK heat-flux speeds |
| `test_gmodes` | **g-modes** (relativistic Cowling + buoyancy): barotropic⇒none; stratified⇒real g-tower below f; converged |
| `test_polar_gr` | **full-GR polar f-mode** (Lindblom–Detweiler + Zerilli): f<Cowling, finite GW damping τ, matches Andersson–Kokkotas to ≈8–10%; + polar w-mode |
| `test_shum_frame` | Shum BDNK **fluid-frame** characteristic structure: production-frame c₊=√3c_s, c₋≈0.0183c_s; well-posed 0<q̂<ŝ (eq 71); c₊≤1 along the M=1.4 star |
| `test_general_eos` | **general/realistic EOS** (Read 2009 piecewise polytrope SLy/APR4/H4/MS1): TOV M_max & R_1.4 vs published to ≤1.1%; causal to ~10¹⁵ g/cm³; SLy 1.4 M☉ f-mode 2.41 kHz |
| `test_slow_rotation` | **rotation** (Hartle slow-rotation moment of inertia I): Ī–C vs Breu–Rezzolla universal relation ≤1.6% for SLy/APR4/H4/MS1; I normalization-independent; SLy 1.4 M☉ I=1.36×10⁴⁵ g cm²; central drag ω(0)/Ω≈0.45 |
| `test_rmodes` | **r-modes + CFS instability window** (LOM98): l=m=2 freqs σ=⅔Ω/ω=−4⁄3Ω; n=1 Lane–Emden τ_GW/τ_sv to ≤0.2%; window min Ω_c/Ω_K≈0.057 at T~2×10⁹ K; causal ζ_eff(ωτ) widens high-T edge vs Navier–Stokes |
| `test_convergence` | **convergence + numerical-vs-physical**: RMode 2nd-order (Richardson p=2); polar-viscous stable in its Nr≲200 window; **SphBDNK viscous damping = MIXED** — ν_mom dissipative & κ_Q reactive (physical signs), but the time-domain rate is σ_ko-/resolution-contaminated (use the converged eigen-solvers for quantitative damping). Fig `convergence_viscous.png` |
| `test_frame_independence` | **BDNK frame-independence** (the headline gauge check): physical observables invariant under frame change at fixed (η,ζ,κ_Q) — dispersion sound speed/damping flat to ~0.003% while non-hydro gaps move ∝1/τ; axial w-mode flat to ≤0.12%; the polar ~5–8% is a **solver offset, not gauge**. Fig `frame_independence.png` |
| `test_cross_method` | **cross-method damping** (eigenvalue vs dissipation-integral vs time-domain): polar **p1 agrees to 85%** (rigorous independent check); polar f-mode contaminated (~3.4×, semi-quant); axial w-mode opposite-sign (spacetime mode, fluid-dissipation inapplicable — physical); time-domain sign/trend confirmed. Fig `cross_method_damping.png` |
| `test_universal_relations` | **universal relations across sequences** (not single stars): AK1998 f-mode holds EOS-independently (full-GR within ~6–14%, 1/τ within ~17%; Cowling collapses 3 EOS to one curve, RMS 5.8%); Breu–Rezzolla Ī–C across full sequence (4 EOS, ≤5.6%); r-mode window near-universal. Fig `universal_relations.png` |
| `test_entropy` | **entropy / second law**: BDNK entropy production (2ησ²+ζθ²+κ_Q-heat)/T positive-definite; shear **energy-entropy consistency ∫σ_S·T dV = 2γE to 100.0000%**; κ_Q "reactive" reconciled (heat entropy κ_Q(∂δT)²/T²>0, but δT sector un-evolved in Cowling — no 2nd-law violation). Fig `entropy_second_law.png` |
| `test_causality_map` | **causality phase diagram + reconciliation**: `Causality.jl` biquadratic vs Kovtun/Shum = **documented convention split, not a bug** (char-speed ≡ Kovtun large-k slope to 1e-16; only frame ratios ŝ,â map); production frame causal+well-posed for the soft Shum star at every radius, but **c₊=√3·cs ⇒ superluminal in the core of realistic stiff EOS** (SLy central cs>1/√3). Fig `causality_phase_diagram.png` |
| `test_secondary` | **secondary robustness** (5): linearity E∝A² (solver fully linear); long-term stable (inviscid growth *saturates*, no blow-up, viscous decays); f-mode BC-reflection-insensitive (0.02% under reflecting node); ζ_eff(ωτ) ≡ Israel–Stewart bulk form; Newtonian limit — Cowling f→√(GM/R³) ✓, full-GR f-mode low-C ill-conditioned (excluded). Fig `secondary_tests.png` |
| `test_noise_excitation` | **Gaussian-noise mode excitation**: a valid+causal band-limited random pluck (`seed_sphbdnk_noise!`) excites the physical ℓ=2/4 f/p spectrum — peaks **seed-independent** in frequency (scatter <0.04 kHz, amplitudes vary 60×), matching the Cowling eigensolver to ~2–8%; viscosity damps the peaks. Fig `noise_excitation.png` |
| `test_fvcartesian` | **FV flux-conservative Cartesian + atmosphere** (`fvcartesian/`): Valencia GRHydro, HLL + MinMod, well-balanced. **1D radial: works** — exactly well-balanced static star + grid-converged radial mode 3.08 kHz. **3D Cartesian: FV-insufficient** — static star well-balanced, but the ℓ=2-perturbed staircased surface drives a secular instability that *worsens with resolution* (needs DG, not refinement). Fig `fv_cartesian.png` |
| `test_dg` | **nodal DG (RKDG) fallback** (`dg/`, 1+1D shock/radial, 2+1D axisymmetric): **shocks captured** — SR shock tube vs exact Riemann (front/state match, no oscillations, positivity, convergent); radial star well-balanced (f=3.2 kHz, ~FV); **staircase cured** — 2D ℓ=2 drift *improves* with refinement (2.29→0.92%, inverse of FV). f-mode frequency not extracted in 2D. Fig `dg_shock_capturing.png` |
| `test_dg3d` | **full 3+1D nodal DG** (`dg/DGCart3D.jl`): 3D tensor-product RKDG GRHydro. 3D static well-balanced (genuine lake-at-rest cancel, drift 2e-13); 3D shock (axis+diagonal, positivity); **staircase cured in 3D** (ℓ=2 drift improves ~200× with refinement, inverts the 3D-FV failure); **non-axisymmetric m≠0 captured** (Y₂₂ proper-projection genuine, oscillates) — the capability 2D-axisymmetric couldn't do. f-mode freq compute-limited (foundation, not production). Fig `dg_3d.png` |
| `test_tidal` | **tidal deformability + dynamical tides** (`TidalDeformability.jl`): static ℓ=2 Love k₂/Λ (Hinderer/PPL y(r) ODE) — Λ(1.4) SLy 297/APR4 249/H4 860 within ~4% of published; **I-Love-Q completed** with the moment of inertia (Yagi–Yunes, 0.68%); GW170817 (SLy/APR4 <800 consistent, H4/MS1 disfavored). Dynamical Λ_eff(ω): static limit exact, f-mode resonance enhancement. Fig `tidal_deformability.png` |
| `test_doubly_diffusive` | **new viscous+heat modes** (`DoublyDiffusive.jl`): GSF rotational + thermohaline/semiconvective doubly-diffusive instabilities under differential rotation. Ideal limits reduce **exactly** to Solberg–Høiland (6.9e-18) and Brunt–Väisälä; doubly-diffusive hallmark (growth→0 as diffusivity→0). Honest: known families (GSF 1967, Stern 1960) realized in BDNK; **causal effect channel-dependent** — thermohaline gets a high-k UV cutoff, GSF only an O(sτ) growth-slowdown (boundary τ-invariant). Fig `doubly_diffusive_modes.png` |
| `test_dyngr` | **1+1D DYNAMICAL GR — Stage 2** (`dyngr/DynGR1D.jl`): the **first non-Cowling engine** — radial Valencia hydro + constrained metric (Hamiltonian constraint + polar slicing) **re-solved each step from the matter**. Metric genuinely evolves (freeze→Δ=0 exactly, dynamical→moves; decisive); TOV stationary (drift 2.5e-11, no cheat); **collapse-to-BH** (lapse α_c→0.001, 2m/r→0.95, ρ_c×4.8; non-collapse control). Radial mode **independently validated**: the dynamical engine (2.02 kHz) matches the corrected Chandrasekhar SL eigensolver (2.124 kHz) and an independent GHZ shooting solver (2.124 kHz, agree to 6 digits) within ~5%, and the published Kokkotas–Ruoff n=1 κ=100 polytrope band (2.1–2.3 kHz, our exact EOS family); sign-of-shift resolved (ω²_GR<ω²_Cowling, GR softens — the prior eigensolver had a P,Q↔W exponent swap, now fixed). F²→0 at M_max behavioral. Figs `dyngr_1d.png`, `dyngr_radial_validation.png` |

---

## `repro/` — reference-paper reproductions

Self-contained scripts that `include` the package, implement a paper's method on
the API, run, and print a numeric comparison against the published target.
Covers Stage 1B axial QNMs, Stage 1C nonlinear Cowling BDNK (Shum), PMP viscous
1D, and the analytic dispersion/causality fleet. Each `.jl` has paired `.txt`
data dumps. See **[repro/README.md](repro/README.md)** for the full module table,
achieved accuracies, and inter-module dependencies.

Run convention (from the package root so relative `include`s resolve):
```bash
cd code/BDNKStar
JULIA_NUM_THREADS=4 julia --project=. repro/<module>.jl
```

---

## `viz/` — figure generators → `figures/`

One script per figure, in a separate environment ([viz/Project.toml](viz/Project.toml),
`CairoMakie` only). Each activates its own env, includes `../src/BDNKStar.jl`,
and renders a `.png` into `figures/`. Naming mirrors the figures (e.g.
[viz/tov_reproduction.jl](viz/tov_reproduction.jl) → `figures/tov_reproduction.png`,
`viz/kovtun*.jl` → the Kovtun dispersion panels, `viz/shum_*.jl` → the Shum
nonlinear-Cowling QNM/decay figures, `viz/coverage_map.jl` →
`figures/coverage_map.png` reproduction-coverage overview).

```bash
julia --project=code/BDNKStar/viz code/BDNKStar/viz/<name>.jl
```

---

## Reproduction stages (roadmap)

`STEP 0` shared EOS + primitive recovery (this trunk) → `1A` radial linear
(Caballero–Yunes 2506.09149) → `1B` axial linear (Redondo-Yuste 2411.16841 /
Bussières) → `1C` nonlinear Cowling (Shum 2509.15303) → `2` 1+1D dynamical GR +
realistic EOS + collapse → `3` 3+1D Cowling non-radial modes → `4` production.
