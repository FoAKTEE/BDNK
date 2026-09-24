# Plan — GR-coupled 1+1D BDNK spherical collapse in AthenaK

*DAG node `s2.gr_coupling` (Fork A), 2026-09-20. Decision `{{GR_BACKEND}}` = reuse an NR
backend rather than hand-roll the Einstein evolution; this plan makes that concrete.*

## 0. What "1+1D in AthenaK" can and cannot mean

AthenaK is a Cartesian, Kokkos code. Its dynamical-spacetime GRMHD (`src/dyn_grmhd`,
`DynGRMHDPS<EOSPolicy, ErrorPolicy>`) evolves the Valencia variables on a Z4c spacetime
(`src/z4c`, puncture gauge: 1+log slicing + Gamma-driver shift), with the matter coupled back
through `SetTmunu` and the geometric sources through `AddCoordTerms`. There is no curvilinear
or 1D-spherical GR path. The unmagnetized TOV problem (`src/pgen/dyn_grmhd/dyngr_tov.cpp`,
`inputs/dyn_grmhd/whisky_tov.athinput`) already does what we need for the ideal fluid: an
octant [0, L]³ with reflecting inner faces, K=100 Γ=2 polytrope, and a `v_pert` radial kick
(v_r = ½ v_pert (3x − x³), x = r/R); a `<z4c>` block makes the spacetime dynamical, an `<adm>`
block alone freezes it (Cowling, the way `dyn_grmhd.cpp:231` reads it).

So the honest statement of the task is:

- **the genuine 1+1D reference is ours**: `code/BDNKStar/src/dyngr/DynGR1D.jl` — areal gauge,
  polar slicing, constrained metric solve every substage; its collapse test (`test_dyngr.jl`
  §4: ε_c = 0.003, inward kick A = −0.03 → α_c → 0, 2m/r → 1) is the ideal-fluid baseline;
- **the AthenaK run is the spherically-symmetric problem on the octant grid** with a
  free-evolution (Z4c) spacetime — a different gauge, a different formulation, the same
  physics. That is the point: an independent code with a real NR backend against which the
  1+1D engine, and later the 3D non-spherical extension, are checked.

A true 1+1D reduction *inside* AthenaK (spherical coordinates, radial ADM solve) would be a
separate code project, not a use of AthenaK; it is not in this plan.

## 1. Phases and gates

### Phase 0 — build, run the stock TOV, fix the scale of the problem (≈ 1 week)

1. **Machine.** This laptop has no CMake/compiler/container (VALIDATION §7.11) and the
   checkout in `ref-code/athenak` (commit c5a0d7f, Kokkos 4.7.02, a `build-cpu/` from a Linux
   host) is gitignored. Options, in order: WSL2 Ubuntu on the laptop for the 64³ development
   runs (8 cores; Kokkos OpenMP); the Caltech cluster with a CUDA build for the 128³ / SMR
   production matrix. Fork `athenak` under `FoAKTEE/` on a branch `bdnk`; the inputs, analysis
   scripts and validation record live in **this** repo under `code/athenak-bdnk/`.
2. **Build** `cmake -D PROBLEM=dyngr_tov -D Kokkos_ENABLE_OPENMP=ON` (CUDA variant on the
   cluster). Run `whisky_tov.athinput` as shipped (64³ octant over [0, 102.4]³, ρ_c = 1.28×10⁻³,
   v_pert = −0.024, `<adm>` only = Cowling). Record: cost per step and per M⊙, ρ_c(t), the
   radial spectrum (F = 2.686, H₁ = 4.550 kHz for this star — our `radial_cowling_spectrum`),
   M_b conservation, atmosphere behaviour. This is also the AthenaK counterpart of the
   DG–FD hybrid star of VALIDATION §7.12–7.13 (same star, same kick family) — one figure.
3. **Turn the spacetime on**: add the `<z4c>` block (puncture-gauge defaults from
   `inputs/z4c/onepuncture/z4c_onepuncture_ahf.athinput`: `lapse_oplog=2`, `shift_eta=2`,
   `diss=0.02`, `damp_kappa1/2`), same star, no kick. Gate: the star stays static (ρ_c drift
   at the 10⁻³ level over 1000 M⊙, Hamiltonian constraint bounded) and the full-GR F-mode
   appears **below** the Cowling one — our `DynGR1D` gives F_GR = 2.12 kHz for the ε_c = 0.0015
   star (Kokkotas–Ruoff band 2.0–2.45 kHz), so run that star too.

### Phase 1 — ideal-fluid collapse: AthenaK vs the 1+1D engine (≈ 1–2 weeks)

The cross-code gate before any BDNK physics.

> **Status 2026-09-22 — the blocker found here has been fixed.** Setting the comparison up showed
> that `DynGR1D`'s raw operator carried a 5–16% *resolution-independent* static momentum imbalance;
> it is now well-balanced by construction (VALIDATION.md §7.15): the residual converges at second
> order, the engine reproduces the Chandrasekhar full-GR eigenvalue to **+0.32%** and the Cowling
> one to −0.26% (it was −4.9% and +15%), the two `wellbalanced` settings agree on every collapse
> number, and the outward-kicked unstable star **migrates** to the stable branch instead of
> dispersing. The 1+1D reference for this phase is therefore usable again. One recorded result was
> overturned in the process: the old "collapse" of the ε_c = 0.003 star was an artifact of the
> imbalance — that star is on the stable branch and merely oscillates — so the collapse benchmark
> is now the unstable-branch star of Font et al. 2002.
>
> Two practical points, both measured: an **unperturbed** star on the unstable branch is decided by
> the sign of the truncation error, not by physics, so both migration and collapse carry an explicit
> signed kick; and the kick must sit above the code's own noise (AthenaK's Z4c star breathes at
> ±0.8% in ρ_max, so the migration run uses v_pert = +0.05, the collapse runs −0.01 and −0.03).
> `p_pert` in AthenaK's `dyngr_tov.cpp` is dead code — a local `Real p_pert = 0.` shadows the
> parameter — so `v_pert` is the only trigger.

> **Status 2026-09-23 — the p1a comparison was not comparing the same initial data.**
> AthenaK's `dyngr_tov` pgen writes the Valencia primitive `w0(IVX+a) = vr*x_a/r` on an
> ISOTROPIC grid whose spatial metric is g_ij = psi^4 delta_ij with psi^4 = (r_schw/r_iso)^2,
> and the Lorentz factor is built as W = sqrt(1 + g_ij u^i u^j), so that primitive is u^i = W v^i
> in COORDINATE indices and the orthonormal speed it carries is **psi^2 * vr**. `DynGR1D`'s
> `seed_dyngr_velocity!` sets the orthonormal speed directly. The same `v_pert` is therefore a
> physically STRONGER kick in AthenaK, by psi^2 = 1.69 at the centre falling to 1.26 at the
> surface — least-squares equivalent **1.41x**. AthenaK's v_pert = -0.030 acts like DynGR1D
> A = -0.0423. (The isotropic radius integrated from our TOV solution reproduces AthenaK's
> printed R_iso = 6.52823 to 6.52820, so the transformation is confirmed.)
>
> That matters here because this star is AT A BIFURCATION: DynGR1D's collapse threshold for it
> lies between A = -0.042 (rho x1.57, bounces) and A = -0.050 (rho x5.08, alpha -> 0.001), i.e.
> exactly where AthenaK's effective kick lands. Correcting the seed closes part of the gap —
> DynGR1D with the psi^2-weighted profile reaches rho x1.62 instead of x1.31 — but not all of
> it: AthenaK still collapses (rho x25.6, alpha_min 0.040) where DynGR1D bounces.
>
> Two consequences for this phase. (a) **The eps_c = 0.003 star cannot serve as the cross-code
> gate.** A setup whose verdict flips under a 40% change in kick amplitude tests the bifurcation,
> not the codes. The Font unstable star does serve: it is well past the turning point, both codes
> collapse it, and their central proper times agree (AthenaK 10.43 at alpha < 0.05, DynGR1D 10.01
> at the horizon). (b) Any future comparison must match the ORTHONORMAL velocity, not the input
> parameter — either seed DynGR1D with the psi^2(r) weight or scale AthenaK's `v_pert` by 1/1.41.
>
> Separately: refinement is NOT the explanation. One extra static level at the centre
> (dx 0.25 -> 0.125) reproduces the unrefined AthenaK run to four significant figures in both
> rho_max and alpha_min through t = 30, so AthenaK's answer is resolution-converged.

1. Two collapse setups, both already in our hands: (a) the `DynGR1D` test star ε_c = 0.003 with
   A = −0.03 (`seed_dyngr_velocity!(profile=:linear)` ↔ AthenaK `v_pert` — match the profile,
   the AthenaK one is cubic); (b) the Font et al. 2002 unstable-branch star (ρ_c = 8.0×10⁻³,
   K=100, Γ=2; migration without a kick, collapse with an inward kick — confirm the numbers
   against `ref-paper/sources` before running).
2. **Compare gauge-invariantly.** `DynGR1D` is in polar slicing/areal radius, AthenaK in the
   puncture gauge: α_c(t) and ρ_c(t) are not comparable in coordinate time. Compare ρ_c against
   the central proper time τ_c = ∫α_c dt, the apparent-horizon formation (AthenaK: `z4c`
   fastflow/`horizon_dump`; ours: max 2m/r → 1 on the polar-sliced grid) by AH mass and by
   the baryon mass inside it, and the fraction of M_b outside the horizon at late times.
   Gate: AH mass to a few %, collapse-or-migration verdict identical for both setups.
> **Result 2026-09-23 — the Font MIGRATION run completed and agrees qualitatively, not
> quantitatively.** It is the only dynamical-spacetime AthenaK run in this phase that finished:
> exit 0, the full 800 M_sun, zero primitive-solve errors, where both collapse runs lost the
> 3-metric at the centre. The star migrates to the stable branch in both codes, so the VERDICT
> gate is met on this setup. On the gauge-invariant clock, over the common window tau <= 396:
>
> | | AthenaK | DynGR1D | diff |
> |---|---|---|---|
> | min rho_c | 2.823e-4 | 3.041e-4 | -7.2% |
> | min / max alpha | 0.2785 / 0.8423 | 0.2725 / 0.8145 | +2.2% / +3.4% |
> | dominant radial period (periodogram) | 121.3 | 108.0 | **+12.3%** |
> | mean rho_c over the window | 1.215e-3 | 1.545e-3 | -21.4% |
>
> The envelope agrees at a few percent; the PERIOD does not, and the 21% gap in the windowed mean
> is a consequence of that phase drift rather than an independent discrepancy.
>
> **The simplest reading is that AthenaK's 64^3 octant is not converged for this problem**, and
> every 10-20% discrepancy in this phase is its discretisation error: its initial central density
> is already 3.0% low at that resolution (0.8% with one refinement level), and the Font collapse
> pair moved its observables by 11-19% between exactly those two resolutions. No convergence
> information exists for the migration run — it was done at one resolution only.
>
> **And the cause of the migration discrepancy is identifiable: the run loses 13.3% of its
> baryon mass, because THE BOX IS TOO SMALL.** (An earlier note in this file called it surface
> shedding into the atmosphere; the density slices disprove that and it is corrected here.) The
> x3 = 0 slices show the expanded envelope reaching the outermost cell layer with rho of 1e-6 to
> 7e-6, four decades above the 1e-10 atmosphere, and `ox_bc = diode` is outflow-only:
>
> | t | rho at the outer face | radius where rho < 1e-6 |
> |---|---|---|
> | 0 | 1.0e-10 (atmosphere) | 4.38 |
> | 50 | 1.0e-10 | 10.88 |
> | 100 | 6.8e-06 | 15.88 (clipped by the box) |
> | 300 | 2.3e-06 | 15.88 |
> | 600 | 1.2e-06 | 15.88 |
>
> The timing settles it: the loss is EXACTLY zero through t = 75, across the whole first
> excursion down to 4% of the initial central density (its minimum is at t = 92) while the
> envelope is still at r ~ 11; it starts at t ~ 100, precisely when the envelope first touches
> the boundary, and then accumulates -- 0.45% at t = 100, 0.8% at 200, 7.6% at 300, 12.2% at 600,
> 13.3% at 800. The 1+1D counterpart conserves M_b to ~1e-4 on the same problem because its grid
> runs to 4R.
>
> A star that has lost 13% of its mass oscillates at a different period, and the measured period
> difference is 12.3%; those two being that close makes the outflow the leading candidate for the
> whole quantitative discrepancy on this setup, ahead of gauge or formulation.
>
> **Fix under test (`p1b_migration_font_box32`):** domain doubled to [0,32] with one static level
> put back over the original [0,16], so the star keeps exactly its former cell size (0.25) and
> the former timestep, at 1.875x the cells rather than the 8x a uniform 128^3 over the larger box
> would cost. Verified controlled: identical mesh resolution on the star, identical dt, and
> rho_c(0) and alpha_min(0) bit-for-bit identical to the [0,16] run.
>
> **Result 2026-09-24 — the box fix works; its refinement interface then fails at t = 221.5.**
> `p1b_migration_font_box32` (domain [0,32], one static level back over [0,16], identical star,
> cell size and dt) was run to a target of t = 300. Baryon mass, at matched times:
>
> | t | box32 | original [0,16] |
> |---|---|---|
> | 100 | -1.3e-6 | -4.5e-3 |
> | 150 | -4.7e-6 | -8.1e-3 |
> | 200 | -5.5e-6 | -8.0e-3 |
>
> The loss is gone — a factor ~1450 at t = 200 — which confirms the diagnosis above: the mass left
> through the diode boundary. The run is also bit-deterministic at fixed thread count: 667
> samples identical to an earlier copy killed at t = 166.5. The density trajectory tracks the
> original to 1-1.5% through t ~ 190, then drifts in phase (9-13% by t = 210-220), as it should
> for a star that has kept mass the original lost.
>
> **Then a NEW failure, unrelated to the star: a Z4c blow-up at the static-refinement interface.**
> First NaN just after t = 219, last physical sample t = 221.5, at (9.375, 9.375, 15.875) — the
> last fine cell before the fine/coarse boundary at z = 16, not the centre. The spacetime is what
> failed: lapse -2.35, psi^4 = -25.0, g_zz = -417, K_zz ~ -1.2e13, with beta^z = 3.67 while
> beta^x, beta^y ~ 1e-3. It was sudden and local: the global H and C norms were flat right up to
> it (H ~ 5e-5, C ~ 2e-3). The original [0,16] run had no refinement interface and ran cleanly
> to t = 800, so enlarging the box through SMR traded the boundary outflow for an interface
> instability.
>
> Consequence: the t = 200-300 window, where the original lost most of its mass (0.8% -> 7.6%),
> was NOT reached. The conclusion that the loss is fixed does not depend on it — that loss needs
> the envelope to reach the domain boundary, which box32 removes by construction — but the
> migration PERIOD with the mass conserved, the comparison the fix was for, is still unmeasured.
> Options: a uniform [0,32] at 128^3 (no interface, 8x the cost, cluster-scale); a uniform
> [0,32] at 64^3 (no interface, same cost, but half the resolution on the star); or stabilising
> the interface (more Kreiss-Oliger dissipation, stronger shift damping, or moving it outward).
>
> **In progress:** the uniform 64^3 control, `p1b_migration_font_box32_uni64`, started
> 2026-09-24 15:10, to t = 300. It answers whether the interface was the culprit (does it
> survive past t = 221.5?) and whether the outflow diagnosis holds with no interface at all (is
> mass conserved?). It cannot give the period: at dx = 0.5 its initial central density is
> already 11% below exact (7.099e-3 vs 7.993e-3, against 3% at dx = 0.25), so its resolution
> differs from both runs it is compared with. Result to be recorded here.
>
> **Phase-1 status: the VERDICT half passes on both Font setups; the QUANTITATIVE half fails.**
> Nothing in this phase is converged to the few percent the gate asks for, and the apparent
> exceptions (the coarse collapse run's 4% proper-time agreement) did not survive refinement.

> **Result 2026-09-23 — the Font resolution pair does NOT converge; this gate FAILS.**
> Pair: the 64^3 octant run (dx_centre = 0.25) and the same run with one static level on the
> centre (dx_centre = 0.125), identical star and kick. The observable is the central proper time
> tau = int(alpha_min dt) at which the lapse collapses:
>
> | marker | tau coarse | tau refined | shift |
> |---|---|---|---|
> | alpha_min < 0.20 | 7.157 | 6.330 | -11.6% |
> | alpha_min < 0.15 | 8.299 | 7.385 | -11.0% |
> | alpha_min < 0.10 | 8.976 | 7.972 | -11.2% |
> | alpha_min < 0.05 | 10.429 | 8.393 | -19.5% |
>
> An 11% shift, consistent in sign across three markers and worse at the fourth, is not a
> converged quantity, so no Richardson extrapolation is meaningful from two points. Peak
> rho/rho_c(0) moves 8.74 -> 13.72 and min alpha 0.0456 -> 0.0222, neither converged either.
>
> **The coarse run's 4% agreement with DynGR1D (10.43 vs 10.01) was therefore fortuitous.**
> Refining moves AthenaK AWAY from the 1+1D engine, to 8.39, i.e. 16% below it. The verdicts
> still agree — both codes collapse this star — but no quantitative observable of the collapse is
> converged in AthenaK at these resolutions, so the quantitative half of the Phase-1 gate ("AH
> mass to a few %") cannot be claimed.
>
> The breakdown itself is resolution-INDEPENDENT: the 3-metric loses positive-definiteness at the
> innermost cells at t = 65.75 (coarse) and t = 64.50 (refined). That the failure time barely
> moves while the physics before it does not converge points at something structural in evolving
> the forming puncture on this octant setup, not at cell size. Before this gate can be retried:
> a horizon finder (needs a full box, not an octant), and either excision/puncture handling or
> more levels concentrated at the centre rather than one level over the whole star.

> **Result 2026-09-23 — the AH-mass gate is out of reach on this machine, and the same number
> explains the breakdowns.** The horizon of the M_AH = 1.2734 remnant sits at isotropic radius
> M/2 = 0.637. Cells across that radius:
>
> | central dx | 0.5 | 0.25 (the octant) | 0.125 | 0.0625 | 0.03125 |
> |---|---|---|---|---|---|
> | cells | 1.3 | **2.5** | 5.1 | 10.2 | 20.4 |
>
> So every run in this phase has been resolving the forming horizon with about two and a half
> cells. That is why they lose positive-definiteness of the 3-metric at the innermost cells just
> as the lapse collapses, and it is why ONE extra level did not help: it took the horizon from
> 2.5 cells to 5.1, still hopeless, and the breakdown time duly barely moved (65.75 -> 64.50).
> The breakdown and the missing horizon mass are the same defect, and both need 3-4 more levels,
> not one.
>
> Two further constraints. `src/z4c/fastflow.cpp` has no symmetry or domain handling, so 7/8 of
> the horizon surface lies outside an octant: the finder needs a FULL box, 8x the volume. And
> AthenaK does not subcycle, so every level halves the global timestep.
>
> **DEFERRED TO THE CLUSTER, with the run configured and smoke-tested:
> `code/athenak-bdnk/inputs/p1b_collapse_font_ahf_cluster.athinput`.** Full box [-16,16], root
> 64^3 with 16^3 blocks, four nested centre levels reaching dx = 0.03125 (20 cells across the
> horizon radius), diode inner faces instead of reflecting, `<fastflow>` searching once per
> M_sun from t = 30, to t = 80.
>
> Validated locally to 2 cycles: 288 MeshBlocks and 4 physical levels as designed, TOV data
> identical to the octant runs (R = 5.8381, M = 1.44759), and fastflow initialising and writing
> `tov.horizon_summary_0.txt`. One configuration bug was caught and fixed in the process:
> `use_puncture_0 = 0` copied from the onepuncture example means "track puncture index 0" and
> aborts with `punc = 0 > npunct = 0`; with no punctures it must be -1.
>
> Cost measured on that smoke test rather than estimated: 2.07e4 zone-cycles/s on 2 threads,
> and t = 80 needs 10240 cycles x 1.18M cells = **1.21e10 zone-cycles — 6.8 days on 2 threads,
> 1.7 to 2.3 days on 8.** That is risk item 6 arriving with a number.
>
> So: the verdict half of Phase 1 is met on both Font setups and stands; the AH-mass half is
> deferred, ready to submit, and is not attempted on this machine.

3. Resolution: 64³ and 128³ octant, plus one SMR level on the star; Richardson on the AH
   formation proper time. Record the cost table (this sets the size of Phase 3).

### Phase 2 — BDNK in AthenaK's `dyn_grmhd` (≈ 4–6 weeks; the real work)

**Formulation** (the one already fixed in `progress/understanding/SYNTHESIS.md` §1–2 and
implemented in `Recovery.jl` / `ConformalBDNK.jl` / `shum_core.jl`): T^{ab} = (ε+A)u^a u^b +
(P+Π)Δ^{ab} + Q^a u^b + Q^b u^a − 2ησ^{ab}, barotropic trunk (β_ε = τ_Q c_s², β_n = 0), Shum
frame (ŝ, â, q̂) = (1, 1, 0.999) with (η̂, ζ̂), τ_ε = V̂L etc. The corrections contain **time
derivatives of the primitives**, so the conserved variables are not an algebraic function of
the primitives: recovery is the PMP gradient-frozen linear solve (SYNTHESIS §2.2), which for a
barotrope is a 2×2 (radial) / 4×4 (3D) linear system, closed-form conformal, numerical
Jacobian for a general EOS (the highest-risk trunk item — do barotropic first).

**Where it enters AthenaK** (a new module `src/bdnk/` next to `dyn_grmhd/`, its own task
list, reusable — the CCM rule of `progress/CLAUDE.md`):

| piece | AthenaK hook | what changes |
|---|---|---|
| spatial gradients of (ln ε, u^i) | new task `BDNK_Gradients` after `MHD_RecvU` | 4th-order FD like Z4c's derivatives; stored per stage |
| time derivatives (ε̇, u̇^i) | new task after cons2prim | the linear solve from (q − q_PF(p, frozen ∂p)); first stage of each step uses the previous step's values (Chabanov–Rezzolla-style backward difference) as the guess; PF fallback below Δ_visc (`{{TOL}}`: off for benchmarks) |
| conserved variables and fluxes | `CalcFluxes` / rsolvers (`hlle_dyn_grmhd.hpp`) | face states carry A, Π, Q^i, σ^{ij}; HLL wave-speed estimates use the BDNK characteristic speeds (c₊ = √3 c_s in the production frame — bigger than the acoustic speed, so both the Riemann solver and `newdt` must know them) |
| Einstein-equation source | `SetTmunu` (`z4c/tmunu.cpp`) | the full T^{ab} including corrections, i.e. the reconstructed time derivatives feed the Z4c matter terms — this is the "GR-coupled BDNK" content and does not exist anywhere else in the toolkit |
| geometric sources | `AddCoordTerms` | T^{ab}∂γ/∂α terms with the corrected T^{ab} |
| floors / limiting | `reset_floor` error policy | keep p + Π ≥ 0.1 p (Chabanov–Rezzolla σ = −0.9 rule), |A| ≤ ε; a per-cell causality flag (port of `Causality.jl`'s biquadratic: real, non-negative, subluminal speeds) and the reduction-constraint monitor (Fantini–Rubio: the first-order constraints of the auxiliary variables grow with no damping scheme available — record, do not hide) |
| time step | `newdt` | explicit SSP-RK3 with the BDNK speeds; IMEX only if τ_ε < few Δt (Shum's frame choice L ~ R keeps it explicit; check the ratio in the collapse phase, where gradients steepen and L should follow the local scale) |

**Verification ladder for the module** (each an existing number in
`progress/reproduction/LEDGER.md` or in the Julia tests, so no new targets are invented):

1. flat space, `<adm>` Minkowski (`minkowski=true` in the TOV pgen or a slab): the steady
   shock (Pandya/PMP: e_L=1, v_L=0.8 → e_R = 4.4074, v_R = 0.41667, η₀ = 0.2) and the smooth
   Gaussian of `ConformalEvolution` against `1D_conformal_bdnk` — the RH half of the
   `s2.gr_coupling` gate; run as a 3D slab (nx2 = nx3 = 4);
2. linear sound-wave damping vs Kovtun (Im ω = −k²Γ_s/2 at low k) — the `s1a`-level check of
   the shear/bulk normalization in 3D;
3. Cowling star (`<adm>` only) with the Shum frame vs our R5 reproduction: F = 2.69 kHz,
   H₁ ≈ 4.55, decay 1/τ in the 0.0011–0.0021 M⊙⁻¹ band (`shum_qnm_production.jl`, VALIDATION
   §3) and vs `SphBDNK`'s frame-independence result (γ_visc invariant to 0.045% across the
   three Clarisse frames);
4. dynamical spacetime, stable star, η̂ = ζ̂ = 0: byte-level agreement with Phase 0/1 (the
   module must be a no-op at zero dissipation — a regression test).

### Phase 3 — the BDNK collapse matrix (≈ 2–3 weeks of runs)

Setups (a) and (b) of Phase 1, with:

- η̂, ζ̂ ∈ {0, 0.02, 0.05, 0.1} (the Shum stable window τ_ε = (4/3)η̂ + ζ̂ ≲ 0.1), bulk-only and
  shear-only rows (bulk is what the IS literature says matters for collapse);
- two frames at fixed (η, ζ) — the production frame and one more from `{{FRAME_SET}}` — so
  that the observables below are shown frame-independent (extends `test_frame_independence`
  to the nonlinear GR regime; a frame-dependent AH time would mean the module is wrong);
- resolutions 64³ / 128³ octant + SMR, one Richardson triple per physics point.

Observables: ρ_c(τ_c), α_c, AH formation proper time and mass, M_b outside the horizon, the
BDNK entropy production ∫(2ησ² + ζθ²)/T (the `test_entropy` machinery), max |A|/ε and
|Π|/p (BDNK's Reynolds-number analogues), the causality-flag and reduction-constraint
histories, Z4c Hamiltonian/Θ norms.

Physics questions the matrix answers (the ones `s2.gr_coupling` → `s2.is_contrast` were
created for): does BDNK bulk viscosity delay or prevent the collapse the way IS does in
Chabanov–Rezzolla (their high-ζ migration does not collapse and the shock stalls at ≈ 60 km);
how the AH time scales with ζ̂ and η̂; whether any point in the window trips the causality flag
before the horizon forms (gradients steepen, the frame relaxation scale does not).

### Phase 4 — cross-code closure and write-up (≈ 1 week)

- Extend `DynGR1D.jl` with the same constitutive terms (its RHS is 400 lines; the recovery is
  `cons2prim_bdnk_barotropic` + the linear solve of `shum_core.jl`) so the genuine 1+1D
  engine runs the same matrix on the laptop in minutes. The AthenaK ↔ DynGR1D comparison of
  the AH proper time and mass per physics point is the closure of `s2.gr_coupling`.
- VALIDATION.md section (§8, "AthenaK"), DAG update (`s2.gr_coupling` → three nodes:
  `s2a.athenak_baseline`, `s2b.bdnk_module`, `s2c.collapse_matrix`; `s2.is_contrast` then has
  its numbers to contrast against), LEDGER rows, figures (VLM-inspected per the project rule).

## 2. Risks, and what is decided in advance

1. **Cartesian vs spherical**: the octant grid breaks spherical symmetry at the (Δx)² level and
   the atmosphere/floors act on the surface; compare invariants only (AH mass, proper times,
   M_b), never coordinate profiles. The same discipline as VALIDATION §7.10–7.13.
2. **The general-EOS Jacobian** of the BDNK recovery is not in any source in closed form.
   Barotropic Γ=2 first; a tabulated/hybrid EOS only after the barotropic matrix is done.
3. **Stiffness near the horizon**: τ_ε ∝ L with L fixed by the initial star; as the collapse
   steepens gradients the explicit step is set by the BDNK speeds, not by τ — measure the
   ratio Δt/τ_ε through the run; switch to IMEX (`s1c.imex` machinery) only if it exceeds ~1.
4. **Reduction-constraint growth** (Fantini–Rubio 2506.06430): no damping scheme exists; the
   run records the constraint norm and the analysis reports at what time it becomes O(1)
   relative to the fields. If that happens before the horizon forms the result is "not
   trustworthy after t_*", not a fabricated collapse.
5. **AthenaK's FOFC and `reset_floor`** were written for an ideal fluid: a large negative Π
   looks like a negative pressure to them. The Chabanov–Rezzolla limiter (Π ≥ −0.9p) is applied
   inside the recovery, before the floors see the state, and every limiter hit is counted.
6. **Cost**: 128³ octant with Z4c is a GPU job (hours per run); the laptop does 64³ development
   runs and every Julia cross-check. Phase 0 measures the real numbers before the matrix is
   sized.

## 3. Immediate next actions

1. Set up WSL2 + build AthenaK (`dyngr_tov`), run `whisky_tov` as shipped, log cost and the
   Cowling spectrum → `code/athenak-bdnk/README.md` first entry.
2. Fork AthenaK to `FoAKTEE/athenak`, branch `bdnk`; add the Z4c-on TOV input and the two
   collapse inputs (`tov_collapse_eps0003.athinput`, `tov_migration_font2002.athinput`).
3. Write `code/athenak-bdnk/analysis/ah_compare.py` (AH time/mass vs the `DynGR1D` record) so
   Phase 1's gate is scripted before the runs exist.
4. Start `src/bdnk/` with the barotropic recovery ported from `Recovery.jl` + `shum_core.jl`
   and its unit test against the Julia round-trip numbers (STEP-0 gate: ≤ 10⁻¹⁰).
