# Source: jl-NeutronStarOscillations (Keeble & Redondo-Yuste 2026)

Local Julia package: `/data/haiyangw/claude/BDNK/ref-code/NeutronStarOscillations/lennoxkeeble-NeutronStarOscillations.jl-786dbc9`

Julia code for **linearized radial perturbations of spherically-symmetric neutron stars**, in the
frequency and time domains, for **perfect fluid (PF)**, **Eckart**, and **BDNK** first-order viscous
hydrodynamics. Cold polytropic EOS only: `p = κ ε^(1+1/n)`. This is the single most build-critical
in-language reference for the Julia BDNK NS reproduction.

Tested on Julia 1.12.0 (compat: julia ≥ 1.10.10). Deps: DifferentialEquations, NonlinearSolve, Arpack,
Dierckx (splines), HDF5, StaticArrays, SparseArrays, FFTW, CairoMakie.

---

## 0. Units & constants (`src/Constants.jl`)

Geometrized-CGS hybrid: everything internal is in **km** (length) and **km^-2** (energy density,
pressure). Key conversions:
- `c = 299792458` (m/s); `sec_to_km = c*1e-3`; `kHz_to_km = 1e3/sec_to_km`
- `G = 6.6743015e-11`; `M_sun = 1.988416e30` kg; `kg_to_m = G/c^2`
- `gram_per_cm3_to_km_minus2 = kg_to_m*1e-3/kg_to_g/cm_to_km^3` — converts ε [g/cm^3] → [km^-2]
- `Msun_to_km = M_sun*kg_to_m*1e-3`
- `MeV_per_fm3_to_km_minus2`, `dyne_per_cm2_to_km_minus2`, `MeV_to_km`, etc.
Frequencies output as `ω²/(2π·kHz_to_km)²` → kHz².

## 1. Star object / EOS (`src/NeutronStarOscillations.jl`)

`mutable struct Star` (l.10) holds εc_SI [km^-2], pc_SI, κ, n, EOS functions, viscous params
(η, ζ, τε, τP, τQ, L), numerics, ID functions.

EOS closures built in `Star(...)` constructor (l.127, l.164-171):
```julia
γ = 1.0 + 1.0/n
EOS_ε(p)  = (p/κ)^(1/γ)          # ε(p)
EOS_p(ε)  = κ * ε^γ              # p(ε)
dp_dε(ε)  = γ*κ*ε^(γ-1)          # = cs²  (sound speed squared)
d2p_dε2(ε)= γ*(γ-1)*κ*ε^(γ-2)
εc_SI = εc_cgs * gram_per_cm3_to_km_minus2
pc_SI = EOS_p(εc_SI)
```
**Sound speed: `cs = sqrt(dp_dε(ε)) = sqrt(γκ ε^(γ-1))`.** No separate primitive recovery (analytic
polytrope, ε↔p invertible). `cs_prime` along the star uses chain rule:
`cs_prime = d2p_dε2(ε)·p'(r)/(2 cs³)` (`cs_prime_func`, FrequencyDomain l.165, BDNKInitialData l.15)
with TOV `p'(r) = -0.5*((-1 + r/(r-2m) + 8πr³p/(r-2m))*(p+ε))/r`.

`ptol < 0` ⇒ termination pressure = `|ptol|·pc_SI` (fractional). Example uses `ptol = -1e-6`,
`ptol_TD = -1e-3`.

## 2. TOV background (`src/TOV.jl`)

Two solvers. Units: m,r in km; p,ε in km^-2. Equations (geometrized, `c=G=1`, factor `4π`):
```
m'  = 4π r² ε
ν'  = (2m + 8π r³ p)/(r² - 2 r m)          # metric: ds² = -e^ν dt² + e^λ dr² + ...
p'  = -((m + 4π r³ p)(p+ε))/(r(r-2m))
```
(`TOV.Explicit.m_prime/nu0_prime/p_prime`, l.187-189). ε from EOS each step.

- **Explicit RK4** (`TOV.Explicit.solve`, l.201): classic RK4 in r, step h. Center regularity enforced
  by `regularity`/`shift` flags only at first step. Stops when `p ≤ ptol`. ν shifted at end to match
  Schwarzschild exterior: `ν += -ν[end] + log(1 - 2m[end]/r[end])`.
- **Implicit Crank–Nicholson** (`TOV.Implicit.solve`, l.442): 2nd-order, Newton iteration each step
  (`F1..F4` residuals l.358-364, analytic Jacobian-inverse `A11..A44` l.366-396, `newton_iterate!`).
  Vars u=[m,p,ε,ν]. `TOV_iter_tol=1e-15`, `TOV_max_iter=10`, `TOV_initial_r=1e-15`.

Frequency/time-domain code uses **Explicit RK4** at `h_TOV=1e-4` km then spline-interpolates
(`Dierckx.Spline1D`, order 5, s=0, bc="error") onto the computational grid.

## 3. BDNK characteristic speeds & causality

### Characteristic speeds (`src/BDNKCharacteristicSpeeds.jl`)
Radial char. speeds c± solve `Λ2 c⁴ - 2Λ1 c² + Λ0 = 0`:
```
c²± = (Λ1 ± sqrt(Λ1² - Λ0)) / Λ2
Λ0 = (4 L⁴(3ζ+4η)⁴ (τP-1) τQ² τε (p+ε)² cs⁴)/81
Λ1 = (L²(3ζ+4η)²(τε + τQ(τP+τε))(p+ε) cs²)/9
Λ2 = (2 L²(3ζ+4η)² τQ τε (p+ε))/9
```
(`BDNKCharacteristicSpeeds.compute`, `Λ0/Λ1/Λ2`, `c_plus_sq/c_minus_sq`). cs²=dp/dε.

### Causality check (`NeutronStarOscillations.jl` l.869-888, `check_causal_frame`)
Upper bounds on cs²(center) for causal frame:
- τP < 1: **always violated** (need τP ≥ 1).
- τP == 1: `cs²_max1 = τQ τε/(τQ + τε + τQ τε)` (`cs2_max1`).
- τP > 1: `cs²_max2 = (τε + τQ(τP+τε) - sqrt(-4(τP-1)τQ²τε + (τε+τQ(τP+τε))²))/(2(τP-1)τQ)` (`cs2_max2`).
Checks `dp_dε(εc_SI) ≤ cs²_max`. With example (τε=15,τP=1.5,τQ=20): defines a max central sound speed.

## 4. BDNK linearized perturbation system

### Variable set (full, non-Cowling, time domain)
Perturbations (functions of t,r): **u1 = δu** (radial 4-velocity pert.), **u2 = δε** (energy density
pert.), **u3 = δλ** (metric pert., grr). Auxiliary first-order vars: **w_i = ∂_r u_i**, **v_i = ∂_t u_i**
(v1=∂t δu, v2=∂t δε, v3=∂t δλ). v3 is eliminated algebraically.

The system = 2 wave-like equations for (δu, δε) + 1 **constraint** equation for δλ (no time deriv. in
first-order form). Authors solve the **constrained** form (more stable) and monitor independent residuals
of both constrained and original wave-like systems.

### PDE coefficients (`src/BDNKDiffEqCoefficients.jl`, `struct CNSystem`)
At each radial grid point the discretized equations use coefficient arrays evaluated on the TOV
background: A1..A14 (Eq1, the δλ/EFE constraint), B1..B14 (Eq2), C1..C8 (Eq3), plus center coeffs
D1..D4, E1..E6, F1..F3, G1..G3. Each is an explicit closed-form function of
`(m,p,ε,ν,cs,cs',cs'',r,η,ζ,τε,τP,τQ,L)`. They are long but fully analytic — directly portable.

Representative simple ones:
- `D1 = -8Lπ(3ζ+4η)τε(p+ε)/(3 e^{ν/2})`, `D2 = -8π/3`, `D3 = 1`, `D4 = -8Lπ(3ζ+4η)τε/(9 e^{ν/2})`.
- `E3 = 2L(3ζ+4η)τQ cs²`, `E6 = L(3ζ+4η)τε/(3 e^ν)`,
  `E4 = L(3ζ+4η)(τQ+τε)(p+ε)/e^ν`.
- `G1 = -8`, `G3 = -2/(p+ε)`, `G2 = 3 e^{ν/2}/((3Lπζτε+4Lπητε)(p+ε))`.
- `C1` (Eq3, the δλ constraint) is moderate length — see file l.286.
- A1, B1 are very long (full EFE / momentum constraint). All keyed off `(3ζ+4η)`, `(p+ε)`, `cs²`,
  `e^{ν/2}` and TOV combinations `(m+4πr³p)`, `(r-2m)`.

### Discretization of the equations (`src/BDNKLinearSystem.jl`)
Crank–Nicholson in time (couples levels n, n+1 at half step n+1/2; functions averaged
`(np1+n)/2`, time deriv `(np1-n)/k`). Spatial derivatives 2nd-order central interior, one-sided at
surface. Per interior point j, three nontrivial equations:
- `InteriorEqs.Eq1` (A-coeffs; δλ/EFE), `Eq2` (B-coeffs; momentum/δu wave), `Eq3` (C-coeffs; δε wave).
- Plus 2 "trivial" time-reduction CN eqs: `u1np1 - u1n = k(v1np1+v1n)/2`, same for u2.
Center point (`module Center`) has special `Eq` using forward differences + the E-coeffs that constrains
(u2, v2) at r=0. Surface (`SurfaceEqs`) uses backward differences (j-3..j stencils).

`v3` (=∂t δλ) is reconstructed from the others via `Forward_v3/Interior_v3/Backward_v3`
(`compute_v3!`, l.10-31) — algebraic, no extra evolution.

### Linear solve (time stepping, `src/TimeDomain.jl` l.1111-1158)
System is **linear** ⇒ Jacobian assembled ONCE:
```
compute_jacobian!(jacobian, CNSystem, h, k, nPointsSpace)   # block layout, l.116
sparse_jac = sparse(jacobian); jac_fact = factorize(sparse_jac)   # LU once
# each step:
compute_RHS!(linear_system_RHS, ...dummy np1..., u..n, ..., h, k)
ldiv!(sol, jac_fact, linear_system_RHS)                     # reuse factorization
```
Unknown vector ordering: u1(j=2..N), u2(j=1..N), u3(j=2..N), v1(j=2..N), v2(j=1..N); size `5N-3`
(u1,u3,v1 omitted at center j=1 by regularity). Solution unpacked to u1_np1[2:end], u2_np1, u3_np1[2:end],
v1_np1[2:end], v2_np1. Then KO dissipation applied, v3 recomputed, independent residuals computed.

## 5. Time-domain Cowling scheme (`src/CowlingTimeDomain.jl`, `src/CowlingTDEqs.jl`)

Cowling (fixed metric) BDNK uses **RK4** explicit (not CN) — `CowlingTimeDomain.solve` l.1092 has
k1..k4 arrays. 6 perturbation vars carried: u1=δu, u2=δε, w1=∂rδu, w2=∂rδε, v1=∂tδu, v2=∂tδε.
RHS time-derivative equations:
- `CowlingTDEqs.RK4.BDNK.dv1_dt` (l.19) — the δu evolution, very long, `/(4Lr³(3ζ+4η)τQ(p+ε)cs⁴)`.
- `CowlingTDEqs.RK4.BDNK.dv2_dt` (l.21) — the δε evolution, `/(4r³τε)`, contains `cs_prime_prime`.
- `du1_dt=v1, du2_dt=v2, dw1_dt=∂r v1, dw2_dt=∂r v2` (trivial reductions).
- Surface variants `RK4Surface.BDNK.*` use backward stencils.
- Center: `CowlingTDEqs.CenterEqs.BDNK.compute_all` (l.49) with forward FD (`FD_1/2/3_deriv`); enforces
  `dt_u1_0=0, dt_w2_0=0, dt_v1_0=0, dt_u2_0=v2, dt_w1_0=∂r v1`, and an explicit `dt_v2_0(...)` (l.79).

Time step: `k = CFL·h` (l.1103). Background TOV at `h_TOV=1e-4`, spline-interpolated, then **downsampled
to grid spacing h**. Λ0,Λ1,Λ2 computed per point and saved.

PF Cowling pulsation ODE in second-order form is in `FrequencyDomain.PerfectFluidCowling.Matrix.A0/A1/A2`
(see §6); the PF/Eckart non-Cowling time domain (`TimeDomain.jl`) integrates ξ via RK4 (`RK4Eq`,
`RK4EqSurface` l.390-391).

### Kreiss–Oliger dissipation (`kreiss_oliger`, CowlingTimeDomain l.1056 / TimeDomain l.529)
5-point: `coef*(u_{j+2} - 4u_{j+1} + 6u_j - 4u_{j-1} + u_{j-2})/16`, subtracted from advanced step,
applied only on interior `3:N-2`. Example `KO = 0.2`.

## 6. Frequency domain (`src/FrequencyDomain.jl`) — PF & Eckart only (no BDNK FD)

Strategy (docstring l.4-6): **matrix method** to locate eigenvalues (good initial guess) → **shooting
method** to refine, driving boundary-condition residual (Lagrangian pressure pert. = 0 at surface) to 0.

### Cowling PF pulsation equation (`PerfectFluidCowling.Matrix`, l.148-151)
Write as `A2(r) ξ'' + A1(r) ξ' + A0(r) ξ = ω² ξ`:
```
A2 = -(e^ν (r-2m) cs²)/r
A1 = e^ν(4πr³p - 2r cs((1+2πr²ε)cs + r cs') + m(1+5cs²+4r cs cs'))/r²
A0 = 2e^ν( m²(1+5cs(cs-2r cs')) + r m(-1-8πr²p + cs((-5+8πr²ε)cs + r(9+8πr²ε)cs'))
          + r²(cs² + 2r(πr(p+ε+8πr²p ε - 8πr²ε²cs²) - (1+2πr²ε)cs cs')) )/(r³(r-2m))
```
ξ(0)=0 (regularity, center row dropped). cs=sqrt(dp_dε(ε)).

### Matrix method (`fill_matrix!` l.291, `get_eigensystem` l.191)
Discretize at interior nodes with central differences:
`Ai_minus1 ξ_{i-1} + Aii ξ_i + Ai_plus1 ξ_{i+1} = ω² ξ_i` where (note A3≡A2 stencil var naming):
```
Ai_minus1 = (2A3 - A2 h)/(2h²),  Aii = A1 - 2A3/h²,  Ai_plus1 = (2A3 + A2 h)/(2h²)
```
Surface row enforces Δp=0 (`An_minus1`, `Ann` l.160-161). Build matrix A, `eigen(A)`, sort, take lowest
`N_eigvals`; eigenvalues returned as `ω²/(2π kHz_to_km)²` (kHz²). Eigenvectors prepended with ξ(0)=0 and
normalized so `|real(ξ(R))| = 1 km` (`normalize_eigvecs!` l.128).

### Shooting method (`Shoot`, l.330-420)
For each matrix eigenvalue ω_init, root-find ω so the surface BC `ξ'_BC(ω)=0`:
- `NonlinearSolve.RobustMultiNewton()`, `AbsTerminationMode()`, reltol/abstol `1e-12`, maxiter 200.
- Inner integration: RK4 of the 1st-order system `ξ'=X`, `X' = (A0 ξ + A1 X - ω² ξ)/A2`
  (`X_prime` l.447, `integrate`/`fast_integrate`), from center (shift 1e-18) to surface.
- Returns refined ω (→ kHz²), eigenvector, residual, retcode; saved to HDF5.

`N_eigvals = 5` (example). mode 0 = fundamental, 1 = first overtone, etc.

## 7. BDNK initial data (`src/BDNKInitialData.jl`)

User specifies δu, δε and their r- and t-derivatives at t=0 (`δu_ID, δu_dr_ID, δu_dt_ID, δε_ID,
δε_dr_ID, δε_dt_ID`). The metric perturbation **δλ (u3) is solved from a lower-order constraint ODE**:
`du3_dr(u1,u2,u3,w1,w2,v1,v2,m,p,ε,ν,cs,cs',r,η,ζ,τε,τP,τQ,L)` (l.18) integrated by RK4 from r=0 with
δλ(0)=0 (`integrate` l.67, `compute_initial_data` l.20). RK4 uses background at half-steps ⇒ pulsation
solved at steps 2h while TOV is at h. Regularity at center enforced by zeroing k1 at first step.

Built-in ID options (params.jl): Gaussian `A/exp((r-r0)²/w²)` (amplitude/center/width), or PF/Eckart
eigenvector as ID. Regularity requires `δu(0)=0` (checked, error if `|δu(0)|>1e-16`).

## 8. Independent residuals (`src/BDNKIndependentResidual.jl`)

Distinct (leapfrog n-1,n,n+1) re-discretization of the EOM to verify the CN solver converges to the
true PDEs. Two IRs for constrained system (δu, δε wave eqs), three for original wave-like system.
`ConstrainedSystem.{Center,Interior,Surface}.Eq1/Eq2`, `WaveSystem.*`. Uses the same A/B/C/E CNSystem
coefficients. Names exported: `["Constrained_IR1","Constained_IR2","Wave_IR1","Wave_IR2","Wave_IR3"]`.

## 9. Concrete example parameters (`src/Examples/params.jl`)

Stellar: `eps_central = 5.5e15` g/cm³ (precompile workload uses 3.0e15, n=0.8, κ=700);
**n=1.0, κ=100.0** [km^(-2/n)] (⇒ γ=2). Viscous: **η=1e-2, ζ=1e-2** (dimensionless), relaxation
**τε=15.0, τP=1.5, τQ=20.0** (dimensionless), **L=1.0 km** (viscous length scale).
Numerics: `ptol=-1e-6`, `ptol_TD=-1e-3`, `NL_reltol=NL_abstol=1e-12`, `NL_maxiter=200`, `N_eigvals=5`,
`TOV_iter_tol=1e-15`, `TOV_max_iter=10`, `TOV_initial_r=1e-15`, `TOV_max_steps=1e11`.
Time domain: `KO=0.2`, `CFL=0.1`, `h_save=4e-2` km, `total_time_ms=0.05`, `dt_save_ms=total/100`,
`save_every=50`. Gaussian ID: amplitude 1.0, center 5.0 km, width 0.5 km.
PF star ⇐ all viscous params 0; Eckart star ⇐ τε=τP=τQ=0 (η,ζ kept); BDNK ⇐ all set.

## Build pointers (file:function)

- EOS / Star / causality: `NeutronStarOscillations.jl:Star` (l.127), `:check_causal_frame`,
  `:cs2_max1/cs2_max2`.
- TOV: `TOV.jl:Explicit.solve` (RK4), `:Implicit.solve` (CN+Newton), `:Explicit.{m_prime,nu0_prime,p_prime}`.
- Char. speeds: `BDNKCharacteristicSpeeds.jl:compute,Λ0,Λ1,Λ2,c_plus_sq,c_minus_sq`.
- BDNK PDE coeffs: `BDNKDiffEqCoefficients.jl:CNSystem` (A1..A14,B1..B14,C1..C8,D1..D4,E1..E6,F1..F3,G1..G3).
- BDNK CN discretization: `BDNKLinearSystem.jl:InteriorEqs.{Eq1,Eq2,Eq3}`, `SurfaceEqs.*`, `Center.Eq`,
  `Jacobian.compute_jacobian!`, `RHS.compute_RHS!`, `compute_v3!`.
- BDNK time stepping: `TimeDomain.jl:solve` (full, l.~940-1180; LU-once + ldiv!).
- BDNK Cowling RHS: `CowlingTDEqs.jl:RK4.BDNK.{dv1_dt,dv2_dt}`, `RK4Surface.BDNK.*`,
  `CenterEqs.BDNK.compute_all`; driver `CowlingTimeDomain.jl:solve` (l.1092).
- Frequency domain: `FrequencyDomain.jl:PerfectFluidCowling.Matrix.{A0,A1,A2,fill_matrix!,get_eigensystem}`,
  `.Shoot.{iterate_frequencies,integrate,X_prime}`.
- Initial data: `BDNKInitialData.jl:compute_initial_data,integrate,du3_dr`.
- Independent residuals: `BDNKIndependentResidual.jl:ConstrainedSystem.*,WaveSystem.*`.
- KO dissipation: `*TimeDomain.jl:kreiss_oliger`.

## Open questions / gaps
- No BDNK frequency-domain solver (FD is PF + Eckart only); BDNK is time-domain only.
- arXiv id not yet published (README placeholder "Keeble:2026"); no benchmark frequency table in repo.
- General BDNK stress-tensor / constitutive relations are NOT written out in covariant form in the code;
  only the already-linearized, already-radially-reduced coefficient functions appear. The mapping from
  the general BDNK first-order transport coefficients to (η,ζ,τε,τP,τQ,L) lives in the (unpublished) paper.
- d2p_dε2 at center handled via cs_prime[1]=0 and FD for cs'' (`FiniteDiffOrder4`).
