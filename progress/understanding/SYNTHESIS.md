# BDNK Neutron-Star Code — Build Specification & Project DAG

**Lead-architect synthesis** reconciling Kovtun (1907.08191), Bemfica-Disconzi-Noronha (2009.11388),
Pandya-Most-Pretorius (2209.09265 / 2201.12317), the conformal C reference
(`ref-code/1D_conformal_bdnk`), the Julia `NeutronStarOscillations.jl` package
(Keeble & Redondo-Yuste), the Shum-et-al nonlinear Cowling paper (2509.15303), the
Caballero-Yunes radial-mode paper (2506.09149), the Redondo-Yuste / Bussières axial
papers (2411.16841, 2604.13208), the frame-constraints notebook, and the IS baseline
(Chabanov-Rezzolla 2311.13027).

Conventions: **mostly-plus signature `(-,+,+,+)`**, `u^a u_a = -1`,
geometrized units `c = G = 1`, `8π` in Einstein/TOV. Greek/latin indices interchangeable.

---

## 1. Canonical BDNK formalism (reconciled notation)

### 1.1 Stress tensor and currents

Kovtun's hydrodynamic-frame decomposition (1907.08191 eq 2.2), specialized to the
Bemfica-Disconzi-Noronha **working frame** with `J^a = n u^a` (no gradient correction to
the particle current — adopted by every numerical source here):

```
T^{ab} = (ε + A) u^a u^b + (P + Π) Δ^{ab} + Q^a u^b + Q^b u^a − 2 η σ^{ab}
J^a    = n u^a
Δ^{ab} = g^{ab} + u^a u^b           (spatial projector, Δ^{ab} u_b = 0)
ρ      ≡ ε + P                       (enthalpy density)
c_s²   ≡ (∂P/∂ε)_s̄                   (adiabatic sound speed²; = γP/ρ for ideal gas)
```

`ε` is the equilibrium energy density (the **defining scalar of the frame**, NOT `E = u_a u_b T^{ab}`),
`P = P(ε,n)` the equilibrium pressure, `n` baryon number, `u^a` the (frame) 4-velocity.

### 1.2 First-order constitutive corrections

This is the single canonical set; every source maps onto it (cross-walk in §1.4):

```
A   = τ_ε [ u^c ∇_c ε + ρ ∇_c u^c ]                                          (E-correction)
Π   = −ζ ∇_c u^c + τ_P [ u^c ∇_c ε + ρ ∇_c u^c ]                             (bulk-like)
Q^a = τ_Q ρ u^c ∇_c u^a + β_ε Δ^{ac} ∇_c ε + β_n Δ^{ac} ∇_c n               (heat flux)
σ^{ab} = (1/2)( Δ^{ac}Δ^{bd} + Δ^{ad}Δ^{bc} − (2/3) Δ^{ab}Δ^{cd} ) ∇_c u_d   (TT shear)
```

The 8 frame coefficients: `{ τ_ε, τ_P, τ_Q, η, ζ, σ, β_ε, β_n }`. They are *first-order*
("BDNK frame") coefficients — NOT Israel-Stewart second-order relaxation times; there is no
extra evolution equation for `A, Π, Q`.

Heat-flux coefficients fixed by EOS + conductivity `σ` (BDN 2009.11388 eq 1062;
PMP 2209.09265 eq 18-19):

```
β_ε = τ_Q (∂P/∂ε)_n + σ (T ρ / n) (∂(μ/T)/∂ε)_n
β_n = τ_Q (∂P/∂n)_ε + σ (T ρ / n) (∂(μ/T)/∂n)_ε
```

For a **barotropic / single-fluid** treatment (the NS trunk; Redondo-Yuste 2411.16841,
Bussières 2604.13208, Shum 2509.15303): drop `n` and `σ`, so `Q^a = τ_Q [ ρ u^c∇_c u^a + c_s² Δ^{ac} ∇_c ε ]`,
i.e. `β_ε = τ_Q c_s²`, `β_n = 0`. **This is the formalism Stages 1A/1B/1C actually use.**

### 1.3 Ideal-gas ("gamma-law") microphysics (PMP 2209.09265 eq 25-37)

Used at STEP 0 for general-EOS unit tests and Stage-2 collapse:

```
P(ε,n) = (Γ−1) m n e,   ε = m n (1+e),   e = specific internal energy,   Γ ∈ (1,2)
c_s² = ΓP/ρ
(∂P/∂ε)_n = Γ−1,   (∂P/∂n)_ε = −(Γ−1) m
κ_ε = −(Γ−1) ε ρ² / (n² P)        [ = (ρ²T/n)(∂(μ/T)/∂ε)_n ]
κ_n =  (ρ/(n²P))[ (Γ−1)ε² + P² ]
κ_s = κ_ε + κ_n
ω   ≡ κ_s/κ_ε = m n P /(ε ρ),     α ≡ (∂P/∂ε)_n / c_s² = (Γ−1)/c_s²
```

### 1.4 Cross-source notation cross-walk

| Concept | Kovtun 1907.08191 | BDN 2009.11388 | PMP 2209.09265 / C-code | Shum 2509.15303 | This synthesis |
|---|---|---|---|---|---|
| E-correction | `E−ε` (ε1 Ṫ/T+ε2 θ) | `A` | `A` / `compute_A` | `A` | **`A`** |
| bulk-like | `P−p` (π1,π2) | `Π` | `Π` | `Π` | **`Π`** |
| heat flux | `Q^a` (θ1,θ2) | `Q^a` | `Q^a`/`Qx` | `Q^μ` | **`Q^a`** |
| shear | `−η σ` | `−2η σ` | `−2η σ`/`m2sxx` | `−2η σ` | **`−2η σ^{ab}`** |
| relax. (energy) | `ε1`→ — | `τ_ε` | `τ_ε`=(3/4ε)χ₀e^{3/4} | `τ_eps` | **`τ_ε`** |
| relax. (bulk) | `π1` | `τ_P` | `τ_P=τ_ε/3` (conf.) | `τ_p` | **`τ_P`** |
| relax. (heat) | `θ`=θ1=θ2 | `τ_Q` | `τ_Q`=(3/4ε)λ₀e^{3/4} | `τ_Q` | **`τ_Q`** |

Frame freedom: Kovtun's frame-invariants `f_i = π_i − (∂p/∂ε)ε_i − (∂p/∂n)ν_i`, `ℓ_i = γ_i − (n/ρ)θ_i`, and `η`
are the only genuine physics; everything else is a frame choice. The C reference and PMP fix the
conformal **luminal frame** `(χ₀,λ₀)=(25/4, 25/7)·η₀` (max characteristic speed = c). The NS papers
fix frames via dimensionless hats (§3).

---

## 2. Primitive-recovery algorithm (the STEP-0 trunk module)

The dissipative pieces contain **time derivatives of primitives** (through `u^c∇_c`). So there is no
algebraic 2D root-find for `(ε,u)` per cell. Instead **evolve the primitives by integrating their time
derivatives** (method of lines on primitives), reconstructing `(ε̇, u̇)` each substage by a
**gradient-frozen shifted linear solve** with a perfect-fluid fallback. This is exactly the structure of
`solver.c:compute_xiD / compute_uxD`; PMP 2209.09265 §IV and 2201.12317 prove the solve is *linear*
(hence always invertible for physical states) because `T^{ab}` is linear in the time-derivative primitives.

### 2.1 Variables

- Conserved (evolved by flux divergence): `q = (T^{tt}, T^{ti}, N^t = γ n)` (densitized by `√γ` in curved space).
- Primitives level 0 (trivially evolved): `p0 = (ξ = ln ε, u^i, [n])`. Evolve `ln ε` (not `ε`) for positivity/stability (PMP, C-code).
- Primitives level 1 (reconstructed): `p1 = (ξ̇, u̇^i, [ṅ])` = the time derivatives needed by the constitutive relations.
- Frozen spatial gradients: `∂_x ξ`, `∂_x u^i` (and curved-space spatial covariant pieces) computed once per substage and held fixed during the inversion.

### 2.2 Per-substage algorithm (gradient-frozen shifted Newton + PF fallback)

```
INPUT  : p0 = (ξ, u^i) at cell, conserved q = (T00, T0i), frozen spatial grads (ξ', u^i')
PARAMS : frame coeffs as functions of (ε,n) ; tolerance Δ_visc (TOL)
OUTPUT : p1 = (ξ̇, u̇^i)

1. Build the perfect-fluid (no-time-derivative) stress with the FROZEN spatial gradients:
   T00_PF = T^{tt}|_{ξ̇=u̇=0},   T0i_PF = T^{ti}|_{ξ̇=u̇=0}.
   (Conformal closed forms: solver.c TttPF/TtxPF; general EOS: assemble from §1.2 with ε̇,u̇ → 0.)

2. SHIFT (deficit) = perfect-fluid-vs-actual mismatch, normalized by leading viscosity scale:
   d_tt = (T00_PF − T00)/η₀ ,   d_ti = (T0i_PF − T0i)/η₀ .
   (η₀ = the overall viscosity scale; if η₀ == 0 set deficits 0 to avoid NaN.)

3. FALLBACK test: if |η₀ · d_tt| < Δ_visc AND |η₀ · d_ti| < Δ_visc:
       → viscous correction below numerical-viscosity floor; set d=0 and do the
         PERFECT-FLUID recovery for p0 directly (do NOT use trivial evolution), VISC=0.
         General EOS PF recovery = small 1D/2D Newton root-find (see §2.3).
         Conformal closed form: ε = −T00 + √(4 T00² − 3 T01²),
                                u^x = 3 T01 / √((3T00+ε)² − (3T01)²).
   else VISC=1, continue with the linear BDNK solve.

4. LINEAR BDNK solve. T^{ab} is linear in p1:  q = q0(p0) + η₀ [ J(p0)·p1 + b(p0, frozen grads) ].
   The 2×2 (or 3×3 with n) Jacobian J = ∂(T00,T0i)/∂(ξ̇,u̇^i) is assembled analytically.
   Solve  p1 = J⁻¹ · [ (1/η₀)(q − q0) − b ]  =  (shifted form)  J⁻¹ · d  +  p1_PF,
   where p1_PF = pure spatial-advection part (the C-code's first terms in xiD/uxD that
   survive at zero dissipation). Explicit conformal closed form: solver.c compute_xiD/uxD
   with common denominator DEN = e^{0.75}(9 ch l + 12 ch(l−1)U² + 4(ch(l−3)−l)U⁴),
   ch=χ₀/η₀, l=λ₀/η₀, U=u^x. For general EOS, J⁻¹ is a 2×2 explicit inverse (det = DEN/η₀).

5. RETURN p1 = (ξ̇, u̇^i). Trivially integrate p0: ξ += dt·ξ̇, u^i += dt·u̇^i (MoL substage).
```

**"Shifted Newton" clarification.** For the *conformal* EOS the inversion is fully closed-form (no
iteration). For a *general/tabulated NS EOS* the Jacobian `J` is built numerically (or from the symbolic
EOS), and the "Newton" is a **single Newton step on a linear system** (one solve, no iteration loop)
because `T^{ab}` is exactly linear in `p1`. The only place a genuine iterative Newton appears is the PF
**fallback** root-find in step 3 for a general EOS (`ε = −T00 + √(…)` has no closed form off-conformal).
Recommended fallback solver: 1D Newton on pressure/`W` à la Kastaun (cf. Chabanov-Rezzolla 2311.13027
`f(z)=z−r/h'(z)`), bracketed and clamped.

### 2.3 General-EOS PF fallback (1D Newton)

Given `(T00, T0i)` and EOS `P(ε)` (barotropic) or `P(ε,n)`:
- form `S² = γ^{ij} T0i T0j`, guess `W`, define `f(W) = …` (standard GRHD cons2prim);
- 1D Newton with analytic `dP/dε = c_s²`; reltol `1e-12`, ≤ 50 iters;
- on non-convergence or unphysical state → atmosphere reset (Stage 1C/2).

### 2.4 NeutronStarOscillations.jl contrast (linear benchmark path)

For the **linear** Stages 1A/1B there is *no* cons2prim at all: the evolved variables are the
perturbation primitives themselves (`δu, δε, δλ`), so `con2prim` is the identity. The Julia package
assembles the linear operator once (`Jacobian.compute_jacobian!`), `factorize()` LU once, `ldiv!` each
step; `v3 = ∂_t δλ` is reconstructed algebraically (`compute_v3!`). STEP-0's `bdnk_recovery` module is
therefore exercised *only* by the nonlinear Stages 1C/2 — but it is built and unit-tested at STEP 0
against the conformal closed form so the nonlinear stages inherit a verified trunk.

---

## 3. Causality / stability inequalities (STEP-0 monitor flag)

Monitor pointwise (per cell, every substage); raise the **CAUSALITY FLAG** if any fails.

### 3.1 Full general-frame BDN conditions (2009.11388 Theorem I; PMP eq 20)

Assumption (A1): `ρ, τ_ε, τ_Q, τ_P > 0` and `η, ζ, σ ≥ 0`. Define

```
A = ρ τ_ε τ_Q
B = −τ_ε ( ρ c_s² τ_Q + ζ + 4η/3 + σ κ_s ) − ρ τ_P τ_Q
C = τ_P ( ρ c_s² τ_Q + σ κ_s ) − β_ε ( ζ + 4η/3 )
D = ρ c_s² (τ_ε + τ_Q) + ζ + 4η/3 + σ κ_ε
E = σ ( (∂P/∂ε)_n κ_s − c_s² κ_ε )
```

**Causality** (necessary & sufficient):
```
(a)  ρ τ_Q > η
(b)  B² ≥ 4 A C ≥ 0
(d)  2A > −B ≥ 0
(e)  A + σ κ_s τ_P > −B − C        (equivalently the (e) form in BDN)
```
**Linear stability** (strict, with `η>0`): the five `stabA1..stabE` inequalities of BDN Sec VII (D,E,B,C
above). Use the **reduced hat** form (frame-constraints notebook In[18]) for cheap runtime checks once a
frame is fixed.

### 3.2 Ideal-gas reduced bounds (PMP eq 44, 71; frame notebook headline)

With the frame ansatz `τ_ε=τ_Q=L V̂ τ̂`, `τ_P=2(Γ−1)L V̂`, `η=ρc_s²L η̂`, `ζ=ρc_s²L ζ̂`,
`σ=(V̂Lρc_s²/(−κ_ε))σ̂`, `V̂=(4η/3+ζ)/(ρc_s²L)`:

```
STABILITY (headline):   σ̂ ≤ 1/3
CAUSALITY (headline):   τ̂ ≥ [ (Γ−1)(2−c_s²) + c_s² ] / (1 − c_s²)      (use > to strictly exclude c₊=1)
2nd-law param box:      0 < c_s² < 1,  α ≥ 1,  0 < ω < 3−2√2 ≈ 0.2
```
As `c_s² → 1` the ansatz forces `τ̂ → ∞` (stiff-EOS caveat: a different frame may be needed near
`c_s²≈1`).

### 3.3 Shum-et-al spherically-symmetric Cowling frame (2509.15303 eq 67-71)

Parametrize by `(s_hat, a_hat, q_hat)` (frame) + `(η̂, ζ̂)` (magnitude), `V̂=(4/3)η̂+ζ̂`:
```
τ_ε = V̂ L,  τ_P = ŝ c_s² L V̂,  τ_Q = â L V̂,  β_ε = c_s² â V̂ L,  η = q̂ L c_s² ρ η̂,  ζ = q̂ L c_s² ρ ζ̂
characteristic speeds:
  c₀  = c_s √( q̂ η̂ / (â V̂) )
  c±  = c_s √{ [ â(1+ŝ) + q̂ ± √( q̂² + â²(4q̂+(ŝ−1)²) + 2 â q̂ (1+ŝ) ) ] / (2â) }
Well-posedness + linear stability:   0 < q̂ < ŝ
Causality:   q̂ < [(1−c_s²)/c_s²] · [(1−ŝ c_s²)/(c_s² + 1/â)]   AND   ŝ < 1/c_s²
```
Production frame used: `(ŝ,â,q̂) = (1, 1, 0.999)`.

### 3.4 NeutronStarOscillations.jl radial check (the runtime `check_causal_frame`)

Radial characteristic speeds `c²± = (Λ1 ± √(Λ1²−Λ0))/Λ2` with
```
Λ0 = 4 L⁴ (3ζ+4η)⁴ (τ_P−1) τ_Q² τ_ε (p+ε)² c_s⁴ / 81
Λ1 = L² (3ζ+4η)² (τ_ε + τ_Q(τ_P+τ_ε)) (p+ε) c_s² / 9
Λ2 = 2 L² (3ζ+4η)² τ_Q τ_ε (p+ε) / 9
```
Verdict: require `τ_P ≥ 1`; if `τ_P==1` then `c_s²(0) ≤ τ_Q τ_ε/(τ_Q+τ_ε+τ_Q τ_ε)`; if `τ_P>1` then
`c_s²(0) ≤ ( τ_ε+τ_Q(τ_P+τ_ε) − √(−4(τ_P−1)τ_Q²τ_ε+(τ_ε+τ_Q(τ_P+τ_ε))²) ) / (2(τ_P−1)τ_Q)`.
`τ_P < 1 ⇒ acausal`. (These are the *dimensionless* `τ` of that package's convention; not identical to
the dimensionful `τ_ε` above — keep conventions separate.)

### 3.5 Caballero-Yunes heat-conduction stability criterion (2506.09149) — EXPLICIT

For **radial** modes with heat conduction `κ` (Eckart / BDNK / MIS all share it), the heat operator
`𝔾` is built from the eigenfunction-weighted integral of the factor `(c_s² − c_n²)`:

```
c_s² ≡ (∂p/∂ε)_s   (adiabatic),   c_n² ≡ (∂p/∂ε)_n   (fixed baryon number)

NECESSARY heat-conduction radial-stability criterion:    c_s² − c_n² ≥ 0   (pointwise over the star)
```
i.e. adiabatic sound speed ≥ fixed-`n` sound speed. If `c_n² > c_s²` anywhere, heat conduction can drive
a radial instability. The precise statement is `⟨φ|𝔾|φ⟩` non-destabilizing via
`Y[φ] = −(κT ρ/n) ∂_r[ (ρ/(nT))(c_s²−c_n²) e^{−Φ} ∂_r φ / r² ]`; `c_s² ≥ c_n²` locally guarantees it.
Bulk+shear operator `𝔽` is positive-definite ⇒ **unconditional** radial stability to `ζ, η`.

**STEP-0 causality flag (one line):** raise flag unless (3.1 a,b,d,e) hold pointwise AND the active-frame
reduced bound of (3.2)/(3.3) holds AND (for any heat-conducting run) `c_s² − c_n² ≥ 0`.

---

## 4. Benchmark table (per-stage acceptance gates)

See the structured `benchmark_table`. Headline numbers:

- **STEP0 EOS:** polytrope `p=κε^γ`, `γ=2`, `κ=100`, exact round-trip & `c_s²=γκε^{γ-1}`; ideal-gas
  `c_s²=ΓP/ρ`, `(∂P/∂ε)_n=Γ−1`. Tol exact (machine).
- **STEP0 causality:** for `(η=ζ=1e-2, τ_ε=15, τ_P=1.5, τ_Q=20, L=1)` reproduce `Λ0,Λ1,Λ2` and
  `c_s²_max2` at center; for conformal `(χ₀,λ₀)=(25/4,25/7)η₀` confirm `c±=1` (luminal).
- **STEP0 recovery:** conformal closed-form `ε=−T00+√(4T00²−3T01²)` matched by the linear BDNK solve in
  the inviscid limit to `1e-10`.
- **S1A:** TOV `M,R` for `n=1,κ=100,ε_c=5.5e15 g/cm³`; ideal radial `ω²` first 5 modes; heat criterion
  toy-EOS instability when `c_n²>c_s²`.
- **S1B:** EOS1 `κ=100 km², n=1 → M=1.27 M☉, R=8.86 km`; `ℓ=2` w-mode `(10.4884 kHz, 29.5870 µs)` at
  `η_c=3e29`, drifting to `(10.0898 kHz, 30.8857 µs)` at `η_c=1e31` (frame A1, Table 2); η-mode
  `Im ω → 0` as `η_c→0`; avoided crossing onset `η_c ≳ 1e30`.
- **S1C:** TOV `ρ0,c=0.00128 M☉⁻², M_T=1.4 M☉`, EOS eq.53 `Γ=2,κ=100`; QNM `F=2.69 kHz`, `H1≈4.55 kHz`,
  `H2≈6.36 kHz`; decay rates `1/τ ∈ {0.00157,0.00150,0.00215,0.00182} M☉⁻¹`; stable window
  `τ_ε=(4/3)η̂+ζ̂ ≲ 0.1`.
- **S2:** PMP Bjorken ODE 4th-order conv `Q_N→16` (`τ=1→20, n0=0.1, Γ=4/3`); steady shock
  `e_L=1,v_L=0.8 → e_R=4.4074, v_R=0.41667, η₀=0.2`; IS-contrast migration (2311.13027): no BH collapse,
  shock stalls ~60 km at high ζ.

---

## 5. Project build DAG

See structured `dag_nodes`. Trunk = STEP 0 (`step0.*`). Linear benchmarks Stage 1A/1B branch off the EOS+TOV
sub-trunk but do NOT need `bdnk_recovery` (linear, identity con2prim). Nonlinear Stage 1C is the first real
consumer of `bdnk_recovery` + HRSC + IMEX. Stage 2 (fork A: dynamical GR + collapse) and Stage 3 (fork B:
3+1D Cowling non-radial) both fork from a validated 1C core. Stage 4 = production. **Blocking next nodes:**
`step0.eos`, `step0.con2prim_ideal`, `step0.causality` (all immediately buildable in parallel).

---

## 6. Human decision points

See structured `human_decision_points`:
- `{{EOS_TABLE}}` — default the Shum analytic `p(ε)=[1+2εκ−√(1+4εκ)]/(2κ)` (Γ=2,κ=100) for the 1C
  reproduction trunk; tabulated/SLy deferred to Stage 2/4.
- `{{TRANSPORT}}` — parametrized hats `(η̂,ζ̂,τ̂,σ̂,ŝ,â,q̂)` for benchmarks; physical `η_c [g cm⁻¹ s⁻¹]`
  only at Stage 1B (to match Bussières η_c) and Stage 4.
- `{{FRAME_SET}}` — Shum `(ŝ,â,q̂)=(1,1,0.999)` for 1C; conformal luminal `(25/4,25/7)` for STEP-0
  recovery unit test; Bussières A1/A2/B1/B2 for 1B.
- `{{TOL}}` — `Δ_visc` (PF fallback): default **off** (`TOL<0`, always BDNK solve) for benchmarks; turn on
  + tune down with resolution only in production.
- `{{GR_BACKEND}}` — Cowling (frozen metric) for 1A/1B/1C/3; full dynamical 1+1D GR (harmonic or
  BSSN/Z4c-like) only at Stage 2; recommend reusing an existing NR backend (AthenaK-style) rather than
  hand-rolling.

---

## 7. Risks

- **General-EOS shifted inversion is not in any source in closed form** (only conformal). The symbolic
  2×2 Jacobian for a tabulated NS EOS must be re-derived (Mathematica) or built numerically — highest-risk
  trunk item.
- **Stiff EOS `c_s²→1`** breaks the PMP frame ansatz (`τ̂→∞`); the realistic-EOS Stage-2 may need a
  different frame.
- **Constraint non-propagation** (Fantini-Rubio 2506.06430): the BDNK first-order *reduction* constraints
  grow exponentially in flat-space; relevant if Stage 2 uses an auxiliary-variable recast — no
  constraint-damping scheme exists yet.
- **No shipped numeric benchmark for the Julia linear package** (Keeble citation placeholder); 1A/1B mode
  tables must be cross-checked against Caballero-Yunes (analytic) and Bussières (Table 2) instead.
- **Sign/convention drift** across sources (4π vs 8π, e^ν vs e^{−2Φ}, mostly-plus); pin once in STEP-0 and
  unit-test the TOV exterior `e^ν → 1−2M/R`.

---

## 8. Recommended build order

1. `step0.eos` + `step0.con2prim_ideal` + `step0.causality` (parallel, all blocking).
2. `step0.bdnk_recovery` (conformal closed-form first; then general-EOS linear-solve scaffold).
3. `s1a.tov_background` (shared by 1A/1B/1C).
4. `s1a.radial_eig` + `s1a.heat_criterion` (linear, fast, validates EOS+TOV+causality plumbing).
5. `s1b.axial_wave_eqs` + `s1b.qnm_freqdomain` (reuse TOV; Leaver exterior + interior shooting).
6. `s1c.hrsc_core` → `s1c.imex` → `s1c.qnm_extract` (first nonlinear consumer of `bdnk_recovery`).
7. Fork: `s2.*` (dynamical GR + collapse) and `s3.*` (3+1D Cowling non-radial) in parallel.
8. `s4.production`.
