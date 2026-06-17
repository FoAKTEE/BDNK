# Fantini & Rubio (2025) — Constraint evolution in first-order viscous relativistic fluids

- arXiv: **2506.06430**; PRD **112**, 063038 (2025), published 22 Sep 2025.
- Authors: Delfina Fantini (UNC, Córdoba), Marcelo E. Rubio (GSSI).
- Setting: **conformal** BDNK viscous fluid, **flat spacetime**, plane (1+1) symmetry. NOT coupled to dynamical GR in the numerics (Minkowski).

## One-line takeaway (IMPORTANT for our project)
This paper is a **negative result / cautionary note**: it shows that the differential constraints arising from the BDNK *first-order reduction* (recasting the 2nd-derivative BDNK system into a 1st-order system via auxiliary variables) do **NOT** propagate correctly — constraint norms grow **exponentially** even though the fluid fields evolve smoothly. The authors propose **NO constraint-damping scheme**; they pose finding a constraint-preserving reduction as an open problem. So for Stage-2 (constraint damping) this source primarily warns against a naive auxiliary-variable reduction and motivates building damping, but gives no ready damping prescription.

## (1) BDNK conformal stress tensor & constitutive relations
T^{ab} = (ε + A)(u^a u^b + (1/3)Π^{ab}) + 2 u^{(a} Q^{b)} − η σ^{ab}   (Eq. 1)

- A = 3χ( θ^{-1} u^c ∇_c θ + (1/3) ∇_c u^c )   (Eq. 2)
- Q_a = λ( θ^{-1} Π_a^c ∇_c θ + u^c ∇_c u_a )   (Eq. 3)
- σ^{ab} = Π_a^c∇_c u^b + Π_b^c∇_c u^a − (2/3)Π^{ab}∇_c u^c   (Eq. 4)
- Π^{ab} = g^{ab} + u^a u^b (projector), X ≡ log θ.

Transport coefficients tied to shear viscosity η(θ) > 0 (analytic):
χ = a₁ η,  λ = a₂ η.

## (2) Equation of state
Pure radiation / conformal: ε = ε₀ θ⁴,  ε₀ > 0 const, θ > 0 temperature.   (Eq. 5)

## (3) Causality / stability inequalities on coefficients
a₁ > 4,  a₂ ≥ 3 a₁ / (a₁ − 1).
These give linear stability + causality (BDNK conditions) but — the central point — do NOT guarantee constraint propagation.

## (4) First-order reduction: auxiliary variables & evolution system
Auxiliary (gradient) variables (Eq. 7):
S_a^b = Π_a^c ∇_c u^b,   S^a = u^c ∇_c u^a.

Evolution variables: {A, Q^a, S^a, S_a^b, θ, u^a}.

Evolution system (Eqs. 8–13):
- u^c∇_c A + ∇_c Q^c + r₁ = 0                                   (Eq. 8)
- Π^{ac}∇_c A + 3 u^c∇_c Q^a + B_a^{cde}∇_e S^c_d + r₂ = 0        (Eq. 9)
- −(1/χ)Π^{ac}∇_c A + (3/λ)u^c∇_c Q^a − 3u^c∇_c S^a + Π^{ac}∇_c S_d^d + r₃ = 0  (Eq. 10)
- u^c∇_c S_a^b − Π_a^d∇_d S^b + r₄ = 0                           (Eq. 11)
- (1/θ)u^c∇_c θ + (1/3)∇_c u^c + r₅ = 0                          (Eq. 12)
- (1/θ)Π^{ac}∇_c θ + u^c∇_c u^a + r₆ = 0                         (Eq. 13)

Source terms (algebraic, Appendix A, Eqs. 41–46), partial:
r₁ = (4/3)A S + Q^a s_a − (1/2)η σ² + (4ε₀/3χ) A θ⁴
r₅ = −A/(3χ),   r₆^a = −Q^a/λ.

## (5) Constraints
Algebraic (from auxiliary defs): S_a^b = Π_a^c∇_c u^b,  S^a = u^c∇_c u^a.

Differential constraints (must hold if reduction = original BDNK):
- C₁^a ≡ Q^a − λ Π^{ab}∇_b X = 0          (Eq. 19; X = log θ)
- C₂^a ≡ (4ε/3λ) Q^a + ((λ+χ)/3χ) Π^{ab}∇_b A   (from Eq. 15)

## (6) Constraint propagation analysis (CORE RESULT)
For uniform-velocity configs (∇_a u^b = 0):
- C₁^a propagates: Ċ₁^a = 0 along u^a.
- C₂^a FAILS. Conservation of C₂ would demand (constraint-conserving eqn):
  3χ Ẍ + 4ε Ẋ + λ ε^{-1} Π^{ab}∇_a ε ∇_b X − (3λ/4ε)(λ+χ) ΔẊ = 0
  but the actual evolution gives:
  3χ Ẍ + 4ε X + λ ΔX = 0.
  These are incompatible ⇒ "constraint propagation fails."

There is NO explicit closed homogeneous strongly-hyperbolic constraint-subsystem derived (despite the abstract framing); the analysis is by direct substitution and numerics. BDNK's strong-hyperbolicity (diagonalizable principal part, real eigenvalues) is questioned because the constraints are not conserved.

## (7) Plane-symmetric flat-spacetime equations actually evolved
Quasilinear 1+1 form:  ∂_t Φ + M(Φ) ∂_x Φ + Ψ(Φ) = 0   (Eq. 24)
Φ = { A, Q^1, S^1, S^{11}, θ, u^1 }.
- (u^1, θ) decouple into a 2×2 subsystem (Eq. 32), principal matrix has γ/(3+2(u^1)^2) terms (γ = Lorentz factor √(1+(u^1)²)).
- Remaining 4 eqs (Eqs. 8–11 projected) → 4×4 linear solve for time derivatives; explicit in Appendix B.

Constraint relations in planar symmetry (Eqs. 35–36), used to define monitored constraints:
S^1 = (u^1)² ∂_x u^1 / [(3+2(u^1)²)√(1+(u^1)²)]
      − 3 u^1 √(1+(u^1)²) ∂_x θ / [θ(3+2(u^1)²)]
      − u^1(−3χ Q^1 + λ A u^1)√(1+(u^1)²) / [λχ(3+2(u^1)²)]
S^{11} = 3(1+(u^1)²) ∂_x u^1 / (3+2(u^1)²)
      − 3 u^1(1+(u^1)²) ∂_x θ / [θ(3+2(u^1)²)]
      + u^1(3χ Q^1 + λ A u^1)(1+(u^1)²) / [λχ(3+2(u^1)²)]

Monitored constraint quantities (Eqs. 38–39):
C₁ ≡ S^1 − S^1(Φ,∂_xΦ),   C₂ ≡ S^{11} − S^{11}(Φ,∂_xΦ).

## (8) Numerical method
- Spatial: centered 2nd-order finite differences with summation-by-parts (SBP).
- Dissipation: 6th-order Kreiss–Oliger, σ_diss = 0.25 (Eq. 37).
- Time: explicit 4th-order Runge–Kutta (RK4).
- Domain: x ∈ [−80, 80], N = 10001 points, Δx = 0.016, Δt = 0.004 (CFL ≈ 0.25).
- BCs: periodic (domain large enough to avoid boundary reflections during run).
- Final time t_f = 80.
- Convergence factor Q(t) ≈ 4; observed orders p≈2 (space), p≈4 (time).

## (9) Numerical benchmark values
Gaussian initial data (Table 1), σ = 5, x_c = 0, form f_0 + f_1·exp(−(x−x_c)²/σ²):
- v (=u^1 boost / velocity): f_0 = 0.3, f_1 = 0.1
- θ: f_0 = 300, f_1 = 50
- A: f_0 = 0, f_1 = 0
- Q^1: f_0 = 0, f_1 = 0.5

Key result (Fig. 2): L² and L∞ norms of C₁, C₂ grow exponentially from early times; C₂ grows faster, ~4 orders of magnitude, while fields stay smooth.

## Relevance to our BDNK Julia build (Stage 2 = constraint damping)
- WARNING source: a naive auxiliary-variable first-order reduction of BDNK introduces differential constraints that do NOT propagate and blow up exponentially.
- Implication: if we use such a reduction we MUST add constraint damping (Z4c/GBSSN-style) or use a reduction that closes without spurious constraints. The paper gives none — open problem.
- Units: conformal, dimensionless θ~300 scale; flat spacetime, geometrized; no NS-specific scales here.

## Open questions
- No explicit constraint-subsystem matrix/eigenvalues given — abstract framing of "strongly hyperbolic constraint system" not matched by an explicit damping-ready form in the extracted text; would need to read full PDF Sec. III/IV to confirm whether a homogeneous constraint evolution system is written.
- Full r₂, r₃, r₄ and B_a^{cde} tensor not extracted (Appendix A/B); needed if we replicate their exact reduction.
- Whether a₁, a₂ choices used in the runs are stated (e.g., a₁, a₂ numeric values) not captured.
