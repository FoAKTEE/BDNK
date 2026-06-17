# Keeble & Pretorius 2025 — First-Order Viscous Relativistic Hydrodynamics on the Two-Sphere

- arXiv: **2508.20998** (v2). Published PRD **112, 124034** (2025), APS DOI 10.1103/d4wd-zj7w.
- Authors: Lennox S. Keeble, Frans Pretorius (Princeton).
- One line: BDNK first-order viscous hydro for a **4D conformal fluid in Minkowski**, with the fluid dynamically constrained to the surface of a geometric 2-sphere, solved by 4th-order finite difference in cubed-sphere coordinates.
- Relevance to our project: **Stage 3 angular-derivative gradient terms.** The 2-sphere restriction gives explicit angular (theta, phi) covariant-derivative structure for the BDNK gradient terms (acceleration, expansion, shear, vorticity) and a concrete cubed-sphere FD implementation to mirror.

## 1. BDNK stress-energy tensor (general frame), Eq. (6a), (7), (8)

```
T^{munu} = (e + A) u^mu u^nu + (P + Pi) Delta^{munu} + Q^mu u^nu + u^mu Q^nu - 2 eta s^{munu}      (6a)
```
Out-of-equilibrium corrections:
```
A   = tau_e [ u^l ∇_l e + (e+P) ∇_l u^l ]                                   (7a)
Pi  = -zeta ∇_l u^l + tau_P [ u^l ∇_l e + (e+P) ∇_l u^l ]                    (7b)
Q^mu= tau_Q (e+P) u^l ∇_l u^mu + beta_e Delta^{ml} ∇_l e + beta_n Delta^{ml} ∇_l n   (8a)
s^{munu} = (1/2) Delta^{mr} Delta^{ns} ( ∇_r u_s + ∇_s u_r - (2/3) g_{rs} ∇_l u^l )   (8b)
```
- Delta^{munu} = g^{munu} + u^mu u^nu (spatial projector). s^{munu} is the shear (traceless, transverse).
- General-frame transport coefficients: tau_e, tau_P, tau_Q (relaxation-time-like, "frame" coefficients), beta_e, beta_n, plus equilibrium eta (shear), zeta (bulk).

## 2. Conformal equation of state, Eq. (10), (11), (13)

Zero chemical potential conformal fluid:
```
P(e,n) = e/3
tau_P = tau_e / 3
zeta  = 0
beta_e = tau_Q / 3 ,  beta_n = 0                                            (10a-b)
```
Coefficient parametrization (relaxation times in terms of chi, lambda):
```
tau_e = 3 chi / (4 e) ,   tau_Q = 3 lambda / (4 e)                          (11)
```
Conformal dimensional scaling (e ∝ T^4; {eta,chi,lambda} ∝ T^3):
```
eta    = eta0   e^{3/4}
chi    = chi0   e^{3/4}
lambda = lambda0 e^{3/4}                                                    (13)
```

## 3. Causality / stability / well-posedness inequalities, Eq. (12)

```
eta > 0
chi = a1 * eta
lambda >= 3 eta a1 / (a1 - 1)
a1 >= 4
```
**Frame B** choice used in runs (max characteristic speed = speed of light), Eq. (14):
```
(lambda0, chi0) = (25 eta0 / 7 ,  25 eta0 / 4)
```
(corresponds to a1 = 25/4 = 6.25 since chi0/eta0 = 25/4.)

## 4. Geometry & evolved (primitive) variables

Line element (flat 4D, fluid pinned to S^2 of radius R), Eq. (17):
```
ds^2 = -dt^2 + dr^2 + R^2 dtheta^2 + R^2 sin^2(theta) dphi^2
```
Velocity ansatz (no radial flow), Eq. (18):
```
u^mu = [u^t, 0, u^theta, u^phi]^T
```
Normalization, Eq. (19):
```
u^t = sqrt( 1 + R^2 (u^theta)^2 + R^2 sin^2(theta) (u^phi)^2 )
```
**Evolved primitive fields: {e, n, u^theta, u^phi}.** No conservative->primitive (cons2prim) recovery step: the BDNK and perfect-fluid densities are written directly in terms of the primitives (Sec. III.1). u^t is algebraically reconstructed from the normalization each step.

Conservation laws solved, Eq. (1):
```
∇_mu T^{munu} = 0 ,   ∇_mu J^mu = 0
```
Vorticity 2-form, Eq. (24): omega_{munu} = ∇_[mu u_nu]. Scalar vorticity density used for diagnostics, Eq. (26):
```
W^0 = -csc(theta) ∂_phi u^theta + sin(theta) ∂_theta u^phi + 2 cos(theta) u^phi
```

## 5. Numerical method (Sec. III.2, App. D)

- **Cubed-sphere coordinates**: six singularity-free coordinate charts covering S^2 continuously (avoids the polar coordinate singularity of theta,phi).
- **Method of lines.** Spatial: 4th-order-accurate centered FD stencils. Temporal: 4th-order Runge-Kutta (RK4).
- **Kreiss-Oliger dissipation** applied as a filtering operation between time steps; KO strength = 0.5 in all runs.
- Inter-patch boundaries handled by cubic-spline interpolation across patch edges, with Jacobian transforms for the (vector) velocity components.
- CFL: lambda = Delta_t / h (grid spacing h). Runs use lambda in {0.1, 0.2}.
- Convergence: 4th-order convergence demonstrated (Table 2 reports shrinking frequency residuals, e.g. Euler 0.004% -> 0.0005% -> 0.0003% across N=2^3+1, 2^4+1, 2^5+1).

### Run table (Table 1)
| Configuration       | Grid pts | CFL lambda | KO   |
|---------------------|----------|------------|------|
| 1D Gaussian         | 2^14     | 0.1        | 0.5  |
| 2D Perturbations    | 2^12     | 0.2        | 0.5  |
| 2D Gaussian         | 2^14     | 0.2        | 0.5  |
| 2D Kelvin-Helmholtz | 2^16     | 0.1        | 0.5  |

## 6. Initial data

Gaussian pulse, Eq. (41):
```
e(theta,phi) = e0 + A exp(-theta^2 / w^2) ,   u^theta = u^phi = 0
e0 = 0.1, A = 0.4, w = 45 deg
```
Kelvin-Helmholtz, Eq. (42):
```
e      = 1
n      = 1 + (1/2)[ tanh((theta-theta1)/a) - tanh((theta-theta2)/a) ]
u^theta= -A sin(2 phi) [ exp(-(theta-theta1)^2/sigma^2) + exp(-(theta-theta2)^2/sigma^2) ]
u^phi  = (cs/4)[ tanh((theta-theta1)/a) - tanh((theta-theta2)/a) - 1 ]
cs = 1/sqrt(3), a = 0.08, sigma = 0.2, theta1 = 70 deg, theta2 = 110 deg, A in {0, 0.01}
```

## 7. Thermo / viscosity

- KSS bound, Eq. (16): (eta/s)_min = hbar/(4 pi kB) = 1/(4 pi) in natural units.
- s = (e + P)/T, with P = e/3 (Sec. II.2).
- Studied eta/s = (1..20) * (eta/s)_min, i.e. {1,2,3,10,20}/(4 pi). For eta/s >= 10/(4pi) on the Gaussian, steep gradients form and the flow diverges from equilibrium (numerical evidence of singularity formation from smooth data).
- e ∝ T^4; the paper works dimensionlessly (e0/T^4 ≈ 10 quoted dimensionally). No explicit T(e) closed form given in the excerpt.

### Regime-of-validity diagnostics (Sec. II.3)
- Weak-energy-condition probe: u_mu u_nu T^{munu} = e + A.
- Relative first-order correction: |T^{munu}_(1)| / |T^{munu}_(0)| ; values >~ 1 signal viscous terms dominate => out of validity regime.

## Open questions
- Exact closed-form characteristic-speed formula not given in the excerpts (only the frame-B "= c" statement and a1>=4 constraint).
- Explicit T(e) normalization constant and the precise relation eta0 <-> (eta/s) not copied verbatim.
- Detailed cubed-sphere patch overlap/interpolation stencil widths (App. D) not extracted in full.
- Whether n is genuinely advected nontrivially in the conformal (beta_n=0) case, vs only used as a passive tracer in KH.
