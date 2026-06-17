# arXiv:2009.11388 — First-Order General-Relativistic Viscous Fluid Dynamics (BDNK)

**Authors:** Fabio S. Bemfica, Marcelo M. Disconzi, Jorge Noronha
**Ref:** Phys. Rev. X 12, 021044 (2022). Source TeX: `Paper_update_final.tex` (arXiv e-print).
**Role:** Foundational BDNK paper at finite baryon density, in the Eckart-like (particle) frame `J^mu = n u^mu`. Definitive causality + strong-hyperbolicity + stability theorems.

## Conventions (l.300)
- Metric signature **mostly plus** `(-+++)`. Greek 0..3, Latin 1..3.
- Natural units `c = hbar = k_B = 1`. `nabla_mu` = spacetime covariant derivative.
- `u_mu u^mu = -1`, `Delta_{mu nu} = g_{mu nu} + u_mu u_nu` (projector orthogonal to u).
- TSS projector: `Delta^{mu nu}_{ab} = 1/2 (Delta^mu_a Delta^nu_b + Delta^mu_b Delta^nu_a - (2/3) Delta^{mu nu} Delta_{ab})`.
- Shear: `sigma^{mu nu} = Delta^{mu nu a b} nabla_a u_b`.
- Vorticity: `omega_{mu nu} = 1/2 (Delta_mu^l nabla_l u_nu - Delta_nu^l nabla_l u_mu)`.
- `rho = epsilon + P` (enthalpy density). First law: `epsilon + P = T s + mu n`.

## EOM / governing equations
`nabla_mu J^mu = 0`, `nabla_mu T^{mu nu} = 0`, coupled to Einstein `R_{mu nu} - R/2 g_{mu nu} = 8 pi G T_{mu nu}` (eqs 1-2). EOS `P = P(epsilon, n)`.

## General Kovtun-frame decomposition (eqs 5-6, l.974-982)
```
J^mu = N u^mu + J^mu_perp
T^{mu nu} = E u^mu u^nu + P_cal Delta^{mu nu} + u^mu Q^nu + u^nu Q^mu + T^{mu nu}_TSS
```
with `N = -u_mu J^mu`, `E = u_mu u_nu T^{mu nu}`, `P_cal = Delta_{mu nu}T^{mu nu}/3`,
`J^nu_perp = Delta^nu_mu J^mu`, `Q^nu = -u_mu T^{m l} Delta_l^nu`, `T^{mu nu}_TSS = Delta^{mu nu}_{ab} T^{ab}`.

### Most general first-order constitutive relations (eq generaldefKovtun, l.1006-1013)
```
E    = epsilon + e1 (u^a nabla_a T)/T + e2 nabla_a u^a + e3 u^a nabla_a(mu/T)
P_cal= P       + pi1(u^a nabla_a T)/T + pi2 nabla_a u^a + pi3 u^a nabla_a(mu/T)
N    = n       + nu1(u^a nabla_a T)/T + nu2 nabla_a u^a + nu3 u^a nabla_a(mu/T)
Q^mu = th1 (Delta^{mu n} nabla_n T)/T + th2 u^a nabla_a u^mu + th3 Delta^{mu n} nabla_n(mu/T)
J^mu = g1  (Delta^{mu n} nabla_n T)/T + g2  u^a nabla_a u^mu + g3  Delta^{mu n} nabla_n(mu/T)
T^{mu nu}_TSS = -2 eta sigma^{mu nu}
```
Thermodynamic consistency forces `g1 = g2`, `th1 = th2`. 14 scalars + eta; 8 removable by first-order field redefinitions, leaving eta, zeta, sigma + 3 relaxation times.

## THE WORKING FRAME OF THIS PAPER (eq finaltheory, l.1049-1056)
Eckart-type particle frame `J^mu = n u^mu` (so `g_i = nu_i = 0`); uses (epsilon, n, u^mu) as variables.
```
J^mu       = n u^mu
T^{mu nu}  = (epsilon + A) u^mu u^nu + (P + Pi) Delta^{mu nu} - 2 eta sigma^{mu nu}
             + u^mu Q^nu + u^nu Q^mu
A   = tau_eps [ u^l nabla_l epsilon + (epsilon+P) nabla_l u^l ]                         (eq 10c)
Pi  = -zeta nabla_l u^l + tau_P [ u^l nabla_l epsilon + (epsilon+P) nabla_l u^l ]
Q^nu= tau_Q (epsilon+P) u^l nabla_l u^nu + beta_eps Delta^{nu l} nabla_l epsilon
             + beta_n Delta^{nu l} nabla_l n
```
with (eq definebetas, l.1062-1063):
```
beta_eps = tau_Q (dP/d epsilon)_n   + sigma T (eps+P)/n (d(mu/T)/d epsilon)_n
beta_n   = tau_Q (dP/d n)_epsilon   + sigma T (eps+P)/n (d(mu/T)/d n)_epsilon
```
- `tau_eps, tau_P, tau_Q` = relaxation-time-dimension coefficients (the BDNK frame coefficients).
- `eta` shear, `zeta` bulk viscosity, `sigma` charge/baryon conductivity.
- On shell: `A = O(d^2)`, `Pi = -zeta nabla u + O(d^2)`, `Q_nu = sigma T (eps+P)/n Delta nabla(mu/T) + O(d^2)`. Reduces to Eckart at O(d).
- Eckart heat flux relation (l.1037-1039): on shell `(eps+P) u^l nabla_l u^mu + Delta^{mu l}nabla_l P = O(d^2)`, and `dP/(eps+P) = dT/T + nT/(eps+P) d(mu/T)`.

### Explicit second-order EOM (eq EOMcurrent, l.1075-1079)
```
u^l nabla_l n + n nabla_l u^l = 0
u^l nabla_l eps + (eps+P) nabla_l u^l = - u^l nabla_l A - (A+Pi) nabla_l u^l
   - nabla_mu Q^mu - Q^mu u^l nabla_l u_mu + 2 eta sigma_{mu nu} sigma^{mu nu}
(eps+P) u^n nabla_n u^b + Delta^{b l} nabla_l P = -(A+Pi) u^n nabla_n u^b - Delta^{b l} nabla_l Pi
   + Delta^b_l nabla_mu(2 eta sigma^{m l}) - u^l nabla_l Q^b - (4/3) nabla_l u^l Q^b
   - Q_mu sigma^{m b} - Q_mu omega^{m b}
```

## Speed of sound & kappa (l.1153-1160)
```
c_s^2 = (dP/d eps)_sbar = (dP/d eps)_n + (n/rho)(dP/d n)_eps
kappa_s = (rho^2 T/n)[d(mu/T)/d eps]_sbar
        = (rho^2 T/n)[d(mu/T)/d eps]_n + T rho [d(mu/T)/d n]_eps
```
(sbar = equilibrium entropy per particle.)

## THEOREM I — Causality (necessary & sufficient), l.1164-1188
Assumption (A1): `rho, tau_eps, tau_Q, tau_P > 0` and `eta, zeta, sigma >= 0`.
Then causality holds iff (eq C_conditions, l.1179-1184):
```
(a) rho tau_Q > eta
(b) [ tau_eps (rho c_s^2 tau_Q + zeta + 4eta/3 + sigma kappa_s) + rho tau_P tau_Q ]^2
       >= 4 rho tau_eps tau_Q [ tau_P (rho c_s^2 tau_Q + sigma kappa_s) - beta_eps (zeta + 4eta/3) ] >= 0
(d) 2 rho tau_eps tau_Q > tau_eps (rho c_s^2 tau_Q + zeta + 4eta/3 + sigma kappa_s) + rho tau_P tau_Q >= 0
(e) rho tau_eps tau_Q + sigma kappa_s tau_P
       > tau_eps (rho c_s^2 tau_Q + zeta + 4eta/3 + sigma kappa_s)
         + rho tau_P tau_Q (1 - c_s^2) + beta_eps (zeta + 4eta/3)
```
Holds with or without dynamical metric. (Labels a,b,d,e are the paper's; "c" was commented out.)

The characteristic tensor (principal part, l.1146):
`C^{mu a b}_nu = (tau_P rho - zeta - eta/3) Delta^{mu(a} delta^{b)}_nu + (rho tau_Q u^a u^b - eta Delta^{ab}) delta^mu_nu`.

## Strong hyperbolicity / Local well-posedness (Sec V, Thm II, Prop I)
Under the same conditions the full nonlinear GR system is strongly hyperbolic; the first-order reduced system is diagonalizable; LWP in Sobolev `H^{N}`, `N >= 5`. Key tensor (l.1272):
`Pi^{m l a}_nu = -eta(Delta^{m l} delta^a_nu + Delta^{a l} delta^m_nu) + (rho tau_P - zeta + 2eta/3) Delta^{m a} delta^l_nu`.

## Stability theorem (Sec VII, l.1586-1620)
Define
```
D = rho c_s^2 (tau_eps + tau_Q) + zeta + 4eta/3 + sigma kappa_eps
E = sigma [ p'_eps kappa_s - c_s^2 kappa_eps ]
  = sigma T rho [ (dP/d eps)_n (d(mu/T)/d n)_eps - (dP/d n)_eps (d(mu/T)/d eps)_n ]   (>= 0)
kappa_s = kappa_eps + kappa_n,  kappa_eps=(T rho^2/n)[d(mu/T)/d eps]_n,
kappa_n=(T rho)[d(mu/T)/d n]_eps,  p'_eps=(dP/d eps)_n
B = -tau_eps (rho c_s^2 tau_Q + zeta + 4eta/3 + sigma kappa_s) - rho tau_P tau_Q   (|B| = -B > 0)
C = tau_P (rho c_s^2 tau_Q + sigma kappa_s) - beta_eps (zeta + 4eta/3)
```
**Statement:** The system is linearly stable (Cowling approx, flat `eta_{mu nu}`, `delta g = 0`) if it is causal within the **strict** form of C_conditions, with additionally `eta > 0`, and (eq stability, l.1603-1609):
```
(a) (tau_eps+tau_Q)|B| >= tau_eps tau_Q D >= rho c_s^2 tau_eps tau_Q (tau_eps+tau_Q)
(b) (tau_eps+tau_Q)|B| D + rho tau_eps tau_Q (tau_eps+tau_Q) E > tau_eps tau_Q D^2 + rho(tau_eps+tau_Q)^2 C
(c) c_s^2 D - E >= rho c_s^4 (tau_eps+tau_Q)
(d) (tau_eps+tau_Q)[ |B|(c_s^2 D - 2E) + 2 c_s^2 rho tau_eps tau_Q E + C D ]
       > 2 c_s^2 rho (tau_eps+tau_Q)^2 C + tau_eps tau_Q D (c_s^2 D - E)
(e) |B| D [ C(tau_eps+tau_Q) + E tau_eps tau_Q ] + 2 rho tau_eps tau_Q (tau_eps+tau_Q) C E
       > rho C^2 (tau_eps+tau_Q)^2 + tau_eps tau_Q (C D^2 + rho tau_eps tau_Q E^2) + B^2 E (tau_eps+tau_Q)
```
Derived from Routh-Hurwitz on the sound-channel polynomial (l.1656). New general theorem (Sec VI): strong hyperbolicity + causality + stability-in-LRF => stability in any boosted frame.

## Entropy production (l.1100-1102)
`nabla_mu S^mu = 2 eta sigma^2/T + zeta (nabla_mu u^mu)^2 / T + sigma T [Delta nabla(mu/T)]^2 + O(d^3) >= 0` when `eta, zeta, sigma >= 0`.

## Open questions / build notes
- Paper does NOT prescribe tau_eps, tau_P, tau_Q numerically; they are EOS-dependent free coefficients to be chosen within the inequality regions. No benchmark numbers in this paper.
- `kappa_s` vs `kappa_eps` distinction matters: causality uses kappa_s, stability D uses kappa_eps.
- For zero-chemical-potential reduction (sigma=0, n decoupled) the theory collapses to Bemfica-Disconzi-Noronha 2019 / Kovtun 2019 single-component BDNK.
