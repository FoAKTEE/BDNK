# BDNKStar Reproduction — Validation Ledger

Merged exact-target ledger across 9 source papers, organized by DAG node. Each entry lists the concrete numeric target, the conditions under which it holds, the tolerance, and the source (table / equation / text line). Targets sourced from **tables** (highest confidence, directly tabulated numbers) are flagged ⭐ and listed in the "Reproduce first" set at the bottom.

Legend:
- ⭐ = table-sourced numeric target (highest-confidence; reproduce first)
- ◻ = text-line / caption numeric target (medium confidence)
- ▽ = figure-only target (qualitative or order-of-magnitude; read off plot)
- ⊙ = analytic / symbolic target (exact equation, no number to match)

Units: M_sun (geometric G=c=1) unless stated; cgs g cm^-1 s^-1 for viscosities; kHz / µs / s^-1 for physical-unit conversions.

---

## DAG node: `step0.causality` — BDNK causality/stability theorems

### Paper arXiv:2009.11388 (Bemfica-Disconzi-Noronha, PRX 12, 021044) — `step0.causality`
Purely analytic paper (no tables/figures). Reproduction = symbolic verification of the causality + stability inequalities, plus one numeric non-emptiness unit-test.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 0.1 ⊙ | Causality (C-a) shear sub-luminality | `rho*tau_Q > eta` (rho=eps+P) | Theorem I, (A1): tau_eps,tau_Q,tau_P>0; eta,zeta,sigma>=0; nec.+suff. | strict `>` | Eq. C_condition_a (l.1179) |
| 0.2 ⊙ | Causality (C-b) sound-sector reality/non-neg | `[tau_eps(rho c_s^2 tau_Q+zeta+4eta/3+sigma kappa_s)+rho tau_P tau_Q]^2 >= 4 rho tau_eps tau_Q[tau_P(rho c_s^2 tau_Q+sigma kappa_s)-beta_eps(zeta+4eta/3)] >= 0` | Theorem I under (A1) | first `>=0`, second `>=0` | Eq. C_condition_b (l.1180-1181) |
| 0.3 ⊙ | Causality (C-d) c_+<1 & c_->=0 | `2 rho tau_eps tau_Q > tau_eps(rho c_s^2 tau_Q+zeta+4eta/3+sigma kappa_s)+rho tau_P tau_Q >= 0` | Theorem I under (A1) | first strict `>`, second `>=0` | Eq. C_condition_d (l.1183) |
| 0.4 ⊙ | Causality (C-e) c_+<1 (A+B+C>0) | `rho tau_eps tau_Q + sigma kappa_s tau_P > tau_eps(rho c_s^2 tau_Q+zeta+4eta/3+sigma kappa_s)+rho tau_P tau_Q(1-c_s^2)+beta_eps(zeta+4eta/3)` | Theorem I under (A1) | strict `>` | Eq. C_condition_e (l.1184) |
| 0.5 ⊙ | Char. determinant coeffs + speeds | `A=rho tau_eps tau_Q; B=-tau_eps(rho c_s^2 tau_Q+zeta+4eta/3+sigma kappa_s)-rho tau_P tau_Q; C=tau_P(rho c_s^2 tau_Q+sigma kappa_s)-beta_eps(zeta+4eta/3); c_1=eta/(rho tau_Q); c_pm=(-B±sqrt(B^2-4AC))/(2A)`, mult n_1=3,n_pm=1 | harmonic gauge | `0<=c_a<1` strict upper | Eq. det_a/det_b + ABC (l.1803-1815) |
| 0.6 ⊙ | Assumption (A1) | `tau_eps,tau_Q,tau_P>0` and `eta,zeta,sigma>=0` | base assumption | strict for taus, >=0 visc. | (A1) l.1171 |
| 0.7 ⊙ | Stability (S-a) | `(tau_eps+tau_Q)|B| >= tau_eps tau_Q D >= rho c_s^2 tau_eps tau_Q(tau_eps+tau_Q)` | Cowling, causal-strict + eta>0; D=rho c_s^2(tau_eps+tau_Q)+zeta+4eta/3+sigma kappa_eps | `>=` | Eq. stability_a (l.1603) |
| 0.8 ⊙ | Stability (S-b) | `(tau_eps+tau_Q)|B| D + rho tau_eps tau_Q(tau_eps+tau_Q)E > tau_eps tau_Q D^2 + rho(tau_eps+tau_Q)^2 C` | as (S-a); E=sigma[p'_eps kappa_s-c_s^2 kappa_eps]>=0 | strict `>` | Eq. stability_b (l.1605) |
| 0.9 ⊙ | Stability (S-c) | `c_s^2 D - E >= rho c_s^4(tau_eps+tau_Q)` | as (S-a) | `>=` | Eq. stability_c (l.1606) |
| 0.10 ⊙ | Stability (S-d) | `(tau_eps+tau_Q)[|B|(c_s^2 D-2E)+2 c_s^2 rho tau_eps tau_Q E+C D] > 2 c_s^2 rho(tau_eps+tau_Q)^2 C + tau_eps tau_Q D(c_s^2 D-E)` | as (S-a) | strict `>` | Eq. stability_d (l.1607) |
| 0.11 ⊙ | Stability (S-e) | `|B|D[C(tau_eps+tau_Q)+E tau_eps tau_Q]+2 rho tau_eps tau_Q(tau_eps+tau_Q)C E > rho C^2(tau_eps+tau_Q)^2 + tau_eps tau_Q(C D^2+rho tau_eps tau_Q E^2)+B^2 E(tau_eps+tau_Q)` | as (S-a); RH (iv) | strict `>` | Eq. stability_e (l.1608-1609) |
| 0.12 ◻ | **Non-emptiness numeric unit-test** (THE one hard number) | `tau_Q=tau_eps; tau_P=c_s^2 tau_eps; c_s^2=p'_eps=1/2; rho tau_eps=8(zeta+4eta/3); kappa_eps=kappa_s/2=1/4 (kappa_s=1/2)`; verify ALL (C) + (S) hold strictly for `sigma/(zeta+4eta/3)={0, 1/4, 1}` | P=P(eps) barotrope, Cowling; zeta+4eta/3>0, eta>0 | all conditions strict | Sec. "Fulfilling...", l.1712 |
| 0.13 ⊙ | Hydro dispersion modes | shear `omega=-i k^2 eta/(eps+P)`; sound `omega=±c_s k - i k^2 Gamma_s/2`; heat `omega=-i D k^2`; non-hydro `omega=-i/tau_Q`, `-i/tau_eps` | LRF linearization, low-k | leading low-k | l.1658 |

### Paper arXiv:1907.08191 (Kovtun, JHEP 10 (2019) 034) — `step0.causality`
Purely analytic; general-frame stability. Reproduction = Routh-Hurwitz / polynomial-root verification.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 0.14 ⊙ | Shear-channel stability | `theta > eta > 0` | uncharged, general frame, all k, any v0; theta1=theta2 | exact | Eq. (4.7) eq:shear-stability |
| 0.15 ⊙ | Sound stability cond.1 (RH e>0) | `eps2 + pi1 > gamma_s/v_s^2 + v_s^2 eps1`; gamma_s=(4/3)eta+zeta | uncharged at rest, all k | exact | Eq. (4.16) eq:sound-c1 |
| 0.16 ⊙ | Sound stability cond.2 (nonlinear RH) | `(eb1^2/v_s^2)+v_s^2(eb1-eb2)(eb1+thb)^2(eb1-pb1)+(eb1+thb)[2 eb1^2-eb1(eb2+pb1)+(thb+eb2)(thb+pb1)] > 0`; bars = /gamma_s, eb1=v_s^2 eps1/gamma_s | uncharged at rest, all k | exact | Eq. (4.17) eq:sound-c2 |
| 0.17 ⊙ | Large-k sound speed (causality) | `(eb1 thb/v_s^2)c^4 + [eb1(eb2+pb1-eb1-1/v_s^2)-thb(eb2+pb1)-eb2 pb1]c^2 + thb[v_s^2(eb2+pb1-eb1)-1]=0`, require `0<c^2<1` | uncharged at rest, k→∞ | exact | Eq. (4.18) eq:cs0largek |
| 0.18 ⊙ | Moving-fluid gap-stability (necessary) | `v_s^2 eps1 + theta > gamma_s/(1-v_s^2)` | uncharged at v0≠0; excludes Landau-Lifshitz eps1=theta=0 | exact | Eq. (4.22) |
| 0.19 ⊙ | Small-k necessary stability | `gamma_s>0, eps1>0, theta>0` | uncharged at rest, small k; gapped omega=-i(eps0+p0)/(v_s^2 eps1), -i(eps0+p0)/theta | exact | Eq. (4.13) region eq:wgaps-2 |
| 0.20 ⊙ | Conformal (uncharged 3+1) stability | `1 - 3eta/theta - eta/pi1 > 0` AND `pi1 > 4eta` | v_s=1/sqrt3, zeta=0, eps1=3pi1; agrees w/ BDN | exact | Eq. (4.24) eq:cft-stable |
| 0.21 ⊙ | Bulk viscosity (general-frame) | `zeta = v_s^2(pi1 - v_s^2 eps1) - pi2 + v_s^2 eps2` | uncharged, frame-invariant | exact | Footnote fn:zeta; App. (A.16) |
| 0.22 ⊙ | Frame invariants / genuine-coeff count | `f_i=pi_i-(dp/deps)_n eps_i-(dp/dn)_eps nu_i; l_i=gamma_i-[n/(eps+p)]theta_i`; only eta,zeta,sigma genuine (3) | general frame; theta1=theta2, gamma1=gamma2 | exact | Eq. (2.11) + (2.15) |
| 0.23 ⊙ | Boost velocity addition c_v(phi) | `c_v(phi)=v0(1-c0^2)/(1-c0^2 v0^2)cosphi ± [c0/(1-c0^2 v0^2)]sqrt{(1-v0^2)[1-v0^2 c0^2-v0^2(1-c0^2)cos^2 phi]}` | mode omega=±c0|k| boosted to v0 | exact | Eq. (4.3) eq:cvphi |
| 0.24 ▽ | Fig.2 stability-region illustrative params | `eps2=0, pi1/gamma_s=3/v_s^2`; region grows for smaller v_s; origin excluded in causal panel | viz of (4.16)+(4.17)+(4.18) | illustrative | Fig. fig:soundstability-1 |
| 0.25 ▽ | Fig.5 sound-dispersion params | `v0=0.9, v_s=0.5, v_s^2 eps1/gamma_s=3, theta/gamma_s=4, eps2=0, pi1/gamma_s=3/v_s^2` | near stability boundary | illustrative | Fig. 5 |
| 0.26 ▽ | Fig.3/1/4 render params | Fig.3 v0=0.9, theta/eta=2; Fig.1 c0=1/2; Fig.4 v_s=1/2 at v0=0.1, 0.9 | render-only | illustrative | Figs. fig:cvphi etc. |

---

## DAG node: `step0.bdnk_recovery` — ideal-gas BDNK flat-space tests

### Paper arXiv:2209.09265 — `step0.bdnk_recovery`
Gamma-law ideal gas, Gamma=4/3 all figures. Bjorken + steady/dynamical shock + heat-flow.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 0.27 ⊙ | Bjorken inviscid solution | `eps=m n_0 tau^-1[1+e_0 tau^-(Gamma-1)]`; diagnostic `eps_dot+Gamma eps/tau = m n_0(Gamma-1)/tau^2` | Gamma=4/3, m=1, n_0=0.1, eps_0=0.25, eps_dot_0∈{-2,0,2}, tau∈[1,20], Vhat=1/10, sigma_hat=0 | plotting precision | eq:inviscid_bjorken / fig:bjorken |
| 0.28 ⭐ | **Bjorken RK4 conv. factor Q_N** | tau_hat=0.5: `34.8,18.7,16.9`; tau_hat=1: `18.4,16.9,16.3`; tau_hat=2: `16.9,16.3,16.1` (N=2^11); Q_N→16 | eps_dot=-2, RK4, 4th-order centered-FD residual of Bjorken EOM | within ~1 of 16 at finite N | Table table:ODE_conv |
| 0.29 ◻ | Bjorken max char. speed c_+ | tau_hat=0.5: ~1.3 (always superluminal); tau_hat=1: ~1.05 early, ~0.9 late; tau_hat=2: ~0.7 (always subluminal) | Gamma=4/3,m=1,Vhat=1/10,sigma_hat=0 | ~2 sig figs | Sec. Bjorken text + fig caption |
| 0.30 ◻ | Steady-shock right state for left {1,0.8,0.1} | `{eps_R,v_R,n_R} ≈ {4.439, 0.4143, 0.2929}` (RH-solved); task target e_R=4.4074, v_R=5/12 | Gamma=4/3, m=0.1, Vhat=2/15, sigma_hat=0, tau_hat=1.5 | ~1e-3 | fig:shockwave_profile + eq:Rankine_Hugoniot |
| 0.31 ◻ | RH left→right pairs (dynamical shock ICs) | `{1,0.9,1}_L => {11.5174,0.354727,5.44212}_R`; `{1,0.6,1}_L => {1.33795,0.514414,1.25027}_R` | Gamma=4/3, m=0.1, verbatim | 6 sig figs verbatim | eq:shockwave_params |
| 0.32 ⭐ | **Shockwave RK4 conv. factor Q_N** | `15.9, 15.9, 15.9` (N=2^13); Q_N→16 | RK4, 4th-order centered-FD residual of T^{tx}_{,x}=0 | ~15.9, →16 | Table table:ODE_conv |
| 0.33 ▽ | WENO/CWENO PDE conv. factor Q_N | Q_N→4 (2nd-order), plateau ~4 early until boundary (t~80 shock, t~150 heat) | FV Pandya:2022pif, TVD-RK2, CFL=0.1 | figure-read | Fig. fig:conv_plot |
| 0.34 ⊙ | Causality bound on tau_hat (simplified) | `tau_hat >= [(Gamma-1)(2-c_s^2)+c_s^2]/(1-c_s^2)`; use `>` for strict |c_+|<1 | frame eq:hydro_frame | exact | eq:fully_simplified_caus_const |
| 0.35 ⊙ | Char. speeds c_pm^2, c_1^2 | `c_pm^2=(c_s^2/2 tau_hat)(2 alpha-omega sigma_hat+tau_hat+1 ± [...]^{1/2})`; `c_1^2=c_s^2 eta/(V tau_hat)` | frame eq:hydro_frame | exact | eq:cpmsq, eq:c1sq |
| 0.36 ⊙ | Linear stability bound on sigma_hat | `sigma_hat <= 1/3` (full); simple `sigma_hat<=1/2` | frame eq:hydro_frame; 0<omega<3-2sqrt2, alpha>=1 | exact | eq:sigma_bound |
| 0.37 ▽ | Heat-flow nonzero-conductivity test | sigma_hat=0: eps_dot→0; sigma_hat=1/3: eps_dot→finite nonzero | Gamma=4/3, m=0.1, Vhat=2/15, tau_hat=1.5 | qualitative zero vs nonzero | Fig. fig:heat_stationary |

---

## DAG node: `s1c.hrsc_core` — conservative finite-volume BDNK scheme

### Paper arXiv:2201.12317 (Pandya, Most, Pretorius) — `s1c.hrsc_core`
Conformal BDNK (mu=0, c_s^2=1/3), flat-space FV code tests.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 0.38 ⊙ | Frame coeff chi_0 | `(25/4) eta_0 = 6.25 eta_0` | conformal, frame fixing char speeds=1 | exact | eq:frame_coeffs (l.455) |
| 0.39 ⊙ | Frame coeff lambda_0 | `(25/7) eta_0 = 3.571428... eta_0` | same frame | exact | eq:frame_coeffs (l.455) |
| 0.40 ⊙ | Shear viscosity scaling eta | `eta_0 eps^{3/4}` | conformal BDNK | exact | eq:frame_coeffs_def (l.441) |
| 0.41 ⊙ | KT max local speed a | `1` | char speeds exactly unity | exact | eq:KT_flux (l.895) |
| 0.42 ⊙ | WENO5 linear weights d_k | `(3/10, 3/5, 1/10)` | reconstruction to faces | exact | eq:WENO_weights (l.1742) |
| 0.43 ⊙ | CWENO derivative linear weights d_k | `(1/6, 2/3, 1/6)` | 4th-order derivative stencil | exact | l.826 |
| 0.44 ◻ | Steady-state shock left state | `eps_L=1, v_L=0.8` | 1D steady shock, eta_0=0.2, w=10, L=200, lambda=0.1 | exact input | l.1386 |
| 0.45 ⊙ | Steady-shock right state eps_R (RH) | `119/27 = 4.407407...` | from RH with eps_L=1, v_L=0.8 | exact analytic | eq:rankine_hugoniot (l.1360) |
| 0.46 ⊙ | Steady-shock right state v_R (RH) | `5/12 = 0.416666...` | from RH with v_L=0.8 | exact analytic | eq:rankine_hugoniot (l.1360) |
| 0.47 ▽ | Adaptive-PVR error → machine precision | `Delta_eta <~ 1e-7` → error ~1e-15 | eta_0=0.2 steady shock | order-of-mag (figure) | Fig. fig:steady_state |
| 0.48 ⊙/▽ | Conv. factor Q_N smooth limit | `4` (2nd-order); CN residual of grad_c T^{cx}, 1-norm | viscous rotor & KH at N=2^8,2^9,2^10 | →4 as h→0 | l.1827 / Fig. fig:conv_plot |
| 0.49 ▽ | Oblique-shock residual self-conv. scaling | N=2^9 ×4, N=2^10 ×16 (overlap) | 2nd-order at t=220 | 2nd order (×4/doubling) | Fig. fig:conv_plot |
| 0.50 ▽ | Discrete T^{tt} conservation (FV vs FD) | FV ~1e-15; Pandya-FD ~1e-3 (~12 orders worse) | 1D Gaussian, eta_0=0.2 | order-of-mag (figure) | Fig. fig:Tab_cons |
| 0.51 ⭐ | **Courant factor lambda (Table I)** | Gaussian 0.1/0.5; rotor 0.1/0.5; shock-tube 0.1/0.5; oblique 0.1/0.1; steady-shock 0.1/0.5; KH 0.5/0.5 | per-test (used/max) | exact | Table I (table:courant l.1039) |
| 0.52 ⊙ | Boundary ghost cells | 3 outermost per direction; X_{k,j}:=X_{3,j} (k∈[0,2]), X_{N-4,j} (k∈[N-3,N-1]) | ghost / periodic | exact | eq:ghost_cells (l.988) |

---

## DAG node: `s1a.*` — radial perturbations, causal viscous fluids (analytic)

### Paper arXiv:2506.09149 — `s1a.*`
ZERO figures/tables/numeric mode data (verified by grep). All targets analytic.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1a.1 ⊙ | Sufficient high-freq heat-conduction stability | `inf(L/P) + inf(K/W) >= 0` (infima over (0,R)); not sufficient for whole star | small visc O(alpha), kappa>=0, c_s^2∈[0,1], high-freq phi~A exp(ikr) | exact | Eq. eq:condi (l.1039-1041) |
| 1a.2 ⊙ | Nec.+suff. heat-conduction stability | `Re(i omega_j^{(1)}) = -[∫L|d2phi|^2+M|dphi|^2+N|phi|^2]/[∫P|dphi|^2+Q|phi|^2] - [∫K|dphi|^2+H|phi|^2]/[∫W|phi|^2]`; stable ⇔ RHS<0 ∀ eigenvectors | first-order small visc; all PF eigenvectors | exact | Eq. eq:heatint (l.978-989) |
| 1a.3 ⊙ | L, K heat-conduction operators | `L=-kappa((e+p)/n)^2 (e^{-Lam-4Phi}/r^2) c_n^2[c_s^2-c_n^2]`; `K=kappa((e+p)/n)^2 (e^{Lam-2Phi}/r^2)[c_s^2-c_n^2]` | static TOV bg, metric ds^2=-e^{-2Phi}dt^2+e^{2Lam}dr^2+r^2 dOmega^2 | exact | l.1052-1056; App B |
| 1a.4 ⊙ | Bulk+shear ALWAYS-stable | `<phi|F phi> = ∫ zeta(e^{Lam-2Phi}/r^2)|dphi|^2 + (4/3)eta(e^{Lam-2Phi}/r^2)|dphi-(3/r)phi|^2 >= 0` (any visc magnitude when kappa=0) | zeta,eta>=0; boundary term →0 | exact | Eq. eq:perfectsquare (l.866-869) |
| 1a.5 ⊙ | Heat stability sign logic | `K>=0 ⇔ c_n^2<=c_s^2`; `L>=0 ⇔ c_n^2<=0 OR c_n^2>=c_s^2`; suff(i): c_n^2<=0 everywhere; suff(ii): eq:condi | kappa>=0, c_s^2∈[0,1] | exact | l.1058-1071, 1497-1501 |
| 1a.6 ⊙ | BDNK fast-timescale eigenfreqs | `Omega_j^± = (i/2)e^{-Phi}[1/tau_E+1/tau_Q] ± e^{-Phi}sqrt(e^{-2Lam}c_tau^2 lambda_j - (1/4)[1/tau_E-1/tau_Q]^2)`, lambda_j=((j+1/2)pi alpha/R)^2 | fast timescale tau=t/alpha; BCs d_rho psi_E(0)=0, psi_E(R/alpha)=0 | exact | Eq. eq:fasteig (l.1345-1348) |
| 1a.7 ⊙ | Second-sound speed + causality (BDNK) | `c_tau^2 = tau_P/tau_E`; speed `e^{-Lam-Phi}sqrt(tau_P/tau_E)`; causality `tau_P < tau_E`; all fast modes decaying | BDNK fast timescale | exact | c_tau^2 l.1330; constraint l.1332-1333 |
| 1a.8 ⊙ | MIS fast-timescale behavior | no oscillation/rescaling; all fast modes strictly decaying, stable for (tau_0,tau_1,tau_2)>0 | Maxwell-Cattaneo MIS | exact | l.1397-1407, 1456 |
| 1a.9 ⊙ | Mathews ideal-gas instance (illustrative) | `p=(1/3)(e^2-(mn)^2)/e; c_n^2=(1/3)(mn/e)^2>=0; c_s^2-c_n^2=-(2/3)(n/(e+p))(mn/e)m<=0 ⇒ L>=0, K<=0` | particle mass m; sign illustration only | exact | l.1077-1090 |

---

## DAG node: `s1b.*` / `s1b.axial_wave_eqs` — axial oscillations & perturbation formalism

### Paper arXiv:2604.13208 ("Axial Oscillations of Viscous Neutron Stars") — `s1b.*`
Polytrope EOS1/EOS2; w-mode QNM tables (frame A/B). **Richest table set in the program.**

Stellar models (ρ_c=3e15 g/cm³):

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1b.1 ◻ | M (EOS1) | `1.27 M_sun` | EOS1 κ=100 km² n=1, ρ_c=3e15 | 3 sig figs | text after eq:EoS |
| 1b.2 ◻ | R (EOS1) | `8.86 km` | EOS1, ρ_c=3e15 | as printed | text after eq:EoS |
| 1b.3 ⭐ | M/R (EOS1) | `0.21` | EOS1, ρ_c=3e15 | as printed | caption Table tab:modes_frames |
| 1b.4 ◻ | M (EOS2) | `1.54 M_sun` | EOS2 κ=700 km^2.5 n=0.8, ρ_c=3e15 | as printed | text after eq:EoS |
| 1b.5 ◻ | R (EOS2) | `8.78 km` | EOS2, ρ_c=3e15 | as printed | text after eq:EoS |

ℓ=2 fundamental w-mode (f [kHz], τ [µs]), EOS1, ρ_c=3e15, frames [A1, A2, B1, B2] — **Table tab:modes_frames (6 sig figs)**:

| # | η_c [cgs] | A1 | A2 | B1 | B2 | Tol | Source |
|---|-----------|-----|-----|-----|-----|-----|--------|
| 1b.6 ⭐ | 3e29 | (10.4884, 29.5870) | (10.4884, 29.5870) | (10.4868, 29.5894) | (10.4868, 29.5891) | as printed | Table tab:modes_frames |
| 1b.7 ⭐ | 5e29 | (10.4795, 29.6169) | (10.4794, 29.6167) | (10.4769, 29.6194) | (10.4768, 29.6200) | as printed | Table tab:modes_frames |
| 1b.8 ⭐ | 8e29 | (10.4661, 29.6619) | (10.4660, 29.6608) | (10.4622, 29.6639) | (10.4619, 29.6658) | as printed | Table tab:modes_frames |
| 1b.9 ⭐ | 1e30 | (10.4571, 29.6917) | (10.4571, 29.6898) | (10.4523, 29.6938) | (10.4522, 29.6964) | as printed | Table tab:modes_frames |
| 1b.10 ⭐ | 3e30 | (10.3692, 29.9752) | (10.3705, 29.9733) | (10.3564, 29.9658) | (10.3539, 29.9699) | as printed | Table tab:modes_frames |
| 1b.11 ⭐ | 5e30 | (10.2854, 30.2463) | (10.2874, 30.2467) | (10.2783, 30.1107)† | (10.2625, 30.2600) | as printed | Table tab:modes_frames |
| 1b.12 ⭐ | 8e30 | (10.1659, 30.6362) | (10.1687, 30.6394) | (10.1293, 30.5004)† | (10.1443, 30.6579) | as printed | Table tab:modes_frames |
| 1b.13 ⭐ | 1e31 | (10.0898, 30.8857) | (10.0932, 30.8905) | (10.0608, 30.7400)† | (10.1271, 30.8477) | as printed | Table tab:modes_frames |

† B1 shows mode-avoidance dip (η_c >= 5e30).

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1b.14 ◻ | PF (η→0) fundamental ℓ=2 w-mode | `(10.50, 29.54)` (f_0, τ_0) | EOS1, ρ_c=3e15; <0.1% vs Kokkotas priv. comm | sub-percent | text/comment l.543 (commented), benchmark l.542 |
| 1b.15 ◻ | Compactness slope Δf/f_0 vs M/R | `-1.8` per unit M/R | both EOS, η_c=5e29, frame A2; intercept EOS-dep (not given) | 2 sig figs | Eq. l.569 |
| 1b.16 ◻ | Compactness slope Δτ/τ_0 vs M/R | `-5.0` per unit M/R | both EOS, η_c=5e29, frame A2 | as printed | Eq. l.569 |
| 1b.17 ⭐ | Frame table A1 max c_s² | `0.83` | frame A1, A, τ̂=10 | as printed | Table tab:frames |
| 1b.18 ⭐ | Frame table A2 max c_s² | `0.91` | frame A2, A, τ̂=20 | as printed | Table tab:frames |
| 1b.19 ⭐ | Frame table B1 max c_s² | `0.8` | frame B1, B, τ̂=10 | as printed | Table tab:frames |
| 1b.20 ⭐ | Frame table B2 max c_s² | `0.9` | frame B2, B, τ̂=20 | as printed | Table tab:frames |
| 1b.21 ▽ | η-mode onset (figure only) | Im ω→0 as η_c→0; kHz freq, ms-scale damping; mode avoidance, never cross | EOS1, ρ_c=3e15, two frames | figure-only | Fig. fig:eta_modes |
| 1b.22 ▽ | UCO viscous damping bound | `|Im ω_ℓ| ≲ 1e-2` (geometric M=1); independent of ℓ, compactness; freqs unshifted (≤ few %) | const-density, frame B, η_c=1e31 | figure / order-of-mag | Fig. fig:ucos |
| 1b.23 ⊙ | KSS shear-viscosity estimate (UCO) | `5.4e31 (M_sun/M) g cm^-1 s^-1` | BH-mimicker entropy estimate | order-of-mag | unnumbered Eq. l.222 |

### Paper arXiv:2411.16841 ("Perturbations of relativistic dissipative stars") — `s1b.axial_wave_eqs`
Formalism paper; NO numeric tables/figures/QNMs. All targets symbolic master equations.

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1b.24 ⊙ | Axial inviscid master eq (RW) | `-d2psi/dt2 + d2psi/drstar2 - (e^Phi/r^2)(lambda^2-6M/r+4pi r^2(E-P)) psi = 0`; psi=r^-1 e^{Phi/2} k_n | axial/odd, PF limit, ℓ>=2; reduces to RW outside (E=P=0) | exact symbolic | eq:Master_Odd_Inviscid (~l.514) |
| 1b.25 ⊙ | Axial viscous GW master eq (coupled 1/2) | `... = 16 pi eta e^{Phi/2} dpsi/dt - (8 pi e^Phi/r)(2 deta/drstar + eta dPhi/drstar) beta` | axial BDNK shear eta(r); stationary bg mu=U=0 | exact symbolic | eq:Master_Odd_Viscous (~l.532) |
| 1b.26 ⊙ | Axial viscous fluid eq (NOVEL viscous mode, coupled 2/2) | `-tau_Q d2beta/dt2 + (eta/(E+P)) d2beta/drstar2 + eta(b1 dbeta/drstar+b2 dbeta/dt+b3 beta) = c1 d2psi/dtdrstar+c2 dpsi/dt+c3 dpsi/drstar+c4 psi` | b_i,c_i in SM; K=8 pi r e^{(Phi+Lam)/2}(E+P); tau→tau_Q | exact symbolic | eq:Master_Odd_Fluid_Viscous (~l.549) |
| 1b.27 ⊙ | Axial char. speeds | `c_GW=1; c_Viscous=sqrt(eta tau_Q/(E+P))` | principal-part of coupled axial pair | exact | text ~l.560 |
| 1b.28 ⊙ | Regge-Wheeler potential V_RW (axial) | `V_RW=(e^Phi/r^2)(lambda^2-6M/r+4pi r^2(E-P))`, lambda^2=ℓ(ℓ+1) | axial/odd interior; exterior E=P=0 | exact | from eq:Master_Odd_Inviscid (~l.514) |
| 1b.29 ⊙ | Even-sector char. speeds + BDNK bounds | `c2_GW=1 (mult 2); c2_Viscous=eta/(rho tau_Q); c2_pm=(C1±C2)/(2 rho tau_E tau_Q)`; bounds tau_Q,tau_E,tau_P>0, 0<=eta<=rho tau_Q, C1^2>=4 c_s^2 rho tau_E tau_Q^2(rho tau_P-V), 0<=C1+C2<=2 rho tau_E tau_Q | even full viscous, 5 vars; recovers Bemfica:2020zjp | exact | eqs ~l.738-759 |

---

## DAG node: `s1c.*` — full BDNK neutron-star evolution

### Paper arXiv:2509.15303 — `s1c.*`
Closed-form EOS (κ=100, Γ=2); Cowling BDNK spherical evolution; f-mode QNM freq + decay tables. **Primary s1c validation set.**

Stellar config:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1c.1 ◻ | Total grav. mass M_T | `1.4 M_sun` | ρ_{0,c}=0.00128 (eps_c=0.00144), EoS κ=100 Γ=2 | 2 sig figs | text l.597 |
| 1c.2 ◻ | Max c_s in star | `~0.45 c` | NS center, c_s^2<1/3 satisfied | as stated | text l.653 |

QNM frequencies [kHz], Dr=0.002, t_f=8000, Blackman FFT — **Table QNM_freq_table** (±0.01):

| # | Mode | PF | smallSB-F2 | highB-F9 | Tol | Source |
|---|------|-----|-----------|----------|-----|--------|
| 1c.3 ⭐ | F (f-mode) | 2.69 | 2.69 | 2.67 | ±0.01 | Table QNM_freq_table |
| 1c.4 ⭐ | H1 | 4.55 | 4.60 | 4.60 | ±0.01 | Table QNM_freq_table |
| 1c.5 ⭐ | H2 | 6.36 | 6.36 | 6.30 | ±0.01 | Table QNM_freq_table |

f-mode decay rates [M_sun^-1], Dr=0.002 — **Table QNM_decay_table**:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1c.6 ⭐ | 1/tau_l smallSB-F2 | `0.00157` | linear fit | ±0.00001 | Table QNM_decay_table |
| 1c.7 ⭐ | 1/tau_nl smallSB-F2 | `0.00157` | nonlinear fit | — | Table QNM_decay_table |
| 1c.8 ⭐ | 1/tau_l medS-F2 | `0.00150` | linear fit | — | Table QNM_decay_table |
| 1c.9 ⭐ | 1/tau_l highB-F9 | `0.00215` | linear fit | — | Table QNM_decay_table |
| 1c.10 ⭐ | 1/tau_l medSB-F9 | `0.00182` | linear fit | — | Table QNM_decay_table |
| 1c.11 ⭐ | omega_nl (all four cases) | `0.0834 M_sun^-1` (= f=2.71 kHz) | nonlinear fit | — | Table QNM_decay_table |

Decay rate vs Dr [M_sun^-1], Dr=0.0032,0.0028,0.0024,0.0020 — **Table tab:QNM_table**:

| # | Case | Values (Dr decreasing) | Tol | Source |
|---|------|------------------------|-----|--------|
| 1c.12 ⭐ | smallSB-F2 | 0.0019, 0.0018, 0.0017, 0.0016 | as printed | Table tab:QNM_table |
| 1c.13 ⭐ | medS-F2 | 0.0018, 0.0017, 0.0016, 0.0015 | as printed | Table tab:QNM_table |
| 1c.14 ⭐ | highB-F9 | 0.0024, 0.0024, 0.0023, 0.0022 | as printed | Table tab:QNM_table |
| 1c.15 ⭐ | medSB-F9 | 0.0021, 0.0020, 0.0019, 0.0018 | as printed | Table tab:QNM_table |
| 1c.16 ⭐ | PF | 0.00023, 0.00021, 0.00019, 0.00018 | as printed | Table tab:QNM_table |

Continuum extrapolation — **Table tab:QNM_table**:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1c.17 ⭐ | Continuum 1/tau_0 | smallSB-F2=0.0011, medS-F2=0.0010, highB-F9=0.0017, medSB-F9=0.0013 [M_sun^-1] | extrapolation eq, PF→0; highB-F9 uses 3 highest-res | as printed | Table tab:QNM_table |
| 1c.18 ⭐ | Continuum 1/tau_0 (SI) | smallSB-F2=220, medS-F2=200, highB-F9=350, medSB-F9=260 [s^-1] | extrapolated physical units | as printed | Table tab:QNM_table |
| 1c.19 ◻ | Convergence order p | BDNK p=1 (marginal), PF p=0.54 | fit of extrapolation eq; scheme expected order 3 | as printed | text l.759, 794 + Table |
| 1c.20 ⭐ | f-mode freq convergence smallSB-F2 | F=2.69/2.69/2.67, H1=4.60/4.60/4.61, H2=6.36/6.36/6.33 [kHz] | Dr=0.0028/0.002/0.001 | as printed | Table QNM_freq_convergence |

Numerical-method targets:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 1c.21 ◻ | CFL Courant factor | `0.25` (Dt/Dr) | SSP-RK3 | exact | text l.583 |
| 1c.22 ◻ | Box size r_max | `20 M_sun` | staggered grid, outflow BC | exact | text l.586 |
| 1c.23 ◻ | Resolution range | `0.001 - 0.0032 M_sun` | high-res runs | as printed | text l.586 |
| 1c.24 ◻ | Min num. dissipation speed | `0.1 c` | FDOC floor near surface | exact | text l.584 |
| 1c.25 ◻ | Atmosphere threshold rho_{0,atms} | `1e-12 M_sun^-2` (set 1e-13) | atmosphere | as stated | text l.584 |
| 1c.26 ◻ | Stable-window boundary | `tau_eps <~ 0.1` | empirical, high-res Dr=0.001 | empirical | text l.615 |

---

## DAG node: `s2.is_contrast` — Israel-Stewart/MIS bulk-viscosity contrast

### Paper arXiv:2311.13027 (Chabanov & Rezzolla) — `s2.is_contrast`
Hybrid polytrope; MIS/IS bulk viscosity; oscillating-TOV / migration / BNS escalating tests.

Reynolds-number + stellar models:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 2.1 ◻ | R^{-1}_min | `≃ -0.1` | Γ_id=2, κ=100, ρ=0.00128, σ=-0.9; limiter (i) | order-of-mag (≃) | Eq. R^{-1}_min=σp/(e+p) |
| 2.2 ◻ | R^{-1}_max | `≃ 0.8` | Γ_id=2, κ=100, ρ=0.00128; limiter (ii) | ≃ | Eq. R^{-1}_max=(e-p)/(e+p) |
| 2.3 ◻ | Osc-TOV ADM mass | `1.4 M_sun` | κ=100, Γ=2, Γ_th=1.1, ρ_c=1.28e-3 | exact | Sec. num_vis text |
| 2.4 ◻ | Osc-TOV radius R | `14.2 km` | same | ±0.1 km | Sec. num_vis text |
| 2.5 ◻ | Osc-TOV central density | `1.28e-3 M_sun^-2 (≈7.91e14 g/cm³)` | κ=100,Γ=2,Γ_th=1.1 | exact | Sec. num_vis text |
| 2.6 ◻ | Input ζ_h set (osc) | `[0, ~9.42e25, ~1.98e26, ~8.20e26]` g s^-1 cm^-1 | osc TOV M=1.4 | as listed | Sec. num_vis text |

Fit coefficients — **Tables tab:fits / tab:fits_old (verbatim)**:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 2.7 ⭐ | Fit ζ_a/ζ_s/p_h (λ≈3.48e4 m) | zero: ζ_s=4.4e29, p_h=1.79; low: ζ_a=9.89e25, ζ_s=8.87e32, p_h=3.57; med: ζ_a=1.6e26, ζ_s=5.65e32, p_h=3.45; high: ζ_a=5.36e26, ζ_s=2.75e32, p_h=3.23 | M=1.4 TOV, fit eq:numvisfit | verbatim | Table tab:fits |
| 2.8 ⭐ | Fit coeffs (λ≈1.74e4 m) | zero: ζ_s=1.66e33, p_h=3.62; low: ζ_s=3.49e32, p_h=3.23; med: ζ_a=8.87e24, ζ_s=1.04e31, p_h=2.32; high: ζ_a=3.66e26, ζ_s=8.29e31, p_h=2.83 | M=0.7 TOV | verbatim | Table tab:fits_old |
| 2.9 ◻ | Measured vs input ζ agreement | low ≲6% / med ≲20% / high ≲35% | osc test | upper bounds | Sec. num_vis text |
| 2.10 ◻ | Numerical viscosity floor | `≲1e26 g cm^-1 s^-1` | M=1.4 TOV | upper bound | Conclusions text |
| 2.11 ◻ | Small-TOV M / R (App B) | `M=0.7 M_sun, R=7.1 km` | κ=25, Γ_th=2, ρ_c=5.12e-3 | exact / ±0.1 km | App. small_star text |

Migration test:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 2.12 ◻ | Migration initial ρ_c (unstable) | `8e-3 M_sun^-2 (≈4.94e15 g/cm³)` | κ=100, Γ=Γ_th=2 | exact | Sec. mig text |
| 2.13 ◻ | Migration star mass M | `1.447 M_sun` | same | exact | Sec. mig text |
| 2.14 ◻ | Migration bulk viscosities | `[0, 4e28, 2e29, 1e30]` g s^-1 cm^-1 | migration test | exact | Sec. mig text |
| 2.15 ▽ | Migration p_th/p (high-visc) | `up to ~10%` | high-visc migration | figure-read / text bound | Conclusions + Fig. mig_central |
| 2.16 ◻ | Migration |R^{-1}| threshold | `≳ 1e-2` | central inverse Reynolds | ~ | Sec. mig text / Fig. mig_central |
| 2.17 ▽ | Migration 2nd-shock stall radius | `~60 km` | high-visc, 2nd shock | figure-read | Sec. mig text / Fig. mig_shocks |
| 2.18 ◻ | Migration density contour levels | `4.5 × [1e14,1e12,1e10,1e8] g/cm³` | Fig. contours; inner=ρ_h, outer=ρ_l | exact | Fig. mig_shocks |

BNS merger:

| # | Quantity | Target | Conditions | Tol | Source |
|---|----------|--------|------------|-----|--------|
| 2.19 ◻ | BNS ζ0 (high-visc ref) | `1e30 g cm^-1 s^-1` | ζ_h∈[ζ0,ζ0/2,ζ0/5,0] | exact | Sec. binary text |
| 2.20 ◻ | BNS ρ_h, ρ_l, τ_h | ρ_h≈4.52e14, ρ_l≈1.13e12 g/cm³; τ_h≈2.7e-4 ms | BNS setup | ≈ | Sec. binary text |
| 2.21 ◻ | BNS R^{-1} / |Π|/p peak | R^{-1}>~1%; |Π|/p~20% | high-visc, center, persists to ~1.5 ms | ~ | Sec. binary + Fig. bin_central |
| 2.22 ▽ | BNS |Π|/ρh, |Π|/p at ~5 ms | (|Π|/ρh)_c≈1e-5; (|Π|/p)_c≈1e-4 | high-visc center, after ~2-order drop | figure-read | Sec. binary / Fig. bin_central |
| 2.23 ▽ | BNS ejecta suppression | high-visc ejecta ≈20% of zero-visc (factor ~5) | Bernoulli, detector ~517 km | ~ | Abstract + Fig. outflow_tot |
| 2.24 ◻ | BNS snapshot times | bin_temp_pip 4.22 ms; rho_temp_xz 23 ms; bin_bernoulli 7.11 ms; HMNS end ~24 ms | BNS figures | exact | Figs. (various) |
| 2.25 ◻ | BNS detector / GW extraction / contour | detector ~517 km; GW ~740 km @100 Mpc; Mollweide contour 8e-7 M_sun^-2 | BNS outflow + GW | ~ / exact | Figs. outflow/waves |
| 2.26 ◻ | Limiter parameter σ | `-0.9` | Π≥σp limiter (i), -1<=σ<=0 | exact (typical) | Sec. limiting (i) |

---

## ⭐ Reproduce-first set (table-sourced, highest-confidence numeric targets)

Ordered roughly by DAG dependency (step0 frame checks → s1 mode tables → s2 fit tables). These are the directly-tabulated numbers; matching them validates each node before tackling text/figure targets.

**step0.bdnk_recovery** (convergence sanity, gates the scheme)
- 0.28 Bjorken RK4 Q_N → 16 (Table table:ODE_conv)
- 0.32 Shockwave RK4 Q_N = 15.9 → 16 (Table table:ODE_conv)

**s1c.hrsc_core** (FV code config)
- 0.51 Courant factor λ per-test (Table I)

**s1b.\*** (axial w-mode tables — the richest exact set)
- 1b.6–1b.13 ℓ=2 w-mode (f,τ) for η_c=3e29…1e31 across frames A1/A2/B1/B2 (Table tab:modes_frames, 6 sig figs)
- 1b.3 M/R(EOS1)=0.21; 1b.17–1b.20 frame max c_s² {0.83, 0.91, 0.8, 0.9} (Table tab:frames)

**s1c.\*** (full-BDNK NS evolution — primary deliverable)
- 1c.3–1c.5 QNM freqs F/H1/H2 for PF/smallSB-F2/highB-F9 (Table QNM_freq_table, ±0.01)
- 1c.6–1c.11 f-mode decay rates + omega_nl=0.0834 (Table QNM_decay_table)
- 1c.12–1c.16 decay rate vs Dr, all five cases (Table tab:QNM_table)
- 1c.17–1c.18 continuum 1/tau_0 (geometric + SI) (Table tab:QNM_table)
- 1c.20 f-mode freq convergence smallSB-F2 (Table QNM_freq_convergence)

**s2.is_contrast** (IS/MIS bulk-viscosity calibration)
- 2.7 Fit ζ_a/ζ_s/p_h, λ≈3.48e4 m (Table tab:fits, verbatim)
- 2.8 Fit ζ_a/ζ_s/p_h, λ≈1.74e4 m (Table tab:fits_old, verbatim)

> Note: `step0.causality` (2009.11388, 1907.08191), `s1a.*` (2506.09149), and `s1b.axial_wave_eqs` (2411.16841) have **no table/figure numeric targets** — they are validated symbolically against equations. Only 2009.11388's non-emptiness unit-test (0.12) carries concrete numbers, and it is the canonical first analytic check.

---

## Summary counts

- Papers merged: **9**
- DAG nodes: 6 distinct (step0.causality ×2 papers, step0.bdnk_recovery, s1c.hrsc_core, s1a.\*, s1b.\* + s1b.axial_wave_eqs, s1c.\*, s2.is_contrast)
- Total exact targets cataloged: **111**
- Table-sourced (⭐) targets flagged to reproduce first: see set above
