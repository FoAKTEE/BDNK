# Source: c-1D_conformal_bdnk (Alex Pandya, 1D conformal BDNK C solver)

LOCAL C code: `/data/haiyangw/claude/BDNK/ref-code/1D_conformal_bdnk/`
Files: `solver.c` (21 KB, the everything-file), `parameters.h`, `makefile`, `plot_Ttt.py`, `README.md`.
Algorithm reference: **arXiv:2201.12317** (Pandya, Most, Pretorius — "Conservative finite volume scheme for first-order viscous relativistic hydrodynamics"). README author Alex Pandya, dated April 3-4 2024.

Problem solved: **conformal** BDNK in **flat spacetime**, Cartesian, **slab symmetry** (invariance under y,z translations) -> effectively 1+1D. Conformal EOS so `p = eps/3`. Primitive vars are `xi = ln(eps)` and `ux` (the x-component of 4-velocity; `ut = sqrt(1+ux^2)`, signature +/-... with `u^mu u_mu`). Conserved vars evolved: `Ttt`, `Ttx`.

---

## 1. Hydrodynamic frame coefficients (parameters.h:33-37; ref eqs 19-20 of 2201.12317)

Conformal frame: all transport coefficients scale as `eps^{1/4}` (=`(exp xi)^0.25`).

```c
#define EPS0 (10.)
const double eta0    = 1.*pow(EPS0,0.25)/(3.*M_PI);   // = EPS0^0.25 / (3 pi)
const double lambda0 = 25./7.*eta0;                   // lambda0 = (25/7) eta0
const double chi0    = 25./4.*eta0;                   // chi0    = (25/4) eta0
```
So with EPS0=10: `eta0 = 10^0.25/(3 pi) = 1.77828/9.42478 = 0.188682`.
- `eta0`  = shear viscosity scale (eta = eta0 * eps^{1/4})
- `lambda0` = energy-flux relaxation coeff scale (the "lambda" frame coeff; multiplies Qx)
- `chi0`  = energy-density relaxation coeff scale (the "chi" frame coeff; multiplies A)

In the constitutive tensor these enter as `eta0*eps^{0.75}`, `lambda0*eps^{0.75}`, `chi0*eps^{0.75}` (i.e. coefficient ~ eps * eps^{-1/4} because they multiply gradients of xi=ln eps; the `eps^{0.75}=(exp XI)^0.75` factor appears explicitly in compute_A/Qx/m2sxx).

Ratios used in primitive recovery: `ch = chi0/eta0 = 25/4 = 6.25`, `l = lambda0/eta0 = 25/7 = 3.571428...`.

**Causality/stability:** The specific ratios chi0=(25/4)eta0, lambda0=(25/7)eta0 are the canonical causal+stable choice from 2201.12317 for the conformal case. The code does NOT re-derive or check the inequalities at runtime; it just hardcodes these values. (Causality conditions in BDNK conformal: roughly chi0 >= 4 eta0 and lambda0 >= ... ; the chosen values satisfy them. Exact inequality forms are in the paper, not in the code — see open questions.)

---

## 2. Constitutive tensor pieces (solver.c:287-307)

These build the first-order dissipative corrections. Arguments: `XI=xi`, `UX=ux`, `xiP=d_x xi`, `uxP=d_x ux`, `xiD=d_t xi` (xi-dot), `uxD=d_t ux` (ux-dot). `ut = sqrt(1+UX^2)`.

**compute_A** (scalar correction to energy density, the "A" / delta-eps piece, coeff chi0):
```
A = (chi0 * (exp XI)^0.75 * ( 4 UX uxD + 4 sqrt(1+UX^2) uxP
        + 3 (1+UX^2) xiD + 3 UX sqrt(1+UX^2) xiP )) / (4 sqrt(1+UX^2))
```

**compute_Qx** (heat/energy flux x-component, coeff lambda0):
```
Qx = ((exp XI)^0.75 * lambda0 * ( 4 sqrt(1+UX^2) uxD + 4 UX uxP
        + UX sqrt(1+UX^2) xiD + (1+UX^2) xiP )) / 4
```

**compute_m2sxx** ("minus 2 eta sigma_xx" = -2*eta*shear_xx, the shear-stress xx piece, coeff eta0):
```
m2sxx = ( -4 (exp XI)^0.75 * eta0 * sqrt(1+UX^2) * ( UX uxD + sqrt(1+UX^2) uxP ) ) / 3
```

### Full stress tensor assembly (solver.c:309-402)
`eps = exp(xi)`. Other components built from A, Qx, m2sxx via boost relations:
- `Qt = ux*Qx/ut`
- `m2etasigmatx = ux*m2sxx/ut`
- `m2etasigmatt = ux*m2etasigmatx/ut = ux^2*m2sxx/ut^2`

**T_tt** (solver.c:309): `Deltatt = -1 + ut^2`
```
T_tt = eps*(ut^2 + Deltatt/3) + A*(ut^2 + Deltatt/3) + 2*Qt*ut + m2etasigmatt
```
**T_tx** (solver.c:341): `Deltatx = ut*ux`
```
T_tx = eps*(ut*ux + Deltatx/3) + A*(ut*ux + Deltatx/3) + Qt*ux + ut*Qx + m2etasigmatx
```
**T_xx** (solver.c:361): `Deltaxx = 1 + ux^2`
```
T_xx = eps*(ux^2 + Deltaxx/3) + A*(ux^2 + Deltaxx/3) + 2*Qx*ux + m2sxx
```
Structure = `(eps+A)*(u^mu u^nu + Delta^{mu nu}/3) + Q^mu u^nu + Q^nu u^mu + (-2 eta sigma^{mu nu})`.
Conformal so the bulk/trace correction is folded into the eps/3 (pressure = eps/3) plus A.

**Perfect-fluid versions** T_tt0/T_tx0/T_xx0 (solver.c:329,380,392): just drop A,Qx,m2sxx:
```
T_tt0 = eps*(ut^2 + Deltatt/3);  T_tx0 = eps*(ut*ux + Deltatx/3);  T_xx0 = eps*(ux^2 + Deltaxx/3)
```
(T_tt0 used in initial data; the *0 funcs are declared but T_tt/T_tx are what set the conserved ID, with all dissipative args =0 so they reduce to PF anyway.)

---

## 3. Primitive recovery: compute_xiD / compute_uxD (solver.c:404-497) — THE KEY ROUTINE

This is the crux of BDNK first-order: there is NO algebraic cons->prim. Instead you solve for the **time derivatives** `xiD = d_t xi` and `uxD = d_t ux` given the current `(xi, ux, d_x xi, d_x ux)` and the conserved `(Ttt=T00, Ttx=T01)`. Then you time-integrate xi, ux directly (method-of-lines on the primitives). The conserved Ttt,Ttx are evolved by their own conservation-law fluxes; the prim time-derivatives are obtained by *inverting* the constitutive relations T_tt(...,xidot,uxdot)=T00 and T_tx(...,xidot,uxdot)=T01 for (xidot,uxdot).

### Step-by-step (identical preamble in both compute_xiD and compute_uxD):

`ch = chi0/eta0 = 6.25`, `l = lambda0/eta0 = 25/7`.

1. **Compute the perfect-fluid-plus-spatial-gradient pieces** `TttPF`, `TtxPF` — these are what T00,T01 would be if xidot=uxdot=0 (i.e. all time derivatives in the stress tensor set to zero, keeping spatial gradients):
```
TttPF = (9 e + 8 e UX^4 + 6 UX^2 (3 e - 2 e^0.75 eta0 uxP) + 3 e^0.75 eta0 UX^3 xiP)
        / (9 + 6 UX^2)
TtxPF = (UX sqrt(1+UX^2) (8 e UX^2 + 12 (e - e^0.75 eta0 uxP) + 3 e^0.75 eta0 UX xiP))
        / (9 + 6 UX^2)
```
   where `e = exp(XI)`, `e^0.75 = pow(exp(XI),0.75)`.

2. **Form the "deficit" that must be supplied by the time derivatives**, normalized by eta0:
```
if (eta0 != 0):  dtt = (TttPF - T00)/eta0;   dtx = (TtxPF - T01)/eta0;
else:            dtt = 0; dtx = 0;            // guard against NaN
```
   So `dtt`, `dtx` measure how far the actual conserved (T00,T01) are from the no-time-derivative stress. They are the driving terms for xidot,uxdot.

3. **Perfect-fluid fallback test (role of TOL):**
```
if (eta0*fabs(dtt) < TOL || eta0*fabs(dtx) < TOL) {
    dtt = 0; dtx = 0;
    eps = -T00 + sqrt(4*T00^2 - 3*T01^2);          // PF algebraic cons->prim
    xi[n][i] = log(eps);
    ux[n][i] = 3*T01 / sqrt((3*T00+eps)^2 - (3*T01)^2);
    VISC[0][i] = 0;   // (only set in compute_xiD) flag: PF solve used in this cell
} else {
    VISC[0][i] = 1;   // (only in compute_xiD) flag: BDNK solve used
}
```
   - `TOL` is set in parameters.h (`#define TOL (-1)`). Because `eta0*fabs(...)` >= 0 always, `... < TOL` is FALSE whenever TOL<0. **TOL<0 => fallback NEVER triggers => always use BDNK primitive solve.** This is the default/intended production setting (comment lines 92-101).
   - If TOL>0, then in cells where the viscous correction is tiny (`eta0*|dtt| < TOL` or `eta0*|dtx| < TOL`), it overwrites xi,ux in place with the **perfect-fluid algebraic inversion** and zeros the deficits so xidot,uxdot become pure-spatial. This is a stability crutch for under-resolved cells; comment says ideally only used when BDNK solve is unstable at low res.
   - **PF cons->prim formula (conformal, p=eps/3):**
     `eps = -Ttt + sqrt(4 Ttt^2 - 3 Ttx^2)`,  `ux = 3 Ttx / sqrt((3 Ttt + eps)^2 - (3 Ttx)^2)`.
     (Derivable from inverting T_tt0,T_tx0; note discriminant `4 T00^2 - 3 T01^2` must be > 0.)

4. **Return the analytic solution for the time derivative.** Common denominator:
```
DEN = e^0.75 * (9 ch l + 12 ch (l-1) UX^2 + 4 (ch(l-3) - l) UX^4)
```
   **compute_xiD returns** `xiD = d_t xi`:
```
xiD = -2 (2 uxP + UX (1+UX^2) xiP) / ( sqrt(1+UX^2) (3 + 2 UX^2) )
      - 4 sqrt(1+UX^2) (3 l + (-4 + 4 ch + 6 l) UX^2) * dtt / DEN
      + 4 (3 (ch + 2 l) UX + (-4 + 4 ch + 6 l) UX^3) * dtx / DEN
```
   **compute_uxD returns** `uxD = d_t ux`:
```
uxD = - sqrt(1+UX^2) (8 UX uxP + 3 xiP) / (12 + 8 UX^2)
      + 3 UX sqrt(1+UX^2) (4 ch + l + 2 (2 ch + l) UX^2) * dtt / DEN
      - 3 (1+UX^2) (3 ch + 2 (2 ch + l) UX^2) * dtx / DEN
```
   First term in each = the "spatial-advection" part (present even at zero dissipation); the dtt/dtx terms = the viscous inversion that couples back the conserved-vs-PF deficit.

**Interpretation:** This is an explicit 2x2 linear solve of `[dT_tt/d(xidot), dT_tt/d(uxdot); dT_tx/d(xidot), dT_tx/d(uxdot)] [xidot; uxdot] = [T00 - (stuff); T01 - (stuff)]` done by hand (Mathematica-generated), so it is O(1) per cell, no Newton iteration. The denominator DEN is the determinant of that Jacobian times e^0.75; it can vanish/sign-change -> potential breakdown, which is why the TOL/PF fallback exists.

---

## 4. Spatial discretization, flux, reconstruction (solver.c:62-181, 499-560)

### WENO5 reconstruction (solver.c:62-136)
Fifth-order WENO with smoothness indicators beta and the WENO epsilon = `epsW = 1e-3` (parameters.h:85).
- `WENO_reconst_qRx(q,n,i)` returns reconstructed value **at x=i+1/2 from the left/+ stencil shifted by i+=1** (linear weights d0=1/10, d1=3/5, d2=3/10). Comment says "returns qL at x=i+1/2" (left state at interface, despite name R).
- `WENO_reconst_qLx(q,n,i)` returns the other interface state (linear weights d0=3/10,d1=3/5,d2=1/10). Comment "returns qR at x=i+1/2".
- Both use identical Jiang-Shu beta indicators:
```
beta0 = 1/4 (3 q_i - 4 q_{i+1} + q_{i+2})^2 + 13/12 (q_i - 2 q_{i+1} + q_{i+2})^2
beta1 = 1/4 (q_{i+1} - q_{i-1})^2          + 13/12 (q_{i-1} - 2 q_i + q_{i+1})^2
beta2 = 1/4 (q_{i-2} - 4 q_{i-1} + 3 q_i)^2 + 13/12 (q_{i-2} - 2 q_{i-1} + q_i)^2
wt_k = d_k/(epsW+beta_k)^2;  w_k = wt_k/sum
qrec = w0 v0 + w1 v1 + w2 v2
```
   Wrappers `reconst_qLx`,`reconst_qRx` just call these.
   IMPORTANT naming: in flux_x, `reconst_qLx(...,i)` gives the L state and `reconst_qRx(...,i)` the R state at interface i+1/2.

### Centered WENO derivative wDx / Dx (solver.c:139-181)
`wDx` = 4th-order-in-smooth-regions centered WENO finite difference of a primitive array (weights C1=1/6,C2=2/3,C3=1/6, three 2nd-order one-sided/centered stencils blended by the same beta). Used to compute `d_x xi`, `d_x ux` (the xiP,uxP arrays).
`Dx(arr,n,i)` = wDx for interior `1<i<N-2`, else fall back to plain centered `(arr[i+1]-arr[i-1])/(2 dx)`.

### Kurganov-Tadmor (KT) central flux (solver.c:499-540, flux_x)
Conserved flux at interface i+1/2. Reconstruct all of xi,ux,xiD,uxD,xiP,uxP to L and R states, build the physical flux from the FULL viscous stress tensor:
- COMP==0 (Ttt equation): physical flux = `T_tx`, dissipation jump = `T_tt(R)-T_tt(L)`.
- COMP==1 (Ttx equation): physical flux = `T_xx`, dissipation jump = `T_tx(R)-T_tx(L)`.
```
a = 1.0;                                  // KT max wave speed = speed of light (conformal)
flux_iph = 0.5*( iphL + iphR - a*jump_iph )   // Rusanov/LLF-style KT central flux
```
The local max speed is hardwired to c=1 (justified: conformal BDNK characteristic speeds <= c). Then divergence:
```
Ttx_cx(n,i) = (flux_x(n,i,0) - flux_x(n,i-1,0))/dx     // d_x of T_tx  -> RHS of Ttt
Txx_cx(n,i) = (flux_x(n,i,1) - flux_x(n,i-1,1))/dx     // d_x of T_xx  -> RHS of Ttx
```

---

## 5. Time integration: Heun (RK2 / explicit trapezoidal) (solver.c:562-643)

`Heun_solve_system(n)` with 3 time levels TL=3 (n, n+1, n+2). Predictor (Euler) then corrector (trapezoid).

**Predictor (n -> n+1):**
```
Ttt[n+1] = Ttt[n] + dt*( -Ttx_cx(n,i) )
Ttx[n+1] = Ttx[n] + dt*( -Txx_cx(n,i) )
xi[n+1]  = xi[n]  + dt*xiD[n][i]      // (BDNK_PRIM branch)
ux[n+1]  = ux[n]  + dt*uxD[n][i]
```
Then recompute xiP,uxP (spatial derivs) and xiD,uxD (prim time-derivs via compute_xiD/uxD) at n+1.

**Corrector (n -> n+2):**
```
Ttt[n+2] = Ttt[n] + dt/2*( -Ttx_cx(n,i) - Ttx_cx(n+1,i) )
Ttx[n+2] = Ttx[n] + dt/2*( -Txx_cx(n,i) - Txx_cx(n+1,i) )
xi[n+2]  = xi[n]  + dt/2*( xiD[n][i] + xiD[n+1][i] )
ux[n+2]  = ux[n]  + dt/2*( uxD[n][i] + uxD[n+1][i] )
```
Then recompute xiP,uxP and xiD,uxD at n+2. Ghost cells set after every stage. `copy_forward` copies level TL-1 (=n+2) back to level 0 for the next step (main loop calls `Heun_solve_system(0)` each iteration).

**PRIM macro CRITICAL BUILD NOTE:** solver.c uses `#if PRIM == BDNK_PRIM ... #elif PRIM == PF_PRIM`, but **none of PRIM, BDNK_PRIM, PF_PRIM are #defined anywhere** in the repo. In C preprocessor, undefined identifiers in `#if` are replaced by 0, so `PRIM == BDNK_PRIM` becomes `0 == 0` = TRUE. => **The BDNK_PRIM branch is compiled by default** (xi,ux integrated from xiD,uxD). To use the pure perfect-fluid evolution instead you would `#define PRIM PF_PRIM` and `#define BDNK_PRIM 1` style — but as shipped it builds BDNK. (Verified via `gcc -E`.)

---

## 6. Grid, BCs, ID, parameters (parameters.h)

- `N = RES_MULTIPLE*BASE_RESOLUTION + 1` (default 128+1=129 cells); `MAX_TIMESTEP = RES_MULTIPLE*BASE_NUM_TIMESTEP+1` (1025). The +1 keeps cell centers/time levels aligned under resolution doubling for convergence tests.
- `dx = (X_MAX-X_MIN)/(N-1)`; domain X_MIN=-200, X_MAX=200 -> dx = 400/128 = 3.125.
- `dt = CFL*dx`, CFL=0.1 -> dt = 0.3125.
- Ghost cells: 3 on each side (loops run i=3..N-4). BC=GHOST (outflow: copy index 3 into 0,1,2 and index N-4 into N-3,N-2,N-1) or PERIODIC.
- Initial data (default GAUSSIAN): `eps = A*exp(-(x-mean)^2/spread) + const` with A=1, mean=0, spread=25, const=0.1 (const+coeff prevent v->c NaNs); ux=0. `xi=ln(eps)`. Conserved set via T_tt,T_tx with all dissipative args 0 (=PF). Also STEP and SMOOTH_SHOCK (uses PF jump conditions: epsR=(epsL-9 vL^2 epsL)/(3(vL^2-1)), vR=1/(3 vL), tanh/erf profile).
- Output: plaintext, one file per var in `datafiles/` (must pre-exist or segfault). Rows=time, cols=x. Vars: xi,ux,Ttt,Ttx,xiD,uxD,VISC,eps. Save every TS_STEP=10 steps.

---

## 7. Open questions / gaps
- The explicit causality/stability inequalities for (eta0,lambda0,chi0) are NOT in the code; only the canonical ratios (25/7, 25/4) are hardcoded. Need 2201.12317 eqs 19-20 + causality section for the inequality forms used for neutron-star (non-conformal) generalization.
- `PRIM`/`BDNK_PRIM`/`PF_PRIM` are undefined macros -> rely on C's "undefined=0" so BDNK branch compiles. A Julia port should make this an explicit enum/flag, not implicit.
- The big return expressions in compute_xiD/uxD/compute_*PF are Mathematica-generated; their derivation (the 2x2 Jacobian inversion of the constitutive relations) is not shown in code. For a non-conformal NS EOS these must be re-derived.
- `eps[]`, `xiP`/`uxP` global arrays and `T_*0` PF functions are partly vestigial; `flux_xx`,`flux_tx`,`Ttt_pv`,`Ttx_pv`,`res` arrays declared but largely unused.
