# Source: nb-BDNK_frame_constraints (STEP 0, critical)

**File:** `/data/haiyangw/claude/BDNK/ref-code/BDNK_frame_constraints/BDNK_frame_constraints.nb`
(Mathematica 13 notebook, plain-text .nb, 2541 lines). Author: Alex Pandya, Aug 19 2022.

**Purpose:** Symbolically verify that the gamma-law ("Gamma-law") BDNK hydrodynamic frame
of Pandya-Most-Pretorius (PMP, arXiv:2209.09265, PRD 106, 123036) satisfies both the
**causality** and **linear-stability** constraints of Bemfica-Disconzi-Noronha (BDN,
arXiv:2009.11388, PRX 12, 021044), eqs (47)/(48) of BDN. These define the STEP-0
causality monitor.

References embedded in notebook:
- [1] F. S. Bemfica, M. M. Disconzi, J. Noronha, PRX 12, 021044 (2022), 2009.11388.
- [2] A. Pandya, E. R. Most, F. Pretorius, PRD 106, 123036 (2022), 2209.09265.

---

## 1. Shorthand symbols A,B,C,D,E (BDN eq 47-48 building blocks)

Defined in the cell `In[40]` (notebook lines ~630-688). Background-state quantities:
`rho` = energy density (rest-frame), `csSq` = c_s^2 = (dp/de)_s sound speed squared,
`tauE` = tau_epsilon, `tauQ`, `tauP` = relaxation times, `eta` = shear visc,
`V` (bulk-related coeff), `sigma` (conductivity-like coeff), `kappaE`, `kappaS`
(heat-flux coefficients), `alpha` free param, `omega = kappaS/kappaE`.

```
betaE = tauQ*alpha*csSq + (sigma/rho)*kappaE          (* beta_epsilon *)
A = rho*tauE*tauQ
B = -tauE*(rho*csSq*tauQ + V + sigma*kappaS) - rho*tauP*tauQ
C = tauP*(rho*csSq*tauQ + sigma*kappaS) - betaE*V
D = rho*csSq*(tauE + tauQ) + V + sigma*kappaE
E = sigma*(alpha*csSq*kappaS - csSq*kappaE)
```

These are the coefficients of the dispersion relation polynomial; A,B,C are also used
for the characteristic-speed quadratic (see section 6).

---

## 2. Linear stability constraints — BDN eq (48), form (...) > 0

Cell `In[12]` (lines ~184-310). (48a) is split into stabA1, stabA2.

```
stabA1 = (tauE+tauQ)*Abs[B] - tauE*tauQ*D
stabA2 = tauE*tauQ*D - rho*csSq*tauE*tauQ*(tauE+tauQ)
stabB  = (tauE+tauQ)*Abs[B]*D + rho*tauE*tauQ*(tauE+tauQ)*E
         - ( tauE*tauQ*D^2 + rho*(tauE+tauQ)^2*C )
stabC  = csSq*D - E - rho*csSq^2*(tauE+tauQ)
stabD  = (tauE+tauQ)*( Abs[B]*(csSq*D - 2 E) + 2*csSq*rho*tauE*tauQ*E + C*D )
         - ( 2*csSq*rho*(tauE+tauQ)^2*C + tauE*tauQ*D*(csSq*D - E) )
stabE  = Abs[B]*D*( C*(tauE+tauQ) + E*tauE*tauQ )
         + 2*rho*tauE*tauQ*(tauE+tauQ)*C*E
         - ( rho*C^2*(tauE+tauQ)^2 + tauE*tauQ*(C*D^2 + rho*tauE*tauQ*E^2)
             + B^2*E*(tauE+tauQ) )
```
All required `> 0` (stabA1, stabA2, stabC may be `>= 0`).

### 2a. Reduced/rescaled stability constraints (from [2]), cell `In[18]` (lines ~321-416)
Rescaled hat quantities:
```
Bh = B/(rho*csSq*tauE*tauQ)
Ch = C/(rho*csSq^2*tauE*tauQ)
Dh = D/(rho*csSq*(tauE+tauQ))
Eh = E/(rho*csSq^2*(tauE+tauQ))
```
Reduced constraints:
```
rStabA1 = Abs[Bh] - Dh
rStabA2 = Dh - 1
rStabB  = Abs[Bh]*Dh + Eh - Dh^2 - Ch
rStabC  = Dh - Eh - 1
rStabD  = Abs[Bh]*Dh + Eh - Dh^2 - Ch
          - ( 2*Abs[Bh]*Eh + Ch - Dh*Eh - Eh - Ch*Dh )
rStabE  = ( Abs[Bh]*Dh + Eh - Dh^2 - Ch )*Ch
          - ( Eh + Bh^2 - Ch - Abs[Bh]*Dh )*Eh
```
The notebook proves (cell In[28]/In[34]) that `rStab* = quot * stab*` with `quot > 0`
(assuming csSq>0, rho>0, tauE>0, tauQ>0), so reduced ⇒ unreduced. Equivalence holds.

---

## 3. Causality constraints — BDN, form (...) > 0

Cell `In[67]` (lines ~1175-1273). (CAUS B) and (CAUS C) of [1] are each split in two.
```
causA  = rho*tauQ - eta
causB1 = ( tauE*(rho*csSq*tauQ + V + sigma*kappaS) + rho*tauP*tauQ )^2
         - 4*rho*tauE*tauQ*( tauP*(rho*csSq*tauQ + sigma*kappaS) - betaE*V )
causB2 = 4*rho*tauE*tauQ*( tauP*(rho*csSq*tauQ + sigma*kappaS) - betaE*V )
causC1 = 2*rho*tauE*tauQ
         - ( tauE*(rho*csSq*tauQ + V + sigma*kappaS) + rho*tauP*tauQ )
causC2 = tauE*(rho*csSq*tauQ + V + sigma*kappaS) + rho*tauP*tauQ
causD  = rho*tauE*tauQ + sigma*kappaS*tauP
         - ( tauE*(rho*csSq*tauQ + V + sigma*kappaS)
             + rho*tauP*tauQ*(1 - csSq) + betaE*V )
```
All `>= 0` (causA `> 0` strictly when checked at end). Mathematica directly proves
causB1, causB2, causC2 >= 0. causA, causC1, causD need the extra causality assumption
on tau-hat (section 5). causD = 0 means the fastest characteristic speed c+ = 1 (speed
of light): allowed (causal) but excluded by BDN's open-set treatment; replace `>=` with
`>` in the last causAssump entry to forbid it.

---

## 4. Gamma-law (ideal gas) BDNK frame ANSATZ from PMP [2]

Cell `In[50]` (lines ~765-805). This is THE explicit hydrodynamic frame for STEP 0.
Free dimensionless params: `t` = tau-hat, `Vh` = V-hat, `sigmah` = sigma-hat, `alpha`,
`omega`. `L` is a length/timescale, `kappaE` an overall heat-coeff scale.
```
tauE = t*L*Vh
tauQ = t*L*Vh                         (* tau_epsilon = tau_Q *)
tauP = 2*alpha*csSq*L*Vh
V     = Vh*L*rho*csSq
sigma = sigmah*Vh*L*rho*csSq/(-kappaE)
kappaS = omega*kappaE
```
Plus from In[40]: `betaE = tauQ*alpha*csSq + (sigma/rho)*kappaE`.

Shear/bulk in hat form (from causAssump / characteristic-speed cells):
```
eta = etah * rho * L * csSq           (* shear viscosity *)
Vh  = (4/3)*etah + zetahOveretah*etah  (* relation V-hat <-> shear+bulk *)
                                       (* i.e. V-hat = (4/3 + zetah/etah)*etah *)
```
So `Vh >= (4/3) etah` (bulk zeta-hat >= 0). etah = eta-hat, zetah/etah = zeta-hat/eta-hat.

### Frame ansatz simplifies the reduced shorthand (cell In[56], output):
```
Bh = -1 + (-1 - 2 alpha + sigmah*omega)/t
Ch = alpha/t + (sigmah - 2 alpha sigmah omega)/t^2
Dh = 1 + (1 - sigmah)/(2 t)
Eh = -( sigmah*(-1 + alpha*omega) )/(2 t)
```
Only alpha, omega, sigmah appear (after rescaling). NOTE: there is a sign nuance —
`sigma = sigmah*Vh*L*rho*csSq/(-kappaE)`; the notebook treats kappaE so that physical
sigma stays consistent. Treat sigmah, omega, alpha as the controllable knobs.

---

## 5. Parameter ranges / assumptions (the monitor's admissible box)

### assump (cell In[60], lines ~972-1009) — physical ranges + stability requirement:
```
csSq > 0, csSq < 1,
L > 0, rho > 0, t > 0,
omega > 0, omega < 1,
alpha >= 1,
Vh > 0,
alpha*omega > 0, alpha*omega < 1/2,
sigmah >= 0, sigmah <= 1/3        <-- KEY stability bound (PMP App. A)
```
`sigmah <= 1/3` is the nonlinear linear-stability requirement.

### causAssump (cell In[76], lines ~1488-1517) — causality:
```
Vh >= (4/3) etah,
Gm1 > 0, Gm1 < 1,          (* Gm1 = Gamma - 1, Gamma = adiabatic index *)
zetahOveretah >= 0,
etah > 0,
t >= ( Gm1*(2 - csSq) + csSq ) / ( 1 - csSq )   <-- KEY causality bound on tau-hat
```

### The two headline nonlinear frame inequalities (Summary cell, lines ~1799-1828):
- **Stability:**   `sigma-hat <= 1/3`
- **Causality:**   `tau-hat >= [ (Gamma - 1)(2 - c_s^2) + c_s^2 ] / ( 1 - c_s^2 )`

With `>=` the saturating case has c+ = 1 (light speed); use `>` to strictly exclude it.
`totalAssump = Join[assump, causAssump]` makes ALL six stability + all causality
constraints evaluate to True (cell In[102], outputs Out[102..114] all `True`).

---

## 6. Characteristic speeds (squared), cell In[115] (lines ~2216-2384)

```
c1Sq = eta/(rho*tauQ)  = csSq*etah/(t*Vh)          (* shear/transverse mode *)

cpSq, cmSq = (-B +/- Sqrt[B^2 - 4 A C])/(2 A)       (* sound/longitudinal modes *)
```
With the ansatz substituted:
```
c{p,m}Sq = (csSq/(2 t)) * ( 1 + t + 2 alpha - sigmah*omega
            +/- Sqrt[ 1 + t^2 + 4 alpha(1+alpha) - 2 sigmah(2+omega)
                      + t(2 - 2 sigmah omega)
                      + sigmah omega (4 alpha + sigmah omega) ] )
```
For causality all these squared speeds must be in [0,1] (real, subluminal). causD>=0
is exactly the condition c+ <= 1.

---

## STEP-0 monitor recipe (concrete)
Given background (rho, csSq, Gamma) and chosen frame knobs (t=tau-hat, Vh=V-hat,
sigmah, alpha, omega, etah, zetah/etah):
1. Check physical ranges (assump): csSq in (0,1), omega in (0,1), alpha>=1,
   alpha*omega < 1/2, sigmah in [0, 1/3], Vh > 0, Vh >= (4/3)etah, etah>0,
   zetah/etah >= 0, Gamma-1 in (0,1).
2. Stability monitor:  sigmah <= 1/3.
3. Causality monitor:  tau-hat >= ((Gamma-1)(2-csSq)+csSq)/(1-csSq).
4. (Optional full check) form A..E (sec 1) and evaluate causA, causB1/2, causC1/2,
   causD, stabA1..stabE > 0; or equivalently c1Sq, cpSq, cmSq in [0,1].

For an ideal gamma-law gas, c_s^2 is a function of Gamma and the thermodynamic state
(c_s^2 = Gamma p / (e + p) for the relativistic ideal gas / Gamma-law EOS used in PMP).

## OPEN QUESTIONS
- The notebook does NOT spell out the gamma-law EOS relation p(rho,...) or
  c_s^2(Gamma) explicitly; only Gm1=Gamma-1 in (0,1) and csSq in (0,1) appear.
  The closed-form EOS and c_s^2 expression must come from PMP 2209.09265 main text.
- Sign convention of `sigma = sigmah*Vh*L*rho*csSq/(-kappaE)` and the role/scale of
  kappaE / L are not pinned to physical units here; confirm against PMP Sec II/III.
- No numerical benchmark values (no specific tau-hat, sigma-hat numbers) are given;
  the notebook is purely symbolic (all checks return `True`).
