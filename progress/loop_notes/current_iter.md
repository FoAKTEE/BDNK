# Iter 009 — dispatch 10-agent multi-team reproduction workflow

**Paper anchors.** Bussières 2604.13208 (R4 axial QNM); Shum 2509.15303 (R5
nonlinear Cowling); Kovtun/CR/PMP/Pandya (analytic fleet).

**Dispatched (background workflow wu7k5dafo, 10 agents, claude 4.8 ultracode).**
- Fleet (4 parallel, analytic): kovtun_sound, is_contrast (CR IS sound speed),
  pmp_causality (acaus_instab classification), pandya_char (char speeds + cons).
- **R4 axial-QNM team** (3-stage pipeline): axial wave-eqs (Redondo-Yuste 14a/14b)
  → shooting+Leaver QNM solver → adversarial verify vs Bussières (10.4884 kHz,
  29.587 µs).
- **R5 nonlinear-Cowling team** (3-stage pipeline): general-EOS BDNK recovery +
  isotropic coords → SSP-RK3 stellar evolution → FFT QNM extraction + verify vs
  Shum (F=2.69/H1=4.55/H2=6.36 kHz, decay 0.00157 M☉⁻¹).
- Agents write self-contained repro/*.jl, run Julia, return honest matched/gap;
  NO fabrication; NO figure gen (orchestrator does figures on matched results).

**On completion (orchestrator).** Integrate verified modules into the package;
generate + VLM the figures for matched reproductions; append knowledge/error
rows under the DAG nodes (s1b.*, s1c.*, s2.is_contrast); commit per substage.

**Loop status.** active:true (ON). gate=continue. The fleet is the iteration's
work; the loop persists.

**Verifier output (pre-dispatch).** julia test/runtests.jl → all pass (8 solid
nodes); Kovtun×3 + Bjorken figures VLM-matched.
