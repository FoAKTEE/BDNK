# Equation DAG — paper_bdnk-hmns

16 nodes. Legend: ● solid · ◐ preliminary · ○ hypothesis · ✗ blocking · □ future · △ concept-advance. Node badge `k<N>` = knowledge records under the node, `t<N>✗<F>` = trials (F failed). Dashed edge = predecessor outside scope.

```mermaid
flowchart TD
  n_step0_eos["● <b>step0.eos</b> · k2<br/>Add Shum Γ=2 cold polytrope p(ε)=[1+2κε−√(1+4κε)]/2κ (1C …"]:::solid
  n_step0_con2prim_ideal["● <b>step0.con2prim_ideal</b> · k2 t2<br/>STEP 0 validation figure (round-trip vs |v|, tabulated N⁻…"]:::solid
  n_step0_causality["● <b>step0.causality</b> · k1<br/>Causality/stability monitor: full BDN Thm-I (a,b,d,e), re… <br/><i>BDNK characteristic speeds, Causality bounds, Simplified constraints (eq.44), causality (A1)+(a) (1171,1179)…</i>"]:::solid
  n_step0_bdnk_recovery["● <b>step0.bdnk_recovery</b> · k2 t1<br/>Conformal BDNK constitutive + compute_xiD/uxD recovery po…"]:::solid
  n_s1a_tov_background["● <b>s1a.tov_background</b> · k3 t2<br/>Reproduce Shum M_T=1.40016 M☉ (target 1.4; M☉=G=c=1 units…"]:::solid
  n_s1a_radial_eig["◐ <b>s1a.radial_eig</b> · k2 t1<br/>Radial Cowling eigensolver (NSO A0/A1/A2 operator); spect…"]:::preliminary
  n_s1a_heat_criterion["□ <b>s1a.heat_criterion</b> · k1<br/>Heat-conduction stability check: build 𝔽 (positive-defini… <br/><i>Viscous operator F (Eq.64-65), Heat source Y / G, First-order freq shift (Eq.62), Sound speeds (Eq.26,27)</i>"]:::future
  n_s1b_axial_wave_eqs["□ <b>s1b.axial_wave_eqs</b> · k1<br/>Axial (odd-parity) linear sector: assemble the two couple… <br/><i>axial O_l wave eq (eq.14a), axial O_n wave eq (eq.14b), Eqs.17-18 (axial QNM system), Eq.19 (Regge-Wheeler potential)…</i>"]:::future
  n_s1b_qnm_freqdomain["□ <b>s1b.qnm_freqdomain</b> · k1<br/>Axial QNM eigenvalue solver: interior two-seed shooting +… <br/><i>Eqs.17-18 (axial QNM system), Eq.24 (surface regularity BC), Eq.41 (compactness scaling)</i>"]:::future
  n_s1c_hrsc_core["◐ <b>s1c.hrsc_core</b> · k2 t1<br/>Flat-space conformal BDNK evolution engine (WENO5+KT+Heun…"]:::preliminary
  n_s1c_imex["□ <b>s1c.imex</b> · k1<br/>Time integration: SSP-RK3 (Shum) / Heun RK2 (conformal) w… <br/><i>Heun RK2, Kreiss-Oliger dissipation</i>"]:::future
  n_s1c_qnm_extract["□ <b>s1c.qnm_extract</b> · k1<br/>QNM extraction from the nonlinear ringdown: fit ε̃_c(t)=A… <br/><i>QNM fit, Table 2</i>"]:::future
  n_s2_gr_coupling["□ <b>s2.gr_coupling</b> · k1<br/>Fork A: couple BDNK hydro to 1+1D dynamical GR (harmonic … <br/><i>Bjorken EOM (eq.50, Milne), Rankine-Hugoniot (steady shock), Steady-shock erf initial data</i>"]:::future
  n_s2_is_contrast["□ <b>s2.is_contrast</b> · k1<br/>IS (Müller-Israel-Stewart) baseline + BDNK-vs-IS migratio… <br/><i>Eq.30 c_s'^2, Eq.31 causality fix, Eq.9 cons vars</i>"]:::future
  n_s3_cowling_3p1["□ <b>s3.cowling_3p1</b> · k1<br/>Fork B: 3+1D Cowling non-radial BDNK modes; cubed-sphere … <br/><i>T_munu (eq.16), constitutive (eqs.17-19), Velocity ansatz, Vorticity diagnostic</i>"]:::future
  n_s4_production["□ <b>s4.production</b> · k1<br/>Production runs: realistic tabulated EOS, physical transp… <br/><i>Eq.58 EOS, Frame choice (eq.41-43)</i>"]:::future
  n_step0_eos --> n_step0_causality
  n_step0_con2prim_ideal --> n_step0_bdnk_recovery
  n_step0_causality --> n_step0_bdnk_recovery
  n_step0_eos --> n_s1a_tov_background
  n_s1a_tov_background --> n_s1a_radial_eig
  n_step0_causality --> n_s1a_radial_eig
  n_s1a_radial_eig --> n_s1a_heat_criterion
  n_s1a_tov_background --> n_s1b_axial_wave_eqs
  n_step0_causality --> n_s1b_axial_wave_eqs
  n_s1b_axial_wave_eqs --> n_s1b_qnm_freqdomain
  n_step0_bdnk_recovery --> n_s1c_hrsc_core
  n_s1a_tov_background --> n_s1c_hrsc_core
  n_s1c_hrsc_core --> n_s1c_imex
  n_s1c_imex --> n_s1c_qnm_extract
  n_s1c_qnm_extract --> n_s2_gr_coupling
  n_s2_gr_coupling --> n_s2_is_contrast
  n_s1c_qnm_extract --> n_s3_cowling_3p1
  n_s2_is_contrast --> n_s4_production
  n_s3_cowling_3p1 --> n_s4_production
  classDef solid fill:#e6ffed,stroke:#28a745,color:#000;
  classDef preliminary fill:#fff8e1,stroke:#d4a017,color:#000;
  classDef hypothesis fill:#e7f0ff,stroke:#4977c7,color:#000;
  classDef blocking fill:#ffe3e3,stroke:#d33,color:#000;
  classDef future fill:#f2f2f2,stroke:#999,color:#555;
  classDef amended fill:#f0f0f0,stroke:#aaa,color:#888;
```
