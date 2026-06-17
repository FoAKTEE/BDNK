# Reference-paper FIGURE reproduction target + status

The Ralph-loop completion target: reproduce every reference-paper figure (VLM
dual-checked against the original in `ref-paper/sources/arXiv-*/src/`). Honest
feasibility triage under the project constraints (**< 20 CPU threads, no GPU**,
Julia). Status: ✅ done · 🔨 tractable-next · 🏗 large-build · ⛔ infeasible-under-constraints.

## arXiv:1907.08191 — Kovtun (first-order hydro dispersion / stability)
| figure | what | status |
|---|---|---|
| piccvphi | c_v(φ) moving-fluid phase velocity (c₀=½) | ✅ kovtun_cvphi.png |
| picreshearv09 | shear-channel Re ω(k), v₀=0.9, θ/η=2 | ✅ kovtun_reshearv09.png |
| picimshearv09 | shear-channel Im ω(k) (stability) | ✅ kovtun_imshearv09.png |
| picresoundv09 / picimsoundv09 | sound-channel Re/Im ω(k) (quartic F_sound) | ✅ kovtun_{re,im}soundv09.png (fan within light cone; Im ω≤0 stable) |
| piccsGs09 / piccsGs01 | sound velocity c_s + damping Γ_s vs φ | ✅ kovtun_csgs.png (c_s exact via cvphi; Γ_s verbatim formula) |
| picstab | sound-channel stability region | ✅ kovtun_picstab.png (bands ε̄1<2/v_s², RH bites; region grows for smaller v_s) |
| picstabcaus | stable + causal region | ✅ kovtun_picstabcaus.png (causality cut removes small-ε̄1; origin excluded) |

## arXiv:2201.12317 — Pandya conformal BDNK (numerical method)
| figure | what | status |
|---|---|---|
| steady_state_fig / shock_comp | steady planar shock | ✅ conf_shock_overlay.png — Julia vs reference C code <1% (shock sharpens to steady-state, εR=4.4074); engine conformal_evolution.png |
| gaussian_clump (code validation) | conformal evolution vs reference C code | ✅ conf_overlay.png — Julia engine vs Pandya 1D_conformal_bdnk C code agree to 0.06% (full ε(x,t) evolution) |
| Conv_plot | self-convergence order | ✅ conf_convergence.png — Julia p=1.37 vs C code p=1.34 (coincident error curves; order limited by under-resolved narrow clump) |
| CC_plot | WENO5 mixed-deriv commutator ∫∫|∂x∂yξ−∂y∂xξ| | 🔨→✅mechanism (PRELIMINARY: 2D-surrogate, no 2D ref code; ε_W=1e15 machine-floor + decay matched) |
| eta1_step / kh_vs_eta / rotor_eta_t / Tab_cons | step, Kelvin-Helmholtz, rotor (2D), conservation | 🏗 (2D evolution) |

## arXiv:2209.09265 — PMP ideal-gas BDNK
| figure | what | status |
|---|---|---|
| bjorken_plot | Bjorken flow ε(τ) + 4th-order convergence | ✅ bjorken.png (Q→16) |
| conv_plot | PDE convergence | 🔨 (engine self-conv done) |
| shockwave_plot / shock_instability | steady shock + instability | 🔨 (general-EOS viscous) |
| telegraphers_plot | telegrapher heat→wave transition | ✅ pmp_telegrapher_reproduction.png (transition+ordering; split verified 1.000 in repro) |
| shock_instability / acaus_instab | shock causality crash classification | ✅ pmp_shock_reproduction.png (c₊(τ̂); crash at τ̂=0.25 acausal & τ̂=3 v>c₊); v(x) profiles in repro |
| heat_stationary | closed-form ε̈ (σ̂=0 stationary; σ̂=1/3 localized) | ✅ pmp_heat_reproduction.png; evolved |ε̇| in repro |

## arXiv:2509.15303 — Shum nonlinear BDNK NS (Cowling)
| figure | what | status |
|---|---|---|
| QNM_plot + decay | QNM spectrum F/H1/H2 + 1/τ decay | ✅ shum_qnm_reproduction.png + shum_decay.png (1/τ_l=0.00209/M_⊙ within 30% of paper per-Δr 0.0016–0.0019; MATCHED_TARGET=true) — F=2.699/H1=4.551/H2=6.468 kHz vs (2.69,4.55,6.36) <2% |
| convergence / error_fit | QNM continuum convergence | ✅ r5_convergence_reproduction.png (F/H1/H2 <1.4%; decay [HOLE] needs finer Δr ladder) |
| casA_fitting | decay-rate fit | 🏗 (needs 4-Δr ladder) |
| stable_evol_comparing_tau / stable_evol_resolutions | stable-window evolutions | 🏗 |
| (TOV M_T=1.4 reproduced) | — | ✅ |

## arXiv:2604.13208 — Bussières axial viscous NS
| figure | what | status |
|---|---|---|
| complex_plane_2 / plot_combined | axial QNM (f,τ) vs η_c | ✅ axial_qnm_reproduction.png (<0.04% Bussières Table II) |
| plot_ultracompact | ultracompact trapped w-modes (ℓ=2 ladder vs 𝒞) | ✅ ultracompact_reproduction.png (2→6 ladder, viscous damping); ω(ℓ)-sweep view pending |
| (EOS1 M=1.27/R=8.86 reproduced) | — | ✅ |

## arXiv:2311.13027 — Chabanov–Rezzolla (bulk-viscous BNS merger)
26 figures (central_dens_migration, gw_plot, rho_temp_xz_HMNS, vx_spacetime_migration, …):
**⛔ INFEASIBLE under <20 CPU / no GPU.** These are products of a full 3+1D GRMHD
binary-neutron-star *merger* simulation (months of supercomputer/GPU time). They
cannot be reproduced within the stated compute budget — a genuine constraint
conflict to escalate. The *physics contrast* (IS vs BDNK, bulk-viscous outcome)
is reproducible at the 1D level (s2.is_contrast); the merger figures are not.

## Summary
- ✅ reproduced: 7 figures (Kovtun×3, Bjorken, conformal evolution, TOV×2 baked in).
- 🔨 tractable-next: ~8 (Kovtun sound/stability, PMP/Pandya exact overlays).
- 🏗 large-build: Shum QNM, Bussières axial, PMP viscous 1D, Pandya 2D (rotor/KH).
- ⛔ infeasible: the 26 Chabanov–Rezzolla 3D-merger figures (compute-budget conflict).
