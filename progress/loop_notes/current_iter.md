# Iter 001 — bootstrap + STEP 0 trunk

**Paper anchor.** Project bootstrap; STEP 0 (shared EOS + primitive recovery).

**Shipped.**
- Understanding fan-out (14 sources) → reformulated notes + SYNTHESIS.md.
- 16-node project DAG (Mermaid), doubly linked to knowledge + error ledgers.
- Julia pkg `BDNKStar` STEP 0: Units, Numerics(Brent), EOS (PolytropeEnergy,
  ShumPolytrope, IdealGas, TabulatedBarotrope), Recovery (barotropic/general
  con2prim + BDNK gradient-frozen), Transport, Causality.
- STEP 0 validation figure (3 panels, log scale) + VLM-inspected.

**Verifier output.**
- `julia test/runtests.jl` → 189/189 pass (after Shum add).
- round-trip max rel err 4.3e-15 (gate 1e-10); tabulated 1e-5→8e-12 (N=25→800).
- VLM: panel A all EOS ≪ gate; B clean N⁻⁴; C char speeds + τ_P<1 acausal mark.

**Next-3 roadmap.**
1. step0.bdnk_recovery → solid (port conformal compute_xiD/uxD, CrossCheck).
2. s1a.tov_background (TOV RK4 + implicit-CN cross-check).
3. s1a.radial_eig + heat_criterion (CY c_s²−c_n²≥0).

**Research-state delta.** STEP 0 gate accepted; s1a.tov_background now blocking.
