# Iter 012 — dispatch round-2 10-agent workflow (R5 production + PMP viscous)

**Dispatched (background workflow wscn1tor1, 10 agents, claude 4.8 ultracode).**
- **R5 production team** (5 agents): perf-optimize shum_evolve -> resolution
  ladder (Dr=0.04/0.02/0.01, t_f->8000, background runs) -> FFT + damped-sinusoid
  fit + Richardson extrapolation -> match Shum F=2.69/H1=4.55/H2=6.36 kHz +
  decay 0.0011 /M☉ (continuum). The overtones+decay R1 left open.
- **PMP viscous-1D team** (4 agents): general-EOS ideal-gas viscous BDNK 1D core
  -> telegrapher / shock-instability+acaus / heat-stationary figures.
- **R4 ultracompact** (1 agent): plot_ultracompact + eta-mode family + avoided crossing.
- Agents write self-contained repro/*.jl, run Julia (background for long runs),
  honest matched/gap, NO fabrication, NO figure gen.

**On completion (orchestrator).** Independently re-verify; integrate; generate +
VLM figures (R5 QNM_plot, PMP plots); append ledgers; commit per substage.

**Prior round (wu7k5dafo) results — committed:** R4 axial QNM <0.04% (SOLID,
re-verified); fleet 4/4 (kovtun_sound, is_contrast, pmp_causality, pandya_char);
R5 fundamental f_nl=2.70 kHz matched (PRELIM). 11 solid + 4 prelim nodes.

**Loop:** active:true (ON), gate=continue.
