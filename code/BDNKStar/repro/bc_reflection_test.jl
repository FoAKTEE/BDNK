#=
    bc_reflection_test.jl — BOUNDARY-CONDITION / REFLECTION sensitivity of the
    extracted f/p mode frequencies. Data only.

    TWO independent probes:

    (A) SphBDNK time-domain ℓ=2 f-mode. Vary
          (1) the SURFACE ghost treatment in _fill_ghosts! :
                :outflow  F[Nr+2]=F[Nr+1]      (production zero-gradient)
                :dirichlet F[Nr+2]=-F[Nr+1]     (reflecting / fixed pressure node)
                :extrap   F[Nr+2]=2F[Nr+1]-F[Nr] (linear extrapolation)
                :copy2    F[Nr+2]=F[Nr]          (wider zero-gradient stencil)
          (2) the global Kreiss–Oliger dissipation σ_ko (the only "damping layer"
              present): a reflection-contaminated mode shifts with σ_ko / BC.
        The f-mode frequency is read from the periodogram peak of the ℓ=2 quadrupole.

    (B) NonRadialModes eigensolver. Vary the outer matching radius rf=R(1-δ_out)
        and the centre start r0=R·ε0. A reflection/BC-contaminated eigenvalue moves
        with the surface match; a true eigenmode is insensitive.

    RUN:  cd code/BDNKStar && export PATH="$HOME/.local/bin:$PATH" \
          && JULIA_NUM_THREADS=3 julia --project=. repro/bc_reflection_test.jl
=#
using Printf
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
const S    = BDNKStar.SphBDNK
const NRM  = BDNKStar.NonRadialModes
const MSUN = BDNKStar.Units.Msun_to_km

# ---- shared M=1.4 n=1 K=100 polytrope background ---------------------------
eos = ShumPolytrope(100.0)
εc  = 0.00128 + 100*0.00128^2

# =====================================================================
#  (A) SphBDNK time-domain f-mode vs SURFACE BC and KO damping
# =====================================================================
# The production surface ghost is F[Nr+2,j] = F[Nr+1,j] (zero-gradient outflow).
# We swap ONLY that surface line; centre parity & pole parity are kept identical.
function patch_surface!(mode::Symbol)
    surf = mode === :outflow   ? :(F[Nr+2,j] = F[Nr+1,j]) :
           mode === :dirichlet ? :(F[Nr+2,j] = -F[Nr+1,j]) :
           mode === :extrap    ? :(F[Nr+2,j] = 2F[Nr+1,j] - F[Nr,j]) :
           mode === :copy2     ? :(F[Nr+2,j] = F[Nr,j]) :
           error("bad surface mode")
    @eval S begin
        function _fill_ghosts!(st::SphBDNKState, Nr, Nθ)
            @inbounds for jj in 1:Nθ
                j = jj+1
                for (F,par) in ((st.δρ,1),(st.δSr,-1),(st.δSθ,1),(st.δvr,-1),(st.δvθ,1))
                    F[1,j] = par*F[2,j]            # centre parity (unchanged)
                    $surf                          # SURFACE — varied
                end
            end
            @inbounds for ii in 1:Nr
                i = ii+1
                for (F,par) in ((st.δρ,1),(st.δSr,1),(st.δSθ,-1),(st.δvr,1),(st.δvθ,-1))
                    F[i,1] = par*F[i,2]; F[i,Nθ+2] = par*F[i,Nθ+1]   # poles (unchanged)
                end
            end
        end
    end
end

# extract the ℓ=2 quadrupole-moment peak frequency [kHz] over a physical band
function fmode_kHz(s, e; dt, nsteps, sample, νlo, νhi, nν=4000)
    st = S.SphBDNKState(s.grid.Nr, s.grid.Nθ)
    S.seed_sphbdnk_l2!(st, e; A=1e-3)
    ts, q2, en = S.evolve_sphbdnk!(st, e; dt=dt, nsteps=nsteps, sample=sample)
    all(isfinite, q2) || return (NaN, false, NaN, NaN)
    νs = range(νlo, νhi; length=nν)
    P  = S.periodogram(ts, q2, νs)
    νpk = νs[argmax(P)]
    fkHz = S.freq_kHz_cyclic(νpk)
    (fkHz, true, en[end]/en[1], maximum(abs, q2))
end

# build star + run-control common to (A)
Nr, Nθ = 200, 24
s  = build_sphstar(eos, εc; Nr=Nr, Nθ=Nθ)
dt = 0.2 * s.grid.dr
nsteps, sample = 12000, 4
# f-mode band, in CYCLIC geometric freq ν = f[kHz]·Lunit·kHz_to_km. The ℓ=2 f-mode
# of this M=1.4 Cowling star sits near 1.9 kHz; use a band that brackets it but
# excludes the centre/grid content below ~1 kHz and the p-mode comb above ~2.6 kHz.
kHzgeom(f) = f * MSUN * BDNKStar.Units.kHz_to_km
νlo = kHzgeom(1.3)
νhi = kHzgeom(2.5)
patch_surface!(:outflow)
e0 = setup_sphbdnk(s; η̂=0.0, σ_ko=0.02)
f_probe, ok, _, _ = fmode_kHz(s, e0; dt=dt, nsteps=nsteps, sample=sample, νlo=νlo, νhi=νhi)
@printf("[probe] production outflow, σ_ko=0.02  f-mode peak=%.4f kHz  finite=%s\n", f_probe, ok)

println("\n=== (A1) SURFACE BC sweep  (η̂=0, σ_ko=0.02, Nr=$Nr Nθ=$Nθ) ===")
println("surface_BC        f_kHz      finite   E_end/E0   |q2|max")
fA1 = Float64[]
for m in (:outflow, :dirichlet, :extrap, :copy2)
    patch_surface!(m)
    e = setup_sphbdnk(s; η̂=0.0, σ_ko=0.02)
    f, okf, eratio, qmax = fmode_kHz(s, e; dt=dt, nsteps=nsteps, sample=sample, νlo=νlo, νhi=νhi)
    okf && isfinite(f) && push!(fA1, f)
    @printf("%-14s  %8.4f   %-6s   %8.3f   %.3e\n", m, f, okf, eratio, qmax)
end
if length(fA1) ≥ 2
    @printf(">> (A1) f-mode spread across SURFACE BCs: min=%.4f max=%.4f  Δ=%.4f kHz (%.2f%%)\n",
            minimum(fA1), maximum(fA1), maximum(fA1)-minimum(fA1),
            100*(maximum(fA1)-minimum(fA1))/(sum(fA1)/length(fA1)))
end

println("\n=== (A2) KO damping-layer σ_ko sweep  (η̂=0, surface=outflow) ===")
patch_surface!(:outflow)
println(" σ_ko       f_kHz      finite   E_end/E0")
fA2 = Float64[]
for σ in (0.0, 0.01, 0.02, 0.04, 0.08)
    e = setup_sphbdnk(s; η̂=0.0, σ_ko=σ)
    f, okf, eratio, _ = fmode_kHz(s, e; dt=dt, nsteps=nsteps, sample=sample, νlo=νlo, νhi=νhi)
    okf && isfinite(f) && push!(fA2, f)
    @printf("%6.3f    %8.4f   %-6s   %8.3f\n", σ, f, okf, eratio)
end
if length(fA2) ≥ 2
    @printf(">> (A2) f-mode spread across σ_ko: min=%.4f max=%.4f  Δ=%.4f kHz (%.2f%%)\n",
            minimum(fA2), maximum(fA2), maximum(fA2)-minimum(fA2),
            100*(maximum(fA2)-minimum(fA2))/(sum(fA2)/length(fA2)))
end

# restore production surface for any later use
patch_surface!(:outflow)

# =====================================================================
#  (B) NonRadialModes eigensolver vs outer match radius + centre start
# =====================================================================
println("\n=== (B) NonRadialModes f/p1: outer match rf=R(1-δ_out), centre r0=R·ε0 ===")
# nonradial_cowling_spectrum hard-codes r0=R*1e-4, rf=R*(1-1e-3). Re-implement the
# driver with those exposed, reusing the module's internal shooter & background.
function spectrum_bc(eos, εc; l=2, nmodes=2, N=4000, h_tov=2e-4,
                     δ_out=1e-3, ε0=1e-4, ω2lo=5e-4, ω2hi=0.12, nscan=1200)
    star = solve_tov(eos, εc; h=h_tov)
    R = star.R
    rt, mt, νt, et = star.r, star.m, star.ν, star.ε
    εf = εc*1e-9
    bg(r) = (NRM._lin(rt,mt,r), NRM._lin(rt,νt,r), max(NRM._lin(rt,et,r), εf))
    r0 = R*ε0
    rf = R*(1 - δ_out)
    D(ω2) = NRM._shoot(eos, bg, l, ω2, r0, rf, N)[1]
    grid = range(ω2lo, ω2hi; length=nscan)
    ω2s = Float64[]; Dprev = D(first(grid)); ω2prev = first(grid)
    for ω2 in Iterators.drop(grid, 1)
        Dc = D(ω2)
        if isfinite(Dc) && isfinite(Dprev) && Dc*Dprev < 0
            res = BDNKStar.Numerics.brent(D, ω2prev, ω2; xtol=1e-12)
            res.converged && push!(ω2s, res.root)
        end
        Dprev = Dc; ω2prev = ω2
    end
    sort!(ω2s)
    k = min(nmodes, length(ω2s))
    [NRM.freq_kHz_from_omega2(w) for w in ω2s[1:k]], R
end

println("\n-- (B1) outer match radius δ_out  (centre ε0=1e-4 fixed) --")
println(" δ_out       f_kHz     p1_kHz")
fB_f = Float64[]; fB_p = Float64[]
for δ in (3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
    fr, _ = spectrum_bc(eos, εc; l=2, nmodes=2, N=4000, δ_out=δ, ε0=1e-4)
    f  = length(fr)≥1 ? fr[1] : NaN
    p1 = length(fr)≥2 ? fr[2] : NaN
    isfinite(f)  && push!(fB_f, f)
    isfinite(p1) && push!(fB_p, p1)
    @printf("%7.1e   %7.4f   %7.4f\n", δ, f, p1)
end
if length(fB_f) ≥ 2
    @printf(">> (B1) f-mode spread vs δ_out:  Δ=%.4f kHz (%.3f%%)\n",
            maximum(fB_f)-minimum(fB_f), 100*(maximum(fB_f)-minimum(fB_f))/(sum(fB_f)/length(fB_f)))
end
if length(fB_p) ≥ 2
    @printf(">> (B1) p1-mode spread vs δ_out: Δ=%.4f kHz (%.3f%%)\n",
            maximum(fB_p)-minimum(fB_p), 100*(maximum(fB_p)-minimum(fB_p))/(sum(fB_p)/length(fB_p)))
end

println("\n-- (B2) centre start ε0  (outer δ_out=1e-3 fixed) --")
println("  ε0        f_kHz     p1_kHz")
fB2_f = Float64[]; fB2_p = Float64[]
for e0 in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5)
    fr, _ = spectrum_bc(eos, εc; l=2, nmodes=2, N=4000, δ_out=1e-3, ε0=e0)
    f  = length(fr)≥1 ? fr[1] : NaN
    p1 = length(fr)≥2 ? fr[2] : NaN
    isfinite(f)  && push!(fB2_f, f)
    isfinite(p1) && push!(fB2_p, p1)
    @printf("%7.1e   %7.4f   %7.4f\n", e0, f, p1)
end
if length(fB2_f) ≥ 2
    @printf(">> (B2) f-mode spread vs centre ε0:  Δ=%.4f kHz (%.3f%%)\n",
            maximum(fB2_f)-minimum(fB2_f), 100*(maximum(fB2_f)-minimum(fB2_f))/(sum(fB2_f)/length(fB2_f)))
end
if length(fB2_p) ≥ 2
    @printf(">> (B2) p1-mode spread vs centre ε0: Δ=%.4f kHz (%.3f%%)\n",
            maximum(fB2_p)-minimum(fB2_p), 100*(maximum(fB2_p)-minimum(fB2_p))/(sum(fB2_p)/length(fB2_p)))
end

println("\nDONE.")
