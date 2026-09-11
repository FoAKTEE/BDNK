#=
    BDNK viscoresistive relativistic MHD — CAUSALITY / FRONT-VELOCITY figure
    (Lier, Armas, Porth 2026, arXiv:2606.22691;  Figs.1–2 analogues).

    (A) SOUND front velocity² W²(τ_X) and the Eq.16 subluminal WINDOW
        [2D_u, τ_ε−D_u−3D_uτ_ε/(τ_u−D_ε)] in which the larger branch is real and
        W²≤1 (luminal boundary W²=1 at the right edge).
    (B) The ALFVÉN + MAGNETOSONIC front velocities |Re W(θ)| vs propagation angle
        θ on an Orszag–Tang background, for τ_X=0 (anti-diffusive: a complex
        front appears, |Im W|>0) vs τ_X=0.2 (real, bounded), with v_max marked
        — the decisive demonstration that τ_X is required (Eq.B13).
    (C) The boosted-telegrapher benchmark: numerical-vs-analytic profile +
        the N⁻² convergence of the L1 error (Eq.21, Fig.2).

    Run: julia --project=code/BDNKStar/viz code/BDNKStar/viz/bdnk_mhd_causality.jl
=#
using Pkg
Pkg.activate(@__DIR__)
include(joinpath(@__DIR__, "..", "src", "BDNKStar.jl"))
using .BDNKStar
using CairoMakie
using Printf
CairoMakie.activate!(type="png")

outdir = joinpath(@__DIR__, "..", "figures"); isdir(outdir) || mkpath(outdir)

fig = Figure(size=(1500, 470))

# ── Panel A: SOUND W²(τ_X) + Eq.16 subluminal window ─────────────────────────
axA = Axis(fig[1,1], xlabel="τ_X", ylabel="sound front  W²",
           title="(A) Sound front W²(τ_X) + subluminal window (Eq.15/16)")
Du, Dε, τu = 1e-2, 2e-3, 2e-1
lo, hi = sound_subluminal_window(Du, Dε, τu)
τXs = range(0.0, hi*1.25; length=400)
Wp = [real(sound_front_W2(Du, Dε, τu, t; branch=:plus))  for t in τXs]
Wm = [real(sound_front_W2(Du, Dε, τu, t; branch=:minus)) for t in τXs]
Wp_im = [abs(imag(sound_front_W2(Du, Dε, τu, t; branch=:minus))) for t in τXs]
band!(axA, [lo, hi], [-0.2, -0.2], [1.3, 1.3]; color=(:seagreen, 0.12))
hlines!(axA, [1.0]; color=:black, linestyle=:dash, linewidth=1.2, label="luminal W²=1")
lines!(axA, τXs, Wp; color=:crimson,  linewidth=2.5, label="W²₊ (larger branch)")
lines!(axA, τXs, Wm; color=:dodgerblue, linewidth=2.5, label="W²₋ (smaller branch)")
vlines!(axA, [lo]; color=:seagreen, linewidth=1.5, linestyle=:dot)
vlines!(axA, [hi]; color=:seagreen, linewidth=1.5)
text!(axA, lo, 1.15; text="2D_u", color=:seagreen, fontsize=11, align=(:center,:bottom))
text!(axA, hi, 1.15; text="τ_X^max", color=:seagreen, fontsize=11, align=(:center,:bottom))
ylims!(axA, -0.15, 1.3); xlims!(axA, 0.0, hi*1.25)
axislegend(axA; position=:lt, framevisible=true, labelsize=10)

# ── Panel B: Alfvén + magnetosonic |Re W(θ)| on OT, τ_X=0 vs 0.2 ─────────────
axB = Axis(fig[1,2], xlabel="propagation angle  θ  [rad]", ylabel="front velocity  |Re W(θ)|",
           title="(B) Alfvén+magnetosonic branches on Orszag–Tang\nτ_X=0 (unstable, Im W≠0) vs τ_X=0.2 (real)")
b_ot = ot_max_b()
θs = range(0, π/2; length=200)
for (τX, col, lab) in ((0.0, :firebrick, "τ_X=0"), (0.2, :navy, "τ_X=0.2"))
    c = bdnk_coeffs_from_Dmaps(Du=1e-2, Dε=2e-3, rb=1e-2, τu=2e-1, τX=τX,
                               τb=8e-2, ε=30.0, b2=b_ot^2)
    # plot the max |Re W| over all branches at each θ
    reW = Float64[]; imW = Float64[]
    for θ in θs
        Ws = ComplexF64[]
        xp, xm = alfven_x(c, θ)
        push!(Ws, sqrt(complex(xp))); push!(Ws, sqrt(complex(xm)))
        for x in magnetosonic_x(c, θ); push!(Ws, sqrt(complex(x))); end
        push!(reW, maximum(abs ∘ real, Ws))
        push!(imW, maximum(abs ∘ imag, Ws))
    end
    lines!(axB, collect(θs), reW; color=col, linewidth=2.5, label="$lab  |Re W|")
    if maximum(imW) > 1e-6
        lines!(axB, collect(θs), imW; color=col, linewidth=2.0, linestyle=:dash,
               label="$lab  |Im W| (anti-diffusive)")
    end
    r = front_velocity_max(c; nθ=721)
    scatter!(axB, [r.θ_at_vmax], [r.vmax]; color=col, markersize=11, marker=:star5)
    @printf("OT τ_X=%.1f : v_max=%.5f  max|ImW|=%.5f\n", τX, r.vmax, r.max_imW)
end
hlines!(axB, [1.0]; color=:black, linestyle=:dot, linewidth=1.0)
ylims!(axB, 0.0, 1.05)
axislegend(axB; position=:lt, framevisible=true, labelsize=9)

# ── Panel C: telegrapher numerical-vs-analytic + N⁻² convergence ─────────────
gC = fig[1,3] = GridLayout()
axC1 = Axis(gC[1,1], xlabel="x", ylabel="b(t,x)",
            title="(C) Boosted-telegrapher benchmark (Eq.21)")
τb, rpb, k, tend = 0.5, 0.4, 3.0, 1.0
xg, bnum, _ = solve_telegrapher_1d(256, tend, k, τb, rpb)
ban = [telegrapher_analytic(tend, xi, k, τb, rpb) for xi in xg]
lines!(axC1, xg, ban; color=:black, linewidth=2.5, label="analytic")
scatter!(axC1, xg[1:8:end], bnum[1:8:end]; color=:crimson, markersize=6, label="numerical N=256")
axislegend(axC1; position=:rb, framevisible=true, labelsize=9)

axC2 = Axis(gC[2,1], xlabel="N", ylabel="L1 error", xscale=log10, yscale=log10,
            title="2nd-order convergence (∝ N⁻²)")
Ns = [32, 64, 128, 256, 512, 1024]
errs, orders = telegrapher_convergence(Ns, tend, k, τb, rpb)
scatter!(axC2, Ns, errs; color=:dodgerblue, markersize=10, label="L1 error")
ref = errs[1] .* (Ns[1] ./ Ns).^2
lines!(axC2, Ns, ref; color=:black, linestyle=:dash, linewidth=1.5, label="∝ N⁻²")
axislegend(axC2; position=:lb, framevisible=true, labelsize=9)
@printf("telegrapher observed orders: %s\n", string(round.(orders, sigdigits=4)))

Label(fig[0, :], "BDNKStar — reproduce Lier–Armas–Porth (2026, 2606.22691): BDNK viscoresistive MHD causality / front velocities",
      fontsize=14, font=:bold)

outfile = joinpath(outdir, "bdnk_mhd_causality.png")
save(outfile, fig; px_per_unit=2)
println("wrote ", outfile)

# ── console validation table ─────────────────────────────────────────────────
println("\n── KH initial v_max vs paper (Table I) ──")
b_kh = kh_initial_b()
khpar = [(:a,1e-4,5e-5,1e-4,1e-3,5e-4,5e-4,0.9098),
         (:b,1e-4,5e-5,1e-3,1e-3,5e-4,5e-3,0.9098),
         (:c,1e-3,5e-4,1e-4,8e-3,5e-3,5e-4,0.9503),
         (:d,1e-3,5e-4,1e-3,8e-3,5e-3,5e-3,0.9503)]
for (nm,Du,Dε,rb,τu,τX,τb,paper) in khpar
    c = bdnk_coeffs_from_Dmaps(Du=Du,Dε=Dε,rb=rb,τu=τu,τX=τX,τb=τb,ε=1.0,b2=b_kh^2)
    v = front_velocity_max(c; nθ=721).vmax
    @printf("  KH-%s: computed=%.4f  paper=%.4f  Δ=%+.4f\n", nm, v, paper, v-paper)
end
