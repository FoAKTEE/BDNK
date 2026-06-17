# Reproduce Pandya arXiv:2201.12317:
#   (1) CC_plot intent / eq:frame_coeffs (line 458): the PMP luminal conformal
#       frame (chi0,lambda0)=(25/4,25/7)*eta0 "fixes the characteristic speeds
#       to be exactly unity" (= c = 1).
#   (2) Tab_cons (Fig. fig:Tab_cons, sec:constraint_tests, lines 1060-1110):
#       discrete conservation -- the finite-volume conformal-BDNK scheme
#       conserves the spatial integrals of T^{tt}, T^{tx} to machine precision
#       (~1e-15) until the pulse reaches the boundary; periodic BCs make the
#       discrete conservation exact for all time.
#
# Grounding:
#   eq:frame_coeffs_def (lines 440-447): eta = eta0 e^{3/4}; tau_eps = 3/(4e) chi,
#       chi = chi0 e^{3/4}; tau_Q = 3/(4e) lambda, lambda = lambda0 e^{3/4}.
#   conformal conditions (lines 425-428): P=e/3, Pi=A/3, zeta=0,
#       beta_eps = tau_Q/3, beta_n=0, tau_P = tau_eps/3.
#   eq:frame_coeffs (line 457): (chi0,lambda0) = (25/4 eta0, 25/7 eta0).
#   eq:gaussian_ID (lines ~1075): eps = A e^{-x^2/w^2}+delta, u^x=0,
#       A=1, w=25, delta=0.1, L=200, eta0=0.2.

include("/data/haiyangw/claude/BDNK/code/BDNKStar/src/BDNKStar.jl")
using .BDNKStar
using Printf

println("="^72)
println("PART 1: characteristic speeds of the PMP luminal conformal frame")
println("="^72)

using .BDNKStar.ConformalBDNK
using .BDNKStar.ConformalBDNK: T_tt, T_tx, T_xx

# -------------------------------------------------------------------------
# Conformal-BDNK characteristic speeds from the PRINCIPAL SYMBOL of the
# authoritative ported solver (ConformalBDNK.T_tt/T_tx/T_xx; solver.c, Pandya
# arXiv:2201.12317). Conserved q=(T^tt,T^tx); slab conservation laws
#     d_t T^tt + d_x T^tx = 0,   d_t T^tx + d_x T^xx = 0.
# Each T^{ab} is LINEAR in the first-derivative primitives g in
# {xi_x,u_x,xi_t,u_t}:  T^{ab} = P^{ab}(xi,u) + sum_g c^{ab}_g g, with
#   g = d_alpha p_j,  alpha in {x,t},  p=(xi,u).
# A plane wave  p_j ~ exp(i k (x - v t))  sends d_t->(-ikv)=ik*s_t, d_x->ik*s_x
# with s_t=-v, s_x=1, so the SECOND-derivative (principal) terms give a 2x2
# symbol S(v) with
#   S[eq,j] = sum_alpha s_alpha ( s_t c^{(top)}_{(alpha,j)} + s_x c^{(bot)}_{(alpha,j)} )
# (top/bot = the two T-components in each conservation eq). The characteristic
# speeds are the roots of det S(v) = 0 (a quartic in v). eq:frame_coeffs (line
# 457-458) claims the chosen frame fixes the (maximum) speeds to unity.
# -------------------------------------------------------------------------
function conf_char_speeds(fr::ConformalFrame; xi0=0.0, u0=0.0)
    Tfs = (T_tt, T_tx, T_xx)
    # arg indices into (xi,u,xi_x,u_x,xi_t,u_t): 3=xi_x,4=u_x,5=xi_t,6=u_t
    dTdarg(Tf, ai) = begin
        h = 1e-6; base = [xi0, u0, 0.0, 0.0, 0.0, 0.0]
        a = copy(base); a[ai] += h; b = copy(base); b[ai] -= h
        (Tf(fr, a...) - Tf(fr, b...)) / (2h)
    end
    gmap = [(3,:x,1), (4,:x,2), (5,:t,1), (6,:t,2)]   # (argidx, alpha, j)
    c = Dict{Tuple{Int,Symbol,Int},Float64}()
    for comp in 1:3, (ai, al, j) in gmap
        c[(comp, al, j)] = dTdarg(Tfs[comp], ai)
    end
    detS(v) = begin
        st = -v; sx = 1.0; S = zeros(2,2)
        for j in 1:2, (eq, top, bot) in ((1,1,2), (2,2,3))
            s = 0.0
            for (al, sal) in ((:t, st), (:x, sx))
                s += sal * (st*get(c,(top,al,j),0.0) + sx*get(c,(bot,al,j),0.0))
            end
            S[eq, j] = s
        end
        S[1,1]*S[2,2] - S[1,2]*S[2,1]
    end
    roots = Float64[]
    vs = range(-3, 3; length=600001); prev = detS(vs[1])
    for k in 2:length(vs)
        cur = detS(vs[k])
        if prev == 0 || (prev < 0) != (cur < 0)
            a = vs[k-1]; b = vs[k]; fa = detS(a)
            for _ in 1:80
                m = (a+b)/2; fm = detS(m)
                (fa < 0) != (fm < 0) ? (b = m) : (a = m; fa = fm)
            end
            push!(roots, (a+b)/2)
        end
        prev = cur
    end
    return sort(roots)
end

eta0 = 0.2
frL = ConformalFrame(eta0, (25/7)*eta0, (25/4)*eta0)   # (eta0, lambda0, chi0) PMP luminal

maxdev = 0.0
println("\n  e         u        char speeds v                              max|v|")
for xi0 in (log(0.1), 0.0, log(5.0), log(100.0)), u0 in (0.0, 0.3, 0.8)
    r = conf_char_speeds(frL; xi0=xi0, u0=u0)
    mx = maximum(abs.(r))
    global maxdev = max(maxdev, abs(mx - 1.0))
    @printf("  %-8.3g  %.1f    %-40s  %.12f\n",
            exp(xi0), u0, join([@sprintf("%+.4f", x) for x in r], " "), mx)
end
# control: a NON-luminal frame (chi0,lambda0)=(1,1)*eta0 must be superluminal
rbad = conf_char_speeds(ConformalFrame(eta0, eta0, eta0))
println("\n  control  non-luminal (chi0,lambda0)=(1,1)eta0:  v=",
        join([@sprintf("%+.4f", x) for x in rbad], " "), "  max|v|=", maximum(abs.(rbad)))
println("\n  max |v_max - 1| over PMP-luminal samples = ", maxdev,
        "   (eq:frame_coeffs line 458: speeds = unity; a=1, lines 494/911)")
luminal_ok = maxdev < 1e-6

println("\n", "="^72)
println("PART 2: discrete conservation of T^{tt}, T^{tx} (Tab_cons)")
println("="^72)

using .BDNKStar.ConformalEvolution

# Gaussian ID per eq:gaussian_ID (A=1, w=25, delta=0.1, L=200), eta0=0.2.
# The package init_gaussian uses exp(-(x-x0)^2/w); to match exp(-x^2/w^2) we
# pass w_code = w^2 = 625. eta0 enters via the PMP luminal frame at the
# reference energy scale (the frame's overall eta0 ~ e^{1/4}/(3pi)); to set the
# paper's eta0=0.2 we build the ConformalFrame directly with eta0=0.2.
fr = frL   # reuse the PMP luminal frame built in Part 1 (eta0=0.2)

# helper: trapezoid-free simple sum * dx over INTERIOR (non-ghost) cells
function integrals(s)
    N = length(s.x); NG = 3
    Itt = 0.0; Itx = 0.0
    for i in NG+1:N-NG
        Itt += s.Ttt[i]; Itx += s.Ttx[i]
    end
    return Itt*s.dx, Itx*s.dx
end

# ---- (a) PERIODIC smooth test: exact discrete conservation for all time ----
println("\n[a] periodic smooth (cosine-bump) test -- discrete conservation exact")
Nper = 256
xper = collect(range(-200.0, 200.0; length=Nper)); dxp = xper[2]-xper[1]
sp = ConfState(fr, xper, dxp, 0.1*dxp,
               zeros(Nper),zeros(Nper),zeros(Nper),zeros(Nper),
               zeros(Nper),zeros(Nper),zeros(Nper),zeros(Nper), true)
for i in 1:Nper
    # smooth periodic energy bump (period = domain length), u=0
    eps = 0.1 + 0.5*(1 + cos(2π*xper[i]/400.0))^2
    sp.ξ[i] = log(eps); sp.u[i] = 0.0
    sp.Ttt[i] = ConformalBDNK.T_tt(fr, sp.ξ[i], 0.0,0,0,0,0)
    sp.Ttx[i] = ConformalBDNK.T_tx(fr, sp.ξ[i], 0.0,0,0,0,0)
end
# initialize aux (call evolve! 0 steps does not; mirror _update_aux via 1 dummy step? use evolve 0)
Itt0, Itx0 = integrals(sp)
nsteps = 2000
evolve!(sp, nsteps)
Itt1, Itx1 = integrals(sp)
rel_tt_p = abs(Itt1-Itt0)/abs(Itt0)
abs_tx_p = abs(Itx1-Itx0)            # Itx0 ~ 0 (u=0 symmetric), use absolute
@printf("  steps=%d  dt=%.4g  T_final=%.3g\n", nsteps, sp.dt, nsteps*sp.dt)
@printf("  integral T^{tt}:  init=%.15e  final=%.15e\n", Itt0, Itt1)
@printf("  integral T^{tx}:  init=%.3e  final=%.3e\n", Itx0, Itx1)
@printf("  rel drift |dI_tt|/|I_tt| = %.3e\n", rel_tt_p)
@printf("  abs drift |dI_tx|        = %.3e\n", abs_tx_p)
cons_periodic_ok = (rel_tt_p < 1e-12) && (abs_tx_p < 1e-12)
println("  PERIODIC conservation to ~1e-12: ", cons_periodic_ok)

# ---- (b) PAPER Gaussian ID (outflow): conserved until pulse hits boundary ----
println("\n[b] paper Gaussian ID (eq:gaussian_ID, outflow) -- conserved pre-boundary")
sg = init_gaussian(fr; N=257, xmin=-200.0, xmax=200.0, A=1.0, x0=0.0, w=625.0, c=0.1, cfl=0.1)
Itt0g, Itx0g = integrals(sg)
# Evolve a time well short of the pulse front reaching |x|=200. Initial pulse
# support half-width ~48 (eps falls to ~delta there); fronts move at speed<=c=1,
# so the pulse stays interior for t < 200-48 ~ 150. Use ~78 to be safe; the
# paper (lines 1084-1090) states T^tt is conserved to ~1e-15 at such times.
nstepsg = 500
evolve!(sg, nstepsg)
Itt1g, Itx1g = integrals(sg)
rel_tt_g = abs(Itt1g-Itt0g)/abs(Itt0g)
abs_tx_g = abs(Itx1g-Itx0g)
@printf("  steps=%d  dt=%.4g  T_final=%.3g  (boundary reach ~ t>=%.0f)\n",
        nstepsg, sg.dt, nstepsg*sg.dt, 200.0-48.0)
@printf("  integral T^{tt}:  init=%.15e  final=%.15e\n", Itt0g, Itt1g)
@printf("  rel drift |dI_tt|/|I_tt| = %.3e\n", rel_tt_g)
@printf("  abs drift |dI_tx|        = %.3e\n", abs_tx_g)
cons_gauss_ok = rel_tt_g < 1e-12
println("  GAUSSIAN (outflow) pre-boundary conservation to ~1e-12: ", cons_gauss_ok)

println("\n", "="^72)
println("SUMMARY")
println("  luminal char speeds (=1):           ", luminal_ok, "   max|c-1|=", maxdev)
println("  periodic discrete conservation:     ", cons_periodic_ok,
        "   rel_tt=", rel_tt_p, " abs_tx=", abs_tx_p)
println("  Gaussian pre-boundary conservation: ", cons_gauss_ok, "   rel_tt=", rel_tt_g)
println("  ALL PASS: ", luminal_ok && cons_periodic_ok && cons_gauss_ok)
println("="^72)
