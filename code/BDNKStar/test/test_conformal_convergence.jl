using Test
using BDNKStar

# Self-convergence of the conformal BDNK evolution: a smooth Gaussian evolved on
# three aligned grids (N, 2N, 4N) with dt ∝ dx and the same final time. The
# convergence factor Q = ‖f_h−f_{h/2}‖ / ‖f_{h/2}−f_{h/4}‖ → 2^p. In the SMOOTH
# region WENO5 (5th-order space) dominates ⇒ p≈3–5; the LEDGER's Q→4 (2nd order,
# Heun-RK2-limited) is the *with-features* asymptote. We require clear, ordered
# self-convergence between 2nd and 5th order — the documented code requirement.
@testset "Conformal evolution: ordered self-convergence (2nd–5th order)" begin
    fr = pmp_luminal_frame(10.0)
    K = 40
    run(N, steps) = begin
        s = init_gaussian(fr; N=N, xmin=-100.0, xmax=100.0, A=0.3, x0=0.0, w=25.0, c=0.5, cfl=0.1)
        evolve!(s, steps)
        energy_density(s)
    end
    ε1 = run(65,  K)        # h
    ε2 = run(129, 2K)       # h/2  (dt halves ⇒ 2K steps to same time)
    ε3 = run(257, 4K)       # h/4
    I = 4:62                # interior (skip 3 ghosts each side)
    d12 = [ε1[i] - ε2[2i-1] for i in I]
    d23 = [ε2[2i-1] - ε3[4i-3] for i in I]
    Q = sqrt(sum(abs2, d12)) / sqrt(sum(abs2, d23))
    @info "conformal evolution self-convergence factor" Q order=log2(Q)
    @test 3.8 < Q < 36      # ordered convergence, 2nd (Q=4) … 5th (Q=32) order
end
