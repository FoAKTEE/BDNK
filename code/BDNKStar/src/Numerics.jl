#=
    Numerics — small, dependency-free root-finding utilities used by the
    conservative->primitive recovery. A guarded Brent solver (derivative-free,
    bracketing, superlinear) is the workhorse: primitive recovery is a 1-D root
    find in the pressure once the velocity/energy are eliminated, and Brent gives
    machine-precision roots without needing the (messy) analytic Jacobian.
=#
module Numerics

export brent, RootResult

struct RootResult
    root::Float64
    fval::Float64
    iters::Int
    converged::Bool
end

"""
    brent(f, a, b; xtol=1e-14, ftol=0.0, maxiter=200)

Brent's method on a bracket `[a, b]` with `f(a)·f(b) ≤ 0`. Returns a
`RootResult`. Combines bisection (guaranteed) with inverse-quadratic /
secant steps (fast). `xtol` is the absolute tolerance on the root.
"""
function brent(f, a::Float64, b::Float64; xtol::Float64=1e-14, ftol::Float64=0.0,
               maxiter::Int=200)
    fa = f(a); fb = f(b)
    if fa == 0.0
        return RootResult(a, 0.0, 0, true)
    elseif fb == 0.0
        return RootResult(b, 0.0, 0, true)
    end
    if fa * fb > 0.0
        # Not bracketed — report the better endpoint, flagged unconverged.
        return RootResult(abs(fa) < abs(fb) ? a : b, min(abs(fa), abs(fb)), 0, false)
    end
    if abs(fa) < abs(fb)
        a, b = b, a; fa, fb = fb, fa          # ensure |f(b)| ≤ |f(a)|
    end
    c = a; fc = fa
    d = b - a
    mflag = true
    for it in 1:maxiter
        if fb == 0.0 || abs(b - a) < xtol
            return RootResult(b, fb, it, true)
        end
        if fa != fc && fb != fc
            # inverse quadratic interpolation
            s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) +
                c*fa*fb/((fc-fa)*(fc-fb))
        else
            s = b - fb*(b-a)/(fb-fa)            # secant
        end
        cond1 = !((3a+b)/4 ≤ s ≤ b || b ≤ s ≤ (3a+b)/4)
        cond2 = mflag && abs(s-b) ≥ abs(b-c)/2
        cond3 = !mflag && abs(s-b) ≥ abs(c-d)/2
        cond4 = mflag && abs(b-c) < xtol
        cond5 = !mflag && abs(c-d) < xtol
        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (a+b)/2; mflag = true          # bisection fallback
        else
            mflag = false
        end
        fs = f(s)
        d = c; c = b; fc = fb
        if fa*fs < 0.0
            b = s; fb = fs
        else
            a = s; fa = fs
        end
        if abs(fa) < abs(fb)
            a, b = b, a; fa, fb = fb, fa
        end
        if ftol > 0.0 && abs(fs) ≤ ftol
            return RootResult(s, fs, it, true)
        end
    end
    return RootResult(b, fb, maxiter, abs(b-a) < 100*xtol)
end

end # module Numerics
