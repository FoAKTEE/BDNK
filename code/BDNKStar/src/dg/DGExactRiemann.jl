#=
    DGExactRiemann — EXACT special-relativistic Riemann solver for a Γ-law ideal
    gas (Martí & Müller 1994, 2003 living review). Used to validate the Stage-1
    RKDG shock tube: the contact pressure p* is the root of the function matching
    the left- and right-going wave velocity functions; the full self-similar
    solution u(x/t) is then assembled (rarefaction fan / shock / contact).

    For the standard tests (initial velocities zero) the relevant pieces are:
      * shock relations via the relativistic Rankine–Hugoniot jump (Taub adiabat),
        giving the post-shock velocity as a function of p*;
      * rarefaction via the Riemann invariant (integral of dp/(ρh√(...)) along the
        isentrope), giving the velocity behind a self-similar fan.
    We solve p* by bisection on Δv(p*) = v3(p*) − v4(p*) = 0 (3=left state piece,
    4=right state piece) and then sample the solution at requested x/t.
=#
module DGExactRiemann

export exact_riemann_sr, sample_riemann_sr, RiemannSol

# enthalpy / sound speed for Γ-law from (ρ,p)
@inline _h(ρ,p,Γ) = 1.0 + Γ/(Γ-1.0)*p/ρ
@inline _cs(ρ,p,Γ) = sqrt(Γ*p/(ρ*_h(ρ,p,Γ)))

# ---------------------------------------------------------------------------
# Relativistic shock: given a known state (ρa,pa) and post-shock pressure p,
# return the post-shock velocity (in the frame where the unshocked fluid is at
# rest) and post-shock density, via the Taub adiabat + shock jump.
# Martí–Müller (2003) eqs. (for a shock moving into state a).
# ---------------------------------------------------------------------------
# Post-shock density (Taub adiabat) for pressure p shocking into (ρa,pa).
function _shock_density(ρa, pa, Γ, p)
    g = Γ/(Γ-1.0); ha=_h(ρa,pa,Γ)
    # Taub adiabat h² - ha² = (p-pa)(h/ρ + ha/ρa), ρ = g p/(h-1).
    f(h)=begin ρ=g*p/(h-1.0); h^2-ha^2-(p-pa)*(h/ρ+ha/ρa) end
    hlo=1.0+1e-12; hhi=max(ha*100.0,10.0); flo=f(hlo); fhi=f(hhi); it=0
    while flo*fhi>0 && it<200; hhi*=2; fhi=f(hhi); it+=1 end
    for _ in 1:300
        h=0.5*(hlo+hhi); fh=f(h)
        if flo*fh≤0; hhi=h; fhi=fh else hlo=h; flo=fh end
        (hhi-hlo)<1e-15 && break
    end
    h=0.5*(hlo+hhi)
    return g*p/(h-1.0)
end

# Post-shock state shocking into (ρa,pa) at rest. Returns (|Δv|, ρ): the speed of
# the shocked fluid RELATIVE to the (at-rest) pre-shock state. Martí–Müller
# (2003) relative-velocity across a shock:
#   v₁₂ = √[ (p−pa)(e−ea) / ((ea+p)(e+pa)) ]
# with e the TOTAL energy density (e = ρ + p/(Γ−1)).
function _shock_state(ρa, pa, Γ, p)
    ρ = _shock_density(ρa, pa, Γ, p)
    ea = ρa + pa/(Γ-1.0)
    e  = ρ  + p /(Γ-1.0)
    num = (p-pa)*(e-ea)
    den = (ea+p)*(e+pa)
    v = sqrt(max(num/den, 0.0))
    return v, ρ
end

# ---------------------------------------------------------------------------
# Rarefaction: integrate the Riemann invariant from (ρa,pa) down to pressure p,
# returning the velocity (a at rest) and density behind the fan. ODE in p:
#   dv/dp = ± 1/(ρ h c_s γ²(1∓vc_s)... ) — we integrate the invariant in terms of
# the self-similar variable using the standard relativistic isentrope ρ∝p^{1/Γ}.
# We integrate dv = ∓ √(1-v²) /(ρ h cs) dp / (1 - v cs ...) — use the known
# closed form for the velocity through a fan:
#   v(p) via the relativistic Riemann invariant J = (1/2)ln[(1+v)/(1-v)] ∓
#          (1/√(Γ-1)) ln[(√(Γ-1)+cs)/(√(Γ-1)-cs)].
# ---------------------------------------------------------------------------
function _rarefaction_state(ρa, pa, Γ, p; left::Bool)
    # isentrope: p/ρ^Γ const ⇒ ρ = ρa (p/pa)^{1/Γ}. Integrate the relativistic
    # Riemann invariant ODE for the velocity (pre-state at rest, va=0):
    #   d(artanh v)/dp = ∓ cs/( ρ h (1−v²)·??? )  — use the exact simple-wave
    # relation  d[½ln((1+v)/(1−v))] = ± cs/(ρ h) · 1/√(1−v²)·... . The clean,
    # frame-correct integrand (Rezzolla–Zanotti eq. 4.95) is
    #   dv/dp = ± (1/(ρ h W²)) · 1/(cs) ... ; we integrate the unambiguous
    #   d(artanh v) = ± cs/( (ε+p) ) dp   along the isentrope (ε+p = ρ h).
    # Here ε+p = ρ h is the total enthalpy density. Sign: LEFT wave ⇒ v rises as
    # p falls (expansion accelerates fluid forward) ⇒ +.
    ρ  = ρa*(p/pa)^(1.0/Γ)
    csa = _cs(ρa,pa,Γ); cs = _cs(ρ,p,Γ)
    sΓ = sqrt(Γ-1.0)
    # closed-form Γ-law relativistic rarefaction (Martí–Müller 2003 eq. 32-33),
    # pre-state at rest (va=0):
    #   (1+v)/(1-v) = [ ((√(Γ-1)-cs)(√(Γ-1)+csa)) / ((√(Γ-1)+cs)(√(Γ-1)-csa)) ]^{2/√(Γ-1)}
    base = ((sΓ-cs)*(sΓ+csa)) / ((sΓ+cs)*(sΓ-csa))
    R = base^(2.0/sΓ)
    v = (R-1.0)/(R+1.0)
    return v, ρ
end

# |Δv| jump magnitude across a wave (shock or rarefaction) at contact pressure p.
# This is the relativistic relative velocity between the pre-state and the star
# (contact) region. For rarefactions, _rarefaction_state already returns the
# velocity relative to the (at-rest) pre-state.
function _wave_jump(ρa, pa, Γ, p; left::Bool)
    if p > pa
        v, _ = _shock_state(ρa,pa,Γ,p)
        return abs(v)
    else
        v, _ = _rarefaction_state(ρa,pa,Γ,p; left=left)
        return abs(v)
    end
end

# Contact velocity v* implied by the LEFT family at pressure p, with vL=0:
#   the left wave moves the fluid to +x by the jump magnitude.
# (Both initial velocities are zero; the matching is v*_L = v*_R.)
function _wave_velocity(ρa, pa, Γ, p; left::Bool)
    jmp = _wave_jump(ρa,pa,Γ,p; left=left)
    return left ? jmp : -jmp
end

struct RiemannSol
    Γ::Float64
    ρL::Float64; vL::Float64; pL::Float64
    ρR::Float64; vR::Float64; pR::Float64
    pstar::Float64
    vstar::Float64
    x0::Float64
end

"""
    exact_riemann_sr(ρL,vL,pL, ρR,vR,pR; Γ=5/3, x0=0.5) -> RiemannSol

Exact special-relativistic Riemann solution for a Γ-law gas with INITIALLY
ZERO velocities (vL=vR=0; the standard shock-tube/blast tests). Solves the
contact pressure p* by bisection. (Nonzero initial velocities are not boosted
here — the validation tests all use v=0 states.)
"""
function exact_riemann_sr(ρL,vL,pL, ρR,vR,pR; Γ::Float64=5/3, x0::Float64=0.5)
    @assert abs(vL) < 1e-12 && abs(vR) < 1e-12 "exact solver assumes v=0 initial states"
    # both pre-states at rest ⇒ contact velocity from each family is the jump
    # magnitude; match jmp_L(p) = jmp_R(p).
    Δv(p) = _wave_jump(ρL,pL,Γ,p; left=true) - _wave_jump(ρR,pR,Γ,p; left=false)
    plo = min(pL,pR)*(1+1e-12); phi = max(pL,pR)*(1-1e-12)
    flo=Δv(plo); fhi=Δv(phi)
    # If the crossing is outside [min,max] (double shock / double rarefaction),
    # expand the bracket.
    it=0
    while flo*fhi>0 && it<200
        phi*=2; fhi=Δv(phi); it+=1
        phi>1e14 && break
    end
    p=0.5*(plo+phi)
    for _ in 1:300
        p=0.5*(plo+phi); fp=Δv(p)
        if flo*fp ≤ 0; phi=p; fhi=fp else plo=p; flo=fp end
        (phi-plo) < 1e-13*(1+p) && break
    end
    vstar = _wave_jump(ρL,pL,Γ,p; left=true)   # contact velocity (positive ⇒ rightward)
    return RiemannSol(Γ,ρL,vL,pL,ρR,vR,pR,p,vstar,x0)
end

# ---------------------------------------------------------------------------
# Sample the exact solution at coordinate x and time t>0. Returns (ρ,v,p).
# Self-similar in ξ=(x-x0)/t. We classify by the wave fans.
# ---------------------------------------------------------------------------
function sample_riemann_sr(sol::RiemannSol, x::Float64, t::Float64)
    Γ=sol.Γ; ξ = t>0 ? (x-sol.x0)/t : (x<sol.x0 ? -1e9 : 1e9)
    pstar=sol.pstar; vstar=sol.vstar
    # LEFT family
    if pstar > sol.pL
        # left shock: speed Vs (into L). compute from shock state
        vsh, ρsh = _shock_state(sol.ρL,sol.pL,Γ,pstar)
        # shock speed: from mass flux; reuse: Vs solves continuity. Approx by
        # the relativistic shock speed via post/pre states.
        Vs = _left_shock_speed(sol.ρL,sol.pL,Γ,pstar)
        if ξ ≤ Vs
            return (sol.ρL, 0.0, sol.pL)
        end
        ρ3 = ρsh
        # between shock and contact: uniform star-left state
        left_inside = (ξ, )  # fallthrough handled below by contact test
        if ξ ≤ vstar
            return (ρ3, vstar, pstar)
        end
    else
        # left rarefaction fan
        csL = _cs(sol.ρL,sol.pL,Γ)
        ξhead = (0.0 - csL)/(1.0 - 0.0*csL)   # vL=0
        ρ3 = sol.ρL*(pstar/sol.pL)^(1.0/Γ)
        cs3 = _cs(ρ3,pstar,Γ)
        ξtail = (vstar - cs3)/(1.0 - vstar*cs3)
        if ξ ≤ ξhead
            return (sol.ρL, 0.0, sol.pL)
        elseif ξ < ξtail
            return _left_fan_sample(sol.ρL,sol.pL,Γ,ξ)
        elseif ξ ≤ vstar
            return (ρ3, vstar, pstar)
        end
    end
    # RIGHT family (ξ > vstar)
    if pstar > sol.pR
        vsh, ρsh = _shock_state(sol.ρR,sol.pR,Γ,pstar)
        Vs = _right_shock_speed(sol.ρR,sol.pR,Γ,pstar)
        if ξ ≥ Vs
            return (sol.ρR, 0.0, sol.pR)
        else
            return (ρsh, vstar, pstar)
        end
    else
        csR=_cs(sol.ρR,sol.pR,Γ)
        ξhead=(0.0+csR)/(1.0+0.0*csR)
        ρ4=sol.ρR*(pstar/sol.pR)^(1.0/Γ)
        cs4=_cs(ρ4,pstar,Γ)
        ξtail=(vstar+cs4)/(1.0+vstar*cs4)
        if ξ ≥ ξhead
            return (sol.ρR,0.0,sol.pR)
        elseif ξ > ξtail
            return _right_fan_sample(sol.ρR,sol.pR,Γ,ξ)
        else
            return (ρ4, vstar, pstar)
        end
    end
end

# Shock speed in the lab frame (pre-state at rest). From mass flux j across the
# shock: j² = (p-pa)/(Va - V) with V=h/ρ (Taub volume), and the shock velocity
# Vs = ± j/(ρa Wa √(1 + j²/(ρa²))) (pre-state at rest ⇒ Wa=1). Sign: '+' = right.
function _shock_speed_mag(ρa,pa,Γ,p)
    ha=_h(ρa,pa,Γ); ρ=_shock_density(ρa,pa,Γ,p)
    Va=ha/ρa; V=_h(ρ,p,Γ)/ρ
    j2=(p-pa)/(Va-V)              # >0 for a compressive shock (V<Va)
    j2=max(j2,0.0); j=sqrt(j2)
    Da=ρa
    return j/(Da*sqrt(1.0 + j2/Da^2))
end
_left_shock_speed(ρa,pa,Γ,p)  = -_shock_speed_mag(ρa,pa,Γ,p)
_right_shock_speed(ρa,pa,Γ,p) =  _shock_speed_mag(ρa,pa,Γ,p)

# rarefaction fan sample at self-similar speed ξ (left fan: ξ=(v-cs)/(1-vcs)).
# Parametrize by pressure p∈(p*, pa); v(p),cs(p) from the closed form; find p so
# the characteristic speed equals ξ.
function _left_fan_sample(ρa,pa,Γ,ξ)
    function vcs(p)
        v,_=_rarefaction_state(ρa,pa,Γ,p; left=true)
        ρ=ρa*(p/pa)^(1.0/Γ); cs=_cs(ρ,p,Γ)
        return v, cs, ρ
    end
    # ξ as function of p: ξ(p)=(v-cs)/(1-v cs). At p=pa: v=0 ⇒ ξ=-csa (head).
    # As p↓, v↑ ⇒ ξ↑ (tail). Monotone ⇒ bisect on p.
    g(p)=begin v,cs,_=vcs(p); (v-cs)/(1-v*cs)-ξ end
    plo=1e-12*pa; phi=pa
    glo=g(plo); ghi=g(phi)
    p=0.5*(plo+phi)
    for _ in 1:200
        p=sqrt(plo*phi); gp=g(p)
        if glo*gp≤0; phi=p; ghi=gp else plo=p; glo=gp end
        (phi-plo)<1e-14*pa && break
    end
    v,cs,ρ=vcs(p)
    return (ρ, v, p)
end
function _right_fan_sample(ρa,pa,Γ,ξ)
    ρ,v,p=_left_fan_sample(ρa,pa,Γ,-ξ)
    return (ρ,-v,p)
end

end # module DGExactRiemann
