# phase_rig.py: Theorem cp:thm:rig of Section cp:sec:cert (A certificate for one point of the diagram) of the chapter on the Thorin phase
# diagram: the finite part of the statement (0.2, 0.6) in G, run in ball arithmetic (python-flint, Arb), with every truncation bounded
# explicitly and added to the radii. For t = 0.2 it certifies that xi_t has no zero in the closed rectangle Q_T = [0.6, 0.9] x [0, T],
# T = 60, 120 and 160, by the winding number of xi_t around the boundary of Q_T. It also checks the hypotheses of Theorem cp:thm:cert as
# the corrected theorem states them (the Cauchy radius r > 0 free at every step, the factor e^{-2 eta (theta_j - r)} replaced by 1 when
# theta_j <= r, the ball radii written eps_j so that they do not clash with r), and runs the same machinery on a positive control.
# Labels: cp:lem:bottom, cp:lem:right, cp:lem:tilt, cp:lem:deriv, cp:lem:grid, cp:thm:cert, cp:thm:rig; Table cp:tab:K (the constants).
# Method. xi_t(s) = 8 int_0^oo e^{t u^2} Phi(u) cosh(2 w u) du = 4 int_R ..., w = s - 1/2, Phi = Phi_dBN (the kernel of Section
# sec:ch21:flow, Phi(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})). The evaluator is the trapezoid rule on the
# real axis (the Polya route of polya_flint.py with the weight e^{t u^2}): the node values c_j = e^{t u_j^2} Phi(u_j), u_j = j h,
# 0 <= j <= J, are computed once, and each point costs one sum of J + 1 terms, xi_t = sum_j 8 h w_j c_j [cosh(2 w u_j)], w_0 = 1/2.
# The boundary of Q_T is walked counterclockwise, in chunks that run in parallel: the right edge up, the top edge to the left, the left
# edge down; the bottom edge is the real segment, where xi_t > 0 (cp:lem:bottom), so it contributes the argument change 0. At every step
# the step length is 0.8 |xi_t(s_j)|_low / D, D being the least over three tilts eta (0.30, 0.35, 0.38) and three radii r of the derivative
# bound K_eta(t, |x| + r) e^{-2 eta (theta - r)} / r of cp:lem:deriv, and the inequality of cp:thm:cert is then verified in ball arithmetic
# with the exact step length and the exact lowest height of the step; the argument increments Arg(B_{j+1}/B_j) are summed as balls.
# Error control. Everything below is added to the radii or checked by an inequality in ball arithmetic.
#  (1) theta series: the loop stops when the tail bound  |e^{9u}| (2 pi^2 + 3 pi / x) (4 / (pi x))^2 e^{-2} q^{(N+1)^2} / (1 - q^{N+1}),
#      q = e^{-pi x / 2}  (x a lower bound of Re e^{4u}), is below 2^{-prec-30}; it covers both terms of the kernel (the term 3 pi n^2 e^{5u}
#      is absorbed through |e^{-4u}| <= 1/x, which repairs the display of the source); the tail is added as a disc, to the real and the
#      imaginary radius (the source adds it to the real radius only).
#  (2) the cut of the u axis: nodes u_j <= U = J h with U >= 3, and the sum over the omitted nodes plus nothing else is bounded by
#      8 G(U) / lambda, G(u) = 60 exp(t u^2 + (9 + 2 sigma) u - pi e^{4u}) decreasing for u >= 3, -G'/G >= lambda = 4 pi e^{4U} - 2 t U - 9 - 2 sigma.
#  (3) the trapezoid discretisation, by Poisson summation on the strip |Im u| < a, a = 0.30 < pi/8:
#      |error| <= 2 K_a(t, sigma) cosh(2 a theta) / (e^{2 pi a / h} - 1), K_a of cp:lem:tilt, itself bounded as in (4).
#  (4) the constants K_eta(t, x) = 8 int_0^oo e^{t (u^2 - eta^2)} |Phi(u + i eta)| cosh(2 x u) du (the even form): upper bounds by
#      interval cells of width 2^{-k} on [0, 3], plus the tail u >= 3 bounded by 8 G(3) / lambda with the kernel bound
#      |Phi(u + i eta)| <= 60 e^{9u} exp(-pi cos(4 eta) e^{4u}), u >= 3 (the omitted tail of the source, review item M7).
#  (5) the grid step of cp:lem:grid: the inequality  |s_{j+1} - s_j| K_eta(t, x_bar + r) e^{-2 eta (theta_j - r)} / r + eps_j + eps_{j+1} < min |B_j|
#      is checked in ball arithmetic at every step, theta_j = min(Im s_j, Im s_{j+1}), the exponential replaced by 1 when theta_j <= r.
#  (6) the arithmetic: all quantities are balls; the grid points are exact dyadic numbers except for the abscissas 0.1 and 0.4 of sigma,
#      which carry a radius of 2^{-prec} (a ball containing the real number).
# Positive control: the same code on Q' = [0.4, 0.6] x [0, T_c], which contains the zeros of xi_t on the critical line below T_c;
# the certified winding must be their number (13 for T_c = 60). A certified count of sign changes of xi_t on the line (the intermediate
# value theorem on a grid of step 1/4) is reported for each height. The evaluator is compared with Arb's Gauss Legendre integrator at
# a few points (the two routes share only the kernel bound).
# Not checked here: the barrier X of Theorem cp:thm:cert (a hypothesis from Polymath15, in the variable of H_t, height X/2 for xi_t);
# the run certifies the rectangle up to T only.
# Usage: python3 phase_rig.py [--quick] [--out FILE] [--workers N] [--heights 60,120,160] [--control 60] [--kcells K]
# Output: phase_rig.json (default: ../data/phase_rig.json from this file): every certified quantity as a string with its radius, and a
# pass/fail per hypothesis. The full case is for the cluster; --quick runs Q_60 and the control Q'_60 with coarser bounds of K in a few seconds.
import argparse
import json
import math
import multiprocessing as mpr
import os
import sys
import time

from flint import acb, arb, ctx

TT = "0.2"
ETAS = ["0.30", "0.35", "0.38"]
RADS = ["0.25", "0.375", "0.5"]
STRIP = "0.30"
KPREC = 128
EV = None


# helpers
def enc(x):
    m, e = x.mid().man_exp()
    rm, re_ = x.rad().man_exp()
    return (int(m), int(e), int(rm), int(re_))


def dec(t):
    return arb(arb(t[0]) * arb(2) ** t[1], arb(t[2]) * arb(2) ** t[3])


def sball(x, n=30):
    return x.str(n, radius=True)


def sup(x, n=20):
    return x.upper().str(n, radius=False)


def xval(s):
    return sum((arb(p) for p in s.split("+")), arb(0))


def eps(z):
    return (z.real.rad() ** 2 + z.imag.rad() ** 2).sqrt().upper()


def contains_zero(z):
    return bool(z.real.lower() <= 0 and z.real.upper() >= 0 and z.imag.lower() <= 0 and z.imag.upper() >= 0)


# the theta series
def theta_tail(abs_e9, xlow, N):
    """Upper bound for sum_{n > N} |(2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})|, given |e^{9u}| <= abs_e9 and
    0 < xlow <= Re e^{4u}: |term_n| <= |e^{9u}| (2 pi^2 + 3 pi / x) n^4 e^{-pi n^2 x}, n^4 e^{-pi n^2 x / 2} <= (4 / (pi x))^2 e^{-2},
    and n^2 - (N+1)^2 >= 2 (N+1) (n - N - 1) for the geometric sum."""
    pi = arb.pi()
    m = (4 / (pi * xlow)) ** 2 * arb(-2).exp()
    q1 = (-pi * (N + 1) ** 2 * xlow / 2).exp()
    q2 = (-pi * (N + 1) * xlow).exp()
    return abs_e9 * (2 * pi ** 2 + 3 * pi / xlow) * m * q1 / (1 - q2)


def phi(u, prec, cplx):
    """Enclosure of Phi_dBN(u): u an arb >= 0 (cplx False) or an acb (cplx True; indeterminate when the lower bound of Re e^{4u} is
    <= 0.01, which the integrator handles by subdivision). The tail is added as a disc."""
    pi = arb.pi()
    if cplx:
        xlow = ((4 * u.real).exp() * (4 * u.imag).cos()).lower()
        if not (xlow > arb("0.01")):
            return acb(arb("nan"), arb("nan"))
        abs_e9 = (9 * u.real).exp().upper()
    else:
        xlow = (4 * u).exp().lower()
        abs_e9 = (9 * u).exp().upper()
    xl = float(xlow.mid())
    N = max(4, int(math.sqrt(2.0 * (prec + 40) * 0.6931471805599453 / (math.pi * xl))) + 2)
    e4 = (4 * u).exp()
    e5 = (5 * u).exp()
    e9 = (9 * u).exp()
    s = acb(0) if cplx else arb(0)
    for n in range(1, N + 1):
        s += (2 * pi ** 2 * n ** 4 * e9 - 3 * pi * n ** 2 * e5) * (-pi * n * n * e4).exp()
    thr = arb(2) ** (-(prec + 30))
    tail = theta_tail(abs_e9, xlow, N)
    while not (tail < thr):
        N += 1
        s += (2 * pi ** 2 * N ** 4 * e9 - 3 * pi * N ** 2 * e5) * (-pi * N * N * e4).exp()
        tail = theta_tail(abs_e9, xlow, N)
    r = tail.upper()
    return s + (acb(arb(0, r), arb(0, r)) if cplx else arb(0, r))


def gtail(U, t, x, c):
    """8 G(U) / lambda >= 8 int_U^oo G, G(u) = 60 exp(t u^2 + (9 + 2x) u - pi c e^{4u}), for U >= 3 and c a lower bound of cos(4 eta)
    (c = 1 for the real axis): G is decreasing and -G'/G >= lambda there."""
    pi = arb.pi()
    e4U = (4 * U).exp()
    lam = 4 * pi * c * e4U - 2 * t * U - 9 - 2 * x
    if not (U >= 3 and lam > 0 and pi * c * e4U > 1000):
        raise RuntimeError("tail hypotheses fail")
    return 8 * 60 * (t * U * U + (9 + 2 * x) * U - pi * c * e4U).exp() / lam


# the constants K_eta(t, x)
def kjob(args):
    """Cells of width 2^-k on [0, 3]: upper bounds of e^{t (u^2 - eta^2)} |Phi(u + i eta)|, then K upper bounds for the list xs."""
    eta_s, k, xs = args
    ctx.prec = KPREC
    t = arb(TT)
    eta = arb(eta_s)
    w = arb(1) / 2 ** k
    M = 3 * 2 ** k
    A = []
    for i in range(M):
        lo = w * i
        ur = lo.union(lo + w)
        val = (t * (ur * ur - eta * eta)).exp() * abs(phi(acb(ur, eta), KPREC, True))
        up = val.upper()
        if not (up == up):
            raise RuntimeError("indeterminate cell")
        A.append(up)
    c = (4 * eta).cos().lower()
    out = {}
    for xs_ in xs:
        x = xval(xs_)
        tot = arb(0)
        for i in range(M):
            tot += A[i] * (2 * x * (w * (i + 1))).cosh().upper()
        K = 8 * w * tot + gtail(arb(3), t, x, c)
        out[xs_] = enc(K.upper())
    return eta_s, out


def compute_K(k, xs, workers):
    jobs = [(e, k, xs) for e in ETAS]
    if workers > 1:
        with mpr.Pool(min(workers, len(jobs))) as pool:
            res = pool.map(kjob, jobs)
    else:
        res = [kjob(j) for j in jobs]
    return {e: o for e, o in res}


# the evaluator
class Ev:
    def __init__(self, cfg):
        ctx.prec = cfg["prec"]
        self.prec = cfg["prec"]
        self.t = arb(cfg["t"])
        self.h = arb(cfg["h"])
        self.J = cfg["J"]
        self.a = arb(cfg["a"])
        self.pi = arb.pi()
        self.Kd = {k: dec(v) for k, v in cfg["Kdisc"].items()}
        self.u = [self.h * j for j in range(self.J + 1)]
        self.u2 = [2 * u for u in self.u]
        w8 = []
        for j, u in enumerate(self.u):
            c = (self.t * u * u).exp() * phi(u, self.prec, False)
            w8.append(8 * self.h * c * (arb(1) / 2 if j == 0 else 1))
        self.w8 = w8
        self.Ecut = gtail(self.h * self.J, self.t, arb(1), arb(1))
        self.expo = (2 * self.pi * self.a / self.h).exp() - 1
        self.cache = {}

    def tables(self, sp):
        A, B = [], []
        for w8, u2 in zip(self.w8, self.u2):
            x = sp * u2
            A.append(w8 * x.cosh())
            B.append(w8 * x.sinh())
        return A, B

    def xi(self, sp, th, key=None, big=False):
        """Ball enclosure of xi_t(1/2 + sp + i th), sp an arb (|sp| <= 0.4, or <= 1 when big), th an arb >= 0."""
        if key is not None:
            tb = self.cache.get(key)
            if tb is None:
                tb = self.cache[key] = self.tables(sp)
        else:
            tb = self.tables(sp)
        A, B = tb
        re = arb(0)
        im = arb(0)
        for a_, b_, u2 in zip(A, B, self.u2):
            s_, c_ = (th * u2).sin_cos()
            re += a_ * c_
            im += b_ * s_
        K = self.Kd["1.0" if big else "0.4"]
        Ed = 2 * K * (2 * th * self.a).cosh() / self.expo
        E = (Ed + self.Ecut).upper()
        return acb(re + arb(0, E), im + arb(0, E))


def init_worker(cfg):
    global EV
    EV = Ev(cfg)


# the march
def march(spec):
    """Walk one chunk of the boundary (kind V: vertical at abscissa sp, from height p0 to p1; kind H: horizontal at height theta, from
    abscissa p0 to p1) by the grid lemma; verify the inequality of Theorem cp:thm:cert at every step; return the sum of the Arg."""
    ev = EV
    ctx.prec = ev.prec
    combos = []
    for es, rs, kenc in spec["combos"]:
        Kb = dec(kenc)
        combos.append((float(es), float(rs), arb(es), arb(rs), Kb, float(Kb.upper().mid())))
    kind, d, key = spec["kind"], spec["d"], spec.get("key")
    if kind == "V":
        sp = dec(spec["sp"])
        pos, end = arb(spec["p0"]), arb(spec["p1"])

        def mk(p):
            return sp, p
    else:
        th = arb(spec["theta"])
        pos, end = dec(spec["p0"]), dec(spec["p1"])

        def mk(p):
            return p, th

    def best(thm):
        bd = None
        for ef, rf, eta, r, Kb, Kf in combos:
            f = math.exp(-2.0 * ef * (thm - rf)) if thm > rf else 1.0
            Dv = Kf * f / rf
            if bd is None or Dv < bd[0]:
                bd = (Dv, eta, r, Kb)
        return bd

    tot = arb(0)
    steps = fails = incr_viol = tilt_viol = 0
    maxinc = maxratio = 0.0
    big = spec.get("big", False)
    tilt_pairs, seen = [], set()
    for ef, rf, eta, r, Kb, Kf in combos:        # K at x_bar + r >= K at x_bar, for the smallest r of each tilt
        if ef not in seen:
            seen.add(ef)
            tilt_pairs.append((eta, Kb))
    sp0, th0 = mk(pos)
    v = ev.xi(sp0, th0, key, big)
    nev = 1
    while True:
        remf = d * float((end - pos).mid())
        if remf <= 0.0:
            break
        vlow = abs(v).lower()
        if not (vlow > 0):
            fails += 1
            break
        vl = float(vlow.mid())
        thc = float(th0.mid())
        hcap = 0.5
        ok = False
        for attempt in range(60):
            if kind == "V" and d < 0:
                hh0 = min(hcap, 0.8 * vl / best(thc)[0])
                hh = min(hh0, 0.8 * vl / best(max(thc - hh0, 0.0))[0])
            else:
                hh = min(hcap, 0.8 * vl / best(thc)[0])
            hh = max(hh, 1e-9)
            pos2 = end if hh >= remf * (1 - 1e-12) else pos + d * arb(hh)
            sp2, th2 = mk(pos2)
            v2 = ev.xi(sp2, th2, key, big)
            nev += 1
            L = (pos2 - pos) * d
            thmin = th2 if (kind == "V" and d < 0) else th0
            Dv, eta, r, Kb = best(float(thmin.mid()))
            expo = (-2 * eta * (thmin - r)).exp() if (thmin > r) else arb(1)
            lhs = L * Kb * expo / r + eps(v) + eps(v2)
            if lhs < vlow:
                ok = True
                break
            hcap = hh * 0.5
        if not ok:
            fails += 1
            break
        a = (v2 / v).arg()
        tot += a
        maxinc = max(maxinc, float(abs(a).upper().mid()))
        maxratio = max(maxratio, float((lhs / vlow).upper().mid()))
        # consistency checks, necessary conditions that are not part of the proof: the computed values must obey the bounds that the
        # proof uses, |xi(s_{j+1}) - xi(s_j)| <= L D (cp:lem:deriv) and |xi(s)| <= K_eta(t, x) e^{-2 eta theta} (cp:lem:tilt)
        if abs(v2 - v).lower() > L * Kb * expo / r:
            incr_viol += 1
        for etaa, Kt in tilt_pairs:
            if abs(v2).lower() > Kt * (-2 * etaa * th2).exp():
                tilt_viol += 1
        pos, v, th0 = pos2, v2, th2
        steps += 1
    return dict(steps=steps, nev=nev, fails=fails, arg=enc(tot), maxinc=maxinc, maxratio=maxratio, edge=spec["edge"], idx=spec["idx"],
                incr_viol=incr_viol, tilt_viol=tilt_viol)


def run_rect(sp_lo, sp_hi, T, combos, pool, clen=4.0):
    """Counterclockwise walk of [sp_lo, sp_hi] x [0, T] (abscissas are sigma - 1/2): right edge up, top edge, left edge down."""
    nch = max(1, int(math.ceil(T / clen)))
    b = [min(i * clen, float(T)) for i in range(nch + 1)]
    specs = []
    for i in range(nch):
        specs.append(dict(kind="V", d=1, sp=enc(sp_hi), key="R", p0=b[i], p1=b[i + 1], combos=combos, edge="R", idx=i))
    specs.append(dict(kind="H", d=-1, theta=float(T), p0=enc(sp_hi), p1=enc(sp_lo), combos=combos, edge="T", idx=0))
    for i in reversed(range(nch)):
        specs.append(dict(kind="V", d=-1, sp=enc(sp_lo), key="L", p0=b[i + 1], p1=b[i], combos=combos, edge="L", idx=i))
    res = pool.map(march, specs) if pool is not None else [march(s) for s in specs]
    out = {}
    for e in "RTL":
        rs = [r for r in res if r["edge"] == e]
        arg = sum((dec(r["arg"]) for r in rs), arb(0))
        out[e] = dict(steps=sum(r["steps"] for r in rs), arg=arg, fails=sum(r["fails"] for r in rs))
    out["maxinc"] = max(r["maxinc"] for r in res)
    out["maxratio"] = max(r["maxratio"] for r in res)
    out["incr_viol"] = sum(r["incr_viol"] for r in res)
    out["tilt_viol"] = sum(r["tilt_viol"] for r in res)
    out["fails"] = sum(r["fails"] for r in res)
    out["steps"] = sum(r["steps"] for r in res)
    out["nev"] = sum(r["nev"] for r in res)
    return out


def params_for(T, Ka):
    a = float(STRIP)
    need = 2 * T * a + math.pi * T / 4 + math.log(2 * Ka / 5e-4)
    h = math.floor(2 * math.pi * a / need * 1e4) / 1e4
    return dict(h=h, J=int(math.ceil(3.0 / h)) + 1, prec=int(1.15 * T) + 150)


def winding_report(tot, W_expected):
    """tot: ball for the sum of the Arg. Certified when tot lies within pi of 2 pi W (so that W is the winding number)."""
    pi = arb.pi()
    W = int(round(float((tot / (2 * pi)).mid())))
    cert = bool(abs(tot - 2 * pi * W).upper() < pi.lower())
    return W, cert, bool(W == W_expected)


def line_signs(ev, T, step=0.25):
    ctx.prec = ev.prec
    n = int(round(T / step))
    prev = 0
    changes = und = 0
    for k in range(n + 1):
        z = ev.xi(arb(0), arb(k * step), "C")
        re = z.real
        sgn = 1 if re.lower() > 0 else (-1 if re.upper() < 0 else 0)
        if sgn == 0:
            und += 1
            continue
        if prev != 0 and sgn != prev:
            changes += 1
        prev = sgn
    return changes, und


def crosscheck(ev, pts):
    ctx.prec = ev.prec
    out = []
    tol = arb(2) ** (-ev.prec + 20)
    Et = gtail(arb(3), ev.t, arb(1), arb(1))
    for sp_s, th_s in pts:
        sp, th = arb(sp_s), arb(th_s)
        z1 = ev.xi(sp, th, None, True)
        w = acb(sp, th)

        def f(u, analytic):
            return (ev.t * u * u).exp() * phi(u, ev.prec, True) * (2 * w * u).cosh()

        t0 = time.time()
        I = 8 * acb.integral(f, 0, 3, rel_tol=tol, abs_tol=tol) + acb(arb(0, Et.upper()), arb(0, Et.upper()))
        dm = abs(acb(z1.real.mid(), z1.imag.mid()) - acb(I.real.mid(), I.imag.mid())).upper()
        out.append(dict(sigma_minus_half=sp_s, theta=th_s, trapezoid=sball(z1.real, 20) + " + i " + sball(z1.imag, 20),
                        gauss_legendre=sball(I.real, 20) + " + i " + sball(I.imag, 20), consistent=contains_zero(z1 - I),
                        midpoint_difference=dm.str(5, radius=False), trapezoid_radius=z1.real.rad().str(5, radius=False),
                        seconds=round(time.time() - t0, 2)))
    return out


def flow_brackets(ev):
    """The five zeros of Table tab:ch21:flowzeros (printed to 5 decimals): a sign change of xi_t on the critical line at the printed value
    -/+ 6e-6 certifies a zero in that bracket (xi_t is real there)."""
    ctx.prec = ev.prec
    out = []
    for z in ("14.07681", "20.94669", "24.95221", "32.88801", "49.72988"):
        sg = []
        for sgn in (-1, 1):
            r = ev.xi(arb(0), arb(z) + sgn * arb("6e-6"), "C").real
            sg.append(1 if r.lower() > 0 else (-1 if r.upper() < 0 else 0))
        out.append(dict(printed=z, signs_at_minus_plus_6e_6=sg, zero_in_bracket=bool(sg[0] * sg[1] == -1)))
    return out


# main
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="reduced case: height 60 only and the control at 60, coarse K cells (3 * 2^9 cells)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "phase_rig.json"))
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--heights", default=None, help="comma separated, default 60,120,160 (quick: 60)")
    ap.add_argument("--control", type=int, default=None, help="height of the control rectangle (default 60)")
    ap.add_argument("--kcells", type=int, default=None, help="k: cells of width 2^-k for the bounds of K (default 14; quick 9)")
    args = ap.parse_args()
    heights = [int(x) for x in (args.heights or ("60" if args.quick else "60,120,160")).split(",")]
    Tc = args.control or 60
    kk = args.kcells or (9 if args.quick else 14)
    t_all = time.time()
    import flint
    rep = dict(script="phase_rig.py", quick=bool(args.quick), t=TT, tilts=ETAS, radii=RADS, strip_a=STRIP, kcells_log2=kk, workers=args.workers,
               python_flint=flint.__version__,
               statement="Theorem cp:thm:rig: xi_{0.2} has no zero s with 0.6 <= Re s <= 0.9 and 0 <= Im s <= T, for each T in heights below that passes")
    # 1. the constants K_eta(t, x): upper bounds
    xs = ["0.1", "0.4", "0.7", "0.9", "1.0"] + ["%s+%s" % (xb, r) for xb in ("0.4", "0.1") for r in RADS]
    t0 = time.time()
    KE = compute_K(kk, xs, args.workers)
    ctx.prec = KPREC
    rep["K_upper_bounds"] = {e: {x: sup(dec(KE[e][x])) for x in xs} for e in ETAS}
    rep["K_seconds"] = round(time.time() - t0, 1)
    print("K bounds", {e: {x: sup(dec(KE[e][x]), 6) for x in ("0.4", "0.7", "0.9")} for e in ETAS}, "seconds", rep["K_seconds"], flush=True)
    Ka = float(dec(KE[STRIP]["0.4"]).upper().mid())
    # 2. the hypotheses that are analytic: cp:lem:right at t = 0.2
    t = arb(TT)
    right = arb(1) / 2 + (1 - 2 * t).sqrt() / 2
    rep["lemma_right_bound_beta"] = sball(right, 12)
    rep["lemma_right_ok"] = bool(right < arb("0.9"))
    rep["barrier_X"] = "not computed: hypothesis of Theorem cp:thm:cert (Polymath15, variable of H_t, height X/2 of xi_t); the run covers [0.6, 0.9] x [0, T] only"
    rep["heights"] = {}
    allpass = rep["lemma_right_ok"]

    def with_pool(cfg):
        init_worker(cfg)
        if args.workers > 1:
            return mpr.Pool(args.workers, initializer=init_worker, initargs=(cfg,))
        return None

    for T in heights + [("control", Tc)]:
        control = isinstance(T, tuple)
        if control:
            T = T[1]
        P = params_for(T, Ka)
        cfg = dict(t=TT, a=STRIP, Kdisc={"0.4": KE[STRIP]["0.4"], "1.0": KE[STRIP]["1.0"]}, **P)
        t0 = time.time()
        pool = with_pool(cfg)
        ev = EV
        ctx.prec = ev.prec
        blk = dict(T=T, prec_bits=P["prec"], h=P["h"], J=P["J"], nodes=P["J"] + 1, setup_seconds=round(time.time() - t0, 1))
        # error budget at the top of the rectangle
        th_top = arb(T)
        Ed = 2 * ev.Kd["0.4"] * (2 * th_top * ev.a).cosh() / ev.expo
        blk["error_budget"] = dict(trapezoid_discretisation_at_top=sup(Ed, 6), cut_of_the_u_axis=sup(ev.Ecut, 6),
                                   theta_series="tail below 2^-(prec+30) per node, added as a disc")
        sps = (arb("-0.1"), arb("0.1")) if control else (arb("0.1"), arb(2) / 5)
        xbar = "0.1" if control else "0.4"
        combos = [(e, r, KE[e]["%s+%s" % (xbar, r)]) for e in ETAS for r in RADS]
        # bottom edge: xi_t > 0 at the corners (cp:lem:bottom is analytic)
        b0, b1 = ev.xi(sps[0], arb(0), None, False), ev.xi(sps[1], arb(0), None, False)
        blk["bottom_corners"] = [sball(b0.real, 20), sball(b1.real, 20)]
        blk["bottom_ok"] = bool(b0.real.lower() > 0 and b1.real.lower() > 0)
        t0 = time.time()
        R = run_rect(sps[0], sps[1], T, combos, pool)
        if pool is not None:
            pool.close()
            pool.join()
        tot = R["R"]["arg"] + R["T"]["arg"] + R["L"]["arg"]
        Wexp = None
        if control:
            ch, und = line_signs(ev, T)
            Wexp = ch if und == 0 else None
            blk["sign_changes_on_line"] = ch
            blk["undetermined_grid_points"] = und
        else:
            ch, und = line_signs(ev, T)
            blk["sign_changes_on_line_below_T"] = ch
            blk["undetermined_grid_points"] = und
            Wexp = 0
        W, cert, match = winding_report(tot, Wexp if Wexp is not None else 0)
        blk.update(rectangle="[%.2f, %.2f] x [0, %d]" % (0.5 + float(sps[0].mid()), 0.5 + float(sps[1].mid()), T),
                   steps=dict(right_up=R["R"]["steps"], top=R["T"]["steps"], left_down=R["L"]["steps"], total=R["steps"]),
                   evaluations=R["nev"],
                   arg_change=dict(right_up=sball(R["R"]["arg"]), top=sball(R["T"]["arg"]), left_down=sball(R["L"]["arg"]),
                                   bottom="0 (xi_t > 0 on the real axis, cp:lem:bottom)", sum=sball(tot)),
                   left_edge_arg_over_pi=sball(-R["L"]["arg"] / arb.pi(), 12),
                   winding_enclosure=sball(tot / (2 * arb.pi())), winding_nearest_integer=W, winding_certified=cert,
                   max_arg_increment=R["maxinc"], max_lhs_over_min_modulus=R["maxratio"], failed_steps=R["fails"],
                   march_seconds=round(time.time() - t0, 1))
        blk["hypotheses"] = {
            "bottom edge: xi_t > 0 at the corners (cp:lem:bottom)": blk["bottom_ok"],
            "right edge: 1/2 + sqrt(1 - 2t)/2 < 0.9 (cp:lem:right)": rep["lemma_right_ok"],
            "grid inequality of cp:thm:cert at every step (no failed step)": R["fails"] == 0,
            "sum of the Arg within pi of 2 pi W, W the winding number": cert,
            "every Arg increment within (-pi/2, pi/2)": bool(R["maxinc"] < math.pi / 2),
            "consistency: |xi(s_j+1) - xi(s_j)| <= |s_j+1 - s_j| D and |xi| <= K e^{-2 eta theta} at every grid point": bool(
                R["incr_viol"] == 0 and R["tilt_viol"] == 0)}
        blk["consistency_violations"] = dict(increment_bound=R["incr_viol"], tilt_bound=R["tilt_viol"])
        if control:
            blk["control_expected_zeros_on_line"] = Wexp
            blk["hypotheses"]["control: winding equals the number of sign changes on the line"] = bool(cert and match)
        else:
            blk["hypotheses"]["winding number is 0, so no zero in the closed rectangle"] = bool(cert and W == 0)
        ok = all(blk["hypotheses"].values())
        blk["pass"] = ok
        allpass = allpass and ok
        if control:
            rep["control"] = blk
        else:
            rep["heights"][str(T)] = blk
        print("T=%d %s steps %d fails %d winding %s (%s) W=%d cert %s seconds %.1f" % (
            T, "control" if control else "main", R["steps"], R["fails"], sball(tot / (2 * arb.pi()), 8), blk["rectangle"], W, cert, time.time() - t0),
            flush=True)
        # real axis values and the cross-check only once, on the first block
        if "xi_real_axis" not in rep:
            tab = {}
            for sg, sp_s in (("0.5", "0"), ("0.6", "0.1"), ("0.9", "0.4"), ("1.0", "0.5"), ("1.5", "1.0")):
                z = ev.xi(arb(sp_s), arb(0), None, True)
                tab[sg] = sball(z.real, 25)
            rep["xi_real_axis"] = tab
            pts = [("0", "3.5"), ("0.1", str(T // 2)), ("0.4", str(T))] if args.quick else [("0", "14.5"), ("0.1", "30.25"), ("0.4", str(T // 2)), ("0.4", str(T))]
            rep["crosscheck"] = crosscheck(ev, pts)
            rep["crosscheck_ok"] = all(c["consistent"] for c in rep["crosscheck"])
            rep["flow_zeros_table_brackets"] = flow_brackets(ev)
            allpass = allpass and rep["crosscheck_ok"]
        rep["all_pass"] = bool(allpass)
        rep["seconds_so_far"] = round(time.time() - t_all, 1)
        json.dump(rep, open(args.out, "w"), indent=1)
    rep["seconds_total"] = round(time.time() - t_all, 1)
    rep["theorem_rig_holds_at"] = [int(T) for T, b in rep["heights"].items() if b["pass"]]
    json.dump(rep, open(args.out, "w"), indent=1)
    print("all_pass", rep["all_pass"], "seconds", rep["seconds_total"], "->", args.out, flush=True)
    return 0 if rep["all_pass"] else 1


if __name__ == "__main__":
    mpr.freeze_support()
    sys.exit(main())
