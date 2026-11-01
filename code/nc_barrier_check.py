# nc_barrier_check.py: the computed barrier check of the chapter "Computing RH: the barrier method in Thorin form" (ch/nc_barrier.tex),
# Section barrier:sec:check, Table barrier:tab:check and the numbers quoted with Proposition barrier:prop:check.
#
# What is computed. Normalisation of Section ch:dbn: xi_t(s) = 8 H_t(-2i(s - 1/2)) = 8 int_0^oo e^{t u^2} Phi_dBN(u) cosh(2 (s - 1/2) u) du,
# so a zero x + iy of H_t is the zero 1/2 - y/2 + i x/2 of xi_t. The barrier of Polymath15 (Theorem 1.2(iii) of their paper, in the
# simplified form X <= x <= X + 1, y0 <= y <= 1, 0 <= t <= t0) is, after the reflection s -> 1 - s, the rectangle
#     B = [1/2 + y0/2, 1] x [X/2, X/2 + 1/2]   in s,  for every t in [0, t0].
# Here t0 = 0.2, y0 = 0.2, and X/2 = THB = 497.4, chosen by a coarse scan of |zeta| at t = 0 (scan_window below) for a window of height 1/2
# between consecutive zeros on which |zeta| stays large, in the spirit of the choice X = 6e10 + 83952 - 0.5 of Polymath15.
# 1. Evaluator: the trapezoid rule on [0, 3] with step h for the Polya integral (the Polya route of Section ch:polya with the weight
#    e^{t u^2}), at working precision 30 + 0.341 * theta + 15 digits; h from the Poisson summation bound with a = 0.3 (as in phase_cert.py)
#    aiming at 30 correct digits relative to |xi|. Checked at t = 0 against xi(s) = s(s-1)/2 pi^{-s/2} Gamma(s/2) zeta(s) from mpmath, and
#    against a second rule with step 0.8 h, at the four corners of B.
# 2. For each t_i of a grid of [0, t0]: the winding number of xi_{t_i} around the boundary of B by an adaptive march with
#    |xi(s_{j+1}) - xi(s_j)| < |xi(s_j)| / 2 (floating point: the step criterion is tested at the grid points only).
# 3. Between grid times: at every march point s_j of time t_i, the Taylor sum sum_{k=1}^{K} D^k / k! |d_t^k xi_{t_i}(s_j)|, with
#    d_t^k xi_t(s) = 8 int u^{2k} e^{t u^2} Phi cosh(2 w u) du (the backward heat equation d_t xi_t = (1/4) d_w^2 xi_t), is required to be
#    below |xi_{t_i}(s_j)| / 2 for D = t_{i+1} - t_i; K = 8 and the last term is reported. This is Rouche's criterion (B2) of Section
#    sec:ft:persist with the Taylor form that Newton's obstruction (Proposition nw:prop:height) forces; it is computed, not a bound.
# 3a. Also one Taylor step from t = 0 with D = t0 = 0.2 (single_step_from_0).
# 3b. Tail of the Taylor series beyond order KTAY, bounded through the tilted contour u -> u + i eta (Lemma cp:lem:tilt):
#    |sum_{k > K} D^k / k! d_t^k xi_t(s)| <= 8 e^{-2 eta theta} int_0^oo e^{t (u^2 - eta^2)} |Phi(u + i eta)| cosh(2 x u)
#                                            (D |u + i eta|^2)^{K+1} e^{D |u + i eta|^2} / (K+1)! du,
#    evaluated (floating point, Gauss Legendre) at x = 1/2, theta = 497.4, eta = 0.38, t = t0 - D, for D = 0.01 and D = 0.2; the
#    same integral with the factor (D |u+i eta|^2)^{K+1} ... replaced by 1 is K_eta(t, 1/2) of Lemma cp:lem:tilt. The untilted tail
#    (eta = 0), which is the absolute value bound of Proposition nw:prop:height, is reported for comparison.
# 4. Thorin form: the minimum of v_t(sigma, theta) = Re (xi_t'/xi_t)(sigma + i theta) over a 9 x 11 grid of B at each t_i
#    (Proposition cp:prop:thorincert; v_t / pi is the Thorin density of the shifted clock at basepoint sigma when sigma is to the right of
#    every zero, Theorem cp:thm:transfer).
# 5. Control: at t = 0 the function F(s) = xi(s) (w^2 - w0^2)(w^2 - conj(w0)^2), w = s - 1/2, w0 = 1/4 + i (THB + 1/4), which has exactly
#    one zero, 3/4 + i (THB + 1/4), in B; the march must return winding 1.
# Everything is floating point (mpmath); nothing here is certified. Output: data/nc_barrier_check.json.
# Usage: python3 nc_barrier_check.py [--quick] [--workers N]
import argparse, json, math, os, sys, time
import multiprocessing as mpr
import mpmath as mp

T0, Y0 = mp.mpf("0.2"), mp.mpf("0.2")
THB = mp.mpf("497.4")
DIGITS = 30
KTAY = 8
EV = None


def phi_dbn(u):
    u = mp.mpf(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    eps = mp.mpf(10) ** (-(mp.mp.dps + 10))
    s, n = mp.mpf(0), 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        s += term
        if n >= 3 and abs(2 * mp.pi ** 2 * n ** 4 * e9 * mp.exp(-mp.pi * n * n * e4)) < eps * max(abs(s), mp.mpf(10) ** -400):
            return s
        n += 1


class Evaluator:
    """xi_t(s), xi_t'(s) and d_t^k xi_t(s) (k = 1..KTAY) by the trapezoid rule; the node data are computed once."""

    def __init__(self, thmax, hscale=1):
        mp.mp.dps = DIGITS + int(0.341 * float(thmax)) + 15
        a = mp.mpf("0.3")
        need = 2 * thmax * a + mp.pi * thmax / 4 + (DIGITS + 15) * mp.log(10) + mp.log(60)
        self.h = mp.mpf(hscale) * mp.floor(2 * mp.pi * a / need * 10 ** 6) / 10 ** 6
        self.J = int(mp.ceil(3 / self.h)) + 1
        self.u = [j * self.h for j in range(self.J + 1)]
        self.w = [8 * self.h * (mp.mpf(1) / 2 if j == 0 else 1) * phi_dbn(self.u[j]) for j in range(self.J + 1)]
        self.u2 = [x * x for x in self.u]
        self.tcache = {}
        self.scache = {}

    def tw(self, t):
        key = mp.nstr(t, 20)
        if key not in self.tcache:
            self.tcache = {key: [w * mp.exp(t * q) for w, q in zip(self.w, self.u2)]}
        return self.tcache[key]

    def hyp(self, sp):
        key = mp.nstr(sp, 30)
        if key not in self.scache:
            if len(self.scache) > 4:
                self.scache = {}
            self.scache[key] = [(mp.cosh(2 * sp * x), mp.sinh(2 * sp * x)) for x in self.u]
        return self.scache[key]

    def __call__(self, t, s, ktay=0, deriv=False):
        s = mp.mpmathify(s)
        sp, th = mp.re(s) - mp.mpf(1) / 2, mp.im(s)
        W, H = self.tw(t), self.hyp(sp)
        f = mp.mpc(0)
        d = mp.mpc(0)
        tk = [mp.mpc(0)] * (ktay + 1)
        for j in range(self.J + 1):
            c, sn = mp.cos_sin(2 * th * self.u[j])
            ch, sh = H[j]
            g = mp.mpc(ch * c, sh * sn)
            wj = W[j]
            f += wj * g
            if deriv:
                d += wj * 2 * self.u[j] * mp.mpc(sh * c, ch * sn)
            if ktay:
                q = self.u2[j]
                p = wj
                for k in range(1, ktay + 1):
                    p = p * q
                    tk[k] += p * g
        return f, d, tk[1:]


def phi_c(u):
    """Phi_dBN at complex u with Re u >= 0, summed until the explicit tail bound of phase_cert.py is below 10^-(dps+10)."""
    u = mp.mpmathify(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    x = mp.re(e4)
    eps = mp.mpf(10) ** (-(mp.mp.dps + 10))
    s, n = mp.mpc(0), 1
    while True:
        s += (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        if n >= 3:
            q1 = mp.exp(-mp.pi * (n + 1) ** 2 * x / 2)
            q2 = mp.exp(-mp.pi * (n + 1) * x)
            if abs(e9) * (2 * mp.pi ** 2 + 3 * mp.pi / x) * (4 / (mp.pi * x)) ** 2 * mp.exp(-2) * q1 / (1 - q2) < eps:
                return s
        n += 1


def tail_bound(D, eta, t, x="0.5", theta=None, K=KTAY):
    with mp.workdps(30):
        D, eta, t, x = mp.mpf(D), mp.mpf(eta), mp.mpf(t), mp.mpf(x)
        theta = THB if theta is None else mp.mpf(theta)
        X, Wt = mp.gauss_quadrature(16, "legendre")
        tot, kint = mp.mpf(0), mp.mpf(0)
        for lo, hi, cells in ((0, 1.5, 64 * 4), (1.5, 3, 8 * 4)):
            hcell = mp.mpf(hi - lo) / cells
            for i in range(cells):
                a0 = lo + i * hcell
                for xi, wi in zip(X, Wt):
                    u = a0 + hcell * (xi + 1) / 2
                    base = mp.exp(t * (u * u - eta * eta)) * abs(phi_c(mp.mpc(u, eta))) * mp.cosh(2 * x * u)
                    r2 = D * (u * u + eta * eta)
                    tot += wi * hcell / 2 * base * r2 ** (K + 1) * mp.exp(r2) / mp.factorial(K + 1)
                    kint += wi * hcell / 2 * base
        return 8 * mp.exp(-2 * eta * theta) * tot, 8 * kint


def xi_mp(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def scan_window(lo=495, hi=505, step="0.05"):
    """Coarse scan at t = 0: the start th of the window [th, th + 1/2] maximising min |zeta| over sigma in {0.6, 0.8, 1.0} and 5 heights."""
    with mp.workdps(20):
        best = None
        th = mp.mpf(lo)
        while th < hi:
            m = min(abs(mp.zeta(mp.mpc(sg, th + mp.mpf(dd)))) for sg in ("0.6", "0.8", "1.0") for dd in ("0", "0.125", "0.25", "0.375", "0.5"))
            if best is None or m > best[0]:
                best = (m, th)
            th += mp.mpf(step)
        z = [mp.zetazero(n).imag for n in range(262, 272)]
        lo_z = max(g for g in z if g < best[1]); hi_z = min(g for g in z if g > best[1] + mp.mpf(1) / 2)
        return dict(start=mp.nstr(best[1], 6), min_abs_zeta=mp.nstr(best[0], 6), zero_below=mp.nstr(lo_z, 12), zero_above=mp.nstr(hi_z, 12))


def corners():
    s0 = mp.mpf(1) / 2 + Y0 / 2
    return [mp.mpc(s0, THB), mp.mpc(1, THB), mp.mpc(1, THB + mp.mpf(1) / 2), mp.mpc(s0, THB + mp.mpf(1) / 2), mp.mpc(s0, THB)]


def march(F):
    """Counterclockwise boundary of B; returns winding, number of points, the points and values, max |dF|/|F|."""
    C = corners()
    pts, vals = [], []
    tot = mp.mpf(0)
    maxratio = mp.mpf(0)
    for A, Bc in zip(C[:-1], C[1:]):
        s, v = A, F(A)
        pts.append(s); vals.append(v)
        h = mp.mpf("0.02")
        d = (Bc - A) / abs(Bc - A)
        while abs(Bc - s) > mp.mpf(10) ** -30:
            step = min(h, abs(Bc - s))
            s2 = Bc if abs(Bc - s) <= h else s + d * step
            v2 = F(s2)
            if abs(v2 - v) > abs(v) / 2:
                h /= 2
                continue
            tot += mp.arg(v2 / v)
            maxratio = max(maxratio, abs(v2 - v) / abs(v))
            grow = abs(v2 - v) < abs(v) / 10
            s, v = s2, v2
            pts.append(s); vals.append(v)
            if grow:
                h = min(h * mp.mpf("1.5"), mp.mpf("0.1"))
    return tot / (2 * mp.pi), pts, vals, maxratio


def init(thmax):
    global EV
    EV = Evaluator(thmax)


def one_time(args):
    ti, D = args
    ti = mp.mpf(ti); D = mp.mpf(D)
    wnd, pts, vals, mr = march(lambda s: EV(ti, s)[0])
    # Taylor in t between t_i and t_i + D at the march points
    worst, lastrel, minabs = mp.mpf(0), mp.mpf(0), None
    for s, v in zip(pts, vals):
        f, _, tk = EV(ti, s, ktay=KTAY)
        terms = [D ** k / mp.factorial(k) * abs(tk[k - 1]) for k in range(1, KTAY + 1)]
        r = sum(terms) / abs(f)
        worst = max(worst, r)
        lastrel = max(lastrel, terms[-1] / abs(f))
        minabs = abs(f) if minabs is None else min(minabs, abs(f))
    # Thorin density on a grid of B
    vmin, arg = None, None
    for i in range(9):
        sg = mp.mpf(1) / 2 + Y0 / 2 + (mp.mpf(1) / 2 - Y0 / 2) * i / 8
        for j in range(11):
            th = THB + mp.mpf(j) / 20
            f, d, _ = EV(ti, mp.mpc(sg, th), deriv=True)
            v = mp.re(d / f)
            if vmin is None or v < vmin:
                vmin, arg = v, (sg, th)
    return dict(t=mp.nstr(ti, 6), winding=mp.nstr(wnd, 8), npoints=len(pts), max_step_ratio=mp.nstr(mr, 4),
                min_abs_xi=mp.nstr(minabs, 6), taylor_ratio=mp.nstr(worst, 4), taylor_last_term_ratio=mp.nstr(lastrel, 3),
                step_to_next=mp.nstr(D, 4), thorin_min=mp.nstr(vmin, 8), thorin_argmin=[mp.nstr(arg[0], 4), mp.nstr(arg[1], 7)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    t_start = time.time()
    out = dict(t0="0.2", y0="0.2", window_theta=[mp.nstr(THB, 6), mp.nstr(THB + mp.mpf(1) / 2, 6)],
               X=mp.nstr(2 * THB, 6), sigma=[mp.nstr(mp.mpf(1) / 2 + Y0 / 2, 3), "1"], digits_target=DIGITS, taylor_order=KTAY)
    out["scan"] = scan_window()
    thmax = THB + 1
    init(thmax)
    out["h"] = mp.nstr(EV.h, 8); out["nodes"] = EV.J + 1; out["dps"] = mp.mp.dps
    # checks at t = 0 against mpmath and against a second step
    ev2 = Evaluator(thmax, hscale=mp.mpf("0.8"))
    chk = []
    for c in corners()[:4]:
        f = EV(mp.mpf(0), c)[0]
        g = xi_mp(c)
        f2 = ev2(T0, c)[0]
        f1 = EV(T0, c)[0]
        chk.append(dict(s=[mp.nstr(mp.re(c), 3), mp.nstr(mp.im(c), 6)], abs_xi0=mp.nstr(abs(g), 6),
                        rel_vs_mpmath_t0=mp.nstr(abs(f - g) / abs(g), 3), rel_two_steps_t02=mp.nstr(abs(f1 - f2) / abs(f1), 3)))
    out["checks"] = chk
    # control
    w0 = mp.mpc(mp.mpf(1) / 4, THB + mp.mpf(1) / 4)
    ctrl = lambda s: EV(mp.mpf(0), s)[0] * ((s - mp.mpf(1) / 2) ** 2 - w0 ** 2) * ((s - mp.mpf(1) / 2) ** 2 - mp.conj(w0) ** 2)
    wc, pc, _, _ = march(ctrl)
    out["control"] = dict(zero=[mp.nstr(mp.re(w0) + mp.mpf(1) / 2, 3), mp.nstr(mp.im(w0), 7)], winding=mp.nstr(wc, 8), npoints=len(pc))
    # time grid: step 0.01 (quick: 0.05)
    step = mp.mpf("0.05") if a.quick else mp.mpf("0.01")
    n = int(mp.nint(T0 / step))
    grid = [(mp.nstr(i * step, 10), mp.nstr(step if i < n else 0, 10)) for i in range(n + 1)]
    with mpr.Pool(a.workers, initializer=init, initargs=(thmax,)) as pool:
        rows = pool.map(one_time, grid)
    out["rows"] = rows
    # one Taylor step from t = 0 over the whole interval [0, t0]
    r0 = one_time((mp.nstr(0, 5), mp.nstr(T0, 5)))
    out["single_step_from_0"] = dict(D=r0["step_to_next"], taylor_ratio=r0["taylor_ratio"], taylor_last_term_ratio=r0["taylor_last_term_ratio"])
    tails = {}
    for D in ("0.01", "0.2"):
        tb, Kc = tail_bound(D, "0.38", T0 - mp.mpf(D))
        tu, _ = tail_bound(D, "0", T0 - mp.mpf(D))
        tails[D] = dict(tilted_tail_bound=mp.nstr(tb, 4), untilted_tail_bound=mp.nstr(tu * 1, 4), K_eta038_x05=mp.nstr(Kc, 6),
                        ratio_to_min_abs_xi=mp.nstr(tb / min(mp.mpf(r["min_abs_xi"]) for r in rows), 4))
    out["taylor_tails"] = tails
    out["all_windings_zero"] = all(abs(float(r["winding"])) < 1e-6 for r in rows)
    out["max_taylor_ratio"] = max(float(r["taylor_ratio"]) for r in rows[:-1])
    out["seconds"] = round(time.time() - t_start, 1)
    here = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(here, "..", "data", "nc_barrier_check%s.json" % ("_quick" if a.quick else ""))
    json.dump(out, open(p, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
