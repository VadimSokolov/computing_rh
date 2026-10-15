# phase_cert.py: Table cp:tab:K and the floating point march of Section cp:sec:march (A certificate for one point of the diagram) of the
# chapter on the Thorin phase diagram; the floating point part that precedes the rigorous run phase_rig.py (Theorem cp:thm:rig).
# Labels: Table cp:tab:K, Lemma cp:lem:tilt (the constants K_eta(t, x)), Lemma cp:lem:grid (the adaptive step), Theorem cp:thm:cert,
# Theorem thm:ch13:perturb (the control function), Table tab:ch21:flowzeros (the five zeros of xi_0.2 used as a check of the evaluator).
# Method. t = 0.2, xi_t(s) = 8 int_0^oo e^{t u^2} Phi(u) cosh(2 (s - 1/2) u) du, Phi = Phi_dBN (Section sec:ch21:flow), evaluated by the trapezoid
# rule on the real axis (the Polya route of polya_flint.py with the weight e^{t u^2}) in mpmath at 40 digits of output: working precision
# 40 + 0.34 T + 12 digits, step h from the Poisson summation bound 2 K_a cosh(2 a T) / (e^{2 pi a / h} - 1) with a = 0.30. The constants
#   K_eta(t, x) = 4 int_R |e^{t (u + i eta)^2}| |Phi(u + i eta)| e^{2 x u} du = 8 int_0^oo e^{t (u^2 - eta^2)} |Phi(u + i eta)| cosh(2 x u) du
# are computed in the even form, with the theta series of Phi summed at every u >= 0 until its explicit tail bound is below 10^-(dps+10)
# (Phi is even, so u < 0 is never needed); Gauss Legendre with 16 nodes on cells of width 1/42 (u < 1.5, where the terms n = 1, 2 of the series interfere) and 3/16 (u > 1.5), checked against half as many cells. The
# series of the source's part2.py, truncated at 60 terms and used at negative u where Re e^{4(u + i eta)} = e^{4u} cos(4 eta) is tiny,
# is reproduced as K_truncated_series (it gives 32.11 at eta = 0.30, x = 0.4; the correct value is 23.67), to document the error.
# The march: the boundary of Q = [0.6, 0.9] x [0, T] (T = 60) counterclockwise, with the adaptive step that keeps
# |xi(s_{j+1}) - xi(s_j)| < |xi(s_j)| / 2 (halved on failure, grown by 1.5 after a step with a change below one tenth, capped at 0.5); the
# winding is the sum of arg(xi(s_{j+1}) / xi(s_j)) over 2 pi. The zeros of xi_t on the critical line below T are located by sign changes on
# a grid of step 1/4 and bisection. Control: F(s) = xi_t(s) q(s), q(s) = (1 - w^2 / w0^2)(1 - w^2 / conj(w0)^2) / ((1 + w^2 / A^2)(1 + w^2 / B^2)),
# w = s - 1/2, w0 = 1/4 + i (A + B) / 2, which moves the pair of zeros A, B of xi_t (the 4th and 5th, 30.330 and 32.888) to the four zeros
# 1/2 +- 1/4 +- i (A + B) / 2; exactly one of them, 0.75 + 31.609 i, lies in Q_60, so the winding must be 1.
# Error control: this script is floating point (the radius of nothing is tracked); its evaluator is checked against mpmath's quad at three
# points, against the five zeros of Table tab:ch21:flowzeros (to 5 decimals), and the rigorous script phase_rig.py reproduces its values in
# ball arithmetic. Output: phase_cert.json (default ../data/phase_cert.json from this file), values as strings with 30 digits.
# Usage: python3 phase_cert.py [--quick] [--out FILE] [--workers N]; --quick runs T = 30 at 25 digits with the pair 14.077, 20.947.
import argparse
import json
import math
import multiprocessing as mpr
import os
import sys
import time

import mpmath as mp

TT = "0.2"
ETAS = ["0.30", "0.35", "0.38"]
XS = ["0.4", "0.7", "0.9", "1.0"]
FLOW = ["14.07681", "20.94669", "24.95221", "32.88801", "49.72988"]     # Table tab:ch21:flowzeros, t = 0.2
XI = None


def ns(x, n=30):
    return mp.nstr(x, n)


# the kernel Phi_dBN
def phi_dbn(u, cap=None):
    """Phi_dBN(u) = sum_n (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u}) for Re u >= 0 (u -> -u otherwise: Phi is even), summed
    until |e^{9u}| (2 pi^2 + 3 pi / x) (4 / (pi x))^2 e^{-2} q^{(n+1)^2} / (1 - q^{n+1}) < 10^-(dps+10), x = Re e^{4u}, q = e^{-pi x / 2}
    (an upper bound of the rest of the series); cap limits the number of terms (used only to reproduce the source's truncation)."""
    u = mp.mpmathify(u)
    if mp.re(u) < 0:
        u = -u
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    x = mp.re(e4)
    eps = mp.mpf(10) ** (-(mp.mp.dps + 10))
    s = mp.mpf(0)
    n = 1
    while True:
        s += (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        if cap is not None:
            if n >= cap or (n > 2 and abs(s) > 0 and abs((2 * mp.pi ** 2 * n ** 4 * e9) * mp.exp(-mp.pi * n * n * e4)) < eps):
                return s
        elif n >= 3:
            q1 = mp.exp(-mp.pi * (n + 1) ** 2 * x / 2)
            q2 = mp.exp(-mp.pi * (n + 1) * x)
            if abs(e9) * (2 * mp.pi ** 2 + 3 * mp.pi / x) * (4 / (mp.pi * x)) ** 2 * mp.exp(-2) * q1 / (1 - q2) < eps:
                return s
        n += 1


def gl_nodes(n, cells, a, b):
    X, W = mp.gauss_quadrature(n, "legendre")
    h = (b - a) / cells
    out = []
    for i in range(cells):
        lo = a + i * h
        for xi, wi in zip(X, W):
            out.append((lo + h * (xi + 1) / 2, wi * h / 2))
    return out


def K_row(args):
    """K_eta(t, x) for x in XS, by the even form, and the check by a coarser rule; also the source's truncated series for x = 0.4.
    The integrand has a fine oscillation for u < 1.5 (the terms n = 1, 2 of the series interfere, with a phase rotating at 4 pi n^2 e^{4u}
    sin 4 eta), so the cells there are narrow: 64 sc cells of [0, 1.5] and 8 sc cells of [1.5, 3] with 16 nodes (coarse: half as many)."""
    eta_s, sc, dps = args
    with mp.workdps(dps):
        t, eta = mp.mpf(TT), mp.mpf(eta_s)
        res = {}
        for ca, cb, tag in ((64 * sc, 8 * sc, "fine"), (32 * sc, 4 * sc, "coarse")):
            nodes = gl_nodes(16, ca, mp.mpf(0), mp.mpf("1.5")) + gl_nodes(16, cb, mp.mpf("1.5"), mp.mpf(3))
            A = [(u, w, mp.exp(t * (u * u - eta * eta)) * abs(phi_dbn(mp.mpc(u, eta)))) for u, w in nodes]
            res[tag] = {x: 8 * sum(w * a * mp.cosh(2 * mp.mpf(x) * u) for u, w, a in A) for x in XS}
        out = {x: ns(res["fine"][x]) for x in XS}
        out["max_abs_difference_fine_minus_coarse"] = ns(max(abs(res["fine"][x] - res["coarse"][x]) for x in XS), 3)
        # the source's route: full line [-3.2, 3.2], series truncated at 60 terms at every u, weight e^{2 x u}
        nodes = gl_nodes(24, 32, mp.mpf("-3.2"), mp.mpf("3.2"))
        f = lambda u: abs(mp.exp(t * (mp.mpc(u, eta)) ** 2)) * abs(phi_dbn_trunc(mp.mpc(u, eta))) * mp.exp(2 * mp.mpf("0.4") * u)
        out["K_truncated_series_x0.4"] = ns(4 * sum(w * f(u) for u, w in nodes), 10)
        return eta_s, out


def phi_dbn_trunc(u):
    """The series as the source's xit.py sums it, with no evenness: terms n = 1.. until the term is below 10^-(dps+5) (n > 2) or n = 61."""
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    s = mp.mpf(0)
    n = 1
    while True:
        term = (2 * mp.pi ** 2 * n ** 4 * e9 - 3 * mp.pi * n ** 2 * e5) * mp.exp(-mp.pi * n * n * e4)
        s += term
        if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5) and n > 2:
            break
        n += 1
        if n > 60:
            break
    return s


# the evaluator
class Xi:
    def __init__(self, t, T, digits):
        self.digits = digits
        mp.mp.dps = digits + int(T * math.pi / 4 / math.log(10)) + 12
        a = mp.mpf("0.3")
        need = 2 * T * a + mp.pi * T / 4 + (digits + 5) * mp.log(10) + mp.log(2 * 30)
        h = mp.floor(2 * mp.pi * a / need * 10 ** 4) / 10 ** 4
        self.h = h
        self.J = int(mp.ceil(3 / h)) + 1
        self.t = mp.mpf(t)
        self.u2, self.w8 = [], []
        for j in range(self.J + 1):
            u = j * h
            self.u2.append(2 * u)
            self.w8.append(8 * h * (mp.mpf(1) / 2 if j == 0 else 1) * mp.exp(self.t * u * u) * phi_dbn(u))
        self.cache = {}
        self.cs = getattr(mp, "cos_sin", None)

    def tab(self, sp):
        key = ns(sp, 25)
        tb = self.cache.get(key)
        if tb is None:
            tb = [(w * mp.cosh(sp * u2), w * mp.sinh(sp * u2), u2) for w, u2 in zip(self.w8, self.u2)]
            if len(self.cache) < 8:
                self.cache[key] = tb
        return tb

    def __call__(self, s):
        s = mp.mpmathify(s)
        sp, th = mp.re(s) - mp.mpf(1) / 2, mp.im(s)
        re = im = mp.mpf(0)
        for A, B, u2 in self.tab(sp):
            if self.cs is not None:
                c, sn = self.cs(th * u2)
            else:
                c, sn = mp.cos(th * u2), mp.sin(th * u2)
            re += A * c
            im += B * sn
        return mp.mpc(re, im)


def load_zeta_zeros():
    """The ordinates of the zeros of zeta from data/zeros_1700.txt (70 digits), when the file is found; otherwise None."""
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ("zeros_1700.txt", os.path.join(here, "..", "data", "zeros_1700.txt"), os.path.join(here, "zeros_1700.txt")):
        if os.path.exists(p):
            return [mp.mpf(l.strip()) for l in open(p) if l.strip()]
    return None


def init_worker(t, T, digits):
    global XI
    XI = Xi(t, T, digits)


# the march
def march_rect(F, T, eps="0.1", sig1="0.9"):
    """Counterclockwise boundary of [1/2 + eps, sig1] x [0, T]; adaptive step with |F(s2) - F(s)| < |F(s)| / 2; returns the sum of arg."""
    eps, sig1, T = mp.mpf(eps), mp.mpf(sig1), mp.mpf(T)
    corners = [mp.mpc(mp.mpf(1) / 2 + eps, 0), mp.mpc(sig1, 0), mp.mpc(sig1, T), mp.mpc(mp.mpf(1) / 2 + eps, T), mp.mpc(mp.mpf(1) / 2 + eps, 0)]
    nev, maxratio, edges = 0, mp.mpf(0), []
    for A, B in zip(corners[:-1], corners[1:]):
        s, v = A, F(A)
        nev += 1
        h = mp.mpf("0.1")
        d = (B - A) / abs(B - A)
        tot = mp.mpf(0)
        while abs(B - s) > 0:
            step = min(h, abs(B - s))
            s2 = B if abs(B - s) <= h else s + d * step
            v2 = F(s2)
            nev += 1
            if abs(v2 - v) > abs(v) / 2:
                h /= 2
                continue
            tot += mp.arg(v2 / v)
            maxratio = max(maxratio, abs(v2 - v) / abs(v))
            grow = abs(v2 - v) < abs(v) / 10
            s, v = s2, v2
            if grow:
                h = min(h * mp.mpf("1.5"), mp.mpf("0.5"))
        edges.append(tot)
    return dict(winding=sum(edges) / (2 * mp.pi), edges=edges, nev=nev, maxratio=maxratio)


def scan_chunk(args):
    lo, hi = args
    out = []
    step = mp.mpf(1) / 4
    prev = None
    k = int(lo * 4)
    while k * step <= hi:
        th = k * step
        val = mp.re(XI(mp.mpc(mp.mpf(1) / 2, th)))
        if prev is not None and prev[1] * val < 0:
            out.append((prev[0], th))
        prev = (th, val)
        k += 1
    return [(ns(a, 20), ns(b, 20)) for a, b in out]


def refine(br):
    a, b = mp.mpf(br[0]), mp.mpf(br[1])
    f = lambda th: mp.re(XI(mp.mpc(mp.mpf(1) / 2, th)))
    fa = f(a)
    for _ in range(46):
        m = (a + b) / 2
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return ns((a + b) / 2, 25)


def run_march(args):
    kind, T, pair = args
    zs = [mp.mpf(z) for z in pair]
    if kind == "plain":
        return march_rect(lambda s: XI(s), T)
    A, B = zs
    g = (A + B) / 2
    w0 = mp.mpc(mp.mpf(1) / 4, g)

    def q(s):
        w = s - mp.mpf(1) / 2
        return (1 - w * w / w0 ** 2) * (1 - w * w / mp.conj(w0) ** 2) / ((1 + w * w / A ** 2) * (1 + w * w / B ** 2))

    return march_rect(lambda s: XI(s) * q(s), T)


# main
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="reduced case: T = 30, 25 digits, the control pair 14.077, 20.947")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "phase_cert.json"))
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()
    T = 30 if args.quick else 60
    digits = 25 if args.quick else 40
    t_all = time.time()
    rep = dict(script="phase_cert.py", quick=bool(args.quick), t=TT, T=T, digits=digits, workers=args.workers)
    pool = None
    if args.workers > 1:
        pool = mpr.Pool(args.workers, initializer=init_worker, initargs=(TT, T, digits))
    init_worker(TT, T, digits)
    mapf = pool.map if pool is not None else (lambda f, xs: [f(x) for x in xs])
    rep.update(dps_working=mp.mp.dps, trapezoid_h=ns(XI.h, 6), trapezoid_nodes=XI.J + 1)
    # 1. Table cp:tab:K, left part: K_eta(0.2, x)
    t0 = time.time()
    rows = dict(mapf(K_row, [(e, 1 if args.quick else 2, 20 if args.quick else 30) for e in ETAS]))
    rep["K_eta"] = rows
    rep["K_seconds"] = round(time.time() - t0, 1)
    print("K", {e: {x: rows[e][x][:7] for x in XS} for e in ETAS}, rep["K_seconds"], flush=True)
    # 2. right part: xi_0.2 on the real axis, and the evaluator against mpmath's quad
    rep["xi_real"] = {sg: ns(mp.re(XI(mp.mpf(sg))), 30) for sg in ("0.5", "0.6", "0.9", "1.0", "1.5")}
    chk = []
    with mp.workdps(40):
        for s in ("0.5", mp.mpc("0.6", "14"), mp.mpc("0.9", "20")):
            s = mp.mpmathify(s)
            f = lambda u: mp.exp(mp.mpf(TT) * u * u) * phi_dbn(u) * mp.cosh(2 * (s - mp.mpf(1) / 2) * u)
            ref = 8 * mp.quad(f, mp.linspace(0, 3, 13))
            mine = XI(s)
            chk.append(dict(s=str(s), quad=ns(ref, 15), trapezoid=ns(mine, 15), abs_difference=ns(abs(ref - mine), 3)))
    rep["evaluator_vs_quad"] = chk
    # 3. zeros of xi_t on the critical line below T
    t0 = time.time()
    nch = max(1, args.workers * 2)
    cuts = [mp.mpf(int(T * 4 * i / nch)) / 4 for i in range(nch + 1)]     # multiples of the grid step 1/4, so that no bracket is lost
    brs = sum(mapf(scan_chunk, [(max(cuts[i], mp.mpf(1) / 4), cuts[i + 1]) for i in range(nch)]), [])
    zeros = sorted(mp.mpf(z) for z in mapf(refine, brs))
    rep["zeros_on_line"] = [ns(z, 25) for z in zeros]
    rep["number_of_zeros_on_line"] = len(zeros)
    rep["flow_zeros_check"] = [dict(printed=f, nearest=ns(min(zeros, key=lambda z: abs(z - mp.mpf(f))), 10),
                                    abs_difference=ns(min(abs(z - mp.mpf(f)) for z in zeros), 3)) for f in FLOW if mp.mpf(f) < T]
    rep["zeros_seconds"] = round(time.time() - t0, 1)
    gam = load_zeta_zeros()
    if gam is not None:
        below = [g for g in gam if g < T]
        rep["zeta_zeros_below_T"] = len(below)
        if len(below) == len(zeros):
            sh = [z - g for z, g in zip(zeros, below)]
            rep["shift_of_the_zeros_by_the_flow"] = dict(min_abs=ns(min(abs(x) for x in sh), 6), max_abs=ns(max(abs(x) for x in sh), 6),
                                                         mean=ns(sum(sh) / len(sh), 6), all_negative=bool(all(x < 0 for x in sh)))
    print("zeros", len(zeros), [ns(z, 7) for z in zeros], flush=True)
    # 4. the marches
    pair_idx = (0, 1) if args.quick else (3, 4)
    pair = [ns(zeros[pair_idx[0]], 25), ns(zeros[pair_idx[1]], 25)]
    t0 = time.time()
    res = mapf(run_march, [("plain", T, pair), ("control", T, pair)])
    for tag, r in zip(("march_xi", "march_control"), res):
        rep[tag] = dict(winding=ns(r["winding"], 12), evaluations=r["nev"], max_relative_step_change=ns(r["maxratio"], 5),
                        arg_bottom=ns(r["edges"][0], 12), arg_right=ns(r["edges"][1], 12), arg_top=ns(r["edges"][2], 12),
                        arg_left=ns(r["edges"][3], 12), arg_left_up_over_pi=ns(-r["edges"][3] / mp.pi, 12))
    A, B = mp.mpf(pair[0]), mp.mpf(pair[1])
    rep["control"] = dict(pair=pair, g=ns((A + B) / 2, 15), new_zero_in_Q="0.75 + %s i" % ns((A + B) / 2, 8),
                          winding_expected=1, winding_computed=rep["march_control"]["winding"])
    rep["march_seconds"] = round(time.time() - t0, 1)
    print("winding xi", rep["march_xi"]["winding"], rep["march_xi"]["evaluations"], "control", rep["march_control"]["winding"],
          rep["march_control"]["evaluations"], flush=True)
    rep["checks"] = {
        "winding of xi_t is 0 to the working precision": bool(abs(mp.mpf(rep["march_xi"]["winding"])) < mp.mpf(10) ** (-digits + 5)),
        "control winding is 1 to the working precision": bool(abs(mp.mpf(rep["march_control"]["winding"]) - 1) < mp.mpf(10) ** (-digits + 5)),
        "evaluator agrees with quad": all(mp.mpf(c["abs_difference"]) < mp.mpf(10) ** (-digits + 5) for c in chk),
        "zeros of Table tab:ch21:flowzeros reproduced to 5 decimals": all(mp.mpf(c["abs_difference"]) < mp.mpf("1e-5") for c in rep["flow_zeros_check"]),
        "K at eta = 0.30 increasing in x": bool(mp.mpf(rows["0.30"]["0.4"]) < mp.mpf(rows["0.30"]["0.7"]) < mp.mpf(rows["0.30"]["0.9"]))}
    rep["all_pass"] = all(rep["checks"].values())
    rep["seconds_total"] = round(time.time() - t_all, 1)
    if pool is not None:
        pool.close()
    json.dump(rep, open(args.out, "w"), indent=1)
    print("all_pass", rep["all_pass"], "seconds", rep["seconds_total"], "->", args.out, flush=True)
    return 0 if rep["all_pass"] else 1


if __name__ == "__main__":
    mpr.freeze_support()
    sys.exit(main())
