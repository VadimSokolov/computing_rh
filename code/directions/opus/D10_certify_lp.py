#!/usr/bin/env python3
"""D10, Hopper proposal: repeat the check of D10_interval.py in ball arithmetic (python-flint, arb), scalar loops.

Certifies, for eta = 1, q at the height RC, for all boxes of w, xi (u = xi a(w)), qi (x_q = qi a(RC)) and all zeros:
    (t <= 0 or t A_x < D)  and for each zero: (s > A_x) or (mu <= 0) or (mu A_x < D),
with D > 0 known by hand for w < kappa_1.  Same monotone endpoint bounds as D10_interval.py; a box that fails is bisected.
Modes: real (real zeros and the Gaussian term), pair (lower pairs), far (pairs of height in [w,1] at distance >= 4, w >= 0.30).
Without python-flint it runs on a float shim for a smoke test only, and says so: that run certifies nothing.
Usage on Hopper (venv with python-flint): python D10_certify_lp.py mode wmin wmax [slice nslices]
Size: real 1.2e4 boxes x 260 cells, pair 1.2e4 x 130 x 8, far 3e3 x 30 x 8; one core, well under an hour per mode.
"""
import itertools
import math
import sys

try:
    from flint import arb, ctx
    ctx.prec = 160
    BALLS = True
except ImportError:
    BALLS = False

    class arb(float):  # shim: float semantics, for logic tests only
        def __new__(cls, x):
            return float.__new__(cls, x)

        def _w(f):
            return lambda s, *o: arb(f(float(s), *[float(v) for v in o]))

        __add__, __radd__ = _w(lambda a, b: a + b), _w(lambda a, b: b + a)
        __sub__, __rsub__ = _w(lambda a, b: a - b), _w(lambda a, b: b - a)
        __mul__, __rmul__ = _w(lambda a, b: a * b), _w(lambda a, b: b * a)
        __truediv__, __rtruediv__ = _w(lambda a, b: a / b), _w(lambda a, b: b / a)
        __pow__ = _w(lambda a, b: a ** b)
        __neg__ = _w(lambda a: -a)
        atan = _w(math.atan)
        sqrt = _w(math.sqrt)
        upper = _w(lambda a: a)
        lower = _w(lambda a: a)

        @staticmethod
        def pi():
            return arb(math.pi)


def A(x):
    return arb(x)


PI2 = arb.pi() / 2
RC = A(0.9416804581932639)
AC2 = 1 - RC * RC
AC = AC2.sqrt()
DMIN, DMID = 0.47, 4.0


def umax(*vals):
    """Exact upper end of the largest upper end."""
    us = [v.upper() for v in vals]
    best = us[0]
    for x in us[1:]:
        if x >= best:
            best = x
    return best


def lmin(*vals):
    ls = [v.lower() for v in vals]
    best = ls[0]
    for x in ls[1:]:
        if x <= best:
            best = x
    return best


def atan2p(y, x):  # y > 0
    return PI2 - (x / y).atan()


def fmob(x):
    return 2 * x / (1 + x * x)


def common(box):
    w0, w1, x0, x1, q0, q1 = [A(v) for v in box]
    a2lo, a2hi = 1 - w1 * w1, 1 - w0 * w0
    alo, ahi = a2lo.sqrt(), a2hi.sqrt()
    ral, rah = w0 / ahi, w1 / alo
    fl, fh = fmob(x0), fmob(x1)
    Ap_hi = umax(ral * fh, rah * fh).atan()
    rq = RC / AC
    Aq_lo = (rq * fmob(q0)).atan()
    lam_lo, lam_hi = w0 / RC, w1 / RC
    D = (PI2 - Ap_hi) - lam_hi * (PI2 - Aq_lo)
    xs_min = A(0.0) if (box[2] <= 0 <= box[3]) else min(x0 * x0, x1 * x1)
    Ax = 2 * w1 * (1 - xs_min) / (a2lo * (1 + xs_min) ** 2 + 4 * w0 * w0 * xs_min)
    ul, uh = lmin(x0 * alo, x0 * ahi), umax(x1 * alo, x1 * ahi)
    ql, qh = lmin(q0 * AC), umax(q1 * AC)
    return w0, w1, lam_lo, lam_hi, D, Ax, ql - uh, qh - ul


def okc(mu, Ax, D):
    return mu <= 0 or mu * Ax < D


def check_real(box, n1=230, n2=30):
    w0, w1, ll, lh, D, Ax, tl, th = common(box)
    if not (th <= 0 or th * Ax < D):
        return False
    r = (DMID / DMIN) ** (1 / n1)
    for j in range(n1):
        d0, d1 = A(DMIN * r ** j), A(DMIN * r ** (j + 1))
        s_lo = w0 / (d1 * d1 + w0 * w0)
        if s_lo > Ax:
            continue
        EL = (w1 / d0).atan() - ll * atan2p(RC, th + d1)
        ER = -(w0 / d1).atan() + lh * atan2p(RC, d0 - th)
        for E in (EL, ER):
            if not (E <= 0 or E * Ax < D * s_lo):
                return False
    for j in range(n2):
        r0, r1 = A(j / (n2 * DMID)), A((j + 1) / (n2 * DMID))
        s_lo = w0 * r0 * r0 / (1 + w0 * w0 * r0 * r0)
        if s_lo > Ax:
            continue
        trl, trh = lmin(tl * r0, tl * r1), umax(th * r0, th * r1)
        X1 = w1 * r1
        pref = 1 + w1 * w1 * r1 * r1
        DLl, DLh, DRl, DRh = 1 + trl, 1 + trh, 1 - trh, 1 - trl
        Y0, Y1 = RC * r0 / DLh, RC * r1 / DLl
        g = th / (1 + th * r0)
        pl, ph = 1 - Y1 * Y1 / 3, 1 - Y0 * Y0 / 3 + Y0 ** 4 / 5
        val = RC * RC / (3 * DLl * DLl) - w0 * w0 * (A(1) / 3 - X1 * X1 / 5)
        inner = umax(g * ph, g * pl) + umax(r0 * val, r1 * val)
        if not okc(umax(inner * pref, inner), Ax, D):
            return False
        Y0, Y1 = RC * r0 / DRh, RC * r1 / DRl
        g = th / (1 - th * r1)
        pl, ph = 1 - Y1 * Y1 / 3, 1 - Y0 * Y0 / 3 + Y0 ** 4 / 5
        val = w1 * w1 / 3 - RC * RC * (A(1) / 3 - Y1 * Y1 / 5) / (DRh * DRh)
        inner = umax(g * ph, g * pl) + umax(r0 * val, r1 * val)
        if not okc(umax(inner * pref, inner), Ax, D):
            return False
    return True


def check_pair(box, n1=115, n2=15, nth=8):
    w0, w1, ll, lh, D, Ax, tl, th = common(box)
    r = (DMID / DMIN) ** (1 / n1)
    for k in range(nth):
        t0, t1 = A(k / nth), A((k + 1) / nth)
        h0, h1 = w0 * t0, w1 * t1
        m0, m1 = w0 * (1 - t1), w1 * (1 - t0)
        p0, p1 = w0 * (1 + t0), w1 * (1 + t1)
        for j in range(n1):
            d0, d1 = A(DMIN * r ** j), A(DMIN * r ** (j + 1))
            s_lo = 2 * w0 * (d0 * d0 + w0 * w0 * (1 - t1 * t1)) / ((d1 * d1 + m1 * m1) * (d1 * d1 + p1 * p1))
            if s_lo > Ax:
                continue
            EL = (m1 / d0).atan() + (p1 / d0).atan() - ll * (atan2p(RC - h1, th + d1) + atan2p(RC + h0, th + d1))
            ER = -(m0 / d1).atan() - (p0 / d1).atan() + lh * (atan2p(RC + h1, d0 - th) + atan2p(RC - h0, d0 - th))
            for E in (EL, ER):
                if not (E <= 0 or E * Ax < D * s_lo):
                    return False
        for j in range(n2):
            r0, r1 = A(j / (n2 * DMID)), A((j + 1) / (n2 * DMID))
            s_lo = 2 * w0 * r0 * r0 / ((1 + m1 * m1 * r1 * r1) * (1 + p1 * p1 * r1 * r1))
            if s_lo > Ax:
                continue
            trl, trh = lmin(tl * r0, tl * r1), umax(th * r0, th * r1)
            DLl, DLh, DRl, DRh = 1 + trl, 1 + trh, 1 - trh, 1 - trl
            pref = 1 + 4 * w1 * w1 * r1 * r1 + w1 ** 4 * r1 ** 4
            hs0, hs1 = h0 * h0, h1 * h1
            g = th / (1 + th * r0)
            val = (RC * RC + 3 * hs1) / (3 * DLl ** 3) - (w0 * w0 + 3 * hs0) / 3 + (2 * w1) ** 4 * r1 * r1 / 5
            inner = g + umax(r0 * val, r1 * val)
            if not okc(umax(inner * pref, inner), Ax, D):
                return False
            g = th / (1 - th * r1)
            Yc = (RC + h1) ** 4 * r1 * r1 / (5 * DRl * DRl)
            val = (w1 * w1 + 3 * hs1) / 3 - ((RC * RC + 3 * hs0) / 3 - Yc) / DRh ** 3
            inner = g + umax(r0 * val, r1 * val)
            if not okc(umax(inner * pref, inner), Ax, D):
                return False
    return True


def check_far(box, n2=30, nth=8):
    w0, w1, ll, lh, D, Ax, tl, th = common(box)
    for k in range(nth):
        t0, t1 = A(k / nth), A((k + 1) / nth)
        h0, h1 = w0 * (1 - t0) + t0, w1 * (1 - t1) + t1
        hs0, hs1 = h0 * h0, h1 * h1
        for j in range(n2):
            r0, r1 = A(j / (n2 * DMID)), A((j + 1) / (n2 * DMID))
            trl, trh = lmin(tl * r0, tl * r1), umax(th * r0, th * r1)
            DLl, DLh, DRl, DRh = 1 + trl, 1 + trh, 1 - trh, 1 - trl
            pref = (1 + 2 * (w1 * w1 + hs1) * r1 * r1 + hs1 * hs1 * r1 ** 4) / (1 - hs1 * r1 * r1)
            g = th / (1 + th * r0)
            val = ((RC * RC + 3 * hs1) / 3 + A(1e-7)) / DLl ** 3 - (w0 * w0 + 3 * hs0) / 3 + (w1 + h1) ** 5 * r1 * r1 / (10 * w0)
            inner = g + umax(r0 * val, r1 * val)
            if not okc(umax(inner * pref, inner), Ax, D):
                return False
            g = th / (1 - th * r1)
            low = (RC * RC + 3 * hs0) / 3 - (RC + h1) ** 5 * r1 * r1 / (10 * RC * DRl * DRl)
            val = (w1 * w1 + 3 * hs1) / 3 + (h1 - w0) ** 5 * r1 * r1 / (10 * w0) - low / DRh ** 3
            inner = g + umax(r0 * val, r1 * val)
            if not okc(umax(inner * pref, inner), Ax, D):
                return False
    return True


def split(box):
    wid = [(box[1] - box[0]) * 4.0, box[3] - box[2], (box[5] - box[4]) * 1.5]
    k = wid.index(max(wid))
    mid = 0.5 * (box[2 * k] + box[2 * k + 1])
    a, b = list(box), list(box)
    a[2 * k + 1] = mid
    b[2 * k] = mid
    return tuple(a), tuple(b), max(wid)


def main():
    mode, wmin, wmax = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    sl, ns = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (0, 1)
    check = {"real": check_real, "pair": check_pair, "far": check_far}[mode]
    nw = max(1, int((wmax - wmin) / 0.005))
    grid = list(itertools.product(range(nw), range(20), range(8)))[sl::ns]
    n, fails = 0, 0
    for i, j, k in grid:
        stack = [(wmin + (wmax - wmin) * i / nw, wmin + (wmax - wmin) * (i + 1) / nw, -1 + j / 10, -1 + (j + 1) / 10, -1 + k / 4, -1 + (k + 1) / 4)]
        while stack:
            b = stack.pop()
            n += 1
            if check(b):
                continue
            x, y, wd = split(b)
            if wd < 2e-5:
                print("FAILED box", b, flush=True)
                fails += 1
                continue
            stack += [x, y]
    print(("ball arithmetic" if BALLS else "FLOAT SHIM, NOT CERTIFIED"), mode, f"w in [{wmin},{wmax}] slice {sl}/{ns}: boxes {n}, failed {fails}", flush=True)


if __name__ == "__main__":
    main()
