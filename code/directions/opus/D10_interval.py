#!/usr/bin/env python3
"""D10: interval check of inequality (I) of D10.md, Section 3 (one pair, eta = 1, q at the height rc).

For every w in [wmin,wmax], xi in [-1,1] (u = xi a(w)), xq_ratio in [-1,1] (x_q = xq_ratio a(rc)) it verifies
  D > 0,   (x_q-u) A_x < D,   and for every zero at distance d (left or right): s > A_x or E A_x < D s,
which contradicts a rising curve at one level through the heights w and rc.  Modes: real (real zeros, multiplicity free),
pair (lower pairs h = theta w, theta in [0,1]).  Distances d < 0.47 are excluded by s > A_x (for w <= kappa_1, see the text).
D > 0 is not tested: it is proved by hand for w < kappa_1, and a box with D_lo <= 0 passes only if E <= 0 or s > A_x on it.
Mode far: pairs of height h in [w,1] with d >= 4 (use with w >= 0.30).  Mode test: enclosures against exact values at random points.
Double precision, numpy, with outward widening of about 5 ulp per elementary operation.  This is an interval check, NOT
ball arithmetic and NOT a certificate; D10_certify_lp.py repeats it on the leaves in arb.
Usage: timeout 590 python3 -I D10_interval.py mode wmin wmax [leaves.npy]
"""
import sys
import time
import numpy as np

EPS = 1.2e-15
RC = 0.9416804581932639  # any fixed height works; this double is used as exact
AC2 = 1 - RC * RC
AC_LO, AC_HI = np.sqrt(AC2) * (1 - 2 * EPS), np.sqrt(AC2) * (1 + 2 * EPS)
PI2_LO, PI2_HI = np.pi / 2 * (1 - EPS), np.pi / 2 * (1 + EPS)
DMIN, DMID = 0.47, 4.0
DEBUG = None


def dn(x):
    return x - np.abs(x) * EPS - 1e-300


def up(x):
    return x + np.abs(x) * EPS + 1e-300


def f_mob(x):
    return 2 * x / (1 + x * x)


def common(b):
    w0, w1, x0, x1, q0, q1 = b
    a2lo, a2hi = dn(1 - up(w1 * w1)), up(1 - dn(w0 * w0))
    alo, ahi = dn(np.sqrt(a2lo)), up(np.sqrt(a2hi))
    ral, rah = dn(w0 / ahi), up(w1 / alo)
    fl, fh = dn(f_mob(x0)), up(f_mob(x1))
    Ap_hi = up(np.arctan(up(np.maximum(ral * fh, rah * fh))))
    gl, gh = dn(f_mob(q0)), up(f_mob(q1))
    rql, rqh = dn(RC / AC_HI), up(RC / AC_LO)
    Aq_lo = dn(np.arctan(dn(np.minimum(rql * gl, rqh * gl))))
    lam_lo, lam_hi = dn(w0 / RC), up(w1 / RC)
    D = dn(dn(PI2_LO - Ap_hi) - up(lam_hi * up(PI2_HI - Aq_lo)))
    xs_min = np.where((x0 <= 0) & (x1 >= 0), 0.0, np.minimum(x0 * x0, x1 * x1))
    xs_min = dn(xs_min).clip(0)
    den = dn(a2lo * (1 + xs_min) ** 2 + 4 * w0 * w0 * xs_min)
    Ax = up(2 * w1 * up(1 - xs_min) / den)
    ul, uh = np.minimum(x0 * alo, x0 * ahi), np.maximum(x1 * alo, x1 * ahi)
    ql, qh = np.minimum(q0 * AC_LO, q0 * AC_HI), np.maximum(q1 * AC_LO, q1 * AC_HI)
    tl, th = dn(ql - uh), up(qh - ul)
    return w0, w1, lam_lo, lam_hi, D, Ax, tl, th


def cells(n1=230, n2=30):
    e = np.geomspace(DMIN, DMID, n1 + 1)
    r = np.linspace(0, 1 / DMID, n2 + 1)
    r[-1] = up(r[-1])
    return e[:-1], e[1:], r[:-1], r[1:]


def check_real(b):
    w0, w1, ll, lh, D, Ax, tl, th = common(b)
    ok = (th <= 0) | (up(th * Ax) < D)
    dbg = {'common': (D, Ax, tl, th)} if DEBUG is not None else None
    d0, d1, r0, r1 = cells()
    w0_, w1_, ll_, lh_, D_, Ax_, tl_, th_ = [v[:, None] for v in (w0, w1, ll, lh, D, Ax, tl, th)]
    # direct cells
    s_lo = dn(w0_ / up(d1 * d1 + w0_ * w0_))
    EL = dn(0.0 + up(np.arctan(up(w1_ / d0)))) - dn(ll_ * dn(np.arctan2(RC, up(th_ + d1))))
    ER = -dn(np.arctan(dn(w0_ / d1))) + up(lh_ * up(np.arctan2(RC, dn(d0 - th_))))
    EL, ER = up(EL), up(ER)
    okd = (s_lo > Ax_) | (((EL <= 0) | (up(EL * Ax_) < dn(D_ * s_lo))) & ((ER <= 0) | (up(ER * Ax_) < dn(D_ * s_lo))))
    ok &= okd.all(axis=1)
    if dbg is not None:
        dbg['direct'] = (d0, d1, s_lo, EL, ER)
    # far cells rho = 1/d in [r0,r1]
    trl = np.where(tl_ < 0, tl_ * r1, tl_ * r0)
    trh = np.where(th_ > 0, th_ * r1, th_ * r0)
    DLl, DLh, DRl, DRh = dn(1 + dn(trl)), up(1 + up(trh)), dn(1 - up(trh)), up(1 - dn(trl))
    s_lo = dn(w0_ * r0 * r0 / up(1 + w0_ * w0_ * r0 * r0))
    X1 = w1_ * r1
    pref = up(1 + w1_ * w1_ * r1 * r1)
    # left
    Y0, Y1 = RC * r0 / DLh, RC * r1 / DLl
    g = up(th_ / dn(1 + th_ * r0))
    pl, ph = dn(1 - Y1 * Y1 / 3), up(1 - Y0 * Y0 / 3 + Y0 ** 4 / 5)
    i1 = np.where(g >= 0, g * ph, g * pl)
    val = up(RC * RC / (3 * DLl * DLl)) - dn(w0_ * w0_ * (1 / 3 - X1 * X1 / 5))
    inner = up(i1 + np.maximum(r0 * val, r1 * val))
    mu = np.where(inner >= 0, inner * pref, inner)
    okL = (s_lo > Ax_) | (mu <= 0) | (up(mu * Ax_) < D_)
    muL = mu
    # right
    Y0, Y1 = RC * r0 / DRh, RC * r1 / DRl
    g = up(th_ / dn(1 - th_ * r1))
    pl, ph = dn(1 - Y1 * Y1 / 3), up(1 - Y0 * Y0 / 3 + Y0 ** 4 / 5)
    i1 = np.where(g >= 0, g * ph, g * pl)
    val = up(w1_ * w1_ / 3) - dn(RC * RC * (1 / 3 - Y1 * Y1 / 5) / (DRh * DRh))
    inner = up(i1 + np.maximum(r0 * val, r1 * val))
    mu = np.where(inner >= 0, inner * pref, inner)
    okR = (s_lo > Ax_) | (mu <= 0) | (up(mu * Ax_) < D_)
    ok &= (okL & okR).all(axis=1)
    if dbg is not None:
        dbg['far'] = (r0, r1, s_lo, muL, mu)
        DEBUG.update(dbg)
    return ok


def check_pair(b, nth=8):
    w0, w1, ll, lh, D, Ax, tl, th = common(b)
    ok = np.ones_like(w0, dtype=bool)
    d0, d1, r0, r1 = cells(115, 15)
    dbgs = []
    w0_, w1_, ll_, lh_, D_, Ax_, tl_, th_ = [v[:, None] for v in (w0, w1, ll, lh, D, Ax, tl, th)]
    te = np.linspace(0, 1, nth + 1)
    for t0, t1 in zip(te[:-1], te[1:]):
        h0, h1 = dn(w0_ * t0), up(w1_ * t1)
        m0, m1 = dn(w0_ * (1 - t1)).clip(0), up(w1_ * (1 - t0))  # w-h
        p0, p1 = dn(w0_ * (1 + t0)), up(w1_ * (1 + t1))          # w+h
        den = up(d1 * d1 + m1 * m1) * up(d1 * d1 + p1 * p1)
        s_lo = dn(2 * w0_ * dn(d0 * d0 + w0_ * w0_ * (1 - t1 * t1)) / up(den))
        EL = up(np.arctan(up(m1 / d0))) + up(np.arctan(up(p1 / d0)))
        EL = up(EL) - dn(ll_ * dn(dn(np.arctan2(dn(RC - h1), up(th_ + d1))) + dn(np.arctan2(dn(RC + h0), up(th_ + d1)))))
        ER = -dn(dn(np.arctan(dn(m0 / d1))) + dn(np.arctan(dn(p0 / d1))))
        ER = ER + up(lh_ * up(up(np.arctan2(up(RC - h0), dn(d0 - th_))) + up(np.arctan2(up(RC + h1), dn(d0 - th_)))))
        EL, ER = up(EL), up(ER)
        okd = (s_lo > Ax_) | (((EL <= 0) | (up(EL * Ax_) < dn(D_ * s_lo))) & ((ER <= 0) | (up(ER * Ax_) < dn(D_ * s_lo))))
        ok &= okd.all(axis=1)
        # far: mu <= R [t/(1+t rho) + rho(((rc^2+3h^2)/3)/D^3 - (w^2+3h^2)/3 + (2w)^4 rho^2/5)], and the mirror
        trl = np.where(tl_ < 0, tl_ * r1, tl_ * r0)
        trh = np.where(th_ > 0, th_ * r1, th_ * r0)
        DLl, DLh, DRl, DRh = dn(1 + dn(trl)), up(1 + up(trh)), dn(1 - up(trh)), up(1 - dn(trl))
        pref = up(1 + 4 * w1_ * w1_ * r1 * r1 + w1_ ** 4 * r1 ** 4)
        sd_lo, EL_d, ER_d = s_lo, EL, ER
        s_lo = dn(2 * w0_ * r0 * r0 / up(up(1 + m1 * m1 * r1 * r1) * up(1 + p1 * p1 * r1 * r1)))
        hs0, hs1 = dn(h0 * h0), up(h1 * h1)
        g = up(th_ / dn(1 + th_ * r0))
        val = up((RC * RC + 3 * hs1) / (3 * DLl ** 3)) - dn((w0_ * w0_ + 3 * hs0) / 3) + up((2 * w1_) ** 4 * r1 * r1 / 5)
        inner = up(g + np.maximum(r0 * val, r1 * val))
        mu = np.where(inner >= 0, inner * pref, inner)
        okL = (s_lo > Ax_) | (mu <= 0) | (up(mu * Ax_) < D_)
        muL = mu
        g = up(th_ / dn(1 - th_ * r1))
        Yc = up((RC + h1) ** 4 * r1 * r1 / (5 * DRl * DRl))
        val = up((w1_ * w1_ + 3 * hs1) / 3) - dn(((RC * RC + 3 * hs0) / 3 - Yc) / DRh ** 3)
        inner = up(g + np.maximum(r0 * val, r1 * val))
        mu = np.where(inner >= 0, inner * pref, inner)
        okR = (s_lo > Ax_) | (mu <= 0) | (up(mu * Ax_) < D_)
        ok &= (okL & okR).all(axis=1)
        if DEBUG is not None:
            dbgs.append((t0, t1, sd_lo, EL_d, ER_d, s_lo, muL, mu))
    if DEBUG is not None:
        DEBUG.update({'common': (D, Ax, tl, th), 'pair': (d0, d1, r0, r1, dbgs)})
    return ok


def check_far(b, nth=8):
    w0, w1, ll, lh, D, Ax, tl, th = common(b)
    ok = np.ones_like(w0, dtype=bool)
    _, _, r0, r1 = cells(1, 30)
    dbgs = []
    w0_, w1_, D_, Ax_, tl_, th_ = [v[:, None] for v in (w0, w1, D, Ax, tl, th)]
    te = np.linspace(0, 1, nth + 1)
    for t0, t1 in zip(te[:-1], te[1:]):
        h0, h1 = dn(w0_ * (1 - t0) + t0), up(w1_ * (1 - t1) + t1)
        hs0, hs1 = dn(h0 * h0), up(h1 * h1)
        trl = np.where(tl_ < 0, tl_ * r1, tl_ * r0)
        trh = np.where(th_ > 0, th_ * r1, th_ * r0)
        DLl, DLh, DRl, DRh = dn(1 + dn(trl)), up(1 + up(trh)), dn(1 - up(trh)), up(1 - dn(trl))
        pref = up((1 + 2 * (w1_ * w1_ + hs1) * r1 * r1 + hs1 * hs1 * r1 ** 4) / dn(1 - hs1 * r1 * r1))
        wh1 = up(w1_ + h1)
        g = up(th_ / dn(1 + th_ * r0))
        val = up((RC * RC + 3 * hs1) / 3 + 1e-7) / DLl ** 3
        val = where_pos_div(up((RC * RC + 3 * hs1) / 3 + 1e-7), DLl, DLh)
        val = up(val) - dn((w0_ * w0_ + 3 * hs0) / 3) + up(wh1 ** 5 * r1 * r1 / (10 * w0_))
        inner = up(g + np.maximum(r0 * val, r1 * val))
        muL = np.where(inner >= 0, inner * pref, inner)
        okL = (muL <= 0) | (up(muL * Ax_) < D_)
        g = up(th_ / dn(1 - th_ * r1))
        low = dn((RC * RC + 3 * hs0) / 3) - up((RC + h1) ** 5 * r1 * r1 / (10 * RC * DRl * DRl))
        val = up((w1_ * w1_ + 3 * hs1) / 3) + up((h1 - w0_) ** 5 * r1 * r1 / (10 * w0_)) - dn(low / DRh ** 3)
        inner = up(g + np.maximum(r0 * val, r1 * val))
        muR = np.where(inner >= 0, inner * pref, inner)
        okR = (muR <= 0) | (up(muR * Ax_) < D_)
        ok &= (okL & okR).all(axis=1)
        if DEBUG is not None:
            dbgs.append((t0, t1, muL, muR))
    if DEBUG is not None:
        DEBUG.update({'common': (D, Ax, tl, th), 'farp': (r0, r1, dbgs)})
    return ok


def where_pos_div(x, Dl, Dh):
    return x / Dl ** 3


def test(N=4000):
    global DEBUG
    DEBUG = {}
    rng = np.random.default_rng(5)
    viol = {k: 0 for k in ('D', 'Ax', 't', 's', 'EL', 'ER', 'sf', 'muL', 'muR', 'p_s', 'p_EL', 'p_ER', 'p_sf', 'p_muL', 'p_muR', 'f_muL', 'f_muR')}
    wid = np.array([0.005, 0.1, 0.25])
    for wlo, whi, fn, kind in ((0.02, 0.3925, check_real, 'r'), (0.02, 0.3925, check_pair, 'p'), (0.30, 0.3925, check_far, 'f')):
        lo = rng.uniform([wlo, -1, -1], [whi - wid[0], 1 - wid[1], 1 - wid[2]], (N, 3))
        b = (lo[:, 0], lo[:, 0] + wid[0], lo[:, 1], lo[:, 1] + wid[1], lo[:, 2], lo[:, 2] + wid[2])
        fn(b)
        D, Ax, tl, th = DEBUG['common']
        w = lo[:, 0] + rng.uniform(0, 1, N) * wid[0]
        xi = lo[:, 1] + rng.uniform(0, 1, N) * wid[1]
        qi = lo[:, 2] + rng.uniform(0, 1, N) * wid[2]
        a = np.sqrt(1 - w * w)
        u, xq = xi * a, qi * np.sqrt(AC2)
        lam = w / RC
        Ap, Aq = np.angle((u + 1j * w) ** 2 + 1), np.angle((xq + 1j * RC) ** 2 + 1)
        zp = u + 1j * w
        Ax_t = (2 * zp / (zp * zp + 1)).imag
        t = xq - u
        viol['D'] += np.sum(D > (np.pi / 2 - Ap) - lam * (np.pi / 2 - Aq))
        viol['Ax'] += np.sum(Ax < Ax_t)
        viol['t'] += np.sum((tl > t) | (th < t))
        I = np.arange(N)

        def kap(d, y, h):
            return np.arctan((y + h) / d) + np.arctan((y - h) / d)

        def KK(y, d):
            return y / (d * d + y * y)

        if kind == 'r':
            d0, d1, s_lo, EL, ER = DEBUG['direct']
            c = rng.integers(0, len(d0), N)
            d = d0[c] + rng.uniform(0, 1, N) * (d1[c] - d0[c])
            s = KK(w, d)
            E_L = np.arctan2(w, d) - lam * np.arctan2(RC, t + d)
            E_R = np.arctan2(w, -d) - np.pi - lam * (np.arctan2(RC, t - d) - np.pi)
            viol['s'] += np.sum(s_lo[I, c] > s)
            viol['EL'] += np.sum(EL[I, c] < E_L)
            viol['ER'] += np.sum(ER[I, c] < E_R)
            r0, r1, sf, muL, muR = DEBUG['far']
            c = rng.integers(0, len(r0), N)
            rho = r0[c] + rng.uniform(0, 1, N) * (r1[c] - r0[c])
            d = 1 / rho
            s = KK(w, d)
            E_L = np.arctan2(w, d) - lam * np.arctan2(RC, t + d)
            E_R = np.arctan2(w, -d) - np.pi - lam * (np.arctan2(RC, t - d) - np.pi)
            viol['sf'] += np.sum(sf[I, c] > s)
            viol['muL'] += np.sum(muL[I, c] < E_L / s)
            viol['muR'] += np.sum(muR[I, c] < E_R / s)
        if kind == 'p':
            d0, d1, r0, r1, dbgs = DEBUG['pair']
            for (t0, t1, sd, EL, ER, sf, muL, muR) in dbgs:
                th_ = t0 + rng.uniform(0, 1, N) * (t1 - t0)
                h = w * th_
                c = rng.integers(0, len(d0), N)
                d = d0[c] + rng.uniform(0, 1, N) * (d1[c] - d0[c])
                s = KK(w + h, d) + KK(w - h, d)
                E_L = kap(d, w, h) - lam * (np.arctan2(RC - h, t + d) + np.arctan2(RC + h, t + d))
                E_R = -kap(d, w, h) - lam * (np.arctan2(RC - h, t - d) + np.arctan2(RC + h, t - d) - 2 * np.pi)
                viol['p_s'] += np.sum(sd[I, c] > s)
                viol['p_EL'] += np.sum(EL[I, c] < E_L)
                viol['p_ER'] += np.sum(ER[I, c] < E_R)
                c = rng.integers(0, len(r0), N)
                rho = r0[c] + rng.uniform(0, 1, N) * (r1[c] - r0[c])
                d = 1 / rho
                s = KK(w + h, d) + KK(w - h, d)
                E_L = kap(d, w, h) - lam * kap(t + d, RC, h)
                E_R = -kap(d, w, h) + lam * kap(d - t, RC, h)
                viol['p_sf'] += np.sum(sf[I, c] > s)
                viol['p_muL'] += np.sum(muL[I, c] < E_L / s)
                viol['p_muR'] += np.sum(muR[I, c] < E_R / s)
        if kind == 'f':
            r0, r1, dbgs = DEBUG['farp']
            for (t0, t1, muL, muR) in dbgs:
                th_ = t0 + rng.uniform(0, 1, N) * (t1 - t0)
                h = w * (1 - th_) + th_
                c = rng.integers(0, len(r0), N)
                rho = r0[c] + rng.uniform(0, 1, N) * (r1[c] - r0[c])
                d = 1 / rho
                s = KK(w + h, d) + KK(w - h, d)
                E_L = kap(d, w, h) - lam * kap(t + d, RC, h)
                E_R = -kap(d, w, h) + lam * kap(d - t, RC, h)
                viol['f_muL'] += np.sum(muL[I, c] < E_L / s)
                viol['f_muR'] += np.sum(muR[I, c] < E_R / s)
    print("enclosure violations at random points (all must be 0):", viol)
    DEBUG = None


def main():
    if sys.argv[1] == "test":
        test()
        return
    mode, wmin, wmax = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    save = sys.argv[4] if len(sys.argv) > 4 else None
    check = {"real": check_real, "pair": check_pair, "far": check_far}[mode]
    M = {"real": 1500, "pair": 200, "far": 400}[mode]
    we = np.linspace(wmin, wmax, max(2, int((wmax - wmin) / 0.005)) + 1)
    xe = np.linspace(-1, 1, 21)
    qe = np.linspace(-1, 1, 9)
    grid = np.array([(a, b, c, d, e, f) for a, b in zip(we[:-1], we[1:]) for c, d in zip(xe[:-1], xe[1:]) for e, f in zip(qe[:-1], qe[1:])])
    stack = [grid]
    leaves, nbox, t0, smallest = [], 0, time.time(), 1.0
    weights = np.array([4.0, 1.0, 1.5])
    while stack:
        cur = stack.pop()
        if len(cur) > M:
            stack.append(cur[M:])
            cur = cur[:M]
        ok = check(tuple(cur.T))
        nbox += len(cur)
        if save:
            leaves.append(cur[ok])
        bad = cur[~ok]
        if len(bad):
            widths = (bad[:, 1::2] - bad[:, 0::2]) * weights
            k = np.argmax(widths, axis=1)
            smallest = min(smallest, widths.max(axis=1).min())
            if widths.max(axis=1).min() < 2e-5:
                i = np.argmin(widths.max(axis=1))
                print("STOP: box too small and still failing:", bad[i], flush=True)
                break
            idx = np.arange(len(bad))
            mid = 0.5 * (bad[idx, 2 * k] + bad[idx, 2 * k + 1])
            left, right = bad.copy(), bad.copy()
            left[idx, 2 * k + 1] = mid
            right[idx, 2 * k] = mid
            stack.append(np.vstack((left, right)))
        if nbox % (50 * M) < M:
            print(f"  boxes {nbox}, stack {sum(len(s) for s in stack)}, smallest weighted width {smallest:.2e}, {time.time()-t0:.0f} s", flush=True)
        if time.time() - t0 > 560:
            print("TIME LIMIT; stack", sum(len(s) for s in stack), "w range left", min(s[:, 0].min() for s in stack), max(s[:, 1].max() for s in stack), flush=True)
            break
    else:
        print(f"PASSED {mode} on w in [{wmin},{wmax}]: boxes evaluated {nbox}, smallest weighted width {smallest:.2e}, {time.time()-t0:.0f} s", flush=True)
        if save:
            np.save(save, np.vstack(leaves))
            print("leaves", sum(len(x) for x in leaves), "saved to", save)


if __name__ == "__main__":
    main()
