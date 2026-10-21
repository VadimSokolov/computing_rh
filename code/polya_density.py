# polya_density.py: the density f_X of X = Z sqrt(2Y) (Y the reciprocal Xi law, Z standard normal) and the numerics of Sections pol:sec:pf and "Problem pol:op:pf2: f_X is log concave" (Theorem pol:thm:pf2, Lemma pol:lem:residue, Theorem pol:thm:approxpf, Proposition pol:prop:hankel) in Chapter 16 (Polya and the Zeros of Entire Functions). Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were not read, run or copied).
# Notation. phi(t^2) = xi(1/2)/xi(1/2 + t) is the characteristic function of X, f_X^(j)(x) = (1/pi) int_0^inf phi(t^2) d^j/dx^j cos(tx) dt, N = f f'' - f'^2, gamma_k are the ordinates of the zeros (data/zeros_1700.txt, 68 digits), c_k = -Xi(0)/Xi'(gamma_k) the residues (Xi'(gamma) = i xi'(1/2 + i gamma) = i [s(s-1)/2 pi^(-s/2) Gamma(s/2)] zeta'(s) at s = 1/2 + i gamma), m_j = (1/pi) int_0^inf t^j phi(t^2) dt, and C_j(T) = (1/2 pi) int |u - iT|^j xi(1/2)/|xi(1/2 + u - iT)| du.
# Checks. Total mass, variance (against 2 a_1 = xi''(1/2)/xi(1/2)), second differences of log f_X on [0, 1.6] at step 0.02, the normalised 3 by 3 determinants of random ordered triples in [-0.8, 0.8] and the tail slope on [1.2, 1.6]. The residues c_1 = 359.524, c_2 = -28007.2, c_1 c_2 (gamma_2 - gamma_1)^2 = -4.776e8, the constants C_0, C_1, C_2 = 5.99e4, 1.38e6, 3.20e7 at T = 23 and the same constants at T = (gamma_12 + gamma_13)/2, a numerical test of the bound |E_T^(j)(x)| <= C_j(T) e^(-Tx) against direct quadrature, m_0..m_4 = 1.899, 7.321, 44.94, 355.8, 3361. Theorem pol:thm:pf2: part (A) on [0, 0.95] (the chapter has 1087 cells), part (B) on [0.95, 2.1] with twelve zeros (26 cells, sum -2.18e8 and remainder 7.0e6 at x = 0.95), part (C) from 2.1 with T = 23 (4.687e8 < 4.776e8), and -(log f_X)'' at x = 0.95 (about 2.9) and at x = 4 (4e-9), with the series of the figure on 400 points of [0, 4]. Theorem pol:thm:approxpf: S <= 1.533e-12, eps_0, eps_1, eps_2 for K = 300 (1.57e-10, 1.3e-9, 1.3e-8), f_0(0) <= 1.976, the minor bounds -1.1e-8 and -1.2e-7 recomputed from eps_0, and the scan of 6050 Toeplitz type minors for n = 3 and n = 4 (smallest 1.6e-10 and 3.5e-19). The Hankel quantities -H_3/f^3 and H_4/f^4 on [0, 4] at step 0.02 (smallest 1.05e-21 and 2.08e-39 at x = 4) and the Cauchy Binet formula of Proposition pol:prop:hankel with its sign claim on random exponential sums.
# Method. f_X^(j) is computed from composite Gauss Legendre panels of width 10 (96 nodes each) on t in [0, 170] (phi(t^2) < 1e-89 beyond) at 80 digits, the values phi(t_i^2) computed once from the closed form of xi with mpmath; the same nodes serve every x, for |x| <= 1.6 by the chapter's rule and, as a check, up to 4. Beyond 1.6 the residue expansion f^(j)(x) = sum c_k (-gamma_k)^j e^(-gamma_k x) is used with 14 zeros (the chapter) and 40 zeros (check). C_j(T) is (1/pi) int_0^inf of an even integrand, by Gauss Legendre panels of width 1 with 24 nodes until the panel contribution is below 1e-25 of the total. Certificate (A) is implemented as the text describes it, with cells chosen greedily from the left (the largest h with h B(h) <= theta |N(x_i)|, B the bound on sup |N'| built from |f^(j)(x_i)| + h m_(j+1); theta = 0.99, 0.9, 0.5 are all reported because the chapter does not say how its cells were chosen). Certificate (B) writes N e^((gamma_1 + gamma_2)x) = S(x) + remainder with S = sum over k < l <= 12 of c_k c_l (gamma_l - gamma_k)^2 e^(-(gamma_k + gamma_l - gamma_1 - gamma_2) x), bounds the remainder by e^((gamma_1 + gamma_2) x) [e^(-Tx) (P C_2 + P'' C_0 + 2 P' C_1) + e^(-2Tx) (C_0 C_2 + C_1^2)] with P^(j) replaced by its majorant sum |c_k| gamma_k^j e^(-gamma_k x), and on a cell [a, b] bounds S above termwise (positive terms at a, negative terms at b) and the remainder by its value at a. Certificate (C) is the same remainder bound with two terms and T = 23. The approximation constants are the integrals of the theorem by Gauss Legendre panels of width 4. Hankel determinants and minors are determinants in mpmath at 80 digits.
# Error control. (a) Quadrature against the residue expansion with 12, 14 and 40 zeros on [0.95, 1.6] for j = 0..3 (relative differences recorded). (b) The bound of Lemma pol:lem:residue tested at 13 abscissae for T = 23 and T = (gamma_12 + gamma_13)/2 (ratio of the actual error to C_j e^(-Tx)). (c) N e^((gamma_1 + gamma_2)x) from quadrature against S(x) and its remainder bound, and against c_1 c_2 (gamma_2 - gamma_1)^2 with the bound of (C). (d) The total mass must be 1 and the variance 2 a_1, with the integral of f_X over [0, 1.6] done by 96 point Gauss Legendre in x on the quadrature values and the tail from the residues. (e) The scan and the Hankel grid are repeated with 40 residues instead of 14 beyond 1.6, and the Hankel quantities from the residues are compared with those from quadrature at x = 1.62 to 2.0. (f) The Cauchy Binet formula is checked on a random exponential sum at 60 digits and its sign claim at 200 digits for n up to 8. Nothing here is interval arithmetic: the certificate parts are floating point at 30 to 80 digits with explicit analytic bounds, as in the chapter.
# Output: data/polya_density.json (values as strings: constants, the three certificate parts, the scan, the Hankel grid, the series of the figure) and, beside it, the log data/polya_density.log.
# Usage: python3 polya_density.py [--quick] [--out FILE] [--zeros FILE] [--workers W]; --quick uses 40 digits, 11 panels of 48 nodes, 16 residues, three values of d in the scan and 40 points of the figure series (about 6 s in one process, writes to /tmp unless --out is given); the full case is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -t 120 -g polya_density.json -g polya_density.log code/polya_density.py data/zeros_1700.txt -- --out polya_density.json
import argparse
import json
import math
import pickle
import random
import sys
import tempfile
import time
from itertools import combinations
from multiprocessing import get_context
from pathlib import Path

from mpmath import mp, mpf, mpc
from mpmath.calculus.quadrature import GaussLegendre

T0 = time.time()
LOG = []
G = {}
WORKERS = 1
_GL = {}


def log(*a):
    s = '[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(s, flush=True)
    LOG.append(s)


def S(x, n=None):
    return mp.nstr(x, n if n else mp.dps + 5)


def fs(x, n=12):
    return mp.nstr(x, n)


def X(s):
    if '/' in s:
        p, q = s.split('/')
        return mpf(p) / mpf(q)
    return mpf(s)


def find_zeros(arg):
    here = Path(__file__).resolve().parent
    for c in ([arg] if arg else []) + [here / 'zeros_1700.txt', Path('data/zeros_1700.txt'), here.parent / 'data' / 'zeros_1700.txt']:
        if c and Path(c).is_file():
            return Path(c)
    raise SystemExit('zeros_1700.txt not found; pass --zeros')


def gl_nodes(deg):
    key = (deg, mp.prec)
    if key not in _GL:
        _GL[key] = GaussLegendre(mp).calc_nodes(deg, mp.prec)
    return _GL[key]


def xi_s(s):
    if s == 1:
        return mpf(1) / 2
    return s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def XI0():
    return xi_s(mpf(1) / 2)


def pmap(fn, items):
    """Map over small picklable tasks: a fork pool (workers > 1) or the same code in this process with a pickle round trip."""
    items = list(items)
    if WORKERS > 1 and len(items) > 1:
        with get_context('fork').Pool(min(WORKERS, len(items))) as pool:
            return pool.map(fn, items, chunksize=1)
    saved = mp.dps
    out = [pickle.loads(pickle.dumps(fn(pickle.loads(pickle.dumps(it))))) for it in items]
    mp.dps = saved
    return out


def panel_points(npan, wd, deg):
    nd = gl_nodes(deg)
    h = mpf(wd) / 2
    ts, ws = [], []
    for p in range(npan):
        c = (2 * p + 1) * h
        for x, w in nd:
            ts.append(c + h * x)
            ws.append(h * w)
    return ts, ws


def _t_phi(arg):
    dps, tstr = arg
    mp.dps = dps
    x0 = XI0()
    half = mpf(1) / 2
    return [S(x0 / xi_s(half + mpf(t))) for t in tstr]


def build_table(dps, npan, wd, deg, jmax):
    """Nodes t_i and the weights W_i t_i^j with W_i = w_i phi(t_i^2)/pi, so that f^(j)(x) = sum W_i t_i^j (d^j/dx^j cos)(t_i x)/t_i^j."""
    mp.dps = dps
    ts, ws = panel_points(npan, wd, deg)
    nch = WORKERS * 3 if WORKERS > 1 else 1
    idx = [list(range(i, len(ts), nch)) for i in range(nch)]
    tstr = [S(t, dps + 10) for t in ts]
    res = pmap(_t_phi, [(dps, [tstr[k] for k in ii]) for ii in idx])
    mp.dps = dps
    ph = [None] * len(ts)
    for ii, r in zip(idx, res):
        for k, v in zip(ii, r):
            ph[k] = mpf(v)
    W = [w * p / mp.pi for w, p in zip(ws, ph)]
    Wt = [W]
    for j in range(jmax):
        Wt.append([a * t for a, t in zip(Wt[-1], ts)])
    return dict(t=ts, Wt=Wt, n=len(ts), tmax=npan * wd)


def fder(Q, x, jmax):
    cs = [mp.cos_sin(t * x) for t in Q['t']]
    cc = [a for a, b in cs]
    ss = [b for a, b in cs]
    out = []
    for j in range(jmax + 1):
        r = j % 4
        v = mp.fdot(Q['Wt'][j], cc if r % 2 == 0 else ss)
        out.append(v if r in (0, 3) else -v)
    return out


def f0(Q, x):
    return mp.fdot(Q['Wt'][0], [mp.cos(t * x) for t in Q['t']])


def rder(c, g, K, x, jmax):
    e = [mp.exp(-g[k] * x) for k in range(K)]
    return [mp.fsum(c[k] * (-g[k]) ** j * e[k] for k in range(K)) for j in range(jmax + 1)]


def _t_res(arg):
    dps, gs = arg
    mp.dps = dps
    x0 = XI0()
    out = []
    for gstr in gs:
        s = mpc(mpf(1) / 2, mpf(gstr))
        A = s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)
        Xip = mpc(0, 1) * A * mp.zeta(s, 1, 1)
        out.append((S(-x0 / Xip.real), mp.nstr(abs(Xip.imag) / abs(Xip.real), 5)))
    return out


def residues(zs, K, dps, nch):
    idx = [list(range(i, K, nch)) for i in range(nch)]
    res = pmap(_t_res, [(dps, [zs[k] for k in ii]) for ii in idx])
    mp.dps = dps
    c = [None] * K
    im = mpf(0)
    for ii, r in zip(idx, res):
        for k, (cs, ims) in zip(ii, r):
            c[k] = mpf(cs)
            im = max(im, mpf(ims))
    return c, im


def _t_grid(arg):
    which, dps, xs, jmax = arg
    mp.dps = dps
    Q = G[which]
    return [[S(v) for v in fder(Q, X(s), jmax)] for s in xs]


def quad_derivs(which, xs, jmax, per=3):
    dps = G['DH'] if which == 'H' else G['DL']
    chunks = [xs[i:i + per] for i in range(0, len(xs), per)]
    res = pmap(_t_grid, [(which, dps, ch, jmax) for ch in chunks])
    mp.dps = dps
    return [[mpf(v) for v in row] for r in res for row in r]


def _t_C(arg):
    Ts, dps, p0, p1, wd, deg = arg
    mp.dps = dps
    T = mpf(Ts)
    x0 = XI0()
    nd = gl_nodes(deg)
    h = mpf(wd) / 2
    half = mpf(1) / 2
    I = [mpf(0)] * 3
    for p in range(p0, p1):
        c = (2 * p + 1) * h
        for x, w in nd:
            u = c + h * x
            ab = abs(xi_s(mpc(half + u, T)))
            base = x0 / ab * h * w
            r2 = u * u + T * T
            I[0] += base
            I[1] += base * mp.sqrt(r2)
            I[2] += base * r2
    return [S(v) for v in I]


def Cconst(T, P):
    dps = P['DC']
    mp.dps = dps
    ntask = max(WORKERS, 2)
    cs = 8 if WORKERS > 1 else 5
    total = [mpf(0)] * 3
    p = 0
    while True:
        res = pmap(_t_C, [(S(T), dps, p + i * cs, p + (i + 1) * cs, P['Cwd'], P['Cdeg']) for i in range(ntask)])
        mp.dps = dps
        last = None
        for r in res:
            last = [mpf(v) for v in r]
            for j in range(3):
                total[j] += last[j]
        p += ntask * cs
        if last[0] <= P['Ctol'] * total[0]:
            break
        if p > 3000:
            raise RuntimeError('C_j integral did not converge')
    return [t / mp.pi for t in total], p * P['Cwd']


def _t_certA(arg):
    theta_s, keep = arg
    dps = G['DL']
    mp.dps = dps
    theta = mpf(theta_s)
    Q = G['L']
    m = [mp.fsum(Q['Wt'][j]) for j in range(5)]
    xend = mpf('0.95')
    x = mpf(0)
    cells = []
    hs = []
    margins = []
    N0 = None
    while x < xend:
        d = fder(Q, x, 3)
        N = d[0] * d[2] - d[1] ** 2
        if N0 is None:
            N0 = N
        if N >= 0:
            return dict(ok=False, x=S(x, 10), N=S(N, 6))
        a0, a1, a2, a3 = [abs(v) for v in d]
        aN = -N

        def B(h):
            return (a0 + h * m[1]) * (a3 + h * m[4]) + (a1 + h * m[2]) * (a2 + h * m[3])
        hm = xend - x
        if hm * B(hm) <= theta * aN:
            h = hm
        else:
            lo, hi = mpf(0), hm
            for _ in range(55):
                mid = (lo + hi) / 2
                if mid * B(mid) <= theta * aN:
                    lo = mid
                else:
                    hi = mid
            h = lo
        if h < mpf(10) ** (-14):
            return dict(ok=False, x=S(x, 10), h=S(h, 6))
        margins.append((aN - h * B(h)) / aN)
        hs.append((h, x))
        if keep:
            cells.append((S(x, 10), S(h, 8), S(N, 6)))
        x += h
    hmin = min(hs)
    hmax = max(hs)
    return dict(ok=True, theta=theta_s, ncells=len(hs), h_min=S(hmin[0], 6), h_min_at=S(hmin[1], 8), h_max=S(hmax[0], 6), h_max_at=S(hmax[1], 8),
                min_relative_margin=S(min(margins), 6), N_at_0=S(N0, 8), N_last_start=S(N, 6), x_last_start=S(x, 8), m=[S(v, 12) for v in m],
                cells=cells)


def cert_B(c, g, C, TB, dps=50):
    mp.dps = dps
    K = 12
    g1, g2 = g[0], g[1]
    C0, C1, C2 = C
    terms = [(c[k] * c[l] * (g[l] - g[k]) ** 2, g[k] + g[l] - g1 - g2) for k in range(K) for l in range(k + 1, K)]

    def Sval(x):
        return mp.fsum(A * mp.exp(-lam * x) for A, lam in terms)

    def Sup(a, b):
        return mp.fsum((A * mp.exp(-lam * a) if A > 0 else A * mp.exp(-lam * b)) for A, lam in terms)

    def Rm(x):
        e = [mp.exp(-g[k] * x) for k in range(K)]
        P0 = mp.fsum(abs(c[k]) * e[k] for k in range(K))
        P1 = mp.fsum(abs(c[k]) * g[k] * e[k] for k in range(K))
        P2 = mp.fsum(abs(c[k]) * g[k] ** 2 * e[k] for k in range(K))
        return mp.exp((g1 + g2 - TB) * x) * (P0 * C2 + P2 * C0 + 2 * P1 * C1) + mp.exp((g1 + g2 - 2 * TB) * x) * (C0 * C2 + C1 ** 2)

    def cell_ok(a, b):
        return Sup(a, b) + Rm(a) < 0
    xa, xb = mpf('0.95'), mpf('2.1')
    cells = []
    a = xa
    ok = True
    while a < xb:
        if cell_ok(a, xb):
            b = xb
        else:
            if not cell_ok(a, a):
                ok = False
                break
            lo, hi = a, xb
            for _ in range(80):
                mid = (lo + hi) / 2
                if cell_ok(a, mid):
                    lo = mid
                else:
                    hi = mid
            b = lo
        if b <= a or len(cells) > 100000:
            ok = False
            break
        cells.append((a, b))
        a = b
    M = 1
    while M <= 2000:
        if all(cell_ok(xa + i * (xb - xa) / M, xa + (i + 1) * (xb - xa) / M) for i in range(M)):
            break
        M += 1
    else:
        M = None
    pass26 = all(cell_ok(xa + i * (xb - xa) / 26, xa + (i + 1) * (xb - xa) / 26) for i in range(26))
    grid = [xa + (xb - xa) * i / 1000 for i in range(1001)]
    worst = max((Sval(x) + Rm(x), x) for x in grid)
    A12 = c[0] * c[1] * (g2 - g1) ** 2
    return dict(ok=ok, greedy_cells=len(cells), uniform_cells_needed=M, chapter_26_uniform_cells_pass=pass26, S_at_095=Sval(xa), Rm_at_095=Rm(xa),
                S_at_21=Sval(xb), Rm_at_21=Rm(xb), A12=A12, worst_S_plus_R=worst[0], worst_at=worst[1], greedy_first_cells=[(S(p, 6), S(q, 6)) for p, q in cells[:5]],
                Sval=Sval, Rm=Rm)


def cert_C(c, g, C, T=23, dps=50):
    mp.dps = dps
    g1, g2 = g[0], g[1]
    C0, C1, C2 = C
    A12 = c[0] * c[1] * (g2 - g1) ** 2

    def maj(x):
        e = [mp.exp(-g[k] * x) for k in range(2)]
        P0 = mp.fsum(abs(c[k]) * e[k] for k in range(2))
        P1 = mp.fsum(abs(c[k]) * g[k] * e[k] for k in range(2))
        P2 = mp.fsum(abs(c[k]) * g[k] ** 2 * e[k] for k in range(2))
        return mp.exp((g1 + g2 - T) * x) * (P0 * C2 + P2 * C0 + 2 * P1 * C1) + mp.exp((g1 + g2 - 2 * T) * x) * (C0 * C2 + C1 ** 2)
    lo, hi = mpf('0.1'), mpf(10)
    for _ in range(100):
        mid = (lo + hi) / 2
        if maj(mid) > abs(A12):
            lo = mid
        else:
            hi = mid
    return dict(A12=A12, maj_at_21=maj(mpf('2.1')), x_where_maj_equals_A12=hi, maj=maj)


def hankel_pair(D):
    H3 = mp.matrix(3, 3)
    H4 = mp.matrix(4, 4)
    for i in range(3):
        for j in range(3):
            H3[i, j] = D[i + j]
    for i in range(4):
        for j in range(4):
            H4[i, j] = D[i + j]
    return -mp.det(H3) / D[0] ** 3, mp.det(H4) / D[0] ** 4


def _t_rand(arg):
    dps, seed, ntri, lo = arg
    mp.dps = dps
    Q = G['L']
    rng = random.Random(seed)
    best = None
    nbad = 0
    for _ in range(ntri):
        xs = sorted(rng.uniform(-lo, lo) for _ in range(3))
        ys = sorted(rng.uniform(-lo, lo) for _ in range(3))
        A = mp.matrix(3, 3)
        for i in range(3):
            for j in range(3):
                A[i, j] = f0(Q, abs(mpf(xs[i]) - mpf(ys[j])))
        nm = mp.det(A) / (A[0, 0] * A[1, 1] * A[2, 2])
        if nm <= 0:
            nbad += 1
        if best is None or nm < best[0]:
            best = (nm, xs, ys)
    return (S(best[0], 10), nbad, xs_str(best[1]), xs_str(best[2]))


def xs_str(v):
    return [repr(t) for t in v]


def _t_scan(arg):
    d, Fstr, dps, nlist, emax, base = arg
    mp.dps = dps
    F = [mpf(s) for s in Fstr]
    out = {}
    for n in nlist:
        best = None
        nbad = 0
        cnt = 0
        for dx in (d, d + 1):
            for e in range(-emax, emax + 1):
                A = mp.matrix(n, n)
                for i in range(n):
                    for j in range(n):
                        A[i, j] = F[abs((i + base) * dx - (j + base) * d - e)]
                dg = mpf(1)
                for i in range(n):
                    dg *= A[i, i]
                nm = mp.det(A) / dg
                cnt += 1
                if nm <= 0:
                    nbad += 1
                if best is None or nm < best[0]:
                    best = (nm, dx, e)
        out[n] = (S(best[0], 12), best[1], best[2], nbad, cnt)
    return out


def eps_constants(gam, K, P, dps):
    mp.dps = dps
    H = mpf(3) * 10 ** 12
    Sb = (mp.log(H / (2 * mp.pi)) + 2) / (2 * mp.pi * H)
    ig = [1 / g ** 2 for g in gam[:K]]
    nd = gl_nodes(4)
    h = mpf(2)
    I = [mpf(0)] * 3
    F0 = mpf(0)
    for p in range(P['epan']):
        c = (2 * p + 1) * h
        for x, w in nd:
            t = c + h * x
            t2 = t * t
            Pi = mp.exp(-mp.fsum(mp.log(1 + t2 * q) for q in ig))
            wt = h * w * Pi
            F0 += wt
            ex = 2 * Sb * t2 + Sb ** 2 * t2 * t2
            for j in range(3):
                I[j] += wt * t ** j * ex
    return Sb, [v / mp.pi for v in I], F0 / mp.pi


def cauchy_binet_checks():
    out = {}
    mp.dps = 60
    rng = random.Random(20261016)
    m = 7
    gam = sorted(mpf(rng.uniform(0.5, 8.0)) for _ in range(m))
    c = [mpf(rng.choice([-1, 1]) * rng.uniform(0.2, 5.0)) for _ in range(m)]
    x = mpf('0.37')
    rows = {}
    for n in (3, 4):
        D = [mp.fsum(c[k] * (-gam[k]) ** r * mp.exp(-gam[k] * x) for k in range(m)) for r in range(2 * n - 1)]
        H = mp.det(mp.matrix([[D[i + j] for j in range(n)] for i in range(n)]))
        CB = mp.fsum(mp.fprod(c[k] for k in Sset) * mp.fprod((gam[b] - gam[a]) ** 2 for a, b in combinations(Sset, 2)) * mp.exp(-x * mp.fsum(gam[k] for k in Sset))
                     for Sset in combinations(range(m), n))
        rows[n] = dict(H=mp.nstr(H, 20), cauchy_binet=mp.nstr(CB, 20), relative_difference=mp.nstr(abs(H - CB) / abs(H), 5))
    out['random_exponential_sum'] = dict(terms=m, x='0.37', rows=rows)
    mp.dps = 200
    m2 = 10
    gam2 = [mpf(2 * k + 1) for k in range(m2)]
    c2 = [(-1) ** k * (1 + mpf(rng.random())) for k in range(m2)]
    x2 = mpf(8)
    signs = {}
    for n in range(1, 9):
        D = [mp.fsum(c2[k] * (-gam2[k]) ** r * mp.exp(-gam2[k] * x2) for k in range(m2)) for r in range(2 * n - 1)]
        H = mp.det(mp.matrix([[D[i + j] for j in range(n)] for i in range(n)]))
        lead = mp.fprod(c2[k] for k in range(n)) * mp.fprod((gam2[b] - gam2[a]) ** 2 for a, b in combinations(range(n), 2)) * mp.exp(-x2 * mp.fsum(gam2[k] for k in range(n)))
        signs[n] = dict(sign_H=int(mp.sign(H)), sign_leading_term=int(mp.sign(lead)), claimed=(-1) ** (n * (n - 1) // 2), ratio_H_over_leading=mp.nstr(H / lead, 10),
                        floor_n_half_sign=(-1) ** (n // 2))
    out['sign_claim'] = dict(x=8, gamma=[int(v) for v in gam2], alternating_residues=True, rows=signs)
    return out


def main():
    global WORKERS
    ap = argparse.ArgumentParser(description='density of X = Z sqrt(2Y) and the log concavity and PF_n numerics of Chapter 16')
    ap.add_argument('--quick', action='store_true', help='small case in one process; writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/polya_density.json; quick: a /tmp file)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    quick = a.quick
    WORKERS = 1 if (quick or not sys.platform.startswith('linux')) else max(1, a.workers)
    out = Path(a.out or (str(Path(tempfile.gettempdir()) / 'polya_density_quick.json') if quick else 'data/polya_density.json'))
    if quick:
        P = dict(DH=40, nH=11, gH=5, DL=20, nL=10, gL=4, DC=15, KR=16, Cwd=2, Cdeg=3, Ctol=mpf('1e-9'), fig=40, gx=4, dmax=3, seeds=2, nlarge=0, epan=20,
                 thetas=[('0.99', True)])
    else:
        P = dict(DH=80, nH=17, gH=6, DL=30, nL=12, gL=5, DC=30, KR=60, Cwd=1, Cdeg=4, Ctol=mpf('1e-25'), fig=400, gx=6, dmax=25, seeds=20, nlarge=3000, epan=30,
                 thetas=[('0.99', True), ('0.9', False), ('0.5', False)])
    DH, DL = P['DH'], P['DL']
    K12, K14, K40 = 12, 14, min(40, P['KR'])
    R = dict(script='polya_density.py', mode='quick' if quick else 'full')
    R['parameters'] = {k: (str(v) if not isinstance(v, (int, str, list)) else v) for k, v in P.items()}
    R['parameters'].update(workers=WORKERS, residues_used_beyond_1_6=dict(chapter=K14, check=K40))
    zf = find_zeros(a.zeros)
    zs = [ln.strip() for ln in zf.read_text().split() if ln.strip()][:1700]
    mp.dps = DH
    G['DH'], G['DL'] = DH, DL
    G['H'] = build_table(DH, P['nH'], 10, P['gH'], 6)
    log('high table: %d nodes on t in [0, %d] at %d digits' % (G['H']['n'], G['H']['tmax'], DH))
    G['L'] = build_table(DL, P['nL'], 10, P['gL'], 4)
    log('low table: %d nodes on t in [0, %d] at %d digits' % (G['L']['n'], G['L']['tmax'], DL))
    mp.dps = DH
    gam = [mpf(s) for s in zs]
    x0 = XI0()
    log('xi(1/2) = %s (chapter 0.497120778188314)' % fs(x0, 20))
    R['xi_half'] = S(x0, 40)

    # residues and the constants of the text
    c, imr = residues(zs, P['KR'], DH, max(1, P['KR'] // 4))
    mp.dps = DH
    A12 = c[0] * c[1] * (gam[1] - gam[0]) ** 2
    log('residues: %d zeros; largest |Im Xi\'|/|Re Xi\'| %s' % (P['KR'], fs(imr, 3)))
    log('c_1 = %s (chapter 359.524), c_2 = %s (chapter -28007.2), c_3 = %s, c_4 = %s' % (fs(c[0], 10), fs(c[1], 10), fs(c[2], 8), fs(c[3], 8)))
    log('c_1 c_2 (gamma_2 - gamma_1)^2 = %s (chapter -4.776e8); signs alternate: %s' % (fs(A12, 10), all((c[k] > 0) == (k % 2 == 0) for k in range(P['KR']))))
    R['residues'] = dict(K=P['KR'], gamma_1_to_14=[S(g, 30) for g in gam[:14]], c_1_to_14=[S(v, 30) for v in c[:14]], c_all=[S(v, 25) for v in c],
                         max_imag_over_real=fs(imr, 3), A12=S(A12, 25), chapter=dict(c1='359.524', c2='-28007.2', A12='-4.776e8'),
                         signs_alternate=all((c[k] > 0) == (k % 2 == 0) for k in range(P['KR'])))
    ks = [k for k in (1, 2, 5, 10, 20, 30, 40, 60) if k <= P['KR']]
    growth = [(k, fs(mp.log(abs(c[k - 1])) - mp.pi * gam[k - 1] / 4, 8), fs(mp.log(gam[k - 1]), 6), fs(mp.log(abs(c[k - 1])) / gam[k - 1], 6)) for k in ks]
    log('growth of the residues: log|c_k| - pi gamma_k/4 (compare log gamma_k) and log|c_k|/gamma_k (pi/4 = %s): %s' % (
        fs(mp.pi / 4, 6), '; '.join('k=%d: %s (%s), %s' % t for t in growth)))
    ph160 = x0 / xi_s(mpf(1) / 2 + 160)
    ph170 = x0 / xi_s(mpf(1) / 2 + 170)
    log('phi(t^2) at t = 160: %s (chapter: below 1e-80), at t = 170: %s (end of the high table at t = %d)' % (fs(ph160, 5), fs(ph170, 5), G['H']['tmax']))
    R['residue_growth'] = dict(rows=[dict(k=t[0], log_abs_c_minus_pi_gamma_over_4=t[1], log_gamma=t[2], log_abs_c_over_gamma=t[3]) for t in growth], pi_over_4=S(mp.pi / 4, 10))
    R['phi_decay'] = dict(phi_at_160=S(ph160, 8), phi_at_170=S(ph170, 8), chapter='below 1e-80 at t = 160')
    mj = [mp.fsum(G['H']['Wt'][j]) for j in range(7)]
    log('m_0..m_4 = %s (chapter 1.899, 7.321, 44.94, 355.8, 3361); f_X(0) = m_0, f_X\'\'(0) = -m_2' % ', '.join(fs(v, 5) for v in mj[:5]))
    R['moments'] = dict(m_0_to_6=[S(v, 25) for v in mj], chapter=['1.899', '7.321', '44.94', '355.8', '3361'])
    TB = (gam[11] + gam[12]) / 2
    log('T_B = (gamma_12 + gamma_13)/2 = %s (chapter 57.90)' % fs(TB, 10))

    # quadrature against the residue expansion, the Lemma residue bound
    XL = ['3/10', '1/2', '7/10', '19/20', '1/1', '11/10', '6/5', '13/10', '7/5', '3/2', '8/5', '2/1', '21/10']
    xl = [X(s) for s in XL]
    DQ = quad_derivs('H', XL, 3, per=2)
    mp.dps = DH
    cross = {}
    for K in (K12, K14, K40):
        rows = []
        for xv, D in zip(xl, DQ):
            if xv < mpf('0.95') or xv > mpf('1.6'):
                continue
            Dr = rder(c, gam, K, xv, 3)
            rows.append((fs(xv, 4), [fs(abs(D[j] - Dr[j]) / abs(D[j]), 4) for j in range(4)]))
        cross[str(K)] = rows
        log('quadrature vs %d residues, relative difference of f, f\', f\'\', f\'\'\' at x = 0.95: %s; at x = 1.6: %s' % (K, ' '.join(rows[0][1]), ' '.join(rows[-1][1])))
    R['cross_check_quadrature_vs_residues'] = cross
    C23, urange = Cconst(mpf(23), P)
    log('C_0, C_1, C_2 at T = 23: %s (chapter 5.99e4, 1.38e6, 3.20e7); integral to u = %d' % (', '.join(fs(v, 6) for v in C23), urange))
    CB_, urange = Cconst(TB, P)
    log('C_0, C_1, C_2 at T = T_B: %s; integral to u = %d' % (', '.join(fs(v, 6) for v in CB_), urange))
    mp.dps = DC = P['DC']
    XiT = abs(xi_s(mpc(mpf(1) / 2, mpf(23))))
    mp.dps = DH
    R['C'] = dict(T23=dict(C=[S(v, 15) for v in C23], chapter=['5.99e4', '1.38e6', '3.20e7']), TB=dict(T=S(TB, 20), C=[S(v, 15) for v in CB_]), abs_Xi_23=fs(XiT, 8),
                  digits=DC, tolerance=str(P['Ctol']))
    lem = {}
    for name, T, K, CC in (('T23', mpf(23), 2, C23), ('TB', TB, 12, CB_)):
        worst = [mpf(0)] * 3
        rows = []
        for xv, D in zip(xl, DQ):
            Dr = rder(c, gam, K, xv, 2)
            rr = [abs(D[j] - Dr[j]) / (CC[j] * mp.exp(-T * xv)) for j in range(3)]
            worst = [max(w, r) for w, r in zip(worst, rr)]
            rows.append((fs(xv, 4), [fs(r, 4) for r in rr]))
        lem[name] = dict(max_ratio=[fs(w, 5) for w in worst], rows=rows)
        log('Lemma residue bound, T = %s with %d residues: largest |E^(j)|/(C_j e^(-Tx)) over 13 abscissae: %s (must be <= 1)' % (fs(T, 6), K, ' '.join(fs(w, 4) for w in worst)))
    R['lemma_residue_check'] = lem

    # the certificate
    resA = pmap(_t_certA, P['thetas'])
    mp.dps = DH
    for r in resA:
        if not r['ok']:
            log('part (A) FAILED', r)
        else:
            log('part (A) theta = %s: %d cells (chapter 1087), h from %s (at x = %s) to %s, smallest relative margin %s, N(0) = %s, N at last cell start %s' % (
                r['theta'], r['ncells'], r['h_min'], r['h_min_at'], r['h_max'], r['min_relative_margin'], r['N_at_0'], r['N_last_start']))
    R['certificate_A'] = dict(runs=[{k: v for k, v in r.items() if k != 'cells'} for r in resA], chapter_cells=1087,
                              first_run_cells_x_h_N=[r['cells'] for r in resA if r['ok'] and r['cells']][:1])
    cb = cert_B(c, gam, CB_, TB)
    log('part (B): greedy cells %d, equal cells needed %s, 26 equal cells pass: %s (chapter 26)' % (cb['greedy_cells'], cb['uniform_cells_needed'], cb['chapter_26_uniform_cells_pass']))
    log('part (B): S(0.95) = %s (chapter -2.18e8), remainder bound at 0.95 = %s (chapter 7.0e6); S(2.1) = %s, remainder at 2.1 = %s; max of S + R on [0.95, 2.1] = %s at x = %s' % (
        fs(cb['S_at_095'], 6), fs(cb['Rm_at_095'], 6), fs(cb['S_at_21'], 6), fs(cb['Rm_at_21'], 6), fs(cb['worst_S_plus_R'], 6), fs(cb['worst_at'], 6)))
    cc = cert_C(c, gam, C23)
    log('part (C): remainder majorant at x = 2.1 is %s (chapter 4.687e8), |c_1 c_2 (gamma_2 - gamma_1)^2| = %s (chapter 4.776e8); the majorant falls below it from x = %s' % (
        fs(cc['maj_at_21'], 6), fs(abs(cc['A12']), 6), fs(cc['x_where_maj_equals_A12'], 6)))
    R['certificate_B'] = dict(ok=cb['ok'], greedy_cells=cb['greedy_cells'], uniform_cells_needed=cb['uniform_cells_needed'], chapter_26_uniform_cells_pass=cb['chapter_26_uniform_cells_pass'],
                              S_at_095=S(cb['S_at_095'], 10), remainder_at_095=S(cb['Rm_at_095'], 10), S_at_21=S(cb['S_at_21'], 10), remainder_at_21=S(cb['Rm_at_21'], 10),
                              A12=S(cb['A12'], 10), max_S_plus_R=S(cb['worst_S_plus_R'], 10), max_at=S(cb['worst_at'], 8), first_greedy_cells=cb['greedy_first_cells'],
                              chapter=dict(cells=26, S_at_095='-2.18e8', remainder='7.0e6'))
    R['certificate_C'] = dict(majorant_at_21=S(cc['maj_at_21'], 10), A12=S(cc['A12'], 10), x_where_majorant_equals_abs_A12=S(cc['x_where_maj_equals_A12'], 8),
                              chapter=dict(majorant='4.687e8', A12='-4.776e8'))
    mp.dps = DH
    # direct check of N e^{(g1+g2)x} against S(x) + remainder and against A12
    chk = []
    for xs_, xv, D in zip(XL, xl, DQ):
        Ne = (D[0] * D[2] - D[1] ** 2) * mp.exp((gam[0] + gam[1]) * xv)
        mp.dps = 50
        dB = abs(Ne - cb['Sval'](xv)) / cb['Rm'](xv)
        dC = abs(Ne - A12) / cc['maj'](xv)
        mp.dps = DH
        chk.append(dict(x=xs_, N_exp=fs(Ne, 8), S=fs(cb['Sval'](xv), 8), ratio_B=fs(dB, 4), ratio_C=fs(dC, 4)))
    R['certificate_direct_check'] = chk
    log('N e^((gamma_1+gamma_2)x) minus S(x), over the (B) remainder bound, at x = 0.95: %s, at x = 1.6: %s; minus A12 over the (C) majorant at x = 1.6: %s, at x = 2.1: %s' % (
        chk[3]['ratio_B'], chk[10]['ratio_B'], chk[10]['ratio_C'], chk[12]['ratio_C']))

    # -(log f)'' at 0.95 and 4 and the series of the figure
    mp.dps = DH
    Dx = quad_derivs('H', ['19/20', '4/1'], 2, per=1)
    mp.dps = DH

    def L2(D):
        return (D[1] ** 2 - D[0] * D[2]) / D[0] ** 2
    asym = lambda xv: abs(c[1]) / c[0] * (gam[1] - gam[0]) ** 2 * mp.exp(-(gam[1] - gam[0]) * xv)
    lq95 = L2(Dx[0])
    lr95 = L2(rder(c, gam, K12, mpf('0.95'), 2))
    lq4 = L2(Dx[1])
    lr4 = L2(rder(c, gam, K12, mpf(4), 2))
    log('-(log f)\'\' at x = 0.95: quadrature %s, 12 residues %s (chapter about 2.9); at x = 4: quadrature %s, 12 residues %s, leading term %s (chapter 4e-9)' % (
        fs(lq95, 8), fs(lr95, 8), fs(lq4, 8), fs(lr4, 8), fs(asym(mpf(4)), 8)))
    npts = P['fig']
    xsf = ['%d/%d' % (4 * m, npts - 1) for m in range(npts)]
    Df = quad_derivs('H', xsf, 2, per=5)
    mp.dps = DH
    qser = [L2(D) for D in Df]
    rser = [L2(rder(c, gam, K12, X(s), 2)) if X(s) >= mpf('0.95') else None for s in xsf]
    dif = [abs(q - r) / q for q, r in zip(qser, rser) if r is not None]
    log('figure series: %d points; minimum of -(log f)\'\' over the series %s at x = %s; positive everywhere: %s; largest relative gap quadrature vs 12 residues on [0.95, 4]: %s' % (
        npts, fs(min(qser), 6), fs(X(xsf[qser.index(min(qser))]), 5), all(q > 0 for q in qser), fs(max(dif), 4)))
    R['log_concavity'] = dict(L2_at_095=dict(quadrature=S(lq95, 15), residues12=S(lr95, 15), chapter='about 2.9'),
                              L2_at_4=dict(quadrature=S(lq4, 15), residues12=S(lr4, 15), leading_term=S(asym(mpf(4)), 15), chapter='4e-9'),
                              L2_at_095_leading_term=S(asym(mpf('0.95')), 10),
                              series=dict(points=npts, x=[S(X(s), 8) for s in xsf], quadrature=[S(q, 10) for q in qser],
                                          note='x = 4 m/(points - 1) for m = 0..points - 1; the figure draws quadrature (solid) for x <= 0.95 and residues12 (dashed) for x >= 0.95; both columns are recorded where they exist',
                                          solid_quadrature_x_le_0_95=[S(q, 10) if X(s) <= mpf('0.95') else None for q, s in zip(qser, xsf)],
                                          residues12=[None if r is None else S(r, 10) for r in rser], min_over_series=S(min(qser), 10),
                                          all_positive=all(q > 0 for q in qser), max_relative_gap_on_0_95_to_4=fs(max(dif), 5)))

    # the grid x = 0.02 m for m = 0..100 (quadrature), used for the second differences, the scan table and the Hankel grid
    mp.dps = DH
    Dg = quad_derivs('H', ['%d/50' % m for m in range(101)], 6, per=2)
    mp.dps = DH
    F = [D[0] for D in Dg]
    lg = [mp.log(F[m]) for m in range(81)]
    sd = [lg[m + 1] - 2 * lg[m] + lg[m - 1] for m in range(1, 80)]
    sd0 = 2 * (lg[1] - lg[0])
    mx = max(range(79), key=lambda i: sd[i])
    log('second differences of log f_X, x = 0.02 m, m = 1..79: all negative %s; largest %s at x = %s (chapter -2.6e-5); at x = 0: %s' % (
        all(v < 0 for v in sd), fs(sd[mx], 6), fs(mpf(mx + 1) / 50, 4), fs(sd0, 6)))
    slope_sec = (lg[80] - lg[60]) / mpf('0.4')
    xm = [mpf(m) / 50 for m in range(60, 81)]
    ym = [lg[m] for m in range(60, 81)]
    xb_, yb_ = mp.fsum(xm) / len(xm), mp.fsum(ym) / len(ym)
    slope_ls = mp.fsum((u - xb_) * (v - yb_) for u, v in zip(xm, ym)) / mp.fsum((u - xb_) ** 2 for u in xm)
    log('tail slope of log f_X on [1.2, 1.6]: secant %s, least squares %s (chapter -14.09); -gamma_1 = %s' % (fs(slope_sec, 8), fs(slope_ls, 8), fs(-gam[0], 8)))
    R['second_differences'] = dict(all_negative_interior=all(v < 0 for v in sd), largest=S(sd[mx], 10), largest_at_x=S(mpf(mx + 1) / 50, 6), at_zero=S(sd0, 10),
                                   values=[S(v, 8) for v in sd], chapter_largest='-2.6e-5')
    R['tail_slope'] = dict(secant_1_2_to_1_6=S(slope_sec, 12), least_squares=S(slope_ls, 12), minus_gamma_1=S(-gam[0], 12), chapter='-14.09')

    # total mass and variance
    nd = gl_nodes(P['gx'])
    ha = mpf('0.8')
    xg = [ha * (1 + x) for x, w in nd]
    wg = [ha * w for x, w in nd]
    Df0 = quad_derivs('H', [S(v, DH + 10) for v in xg], 0, per=4)
    mp.dps = DH
    fv = [D[0] for D in Df0]
    aa = mpf('1.6')
    I0 = mp.fsum(w * v for w, v in zip(wg, fv))
    I2 = mp.fsum(w * v * x * x for w, v, x in zip(wg, fv, xg))
    t0 = mp.fsum(c[k] * mp.exp(-gam[k] * aa) / gam[k] for k in range(K40))
    t2 = mp.fsum(c[k] * mp.exp(-gam[k] * aa) * (aa ** 2 / gam[k] + 2 * aa / gam[k] ** 2 + 2 / gam[k] ** 3) for k in range(K40))
    mass = 2 * (I0 + t0)
    var = 2 * (I2 + t2)
    a1 = mp.diff(xi_s, mpf(1) / 2, 2) / (2 * x0)
    log('total mass %s (chapter 1.0000000); variance %s (chapter 0.046209986), 2 a_1 from xi\'\'(1/2)/xi(1/2) = %s, difference %s' % (fs(mass, 15), fs(var, 15), fs(2 * a1, 15), fs(var - 2 * a1, 4)))
    R['mass_variance'] = dict(mass=S(mass, 30), mass_minus_1=fs(mass - 1, 5), variance=S(var, 30), two_a1=S(2 * a1, 30), difference=fs(var - 2 * a1, 5),
                              integral_0_to_1_6=S(I0, 25), tail_residues=K40, chapter=dict(mass='1.0000000', variance='0.046209986', two_a1='0.0462099862'))

    # random ordered triples
    tasks = [(DL, s, 60, 0.8) for s in range(P['seeds'])]
    rr = pmap(_t_rand, tasks)
    mp.dps = DH
    mins = [mpf(r[0]) for r in rr]
    srt = sorted(mins)
    log('normalised 3x3 determinants of 60 random ordered triples in [-0.8, 0.8]: smallest over %d seeds: %s; per seed minimum: lowest %s, median %s, highest %s; seed 0: %s (chapter 1.2e-3); non positive: %d' % (
        P['seeds'], fs(srt[0], 5), fs(srt[0], 5), fs(srt[len(srt) // 2], 5), fs(srt[-1], 5), fs(mins[0], 5), sum(r[1] for r in rr)))
    big = None
    if P['nlarge']:
        rl = pmap(_t_rand, [(DL, 1000 + s, P['nlarge'] // 12, 0.8) for s in range(12)])
        mp.dps = DH
        bm = min(rl, key=lambda r: mpf(r[0]))
        big = dict(triples=12 * (P['nlarge'] // 12), smallest=bm[0], at_x=bm[2], at_y=bm[3], non_positive=sum(r[1] for r in rl))
        log('normalised 3x3 determinants of %d random ordered triples: smallest %s, non positive: %d' % (big['triples'], big['smallest'], big['non_positive']))
    R['random_triples'] = dict(per_seed_minimum_of_60=[S(v, 8) for v in mins], seed0=S(mins[0], 8), lowest=S(srt[0], 8), median=S(srt[len(srt) // 2], 8), highest=S(srt[-1], 8),
                               non_positive=sum(r[1] for r in rr), large_sample=big, chapter='1.2e-3 for one unspecified sample of 60')

    # approximation constants
    Sb, eps, f0b = eps_constants(gam, 300, P, 30)
    mp.dps = 30
    bounds = {n: n * mp.factorial(n) * (2 * eps[0] + f0b) ** (n - 1) * eps[0] for n in (3, 4)}
    log('S bound %s (chapter 1.533e-12); eps_0, eps_1, eps_2 = %s (chapter 1.57e-10, 1.3e-9, 1.3e-8); f_0(0) <= %s (chapter 1.976)' % (fs(Sb, 6), ', '.join(fs(v, 4) for v in eps), fs(f0b, 6)))
    log('minor lower bounds -n n! (2 eps_0 + f_0(0))^(n-1) eps_0: n = 3: -%s (chapter -1.1e-8); n = 4: -%s (chapter -1.2e-7)' % (fs(bounds[3], 4), fs(bounds[4], 4)))
    R['approximation'] = dict(S=S(Sb, 10), eps=[S(v, 10) for v in eps], f0_bound=S(f0b, 10), minor_bound_n3=S(-bounds[3], 6), minor_bound_n4=S(-bounds[4], 6), K=300,
                              chapter=dict(S='1.533e-12', eps=['1.57e-10', '1.3e-9', '1.3e-8'], f0='1.976', n3='-1.1e-8', n4='-1.2e-7'))
    mp.dps = DH

    # the scan of Toeplitz type minors
    Ftab = {}
    for name, K in (('K40', K40), ('K14', K14)):
        Ft = []
        for m in range(141):
            Ft.append(F[m] if m <= 80 else rder(c, gam, K, mpf(m) / 50, 0)[0])
        Ftab[name] = Ft
    gap = max(abs(Ftab['K14'][m] - Ftab['K40'][m]) / Ftab['K40'][m] for m in range(81, 141))
    gapq = max(abs(Ftab['K40'][m] - F[m]) / F[m] for m in range(81, 101))
    log('table of f_X(0.02 m): residues with 14 against 40 zeros for m = 81..140: largest relative difference %s; 40 residues against quadrature for m = 81..100: %s' % (fs(gap, 4), fs(gapq, 4)))
    scan = {}
    dmax = P['dmax']
    for name, base, K in (('main_K40_index0', 0, 'K40'), ('chapter_K14_index0', 0, 'K14'), ('K40_index1', 1, 'K40')):
        Fs = [S(v) for v in Ftab[K]]
        rs = pmap(_t_scan, [(d, Fs, DH, (3, 4), 60, base) for d in range(1, dmax + 1)])
        mp.dps = DH
        sc = {}
        for n in (3, 4):
            best = min(((mpf(r[n][0]), d + 1, r[n][1], r[n][2]) for d, r in enumerate(rs)), key=lambda t: t[0])
            sc[str(n)] = dict(smallest=S(best[0], 10), d=best[1], dx=best[2], e=best[3], non_positive=sum(r[n][3] for r in rs), minors=sum(r[n][4] for r in rs),
                              per_d_minimum=[mp.nstr(mpf(r[n][0]), 6) for r in rs])
        scan[name] = sc
        log('scan %s: n = 3: %d minors, smallest %s at d = %d, dx = %d, e = %d, non positive %d; n = 4: %d minors, smallest %s at d = %d, dx = %d, e = %d, non positive %d' % (
            name, sc['3']['minors'], sc['3']['smallest'], sc['3']['d'], sc['3']['dx'], sc['3']['e'], sc['3']['non_positive'], sc['4']['minors'], sc['4']['smallest'],
            sc['4']['d'], sc['4']['dx'], sc['4']['e'], sc['4']['non_positive']))
    R['scan'] = dict(runs=scan, chapter=dict(minors=6050, smallest_n3='1.6e-10', smallest_n4='3.5e-19'),
                     table_gap_K14_vs_K40=fs(gap, 5), table_gap_K40_vs_quadrature=fs(gapq, 5))

    # Hankel grid
    hk = {}
    for name, K in (('chapter_K14', K14), ('check_K40', K40)):
        v3, v4 = [], []
        for m in range(201):
            if m <= 80:
                D = Dg[m]
            else:
                D = rder(c, gam, K, mpf(m) / 50, 6)
            h3, h4 = hankel_pair(D)
            v3.append(h3)
            v4.append(h4)
        i3 = min(range(201), key=lambda i: v3[i])
        i4 = min(range(201), key=lambda i: v4[i])
        hk[name] = dict(v3=v3, v4=v4, i3=i3, i4=i4)
        log('Hankel grid (%s): -H3/f^3 positive everywhere %s, smallest %s at x = %s; H4/f^4 positive everywhere %s, smallest %s at x = %s (chapter 1.05e-21 and 2.08e-39 at x = 4)' % (
            name, all(v > 0 for v in v3), fs(v3[i3], 6), fs(mpf(i3) / 50, 4), all(v > 0 for v in v4), fs(v4[i4], 6), fs(mpf(i4) / 50, 4)))
    qa = []
    for m in range(81, 101):
        h3, h4 = hankel_pair(Dg[m])
        qa.append((m, h3, h4))
    d14, d40 = [], []
    for m, h3, h4 in qa:
        for K, dl in ((K14, d14), (K40, d40)):
            r3, r4 = hankel_pair(rder(c, gam, K, mpf(m) / 50, 6))
            dl.append(max(abs(r3 - h3) / h3, abs(r4 - h4) / h4))
    log('Hankel at x = 1.62..2.0: relative difference of the quantities, residues against quadrature: 14 zeros %s, 40 zeros %s' % (fs(max(d14), 4), fs(max(d40), 4)))
    R['hankel'] = dict(grid_step='0.02', range='[0,4]', chapter=dict(smallest_minus_H3_over_f3='1.05e-21', smallest_H4_over_f4='2.08e-39', at='x = 4'),
                       runs={name: dict(minus_H3_over_f3_min=S(v['v3'][v['i3']], 10), at_x3=S(mpf(v['i3']) / 50, 6), H4_over_f4_min=S(v['v4'][v['i4']], 10), at_x4=S(mpf(v['i4']) / 50, 6),
                                        all_positive_3=all(t > 0 for t in v['v3']), all_positive_4=all(t > 0 for t in v['v4']),
                                        minus_H3_over_f3=[S(t, 8) for t in v['v3']], H4_over_f4=[S(t, 8) for t in v['v4']]) for name, v in hk.items()},
                       quadrature_vs_residues_x_1_62_to_2_0=dict(K14=fs(max(d14), 5), K40=fs(max(d40), 5)))
    cbk = cauchy_binet_checks()
    mp.dps = DH
    for n, rw in cbk['random_exponential_sum']['rows'].items():
        log('Cauchy Binet on a random exponential sum of %d terms, n = %d: relative difference %s' % (cbk['random_exponential_sum']['terms'], n, rw['relative_difference']))
    log('sign claim at x = 8 with alternating residues: sign H_n = (-1)^(n(n-1)/2) for n = 1..8: %s; equals the sign of the term S = {1..n}: %s' % (
        all(r['sign_H'] == r['claimed'] for r in cbk['sign_claim']['rows'].values()), all(r['sign_H'] == r['sign_leading_term'] for r in cbk['sign_claim']['rows'].values())))
    R['cauchy_binet'] = cbk
    R['timing_seconds'] = fs(mpf(time.time() - T0), 5)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(R, indent=1) + '\n')
    log('wrote', out)
    out.with_suffix('.log').write_text('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
