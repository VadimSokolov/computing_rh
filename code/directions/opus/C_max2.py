#!/usr/bin/env python3
"""Track C, script 9. The exactly symmetric recentred maximum that agrees with the Polya law on [-a,a].

Given n >= 2, a is defined by H(-a) = 2^{-n}. Q = H on [-a,a], and beyond a by Q(w+2a) = 1-(1-Q(w)^{1/n})^n. Then F(x) =
Q(x-a)^{1/n} is the distribution function of a symmetric spin Y, and max(Y_1,...,Y_n)-a has distribution function Q,
which is symmetric. TV(Q,H) <= 2^{1-n}. Are Q and Y ferromagnetic? Count zeros strictly above the real axis.

Run: timeout 500 python3 -I C_max2.py   (computed in double precision, not certified)
"""
import numpy as np
from scipy import optimize, special

XI0 = 0.4971207781883141


def p_density(x):
    x = np.abs(np.asarray(x, dtype=float))
    s = np.zeros_like(x)
    for n in range(1, 9):
        y = np.pi * n**2 * np.exp(2 * x)
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * x) - 3 * np.pi * n**2 * np.exp(2.5 * x)) * np.exp(-y)
    return 2 * s / XI0


def GH(x):
    """P(X > x) for x >= 0, closed form: sum over n of (pi n^2)^(-1/4) [2 v^(5/4) e^(-v) - (1/2) Gamma(5/4, v)] / xi(1/2), v = pi n^2 e^(2x)."""
    x = np.asarray(x, dtype=float)
    s = np.zeros_like(x)
    for n in range(1, 9):
        v = np.pi * n**2 * np.exp(2 * x)
        s += (np.pi * n**2) ** (-0.25) * (2 * v**1.25 * np.exp(-v) - 0.5 * special.gamma(1.25) * special.gammaincc(1.25, v))
    return s / XI0


def next_level(Q, G, q, n):
    """(Q,G,q) at u -> (Q,G,q) at u + 2a."""
    logQ = np.where(G < 0.5, np.log1p(-G), np.log(Q))
    comp = -np.expm1(logQ / n)  # 1 - Q^(1/n)
    Gn = comp**n
    qn = comp ** (n - 1) * np.exp((1.0 / n - 1) * logQ) * q
    return 1 - Gn, Gn, qn


def Qq(w, a, n):
    """Distribution function, complement and density of Q at w >= -a. All of w in one interval [(2j-1)a,(2j+1)a]."""
    j = int(np.ceil((np.max(w) - a) / (2 * a) - 1e-12)) if np.max(w) > a else 0
    u = w - 2 * j * a
    Gb = np.where(u >= 0, GH(np.abs(u)), 1 - GH(np.abs(u)))
    Qb = np.where(u >= 0, 1 - GH(np.abs(u)), GH(np.abs(u)))
    qb = p_density(u)
    for _ in range(j):
        Qb, Gb, qb = next_level(Qb, Gb, qb, n)
    return Qb, Gb, qb


def panel_nodes(lo, hi, width=0.05, order=24):
    k = max(1, int(np.ceil((hi - lo) / width)))
    x, w = np.polynomial.legendre.leggauss(order)
    h = (hi - lo) / k
    mids = lo + (np.arange(k) + 0.5) * h
    return (mids[:, None] + x[None, :] * h / 2).ravel(), np.tile(w * h / 2, k)


def winding_up(fun, x0, x1, y0, y1, n):
    pts = np.concatenate([
        x0 + (x1 - x0) * np.arange(n) / n + 1j * y0,
        x1 + 1j * (y0 + (y1 - y0) * np.arange(n) / n),
        x1 - (x1 - x0) * np.arange(n) / n + 1j * y1,
        x0 + 1j * (y1 - (y1 - y0) * np.arange(n) / n),
    ])
    pts = np.append(pts, pts[0])
    vals = fun(pts)
    return np.sum(np.angle(vals[1:] / vals[:-1])) / (2 * np.pi)


def laws(n):
    a = optimize.brentq(lambda s: GH(np.array([s]))[0] - 2.0 ** (-n), 1e-9, 2.5, xtol=1e-14)
    # the law Q on w >= 0 and the spin Y on x >= 0
    xs_q, ws_q, q_vals, tv = [], [], [], 0.0
    xs_y, ws_y, y_vals = [], [], []
    wmax = 3.2
    # Q: [0,a], then [(2j-1)a,(2j+1)a]
    edges = [0.0, a]
    while edges[-1] < wmax:
        edges.append(edges[-1] + 2 * a)
    for lo, hi in zip(edges[:-1], edges[1:]):
        w, wt = panel_nodes(lo, hi)
        Q, G, q = Qq(w, a, n)
        xs_q.append(w)
        ws_q.append(wt)
        q_vals.append(q)
        tv += np.sum(wt * np.abs(q - p_density(w)))
    # Y on x >= 0: w = x - a runs through [-a,a], [a,3a], ...
    e2 = [-a, a]
    while e2[-1] < wmax:
        e2.append(e2[-1] + 2 * a)
    for lo, hi in zip(e2[:-1], e2[1:]):
        w, wt = panel_nodes(lo, hi)
        Q, G, q = Qq(w, a, n)
        logQ = np.where(G < 0.5, np.log1p(-G), np.log(Q))
        y = np.exp((1.0 / n - 1) * logQ) * q / n
        xs_y.append(w + a)
        ws_y.append(wt)
        y_vals.append(y)
    return a, tv, [np.concatenate(v) for v in (xs_q, ws_q, q_vals)], [np.concatenate(v) for v in (xs_y, ws_y, y_vals)]


def kink(n):
    """J_0 = q'(a+) - q'(a-) = -((n-1)/n) 2^(n+1) p(a)^2 - 2 p'(a), with the ratio R = |p'(a)| G(a) / p(a)^2 against (n-1)/n."""
    a = optimize.brentq(lambda s: GH(np.array([s]))[0] - 2.0 ** (-n), 1e-9, 2.5, xtol=1e-14)
    h = 1e-3
    xs = a + h * np.array([-2.0, -1.0, 1.0, 2.0])
    v = p_density(xs)
    dp = (v[0] - 8 * v[1] + 8 * v[2] - v[3]) / (12 * h)
    pa = p_density(np.array([a]))[0]
    J0 = -((n - 1) / n) * 2.0 ** (n + 1) * pa**2 - 2 * dp
    return a, J0, abs(dp) * 2.0 ** (-n) / pa**2


def far_counts(data, T, Ys):
    xs, ws, d = data
    keep = (d > 1e-300) & (xs < 2.8)
    xs, ws, d = xs[keep], ws[keep], d[keep]
    phi = lambda th: np.array([2 * np.sum(ws * d * np.cos(t * xs)) for t in np.atleast_1d(th)])
    return [winding_up(phi, 0.25, T, 0.25, Y, 1000) for Y in Ys]


def count(data, T, Y):
    xs, ws, d = data
    keep = d > 1e-300
    xs, ws, d = xs[keep], ws[keep], d[keep]
    phi = lambda th: np.array([2 * np.sum(ws * d * np.cos(t * xs)) for t in np.atleast_1d(th)])
    grid = np.arange(0.05, T, 0.01)
    v = phi(grid)
    sc = int(np.sum(v[:-1] * v[1:] < 0))
    up = winding_up(phi, 0.25, T, 0.25, Y, 900)
    return 2 * np.sum(ws * d), sc, up


def main():
    print("GH(0) =", GH(np.array([0.0]))[0], "(must be 1/2)")
    print("  n      a_n     2^(1-n)     TV(Q,p)  | Q: mass, real sign changes, zeros above the axis | spin Y: mass, sign changes, zeros above the axis  (box (0.25,40)x(0.25,3))")
    for n in (2, 3, 4, 6, 8, 10, 15, 20, 30):
        a, tv, dq, dy = laws(n)
        mq, sq, uq = count(dq, 40.0, 3.0)
        my, sy, uy = count(dy, 40.0, 3.0)
        print(f"  {n:2d}  {a:.5f}  {2.0**(1-n):.3e}  {tv:.4e} | {mq:.8f} {sq:3d} {uq:6.2f} | {my:.8f} {sy:3d} {uy:6.2f}", flush=True)
    print("first kink of Q_n: n, a_n, J_0, R = |p'(a)| G(a) / p(a)^2, (n-1)/n")
    for n in range(2, 31):
        a, J0, R = kink(n)
        print(f"  {n:2d} {a:.5f} {J0:+.4e} {R:.6f} {(n-1)/n:.6f}", flush=True)
    print("zeros of phi_Q in (0.25,45) x (0.25,Y), Y = 10, 30, 100, 200: nonzero means zeros off the real axis")
    for n in (2, 3, 4, 8, 15, 20):
        a, tv, dq, dy = laws(n)
        print(f"  n={n:2d}: Q {np.round(far_counts(dq, 45.0, (10.0, 30.0, 100.0, 200.0)), 2)}; spin Y {np.round(far_counts(dy, 45.0, (10.0, 30.0, 100.0, 200.0)), 2)}", flush=True)


if __name__ == "__main__":
    main()
