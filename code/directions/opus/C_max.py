#!/usr/bin/env python3
"""Track C, script 4. Maxima, medians and symmetrised maxima of ferromagnetic spins.

(a) max(U1,U2): closed form, list of zeros, no imaginary zero of the moment generating function
(b) maxima of sums of uniforms: zero counts by the argument principle
(c) medians: uniform (symmetric beta, ferromagnetic) against the triangle (not ferromagnetic)
(d) symmetrised maxima of n triangles, real zeros against the box count
(e) the Polya law as a recentred maximum of n iid symmetric spins: necessary conditions and floor
(f) best total variation of a recentred maximum of n iid Gaussian spins to the Polya law

Run: timeout 500 python3 -I C_max.py
Status: computed in double precision (mpmath for the closed form), not certified.
"""
import numpy as np
import mpmath as mp
from scipy import optimize, special

XI0 = 0.4971207781883141
mp.mp.dps = 30


def p_density(x):
    x = np.abs(np.asarray(x, dtype=float))
    s = np.zeros_like(x)
    for n in range(1, 8):
        y = np.pi * n**2 * np.exp(2 * x)
        s += (2 * np.pi**2 * n**4 * np.exp(4.5 * x) - 3 * np.pi * n**2 * np.exp(2.5 * x)) * np.exp(-y)
    return 2 * s / XI0


def uniform_sum_density(ncomp, xs):
    """Density and distribution function of U_1+...+U_ncomp on the grid xs (uniform step), by the spline formula."""
    # f is the B spline of order ncomp, scaled: U on [-1,1]
    f = np.zeros_like(xs)
    F = np.zeros_like(xs)
    n = ncomp
    for k in range(n + 1):
        c = (-1) ** k * special.comb(n, k)
        t = np.maximum(xs + n - 2 * k, 0.0)
        f += c * t ** (n - 1)
        F += c * t**n
    f = f / (2**n * special.factorial(n - 1))
    F = F / (2**n * special.factorial(n))
    return f, np.clip(F, 0, 1)


def gl_grid(a, b, panels, order=32):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    h = (b - a) / panels
    mids = a + (np.arange(panels) + 0.5) * h
    xs = (mids[:, None] + nodes[None, :] * h / 2).ravel()
    ws = np.tile(weights * h / 2, panels)
    return xs, ws


def winding(fun, x0, x1, Y, n=2000):
    pts = np.concatenate([
        x0 + (x1 - x0) * np.arange(n) / n - 1j * Y,
        x1 + 1j * (-Y + 2 * Y * np.arange(n) / n),
        x1 - (x1 - x0) * np.arange(n) / n + 1j * Y,
        x0 + 1j * (Y - 2 * Y * np.arange(n) / n),
    ])
    pts = np.append(pts, pts[0])
    vals = fun(pts)
    return np.sum(np.angle(vals[1:] / vals[:-1])) / (2 * np.pi)


def winding_up(fun, x0, x1, y0, y1, n=2000):
    pts = np.concatenate([
        x0 + (x1 - x0) * np.arange(n) / n + 1j * y0,
        x1 + 1j * (y0 + (y1 - y0) * np.arange(n) / n),
        x1 - (x1 - x0) * np.arange(n) / n + 1j * y1,
        x0 + 1j * (y1 - (y1 - y0) * np.arange(n) / n),
    ])
    pts = np.append(pts, pts[0])
    vals = fun(pts)
    return np.sum(np.angle(vals[1:] / vals[:-1])) / (2 * np.pi)


def part_a():
    print("(a) M=max(U1,U2), density (1+x)/2 on [-1,1]: L(z) = (z e^z - sinh z)/z^2")
    f = lambda z: z * mp.exp(z) - mp.sinh(z)
    zs = []
    for k in range(1, 9):
        # fixed point of e^{2z} = -1/(2z-1) on the sheet k, then Newton
        z = mp.mpc(-1, (k + 0.5) * mp.pi)
        for _ in range(200):
            z = (mp.log(mp.mpf(1) / (2 * z - 1)) + 1j * mp.pi * (2 * k + 1)) / 2
        assert abs(f(z)) < 1e-20 * abs(z)
        zs.append(z)
    print("  zeros z with Im z>0:", [mp.nstr(z, 7) for z in zs])
    print("  Re z < 0 and Re z ~ -(1/2) log(2|z|): the zeros leave the imaginary axis logarithmically")
    # exact check on phi at theta: Re phi = sin(th)/th, Im phi = (sin(th)-th cos(th))/th^2
    th = np.linspace(0.01, 200, 200000)
    mod = np.sqrt((np.sin(th) / th) ** 2 + ((np.sin(th) - th * np.cos(th)) / th**2) ** 2)
    print(f"  min |phi(theta)| * theta on (0,200] = {np.min(mod*th):.6f} (no real zero)")


def part_b():
    print("(b) maxima of sums of uniforms: density n F^(n-1) f; count zeros of the transform in a box, and min|phi| on the real line")
    for ncomp, n in ((2, 2), (3, 2), (2, 3), (3, 3), (6, 2)):
        c = float(ncomp)
        panels = 16 * ncomp
        pts = []
        ws_all = []
        for i in range(2 * ncomp):
            xs, ws = gl_grid(-c + i * 1.0, -c + (i + 1) * 1.0, 8)
            pts.append(xs)
            ws_all.append(ws)
        xs = np.concatenate(pts)
        ws = np.concatenate(ws_all)
        f, F = uniform_sum_density(ncomp, xs)
        dens = n * F ** (n - 1) * f
        mean = np.sum(ws * dens * xs)
        phi = lambda th: np.array([np.sum(ws * dens * np.exp(1j * t * (xs - mean))) for t in np.atleast_1d(th)])
        realg = np.arange(0.05, 60, 0.01)
        mn = np.min(np.abs(phi(realg)))
        w = winding(phi, 0.25, 60.25, 8.0, n=1500)
        print(f"  S=sum of {ncomp} uniforms, max of {n}: min |phi| on (0,60) = {mn:.3e}; zeros in (0.25,60.25)x(-8,8): {w:.2f}")


def part_c():
    print("(c) medians of three: uniform gives (3/4)(1-x^2), symmetric beta with zeros tan th = th (all real). Triangle:")
    xs, ws = gl_grid(0, 2, 40)
    f, F = uniform_sum_density(2, xs)
    dens = 6 * F * (1 - F) * f
    phi = lambda th: np.array([2 * np.sum(ws * dens * np.cos(t * xs)) for t in np.atleast_1d(th)])
    print(f"  mass check {2*np.sum(ws*dens):.12f}; theta^2 phi(theta) at 200, 400, 800 (limit 3/4 by the kink at 0):",
          np.round(phi(np.array([200.0, 400.0, 800.0])) * np.array([200.0, 400.0, 800.0]) ** 2, 5), flush=True)
    grid = np.arange(0.05, 200, 0.005)
    v = phi(grid)
    sc = np.where(v[:-1] * v[1:] < 0)[0]
    print(f"  real sign changes on (0,200): {len(sc)}, at {np.round(grid[sc][:12],3)}")
    fc = lambda z: np.array([2 * np.sum(ws * dens * np.cos(t * xs)) for t in z])
    for box in ((0.25, 40.25, 10.0), (40.25, 80.25, 10.0)):
        print(f"  zeros in ({box[0]},{box[1]})x(-{box[2]},{box[2]}): {winding(fc, box[0], box[1], box[2], n=1000):.2f}")


def part_d():
    print("(d) symmetrised maxima: density (n/2)(F^(n-1)+(1-F)^(n-1)) f; n=2 returns the spin. Sign changes count odd zeros only")
    for ncomp in (2, 3):
        xs, ws = gl_grid(0, ncomp, 20 * ncomp)
        f, F = uniform_sum_density(ncomp, xs)
        for n in (2, 3, 4, 6, 10):
            dens = (n / 2) * (F ** (n - 1) + (1 - F) ** (n - 1)) * f
            phi = lambda th: np.array([2 * np.sum(ws * dens * np.cos(t * xs)) for t in np.atleast_1d(th)])
            grid = np.arange(0.05, 60.25, 0.005)
            v = phi(grid)
            sc = int(np.sum(v[:-1] * v[1:] < 0))
            w = winding(phi, 0.25, 60.25, 6.0, n=1200)
            up = winding_up(phi, 0.25, 60.25, 0.25, 6.0, n=1200)
            print(f"  spin = sum of {ncomp} uniforms, n={n:2d}: mass {2*np.sum(ws*dens):.10f}; real sign changes on (0.05,60.25): {sc}; zeros in (0.25,60.25)x(-6,6): {w:.2f}; above the line Im = 1/4, (0.25,60.25)x(0.25,6): {up:.2f}", flush=True)


def part_e():
    print("(e) Polya law as a recentred maximum of n iid symmetric spins: a must satisfy e^{4a}=n and H(-a)=2^{-n}")
    xs = np.linspace(-2.6, 2.6, 520001)
    dx = xs[1] - xs[0]
    pv = p_density(xs)
    H = np.cumsum(pv) * dx - 0.5 * pv * dx
    H = H / H[-1]
    Hf = lambda x: np.interp(x, xs, H)
    print("   n   a_tail=log(n)/4   a_quart: H(-a)=2^-n   H(-a_tail)   2^-n    floor(n) on eps by the symmetry defect")
    for n in (2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 30, 50):
        at = 0.25 * np.log(n)
        target = 2.0 ** (-n)
        # H(-a)=target: use the tail series for large n
        if target > 1e-12:
            aq = optimize.brentq(lambda a: Hf(-a) - target, 0, 2.9)
        else:
            aq = float("nan")
        # floor: minimise over a the sup over y of the first order defect
        best = np.inf
        for a in np.linspace(max(0.0, min(at, aq if aq == aq else at) - 0.25), max(at, aq if aq == aq else at) + 0.25, 201):
            ys = np.linspace(0, 1.5, 3001)
            h1, h2 = Hf(ys - a), Hf(-ys - a)
            good = (h1 > 1e-280) & (h2 > 1e-280)
            if not good.any():
                continue
            A = np.abs(h1[good] ** (1 / n) + h2[good] ** (1 / n) - 1)
            den = (h1[good] ** (1 / n - 1) + h2[good] ** (1 / n - 1)) / n
            best = min(best, np.max(A / den))
        print(f"  {n:3d}  {at:10.5f}  {aq:14.5f}   {Hf(-at):.4e}  {target:.3e}  {best:.3e}")


def part_f():
    print("(f) best TV of a recentred max of n iid N(0,s^2) spins to the Polya law")
    xs = np.linspace(-2.5, 2.5, 50001)
    dx = xs[1] - xs[0]
    pv = p_density(xs)
    for n in (2, 3, 5, 10, 30):
        def tv(par):
            s, a = abs(par[0]), par[1]
            u = (xs + a) / s
            d = n * special.ndtr(u) ** (n - 1) * np.exp(-u * u / 2) / (s * np.sqrt(2 * np.pi))
            return 0.5 * np.sum(np.abs(d - pv)) * dx
        r = optimize.minimize(tv, [0.215, 0.1 * np.log(n)], method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9})
        print(f"  n={n:3d}: TV={r.fun:.4e}, s={abs(r.x[0]):.4f}, a={r.x[1]:.4f}")


if __name__ == "__main__":
    import traceback
    for part in (part_a, part_b, part_c, part_d, part_e, part_f):
        try:
            part()
        except Exception:
            traceback.print_exc()
