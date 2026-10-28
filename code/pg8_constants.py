# pg8_constants.py: the constants and the shared quadrature library for note X8 (Chapter 22, Section mi:sec:mixture, "The Polya mixture of
# psi_M"): the residue r_0 of Proposition mi:prop:polar (placeholder [N8.1]) and the checks of Proposition mi:prop:meixner (the identity
# mi:eq:nine, U(0), mass, variance and mean of the clock of q_phi). Added for the book in October 2026, written from the mathematics (the
# collaborator's scripts were not read, run or copied). The module also holds the library used by pg8_meixner.py ([N8.3]), pg8_E.py ([N8.4])
# and pg8_main.py ([N8.5]); those scripts import it, so it must sit beside them (on Hopper pass it as an input file).
# Conventions. Book units: Phi(u) = sum_{n>=1} 2 (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}), analytic and even in |Im u| < pi/4,
# Xi(z) = int Phi(u) e^{izu} du; c_0 = cos(2 pi/9), c_1 = cos(pi/9); r_0 = sin(2 pi/9) Phi(i pi/9) / (18 pi^2).
# Method for r_0. Phi(i pi/9) is the theta series at e^{2u} = e^{2 pi i/9}, whose n-th term is at most 2 (2 pi^2 + 3 pi) n^4 exp(-pi c_0 n^2)
# (mi4.tex, after Theorem mi:thm:tail); the tail n >= N is bounded in ball arithmetic by the first omitted bound times 1/(1 - rho), rho the
# ratio of consecutive bounds at n = N, so r_0 is CERTIFIED (arb ball, radius reported). Phi(i y) is real, so the real part of the partial sum
# is taken and the size of its imaginary part is logged.
# Library (used by the other pg8 scripts). The transform F(z) = int K(u) e^{izu} du of an even kernel K, analytic and decaying in the strip
# |Im u| < pi/4 and real on R, is evaluated on the line Im u = tau = 0.77 (the poles of the multipliers are at height >= pi/2 > tau, and
# the removable singularity of the kernel of Xi^E at i pi/9 lies between the lines): with u = v + i tau and G(v) = K(v + i tau), G(-v) =
# conj G(v), and the trapezoid rule of step 1/H on |v| <= V gives F_s(z) = e^{z tau} F(z) = (1/H)[G_0 + sum_{j>=1} (G_j e^{i z j/H} + conj(G_j)
# e^{-i z j/H})]. F_s has the zeros of F. Error control: the aliasing error is of the size of F_s(z +- 2 pi H) (a factor exp(-0.015 (2 pi H - |z|))
# smaller than the data: 1e-39 at H = 1000, |z| = 750) and the truncation error of the size of K at v = V (below 1e-100 at V = 4.2); both are
# ESTIMATED by a second run at a different step and cut (compare()), not bounded by a theorem, so statements derived from F_s are marked
# computed unless the script says otherwise. Ball arithmetic (python-flint, arb) is used for every evaluation, so the rounding error is
# bounded; the argument principle uses the Lipschitz bound of F_s from the absolute values of the terms (see track()).
# Output: data/pg8_constants.json and data/pg8_constants.log. Usage: python3 pg8_constants.py [--quick] [--out FILE]; ten seconds; the
# library needs python-flint (the book venv on Hopper, or any venv with python-flint 0.9) and runs as `timeout 30 python3 -I` on the laptop.
import argparse
import json
import math
import os
import sys
import time

import mpmath as mp

try:
    import numpy as np
    from flint import acb, arb, ctx
except ImportError:  # the constants part needs only mpmath, the library needs flint and numpy
    np = None
    acb = arb = ctx = None

TAU_NUM, TAU_DEN = 77, 100  # tau = 0.77


class Log:
    def __init__(self, path=None):
        self.f = open(path, "w") if path else None

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        if self.f:
            self.f.write(s + "\n")
            self.f.flush()


def S(x, n=30):
    """arb/acb to a string of n digits, without the radius."""
    return x.str(n, radius=False)


# ---------------------------------------------------------------- kernels (functions of an acb u with Re u >= 0)
def phi_series(u):
    """Phi(u) for Re u >= 0, |Im u| < pi/4, the theta series cut where pi n^2 Re e^{2u} exceeds the working precision."""
    pi = arb.pi()
    e2 = (2 * u).exp()
    a = (u * 9 / 2).exp()
    b = (u * 5 / 2).exp()
    lim = ctx.prec * 0.70 + 30
    s = acb(0)
    n = 1
    while True:
        x = pi * n * n * e2
        s += 2 * (2 * pi * pi * n ** 4 * a - 3 * pi * n * n * b) * (-x).exp()
        if n >= 2 and float(x.real.mid()) > lim + 4 * math.log(n):
            return s
        n += 1


def kern_meixner(phi):
    """Phi(u) (1 + cos phi) / (cosh u + cos phi): the kernel of q_phi * Xi."""
    c = phi.cos()
    return lambda u: phi_series(u) * (1 + c) / (u.cosh() + c)


def kern_E():
    """cosh(9u/2) exp(-2 pi (cosh 2u - c_0)) / (cosh 2u - c_0): the kernel of calE_{2 pi/9} (a_0 = 2 pi)."""
    c0 = (2 * arb.pi() / 9).cos()
    tp = 2 * arb.pi()

    def f(u):
        ch = (2 * u).cosh()
        return (u * 9 / 2).cosh() * (-tp * (ch - c0)).exp() / (ch - c0)

    return f


def kern_main(r0, N, with_E=True):
    """8 pi^2 cosh(9u/2) [sum_{n<=N} n^4 exp(-2 pi n^2 cosh 2u) + r_0 exp(-2 pi (cosh 2u - c_0))/(cosh 2u - c_0)]: xi^(inf) + Xi^E."""
    c0 = (2 * arb.pi() / 9).cos()
    tp = 2 * arb.pi()
    p8 = 8 * arb.pi() ** 2

    def f(u):
        ch = (2 * u).cosh()
        s = acb(0)
        for n in range(1, N + 1):
            s += n ** 4 * (-tp * n * n * ch).exp()
        if with_E:
            s += r0 * (-tp * (ch - c0)).exp() / (ch - c0)
        return p8 * (u * 9 / 2).cosh() * s

    return f


# ---------------------------------------------------------------- the quadrature object
class Quad:
    def __init__(self, fn, H, V):
        tau = arb(TAU_NUM) / TAU_DEN
        self.H, self.J = H, int(math.ceil(V * H))
        self.h = arb(1) / H
        self.G = [fn(acb(arb(j) / H, tau)) for j in range(self.J + 1)]
        self.Gc = [g.conjugate() for g in self.G]
        self.np = np.array([complex(g) for g in self.G])
        self.v = np.arange(self.J + 1) / H
        self.absG = np.abs(self.np)

    def eval(self, z, deriv=False):
        """F_s(z) (and F_s'(z)) as acb balls, z an acb. Each exponential is computed directly: a recurrence w^j in rectangular
        complex balls loses a factor sqrt(2) per step on the unit circle."""
        h = self.h
        iz = acb(0, 1) * z
        s = self.G[0]
        ds = acb(0)
        for j in range(1, self.J + 1):
            vj = arb(j) / self.H
            a = self.G[j] * (iz * vj).exp()
            b = self.Gc[j] * (-iz * vj).exp()
            s += a + b
            if deriv:
                ds += vj * (a - b)
        if deriv:
            return h * s, acb(0, 1) * h * ds
        return h * s

    def taylor(self, z, K):
        """[F_s^(k)(z)/k! for k = 0..K] as acb balls."""
        h = self.h
        iz = acb(0, 1) * z
        sums = [acb(0) for _ in range(K + 1)]
        sums[0] = self.G[0]
        for j in range(1, self.J + 1):
            vj = arb(j) / self.H
            a = self.G[j] * (iz * vj).exp()
            b = self.Gc[j] * (-iz * vj).exp()
            sums[0] += a + b
            pw = acb(1)
            iv = acb(0, 1) * vj
            for k in range(1, K + 1):
                pw *= iv
                sums[k] += pw * (a + (-1) ** k * b)
        fact = 1
        out = []
        for k in range(K + 1):
            if k > 1:
                fact *= k
            out.append(h * sums[k] / fact)
        return out

    def bound(self, Y, k=1):
        """sup over |Im z| <= Y of |F_s^(k)(z)| <= (1/H) sum_j v_j^k |G_j| 2 cosh(Y v_j) (absolute values of the terms), as a float."""
        v = self.v[1:]
        return float(np.sum(v ** k * self.absG[1:] * 2 * np.cosh(Y * v)) / self.H) * 1.0000001

    def scan(self, t):
        """F_s(t) on a float array t in double precision (no ball arithmetic), block-wise."""
        out = np.empty(len(t))
        v = self.v[1:]
        Gp = self.np[1:]
        for i in range(0, len(t), 300):
            tb = t[i:i + 300]
            E = np.exp(1j * np.outer(tb, v))
            out[i:i + 300] = (self.np[0].real + 2 * (E @ Gp).real) / self.H
        return out


def scan_many(quads, t):
    """F_s(t) for several quads with the same nodes in one pass (shared exponentials)."""
    v = quads[0].v[1:]
    Gp = np.stack([q.np[1:] for q in quads], axis=1)
    G0 = np.array([q.np[0].real for q in quads])
    out = np.empty((len(t), len(quads)))
    for i in range(0, len(t), 300):
        E = np.exp(1j * np.outer(t[i:i + 300], v))
        out[i:i + 300] = (G0 + 2 * (E @ Gp).real) / quads[0].H
    return out


def sign_changes(t, f):
    idx = np.where(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    return [(t[i], t[i + 1]) for i in idx]


def refine_real(q, a, b, tol=None):
    """Zero of F_s in [a, b] (floats, sign change) by safeguarded Newton in ball arithmetic; returns an exact arb (the midpoint)."""
    tol = tol or arb(2) ** (-(ctx.prec - 60))
    lo, hi = arb(a), arb(b)
    flo = q.eval(acb(lo)).real
    t = (lo + hi) / 2
    for _ in range(30):
        f, df = q.eval(acb(t), True)
        f, df = f.real, df.real
        if (f > 0) == (flo > 0):
            lo = t
        else:
            hi = t
        tn = t - f / df
        tn = tn.mid()
        if not (lo < tn < hi):
            tn = ((lo + hi) / 2).mid()
        if abs(tn - t) < tol:
            return tn
        t = tn
    return t


def cert_real(q, t, delta):
    """True if F_s(t - delta) and F_s(t + delta) are nonzero balls of opposite signs (a zero in between, certified for F_s)."""
    fa = q.eval(acb(t - delta)).real
    fb = q.eval(acb(t + delta)).real
    return (fa > 0 and fb < 0) or (fa < 0 and fb > 0)


def newton_c(q, z, tol=1e-28, maxit=40):
    z = acb(z)
    for _ in range(maxit):
        f, df = q.eval(z, True)
        dz = f / df
        z = (z - dz).mid() if hasattr(z, "mid") else z - dz
        if float(abs(dz).mid()) < tol:
            break
    return z


def krawczyk(q, z, r, qerr):
    """Interval Newton in the complex plane: X = z +- r (a box), N(X) = m - (F_s(m) + e)/F_s'(X), |e| <= qerr the quadrature
    allowance. Returns (True, N(X)) if N(X) lies in the interior of X, which proves one simple zero of the (discretised) F_s in X."""
    m = z.mid()
    X = acb(arb(m.real.mid(), r), arb(m.imag.mid(), r))
    f = q.eval(m) + acb(arb(0, qerr), arb(0, qerr))
    _, dX = q.eval(X, True)
    Nx = m - f / dX
    ok = (Nx.real.lower() > X.real.lower() and Nx.real.upper() < X.real.upper()
          and Nx.imag.lower() > X.imag.lower() and Nx.imag.upper() < X.imag.upper())
    return bool(ok), Nx


def enclose(q, z, qerr, radii=(1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4)):
    """krawczyk() with the smallest radius of the list for which the interval Newton operator maps the box into its interior."""
    for r in radii:
        ok, Nx = krawczyk(q, z, arb(r), qerr)
        if ok:
            return True, r, Nx
    return False, radii[-1], Nx


# ---------------------------------------------------------------- argument principle with a Lipschitz test on every step
_Q = None  # set before a Pool is created (fork)


def track(task):
    """Change of arg F_s along the segment z0 -> z1 (python complex) in adaptive steps with a Taylor test. At z_k the balls
    c_i = F_s^(i)(z_k)/i! (i <= K) are computed; a step of length L is accepted if sum_{i=1..K} |c_i| L^i + L^(K+1) B_(K+1)/(K+1)! <= SAFE |F_s(z_k)|_lower,
    where B_(K+1) = (1/H) sum v^(K+1) |G| 2 cosh(Y v) bounds sup |F_s^(K+1)| on |Im z| <= Y (absolute values of the terms). Then
    |F_s(z) - F_s(z_k)| < |F_s(z_k)| on the step: no zero, and the change of arg is arg(F_s(z_(k+1)) / F_s(z_k)) in (-pi/6, pi/6).
    Returns (change, steps, smallest accepted step, max radius of the arg balls)."""
    z0, z1, Lmax, SAFE, K = task
    q = _Q
    d = z1 - z0
    total = abs(d)
    z = z0
    cs = q.taylor(acb(z.real, z.imag), K)
    pos, tot, steps, minL, maxrad = 0.0, 0.0, 0, 1e300, 0.0
    fk1 = math.factorial(K + 1)
    while pos < total:
        m = float(cs[0].abs_lower().mid())
        cabs = [float(c.abs_upper().mid()) for c in cs]
        L = min(Lmax, total - pos)
        while True:
            Y = min(max(abs(z.imag), abs(z1.imag)), abs(z.imag) + L)
            R = L ** (K + 1) * q.bound(Y, K + 1) / fk1
            if sum(cabs[i] * L ** i for i in range(1, K + 1)) + R <= SAFE * m:
                break
            L /= 1.4
            if L < 1e-9:
                raise RuntimeError("step collapsed at %r (|F| lower %g)" % (z, m))
        pos2 = pos + L
        z2 = z1 if pos2 >= total * (1 - 1e-12) else z0 + d * (pos2 / total)
        if abs(z2 - z) > L * 1.0000001:  # rounding of the point: the test was for length L
            raise RuntimeError("step longer than tested")
        cs2 = q.taylor(acb(z2.real, z2.imag), K)
        a = (cs2[0] / cs[0]).arg()
        rad = float(a.rad())
        maxrad = max(maxrad, rad)
        if rad > 1e-6:
            raise RuntimeError("ball too wide at %r" % (z2,))
        tot += float(a.mid())
        minL = min(minL, L)
        z, cs, pos = z2, cs2, (total if z2 == z1 else pos2)
        steps += 1
    return tot, steps, minL, maxrad


def argument_count(q, x0, x1, Y, cpus, Lmax=0.25, SAFE=0.5, K=6, chunk=10.0, log=print):
    """Number of zeros of F_s in the rectangle [x0, x1] x [-Y, Y] by the argument principle along its boundary (counterclockwise)."""
    global _Q
    _Q = q
    V = [complex(x0, -Y), complex(x1, -Y), complex(x1, Y), complex(x0, Y), complex(x0, -Y)]
    tasks = []
    for a, b in zip(V[:-1], V[1:]):
        n = max(1, int(math.ceil(abs(b - a) / chunk)))
        pts = [a + (b - a) * (k / n) for k in range(n)] + [b]
        pts = [complex(round(p.real, 9), round(p.imag, 9)) for p in pts]
        pts[0], pts[-1] = a, b
        tasks += [(pts[k], pts[k + 1], Lmax, SAFE, K) for k in range(n)]
    t0 = time.time()
    if cpus > 1:
        import multiprocessing as mpr
        with mpr.get_context("fork").Pool(cpus) as pool:
            res = pool.map(track, tasks, chunksize=1)
    else:
        res = [track(t) for t in tasks]
    tot = sum(r[0] for r in res)
    n = tot / (2 * math.pi)
    steps = sum(r[1] for r in res)
    log("  argument principle: total change of arg / 2pi = %.10f, steps %d, smallest step %.3g, max arg radius %.2g, %.1f s"
        % (n, steps, min(r[2] for r in res), max(r[3] for r in res), time.time() - t0))
    if abs(n - round(n)) > 1e-6:
        raise RuntimeError("winding number not an integer: %r" % n)
    return int(round(n)), steps


def compare(q1, q2, pts):
    """Largest relative difference of F_s between two quadratures (different step and cut) at the real points pts (list of float)."""
    worst = 0.0
    for t in pts:
        a = q1.eval(acb(t)).real.mid()
        b = q2.eval(acb(t)).real.mid()
        worst = max(worst, float(abs(a - b) / max(abs(a), arb("1e-300"))))
    return worst


def grid_minima(q, xs, ys):
    """Local minima of |F_s| on the grid xs (increasing) x ys (increasing, all > 0): list of (x, y). Double precision."""
    v = q.v[1:]
    Gp = q.np[1:]
    A = np.empty((len(ys), len(xs)))
    for i in range(0, len(xs), 400):
        xb = xs[i:i + 400]
        E = np.exp(1j * np.outer(xb, v))
        Ec = E.conj()
        for k, y in enumerate(ys):
            val = q.np[0] + E @ (Gp * np.exp(-y * v)) + Ec @ (Gp.conj() * np.exp(y * v))
            A[k, i:i + 400] = np.abs(val)
    out = []
    ny, nx = A.shape
    for k in range(ny):
        for i in range(1, nx - 1):
            a = A[k, i]
            if a <= A[k, i - 1] and a <= A[k, i + 1] and (k == 0 or a <= A[k - 1, i - 1:i + 2].min()) and (k == ny - 1 or a <= A[k + 1, i - 1:i + 2].min()):
                out.append((float(xs[i]), float(ys[k])))
    return out


def find_nonreal(q, xs, ys, real_zeros, tolre=0.2):
    """Nonreal zeros (Im > 1e-6, Re > 0) of F_s found by Newton from the grid minima; returns a sorted list of acb."""
    cands = grid_minima(q, xs, ys)
    zs = []
    rz = np.array([float(r) for r in real_zeros]) if len(real_zeros) else np.array([])
    for (x, y) in cands:
        if y == ys[0] and len(rz) and np.min(np.abs(rz - x)) < tolre:
            continue
        try:
            z = newton_c(q, acb(x, y))
        except (ZeroDivisionError, ValueError, ArithmeticError):
            continue
        zr, zi = float(z.real.mid()), float(z.imag.mid())
        if abs(zi) > 1e-6 and zr > 0 and abs(zi) < 3.5 and zi > 0:
            if not any(abs(zr - float(w.real.mid())) < 1e-8 and abs(zi - float(w.imag.mid())) < 1e-8 for w in zs):
                zs.append(z)
    zs.sort(key=lambda w: float(w.real.mid()))
    return zs


def safe_edge(x, roots, step, margin=0.06):
    """Move x by multiples of step until no real zero of the list lies within margin (the boundary must not pass near a zero)."""
    rz = [float(r) for r in roots]
    while any(abs(r - x) < margin for r in rz):
        x += step
    return x


def analyse(name, q, q2, tmax, box, cpus, log, quick=False, dt=0.01, xstep=0.05, ap=True, nonreal=True):
    """The common analysis of a transform F_s: real zeros in (0, tmax] by sign changes of F_s on a grid of step dt (double precision) refined by
    Newton in ball arithmetic, each certified by the signs of F_s at root -+ 1e-12 (balls); comparison with the second quadrature q2;
    argument principle on box = (x0, x1, Y) (x0, x1 moved off the real zeros); nonreal zeros from grid minima plus Newton.
    Returns a dict (numbers as strings)."""
    t0 = time.time()
    res = {"name": name}
    t = np.arange(0.005, tmax + 0.5, dt)
    f = q.scan(t)
    sc = sign_changes(t, f)
    log("[%s] grid scan: %d sign changes up to t = %g (double precision), %.1f s" % (name, len(sc), t[-1], time.time() - t0))
    res["scan_sign_changes_total"] = len(sc)
    roots = []
    ncert = 0
    if not quick:
        delta = arb(10) ** -12
        for (a, b) in sc:
            r = refine_real(q, a, b)
            roots.append(r)
            ncert += bool(cert_real(q, r, delta))
        log("[%s] %d zeros refined to ball arithmetic, %d certified as sign changes of F_s at root -+ 1e-12, %.1f s" % (name, len(roots), ncert, time.time() - t0))
    else:
        roots = [arb((a + b) / 2) for (a, b) in sc]
        g0 = np.diff([float(r) for r in roots])
        for k in np.argsort(g0)[:8]:  # quick mode: only the pairs with the smallest gaps are refined (the others stay grid midpoints, +-0.005)
            for i in (int(k), int(k) + 1):
                roots[i] = refine_real(q, sc[i][0], sc[i][1])
    gaps = [float(roots[i + 1] - roots[i]) for i in range(len(roots) - 1)]
    res["roots"] = [S(r, 30) for r in roots]
    res["n_roots_certified"] = ncert
    if gaps:
        k = int(np.argmin(gaps))
        res["least_gap"] = {"gap": "%.12g" % gaps[k], "lower_zero": S(roots[k], 20), "upper_zero": S(roots[k + 1], 20), "index": k + 1}
        log("[%s] least gap between consecutive real zeros: %.10g between %s and %s (zeros number %d and %d)" % (name, gaps[k], S(roots[k], 15), S(roots[k + 1], 15), k + 1, k + 2))
    # comparison with the second quadrature (step and cut differ), at a sample of zeros
    if not quick and q2 is not None and roots:
        sample = [float((roots[i] + roots[i + 1]) / 2) for i in range(0, len(roots) - 1, max(1, len(roots) // 12))]
        d = compare(q, q2, sample)
        res["quadrature_rel_difference_at_zeros_sample"] = "%.3g" % d
        log("[%s] relative difference of F_s between the two quadratures (H=%d,V=%g and H=%d,V=%g) at %d sample points between zeros: %.3g" % (name, q.H, q.J / q.H, q2.H, q2.J / q2.H, len(sample), d))
        # and the difference of the zeros themselves
        zz = [refine_real(q2, float(roots[i]) - 1e-6, float(roots[i]) + 1e-6) for i in range(0, len(roots), max(1, len(roots) // 12))]
        dz = max(float(abs(zz[k] - roots[i])) for k, i in enumerate(range(0, len(roots), max(1, len(roots) // 12))))
        res["quadrature_zero_difference_sample"] = "%.3g" % dz
        log("[%s] largest difference of the zero itself between the two quadratures (sample): %.3g" % (name, dz))
    if ap and box is not None:
        x0, x1, Y = box
        if roots:
            x0 = safe_edge(x0, roots, 0.1)
            x1 = safe_edge(x1, roots, -0.1)
        nreal = sum(1 for r in roots if x0 < float(r) < x1)
        n_ap, steps = argument_count(q, x0, x1, Y, cpus, log=log)
        res["box"] = [x0, x1, Y]
        res["real_zeros_in_box"] = nreal
        res["zeros_in_box_argument_principle"] = n_ap
        res["ap_steps"] = steps
        log("[%s] box [%g,%g]x[-%g,%g]: argument principle %d zeros, %d real (sign changes), so %d nonreal" % (name, x0, x1, Y, Y, n_ap, nreal, n_ap - nreal))
        if nonreal and n_ap != nreal:
            xs = np.arange(x0 + 0.02, x1, xstep)
            ys = np.array([0.002, 0.005, 0.01, 0.02, 0.04, 0.07] + list(np.arange(0.1, Y + 0.05, 0.1)))
            zs = find_nonreal(q, xs, ys, roots)
            res["nonreal_expected_upper"] = (n_ap - nreal) // 2
            res["nonreal_found_upper"] = len(zs)
            res["nonreal_zeros_upper"] = [[S(z.real, 20), S(z.imag, 20)] for z in zs]
            log("[%s] nonreal zeros in the upper half plane found: %d (argument principle expects %d)" % (name, len(zs), (n_ap - nreal) // 2))
            for z in zs:
                log("    %s + %s i" % (S(z.real, 16), S(z.imag, 16)))
    elif nonreal and box is not None and not ap:
        x0, x1, Y = box
        xs = np.arange(x0 + 0.02, x1, xstep)
        ys = np.array([0.002, 0.005, 0.01, 0.02, 0.04, 0.07] + list(np.arange(0.1, Y + 0.05, 0.1)))
        zs = find_nonreal(q, xs, ys, roots)
        res["nonreal_found_upper"] = len(zs)
        res["nonreal_zeros_upper"] = [[S(z.real, 20), S(z.imag, 20)] for z in zs]
        log("[%s] (no argument principle) nonreal zeros in the upper half plane found in [%g,%g]: %d" % (name, x0, x1, len(zs)))
        for z in zs:
            log("    %s + %s i" % (S(z.real, 16), S(z.imag, 16)))
    res["seconds"] = round(time.time() - t0, 1)
    return res


def load_ordinates(path, n=1700):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.strip():
                out.append(mp.mpf(line.split()[0]))
    return out[:n]


# ---------------------------------------------------------------- N8.1
def phi_imag_axis(y, N, prec):
    """Real part of Phi(i y) from the theta series with n < N, plus the rigorous tail bound; arb. Returns (ball, imag_partial, tail)."""
    ctx.prec = prec
    pi = arb.pi()
    u = acb(0, y)
    e2 = (2 * u).exp()
    a = (u * 9 / 2).exp()
    b = (u * 5 / 2).exp()
    s = acb(0)
    for n in range(1, N):
        s += 2 * (2 * pi * pi * n ** 4 * a - 3 * pi * n * n * b) * (-(pi * n * n * e2)).exp()
    c = (2 * y).cos()
    bn = lambda n: 2 * (2 * pi * pi + 3 * pi) * arb(n) ** 4 * (-pi * c * n * n).exp()
    rho = (arb(N + 1) / N) ** 4 * (-pi * c * (2 * N + 1)).exp()
    tail = bn(N) / (1 - rho)
    return s.real + arb(0, tail.upper()), s.imag, tail


def main_constants(args, log):
    mp.mp.dps = 60
    res = {}
    t0 = time.time()
    # --- r_0, certified
    prec = 260
    ctx.prec = prec
    y = arb.pi() / 9
    N = 14
    ph, im, tail = phi_imag_axis(y, N, prec)
    ph2, im2, tail2 = phi_imag_axis(y, N + 4, prec)
    s29 = (2 * arb.pi() / 9).sin()
    r0 = s29 * ph / (18 * arb.pi() ** 2)
    C92 = 4 * arb.pi() / 9 * ph
    log("Phi(i pi/9)        =", ph.str(55), "(ball, includes the tail bound %s)" % tail.str(3))
    log("imaginary part of the partial sum (n<%d) = %s (the full sum is real)" % (N, im.str(5)))
    log("Phi(i pi/9) with N=%d terms agrees to %s" % (N + 4, (ph - ph2).str(3)))
    log("C_{9/2} = 4 pi/9 Phi(i pi/9) =", C92.str(30, radius=False), "(chapter: 3.5527303732)")
    log("r_0 = sin(2pi/9) Phi(i pi/9)/(18 pi^2) =", r0.str(55))
    r0chk = s29 / (8 * arb.pi() ** 3) * C92
    log("r_0 from sin(2pi/9)/(8 pi^3) C_{9/2}: difference %s" % ((r0 - r0chk).str(3)))
    log("4 pi^2 r_0 =", (4 * arb.pi() ** 2 * r0).str(12, radius=False), "(orientation 0.36345 in the reviews)")
    res["r0"] = r0.str(50, radius=False)
    res["r0_ball"] = r0.str(55)
    res["r0_radius"] = r0.rad().str(5)
    res["Phi_i_pi_9"] = ph.str(50, radius=False)
    res["C_9_2"] = C92.str(30, radius=False)
    res["r0_status"] = "certified (arb ball with rigorous tail bound)"
    ctx.prec = 200

    # --- checks in mpmath
    mp.mp.dps = 50
    c1 = mp.cos(mp.pi / 9)
    c0 = mp.cos(2 * mp.pi / 9)
    log("c_1 = cos(pi/9) = %s; c_0 = cos(2pi/9) = %s; c_0 = 2 c_1^2 - 1 : %s" % (mp.nstr(c1, 20), mp.nstr(c0, 20), mp.nstr(c0 - (2 * c1 ** 2 - 1), 3)))
    four = lambda u: (mp.cosh(7 * u / 2) + 2 * c1 * mp.cosh(5 * u / 2) + (1 + 2 * c0) * mp.cosh(3 * u / 2) + (1 + 2 * c1) * mp.cosh(u / 2))
    res["identity_nine"] = {}
    for u in [mp.mpf(0), mp.mpf("0.2"), mp.mpf("0.74"), mp.mpf(2)]:
        if u == 0:
            lhs = 1 / (2 * (1 - c1))
            lhs2 = 1 / (1 - c0)
            rhs2 = four(u) / (1 + c1)
        else:
            lhs = mp.cosh(9 * u / 2) / (2 * (mp.cosh(u) - c1))
            lhs2 = mp.cosh(9 * u / 2) / (mp.cosh(2 * u) - c0)
            rhs2 = four(u) / (mp.cosh(u) + c1)
        d1 = abs(lhs - four(u)) / lhs
        d2 = abs(lhs2 - rhs2) / lhs2
        log("mi:eq:nine at u = %s: U = %s, relative difference to the four term sum %s; second identity %s" % (mp.nstr(u, 3), mp.nstr(lhs, 20), mp.nstr(d1, 3), mp.nstr(d2, 3)))
        res["identity_nine"][str(u)] = [mp.nstr(d1, 3), mp.nstr(d2, 3)]
    log("U(0): 1/(2(1-c_1)) = %s, four term sum = (1+2c_1)^2-type value %s" % (mp.nstr(1 / (2 * (1 - c1)), 25), mp.nstr(four(0), 25)))
    res["U0"] = [mp.nstr(1 / (2 * (1 - c1)), 25), mp.nstr(four(0), 25)]

    res["q_phi"] = {}
    for name, phi in [("0", mp.mpf(0)), ("pi/9", mp.pi / 9), ("pi/4", mp.pi / 4), ("pi/2", mp.pi / 2)]:
        cp = mp.cos(phi)
        if phi == 0:
            q = lambda x: 2 * x / mp.sinh(mp.pi * x) if x != 0 else 2 / mp.pi
        else:
            q = lambda x, phi=phi, cp=cp: (1 + cp) * mp.sinh(phi * x) / (mp.sin(phi) * mp.sinh(mp.pi * x)) if x != 0 else (1 + cp) * phi / (mp.sin(phi) * mp.pi)
        mass = 2 * mp.quad(q, [0, 2, 10, 40, 120])
        var = 2 * mp.quad(lambda x: x * x * q(x), [0, 2, 10, 40, 120])
        mean_sig = mp.nsum(lambda k: 1 / ((2 * k + 1) * mp.pi - phi) ** 2 + 1 / ((2 * k + 1) * mp.pi + phi) ** 2, [0, mp.inf])
        # characteristic function at u = 0.7 against (1+cos)/(cosh u + cos)
        cf = 2 * mp.quad(lambda x: q(x) * mp.cos(mp.mpf("0.7") * x), [0, 2, 10, 40, 120])
        cf0 = (1 + cp) / (mp.cosh(mp.mpf("0.7")) + cp)
        log("q_phi, phi = %s: mass %s; variance %s vs 1/(1+cos) %s; E varsigma = sum %s vs 1/(2(1+cos)) %s; cf(0.7) %s vs %s"
            % (name, mp.nstr(mass, 20), mp.nstr(var, 20), mp.nstr(1 / (1 + cp), 20), mp.nstr(mean_sig, 20), mp.nstr(1 / (2 * (1 + cp)), 20), mp.nstr(cf, 15), mp.nstr(cf0, 15)))
        res["q_phi"][name] = {"mass": mp.nstr(mass, 25), "variance": mp.nstr(var, 25), "variance_expected": mp.nstr(1 / (1 + cp), 25),
                              "E_varsigma": mp.nstr(mean_sig, 25), "E_varsigma_expected": mp.nstr(1 / (2 * (1 + cp)), 25),
                              "cf_0.7": mp.nstr(cf, 20), "cf_0.7_expected": mp.nstr(cf0, 20)}
    log("time %.1f s" % (time.time() - t0))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or (os.path.join(here, "..", "data", "pg8_constants.json") if not a.quick else "/tmp/pg8_constants.json")
    log = Log(os.path.splitext(out)[0] + ".log")
    res = main_constants(a, log)
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    log("wrote", out)
