# yo_matrix.py: the prime side matrix of Weil's damped form on a window [-a, a], shared by the scripts of Chapter 17 (Yoshida's Hermitian Forms, Weil Positivity and the Thorin Program), of the deconvolution chapter (Legendre basis, damping y = alpha - 1/2) and of Zhu's chapter (cos and sin modes, 60 digits, even and odd sectors).
# Added for the book in October 2026, written from the mathematics (the chapter authors' own scripts were neither read nor run).
# Mathematics. For an even g with Fourier transform ghat(t) = int g(x) e^(itx) dx, Weil's functional is W(g) = ghat(i/2) + ghat(-i/2) - (gamma_E + log pi) g(0) + int_0^oo 2 (e^(-2x) g(0) - e^(-x/2) g(x)) / (1 - e^(-2x)) dx - 2 sum_n Lambda(n) n^(-1/2) g(log n), applied to the even part of g. For h supported in [-a, a] and a damping y = alpha - 1/2 >= 0 the damped form is Q_y(h) = W(e^(-y|x|) (h * h~)) with h~(x) = conj(h(-x)); on a basis phi_i the matrix is M_ij = W(e^(-y|x|) (phi_i * phi_j~)). Since g = e^(-y|x|) (phi_i * phi_j~) has support in [-2a, 2a], only the prime powers q <= e^(2a) enter, with weight Lambda(q) q^(-alpha) g0(log q), and the integral to infinity is the integral over [0, 2a] plus the closed form tail -g(0) log(1 - e^(-4a)).
# Under RH (and for alpha >= 1 unconditionally) Q_y(h) = (1/pi) int |hhat(theta)|^2 Re xi'/xi(alpha + i theta) dtheta, so the least eigenvalue is at least 2 xi'(alpha)/xi(alpha); at alpha = 1/2 it collapses (below 1e-15 from a = 0.8 on).
# Bases (orthonormal in L^2[-a, a]). 'fourier': chi_n(u) = (2a)^(-1/2) e^(i pi n u/a), |n| <= size, a matrix of order 2 size + 1 indexed n = -size..size. 'legendre': phi_j(u) = sqrt((2j+1)/(2a)) P_j(u/a), j = 0..size-1. 'cossin' (real sector bases): even sector psi_0 = (2a)^(-1/2), psi_j = a^(-1/2) cos(pi j u/a), j = 0..size-1 (size functions, highest order size-1); odd sector phi_j = a^(-1/2) sin(pi j u/a), j = 1..size (size functions, highest order size). So 'size cos modes and size sin modes up to order K' means size = K + 1 even and size = K odd.
# The matrix is real symmetric in all three bases: the even part of chi_n * chi_m~ is real and symmetric in (n, m), and the form commutes with u -> -u, so it splits into an even and an odd sector. For the Fourier basis the entries are closed forms in one dimensional integrals: for n != m, M_nm = (-1)^(n-m) (T(m) - T(n)) / (2 pi (n - m)) with T(m) = 4 int_0^(2a) e^(-yx) cosh(x/2) sin(pi m x/a) dx - 2 int_0^(2a) e^(-(y+1/2)x) sin(pi m x/a)/(1 - e^(-2x)) dx - 2 sum_q Lambda(q) q^(-alpha) sin(pi m log q/a), T(-m) = -T(m) (a Loewner matrix, a divided difference of the odd function T); the diagonal entries are quadratures of the same kind. The cos and sin sector blocks are M_jl +- M_j,-l (with sqrt 2 M_0l in the first row of the even block). The Legendre entries come from the exact correlations (Gauss Legendre on the overlap [x - a, a], nt = size + 6 nodes, exact for polynomials).
# Quadrature. The integrals over [0, 2a] use Gauss Legendre with `nodes` nodes (default: enough for the highest frequency pi size/a and for the poles of 1/(1 - e^(-2x)) at x = +-i pi, with a margin; the removable singularity at x = 0 is never hit since the nodes are interior, and the cancellation near 0 costs about 5 digits of the working precision, below the 15 guard digits added in extended precision).
# Arithmetic. dps <= 15: numpy float64 (the result is an ndarray). dps > 15: mpmath at dps + 15 digits (the result is an mpmath matrix). At alpha = 1/2 the least eigenvalue is below 1e-15 from a = 0.8 on and needs dps >= 50 (more for larger a: check by running at two values of dps). For y > 0 double precision is enough: the floor is 2 xi'(alpha)/xi(alpha) >= 0.0046 for alpha >= 0.55, the entries are accurate to about 1e-13 and the least eigenvalue to about 1e-12 (tested here against 30 digits, see --selftest).
# INTERFACE (stable; import with sys.path.insert(0, directory of this file); import yo_matrix as ym). Common arguments: a = half length of the window; the damping is given by alpha (y = alpha - 1/2, default alpha = 0.5) or by y (not both); a, alpha, y may be float, int, str or mpmath numbers (a float is read as its shortest decimal, so 0.35 means 0.35 exactly in extended precision).
#   ym.build_matrix(a, alpha=None, y=None, basis='fourier', size=12, sector=None, dps=15, nodes=None) -> ndarray or mpmath matrix. sector None gives the full matrix for 'fourier' and 'legendre'; 'even' or 'odd' gives the sector block (for 'legendre' the rows and columns with j even or odd; for 'fourier' the cos and sin block with jmax = size, of order size + 1 and size; 'cossin' requires a sector).
#   ym.spectrum(same arguments) -> ascending eigenvalues (ndarray of floats, or a list of mpmath numbers), computed sector by sector when sector is None.
#   ym.least_eigenvalue(same arguments) -> the least eigenvalue (float or mpmath number).
#   ym.weil_functional(g, X, alpha=None, y=None, dps=15, nodes=None) -> W(e^(-y|x|) g) for an even scalar function g (called as g(x) with x >= 0, a float or an mpmath number), with g negligible beyond X (X = 2a for a compactly supported autocorrelation): the prime powers q <= e^X, the integral over [0, X] and the tail -g(0) log(1 - e^(-2X)).
#   ym.symbol_exact(alpha, tau, dps=20) -> Re xi'/xi(alpha + i tau) (mpmath); ym.b_alpha(alpha, dps=30) -> xi'(alpha)/xi(alpha); ym.prime_powers(L) -> (q list, p list) of the prime powers q = p^k <= e^L.
# Usage: timeout 60 python3 -I code/yo_matrix.py --selftest [--json FILE] [--log FILE] runs the unit tests (under 30 s on one core) and writes data/yo_matrix.json and data/yo_matrix.log (or ./yo_matrix.json and ./yo_matrix.log outside the repository).
import argparse
import cmath
import contextlib
import json
import math
import os
import sys
import time

import mpmath
import numpy as np

GUARD = 15
EULER = 0.5772156649015328606
_CACHE = {}


def _gl_refine(x0, n, wp):
    """Gauss Legendre nodes and weights at working precision wp digits, by Newton's method from double precision nodes (inside workdps)."""
    iters = int(math.ceil(math.log2(wp / 12.0))) + 1
    xs = [None] * n
    ws = [None] * n
    one = mpmath.mpf(1)
    for i in range(n // 2, n):
        x = mpmath.mpf(0) if (n % 2 == 1 and i == n // 2) else mpmath.mpf(float(x0[i]))
        for it in range(iters + 1):
            p0, p1 = one, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            if it < iters:
                x = x - p1 / dp
        xs[i], ws[i] = x, 2 / ((1 - x * x) * dp * dp)
        xs[n - 1 - i], ws[n - 1 - i] = -x, ws[i]
    sx = np.empty(n, dtype=object)
    sw = np.empty(n, dtype=object)
    for i in range(n):
        sx[i], sw[i] = xs[i], ws[i]
    return sx, sw


class Arith:
    """Backend: numpy float64 for dps <= 15, mpmath at dps + GUARD digits (numpy object arrays) above."""

    def __init__(self, dps):
        self.dps = int(dps)
        self.ext = self.dps > 15
        if self.ext:
            self.wp = self.dps + GUARD
            uf = lambda f: np.frompyfunc(f, 1, 1)
            self.exp, self.expm1, self.sin, self.cos = uf(mpmath.exp), uf(mpmath.expm1), uf(mpmath.sin), uf(mpmath.cos)
            self.cosh, self.log, self.sqrt = uf(mpmath.cosh), uf(mpmath.log), uf(mpmath.sqrt)
        else:
            self.exp, self.expm1, self.sin, self.cos = np.exp, np.expm1, np.sin, np.cos
            self.cosh, self.log, self.sqrt = np.cosh, np.log, np.sqrt

    def ctx(self):
        return mpmath.workdps(self.wp) if self.ext else contextlib.nullcontext()

    def num(self, v):
        if not self.ext:
            return float(v)
        if isinstance(v, mpmath.mpf):
            return +v
        if isinstance(v, float):
            return mpmath.mpf(repr(v))
        return mpmath.mpf(v)

    def vec(self, vals):
        if not self.ext:
            return np.array([float(v) for v in vals])
        out = np.empty(len(vals), dtype=object)
        for i, v in enumerate(vals):
            out[i] = self.num(v)
        return out

    def zeros(self, shape, cplx=False):
        if not self.ext:
            return np.zeros(shape, dtype=complex if cplx else float)
        out = np.empty(shape, dtype=object)
        out.fill(mpmath.mpf(0))
        return out

    def eye(self, n, cplx=False):
        out = self.zeros((n, n), cplx)
        for i in range(n):
            out[i, i] = mpmath.mpf(1) if self.ext else 1.0
        return out

    @property
    def pi(self):
        return +mpmath.pi if self.ext else math.pi

    @property
    def euler(self):
        return +mpmath.euler if self.ext else EULER

    def gl(self, n):
        key = (n, self.wp if self.ext else 0)
        if key not in _CACHE:
            x0, w0 = np.polynomial.legendre.leggauss(n)
            _CACHE[key] = _gl_refine(x0, n, self.wp) if self.ext else (x0, w0)
        return _CACHE[key]


def prime_powers(L):
    """The prime powers q = p^k with log q <= L, ascending, as two integer lists (q, p)."""
    n = int(math.exp(L)) + 1
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    out = []
    for p in np.nonzero(sieve)[0]:
        q = int(p)
        while q <= n:
            if math.log(q) <= L:
                out.append((q, int(p)))
            q *= int(p)
    out.sort()
    return [q for q, _ in out], [p for _, p in out]


def _damping(ar, alpha, y):
    if alpha is not None and y is not None:
        raise ValueError('give alpha or y, not both')
    half = ar.num('0.5')
    if y is None:
        al = half if alpha is None else ar.num(alpha)
        return al, al - half
    yy = ar.num(y)
    return yy + half, yy


def _nq_default(kind, a, size, dps):
    osc = (int(math.ceil(math.pi * size)) + 20) if kind == 'fourier' else (2 * size + 20)
    z = complex(-1.0, math.pi / float(a))
    r = max(abs(z + cmath.sqrt(z * z - 1)), abs(z - cmath.sqrt(z * z - 1)))
    pole = int(math.ceil((dps + 5) * 2.302585 / (2 * math.log(r)))) + 10
    return max(60, osc + pole)


def _quad_data(ar, X, y, nq):
    """Nodes x in (0, X), weights, and the factors e^(-yx), 1 - e^(-2x), e^(-(y+1/2)x), e^(-2x), cosh(x/2)."""
    s, w = ar.gl(nq)
    x = (s + 1) * (X / 2)
    wx = w * (X / 2)
    e1 = ar.exp(x * (-y))
    om = -ar.expm1(x * (-2))
    ea = e1 * ar.exp(x * (-0.5))
    e2 = ar.exp(x * (-2))
    ch = ar.cosh(x / 2)
    return x, wx, e1, om, ea, e2, ch


def _prime_data(ar, X, alpha):
    qs, ps = prime_powers(float(X))
    if not qs:
        return None, None
    lq = ar.log(ar.vec(qs))
    lam = ar.log(ar.vec(ps))
    return lq, lam * ar.exp(lq * (-alpha))


def _assemble(ar, X, y, alpha, corr, g0, nq, cplx=False):
    """W(e^(-y|x|) G) for a matrix valued even function G supported in [-X, X]: corr(x) is G(x) at x >= 0, g0 = G(0)."""
    n = g0.shape[0]
    x, wx, e1, om, ea, e2, ch = _quad_data(ar, X, y, nq)
    pk = wx * e1 * ch * 4 - wx * ea * 2 / om
    rk = wx * e2 * 2 / om
    acc = ar.zeros((n, n), cplx)
    for k in range(nq):
        acc = acc + corr(x[k]) * pk[k] + g0 * rk[k]
    tail = ar.log(-ar.expm1(X * (-2)))
    acc = acc - g0 * (ar.euler + ar.log(ar.pi) + tail)
    lq, pw = _prime_data(ar, X, alpha)
    if lq is not None:
        for k in range(len(lq)):
            acc = acc - corr(lq[k]) * (pw[k] * 2)
    return acc


def weil_functional(g, X, alpha=None, y=None, dps=15, nodes=None):
    """W(e^(-y|x|) g) for an even scalar g negligible beyond X (see the header)."""
    ar = Arith(dps)
    with ar.ctx():
        alpha_, y_ = _damping(ar, alpha, y)
        X_ = ar.num(X)

        def mat(v):
            out = ar.zeros((1, 1))
            out[0, 0] = v
            return out

        nq = nodes or max(200, _nq_default('legendre', X_ / 2, 30, dps))
        res = _assemble(ar, X_, y_, alpha_, lambda x: mat(g(x)), mat(g(ar.num(0))), nq)
        return res[0, 0]


def _fourier_vectors(ar, a, y, alpha, N, nq):
    """T(m) (m = 0..N) and the diagonal D(n) (n = 0..N) of the Fourier basis matrix (see the header)."""
    x, wx, e1, om, ea, e2, ch = _quad_data(ar, 2 * a, y, nq)
    idx = ar.vec(list(range(N + 1)))
    th = np.multiply.outer(idx, x * (ar.pi / a))
    S, C = ar.sin(th), ar.cos(th)
    tri = 1 - x / (2 * a)
    lq, pw = _prime_data(ar, 2 * a, alpha)
    if lq is not None:
        thp = np.multiply.outer(idx, lq * (ar.pi / a))
        Ps = ar.sin(thp).dot(pw)
        Pc = ar.cos(thp).dot(pw * (1 - lq / (2 * a)))
    else:
        Ps = Pc = ar.zeros(N + 1)
    T = S.dot(wx * e1 * ch) * 4 - S.dot(wx * ea / om) * 2 - Ps * 2
    arch = (e2[None, :] - (ea * tri)[None, :] * C) * 2 / om[None, :]
    const = ar.euler + ar.log(ar.pi) + ar.log(-ar.expm1(a * (-4)))
    D = C.dot(wx * e1 * ch * tri) * 4 + arch.dot(wx) - Pc * 2 - const
    return T, D


def _block_np(T, D, kind, jmax):
    """Double precision, vectorised: the full matrix (kind None, n, m = -jmax..jmax) or the cos (even) and sin (odd) sector block."""
    def tf(i):
        return np.sign(i) * T[np.abs(i)]

    def M(n, m):
        k = n - m
        ks = np.where(k == 0, 1, k)
        off = np.where(k % 2 == 1, -1.0, 1.0) * (tf(m) - tf(n)) / (2 * math.pi * ks)
        return np.where(k == 0, D[np.abs(n)], off)

    if kind is None:
        r = np.arange(-jmax, jmax + 1)
        return M(r[:, None], r[None, :])
    if kind == 'even':
        j = np.arange(0, jmax + 1)
        B = M(j[:, None], j[None, :]) + M(j[:, None], -j[None, :])
        B[0, 0] /= 2
        B[0, 1:] /= math.sqrt(2)
        B[1:, 0] /= math.sqrt(2)
        return B
    j = np.arange(1, jmax + 1)
    return M(j[:, None], j[None, :]) - M(j[:, None], -j[None, :])


def _fourier_entry(T, D, pi):
    def tof(n):
        return T[n] if n >= 0 else -T[-n]

    def entry(n, m):
        if n == m:
            return D[abs(n)]
        k = n - m
        v = (tof(m) - tof(n)) / (2 * pi * k)
        return -v if (k & 1) else v

    return entry


def _fourier_block(entry, kind, jmax, s2):
    if kind == 'even':
        idx = range(0, jmax + 1)

        def f(j, l):
            if j == 0 and l == 0:
                return entry(0, 0)
            if j == 0 or l == 0:
                return s2 * entry(j, l)
            return entry(j, l) + entry(j, -l)
    else:
        idx = range(1, jmax + 1)

        def f(j, l):
            return entry(j, l) - entry(j, -l)
    return [[f(j, l) for l in idx] for j in idx]


def _legendre_corr(ar, a, M, nt):
    st, wt = ar.gl(nt)
    nf = ar.sqrt(ar.vec([2 * j + 1 for j in range(M)]) / (2 * a))

    def phi(t):
        z = t / a
        P = ar.zeros((M, len(z)))
        P[0] = z * 0 + 1
        if M > 1:
            P[1] = z
        for k in range(1, M - 1):
            P[k + 1] = (z * P[k] * (2 * k + 1) - P[k - 1] * k) / (k + 1)
        return P * nf[:, None]

    def corr(x):
        half = a - x / 2
        t = st * half + x / 2
        g = phi(t).dot(phi(t - x).T * wt[:, None]) * half
        return (g + g.T) / 2

    return corr


def _build(a, alpha, y, basis, size, sector, ar, nodes):
    """Internal: the matrix as a nested list (fourier, cossin) or ndarray (legendre) in the backend arithmetic."""
    alpha_, y_ = _damping(ar, alpha, y)
    a_ = ar.num(a)
    if basis in ('fourier', 'cossin'):
        if basis == 'cossin':
            if sector not in ('even', 'odd'):
                raise ValueError("basis 'cossin' needs sector 'even' or 'odd'")
            jmax = size - 1 if sector == 'even' else size
        else:
            jmax = size
        nq = nodes or _nq_default('fourier', a_, jmax, ar.dps)
        T, D = _fourier_vectors(ar, a_, y_, alpha_, jmax, nq)
        if not ar.ext:
            return _block_np(T, D, sector if (basis == 'cossin' or sector is not None) else None, jmax)
        entry = _fourier_entry(list(T), list(D), ar.pi)
        if basis == 'fourier' and sector is None:
            rng = range(-jmax, jmax + 1)
            return [[entry(n, m) for m in rng] for n in rng]
        return _fourier_block(entry, sector, jmax, ar.sqrt(ar.num(2)))
    if basis == 'legendre':
        nq = nodes or _nq_default('legendre', a_, size, ar.dps)
        corr = _legendre_corr(ar, a_, size, size + 6)
        B = _assemble(ar, 2 * a_, y_, alpha_, corr, ar.eye(size), nq)
        if sector is not None:
            keep = [j for j in range(size) if j % 2 == (0 if sector == 'even' else 1)]
            B = B[np.ix_(keep, keep)]
        return B.tolist()
    raise ValueError('basis must be fourier, legendre or cossin')


def _export(ar, B):
    if ar.ext:
        return mpmath.matrix(B)
    return np.array(B, dtype=float)


def build_matrix(a, alpha=None, y=None, basis='fourier', size=12, sector=None, dps=15, nodes=None):
    ar = Arith(dps)
    with ar.ctx():
        return _export(ar, _build(a, alpha, y, basis, size, sector, ar, nodes))


def _eigs(ar, B):
    if not ar.ext:
        return list(np.linalg.eigvalsh(np.array(B, dtype=float)))
    E = mpmath.eigsy(mpmath.matrix(B), eigvals_only=True)
    return sorted(E[i] for i in range(len(E)))


def spectrum(a, alpha=None, y=None, basis='fourier', size=12, sector=None, dps=15, nodes=None):
    ar = Arith(dps)
    with ar.ctx():
        if sector is None and basis == 'fourier':
            blocks = [_build(a, alpha, y, 'fourier', size, s, ar, nodes) for s in ('even', 'odd')]
        elif sector is None and basis == 'cossin':
            blocks = [_build(a, alpha, y, 'cossin', size, s, ar, nodes) for s in ('even', 'odd')]
        elif sector is None and basis == 'legendre':
            B = _build(a, alpha, y, 'legendre', size, None, ar, nodes)
            blocks = [[[B[i][j] for j in range(par, size, 2)] for i in range(par, size, 2)] for par in (0, 1)]
        else:
            blocks = [_build(a, alpha, y, basis, size, sector, ar, nodes)]
        ev = sorted(sum((_eigs(ar, B) for B in blocks if len(B) > 0), []))
    return np.array(ev) if not ar.ext else ev


def least_eigenvalue(a, alpha=None, y=None, basis='fourier', size=12, sector=None, dps=15, nodes=None):
    return spectrum(a, alpha, y, basis, size, sector, dps, nodes)[0]


def symbol_exact(alpha, tau, dps=20):
    """Re xi'/xi(alpha + i tau) with xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s): 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'/zeta(s)."""
    with mpmath.workdps(dps):
        s = mpmath.mpc(alpha, tau)
        v = 1 / s + 1 / (s - 1) - mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s / 2) / 2 + mpmath.zeta(s, derivative=1) / mpmath.zeta(s)
        return v.real


def b_alpha(alpha, dps=30):
    """xi'(alpha)/xi(alpha); at alpha = 1 the closed form 1 + gamma_E/2 - log(pi)/2 - log 2 (the poles of 1/(s-1) and zeta'/zeta cancel)."""
    with mpmath.workdps(dps + 10):
        if alpha == 1:
            v = 1 + mpmath.euler / 2 - mpmath.log(mpmath.pi) / 2 - mpmath.log(2)
        else:
            s = mpmath.mpf(alpha)
            v = 1 / s + 1 / (s - 1) - mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s / 2) / 2 + mpmath.zeta(s, derivative=1) / mpmath.zeta(s)
        return +v


# ---------------------------------------------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------------------------------------------

def _generic_fourier(a, y, N, nq, nt):
    """Independent route for the Fourier basis: the correlations chi_n * chi_m~ at x and -x by quadrature (complex), then the same assembly."""
    ar = Arith(15)
    s, w = np.polynomial.legendre.leggauss(nt)
    ns = np.arange(-N, N + 1)

    def F(xx):
        lo, hi = max(-a, xx - a), min(a, a + xx)
        t = (lo + hi) / 2 + (hi - lo) / 2 * s
        E1 = np.exp(1j * math.pi * np.outer(ns, t) / a) / math.sqrt(2 * a)
        E2 = np.exp(1j * math.pi * np.outer(ns, t - xx) / a) / math.sqrt(2 * a)
        return (E1 * w).dot(E2.conj().T) * (hi - lo) / 2

    def corr(x):
        return (F(float(x)) + F(-float(x))) / 2

    g0 = np.eye(2 * N + 1, dtype=complex)
    return _assemble(ar, 2 * a, y, y + 0.5, corr, g0, nq, cplx=True)


def _panels(lo, hi, width, nn):
    s, w = np.polynomial.legendre.leggauss(nn)
    edges = np.linspace(lo, hi, int(round((hi - lo) / width)) + 1)
    t = np.concatenate([(e0 + e1) / 2 + (e1 - e0) / 2 * s for e0, e1 in zip(edges[:-1], edges[1:])])
    wt = np.concatenate([(e1 - e0) / 2 * w for e0, e1 in zip(edges[:-1], edges[1:])])
    return t, wt


def _hann_ft(theta, a):
    """FT of the vectors h_even = (1 + cos(pi x/a))^2 = 3/2 + 2 cos + cos(2.)/2 and h_odd = sin(pi x/a) + sin(2 pi x/a)/2 on [-a, a]."""
    b = math.pi / a
    ev = 2 * np.sin(theta * a) * (1.5 / theta - 2 * theta / (theta ** 2 - b ** 2) + 0.5 * theta / (theta ** 2 - (2 * b) ** 2))
    od = np.sin(theta * a) * (2 * b / (b ** 2 - theta ** 2) - 2 * b / (4 * b ** 2 - theta ** 2))
    return ev, od


def selftest(json_path=None, log=print):
    t_start = time.time()
    results = []

    def record(name, ok, **info):
        info = {k: (v if isinstance(v, (str, bool)) else repr(float(v)) if isinstance(v, (float, np.floating)) else v) for k, v in info.items()}
        results.append(dict(name=name, passed=bool(ok), **info))
        log('%-4s %-34s %s' % ('PASS' if ok else 'FAIL', name, ' '.join('%s=%s' % (k, v) for k, v in info.items())))

    # 1. Gauss Legendre nodes at 40 digits integrate x^20 exactly (12 and 13 nodes).
    ar = Arith(40)
    with ar.ctx():
        errs = []
        for n in (12, 13):
            s, w = ar.gl(n)
            errs.append(abs(sum(w * s ** 20) - mpmath.mpf(2) / 21))
        record('gauss_legendre_40_digits', max(errs) < mpmath.mpf(10) ** -35, max_error=float(max(errs)))

    # 2. Prime side against zero side for the Gaussian g(x) = exp(-(x/c)^2), c = 1/2: 2 sum_{gamma>0} ghat(gamma), ghat(t) = c sqrt(pi) exp(-c^2 t^2/4).
    c = 0.5
    with mpmath.workdps(35):
        zer = [mpmath.zetazero(k).imag for k in range(1, 7)]
        zero_side = 2 * sum(c * mpmath.sqrt(mpmath.pi) * mpmath.exp(-c * c * g * g / 4) for g in zer)
    w_d = weil_functional(lambda x: math.exp(-(x / c) ** 2), 12.0, y=0.0, dps=15, nodes=240)
    record('weil_gaussian_double', abs(w_d - float(zero_side)) < 1e-13, prime_side=w_d, zero_side=float(zero_side), diff=abs(w_d - float(zero_side)))
    w_e = weil_functional(lambda x: mpmath.exp(-(x / mpmath.mpf(c)) ** 2), 12, y=0, dps=30, nodes=260)
    record('weil_gaussian_30_digits', abs(w_e - zero_side) < mpmath.mpf(10) ** -24, prime_side=mpmath.nstr(w_e, 28), zero_side=mpmath.nstr(zero_side, 28), diff=float(abs(w_e - zero_side)))

    # 3. The closed form Fourier matrix against the independent quadrature route (complex correlations), real symmetric.
    a, al = 0.7, 0.6
    B = build_matrix(a, al, basis='fourier', size=3)
    G = _generic_fourier(a, al - 0.5, 3, 160, 80)
    record('fourier_closed_form_vs_quadrature', np.abs(G.real - B).max() < 1e-12, max_diff=np.abs(G.real - B).max())
    record('hermitian_real_symmetric', np.abs(G.imag).max() < 1e-12 and np.abs(B - B.T).max() == 0 and np.abs(G - G.conj().T).max() < 1e-12,
           max_imag=np.abs(G.imag).max(), max_asymmetry_generic=np.abs(G - G.conj().T).max())

    # 4. The three bases give the same least eigenvalue (algebraic convergence in the size): a = 0.8, alpha = 0.6.
    a, al = 0.8, 0.6
    lf = [least_eigenvalue(a, al, basis='fourier', size=n) for n in (12, 24, 48)]
    lc = min(least_eigenvalue(a, al, basis='cossin', size=49, sector='even'), least_eigenvalue(a, al, basis='cossin', size=48, sector='odd'))
    ll = [least_eigenvalue(a, al, basis='legendre', size=m) for m in (16, 24, 40)]
    record('three_bases_agree', abs(lf[2] - ll[2]) < 2e-6 and abs(lf[2] - lc) < 1e-12,
           fourier_12_24_48='%.10f %.10f %.10f' % tuple(lf), cossin_even_odd_min_48='%.10f' % lc, legendre_16_24_40='%.10f %.10f %.10f' % tuple(ll),
           fourier48_minus_legendre40=lf[2] - ll[2])
    # nodes check: 1.5 times the nodes changes nothing
    n0 = _nq_default('fourier', a, 24, 15)
    d1 = abs(least_eigenvalue(a, al, basis='fourier', size=24, nodes=int(1.5 * n0)) - lf[1])
    n1 = _nq_default('legendre', a, 24, 15)
    d2 = abs(least_eigenvalue(a, al, basis='legendre', size=24, nodes=int(1.5 * n1)) - ll[1])
    record('quadrature_nodes_converged', max(d1, d2) < 1e-13, fourier_change=d1, legendre_change=d2, nodes_fourier=n0, nodes_legendre=n1)

    # 5. Multiplier identity: Q_y(h) = (1/pi) int |hhat|^2 Re xi'/xi(alpha + i theta) dtheta for trigonometric vectors vanishing to high order at +-a.
    a = 1.0
    ce = np.array([1.5 * math.sqrt(2 * a), 2 * math.sqrt(a), 0.5 * math.sqrt(a)])
    co = np.array([math.sqrt(a), 0.5 * math.sqrt(a)])
    cf = np.zeros(5)
    cf[2] = ce[0]
    cf[1] = cf[3] = math.sqrt(2 * a)
    cf[0] = cf[4] = math.sqrt(2 * a) / 4
    xs, ws = np.polynomial.legendre.leggauss(60)
    hann2 = lambda u: (1 + np.cos(math.pi * u / a)) ** 2
    pj = np.array([math.sqrt((2 * j + 1) / (2 * a)) * np.polynomial.legendre.Legendre.basis(j)(xs) for j in range(20)])
    cl = (pj * ws * a).dot(hann2(a * xs))
    for al, width, nn, top in ((1.0, 0.5, 10, 60.0), (0.75, 0.25, 10, 60.0), (0.6, 0.2, 12, 40.0)):
        t0 = time.time()
        th, wt = _panels(0.0, top, width, nn)
        with mpmath.workdps(15):
            v = np.array([float(symbol_exact(al, float(t), 15)) for t in th])
        ev, od = _hann_ft(th, a)
        ref_e = 2 / math.pi * np.sum(wt * ev ** 2 * v)
        ref_o = 2 / math.pi * np.sum(wt * od ** 2 * v)
        be = build_matrix(a, al, basis='cossin', size=3, sector='even')
        bo = build_matrix(a, al, basis='cossin', size=2, sector='odd')
        bf = build_matrix(a, al, basis='fourier', size=2)
        bl = build_matrix(a, al, basis='legendre', size=20)
        q_e, q_o, q_f, q_l = ce @ be @ ce, co @ bo @ co, cf @ bf @ cf, cl @ bl @ cl
        tol = 1e-7
        record('multiplier_identity_alpha_%g' % al, max(abs(q_e - ref_e), abs(q_o - ref_o), abs(q_f - ref_e), abs(q_l - ref_e)) < tol,
               even_integral=ref_e, even_cossin=q_e, even_fourier=q_f, even_legendre=q_l, odd_integral=ref_o, odd_cossin=q_o,
               symbol_evaluations=len(th), seconds=time.time() - t0)

    # 6. Double precision against 30 digits for y > 0: a = 1.5, alpha = 0.55, N = 18, and Legendre 14.
    d_f = least_eigenvalue(1.5, 0.55, basis='fourier', size=18)
    e_f = least_eigenvalue(1.5, 0.55, basis='fourier', size=18, dps=30)
    d_l = least_eigenvalue(1.5, 0.55, basis='legendre', size=14)
    e_l = least_eigenvalue(1.5, 0.55, basis='legendre', size=14, dps=30)
    err = max(abs(d_f - float(e_f)), abs(d_l - float(e_l)))
    record('double_vs_30_digits_alpha_0.55', err < 1e-12, fourier_18_double=d_f, fourier_18_30digits=mpmath.nstr(e_f, 20), legendre_14_double=d_l,
           legendre_14_30digits=mpmath.nstr(e_l, 20), max_diff=err)

    # 7. Floor: b_alpha at 1 against the closed form and 2 b_alpha for the printed values.
    b = [float(b_alpha(al)) for al in (0.55, 0.6, 0.75, 1)]
    record('floor_row_2b', all(abs(2 * x - r) < 6e-8 for x, r in zip(b, (0.0046210, 0.0092418, 0.0231027, 0.0461914))), two_b=' '.join('%.7f' % (2 * x) for x in b))

    ok = all(r['passed'] for r in results)
    log('selftest: %s in %.1f s' % ('all passed' if ok else 'FAILED', time.time() - t_start))
    out = dict(script='code/yo_matrix.py --selftest', seconds=time.time() - t_start, all_passed=ok, tests=results,
               versions=dict(numpy=np.__version__, mpmath=mpmath.__version__, python=sys.version.split()[0]))
    if json_path:
        with open(json_path, 'w') as fh:
            json.dump(out, fh, indent=1)
    return ok


def _default_paths(stem):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    if os.path.basename(here) == 'code' and os.path.isdir(os.path.join(root, 'data')):
        return os.path.join(root, 'data', stem)
    return os.path.join(os.getcwd(), stem)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--json', default=None)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    if args.selftest:
        jpath = args.json or _default_paths('yo_matrix.json')
        lpath = args.log or _default_paths('yo_matrix.log')
        with open(lpath, 'w') as lf_:
            def say(msg):
                print(msg)
                lf_.write(msg + '\n')
                lf_.flush()
            sys.exit(0 if selftest(jpath, say) else 1)
    ap.print_help()
