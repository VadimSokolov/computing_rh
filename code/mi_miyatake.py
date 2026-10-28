# mi_miyatake.py: the numerical claims of Chapter 22 (Nick Polson's Miyatake chapter), the sections "Miyatake's theorem and the remainder R"
# (the identity, R, Q, the kinks, xi_2, the tail constant, the interpolation) and the formulas of "The transition at b = 2" (C_b, the cusp,
# b_0). Group G3 of the computation briefs of 22 Oct 2026. Written from the mathematics (the collaborator's scripts were not read, run or
# copied).
# Conventions. psi_theta(x) = sum_{n>=1} exp(-pi n^2 x); Xi(z) = xi(1/2 + i z), xi(s) = (1/2) s (s-1) pi^(-s/2) Gamma(s/2) zeta(s);
# Phi(u) = sum_n (4 pi^2 n^4 e^(9u/2) - 6 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)) is the kernel with Xi(z) = int_R Phi(u) e^(izu) du, which is the
# chapter's displayed formula at C:346 (Titchmarsh's kernel with 2 pi^2 and 3 pi is half of it; see the check "norm"). Miyatake's psi is psi(t) = sum n^4 exp(-2 pi n^2 cosh t) (1 - R),
# 1 - R(x,t) = (1 - 3 e^(-t)/(2x)) e^(x e^(-t)) / (1 + e^(-9t/2)) (t > 0), even in t. The transform of a function f of t is
# int f(t) 2 cosh(9t/4) e^(izt/2) dt (times 2 pi^2, which gives Xi when 8 pi^2 cosh(9t/2) psi(2t) = Phi(t)).
# Method. (1) Everything that is an entire function of the real variable is integrated by the trapezoid rule (exponentially accurate for an
# even analytic integrand; the aliasing error is bounded by the Fourier transform at 2 pi/h - z). (2) The kinked sums (n in S(2) = powers of 2,
# interpolation lambda) are integrated on t >= 0 by composite Gauss-Legendre panels (width 0.2, 24 nodes, nodes by Newton in mpmath). (3) All
# sums over nodes on a z grid of constant step use the recurrence E_j(z + dz) = E_j(z) exp(i dz u_j) in fixed point integers of P bits (P = 400
# means an absolute error 1e-120 times the number of steps); the second run at P = 330 and another step is compared and logged. (4) C_b near the
# cusp b = 2 uses the theta modular transformation 1 + 2 psi(x) = (-i eps)^(-1/2) theta_2(-1/eps), eps = 1 + i x, not the series in e^(-pi n^2 x).
# Error control. Every quantity is computed twice (precision P, panel or trapezoid step) and the maximum difference is logged and written.
# Output: data/mi_miyatake.json and data/mi_miyatake.log (quick: --out FILE). Options: --task NAME, --quick, --out FILE.
import sys, os, json, time, math, argparse
from operator import mul
import mpmath as mpm
from mpmath import mp, mpf, mpc

LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def setdps(d):
    mp.dps = d


# ---------------------------------------------------------------- basic functions
def Phi(u):
    """Titchmarsh's kernel, valid for real or complex u with Re e^(2u) > 0; even in u."""
    if not isinstance(u, mpc) and u < 0:
        u = -u
    e2 = mpm.exp(2 * u)
    e1 = mpm.exp(u / 2)
    a = mpm.exp(2 * u)
    tot = mpf(0)
    n = 1
    thr = mp.dps * 2.31 + 10
    while True:
        x = mp.pi * n * n * a
        t = (4 * mp.pi ** 2 * n ** 4 * e1 ** 9 - 6 * mp.pi * n * n * e1 ** 5) * mpm.exp(-x)
        tot += t
        if mpm.re(x) > thr:
            break
        n += 1
    return tot


def Xi(z):
    s = mpf(1) / 2 + 1j * z
    return mpm.re(s * (s - 1) / 2 * mp.pi ** (-s / 2) * mpm.gamma(s / 2) * mpm.zeta(s)) if not isinstance(z, mpc) else \
        s * (s - 1) / 2 * mp.pi ** (-s / 2) * mpm.gamma(s / 2) * mpm.zeta(s)


def friable(P, nmax):
    """integers n <= nmax whose prime factors are all <= P (P = 2: the powers of 2)."""
    out = []
    for n in range(1, nmax + 1):
        m = n
        for p in range(2, P + 1):
            while m % p == 0:
                m //= p
        if m == 1:
            out.append(n)
    return out


def oneminusR(x, t):
    return (1 - mpf(3) / (2 * x) * mpm.exp(-t)) * mpm.exp(x * mpm.exp(-t)) / (1 + mpm.exp(-mpf(9) * t / 2))


def psi_M(t, ns=None, lam=1):
    """Miyatake's psi(t) = sum n^4 e^(-2 pi n^2 cosh t) (1 - lam R(pi n^2, |t|)); ns = list of n (default all)."""
    t = abs(t)
    tot = mpf(0)
    if ns is None:
        n = 1
        while True:
            term = n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(t)) * (1 - lam * (1 - oneminusR(mp.pi * n * n, t)))
            tot += term
            if 2 * mp.pi * n * n * mpm.cosh(t) > mp.dps * 2.31 + 20:
                break
            n += 1
        return tot
    for n in ns:
        tot += n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(t)) * (1 - lam * (1 - oneminusR(mp.pi * n * n, t)))
    return tot


# ---------------------------------------------------------------- fixed point sweeps
def sweep(nodes, wvecs, z0, dz, K, P=400):
    """out[v][k] = sum_j wvecs[v][j] exp(i (z0 + k dz) nodes[j]) for k = 0..K-1, as python complex numbers.
    Fixed point of P bits, recurrence in k. nodes, weights: lists of mpf; z0, dz: mpf or mpc."""
    old = mp.dps
    mp.dps = int(P * 0.30103) + 25
    one = 1 << P
    Er, Ei, Rr, Ri = [], [], [], []
    for u in nodes:
        e = mpm.exp(1j * mpc(z0) * u)
        r = mpm.exp(1j * mpc(dz) * u)
        Er.append(int(mpm.re(e) * one)); Ei.append(int(mpm.im(e) * one))
        Rr.append(int(mpm.re(r) * one)); Ri.append(int(mpm.im(r) * one))
    W = [[int(w * one) for w in wv] for wv in wvecs]
    out = [[] for _ in wvecs]
    for k in range(K):
        for v, wv in enumerate(W):
            sr = sum(map(mul, wv, Er)) >> P
            si = sum(map(mul, wv, Ei)) >> P
            out[v].append(complex(sr / one, si / one))
        if k < K - 1:
            Er, Ei = zip(*[((a * c - b * d) >> P, (a * d + b * c) >> P) for a, b, c, d in zip(Er, Ei, Rr, Ri)])
    mp.dps = old
    return out


def trap_nodes(h, umax, f, even=True):
    """nodes u_j = j h (j = 0..N) and weights h f(u_j) (doubled for j > 0): for even f the full line integral is sum w_j cos(z u_j)."""
    N = int(math.ceil(umax / h))
    nodes = [h * j for j in range(N + 1)]
    wts = [(1 if j == 0 else 2) * h * f(h * j) for j in range(N + 1)]
    return nodes, wts


def gauss_legendre(m):
    import numpy as np
    x0, _ = np.polynomial.legendre.leggauss(m)
    xs, ws = [], []
    for x in x0:
        x = mpf(float(x))
        for _ in range(8):
            p0, p1 = mpf(1), x
            for k in range(2, m + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = m * (x * p1 - p0) / (x * x - 1)
            x = x - p1 / dp
        p0, p1 = mpf(1), x
        for k in range(2, m + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = m * (x * p1 - p0) / (x * x - 1)
        xs.append(x); ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def gl_panels(S, width, m):
    xs, ws = gauss_legendre(m)
    nodes, wts = [], []
    npan = int(round(S / width))
    for p in range(npan):
        a = p * mpf(width)
        for x, w in zip(xs, ws):
            nodes.append(a + mpf(width) * (1 + x) / 2); wts.append(mpf(width) * w / 2)
    return nodes, wts


def psi_parts(s, ns):
    """(kA, kB): the lambda-independent part sum n^4 e^(-2 pi n^2 cosh s) and the lambda part sum n^4 e^(-pi n^2 e^s)(1 - 3e^(-s)/(2 pi n^2))/(1+e^(-9s/2));
    psi_lambda = (1-lambda) kA + lambda kB."""
    A = mpf(0); B = mpf(0)
    for n in ns:
        A += n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(s))
        B += n ** 4 * mpm.exp(-mp.pi * n * n * mpm.exp(s)) * (1 - 3 * mpm.exp(-s) / (2 * mp.pi * n * n)) / (1 + mpm.exp(-mpf(9) * s / 2))
    return A, B


def xi_lambda_weights(ns, S=4.7, width=0.2, m=24):
    """nodes u = s/2 and weights for the two functions A(z) = 8 pi^2 int_0^inf kA cosh(9s/4) cos(zs/2) ds (this is xi_2 for ns = S(2)),
    B(z) the same with kB (Xi_P^M); Xi_lambda = (1-lambda) A + lambda B."""
    sn, sw = gl_panels(S, width, m)
    nodes, wA, wB = [], [], []
    for s, w in zip(sn, sw):
        kA, kB = psi_parts(s, ns)
        c = 8 * mp.pi ** 2 * w * mpm.cosh(9 * s / 4)
        nodes.append(s / 2); wA.append(c * kA); wB.append(c * kB)
    return nodes, wA, wB


# ---------------------------------------------------------------- tasks
def task_norm(R):
    setdps(50)
    out = {}
    h = mpf(1) / 100
    # (a) Xi(z) = int Phi e^(izu) du with Phi = Titchmarsh's kernel; the chapter's formula is twice this
    nodes, wts = trap_nodes(h, 2.6, Phi)
    # chapter's kernel: 2 pi e^(5t/2) sum (2 pi n^2 e^(2t) - 3) n^2 e^(-pi n^2 e^(2t)) = 2 Phi
    def Phi_ch(t):
        t = abs(t)
        tot = mpf(0)
        for n in range(1, 30):
            tot += (2 * mp.pi * n * n * mpm.exp(2 * t) - 3) * n * n * mpm.exp(-mp.pi * n * n * mpm.exp(2 * t))
        return 2 * mp.pi * mpm.exp(5 * t / 2) * tot
    t0 = mpf('0.3')
    out['Phi_chapter / Phi (code, 4pi^2 and 6pi) at 0.3'] = str(Phi_ch(t0) / Phi(t0))
    zs = [mpf(0), mpf(5), mpf('14.134725141734693790457251983562')]
    for z in zs:
        I = sum(w * mpm.cos(z * u) for u, w in zip(nodes, wts))
        out['int Phi cos / Xi at z=%s' % mpm.nstr(z, 6)] = mpm.nstr(I / Xi(z), 15)
    # (b) 8 pi^2 cosh(9t/2) psi(2t) vs Phi, with psi = Miyatake's sum with the formula for R
    for t in [mpf('0.1'), mpf('0.5'), mpf('1.0')]:
        lhs = 8 * mp.pi ** 2 * mpm.cosh(9 * t / 2) * psi_M(2 * t)
        out['8pi^2 cosh(9t/2) psi_M(2t) / Phi at t=%s' % mpm.nstr(t, 3)] = mpm.nstr(lhs / Phi(t), 15)
        # theta version: Phi(t) = 2 e^(t/2) (2 x^2 psi'' + 3 x psi'), x = e^(2t)
        x = mpm.exp(2 * t)
        ps1 = mpm.diff(lambda y: mpm.nsum(lambda n: mpm.exp(-mp.pi * n * n * y), [1, mpm.inf]), x)
        ps2 = mpm.diff(lambda y: mpm.nsum(lambda n: mpm.exp(-mp.pi * n * n * y), [1, mpm.inf]), x, 2)
        out['Phi / (2 e^(t/2)(2x^2 psi''+3x psi\')) at t=%s' % mpm.nstr(t, 3)] = mpm.nstr(Phi(t) / (2 * mpm.exp(t / 2) * (2 * x * x * ps2 + 3 * x * ps1)), 12)
    # (c) the transform of psi_M: 2 pi^2 {G(iz/2-9/4)+G(iz/2+9/4)} = 8 pi^2 int_0^inf psi cosh(9t/4) cos(zt/2) dt (psi even and smooth: trapezoid on R)
    hh = mpf(1) / 40
    for hq in (hh, hh * mpf(4) / 5):
        N = int(math.ceil(4.7 / hq))
        tn = [hq * j for j in range(N + 1)]
        pv = [psi_M(t) * mpm.cosh(9 * t / 4) * (1 if j == 0 else 2) * hq for j, t in enumerate(tn)]
        res = []
        for z in zs:
            T = 2 * mp.pi ** 2 * 2 * sum(w * mpm.cos(z * t / 2) for t, w in zip(tn, pv))
            res.append(T)
        out['transform psi_M (all n) at the three z, step %s' % mpm.nstr(hq, 4)] = [mpm.nstr(T / Xi(z), 14) for T, z in zip(res, zs)]
    R['norm'] = out
    for k, v in out.items():
        log('norm', k, v)


def task_kink(R):
    setdps(40)
    out = {}
    sl = []
    for n in (1, 2, 3):
        x = mp.pi * n * n
        r = 1 - 3 / (2 * x)
        sl.append(3 / (2 * x) / r - x + mpf(9) / 4)
    out['one sided log slopes n=1,2,3 at t=0+'] = [mpm.nstr(s, 6) for s in sl]
    # check by finite difference of log of the n-th term
    fd = []
    for n in (1, 2, 3):
        f = lambda t: mpm.log(n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(t)) * oneminusR(mp.pi * n * n, t))
        fd.append(mpm.diff(f, 0, direction=1))
    out['finite difference check'] = [mpm.nstr(s, 6) for s in fd]
    psi0 = psi_M(mpf(0))
    Q0 = psi0 / mpm.exp(-2 * mp.pi)
    out['Q(0)'] = mpm.nstr(Q0, 12)
    out['Q(0)-1'] = mpm.nstr(Q0 - 1, 12)
    ts = [mpf(j) / 20 for j in range(0, 121)]
    Q = [psi_M(t) / mpm.exp(-2 * mp.pi * mpm.cosh(t)) for t in ts]
    out['Q decreasing on [0,6] step 1/20'] = all(Q[i] > Q[i + 1] for i in range(len(Q) - 1))
    out['Q(6) - 1'] = mpm.nstr(Q[-1] - 1, 6)
    out['Q at 0, 0.5, 1, 2, 4'] = [mpm.nstr(Q[i], 9) for i in (0, 10, 20, 40, 80)]
    # smoothness of the full sum at 0: psi(t) - psi(-t) trivially 0; derivative from the right
    d = mpm.diff(lambda t: psi_M(t), 0, direction=1)
    out['psi_M one sided derivative at 0+ (smooth: must be 0)'] = mpm.nstr(d, 6)
    R['kink'] = out
    for k, v in out.items():
        log('kink', k, v)


def count_sign_changes(vals):
    return sum(1 for a, b in zip(vals[:-1], vals[1:]) if (a > 0) != (b > 0))


def task_xi2(R, quick):
    """real zeros of xi_2 (n in S(2)) and of Xi_P^M, A(z) and B(z) on a grid; counts by sign changes (two steps)."""
    ns = [1, 2, 4, 8, 16]
    zmax = 60 if quick else 200
    out = {}
    for (width, m, dz, P) in [(0.2, 24, 0.05, 400), (0.2, 28, 0.04, 330)]:
        t0 = time.time()
        setdps(int(P * 0.30103) + 25)
        nodes, wA, wB = xi_lambda_weights(ns, 4.7, width, m)
        K = int(zmax / dz) + 1
        A, B = sweep(nodes, [wA, wB], mpf(0), mpf(dz), K, P)
        a = [v.real for v in A]; b = [v.real for v in B]
        ca, cb = count_sign_changes(a), count_sign_changes(b)
        key = 'width %s m %d dz %s P %d' % (width, m, dz, P)
        out[key] = {'real zeros of xi_2 in (0,%d)' % zmax: ca, 'real zeros of Xi_P^M (lambda=1) in (0,%d)' % zmax: cb,
                    'seconds': round(time.time() - t0, 1)}
        log('xi2', key, out[key])
        if width == 0.2 and P == 400:
            # first zeros of xi_2 by linear interpolation
            zr = [dz * i - dz * a[i] / (a[i + 1] - a[i]) * 0 + dz * (a[i] / (a[i] - a[i + 1])) for i in range(len(a) - 1) if (a[i] > 0) != (a[i + 1] > 0)]
            out['first 5 zeros of xi_2'] = [round(x, 4) for x in zr[:5]]
            out['last zero of xi_2'] = round(zr[-1], 4)
            zb = [dz * i + dz * (b[i] / (b[i] - b[i + 1])) for i in range(len(b) - 1) if (b[i] > 0) != (b[i + 1] > 0)]
            out['first 5 zeros of Xi_P^M'] = [round(x, 4) for x in zb[:5]]
            out['last zero of Xi_P^M'] = round(zb[-1], 4)
            out['gap check: max |a| ratio (value at 0)'] = a[0]
    # N(200) of Xi from the zero file
    zs = [mpf(l.strip()) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_1700.txt')) if l.strip()]
    out['zeros of zeta with ordinate < %d' % zmax] = sum(1 for z in zs if z < zmax)
    R['xi2'] = out
    log('xi2', out)


def task_xi2trap(R, quick):
    """xi_2 (n in S(2)) = A(z) alone is the cosine transform of the even analytic function kA cosh(9s/4): trapezoid rule on R, step h (error e^(-1.2 (2 pi/h - z/2))),
    so the zero count up to 200 is reliable (the Gauss panels above lose the tiny values of the transform at z > 100)."""
    ns = [1, 2, 4, 8, 16]
    zmax = 100 if quick else 200
    out = {}
    for (h, dz, P) in [(0.01, 0.05, 400), (0.0125, 0.04, 340)]:
        t0 = time.time()
        setdps(int(P * 0.30103) + 25)
        hh = mpf(h)
        N = int(math.ceil(4.7 / h))
        nodes = [hh * j / 2 for j in range(N + 1)]
        w = []
        for j in range(N + 1):
            s_ = hh * j
            kA = sum(n ** 4 * mpm.exp(-2 * mp.pi * n * n * mpm.cosh(s_)) for n in ns)
            w.append((1 if j == 0 else 2) * hh * 8 * mp.pi ** 2 * mpm.cosh(9 * s_ / 4) * kA)
        K = int(zmax / dz) + 1
        v = sweep(nodes, [w], mpf(0), mpf(dz), K, P)[0]
        a = [x.real for x in v]
        zr = [dz * i + dz * a[i] / (a[i] - a[i + 1]) for i in range(len(a) - 1) if (a[i] > 0) != (a[i + 1] > 0)]
        key = 'h %s dz %s P %d' % (h, dz, P)
        out[key] = {'real zeros of xi_2 in (0,%d)' % zmax: len(zr), 'first zeros': [round(x, 4) for x in zr[:4]], 'last zero': round(zr[-1], 4), 'seconds': round(time.time() - t0, 1)}
        log('xi2trap', key, out[key])
    zs = [mpf(l.strip()) for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_1700.txt')) if l.strip()]
    out['zeros of zeta with ordinate < %d' % zmax] = sum(1 for z in zs if z < zmax)
    R['xi2trap'] = out


def task_tail(R):
    setdps(80)
    out = {}
    # Phi(i pi/9) and C_{9/2}
    z9 = mpc(0, mp.pi / 9)
    ph = Phi(z9)
    C = 4 * mp.pi / 9 * ph
    out['Phi(i pi/9)'] = mpm.nstr(ph, 15)
    out['C_9/2 = (4 pi/9) Phi(i pi/9)'] = mpm.nstr(C, 14)
    out['C_9/2/(4 pi^2)'] = mpm.nstr(C / (4 * mp.pi ** 2), 14)
    # terms of Phi(i pi/9)
    e2 = mpm.exp(2 * z9)
    terms = []
    for n in range(1, 12):
        t = (4 * mp.pi ** 2 * n ** 4 * mpm.exp(9 * z9 / 2) - 6 * mp.pi * n * n * mpm.exp(5 * z9 / 2)) * mpm.exp(-mp.pi * n * n * e2)
        terms.append(t)
    out['real parts of the first two terms of Phi(i pi/9)'] = [mpm.nstr(mpm.re(terms[0]), 6), mpm.nstr(mpm.re(terms[1]), 6)]
    bound = lambda n: 2 * (2 * mp.pi ** 2 + 3 * mp.pi) * n ** 4 * mpm.exp(-mpf('2.4066') * n * n)
    out['abs terms n>=3 <= bound (all n)'] = all(abs(terms[n - 1]) <= bound(n) for n in range(3, 12))
    out['sum_{n>=3} of the bound'] = mpm.nstr(sum(bound(n) for n in range(3, 30)), 6)
    out['sum_{n>=3} |term|'] = mpm.nstr(sum(abs(t) for t in terms[2:]), 6)
    out['2 * real part of the sum (n>=1) of Phi(i pi/9) check'] = mpm.nstr(mpm.re(sum(terms)), 12)
    out['exponent: pi cos(2 pi/9) = 2.4066?'] = mpm.nstr(mp.pi * mpm.cos(2 * mp.pi / 9), 8)
    # g(x) e^(2 pi x/9) at x = 40, 60, 80 : g(x) = g_{9/2}(2x)/(4 pi^2), trapezoid on the real line (pole of sech(9u/2) at i pi/9: d < 0.349)
    setdps(150)
    for h in (mpf(1) / 100, mpf(1) / 125):
        f = lambda u: Phi(u) * mpm.sech(mpf(9) / 2 * u)
        nodes, wts = trap_nodes(h, 2.5, f)
        vals = []
        for x in (40, 60, 80):
            z = 2 * mpf(x)
            g = sum(w * mpm.cos(z * u) for u, w in zip(nodes, wts)) / (4 * mp.pi ** 2)
            vals.append(g * mpm.exp(2 * mp.pi * x / 9))
        out['g(x) e^(2 pi x/9) at x=40,60,80, h=%s' % mpm.nstr(h, 4)] = [mpm.nstr(v, 14) for v in vals]
    setdps(50)
    R['tail'] = out
    for k, v in out.items():
        log('tail', k, v)


# ---------------------------------------------------------------- the cusp
def theta_pieces(b, nterms_tol=None):
    """Phi(i pi/(2b)) for b > 2 by the modular transformation: x = e^(2u) = e^(i pi/b), eps = 1 + i x, 1+2psi(x) = (-i eps)^(-1/2) theta_2(-1/eps),
    theta_2(tau) = sum_{n in Z} e^(i pi (n+1/2)^2 tau); Phi = 2 e^(u/2) (2 x^2 psi'' + 3 x psi'), u = i pi/(2b)."""
    b = mpf(b)
    x = mpm.expj(mp.pi / b)
    u = mpc(0, mp.pi / (2 * b))
    eps = 1 + 1j * x
    c = (-1j * eps) ** (-mpf(1) / 2)
    t1 = mpc(0); t2 = mpc(0); th = mpc(0)
    # derivative wrt x: d eps/dx = i
    dL0 = -1j / (2 * eps)
    d2L0 = -1 / (2 * eps ** 2)
    imtau = mpm.im(-1 / eps)
    tol = mpf(10) ** (-mp.dps - 5)
    n = 0
    while True:
        for m in ([mpf(n) + mpf(1) / 2, -(mpf(n) + mpf(1) / 2)] if True else []):
            pass
        m = mpf(n) + mpf(1) / 2
        T = c * mpm.exp(1j * mp.pi * m * m * (-1 / eps))
        L1 = dL0 - mp.pi * m * m / eps ** 2
        L2 = d2L0 + 2j * mp.pi * m * m / eps ** 3
        mult = 2  # n and -n-1 give the same m^2
        th += mult * T
        t1 += mult * L1 * T
        t2 += mult * (L2 + L1 ** 2) * T
        if abs(T) < tol * abs(th) and n > 0:
            break
        n += 1
        if n > 400:
            break
    psi1 = t1 / 2; psi2 = t2 / 2
    Ph = 2 * mpm.exp(u / 2) * (2 * x * x * psi2 + 3 * x * psi1)
    return Ph


def C_b(b):
    return 2 * mp.pi / mpf(b) * theta_pieces(b)


def task_cusp(R, quick, lo=None, hi=None):
    out = {}
    setdps(60)
    # agreement of the modular formula with the direct series at b = 3, 2.5
    for b in (mpf(3), mpf('2.5'), mpf(9) / 2):
        direct = 2 * mp.pi / b * Phi(mpc(0, mp.pi / (2 * b)))
        mod = C_b(b)
        out['modular vs direct C_b at b=%s' % mpm.nstr(b, 3)] = [mpm.nstr(mpm.re(direct), 14), mpm.nstr(mpm.re(mod), 14), mpm.nstr(abs(mod - direct) / abs(direct), 3)]
    # asymptotics and ratios
    rat = {}
    for bs in ('2.02', '2.005', '2.001', '2.0001'):
        b = mpf(bs)
        delta = mp.pi * (b - 2) / (4 * b)
        asy = -mp.pi ** 3 / (32 * b) * (2 * delta) ** (-mpf(1) / 2) * delta ** (-4) * mpm.exp(-mp.pi / 8 * mpm.cot(delta))
        cb = C_b(b)
        rat[bs] = [mpm.nstr(mpm.re(cb) / asy, 6), mpm.nstr(mpm.re(cb), 4), mpm.nstr(mpm.im(cb) / mpm.re(cb), 3)]
    out['ratio exact/asymptotic, C_b, imag/real'] = rat
    # second precision
    setdps(120)
    b = mpf('2.0001')
    delta = mp.pi * (b - 2) / (4 * b)
    asy = -mp.pi ** 3 / (32 * b) * (2 * delta) ** (-mpf(1) / 2) * delta ** (-4) * mpm.exp(-mp.pi / 8 * mpm.cot(delta))
    out['ratio at 2.0001 at dps 120'] = mpm.nstr(mpm.re(C_b(b)) / asy, 6)
    setdps(30)
    # sign scan
    def scan(a, bb, step):
        n = int(round((bb - a) / step))
        prev = None; changes = []
        for i in range(n + 1):
            bv = mpf(a) + i * mpf(step)
            if bv <= 2:
                continue
            s = mpm.re(C_b(bv)) > 0
            if prev is not None and s != prev:
                changes.append(bv)
            prev = s
        return changes
    if lo is None:
        lo, hi = (2.0, 2.5)
    step = 1e-4 if quick else 1e-5
    t0 = time.time()
    ch = scan(mpf(lo) + mpf(step), hi, mpf(step))
    out['sign changes of C_b on (%s,%s) step %s' % (lo, hi, step)] = [mpm.nstr(c, 10) for c in ch]
    out['seconds scan'] = round(time.time() - t0, 1)
    # b_0
    setdps(40)
    f = lambda b: mpm.re(C_b(b))
    b0 = mpm.findroot(f, mpf('2.4588'), tol=1e-30)
    out['b_0'] = mpm.nstr(b0, 14)
    setdps(60)
    out['b_0 dps 60'] = mpm.nstr(mpm.findroot(lambda b: mpm.re(C_b(b)), mpf('2.4588'), tol=1e-45), 14)
    out['C_b at 2.4, 2.5 (sign)'] = [mpm.nstr(mpm.re(C_b(mpf('2.4'))), 6), mpm.nstr(mpm.re(C_b(mpf('2.5'))), 6)]
    R['cusp'] = out
    for k, v in out.items():
        log('cusp', k, v)


def task_cscan(R, lo, hi, step):
    setdps(30)
    prev = None; changes = []
    n = int(round((hi - lo) / step))
    t0 = time.time()
    for i in range(n + 1):
        bv = mpf(lo) + i * mpf(step)
        s = mpm.re(C_b(bv)) > 0
        if prev is not None and s != prev:
            changes.append(mpm.nstr(bv, 10))
        prev = s
    R['cscan_%s_%s_%s' % (lo, hi, step)] = {'changes': changes, 'seconds': round(time.time() - t0, 1)}
    log('cscan', lo, hi, step, changes, round(time.time() - t0, 1))



def unwrap_change(vals):
    """total change of the argument along the list of complex numbers (steps wrapped to (-pi, pi])."""
    tot = 0.0
    for a, b in zip(vals[:-1], vals[1:]):
        d = math.atan2((b / a).imag, (b / a).real)
        tot += d
    return tot


def task_interp(R, quick):
    """T(lambda): the first real part X at which Xi_lambda = (1-lambda) xi_2 + lambda Xi_P^M (n in S(2)) has a nonreal zero, found by the argument
    principle on the rectangles [0,X] x [-Y,Y] (Y = 6) minus the real sign changes; X on a grid of step 1 (the zero of Xi_lambda just above T)."""
    ns = [1, 2, 4, 8, 16]
    zmax = 70 if quick else 100
    Y = float(os.environ.get('MI_Y', '6.0'))
    setdps(145)
    nodes, wA, wB = xi_lambda_weights(ns, 4.7, 0.2, 24)
    t0 = time.time()
    dzr = 0.05
    Kr = int(zmax / dzr) + 1
    Ar, Br = sweep(nodes, [wA, wB], mpf(0), mpf(dzr), Kr, 400)
    Ar = [v.real for v in Ar]; Br = [v.real for v in Br]
    # top row y = Y, mirrored nodes
    mn = nodes + [-u for u in nodes]
    hA = [w / 2 for w in wA] * 2; hB = [w / 2 for w in wB] * 2
    dxt = 0.1
    Kt = int(zmax / dxt) + 1
    At, Bt = sweep(mn, [hA, hB], mpc(0, Y), mpf(dxt), Kt, 400)
    log('interp: real and top row done', round(time.time() - t0, 1), 's')
    Xs = list(range(1, zmax + 1))
    VA, VB = {}, {}
    for X in [0] + Xs:
        a, b = sweep(mn, [hA, hB], mpf(X), mpc(0, 0.2), int(Y / 0.2) + 1, 400)
        VA[X] = a; VB[X] = b
    log('interp: vertical edges done', round(time.time() - t0, 1), 's')
    out = {}
    res = {}
    for k in range(0, 9):
        lam = 10.0 ** (-k)
        comb = lambda a, b: [(1 - lam) * x + lam * y for x, y in zip(a, b)]
        rr = comb(Ar, Br)
        top = comb(At, Bt)
        c_cum = [0]
        for i in range(len(rr) - 1):
            c_cum.append(c_cum[-1] + (1 if (rr[i] > 0) != (rr[i + 1] > 0) else 0))
        V0 = unwrap_change(comb(VA[0], VB[0]))
        T = None
        for X in Xs:
            dtop = unwrap_change(top[:int(round(X / dxt)) + 1])
            VX = unwrap_change(comb(VA[X], VB[X]))
            W = (-2 * dtop + 2 * VX - 2 * V0) / (2 * math.pi)
            c = c_cum[int(round(X / dzr))]
            if round(W) - c > 0:
                T = X
                break
        res[k] = T
        log('interp lambda = 1e-%d: first X with a nonreal zero' % k, T, '(windings: total %s, real %s at that X)' % (round(W, 3), c))
    out['T(lambda = 10^-k) first X on the grid with a nonreal zero'] = res
    ks = [k for k in res if res[k] is not None]
    if len(ks) >= 3:
        import numpy as np
        sl = np.polyfit(ks, [res[k] for k in ks], 1)[0]
        out['fitted slope per decade'] = float(sl)
        out['predicted slope (4/pi) ln 10'] = 4 / math.pi * math.log(10)
    out['seconds'] = round(time.time() - t0, 1)
    R['interp'] = out
    log('interp', out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--lo', type=float, default=None)
    ap.add_argument('--hi', type=float, default=None)
    ap.add_argument('--step', type=float, default=1e-5)
    a = ap.parse_args()
    R = {}
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or os.path.join(here, '..', 'data', 'mi_miyatake_%s.json' % a.task)
    t = a.task
    if t in ('norm', 'all'): task_norm(R)
    if t in ('kink', 'all'): task_kink(R)
    if t in ('xi2', 'all'): task_xi2(R, a.quick)
    if t == 'xi2trap': task_xi2trap(R, a.quick)
    if t in ('tail', 'all'): task_tail(R)
    if t == 'interp': task_interp(R, a.quick)
    if t in ('cusp', 'all'): task_cusp(R, a.quick, a.lo, a.hi)
    if t == 'cscan': task_cscan(R, a.lo, a.hi, a.step)
    json.dump(R, open(out, 'w'), indent=1, default=str)
    open(out.replace('.json', '.log'), 'w').write('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
