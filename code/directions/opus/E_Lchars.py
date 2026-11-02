#!/usr/bin/env python3
# Track E. GHS and the EMN condition across real primitive Dirichlet L functions, against zeta and DH (control).
# For each law: density positive on a grid?, sign of V''''(0) (EMN needs V'''>=0 near 0), and
# kappa_3(h) = (d/dh)^3 log Lambda(1/2+h) on a grid (GHS needs <=0 for h>0). Values at grid points are balls.
# Run: timeout 500 python3 -I E_Lchars.py
import sys, math
sys.path.insert(0, '/tmp/e_lib')
from flint import arb, acb, arb_series, acb_series, ctx
ctx.prec = 128
ctx.cap = 8
pi = arb.pi()
sq5 = arb(5).sqrt()
kappa = ((10 - 2 * sq5).sqrt() - 2) / (sq5 - 1)


def kron(D, n):
    def leg(D, p):
        if p == 2:
            if D % 2 == 0:
                return 0
            return 1 if D % 8 in (1, 7) else -1
        v = pow(D % p, (p - 1) // 2, p)
        return -1 if v == p - 1 else v
    r = 1
    m = n
    p = 2
    while m > 1:
        if p * p > m:
            p = m
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            r *= leg(D, p)**e
            if r == 0:
                return 0
        p += 1
    return r


def make(name, q, coef, odd):
    return dict(name=name, q=q, coef=[coef(j) for j in range(q + 1)], odd=odd)


LAWS = [make('DH q=5', 5, lambda n: [arb(0), arb(1), kappa, -kappa, arb(-1), arb(0)][n], True)]
for D in (-3, -4, -7, -8, -11, -15, -19, -20, -23, -24):
    LAWS.append(make('D=%d' % D, -D, lambda n, D=D: (kron(D, n) if n % (-D) and math.gcd(n, -D) == 1 else 0), True))
for D in (5, 8, 12, 13, 17, 21, 24):
    LAWS.append(make('D=%d' % D, D, lambda n, D=D: (kron(D, n) if n % D and math.gcd(n, D) == 1 else 0), False))


def logLam(law, c):
    """log Lambda(s) series at s=c+h, real c>=1/2 exact. L(s)=q^{-s} sum_j a_j zeta(s,j/q); no pole since sum a_j=0."""
    q = law['q']
    s = acb_series([acb(c), acb(1)], prec=6)
    tot = None
    for j in range(1, q):
        a = law['coef'][j]
        if a == 0:
            continue
        z = s.zeta(acb(arb(j) / q)) * a
        tot = z if tot is None else tot + z
    L = tot * ((-s * arb(q).log()).exp())
    shift = 1 if law['odd'] else 0
    g = ((s + shift) * 0.5).lgamma() + s * (arb(q) / pi).log() * 0.5
    return g + L.log()


def logxi(c):
    s = acb_series([acb(c), acb(1)], prec=6)
    if c < 2:
        r = 1 - s
        return s.log() - r * pi.log() * 0.5 + ((3 - s) * 0.5).lgamma() + ((-1) * r.zeta()).log()
    return s.log() + (s - 1).log() + (s * 0.5).lgamma() - s * pi.log() * 0.5 + s.zeta().log()


def dens_series(law, t0):
    q = law['q']
    T = arb_series([arb(t0), arb(1)], prec=6)
    w = (2 * T).exp()
    N = int(16 * math.sqrt(q)) + 10
    tot = None
    for n in range(1, N):
        a = law['coef'][n % q]
        if a == 0:
            continue
        if law['odd']:
            term = a * n * ((-pi * n**2 * w * (arb(1) / q)).exp())
        else:
            term = a * ((-pi * n**2 * w * (arb(1) / q)).exp())
        tot = term if tot is None else tot + term
    return tot * ((T * (arb(3) / 2 if law['odd'] else arb(1) / 2)).exp())


def report(name, kfun, dens):
    grid = [arb(j) / 16 for j in range(1, 241)]
    mx = None
    arg = None
    for h in grid:
        v = kfun(arb(0.5) + h)
        v3 = v[3] * 6
        u = float(v3.real.mid())
        if mx is None or u > mx[0]:
            mx = (u, float(h), v3.real)
    ghs = 'pass' if mx[2].upper() < 0 else ('FAIL' if mx[2].lower() > 0 else 'undecided')
    s = name.ljust(10) + ' GHS on grid (0,15]: %-4s max kappa_3 %s at h=%.4f' % (ghs, mx[2].str(5), mx[1])
    if dens is not None:
        mn = min(float(dens(j / 32.)[0].mid()) for j in range(0, 97))
        d = dens(0.0)
        L = d.log()
        v2 = -2 * L[2]
        v4 = -24 * L[4]
        s += " | min density on [0,3] %.3e; V''(0)=%.4f V''''(0)=%.4f" % (mn, float(v2), float(v4))
    print(s, flush=True)


report('zeta', logxi, None)
for law in LAWS:
    report(law['name'], lambda c, law=law: logLam(law, c), lambda t, law=law: dens_series(law, t))
