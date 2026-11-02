#!/usr/bin/env python3
# Track E. Ball arithmetic tests of the Griffiths, Hurst and Sherman inequality (GHS) and of the
# Ellis, Monroe and Newman condition (V' convex) for the Polya law and the Davenport and Heilbronn law.
# Run: timeout 300 python3 -I E_ghs_emn.py
import sys
sys.path.insert(0, '/tmp/e_lib')
from flint import arb, acb, arb_series, acb_series, ctx
ctx.prec = 128
ctx.cap = 8
pi = arb.pi()
sq5 = arb(5).sqrt()
kappa = ((10 - 2 * sq5).sqrt() - 2) / (sq5 - 1)
A = [1, kappa, -kappa, -1]
fac = [1, 1, 2, 6, 24]


def cell(a, b):
    """ball containing [a,b], a,b exact dyadic"""
    a = arb(a)
    b = arb(b)
    return (a + b) / 2 + (b - a) / 2 * arb(0, 1)


def ser(center):
    return acb_series([acb(center), acb(1)], prec=7)


def logxi_series(c):
    """log xi(s) as a series in h at s=c+h, c real. Direct for c>=2, reflected for c<2."""
    s = ser(c)
    if c.lower() >= 2:
        return s.log() + (s - 1).log() + (s * 0.5).lgamma() - s * pi.log() * 0.5 + s.zeta().log() - acb(2).log()
    # xi(s) = s pi^{-(1-s)/2} Gamma((3-s)/2) (-zeta(1-s))
    r = 1 - s
    return s.log() - r * pi.log() * 0.5 + ((3 - s) * 0.5).lgamma() + ((-1) * r.zeta()).log()


def f_dh(s):
    tot = None
    for j, a in enumerate(A):
        z = s.zeta(acb(arb(j + 1) / 5)) * a
        tot = z if tot is None else tot + z
    return tot * ((-s * arb(5).log()).exp())


def logLam_series(c):
    """Lambda(s) = (5/pi)^{s/2} Gamma((s+1)/2) f(s); uses Lambda(s)=Lambda(1-s) for c<1.5."""
    s = ser(c)
    if c.lower() >= 1.5:
        return s * (arb(5) / pi).log() * 0.5 + ((s + 1) * 0.5).lgamma() + f_dh(s).log()
    r = 1 - s
    return r * (arb(5) / pi).log() * 0.5 + ((r + 1) * 0.5).lgamma() + f_dh(r).log()


def deriv(series, k):
    v = series[k] * arb(1)
    for i in range(2, k + 1):
        v = v * i
    return v


def grid_ghs(name, logf, Hmax, n=240):
    """kappa_3(h) at exact grid points, each value a certified ball."""
    mx = None
    arg = None
    neg = 0
    for j in range(1, n + 1):
        h = arb(j) * arb(Hmax) / n
        v = deriv(logf(arb(0.5) + h), 3).real
        if v.upper() < 0:
            neg += 1
        u = float(v.upper())
        if mx is None or u > mx:
            mx = u
            arg = float(h)
    print(name, 'kappa_3 at %d grid points in (0,%g]: certified negative at %d; largest value %.6g at h=%.4f' % (n, Hmax, neg, mx, arg), flush=True)


def theta_P(tb):
    """Phi(t)=exp(-pi e^{2t}) G(t): returns (c,G) with c=pi. Tail n>=7: relative radius 1e-50 (see E.md)."""
    T = arb_series([tb, arb(1)], prec=6)
    w = (2 * T).exp()
    tot = None
    for n in range(1, 7):
        term = (2 * pi**2 * n**4 * (T * (arb(9) / 2)).exp() - 3 * pi * n**2 * (T * (arb(5) / 2)).exp()) * ((-pi * (n**2 - 1) * w).exp())
        tot = term if tot is None else tot + term
    return pi, tot


def theta_DH(tb):
    """e^{3t/2} theta(e^{2t}) = exp(-(pi/5) e^{2t}) G(t): returns (c,G) with c=pi/5."""
    T = arb_series([tb, arb(1)], prec=6)
    w = (2 * T).exp()
    tot = None
    for n in range(1, 21):
        r = n % 5
        if r == 0:
            continue
        a = [arb(0), arb(1), kappa, -kappa, arb(-1)][r]
        term = a * n * ((-pi * (n**2 - 1) * w * (arb(1) / 5)).exp())
        tot = term if tot is None else tot + term
    return pi / 5, tot * ((T * (arb(3) / 2)).exp())


def V_derivs(dens, tb):
    # V=c e^{2t} - log G: returns the first four derivatives
    c, G = dens(tb)
    L = G.log()
    e = arb(0, 1e-50)
    E = (2 * tb).exp()
    V1 = 2 * c * E - L[1] + e
    V2 = 4 * c * E - 2 * L[2] + e
    V3 = 8 * c * E - 6 * L[3] + e
    V4 = 16 * c * E - 24 * L[4] + e
    return V1, V2, V3, V4


def adaptive_pos(fn, a, b, minw=2.0**-14):
    v = fn(cell(a, b))
    if v.lower() > 0:
        return 1, 0, float(v.lower())
    if b - a < minw:
        print('  V undecided on [%.8f,%.8f]: %s' % (a, b, v.str(6)))
        return 1, 1, float(v.mid())
    m = (a + b) / 2
    c1, u1, w1 = adaptive_pos(fn, a, m, minw)
    c2, u2, w2 = adaptive_pos(fn, m, b, minw)
    return c1 + c2, u1 + u2, min(w1, w2)


def scan_emn(name, dens, T, t0=1 / 16.):
    V = V_derivs(dens, arb(0))
    print(name, "V''(0)=%s V''''(0)=%s" % (V[1].str(10), V[3].str(10)))
    c, u, w = adaptive_pos(lambda x: V_derivs(dens, x)[3], 0.0, t0)
    print(name, "V''''>0 on [0,%g]: cells %d undecided %d, min lower bound %.5g" % (t0, c, u, w))
    tot = [0, 0]
    mn = 1e99
    a = t0
    while a < T:
        b = min(a + 1 / 16., T)
        c, u, w = adaptive_pos(lambda x: V_derivs(dens, x)[2], a, b)
        tot[0] += c
        tot[1] += u
        mn = min(mn, w)
        a = b
    print(name, "V'''>0 on [%g,%g]: cells %d undecided %d, min lower bound %.5g" % (t0, T, tot[0], tot[1], mn), flush=True)


print('kappa=', kappa.str(12))
print('Polya: kappa_2(0), kappa_4(0), kappa_6(0):', [deriv(logxi_series(arb(0.5)), k).real.str(12) for k in (2, 4, 6)])
print('DH   : kappa_2(0), kappa_4(0), kappa_6(0):', [deriv(logLam_series(arb(0.5)), k).real.str(12) for k in (2, 4, 6)])
grid_ghs('Polya', logxi_series, 9.6)
grid_ghs('DH', logLam_series, 12.6)
scan_emn('Polya', theta_P, 1.5)
scan_emn('DH', theta_DH, 2.5)
