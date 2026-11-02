"""Theorem gramexact repaired: analytic bound for t>=34 and ball arithmetic certificate on [0,g_4]. Run: timeout 300 python3 -I N6_c.py"""
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from N6_common import *
ctx.prec = 256
LOGPI = arb.pi().log()
def theta_b(t):
    t = acb(t)
    return (acb(arb(1) / 4) + I_ * t / 2).lgamma().imag - t.real * LOGPI / 2
def thetap_b(t):
    return (acb(arb(1) / 4, arb(t) / 2)).digamma().real / 2 - LOGPI / 2
def gram(k, x0):
    x = arb(x0)
    for _ in range(40):
        x = arb((x - (theta_b(x) - k * arb.pi()) / thetap_b(x)).mid())
    return x
G = {}
for k, x0 in ((-1, 9.667), (0, 17.846), (1, 23.17), (2, 27.67), (3, 31.718), (4, 35.467), (5, 39.0)):
    G[k] = gram(k, x0)
    xl, xr = float(G[k]) - 1e-12, float(G[k]) + 1e-12
    ok = ((theta_b(arb(xl)) - k * arb.pi()) < 0) and ((theta_b(arb(xr)) - k * arb.pi()) > 0)
    v = psi_ball(acb(arb(float(G[k]), 2e-12)), 1)
    print('g_%d=%.10f bracket ok %s; sign of psi(2g_k) on the ball: %s' % (k, float(G[k]), ok, '+' if v.real > 0 else ('-' if v.real < 0 else '?')))
f1 = lambda t: psi_ball(t, 1)
c = 0.5
pi_ = arb.pi()
def integrand(u, analytic):
    u = acb(u)
    return 4 * u * (2 * pi_ * pi_ * (9 * u).cosh() - 3 * pi_ * (5 * u).cosh()) * (-2 * pi_ * (4 * u).cosh()).exp() * (2 * c * u).cosh()
L = acb.integral(integrand, 0, 4).real
Lup = arb(float(L.upper()) * 1.0001 + 1e-20)
print('Lipschitz constant of psi(2t) on |Im t|<=%.2f: L<=%.6f' % (c, float(Lup)))
def walk(fun, p, q, Lup):
    """Rigorous arg change along [p,q]: from t_j step r_j=|v_j|/(2L), so psi stays in the disc |v-v_j|<=|v_j|/2 on the step."""
    total = arb(0); n = 0
    d = abs(q - p); u = (q - p) / d; pos = 0.0
    v = fun(acb(p.real, p.imag))
    while pos < d:
        r = float((abs(v) / (2 * Lup)).lower()) * 0.99
        if not r > 0:
            raise RuntimeError('zero near boundary')
        pos2 = min(d, pos + r)
        t2 = p + u * pos2
        v2 = fun(acb(t2.real, t2.imag))
        total = total + (v2 / v).arg()
        v, pos = v2, pos2; n += 1
    return total, n
a = float(G[-1])
total = arb(0); ns = 0
cs = [complex(-a, -c), complex(a, -c), complex(a, c), complex(-a, c), complex(-a, -c)]
for u_, v_ in zip(cs[:-1], cs[1:]):
    tt, n = walk(f1, u_, v_, Lup); total += tt; ns += n
print('winding [-g_-1,g_-1]x[-1/2,1/2] by the Lipschitz walk: %s, %d steps' % ((total / (2 * pi_)).str(8), ns))
w, ok, ns = cert_winding(f1, float(G[-1]), float(G[0]), -c, c, 0.05, 22)
print('winding [g_-1,g_0]x[-1/2,1/2]: %s (ok %s, %d pieces)' % (w.str(8) if w is not None else 'NA', ok, ns))
for k in range(0, 5):
    w, ok, ns = cert_winding(f1, float(G[k]), float(G[k + 1]), -c, c, 0.05, 22)
    print('winding [g_%d,g_%d]x[-1/2,1/2]: %s (ok %s, %d pieces)' % (k, k + 1, w.str(8) if w is not None else 'NA', ok, ns))
print('Analytic bound for real t')
pi = arb.pi(); PI2 = pi * pi
def bounds(t):
    y = t / 2.0
    A = float(PI2) / y
    D = math.hypot(0.25, y)
    Em = (y + 1.5) / D * (math.exp(A) - 1)
    dE = (0.5 * ((A - 1) * math.exp(A) + 1) + 0.75 * (A / y) * math.exp(A)) / D + 0.5 * Em / D
    return Em, dE
def refl_ball(t):
    t = acb(t); y = I_ * t / 2
    n1 = acb(arb(9) / 4) + y; n2 = acb(arb(5) / 4) + y
    z = acb(2 * pi)
    r = 2 * PI2 * (-pi / (2 * (pi * n1).sin())) * z.bessel_i(n1) - 3 * pi * (-pi / (2 * (pi * n2).sin())) * z.bessel_i(n2)
    return r / (gamma_inf_ball(acb(arb(1) / 2) + I_ * t) / 2)
for t in (32, 33, 34, 36, 40, 60, 100, 700):
    Em, dE = bounds(t)
    R = float(abs(refl_ball(acb(t))).upper())
    dsk = acb(arb(t, 1), arb(0, 1))
    Rd = float(abs(refl_ball(dsk)).upper())
    tot = Em + R
    tp = float(thetap_b(arb(t)))
    print(' t=%-4d |E_m|<=%.4f |R|<=%.2e sup_disc|R|<=%.2e  |E|<=%.4f |E_m\'|<=%.4f |delta\'|<=%s theta\'=%.4f' % (t, Em, R, Rd, tot, dE, ('%.4f' % ((dE + Rd) / (1 - tot))) if tot < 1 else 'n/a', tp))
