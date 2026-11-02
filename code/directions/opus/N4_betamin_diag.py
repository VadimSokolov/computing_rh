#!/usr/bin/env python3
"""N4: (a) why the extended beta_min table of Note 1 is wrong, tested on my own reconstruction of its sampling rule (root bound B, v=u^2 sampled at B*i/3000 and B*10^(-i/200));
(b) Arb certified zeros up to gamma_8000 and the two block certificates.
Run: timeout 580 python3 -I N4_betamin_diag.py > N4_betamin_diag.log. Arb for the signs and the zeros."""
import sys, os, time, bisect
sys.dont_write_bytecode = True
sys.path.insert(0, '/tmp/d1lib')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_poly, ctx
t0 = time.time()
ctx.prec = 6000
zs = acb.zeta_zeros(1, 320)
gam = [z.imag for z in zs]
print('zeros 320 at 6000 bits', round(time.time() - t0, 1), flush=True)
def coeffs(N, beta):
    s2 = 2 * arb(beta)
    P = [arb(1)]
    for j in range(N):
        a2 = (2 * gam[j]) ** 2
        new = [arb(0)] * (len(P) + 2)
        for k, c in enumerate(P):
            new[k] = new[k] + c - c / (a2 * s2)
            new[k + 2] = new[k + 2] + c / (a2 * s2 ** 2)
            if k >= 1:
                new[k] = new[k] - 2 * k * c / (a2 * s2)
            if k >= 2:
                new[k - 2] = new[k - 2] + k * (k - 1) * c / a2
        P = new
    return P[0::2]
cases = [(40, 5.50e-4), (60, 4.91e-4), (80, 4.47e-4), (120, 3.76e-4), (160, 2.09e-4), (200, 1.17e-4), (300, 0.450e-4)]
print('N beta_note log10B u_min_sampled(=sqrt(B*10^-5.995)) neg_on_sampling_rule u_first_neg(fine grid in u) Q there sign_certain')
for N, b in cases:
    c = coeffs(N, b)
    poly = arb_poly(c)
    ratios = [abs(x / c[-1]) for x in c[:-1]]
    B = arb(1)
    mx = ratios[0]
    for r in ratios[1:]:
        if float((r.log()).mid()) > float((mx.log()).mid()):
            mx = r
    B = 1 + mx
    lgB = float((B.log() / arb(10).log()).mid())
    vmin = B * arb(10) ** (-arb(1199) / 200)
    umin = vmin.sqrt()
    neg = 0
    amb = 0
    for i in range(0, 3001):
        v = B * arb(i) / 3000
        val = poly(v)
        if val < 0:
            neg += 1
        elif not val > 0:
            amb += 1
    for i in range(1, 1200):
        v = B * arb(10) ** (-arb(i) / 200)
        val = poly(v)
        if val < 0:
            neg += 1
        elif not val > 0:
            amb += 1
    first = None
    for i in range(0, 6001):
        u = arb(i) / 6000 * 3
        val = poly(u * u)
        if val < 0:
            first = (float(u.mid()), float(val.mid()), True)
            break
    print(N, b, round(lgB, 2), '%.3e' % float(umin.mid()), neg, 'ambiguous', amb, first, 'time', round(time.time() - t0, 1), flush=True)
print()
ctx.prec = 128
t1 = time.time()
all8 = acb.zeta_zeros(1, 8000)
G = [z.imag for z in all8]
print('Arb zeros 1..8000 at 128 bits', round(time.time() - t1, 1), 's; gamma_8000', G[-1].str(25), 'max radius', max(float(g.rad()) for g in G), flush=True)
F = [float(g.mid()) for g in G]
assert all(G[i] < G[i + 1] for i in range(len(G) - 1))
def closest(t):
    x = float(t.mid())
    i = bisect.bisect_left(F, x)
    best = None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(G):
            d = abs(G[j] - t)
            if best is None or float(d.mid()) < float(best[0].mid()):
                best = (d, j)
    return best
g1 = G[0]
res1 = []
for a in range(1, len(G)):
    t = 2 * G[a] - g1
    if t < G[-1]:
        d, j = closest(t)
        res1.append((float(d.mid()), a + 1, j + 1, d > arb(1) / 10 ** 6))
res1.sort()
print('three term test from gamma_1, a up to', max(r[1] for r in res1), 'tested', len(res1), 'smallest', res1[:3], 'all certified nonzero', all(r[3] for r in res1))
top = max(r[1] for r in res1)
print('K0 reach', int((F[top - 1] / F[0] - 1) / 2), 'gamma_a max', F[top - 1])
res2 = []
for a in range(len(G)):
    t = G[a] - 2 * g1
    if t >= g1 - arb(1) / 10 ** 6:
        d, j = closest(t)
        res2.append((float(d.mid()), a + 1, j + 1, d > arb(1) / 10 ** 6))
res2.sort()
print('difference 2 gamma_1 test tested', len(res2), 'smallest', res2[:3], 'all certified nonzero', all(r[3] for r in res2))
print('time', round(time.time() - t0, 1))
