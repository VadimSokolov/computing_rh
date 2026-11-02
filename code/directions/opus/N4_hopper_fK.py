#!/usr/bin/env python3
"""Hopper proposal for N4, not run at full size. Usage: python3 -I N4_hopper_fK.py Kmax out.txt
Zeros by Arb, then v1(K)/2beta_K (double precision), a one point Arb evaluation of log|phi_K| at the first K where the ratio passes 1, and the reflection certificate of S5."""
import sys, math, time
sys.dont_write_bytecode = True
sys.path.insert(0, '/tmp/d1lib')
import numpy as np
from flint import arb, acb, acb_series, ctx
t0 = time.time()
Kmax = int(sys.argv[1])
out = open(sys.argv[2], 'w')
def say(*a):
    print(*a, file=out, flush=True)
ctx.prec = 128
G = []
for start in range(1, Kmax + 1, 2000):
    G += [z.imag for z in acb.zeta_zeros(start, min(2000, Kmax - start + 1))]
    say('zeros', len(G), 'time', round(time.time() - t0, 1))
F = np.array([float(g.mid()) for g in G])
z = acb_series([acb(0.5), acb(1)], 3).zeta()
a0, a1, a2 = z[0], z[1], 2 * z[2]
logz2 = a2 / a0 - (a1 / a0) ** 2
trig = arb.pi() ** 2 + 8 * arb.const_catalan()
S = ((-8 + trig / 4 + logz2.real) / 2)
say('sum gamma^-2 (ball)', S.str(20))
Sf = float(S.mid())
inv2 = 1 / F ** 2
cum = np.cumsum(inv2)
rs = np.linspace(1.02, 1.5, 241)
grid = sorted(set(int(k) for k in np.geomspace(min(1000, Kmax), Kmax, 40)))
first = None
for K in grid:
    T = F[K - 1]
    betaK = Sf - cum[K - 1]
    th = rs * T
    vals = np.array([np.sum(np.log(np.abs(1 - t * t * inv2[:K]))) for t in th])
    sup = vals / th ** 2
    i = int(np.argmax(sup))
    ratio = sup[i] / betaK
    say('K', K, 'T %.3f' % T, 'ratio v1/2beta_K %.5f' % ratio, 'r* %.3f' % rs[i], 'log|phi_K| %.2f' % (vals[i] - betaK * th[i] ** 2))
    if ratio > 1 and first is None:
        first = (K, float(th[i]))
label = 'first K with ratio>1:'
if first is None:
    first = (grid[-1], float(th[i]))
    label = 'no K with ratio>1; test evaluation at the last K:'
if first is not None:
    K, t = first
    t = arb(t)
    betaK = S - sum((1 / (g * g) for g in G[:K]), arb(0))
    lg = -betaK * t * t
    for g in G[:K]:
        lg = lg + abs(1 - t * t / (g * g)).log()
    say(label, K, 'theta', float(t.mid()), 'log|phi_K| ball', lg.str(12), 'certified positive', lg > 0)
g1 = F[0]
idx = np.nonzero(F - 2 * g1 >= g1 - 1e-6)[0]
t = F[idx] - 2 * g1
j = np.searchsorted(F, t)
d = np.minimum(np.abs(F[np.clip(j, 0, len(F) - 1)] - t), np.abs(F[np.clip(j - 1, 0, len(F) - 1)] - t))
m = int(np.argmin(d))
say('reflection test: tested', len(idx), 'smallest |gamma_a-gamma_c-2gamma_1| %.3e' % d[m], 'a', int(idx[m]) + 1, 'K0 reach', int((F[-1] / g1 - 1) / 2))
say('time', round(time.time() - t0, 1))
