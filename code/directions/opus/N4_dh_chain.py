#!/usr/bin/env python3
"""N4: the chain of Proposition synth(iii) run on the Davenport and Heilbronn law. Positive syntheses in height order at beta=Var/4, through the off line pairs, converge weakly to V_DH+sqrt(2 beta)Z, whose characteristic function has nonreal zeros.
Run: timeout 300 python3 -I N4_dh_chain.py > N4_dh_chain.log. Computed (double precision against mpmath), not certified."""
import sys, os, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
import N4_common as C
t0 = time.time()
mp.mp.dps = 20
VarDH = 0.049034788
beta = VarDH / 4
d = C.dh_json()
items = sorted([(float(s), 0.0) for s in d['online']] + [(float(t), float(s) - 0.5) for s, t in d['offline']])
h = 0.05
x = np.arange(0, 64 + h / 2, h)
L0 = C.dh_Lam(mp.mpf(1) / 2)
target = np.array([float(mp.re(C.dh_Lam(mp.mpf(1) / 2 + 0.5j * mp.mpf(float(v))) / L0)) for v in x]) * np.exp(-beta * x * x)
print('target char fn computed', round(time.time() - t0, 1), 'at x=64 %.2e' % abs(target[-1]), flush=True)
u = np.arange(0, 3 + 0.001, 0.002)
C_ = np.cos(np.outer(u, x))
w = np.full_like(x, h / np.pi); w[0] /= 2; w[-1] /= 2
inv = lambda F: C_ @ (w * F)
ft = inv(target)
print('target density: mass %.8f, min %.3e, Var %.6f (expected %.6f)' % (2 * np.trapezoid(ft, u), ft.min(), 2 * np.trapezoid(u * u * ft, u), VarDH + 2 * beta), flush=True)
F = np.exp(-beta * x * x)
K = 0; s = 0.0
want = {5, 43, 45, 49, 59, 130}
print('factors K, pairs included, budget Sum Re 1/4gamma^2 (needs <= 2 beta = %.6f), min f_K, TV(f_K, target), sup|phi_K|' % (2 * beta))
pairs = 0
for a, e in items:
    if e == 0.0:
        F = F * (1 - x * x / (4 * a * a)); K += 1; s += 1 / (4 * a * a)
    else:
        g2 = (a + 1j * e) ** 2
        F = F * np.abs(1 - x * x / (4 * g2)) ** 2; K += 2; pairs += 1; s += 2 * (1 / (4 * g2)).real
    if K in want:
        f = inv(F)
        print(K, pairs, '%.6f' % s, '%.3e' % f.min(), '%.5f' % np.trapezoid(np.abs(f - ft), u), '%.6f' % np.abs(F).max(), flush=True)
print('time', round(time.time() - t0, 1))
