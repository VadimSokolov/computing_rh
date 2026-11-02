#!/usr/bin/env python3
"""N4, Note 3 numbers: Gumbel form, the spline coefficients, the dual density, PF minors for zeta and for the DH control.
Run: timeout 580 python3 -I N4_note3.py > N4_note3.log. Computed (mpmath, 80 digits), not certified."""
import sys, os, math, time, random
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import numpy as np
import N4_common as C
t0 = time.time()
mp.mp.dps = 40
half = mp.mpf(1) / 2
print('A. Gumbel form')
def k(y):
    return (2 * mp.exp(9 * y / 4) - 3 * mp.exp(5 * y / 4)) * mp.exp(-mp.exp(y))
def form(u, N=40):
    return mp.pi ** (-mp.mpf(1) / 4) * mp.fsum(mp.mpf(n) ** (-half) * k(4 * u + mp.log(mp.pi * n * n)) for n in range(1, N + 1))
for u in (mp.mpf('0.25'), mp.mpf('-0.25'), mp.mpf(1)):
    print('u', u, 'kernel form', mp.nstr(form(u), 18), 'Phi', mp.nstr(C.Phi_raw(abs(u)), 18))
for z in (half, 2, 5):
    I = mp.quad(lambda y: mp.exp(z * y / 4) * k(y), [-mp.inf, -5, 0, 2, 4])
    print('Mellin z', z, mp.nstr(I, 16), mp.nstr((z - 1) / 2 * mp.gamma(mp.mpf(5) / 4 + mp.mpf(z) / 4), 16))
for z in (3, 5):
    I = 2 * mp.quad(lambda u: mp.cosh(z * u) * C.Phi_raw(u, 20), [0, 0.25, 0.5, 1, 2])
    s = half + mp.mpf(z) / 2
    print('theorem factor z', z, mp.nstr(I, 14), mp.nstr(C.xi(s) / 4, 14))
print('\nB. spline coefficients')
zs = C.zero_strings(50)
xh = C.xi(half)
cs = []
for j in range(5):
    g = mp.mpf(zs[j])
    d = mp.diff(lambda t: C.Xi(t), g)
    cs.append(2 * xh / abs(d))
    print('c_%d' % (j + 1), mp.nstr(cs[-1], 8))
print('\nC. dual density')
mp.mp.dps = 80
half = mp.mpf(1) / 2
xh = C.xi(half)
h = mp.mpf(1) / 2
W = [j * h for j in range(0, 801)]
G = [xh / C.xi(half + w / 2) for w in W]
print('G at w=400', mp.nstr(G[-1], 3), 'time', round(time.time() - t0, 1), flush=True)
def pL(x, G=G, W=W, h=h, der=0):
    s = mp.mpf(0)
    for j in range(1, len(W)):
        if der == 0:
            s += G[j] * mp.cos(W[j] * x)
        elif der == 1:
            s += -G[j] * W[j] * mp.sin(W[j] * x)
        else:
            s += -G[j] * W[j] ** 2 * mp.cos(W[j] * x)
    if der == 0:
        s += G[0] / 2
    return s * h / mp.pi
gs = [mp.mpf(z) for z in zs]
for x in ('0.5', '0.55', '0.6', '0.8'):
    x = mp.mpf(x)
    spl = [mp.fsum((-1) ** j * cs[j] * mp.exp(-2 * gs[j] * x) for j in range(n)) for n in (1, 3, 5)]
    print('x', x, 'p_L', mp.nstr(pL(x), 7), '1,3,5 term spline', [mp.nstr(v, 7) for v in spl], flush=True)
print('p_L(0)', mp.nstr(pL(mp.mpf(0)), 8))
def logconc(xs, G, W, h):
    out = []
    for x in xs:
        x = mp.mpf(x)
        p = pL(x, G, W, h); d1 = pL(x, G, W, h, 1); d2 = pL(x, G, W, h, 2)
        out.append((float((d2 * p - d1 * d1) / p ** 2), float(x), float(p)))
    return out
lc = logconc([i / 100 for i in range(0, 61, 3)], G, W, h)
print('zeta: max (log p_L)\'\' on [0,0.6]', max(lc)[:2], 'at 0', lc[0][0], 'min p', min(a[2] for a in lc))
lc = logconc([0.6, 0.8, 1.0], G, W, h)
print('zeta: (log p_L)\'\' at 0.6,0.8,1.0', [a[0] for a in lc])
def lattice(G, W, h, step, nmax):
    arr = [mp.mpf(0)] * (nmax + 1)
    for j in range(len(W)):
        c1 = mp.cos(W[j] * step)
        a, b = mp.mpf(1), c1
        wt = G[j] if j > 0 else G[j] / 2
        arr[0] += wt
        arr[1] += wt * c1
        for n in range(2, nmax + 1):
            a, b = b, 2 * c1 * b - a
            arr[n] += wt * b
    return [v * h / mp.pi for v in arr]
def minors(arr, step, R, orders, ncfg, seed):
    rnd = random.Random(seed)
    M = int(round(R / step))
    out = {}
    for r in orders:
        worst = None; neg = 0
        for _ in range(ncfg):
            xs = sorted(rnd.sample(range(-M, M + 1), r))
            ys = sorted(rnd.sample(range(-M, M + 1), r))
            A = mp.matrix(r, r)
            nrm = mp.mpf(1)
            for i in range(r):
                for j in range(r):
                    A[i, j] = arr[abs(xs[i] - ys[j])]
                nrm *= arr[abs(xs[i] - ys[i])]
            v = mp.det(A) / nrm
            if v < 0:
                neg += 1
            if worst is None or v < worst:
                worst = v
        out[r] = (mp.nstr(worst, 4), neg)
        print('order', r, 'smallest normalised minor', out[r][0], 'negatives', neg, 'of', ncfg, 'time', round(time.time() - t0, 1), flush=True)
    return out
print('\nD. PF minors, zeta, ordered points on a lattice of step 0.0005 in [-0.25,0.25]')
arr = lattice(G, W, h, mp.mpf('0.0005'), 1000)
print('lattice done', round(time.time() - t0, 1), flush=True)
minors(arr, 0.0005, 0.25, [2, 3, 4, 5, 6], 2000, 1)
print('\nE. DH dual density')
mp.mp.dps = 60
half = mp.mpf(1) / 2
Lh = C.dh_Lam(half)
hd = mp.mpf(1) / 4
Wd = [j * hd for j in range(0, 1601)]
kp = mp.mpf(C.KAPPA)
f1 = -(mp.digamma(mp.mpf(1) / 5) + kp * mp.digamma(mp.mpf(2) / 5) - kp * mp.digamma(mp.mpf(3) / 5) - mp.digamma(mp.mpf(4) / 5)) / 5
def Lamd(s):
    return mp.sqrt(5 / mp.pi) * f1 if s == 1 else C.dh_Lam(s)
Gd = [Lh / Lamd(half + w / 2) for w in Wd]
print('f(1)', mp.nstr(f1, 12))
print('Lam(1/2)', mp.nstr(Lh, 12), 'G at 400', mp.nstr(mp.re(Gd[-1]), 3), 'all real positive', all(mp.re(g) > 0 for g in Gd), 'time', round(time.time() - t0, 1), flush=True)
Gd = [mp.re(g) for g in Gd]
lc = logconc([i / 50 for i in range(0, 81, 4)], Gd, Wd, hd)
print('DH: max (log p_L)\'\' on [0,1.6]', max(lc)[:2], 'min p', min(a[2] for a in lc), flush=True)
arrd = lattice(Gd, Wd, hd, mp.mpf('0.002'), 700)
print('lattice done', round(time.time() - t0, 1), flush=True)
print('DH minors on [-0.7,0.7], step 0.002')
minors(arrd, 0.002, 0.7, [2, 3, 4, 5], 1500, 2)
print('\nF. Truncated splines for DH in height order: uniform convergence with no hypothesis, though the stages after the pair are not PF')
d = C.dh_json()
on = sorted(float(s) for s in d['online'])
off = [(float(t), float(s) - 0.5) for s, t in d['offline']]
items = sorted([(g, 0.0) for g in on] + off)
Wf = np.array([float(w) for w in Wd])
Gf = np.array([float(g) for g in Gd])
prod = np.ones_like(Wf)
xs = np.array([0.0, 0.1, 0.2, 0.4, 0.8])
def invert(y):
    return np.array([(np.sum(y * np.cos(Wf * x)) - 0.5 * y[0] - 0.5 * y[-1] * np.cos(Wf[-1] * x)) * 0.25 / math.pi for x in xs])
target = invert(Gf)
print('p_L^DH at', xs, target)
cnt = 0
for a, e in items:
    cnt += 1
    if e == 0.0:
        prod = prod / (1 + Wf ** 2 / (4 * a * a))
    else:
        g2 = (a + 1j * e) ** 2
        prod = prod / np.abs(1 + Wf ** 2 / (4 * g2)) ** 2
    if cnt in (5, 10, 20, 43, 44, 60, 126):
        print('factors', cnt, 'includes pair at 85.7:', a >= 85.69, 'max|p_N-p_L| on those x %.3e' % np.max(np.abs(invert(prod) - target)), flush=True)
print('time', round(time.time() - t0, 1))
