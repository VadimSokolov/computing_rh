#!/usr/bin/env python3
"""N4, Note 1 numbers other than the synthesis thresholds. Run: timeout 580 python3 -I N4_note1.py > N4_note1.log. Computed (mpmath, numpy), not certified."""
import sys, os, math, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
import N4_common as C
mp.mp.dps = 50
half = mp.mpf(1) / 2
print('A. normalisation')
xh = C.xi(half)
print('xi(1/2)', mp.nstr(xh, 12), 'H0(0)=xi/4', mp.nstr(xh / 4, 10))
h = mp.mpf(1) / 128
us = [j * h for j in range(0, 225)]
ph = [C.Phi_raw(u) for u in us]
def H0(z):
    return 2 * h * (ph[0] / 2 + mp.fsum(ph[j] * mp.cos(z * us[j]) for j in range(1, len(us))))
for z in (0, 5, 20):
    a = 4 * H0(z); b = C.Xi(mp.mpf(z) / 2)
    print('z', z, '4H0', mp.nstr(a, 18), 'Xi(z/2)', mp.nstr(b, 18), 'rel', mp.nstr(abs(a - b) / abs(b), 3))
print('\nB. tilts')
hh = 1 / 4096.0
u = np.arange(0, 1.6 + hh / 2, hh)
Ph = C.Phi_np(u)
tz = lambda f: np.trapezoid(f, u)
mass = 2 * tz(Ph)
p = Ph / mass
EV2 = 2 * tz(u * u * p)
print('mass', mass, 'Var V', EV2, 'const=(1/2)E|V^2-EV^2|', 0.5 * 2 * tz(np.abs(u * u - EV2) * p))
for t in (0.01, 0.05, 0.10, 0.20, 0.22, 0.50):
    pt = p * np.exp(t * u * u); pt /= 2 * tz(pt)
    tv = tz(np.abs(pt - p))
    d = pt - p
    tail = np.concatenate([np.cumsum(((d[1:] + d[:-1]) / 2 * hh)[::-1])[::-1], [0.0]])
    print('t', t, 'TV %.4e' % tv, 'dK %.4e' % np.max(np.abs(tail)), 'Var %.7f' % (2 * tz(u * u * pt)), 'TV/t %.7f' % (tv / t))
print('\nC. Xi at multiples of gamma_1')
zs = C.zero_strings(8000)
gamm = [mp.mpf(s) for s in zs]
g1 = gamm[0]
print('Xi(2g1)', mp.nstr(C.Xi(2 * g1), 6), 'Xi(3g1)', mp.nstr(C.Xi(3 * g1), 6), 'bracket', float(gamm[2]), float(gamm[3]), float(gamm[6]), float(gamm[7]))
print('\nD. cumulants')
mp.mp.dps = 60
h2 = mp.mpf(1) / 128
pv = [C.Phi_raw(j * h2) for j in range(0, 225)]
def mom(k):
    return 2 * h2 * mp.fsum((j * h2) ** (2 * k) * pv[j] for j in range(1, 225)) + (2 * h2 * pv[0] / 2 if k == 0 else 0)
m0 = mom(0)
M = {i: mom(i) / m0 for i in range(1, 6)}
mm = [mp.mpf(1)] + [0] * 10
for i in range(1, 6):
    mm[2 * i] = M[i]
kap = [mp.mpf(0)] * 11
for n in range(1, 11):
    kap[n] = mm[n] - mp.fsum(mp.binomial(n - 1, k - 1) * kap[k] * mm[n - k] for k in range(1, n))
gf = gamm[:400]
for k in range(1, 6):
    zsum = (-1) ** (k + 1) * mp.factorial(2 * k) / k * mp.fsum((2 * g) ** (-2 * k) for g in gf)
    print('kappa_%d' % (2 * k), mp.nstr(kap[2 * k], 8), 'zeros400', mp.nstr(zsum, 8))
T = gf[-1]
print('kappa2 zeros400+smooth tail', mp.nstr(mp.fsum(1 / (2 * g * g) for g in gf) + (mp.log(T / (2 * mp.pi)) + 1) / (4 * mp.pi * T), 8), 'zeros 8000 raw', mp.nstr(mp.fsum(1 / (2 * g * g) for g in gamm), 9))
print('EV4', mp.nstr(mm[4], 6), '3 sigma^4', mp.nstr(3 * mm[2] ** 2, 6), 'Var', mp.nstr(mm[2], 10))
print('\nE. K third derivative')
mp.mp.dps = 30
def K3(z):
    return mp.diff(lambda s: mp.log(C.xi(s)), half + mp.mpf(z) / 2, 3) / 8
for z in (1, 12, 80):
    print('K3(%s)' % z, mp.nstr(K3(z), 4))
vals = [(float(K3(0.05 + (60 - 0.05) * i / 239)), 0.05 + (60 - 0.05) * i / 239) for i in range(240)]
print('grid 240 in [0.05,60]: max', max(vals), 'min', min(vals))
s = 11.0
print('bound at s=11: 2/s^3+2/(s-1)^3-1/(2s^2)=', 2 / s ** 3 + 2 / (s - 1) ** 3 - 1 / (2 * s * s), ' fable bound s=9 version:', 1 / 729 + 2 / 512 - 1 / 162)
print('\nF. certificate')
mp.mp.dps = 60
gfl = np.array([float(g) for g in gamm])
def nearest(t):
    i = int(np.searchsorted(gfl, float(t)))
    c = [j for j in (i - 1, i, i + 1) if 0 <= j < len(gamm)]
    j = min(c, key=lambda q: abs(gamm[q] - t))
    return j, abs(gamm[j] - t)
idx = [a for a in range(len(gamm)) if 2 * g1 < gamm[a] <= 41 * g1]
print('ordinates in (2g1,41g1]:', len(idx), 'indices', idx[0] + 1, idx[-1] + 1)
res = sorted(((nearest(2 * gamm[a] - g1)[1], a + 1, nearest(2 * gamm[a] - g1)[0] + 1) for a in idx), key=lambda r: r[0])
for r in res[:3]:
    print('dist', mp.nstr(r[0], 8), 'a', r[1], 'b', r[2], 't', mp.nstr(2 * gamm[r[1] - 1] - g1, 14), 'gamma_b', mp.nstr(gamm[r[2] - 1], 14))
allidx = [a for a in range(1, len(gamm)) if 2 * gamm[a] - g1 <= gamm[-1]]
res = sorted(((nearest(2 * gamm[a] - g1)[1], a + 1, nearest(2 * gamm[a] - g1)[0] + 1) for a in allidx), key=lambda r: r[0])
print('three term progressions from gamma_1 up to gamma_8000: tested', len(allidx), 'smallest', [(mp.nstr(r[0], 5), r[1], r[2]) for r in res[:4]], 'K0 reach', int((float(gamm[allidx[-1]]) / float(g1) - 1) / 2))
refl = [a for a in range(len(gamm)) if gamm[a] - 2 * g1 >= g1]
res = sorted(((nearest(gamm[a] - 2 * g1)[1], a + 1, nearest(gamm[a] - 2 * g1)[0] + 1) for a in refl), key=lambda r: r[0])
print('reflection test gamma_a-2 gamma_1 ordinate? tested', len(refl), 'smallest', [(mp.nstr(r[0], 5), r[1], r[2]) for r in res[:4]], 'K0 reach', int((float(gamm[-1]) / float(g1) - 1) / 2))
print('\nG. Monte Carlo of Proposition (Williams)')
rng = np.random.default_rng(20261101)
wn = 1 / (np.pi * np.arange(1, 401) ** 2.0)
tail = 2 * (np.pi ** 2 / 6 - np.sum(1 / np.arange(1, 401) ** 2.0)) / np.pi
Rv = []
for _ in range(10):
    Rv.append(rng.gamma(2.0, 1.0, size=(20000, 400)) @ wn + tail)
R2 = np.concatenate(Rv)
w = R2 ** 0.25
print('E R^{1/2}', w.mean(), '+-', w.std() / math.sqrt(len(w)), '2 xi(1/2)', float(2 * C.xi(half)))
L = 0.25 * np.log(R2)
print('size biased mean, var of (1/4)log R^2', np.sum(w * L) / np.sum(w), np.sum(w * L * L) / np.sum(w) - (np.sum(w * L) / np.sum(w)) ** 2)
