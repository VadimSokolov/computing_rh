#!/usr/bin/env python3
"""N4, Note 4 numbers: potential of V, the DH control, the cosh spin (ball cover for U_*'''), trig polynomial, harmonic uniform sums.
Run: timeout 580 python3 -I N4_note4.py > N4_note4.log. mpmath, numpy; python-flint for the ball cover."""
import sys, os, math, time
sys.dont_write_bytecode = True
sys.path.insert(0, '/tmp/d1lib')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
from scipy import special
import N4_common as C
from flint import arb, ctx
t0 = time.time()
mp.mp.dps = 60
print('A. potential of V')
U = lambda u: -mp.log(C.Phi_raw(u))
ty = mp.taylor(U, 0, 4)
print('U2(0)', mp.nstr(2 * ty[2], 8), 'U4(0)', mp.nstr(24 * ty[4], 8))
mp.mp.dps = 40
vals = [(mp.diff(U, 0.02 * j, 3), 0.02 * j) for j in range(1, 101)]
print('U3 grid (0,2], 100 points: min', mp.nstr(min(vals)[0], 5), 'at', min(vals)[1], 'at u=2', mp.nstr(vals[-1][0], 6), '64 pi e^8', mp.nstr(64 * mp.pi * mp.exp(8), 6))
print('tail check, U3 - 64 pi e^{4u} at u=1,2:', mp.nstr(mp.diff(U, 1, 3) - 64 * mp.pi * mp.exp(4), 4), mp.nstr(mp.diff(U, 2, 3) - 64 * mp.pi * mp.exp(8), 4))
print('\nB. DH control, in the V scale: density proportional to e^{3u} theta(e^{4u})')
k = mp.mpf(C.KAPPA)
cyc = [1, k, -k, -1, 0]
def thd(x):
    return mp.fsum(cyc[(n - 1) % 5] * n * mp.exp(-mp.pi * n * n * x / 5) for n in range(1, 61))
def pd(u):
    return mp.exp(3 * u) * thd(mp.exp(4 * u))
print('even?', mp.nstr(pd(mp.mpf('0.3')), 15), mp.nstr(pd(mp.mpf('-0.3')), 15))
Ud = lambda u: -mp.log(pd(u))
tyd = mp.taylor(Ud, 0, 4)
print('U_DH2(0)', mp.nstr(2 * tyd[2], 8), 'U_DH4(0)', mp.nstr(24 * tyd[4], 8), '(opus/E in X scale times 4, 16: %s %s)' % (4 * 2.553638949, 16 * 31.61921421))
vals = [(mp.diff(Ud, 0.02 * j, 3), 0.02 * j) for j in range(1, 101)]
print('U_DH3 grid (0,2]: min', mp.nstr(min(vals)[0], 5), 'at', min(vals)[1], 'all positive', all(v[0] > 0 for v in vals))
h = mp.mpf(1) / 64
us = [j * h for j in range(0, 193)]
pv = [pd(u) for u in us]
def momd(kk):
    return 2 * h * mp.fsum(us[j] ** (2 * kk) * pv[j] for j in range(1, len(us))) + (2 * h * pv[0] / 2 if kk == 0 else 0)
m0 = momd(0)
mm = [mp.mpf(1)] + [mp.mpf(0)] * 12
for i in range(1, 7):
    mm[2 * i] = momd(i) / m0
kap = [mp.mpf(0)] * 13
for n in range(1, 13):
    kap[n] = mm[n] - mp.fsum(mp.binomial(n - 1, j - 1) * kap[j] * mm[n - j] for j in range(1, n))
print('DH cumulants 2..12', [mp.nstr(kap[2 * i], 4) for i in range(1, 7)], 'EV4', mp.nstr(mm[4], 5), '3s^4', mp.nstr(3 * mm[2] ** 2, 5), 'Var', mp.nstr(mm[2], 8))
mp.mp.dps = 25
half = mp.mpf(1) / 2
def K3d(z):
    return mp.diff(lambda s: mp.log(C.dh_Lam(s)), half + mp.mpf(z) / 2, 3) / 8
gv = [(K3d(0.05 + (60 - 0.05) * i / 119), 0.05 + (60 - 0.05) * i / 119) for i in range(120)]
print('K3_DH on 120 points in [0.05,60]: max', mp.nstr(max(gv)[0], 4), 'at', max(gv)[1], 'all negative', all(v[0] < 0 for v in gv), 'time', round(time.time() - t0, 1), flush=True)
print('\nC. cosh spin')
mp.mp.dps = 40
A = 2 * mp.pi ** 2; B = 3 * mp.pi
ms = lambda u: A * mp.cosh(9 * u) - B * mp.cosh(5 * u)
Us = lambda u: 2 * mp.pi * mp.cosh(4 * u) - mp.log(ms(u))
tys = mp.taylor(Us, 0, 4)
print('U*2(0)', mp.nstr(2 * tys[2], 6), 'U*4(0)', mp.nstr(24 * tys[4], 7))
vals = [(mp.diff(Us, 0.02 * j, 3), 0.02 * j) for j in range(1, 101)]
print('U*3 grid (0,2]: min', mp.nstr(min(vals)[0], 5), 'at', min(vals)[1], 'all positive', all(v[0] > 0 for v in vals))
print('trig polynomial sign at k pi/9:', [int(mp.sign(A * mp.cos(9 * j * mp.pi / 9) - B * mp.cos(5 * j * mp.pi / 9))) for j in range(19)])
hh = 1 / 4096.0
uu = np.arange(0, 1.6 + hh / 2, hh)
Ph = C.Phi_np(uu)
tz = lambda f: np.trapezoid(f, uu)
p = Ph / (2 * tz(Ph))
def Q(M, u=uu):
    out = np.zeros_like(u)
    for n in range(1, M + 1):
        out += (2 * np.pi ** 2 * n ** 4 * np.cosh(9 * u) - 3 * np.pi * n ** 2 * np.cosh(5 * u)) * np.exp(-2 * np.pi * n ** 2 * np.cosh(4 * u))
    return out
VarV = 2 * tz(uu * uu * p)
for M in (1, 2, 3, 6):
    q = Q(M); q /= 2 * tz(q)
    v = 2 * tz(uu * uu * q)
    c = math.sqrt(VarV / v)
    qs = np.interp(uu / c, uu, q) / c
    print('Q_%d: Var %.6f TV to p %.4f, after matching variances %.4f, min %.2e' % (M, v, tz(np.abs(q - p)), tz(np.abs(qs - p)), q.min()))
f1 = (2 * np.pi ** 2 * np.exp(9 * uu) - 3 * np.pi * np.exp(5 * uu)) * np.exp(-np.pi * np.exp(4 * uu))
f1m = (2 * np.pi ** 2 * np.exp(-9 * uu) - 3 * np.pi * np.exp(-5 * uu)) * np.exp(-np.pi * np.exp(-4 * uu))
print('true symmetrisation of the n=1 term at u=1.5: value', 0.5 * (f1[int(1.5 / hh)] + f1m[int(1.5 / hh)]), 'negative tail, not doubly exponential; rho_* there', Q(1)[int(1.5 / hh)])
zz = np.arange(0.5, 120, 0.01)
qq = Q(1); qq = qq / (2 * tz(qq))
F = np.array([2 * tz(qq * np.cos(z * uu)) for z in zz])
sc = zz[:-1][np.sign(F[:-1]) != np.sign(F[1:])]
print('transform of rho_*: sign changes in (0.5,120):', len(sc), 'first', np.round(sc[:8], 2), 'time', round(time.time() - t0, 1), flush=True)
print('\nD. ball cover of the GHS condition for rho_* (Arb)')
ctx.prec = 128
Aa = 2 * arb.pi() ** 2; Bb = 3 * arb.pi()
def ders(u):
    c9, s9 = (9 * u).cosh(), (9 * u).sinh()
    c5, s5 = (5 * u).cosh(), (5 * u).sinh()
    m0 = Aa * c9 - Bb * c5
    m1 = 9 * Aa * s9 - 5 * Bb * s5
    m2 = 81 * Aa * c9 - 25 * Bb * c5
    m3 = 729 * Aa * s9 - 125 * Bb * s5
    m4 = 6561 * Aa * c9 - 625 * Bb * c5
    L3 = m3 / m0 - 3 * m1 * m2 / m0 ** 2 + 2 * m1 ** 3 / m0 ** 3
    L4 = m4 / m0 - (4 * m1 * m3 + 3 * m2 ** 2) / m0 ** 2 + 12 * m1 ** 2 * m2 / m0 ** 3 - 6 * m1 ** 4 / m0 ** 4
    return 128 * arb.pi() * (4 * u).sinh() - L3, 512 * arb.pi() * (4 * u).cosh() - L4
def cover(a, b, which, depth=0):
    mid = (a + b) / 2; rad = (b - a) / 2
    ball = arb(mid, rad)
    v = ders(ball)[which]
    if v > 0:
        return 1, float(v.lower())
    if depth > 40:
        raise RuntimeError('cover failed at %g' % a)
    n1, l1 = cover(a, mid, which, depth + 1)
    n2, l2 = cover(mid, b, which, depth + 1)
    return n1 + n2, min(l1, l2)
print('U*4>0 on [0,1/64]: cells, min lower bound', cover(0.0, 1 / 64, 1))
print('U*3>0 on [1/64,2]: cells, min lower bound', cover(1 / 64, 2.0, 0))
c9 = 18.0 ** 3 * (1 + 3 / (2 * math.pi)) * 2 * math.exp(-8)
print('tail u>=2: h=1+e^{-18u}-(3/2pi)(e^{-4u}+e^{-14u}), |h^(j)-delta| <= 18^j*2*(1+3/2pi)*e^{-8}; third derivative of log h at most about', 3 * c9, 'against 128 pi sinh 8 =', 128 * math.pi * math.sinh(8))
print('\nE. Doubly exponential tails inside the Griffiths Simon class: S_b=sum U_k/(2(k+b)), X scale, b*=1.25839697')
bs = 1.25839697
print('gamma_E-psi(1+b*)', np.euler_gamma - special.digamma(1 + bs))
mp.mp.dps = 25
xh = C.xi(mp.mpf(1) / 2)
for s in (10.0, 100.0, 1000.0, 10000.0):
    kk = np.arange(1, 10 ** 7 + 1, dtype=float)
    w = s / (2 * (kk + bs))
    big = w > 1e-3
    val = np.empty_like(w)
    val[big] = w[big] + np.log1p(-np.exp(-2 * w[big])) - np.log(2 * w[big])
    val[~big] = w[~big] ** 2 / 6
    ls = np.sum(val) + s * s / (24 * (1e7 + bs))
    mt = s / 2 * math.log(s / (2 * math.pi * math.e))
    lx = float(mp.log(C.xi(mp.mpf(1) / 2 + s) / xh))
    print('s', s, 'log L_S %.3f' % ls, 'log L_X %.3f' % lx, 'main term (s/2)log(s/2 pi e) %.3f' % mt, 'differences %.3f %.3f' % (ls - mt, lx - mt), flush=True)
print('time', round(time.time() - t0, 1))
