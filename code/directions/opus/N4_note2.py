#!/usr/bin/env python3
"""N4, Note 2 numbers. Run: timeout 580 python3 -I N4_note2.py > N4_note2.log. Computed (mpmath, scipy), not certified."""
import sys, os, math, time
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, mpmath as mp
from scipy import special
import N4_common as C
t0 = time.time()
gam = C.zeros_float(8000)
def prime_powers(X):
    s = np.ones(X + 1, dtype=bool); s[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    ps = np.nonzero(s)[0]
    n = []; Lam = []
    for p in ps:
        lp = math.log(p)
        q = int(p)
        while q <= X:
            n.append(q); Lam.append(lp); q *= int(p)
    n = np.array(n, dtype=float); Lam = np.array(Lam)
    o = np.argsort(n)
    return n[o], Lam[o]
print('A. heat trace table')
mp.mp.dps = 25
pn, pL = prime_powers(100000)
def arch(b):
    b = mp.mpf(b)
    f = lambda r: mp.exp(-b * r * r) * mp.re(mp.digamma(mp.mpf(1) / 4 + 0.5j * r))
    R = mp.sqrt(70 / b)
    pts = [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    pts = [p for p in pts if p < R] + [R]
    I = mp.quad(f, pts) / mp.pi
    return 2 * mp.exp(b / 4) - mp.log(mp.pi) / mp.sqrt(4 * mp.pi * b) + I
for b in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05):
    A = arch(b)
    pr = -2 / mp.sqrt(4 * mp.pi * b) * mp.fsum(mp.mpf(float(l)) / mp.sqrt(int(n)) * mp.exp(-mp.log(int(n)) ** 2 / (4 * b)) for n, l in zip(pn[:3000], pL[:3000]))
    Z = 2 * mp.fsum(mp.exp(-b * mp.mpf(float(g)) ** 2) for g in gam[:1500])
    print('b', b, 'zeros', mp.nstr(Z, 16), 'arch', mp.nstr(A, 16), 'prime', mp.nstr(pr, 5), 'diff', mp.nstr(Z - A - pr, 3), flush=True)
print('\nB. K(z) from primes up to 2e5')
pn2, pL2 = prime_powers(200000)
Ln = np.log(pn2); wgt = pL2 / np.sqrt(pn2)
tau = np.arange(-30, 40.0001, 0.02)
rr = np.exp(tau)
rp = np.real(special.digamma(0.25 + 0.5j * rr))
def Afun(b):
    return 2 * np.exp(b / 4) - math.log(math.pi) / math.sqrt(4 * math.pi * b) + np.sum(np.exp(-b * rr * rr) * rp * rr) * 0.02 / math.pi
def Kcalc(lo, step=0.05):
    n = int(round((math.log(0.3) - lo) / step)) + 1
    lb = np.linspace(lo, math.log(0.3), n)
    bs = np.exp(lb)
    Aa = np.array([Afun(b) for b in bs])
    Pr = np.array([np.sum(wgt * np.exp(-Ln * Ln / (4 * b))) / math.sqrt(4 * math.pi * b) for b in bs])
    return lb, bs, 0.5 * Aa, Pr
lb, bs, HA, Pr = Kcalc(-40.0)
lb2, bs2, HA2, Pr2 = Kcalc(-60.0)
mp.mp.dps = 25
half = mp.mpf(1) / 2
xh = C.xi(half)
def tr(y, lb):
    return (np.sum(y) - 0.5 * (y[0] + y[-1])) * (lb[1] - lb[0])
for z in (0.5, 3, 6, 10, 20):
    c = z * z / 4
    full = tr(-np.expm1(-bs * c) * (HA - Pr), lb)
    ar = tr(-np.expm1(-bs * c) * HA, lb)
    full2 = tr(-np.expm1(-bs2 * c) * (HA2 - Pr2), lb2)
    ex = mp.log(C.xi(half + mp.mpf(z) / 2) / xh)
    print('z', z, 'primes+arch[e^-40,.3] %.8f' % full, 'arch only %.8f' % ar, 'primes+arch[e^-60,.3] %.8f' % full2, 'exact', mp.nstr(ex, 9), 'ratio arch/exact %.3f' % (ar / float(ex)), flush=True)
print('\nC. Euler half plane; exact K at z=3,4,6 against zeros 8000 + tail')
Var = 0.0115524965577
for z in (3, 4, 6):
    ex = mp.log(C.xi(half + mp.mpf(z) / 2) / xh)
    zs = np.sum(np.log1p(z * z / (4 * gam ** 2)))
    tail = z * z / 4 * (2 * Var - np.sum(1 / gam ** 2))
    print('z', z, 'exact', mp.nstr(ex, 9), 'zeros8000+tail %.7f' % (zs + tail))
print('\nD. Spectrogram, sigma=2, t in [5,60], 1101 points')
rg = np.arange(-12, 76.0001, 0.01)
rv = np.real(special.digamma(0.25 + 0.5j * rg))
ts = np.linspace(5, 60, 1101)
def spec(X):
    n, L = prime_powers(X)
    Ln = np.log(n); w = L / np.sqrt(n) * np.exp(-Ln * Ln / 8)
    out = np.zeros_like(ts)
    for i, t in enumerate(ts):
        arch_i = np.sum(np.exp(-2 * (rg - t) ** 2) * rv) * 0.01 / math.pi
        hi = 2 * (np.exp(-2 * ((0.5j - t) ** 2)) + np.exp(-2 * ((0.5j + t) ** 2))).real
        pr = 2 / math.sqrt(2 * math.pi) * np.sum(w * np.cos(t * Ln))
        out[i] = hi - math.log(math.pi) / math.sqrt(2 * math.pi) + arch_i - pr
    return out
lhs = np.array([2 * np.sum(np.exp(-2 * (gam[:200] - t) ** 2) + np.exp(-2 * (gam[:200] + t) ** 2)) for t in ts])
for X in (10 ** 6, 5 * 10 ** 6):
    R = spec(X)
    d = np.abs(R - lhs)
    print('X', X, 'max|RHS-LHS|', d.max(), 'half of it', d.max() / 2, 'at t', ts[int(np.argmax(d))], flush=True)
offs = []
for g in gam[:13]:
    m = np.abs(ts - g) < 0.5
    tt = ts[m][int(np.argmax(R[m]))]
    offs.append(abs(tt - g))
print('peak offsets', np.round(offs, 4), 'max', max(offs))
print('\nE. Truncation table at z=3')
c = 2.25
g1 = gam[0]
n3, L3 = prime_powers(200000)
Ln3 = np.log(n3); w3 = L3 / np.sqrt(n3)
exact3 = float(mp.log(C.xi(half + mp.mpf(3) / 2) / xh))
for X in (2, 3, 5, 10, 20, 50, 100):
    bX = math.log(X) / (2 * g1)
    cut = np.sum(special.exp1(bX * gam ** 2) - special.exp1(bX * (gam ** 2 + c)))
    lbb = np.linspace(-14, math.log(bX), 3001)
    bb = np.exp(lbb)
    sel = n3 > X
    val = np.array([(1 - np.exp(-c * b)) / math.sqrt(4 * math.pi * b) * np.sum(w3[sel] * np.exp(-Ln3[sel] ** 2 / (4 * b))) for b in bb])
    miss = (np.sum(val) - 0.5 * (val[0] + val[-1])) * (lbb[1] - lbb[0])
    print('X', X, 'b_X %.4f' % bX, 'cut %.2e' % cut, 'missing %.2e' % miss, 'K_X-K %.2e' % (miss - cut), 'X^-g1/2 %.2e' % X ** (-g1 / 2), 'ratio %.4f' % (cut / X ** (-g1 / 2)), 'K(3)=%.8f' % exact3, flush=True)
print('\nF. Identities')
mp.mp.dps = 20
La, cc = mp.mpf('1.7'), mp.mpf('0.9')
lhsv = mp.quad(lambda b: (1 - mp.exp(-b * cc)) * mp.exp(-La ** 2 / (4 * b)) * b ** (-1.5), [0, 1, 10, mp.inf]) / mp.sqrt(4 * mp.pi)
print('prime Frullani', mp.nstr(lhsv, 12), mp.nstr((1 - mp.exp(-La * mp.sqrt(cc))) / La, 12))
x_, g_ = mp.mpf('0.3'), mp.mpf(7) + 0.2j
print('subordination', mp.nstr(mp.quad(lambda b: mp.exp(-b * g_ ** 2) * x_ / mp.sqrt(mp.pi) * b ** (-1.5) * mp.exp(-x_ ** 2 / b), [0, 0.1, 1, mp.inf]), 12), mp.nstr(mp.exp(-2 * g_ * x_), 12))
print('1/phi at x=5,10,20', [mp.nstr(xh / C.xi(half + 0.5j * x), 5) for x in (5, 10, 20)])
print('\nG. Finite order monotonicity for the DH control: first order k with (-1)^k W1^(k)(x)<0, zeros to height 200')
d = C.dh_json()
on = np.array([float(s) for s in d['online']])
off = np.array([complex(float(t), float(s) - 0.5) for s, t in d['offline']])
gc = np.concatenate([on + 0j, off, np.conj(off)])
lg = np.log(gc)
def minsign(k, tgrid):
    xs = k / (2 * tgrid)
    best = np.inf; at = None
    for x in xs:
        ex = k * lg - 2 * gc * x
        ref = np.max(ex.real)
        s = np.sum(np.exp(ex - ref)).real
        if s < best:
            best = s; at = x
    return best, at
a0 = off[0].real
found = None
for k in [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 15000, 20000, 25000, 30000, 40000, 50000, 70000, 100000, 150000]:
    tg = np.linspace(a0 - 4, a0 + 4, 1601)
    m, at = minsign(k, tg)
    print('k', k, 'min normalised sum near the pair %.4e' % m, 'at x %.3f' % at, flush=True)
    if m < 0 and found is None:
        found = k
print('first negative among tested k:', found)
lo_k, hi_k = 30000, 40000
while hi_k - lo_k > 100:
    mid = (lo_k + hi_k) // 2
    if minsign(mid, np.linspace(a0 - 2, a0 + 2, 1601))[0] < 0:
        hi_k = mid
    else:
        lo_k = mid
print('DH first negative order near the pair at 85.7 is in', (lo_k, hi_k), 'k/a0^2 =', hi_k / a0 ** 2)
for k in (50, 90, 200):
    print('DH sign bound: at order', k, 'all terms positive for x <', (math.pi / 2 - k * off[0].imag / off[0].real) / (2 * off[0].imag))
print('time', round(time.time() - t0, 1))
