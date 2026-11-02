#!/usr/bin/env python3
"""D2: the saturation height b(N,L) as a detection law: the Davenport and Heilbronn form (conductor q=5) and planted quadruples.
b solves (4L/pi) phi = log(q b/(4 pi)), sin phi = T_N/b; for T_N >= (4 pi/q) e^{2L} it is b = (4 pi/q) e^{2L}.
Run: timeout 60 python3 -I D2_detection.py   Status: floating point, a heuristic reading of the law against the book's certified signs.
"""
import json, math
ROOT = '.'


def b_of(N, L, q=1.0):
    T = math.pi * N / L
    bmax = 4 * math.pi / q * math.exp(2 * L)
    if T >= bmax:
        return bmax
    lo, hi = 1e-12, math.pi / 2
    g = lambda p: (4 * L / math.pi) * p - math.log(q * T / (4 * math.pi * math.sin(p)))
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (lo, m) if g(m) > 0 else (m, hi)
    return T / math.sin(0.5 * (lo + hi))


gamma_off = 85.6993
print('DH: X, b_infty, least N with b(N)>=85.70, cells of wl_dh scans: last positive N, first negative N')
scan = json.load(open(f'{ROOT}/data/wl_dh_scan1.json')) + json.load(open(f'{ROOT}/data/wl_dh_scan2.json'))
for X in (30, 50, 100, 200, 500, 1000, 2000):
    L = 0.5 * math.log(X)
    binf = 4 * math.pi / 5 * X
    Nd = next((N for N in range(1, 400) if b_of(N, L, 5.0) >= gamma_off), None)
    cells = sorted((e['N'], e.get('certified_positive'), e.get('certified_negative')) for e in scan if e.get('X') == X and e.get('parity', 'even') == 'even')
    pos = [c[0] for c in cells if c[1]]
    neg = [c[0] for c in cells if c[2]]
    print(f'  {X:5d} {binf:8.1f} {Nd}  lastpos {max(pos) if pos else None} firstneg {min(neg) if neg else None}  Nmax {cells[-1][0] if cells else None}')
print('planted quadruple: gamma, (1/2)log(gamma/2pi), a with b_infty(a)=gamma i.e. (1/2)log(gamma/4pi)')
for g in (51.37, 199.64, 1000.57, 1999.99):
    print(f'  {g:8.2f} {0.5 * math.log(g / (2 * math.pi)):6.3f} {0.5 * math.log(g / (4 * math.pi)):6.3f}')
gam = [float(x) for x in open(f'{ROOT}/data/zeros_8000.txt') if x.strip()]
T = gam[-1]
s = sum(1 / g ** 2 for g in gam)
tail = (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T) + 2 * 3.4 / T ** 2
print('S_0: first 8000', s, 'tail bound', tail, 'total', s + tail, 'top', T)
print('k<=Psi(gamma_k/2pi) for k=1,2:', [(k, gam[k - 1] / (2 * math.pi) * math.log(gam[k - 1] / (2 * math.pi))) for k in (1, 2)])
print('c(X) of the first lemma: X, w, c')
for X in (3, 5, 13, 50, 200, 1000, 1e4, 1e6, 1e12):
    two = math.log(X)
    w = 1.0
    for _ in range(200):
        w = two - math.log(w)
    print(f'  {X:10.0f} {w:7.3f} {4 * math.log(1 + math.log(w) / w):6.3f}   rate of the law at N=X: {(math.pi / (0.5 * two)) * (1 + math.log(2 * two)):6.3f}')
