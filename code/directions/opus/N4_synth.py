#!/usr/bin/env python3
"""N4, Note 1 section 'synthesis': thresholds beta_min(N), the bound |char fn|<=1, positivity at beta=Var/4, the DH control, TV table.
Run: timeout 590 python3 -I N4_synth.py > N4_synth.log
Ball arithmetic (Arb through python-flint) for Q values at sample points; numpy double precision for the sup bound and the TV table.
Scale: the book's v and the note's beta are related by v = 8 beta (X=2V)."""
import sys, os, time, math
sys.dont_write_bytecode = True
sys.path.insert(0, '/tmp/d1lib')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from flint import arb, acb, ctx, arb_poly, fmpz_poly, fmpz
import N4_common as C

t0 = time.time()
NOTE = {2: 6.19, 4: 6.67, 6: 6.85, 8: 6.88, 10: 6.76, 12: 6.70, 16: 6.50, 20: 6.29, 24: 6.10, 30: 5.88, 40: 5.50, 60: 4.91, 80: 4.47, 120: 3.76, 160: 2.09, 200: 1.17, 300: 0.450}
BOOK = {2: .004955, 3: .005416, 5: .005527, 10: .005441, 16: .005222, 20: .005065, 50: .004146, 100: .003303, 150: .002816, 200: .002487, 300: .002059}
KS = [2, 4, 6, 8, 10, 12, 16, 20, 24, 30, 40, 50, 60, 80, 100, 120, 150, 160, 200, 300]
KMAX = 300

def he_table(kmax):
    y = fmpz_poly([0, 1]); Hm = fmpz_poly([1]); H = y
    tab = [[fmpz(1)]]
    for n in range(1, 2 * kmax):
        Hn = y * H - n * Hm
        Hm, H = H, Hn
        if (n + 1) % 2 == 0:
            c = H.coeffs()
            tab.append([c[i] for i in range(0, len(c), 2)])
    return tab

TAB = he_table(KMAX)
print('He table', round(time.time() - t0, 1), flush=True)
ctx.prec = 5200
ZZ = [z.imag for z in acb.zeta_zeros(1, KMAX)]
print('Arb zeros to', KMAX, round(time.time() - t0, 1), flush=True)

def elem(factors):
    E = [arb(1)]
    for A, B in factors:
        n = len(E)
        if B is None:
            new = E + [arb(0)]
            for j in range(1, n + 1):
                new[j] = new[j] + A * E[j - 1]
        else:
            new = E + [arb(0), arb(0)]
            for j in range(1, n + 1):
                new[j] = new[j] + A * E[j - 1]
            for j in range(2, n + 2):
                new[j] = new[j] + B * E[j - 2]
        E = new
    return E

def Qpoly(E, v):
    K = len(E) - 1
    co = [arb(0) for _ in range(K + 1)]
    w = arb(1); vi = 1 / arb(v)
    for j in range(K + 1):
        wj = E[j] * w
        for i, c in enumerate(TAB[j]):
            co[i] += wj * arb(c)
        w *= vi
    return arb_poly(co)

UNCERTAIN = [0]
def scan(p, K, dy, hint=None, full_min=False):
    ymax = math.sqrt(8 * K + 2)
    ys = np.arange(0, ymax + dy, dy)
    if hint is not None:
        order = sorted(ys, key=lambda a: abs(a - hint))
    else:
        order = ys
    best = None
    for y in order:
        y = float(y)
        val = p(arb(y) * arb(y))
        if val < 0:
            if not full_min:
                return (y, val)
            if best is None or val < best[1]:
                best = (y, val)
        elif not val > 0:
            UNCERTAIN[0] += 1
    return best

def threshold(E, K, hi, steps=18, dy=0.02):
    lo = hi / 2000.0
    hint = None
    r = scan(Qpoly(E, lo), K, dy)
    if r is None:
        return None
    hint = r[0]
    for _ in range(steps):
        mid = math.sqrt(lo * hi)
        r = scan(Qpoly(E, mid), K, dy, hint)
        if r is None:
            hi = mid
        else:
            lo = mid; hint = r[0]
    while True:
        r = scan(Qpoly(E, hi), K, dy / 4.0, hint)
        if r is None:
            break
        lo = hi; hi *= 1.0005; hint = r[0]
    return lo, hi, hint

print('\nA. Threshold v*(K) (grid in y, step 0.02 then 0.005; ball evaluation), beta_min = v*/8. note = 1e4 beta from Note 1 tables; book = fm:tab:fK.', flush=True)
print('K, s_K, v_lo, v_hi, 1e4*beta_true, 1e4*beta_note, note/true, book v*, y of last negative value')
RESULT = {}
for K in KS:
    ctx.prec = 12 * K + 800
    fac = [(1 / (ZZ[k] * ZZ[k]), None) for k in range(K)]
    E = elem(fac)
    sK = float(sum((1 / (ZZ[k] * ZZ[k]) for k in range(K)), arb(0)).mid())
    r = threshold(E, K, sK * 1.0000001)
    lo, hi, hint = r
    RESULT[K] = hi
    print(K, '%.6f' % sK, '%.7f' % lo, '%.7f' % hi, '%.3f' % (hi / 8 * 1e4), NOTE.get(K, ''), ('%.3f' % (NOTE[K] / (hi / 8 * 1e4))) if K in NOTE else '', BOOK.get(K, ''), hint, 't=%.0f' % (time.time() - t0), flush=True)
print('uncertain ball signs so far', UNCERTAIN[0])

print('\nB. Q at the note\'s own beta, certified negative values (Arb)')
for K in [10, 20, 120, 160, 200, 300]:
    ctx.prec = 12 * K + 800
    E = elem([(1 / (ZZ[k] * ZZ[k]), None) for k in range(K)])
    v = 8 * NOTE[K] * 1e-4
    r = scan(Qpoly(E, v), K, 0.02, None, True)
    print('N', K, 'beta_note 1e4', NOTE[K], 'most negative Q on grid:', None if r is None else (r[0], str(r[1])[:40]), flush=True)

gam = C.zeros_float(8000)
print('\nC. Necessary condition |char fn|<=1: beta >= beta1(N) = sup_x x^-2 sum_{j<=N} log|1-x^2/(4 gamma_j^2)| (double precision)')
print('N, beta1*1e4, argmax x, 2 gamma_N, 1e4 beta_true, 1e4 beta_note, log10 sup|char fn| at beta_note')
for N in KS:
    x = np.arange(0.5, 8 * gam[N - 1], 0.01)
    S = np.zeros_like(x)
    for g in gam[:N]:
        S += np.log(np.abs(1 - x * x / (4 * g * g)) + 1e-300)
    r = S / (x * x); i = int(np.argmax(r))
    bn = NOTE.get(N)
    sup = (np.max(S - bn * 1e-4 * x * x) / math.log(10)) if bn else float('nan')
    print(N, '%.4f' % (r[i] * 1e4), '%.2f' % x[i], '%.2f' % (2 * gam[N - 1]), '%.3f' % (RESULT[N] / 8 * 1e4), bn if bn else '', '%.2f' % sup, flush=True)

print('\nD. DH control. Factors in height order, an off line pair is one factor 1+A tau+B tau^2 (K counts 2).')
d = C.dh_json()
items = []
ctx.prec = 2400
for s in d['online']:
    items.append((float(s), 'on', arb(s), None))
for sig, t in d['offline']:
    items.append((float(t), 'off', arb(t), arb(sig) - arb(1) / 2))
items.sort(key=lambda a: a[0])
print('first off line heights', [round(i[0], 3) for i in items if i[1] == 'off'], 'online below 85.7:', sum(1 for i in items if i[1] == 'on' and i[0] < 85.7))
def dh_factors(items):
    fac = []
    for h, typ, a, eta in items:
        if typ == 'on':
            fac.append((1 / (a * a), None))
        else:
            den = (a * a + eta * eta) ** 2
            A = 2 * (a * a - eta * eta) / den
            B = 1 / den
            fac.append((A, B))
    return fac
idx_first = [i for i, it in enumerate(items) if it[1] == 'off'][0]
print('factors, K, off line pairs included, s_K(with Re budget), v*_DH lo, hi, 1e4 beta, Q(min) at v=s_K positive?')
for nf in [5, 10, 20, idx_first, idx_first + 1, idx_first + 5, idx_first + 15, len(items)]:
    its = items[:nf]
    fac = dh_factors(its)
    K = sum(1 if b is None else 2 for a, b in fac)
    ctx.prec = 14 * K + 800
    fac = dh_factors(its)
    E = elem(fac)
    sK = sum((a for a, b in fac), arb(0))
    sKf = float(sK.mid())
    r0 = scan(Qpoly(E, sKf), K, 0.02, None, True)
    r = threshold(E, K, sKf * 1.0000001, steps=14)
    print(nf, K, sum(1 for i in its if i[1] == 'off'), '%.6f' % sKf, '%.7f' % r[0], '%.7f' % r[1], '%.3f' % (r[1] / 8 * 1e4), 'negative at s_K:' , r0 is not None, 't=%.0f' % (time.time() - t0), flush=True)
print('uncertain ball signs', UNCERTAIN[0])

print('\nE. TV at beta=1e-3 against V+sqrt(2 beta)Z (double precision, trapezoid)')
Vx = 0.0115525  # refined below
Var = float(C.mp.mpf('0.0115524965577'))  # half of 0.0231049931154 (opus/E section 4.9)
tau = Var / 2 - np.sum(1 / (4 * gam ** 2))
x = np.arange(0, 500.0001, 0.02)
beta = 1e-3
gauss = np.exp(-beta * x * x)
prod = np.ones_like(x)
Ns = [2, 10, 40, 100, 300, 1000, 3000]
partial = {}
for j, g in enumerate(gam, 1):
    prod *= (1 - x * x / (4 * g * g))
    if j in Ns:
        partial[j] = prod.copy()
true = prod * np.exp(-tau * x * x) * gauss
u = np.arange(0, 1.2001, 0.002)
CM = np.cos(np.outer(u, x))
wt = np.full_like(x, 0.02); wt[0] = 0.01
g_true = CM @ (true * wt) / np.pi
print('tail tau', tau, 'mass check', 2 * np.trapezoid(g_true, u))
for N in Ns:
    f = CM @ (partial[N] * gauss * wt) / np.pi
    tv = np.trapezoid(np.abs(f - g_true), u)
    tv2 = 0.242 * 2 * (Var / 2 - np.sum(1 / (4 * gam[:N] ** 2))) / (Var + 2 * beta)
    print('N', N, 'TV %.4f' % tv, 'min f %.3e' % f.min(), 'TV*N/log(N)^2 %.3f' % (tv * N / math.log(N) ** 2), 'gaussian estimate %.4f' % tv2, flush=True)

print('\nF. Hermite form against the Fourier integral, N=10, beta=1e-3')
ctx.prec = 2000
E10 = elem([(1 / (ZZ[k] * ZZ[k]), None) for k in range(10)])
v = 8e-3
p10 = Qpoly(E10, v)
for uu in (0.1, 0.3):
    xx = 2 * uu
    y = xx / math.sqrt(v)
    fX = math.exp(-y * y / 2) / math.sqrt(2 * math.pi * v) * float(p10(arb(y) * arb(y)).mid())
    fV = 2 * fX  # density of V at u from the density of X=2V at 2u
    fFT = np.sum(np.cos(x * uu) * partial_10 * gauss * wt) / np.pi if False else None
    print('u', uu, 'Hermite', fV)
print('Fourier', [float((np.cos(x * uu) * P10 * gauss * wt).sum() / np.pi) for uu in (0.1, 0.3) for P10 in [np.prod([1 - x * x / (4 * g * g) for g in gam[:10]], axis=0)]])

print('\nG. The bound against the variance preserving width: v1(K)=8 beta1(K) and 2 beta_K=2 sum_{k>K} gamma_k^-2 (book fm:thm:zc, f_K=f_{K,2 beta_K})')
P1 = 0.0231049931154
print('K, T=gamma_K, v1, 2beta_K, ratio, sup log|phi_K| at v=2beta_K, argmax theta/T')
for K in [100, 300, 400, 800, 2000, 4000, 8000]:
    th = np.arange(1.0, 3 * gam[K - 1], 0.02 if K <= 800 else 0.1)
    S = np.zeros_like(th)
    for g in gam[:K]:
        S += np.log(np.abs(1 - th * th / (g * g)) + 1e-300)
    r = S / (th * th); i = int(np.argmax(r))
    bK = P1 - np.sum(1 / gam[:K] ** 2)
    sup = np.max(S - bK * th * th); j = int(np.argmax(S - bK * th * th))
    print(K, '%.3f' % gam[K - 1], '%.6f' % (2 * r[i]), '%.6f' % (2 * bK), '%.4f' % (r[i] / bK), '%.3f' % sup, '%.3f' % (th[j] / gam[K - 1]), flush=True)
from scipy import integrate
def Lfun(T, th):
    head = np.sum(np.log(np.abs(1 - th * th / gam ** 2)))
    if T <= gam[-1]:
        return None
    val, err = integrate.quad(lambda lt: (1 / (2 * np.pi)) * math.log(math.exp(lt) / (2 * np.pi)) * math.log(th * th / math.exp(2 * lt) - 1) * math.exp(lt), math.log(gam[-1]), math.log(T), limit=400)
    return head + val
print('extrapolation with the smooth density beyond gamma_8000: T, K=N(T), v1/2beta_K, sup log|phi_K| at 2beta_K, argmax r')
for T in [1e4, 2e4, 4e4, 6e4, 8e4, 1e5, 2e5, 5e5, 1e6, 1e7]:
    l = math.log(T / (2 * np.pi))
    bK = (1 / (2 * np.pi)) * (l + 1) / T
    best = None; bestr = None
    for r in np.arange(1.02, 1.8, 0.02):
        th = r * T
        L = Lfun(T, th)
        a = L - bK * th * th
        q = L / (th * th)
        if best is None or a > best[0]:
            best = (a, r)
        if bestr is None or q > bestr:
            bestr = q
    print('%.3g' % T, '%.4g' % ((T / (2 * np.pi)) * (l - 1) + 7 / 8), '%.4f' % (bestr / bK), '%.2f' % best[0], '%.2f' % best[1], flush=True)
print('total time', round(time.time() - t0, 1))
