# wp5_analyze.py: fits of the margin law from data/wp5_ladder.json and the comparison with the zeta ladder data/wl_grid.json (run from data/).
import json, math, re, collections
import numpy as np
import mpmath as mp
R = json.load(open('wp5_ladder.json'))
taus = [2, 3, 4, 6, 8, 12, 16, 24, 32]
gen = [r for r in R if r['label'].startswith('genuine') and not r['neg'] and r['stable'] and r.get('lo_frac') is not None]
def key(lab):
    g = 2 if 'y2' in lab else int(re.search(r'g(\d)', lab).group(1)); q = int(re.search(r'q(\d+)', lab).group(1)); return g, q
get = {(r['label'], r['tau'], r['N']): r['log10eps'] for r in gen}
pts = []
for lab in sorted(set(r['label'] for r in gen)):
    g, q = key(lab)
    for t in taus:
        a, b = get.get((lab, t, 16)), get.get((lab, t, 32))
        if a is None or b is None or a > -3: continue
        pts.append((g, q, t, t * math.log(q) / 2, -(b - a) / 16, -(get[(lab, t, 24)] - a) / 8 if (lab, t, 24) in get else None))
print('points (genus, q, tau, L, c16-32):', len(pts))
for g in (1, 2, 3):
    for q in (3, 5, 13, 49):
        row = ['%g:%.2f' % (t, c) for (gg, qq, t, L, c, c2) in pts if gg == g and qq == q]
        print('g%d q%2d' % (g, q), ' '.join(row))
def fit(xs, ys):
    A = np.vstack([xs, np.ones_like(xs)]).T; co = np.linalg.lstsq(A, ys, rcond=None)[0]; return co, np.max(np.abs(A @ co - ys))
sel = [p for p in pts if p[2] / p[0] >= 4]
y = np.array([p[4] for p in sel])
for name, x in (('log10(tau/g)', np.array([math.log10(p[2] / p[0]) for p in sel])), ('log10(L)', np.array([math.log10(p[3]) for p in sel])),
                ('log10(L/g)', np.array([math.log10(p[3] / p[0]) for p in sel]))):
    co, mr = fit(x, y); print('c = a %s + b: a=%.3f b=%.3f max resid %.3f (n=%d)' % (name, co[0], co[1], mr, len(sel)))
# the collapse: c is a function f(x) of x = tau/g = 2L/(g ell) = 2L/(pi rho), rho = g ell/pi the density of zeros per unit r; table of f from all genuine points (mean per x)
fx = collections.defaultdict(list)
for (g, q, t, L, c, c2) in pts:
    if not (g == 3 and q <= 5): fx[round(t / g, 4)].append(c)  # genus 3 over F_3, F_5 has exact kernels at some tau (see checks file)
for (g, tau, c) in [(1, 2.25, .296), (1, 2.5, .537), (2, 4.5, .300), (2, 5.0, .538), (3, 7.5, .539)]: pass
xs_ = sorted(fx); fm = [float(np.mean(fx[x])) for x in xs_]; fs = [float(np.ptp(fx[x])) for x in xs_]
print('\nf(x) = c at x = tau/g (mean over curves, spread):', ' '.join('%.3g:%.2f(%.2f)' % (x, m, s) for x, m, s in zip(xs_, fm, fs)))
f = lambda x: 0.0 if x <= 2 else float(np.interp(math.log(x), np.log(xs_ + [200.0]), fm + [fm[-1] + 4 * math.log10(200 / xs_[-1])]))
print('zeta (wl_grid): observed c over N windows against f(2L/(pi rho(T_N))), T_N = pi N / L, rho(T) = log(T/2pi)/(2pi)')
Z = json.load(open('wl_grid.json'))
by = collections.defaultdict(dict)
for r in Z:
    if 'eps_up' in r: by[r['X']][r['N']] = float(mp.log10(mp.mpf(r['eps_up'])))
for X in (30, 50, 100, 200, 500, 1000, 2000, 5000, 10000):
    L = math.log(X) / 2; ns = sorted(by[X]); row = []
    for i in range(len(ns) - 1):
        n1, n2 = ns[i], ns[i + 1]; nm = (n1 + n2) / 2; T = math.pi * nm / L
        if T < 20: continue
        obs = -(by[X][n2] - by[X][n1]) / (n2 - n1); rho = math.log(T / (2 * math.pi)) / (2 * math.pi); x = 2 * L / (math.pi * rho)
        row.append('N%g-%g: %.2f vs %.2f' % (n1, n2, obs, f(x)))
    print('X=%g L=%.2f' % (X, L), '; '.join(row))
# the fakes: first tau at which the least eigenvalue is certified negative, per N
print('\nfakes: first tau with certified eps<0, for N = 2, 8, 32; and log10|eps| at (tau=32, N=32)')
fk = [r for r in R if r['label'].startswith('fake')]
for lab in sorted(set(r['label'] for r in fk)):
    out = []
    for N in (2, 8, 32):
        t = [r['tau'] for r in fk if r['label'] == lab and r['N'] == N and r['cert_neg']]
        out.append(min(t) if t else None)
    e = [r for r in fk if r['label'] == lab and r['tau'] == 32 and r['N'] == 32][0]
    print(lab, out, e['eps'])
