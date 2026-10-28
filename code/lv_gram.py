# lv_gram.py: numerical claims of Chapter 29, part 2 (open problems chapter, Section "Problem 2" and "Numerical checks"): the optimal distance
# d_N(c) of e^-x from span{sech(n x), n <= N} in H_c (table tab:dN: c = 0.60, 0.75, 1.00 and N = 4, 8, 16, 32, with max |a_n| of the optimal
# coefficients) and the sharpness experiment, squares {n^2} against integers {n} at c = 3/2 and c = 1/2 (N = 4, 8, 16, 24). Added for the book
# in October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. H_c = L^2((0,inf), x^(2c-1) dx). d^2 = ||e^-x||_c^2 - b^T G^(-1) b, ||e^-x||_c^2 = Gamma(2c)/2^(2c), a = G^(-1) b,
# G_jk = <sech(f_j x), sech(f_k x)>_c, b_j = <e^-x, sech(f_j x)>_c.
# Method. Every integral is the trapezoid rule in u = log x (integrand x^(2c) F(e^u), analytic in the strip |Im u| < pi/2, the poles of sech
# at f x = i pi/2, so the rule converges like exp(-pi^2/h)), truncated below at e^(2c u) = 10^-(dps - 15) and above where sech(x) < 10^-(dps - 15).
# The Gram matrix of the nested sets is the leading block of that of the largest set. The system is solved at dps digits by LU.
# Error control: every run is repeated at dps + 30 with a smaller step h, and the maximum relative difference of d_N and of max |a_n| is logged.
# Output: data/lv_gram.json and data/lv_gram.log. Usage: python3 lv_gram.py [--quick] [--out FILE] [--dps D]; the full case takes some minutes
# (cluster: bash misc/tools/hopper_run.sh -c 1 -t 120 -g lv_gram.json -g lv_gram.log code/lv_gram.py -- --out lv_gram.json).
import argparse, json, math, time
import mpmath as mp
ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/lv_gram.json'); ap.add_argument('--dps', type=int, default=0)
a = ap.parse_args()
logf = open(a.out.replace('.json', '.log'), 'w')
def log(*x):
    s = ' '.join(str(i) for i in x); print(s); logf.write(s + '\n'); logf.flush()
def solve(freqs, c, Ns, dps, hfac):
    mp.mp.dps = dps; eps = dps - 15
    c = mp.mpf(c); ulo = -mp.log(10) * eps / (2 * c)
    uhi = mp.log(mp.log(10) * eps / 2 / min(freqs)) + mp.mpf('0.2')
    h = mp.pi**2 / (mp.log(10) * eps) * hfac; m = int(mp.ceil((uhi - ulo) / h)); h = (uhi - ulo) / m
    us = [ulo + k * h for k in range(m + 1)]; xs = [mp.e**u for u in us]
    wt = [h * x**(2 * c) for x in xs]; wt[0] /= 2; wt[-1] /= 2
    ex = [mp.exp(-x) for x in xs]
    S = []
    for f in freqs:
        row = []
        for x in xs:
            e = mp.exp(-f * x); row.append(2 * e / (1 + e * e))
        S.append(row)
    nf = len(freqs); G = mp.matrix(nf, nf); b = mp.matrix(nf, 1)
    for j in range(nf):
        Sj = [wt[k] * S[j][k] for k in range(m + 1)]
        b[j] = mp.fsum(Sj[k] * ex[k] for k in range(m + 1))
        for l in range(j + 1):
            G[j, l] = G[l, j] = mp.fsum(Sj[k] * S[l][k] for k in range(m + 1))
    nrm = mp.gamma(2 * c) / 2**(2 * c); out = {}
    for N in Ns:
        GN = G[:N, :N]; bN = b[:N, :]
        av = mp.lu_solve(GN, bN); d2 = nrm - sum(bN[i] * av[i] for i in range(N))
        out[N] = (mp.sqrt(d2) if d2 > 0 else mp.mpf(-1), max(abs(av[i]) for i in range(N)))
    return out, m
res = {}
cases = []
Nt = [4, 8, 16, 32] if not a.quick else [4, 8]
for c in ['0.6', '0.75', '1']:
    cases.append(('dN', c, list(range(1, Nt[-1] + 1)), Nt))
Ns2 = [4, 8, 16, 24] if not a.quick else [4, 8]
for c in ['1.5', '0.5']:
    cases.append(('sq', c, [k * k for k in range(1, Ns2[-1] + 1)], Ns2))
    cases.append(('int', c, list(range(1, Ns2[-1] + 1)), Ns2))
dps0 = {'dN': 90, 'sq': 130, 'int': 90}
if a.quick: dps0 = {'dN': 50, 'sq': 60, 'int': 50}
for kind, c, fr, Ns in cases:
    dps = a.dps or dps0[kind]; t0 = time.time()
    o1, m1 = solve([mp.mpf(f) for f in fr], c, Ns, dps, 1.0)
    o2, m2 = solve([mp.mpf(f) for f in fr], c, Ns, dps + 30, 0.7)
    mp.mp.dps = 30; md = 0
    for N in Ns:
        d1, a1 = o1[N]; d2_, a2 = o2[N]
        md = max(md, abs(d1 - d2_) / abs(d2_), abs(a1 - a2) / a2)
        log('%s c=%s N=%d d_N=%s max|a|=%s (points %d/%d, %.0fs)' % (kind, c, N, mp.nstr(d2_, 6), mp.nstr(a2, 4), m1, m2, time.time() - t0))
        res[f'{kind},{c},{N}'] = [mp.nstr(d2_, 12), mp.nstr(a2, 8)]
    log('%s c=%s: max relative difference between the two runs (dps %d vs %d): %s' % (kind, c, dps, dps + 30, mp.nstr(md, 3)))
json.dump(res, open(a.out, 'w'), indent=1)
