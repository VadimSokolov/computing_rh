# ht_euler32.py: the truncated Euler hitting times of Nick Polson's residue note (R:595 to 610, tab:euler32, Theorem euler32) and of Chapter
# 24 (H:318 to 337, the rate rho_4(P), rho_5(P)): log E exp(-s Z_P) at alpha = 1/2 for P = 1e3 ... 1e7, s = 0.5, 1, 2, the limits (2s)^(3/2)/6,
# the cumulants kappa_j(P), rho_j(P) = kappa_j/kappa_3^(j/3), the expansion (2s)^(3/2)/6 - rho_4 (2s)^2/24 + rho_5 (2s)^(5/2)/120 against the
# exact value, and the asymptotic 2^(-1/3) P^(-1/6) (log P)^(1/3) of rho_4.
# Conventions. L_P = sum_{p<=P} (log p) G_p with G_p geometric, P(G_p = m) = (1 - q_p) q_p^m, q_p = p^(-1/2). K_P(w) = log E exp(-w (L_P - E L_P))
# = sum_p [log(1 - q_p) - log(1 - q_p p^(-w)) + w (log p) q_p/(1 - q_p)]. Y_P has Laplace exponent -log E exp(-u Y_P) = K_P(sqrt(2u)); Z_P =
# (Y_P - kappa_2)/kappa_3^(2/3), c = kappa_3^(2/3), so log E exp(-s Z_P) = s kappa_2/c - K_P(sqrt(2 s/c)) (the j = 2 term of K_P cancels the
# centring). kappa_j(P) = sum_p (log p)^j Li_(1-j)(q_p) with Li_(-1)(q) = q/(1-q)^2, Li_(-2) = q(1+q)/(1-q)^3, Li_(-3) = q(1+4q+q^2)/(1-q)^4,
# Li_(-4) = q(1+11q+11q^2+q^3)/(1-q)^5 (the series sum_m m^(j-1) q^m).
# Method. Exact finite products over the primes up to P (numpy sieve), mpmath at dps 30 and dps 45; the maximum difference of the printed
# values between the two precisions is logged and written. The cancellation between the terms of K_P is about 1e2 at P = 1e7, so dps 30 is ample.
# Options: --quick (P up to 1e5, dps 20/30), --out FILE.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse
import numpy as np
from mpmath import mp, mpf, log, sqrt, nstr, power

ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/ht_euler32.json')
a = ap.parse_args()
LOG = open(a.out.replace('.json', '.log'), 'w')
def lg(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
Ps = [10 ** 3, 10 ** 4, 10 ** 5] if a.quick else [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
PM = max(Ps)
sieve = np.ones(PM + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(PM ** 0.5) + 1):
    if sieve[i]: sieve[i * i::i] = False
primes = np.nonzero(sieve)[0]
lg('primes up to %d: %d' % (PM, len(primes)))
Ss = [mpf('0.5'), mpf(1), mpf(2)]

def run(dps):
    mp.dps = dps
    out = {}
    # cumulative cumulants
    k = {j: mpf(0) for j in (2, 3, 4, 5)}; EL = mpf(0)
    cum = {}
    idx = 0
    for P in Ps:
        top = int(np.searchsorted(primes, P, side='right'))
        for p in primes[idx:top]:
            p = int(p); lp = log(p); q = 1 / sqrt(p); om = 1 - q
            EL += lp * q / om
            k[2] += lp ** 2 * q / om ** 2
            k[3] += lp ** 3 * q * (1 + q) / om ** 3
            k[4] += lp ** 4 * q * (1 + 4 * q + q * q) / om ** 4
            k[5] += lp ** 5 * q * (1 + 11 * q + 11 * q ** 2 + q ** 3) / om ** 5
        idx = top
        cum[P] = (dict(k), EL)
    for P in Ps:
        top = int(np.searchsorted(primes, P, side='right'))
        kk, EL = cum[P]
        c = kk[3] ** (mpf(2) / 3)
        ws = [sqrt(2 * s / c) for s in Ss]
        Ksum = [mpf(0) for _ in Ss]
        for p in primes[:top]:
            p = int(p); lp = log(p); q = 1 / sqrt(p); om = 1 - q
            l1 = log(om); eq = q * lp / om
            for t, w in enumerate(ws):
                Ksum[t] += l1 - log(1 - q * power(p, -w)) + w * eq
        vals = [s * kk[2] / c - Ksum[t] for t, s in enumerate(Ss)]
        rho4 = kk[4] / kk[3] ** (mpf(4) / 3); rho5 = kk[5] / kk[3] ** (mpf(5) / 3)
        asym = 2 ** (-mpf(1) / 3) * mpf(P) ** (-mpf(1) / 6) * log(P) ** (mpf(1) / 3)
        exp2 = [(2 * s) ** mpf('1.5') / 6 - rho4 * (2 * s) ** 2 / 24 for s in Ss]
        exp3 = [exp2[t] + rho5 * (2 * s) ** mpf('2.5') / 120 for t, s in enumerate(Ss)]
        out[P] = dict(vals=vals, rho4=rho4, rho5=rho5, asym=asym, exp2=exp2, exp3=exp3, k=kk)
    return out
t0 = time.time()
r1 = run(20 if a.quick else 30); r2 = run(30 if a.quick else 45)
mp.dps = 30
dmax = max(abs(r1[P]['vals'][t] - r2[P]['vals'][t]) for P in Ps for t in range(3))
dr = max(abs(r1[P]['rho4'] / r2[P]['rho4'] - 1) for P in Ps)
lg('max abs difference of the printed values between the two precisions: %s ; of rho_4 (rel): %s ; time %.0f s' % (nstr(dmax, 3), nstr(dr, 3), time.time() - t0))
R = {'precision_diff': nstr(dmax, 3)}
lim = [(2 * s) ** mpf('1.5') / 6 for s in Ss]
lg('limits (2s)^(3/2)/6 at s = 0.5, 1, 2:', [nstr(x, 6) for x in lim])
for t, s in enumerate(Ss):
    lg('tab:euler32 s=%s: exact log E exp(-s Z_P) at P = %s: %s' % (nstr(s, 2), [('1e%d' % round(np.log10(P))) for P in Ps], [nstr(r2[P]['vals'][t], 5) for P in Ps]))
    R['exact_s%s' % nstr(s, 2)] = [nstr(r2[P]['vals'][t], 6) for P in Ps]
for P in Ps:
    r = r2[P]
    lg('P=%d: rho_4 = %s (chapter H: 0.54, 0.37, 0.26, 0.19, 0.13), asymptotic 2^(-1/3) P^(-1/6) (log P)^(1/3) = %s, rho_5 = %s; kappa_2, kappa_3 = %s, %s' % (P, nstr(r['rho4'], 5), nstr(r['asym'], 5), nstr(r['rho5'], 5), nstr(r['k'][2], 8), nstr(r['k'][3], 8)))
    lg('   expansion with rho_4 only: %s ; with rho_4, rho_5: %s ; exact: %s ; error of the two term expansion %s, of the three term %s' % ([nstr(x, 5) for x in r['exp2']], [nstr(x, 5) for x in r['exp3']], [nstr(x, 5) for x in r['vals']], [nstr(r['exp2'][t] - r['vals'][t], 2) for t in range(3)], [nstr(r['exp3'][t] - r['vals'][t], 2) for t in range(3)]))
    R['P%d' % P] = dict(rho4=nstr(r['rho4'], 5), rho5=nstr(r['rho5'], 5), asym=nstr(r['asym'], 5), exp3=[nstr(x, 6) for x in r['exp3']], exact=[nstr(x, 6) for x in r['vals']], kappa2=nstr(r['k'][2], 10), kappa3=nstr(r['k'][3], 10))
json.dump(R, open(a.out, 'w'), indent=1)
