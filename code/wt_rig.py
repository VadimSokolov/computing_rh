# wt_rig.py: the numerical claims of Chapter 26 (Weil Positivity as the Thorin Condition) about the prime point process: the maximum of
# |psi(y) - y| / sqrt(y) on (1e3, 1e7] for the primes and for a Cramer model, tab:rig (variance of psi(x+h) - psi(x) - h over 2e5 uniform x in
# (1e6, 1e7), h = delta x, with the GUE prediction log(1/delta) + 1 - gamma_E - log 2 pi and the Poisson level log x - 1), tab:energy (mean of
# r(u)^2 on four windows of u = log y, r = (psi(e^u) - e^u) e^(-u/2), by exact integration between integers) and the constant
# sum_rho |rho|^(-2) = 2 + gamma_E - log 4 pi. Added for the book in October 2026, written from the mathematics (the collaborator's scripts
# were not read, run or copied).
# Conventions. psi(y) = sum_{n <= y} Lambda(n) (a right continuous step function); the Cramer model keeps each n >= 2 independently with
# probability 1/log n and gives it the weight log n. On [n, n+1) psi is constant psi_n, so int_n^(n+1) r^2 du = int (psi_n - y)^2 / y^2 dy
# is exact with the antiderivative -psi_n^2 / y - 2 psi_n log y + y (partial intervals at the window ends). The sup of |psi(y) - y| / sqrt(y)
# on [n, n+1) is attained at y = n or as y -> (n+1)^-, so the maximum over the integers of the values psi(n) - n and psi(n-1) - n (the left
# limit at n), divided by sqrt(n), is the supremum over real y.
# Variance estimators. Two, because "Var / h" with h = delta x varying along the sample is ambiguous: m1 = mean of D^2 / h and m2 = mean of
# D^2 over mean of h, D = psi(x+h) - psi(x) - h (the mean of D is also logged; the variance is about the mean). Standard errors are logged
# (batch means of 20 batches).
# Error control. The prime sums are exact up to float rounding (long double cumulative sums); the sampling error is the logged standard
# error; the sample is repeated with a second seed and the difference is logged. The Cramer model uses seed 20261023 for the tables and
# seeds 1 to 20 for the maximum.
# Output: data/wt_rig.json (numbers as strings) and the log data/wt_rig.log.
# Usage: python3 wt_rig.py [--quick] [--out FILE]; --quick uses N = 2e6 for everything (about 10 seconds) and writes to /tmp unless --out
# is given; the full case (N = 2.1e7) is a cluster job: bash misc/tools/hopper_run.sh -c 1 -t 60 -m 16G -g wt_rig.json -g wt_rig.log
# code/wt_rig.py -- --out wt_rig.json
import argparse
import json
import math
import time

import numpy as np

LOGF = None
EG = 0.5772156649015328606


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF:
        LOGF.write(s + '\n')
        LOGF.flush()


def lam_array(N):
    sieve = np.ones(N + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    lam = np.zeros(N + 1, dtype=np.longdouble)
    primes = np.nonzero(sieve)[0]
    lam[primes] = np.log(primes.astype(np.longdouble))
    for p in primes[primes <= int(N ** 0.5)]:
        q = int(p) * int(p)
        lp = np.log(np.longdouble(p))
        while q <= N:
            lam[q] = lp
            q *= int(p)
    return lam


def cramer_array(N, seed):
    rng = np.random.default_rng(seed)
    n = np.arange(N + 1, dtype=np.float64)
    n[:2] = 2.0
    keep = rng.random(N + 1) < 1.0 / np.log(n)
    keep[:2] = False
    lam = np.zeros(N + 1, dtype=np.longdouble)
    idx = np.nonzero(keep)[0]
    lam[idx] = np.log(idx.astype(np.longdouble))
    return lam


def psi_at(cum, y):
    return cum[np.floor(y).astype(np.int64)]


def rigidity(cum, deltas, rng, ns, xlo=1e6, xhi=1e7):
    out = {}
    x = rng.uniform(xlo, xhi, ns)
    for d in deltas:
        h = d * x
        D = (psi_at(cum, x + h) - psi_at(cum, x) - h).astype(np.float64)
        q = D * D / h
        # batch means for the standard error
        b = np.array_split(q, 20)
        bm = np.array([v.mean() for v in b])
        m2 = (D * D).mean() / h.mean()
        mu = D.mean()
        out[d] = dict(m1=float(q.mean()), se=float(bm.std(ddof=1) / math.sqrt(20)), m2=float(m2), meanD=float(mu))
    return out


def energy(cum, a, b, N):
    lo_y, hi_y = math.exp(a), math.exp(b)
    n0, n1 = int(math.floor(lo_y)), int(math.floor(hi_y))
    n = np.arange(n0, n1 + 1)
    ylo = np.maximum(n.astype(np.longdouble), np.longdouble(lo_y))
    yhi = np.minimum((n + 1).astype(np.longdouble), np.longdouble(hi_y))
    ps = cum[n]
    def F(y):
        return -ps * ps / y - 2 * ps * np.log(y) + y
    tot = (F(yhi) - F(ylo)).sum()
    return float(tot / (b - a))


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or ('/tmp/wt_rig.json' if a.quick else 'data/wt_rig.json')
    LOGF = open(out.rsplit('.', 1)[0] + '.log', 'w')
    N = 2_000_000 if a.quick else 21_000_000
    ns = 20000 if a.quick else 200000
    xmax_rig = 1e6 if a.quick else 1.1e7
    t = time.time()
    lam = lam_array(N)
    cum = np.cumsum(lam)
    say('sieve to', N, 'time', round(time.time() - t, 1), ' psi(N) =', float(cum[N]), ' psi(N)/N =', float(cum[N] / N))
    res = {}
    # constant
    c = 2 + EG - math.log(4 * math.pi)
    res['sum_rho_inv2'] = c
    say('2 + gamma_E - log 4 pi =', c)
    # max |psi - y| / sqrt y on (1e3, Ymax]
    Ymax = 10 ** 7 if not a.quick else 10 ** 6
    n = np.arange(1001, Ymax + 1)
    nn = n.astype(np.longdouble)
    v_at = np.abs(cum[n] - nn) / np.sqrt(nn)
    v_left = np.abs(cum[n - 1] - nn) / np.sqrt(nn)
    mx_int = float(v_at.max())
    mx_real = float(max(v_at.max(), v_left.max()))
    say('max |psi(n) - n|/sqrt n over integers in (1e3, %d]:' % Ymax, mx_int, ' sup over real y (left limits included):', mx_real,
        ' at n =', int(n[np.argmax(np.maximum(v_at, v_left))]))
    res['max_primes_integers'] = mx_int
    res['max_primes_real'] = mx_real
    # Cramer maxima
    mxs = []
    for seed in range(1, 21):
        cl = cramer_array(Ymax + 2, seed)
        cc = np.cumsum(cl)
        n = np.arange(1001, Ymax + 1)
        nn = n.astype(np.longdouble)
        m = float(max((np.abs(cc[n] - nn) / np.sqrt(nn)).max(), (np.abs(cc[n - 1] - nn) / np.sqrt(nn)).max()))
        mxs.append(m)
    mxs = np.array(mxs)
    say('Cramer maximum over 20 seeds: min %.2f median %.2f max %.2f mean %.2f' % (mxs.min(), np.median(mxs), mxs.max(), mxs.mean()))
    res['cramer_max'] = dict(min=float(mxs.min()), median=float(np.median(mxs)), max=float(mxs.max()), mean=float(mxs.mean()), values=[float(v) for v in mxs])
    # rigidity
    deltas = [1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 1e-1]
    xlo, xhi = (1e5, 1e6) if a.quick else (1e6, 1e7)
    rig = {}
    for seed in [1, 2]:
        rng = np.random.default_rng(seed)
        rig['primes_seed%d' % seed] = rigidity(cum, deltas, rng, ns, xlo, xhi)
    cl = cramer_array(N, 20261023)
    cc = np.cumsum(cl)
    rng = np.random.default_rng(1)
    rig['cramer'] = rigidity(cc, deltas, rng, ns, xlo, xhi)
    meanlogx = (1e7 * math.log(1e7) - 1e7 - 1e6 * math.log(1e6) + 1e6) / 9e6
    say('Poisson level mean(log x) - 1 over (1e6, 1e7):', meanlogx - 1)
    res['poisson'] = meanlogx - 1
    res['rig'] = {}
    for d in deltas:
        gue = math.log(1 / d) + 1 - EG - math.log(2 * math.pi)
        r1, r2, rc = rig['primes_seed1'][d], rig['primes_seed2'][d], rig['cramer'][d]
        say('delta %g: primes m1 %.3f (se %.3f; seed 2: %.3f)  m2 %.3f | GUE %.3f | Cramer m1 %.2f (se %.2f) m2 %.2f | meanD primes %.3f' %
            (d, r1['m1'], r1['se'], r2['m1'], r1['m2'], gue, rc['m1'], rc['se'], rc['m2'], r1['meanD']))
        res['rig'][str(d)] = dict(primes=r1, primes2=r2, cramer=rc, gue=gue)
    # energy
    res['energy'] = {}
    for (aa, bb) in [(5, 10), (8, 12), (10, 14), (12, 16.8)]:
        if math.exp(bb) > N - 1:
            continue
        ep = energy(cum, aa, bb, N)
        ec = energy(cc, aa, bb, N)
        say('window [%g, %g]: mean r^2 primes %.4f, Cramer %.3f' % (aa, bb, ep, ec))
        res['energy']['%g,%g' % (aa, bb)] = dict(primes=ep, cramer=ec)
    json.dump(res, open(out, 'w'), indent=1)
    say('written', out, 'time', round(time.time() - t, 1))


if __name__ == '__main__':
    main()
