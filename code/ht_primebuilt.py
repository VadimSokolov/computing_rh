# ht_primebuilt.py: the numerical claims of the prime built hitting times of Nick Polson's hitting times chapter (chapter lines H:657 to H:678, Theorem primebuilt):
# the density d_X(t) = pi^(-1) [ Re A(beta + i t) - Re sum_{n <= X} Lambda(n) n^(-beta - i t) (1 - n/X)^3 + Re X^(1 - beta - i t) W(1 - beta - i t) ], W(z) = 3! Gamma(z)/Gamma(z + 4) = 6/(z (z+1) (z+2) (z+3)),
# A(s) = xi'/xi(s) - zeta'/zeta(s) = 1/s + 1/(s-1) - (1/2) log pi + (1/2) psi(s/2), compared on 0 < t < 2000 with the exact Thorin density pi^(-1) Re xi'/xi(beta + i t).
# Claims: errors 0.26, 0.079, 0.027, 0.009 (beta = 1) and 0.62, 0.32, 0.15, 0.086 (beta = 0.8) for X = 1e2 .. 1e5, decay like X^(1/2 - beta), nonnegativity for X >= 1e3 at beta = 1, 0.9, 0.8, 0.7,
# minima 0.0073 against 0.0074 at beta = 1, and -0.16 near t = 1164 at beta = 0.6 (exact minimum 0.0015), and (by the same grid) the remark that the sum over the zeros is the whole error.
# Conventions. Error = max over the grid of |d_X - exact|; the grid is t = 0 (1e-9), step DT up to 2000 (the exact density is evaluated by mpmath at 15 digits, zeta'/zeta; the prime sums in double precision, Lambda over all prime powers n <= X,
# the n^(-i t) factors shared between the betas). The minimum of each density is taken over the grid and refined by a golden search around the grid minimum for the exact density and for d_X.
# Error control: the whole run is repeated on a second grid (step DT/2 on a subset of t, 0 < t < 400) and the changes of the errors and minima are logged; the exact density at 15 digits against 25 digits on 200 random points.
# Options: --quick (t < 200, step 0.1, beta 1 and 0.8, X 1e2 .. 1e3), --step DT (default 0.02), --out BASE (writes BASE.json and BASE.log). Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse, json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.special import digamma

HERE = Path(__file__).resolve().parent
T0 = time.time()
LOG = None


def log(*a):
    line = '[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(line, flush=True)
    if LOG:
        LOG.write(line + '\n'); LOG.flush()


def exact_chunk(args):
    beta, ts, dps = args
    mp.mp.dps = dps
    out = []
    for t in ts:
        s = mp.mpc(beta, t)
        v = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)
        out.append(float(mp.re(v) / mp.pi))
    return out


def exact(pool, beta, ts, dps=15):
    nch = 4 * pool._processes
    parts = pool.map(exact_chunk, [(beta, ts[i::nch], dps) for i in range(nch)])
    out = np.zeros(len(ts))
    for i in range(nch):
        out[i::nch] = parts[i]
    return out


def lambda_array(N):
    lam = np.zeros(N + 1)
    sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
    for p in range(2, int(N ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    for p in np.nonzero(sieve)[0]:
        q = int(p)
        while q <= N:
            lam[q] = math.log(p); q *= int(p)
    return lam


def dX(ts, betas, Xs, lam):
    N = int(max(Xs))
    n = np.nonzero(lam[:N + 1])[0]
    L = lam[n]; logn = np.log(n.astype(float))
    res = {(b, X): np.zeros(len(ts)) for b in betas for X in Xs}
    wt = {}
    for b in betas:
        for X in Xs:
            m = n <= X
            wt[(b, X)] = (m, L[m] * np.exp(-b * logn[m]) * (1 - n[m] / X) ** 3)
    CH = 400
    for i0 in range(0, len(ts), CH):
        t = ts[i0:i0 + CH]
        E = np.exp(-1j * np.outer(t, logn))
        for b in betas:
            s = b + 1j * t
            A = 1 / s + 1 / (s - 1) - math.log(math.pi) / 2 + digamma(s / 2) / 2
            for X in Xs:
                m, w = wt[(b, X)]
                P = E[:, m] @ w
                z = 1 - b - 1j * t
                W = 6 / (z * (z + 1) * (z + 2) * (z + 3))
                main = np.exp(z * math.log(X)) * W
                res[(b, X)][i0:i0 + CH] = (A.real - P.real + main.real) / math.pi
    return res


def golden_min(f, lo, hi):
    gr = (math.sqrt(5) - 1) / 2
    c = hi - gr * (hi - lo); d = lo + gr * (hi - lo); fc, fd = f(c), f(d)
    for _ in range(40):
        if fc < fd:
            hi, d, fd = d, c, fc; c = hi - gr * (hi - lo); fc = f(c)
        else:
            lo, c, fc = c, d, fd; d = lo + gr * (hi - lo); fd = f(d)
    return min(fc, fd), (lo + hi) / 2


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--step', type=float, default=0.02)
    ap.add_argument('--cpus', type=int, default=8)
    ap.add_argument('--out', default=str(HERE.parent / 'data' / 'ht_primebuilt'))
    a = ap.parse_args()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    LOG = open(a.out + '.log', 'w')
    if a.quick:
        tmax, step, betas, Xs = 200.0, 0.1, [1.0, 0.8], [100, 1000]
    else:
        tmax, step, betas, Xs = 2000.0, a.step, [1.0, 0.9, 0.8, 0.7, 0.6], [100, 1000, 10000, 100000]
    n = int(round(tmax / step))
    ts = np.array([1e-9] + [i * step for i in range(1, n + 1)])
    lam = lambda_array(max(Xs))
    R = {}
    with Pool(a.cpus) as pool:
        ex = {b: exact(pool, b, list(ts)) for b in betas}
        log('exact densities done on %d points' % len(ts))
        d = dX(ts, betas, Xs, lam)
        log('d_X done')
        for b in betas:
            errs = []
            for X in Xs:
                e = np.abs(d[(b, X)] - ex[b])
                i = int(np.argmax(e))
                mn = float(d[(b, X)].min()); imn = int(np.argmin(d[(b, X)]))
                R['%g_%d' % (b, X)] = dict(max_err=float(e[i]), at_t=float(ts[i]), min_dX=mn, argmin_t=float(ts[imn]), nonneg=bool(mn >= 0))
                errs.append(float(e[i]))
                log('beta=%g X=%d: max |d_X - exact| = %.4f at t=%.2f; min d_X = %.5f at t=%.2f; nonnegative on the grid: %s' % (b, X, e[i], ts[i], mn, ts[imn], mn >= 0))
            sl = [math.log(errs[k + 1] / errs[k]) / math.log(Xs[k + 1] / Xs[k]) for k in range(len(Xs) - 1)]
            mnx = float(ex[b].min()); imx = int(np.argmin(ex[b]))
            f = lambda t: exact_chunk((b, [t], 15))[0]
            gm, gt = golden_min(f, ts[max(imx - 1, 0)], ts[min(imx + 1, len(ts) - 1)])
            R['%g_exact' % b] = dict(min_grid=mnx, argmin_grid=float(ts[imx]), min_refined=gm, argmin_refined=gt, log_slopes_of_error_in_X=sl)
            log('beta=%g: exact minimum on the grid %.5f at t=%.2f (refined %.5f at %.3f); local slopes of log error against log X: %s (X^(1/2-beta) gives %.2f)' % (b, mnx, ts[imx], gm, gt, ['%.2f' % s for s in sl], 0.5 - b))
        # checks: finer subgrid and higher precision
        t2 = np.array([i * step / 2 + step / 4 for i in range(int(400 / step * 2))])
        ex2 = {b: exact(pool, b, list(t2)) for b in betas[:2]}
        d2 = dX(t2, betas[:2], Xs, lam)
        for b in betas[:2]:
            for X in Xs:
                e = float(np.max(np.abs(d2[(b, X)] - ex2[b])))
                R['fine_%g_%d' % (b, X)] = e
                log('check on the offset grid (0 < t < 400, step %.3g, shifted): beta=%g X=%d max error %.4f (compare max over the main grid in t < 400: %.4f)' % (step / 2, b, X, e, float(np.max(np.abs(d[(b, X)] - ex[b])[ts < 400]))))
        rng = np.random.default_rng(1)
        tt = list(rng.uniform(1, tmax, 100))
        a15 = np.array(exact(pool, 0.8, tt, 15)); a25 = np.array(exact(pool, 0.8, tt, 25))
        R['exact_15_vs_25'] = float(np.max(np.abs(a15 - a25)))
        log('exact density at 15 digits against 25 digits on 100 random points: max difference %.2e' % R['exact_15_vs_25'])
    json.dump(R, open(a.out + '.json', 'w'), indent=1)
    log('done')


if __name__ == '__main__':
    main()
