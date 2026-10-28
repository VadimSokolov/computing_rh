# wt_heat.py: the numerical claims of Chapter 26 (Weil Positivity as the Thorin Condition) about the heat trace Theta(x) = A(x) - P(x):
# tab:check (Theta from 200 zeros against A - P), tab:regimes (A, P, Theta, digits cancelled, share of the prime 2), the prime cutoff of
# the Prony section (2767, 403 primes, 438 prime powers, the cutoff table at x = 0.062 and 0.11, the 110 digit claim) and tab:primezeros
# (Prony with 16 samples x = 0.02 + 0.006 j and eight exponentials). Added for the book in October 2026, written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# Conventions. g_x(u) = (4 pi x)^(-1/2) exp(-u^2/4x); A(x) = e^(x/4) - (1/2) g_x(0) log(pi) + (1/4pi) int_R e^(-x r^2) Re psi(1/4 + i r/2) dr;
# P(x) = sum_{n >= 2} Lambda(n) n^(-1/2) g_x(log n) over the prime powers n; Theta(x) = sum_{gamma > 0} exp(-gamma^2 x) = A - P.
# Method. The integrand of A is even and analytic in |Im r| < 1/2 (poles of psi at r = +-i/2), so the trapezoid rule of step D on R has
# error of order exp(-pi/D) e^(x/4); every run computes the nodes at step D/2 and returns both the rule of step D (every other node) and of
# step D/2, and logs their difference. The prime sum is finite in practice: the terms with (log n)^2/4x > (dps + 12) ln 10 are below the
# working precision; the cutoff used for each sample is logged, and the effect of a looser cutoff (n <= 2767) is computed separately as
# the exact tail sum_{n > 2767} (to the rule's cutoff). Prony: theta_j = Theta(x_0 + j h), 8 exponentials: solve the 8 x 8 Hankel system for
# the coefficients of the characteristic polynomial, take its roots z_k = exp(-gamma_k^2 h) with mpmath polyroots (extra precision), and
# gamma_k = sqrt(-log z_k / h). Reference zeros: data/zeros_1700.txt (68 digits).
# Error control. dps and D are each changed (second run at dps2, D2) and the largest difference of Theta over the samples is logged; Theta
# from A - P is also compared with the sum over the first 1700 zeros (they agree to the 68 digits of the zeros file).
# Output: data/wt_heat.json (numbers as strings) and the log data/wt_heat.log.
# Usage: python3 wt_heat.py [--quick] [--out FILE] [--zeros FILE] [--procs N]; --quick runs at 30 digits in about 20 seconds on one
# core and writes to /tmp unless --out is given; the full case (110 and 130 digits) is a cluster job:
# bash misc/tools/hopper_run.sh -c 8 -t 120 -g wt_heat.json -g wt_heat.log code/wt_heat.py data/zeros_1700.txt -- --out wt_heat.json
import argparse
import json
import math
import sys
import time
from multiprocessing import Pool

from mpmath import mp, mpf, mpc, exp, log, sqrt, pi, euler, digamma, matrix, lu_solve, polyroots, floor, nstr

LOGF = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF:
        LOGF.write(s + '\n')
        LOGF.flush()


def repsi(args):
    r, dps = args
    mp.dps = dps
    return digamma(mpc(mpf(1) / 4, mpf(r) / 2)).real


def prime_powers(N):
    sieve = bytearray([1]) * (N + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    pp = []
    for p in range(2, N + 1):
        if sieve[p]:
            q = p
            while q <= N:
                pp.append((q, p))
                q *= p
    pp.sort()
    return pp


class Heat:
    def __init__(self, dps, D, xs, procs, Nmax):
        mp.dps = dps
        self.dps, self.D, self.xs = dps, mpf(D), [mpf(x) for x in xs]
        xmin = min(self.xs)
        R = sqrt((dps + 12) * mp.ln(10) / xmin)
        self.M = int(R / (self.D / 2)) + 2  # nodes j * D/2, j = 0..M
        t = time.time()
        with Pool(procs) as pool:
            vals = pool.map(repsi, [((self.D / 2) * j, dps) for j in range(self.M + 1)], chunksize=64)
        mp.dps = dps
        self.vals = vals
        say('  nodes', self.M + 1, 'dps', dps, 'D', nstr(self.D, 5), 'time', round(time.time() - t, 1))
        self.pp = prime_powers(Nmax)

    def A(self, x, fine=True):
        mp.dps = self.dps
        step = 1 if fine else 2
        h = self.D / 2 * step
        # int_R = 2 h (f(0)/2 + sum_{j>=1} f(jh))
        s0 = self.vals[0] / 2
        tot = s0 + sum(exp(-x * ((self.D / 2) * j) ** 2) * self.vals[j] for j in range(step, self.M + 1, step))
        I = 2 * h * tot
        g0 = 1 / sqrt(4 * pi * x)
        return exp(x / 4) - g0 * log(pi) / 2 + I / (4 * pi)

    def P(self, x, nmax=None):
        mp.dps = self.dps
        cut = exp(sqrt(4 * x * (self.dps + 12) * mp.ln(10)))
        if nmax is not None:
            cut = min(cut, mpf(nmax))
        s = mpf(0)
        g0 = 1 / sqrt(4 * pi * x)
        for q, p in self.pp:
            if q > cut:
                break
            u = log(q)
            s += log(p) / sqrt(q) * g0 * exp(-u * u / (4 * x))
        return s, int(cut)

    def P_tail(self, x, n0):
        mp.dps = self.dps
        cut = exp(sqrt(4 * x * (self.dps + 12) * mp.ln(10)))
        s = mpf(0)
        g0 = 1 / sqrt(4 * pi * x)
        for q, p in self.pp:
            if q <= n0:
                continue
            if q > cut:
                break
            u = log(q)
            s += log(p) / sqrt(q) * g0 * exp(-u * u / (4 * x))
        return s


def load_zeros(path, n):
    z = []
    for line in open(path):
        line = line.strip()
        if line:
            z.append(mpf(line.split()[-1]))
        if len(z) >= n:
            break
    return z


def prony(theta, K, h):
    n = len(theta)
    assert n == 2 * K
    Hm = matrix(K, K)
    rhs = matrix(K, 1)
    for j in range(K):
        for i in range(K):
            Hm[j, i] = theta[j + i]
        rhs[j] = -theta[j + K]
    a = lu_solve(Hm, rhs)
    coeffs = [mpf(1)] + [a[K - 1 - i] for i in range(K)]  # z^K + a_{K-1} z^{K-1} + ... + a_0
    roots = polyroots(coeffs, maxsteps=2000, extraprec=4 * mp.prec)
    roots = sorted(roots, key=lambda z: -abs(z))
    return [sqrt(-log(z.real if abs(z.imag) < mpf(10) ** (-mp.dps // 2) else z) / h) for z in roots]


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--zeros', default='data/zeros_1700.txt')
    ap.add_argument('--procs', type=int, default=1)
    a = ap.parse_args()
    out = a.out or ('/tmp/wt_heat.json' if a.quick else 'data/wt_heat.json')
    LOGF = open(out.rsplit('.', 1)[0] + '.log', 'w')
    res = {}
    if a.quick:
        runs = [(30, '0.04'), (34, '0.035')]
        dd = 30
    else:
        runs = [(110, '0.0115'), (130, '0.0097')]
        dd = 110
    mp.dps = 140
    zeros = load_zeros(a.zeros, 1700)
    say('wt_heat: runs (dps, D) =', runs, ' zeros read:', len(zeros))
    xs_check = ['0.005', '0.02', '0.05', '0.1', '0.2']
    xs_reg = ['0.002', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2', '0.3']
    xs_prony = [mpf('0.02') + mpf('0.006') * j for j in range(16)]
    allx = sorted(set([mpf(x) for x in xs_check + xs_reg] + xs_prony))
    # in the quick case only the Prony grid at 30 digits is too ill conditioned for 8 exponentials; quick uses K = 4 below
    heats = []
    for dps, D in runs:
        say('run dps', dps, 'D', D)
        Nmax = 80000 if not a.quick else 3000
        xx = allx if not a.quick else [mpf(x) for x in ['0.02', '0.05', '0.1']] + xs_prony[:8]
        H = Heat(dps, D, xx, a.procs, Nmax)
        heats.append((H, xx))
    H, xx = heats[0]
    H2, _ = heats[1]
    mp.dps = runs[0][0]
    # tab:check and tab:regimes
    say('\n== tab:check and tab:regimes (run 1; coarse vs fine step and run 2 differences)')
    table = {}
    for x in xx:
        mp.dps = runs[0][0]
        Af, Ac = H.A(x, True), H.A(x, False)
        Pv, cut = H.P(x)
        th = Af - Pv
        mp.dps = runs[1][0]
        A2 = H2.A(x, True)
        P2, cut2 = H2.P(x)
        th2 = A2 - P2
        mp.dps = runs[0][0]
        thz = sum(exp(-g * g * x) for g in zeros)
        p2 = log(2) / sqrt(2) * (1 / sqrt(4 * pi * x)) * exp(-log(2) ** 2 / (4 * x))
        row = dict(A=nstr(Af, 25), P=nstr(Pv, 25), Theta=nstr(th, 30), Theta_zeros1700=nstr(thz, 30), cut=cut,
                   dA_step=nstr(abs(Af - Ac), 3), dTheta_run=nstr(abs(th - th2), 3), dTheta_zeros=nstr(abs(th - thz), 3),
                   digits_cancelled=nstr(log(Af / th, 10), 5), share2=nstr(p2 / Pv, 5))
        table[nstr(x, 6)] = row
        say(nstr(x, 6), row)
    res['table'] = table
    # cutoff table
    say('\n== cutoff of the prime sum at dps', runs[0][0])
    cuts = {}
    for j in [0, 7, 15]:
        x = xs_prony[j]
        mp.dps = runs[0][0]
        cut = exp(sqrt(4 * x * dd * mp.ln(10)))
        t = H.P_tail(x, 2767)
        th = H.A(x) - H.P(x)[0]
        cuts[nstr(x, 6)] = dict(cut_110_rule=nstr(cut, 8), tail_after_2767=nstr(t, 5), Theta=nstr(th, 8),
                                digits_kept=nstr(-log(t / th, 10), 5), digits_absolute=nstr(-log(t, 10), 5))
        say(nstr(x, 6), cuts[nstr(x, 6)])
    res['cutoffs'] = cuts
    pp = prime_powers(2767)
    primes = [q for q, p in pp if q == p]
    res['counts_2767'] = dict(primes=len(primes), prime_powers=len(pp), last_prime=primes[-1])
    say('n <= 2767: primes', len(primes), 'prime powers', len(pp), 'largest prime', primes[-1])
    # Prony
    K = 8 if not a.quick else 4
    xp = xs_prony[:2 * K]
    h = mpf('0.006')
    res['prony'] = {}
    for label, nmax in [('full cutoff', None), ('n <= 2767', 2767)]:
        mp.dps = runs[0][0]
        theta = [H.A(x) - H.P(x, nmax)[0] for x in xp]
        mp.dps = runs[0][0] + 20
        gam = prony(theta, K, h)
        mp.dps = runs[0][0]
        say('\n== Prony,', label, ' K =', K, 'samples', 2 * K)
        rows = []
        for k in range(K):
            err = abs(gam[k] - zeros[k])
            rows.append(dict(k=k + 1, est=nstr(gam[k], 22), zero=nstr(zeros[k], 22), err=nstr(err, 3)))
            say(rows[-1])
        res['prony'][label] = rows
    json.dump(res, open(out, 'w'), indent=1)
    say('written', out)


if __name__ == '__main__':
    main()
