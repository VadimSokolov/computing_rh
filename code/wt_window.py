# wt_window.py: the numerical claims of Chapter 26 (Weil Positivity as the Thorin Condition) about the archimedean positivity window,
# Sections sec:regimes (the window, tab:window), sec:work (the scan, the sign structure, tab:pbp, the law log q / 2 gamma_1) and sec:round2
# (the window constants beta, delta, the second sign changes, tab:cert). Added for the book in October 2026, written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# Conventions. (-1)^n A^(n)(x) = (-1/4)^n e^(x/4) + (1/2pi) int_0^infty r^(2n) e^(-x r^2) (Re psi(1/4 + i r/2) - log pi) dr
# (eq:Aspec with lambda = r^2); g_x(u) = (4 pi x)^(-1/2) exp(-u^2/4x), d^n g_x/dx^n = g_x p_n(1/x) with p_0 = 1, p_(n+1)(y) = p_n(y)(-y/2 +
# c y^2) - y^2 p_n'(y), c = u^2/4; P_S(x) = sum over the prime powers m of the primes in S of Lambda(m) m^(-1/2) g_x(log m).
# x_n(S) is the first point (scan upward from 0.004 in steps of 0.002, then bisection to 1e-14) where (-1)^n (A - P_S)^(n) is not positive.
# Method. The integrand is even and analytic in |Im r| < 1/2, so the trapezoid rule of step D has error of order exp(-pi/D); each run
# evaluates with the nodes of step D and step D/2 (the latter is the value used) and a second run changes dps and D: the largest
# differences of every table value are logged. Table tab:cert: the lower bound of Proposition prop:cert (2) with the Rosser Schoenfeld
# constant 1.03883, on the grid 0.01 * 1.08^k, together with a check on a fine grid of u that f(y) = y^(-1/2) d^2 g_x(log y) is
# nonnegative and decreasing on [N_0, infinity); a grid point passes if the bound is positive and f passes. It is a floating point
# check, not an interval certificate.
# Output: data/wt_window.json (numbers as strings) and the log data/wt_window.log.
# Usage: python3 wt_window.py [--quick] [--out FILE]; --quick uses 22 digits, step 0.08 and a short scan (about 25 seconds) and writes to
# /tmp unless --out is given; the full case (35 digits, step 0.04, and a run at 45 digits, step 0.03) is a cluster job:
# bash misc/tools/hopper_run.sh -c 2 -t 120 -g wt_window.json -g wt_window.log code/wt_window.py -- --out wt_window.json
import argparse
import json
import sys
import time

from mpmath import mp, mpf, mpc, exp, log, sqrt, pi, digamma, gamma, findroot, diff, nstr, floor

LOGF = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF:
        LOGF.write(s + '\n')
        LOGF.flush()


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


def polyN(n, c):
    # coefficients (in y = 1/x) of p_n for g_x(u) with c = u^2/4
    p = [mpf(1)]
    for _ in range(n):
        q = [mpf(0)] * (len(p) + 2)
        for i, a in enumerate(p):
            q[i + 1] += -a / 2
            q[i + 2] += c * a
            q[i + 1] += 0
        # - y^2 p'(y)
        for i, a in enumerate(p):
            if i >= 1:
                q[i + 1] -= i * a
        p = q
    return p


def gder(n, x, u):
    c = u * u / 4
    g = exp(-c / x) / sqrt(4 * pi * x)
    y = 1 / x
    p = polyN(n, c)
    return g * sum(a * y ** i for i, a in enumerate(p))


class Arch:
    def __init__(self, dps, D, xmin):
        mp.dps = dps
        self.dps = dps
        self.D = mpf(D)
        self.xmin = mpf(xmin)
        R = sqrt(mpf(140) / self.xmin)
        self.M = int(R / (self.D / 2)) + 2
        lp = log(pi)
        self.val = [digamma(mpc(mpf(1) / 4, (self.D / 2) * j / 2)).real - lp for j in range(self.M + 1)]
        self.r2 = [((self.D / 2) * j) ** 2 for j in range(self.M + 1)]
        self.cache = {}

    def f(self, n, x, fine=True):
        key = (n, x, fine)
        if key in self.cache:
            return self.cache[key]
        mp.dps = self.dps
        step = 1 if fine else 2
        h = self.D / 2 * step
        Rmax2 = mpf(135) / x
        s = mpf(0)
        for j in range(step, self.M + 1, step):
            r2 = self.r2[j]
            if r2 > Rmax2:
                break
            s += r2 ** n * exp(-x * r2) * self.val[j]
        v = (-mpf(1) / 4) ** n * exp(x / 4) + h * s / (2 * pi) * 1  # int_0^inf = h * sum
        # (1/2pi) * int_0^inf  with int_0^inf = h*(f(0)/2 + sum), f(0) = 0 for n >= 1
        self.cache[key] = v
        return v


def Pder(n, x, S_pp):
    s = mpf(0)
    for q, p in S_pp:
        u = log(q)
        s += log(p) / sqrt(q) * gder(n, x, u)
    return s


def first_fail(fun, x0, step, xmax, tol=mpf(10) ** (-14)):
    x = x0
    prev = x
    while x <= xmax:
        if fun(x) <= 0:
            lo, hi = prev, x
            while hi - lo > tol:
                mid = (lo + hi) / 2
                if fun(mid) > 0:
                    lo = mid
                else:
                    hi = mid
            return (lo + hi) / 2
        prev = x
        x += step
    return None


def sign_changes(fun, x0, step, xmax, tol=mpf(10) ** (-12)):
    out = []
    x = x0
    s0 = fun(x) > 0
    prev = x
    while x <= xmax:
        x += step
        s1 = fun(x) > 0
        if s1 != s0:
            lo, hi = prev, x
            while hi - lo > tol:
                mid = (lo + hi) / 2
                if (fun(mid) > 0) == s0:
                    lo = mid
                else:
                    hi = mid
            out.append((lo + hi) / 2)
            s0 = s1
        prev = x
    return out


def run(dps, D, quick):
    mp.dps = dps
    out = {}
    x0 = mpf('0.004') if not quick else mpf('0.012')
    A = Arch(dps, D, x0)
    lp = log(pi)
    # r*, lambda*
    rstar = findroot(lambda r: digamma(mpc(mpf(1) / 4, r / 2)).real - lp, mpf('6.29'))
    lam = rstar ** 2
    out['rstar'] = nstr(rstar, 15)
    out['lambda_star'] = nstr(lam, 15)
    say('r* =', nstr(rstar, 15), ' lambda* =', nstr(lam, 15), ' 4pi^2 =', nstr(4 * pi ** 2, 10), ' psi(1/4) =', nstr(digamma(mpf(1) / 4), 8))
    # window x_n and law
    xn = {}
    for n in range(1, 7):
        fn = lambda x, n=n: A.f(n, x, True)
        xn[n] = first_fail(fn, x0, mpf('0.002') if not quick else mpf('0.004'), mpf(1))
        fc = lambda x, n=n: A.f(n, x, False)
        xc = first_fail(fc, x0, mpf('0.002') if not quick else mpf('0.004'), mpf(1))
        say('x_%d = %s  (n/r*^2 = %s)  step D vs D/2 difference %s' % (n, nstr(xn[n], 12), nstr(n / lam, 8), nstr(abs(xn[n] - xc), 3)))
    out['x_n'] = {str(n): nstr(xn[n], 14) for n in xn}
    out['n_over_lambda'] = {str(n): nstr(n / lam, 10) for n in xn}
    x1 = xn[1]
    say('x_1 =', nstr(x1, 8), '1/(4pi^2) =', nstr(1 / (4 * pi ** 2), 8), ' exp(-(log2)^2/4x_1) =', nstr(exp(-log(2) ** 2 / (4 * x1)), 5),
        ' (log2)^2/4 =', nstr(log(2) ** 2 / 4, 5), ' gamma_1^2 = 199.79 ; (gamma_1^2+1/4)/ln10 =', nstr((mpf('14.134725141734693790457') ** 2 + mpf(1) / 4) / mp.ln(10), 6))
    out['x1_extra'] = dict(inv4pi2=nstr(1 / (4 * pi ** 2), 8), frechet=nstr(exp(-log(2) ** 2 / (4 * x1)), 6),
                           digits_per_x=nstr((mpf('14.134725141734693790457') ** 2 + mpf(1) / 4) / mp.ln(10), 6),
                           gamma1sq=nstr(mpf('14.134725141734693790457') ** 2, 7))
    # offsets and delta, beta
    w = lambda l: (digamma(mpc(mpf(1) / 4, sqrt(l) / 2)).real - lp) / (4 * pi * sqrt(l))
    w1 = diff(w, lam, 1)
    w2 = diff(w, lam, 2)
    beta = (1 + lam * w2 / (2 * w1)) / lam
    c = (digamma(mpf(1) / 4) - lp) / (4 * pi)
    say('beta =', nstr(beta, 8), ' c =', nstr(c, 8), " w'(lam*) =", nstr(w1, 8), " w''(lam*) =", nstr(w2, 8))
    out['beta'] = nstr(beta, 8)
    out['c'] = nstr(c, 8)
    off = {n: xn[n] - n / lam for n in xn}
    say('offsets x_n - n/lambda* (1e-4):', [nstr(off[n] * 10 ** 4, 4) for n in range(1, 7)])
    say('n (x_n - n/lambda* - beta) (1e-3):', [nstr(n * (off[n] - beta) * 10 ** 3, 4) for n in range(1, 7)])
    out['offsets_1e4'] = [nstr(off[n] * 10 ** 4, 5) for n in range(1, 7)]
    out['delta_1e3'] = [nstr(n * (off[n] - beta) * 10 ** 3, 5) for n in range(1, 7)]
    # scans and second sign changes
    xmax = mpf(7) if not quick else mpf(1.5)
    sec = {}
    for n in [1, 2, 3, 4, 5, 6]:
        fn = lambda x, n=n: A.f(n, x, True)
        st = mpf('0.04') if not quick else mpf('0.1')
        ch = sign_changes(fn, x0 + mpf('0.001'), mpf('0.001') if n == 1 else mpf('0.002'), mpf('0.3'))
        ch2 = sign_changes(fn, mpf('0.3'), st, xmax)
        sec[n] = [nstr(z, 9) for z in ch + ch2]
        say('sign changes of (-1)^%d A^(%d) on [0.005, %s]: %s' % (n, n, nstr(xmax, 3), sec[n]))
    out['sign_changes'] = sec
    # balance for second sign change
    bal = {}
    for n in [2, 4, 6]:
        r = findroot(lambda x: 4 ** (-n) * exp(x / 4) - abs(c) * gamma(n + mpf(1) / 2) * x ** (-n - mpf(1) / 2), mpf(2) * n / 3 + 1)
        bal[n] = nstr(r, 6)
        say('balance n = %d: %s' % (n, bal[n]))
    out['balance'] = bal
    # tab:pbp
    laws = {}
    gam1 = mpf('14.134725141734693790457')
    sets = [('none', 1), ('2', 2), ('2,3', 3), ('<=5', 5), ('<=7', 7), ('<=13', 13), ('<=31', 31), ('<=97', 97)]
    allpp = prime_powers(400)
    # powers of the primes in S: all powers up to 1e6
    def powers(pmax):
        res = []
        for q, p in prime_powers(pmax):
            pass
        return res
    pbp = {}
    for label, pmax in sets:
        S_primes = [p for p in range(2, pmax + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
        if label == '2':
            S_primes = [2]
        elif label == '2,3':
            S_primes = [2, 3]
        if label == 'none':
            S_primes = []
        spp = []
        for p in S_primes:
            q = p
            while q <= 10 ** 5:
                spp.append((q, p))
                q *= p
        # first omitted prime power: smallest prime power not generated by S
        nextq = None
        for q, p in allpp:
            if p not in S_primes:
                nextq = q
                break
        row = {'q': int(nextq), 'law': nstr(log(nextq) / (2 * gam1), 6)}
        for n in [1, 2, 3]:
            def fun(x, n=n, spp=spp):
                return A.f(n, x, True) - (-1) ** n * Pder(n, x, spp)
            xx = first_fail(fun, x0, mpf('0.004') if quick else mpf('0.002'), mpf('0.8') if not quick else mpf('0.25'), tol=mpf(10) ** (-10))
            row['x%d' % n] = nstr(xx, 6) if xx is not None else '> 0.8'
        say('S =', label, row)
        pbp[label] = row
        if row['x1'][0] != '>':
            e = abs(mpf(row['x1']) - log(nextq) / (2 * gam1)) / mpf(row['x1']) * 100
            row['err_pct'] = nstr(e, 3)
            row['err_pct_rel_law'] = nstr(abs(mpf(row['x1']) - log(nextq) / (2 * gam1)) / (log(nextq) / (2 * gam1)) * 100, 3)
            say('   law error (percent of x_1(S) / of the law):', row['err_pct'], row['err_pct_rel_law'])
    out['pbp'] = pbp
    # tab:cert
    cert = {}
    for N0 in ([10, 100, 1000] if not quick else [10]):
        pp = prime_powers(N0)
        grid = [mpf('0.01') * mpf('1.08') ** k for k in range(0, 60)]
        best = None
        failed = None
        for x in grid:
            if x > mpf('0.4') or (quick and x > mpf('0.1')):
                break
            lead = A.f(2, x, True) - Pder(2, x, pp)
            u0 = log(N0)
            fN0 = N0 ** (-mpf(1) / 2) * gder(2, x, u0)
            tail = mp.quad(lambda u: exp(u / 2) * gder(2, x, u), [u0, u0 + 2, u0 + 6, u0 + 20, mp.inf])
            Lb = lead - mpf('1.03883') * (N0 * fN0 + tail)
            # f nonnegative decreasing on [N0, infinity): check on a grid of u
            ok = True
            prev = None
            u = u0
            while u < u0 + (30 if not quick else 8):
                fv = exp(-u / 2) * gder(2, x, u)
                if fv < 0 or (prev is not None and fv > prev):
                    ok = False
                    break
                prev = fv
                u += mpf('0.01')
            if Lb > 0 and ok:
                best = x
            else:
                failed = (x, Lb > 0, ok)
                break
        cert[str(N0)] = dict(certified_up_to=nstr(best, 5) if best is not None else None,
                             first_failed_grid_point=nstr(failed[0], 5) if failed else None,
                             bound_positive_at_failure=bool(failed[1]) if failed else None, f_ok_at_failure=bool(failed[2]) if failed else None)
        say('N0 =', N0, cert[str(N0)])
    out['cert'] = cert
    return out


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or ('/tmp/wt_window.json' if a.quick else 'data/wt_window.json')
    LOGF = open(out.rsplit('.', 1)[0] + '.log', 'w')
    runs = [(22, '0.08')] if a.quick else [(35, '0.04'), (45, '0.03')]
    res = {}
    for dps, D in runs:
        say('===== run dps', dps, 'D', D)
        t = time.time()
        res['dps%d' % dps] = run(dps, D, a.quick)
        say('time', round(time.time() - t, 1))
    if len(runs) == 2:
        r1, r2 = res['dps35'], res['dps45']
        d = max(abs(mpf(r1['x_n'][k]) - mpf(r2['x_n'][k])) for k in r1['x_n'])
        say('max difference of x_n between the runs:', nstr(d, 3))
        res['diff_x_n'] = nstr(d, 3)
        for k in r1['pbp']:
            for n in ['x1', 'x2', 'x3']:
                if r1['pbp'][k][n][0] != '>' and r2['pbp'][k][n][0] != '>':
                    dd = abs(mpf(r1['pbp'][k][n]) - mpf(r2['pbp'][k][n]))
                    if dd > mpf(10) ** (-5):
                        say('pbp run difference', k, n, nstr(dd, 3))
    json.dump(res, open(out, 'w'), indent=1)
    say('written', out)


if __name__ == '__main__':
    main()
