# code/hx_cramer.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N6 and decision R6, Proposition sm:prop:cramer): the primes as poles of the Poisson trace K_{1/2}, written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definitions. Under RH K_{1/2}(x) = sum_k e^{-gamma_k x} over the positive ordinates gamma_k (K_{1/2}(x) = sum_k e^{-x sqrt(a_k)}, a_k = gamma_k^2, Section sm:sec:poisson). F_n(eps) = eps K_{1/2}(eps + i log n) = eps sum_k e^{-gamma_k eps} n^(-i gamma_k). The explicit formula predicts the pole -Lambda(n)/(2 pi sqrt n) at x = i log n (prime part of K_{1/2}: -(1/pi) sum Lambda(n) n^(-1/2) x/(x^2 + log^2 n)), so F_n(0+) = -Lambda(n)/(2 pi sqrt n) for the prime powers n and 0 for the others (control n = 6).
# Method. All ordinates up to the height H = 10200 from Arb (python-flint acb.zeta_zeros, each zero a ball whose midpoint is used; the count N(H) is printed from the list, N(10^4) = 10142 is checked: the note's 10260 is wrong), 20 digits (the sum is converged at eps = 0.005 to e^{-gamma eps} <= e^{-51}, and the remainder is bounded by the integral of the density, printed), for n = 2, 3, 4, 5, 7, 8, 9 and 6 and eps = 0.05, 0.02, 0.01, 0.005, real and imaginary parts. Linear extrapolation to eps = 0: primary the line through the two smallest eps (0.01 and 0.005), 2 F(0.005) - F(0.01); the least squares line through all four is logged as well. [N6b] is the largest absolute difference (real and imaginary parts, primary extrapolation) from -Lambda(n)/(2 pi sqrt n) over the eight n; [N6a] is the height H of the list. Without python-flint a list of zeros from mpmath zetazero (slow) is used up to the --nzeros of the options.
# Error control: the list is compared with data/zeros_8000.txt (first 8000 zeros, relative 1e-15), the sums are recomputed at 30 digits for one (n, eps), and the tail of the sum beyond H is bounded by eps int_H^oo e^{-x eps} log(x/2 pi)/(2 pi) dx.
# Output: hx_cramer.json and hx_cramer.log; python3 hx_cramer.py [--quick] [--out FILE] [--H 10200]. Hopper: about a minute. Quick: H = 600.
import argparse
import json
import os
import time

import mpmath as mp

T0 = time.time()
LOGF = None


def say(*a):
    msg = '[%6.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(msg, flush=True)
    if LOGF:
        with open(LOGF, 'a') as f:
            f.write(msg + '\n')


def zeros_upto(H):
    try:
        from flint import acb, ctx
    except Exception:
        zf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_8000.txt')
        assert H <= 8000, 'python-flint is needed for H > 8000'
        return [x for x in (mp.mpf(t) for t in open(zf).read().split()) if x <= H]
    ctx.prec = 100
    # N(H) is about H/(2 pi) log(H/(2 pi e)); take a safe upper guess and cut
    Nest = int(H / 6.2831853 * max(1.0, float(mp.log(H / 17.0))) * 1.1) + 100
    zs = acb.zeta_zeros(1, Nest)
    out = [mp.mpf(z.imag.mid().str(30, radius=False)) for z in zs]
    return [x for x in out if x <= H]


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='hx_cramer.json')
    ap.add_argument('--H', type=float, default=10200.0)
    a = ap.parse_args()
    LOGF = os.path.splitext(a.out)[0] + '.log'
    open(LOGF, 'w').close()
    H = 600.0 if a.quick else a.H
    mp.mp.dps = 30
    g = zeros_upto(H)
    say('zeros up to H = %g: N(H) = %d, gamma_1 = %s, last = %s' % (H, len(g), mp.nstr(g[0], 15), mp.nstr(g[-1], 15)))
    res = {'script': 'hx_cramer.py', 'H': H, 'N_H': len(g)}
    n4 = sum(1 for x in g if x <= 1e4)
    say('N(10^4) from the list: %d (reconcile-P: 10142; the note: 10260)' % n4)
    res['N_1e4'] = n4
    zf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_8000.txt')
    if os.path.exists(zf):
        ref = [mp.mpf(x) for x in open(zf).read().split()]
        m = min(len(ref), len(g))
        d = max(abs(ref[i] - g[i]) / ref[i] for i in range(m))
        say('first %d zeros against data/zeros_8000.txt: max relative difference %s' % (m, mp.nstr(d, 3)))
        res['list_vs_zeros8000'] = mp.nstr(d, 4)
    Hh = g[-1]
    mp.mp.dps = 20
    gam = [mp.mpf(x) for x in g]
    epss = [mp.mpf(x) for x in ('0.05', '0.02', '0.01', '0.005')]
    primes = {2: mp.log(2), 3: mp.log(3), 4: mp.log(2), 5: mp.log(5), 7: mp.log(7), 8: mp.log(2), 9: mp.log(3), 6: mp.mpf(0)}
    # tail beyond H for eps = 0.005: bound
    eps0 = epss[-1]
    tailb = eps0 * mp.quad(lambda x: mp.exp(-x * eps0) * mp.log(x / (2 * mp.pi)) / (2 * mp.pi), [Hh, mp.inf])
    say('tail bound of the sum beyond H at eps = 0.005: %s' % mp.nstr(tailb, 3))
    res['tail_bound_eps0.005'] = mp.nstr(tailb, 4)
    res['table'] = {}
    worst = mp.mpf(0)
    worst_all = mp.mpf(0)
    for n, Lam in primes.items():
        pred = -Lam / (2 * mp.pi * mp.sqrt(n))
        ln = mp.log(n)
        F = []
        for e in epss:
            s = mp.mpc(0)
            for x in gam:
                s += mp.exp(-x * e) * mp.expjpi(-x * ln / mp.pi)
            F.append(e * s)
        ex = 2 * F[3] - F[2]
        # least squares line through (eps, F)
        xs = [float(e) for e in epss]
        mx = sum(xs) / 4
        ls = []
        for part in (lambda z: z.real, lambda z: z.imag):
            ys = [float(part(f)) for f in F]
            my = sum(ys) / 4
            b = sum((xs[i] - mx) * (ys[i] - my) for i in range(4)) / sum((xs[i] - mx) ** 2 for i in range(4))
            ls.append(my - b * mx)
        d_re, d_im = abs(ex.real - pred), abs(ex.imag)
        d = max(d_re, d_im)
        worst_all = max(worst_all, d)
        if n != 6:
            worst = max(worst, d)
        say('n = %d: predicted %s; F(0.05..0.005) re %s; im %s; extrapolated %s%+sj; error %s (least squares: %.3e%+.3ej)' % (
            n, mp.nstr(pred, 10), [mp.nstr(f.real, 8) for f in F], [mp.nstr(f.imag, 5) for f in F], mp.nstr(ex.real, 10), mp.nstr(ex.imag, 6), mp.nstr(d, 4), ls[0], ls[1]))
        res['table'][str(n)] = {'predicted': mp.nstr(pred, 12), 'F_re': [mp.nstr(f.real, 20) for f in F], 'F_im': [mp.nstr(f.imag, 20) for f in F],
                                'extrap_re': mp.nstr(ex.real, 12), 'extrap_im': mp.nstr(ex.imag, 8), 'error': mp.nstr(d, 6), 'lsq_re': ls[0], 'lsq_im': ls[1]}
    res['N6a'] = H
    res['N6b'] = mp.nstr(worst, 4)
    res['N6b_including_n6'] = mp.nstr(worst_all, 4)
    say('N6a: H = %g, N(H) = %d. N6b: largest |extrapolated - predicted| over n = 2,3,4,5,7,8,9: %s (including the control n = 6: %s)' % (H, len(g), mp.nstr(worst, 4), mp.nstr(worst_all, 4)))
    # precision check: n = 2, eps = 0.005 at 30 digits
    mp.mp.dps = 30
    e = mp.mpf('0.005')
    s = mp.mpc(0)
    for x in g:
        s += mp.exp(-x * e) * mp.expjpi(-x * mp.log(2) / mp.pi)
    f30 = e * s
    f20 = mp.mpc(res['table']['2']['F_re'][3]) + 1j * mp.mpf(res['table']['2']['F_im'][3])
    say('precision check n = 2, eps = 0.005: 30 digits %s, 20 digits %s, difference %s' % (mp.nstr(f30, 15), mp.nstr(f20, 15), mp.nstr(abs(f30 - f20), 3)))
    res['precision_check_diff'] = mp.nstr(abs(f30 - f20), 4)
    res['seconds'] = round(time.time() - T0, 1)
    json.dump(res, open(a.out, 'w'), indent=1)
    say('wrote', a.out)


if __name__ == '__main__':
    main()
