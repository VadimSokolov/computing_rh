# code/hx_half.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N5 and decision R4, Proposition sm:prop:half): the cumulants of the positive half A of the reciprocal law under RH, written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definitions. gamma_k are the positive ordinates (the list data/zeros_8000.txt that xg_tilt.py uses). A = sum_k (E_k - 1)/gamma_k has cumulants kappa_m = (m-1)! sum_k gamma_k^(-m); E T = sum gamma^-2 (= Var A), Var T = sum gamma^-4. sd A = (E T)^(1/2); excess kurtosis of A = kappa_4/kappa_2^2 = 6 Var T/(E T)^2; skewness = kappa_3/kappa_2^(3/2) = 2 sum gamma^-3/(E T)^(3/2).
# Method. S_m = sum_{k <= 8000} gamma_k^-m at 50 digits (the list reaches H = gamma_8000 and N(H) = 8000 by the list itself); the tail sum_{gamma > H} gamma^-m = int_{(H, oo)} x^-m dN(x) = H^-m (M(H) - 8000) + int_H^oo x^-m M'(x) dx + m int_H^oo x^(-m-1) E(x) dx, where N = M + E, M(x) = (x/2 pi) log(x/(2 pi e)) + 7/8 (so M' = log(x/2 pi)/(2 pi)) and |E(x)| <= B(x) = 0.1038 log x + 0.2573 log log x + 9.3675 for x >= e (Hasanalizade, Shen and Wong 2022, J. Number Theory 235, as quoted in reconcile-P; to be checked against the paper). The first two terms are exact (closed forms), the last is bounded by m int_H^oo x^(-m-1) B(x) dx (mpmath quad). Rechecks: E T and Var T from the same list and tail against the book's p_1 = 0.0231049931154189707889 and p_2 = 0.0000371725992852696862 (Section sec:ch18:moments; the difference is the tail error, of the size of m int x^(-m-1) B); then sd A, excess kurtosis and skewness from the book's p_1, p_2 and from the list plus tail; the lattice control gamma_k = k (skewness 12 sqrt6 zeta(3)/pi^3, excess kurtosis 12/5).
# Output: hx_half.json and hx_half.log; python3 hx_half.py [--out FILE] (about 1 s). The number [N5] is the skewness to four digits.
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


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='hx_half.json')
    ap.add_argument('--zeros', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_8000.txt'))
    ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    LOGF = os.path.splitext(a.out)[0] + '.log'
    open(LOGF, 'w').close()
    mp.mp.dps = 50
    g = [mp.mpf(x) for x in open(a.zeros).read().split()]
    assert all(g[i] < g[i + 1] for i in range(len(g) - 1))
    H, n = g[-1], len(g)
    say('zero list: %d ordinates, gamma_1 = %s, H = gamma_%d = %s' % (n, mp.nstr(g[0], 15), n, mp.nstr(H, 20)))
    pi = mp.pi
    M = lambda x: x / (2 * pi) * mp.log(x / (2 * pi * mp.e)) + mp.mpf(7) / 8
    B = lambda x: mp.mpf('0.1038') * mp.log(x) + mp.mpf('0.2573') * mp.log(mp.log(x)) + mp.mpf('9.3675')
    say('M(H) = %s, N(H) - M(H) = %s, bound B(H) = %s' % (mp.nstr(M(H), 12), mp.nstr(n - M(H), 6), mp.nstr(B(H), 6)))

    def tailmain(m):
        # int_H^oo x^-m log(x/2 pi)/(2 pi) dx
        L = mp.log(H / (2 * pi))
        return (L / ((m - 1) * H ** (m - 1)) + 1 / ((m - 1) ** 2 * H ** (m - 1))) / (2 * pi)

    res = {'script': 'hx_half.py', 'H': mp.nstr(H, 20), 'count': n, 'N_H_minus_M_H': mp.nstr(n - M(H), 8)}
    S, tail, err = {}, {}, {}
    for m in (2, 3, 4):
        S[m] = sum(x ** (-m) for x in g)
        tail[m] = H ** (-m) * (M(H) - n) + tailmain(m)
        err[m] = m * mp.quad(lambda x: x ** (-m - 1) * B(x), [H, 10 * H, 1000 * H, mp.inf])
        say('m = %d: list sum %s, tail %s (error bound %s), total %s' % (m, mp.nstr(S[m], 22), mp.nstr(tail[m], 8), mp.nstr(err[m], 3), mp.nstr(S[m] + tail[m], 22)))
        res['sum_gamma^-%d' % m] = {'list': mp.nstr(S[m], 25), 'tail': mp.nstr(tail[m], 10), 'tail_error_bound': mp.nstr(err[m], 4), 'total': mp.nstr(S[m] + tail[m], 25)}
    p1 = mp.mpf('0.0231049931154189707889')
    p2 = mp.mpf('0.0000371725992852696862')
    say('E T: book %s, list + tail %s, difference %s (tail error bound %s)' % (mp.nstr(p1, 15), mp.nstr(S[2] + tail[2], 15), mp.nstr(S[2] + tail[2] - p1, 4), mp.nstr(err[2], 3)))
    say('Var T: book %s, list + tail %s, difference %s (tail error bound %s)' % (mp.nstr(p2, 15), mp.nstr(S[4] + tail[4], 15), mp.nstr(S[4] + tail[4] - p2, 4), mp.nstr(err[4], 3)))
    S3 = S[3] + tail[3]
    S3lo, S3hi = S3 - err[3], S3 + err[3]
    for name, ET, VT in (('book p1, p2', p1, p2), ('list + tail', S[2] + tail[2], S[4] + tail[4])):
        sd = mp.sqrt(ET)
        kurt = 6 * VT / ET ** 2
        skew = 2 * S3 / ET ** mp.mpf('1.5')
        say('%s: sd A = %s, excess kurtosis = %s, skewness = %s (range from the tail bound %s to %s)' % (name, mp.nstr(sd, 8), mp.nstr(kurt, 8), mp.nstr(skew, 8), mp.nstr(2 * S3lo / ET ** mp.mpf('1.5'), 8), mp.nstr(2 * S3hi / ET ** mp.mpf('1.5'), 8)))
        res[name] = {'sd': mp.nstr(sd, 10), 'excess_kurtosis': mp.nstr(kurt, 10), 'skewness': mp.nstr(skew, 10), 'skewness_low': mp.nstr(2 * S3lo / ET ** mp.mpf('1.5'), 10), 'skewness_high': mp.nstr(2 * S3hi / ET ** mp.mpf('1.5'), 10)}
    res['N5'] = mp.nstr(2 * S3 / p1 ** mp.mpf('1.5'), 4)
    say('N5: skewness = 2 sum gamma^-3 / (E T)^(3/2) = %s (four digits; the note has 0.415)' % res['N5'])
    sk = 12 * mp.sqrt(6) * mp.zeta(3) / pi ** 3
    say('lattice control gamma_k = k: skewness %s, excess kurtosis %s' % (mp.nstr(sk, 6), mp.nstr(mp.mpf(12) / 5, 6)))
    res['lattice'] = {'skewness': mp.nstr(sk, 8), 'excess_kurtosis': '2.4'}
    res['seconds'] = round(time.time() - T0, 2)
    json.dump(res, open(a.out, 'w'), indent=1)
    say('wrote', a.out)


if __name__ == '__main__':
    main()
