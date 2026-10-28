# cf_info.py: the information numbers of Section "Information rates" of the source note L (Chapter 34, prefix cf; Remark on information of the
# filter, Table tab:info, Proposition cf:prop:detect) and the BBLM coefficient. Group CF-XI, tasks T7 and T9. Written from the mathematics
# (the collaborator's scripts were not read, run or copied).
# Model of the source (L:884-889): one zero with free real part, level a = sigma - 1/2, signal Q_{sigma-beta}(t-gamma) = P_a(t-gamma) with
# a = sigma - beta; d/dbeta of it is -d/da P_a, with int (d_a P_a)^2 dx = 1/(4 pi a^3) (closed form checked here by mp.quad at 30 digits, and the
# derivative formula (x^2-a^2)/(pi (a^2+x^2)^2) by mp.diff). Gaussian errors of sd eps on a grid of spacing Dt: the Fisher information of the
# levels a_1..a_k is (4 pi eps^2 Dt)^(-1) S3_k, S3_k = sum a_l^(-3), sd(beta) = eps sqrt(4 pi Dt) S3_k^(-1/2); here eps = 0.01, Dt = 0.04 (the
# 401 points of [78,94]). Ratios of the sd after level 1 to that after level k: (S3_k/S3_1)^(1/2) (free real parts) and (S5_k/S5_1)^(1/4) (pairs,
# Proposition detect: the pair 1/2 +- delta against a double zero differs by eta d_a^2 P_a, eta = delta^2, with int (d_a^2 P_a)^2 = 3/(4 pi a^5)).
# Pair bound (L1 M6 of review-L1): the Fisher information of eta is (4 pi eps^2 Dt)^(-1) 3 S5_k, sd(eta) = its inverse square root, and the one sided
# 95% bound delta <= sqrt(1.6449 sd(eta)), a second order Gaussian approximation. The maximum over m zeros of independent half normal errors has
# quantile q: sd * z with Phi(z) = (1 + q^(1/m))/2. The levels of Table tab:info are sigma = 1, 0.9, 0.8, 0.7, 0.65, 0.6. T9: the coefficient
# 1/sqrt(12 Lambda_0), Lambda_0 = 1.57314 (the editors check the value against the BBLM paper) and the sizes log(t/2 pi)/sqrt(12 Lambda_0).
# Error control: all closed forms are exact rationals times pi or evaluated at 40 digits; the discrete grid sum of (d_a P_a)^2 over the 401 points
# of the source's grid is compared with the integral (the window truncation costs 2e-4 at a = 0.5 and the trapezoid aliasing, a fourth order pole at x = +-ia, about 4e-4 at a = 0.1).
# Usage: python3 -I cf_info.py [--quick] [--out FILE] (--quick has no effect: the run takes about 2 s); writes data/cf_info.json and cf_info.log.
import argparse
import json
import os
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 40
LEV = {'1': Fraction(1, 2), '0.9': Fraction(2, 5), '0.8': Fraction(3, 10), '0.7': Fraction(1, 5), '0.65': Fraction(3, 20), '0.6': Fraction(1, 10)}
EPS = mp.mpf('0.01')
DT = mp.mpf('0.04')


def P(a, x):
    return a / (mp.pi * (a * a + x * x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or os.path.join(here, '..', 'data', 'cf_info.json')
    logf = os.path.splitext(out)[0] + '.log'
    lines = []

    def log(s):
        print(s)
        lines.append(s)

    res = {}
    # closed forms
    cf = {}
    for av in ('0.1', '0.3', '0.5'):
        aa = mp.mpf(av)
        i1 = 2 * mp.quad(lambda x: ((x * x - aa * aa) / (mp.pi * (aa * aa + x * x) ** 2)) ** 2, [0, aa, 10 * aa, mp.inf])
        i2 = 2 * mp.quad(lambda x: (2 * aa * (aa * aa - 3 * x * x) / (mp.pi * (aa * aa + x * x) ** 3)) ** 2, [0, aa, 10 * aa, mp.inf])
        dP = mp.diff(lambda y: P(y, mp.mpf('0.37')), aa)
        d2P = mp.diff(lambda y: P(y, mp.mpf('0.37')), aa, 2)
        e1 = i1 - 1 / (4 * mp.pi * aa ** 3)
        e2 = i2 - 3 / (4 * mp.pi * aa ** 5)
        f1 = dP - (mp.mpf('0.37') ** 2 - aa * aa) / (mp.pi * (aa * aa + mp.mpf('0.37') ** 2) ** 2)
        f2 = d2P - 2 * aa * (aa * aa - 3 * mp.mpf('0.37') ** 2) / (mp.pi * (aa * aa + mp.mpf('0.37') ** 2) ** 3)
        log('a = %s: int (d_a P)^2 - 1/(4 pi a^3) = %.1e; int (d_a^2 P)^2 - 3/(4 pi a^5) = %.1e; derivative formulas %.1e, %.1e' % (av, e1, e2, f1, f2))
        cf[av] = {'e1': float(e1), 'e2': float(e2), 'd1': float(f1), 'd2': float(f2)}
    res['closed_forms'] = cf
    # discrete grid sum for a zero at the middle of [78,94] (gamma = 86): 401 points
    grid = {}
    for lv, av in LEV.items():
        af = mp.mpf(av.numerator) / av.denominator
        s = sum(((78 + i * DT - 86) ** 2 - af * af) ** 2 / (mp.pi ** 2 * (af * af + (78 + i * DT - 86) ** 2) ** 4) for i in range(401)) * DT
        grid[lv] = float(s / (1 / (4 * mp.pi * af ** 3)) - 1)
    log('relative difference of the grid sum of (d_a P_a)^2 Dt over 401 points to the integral, by level: ' + ', '.join('%s: %.1e' % kv for kv in grid.items()))
    res['grid_sum_vs_integral'] = grid
    names = list(LEV)
    S3 = []
    S5 = []
    s3 = Fraction(0)
    s5 = Fraction(0)
    for nme in names:
        s3 += LEV[nme] ** -3
        s5 += LEV[nme] ** -5
        S3.append(s3)
        S5.append(s5)
    log('S3 (sum a^-3) by level: ' + ', '.join(str(float(x)) for x in S3))
    log('S5 (sum a^-5) by level: ' + ', '.join(str(float(x)) for x in S5))
    free = [float(mp.sqrt(mp.mpf(S3[k].numerator) / S3[k].denominator / (mp.mpf(S3[0].numerator) / S3[0].denominator))) for k in range(6)]
    pair = [float((mp.mpf(S5[k].numerator) / S5[k].denominator / (mp.mpf(S5[0].numerator) / S5[0].denominator)) ** mp.mpf('0.25')) for k in range(6)]
    log('ratios of sd after level 1 to sd after level k, k = 2..6: free real parts ' + ', '.join('%.2f' % x for x in free[1:]) + '   (source 1.72, 2.75, 4.82, 7.76, 13.61)')
    log('   pairs ' + ', '.join('%.2f' % x for x in pair[1:]) + '   (source 1.42, 2.03, 3.27, 4.79, 7.77)')
    log('   unrounded free: ' + ', '.join('%.5f' % x for x in free[1:]) + '; pairs: ' + ', '.join('%.5f' % x for x in pair[1:]))
    res['ratios_free'] = free[1:]
    res['ratios_pair'] = pair[1:]
    sd = [EPS * mp.sqrt(4 * mp.pi * DT) / mp.sqrt(mp.mpf(x.numerator) / x.denominator) for x in S3]
    log('predicted sd of beta (x1e-4) by level 1, 0.9, 0.8, 0.7, 0.65, 0.6: ' + ', '.join('%.3f' % float(x * 10 ** 4) for x in sd) + '   (source 25.1, 14.6, 9.1, 5.2, 3.2, 1.84)')
    res['sd_beta'] = [float(x) for x in sd]
    hw = 1.959964 * sd[5]
    log('95%% half width after all six levels: %.3e (1.96 sd), source 3.6e-4' % float(hw))
    res['half_width_95'] = float(hw)

    def zq(q, m):
        return mp.sqrt(2) * mp.erfinv(q ** (mp.mpf(1) / m))   # Phi(z) = (1+p)/2 with p = q^(1/m): z = sqrt2 erfinv(p)

    for lab, k in (('after all six levels', 5), ('after sigma = 1 alone', 0)):
        log('max over four zeros, %s: sd %.3e, median %.3e, 95%% quantile %.3e' % (lab, float(sd[k]), float(sd[k] * zq(mp.mpf('0.5'), 4)), float(sd[k] * zq(mp.mpf('0.95'), 4))))
        res['max4_' + str(k)] = [float(sd[k]), float(sd[k] * zq(mp.mpf('0.5'), 4)), float(sd[k] * zq(mp.mpf('0.95'), 4))]
    z95 = mp.sqrt(2) * mp.erfinv(mp.mpf('0.9'))     # one sided 95%: Phi(z) = 0.95
    for lab, k in (('after all six levels', 5), ('after sigma = 1 alone', 0)):
        s5v = mp.mpf(S5[k].numerator) / S5[k].denominator
        sdeta = 1 / mp.sqrt(3 * s5v / (4 * mp.pi * EPS ** 2 * DT))
        bound = mp.sqrt(z95 * sdeta)
        log('pair bound %s: sd(eta) = %.4e, delta <= sqrt(%.4f sd(eta)) = %.4e (source 4.4e-3 / 0.035)' % (lab, float(sdeta), float(z95), float(bound)))
        res['pair_bound_' + str(k)] = float(bound)
    log('share of the levels 1 and 0.9 (sigma >= 0.9): %.4f of S3 and %.5f of S5 (source 1.6%% and 0.1%%); level 1 alone: %.4f and %.5f'
        % (float(S3[1] / S3[5]), float(S5[1] / S5[5]), float(S3[0] / S3[5]), float(S5[0] / S5[5])))
    log('filter against prediction at the first level: 27.1/25.1 = %.3f' % (27.1 / 25.1))
    # T9
    lam = mp.mpf('1.57314')
    kc = 1 / mp.sqrt(12 * lam)
    log('T9: 1/sqrt(12 Lambda_0) = %.5f (source 0.23016); 1/sqrt(12) = %.4f (the printed value 0.2887)' % (float(kc), float(1 / mp.sqrt(12))))
    sizes = {}
    for t in ('1e6', '1e7', '1e9', '1e10'):
        sizes[t] = float(mp.log(mp.mpf(t) / (2 * mp.pi)) * kc)
    log('T9: sizes log(t/2 pi)/sqrt(12 Lambda_0): ' + ', '.join('%s: %.3f' % kv for kv in sizes.items()) + '   (reconcile: 2.757 at 1e6, 4.877 at 1e10)')
    res['T9'] = {'coefficient': float(kc), 'sizes': sizes}
    json.dump(res, open(out, 'w'), indent=1)
    open(logf, 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
