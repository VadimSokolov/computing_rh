# maxwald_bondesson.py: the numbers of Theorem mw:thm:bondesson of Section sec:maxwald (ch/maxwald3.tex), "Bondesson classes
# containing GGC", and of the Remark that follows it. Added for the book in October 2026 from the checks of the authors' chapter,
# misc/reviews/nick-maxwald/ (review-3 items M11, M13, M15). With R(lambda) = xi(1/2)/xi(1/2 + sqrt(lambda)), w = v + 1/v,
#  order_ii:     the order pi H/1.01 - 1 of part (ii) for H = 3e12, against 9e12, and the angle bound tan(phi) <=
#                gamma/(gamma^2 - 1/4) < 1.01/gamma for gamma >= 10 (the largest ratio gamma^2/(gamma^2 - 1/4) is at gamma = 10);
#  remark_orders: pi H, the printed trade-off 1e-3 H^2/(2 log(H/gamma_1)) (review-3 M13: a square root is missing) and the order k
#                where the two regimes x >= (52.2 k + 100)/H^2 and x <= 1e-3/k of the derivative ladder meet; the remark prints
#                instead the order up to which the proof of Theorem rz:thm:nmono (ch/rz2.tex) works, from its two branches
#                y <= 0.05/n and 0.05 (H^2 - 1/4)/n > 2 n log(H + 1) + 6, which is (0.05 H^2/(2 log(H + 1)))^(1/2) to leading order;
#  model_pair:   the pair factor 1/((w - w0)(w - conj w0)) of a zero beta + i gamma off the line, w0 = -(y/u + u/y),
#                y = gamma^2 - (beta - 1/2)^2 - 2 i (beta - 1/2) gamma: the first order j at which
#                (-1)^j d^j/dw^j of it, = j! Im((w - w0)^(-j-1))/Im(w0), is not positive, for beta = 0.6, gamma = 100, w = 2 and a
#                range of u (review-3 M15 expects 1571 at u = 1 and 1570 as u tends to 0, against pi gamma/(2(beta - 1/2)) = 1570.8),
#                computed with the complex power at 60 digits and checked against pi/theta, theta = |arg(w - w0)|;
#  hcm_signs:    (-1)^j d^j/dw^j [R(uv) R(u/v)] for j <= 6, u in {0.5, 5, 50, 400}, w in {2.01, 3, 10}, from xi directly (60 digits);
#  hadamard:     R(uv) R(u/v) against the product over the first 1700 zeros of (1 + u w/gamma^2 + u^2/gamma^4)^(-1) with a smooth
#                tail from the density log(x/2 pi)/(2 pi): exponent -1 for zeros on the line (review-3 M11), against the squares
#                that the chapter's display writes.
# Run from a folder holding zeros_1700.txt: python3 maxwald_bondesson.py; writes maxwald_bondesson.json there.
import json
import os
import time

import mpmath as mp

mp.mp.dps = 60
HALF = mp.mpf(1) / 2


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


XI_HALF = xi(HALF)


def R(lam):
    return XI_HALF / xi(HALF + mp.sqrt(lam))


def pair(u, w):
    v = (w + mp.sqrt(w * w - 4)) / 2
    return R(u * v) * R(u / v)


def first_failing(beta, gamma, u, w, jmax=6000):
    b = mp.mpf(beta) - HALF
    g = mp.mpf(gamma)
    y = g ** 2 - b ** 2 - 2j * b * g
    w0 = -(y / u + u / y)
    z = mp.mpf(w) - w0
    theta = abs(mp.arg(z))
    out = {'u': float(u), 'w': float(w), 'Re_w0': float(w0.real), 'Im_w0': float(w0.imag),
           'pi_over_theta': float(mp.pi / theta) if theta != 0 else None}
    if abs(w0.imag) < mp.mpf(10) ** -40:
        out['first_failing_j'] = None
        return out
    zi = 1 / z
    P = zi
    out['first_failing_j'] = None
    for j in range(jmax + 1):
        if P.imag / w0.imag <= 0:
            out['first_failing_j'] = j
            break
        P *= zi
    return out


def main():
    t0 = time.time()
    res = {'working_precision': 'mpmath, 60 digits'}
    H = mp.mpf('3e12')
    order = mp.pi * H / mp.mpf('1.01') - 1
    g10 = mp.mpf(10) ** 2 / (mp.mpf(10) ** 2 - mp.mpf(1) / 4)
    res['order_ii'] = {'H': float(H), 'pi_H_over_1_01_minus_1': float(order), 'exceeds_9e12': bool(order > mp.mpf('9e12')),
                       'ratio_to_9e12': float(order / mp.mpf('9e12')),
                       'tan_phi_bound_at_gamma_10_gamma_over_gamma2_minus_quarter': float(10 / (100 - mp.mpf(1) / 4)),
                       'one_point_01_over_gamma_at_10': 0.101,
                       'max_ratio_gamma2_over_gamma2_minus_quarter_for_gamma_ge_10': float(g10),
                       'gamma_above_which_the_ratio_is_below_1_01': float(mp.sqrt(mp.mpf(1) / 4 / (1 - 1 / mp.mpf('1.01'))))}
    g1 = mp.mpf(open(next(p for p in ('zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt') if os.path.exists(p))).readline())
    L = 2 * mp.log(H / g1)
    printed = mp.mpf('1e-3') * H ** 2 / L
    k_root = (-100 + mp.sqrt(10 ** 4 + 4 * mp.mpf('52.2') * mp.mpf('1e-3') * H ** 2)) / (2 * mp.mpf('52.2'))
    res['remark_orders'] = {'pi_H': float(mp.pi * H), 'two_log_H_over_gamma_1': float(L),
                            'printed_1e-3_H2_over_two_log': float(printed), 'with_square_root': float(mp.sqrt(printed)),
                            'overlap_k_of_52_2k_plus_100_equal_1e-3_H2_over_k': float(k_root),
                            'pi_H_exceeds_order_with_root': bool(mp.pi * H > mp.sqrt(printed)),
                            'pi_H_exceeds_printed': bool(mp.pi * H > printed),
                            'nmono_order_leading_sqrt_0_05_H2_over_two_log_H_plus_1':
                                float(mp.sqrt(mp.mpf('0.05') * H ** 2 / (2 * mp.log(H + 1)))),
                            'nmono_order_exact_root_of_0_05_H2_minus_quarter_eq_2n2_log_H_plus_1_plus_6n':
                                float((-6 + mp.sqrt(36 + 8 * mp.log(H + 1) * mp.mpf('0.05') * (H ** 2 - mp.mpf(1) / 4)))
                                      / (4 * mp.log(H + 1)))}
    beta, gamma = 0.6, 100
    b = mp.mpf(beta) - HALF
    y = mp.mpf(gamma) ** 2 - b ** 2 - 2j * b * mp.mpf(gamma)
    phi = abs(mp.arg(y))
    mpair = {'beta': beta, 'gamma': gamma, 'y': [float(y.real), float(y.imag)], 'phi': float(phi),
             'pi_over_phi': float(mp.pi / phi), 'bound_pi_gamma_over_2_beta_minus_half': float(mp.pi * gamma / (2 * b)),
             'abs_y': float(abs(y))}
    byu = {}
    for u in (1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 0.5, 1, 2, 5, 10, 100, 1e3, 1e4, 1e5, 1e8):
        byu[f'u={u:g}'] = first_failing(beta, gamma, mp.mpf(u), 2)
    mpair['first_failing_order_w2_by_u'] = byu
    mpair['first_failing_order_u1_by_w'] = {f'w={w}': first_failing(beta, gamma, mp.mpf(1), w) for w in (2, 2.01, 3, 10)}
    mpair['u1_w2_first_failing_j'] = byu['u=1']['first_failing_j']
    mpair['u_to_0_first_failing_j'] = byu['u=1e-12']['first_failing_j']
    res['model_pair'] = mpair
    rows = {}
    allpos = True
    for u in (0.5, 5, 50, 400):
        for w in (2.01, 3, 10):
            d = mp.taylor(lambda x: pair(mp.mpf(u), x), mp.mpf(w), 6)
            vals = [(-1) ** j * d[j] * mp.factorial(j) for j in range(7)]
            ok = all(v > 0 for v in vals)
            allpos = allpos and ok
            rows[f'u={u:g},w={w:g}'] = {'values_j0_to_6': [mp.nstr(v, 6) for v in vals], 'all_positive': bool(ok)}
    res['hcm_signs'] = {'all_positive': bool(allpos), 'rows': rows}
    zs = [mp.mpf(x) for x in open(next(p for p in ('zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt') if os.path.exists(p))).read().split()]
    K = len(zs)
    Hc = zs[-1] + mp.pi / mp.log(zs[-1] / (2 * mp.pi))
    had = {}
    for u, w in ((5, 3), (50, 2.01), (0.5, 10)):
        u, w = mp.mpf(u), mp.mpf(w)
        logp = -sum(mp.log(1 + u * w / g ** 2 + u ** 2 / g ** 4) for g in zs)
        tail = -mp.quad(lambda x: mp.log(x / (2 * mp.pi)) / (2 * mp.pi) * mp.log(1 + u * w / x ** 2 + u ** 2 / x ** 4),
                        [Hc, 10 * Hc, mp.inf])
        direct = pair(u, w)
        prod1 = mp.exp(logp + tail)
        had[f'u={float(u):g},w={float(w):g}'] = {'direct_from_xi': mp.nstr(direct, 12), 'product_exponent_minus_1': mp.nstr(prod1, 12),
                                                 'relative_difference_exponent_minus_1': float(abs(prod1 / direct - 1)),
                                                 'product_exponent_minus_2': mp.nstr(prod1 ** 2, 12),
                                                 'relative_difference_exponent_minus_2': float(abs(prod1 ** 2 / direct - 1)),
                                                 'zeros_used': K}
    res['hadamard'] = had
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(res, open('maxwald_bondesson.json', 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'hcm_signs'}, indent=1)[:6000])
    print('hcm_signs all positive:', res['hcm_signs']['all_positive'])


if __name__ == '__main__':
    main()
