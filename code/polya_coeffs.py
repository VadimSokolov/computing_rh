# polya_coeffs.py: the Taylor coefficients of F(s) = xi(1/2 + sqrt s)/xi(1/2), the Turan ratios and the Jensen polynomials of Section pol:sec:lp (Proposition pol:prop:lp and its numerics paragraph) and the shifted Jensen polynomials of the derivative problem (Problem pol:op:deriv), in Chapter 16 (Polya and the Zeros of Entire Functions). Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were not read, run or copied).
# Notation. The authors' draft of the chapter wrote gamma_k = k! a_k for the Polya Schur sequence and gamma_k for the ordinates of the zeros; here the sequence is called g_k = k! a_k and gamma_k is always an ordinate (data/zeros_1700.txt, 68 digits). F(s) = sum a_k s^k, a_k = c_k/c_0, and xi(1/2 + w) = sum c_k w^(2k) with c_k = (2/(2k)!) int_0^inf Phi(u) u^(2k) du, Phi(u) = 2 sum_{n>=1} (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)).
# Checks. (1) Xi(0) = xi(1/2) = 0.497120778188314 from the kernel, against the closed form -(1/8) pi^(-1/4) Gamma(1/4) zeta(1/2). (2) a_1 = 0.0231049931, a_2 = 2.48334054e-4, a_3 = 1.67435263e-6, a_4 = 8.03069742e-9. (3) The Turan ratios g_k^2/(g_(k-1) g_(k+1)) for k = 1..6 (1.0748, 1.0627, 1.0543, 1.0481, 1.0432, 1.0393) and their minimum for k <= 40 (above 1.0109). (4) Every Jensen polynomial J_d(x) = sum_{k<=d} binom(d,k) g_k x^k, 2 <= d <= 30, has only real negative zeros at 120 digits, and 30 times the largest zero of J_30 (chapter: -174.81) against -gamma_1^2 (chapter: -199.79). (5) Every shifted Jensen polynomial J^{d,n}(x) = sum_{j<=d} binom(d,j) g_(n+j) x^j with 2 <= d <= 15, 0 <= n <= 25, d + n <= 40 has only real negative zeros at 150 digits. (6) The zeros of xi'(1/2 + w) with 1 < Im w < 100 in the box |Re w| <= 0.6 (chapter, derivative problem: 29 of them, all on the imaginary axis): the extrema of Xi(t) on [1, 100] and the winding number of xi' around the box.
# Method. Composite Gauss Legendre panels on [0, 3] for the kernel (the kernel is below 1e-250 beyond 2.6), the moments u^(2k) taken at the same nodes, k up to 120. The real rootedness of a Jensen polynomial is decided by Sturm's theorem in exact integer arithmetic (python integers) for the polynomial whose coefficients are binom(d,k) g_k with g_k first rounded to 120 (150 for the shifted ones) significant digits and read as exact dyadic rationals: the number of distinct real roots is V(-inf) - V(+inf) of the Sturm sequence, it equals the degree exactly when all roots are real and simple, and all coefficients are positive so every real root is negative. The root nearest to the origin is then found by a geometric scan (ratio 1.02) and bisection at the working precision. For (6), xi'(s) = xi(s) (1/s + 1/(s - 1) - (log pi)/2 + psi(s/2)/2 + zeta'(s)/zeta(s)) and xi'(s) = -xi'(1 - s) for Re s < 1/2; the zeros on the axis are the extrema of Xi(t) = xi(1/2 + it), located on a grid of step 0.01 by sign changes of differences and refined by bisection on Xi'(t) = Re(i xi'(1/2 + it)), and the zeros inside the box are counted by the winding number of xi' along its boundary (steps of 0.05, the largest phase step is recorded).
# Error control. (a) Two quadrature layouts (panels of width 0.1 with 96 nodes and of width 0.2 with 192 nodes) are compared for every c_k; the largest relative difference is recorded. (b) The kernel value of xi(1/2) is compared with the closed form at the full working precision. (c) a_1 is compared with the sum of 1/gamma_k^2 over the 1700 tabulated zeros plus the density tail (log(gamma/2 pi) + 1)/(2 pi gamma) (about 9 digits expected), and F(-gamma_k^2) for k = 1, 2, 3 is evaluated from the Taylor series (it must vanish up to the truncation of the series at k = 120). (d) The Sturm routine is run on polynomials with known root counts before the real work, and for (6) the number of zeros of Xi on the axis must be 29 and the winding number must be an integer equal to the number of extrema. Nothing here is interval arithmetic: the digits are those of floating point at 175 digits, and the rounding of g_k to 120 or 150 digits is the stated precision of the Jensen test.
# Output: data/polya_coeffs.json (values as strings; a_k to 140 digits for k <= 120, the Turan ratios, per degree the Sturm count and the root nearest to 0, the shifted counts, the checks) and, beside it, the log data/polya_coeffs.log.
# Usage: python3 polya_coeffs.py [--quick] [--out FILE] [--zeros FILE]; --quick uses 70 digits, one quadrature layout pair of 48 and 96 node panels, k up to 60, rounding to 50 and 60 digits and the xi' box only up to height 25 (about 5 s, writes to /tmp unless --out is given); the full case (175 digits) is meant for the Hopper cluster: bash misc/tools/hopper_run.sh -g polya_coeffs.json -g polya_coeffs.log code/polya_coeffs.py data/zeros_1700.txt -- --out polya_coeffs.json
import argparse
import json
import math
import sys
import tempfile
import time
from math import comb
from pathlib import Path

from mpmath import mp, mpf, mpc
from mpmath.calculus.quadrature import GaussLegendre

T0 = time.time()
LOG = []


def log(*a):
    s = '[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(s, flush=True)
    LOG.append(s)


def fs(x, n=30):
    return mp.nstr(x, n)


def find_zeros(arg):
    here = Path(__file__).resolve().parent
    for c in ([arg] if arg else []) + [here / 'zeros_1700.txt', Path('data/zeros_1700.txt'), here.parent / 'data' / 'zeros_1700.txt']:
        if c and Path(c).is_file():
            return Path(c)
    raise SystemExit('zeros_1700.txt not found; pass --zeros')


def kernel_table(dps, npan, inv, deg):
    """Nodes u, weights and the kernel Phi(u) on [0, npan/inv] in panels of width 1/inv with 3*2^(deg-1) Gauss Legendre nodes each."""
    mp.dps = dps
    nodes = GaussLegendre(mp).calc_nodes(deg, mp.prec)
    h = mpf(1) / (2 * inv)
    pi = mp.pi
    pi2 = pi * pi
    cutoff = (dps + 80) * mp.log(10)
    us, ws, ph = [], [], []
    for p in range(npan):
        c = (2 * p + 1) * h
        for x, w in nodes:
            u = c + h * x
            e2 = mp.exp(2 * u)
            a = mp.exp(mpf(9) * u / 2)
            b = mp.exp(mpf(5) * u / 2)
            s = mpf(0)
            n = 1
            while pi * n * n * e2 < cutoff:
                s += (2 * pi2 * n ** 4 * a - 3 * pi * n * n * b) * mp.exp(-pi * n * n * e2)
                n += 1
            us.append(u)
            ws.append(h * w)
            ph.append(2 * s)
    return us, ws, ph


def coefficients(us, ws, ph, kmax):
    """c_k = (2/(2k)!) int Phi(u) u^(2k) du for k = 0..kmax."""
    u2 = [u * u for u in us]
    t = [w * p for w, p in zip(ws, ph)]
    c = []
    for k in range(kmax + 1):
        c.append(2 * mp.fsum(t) / mp.factorial(2 * k))
        t = [a * b for a, b in zip(t, u2)]
    return c


def rnd(x, nd):
    return mpf(mp.nstr(x, nd))


def poly_ints(binoms, coef):
    """Exact integer coefficients (lowest degree first) of sum binoms[k]*coef[k]*x^k, coef[k] mpf read as exact dyadic rationals."""
    parts = []
    for b, c in zip(binoms, coef):
        s, man, ex, bc = c._mpf_
        parts.append((b * (-man if s else man), ex))
    emin = min(e for v, e in parts if v != 0)
    return [v << (e - emin) if v != 0 else 0 for v, e in parts]


def _trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def _prem(a, b):
    """A positive multiple of the remainder of a on division by b (integer arithmetic, positive scalings only)."""
    a = list(a)
    db = len(b) - 1
    lb = b[-1]
    s = 1 if lb > 0 else -1
    while a and len(a) - 1 >= db:
        la = a[-1]
        k = len(a) - 1 - db
        a = [abs(lb) * c for c in a]
        for i, c in enumerate(b):
            a[i + k] -= s * la * c
        _trim(a)
        g = 0
        for c in a:
            g = math.gcd(g, c)
        if g > 1:
            a = [c // g for c in a]
    return a


def sturm_count(p):
    """Number of distinct real roots of the integer polynomial p (lowest degree first), by Sturm's theorem."""
    p = _trim(list(p))
    seq = [p, [i * p[i] for i in range(1, len(p))]]
    while True:
        r = _prem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(sg):
        return sum(1 for i in range(len(sg) - 1) if sg[i] != sg[i + 1])
    sp = [1 if q[-1] > 0 else -1 for q in seq]
    sm = [(1 if q[-1] > 0 else -1) * (-1) ** (len(q) - 1) for q in seq]
    return changes(sm) - changes(sp)


def closest_root(coef, dps):
    """The root of the polynomial (mpf coefficients, lowest first, positive constant term) nearest to 0, assumed negative."""
    mp.dps = dps

    def J(x):
        v = mpf(0)
        for c in reversed(coef):
            v = v * x + c
        return v
    t = mpf(1) / 10
    while J(-t) <= 0:
        t /= 2
    ratio = mpf('1.02')
    while t < mpf(10) ** 12:
        t2 = t * ratio
        if J(-t2) <= 0:
            break
        t = t2
    else:
        return None
    lo, hi = t, t2
    for _ in range(int(dps * 3.4) + 10):
        mid = (lo + hi) / 2
        if J(-mid) > 0:
            lo = mid
        else:
            hi = mid
    return -(lo + hi) / 2


def xi_s(s):
    return s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def xi_prime(s):
    if s.real < mpf(1) / 2:
        return -xi_prime(1 - s)
    return xi_s(s) * (1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, 1, 1) / mp.zeta(s))


def xi_prime_zeros(Hh, wd, step_axis, step_arg, dps):
    """Zeros of xi'(1/2 + w) with 1 < Im w < Hh in |Re w| <= wd: extrema of Xi on the axis and the winding number of xi' around the box."""
    mp.dps = dps
    half = mpf(1) / 2
    Hh = mpf(Hh)
    I = mpc(0, 1)
    na = int((Hh + 3 - 1) / mpf(step_axis))
    ts = [mpf(1) + k * mpf(step_axis) for k in range(na + 1)]
    Xi = [xi_s(mpc(half, t)).real for t in ts]
    nzero = sum(1 for k in range(na) if Xi[k] * Xi[k + 1] < 0 and ts[k] < Hh)
    d = [Xi[k + 1] - Xi[k] for k in range(na)]
    ext = []
    for k in range(1, na):
        if d[k - 1] * d[k] < 0:
            lo, hi = ts[k - 1], ts[k + 1]
            glo = (I * xi_prime(mpc(half, lo))).real
            for _ in range(45):
                mid = (lo + hi) / 2
                gm = (I * xi_prime(mpc(half, mid))).real
                if gm * glo > 0:
                    lo, glo = mid, gm
                else:
                    hi = mid
            ext.append((lo + hi) / 2)
    inside = [t for t in ext if t < Hh]
    wdm = mpf(wd)
    nv = int((Hh - 1) / mpf(step_arg))
    nh = int(2 * wdm / mpf(step_arg))
    path = [mpc(-wdm + k * (2 * wdm / nh), 1) for k in range(nh)]
    path += [mpc(wdm, 1 + k * (Hh - 1) / nv) for k in range(nv)]
    path += [mpc(wdm - k * (2 * wdm / nh), Hh) for k in range(nh)]
    path += [mpc(-wdm, Hh - k * (Hh - 1) / nv) for k in range(nv)]
    path.append(path[0])
    vals = [xi_prime(half + w) for w in path]
    tot = mpf(0)
    mx = mpf(0)
    for u, v in zip(vals[:-1], vals[1:]):
        inc = mp.arg(v / u)
        tot += inc
        mx = max(mx, abs(inc))
    return dict(height=int(Hh), half_width=fs(wdm, 4), xi_zeros_on_axis_1_to_H=nzero, extrema_count_1_to_H=len(inside), extrema=[fs(t, 10) for t in inside],
                first_extremum_beyond_H=fs(next((t for t in ext if t >= Hh), mpf(0)), 10), winding_number=fs(tot / (2 * mp.pi), 12), largest_phase_step=fs(mx, 5),
                boundary_points=len(path), chapter='29 zeros with 1 < Im w < 100 in |Re w| <= 0.6, all on the imaginary axis')


def main():
    ap = argparse.ArgumentParser(description='coefficients, Turan ratios and Jensen polynomials of F(s) = xi(1/2 + sqrt s)/xi(1/2)')
    ap.add_argument('--quick', action='store_true', help='70 digits, k up to 60; writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/polya_coeffs.json; quick: a /tmp file)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    a = ap.parse_args()
    quick = a.quick
    out = Path(a.out or (str(Path(tempfile.gettempdir()) / 'polya_coeffs_quick.json') if quick else 'data/polya_coeffs.json'))
    dps = 70 if quick else 175
    kmax = 60 if quick else 120
    layA, layB = ((15, 5, 5), (15, 5, 6)) if quick else ((30, 10, 6), (15, 5, 7))
    ndj, nds = (50, 60) if quick else (120, 150)
    R = dict(script='polya_coeffs.py', mode='quick' if quick else 'full')
    R['parameters'] = dict(dps=dps, kmax=kmax, layout_A=dict(panels=layA[0], inv_width=layA[1], nodes=3 * 2 ** (layA[2] - 1)),
                           layout_B=dict(panels=layB[0], inv_width=layB[1], nodes=3 * 2 ** (layB[2] - 1)),
                           jensen_digits=ndj, shifted_digits=nds)

    # (d) the Sturm routine on known cases
    tests = [([6, 11, 6, 1], 3), ([1, 0, 1], 0), ([1, 2, 1], 1), ([-1, 0, 0, 1], 1), ([2, -3, 0, 1], 2), ([24, 50, 35, 10, 1], 4), ([1, 1, 1, 1, 1], 0)]
    for p, want in tests:
        assert sturm_count(p) == want, (p, sturm_count(p), want)
    log('Sturm self test passed on', len(tests), 'polynomials')

    # (1) the coefficients from two quadrature layouts
    uA = kernel_table(dps, *layA)
    cA = coefficients(*uA, kmax)
    log('layout A done: nodes', len(uA[0]))
    uB = kernel_table(dps, *layB)
    cB = coefficients(*uB, kmax)
    log('layout B done: nodes', len(uB[0]))
    mp.dps = dps
    c = cB
    maxrel = max(abs(x - y) / abs(y) for x, y in zip(cA, cB))
    digits_agree = int(-mp.log10(maxrel)) if maxrel > 0 else dps
    log('largest relative difference of the two layouts over c_0..c_%d: %s (%d digits agree)' % (kmax, fs(maxrel, 5), digits_agree))
    xi0_closed = -mp.pi ** (-mpf(1) / 4) * mp.gamma(mpf(1) / 4) * mp.zeta(mpf(1) / 2) / 8
    d0 = abs(c[0] - xi0_closed) / xi0_closed
    log('xi(1/2) kernel   ', fs(c[0], 40))
    log('xi(1/2) closed   ', fs(xi0_closed, 40), ' relative difference', fs(d0, 5))
    ak = [x / c[0] for x in c]
    for k in range(1, 5):
        log('a_%d = %s' % (k, fs(ak[k], 20)))
    R['xi_half'] = dict(kernel=fs(c[0], 140), closed_form=fs(xi0_closed, 140), relative_difference=fs(d0, 5), chapter='0.497120778188314')
    R['layouts'] = dict(max_relative_difference=fs(maxrel, 5), digits_agree=digits_agree)
    R['a_k'] = [fs(x, 140) for x in ak]
    R['a_1_to_4'] = dict(computed=[fs(ak[k], 12) for k in range(1, 5)], chapter=['0.0231049931', '2.48334054e-4', '1.67435263e-6', '8.03069742e-9'])

    # (c) comparisons with the zeros
    zf = find_zeros(a.zeros)
    lines = [ln.strip() for ln in zf.read_text().split() if ln.strip()]
    gam = [mpf(s) for s in lines[:1700]]
    s1 = mp.fsum(1 / g ** 2 for g in gam)
    gN = gam[-1]
    tail = (mp.log(gN / (2 * mp.pi)) + 1) / (2 * mp.pi * gN)
    log('a_1 = %s; sum 1/gamma^2 over 1700 zeros %s, plus tail estimate %s, difference %s' % (fs(ak[1], 15), fs(s1, 15), fs(s1 + tail, 15), fs(ak[1] - s1 - tail, 4)))
    Fz = {}
    for k in (1, 2, 3):
        z = -gam[k - 1] ** 2
        terms = [ak[m] * z ** m for m in range(kmax + 1)]
        Fz[k] = (fs(mp.fsum(terms), 5), fs(max(abs(t) for t in terms), 5), fs(abs(terms[-1]), 5))
        log('F(-gamma_%d^2) from %d Taylor terms: %s (largest term %s, last term %s)' % (k, kmax + 1, Fz[k][0], Fz[k][1], Fz[k][2]))
    R['zeros_checks'] = dict(sum_inverse_square_1700=fs(s1, 20), tail_estimate=fs(tail, 8), a1_minus_sum_minus_tail=fs(ak[1] - s1 - tail, 4),
                             F_at_minus_gamma_squared={str(k): dict(value=v[0], largest_term=v[1], last_term=v[2]) for k, v in Fz.items()},
                             gamma_1=fs(gam[0], 20), minus_gamma_1_squared=fs(-gam[0] ** 2, 20))

    # (3) Turan ratios
    g = [mp.factorial(k) * ak[k] for k in range(kmax + 1)]
    ratios = [g[k] ** 2 / (g[k - 1] * g[k + 1]) for k in range(1, 41)]
    kmin = min(range(40), key=lambda i: ratios[i])
    log('Turan ratios k = 1..6:', ' '.join(fs(r, 8) for r in ratios[:6]))
    log('minimum over k = 1..40: %s at k = %d' % (fs(ratios[kmin], 10), kmin + 1))
    R['turan'] = dict(ratios_k1_to_40=[fs(r, 25) for r in ratios], first_six=[fs(r, 8) for r in ratios[:6]], minimum_k_le_40=fs(ratios[kmin], 12), argmin=kmin + 1,
                      decreasing=all(ratios[i + 1] < ratios[i] for i in range(39)), chapter_first_six=['1.0748', '1.0627', '1.0543', '1.0481', '1.0432', '1.0393'],
                      chapter_minimum='1.0109')

    # (4) Jensen polynomials
    jres = {}
    bad = []
    for d in range(2, 31):
        coef = [rnd(g[k], ndj) for k in range(d + 1)]
        ints = poly_ints([comb(d, k) for k in range(d + 1)], coef)
        cnt = sturm_count(ints)
        root = closest_root([comb(d, k) * coef[k] for k in range(d + 1)], ndj + 10)
        jres[d] = dict(distinct_real_roots=cnt, real_rooted_negative=(cnt == d), root_nearest_zero=fs(root, 40), d_times_root=fs(d * root, 25))
        if cnt != d:
            bad.append(d)
    mp.dps = dps
    log('Jensen polynomials J_d, 2 <= d <= 30, coefficients rounded to %d digits: real rooted for all d: %s (failures %s)' % (ndj, not bad, bad))
    log('d times the root nearest to 0: d = 2 %s, d = 10 %s, d = 20 %s, d = 30 %s; -gamma_1^2 = %s' % (
        jres[2]['d_times_root'][:12], jres[10]['d_times_root'][:12], jres[20]['d_times_root'][:12], jres[30]['d_times_root'][:12], fs(-gam[0] ** 2, 10)))
    R['jensen'] = dict(digits=ndj, all_real_rooted=not bad, failures=bad, per_degree={str(d): v for d, v in jres.items()},
                       J30_largest_zero_times_30=jres[30]['d_times_root'], chapter_J30_times_30='-174.81', chapter_minus_gamma1_sq='-199.79')

    # (5) shifted Jensen polynomials
    sres = []
    sbad = []
    for d in range(2, 16):
        for n in range(0, 26):
            if d + n > 40:
                continue
            coef = [rnd(g[n + j], nds) for j in range(d + 1)]
            cnt = sturm_count(poly_ints([comb(d, j) for j in range(d + 1)], coef))
            sres.append((d, n, cnt))
            if cnt != d:
                sbad.append((d, n, cnt))
    log('shifted Jensen polynomials J^{d,n}: %d polynomials (2 <= d <= 15, 0 <= n <= 25, d + n <= 40) at %d digits, real rooted: %d, failures %s' % (
        len(sres), nds, sum(1 for t in sres if t[2] == t[0]), sbad))
    R['shifted_jensen'] = dict(digits=nds, count=len(sres), real_rooted=sum(1 for t in sres if t[2] == t[0]), failures=sbad,
                               note='each polynomial: number of distinct real roots by Sturm equals the degree; all coefficients positive, so all roots are negative')
    # (6) the zeros of xi'(1/2 + w) in the box
    hh, stepa, stepw, dpsx = (25, 0.05, 0.05, 20) if quick else (100, 0.01, 0.05, 25)
    zx = xi_prime_zeros(hh, '0.6', stepa, stepw, dpsx)
    mp.dps = dps
    log('xi\'(1/2 + w), 1 < Im w < %d, |Re w| <= 0.6: zeros of Xi on the axis %d, extrema of Xi (= zeros of xi\' on the axis) %d, winding number of xi\' around the box %s (largest phase step %s); next extremum above the box at %s (chapter: %s)' % (
        hh, zx['xi_zeros_on_axis_1_to_H'], zx['extrema_count_1_to_H'], zx['winding_number'], zx['largest_phase_step'], zx['first_extremum_beyond_H'], '29 for H = 100'))
    log('last three extrema inside the box: %s' % ', '.join(zx['extrema'][-3:]))
    R['xi_prime_zeros'] = zx
    R['timing_seconds'] = fs(mpf(time.time() - T0), 5)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(R, indent=1) + '\n')
    log('wrote', out)
    logp = out.with_suffix('.log')
    logp.write_text('\n'.join(LOG) + '\n')


if __name__ == '__main__':
    main()
