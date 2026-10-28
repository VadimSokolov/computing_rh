# code/hx_far.py (October 2026, from misc/reviews/nick-b4/reconcile-P.md, task N3, Corollary cp:cor:minflow and Problem cp:prob:rvm (a), decision R11): the far line hypothesis v_{0.2}(10, theta) >= (19/2) V_{0.2} for theta >= theta_1(0.2), written from the mathematics (the collaborator's scripts, notes and pickles were not read, run or copied).
# Definitions (Section cp:sec:minflow, Section sec:xg:tilt). v_t(alpha, theta) = Re (xi_t'/xi_t)(alpha + i theta) = Re[l + (t/2) l l' + zeta_t'/zeta_t](alpha + i theta) by xi_t = g e^{t l^2/4} zeta_t (Lemma cp:lem:shift); V_t = (log xi_t)''(1/2) = 2 E T_t = M_1/M_0 and sum_k theta_k(t)^(-4) = Var T_t = (M_1/M_0)^2/4 - (M_2/M_0)/12 with M_k/M_0 = 4^k m_k/m_0, m_k = int_0^oo u^(2k) e^{t u^2} Phi_dBN(u) du (Newton identities, Section sec:ch21:hankel; the same formulas as xg_tilt.py). Then theta_1(t) >= (sum theta_k^(-4))^(-1/4).
# (a) Arb enclosures of V_{0.2} and Var T_{0.2}: m_k by the rigorous integration of python-flint (acb.integral) over [0, 1.8], the series of Phi_dBN cut at n = 14 with the tail added to the ball radius from the lower bound of Re e^{4u} on the evaluation ball (tail <= (2 pi^2 |e^{9u}| + 3 pi |e^{5u}|) 1.01 15^4 e^{-pi L 15^2}, valid for L >= 0.3 where the terms n >= 14 decrease geometrically; the ball is made infinite if L < 0.3), and the part u > 1.8 bounded by 1e-800 (there e^{t u^2} u^4 Phi(u) <= e^{-pi e^{4u}/2}, whose integral is below e^{-2100}). Needs python-flint; without it only a non rigorous mpmath quadrature is logged.
# (b) The scan: v_{0.2}(10, theta) with zeta_t from the Gauss Hermite rule of hx_shift.py (K = 40 and 60 nodes, 40 digits), zeta_t' and l' by central differences with step 1e-8 (zeta_t at the mean of the two shifted values; the second step 2e-8 is compared on a sample), on theta = 12.7 + 0.05 j up to 10^4 (grid A, [N3a] = 10^4) and on a geometric grid of 200 points from 12.7 to 10^8 (grid B). The margin is v - (19/2) V_up with V_up the upper end of the enclosure of (a): its least value on the grid A and its location ([N3b], [N3c]) for K = 60, also restricted to theta >= theta_1 = 14.0768 (the first zero of xi_{0.2}), the least on grid B, the largest difference K = 40 against 60, and the three terms Re l, (t/2) Re(l l') and Re zeta_t'/zeta_t at theta = 14.0768 and 12.7. If the margin is not positive on the range this is reported and nothing is adjusted.
# (c) The largest |zeta_{0.2}(w) - 1| on the circles |w - 10 - i theta| = 5 at theta = 12.7, 50, 1000, 10^6 (256 equally spaced points each, K = 60, and K = 40 for the comparison) ([N3d]); the lower bound of Re l(w) on Re w >= 5, (1/2) psi(5/2) - (1/2) log pi, the mean 5 + (t/2) of it and the standard deviation sqrt(t) (the Gaussian P(Z <= -9.418) of R11 is computed from them).
# Output: hx_far.json and hx_far.log; python3 hx_far.py [--quick] [--out FILE] [--procs N]. Hopper: about 20 minutes on 48 cores. Quick: grid A to theta = 200, K = 40 only, a coarse circle, enclosures at 30 digits.
import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hx_shift as hs

T = '0.2'
ALPHA = '10'
EV = {}


def ev(dps=40):
    if dps not in EV:
        EV[dps] = hs.Ev(dps)
    return EV[dps]


def vcalc(theta, K, delta='1e-8', alpha=ALPHA):
    e = ev()
    e.act()
    s = e.cx(alpha, repr(float(theta)))
    d = e.real(delta)
    tt = e.real(T)
    sp, sm = s + d, s - d
    zp, zm = e.zt(sp, T, K), e.zt(sm, T, K)
    l = e.ell(s)
    dl = (e.ell(sp) - e.ell(sm)) / (2 * d)
    dz = (zp - zm) / (2 * d)
    z0 = (zp + zm) / 2
    t1 = e.to_c(l).real
    t2 = float(e.to_mp(tt / 2 * l * dl).real)
    t3 = e.to_c(dz / z0).real
    return (t1 + t2 + t3, t1, t2, t3)


def vchunk(args):
    ths, Ks = args
    return [[vcalc(th, K) for K in Ks] for th in ths]


def circ(args):
    theta, K, npts = args
    e = ev()
    best = 0.0
    arg = 0.0
    for j in range(npts):
        ph = 2 * math.pi * j / npts
        re = 10 + 5 * math.cos(ph)
        im = theta + 5 * math.sin(ph)
        w = e.cx(repr(re), repr(im))
        z = e.to_c(e.zt(w, T, K))
        d = abs(z - 1)
        if d > best:
            best, arg = d, ph
    return best, arg


def enclosures(quick):
    from flint import acb, arb, ctx
    ctx.prec = 200
    N = 14
    t = arb('0.2')
    pi = arb.pi()
    U = arb('1.8')

    def phi(x, analytic):
        e4, e5, e9 = (4 * x).exp(), (5 * x).exp(), (9 * x).exp()
        tot = acb(0)
        for n in range(1, N + 1):
            tot += (2 * pi ** 2 * n ** 4 * e9 - 3 * pi * n ** 2 * e5) * (-pi * n * n * e4).exp()
        L = e4.real.lower()
        if L < 0.3:
            return acb(arb(0, float('inf')), arb(0, float('inf')))
        Lb = arb(L)
        tail = (2 * pi ** 2 * abs(e9) + 3 * pi * abs(e5)) * arb('1.01') * arb(N + 1) ** 4 * (-pi * Lb * (N + 1) ** 2).exp()
        return tot + acb(arb(0, tail.upper()), arb(0, tail.upper()))

    ms = []
    for k in range(3):
        def f(x, analytic, k=k):
            return x ** (2 * k) * (t * x * x).exp() * phi(x, analytic)
        val = acb.integral(f, acb(0), acb(U), rel_tol=2 ** -150, abs_tol=2 ** -300)
        ms.append(val.real + arb(0, arb('1e-800')))
    m0, m1, m2 = ms
    V = 4 * m1 / m0
    VarT = V * V / 4 - arb(4) / 3 * m2 / m0
    ET = V / 2
    theta1_lb = arb(VarT.upper()) ** (arb(-1) / 4)
    return {'m0': m0, 'm1': m1, 'm2': m2, 'V': V, 'VarT': VarT, 'ET': ET, 'theta1_lower': theta1_lb}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='hx_far.json')
    ap.add_argument('--procs', type=int, default=0)
    a = ap.parse_args()
    hs.LOGF = os.path.splitext(a.out)[0] + '.log'
    open(hs.LOGF, 'w').close()
    say = hs.say
    procs = a.procs or (os.cpu_count() or 4)
    res = {'script': 'hx_far.py', 'quick': a.quick}
    say('hx_far.py, backend', 'flint' if hs.HAVE_FLINT else 'mpmath', 'quick' if a.quick else 'full', procs, 'processes')
    # (a)
    if hs.HAVE_FLINT:
        en = enclosures(a.quick)
        for k in ('m0', 'V', 'VarT', 'theta1_lower'):
            say('(a)', k, '=', en[k].str(30))
        V_up = float(en['V'].upper())
        Vs = en['V'].upper().str(30, radius=False)
        res['enclosures'] = {k: v.str(30) for k, v in en.items()}
        res['V_upper'] = Vs
        res['VarT_upper'] = en['VarT'].upper().str(30, radius=False)
        res['theta1_lower_bound'] = en['theta1_lower'].lower().str(15, radius=False)
        res['V_upper_float'] = V_up
        th1lb = float(en['theta1_lower'].lower())
        say('(a) V_0.2 <= %.12g, sum theta_k^-4 <= %.8g, theta_1(0.2) >= %.6f; (19/2)^2 = 90.25 <= (3/5) theta_lb^2 = %.4f: %s' % (V_up, float(en['VarT'].upper()), th1lb, 0.6 * th1lb ** 2, 90.25 <= 0.6 * th1lb ** 2))
        res['condition_95_ok_from_lower_bound'] = bool(90.25 <= 0.6 * th1lb ** 2)
    else:
        mp.mp.dps = 40
        t = mp.mpf('0.2')
        m = [mp.quad(lambda u, k=k: u ** (2 * k) * mp.exp(t * u * u) * hs.phi_dbn(u, 40), mp.linspace(0, 1.8, 10)) for k in range(3)]
        V = 4 * m[1] / m[0]
        VarT = V * V / 4 - mp.mpf(4) / 3 * m[2] / m[0]
        V_up = float(V)
        say('(a) not rigorous (no python-flint): V =', mp.nstr(V, 20), 'Var T =', mp.nstr(VarT, 15), 'theta_1 >=', mp.nstr(VarT ** (-0.25), 8))
        res['V_upper'] = mp.nstr(V, 20)
    Vup = V_up
    # (b)
    end = 200.0 if a.quick else 1.0e4
    n = int(round((end - 12.7) / 0.05))
    thA = [12.7 + 0.05 * j for j in range(n + 1)]
    thB = [12.7 * (1e8 / 12.7) ** (j / 199.0) for j in range(200)]
    Ks = (40,) if a.quick else (40, 60)
    cs = 100
    with Pool(procs) as pool:
        outA = []
        for r in pool.imap(vchunk, [(thA[i:i + cs], Ks) for i in range(0, len(thA), cs)]):
            outA.extend(r)
        say('grid A done:', len(thA), 'points')
        outB = list(pool.imap(vchunk, [([th], Ks) for th in thB]))
        outB = [r[0] for r in outB]
        say('grid B done:', len(thB), 'points')
        mA = 9.5 * Vup
        Kp = len(Ks) - 1

        def margin_stats(th, out, lo=0.0):
            best = None
            for x, o in zip(th, out):
                if x < lo:
                    continue
                m = o[Kp][0] - mA
                if best is None or m < best[0]:
                    best = (m, x)
            return best
        bA = margin_stats(thA, outA)
        bA1 = margin_stats(thA, outA, 14.0768)
        bB = margin_stats(thB, outB)
        say('(b) 9.5 V_up = %.12f' % mA)
        say('(b) N3b/N3c: least margin on grid A (K=%d): %.8f at theta = %.2f; restricted to theta >= 14.0768: %.8f at %.2f' % (Ks[Kp], bA[0], bA[1], bA1[0], bA1[1]))
        say('(b) least margin on grid B (200 points to 1e8): %.8f at theta = %.6g' % bB)
        dK = max(abs(o[0][0] - o[Kp][0]) for o in outA + outB) if len(Ks) > 1 else None
        say('(b) largest |v(K=40) - v(K=60)| over both grids:', dK)
        res.update({'N3a': end, 'margin_gridA': bA[0], 'margin_gridA_at': bA[1], 'margin_gridA_theta_ge_theta1': bA1[0], 'margin_gridA_theta_ge_theta1_at': bA1[1],
                    'margin_gridB': bB[0], 'margin_gridB_at': bB[1], 'max_K40_vs_K60': dK, 'nine_half_V_up': mA, 'margin_positive_on_range': bool(min(bA1[0], bB[0]) > 0 and bA[0] > 0)})
        # step check: v with step 2e-8 on a sample
        sample = thA[::max(1, len(thA) // 200)]
        dd = 0.0
        for th in sample:
            dd = max(dd, abs(vcalc(th, Ks[Kp], '2e-8')[0] - vcalc(th, Ks[Kp])[0]))
        say('(b) step check: largest |v(delta=2e-8) - v(delta=1e-8)| on %d points: %.3e' % (len(sample), dd))
        res['step_check'] = dd
        # three terms
        terms = {}
        for th in (14.0768, 12.7, 1000.0, 10000.0):
            r = vcalc(th, Ks[Kp])
            terms[str(th)] = {'v': r[0], 'Re_l': r[1], 'half_t_Re_l_dl': r[2], 'Re_dz_over_z': r[3]}
            say('(b) theta = %g: v = %.8f, Re l = %.6f, (t/2) Re(l l\') = %.6f, Re zeta_t\'/zeta_t = %.6f' % (th, r[0], r[1], r[2], r[3]))
        res['terms'] = terms
        # (c)
        npts = 16 if a.quick else 256
        circles = {}
        N3d = 0.0
        for th in ([12.7, 50.0, 1000.0] if a.quick else [12.7, 50.0, 1000.0, 1.0e6]):
            jobs = [(th, K, npts) for K in Ks]
            rr = pool.map(circ, jobs)
            circles[str(th)] = {'K%d' % K: {'max_abs_zeta_t_minus_1': r[0], 'at_phase': r[1]} for K, r in zip(Ks, rr)}
            N3d = max(N3d, rr[-1][0])
            say('(c) theta = %g: max |zeta_0.2 - 1| on the circle: %s' % (th, ['%.6e' % r[0] for r in rr]))
        res['circles'] = circles
        res['N3d'] = N3d
    mp.mp.dps = 30
    lower = mp.digamma(mp.mpf(5) / 2) / 2 - mp.log(mp.pi) / 2
    mean = 5 + mp.mpf('0.1') * lower
    # the real part of the shifted point is Re w + (t/2) Re l + sqrt(t/2) Z (variance t/2 = 0.1)
    sd = mp.sqrt(mp.mpf('0.1'))
    zq = (2 - (5 + mp.mpf('0.1') * lower)) / sd
    p = mp.ncdf(zq)
    say('(c) Re l(w) >= %s on Re w >= 5; the shifted real part reaches Re <= 2 when Z <= %s: probability %s' % (mp.nstr(lower, 8), mp.nstr(zq, 8), mp.nstr(p, 5)))
    res['Re_l_lower'] = mp.nstr(lower, 10)
    res['z_threshold'] = mp.nstr(zq, 10)
    res['prob_Z_below'] = mp.nstr(p, 6)
    res['seconds'] = round(time.time() - hs.T0, 1)
    json.dump(res, open(a.out, 'w'), indent=1, default=str)
    say('wrote', a.out)


if __name__ == '__main__':
    main()
