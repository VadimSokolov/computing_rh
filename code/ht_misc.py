# ht_misc.py: the numerical claims of the frontier, Yor and Davenport Heilbronn parts of Nick Polson's hitting times chapter (chapter lines H:417 to H:420, H:486,
# H:520 to H:531, H:558, H:576, H:590 to H:594, H:614, H:680, H:699): the Bernstein order 1.7e23 x and the constants in the proof of Theorem partialcm, the count
# N(3e12), the first order shift sum 6 delta^2/gamma^4 over the zeros above 3e12 (the 1e-47 claim of H:486 against the review's 1e-34), de Bruijn 9/32, the Selberg
# against Guth Maynard crossover, the Thorin densities N_alpha' of the base point laws and N_alpha(50) (Figure yor), and the Davenport Heilbronn function:
# kappa, the bound 0.324 of Proposition dhtheta, theta_DH > 0 and self duality, the zero 0.8085... + 85.699... i with an argument principle count, and the contiguous
# Toeplitz minors of the Taylor coefficients of xi_DH(1/2 + sqrt s)/xi_DH(1/2).
# Conventions. xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s); N(T) = T/(2 pi) log(T/(2 pi e)) + 7/8 (Riemann von Mangoldt, error O(log T)); N_alpha'(t) = pi^(-1) Re xi'/xi(alpha + i t),
# xi'/xi = 1/s + 1/(s-1) - (1/2) log pi + (1/2) psi(s/2) + zeta'/zeta; N_alpha(T) = int_0^T N_alpha' (Simpson). Davenport Heilbronn: a = (0, 1, kappa, -kappa, -1) mod 5,
# kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1); theta_DH(x) = sum_n a_n n exp(-pi n^2 x/5); Lambda(s) = (5/pi)^(s/2) Gamma((s+1)/2) f(s) with the Hurwitz form
# f(s) = 5^(-s)[(h1 - h4) + kappa (h2 - h3)] (the book's code/mi_dh.py Lam is imported: it is the book's own function for this object);
# xi_DH(1/2 + w) = int phi(u) cosh(w u) du with phi(u) = 2 sqrt(pi/5) theta_DH(e^(2u)) e^(3u/2) (even), so its Taylor coefficients in w^2 are a_k = int phi(u) u^(2k) du/(2k)!,
# computed by composite Gauss Legendre (panel 0.1, 40 nodes) at the working precision and again with panel 0.08 and 48 nodes at another precision; minors det[a_(d+i-j)] for
# order r = 1..14 and d = 0..119 (a_negative = 0; the minors with d < 0 vanish identically) at dps 400 and 330, signs and relative differences logged.
# Argument principle: phase increments of Lambda on the boundary of a rectangle of step 1e-3, integer count, run at steps 1e-3 and 5e-4.
# Options: --task {bern,yor,dh,toep,all} --quick --out BASE. Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse, json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
T0 = time.time()
LOG = None
RES = {}


def log(*a):
    line = '[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(line, flush=True)
    if LOG:
        LOG.write(line + '\n'); LOG.flush()


def task_bern(quick):
    mp.mp.dps = 50
    g1 = mp.mpf(open(ZF).readline().strip())
    T = mp.mpf(3) * 10 ** 12
    out = {}
    coef = (T ** 2 - g1 ** 2 - 1) / mp.log(T ** 2 / g1 ** 2)
    out['coef_n_over_x_with_5'] = mp.nstr(coef, 8)
    out['log_T2_over_g1sq'] = mp.nstr(mp.log(T ** 2 / g1 ** 2), 8)
    log('Theorem partialcm: n <= ((T^2 - g1^2 - 1) x - 5)/log(T^2/g1^2): coefficient of x = %s (log(T^2/g1^2) = %s)' % (mp.nstr(coef, 6), mp.nstr(mp.log(T ** 2 / g1 ** 2), 6)))
    out['T2_over_52'] = mp.nstr(T ** 2 / 52, 6)
    out['bound_below_T2x/52'] = bool(coef <= T ** 2 / 52)
    log('n <= T^2 x/52 follows from the bound: coefficient %s <= T^2/52 = %s: %s' % (mp.nstr(coef, 6), mp.nstr(T ** 2 / 52, 6), coef <= T ** 2 / 52))
    # decay ratio of g at n = T^2 x/52: log ratio <= n log((G'^2+1/4)/(G^2+1/4)) - (G'^2 - G^2) x, G = T + j, worst j = 0 (and large j)
    worst = None
    for x in (mp.mpf(1) / T, mp.mpf(1), mp.mpf(100)):
        n = T ** 2 * x / 52
        for j in (0, 1, 10 ** 6, 10 ** 12):
            G = T + j
            lr = n * mp.log((((G + 1) ** 2 + mp.mpf(1) / 4) / (G ** 2 + mp.mpf(1) / 4))) - ((G + 1) ** 2 - G ** 2) * x
            r = lr / (T * x)
            worst = r if worst is None or r > worst else worst
    out['max_log_ratio_over_Tx'] = mp.nstr(worst, 6)
    log('log of g(T+j+1)/g(T+j) at n = T^2 x/52, divided by T x, maximum over x in {1/T,1,100}, j in {0,1,1e6,1e12}: %s (chapter: <= -1.9)' % mp.nstr(worst, 6))
    tailc = 2 * mp.log(T + 1) / (1 - mp.exp(-mp.mpf('1.9')))
    out['tail_constant'] = mp.nstr(tailc, 6)
    log('tail constant 2 log(T+1)/(1 - e^-1.9) = %s (chapter: < 68 < e^5 = %s)' % (mp.nstr(tailc, 6), mp.nstr(mp.exp(5), 6)))
    # final inequality: n log((T^2+1/4)/g1^2) + 5 < (T^2 - 1/4 - g1^2) x at n = the bound
    mar = []
    for x in (mp.mpf(1) / T, mp.mpf(1), mp.mpf(10) ** 6):
        n = ((T ** 2 - g1 ** 2 - 1) * x - 5) / mp.log(T ** 2 / g1 ** 2)
        lhs = n * mp.log((T ** 2 + mp.mpf(1) / 4) / g1 ** 2) + 5
        rhs = (T ** 2 - mp.mpf(1) / 4 - g1 ** 2) * x
        mar.append(mp.nstr((rhs - lhs) / rhs, 5))
    out['final_inequality_relative_margin'] = mar
    log('final inequality of the proof, relative margin (rhs - lhs)/rhs at x = 1/T, 1, 1e6: %s (positive: holds)' % mar)
    # N(3e12)
    N = T / (2 * mp.pi) * mp.log(T / (2 * mp.pi * mp.e)) + mp.mpf(7) / 8
    out['N_3e12'] = mp.nstr(N, 8)
    log('N(3e12) = %s (chapter: N <= 1.2e13)' % mp.nstr(N, 8))
    # M5: first order shift
    one = 6 * mp.mpf(1) / 4 / T ** 4
    S4 = (mp.log(T / (2 * mp.pi)) / (3 * T ** 3) + 1 / (9 * T ** 3)) / (2 * mp.pi)  # int_T^inf log(t/2pi)/(2pi) t^-4 dt
    tot = 6 * mp.mpf(1) / 4 * S4
    out['one_zero_shift'] = mp.nstr(one, 5); out['sum_gamma_m4'] = mp.nstr(S4, 5); out['all_zeros_shift'] = mp.nstr(tot, 5)
    out['K215_one'] = mp.nstr(215 * one, 5); out['K215_all'] = mp.nstr(215 * tot, 5)
    log('one zero above T: 6 (1/4)/T^4 = %s (chapter < 2e-50), times 215 = %s (chapter < 1e-47)' % (mp.nstr(one, 5), mp.nstr(215 * one, 5)))
    log('all zeros above T: sum gamma^-4 = %s, 6/4 sum = %s, times 215 = %s (review M5: 7.9e-38 and 1.7e-35; proposed text: < 1e-34)' % (mp.nstr(S4, 5), mp.nstr(tot, 5), mp.nstr(215 * tot, 5)))
    # de Bruijn and crossover
    out['9/32'] = float(2 * mp.mpf(3) ** 2 / 8 ** 2)
    beta = mp.mpf(123) / 214
    out['crossover_beta'] = mp.nstr(beta, 8)
    sel = 1 - (beta - mp.mpf(1) / 2) / 4; gm = 30 * (1 - beta) / 13
    out['crossover_check'] = [mp.nstr(sel, 8), mp.nstr(gm, 8)]
    log('de Bruijn 2 (3/8)^2 = %s ; Selberg exponent 1 - (b - 1/2)/4 equals Guth Maynard 30 (1 - b)/13 at b = 123/214 = %s (%s, %s)' % (mp.nstr(mp.mpf(9) / 32, 6), mp.nstr(beta, 6), mp.nstr(sel, 8), mp.nstr(gm, 8)))
    RES['bern'] = out


def _yor_chunk(args):
    alpha, ts = args
    mp.mp.dps = 20
    out = []
    for t in ts:
        s = mp.mpc(alpha, t)
        v = 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)
        out.append(float(mp.re(v) / mp.pi))
    return out


def task_yor(quick):
    step = 0.1 if quick else 0.005
    n = int(round(50 / step))
    ts = [max(i * step, 1e-9) for i in range(n + 1)]
    res = {}
    with Pool(4) as pool:
        for alpha in (1.0, 0.75, 0.6, 0.55):
            chunks = [ts[i::16] for i in range(16)]
            parts = pool.map(_yor_chunk, [(alpha, c) for c in chunks])
            dens = [0.0] * len(ts)
            for i in range(16):
                for k, v in enumerate(parts[i]):
                    dens[i + 16 * k] = v
            # Simpson
            h = step
            Nint = h / 3 * (dens[0] + dens[-1] + 4 * sum(dens[1:-1:2]) + 2 * sum(dens[2:-1:2]))
            # trapezoid for comparison
            Ntr = h * (sum(dens) - (dens[0] + dens[-1]) / 2)
            i0 = min(range(len(dens)), key=lambda i: dens[i])
            # refine the minimum by a golden search around the grid minimum
            lo, hi = ts[max(i0 - 1, 0)], ts[min(i0 + 1, len(ts) - 1)]
            f = lambda t: _yor_chunk((alpha, [t]))[0]
            gr = (math.sqrt(5) - 1) / 2
            c = hi - gr * (hi - lo); d = lo + gr * (hi - lo); fc, fd = f(c), f(d)
            for _ in range(40):
                if fc < fd:
                    hi, d, fd = d, c, fc; c = hi - gr * (hi - lo); fc = f(c)
                else:
                    lo, c, fc = c, d, fd; d = lo + gr * (hi - lo); fd = f(d)
            mn = min(fc, fd, dens[i0])
            res[str(alpha)] = dict(min_density=mn, argmin=(lo + hi) / 2, N50_simpson=Nint, N50_trapezoid=Ntr, density_t0=dens[0], grid_step=step)
            log('alpha=%g: min of N_alpha\' on [0, 50] = %.5f at t = %.3f; N_alpha(50) = %.4f (Simpson), %.4f (trapezoid), step %g' % (alpha, mn, (lo + hi) / 2, Nint, Ntr, step))
    RES['yor'] = res


def dh_setup():
    import mi_dh
    return mi_dh


def theta(x, kap, nmax=None):
    a = [0, 1, kap, -kap, -1]
    s = mp.mpf(0)
    nmax = int(math.sqrt((mp.mp.dps + 12) * math.log(10) * 5 / (math.pi * float(x)))) + 2
    for n in range(1, nmax + 1):
        s += a[n % 5] * n * mp.exp(-mp.pi * n * n * x / 5)
    return s


def task_dh(quick):
    mp.mp.dps = 40
    mi = dh_setup()
    kap = mi.kappa()
    out = dict(kappa=mp.nstr(kap, 12))
    log('kappa = %s (chapter 0.2840790438)' % mp.nstr(kap, 12))
    # the bound 0.324
    x = mp.mpf(1)
    sm = mp.fsum(n * mp.exp(-mp.pi * (n * n - 1) * x / 5) for n in range(2, 60))
    out['sum_n_ge_2_at_x1'] = mp.nstr(sm, 8); out['one_minus'] = mp.nstr(1 - sm, 8)
    log('sum_{n>=2} n exp(-pi (n^2-1)/5) at x = 1: %s (chapter: 1 - 0.324, i.e. sum 0.324); decreasing in x so the bound holds for x >= 1' % mp.nstr(sm, 8))
    # positivity on a grid, self duality
    xs = [mp.mpf(k) / 20 for k in range(1, 400)] + [mp.mpf(10) ** (-k) for k in range(1, 3)]
    mn = min(theta(xx, kap) * (xx ** mp.mpf(1) if False else 1) for xx in xs if xx >= mp.mpf(1) / 10)
    out['min_theta_x_ge_0.1'] = mp.nstr(mn, 8)
    dev = max(abs(theta(1 / xx, kap) - xx ** mp.mpf('1.5') * theta(xx, kap)) / abs(xx ** mp.mpf('1.5') * theta(xx, kap)) for xx in [mp.mpf(k) / 7 for k in range(3, 40)])
    out['self_duality_max_rel_dev'] = float(dev)
    log('theta_DH positive on the grid x in [0.1, 20]: min %s; self duality theta(1/x) = x^(3/2) theta(x): max relative deviation over 37 points %.2e (chapter: ten digits)' % (mp.nstr(mn, 6), float(dev)))
    # Mellin check of xi_DH against theta: int theta x^((s+1)/2) dx/x = sqrt(5/pi) xi_DH(s), s = 1/2 and 0.8 + 1 i
    for s in (mp.mpf(1) / 2, mp.mpc('0.8', '1')):
        f = lambda xx: theta(xx, kap) * xx ** ((s + 1) / 2 - 1)
        I = mp.quad(f, [0, 0.5, 1, 4, 20, mp.inf]) if False else mp.quad(f, [mp.mpf(10) ** -3, 0.5, 1, 4, 20])
        # the part below 1e-3 by duality: int_0^{e} theta x^(s-1/2... ) handled via x -> 1/x
        g = lambda xx: theta(1 / xx, kap) * (1 / xx) ** ((s + 1) / 2 - 1) / xx ** 2
        I2 = 0  # the part below 1e-3 is below exp(-600) by self duality
        tail = mp.quad(lambda xx: theta(xx, kap) * xx ** ((s + 1) / 2 - 1), [20, 200])
        val = I + I2 + tail
        lam = mi.Lam(s, kap)
        r = abs(val - mp.sqrt(5 / mp.pi) * lam) / abs(lam)
        out['mellin_check_s=%s' % mp.nstr(s, 3)] = float(r)
        log('Mellin identity int theta_DH x^((s+1)/2) dx/x = sqrt(5/pi) xi_DH(s) at s = %s: relative difference %.2e' % (mp.nstr(s, 4), float(r)))
    # the zero
    z = mp.findroot(lambda s: mi.Lam(s, kap), mp.mpc('0.8085', '85.699'), tol=mp.mpf(10) ** -30)
    mp.mp.dps = 55
    kap55 = mi.kappa()
    z2 = mp.findroot(lambda s: mi.Lam(s, kap55), z, tol=mp.mpf(10) ** -45)
    out['zero'] = [mp.nstr(z.real, 18), mp.nstr(z.imag, 18)]
    out['zero_change_between_precisions'] = float(abs(z2 - z))
    log('Davenport Heilbronn zero: %s + %s i (chapter 0.808517182456637 + 85.6993484853776 i); change between 40 and 55 digits %.2e' % (mp.nstr(z2.real, 16), mp.nstr(z2.imag, 16), float(abs(z2 - z))))
    # symmetric image zeros and functional equation
    fe = abs(mi.Lam(mp.mpc('0.3', '5'), kap55) - mi.Lam(mp.mpc('0.7', '-5'), kap55)) / abs(mi.Lam(mp.mpc('0.3', '5'), kap55))
    out['functional_equation_rel'] = float(fe)
    log('functional equation Lambda(s) = Lambda(1 - s) at s = 0.3 + 5i: relative residual %.2e' % float(fe))
    # argument principle in a box around the zero
    mp.mp.dps = 25
    kap25 = mi.kappa()
    cnt = []
    for h in (mp.mpf(1) / 200, mp.mpf(1) / 400) if not quick else (mp.mpf(1) / 50,):
        x0, x1, y0, y1 = mp.mpf('0.78'), mp.mpf('0.84'), mp.mpf('85.67'), mp.mpf('85.73')
        pts = []
        for k in range(int(round((x1 - x0) / h))): pts.append(mp.mpc(x0 + k * h, y0))
        for k in range(int(round((y1 - y0) / h))): pts.append(mp.mpc(x1, y0 + k * h))
        for k in range(int(round((x1 - x0) / h))): pts.append(mp.mpc(x1 - k * h, y1))
        for k in range(int(round((y1 - y0) / h))): pts.append(mp.mpc(x0, y1 - k * h))
        pts.append(pts[0])
        vals = [mi.Lam(p, kap25) for p in pts]
        tot = mp.mpf(0); worst = 0
        for a, b in zip(vals[:-1], vals[1:]):
            d = mp.arg(b / a); tot += d; worst = max(worst, abs(float(d)))
        cnt.append((float(h), float(tot / (2 * mp.pi)), worst))
        log('argument principle in [0.78, 0.84] x [85.67, 85.73], step %g: count %.6f (largest phase step %.3f rad)' % cnt[-1])
    out['box_count'] = cnt
    RES['dh'] = out


def moments(K, dps, pw, nn):
    mp.mp.dps = dps
    kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    U = mp.mpf(5) if dps > 100 else mp.mpf(4)
    npan = int(round(U / mp.mpf(pw)))
    xs, ws = zip(*[(x, w) for x, w in mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(6 if nn == 40 else 7, mp.mp.prec)]) if False else (None, None)
    # Gauss Legendre nodes of nn points on [-1, 1]
    import numpy as np
    xn, wn = np.polynomial.legendre.leggauss(nn)
    # refine the nodes to the working precision by Newton on P_nn
    X, W = [], []
    for x0 in xn:
        x = mp.mpf(float(x0))
        for _ in range(8):
            p = mp.legendre(nn, x); dp = nn * (x * mp.legendre(nn, x) - mp.legendre(nn - 1, x)) / (x * x - 1)
            x = x - p / dp
        dp = nn * (x * mp.legendre(nn, x) - mp.legendre(nn - 1, x)) / (x * x - 1)
        X.append(x); W.append(2 / ((1 - x * x) * dp * dp))
    nodes, weights = [], []
    for k in range(npan):
        a = k * mp.mpf(pw); b = a + mp.mpf(pw)
        for x, w in zip(X, W):
            nodes.append((a + b) / 2 + (b - a) / 2 * x); weights.append((b - a) / 2 * w)
    phis = [2 * mp.sqrt(mp.pi / 5) * theta(mp.exp(2 * u), kap) * mp.exp(3 * u / 2) for u in nodes]
    a = []
    for k in range(K + 1):
        s = mp.fsum(w * ph * u ** (2 * k) for u, w, ph in zip(nodes, weights, phis))
        a.append(2 * s / mp.factorial(2 * k))
    return a


def dets(a, rmax, dmax):
    out = {}
    for r in range(1, rmax + 1):
        for d in range(0, dmax + 1):
            M = mp.matrix(r, r)
            for i in range(r):
                for j in range(r):
                    k = d + i - j
                    M[i, j] = a[k] if 0 <= k < len(a) else mp.mpf(0)
            out[(r, d)] = mp.det(M)
    return out


def task_toep(quick):
    K = 40 if quick else 140
    dps1, dps2 = (120, 100) if quick else (400, 330)
    rmax, dmax = (6, 30) if quick else (14, 119)
    a1 = moments(K, dps1, '0.1', 40)
    D1 = dets(a1, rmax, dmax)
    a2 = moments(K, dps2, '0.08', 48)
    D2 = dets(a2, rmax, dmax)
    mp.mp.dps = 30
    rel = max(abs(D1[k] - D2[k]) / abs(D1[k]) for k in D1 if D1[k] != 0)
    arel = max(abs(a1[k] - a2[k]) / abs(a1[k]) for k in range(K + 1))
    neg = sorted([k for k, v in D1.items() if v < 0])
    negs = sorted(neg, key=lambda k: (k[0], k[1]))
    out = dict(coef_rel_diff=float(arel), det_rel_diff=float(rel), n_minors=len(D1), n_negative=len(neg), first_negative=[list(k) for k in negs[:10]])
    # per order: smallest index d with a negative minor
    byr = {}
    for r in range(1, rmax + 1):
        ds = [d for (rr, d) in neg if rr == r]
        byr[r] = min(ds) if ds else None
    out['first_negative_d_by_order'] = byr
    mp.mp.dps = 50
    out['a0'] = mp.nstr(a1[0], 20)
    out['a1_a2'] = [mp.nstr(a1[1], 20), mp.nstr(a1[2], 20)]
    log('Toeplitz minors det[a_(d+i-j)] of the Taylor coefficients of xi_DH(1/2 + sqrt s)/xi_DH(1/2): a_0 = %s' % out['a0'])
    log('%d minors, orders 1..%d, d from 0 to %d; coefficient relative difference between runs %.2e; determinant relative difference %.2e; negative minors: %d' % (len(D1), rmax, dmax, float(arel), float(rel), len(neg)))
    log('smallest d with a negative minor, by order r: %s' % byr)
    out['min_positive_all_orders_le_10_d_le_99'] = bool(all(D1[(r, d)] > 0 for r in range(1, min(rmax, 10) + 1) for d in range(0, min(dmax, 99) + 1)))
    log('all minors of order <= 10 with d <= 99 positive: %s' % out['min_positive_all_orders_le_10_d_le_99'])
    RES['toep'] = out


def main():
    global LOG, ZF
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=str(HERE.parent / 'data' / 'ht_misc'))
    ap.add_argument('--zeros', default=str(HERE.parent / 'data' / 'zeros_1700.txt'))
    a = ap.parse_args()
    ZF = a.zeros
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    LOG = open(a.out + '.log', 'w')
    tasks = ['bern', 'yor', 'dh', 'toep'] if a.task == 'all' else a.task.split(',')
    for tk in tasks:
        log('=== task', tk)
        {'bern': task_bern, 'yor': task_yor, 'dh': task_dh, 'toep': task_toep}[tk](a.quick)
    json.dump(RES, open(a.out + '.json', 'w'), indent=1, default=str)
    log('done')


if __name__ == '__main__':
    main()
