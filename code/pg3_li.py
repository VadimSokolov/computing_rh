# pg3_li.py: Li's coefficients from the heat trace and the tilted clock (Proposition prop:ch20:heat and the paragraph after it, Chapter 20,
# merged from note X3; placeholders N3a, N3b, N3c). Added for the book in October 2026, written from the mathematics (the collaborator's
# scripts were not read, run or copied).
# Conventions. xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2; ell_n(t) = sum_{k<n} (-1)^k n C(n+k,2k+1) t^k/(k+1)! (the coefficient of
# y^(n-1) in q'(y) exp(-t q(y)), q = y/(1-y)^2); phi_n(z) = 2 - r^n - r^(-n), r = (2 sqrt z + i)/(2 sqrt z - i); T~ the clock tilted at u = 1/4.
# Count: every part uses the same 100 cumulants kappa_j, j <= 100 (the quick run 12), stored in the json; N3a uses kappa_2.
# (a) kappa_j = (-1)^(j+1) j! [v^j] F(v), F(v) = log xi(1/2 + (1/4+v)^(1/2)) (even in the root, analytic and zero free for |v| < gamma_1^2 + 1/4
#     = 199.79+1/4), by the trapezoid rule on |v| = rad with NP points (the aliasing error is (rad/200)^NP times the size of F, below 1e-300 for
#     rad = 100, NP = 1024), the branch of log by unwrapping the argument along the circle. Run at two radii and precisions; the maximum relative
#     difference of the kappa_j, j <= 100, is the error of the table. Check: kappa_1 = lambda_1 and Var T~ = kappa_2 = 4 lambda_1 - lambda_2 (N3a).
# (b) lambda_n = sum_j (-1)^(j+1) n/j! C(n+j-1,2j-1) kappa_j for n in {1,2,3,5,10,20,50,100} against tab:ch20:li (20 significant digits printed):
#     every row must agree within half a unit of the last printed digit; N3b is the number of agreeing digits (the printed 20).
# (c) lambda_n = int_0^inf W(t) e^(-t/4) ell_n(t) dt for n <= 4, with W the arithmetic formula for t <= t0 = 0.004 (polar + archimedean +
#     prime terms of eq:ch2:Warith, in the form of ch12.tex:22), integrated exactly in t against ell_n (incomplete gamma functions) and then
#     in u (mp.quad), and W = sum_{gamma} e^(-gamma^2 t) over the first 8000 zeros of data/zeros_8000.txt for t >= t0, integrated exactly.
#     (Zeros with gamma^2 t0 > 160 contribute less than e^(-160); the zeros above the 8000th contribute nothing in either part beyond that.)
#     W(t0) from the arithmetic formula against the zero sum is logged. N3c is the number of agreeing digits against (b).
# (d) sympy: ell_n for n <= 8; the identity of (i) as rational functions for n <= 8; the composition identity of (ii) for n <= 8; the exact
#     value (4 - 2824 e^-5)/25 of int_0^20 e^(-t/4) ell_3(t) (1-t/20)^2 dt.
# Options: --quick (dps 50, NP = 128, kappa_j for j <= 12, one radius, part (c) for n <= 2 at 20 digits; seconds on a laptop), --out FILE.
import sys, json, argparse, os, time
import mpmath as mp

TABLE = {1: '0.023095708966121033814', 2: '0.092345735228046670386', 3: '0.20763892055432480379', 5: '0.57554271446117745243',
         10: '2.2793393631931577437', 20: '8.7692768720932151296', 50: '43.531096488374019299', 100: '118.60377537679132993'}

def zfile(n):
    for d in ('.', '../data', 'data'):
        f = os.path.join(d, 'zeros_8000.txt')
        if os.path.exists(f):
            return [l.strip() for l in open(f).read().split('\n') if l.strip()][:n]
    raise SystemExit('zeros_8000.txt not found')

def logxi(s):
    return mp.log(s * (s - 1) / 2) - s / 2 * mp.log(mp.pi) + mp.loggamma(s / 2) + mp.log(mp.zeta(s))

def kappas(dps, rad, NP, jmax):
    mp.mp.dps = dps
    rad = mp.mpf(rad)
    half = NP // 2
    vals = []
    for k in range(half + 1):
        v = rad * mp.expjpi(mp.mpf(2 * k) / NP)
        s = mp.mpf(1) / 2 + mp.sqrt(mp.mpf(1) / 4 + v)
        # continuous log: log xi = log|xi| + i arg with the argument unwrapped (log of the ratio to the previous value)
        x = s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
        vals.append(x)
    lg = [mp.log(vals[0])]
    for k in range(1, half + 1):
        lg.append(lg[-1] + mp.log(vals[k] / vals[k - 1]))      # increments are small (|arg| < pi), so the principal log is continuous
    F = lg + [mp.conj(lg[NP - k]) for k in range(half + 1, NP)]
    c = []
    for j in range(jmax + 1):
        c.append(mp.fsum(F[k] * mp.expjpi(-mp.mpf(2 * j * k) / NP) for k in range(NP)) / NP / rad ** j)
    kap = [None] + [(-1) ** (j + 1) * mp.factorial(j) * mp.re(c[j]) for j in range(1, jmax + 1)]
    imagmax = max(abs(mp.im(c[j])) * mp.factorial(j) for j in range(1, jmax + 1))
    return kap, imagmax

def lam_from_kappa(n, kap):
    return mp.fsum((-1) ** (j + 1) * n / mp.factorial(j) * mp.binomial(n + j - 1, 2 * j - 1) * kap[j] for j in range(1, n + 1))

def ell_coeffs(n):
    return [(-1) ** k * n * mp.binomial(n + k, 2 * k + 1) / mp.factorial(k + 1) for k in range(n)]

def Lambda_vm(m):
    # von Mangoldt
    for p in range(2, m + 1):
        if m % p == 0:
            q = m
            while q % p == 0:
                q //= p
            return mp.log(p) if q == 1 else mp.mpf(0)
    return mp.mpf(0)

def part_c(ns, dps, t0, log):
    mp.mp.dps = dps
    t0 = mp.mpf(t0)
    zs = [mp.mpf(z) for z in zfile(8000)]
    res = {}
    # W check at t0
    def theta_p(u):
        return (mp.re(mp.digamma(mp.mpf(1) / 4 + 1j * u / 2)) - mp.log(mp.pi)) / 2
    def W_arith(t):
        arch = mp.quad(lambda u: mp.exp(-u * u * t) * theta_p(u), [0, 5, 20, 60, mp.sqrt(2.5 * dps * mp.log(10) / t) + 5]) / mp.pi
        pr = mp.fsum(Lambda_vm(m) / mp.sqrt(m) * mp.exp(-mp.log(m) ** 2 / (4 * t)) for m in range(2, 200)) / mp.sqrt(4 * mp.pi * t)
        return mp.exp(t / 4) + arch - pr
    wa = W_arith(t0); wz = mp.fsum(mp.exp(-g * g * t0) for g in zs)
    log('W(t0 = %s): arithmetic %s, zero sum %s, difference %s' % (mp.nstr(t0, 3), mp.nstr(wa, 25), mp.nstr(wz, 25), mp.nstr(wa - wz, 3)))
    res['W_t0_diff'] = mp.nstr(wa - wz, 3)
    zsel = [g for g in zs if g * g * t0 < 160]
    for n in ns:
        a = ell_coeffs(n)
        # zeros part: sum_gamma int_t0^inf e^(-c t) ell_n dt, c = gamma^2 + 1/4
        zp = mp.mpf(0)
        for g in zsel:
            c = g * g + mp.mpf(1) / 4
            zp += mp.fsum(a[k] * mp.gammainc(k + 1, c * t0, mp.inf) / c ** (k + 1) for k in range(n))
        polar = mp.fsum(a[k] * t0 ** (k + 1) / (k + 1) for k in range(n))
        def inner(u):
            c = u * u + mp.mpf(1) / 4
            return mp.fsum(a[k] * mp.gammainc(k + 1, 0, c * t0) / c ** (k + 1) for k in range(n))
        arch = mp.quad(lambda u: theta_p(u) * inner(u), [0, 1, 5, 20, 100, 400, 3000, mp.inf]) / mp.pi
        prim = mp.mpf(0)
        for m in range(2, 60):
            L = Lambda_vm(m)
            if L != 0:
                b2 = mp.log(m) ** 2 / 4
                prim += L / mp.sqrt(m) * mp.quad(lambda t: t ** (-mp.mpf(1) / 2) * mp.exp(-b2 / t - t / 4) * mp.polyval(a[::-1], t), [0, t0 / 4, t0])
        prim /= mp.sqrt(4 * mp.pi)
        lam = zp + polar + arch - prim
        res[n] = dict(lam=mp.nstr(lam, 30), zeros=mp.nstr(zp, 20), polar=mp.nstr(polar, 20), arch=mp.nstr(arch, 20), prime=mp.nstr(prim, 8))
        log('n=%d (c): lambda = %s  [zeros %s, polar %s, archimedean %s, primes %s]' % (n, mp.nstr(lam, 25), mp.nstr(zp, 12), mp.nstr(polar, 12), mp.nstr(arch, 12), mp.nstr(prim, 5)))
    return res

def part_d(log):
    import sympy as sp
    t, y, z = sp.symbols('t y z', positive=True)
    q = y / (1 - y) ** 2
    ok = True
    for n in range(1, 9):
        ell = sum((-1) ** k * n * sp.binomial(n + k, 2 * k + 1) * t ** k / sp.factorial(k + 1) for k in range(n))
        gen = sp.series(sp.diff(q, y) * sp.exp(-t * q), y, 0, n).removeO()
        ok &= sp.simplify(sp.expand(gen).coeff(y, n - 1) - ell) == 0
        # (i): int_0^inf e^{-(z+1/4)t} ell dt = phi_n(z)
        phi = sum((-1) ** k * n * sp.binomial(n + k, 2 * k + 1) * sp.factorial(k) / (z + sp.Rational(1, 4)) ** (k + 1) / sp.factorial(k + 1) for k in range(n))
        s = sp.sqrt(z)
        r = (2 * s + sp.I) / (2 * s - sp.I)
        ok &= sp.simplify(phi - (2 - r ** n - r ** (-n))) == 0
        # (ii) composition: n [y^n] sum_j (-1)^(j+1) kappa_j q^j/j! has coefficient n/j! C(n+j-1,2j-1) of (-1)^(j+1) kappa_j
        kap = sp.symbols('k1:%d' % (n + 1))
        ser = sp.series(sum((-1) ** (j + 1) * kap[j - 1] * q ** j / sp.factorial(j) for j in range(1, n + 1)), y, 0, n + 1).removeO()
        rhs = sum((-1) ** (j + 1) * sp.Rational(n, 1) / sp.factorial(j) * sp.binomial(n + j - 1, 2 * j - 1) * kap[j - 1] for j in range(1, n + 1))
        ok &= sp.expand(n * sp.expand(ser).coeff(y, n) - rhs) == 0
    log('(d) ell_n generating function, identity (i) as rational functions and the composition identity (ii), n = 1..8: %s' % ok)
    val = sp.integrate(sp.exp(-t / 4) * sum((-1) ** k * 3 * sp.binomial(3 + k, 2 * k + 1) * t ** k / sp.factorial(k + 1) for k in range(3)) * (1 - t / 20) ** 2, (t, 0, 20))
    ex = (4 - 2824 * sp.exp(-5)) / 25
    log('(d) R6: integral = %s, (4 - 2824 e^-5)/25 = %s, equal: %s, numerical value %s' % (sp.simplify(val), sp.simplify(ex), sp.simplify(val - ex) == 0, sp.N(ex, 12)))
    return dict(identities=bool(ok), r6_exact=bool(sp.simplify(val - ex) == 0), r6_value=str(sp.N(ex, 15)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='pg3_li.json')
    a = ap.parse_args()
    L = []
    def log(s):
        print(s); sys.stdout.flush(); L.append(s)
    res = {}
    if a.quick:
        runs = [(50, 100, 128, 12)]
        ns = [1, 2, 3, 5, 10]
    else:
        runs = [(180, 100, 1024, 100), (150, 80, 1024, 100)]
        ns = [1, 2, 3, 5, 10, 20, 50, 100]
    K = []
    for (dps, rad, NP, jm) in runs:
        t0 = time.time()
        kap, im = kappas(dps, rad, NP, jm)
        K.append(kap)
        log('kappa_j, j <= %d: dps %d, radius %d, %d points, %.0f s, largest imaginary part of j! c_j %s' % (jm, dps, rad, NP, time.time() - t0, mp.nstr(im, 3)))
    mp.mp.dps = 50
    rel = max(abs(K[0][j] - K[-1][j]) / abs(K[0][j]) for j in range(1, len(K[0]))) if len(K) > 1 else None
    if rel is not None:
        log('largest relative difference of kappa_j between the two runs: %s' % mp.nstr(rel, 3))
        res['kappa_rel_diff'] = mp.nstr(rel, 3)
    kap = K[0]
    mp.mp.dps = 50
    lam1 = lam_from_kappa(1, kap)
    var = 4 * lam1 - lam_from_kappa(2, kap)
    log('(a) kappa_1 = %s (lambda_1), kappa_2 = Var T~ = 4 lambda_1 - lambda_2 = %s ; kappa_2 directly %s' % (mp.nstr(kap[1], 20), mp.nstr(var, 15), mp.nstr(kap[2], 15)))
    res['kappa'] = {j: mp.nstr(kap[j], 30) for j in range(1, len(kap))}
    res['N3a'] = mp.nstr(var, 9)
    lams = {}
    allok = True
    mind = None
    for n in ns:
        if n > len(kap) - 1:
            continue
        l = lam_from_kappa(n, kap)
        lams[n] = l
        t = mp.mpf(TABLE[n]); dec = len(TABLE[n].split('.')[1])
        d = abs(l - t)
        ok = d <= mp.mpf(10) ** (-dec) / 2 * (1 + mp.mpf(10) ** -10)
        allok &= ok
        rd = d / t
        mind = rd if mind is None or rd > mind else mind
        log('(b) n=%d: lambda_n = %s, table %s, difference %s, within half a unit of the last printed digit: %s' % (n, mp.nstr(l, 24), TABLE[n], mp.nstr(d, 3), ok))
    log('(b) all rows agree to the printed digits: %s (largest relative difference %s)' % (allok, mp.nstr(mind, 3)))
    res['N3b_all_rows_agree_printed_digits'] = bool(allok)
    res['lam_from_kappa'] = {n: mp.nstr(l, 30) for n, l in lams.items()}
    cns = [1, 2] if a.quick else [1, 2, 3, 4]
    pc = part_c(cns, 20 if a.quick else 45, '0.004', log)
    mp.mp.dps = 50
    cmp_ = {}
    for n in cns:
        l = lam_from_kappa(n, kap)
        c = mp.mpf(pc[n]['lam'])
        cmp_[n] = mp.nstr(abs(c - l) / l, 3)
        log('(c) n=%d: lambda from W %s, from kappa %s, relative difference %s' % (n, mp.nstr(c, 22), mp.nstr(l, 22), cmp_[n]))
    res['c'] = pc; res['c_rel_diff'] = cmp_
    res['d'] = part_d(log)
    json.dump(res, open(a.out, 'w'), indent=1, default=str)
    open(a.out.replace('.json', '.log'), 'w').write('\n'.join(L) + '\n')
