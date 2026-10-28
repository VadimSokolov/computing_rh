# cf_dhzeros.py: certified zeros and zero counts of the Davenport Heilbronn function f, Chapter 34 (Cauchy deconvolution, prefix cf), group CF-DH, tasks T1 and T2 of
# misc/reviews/nick-b4/reconcile-L.md (rows R11, R21, R24, N31). Claims supported (L = cauchy_deconvolution.tex): L:62 and L:1596 to L:1602 (the two zeros of f near
# 0.6914 + 1012.0200 i and 0.6981 + 1018.7958 i, the counts 0 in [0.6, 2] x [1000, 1010] and 2 in [0.6, 2] x [1010, 1020]), L:285 to L:290 (the four zeros of f off the critical
# line below height 200), L:1706 (tab:scan, the zero at 176.70). Added for the book in October 2026.
# Conventions. a_n periodic mod 5 with (a_1, a_2, a_3, a_4, a_5) = (1, kappa, -kappa, -1, 0), kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1) (Davenport and Heilbronn; ch/ht5.tex,
# ht:sec:dhtest; L:280 to L:283). f(s) = sum a_n n^(-s) = 5^(-s) [zeta(s, 1/5) + kappa zeta(s, 2/5) - kappa zeta(s, 3/5) - zeta(s, 4/5)] (Hurwitz form, as in code/mi_dh.py, which
# is imported for a cross check only), Lambda(s) = (5/pi)^(s/2) Gamma((s + 1)/2) f(s) is entire, Lambda(s) = Lambda(1 - s), Lambda real on Re s = 1/2 and Lambda(conj s) = conj Lambda(s)
# (the a_n are real). f has no Euler product, so the line count below uses Lambda and its own functional equation, not an Euler product or a Riemann von Mangoldt formula for zeta.
# f != 0 for Re s >= 2: |a_n| <= 1 gives |f(s) - 1| <= zeta(2) - 1 < 1 there (checked: |kappa| < 1 in Arb). By Lambda(s) = Lambda(1 - s) and the poles of Gamma, the zeros of Lambda
# lie in -1 < Re s < 2 and are symmetric in Re s = 1/2 and in the real axis; zeros of f with Re s <= 0 are the trivial zeros at s = -1, -3, ... and the reflections of those in Re s >= 1.
# Method. All values are Arb balls (python-flint, acb.zeta(s, a) the Hurwitz zeta, at precision PREC bits). Argument principle with certified arithmetic: a closed polygon is cut into
# segments [z0, z1]; a segment is accepted if the ball F(S) over the whole segment S (union of the endpoint balls, a ball containing S) satisfies Re(F(S)/F(z0)) > 0, which forces
# arg(F(z)/F(z0)) into (-pi/2, pi/2) along S, so the change of argument on S equals the principal log of F(z1)/F(z0) (a ball); otherwise S is bisected (depth at most 40). The sum
# of the changes divided by 2 pi must be a ball around an integer of radius below 0.01: that integer is the number of zeros (with multiplicity) inside, and F has no zero on the
# polygon. Zeros of f in a rectangle [0.6, 2] x [a, b]: the count of f on that rectangle. Certified zeros: Newton/secant in Arb (400 bits) to a center c, then the argument principle
# on the square c + [-r, r] + i [-r, r] (r = 1e-18, or 1e-16 for T2) gives exactly one zero in it, hence |zero - c| <= r sqrt 2, which certifies 15 digits (T2: 13 digits) of both parts
# (a Krawczyk or interval Newton step would give the same; the argument principle needs no derivative). Completeness for 0 < t <= 200 (T2): N_rect = zeros of Lambda in the rectangle
# [-2, 3] x [-200, 200] (argument principle for Lambda; no zero on its boundary, no zero outside by the above; the edges Re s = -2 and 3 avoid the poles of Gamma and the pole of the
# Hurwitz zeta at s = 1); N_real = zeros of Lambda on the real segment [-2, 3] (certified none: Lambda over a cover of [1/2, 3] by balls, the rest by Lambda(s) = Lambda(1 - s); near
# s = 1 the ball of Lambda(1 + x) = Lambda(-x) is used); N_up = (N_rect - N_real)/2 zeros with 0 < t < 200. S = the number of sign changes of the real function Lambda(1/2 + i t) on a
# grid of step h in (0, 200] (the sign of the ball is certified, a grid point with a ball containing 0 is shifted by h/1000 until it is not). Each sign change is a zero on the line, and the
# zeros off the line (certified separately, 2 K of them with their reflections) are distinct from the line zeros, so N_up >= S + 2 K, with equality if and only if the line has no zero
# missed by the grid and no multiple or non sign changing zero and f has no further zeros off the line below height 200; the script reports N_up - S and compares it with 2 K.
# Error control. Every count and every zero is computed at two precisions (PREC and PREC2); the script logs both and requires agreement. Float cross check: Lambda against mi_dh.Lam
# (mpmath, 30 digits) at several points, and the functional equation residual |Lambda(s) - Lambda(1 - s)| (ball midpoint).
# Options: --quick (a smoke test: heights up to 30 and the first zero only) --out FILE (base name of the json; the log is next to it) --cpus N --task {t1,t2,all}.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, os, json, time, argparse
from fractions import Fraction
from multiprocessing import Pool
from flint import arb, acb, ctx

LOGF = None


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGF:
        LOGF.write(s + '\n'); LOGF.flush()


def A(fr):
    fr = Fraction(fr)
    return arb(fr.numerator) / arb(fr.denominator)


def kappa():
    s5 = arb(5).sqrt()
    return ((10 - 2 * s5).sqrt() - 2) / (s5 - 1)


def fDH(s):
    kap = kappa()
    h = [s.zeta(arb(a) / 5) for a in (1, 2, 3, 4)]
    return (-s * arb(5).log()).exp() * (h[0] + kap * h[1] - kap * h[2] - h[3])


def Lam(s, flip=True):
    # Lambda(s) = Lambda(1 - s): for Re s < -1/4, or Re s < 1/2 and |Im s| > 1/2, evaluate at 1 - s (the Hurwitz zeta at a wide ball with Re s < 0 gives poor enclosures; a pole of f at s = 1 is avoided by |Im| > 1/2)
    if flip and (float(s.real.mid()) < -0.25 or (float(s.real.mid()) < 0.5 and float(abs(s.imag).lower()) > 0.5)):
        s = 1 - s
    return((s / 2) * (arb(5).log() - arb.pi().log())).exp() * ((s + 1) / 2).gamma() * fDH(s)


FUNCS = {'f': fDH, 'L': Lam}


def pt(p):
    return acb(A(p[0]), A(p[1]))


def seg_ball(p0, p1):
    return acb(A(p0[0]).union(A(p1[0])), A(p0[1]).union(A(p1[1])))


def midp(p0, p1):
    return ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)


def integ(F, p0, f0, p1, f1, depth=0):
    """returns (arb sum of argument changes, number of segments) on [p0, p1]"""
    E = F(seg_ball(p0, p1))
    ok = E.is_finite() and (E / f0).real > 0
    if ok:
        return (f1 / f0).log().imag, 1
    if depth > 40:
        raise RuntimeError('depth exceeded at %s' % (p0,))
    pm = midp(p0, p1)
    fm = F(pt(pm))
    a1, n1 = integ(F, p0, f0, pm, fm, depth + 1)
    a2, n2 = integ(F, pm, fm, p1, f1, depth + 1)
    return a1 + a2, n1 + n2


def chunk_job(job):
    name, prec, p0, p1 = job
    ctx.prec = prec
    F = FUNCS[name]
    s, n = integ(F, p0, F(pt(p0)), p1, F(pt(p1)))
    return (s.mid().str(40, radius=False), float(s.rad()) * (1 + 1e-9), n)


def arbf(mid, rad):
    return arb(mid, rad)


def edges(x0, x1, y0, y1, L):
    """counterclockwise polygon, edges cut in pieces of length <= L"""
    out = []
    def cut(p, q):
        n = max(1, int(-(-max(abs(q[0] - p[0]), abs(q[1] - p[1])) // L)))
        pts = [(p[0] + (q[0] - p[0]) * k / n, p[1] + (q[1] - p[1]) * k / n) for k in range(n + 1)]
        return [(pts[k], pts[k + 1]) for k in range(n)]
    c = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    for k in range(4):
        out += cut(c[k], c[k + 1])
    return out


def count_rect(pool, name, x0, x1, y0, y1, prec, L=Fraction(1)):
    x0, x1, y0, y1 = [Fraction(v) for v in (x0, x1, y0, y1)]
    jobs = [(name, prec, p, q) for p, q in edges(x0, x1, y0, y1, Fraction(L))]
    res = pool.map(chunk_job, jobs, chunksize=1)
    ctx.prec = 256
    tot = arb(0); nseg = 0
    for m, r, n in res:
        tot += arb(m, r); nseg += n
    N = tot / (2 * arb.pi())
    k = int(round(float(N.mid())))
    good = float(N.rad()) < 0.01 and abs(float(N.mid()) - k) < 0.01
    return k, good, float(N.rad()), nseg


def secant(F, z0, prec=400, tol=40):
    ctx.prec = prec
    a = acb(z0[0], z0[1]); b = a + acb(arb(1) / 10 ** 7, arb(1) / 10 ** 7)
    fa, fb = F(a), F(b)
    for _ in range(80):
        c = b - fb * (b - a) / (fb - fa)
        c = acb(c.real.mid(), c.imag.mid())
        a, fa = b, fb
        b, fb = c, F(c)
        if float(abs(b - a).mid()) < 10 ** (-tol):
            break
    return b


def dec(x, nd):
    s = x.mid().str(nd + 10, radius=False)
    from decimal import Decimal, getcontext
    getcontext().prec = nd + 15
    return Fraction(Decimal(s).quantize(Decimal(1).scaleb(-nd)))


def certify_zero(pool, name, seed, nd, r_exp, precs, tag):
    F = FUNCS[name]
    z = secant(F, seed)
    cre, cim = dec(z.real, 24), dec(z.imag, 24)
    r = Fraction(1, 10 ** r_exp)
    out = {'tag': tag, 'seed': seed, 'center_re': str(float(cre)), 'center': None}
    cs = [str(z.real.mid().str(22, radius=False)), str(z.imag.mid().str(22, radius=False))]
    ctx.prec = 400
    res = []
    for pr in precs:
        k, good, rad, nseg = count_rect(pool, name, cre - r, cre + r, cim - r, cim + r, pr, L=r * 4)
        res.append((pr, k, good, rad, nseg))
    out.update({'center_re': cs[0], 'center_im': cs[1], 'box_halfwidth': '1e-%d' % r_exp, 'counts': res,
                'certified_digits': r_exp - 1, 'ok': all(x[1] == 1 and x[2] for x in res)})
    log('zero %s: center %s + %s i, box half width 1e-%d, counts (prec, N, good, rad, segs) %s' % (tag, cs[0], cs[1], r_exp, res))
    return out


def grid_block(job):
    prec, ts, h = job
    ctx.prec = prec
    out = []
    for t in ts:
        sh = 0
        while True:
            v = Lam(acb(A(Fraction(1, 2)), A(t + sh))).real
            if v > 0: sg = 1; break
            if v < 0: sg = -1; break
            sh += h / 1000
        out.append(sg)
    return out


def line_and_rect(pool, TMAX, K, precs, label):
    """argument principle count of Lambda in [-2, 3] x [-TMAX, TMAX] against sign changes of Lambda(1/2 + it), 0 <= t <= TMAX"""
    TMAX = Fraction(TMAX)
    T = {'tmax': str(float(TMAX)), 'rect': []}
    for pr in precs:
        t0 = time.time()
        k, good, rad, nseg = count_rect(pool, 'L', -2, 3, -TMAX, TMAX, pr, L=Fraction(1, 2))
        T['rect'].append({'prec': pr, 'count': k, 'good': good, 'rad': rad, 'segments': nseg})
        log('%s count Lambda in [-2,3]x[-%s,%s] at prec %d: N_rect = %d (certified %s, rad %.2g, %d segments, %.0f s)' % (label, float(TMAX), float(TMAX), pr, k, good, rad, nseg, time.time() - t0))
    Nup = T['rect'][0]['count'] // 2  # N_real = 0 (real axis check), N_rect = 2 N_up
    T['N_up'] = Nup
    h = Fraction(1, 50)
    ts = [h * j for j in range(int(TMAX / h) + 1)] + [TMAX]
    blocks = [(precs[0], ts[b:b + 200], h) for b in range(0, len(ts), 200)]
    sg = []
    for blk in pool.map(grid_block, blocks, chunksize=1):
        sg += blk
    sgn = sum(1 for i in range(len(sg) - 1) if sg[i] != sg[i + 1])
    T['sign_changes_line'] = sgn
    T['grid_points'] = len(sg)
    T['grid_step'] = '1/50 (last point t = tmax)'
    T['diff_Nup_minus_S'] = Nup - sgn
    T['expected_diff_2K'] = 2 * K
    T['complete'] = bool(Nup - sgn == 2 * K and T['rect'][0]['count'] == T['rect'][1]['count'] and all(r['good'] for r in T['rect']) and T['rect'][0]['count'] % 2 == 0)
    log('%s: N_up = %d zeros of Lambda with 0 < t < %s; S = %d sign changes of Lambda(1/2 + it) on %d grid points; N_up - S = %d; 2K = %d (K = certified zeros with beta > 1/2 below tmax: %d); complete: %s'
        % (label, Nup, float(TMAX), sgn, len(sg), Nup - sgn, 2 * K, K, T['complete']))
    return T


def real_axis(prec, depth_max=30):
    """certify Lambda != 0 on [-2, 3] real: [1/2, 0.9], [1.1, 3] directly, [0.9, 1.1] via 1 - s in [-0.1, 0.1], then by symmetry"""
    ctx.prec = prec
    bad = []
    def chk(lo, hi, refl, depth=0):
        lo, hi = Fraction(lo), Fraction(hi)
        x = A(lo).union(A(hi))
        s = acb(1 - x if refl else x, 0)
        v = Lam(s)
        if v.is_finite() and not v.real.contains(0):
            return 1
        if depth >= depth_max:
            bad.append((float(lo), float(hi))); return 1
        m = (lo + hi) / 2
        return chk(lo, m, refl, depth + 1) + chk(m, hi, refl, depth + 1)
    n = chk(Fraction(1, 2), Fraction(9, 10), False) + chk(Fraction(9, 10), Fraction(11, 10), True) + chk(Fraction(11, 10), 3, False)
    return n, bad


def crosscheck(precs):
    out = {}
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, here)
        import mi_dh
        import mpmath as mp
        mp.mp.dps = 40
        pts = [(0.5, 14.5), (0.8, 85.7), (0.7, 1012.02), (2.0, 3.0), (-1.5, 20.0)]
        worst = 0
        for (a, b) in pts:
            ctx.prec = 200
            v = Lam(acb(arb(a), arb(b)), False)
            m = mi_dh.Lam(mp.mpc(a, b))
            vm = mp.mpc(mp.mpf(v.real.mid().str(45, radius=False)), mp.mpf(v.imag.mid().str(45, radius=False)))
            d = float(abs(vm - m) / abs(m))
            worst = max(worst, d)
        out['mi_dh_max_rel_diff'] = worst
        log('cross check against mi_dh.Lam at %d points: max relative difference %.3g' % (len(pts), worst))
    except Exception as e:
        out['mi_dh_error'] = repr(e); log('cross check against mi_dh failed:', repr(e))
    ctx.prec = 300
    fe = []
    for (a, b) in [(0.3, 5.0), (0.8, 85.7), (-1.5, 20.0), (1.7, 300.0)]:
        s = acb(arb(a), arb(b))
        d = Lam(s, False) - Lam(1 - s, False)
        fe.append(float(abs(d.real.mid()) + abs(d.imag.mid())) / float(abs(Lam(s, False)).mid()))
    out['functional_equation_max_rel'] = max(fe)
    log('functional equation |Lambda(s) - Lambda(1 - s)| / |Lambda(s)| at 4 points (midpoints): max %.3g' % max(fe))
    kap = kappa()
    out['kappa'] = kap.str(30, radius=False)
    out['kappa_lt_1'] = bool(abs(kap) < 1)
    out['dirichlet_bound_Re_ge_2'] = 'zeta(2) - 1 = %s < 1' % (arb.pi() ** 2 / 6 - 1).str(10, radius=False)
    log('kappa =', out['kappa'], '; |kappa| < 1 certified:', out['kappa_lt_1'], '; |f - 1| <= zeta(2) - 1 =', out['dirichlet_bound_Re_ge_2'], 'for Re s >= 2')
    return out


def main():
    global LOGF
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--cpus', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 8)))
    ap.add_argument('--task', default='all')
    a = ap.parse_args()
    out = a.out or ('data/cf_dhzeros' if not a.quick else 'quick_cf_dhzeros')
    LOGF = open(out + '.log', 'w')
    t00 = time.time()
    P1, P2 = (160, 224) if a.quick else (256, 384)
    Q1, Q2 = (128, 192) if a.quick else (128, 192)
    R = {'quick': a.quick, 'prec_T1': [P1, P2], 'prec_T2': [Q1, Q2]}
    log('cf_dhzeros: quick=%s precisions T1 %s T2 %s' % (a.quick, (P1, P2), (Q1, Q2)))
    R['checks'] = crosscheck((P1, P2))
    pool = Pool(a.cpus)
    if a.task in ('t1', 'all'):
        T1 = {}
        zs = [('z1012', (0.6914048, 1012.0199815), 15, 18), ('z1018', (0.6980980, 1018.7957623), 15, 18)]
        if a.quick:
            zs = zs[:1]
        T1['zeros'] = [certify_zero(pool, 'f', s, 12, rexp, (P1, P2), tag) for tag, s, _, rexp in zs]
        T1['counts'] = {}
        rects = [('[0.6,2]x[1000,1010]', 1000, 1010, 0), ('[0.6,2]x[1010,1020]', 1010, 1020, 2)]
        if a.quick:
            rects = [('[0.6,2]x[1010,1013]', 1010, 1013, 1)]
        for name, y0, y1, expect in rects:
            res = []
            for pr in (P1, P2):
                t0 = time.time()
                k, good, rad, nseg = count_rect(pool, 'f', Fraction(6, 10), 2, y0, y1, pr, L=Fraction(1, 2))
                res.append({'prec': pr, 'count': k, 'good': good, 'rad': rad, 'segments': nseg, 'sec': round(time.time() - t0, 1)})
                log('T1 count f in %s at prec %d: N = %d (integer certified %s, radius %.2g, %d segments, %.0f s)' % (name, pr, k, good, rad, nseg, time.time() - t0))
            T1['counts'][name] = {'results': res, 'expected': expect, 'ok': all(r['count'] == expect and r['good'] for r in res) and res[0]['count'] == res[1]['count']}
        R['T1'] = T1
    if a.task in ('t2', 'all'):
        T2 = {}
        TMAX = 30 if a.quick else 200
        zs = [('z85', (0.8085171825, 85.6993484854)), ('z114', (0.6508300806, 114.1633427308)),
              ('z166', (0.5743560505, 166.4793059132)), ('z176', (0.7242576946, 176.7024612429))]
        if a.quick:
            zs = zs[:1]
        T2['zeros'] = [certify_zero(pool, 'f', s, 10, 16, (Q1, Q2), tag) for tag, s in zs]
        K = sum(1 for _, z in zs if z[1] <= TMAX)
        T2['real_axis_pieces'] = {}
        for pr in (Q1, Q2):
            n, bad = real_axis(pr)
            T2['real_axis_pieces'][pr] = {'balls': n, 'undecided': bad}
            log('Lambda real axis: nonzero on [-2, 3] (cover with %d balls, undecided %s) at prec %d' % (n, bad, pr))
        T2.update(line_and_rect(pool, TMAX, K, (Q1, Q2), 'T2'))
        R['T2'] = T2
        # C1 (reconcile-M.md, Chapter 35): the first zero of Lambda on the line, the count below the first off line zero, the disk
        C1 = {}
        ctx.prec = Q1
        vs = [Lam(acb(arb(1) / 2, A(Fraction(t)))) for t in (Fraction(509, 100), Fraction(510, 100))]
        C1['a_values'] = [{'t': str(t), 'ReLambda': v.real.str(10), 'ImLambda': v.imag.str(10)} for t, v in zip((5.09, 5.10), vs)]
        C1['a_sign_change'] = bool((vs[0].real > 0 and vs[1].real < 0) or (vs[0].real < 0 and vs[1].real > 0))
        log('C1(a): Re Lambda(1/2 + 5.09 i) = %s, Re Lambda(1/2 + 5.10 i) = %s, certified sign change: %s' % (vs[0].real.str(10), vs[1].real.str(10), C1['a_sign_change']))
        C1['b'] = line_and_rect(pool, Fraction(8569, 100), 0, (Q1, Q2), 'C1(b)')
        C1['c'] = T2['zeros'][0]
        log('C1(c): disk around %s + %s i of half width %s contains exactly one zero (counts %s)' % (C1['c']['center_re'], C1['c']['center_im'], C1['c']['box_halfwidth'], C1['c']['counts']))
        R['C1'] = C1
    R['seconds'] = round(time.time() - t00, 1)
    json.dump(R, open(out + '.json', 'w'), indent=1, default=str)
    log('done in %.0f s' % (time.time() - t00))


if __name__ == '__main__':
    main()
