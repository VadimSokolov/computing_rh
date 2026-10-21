# zh_window_eigs.py: Table zh:tab:hp (the least eigenvalues of Weil's form on the windows [-L, L] at 60 digits), the second even eigenvalue, Zhu's reference floors, the zero side check of the normalisation and the damping ladder with the crossovers of Section zh:sec:damp, in Chapter 18 (Zhu's compact windows and the Landau Widom horizon).
# Added for the book in October 2026, written from the mathematics (the chapter author's own scripts were neither read nor run). Uses the shared module code/yo_matrix.py, which builds the prime side matrix of Weil's damped form Q_y(f) = W(e^(-y|x|) f * f~), y = alpha - 1/2, on L^2 orthonormal bases of the window [-a, a]; the half length a of the module is the chapter's L. Basis (the chapter's 'modes up to the order N'): the even sector has the N + 1 functions (2L)^(-1/2) and L^(-1/2) cos(pi j u/L), j = 1..N, the odd sector the N functions L^(-1/2) sin(pi j u/L), j = 1..N (module basis 'cossin' with size N + 1 and N); group Y2 reproduced the eight printed entries at N = 8 and 14 with this convention, the entries at N = 20 and 26 are tested here.
# Statements checked. T: the 16 printed entries of Table zh:tab:hp (C:63 to C:70), the least even and odd eigenvalue at alpha = 1/2 for L = 0.5 and 0.8 and N = 8, 14, 20, 26, and the statements of C:73 and C:77 (the values decrease with N; at L = 0.5 they pass below Zhu's upper bound 1.05e-6; at L = 0.8 the even value lies in [8.9e-18, 2.27e-17] and the odd value in [8.206e-15, 2.347e-14] from N = 20 on; the odd sector lies three orders of magnitude above the even one); the ladder is extended to N = 32, 40, 48, 64, 80, 96 and, in the local chunks, 128 and 160 (larger N is for the cluster). S: the second even eigenvalue at L = 0.8 against Zhu's lower bound 2.085e-12 and his Table 1. R: Zhu's Table 1 reference floors (zhu.txt:897 to 906) against a power law extrapolation of the ladder in N. Z: the normalisation, by the zero side: for the unit vector f of the bottom eigenvalue the partial sums 2 sum_{k <= K} F(gamma_k)^2 over the first K ordinates of data/zeros_1700.txt (F the Fourier transform of f) increase to the eigenvalue (Zhu, zhu.txt:838 to 841, K = 300: 1.767e-16 against 2.5223e-16 and 2.071e-14 against 2.3470e-14); K = 100, 300, 400, 1700, with a tail estimate from the density of the ordinates. D: the damping (C:208 to C:212): the least eigenvalue lambda_min(alpha, L) at alpha - 1/2 = 1e-6 .. 1e-4 for L = 0.5 and 3e-17 .. 1e-15 for L = 0.8 (and further values, among them the printed crossovers), against lambda_min(1/2, L), the floor 2 xi'(alpha)/xi(alpha) (Theorem yo:thm:coercive, under RH) and 4 kappa_1 (alpha - 1/2) with kappa_1 = sum over gamma > 0 of gamma^(-2); the crossover is computed as the damping at which the damped least eigenvalue is twice the undamped one.
# Method. Table and second eigenvalues: ym.spectrum on each sector at dps = 60 (mpmath at 75 digits inside the module), the three lowest eigenvalues of each sector are kept. Convergence in N: the three point power law lambda_N = lambda_inf + C N^(-p) through each triple of consecutive sizes (p from the ratio of the two differences by bisection, lambda_inf then in closed form); it is a fit, not a theorem, and the run of p and lambda_inf over the triples is written so that its drift can be judged. Zero side: the eigenvector by mpmath.eigsy on the matrix of ym.build_matrix (unit Euclidean norm, so unit L^2 norm), F(gamma) in closed form, for the even basis sin(gamma L) times (2L)^(-1/2) 2/gamma for j = 0 and L^(-1/2) (-1)^j 2 gamma/(gamma^2 - k_j^2) for j >= 1, for the odd basis i sin(gamma L) L^(-1/2) (-1)^j 2 k_j/(gamma^2 - k_j^2), with k_j = pi j/L; the sum is 2 sum F(gamma_k)^2 over the first K ordinates (|F|^2 is used, so the factor i is immaterial). The tail estimate is 2 times the integral over t > gamma_K of |F(t)|^2 log(t/2 pi)/(2 pi) dt, with the mean value 1/2 for sin^2(tL) beyond twice the largest k_j and the exact integrand below it (not given when gamma_K is below 1.2 times the largest k_j). Damping: ym.spectrum with y = alpha - 1/2 given as a decimal string; the floor 2 xi'(alpha)/xi(alpha) by ym.b_alpha at alpha = 1/2 + c (the mpmath formula 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'/zeta(s)); kappa_1 = xi'(1/2 + c)/xi(1/2 + c)/(2c) at c = 1e-25 (the relative error is (kappa_2/kappa_1) c^2, about 2e-53, and the cancellation costs 25 of the 70 digits); the crossover by bisection in log c, 11 steps, on the even sector (the odd sector lies far above), with a final evaluation of both sectors.
# Precision and error control. Every eigenvalue is computed at dps = 60 (75 digits inside the module); part prec repeats the least eigenvalue of both sectors at dps + 30 (the module's default number of Gauss Legendre nodes grows with dps) and with 1.5 times the default nodes at dps = 60, and writes the relative changes (the Gauss Legendre error of the archimedean entries is what MINOR 12 of the review asks to control). Further checks written to the JSON: the Rayleigh quotient of the eigenvector against the eigenvalue, the unit norm, the closed form of F against direct quadrature at two ordinates, the monotone decrease in N (the bases are nested, so the least eigenvalue of each sector cannot increase and each value is an upper bound for the window infimum up to the numerical error, which is far below the printed digits), the dependence on N through the ladder, and the agreement with the printed digits (the printed number is compared with ours rounded to the same number of significant digits). Nothing here is a certificate: the values are floating point upper bounds for the infimum of the form over the window, and the zero side sums assume RH for the ordinates used (verified to height 3e12) and the Hadamard type tail for the rest.
# Output: data/zh_window_eigs.json (numbers as strings with 15 to 30 digits, with size and precision; the JSON is updated in place, so partial runs add up) and data/zh_window_eigs.log. Keys: table, fits, prec, zeros, damping, claims, meta. Usage (repository root): timeout 60 python3 -I code/zh_window_eigs.py --quick [--out FILE] (L = 0.5, N = 8 at 40 digits, about 2 s, writes to /tmp unless --out is given); the parts are table, prec, zeros, damping and claims, with --L, --K (list of N), --dps, --zk (list of K for the zero sums), --procs P, --zeros FILE. The local chunks that made data/zh_window_eigs.json (each under 25 s on one core of a laptop, in this order of cost): --part table --L 0.5; --part table --L 0.8 --K 8,14,20,26,32,40,48,64; --part table --L 0.8 --K 80,96,128; --part table --L 0.8 --K 160; --part table --L 0.5 --K 160; --part prec --L 0.5 --K 26,40; --part prec --L 0.8 --K 26,40; --part zeros --L 0.5 --K 26; --part zeros --L 0.8 --K 26,40; --part zeros --L 0.5 --K 40; --part damping --L 0.5 --K 26; --part damping --L 0.8 --K 26; --part damping --L 0.5 --K 40; --part damping --L 0.8 --K 40; --part claims. The complete run is for the Hopper cluster: bash misc/tools/hopper_run.sh -c 8 -t 240 -g zh_window_eigs.json -g zh_window_eigs.log code/zh_window_eigs.py code/yo_matrix.py data/zeros_1700.txt -- --part all --procs 8 --out zh_window_eigs.json
import argparse
import json
import math
import os
import platform
import socket
import sys
import time

import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yo_matrix as ym  # noqa: E402

CHAPTER_N = [8, 14, 20, 26]
LADDER_N = [8, 14, 20, 26, 32, 40, 48, 64, 80, 96]
PRINTED = {0.5: dict(lines=[63, 64, 65, 66], even=['1.336e-6', '1.082e-6', '1.005e-6', '9.743e-7'], odd=['2.51e-4', '2.39e-4', '2.30e-4', '2.24e-4']),
           0.8: dict(lines=[67, 68, 69, 70], even=['5.49e-15', '4.127e-17', '2.651e-17', '2.144e-17'], odd=['2.67e-12', '3.86e-14', '2.333e-14', '2.157e-14'])}
# Zhu (zhu.txt): Table 2 at L = 0.5 (line 1307), Theorem 1.2 and the certified odd interval at L = 0.8, the lower bound of the second even eigenvalue, Table 1 (lines 897 to 906: even, odd, second even), the partial sums at K = 300 (lines 838 to 841).
ZHU = dict(even_up_05='1.05e-6', even_lo_08='8.9e-18', even_up_08='2.27e-17', odd_lo_08='8.206e-15', odd_up_08='2.347e-14', even2_lo_08='2.085e-12')
ZHU_TABLE1 = {0.5: dict(even='9.34e-7', odd='1.94e-4', even2='1.80e-2'), 0.8: dict(even='1.65e-17', odd='1.57e-14', even2='8.38e-12')}
ZHU_K300 = dict(even=('1.767e-16', '2.5223e-16'), odd=('2.071e-14', '2.3470e-14'))
LADDER_C = {0.5: ['1e-6', '3e-6', '1e-5', '3e-5', '1e-4'], 0.8: ['3e-17', '1e-16', '3e-16', '1e-15']}
EXTRA_C = {0.5: ['1e-7', '3e-7', '1.1e-5', '3e-4', '1e-3'], 0.8: ['1e-17', '2.5e-16', '1e-14', '1e-12', '1e-10', '1e-8', '1e-6', '1e-4']}
CHAPTER_CROSSOVER = {0.5: '1.1e-5', 0.8: '2.5e-16'}


class Log:
    def __init__(self, path):
        self.fh = open(path, 'a')

    def __call__(self, msg=''):
        print(msg, flush=True)
        self.fh.write(msg + '\n')
        self.fh.flush()


def default_path(stem):
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    if os.path.basename(here) == 'code' and os.path.isdir(os.path.join(root, 'data')):
        return os.path.join(root, 'data', stem)
    return os.path.join(os.getcwd(), stem)


def find_zeros(given=None):
    here = os.path.dirname(os.path.abspath(__file__))
    for p in ([given] if given else []) + [os.path.join(os.path.dirname(here), 'data', 'zeros_1700.txt'), os.path.join(here, 'zeros_1700.txt'), 'zeros_1700.txt']:
        if p and os.path.exists(p):
            return p
    return None


def read_zeros(path, wp):
    out = []
    with mpmath.workdps(wp):
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    out.append(mpmath.mpf(line.split()[0]))
    return out


def load(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save(path, res):
    with open(path, 'w') as fh:
        json.dump(res, fh, indent=1)


def stamp(res, part, seconds, mode):
    res.setdefault('meta', {})
    res['meta'][part] = dict(seconds=round(seconds, 1), host=socket.gethostname(), mode=mode, date=time.strftime('%Y-%m-%d %H:%M'), python=platform.python_version(),
                             mpmath=mpmath.__version__, slurm_job=os.environ.get('SLURM_JOB_ID'))


def ns(x, n=25):
    return mpmath.nstr(x, n)


def sector_size(K, sector):
    return K + 1 if sector == 'even' else K


def spec(L, c, K, sector, dps, nodes=None):
    """All eigenvalues (ascending, mpmath numbers) of the sector block of order sector_size(K, sector) at the damping c (a decimal string, '0' for none)."""
    return ym.spectrum(str(L), y=c, basis='cossin', size=sector_size(K, sector), sector=sector, dps=dps, nodes=nodes)


def printed_digits(s):
    mant = s.lower().split('e')[0]
    return len(mant.replace('.', '').replace('-', '').lstrip('0'))


def agrees(printed, value):
    n = printed_digits(printed)
    return '%.*e' % (n - 1, float(printed)) == '%.*e' % (n - 1, float(value))


def run_pool(procs, func, tasks, callback):
    if procs > 1 and len(tasks) > 1:
        import multiprocessing
        with multiprocessing.Pool(min(procs, len(tasks))) as pool:
            for r in pool.imap_unordered(func, tasks):
                callback(r)
    else:
        for t in tasks:
            callback(func(t))


# ---------------------------------------------------------------------------------------------------------------------------------
# T, S, R: the table, the second eigenvalues, the convergence in N
# ---------------------------------------------------------------------------------------------------------------------------------

def table_case(arg):
    L, K, dps = arg
    t0 = time.time()
    r = dict(L=L, K=K, dps=dps, even_size=K + 1, odd_size=K)
    for sec in ('even', 'odd'):
        ev = spec(L, '0', K, sec, dps)
        r[sec] = [ns(v, 30) for v in ev[:3]]
    r['seconds'] = round(time.time() - t0, 2)
    return r


def powerlaw_fit(ks, vals):
    """lambda_N = lambda_inf + C N^(-p) through three points (N1 < N2 < N3, decreasing values): (p, lambda_inf) or (None, None)."""
    n1, n2, n3 = ks
    v1, v2, v3 = vals
    d1, d2 = v1 - v2, v2 - v3
    if not (d1 > 0 and d2 > 0):
        return None, None

    def ratio(p):
        return (n1 ** -p - n2 ** -p) / (n2 ** -p - n3 ** -p)
    lo, hi = 1e-3, 40.0
    if not (ratio(lo) <= d1 / d2 <= ratio(hi)):
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ratio(mid) < d1 / d2:
            lo = mid
        else:
            hi = mid
    p = 0.5 * (lo + hi)
    return p, v3 - d2 * n3 ** -p / (n2 ** -p - n3 ** -p)


def part_fits(res, log):
    fits = res.setdefault('fits', {})
    tab = res.get('table', {})
    for L in sorted({r['L'] for r in tab.values()}):
        rows = sorted((r for r in tab.values() if r['L'] == L), key=lambda r: r['K'])
        for sec, idx in (('even', 0), ('odd', 0), ('even', 1)):
            ks = [r['K'] for r in rows]
            vals = [float(mpmath.mpf(r[sec][idx])) for r in rows]
            ent = dict(K=ks, values=[r[sec][idx] for r in rows], monotone_decreasing=all(vals[i + 1] <= vals[i] for i in range(len(vals) - 1)),
                       rel_change_to_next=[ns(mpmath.mpf(rows[i + 1][sec][idx]) / mpmath.mpf(rows[i][sec][idx]) - 1, 4) for i in range(len(rows) - 1)], triples=[])
            for i in range(len(ks) - 2):
                p, inf = powerlaw_fit(ks[i:i + 3], vals[i:i + 3])
                ent['triples'].append(dict(K=ks[i:i + 3], p=None if p is None else round(p, 4), lam_inf=None if inf is None else '%.5e' % inf))
            fits['L%g_%s_%d' % (L, sec, idx + 1)] = ent
            tail = ent['triples'][-1] if ent['triples'] else None
            log('fit L=%g %s eigenvalue %d: K=%s monotone=%s, last triple %s' % (L, sec, idx + 1, ks, ent['monotone_decreasing'], tail))


def part_table(res, log, Ls, Ks, dps, procs, path):
    t0 = time.time()
    tab = res.setdefault('table', {})
    tasks = [(L, K, dps) for L in Ls for K in sorted(Ks, reverse=True)]

    def done(r):
        tab['%g_%d' % (r['L'], r['K'])] = r
        log('table L=%g K=%d dps=%d (sizes %d, %d): even %s %s odd %s %s (%.1f s)' % (r['L'], r['K'], r['dps'], r['even_size'], r['odd_size'], mpmath.nstr(mpmath.mpf(r['even'][0]), 10),
                                                                                  mpmath.nstr(mpmath.mpf(r['even'][1]), 6), mpmath.nstr(mpmath.mpf(r['odd'][0]), 10), mpmath.nstr(mpmath.mpf(r['odd'][1]), 6), r['seconds']))
        save(path, res)
    run_pool(procs, table_case, tasks, done)
    part_fits(res, log)
    stamp(res, 'table', time.time() - t0, 'cluster' if os.environ.get('SLURM_JOB_ID') else 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
# precision and quadrature control
# ---------------------------------------------------------------------------------------------------------------------------------

def prec_case(arg):
    L, K, c, dps = arg
    t0 = time.time()
    r = dict(L=L, K=K, c=c, dps=dps, dps_tested=dps + 30)
    try:
        nq0 = int(ym._nq_default('fourier', mpmath.mpf(str(L)), K, dps))
    except Exception:
        nq0 = None
    for sec in ('even', 'odd'):
        base = spec(L, c, K, sec, dps)[0]
        hi = spec(L, c, K, sec, dps + 30)[0]
        ent = dict(lam=ns(base, 30), rel_diff_dps=ns(abs(hi - base) / base, 3))
        if nq0 is not None:
            more = spec(L, c, K, sec, dps, nodes=int(1.5 * nq0))[0]
            ent.update(nodes_default=nq0, nodes_tested=int(1.5 * nq0), rel_diff_nodes=ns(abs(more - base) / base, 3))
        r[sec] = ent
    r['seconds'] = round(time.time() - t0, 1)
    return r


def part_prec(res, log, Ls, Ks, dps, procs, path):
    t0 = time.time()
    out = res.setdefault('prec', {})
    tasks = [(L, K, c, dps) for L in Ls for K in Ks for c in ('0', '1e-5' if L == 0.5 else '1e-16')]

    def done(r):
        out['L%g_K%d_c%s' % (r['L'], r['K'], r['c'])] = r
        log('prec L=%g K=%d c=%s dps=%d vs %d: even %s (dps %s, nodes %s), odd %s (dps %s, nodes %s) (%.1f s)' % (
            r['L'], r['K'], r['c'], r['dps'], r['dps_tested'], mpmath.nstr(mpmath.mpf(r['even']['lam']), 8), r['even']['rel_diff_dps'], r['even'].get('rel_diff_nodes'),
            mpmath.nstr(mpmath.mpf(r['odd']['lam']), 8), r['odd']['rel_diff_dps'], r['odd'].get('rel_diff_nodes'), r['seconds']))
        save(path, res)
    run_pool(procs, prec_case, tasks, done)
    stamp(res, 'prec', time.time() - t0, 'cluster' if os.environ.get('SLURM_JOB_ID') else 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
# Z: the zero side of the normalisation
# ---------------------------------------------------------------------------------------------------------------------------------

def bottom_vector(L, K, sector, dps):
    """Lowest eigenvalue and unit eigenvector (mpmath, 75 digits) of the sector matrix at alpha = 1/2, with the check of the Rayleigh quotient."""
    M = ym.build_matrix(str(L), y='0', basis='cossin', size=sector_size(K, sector), sector=sector, dps=dps)
    with mpmath.workdps(dps + ym.GUARD):
        E, Q = mpmath.eigsy(M)
        i0 = min(range(len(E)), key=lambda i: E[i])
        vec = [Q[r, i0] for r in range(Q.rows)]
        nrm = mpmath.sqrt(mpmath.fsum(v * v for v in vec))
        vec = [v / nrm for v in vec]
        n = len(vec)
        rq = mpmath.fsum(vec[i] * M[i, j] * vec[j] for i in range(n) for j in range(n))
        lam = E[i0]
        return lam, vec, abs(rq - lam) / lam, abs(mpmath.fsum(v * v for v in vec) - 1)


def make_G(L, sector, vec):
    """G(t) with |F(t)| = |sin(t L) G(t)|, F the Fourier transform of f = sum vec_i basis_i, at the current precision."""
    a = mpmath.mpf(str(L))
    pi = mpmath.pi
    ra, r2a = mpmath.sqrt(a), mpmath.sqrt(2 * a)
    n = len(vec)
    js = list(range(n)) if sector == 'even' else list(range(1, n + 1))
    kk = [pi * j / a for j in js]
    k2 = [x * x for x in kk]
    sg = [1 if j % 2 == 0 else -1 for j in js]

    def G(t):
        t2 = t * t
        if sector == 'even':
            tot = vec[0] * 2 / (t * r2a)
            for i in range(1, n):
                tot += vec[i] * sg[i] * 2 * t / (t2 - k2[i]) / ra
        else:
            tot = mpmath.mpf(0)
            for i in range(n):
                tot += vec[i] * sg[i] * 2 * kk[i] / (t2 - k2[i]) / ra
        return tot
    return G


def closed_form_check(L, sector, vec, gam):
    """Largest relative difference between the closed form of F and direct quadrature of the Fourier integral at two ordinates (30 digits)."""
    with mpmath.workdps(30):
        a = mpmath.mpf(str(L))
        pi = mpmath.pi
        n = len(vec)

        def f(u):
            if sector == 'even':
                tot = vec[0] / mpmath.sqrt(2 * a)
                for j in range(1, n):
                    tot += vec[j] * mpmath.cos(pi * j * u / a) / mpmath.sqrt(a)
            else:
                tot = mpmath.mpf(0)
                for j in range(1, n + 1):
                    tot += vec[j - 1] * mpmath.sin(pi * j * u / a) / mpmath.sqrt(a)
            return tot
        G = make_G(L, sector, vec)
        worst = mpmath.mpf(0)
        for g in (gam[0], gam[4]):
            g = +g
            if sector == 'even':
                fq = 2 * mpmath.quad(lambda u: f(u) * mpmath.cos(g * u), [0, a / 2, a])
            else:
                fq = 2 * mpmath.quad(lambda u: f(u) * mpmath.sin(g * u), [0, a / 2, a])
            fc = mpmath.sin(g * a) * G(g)
            worst = max(worst, abs(fq - fc) / abs(fc))
        return worst


def tail_estimate(L, sector, vec, T, n_modes):
    """2 int_T^inf |F(t)|^2 log(t/2 pi)/(2 pi) dt (mean value 1/2 of sin^2 beyond twice the largest k_j, the exact integrand below it); None if T < 1.2 k_max."""
    with mpmath.workdps(20):
        a = mpmath.mpf(str(L))
        pi = mpmath.pi
        kmax = pi * n_modes / a
        T = mpmath.mpf(T)
        if T < 1.2 * kmax:
            return None
        G = make_G(L, sector, vec)
        T2 = max(T, 2 * kmax)
        total = mpmath.mpf(0)
        if T2 > T:
            pts = [T + (T2 - T) * i / 8 for i in range(9)]
            total += 2 * mpmath.quad(lambda t: (mpmath.sin(t * a) * G(t)) ** 2 * mpmath.log(t / (2 * pi)) / (2 * pi), pts, maxdegree=5)
        total += mpmath.quad(lambda t: G(t) ** 2 * mpmath.log(t / (2 * pi)) / (2 * pi), [T2, 3 * T2, 30 * T2, mpmath.inf], maxdegree=5)
        return total


def part_zeros(res, log, Ls, Ks, dps, zeros_path, zks, path):
    t0 = time.time()
    if zeros_path is None:
        log('zeros: data/zeros_1700.txt not found, part skipped')
        return
    wp = dps + ym.GUARD
    gam = read_zeros(zeros_path, wp)
    out = res.setdefault('zeros', {})
    for L in Ls:
        for K in Ks:
            for sec in ('even', 'odd'):
                t1 = time.time()
                lam, vec, rq_err, norm_err = bottom_vector(L, K, sec, dps)
                with mpmath.workdps(wp):
                    G = make_G(L, sec, vec)
                    a = mpmath.mpf(str(L))
                    F2 = [(mpmath.sin(g * a) * G(g)) ** 2 for g in gam]
                row = dict(L=L, K=K, sector=sec, dps=dps, size=len(vec), n_ordinates=len(gam), lam=ns(lam, 25), rayleigh_minus_eig_rel=ns(rq_err, 3), unit_norm_error=ns(norm_err, 3),
                           closed_form_vs_quadrature_rel=ns(closed_form_check(L, sec, vec, gam), 3), partial={})
                for kk in zks:
                    if kk > len(gam):
                        continue
                    with mpmath.workdps(wp):
                        zs = 2 * mpmath.fsum(F2[:kk])
                        d = lam - zs
                        te = tail_estimate(L, sec, vec, gam[kk - 1], len(vec))
                        row['partial'][str(kk)] = dict(gamma_K=ns(gam[kk - 1], 12), zs=ns(zs, 15), zs_over_lam=ns(zs / lam, 10), lam_minus_zs=ns(d, 8), rel_distance=ns(d / lam, 8),
                                                       tail_estimate=None if te is None else ns(te, 6), distance_over_tail=None if te is None else ns(d / te, 6))
                out['L%g_K%d_%s' % (L, K, sec)] = row
                log('zeros L=%g K=%d %s: lam=%s, Rayleigh rel %s, closed form rel %s (%.1f s)' % (L, K, sec, ns(lam, 12), row['rayleigh_minus_eig_rel'], row['closed_form_vs_quadrature_rel'], time.time() - t1))
                for kk in zks:
                    p = row['partial'].get(str(kk))
                    if p:
                        log('   K=%-5d gamma_K=%-12s 2 sum F^2 = %s, ratio to lam %s, distance %s (relative %s), tail estimate %s, distance/tail %s' % (
                            kk, p['gamma_K'][:9], mpmath.nstr(mpmath.mpf(p['zs']), 8), mpmath.nstr(mpmath.mpf(p['zs_over_lam']), 8), mpmath.nstr(mpmath.mpf(p['lam_minus_zs']), 5),
                            mpmath.nstr(mpmath.mpf(p['rel_distance']), 5), p['tail_estimate'], p['distance_over_tail']))
                save(path, res)
    stamp(res, 'zeros', time.time() - t0, 'cluster' if os.environ.get('SLURM_JOB_ID') else 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
# D: the damping ladder and the crossover
# ---------------------------------------------------------------------------------------------------------------------------------

def kappa1(dps):
    with mpmath.workdps(dps + 30):
        c = mpmath.mpf(10) ** -25
        return ym.b_alpha(mpmath.mpf('0.5') + c, dps=dps) / (2 * c)


def floor_2b(c, dps):
    with mpmath.workdps(dps + 30):
        return 2 * ym.b_alpha(mpmath.mpf('0.5') + mpmath.mpf(c), dps=dps)


def lam_even(L, c, K, dps):
    return spec(L, c, K, 'even', dps)[0]


def find_factor(L, K, dps, lam0, rows, factor=2):
    """The damping c with lam_even(c) = factor lam0 by bisection in log c, bracketed from the ladder rows (11 steps)."""
    pts = sorted((float(mpmath.mpf(r['c'])), float(mpmath.mpf(r['even']) / lam0)) for r in rows)
    lo = max((c for c, q in pts if q < factor), default=None)
    hi = min((c for c, q in pts if q >= factor and (lo is None or c > lo)), default=None)
    steps = 0
    while lo is None and steps < 12:
        c = pts[0][0] / 3 ** (steps + 1)
        if float(lam_even(L, '%.6e' % c, K, dps) / lam0) < factor:
            lo = c
        steps += 1
    steps = 0
    while hi is None and steps < 12:
        c = pts[-1][0] * 3 ** (steps + 1)
        if float(lam_even(L, '%.6e' % c, K, dps) / lam0) >= factor:
            hi = c
        steps += 1
    if lo is None or hi is None:
        return None
    a, b = math.log(lo), math.log(hi)
    for _ in range(11):
        mid = 0.5 * (a + b)
        if float(lam_even(L, '%.8e' % math.exp(mid), K, dps) / lam0) < factor:
            a = mid
        else:
            b = mid
    return '%.4e' % math.exp(0.5 * (a + b))


def damp_case(arg):
    L, K, dps, ladder, extra, search = arg
    t0 = time.time()
    k1 = kappa1(dps)
    base = {sec: spec(L, '0', K, sec, dps)[:2] for sec in ('even', 'odd')}
    lam0 = min(base['even'][0], base['odd'][0])
    rows = []
    for c in ladder + extra:
        secs = ('even', 'odd') if c in ladder else ('even',)
        pr = {sec: spec(L, c, K, sec, dps)[:2] for sec in secs}
        lam = min(v[0] for v in pr.values())
        fl = floor_2b(c, dps)
        with mpmath.workdps(dps + 15):
            cc = mpmath.mpf(c)
            rows.append(dict(c=c, in_ladder=(c in ladder), lam=ns(lam, 25), even=ns(pr['even'][0], 25), odd=ns(pr['odd'][0], 25) if 'odd' in pr else None, sector='even' if lam == pr['even'][0] else 'odd',
                             floor_2b=ns(fl, 15), four_kappa1_c=ns(4 * k1 * cc, 15), lam_over_floor=ns(lam / fl, 8), lam_over_four_kappa1_c=ns(lam / (4 * k1 * cc), 8),
                             lam_over_lam0=ns(lam / lam0, 8), lam_over_max_lam0_4k1c=ns(lam / max(lam0, 4 * k1 * cc), 8)))
    rows.sort(key=lambda r: float(mpmath.mpf(r['c'])))
    out = dict(L=L, K=K, dps=dps, kappa1=ns(k1, 25), lam0=ns(lam0, 25), lam0_even=ns(base['even'][0], 25), lam0_odd=ns(base['odd'][0], 25), lam0_even2=ns(base['even'][1], 25), rows=rows,
               crossover_chapter_formula_from_lam0=ns(lam0 / (4 * k1), 6))
    with mpmath.workdps(dps + 15):
        small = [r for r in rows if r['in_ladder']][0]
        out['first_order_slope_over_4kappa1'] = ns((mpmath.mpf(small['lam']) - lam0) / mpmath.mpf(small['c']) / (4 * k1), 6)
        out['first_order_slope_from_c'] = small['c']
    if search:
        c2 = find_factor(L, K, dps, base['even'][0], rows)
        if c2 is not None:
            fin = {sec: spec(L, c2, K, sec, dps)[0] for sec in ('even', 'odd')}
            out['crossover_factor2'] = dict(c=c2, lam_min=ns(min(fin.values()), 12), ratio_to_lam0=ns(min(fin.values()) / lam0, 8), over_chapter_formula=ns(mpmath.mpf(c2) / (lam0 / (4 * k1)), 6))
        else:
            out['crossover_factor2'] = None
    out['seconds'] = round(time.time() - t0, 1)
    return out


def part_damping(res, log, Ls, Ks, dps, procs, path, quick=False):
    t0 = time.time()
    out = res.setdefault('damping', {})
    tasks = []
    for L in Ls:
        ladder = ['1e-5'] if quick else LADDER_C[L]
        extra = [] if quick else EXTRA_C[L]
        for K in sorted(Ks, reverse=True):
            tasks.append((L, K, dps, ladder, extra, True))

    def done(r):
        out['L%g_K%d' % (r['L'], r['K'])] = r
        log('damping L=%g K=%d dps=%d: kappa1 = %s, lam(1/2) = %s (even %s, odd %s, even2 %s), chapter formula lam0/(4 kappa1) = %s (%.1f s)' % (
            r['L'], r['K'], r['dps'], ns(mpmath.mpf(r['kappa1']), 12), ns(mpmath.mpf(r['lam0']), 8), ns(mpmath.mpf(r['lam0_even']), 8), ns(mpmath.mpf(r['lam0_odd']), 8),
            ns(mpmath.mpf(r['lam0_even2']), 6), r['crossover_chapter_formula_from_lam0'], r['seconds']))
        for w in r['rows']:
            log('   c=%-8s %s lam=%s (%s sector) floor 2b=%s, 4 kappa1 c=%s, lam/floor=%s, lam/(4 kappa1 c)=%s, lam/lam0=%s' % (
                w['c'], '*' if w['in_ladder'] else ' ', ns(mpmath.mpf(w['lam']), 8), w['sector'], ns(mpmath.mpf(w['floor_2b']), 6), ns(mpmath.mpf(w['four_kappa1_c']), 6), ns(mpmath.mpf(w['lam_over_floor']), 6),
                ns(mpmath.mpf(w['lam_over_four_kappa1_c']), 6), ns(mpmath.mpf(w['lam_over_lam0']), 6)))
        log('   first order slope over 4 kappa1 (from c = %s): %s; factor 2 crossover: %s' % (r['first_order_slope_from_c'], r['first_order_slope_over_4kappa1'], r.get('crossover_factor2')))
        save(path, res)
    run_pool(procs, damp_case, tasks, done)
    stamp(res, 'damping', time.time() - t0, 'cluster' if os.environ.get('SLURM_JOB_ID') else 'local')


# ---------------------------------------------------------------------------------------------------------------------------------
# the claims, mechanically
# ---------------------------------------------------------------------------------------------------------------------------------

def part_claims(res, log):
    tab = res.get('table', {})
    fits = res.get('fits', {})
    dm = res.get('damping', {})
    zr = res.get('zeros', {})
    cl = []
    AG, DF, PEND = 'agrees to the printed digits', 'differs', 'PENDING'

    def add(cid, line, text, chap, ours, status, note=''):
        cl.append(dict(id=cid, chapter_line=line, claim=text, chapter_value=chap, ours=ours, status=status, note=note))
        log('claim %-14s C:%-8s %-66s chapter %-18s ours %-26s %s%s' % (cid, line, text, chap, ours, status, (' | ' + note) if note else ''))

    def lam(L, sec, K, idx=0):
        r = tab.get('%g_%d' % (L, K))
        return None if r is None else mpmath.mpf(r[sec][idx])

    def sizes(L):
        return sorted(r['K'] for r in tab.values() if r['L'] == L)

    # T: the sixteen printed entries
    for L in (0.5, 0.8):
        P = PRINTED[L]
        for i, K in enumerate(CHAPTER_N):
            for sec in ('even', 'odd'):
                v = lam(L, sec, K)
                cid, text = 'T%s%g_%d' % (sec[0], L, K), 'L=%g N=%d %s least eigenvalue' % (L, K, sec)
                if v is None:
                    add(cid, P['lines'][i], text, P[sec][i], 'n/a', PEND)
                else:
                    ok = agrees(P[sec][i], v)
                    add(cid, P['lines'][i], text, P[sec][i], ns(v, 8), AG if ok else DF, '' if ok else 'relative difference %s' % ns(abs(v - mpmath.mpf(P[sec][i])) / v, 3))
    # the statements of C:73 and C:77
    for L in (0.5, 0.8):
        for sec in ('even', 'odd'):
            ent = fits.get('L%g_%s_1' % (L, sec))
            ok = ent is not None and len(ent['K']) >= 4
            add('mono_%s%g' % (sec[0], L), '73,77', 'values decrease with N (%s, L=%g)' % (sec, L), 'decrease', ('N=%s' % ent['K']) if ent else 'n/a', (AG if ent['monotone_decreasing'] else DF) if ok else PEND)
    if 26 in sizes(0.5):
        v = lam(0.5, 'even', 26)
        below = [K for K in sizes(0.5) if lam(0.5, 'even', K) < mpmath.mpf(ZHU['even_up_05'])]
        add('Z05', '77', 'L=0.5 even below Zhu upper bound 1.05e-6 (at N=26)', '9.743e-7 < 1.05e-6', 'N=26: %s; first N below: %s' % (ns(v, 6), min(below) if below else 'none'), AG if v < mpmath.mpf(ZHU['even_up_05']) else DF)
        q5 = lam(0.5, 'odd', 26) / v
        add('ratio05', '-', 'odd over even ratio at L=0.5, N=26 (information, no claim in the chapter)', 'none', ns(q5, 5), 'information', 'Zhu Table 1 ratio: 208')
    else:
        add('Z05', '77', 'L=0.5 even below Zhu upper bound 1.05e-6 (at N=26)', '9.743e-7 < 1.05e-6', 'n/a', PEND)
    if 26 in sizes(0.8):
        ks = sizes(0.8)
        lo_e, up_e = mpmath.mpf(ZHU['even_lo_08']), mpmath.mpf(ZHU['even_up_08'])
        lo_o, up_o = mpmath.mpf(ZHU['odd_lo_08']), mpmath.mpf(ZHU['odd_up_08'])
        inside_e = [K for K in ks if lo_e <= lam(0.8, 'even', K) <= up_e]
        add('Z08e', '70,73,77', 'L=0.8 even at N=26 inside [8.9e-18, 2.27e-17]', 'inside', 'N inside: %s' % inside_e, AG if 26 in inside_e else DF)
        mn = min(lam(0.8, 'even', K) for K in ks)
        add('Z08e_lo', '70', 'L=0.8 even at every N above the certified lower bound 8.9e-18', 'above', 'least over N=%s: %s' % (ks, ns(mn, 6)), AG if mn >= lo_e else DF)
        inside_o = [K for K in ks if lo_o <= lam(0.8, 'odd', K) <= up_o]
        need = [K for K in ks if K in CHAPTER_N]
        ok_o = all((K in inside_o) == (K >= 20) for K in need)
        add('Z08o', '73,77', 'L=0.8 odd inside [8.206e-15, 2.347e-14] from N=20 on (not at N=8, 14)', 'N=20, 26 in; 8, 14 out', 'N inside: %s' % inside_o, AG if ok_o else DF)
        q8 = lam(0.8, 'odd', 26) / lam(0.8, 'even', 26)
        add('ratio08', '77', 'odd over even ratio at L=0.8, N=26 ("three orders of magnitude")', '~1e3', ns(q8, 5), AG if 500 < q8 < 5000 else DF, 'Zhu Table 1 ratio: 949')
        v2 = [lam(0.8, 'even', K, 1) for K in ks]
        add('S08', 'review B', 'second even eigenvalue at L=0.8 at least Zhu lower bound 2.085e-12', '>= 2.085e-12', 'least over N: %s; N=%d: %s' % (ns(min(v2), 6), ks[-1], ns(v2[-1], 6)), AG if min(v2) >= mpmath.mpf(ZHU['even2_lo_08']) else DF, 'Zhu Table 1 reference: 8.38e-12')
    else:
        for cid, text in (('Z08e', 'L=0.8 even at N=26 inside [8.9e-18, 2.27e-17]'), ('Z08o', 'L=0.8 odd inside [8.206e-15, 2.347e-14] from N=20 on'), ('ratio08', 'odd over even ratio at L=0.8'), ('S08', 'second even eigenvalue at L=0.8 at least 2.085e-12')):
            add(cid, '70-77', text, '-', 'n/a', PEND)
    # R: Zhu's Table 1 reference floors against the ladder and its power law extrapolation
    for L in (0.5, 0.8):
        for key, fk in (('even', 'L%g_even_1'), ('odd', 'L%g_odd_1'), ('even2', 'L%g_even_2')):
            ent = fits.get(fk % L)
            if ent is None or not ent['K']:
                add('R%g_%s' % (L, key), 'zhu.txt:897', 'Zhu Table 1 reference, L=%g %s' % (L, key), ZHU_TABLE1[L][key], 'n/a', PEND)
                continue
            tr = ent['triples'][-1]['lam_inf'] if ent['triples'] else None
            ref = mpmath.mpf(ZHU_TABLE1[L][key])
            lastv = mpmath.mpf(ent['values'][-1])
            add('R%g_%s' % (L, key), 'zhu.txt:897', 'Zhu Table 1 reference, L=%g %s' % (L, key), ZHU_TABLE1[L][key], 'N=%d: %s; power law: %s' % (ent['K'][-1], ns(lastv, 5), tr), AG if abs(lastv / ref - 1) <= mpmath.mpf('0.005') else PEND,
                'N=%d value over Zhu: %s; power law over Zhu: %s; status AG when the upper bound at the largest N is within 0.5 percent of the reference, else PENDING (larger N, cluster)' % (ent['K'][-1], ns(lastv / ref, 4), ns(mpmath.mpf(tr) / ref, 4) if tr else 'n/a'))
    # Z: the zero side partial sums
    for key in sorted(zr):
        r = zr[key]
        kks = sorted(int(x) for x in r['partial'])
        zs = [mpmath.mpf(r['partial'][str(x)]['zs_over_lam']) for x in kks]
        below = all(z < 1 for z in zs)
        incr = all(zs[i + 1] >= zs[i] for i in range(len(zs) - 1))
        add('Z_%s' % key, '-', 'zero side partial sums below lambda and increasing (%s)' % key, 'converge from below', 'K=%s: ratio %s' % (kks, [ns(z, 6) for z in zs]), AG if (below and incr) else DF,
            'closed form F vs quadrature %s, Rayleigh %s' % (r['closed_form_vs_quadrature_rel'], r['rayleigh_minus_eig_rel']))
    # D: the damping, the bound of M11 and the crossovers of C:212
    k1 = kappa1(60)
    for L, zub, zlo in ((0.5, ZHU['even_up_05'], None), (0.8, ZHU['even_up_08'], ZHU['even_lo_08'])):
        cx = mpmath.mpf(zub) / (4 * k1)
        add('Dar%g_up' % L, '212', 'Zhu upper bound %s over 4 kappa1 (arithmetic, kappa1 = %s)' % (zub, ns(k1, 10)), CHAPTER_CROSSOVER[L], ns(cx, 5), AG if agrees(CHAPTER_CROSSOVER[L], cx) else DF)
        if zlo:
            add('Dar%g_lo' % L, '212', 'Zhu lower bound %s over 4 kappa1 (arithmetic)' % zlo, CHAPTER_CROSSOVER[L], ns(mpmath.mpf(zlo) / (4 * k1), 5), 'information', 'review MINOR 25: the crossover lies between this and the one from the upper bound')
        ent = fits.get('L%g_even_1' % L)
        if ent is not None and ent['triples'] and ent['triples'][-1]['lam_inf']:
            lc = mpmath.mpf(ent['triples'][-1]['lam_inf'])
            add('Dar%g_conv' % L, '212', 'converged lambda (power law over N=%s: %s) over 4 kappa1' % (ent['triples'][-1]['K'], ent['triples'][-1]['lam_inf']), CHAPTER_CROSSOVER[L], ns(lc / (4 * k1), 5),
                AG if agrees(CHAPTER_CROSSOVER[L], lc / (4 * k1)) else DF)
    for L in (0.5, 0.8):
        for K in (26, 40):
            ent = dm.get('L%g_K%d' % (L, K))
            if ent is None:
                for cid, text in (('D_floor', 'M11 second term, lambda_min >= 2 xi\'/xi(alpha)'), ('D_first', 'M11 first term, lambda_min >= lambda_min(1/2, L)'), ('D_cross', 'crossover lambda_min(1/2, L)/(4 kappa1)'), ('D_x2', 'damping at which lambda_min doubles')):
                    add('%s%g_%d' % (cid, L, K), '208-212', text + ' (N=%d)' % K, CHAPTER_CROSSOVER[L] if cid in ('D_cross', 'D_x2') else 'holds', 'n/a', PEND)
                continue
            lam0 = mpmath.mpf(ent['lam0'])
            rows = ent['rows']
            minf = min(mpmath.mpf(w['lam_over_floor']) for w in rows)
            add('D_floor%g_%d' % (L, K), '208,210', 'M11 second term, lambda_min(alpha, L) >= 2 xi\'/xi(alpha), N=%d, %d values of c' % (K, len(rows)), 'holds (RH)', 'least lam/floor: %s' % ns(minf, 5), AG if minf >= 1 else DF)
            min0 = min(mpmath.mpf(w['lam_over_lam0']) for w in rows)
            add('D_first%g_%d' % (L, K), '210', 'M11 first term, lambda_min(alpha, L) >= lambda_min(1/2, L), N=%d' % K, 'holds', 'least lam/lam0: %s' % ns(min0, 5), AG if min0 >= 1 else DF)
            cx = lam0 / (4 * mpmath.mpf(ent['kappa1']))
            add('D_cross%g_%d' % (L, K), '212', 'crossover lambda_min(1/2, L)/(4 kappa1) with our N=%d value %s' % (K, ns(lam0, 5)), CHAPTER_CROSSOVER[L], ns(cx, 4), AG if agrees(CHAPTER_CROSSOVER[L], cx) else DF)
            f2 = ent.get('crossover_factor2')
            add('D_x2_%g_%d' % (L, K), '212', 'damping c at which lambda_min(1/2 + c, L) is twice lambda_min(1/2, L), N=%d' % K, CHAPTER_CROSSOVER[L], f2['c'] if f2 else 'n/a', (AG if agrees(CHAPTER_CROSSOVER[L], f2['c']) else DF) if f2 else PEND,
                'ratio to the chapter formula %s' % f2['over_chapter_formula'] if f2 else '')
    res['claims'] = cl
    res.setdefault('meta', {})['claims_made'] = time.strftime('%Y-%m-%d %H:%M')


# ---------------------------------------------------------------------------------------------------------------------------------
def int_list(s):
    return [int(x) for x in s.split(',')]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Table zh:tab:hp, the zero side check and the damping ladder (see the header)")
    ap.add_argument('--part', default='all', help='table, prec, zeros, damping, claims or all')
    ap.add_argument('--quick', action='store_true', help='L = 0.5, N = 8 at 40 digits (writes to /tmp unless --out is given)')
    ap.add_argument('--out', default=None)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--L', default='0.5,0.8')
    ap.add_argument('--K', default=None, help='list of N (highest mode order); defaults depend on the part')
    ap.add_argument('--dps', type=int, default=60)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--zk', default='100,300,400,1700')
    args = ap.parse_args()
    Ls = [float(x) for x in args.L.split(',')]
    dps = args.dps
    zks = int_list(args.zk)
    if args.quick:
        Ls, dps, zks = [0.5], 40, [100]
        out = args.out or '/tmp/zh_window_eigs_quick.json'
    else:
        out = args.out or default_path('zh_window_eigs.json')
    log = Log(os.path.splitext(out)[0] + '.log')
    log('=== zh_window_eigs.py --part %s%s on %s at %s (L = %s, dps = %d)' % (args.part, ' --quick' if args.quick else '', socket.gethostname(), time.strftime('%Y-%m-%d %H:%M'), Ls, dps))
    res = load(out)
    res['conventions'] = dict(L='window half length a of ym', N='K = highest mode order, even sector size K + 1 (cos j = 0..K), odd sector size K (sin j = 1..K)', alpha='1/2 + c, c = y of the module',
                              basis='cossin', precision='dps %d (module working precision dps + %d)' % (dps, ym.GUARD))
    parts = ['table', 'prec', 'zeros', 'damping', 'claims'] if args.part == 'all' else [args.part]
    for part in parts:
        if part == 'table':
            part_table(res, log, Ls, [8] if args.quick else (int_list(args.K) if args.K else LADDER_N), dps, args.procs, out)
        elif part == 'prec':
            part_prec(res, log, Ls, [8] if args.quick else (int_list(args.K) if args.K else [26]), dps, args.procs, out)
        elif part == 'zeros':
            part_zeros(res, log, Ls, [8] if args.quick else (int_list(args.K) if args.K else [26, 40]), dps, find_zeros(args.zeros), zks, out)
        elif part == 'damping':
            part_damping(res, log, Ls, [8] if args.quick else (int_list(args.K) if args.K else [26, 40]), dps, args.procs, out, quick=args.quick)
        elif part == 'claims':
            part_claims(res, log)
        save(out, res)
    log('=== done')
