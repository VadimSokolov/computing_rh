# dc_ladder.py: the ladder lambda(a, y) of the chapter "Deconvolution as a Damped Weil Criterion" (label prefix dc:, Chapter 20 of the book): Table dc:tab:ladder (Table 1.1), Figure dc:fig:ladder (Figure 1.1) and the prose numbers of Section dc:sec:numerics (source lines 122, 143, 147). Group D2 of the number checks.
# Added for the book in October 2026, written from the mathematics (no script of the chapter's author was read or run). The matrices come from the shared module code/yo_matrix.py (group Y2; Weil's functional on e^(-y|x|) (h * h~), h supported in [-a, a], orthonormal bases), which is used as it is.
# Definition. lambda(a, y) = inf Q_y(h) / ||h||_2^2 over real h supported in [-a, a], Q_y(h) = W(e^(-y|x|) (h * h~)), W Weil's functional; it is the least eigenvalue of the matrix of the form in an orthonormal basis of L^2[-a, a], and every finite basis gives an upper bound. For y > 0 to the right of the zeros lambda(a, y) decreases in a to m(y) = 2 xi'/xi(1/2 + y) (the module's selftest and the ladder below confirm the normalisation: the values approach m(y) from above).
# Statements checked (chapter lines). L1 Table 1.1 (128 to 134): lambda(a, y) for y in {0.5, 0.375, 0.2, 0.1, 0.05} and a in {0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4}, in two bases (Legendre polynomials M = 24, 32, 48, 64, 96, 128, the chapter's method, and the Fourier modes e^(i pi n u/a), |n| <= N, N = 24 to 768), the converged value with its error, the change from M = 32, and whether every value exceeds m(y) and decreases in a (147, "never cross it"). L2 (122): the change from M = 24 to M = 32 is at most 3e-6, and the distance from M = 32 to the converged value. L3 (134, 136, 147, 192): the row y = 0, 0.03336 at a = 0.25 and 9.40e-7 at a = 0.5 with convergence in the size, and the claim that for a at least 0.75 the value is below the double precision floor of about 1e-13 (extended precision at a = 0.75 and 1, the order of the upper bound). L4 (143): "four and a half orders of magnitude between a = 0.25 and a = 0.5, as the prime 2 enters": log10 of the ratio and lambda at a = (1/2) log 2 and a = 0.35. L5 (147): 3.3e-2 and 9.4e-7. L6 Figure 1.1: fig/dc_ladder.pdf from these data, black and white, log scale, dashed lines for m(y), y = 0 for the first two a only. L7: the cross check with Table yo:tab:lam of Chapter 17 (alpha = 1/2 + y, y in {0, 0.05, 0.1, 0.25, 0.5}, a in {0.35, 0.8, 1.5, 2.5, 4}), read from data/yo_damped.json.
# Method. For y > 0 double precision is enough (the module's entries are accurate to about 1e-13, the least eigenvalue to about 1e-12; part check compares with 30 digits). The Legendre values are rigorous upper bounds but their convergence in M is irregular (the increments are not monotone, the minimiser is not smooth at the window ends); the Fourier values converge like N^(-p), p = 1.3 to 1.8, so the converged value is the power law (Richardson) extrapolation v_inf = v3 - (v2 - v3) / ((s3/s2)^p - 1) of the three largest Fourier sizes, p solved from the ratio of the two increments. Its uncertainty is the shift of that estimate when the finest size is dropped (the extrapolation from the previous three sizes), which is conservative (it falls about tenfold per size). The Legendre best value (M = 128) minus the converged value must be positive (every computed value is an upper bound) and is the independent check of the Fourier extrapolation. For y = 0 and a above 0.5 the values are below 1e-13 and need extended precision (mpmath at dps + 15 digits inside the module): part y0x runs the two bases at dps 60 to 80 and repeats the smallest size of each basis at dps + 20; the values there are upper bounds that depend strongly on the size (they are NOT converged), so only their order is reported.
# Error control. (a) part check: double against 30 digits for both bases (y > 0 and y = 0), and 1.5 times the quadrature nodes; (b) two bases, six sizes each, the shift of the extrapolation; (c) every printed value of Table 1.1 is compared digit by digit with the M = 32 value (the chapter's own method) and with the converged value; (d) m(y) = 2 xi'/xi(1/2 + y) at 40 digits by the module (mpmath zeta derivative). Nothing here is a proof: the eigenvalues are floating point, the extrapolation is a model of the convergence.
# Output: data/dc_ladder.json (numbers as strings with the digits trusted; keys check, grid[y=..][a=..], y0[a=..], L4, y0x[a=..], cross[y=..][a=..], summary, meta; a part written again replaces only its own entries, so partial runs add up) and data/dc_ladder.log (appended), fig/dc_ladder.pdf.
# Usage (repository root): timeout 60 python3 -I code/dc_ladder.py --part check|grid|y0|y0x|cross|fig|all [--quick] [--out FILE] [--fig FILE] [--procs P] [--y 0.5,0.2] [--a 0.25,1] [--legendre 24,32] [--fourier 24,48] [--xa 0.75,1] [--xleg 24] [--xfou 12,24] [--xdps 60]. --quick runs a small case of every part in under 30 s in one process and writes to /tmp unless --out is given. The grid of one value of y takes about 12 s locally (--part grid --y 0.5); the complete run, with the extended precision sizes up to M = 96, is meant for the Hopper cluster:
#   bash misc/tools/hopper_run.sh -c 16 -t 240 -g dc_ladder.json -g dc_ladder.log -g dc_ladder.pdf code/dc_ladder.py code/yo_matrix.py data/yo_damped.json -- --part all --procs 16 --fourier 24,48,96,192,384,768,1536 --xleg 24,32,48,64,96 --xfou 12,24,36,48,60 --out dc_ladder.json --fig dc_ladder.pdf
import argparse
import json
import math
import os
import platform
import socket
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import mpmath
import numpy as np

import yo_matrix as ym

mpmath.mp.dps = 50

LEG = [24, 32, 48, 64, 96, 128]
FOU = [24, 48, 96, 192, 384, 768]
A_TABLE = ['0.25', '0.5', '0.75', '1', '1.5', '2', '3', '4']
Y_TABLE = ['0.5', '0.375', '0.2', '0.1', '0.05']
A_CROSS = ['0.35', '0.8', '1.5', '2.5', '4']
Y_CROSS = ['0.05', '0.1', '0.25', '0.5']
LOG2_HALF = '0.34657359027997265'
A_Y0 = ['0.25', '0.3', LOG2_HALF, '0.35', '0.4', '0.45', '0.5']
A_Y0X = ['0.75', '0.8', '1']
DPS_Y0X = {'0.25': 30, '0.5': 30, '0.75': 60, '0.8': 60, '1': 80}
# Table 1.1 as printed in the chapter (lines 128 to 134): m(y) and the eight entries for a = 0.25 ... 4
PRINTED = {
    '0.5': ('0.04619', ['0.08949', '0.04953', '0.04743', '0.04684', '0.04646', '0.04634', '0.04626', '0.04623']),
    '0.375': ('0.03465', ['0.07564', '0.03726', '0.03560', '0.03514', '0.03486', '0.03476', '0.03470', '0.03468']),
    '0.2': ('0.01848', ['0.05605', '0.02001', '0.01901', '0.01875', '0.01859', '0.01854', '0.01851', '0.01850']),
    '0.1': ('0.00924', ['0.04474', '0.01009', '0.00951', '0.00938', '0.00930', '0.00927', '0.00926', '0.00925']),
    '0.05': ('0.00462', ['0.03906', '0.00508', '0.00476', '0.00469', '0.00465', '0.00464', '0.00463', '0.00462']),
}
PRINTED_Y0 = {'0.25': '0.03336', '0.5': '9.40e-7'}


class Log:
    def __init__(self, path):
        self.fh = open(path, 'a') if path else None

    def __call__(self, msg=''):
        print(msg, flush=True)
        if self.fh:
            self.fh.write(msg + '\n')
            self.fh.flush()


def default_path(stem):
    root = os.path.dirname(HERE)
    if os.path.basename(HERE) == 'code' and os.path.isdir(os.path.join(root, 'data')):
        return os.path.join(root, 'data', stem)
    return os.path.join(os.getcwd(), stem)


def fig_default():
    root = os.path.dirname(HERE)
    if os.path.basename(HERE) == 'code' and os.path.isdir(os.path.join(root, 'fig')):
        return os.path.join(root, 'fig', 'dc_ladder.pdf')
    return os.path.join(os.getcwd(), 'dc_ladder.pdf')


def load(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def save(path, res):
    with open(path, 'w') as fh:
        json.dump(res, fh, indent=1)


def stamp(res, part, seconds, mode, **extra):
    meta = res.setdefault('meta', {})
    meta['script'] = 'dc_ladder.py'
    meta[part] = dict(seconds=round(seconds, 1), host=socket.gethostname(), mode=mode, date=time.strftime('%Y-%m-%d %H:%M'),
                      python=platform.python_version(), numpy=np.__version__, mpmath=mpmath.__version__, **extra)


def fnum(v, n=14):
    return mpmath.nstr(mpmath.mpf(v), n)


def sh(s, n=7):
    """A string number shortened to n significant digits (keeps the exponent), '-' stays '-'."""
    return s if s in (None, '-') else mpmath.nstr(mpmath.mpf(s), n)


def parse_list(s, conv=str):
    return [conv(x) for x in s.split(',') if x.strip()] if s else None


def m_of(y):
    """m(y) = 2 xi'/xi(1/2 + y) at 40 digits (y a decimal string); 0 at y = 0."""
    if mpmath.mpf(y) == 0:
        return mpmath.mpf(0)
    with mpmath.workdps(60):
        al = mpmath.mpf('0.5') + mpmath.mpf(y)
        return 2 * ym.b_alpha(al, dps=40)


def eig(a, y, basis, size, dps, nodes=None):
    """The least eigenvalue as an mpf: double precision for dps <= 15, extended precision (mpmath, dps + 15 digits) above."""
    if dps <= 15:
        return mpmath.mpf(float(ym.least_eigenvalue(a, y=y, basis=basis, size=size, nodes=nodes)))
    return +mpmath.mpf(ym.least_eigenvalue(a, y=y, basis=basis, size=size, dps=dps, nodes=nodes))


def series_job(job):
    a, y, basis, size, dps = job
    t = time.time()
    v = eig(a, y, basis, size, dps)
    return job, mpmath.nstr(v, 32), time.time() - t


def run_jobs(func, jobs, procs):
    if procs <= 1 or len(jobs) < 2:
        return [func(j) for j in jobs]
    import multiprocessing
    with multiprocessing.get_context('fork').Pool(procs) as pool:
        return pool.map(func, jobs, chunksize=1)


def rich3(sizes, vals):
    """Power law extrapolation v(s) = v_inf + c s^(-p) from three sizes s1 < s2 < s3 with decreasing values: returns (v_inf, p) or None."""
    s1, s2, s3 = [float(s) for s in sizes]
    v1, v2, v3 = [mpmath.mpf(v) for v in vals]
    d12, d23 = v1 - v2, v2 - v3
    if not (d12 > 0 and d23 > 0):
        return None
    rho = float(d12 / d23)

    def g(p):
        return (s1 ** -p - s2 ** -p) / (s2 ** -p - s3 ** -p)

    lo, hi = 0.05, 12.0
    if not (g(lo) <= rho <= g(hi)):
        return None
    for _ in range(100):
        mid = (lo + hi) / 2
        if g(mid) < rho:
            lo = mid
        else:
            hi = mid
    p = (lo + hi) / 2
    return v3 - d23 / (mpmath.mpf(s3 / s2) ** p - 1), p


def estimate(leg, fou):
    """Converged value and its uncertainty from the Fourier sequence, with the Legendre sequence as the independent check."""
    out = dict(legendre={str(k): fnum(v) for k, v in sorted(leg.items())}, fourier={str(k): fnum(v) for k, v in sorted(fou.items())})
    ns = sorted(fou)
    conv = unc = None
    if len(ns) >= 3:
        t1 = ns[-3:]
        r1 = rich3(t1, [fou[n] for n in t1])
        r0 = rich3(ns[-4:-1], [fou[n] for n in ns[-4:-1]]) if len(ns) >= 4 else None
        if r1 is not None:
            conv = r1[0]
            out['fourier_rich'] = dict(sizes=t1, p='%.3f' % r1[1], estimate=fnum(r1[0]))
            if r0 is not None:
                unc = abs(r1[0] - r0[0])
                out['fourier_rich'].update(previous_sizes=ns[-4:-1], previous_p='%.3f' % r0[1], previous_estimate=fnum(r0[0]))
                out['uncertainty'] = fnum(unc, 3)
            out['converged'] = fnum(conv)
    ms = sorted(leg)
    if len(ms) >= 3:
        rl = rich3(ms[-3:], [leg[k] for k in ms[-3:]])
        if rl is not None:
            out['legendre_rich_last_triple'] = dict(sizes=ms[-3:], p='%.3f' % rl[1], estimate=fnum(rl[0]), note='informational: the Legendre increments are not monotone')
    if ms:
        best = min(leg.values())
        out['legendre_best'] = fnum(best)
        if conv is not None:
            out['legendre_best_minus_converged'] = fnum(best - conv, 4)
    if 24 in leg and 32 in leg:
        out['change_24_to_32'] = fnum(leg[24] - leg[32], 4)
    if 32 in leg and conv is not None:
        out['M32_minus_converged'] = fnum(leg[32] - conv, 4)
    return out, conv, unc


def collect(out):
    vals = {}
    for (a, y, basis, size, dps), s, sec in out:
        vals.setdefault((y, a), {'legendre': {}, 'fourier': {}})[basis][size] = mpmath.mpf(s)
    return vals


def part_check(res, log):
    t0 = time.time()
    chk = {}
    for a, y, basis, size in (('1', '0.2', 'fourier', 24), ('1', '0.2', 'legendre', 16), ('0.25', '0', 'fourier', 12), ('0.5', '0', 'legendre', 12)):
        vd, ve = eig(a, y, basis, size, 15), eig(a, y, basis, size, 30)
        key = '%s a=%s y=%s size=%d' % (basis, a, y, size)
        chk[key] = dict(double=fnum(vd), digits30=fnum(ve, 20), difference=fnum(abs(vd - ve), 3))
        log('check double against 30 digits: %-38s %s %s diff %s' % (key, fnum(vd, 13), fnum(ve, 13), fnum(abs(vd - ve), 3)))
    for basis, size in (('fourier', 24), ('legendre', 24)):
        try:
            n0 = ym._nq_default(basis, 1.0, size, 15)
        except Exception:
            n0 = 200
        v0, v1 = eig('1', '0.2', basis, size, 15), eig('1', '0.2', basis, size, 15, nodes=int(1.5 * n0))
        chk['nodes %s a=1 y=0.2 size=%d' % (basis, size)] = dict(default_nodes=n0, change_with_1_5_times_nodes=fnum(abs(v0 - v1), 3))
        log('check 1.5 times the nodes: %s size %d (%d nodes) change %s' % (basis, size, n0, fnum(abs(v0 - v1), 3)))
    res['check'] = chk
    stamp(res, 'check', time.time() - t0, 'double precision and 30 digits')


def make_table_lines(y, row, log):
    log('y = %s, m(y) = %s (printed %s)' % (y, fnum(m_of(y), 8), PRINTED[y][0]))
    log('   a    printed   M = 32 (chapter method)   converged (+- shift)       M32 - conv   M24 - M32   best Legendre - conv   p')
    for a in A_TABLE:
        e = row.get('a=' + a)
        if not e:
            continue
        i = A_TABLE.index(a)
        log('%5s  %s   %s            %s (%s)   %s   %s   %s   %s' % (a, PRINTED[y][1][i], sh(e['legendre'].get('32', '-'), 9), sh(e.get('converged', '-'), 10), e.get('uncertainty', '-'),
                                                           e.get('M32_minus_converged', '-'), e.get('change_24_to_32', '-'), e.get('legendre_best_minus_converged', '-'),
                                                           e.get('fourier_rich', {}).get('p', '-')))


def part_grid(res, log, ys, avals, legs, fous, procs):
    t0 = time.time()
    jobs = [(a, y, 'legendre', M, 15) for y in ys for a in avals for M in legs] + [(a, y, 'fourier', N, 15) for y in ys for a in avals for N in fous]
    vals = collect(run_jobs(series_job, jobs, procs))
    grid = res.setdefault('grid', {})
    for y in ys:
        row = grid.setdefault('y=' + y, {})
        m = m_of(y)
        for a in avals:
            v = vals[(y, a)]
            e, conv, unc = estimate(v['legendre'], v['fourier'])
            e['m'] = fnum(m)
            allv = list(v['legendre'].values()) + list(v['fourier'].values())
            e['all_computed_above_m'] = bool(min(allv) > m)
            if conv is not None:
                e['converged_above_m'] = bool(conv > m)
                e['excess_over_m'] = fnum(conv - m, 10)
            if a in A_TABLE and y in PRINTED:
                pr = PRINTED[y][1][A_TABLE.index(a)]
                e['printed'] = pr
                if 32 in v['legendre']:
                    e['M32_agrees_with_printed'] = bool('%.5f' % v['legendre'][32] == pr)
                if conv is not None:
                    e['converged_agrees_with_printed'] = bool('%.5f' % conv == pr)
            row['a=' + a] = e
        make_table_lines(y, row, log)
    summarise(res, log)
    stamp(res, 'grid', time.time() - t0, 'double precision', legendre=legs, fourier=fous)


def summarise(res, log):
    g = res.get('grid', {})
    s = dict(entries=0)
    mx24, mx32, mxlb = (0, None), (0, None), (0, None)
    all_above, mono, p32, pcv, bad32, badcv, minrich = True, {}, 0, 0, [], [], None
    for y, row in g.items():
        conv_by_a = []
        for a, e in row.items():
            s['entries'] += 1
            key = (y, a)
            if 'change_24_to_32' in e and abs(float(e['change_24_to_32'])) > mx24[0]:
                mx24 = (abs(float(e['change_24_to_32'])), key)
            if 'M32_minus_converged' in e and float(e['M32_minus_converged']) > mx32[0]:
                mx32 = (float(e['M32_minus_converged']), key)
            if 'legendre_best_minus_converged' in e:
                v = float(e['legendre_best_minus_converged'])
                if v > mxlb[0]:
                    mxlb = (v, key)
                if minrich is None or v < minrich[0]:
                    minrich = (v, key)
            all_above = all_above and e.get('all_computed_above_m', True) and e.get('converged_above_m', True)
            if 'M32_agrees_with_printed' in e:
                if e['M32_agrees_with_printed']:
                    p32 += 1
                else:
                    bad32.append('%s %s' % key)
            if 'converged_agrees_with_printed' in e:
                if e['converged_agrees_with_printed']:
                    pcv += 1
                else:
                    badcv.append('%s %s' % key)
            if 'converged' in e:
                conv_by_a.append((float(a[2:]), float(e['converged'])))
        conv_by_a.sort()
        mono[y] = all(c2 < c1 for (_, c1), (_, c2) in zip(conv_by_a, conv_by_a[1:]))
    s.update(max_abs_change_24_to_32=['%.3e' % mx24[0], mx24[1]], max_M32_minus_converged=['%.3e' % mx32[0], mx32[1]],
             max_legendre_best_minus_converged=['%.3e' % mxlb[0], mxlb[1]], min_legendre_best_minus_converged=(['%.3e' % minrich[0], minrich[1]] if minrich else None),
             every_value_exceeds_m=bool(all_above), converged_decreasing_in_a=mono, M32_entries_agreeing_with_printed=p32, M32_entries_disagreeing=bad32,
             converged_entries_agreeing_with_printed=pcv, converged_entries_disagreeing=badcv)
    res['summary'] = s
    log('summary of the grid (%d entries): max |M24 - M32| %s at %s; max (M32 - converged) %s at %s; best Legendre minus converged between %s and %s' % (
        s['entries'], s['max_abs_change_24_to_32'][0], s['max_abs_change_24_to_32'][1], s['max_M32_minus_converged'][0], s['max_M32_minus_converged'][1],
        s['min_legendre_best_minus_converged'][0] if minrich else '-', s['max_legendre_best_minus_converged'][0]))
    log('every computed and converged value exceeds m(y): %s; converged values decreasing in a: %s; M = 32 agrees with the printed entry in %d, converged agrees in %d (disagree M32: %s; converged: %s)' % (
        all_above, mono, p32, pcv, bad32, badcv))


def part_y0(res, log, avals, legs, fous, procs):
    t0 = time.time()
    jobs = [(a, '0', 'legendre', M, 15) for a in avals for M in legs] + [(a, '0', 'fourier', N, 15) for a in avals for N in fous]
    vals = collect(run_jobs(series_job, jobs, procs))
    y0 = res.setdefault('y0', {})
    for a in avals:
        v = vals[('0', a)]
        e, conv, unc = estimate(v['legendre'], v['fourier'])
        if a in PRINTED_Y0:
            e['printed'] = PRINTED_Y0[a]
        y0['a=' + a] = e
        log('y = 0, a = %-20s M = 24: %s  M = 32: %s  best Legendre: %s  converged: %s (+- %s) p = %s' % (
            a, sh(e['legendre'].get('24', '-')), sh(e['legendre'].get('32', '-')), sh(e.get('legendre_best', '-')), sh(e.get('converged', '-'), 9),
            e.get('uncertainty', '-'), e.get('fourier_rich', {}).get('p', '-')))
    derive_L4(res, log)
    stamp(res, 'y0', time.time() - t0, 'double precision', legendre=legs, fourier=fous)


def value_y0(res, a):
    e = res.get('y0', {}).get('a=' + a)
    if not e:
        return None
    return mpmath.mpf(e['converged']) if 'converged' in e else mpmath.mpf(e['legendre_best'])


def derive_L4(res, log):
    v = {a: value_y0(res, a) for a in ('0.25', '0.5', LOG2_HALF, '0.35')}
    if v['0.25'] is None or v['0.5'] is None:
        return
    out = dict(log10_ratio_025_to_05=mpmath.nstr(mpmath.log10(v['0.25'] / v['0.5']), 6))
    if v[LOG2_HALF] is not None:
        out['lambda_at_half_log2'] = fnum(v[LOG2_HALF], 8)
        out['orders_before_prime_2_enters_025_to_half_log2'] = mpmath.nstr(mpmath.log10(v['0.25'] / v[LOG2_HALF]), 5)
        out['orders_after_half_log2_to_05'] = mpmath.nstr(mpmath.log10(v[LOG2_HALF] / v['0.5']), 5)
    if v['0.35'] is not None:
        out['lambda_at_035'] = fnum(v['0.35'], 8)
        out['chapter17_printed_at_035'] = '1.20e-3'
    series_ = []
    for a in A_Y0:
        x = value_y0(res, a)
        if x is not None:
            series_.append('%s:%s' % (a[:6], mpmath.nstr(x, 6)))
    out['lambda_y0_by_a'] = series_
    res['L4'] = out
    log('L4: log10(lambda(0.25, 0) / lambda(0.5, 0)) = %s; lambda((1/2) log 2) = %s, lambda(0.35) = %s; orders before the prime 2 enters %s, after %s' % (
        out['log10_ratio_025_to_05'], out.get('lambda_at_half_log2'), out.get('lambda_at_035'), out.get('orders_before_prime_2_enters_025_to_half_log2'),
        out.get('orders_after_half_log2_to_05')))


def part_y0x(res, log, alist, legs, fous, xdps, procs):
    t0 = time.time()
    store = res.setdefault('y0x', {})
    for a in alist:
        dps = xdps or DPS_Y0X.get(a, 60)
        jobs = [(a, '0', 'legendre', M, dps) for M in legs] + [(a, '0', 'fourier', N, dps) for N in fous]
        chk = []
        if legs:
            chk.append((a, '0', 'legendre', min(legs), dps + 20))
        if fous:
            chk.append((a, '0', 'fourier', min(fous), dps + 20))
        out = run_jobs(series_job, jobs + chk, procs)
        main = {(j[2], j[3]): (mpmath.mpf(s), sec) for (j, s, sec) in out[:len(jobs)]}
        second = {(j[2], j[3]): mpmath.mpf(s) for (j, s, sec) in out[len(jobs):]}
        e = dict(dps=dps, check_dps=dps + 20, legendre={}, fourier={}, seconds={})
        agree = {}
        for (basis, size), (v, sec) in sorted(main.items()):
            e[basis][str(size)] = fnum(v, 18)
            e['seconds']['%s %d' % (basis, size)] = round(sec, 1)
            if (basis, size) in second:
                d = abs(v - second[(basis, size)]) / abs(v)
                agree['%s %d' % (basis, size)] = mpmath.nstr(d, 3)
        e['relative_difference_dps_vs_dps_plus_20'] = agree
        old = store.get('a=' + a)
        if old and old.get('dps') == dps:
            # a second invocation for the same a and dps adds its sizes to the stored ones
            for key in ('legendre', 'fourier', 'seconds', 'relative_difference_dps_vs_dps_plus_20'):
                merged = dict(old.get(key, {}))
                merged.update(e[key])
                e[key] = merged
        allv = [mpmath.mpf(v) for v in list(e['legendre'].values()) + list(e['fourier'].values())]
        if allv:
            best = min(allv)
            e['best_upper_bound'] = fnum(best, 8)
            e['log10_best_upper_bound'] = mpmath.nstr(mpmath.log10(best), 5)
        store['a=' + a] = e
        log('y = 0, a = %s, dps %d (check %d): Legendre %s; Fourier %s; best upper bound %s; relative differences between the two precisions %s' % (
            a, dps, dps + 20, {k: sh(v, 6) for k, v in e['legendre'].items()}, {k: sh(v, 6) for k, v in e['fourier'].items()}, e.get('best_upper_bound'), agree))
    stamp(res, 'y0x', time.time() - t0, 'extended precision', legendre=legs, fourier=fous)


def find_yo(path):
    cand = [path] if path else []
    cand += [os.path.join(os.path.dirname(HERE), 'data', 'yo_damped.json'), os.path.join(HERE, 'yo_damped.json'), os.path.join(os.getcwd(), 'yo_damped.json')]
    for c in cand:
        if c and os.path.exists(c):
            with open(c) as fh:
                return json.load(fh), c
    return None, None


def lookup(d, keys):
    for k in keys:
        if isinstance(d, dict) and k in d:
            return d[k]
    return None


def part_cross(res, log, ys, avals, legs, fous, procs, yo_path):
    t0 = time.time()
    yo, src = find_yo(yo_path)
    jobs = [(a, y, 'legendre', M, 15) for y in ys for a in avals for M in legs] + [(a, y, 'fourier', N, 15) for y in ys for a in avals for N in fous]
    vals = collect(run_jobs(series_job, jobs, procs))
    cross = res.setdefault('cross', {})
    worst = (0, None)
    worst_m32 = (0, None)
    for y in ys:
        row = cross.setdefault('y=' + y, {})
        alpha = mpmath.mpf('0.5') + mpmath.mpf(y)
        akeys = [mpmath.nstr(alpha, 3), '%g' % float(alpha), '%.1f' % float(alpha), '%.2f' % float(alpha)]
        for a in avals:
            v = vals[(y, a)]
            e, conv, unc = estimate(v['legendre'], v['fourier'])
            e['m'] = fnum(m_of(y))
            if yo is not None:
                ent = lookup(lookup(yo.get('table_comparison', {}), [a, '%g' % float(a), '%.1f' % float(a)]), akeys)
                if ent:
                    e['chapter17_printed'] = ent.get('printed')
                    e['chapter17_finite_section'] = ent.get('lam_N_script_big')
                    e['chapter17_converged'] = ent.get('converged')
                    if ent.get('converged') and conv is not None:
                        d = conv - mpmath.mpf(ent['converged'])
                        e['converged_minus_chapter17_converged'] = fnum(d, 3)
                        if abs(d) > worst[0]:
                            worst = (abs(d), (y, a))
                    if ent.get('lam_N_script_big') and 32 in v['legendre']:
                        d32 = v['legendre'][32] - mpmath.mpf(ent['lam_N_script_big'])
                        e['M32_minus_chapter17_finite_section'] = fnum(d32, 3)
                        if abs(d32) > worst_m32[0]:
                            worst_m32 = (abs(d32), (y, a))
            row['a=' + a] = e
            log('cross y = %-4s a = %-4s M32 %s converged %s (+- %s) chapter 17: finite section %s converged %s printed %s' % (
                y, a, sh(e['legendre'].get('32', '-'), 9), sh(e.get('converged', '-'), 10), e.get('uncertainty', '-'), e.get('chapter17_finite_section', '-'), e.get('chapter17_converged', '-'),
                e.get('chapter17_printed', '-')))
    res.setdefault('cross_summary', {})['source'] = os.path.relpath(src, os.path.dirname(HERE))
    summarise_cross(res, log)
    stamp(res, 'cross', time.time() - t0, 'double precision', legendre=legs, fourier=fous)


def summarise_cross(res, log):
    """Largest differences over every stored entry of the cross check (all invocations)."""
    worst, worst_m32, n = (0.0, None), (0.0, None), 0
    for y, row in res.get('cross', {}).items():
        for a, e in row.items():
            if 'converged_minus_chapter17_converged' in e:
                n += 1
                d = abs(float(e['converged_minus_chapter17_converged']))
                if d > worst[0]:
                    worst = (d, (y, a))
            if 'M32_minus_chapter17_finite_section' in e:
                d = abs(float(e['M32_minus_chapter17_finite_section']))
                if d > worst_m32[0]:
                    worst_m32 = (d, (y, a))
    res.setdefault('cross_summary', {}).update(entries_with_chapter17_converged=n, max_abs_converged_minus_chapter17_converged=['%.3e' % worst[0], worst[1]],
                                                max_abs_M32_minus_chapter17_finite_section=['%.3e' % worst_m32[0], worst_m32[1]])
    log('cross check (%d entries): largest |converged here - converged of Chapter 17| %.3e at %s; largest |M32 - Chapter 17 finite section| %.3e at %s' % (n, worst[0], worst[1], worst_m32[0], worst_m32[1]))


def part_fig(res, log, path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        log('matplotlib is not installed here: move the JSON to data/ and run  timeout 40 python3 -I code/dc_ladder.py --part fig  on the laptop')
        return
    g = res.get('grid', {})
    style = {'0.5': ('o', '0.0'), '0.375': ('s', '0.3'), '0.2': ('^', '0.15'), '0.1': ('v', '0.45'), '0.05': ('D', '0.6')}
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for y in Y_TABLE:
        row = g.get('y=' + y)
        if not row:
            continue
        pts = sorted((float(a[2:]), float(e['converged'])) for a, e in row.items() if 'converged' in e)
        if not pts:
            continue
        mk, col = style[y]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, ms=4.5, lw=1.1, color=col, mfc='white' if y in ('0.375', '0.1') else col, label='y = %s' % y)
        mval = float(mpmath.mpf(next(iter(row.values()))['m']))
        ax.axhline(mval, ls='--', lw=0.7, color=col)
        ax.text(4.55, mval, 'm = %.5f' % mval, va=('top' if y == '0.375' else 'bottom'), ha='right', fontsize=6.5, color='0.2')
    y0 = res.get('y0', {})
    pts = [(float(a), float(mpmath.mpf(y0['a=' + a]['converged']))) for a in ('0.25', '0.5') if 'a=' + a in y0 and 'converged' in y0['a=' + a]]
    if pts:
        ax.plot([p[0] for p in pts], [p[1] for p in pts], marker='*', ms=9, lw=1.6, color='black', label='y = 0 (a = 0.25, 0.5)')
    ax.set_yscale('log')
    ax.set_xlim(0, 4.6)
    ax.set_xlabel('a (half length of the window)')
    ax.set_ylabel(r'$\lambda(a,y)$, converged')
    ax.grid(True, which='major', lw=0.3, color='0.8')
    ax.legend(fontsize=7, loc='lower right', frameon=True, framealpha=1.0, edgecolor='0.6')
    fig.tight_layout()
    fig.savefig(path)
    log('figure written to %s' % path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', default='all', help='check, grid, y0, y0x, cross, crosssum, fig or all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--y', default=None)
    ap.add_argument('--a', default=None)
    ap.add_argument('--legendre', default=None)
    ap.add_argument('--fourier', default=None)
    ap.add_argument('--xa', default=None)
    ap.add_argument('--xleg', default=None)
    ap.add_argument('--xfou', default=None)
    ap.add_argument('--xdps', type=int, default=None)
    ap.add_argument('--yo', default=None, help='path of data/yo_damped.json for the cross check')
    args = ap.parse_args()
    t_all = time.time()
    if args.quick:
        out = args.out or '/tmp/dc_ladder_quick.json'
        figp = args.fig or '/tmp/dc_ladder_quick.pdf'
        legs, fous, ys, avals = [16, 24, 32], [12, 24, 48, 96], ['0.5', '0.1'], ['0.25', '1']
        y0a, xa, xleg, xfou = ['0.25', '0.5'], ['0.75'], [], [8, 12]
        cy, ca = ['0.5'], ['0.8']
    else:
        out = args.out or default_path('dc_ladder.json')
        figp = args.fig or fig_default()
        legs, fous, ys, avals = LEG, FOU, Y_TABLE, A_TABLE
        y0a, xa, xleg, xfou = A_Y0, A_Y0X, [24, 32], [12, 24, 36, 48]
        cy, ca = Y_CROSS, A_CROSS
    legs = parse_list(args.legendre, int) or legs
    fous = parse_list(args.fourier, int) or fous
    ys_g = parse_list(args.y) or ys
    avals_g = parse_list(args.a) or avals
    xa = parse_list(args.xa) or xa
    xleg = (parse_list(args.xleg, int) or []) if args.xleg is not None else xleg
    xfou = (parse_list(args.xfou, int) or []) if args.xfou is not None else xfou
    log = Log(os.path.splitext(out)[0] + '.log')
    log('=== dc_ladder.py %s part %s, %s ===' % (time.strftime('%Y-%m-%d %H:%M'), args.part, 'quick' if args.quick else 'full'))
    res = load(out)
    parts = ['check', 'grid', 'y0', 'y0x', 'cross', 'fig'] if args.part == 'all' else [args.part]
    for part in parts:
        if part == 'check':
            part_check(res, log)
        elif part == 'grid':
            part_grid(res, log, ys_g, avals_g, legs, fous, args.procs)
        elif part == 'y0':
            part_y0(res, log, parse_list(args.a) or y0a, legs, fous, args.procs)
        elif part == 'y0x':
            part_y0x(res, log, xa, xleg, xfou, args.xdps, args.procs)
        elif part == 'cross':
            part_cross(res, log, parse_list(args.y) or cy, parse_list(args.a) or ca, legs, fous, args.procs, args.yo)
        elif part == 'crosssum':
            summarise_cross(res, log)
        elif part == 'fig':
            part_fig(res, log, figp)
        else:
            raise SystemExit('unknown part ' + part)
        save(out, res)
    log('done in %.1f s, written to %s' % (time.time() - t_all, out))


if __name__ == '__main__':
    main()
