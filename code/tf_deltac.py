# tf_deltac.py: Table tf:tab:deltac of Chapter tf (Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction), Section tf:sec:window: the first order visibility distance
# delta_c = sqrt(2/S_0) at the actual zeros of zeta, in two configurations (a quadruple replacing an on line zero, a quadruple inserted at the midpoint of a gap) at six heights, and the statements of the paragraph after the table
# (the ratio delta_c/gap is stable, about 0.45 at midpoints and 0.66 to 0.73 when a zero is replaced; the least values 0.235 and 0.088 at height 74921; the medians fall below 1/2 from height 1419 at midpoints and from 74921 when replacing).
# Added for the book in October 2026, written from the mathematics (the collaborator's three_families.py, deltac and deltarep, was neither read nor run). The review's finding M3 says that both configurations change the number of zeros
# (a quadruple at a midpoint adds two zeros to the gap, one replacing a zero adds one, the book's perturbation, Theorem thm:ch13:perturb, replaces two zeros and keeps the count), so this script also computes the count preserving configuration.
# Definitions. S(theta) = -(log|Xi|)''(theta) = sum over all zeros (theta - theta_rho)^(-2), Xi(theta) = xi(1/2 + i theta); the formula in terms of zeta, zeta', zeta'' and the Gamma factor is in tf_xi.py (S = L2 + zeta''/zeta - (zeta'/zeta)^2 at s = 1/2 + i theta,
# the same as the review's -[(log(theta^2 + 1/4))'' - (1/4) Re psi'(1/4 + i theta/2) + (log|Z|)'']). The first order test at the ordinate g of a quadruple 1/2 +- delta +- i g sees it when S_0(g) - 2/delta^2 < 0, that is delta < delta_c = sqrt(2/S_0(g)),
# S_0 being the contribution of the other zeros (Proposition tf:prop:nearline). For each height n0 in {1, 10, 100, 1000, 10000, 100000}, with gamma_n the n-th ordinate:
#   midpoint: the 20 midpoints m_n = (gamma_n + gamma_{n+1})/2, n = n0 .. n0 + 19, S_0 = S(m_n), gap = gamma_{n+1} - gamma_n.
#   replacement: the 20 zeros gamma_n, n = n0 + 1 .. n0 + 20, S_0(gamma_n) = lim (S(theta) - (theta - gamma_n)^(-2)) in the closed form L2 + zeta'''/(3 zeta') - (zeta''/zeta')^2/4 at 1/2 + i gamma_n (equal to (Xi''^2/4 - Xi' Xi'''/3)/Xi'^2),
#   gap = (gamma_{n+1} - gamma_{n-1})/2 (so the zeros n0 .. n0 + 21 are needed).
#   count preserving (the review's M3): the same 20 midpoints with the two zeros next to the midpoint removed, S_0 = S(m_n) - 8/(gamma_{n+1} - gamma_n)^2 (each is at distance gap/2), gap = gamma_{n+1} - gamma_n; in the lattice S_0 = (pi^2 - 8)/Delta^2, delta_c = 1.034 Delta.
# delta_c = sqrt(2/S_0); per height the median (mean of the two middle values; the upper median is written too) and the least delta_c, and the median over the 20 points of delta_c/gap; the height is gamma_{n0} and the mean gap is the mean local gap of the replacement run.
# Second order variant (not in the chapter, written under the key exact_first_order): the quadruple's own images at -g contribute 2 Re (2g + i delta)^(-2) at theta = g and the mirror of a removed zero is removed too, so delta_c solves
# S_0x + 2 (4 g^2 - delta^2)/(4 g^2 + delta^2)^2 = 2/delta^2, with S_0x = S_0 - (2 gamma)^(-2) (replacement) and S_0x = S_0 - (m + gamma_n)^(-2) - (m + gamma_{n+1})^(-2) (count preserving); only the first row can feel it.
# Zeros: gamma_n for n <= 1700 from data/zeros_1700.txt (68 digits), the others by mpmath's zetazero (22 consecutive zeros at n0 = 10000 and 100000), at the working precision (30 digits).
# Error control. (a) S and S_0 recomputed at twice the precision at three points of the highest height (relative difference printed); (b) the closed form S_0 against its definition, (S(gamma + h) + S(gamma - h))/2 - h^(-2) with h = 1e-4, at three zeros (error O(h^2));
# (c) S from the zeta formula against a numerical second derivative of log|Xi| computed from xi itself (mpmath diff) at two midpoints; (d) mpmath's zetazero against the file at n = 1, 29, 1000 (and 1700 in the full run); (e) the imaginary part of S, which must vanish.
# The chapter's printed values are in the table CLAIM below, used only to print the comparison; nothing is tuned to them. Nothing here is a proof: mpmath is floating point at the stated precision.
# Output: data/tf_deltac.json (numbers as strings with 20 significant digits, of which about 25 are trusted at dps 30; the 20 points of every configuration at every height, the statements, the checks, the comparison with the chapter) and data/tf_deltac.log (the printed output).
# Usage: python3 tf_deltac.py [--quick] [--out FILE] [--workers W] [--zeros FILE] [--dps D] [--heights 1,10,...] [--count 20]; --quick does the four rows n0 = 1, 10, 100, 1000 (all zeros from the file, one process, about 4 s) and writes to /tmp unless --out is given;
# with --quick --heights N0 a single row is tested the same way (one process, 19 s at n0 = 10000 and 23 s at n0 = 100000, the zeros by zetazero at 0.5 s each).
# The full case (six rows, zetazero for n0 = 10000 and 100000, 8 processes) is for the Hopper cluster, run from the repository root:
#   bash misc/tools/hopper_run.sh -c 8 -t 120 -g tf_deltac.json code/tf_deltac.py code/tf_xi.py data/zeros_1700.txt -- --out tf_deltac.json | tee data/tf_deltac.log; mv tf_deltac.json data/
import argparse
import json
import os
import sys
import tempfile
import time
from multiprocessing import get_context

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp  # noqa: E402
import tf_xi as X  # noqa: E402

T0 = time.time()
HEIGHTS = [1, 10, 100, 1000, 10000, 100000]

# The chapter's Table tf:tab:deltac as printed (heights and mean gaps from the review's inventory), as (median, least, ratio) for
# the replacement and the midpoint blocks. Used only to print the comparison.
CLAIM = {
    1: dict(height='14.1', gap='3.177', rep=('2.064', '1.503', '0.716'), mid=('1.339', '0.692', '0.450')),
    10: dict(height='49.8', gap='2.558', rep=('1.825', '1.335', '0.719'), mid=('1.128', '0.589', '0.447')),
    100: dict(height='236.5', gap='1.679', rep=('1.201', '0.951', '0.733'), mid=('0.674', '0.447', '0.454')),
    1000: dict(height='1419.4', gap='1.178', rep=('0.799', '0.654', '0.711'), mid=('0.469', '0.297', '0.452')),
    10000: dict(height='9877.8', gap='0.848', rep=('0.579', '0.419', '0.701'), mid=('0.375', '0.165', '0.444')),
    100000: dict(height='74920.8', gap='0.661', rep=('0.446', '0.235', '0.662'), mid=('0.296', '0.088', '0.437')),
}


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


def init_worker(dps):
    mp.mp.dps = dps


def job_zero(n):
    return n, mp.nstr(mp.zetazero(n).imag, mp.mp.dps + 6)


def job_S(args):
    key, theta = args
    v = X.S_line(mp.mpf(theta))
    return key, mp.nstr(v.real, mp.mp.dps), mp.nstr(abs(v.imag), 3)


def job_S0(args):
    key, gamma = args
    v = X.S0_zero(mp.mpf(gamma))
    return key, mp.nstr(v.real, mp.mp.dps), mp.nstr(abs(v.imag), 3)


def run_pool(func, tasks, workers, dps):
    if workers <= 0:
        init_worker(dps)
        return [func(t) for t in tasks]
    with get_context('fork').Pool(workers, initializer=init_worker, initargs=(dps,)) as pool:
        return list(pool.imap_unordered(func, tasks, chunksize=1))


def delta_exact(S0x, g):
    """delta solving S0x + 2 (4 g^2 - d^2)/(4 g^2 + d^2)^2 = 2/d^2 (the quadruple's images at -g included), by bisection."""
    def f(d):
        return S0x + 2 * (4 * g * g - d * d) / (4 * g * g + d * d) ** 2 - 2 / (d * d)
    d0 = mp.sqrt(2 / S0x)
    lo, hi = d0 / 2, d0 * 2
    for _ in range(110):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def block(points):
    """Statistics of one configuration at one height; points are dicts with S0, delta_c, gap, delta_exact."""
    d = [p['delta_c'] for p in points]
    r = [p['delta_c'] / p['gap'] for p in points]
    de = [p['delta_exact'] for p in points]
    re_ = [p['delta_exact'] / p['gap'] for p in points]
    return dict(median=X.median(d), median_upper=X.upper_median(d), least=min(d), ratio_median=X.median(r),
                ratio_median_upper=X.upper_median(r), ratio_least=min(r), ratio_largest=max(r),
                exact_first_order=dict(median=X.median(de), least=min(de), ratio_median=X.median(re_)), points=points)


def jsonable(x):
    if isinstance(x, dict):
        return {k: jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, mp.mpf):
        return X.fmt(x, 20)
    return x


def main():
    ap = argparse.ArgumentParser(description='Table tf:tab:deltac: the first order visibility distance at the actual zeros (see the header of the file).')
    ap.add_argument('--quick', action='store_true', help='rows n0 = 1, 10, 100, 1000 only, one process, writes to /tmp unless --out')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--dps', type=int, default=30)
    ap.add_argument('--heights', default=None, help='comma separated n0 values')
    ap.add_argument('--count', type=int, default=20)
    a = ap.parse_args()
    heights = [int(x) for x in a.heights.split(',')] if a.heights else (HEIGHTS[:4] if a.quick else HEIGHTS)
    workers = 0 if a.quick else a.workers
    out = a.out or (os.path.join(tempfile.gettempdir(), 'tf_deltac_quick.json') if a.quick else 'tf_deltac.json')
    dps, cnt = a.dps, a.count
    mp.mp.dps = dps
    zfile = X.find_zeros_file(a.zeros)
    table = X.read_zero_table(zfile)
    log('tf_deltac: heights n0 = %s, %d points each, dps %d, %d process(es), zeros file %s (%d zeros)' % (heights, cnt, dps, max(workers, 1), zfile, len(table)))

    # Stage 1: the ordinates gamma_n, n = n0 .. n0 + cnt + 1.
    gam, need = {}, []
    for n0 in heights:
        for n in range(n0, n0 + cnt + 2):
            if n <= len(table):
                gam[n] = mp.mpf(table[n - 1])
            else:
                need.append(n)
    need = sorted(set(need))
    if need:
        log('zetazero for %d ordinates beyond the file (n = %d .. %d)' % (len(need), need[0], need[-1]))
        for n, s in run_pool(job_zero, need, workers, dps):
            gam[n] = mp.mpf(s)
    for n0 in heights:
        for n in range(n0, n0 + cnt + 1):
            assert gam[n + 1] > gam[n], 'zeros out of order at n = %d' % n
    log('ordinates ready; gamma_n0 = ' + ', '.join('%s' % X.fmt(gam[n0], 12) for n0 in heights))

    # Stage 2: S at the midpoints, S_0 at the zeros (largest heights first for the balance of the pool).
    tasks_S, tasks_S0 = [], []
    for n0 in sorted(heights, reverse=True):
        for i in range(cnt):
            n = n0 + i
            tasks_S.append(((n0, i), mp.nstr((gam[n] + gam[n + 1]) / 2, dps + 6)))
        for i in range(cnt):
            tasks_S0.append(((n0, i), mp.nstr(gam[n0 + 1 + i], dps + 6)))
    Sv, S0v, imag_max = {}, {}, mp.mpf(0)
    for key, v, im in run_pool(job_S, tasks_S, workers, dps):
        Sv[key] = mp.mpf(v)
        imag_max = max(imag_max, mp.mpf(im))
    log('S at %d midpoints done' % len(tasks_S))
    for key, v, im in run_pool(job_S0, tasks_S0, workers, dps):
        S0v[key] = mp.mpf(v)
        imag_max = max(imag_max, mp.mpf(im))
    log('S_0 at %d zeros done; largest imaginary part %s' % (len(tasks_S0), mp.nstr(imag_max, 3)))

    # Stage 3: the three configurations.
    rows = {}
    for n0 in heights:
        mid, rep, cp = [], [], []
        for i in range(cnt):
            n = n0 + i
            m = (gam[n] + gam[n + 1]) / 2
            gap = gam[n + 1] - gam[n]
            S = Sv[(n0, i)]
            mid.append(dict(n=n, ordinate=m, gap=gap, S0=S, delta_c=mp.sqrt(2 / S), delta_exact=delta_exact(S, m)))
            S0cp = S - 8 / gap ** 2
            S0cpx = S0cp - 1 / (m + gam[n]) ** 2 - 1 / (m + gam[n + 1]) ** 2
            cp.append(dict(n=n, ordinate=m, gap=gap, S0=S0cp, delta_c=mp.sqrt(2 / S0cp), delta_exact=delta_exact(S0cpx, m)))
        local = []
        for i in range(cnt):
            n = n0 + 1 + i
            g = gam[n]
            gap = (gam[n + 1] - gam[n - 1]) / 2
            local.append(gap)
            S0 = S0v[(n0, i)]
            rep.append(dict(n=n, ordinate=g, gap=gap, S0=S0, delta_c=mp.sqrt(2 / S0), delta_exact=delta_exact(S0 - 1 / (2 * g) ** 2, g)))
        rows[n0] = dict(n0=n0, height=gam[n0], mean_local_gap=sum(local) / cnt, mean_midpoint_gap=(gam[n0 + cnt] - gam[n0]) / cnt,
                        replacement=block(rep), midpoint=block(mid), count_preserving=block(cp))

    # Stage 4: error control.
    checks = {}
    top = max(heights)
    pick = [0] if a.quick else [0, cnt // 2, cnt - 1]
    hp = 2 * dps
    rel_S, rel_S0, rel_lim = mp.mpf(0), mp.mpf(0), mp.mpf(0)
    mp.mp.dps = hp
    for i in pick:
        n = top + i
        gh = [X.ordinate(k, table) if k <= len(table) else mp.zetazero(k).imag for k in (n, n + 1)]
        m = (gh[0] + gh[1]) / 2
        Sh = X.S_line(m).real
        rel_S = max(rel_S, abs(Sh - Sv[(top, i)]) / abs(Sh))
        n = top + 1 + i
        g = X.ordinate(n, table) if n <= len(table) else mp.zetazero(n).imag
        S0h = X.S0_zero(g).real
        rel_S0 = max(rel_S0, abs(S0h - S0v[(top, i)]) / abs(S0h))
        h = mp.mpf(10) ** (-4)
        lim = (X.S_line(g + h).real + X.S_line(g - h).real) / 2 - 1 / h ** 2
        rel_lim = max(rel_lim, abs(lim - S0h) / abs(S0h))
    mp.mp.dps = dps
    checks['precision_S_relative_difference_dps_%d_vs_%d' % (dps, hp)] = X.fmt(rel_S, 3)
    checks['precision_S0_relative_difference_dps_%d_vs_%d' % (dps, hp)] = X.fmt(rel_S0, 3)
    checks['S0_closed_form_vs_limit_h_1e-4_relative_difference'] = X.fmt(rel_lim, 3)
    checks['S0_check_zeros'] = 'n = %s of the highest height %d' % ([top + 1 + i for i in pick], top)
    rel_fd = mp.mpf(0)
    low = min(heights)
    for i in (0, cnt - 1):
        n = low + i
        m = (gam[n] + gam[n + 1]) / 2
        fd = -mp.diff(lambda t: mp.log(abs(X.xi_line(t))), m, 2)
        rel_fd = max(rel_fd, abs(fd - Sv[(low, i)]) / abs(fd))
    checks['S_vs_numerical_second_derivative_of_log_abs_Xi_relative_difference'] = X.fmt(rel_fd, 3)
    zdiff = mp.mpf(0)
    for n in ([1, 29, 1000] + ([] if a.quick else [1700])):
        zdiff = max(zdiff, abs(mp.zetazero(n).imag - mp.mpf(table[n - 1])))
    checks['zetazero_vs_file_largest_difference_n_1_29_1000%s' % ('' if a.quick else '_1700')] = X.fmt(zdiff, 3)
    checks['largest_imaginary_part_of_S'] = mp.nstr(imag_max, 3)
    for k, v in checks.items():
        log('check  %s = %s' % (k, v))

    # Stage 5: print the table and compare with the chapter.
    comp = []

    def compare(name, claim, mine):
        """agrees: the computed value rounds to the printed one; one unit: it is within one unit of the last printed digit (a tie at the next digit); differs: more."""
        dec = len(claim.split('.')[1]) if '.' in claim else 0
        unit = mp.mpf(10) ** (-dec)
        dif = abs(mine - mp.mpf(claim))
        status = 'agrees' if dif <= unit / 2 * (1 + mp.mpf(10) ** (-20)) else ('last digit differs by one unit' if dif < unit else 'differs')
        comp.append(dict(item=name, chapter=claim, computed=X.fmt(mine, 10), status=status, difference_in_units_of_last_printed_digit=X.fmt(dif / unit, 4)))
        return status == 'agrees'

    log('')
    log('%8s %9s %7s | replacing a zero: median least ratio | at a midpoint: median least ratio | count preserving: median least ratio' % ('n0', 'height', 'gap'))
    for n0 in heights:
        r = rows[n0]
        rp, md, cq = r['replacement'], r['midpoint'], r['count_preserving']
        log('%8d %9s %7s |                 %6s %6s %6s |                %6s %6s %6s |                  %6s %6s %6s' % (
            n0, mp.nstr(r['height'], 6), mp.nstr(r['mean_local_gap'], 4),
            mp.nstr(rp['median'], 4), mp.nstr(rp['least'], 4), mp.nstr(rp['ratio_median'], 3),
            mp.nstr(md['median'], 4), mp.nstr(md['least'], 4), mp.nstr(md['ratio_median'], 3),
            mp.nstr(cq['median'], 4), mp.nstr(cq['least'], 4), mp.nstr(cq['ratio_median'], 3)))
        c = CLAIM.get(n0)
        if c:
            compare('n0=%d height' % n0, c['height'], r['height'])
            compare('n0=%d mean local gap' % n0, c['gap'], r['mean_local_gap'])
            for tag, blk, cl in (('replacement', rp, c['rep']), ('midpoint', md, c['mid'])):
                compare('n0=%d %s median delta_c' % (n0, tag), cl[0], blk['median'])
                compare('n0=%d %s least delta_c' % (n0, tag), cl[1], blk['least'])
                compare('n0=%d %s median ratio' % (n0, tag), cl[2], blk['ratio_median'])
    log('')
    log('the upper median (sorted[10]) and the second order variant, for the record:')
    for n0 in heights:
        r = rows[n0]
        for tag in ('replacement', 'midpoint', 'count_preserving'):
            b = r[tag]
            log('  n0 %6d %-16s median %s (upper %s), ratio median %s (upper %s); exact_first_order median %s least %s ratio median %s' % (
                n0, tag, mp.nstr(b['median'], 5), mp.nstr(b['median_upper'], 5), mp.nstr(b['ratio_median'], 4), mp.nstr(b['ratio_median_upper'], 4),
                mp.nstr(b['exact_first_order']['median'], 5), mp.nstr(b['exact_first_order']['least'], 5), mp.nstr(b['exact_first_order']['ratio_median'], 4)))
    log('')
    stm = {}
    rr = [rows[n0]['replacement']['ratio_median'] for n0 in heights]
    rm = [rows[n0]['midpoint']['ratio_median'] for n0 in heights]
    rc = [rows[n0]['count_preserving']['ratio_median'] for n0 in heights]
    stm['replacement_ratio_medians_range'] = [X.fmt(min(rr), 6), X.fmt(max(rr), 6)]
    stm['midpoint_ratio_medians_range'] = [X.fmt(min(rm), 6), X.fmt(max(rm), 6)]
    stm['count_preserving_ratio_medians_range'] = [X.fmt(min(rc), 6), X.fmt(max(rc), 6)]
    stm['lattice_ratios'] = dict(replacement='sqrt(6)/pi = ' + X.fmt(mp.sqrt(6) / mp.pi, 6), midpoint='sqrt(2)/pi = ' + X.fmt(mp.sqrt(2) / mp.pi, 6),
                                 count_preserving='sqrt(2/(pi^2 - 8)) = ' + X.fmt(mp.sqrt(2 / (mp.pi ** 2 - 8)), 6))
    for tag in ('replacement', 'midpoint', 'count_preserving'):
        first = None
        for n0 in heights:
            if rows[n0][tag]['median'] < mp.mpf(1) / 2:
                first = n0
                break
        stm['first_height_median_below_one_half_' + tag] = None if first is None else dict(n0=first, height=X.fmt(rows[first]['height'], 8))
    log('statements of line 203: replacement ratio medians run from %s to %s (chapter: 0.66 to 0.73, lattice %s)' % (
        mp.nstr(min(rr), 4), mp.nstr(max(rr), 4), mp.nstr(mp.sqrt(6) / mp.pi, 4)))
    log('                        midpoint ratio medians run from %s to %s (chapter: about 0.45, lattice %s)' % (mp.nstr(min(rm), 4), mp.nstr(max(rm), 4), mp.nstr(mp.sqrt(2) / mp.pi, 4)))
    log('                        count preserving ratio medians run from %s to %s (lattice %s)' % (mp.nstr(min(rc), 4), mp.nstr(max(rc), 4), mp.nstr(mp.sqrt(2 / (mp.pi ** 2 - 8)), 4)))
    for tag in ('midpoint', 'replacement', 'count_preserving'):
        f = stm['first_height_median_below_one_half_' + tag]
        log('                        first height with median delta_c below 1/2, %s: %s (chapter: %s)' % (
            tag, 'none among the rows computed' if f is None else 'n0 = %d, height %s' % (f['n0'], f['height'][:8]),
            {'midpoint': '1419', 'replacement': '74921', 'count_preserving': 'not in the chapter'}[tag]))
    if 100000 in rows:
        compare('least replacement delta_c at height 74921 (line 203)', '0.235', rows[100000]['replacement']['least'])
        compare('least midpoint delta_c at height 74921 (line 203)', '0.088', rows[100000]['midpoint']['least'])
    nbad = [c for c in comp if c['status'] != 'agrees']
    log('comparison with the chapter: %d numbers compared, %d agree to the printed digits, %d do not' % (len(comp), len(comp) - len(nbad), len(nbad)))
    for c in nbad:
        log('   %s: %s: chapter %s, computed %s (%s units of the last printed digit)' % (c['status'].upper(), c['item'], c['chapter'], c['computed'], mp.nstr(mp.mpf(c['difference_in_units_of_last_printed_digit']), 3)))

    res = dict(script='tf_deltac.py', chapter='tf', table='tf:tab:deltac',
               parameters=dict(dps=dps, heights_n0=heights, points_per_configuration=cnt, workers=workers, quick=bool(a.quick),
                               zeros_file=os.path.basename(zfile), ordinates_by_zetazero=need,
                               definitions='midpoint: S_0 = S(m); replacement: S_0 = lim (S - (theta - gamma)^-2) in closed form; count_preserving: S_0 = S(m) - 8/gap^2; delta_c = sqrt(2/S_0)'),
               rows={str(n0): jsonable(rows[n0]) for n0 in heights}, statements=stm, checks=checks, comparison=comp,
               elapsed_seconds=round(time.time() - T0, 1))
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1)
    log('wrote ' + out)


if __name__ == '__main__':
    main()
