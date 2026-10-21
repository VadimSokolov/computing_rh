# tf_window.py: Table tf:tab:window and Figure tf:fig:window of Chapter tf (Three Families of Generalised Gamma Convolutions and the Deconvolution Obstruction), Section tf:sec:window: what a single basepoint sees of a hypothetical quadruple of
# zeros 1/2 +- delta +- i gamma* inserted at the midpoint gamma* of a gap, at heights 100 (the gap between the 29th and 30th zeros) and 9878 (between the 10000th and 10001st).
# Added for the book in October 2026, written from the mathematics (the collaborator's three_families.py and fig_window.py were neither read nor run).
# Checked: gamma* = 100.0745226 and 9878.218713, gaps 2.48666 and 0.872118, delta_c = sqrt(2/S(gamma*)) = 1.10904 and 0.381193; the levels y = delta k/40 (k = 1 .. 39, and 0.999 delta) at which min over the gap of v_y is negative, for delta in {0.05, 0.2, 0.3, 0.45}
# (chapter: every level at height 100 for every delta, at height 9878 every level for delta <= 0.3 and 0.191 <= y < 0.45, 24 of 40 levels, for delta = 0.45); the first order test S(gamma*) - 2/delta^2 = 3.8873 at height 9878, delta = 0.45; the sign change
# between y = 0.18 and y = 0.19125; the data of the figure (v_y on 121 points of the gap at y = 0.05, 0.15, 0.25, 0.40 for delta = 0.45 at height 9878, minima 0.0575, 0.0706, -0.2602, -4.8939); and the figure itself (--plot), redrawn from this script's own data.
# Reading of the chapter. v_y(theta) is the Thorin density of the basepoint 1/2 + y, v_y(theta) = (1/pi) Re xi'/xi(1/2 + y + i theta) + (1/pi) sum over b = +-delta and g = +-gamma* of (y - b)/((y - b)^2 + (theta - g)^2), the second sum being the four Cauchy
# terms Re 1/(s - rho) of the zeros rho = 1/2 + b + i g of the quadruple, added exactly (the review's formula; it agrees with the chapter, which gives at theta = gamma* the contribution (1/pi)((y - delta)^(-1) + (y + delta)^(-1)) of the two zeros at height gamma*; the mirror
# pair at -gamma* is included here, it is of size 1e-8 at height 9878 and larger at height 100). There is no difference of reading. xi'/xi(s) = 1/s + 1/(s - 1) - (1/2) log pi + (1/2) psi(s/2) + zeta'/zeta(s), by mpmath at 30 digits.
# Method. For each gap (n = 29, 10000), each delta, the 40 levels (and for delta = 0.45 at height 9878 the extra levels 0.18, 0.1895, 0.19), v_y is evaluated at the 41 points theta = gamma* + (j/40) gap, j = -20 .. 20 (the ends are the zeros themselves; the midpoint j = 0 is on the grid);
# the grid minimum is then refined by a golden section search on the two neighbouring grid intervals (25 steps), so that "negative somewhere in the gap" does not depend on the grid; both are written (the chapter's claim is the grid one). For the level at which the sign changes
# the critical level y_c (min over the gap equal to 0) is found by bisection, and so is the critical level at the ordinate theta = gamma* alone (the visibility curve of Proposition tf:prop:curve, the chapter's 0.1895). The first order test is also evaluated over the gap
# (S(theta) plus the quadruple's contribution 2 Re (theta - g + i delta)^(-2) summed over g = +-gamma*, at j = -19 .. 19; S has a pole at the zeros j = +-20), because "the first order test included" in the chapter needs positivity of S over the gap, not at gamma* alone.
# Error control. (a) three values of v recomputed at 15 more digits (relative difference printed); (b) the limit y -> 0 of Proposition tf:prop:nearline, v without the quadruple against (y/pi) S(theta) at y = 1e-3 (relative difference O(y^2)); (c) the Cauchy terms of the quadruple at the ordinate
# against the closed form (1/pi)(2y/(y^2 - delta^2) + the two terms at distance 2 gamma*) of the chapter's proof of Proposition tf:prop:curve; (d) the imaginary part of S (must vanish). The grid is the chapter's, the refinement and the bisection are additions.
# The chapter's printed values are in the tables CLAIM below, used only to print the comparison. Nothing here is a proof: the minima are over a grid with a local refinement, not certified (the chapter's "stays positive across the gap" would need an Arb enclosure and a derivative bound).
# Output: data/tf_window.json (numbers as strings; gaps, S, delta_c, first order margins, the level records of all (gap, delta), their summaries, critical levels, the figure data, checks, comparison) and data/tf_window.log; fig/tf_window.pdf by --plot.
# Usage: python3 tf_window.py [--quick] [--out FILE] [--workers W] [--zeros FILE] [--dps D] [--gaps N,..] [--deltas D,..] [--ks K,..] [--jstep J] [--no-figure]; --quick uses the gap n = 29 only, deltas 0.45 and 1.5 (1.5 is a test of the bisection, not a case of the chapter), 8 levels,
# 11 points, a 31 point figure, one process (about 3 s) and writes to /tmp unless --out is given; the other options restrict a run further, for instance --quick --gaps 10000 --deltas 0.45 --ks 17 --jstep 4 is a one process test of the gap of height 9878 (25 s);
# python3 tf_window.py --plot [--data FILE] [--fig FILE] [--root DIR] draws fig/tf_window.pdf from data/tf_window.json (about 2 s; black and white: line styles and markers, a symmetric logarithmic vertical scale linear within +-0.1).
# The full case (two gaps, 8 processes, a few minutes) is for the Hopper cluster, run from the repository root:
#   bash misc/tools/hopper_run.sh -c 8 -t 120 -g tf_window.json code/tf_window.py code/tf_xi.py data/zeros_1700.txt -- --out tf_window.json | tee data/tf_window.log; mv tf_window.json data/; python3 code/tf_window.py --plot
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
GAPS = [29, 10000]
DELTAS = ['0.05', '0.2', '0.3', '0.45']
FIG_LEVELS = ['0.05', '0.15', '0.25', '0.40']
EXTRA = {(10000, '0.45'): ['0.18', '0.1895', '0.19']}

# The chapter's printed values (and the collaborator's data quoted in the review), used only for the comparison.
CLAIM = {29: dict(gstar='100.0745226', gap='2.48666', delta_c='1.10904'), 10000: dict(gstar='9878.218713', gap='0.872118', delta_c='0.381193')}
CLAIM_FIRST_ORDER = '3.8873'
CLAIM_FIG_MIN = {'0.05': '0.0575', '0.15': '0.0706', '0.25': '-0.2602', '0.40': '-4.8939'}
CLAIM_SIGN = {'0.18': '0.0227', '0.19125': '-0.0046'}


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


def init_worker(dps):
    mp.mp.dps = dps


def nst(x, extra=4):
    return mp.nstr(x, mp.mp.dps + extra)


def v_base(y, theta):
    """(1/pi) Re xi'/xi(1/2 + y + i theta): the Thorin density of the actual zeros at level y."""
    return X.xi_logder(mp.mpc(X.half() + y, theta)).real / mp.pi


def v_quad(y, theta, delta, gstar):
    """The four Cauchy terms of the zeros 1/2 +- delta +- i gstar, divided by pi."""
    tot = mp.mpf(0)
    for b in (delta, -delta):
        for g in (gstar, -gstar):
            tot += (y - b) / ((y - b) ** 2 + (theta - g) ** 2)
    return tot / mp.pi


def v_total(y, theta, delta, gstar):
    return v_base(y, theta) + v_quad(y, theta, delta, gstar)


def S_quad(theta, delta, gstar):
    """The quadruple's contribution to S(theta) = sum (theta - theta_rho)^(-2): 2 Re (theta - g + i delta)^(-2) summed over g = +-gstar."""
    tot = mp.mpf(0)
    for g in (gstar, -gstar):
        x = theta - g
        tot += 2 * (x * x - delta * delta) / (x * x + delta * delta) ** 2
    return tot


def golden_min(f, lo, hi, iters):
    gr = (mp.sqrt(5) - 1) / 2
    a, b = lo, hi
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f1, f2 = f(x1), f(x2)
    for _ in range(iters):
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - gr * (b - a)
            f1 = f(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + gr * (b - a)
            f2 = f(x2)
    return (x1, f1) if f1 < f2 else (x2, f2)


def job_S(args):
    key, theta = args
    v = X.S_line(mp.mpf(theta))
    return key, nst(v.real, 0), mp.nstr(abs(v.imag), 3)


def job_level(args):
    n, d, label, gstar_s, gap_s, y_s, jstep, iters, keep = args
    gstar, gap, delta, y = mp.mpf(gstar_s), mp.mpf(gap_s), mp.mpf(d), mp.mpf(y_s)
    js = list(range(-20, 21, jstep))
    vals = [v_total(y, gstar + gap * j / 40, delta, gstar) for j in js]
    imin = min(range(len(vals)), key=lambda i: vals[i])
    rec = dict(n=n, d=d, label=label, y=y_s, v_min_grid=nst(vals[imin]), j_min=js[imin], v_ordinate=nst(vals[js.index(0)]))
    vmin, off = vals[imin], mp.mpf(0)
    j0 = js[imin]
    if abs(j0) < 20 and iters > 0:
        lo = gstar + gap * max(j0 - jstep, -20) / 40
        hi = gstar + gap * min(j0 + jstep, 20) / 40
        th, vv = golden_min(lambda t: v_total(y, t, delta, gstar), lo, hi, iters)
        if vv < vmin:
            vmin, off = vv, th - gstar
    else:
        off = gstar + gap * j0 / 40 - gstar
    rec['v_min_refined'] = nst(vmin)
    rec['offset_refined'] = nst(off)
    if keep:
        rec['vector'] = [nst(v) for v in vals]
    return rec


def job_fig(args):
    gstar_s, gap_s, d, y_s, chunk, den = args
    gstar, gap, delta, y = mp.mpf(gstar_s), mp.mpf(gap_s), mp.mpf(d), mp.mpf(y_s)
    out = []
    for j in chunk:
        off = gap * j / den
        out.append((j, nst(off), nst(v_total(y, gstar + off, delta, gstar))))
    return y_s, out


def run_pool(func, tasks, workers, dps):
    if workers <= 0:
        init_worker(dps)
        return [func(t) for t in tasks]
    with get_context('fork').Pool(workers, initializer=init_worker, initargs=(dps,)) as pool:
        return list(pool.imap_unordered(func, tasks, chunksize=1))


def bracket_first_negative(recs, key):
    prev = None
    for r in recs:
        if prev is not None and prev[key] >= 0 and r[key] < 0:
            return prev['y'], r['y']
        prev = r
    return None


def summarise(recs, key):
    main = [r for r in recs if not r['label'].startswith('x')]
    neg = [r for r in main if r[key] < 0]
    idx = [main.index(r) for r in neg]
    contiguous = (not idx) or (idx[-1] - idx[0] + 1 == len(idx))
    return dict(levels=len(main), negative=len(neg), all_negative=len(neg) == len(main), none_negative=not neg, contiguous=contiguous,
                first_negative_y=X.fmt(neg[0]['y'], 10) if neg else None, last_negative_y=X.fmt(neg[-1]['y'], 10) if neg else None,
                last_positive_y=X.fmt([r for r in main if r[key] >= 0][-1]['y'], 10) if len(neg) < len(main) and any(r[key] >= 0 for r in main) else None)


def plot(data_path, fig_path):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    with open(data_path) as fh:
        d = json.load(fh)['figure']
    off = np.array([float(x) for x in d['offset']])
    styles = [('-', 'o'), ('--', 's'), ('-.', '^'), (':', 'D')]
    fig, ax = plt.subplots(figsize=(5.0, 2.9))
    ax.axhline(0.0, color='0.7', lw=0.5)
    ax.axvline(0.0, color='0.85', lw=0.5)
    lo, hi = 0.0, 0.0
    for q, ((ls, mk), y) in enumerate(zip(styles, d['levels'])):
        v = np.array([float(x) for x in d['v'][y]])
        lo, hi = min(lo, v.min()), max(hi, v.max())
        ax.plot(off, v, color='black', ls=ls, marker=mk, ms=2.4, markevery=(3 * q, 12), label=r'$\varepsilon=%s$' % y)
    ax.set_yscale('symlog', linthresh=0.1, linscale=0.7)
    ax.set_xlim(off[0], off[-1])
    ax.set_ylim(min(-0.3, 2.0 * lo), max(0.3, 2.0 * hi))
    ticks = [t for t in (-10.0, -1.0, -0.1, 0.0, 0.1, 1.0, 10.0) if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
    ax.set_yticks(ticks)
    ax.set_yticklabels([('%g' % t) for t in ticks])
    ax.set_xlabel(r'$\theta-\gamma^*$   ($\gamma^*=%.4f$, gap $%.4f$)' % (float(d['gstar']), float(d['gap'])))
    ax.set_ylabel(r'$\rho_\varepsilon(\theta)$ (symmetric log scale)')
    ax.legend(loc='upper center', frameon=False, fontsize=6.5, ncol=4)
    fig.tight_layout()
    fig.savefig(fig_path)
    plt.close(fig)
    print('wrote', fig_path, '(levels', d['levels'], ', %d points, delta = %s, gap n = %s)' % (len(off), d['delta'], d['n']))


def main():
    ap = argparse.ArgumentParser(description='Table tf:tab:window and Figure tf:fig:window: what a single basepoint sees of a hypothetical quadruple (see the header of the file).')
    ap.add_argument('--quick', action='store_true', help='gap n = 29 only, a reduced grid, one process, writes to /tmp unless --out')
    ap.add_argument('--out', default=None)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--dps', type=int, default=30)
    ap.add_argument('--gaps', default=None, help='comma separated n of the gaps (n, n + 1); default 29,10000 (29 with --quick)')
    ap.add_argument('--deltas', default=None, help='comma separated delta; default 0.05,0.2,0.3,0.45 (0.45,1.5 with --quick)')
    ap.add_argument('--ks', default=None, help='comma separated k of the levels delta k/40; default 1..39 (a reduced set with --quick)')
    ap.add_argument('--jstep', type=int, default=None, help='stride of j on the 41 point grid; default 1 (4 with --quick)')
    ap.add_argument('--no-figure', action='store_true', help='skip the 4 x 121 point figure data')
    ap.add_argument('--plot', action='store_true', help='draw the figure from the JSON data instead of computing')
    ap.add_argument('--data', default=None)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--root', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    a = ap.parse_args()
    if a.plot:
        plot(a.data or os.path.join(a.root, 'data', 'tf_window.json'), a.fig or os.path.join(a.root, 'fig', 'tf_window.pdf'))
        return
    quick = a.quick
    gaps = [int(x) for x in a.gaps.split(',')] if a.gaps else ([29] if quick else GAPS)
    deltas = a.deltas.split(',') if a.deltas else (['0.45', '1.5'] if quick else DELTAS)
    ks = [int(x) for x in a.ks.split(',')] if a.ks else ([2, 8, 14, 20, 26, 32, 38] if quick else list(range(1, 40)))
    jstep = a.jstep or (4 if quick else 1)
    full = (not quick) and gaps == GAPS and deltas == DELTAS and ks == list(range(1, 40)) and jstep == 1
    iters = 12 if quick else 25
    biters = 8 if quick else 28
    workers = 0 if quick else a.workers
    dps = a.dps
    out = a.out or (os.path.join(tempfile.gettempdir(), 'tf_window_quick.json') if quick else 'tf_window.json')
    mp.mp.dps = dps
    table = X.read_zero_table(X.find_zeros_file(a.zeros))
    log('tf_window: gaps (n, n + 1) for n = %s, deltas %s, %d levels per delta plus the 0.999 level, grid step j/40 with j stride %d, dps %d, %d process(es)' % (gaps, deltas, len(ks), jstep, dps, max(workers, 1)))

    # Stage 1: the gaps, S(gamma*) and delta_c, the first order test over the gap.
    info = {}
    for n in gaps:
        g0, g1 = X.ordinate(n, table), X.ordinate(n + 1, table)
        gstar, gap = (g0 + g1) / 2, g1 - g0
        Sg = X.S_line(gstar)
        info[n] = dict(g0=g0, g1=g1, gstar=gstar, gap=gap, S=Sg.real, S_imag=abs(Sg.imag), delta_c=mp.sqrt(2 / Sg.real))
        log('n = %d: gamma_n = %s, gamma_{n+1} = %s, gamma* = %s, gap = %s, S(gamma*) = %s, delta_c = %s' % (
            n, X.fmt(g0, 12), X.fmt(g1, 12), X.fmt(gstar, 12), X.fmt(gap, 8), X.fmt(Sg.real, 12), X.fmt(info[n]['delta_c'], 10)))
    js_S = [j for j in range(-20, 21, jstep) if abs(j) < 20]
    tasks = [((n, j), nst(info[n]['gstar'] + info[n]['gap'] * j / 40, 6)) for n in gaps for j in js_S]
    Sgap, imag_max = {}, mp.mpf(0)
    for key, v, im in run_pool(job_S, tasks, workers, dps):
        Sgap[key] = mp.mpf(v)
        imag_max = max(imag_max, mp.mpf(im), info[key[0]]['S_imag'])
    first_order = {}
    for n in gaps:
        I = info[n]
        for d in deltas:
            dl = mp.mpf(d)
            across = [(Sgap[(n, j)] + S_quad(I['gstar'] + I['gap'] * j / 40, dl, I['gstar']), j) for j in js_S]
            smin, jmin = min(across)
            first_order[(n, d)] = dict(S_minus_2_over_delta2=I['S'] - 2 / dl ** 2, S_with_quadruple_at_ordinate=I['S'] + S_quad(I['gstar'], dl, I['gstar']),
                                       min_over_gap_of_S_with_quadruple=smin, j_of_min=jmin)
            log('  n = %d, delta = %s: S(gamma*) - 2/delta^2 = %s; with the mirror terms %s; min over j in %d..%d of S + quadruple = %s at j = %d' % (
                n, d, X.fmt(first_order[(n, d)]['S_minus_2_over_delta2'], 10), X.fmt(first_order[(n, d)]['S_with_quadruple_at_ordinate'], 10), js_S[0], js_S[-1], X.fmt(smin, 8), jmin))

    # Stage 2: the level scans.
    tasks = []
    for n in gaps:
        I = info[n]
        for d in deltas:
            dl = mp.mpf(d)
            levels = [('k%d' % k, dl * k / 40) for k in ks] + [('k40', mp.mpf('0.999') * dl)] + [('x' + e, mp.mpf(e)) for e in EXTRA.get((n, d), [])]
            for label, y in levels:
                keep = label.startswith('x') or (n == 10000 and d == '0.45' and label == 'k17')
                tasks.append((n, d, label, nst(I['gstar'], 6), nst(I['gap'], 6), nst(y, 6), jstep, iters, keep))
    tasks.sort(key=lambda t: -t[0])
    log('level scans: %d tasks' % len(tasks))
    recs = {}
    for r in run_pool(job_level, tasks, workers, dps):
        rec = dict(label=r['label'], y=mp.mpf(r['y']), v_min_grid=mp.mpf(r['v_min_grid']), j_min=r['j_min'], v_ordinate=mp.mpf(r['v_ordinate']),
                   v_min_refined=mp.mpf(r['v_min_refined']), offset_refined=mp.mpf(r['offset_refined']), vector=r.get('vector'))
        recs.setdefault((r['n'], r['d']), []).append(rec)
    for k in recs:
        recs[k].sort(key=lambda r: r['y'])
    summ = {k: dict(grid=summarise(v, 'v_min_grid'), refined=summarise(v, 'v_min_refined'), ordinate=summarise(v, 'v_ordinate')) for k, v in recs.items()}

    # Stage 3: the critical levels (bisection) where the sign changes.
    crit = {}
    for (n, d), rl in sorted(recs.items()):
        I = info[n]
        dl, gs, gp = mp.mpf(d), I['gstar'], I['gap']
        h = gp / 40
        entry = {}
        for tag, key in (('gap', 'v_min_refined'), ('ordinate', 'v_ordinate')):
            br = bracket_first_negative(rl, key)
            if br is None:
                entry[tag] = None
                continue
            lo, hi = br
            for _ in range(biters):
                mid = (lo + hi) / 2
                if tag == 'ordinate':
                    f = v_total(mid, gs, dl, gs)
                else:
                    f = min(v_total(mid, gs, dl, gs), golden_min(lambda t: v_total(mid, t, dl, gs), gs - 2 * h, gs + 2 * h, iters if iters > 0 else 10)[1])
                if f >= 0:
                    lo = mid
                else:
                    hi = mid
            entry[tag] = dict(bracket=[X.fmt(br[0], 10), X.fmt(br[1], 10)], critical_level=X.fmt((lo + hi) / 2, 10), bisection_steps=biters)
            log('critical level, n = %d, delta = %s, %s criterion: y_c = %s (grid bracket %s .. %s)' % (n, d, tag, X.fmt((lo + hi) / 2, 8), X.fmt(br[0], 8), X.fmt(br[1], 8)))
        crit[(n, d)] = entry

    # Stage 4: the figure data.
    fig_n = gaps[-1]
    fig = None
    if not a.no_figure:
        I = info[fig_n]
        fig_den = 120
        js_fig = list(range(-60, 61, 4 if quick else 1))
        chunks = [js_fig[i:i + 16] for i in range(0, len(js_fig), 16)]
        ftasks = [(nst(I['gstar'], 6), nst(I['gap'], 6), '0.45', y, ch, fig_den) for y in FIG_LEVELS for ch in chunks]
        figv = {y: [] for y in FIG_LEVELS}
        for y_s, res in run_pool(job_fig, ftasks, workers, dps):
            figv[y_s].extend(res)
        for y in FIG_LEVELS:
            figv[y].sort(key=lambda t: t[0])
        fig = dict(n=fig_n, delta='0.45', gstar=X.fmt(I['gstar'], 14), gap=X.fmt(I['gap'], 10), levels=FIG_LEVELS, points=len(js_fig), points_per_gap='j/120 for j = -60 .. 60 (stride %d here)' % (4 if quick else 1),
                   offset=[t[1] for t in figv[FIG_LEVELS[0]]], v={y: [t[2] for t in figv[y]] for y in FIG_LEVELS},
                   minima={y: X.fmt(min(mp.mpf(t[2]) for t in figv[y]), 12) for y in FIG_LEVELS})
        log('figure data at n = %d, delta = 0.45: minima over %d points: %s' % (fig_n, len(js_fig), ', '.join('y=%s: %s' % (y, mp.nstr(mp.mpf(fig['minima'][y]), 6)) for y in FIG_LEVELS)))

    # Stage 5: error control.
    checks = {}
    top = gaps[-1]
    I = info[top]
    gs, gp = I['gstar'], I['gap']
    dl = mp.mpf('0.45')
    pts = [(mp.mpf('0.25'), gs), (mp.mpf('0.05'), gs + mp.mpf('0.3') * gp), (mp.mpf('0.4'), gs - mp.mpf('0.1') * gp)]
    base = [v_total(y, t, dl, gs) for y, t in pts]
    mp.mp.dps = dps + 15
    hi_vals = [v_total(y, t, dl, gs) for y, t in pts]
    mp.mp.dps = dps
    checks['precision_v_relative_difference_dps_%d_vs_%d' % (dps, dps + 15)] = X.fmt(max(abs(b - c) / abs(c) for b, c in zip(base, hi_vals)), 3)
    y1, t1 = mp.mpf('0.001'), gs + mp.mpf('0.2') * gp
    S1 = X.S_line(t1)
    checks['small_y_limit_relative_difference_of_v_and_y_S_over_pi_at_y_0.001'] = X.fmt(abs(v_base(y1, t1) / (y1 * S1.real / mp.pi) - 1), 3)
    imag_max = max(imag_max, abs(S1.imag))
    y2 = mp.mpf('0.25')
    closed = (2 * y2 / (y2 ** 2 - dl ** 2) + sum((y2 - b) / ((y2 - b) ** 2 + (2 * gs) ** 2) for b in (dl, -dl))) / mp.pi
    checks['Cauchy_terms_at_the_ordinate_vs_closed_form_difference'] = X.fmt(abs(v_quad(y2, gs, dl, gs) - closed), 3)
    checks['largest_imaginary_part_of_S'] = mp.nstr(imag_max, 3)
    for k, v in checks.items():
        log('check  %s = %s' % (k, v))

    def recjson(r):
        d = dict(label=r['label'], y=X.fmt(r['y'], 12), v_min_grid=X.fmt(r['v_min_grid'], 14), j_min=r['j_min'], v_ordinate=X.fmt(r['v_ordinate'], 14),
                 v_min_refined=X.fmt(r['v_min_refined'], 14), offset_of_refined_minimum=X.fmt(r['offset_refined'], 8))
        if r['vector']:
            d['v_on_the_41_point_grid'] = [X.fmt(mp.mpf(x), 14) for x in r['vector']]
        return d

    comp = []

    def write_json():
        res = dict(script='tf_window.py', chapter='tf', table='tf:tab:window', figure_label='tf:fig:window',
                   parameters=dict(dps=dps, gaps_n=gaps, deltas=deltas, levels_per_delta='delta k/40 for k in %s and 0.999 delta' % ('1..39' if not quick else ks), grid='theta = gamma* + (j/40) gap, j = -20 .. 20 (stride %d)' % jstep,
                                   refinement_golden_steps=iters, bisection_steps=biters, workers=workers, quick=bool(quick), extra_levels={'%d,%s' % k: v for k, v in EXTRA.items()},
                                   v_definition='(1/pi) Re xi_prime/xi(1/2 + y + i theta) + (1/pi) sum over b = +-delta, g = +-gamma* of (y - b)/((y - b)^2 + (theta - g)^2)'),
                   gaps={str(n): dict(n=n, gamma_n=X.fmt(i['g0'], 20), gamma_n1=X.fmt(i['g1'], 20), gamma_star=X.fmt(i['gstar'], 20), gap=X.fmt(i['gap'], 15), S_gamma_star=X.fmt(i['S'], 20), delta_c=X.fmt(i['delta_c'], 15)) for n, i in info.items()},
                   first_order={'%d,%s' % k: {kk: (X.fmt(vv, 15) if isinstance(vv, mp.mpf) else vv) for kk, vv in v.items()} for k, v in first_order.items()},
                   levels={'%d,%s' % k: [recjson(r) for r in v] for k, v in recs.items()}, summary={'%d,%s' % k: v for k, v in summ.items()},
                   critical_levels={'%d,%s' % k: v for k, v in crit.items()}, figure=fig, checks=checks, comparison=comp, elapsed_seconds=round(time.time() - T0, 1))
        with open(out, 'w') as fh:
            json.dump(res, fh, indent=1)
        log('wrote ' + out)

    write_json()  # the data first, so that a failure of the comparison below cannot lose them

    # Stage 6: the table and the comparison with the chapter.
    def compare(name, claim, mine):
        dec = len(claim.split('.')[1]) if '.' in claim else 0
        unit = mp.mpf(10) ** (-dec)
        dif = abs(mine - mp.mpf(claim))
        status = 'agrees' if dif <= unit / 2 * (1 + mp.mpf(10) ** (-20)) else ('last digit differs by one unit' if dif < unit else 'differs')
        comp.append(dict(item=name, chapter=claim, computed=X.fmt(mine, 10), status=status, difference_in_units_of_last_printed_digit=X.fmt(dif / unit, 4)))

    def flag(name, chapter, mine, ok):
        comp.append(dict(item=name, chapter=chapter, computed=mine, status='agrees' if ok else 'differs'))

    def cell(s):
        if s['all_negative']:
            return 'all %d y<delta' % s['levels']
        if s['none_negative']:
            return 'none of %d' % s['levels']
        return '%s <= y <= %s (%d of %d)%s' % (mp.nstr(mp.mpf(s['first_negative_y']), 6), mp.nstr(mp.mpf(s['last_negative_y']), 6), s['negative'], s['levels'], '' if s['contiguous'] else ' NOT contiguous')

    log('')
    for n in gaps:
        I = info[n]
        log('Table tf:tab:window, gap n = %d: gamma* = %s, delta_c = %s' % (n, mp.nstr(I['gstar'], 9), mp.nstr(I['delta_c'], 6)))
        for d in deltas:
            s = summ[(n, d)]
            log('   delta = %-5s grid: %-46s refined: %-46s at the ordinate only: %s' % (d, cell(s['grid']), cell(s['refined']), cell(s['ordinate'])))
        if n in CLAIM:
            c = CLAIM[n]
            compare('gap n=%d: gamma*' % n, c['gstar'], I['gstar'])
            compare('gap n=%d: gap' % n, c['gap'], I['gap'])
            compare('gap n=%d: delta_c' % n, c['delta_c'], I['delta_c'])
    if full:
        for d in DELTAS:
            flag('n=29 delta=%s: every level negative' % d, 'all y<delta', cell(summ[(29, d)]['grid']), summ[(29, d)]['grid']['all_negative'])
        for d in ('0.05', '0.2', '0.3'):
            flag('n=10000 delta=%s: every level negative' % d, 'all y<delta', cell(summ[(10000, d)]['grid']), summ[(10000, d)]['grid']['all_negative'])
        g = summ[(10000, '0.45')]['grid']
        compare('n=10000 delta=0.45: first negative level', '0.19125', mp.mpf(g['first_negative_y']) if g['first_negative_y'] else mp.mpf(-1))
        compare('n=10000 delta=0.45: last negative level (0.999 delta)', '0.4496', mp.mpf(g['last_negative_y']) if g['last_negative_y'] else mp.mpf(-1))
        flag('n=10000 delta=0.45: number of negative levels', '24 of 40', '%d of %d' % (g['negative'], g['levels']), g['negative'] == 24 and g['levels'] == 40)
    if (10000, '0.45') in recs:
        compare('n=10000 delta=0.45: S(gamma*) - 2/delta^2', CLAIM_FIRST_ORDER, first_order[(10000, '0.45')]['S_minus_2_over_delta2'])
        rl = {r['label']: r for r in recs[(10000, '0.45')]}
        if 'x0.18' in rl:
            compare('n=10000 delta=0.45: min of v at y = 0.18 (sign change, text)', CLAIM_SIGN['0.18'], rl['x0.18']['v_min_grid'])
        if 'k17' in rl:
            compare('n=10000 delta=0.45: min of v at y = 0.19125 (sign change, text)', CLAIM_SIGN['0.19125'], rl['k17']['v_min_grid'])
        low_neg = [r['y'] for r in recs[(10000, '0.45')] if not r['label'].startswith('x') and r['y'] < mp.mpf('0.19') and r['v_min_refined'] < 0]
        flag('n=10000 delta=0.45: no grid level below 0.19 negative (refined minima)', 'v stays positive across the gap at every grid level below 0.19', 'grid levels below 0.19 negative: %d' % len(low_neg), not low_neg)
        cr = crit.get((10000, '0.45'), {}).get('ordinate')
        if cr:
            compare('n=10000 delta=0.45: critical level at the ordinate (text of tf:prop:curve: delta_det(y) = 0.45 gives y = 0.1895)', '0.1895', mp.mpf(cr['critical_level']))
        flag('n=10000 delta=0.45: first order test positive over the gap (j = -19 .. 19)', 'positive', X.fmt(first_order[(10000, '0.45')]['min_over_gap_of_S_with_quadruple'], 8),
             first_order[(10000, '0.45')]['min_over_gap_of_S_with_quadruple'] > 0)
    if fig is not None and fig_n == 10000:
        for y in FIG_LEVELS:
            compare('figure data, delta=0.45, y=%s: minimum of v over %d points of the gap' % (y, fig['points']), CLAIM_FIG_MIN[y], mp.mpf(fig['minima'][y]))
        flag('figure caption: positive at y = 0.05 and 0.15, negative at 0.25 and 0.4', 'signs + + - -', ' '.join('+' if mp.mpf(fig['minima'][y]) > 0 else '-' for y in FIG_LEVELS),
             [mp.mpf(fig['minima'][y]) > 0 for y in FIG_LEVELS] == [True, True, False, False])
    nbad = [c for c in comp if c['status'] != 'agrees']
    log('comparison with the chapter: %d items compared, %d agree, %d do not' % (len(comp), len(comp) - len(nbad), len(nbad)))
    for c in nbad:
        log('   %s: %s: chapter %s, computed %s' % (c['status'].upper(), c['item'], c['chapter'], c['computed']))

    write_json()


if __name__ == '__main__':
    main()
