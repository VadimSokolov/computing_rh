# phase_scaling.py: the cost of the flow route, Proposition cp:prop:cost and Figure cp:fig:scaling in Section cp:sec:cost (The cost
# of the flow route) of Chapter 14 (The Thorin phase diagram). Added for the book in October 2026, written from the mathematics.
# Data, hard-coded below (13 points: the barrier X and the bound Lambda; RH must be verified to height X/2): Polymath15, Effective
# approximation of heat flow evolution of the Riemann xi function, and a new upper bound for the de Bruijn Newman constant, Res. Math.
# Sci. 6 (2019), arXiv:1904.12438: (a) the theorem Lambda <= 0.22 (t0 = y0 = 0.2) with the barrier X = 6e10 + 83952 - 1/2; (b) Table 1,
# Conditional Lambda results, 12 rows: X = 2e12+129093, 5e12+194858, 2e13+131252, 6e13+123375, 3e14+188911, 2e15+122014, 7e15+68886,
# 6e16+156984, 6e17+88525, 9e18+35785, 2e20+66447, 9e21+70686 for Lambda = 0.21, 0.20, ..., 0.10. Only (a) and the row Lambda = 0.20
# (Platt and Trudgian, who verified RH to 3e12 and state that 2.51e12 suffices) are proved bounds; the other rows are conditional on a
# verification that has not been done. Transcribed from the arXiv PDF by the coordinating session; not re-read by this script's author.
# Fits. (1) Lambda = c / log X by least squares in 1/log X, c = sum Lambda_i/log X_i / sum (log X_i)^-2, over the 13 points and over
# the 12 rows (and, for reference, over the two proved points); (2) log X = a / Lambda + b by ordinary least squares over the 12 rows.
# For each fit: the maximal residual in Lambda and in log10 X; the barrier X = exp(c/delta) (law 1) or exp(a/delta + b) (law 2) for
# delta = 0.15, 0.1, 0.05, 0.01 with the height X/2 and its ratio to the record 3e12; the value of each law at the record (X = 6e12);
# the ratio log X(delta/2) / log X(delta) (2 for law 1: the height is squared; below 2 for law 2). The table gives the heights
# directly down to Lambda = 0.10, so only 0.05 and 0.01 are extrapolations.
# Error control: arithmetic only (mpmath, 30 digits). The data span a factor 1.9 in log X; the fits are order of magnitude statements and
# the output gives the spread between the laws. Output: data/phase_scaling.json (values as strings) and fig/cp_scaling.pdf (grayscale
# safe: line styles and marker fills). Usage: python3 phase_scaling.py [--quick] [--out FILE] [--fig FILE] [--figure] [--workers W];
# --figure redraws the figure from the JSON without computing; --quick writes to /tmp unless --out and --fig are given (it runs in a second).
import argparse
import json
import os
import tempfile
import time
from pathlib import Path

import mpmath as mp

mp.mp.dps = 30
T0 = time.time()
RECORD_HEIGHT = mp.mpf('3e12')          # Platt and Trudgian 2021: RH verified to 3e12
DELTAS = ['0.15', '0.1', '0.05', '0.01']
TABLE1 = [('2e12', 129093, '0.21'), ('5e12', 194858, '0.20'), ('2e13', 131252, '0.19'), ('6e13', 123375, '0.18'),
          ('3e14', 188911, '0.17'), ('2e15', 122014, '0.16'), ('7e15', 68886, '0.15'), ('6e16', 156984, '0.14'),
          ('6e17', 88525, '0.13'), ('9e18', 35785, '0.12'), ('2e20', 66447, '0.11'), ('9e21', 70686, '0.10')]


def fs(x, n=12):
    return mp.nstr(x, n)


def points():
    P = [dict(name='Polymath15 theorem, Lambda <= 0.22', kind='theorem', X=mp.mpf('6e10') + 83952 - mp.mpf('0.5'), bound=mp.mpf('0.22'),
              label='X = 6e10 + 83952 - 1/2, t0 = y0 = 0.2')]
    for lead, off, lam in TABLE1:
        P.append(dict(name='Table 1, Lambda <= ' + lam, kind='PT21_row' if lam == '0.20' else 'table_row', X=mp.mpf(lead) + off,
                      bound=mp.mpf(lam), label='X = %s + %d' % (lead, off)))
    return P


def fit_one(P):
    """c = sum (Lambda_i / L_i) / sum L_i^-2, L_i = log X_i: least squares of Lambda against 1/log X without intercept."""
    u = [1 / mp.log(p['X']) for p in P]
    return sum(p['bound'] * ui for p, ui in zip(P, u)) / sum(ui * ui for ui in u)


def fit_two(P):
    """log X = a / Lambda + b by ordinary least squares (y = log X on x = 1 / Lambda)."""
    x = [1 / p['bound'] for p in P]
    y = [mp.log(p['X']) for p in P]
    n = len(P)
    xb, yb = sum(x) / n, sum(y) / n
    a = sum((xi - xb) * (yi - yb) for xi, yi in zip(x, y)) / sum((xi - xb) ** 2 for xi in x)
    return a, yb - a * xb


def law1(c):
    return dict(lam=lambda lnX: c / lnX, lnx=lambda lam: c / lam)


def law2(a, b):
    return dict(lam=lambda lnX: a / (lnX - b), lnx=lambda lam: a / lam + b)


def resid(P, law):
    rl = max((abs(p['bound'] - law['lam'](mp.log(p['X']))), p['name']) for p in P)
    rx = max((abs(mp.log(p['X']) - law['lnx'](p['bound'])) / mp.log(10), p['name']) for p in P)
    return dict(max_abs_residual_Lambda=fs(rl[0], 6), worst_point_Lambda=rl[1], max_abs_residual_log10X=fs(rx[0], 6), worst_point_log10X=rx[1])


def required(law, d):
    lnX = law['lnx'](mp.mpf(d))
    l10X = lnX / mp.log(10)
    l10H = l10X - mp.log10(2)
    return dict(ln_X=fs(lnX, 10), log10_X=fs(l10X, 8), log10_height_X_over_2=fs(l10H, 8),
                ratio_height_to_record=fs(mp.mpf(10) ** (l10H - mp.log10(RECORD_HEIGHT)), 6))


def compute():
    P = points()
    P13, P12, P2 = P, P[1:], [P[0], P[2]]
    c13, c12, c2 = fit_one(P13), fit_one(P12), fit_one(P2)
    a, b = fit_two(P12)
    a13, b13 = fit_two(P13)
    laws = {'one_parameter_13_points': (law1(c13), dict(c=fs(c13, 12))),
            'one_parameter_12_rows': (law1(c12), dict(c=fs(c12, 12))),
            'one_parameter_two_proved_points': (law1(c2), dict(c=fs(c2, 12))),
            'two_parameter_12_rows': (law2(a, b), dict(a=fs(a, 12), b=fs(b, 12))),
            'two_parameter_13_points': (law2(a13, b13), dict(a=fs(a13, 12), b=fs(b13, 12)))}
    R = {}
    R['source'] = dict(paper='Polymath15, Res. Math. Sci. 6 (2019), arXiv:1904.12438, the theorem Lambda <= 0.22 and Table 1 (Conditional Lambda results)',
                       proved_points='the 0.22 theorem and the row Lambda = 0.20 (Platt and Trudgian, RH verified to 3e12, 2.51e12 suffices); the other rows are conditional',
                       transcription='by the coordinating session from the arXiv PDF; not re-read here (the web tool had no credits)')
    R['points'] = []
    for p in P:
        lnX = mp.log(p['X'])
        R['points'].append(dict(name=p['name'], kind=p['kind'], X=fs(p['X'], 16), X_as_stated=p['label'], log10_X=fs(lnX / mp.log(10), 10),
                                bound=fs(p['bound'], 3), product_Lambda_logX=fs(p['bound'] * lnX, 10), height_X_over_2=fs(p['X'] / 2, 10),
                                ratio_height_to_record=fs(p['X'] / 2 / RECORD_HEIGHT, 6),
                                residual_Lambda_one_parameter_13=fs(p['bound'] - c13 / lnX, 4),
                                residual_log10X_two_parameter_12=fs((lnX - a / p['bound'] - b) / mp.log(10), 4)))
    R['fits'] = {}
    for k, (law, par) in laws.items():
        R['fits'][k] = dict(par, **resid(P12 if '12' in k else (P2 if 'two_proved' in k else P13), law))
    R['required'] = [dict(delta=d, **{k: required(law, d) for k, (law, par) in laws.items()}) for d in DELTAS]
    Xrec = 2 * RECORD_HEIGHT
    R['record'] = dict(height=fs(RECORD_HEIGHT, 4), X=fs(Xrec, 4), law_value={k: fs(law['lam'](mp.log(Xrec)), 6) for k, (law, par) in laws.items()},
                       note='the row Lambda = 0.20 needs X = 5e12 + 194858 (height 2.5e12); the next row needs height 1.0e13 + 65626, '
                            'which is %s times the record' % fs(P[3]['X'] / 2 / RECORD_HEIGHT, 4))
    R['halving'] = {k: [dict(delta=d, ratio_log_X_half_delta_over_log_X_delta=fs(law['lnx'](mp.mpf(d) / 2) / law['lnx'](mp.mpf(d)), 8)) for d in DELTAS[:3]]
                     for k, (law, par) in laws.items()}
    R['halving_note'] = 'law 1: log X(delta/2) = 2 log X(delta), the required X is squared; law 2: the ratio is below 2 (X(delta/2) = X(delta)^2 e^{-b} at delta/2 for b > 0)'
    return R


def draw(R, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    c13 = float(R['fits']['one_parameter_13_points']['c'])
    a, b = float(R['fits']['two_parameter_12_rows']['a']), float(R['fits']['two_parameter_12_rows']['b'])
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    lx = np.linspace(10.0, 52.0, 400)
    lnX = lx * np.log(10)
    ax.plot(lx, c13 / lnX, color='black', lw=1.7, ls='-', label=r'$\Lambda=c/\log X$, $c=%.2f$ (13 points)' % c13)
    ax.plot(lx, a / (lnX - b), color='0.45', lw=1.6, ls='-.', label=r'$\log X=a/\Lambda+b$, $a=%.2f$, $b=%.2f$ (12 rows)' % (a, b))
    for d in (0.1, 0.05):
        ax.axhline(d, color='0.6', lw=0.8, ls='--')
    first = True
    for P in R['points']:
        x, y = float(P['log10_X']), float(P['bound'])
        if P['kind'] == 'theorem':
            ax.plot([x], [y], 'o', color='black', ms=6.5, label='Polymath15 theorem, $\\Lambda\\leq0.22$')
        elif P['kind'] == 'PT21_row':
            ax.plot([x], [y], 's', color='black', ms=6, label='row used by PT21, $\\Lambda\\leq0.2$')
        else:
            ax.plot([x], [y], 'o', mfc='white', mec='black', ms=4.5, mew=1.0, label='rows of Table 1' if first else None)
            first = False
    for lawc, (cx, col) in (('one', (c13, 'black')), ('two', (None, '0.45'))):
        xs = c13 / 0.05 / np.log(10) if lawc == 'one' else (a / 0.05 + b) / np.log(10)
        ax.plot([xs], [0.05], 'D', color=col, ms=4)
        ax.annotate('$10^{%.1f}$' % xs, (xs, 0.05), textcoords='offset points', xytext=(3, 6), fontsize=7, color=col)
    ax.set_xlabel(r'$\log_{10}X$ (RH verified to height $X/2$)')
    ax.set_ylabel(r'bound on $\Lambda$')
    ax.set_ylim(0, 0.24)
    ax.set_xlim(10, 52)
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    fig.tight_layout()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)


def main():
    ap = argparse.ArgumentParser(description='cost of the flow route (cp:prop:cost)')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--fig', default=None)
    ap.add_argument('--figure', action='store_true', help='redraw the figure from --out without computing')
    ap.add_argument('--workers', type=int, default=1, help='unused: the computation takes a second')
    a = ap.parse_args()
    tmp = tempfile.gettempdir()
    out = a.out or (os.path.join(tmp, 'phase_scaling_quick.json') if a.quick else 'data/phase_scaling.json')
    fig = a.fig or (os.path.join(tmp, 'cp_scaling_quick.pdf') if a.quick else 'fig/cp_scaling.pdf')
    if a.figure:
        draw(json.load(open(out)), fig)
        print('drew', fig)
        return
    R = dict(script='phase_scaling.py', mode='quick' if a.quick else 'full', **compute())
    R['seconds'] = round(time.time() - T0, 2)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(R, open(out, 'w'), indent=1)
    try:
        draw(R, fig)
        print('wrote', out, 'and', fig)
    except ImportError:
        print('wrote', out, '(matplotlib is missing: run with --figure where it is installed)')
    for k, v in R['fits'].items():
        print(k, {kk: vv for kk, vv in v.items() if kk in ('c', 'a', 'b', 'max_abs_residual_Lambda', 'max_abs_residual_log10X')})
    for row in R['required']:
        print('delta', row['delta'], {k: v['log10_X'] for k, v in row.items() if k != 'delta'})
    print('law values at the record', R['record']['law_value'])


if __name__ == '__main__':
    main()
