# wl_accuracy.py: Work package 1, "what to fit" item 2: the distance err_k(lambda, N) from the ordinate gamma_k of zeta to the nearest zero of xi_hat (the least eigenvector of the Weil form on the N + 1
# even modes), for every ordinate below pi N / L. `run` computes the cells (floating point at 3.5N + 60 digits and above, general routine of wl_dh.py, no certificate: the signs are in wl_pcert/wl_grid);
# `fit` reads the json and writes the fits (numpy least squares; only entries with -(working digits) + 8 < log10 err < -6 are used).
# Usage (from data/ on Hopper): python3 wl_accuracy.py run --procs 40 --out wl_acc.json ;  laptop: timeout 30 python3 -I code/wl_accuracy.py fit --inp data/wl_acc.json --out data/wl_accuracy_fit.json
import argparse, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def grid():
    return [dict(kind='zeta', parity='even', X=X, N=N, cert=False, roots=True, neig=85) for X in (13, 30, 50, 100, 200, 500, 1000, 2000) for N in (16, 24, 32, 40, 60, 80)]


def fit(inp, out):
    import numpy as np
    import mpmath as mp
    mp.mp.dps = 15
    cells = json.load(open(inp))
    gam = [float(mp.zetazero(k).imag) for k in range(1, 200)]
    rows = []
    for c in cells:
        if 'err_log10' not in c:
            continue
        L = math.log(c['X']) / 2
        eps = float(c['eps']) if c.get('eps') else None
        le = math.log10(eps) if eps and eps > 0 else None
        for k, e in enumerate(c['err_log10'], 1):
            if e is None or not (-c['D'] + 8 < e < -6):
                continue
            rows.append(dict(X=c['X'], N=c['N'], L=L, k=k, g=gam[k - 1], e=e, leps=le, D=c['D'], top=c['top']))
    res = dict(n_entries=len(rows), n_cells=len(cells))

    def lstsq(feat, sel, name):
        S = [r for r in rows if sel(r) and all(f(r) is not None for f in feat)]
        if len(S) < len(feat) + 3:
            return None
        A = np.array([[f(r) for f in feat] for r in S])
        y = np.array([r['e'] for r in S])
        co, *_ = np.linalg.lstsq(A, y, rcond=None)
        r_ = y - A @ co
        return dict(model=name, n=len(S), coef=[float(v) for v in co], rms=float(np.sqrt(np.mean(r_ ** 2))), max_abs=float(np.max(np.abs(r_))), y_rms=float(np.std(y)))
    one = lambda r: 1.0
    models = [
        ('e = a + b N', [one, lambda r: r['N']]),
        ('e = a + b N + c L + d k', [one, lambda r: r['N'], lambda r: r['L'], lambda r: r['k']]),
        ('e = a + b log10(eps)', [one, lambda r: r['leps']]),
        ('e = a + b log10(eps) + c gamma_k', [one, lambda r: r['leps'], lambda r: r['g']]),
        ('e = a + b log10(eps) + c gamma_k + d gamma_k L', [one, lambda r: r['leps'], lambda r: r['g'], lambda r: r['g'] * r['L']]),
        ('e = a + b log10(eps) + c gamma_k + d gamma_k L + f L', [one, lambda r: r['leps'], lambda r: r['g'], lambda r: r['g'] * r['L'], lambda r: r['L']]),
        ('e = a + b log10(eps) + c k + d k L', [one, lambda r: r['leps'], lambda r: r['k'], lambda r: r['k'] * r['L']]),
        ('e = a + b N + c L + d gamma_k + f gamma_k L + h N L', [one, lambda r: r['N'], lambda r: r['L'], lambda r: r['g'], lambda r: r['g'] * r['L'], lambda r: r['N'] * r['L']]),
        ('e = a + b log10(eps) + c L (k = 1 only)', [one, lambda r: r['leps'], lambda r: r['L']]),
    ]
    res['fits_all_k'] = [f for f in (lstsq(ft, lambda r: True, nm) for nm, ft in models) if f]
    res['fits_k1'] = [f for f in (lstsq(ft, lambda r: r['k'] == 1, nm + ' [k=1]') for nm, ft in models if ('gamma' not in nm and ' k' not in nm and 'k =' not in nm.replace('(k = 1 only)', ''))) if f]
    # slopes at fixed X (k = 1): digits gained per added mode, and at fixed N per doubling of lambda
    sl = {}
    for X in sorted(set(r['X'] for r in rows)):
        S = sorted([(r['N'], r['e']) for r in rows if r['X'] == X and r['k'] == 1])
        if len(S) >= 3:
            A = np.array([[1.0, n] for n, _ in S])
            co, *_ = np.linalg.lstsq(A, np.array([e for _, e in S]), rcond=None)
            sl[str(X)] = dict(points=S, digits_per_mode=float(-co[1]), intercept=float(co[0]))
    res['k1_slope_in_N_by_X'] = sl
    dl = {}
    for N in sorted(set(r['N'] for r in rows)):
        S = sorted([(r['L'], r['e']) for r in rows if r['N'] == N and r['k'] == 1])
        if len(S) >= 3:
            A = np.array([[1.0, l / math.log(2)] for l, _ in S])
            co, *_ = np.linalg.lstsq(A, np.array([e for _, e in S]), rcond=None)
            dl[str(N)] = dict(points=[(round(l, 4), e) for l, e in S], digits_per_doubling_of_lambda=float(-co[1]))
    res['k1_slope_in_L_by_N'] = dl
    # dependence on k at fixed (X, N): digits lost per ordinate, and against the distance to the top
    byk = {}
    for c in cells:
        if 'err_log10' in c and c['N'] in (40, 60) and c['X'] in (50, 200, 1000):
            es = [(k, e) for k, e in enumerate(c['err_log10'], 1) if e is not None and -c['D'] + 8 < e < -6]
            if len(es) >= 4:
                A = np.array([[1.0, k] for k, _ in es])
                co, *_ = np.linalg.lstsq(A, np.array([e for _, e in es]), rcond=None)
                byk['%d:%d' % (c['X'], c['N'])] = dict(n_ordinates=len(es), digits_per_ordinate=float(co[1]), e_k1=es[0][1], e_klast=es[-1][1], klast=es[-1][0], top=c['top'])
    res['by_k'] = byk
    res['k1_NL_N_le_40'] = lstsq([one, lambda r: r['N'] * r['L'], lambda r: r['N'], lambda r: r['L']], lambda r: r['k'] == 1 and r['N'] <= 40, 'e_1 = a + b N L + c N + d L (N <= 40)')
    best = [f for f in res['fits_all_k'] if f['model'].startswith('e = a + b log10(eps) + c gamma_k + d gamma_k L') and 'f L' not in f['model']]
    if best:
        b = best[0]
        res['summary'] = 'Accuracy law (fit over all (X, N, k) entries, rms %.1f digits of std %.1f): log10 err_k = %.2f + %.3f log10(eps_N) %+.3f gamma_k %+.3f gamma_k L; the least squares fit for k = 1 alone gives the coefficients listed under [k=1].' % (b['rms'], b['y_rms'], b['coef'][0], b['coef'][1], b['coef'][2], b['coef'][3])
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k.startswith('fits')}, indent=1)[:3000])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['run', 'fit'])
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--inp', default='wl_acc.json')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    if a.mode == 'run':
        from multiprocessing import Pool
        import wl_dh
        tasks = sorted(grid(), key=lambda t: -t['N'] * 10 - t['X'] / 1000)
        with Pool(a.procs) as p:
            res = list(p.imap_unordered(wl_dh._cell, tasks))
        json.dump(res, open(a.out or 'wl_acc.json', 'w'), indent=1)
    else:
        fit(a.inp, a.out or 'wl_accuracy_fit.json')
