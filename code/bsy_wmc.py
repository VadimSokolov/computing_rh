# bsy_wmc.py: the numerical claims of Section bsy:sec:wmc of Chapter 35 (Theorem bsy:thm:wmc, weighted Monte Carlo for integrands flat at the boundary; Lemma bsy:lem:trap in one
# dimension; Remark bsy:rem:flat), group BSY-E. It fills Table bsy:tab:wmc: root mean square error (RMSE) of the Voronoi, trapezoid and Delaunay rules at N independent uniform nodes on
# Q = [-1,1]^d, for N = 1000 * 4^k, with the slopes fitted for N >= 4000 and their standard errors.
# Integrands. The bump g(x) = prod_i b(x_i/0.9), b(u) = exp(-1/(1-u^2)) for |u| < 1 and 0 otherwise (case (a) of the theorem with delta = 0.1; reference integral (0.9 int_{-1}^{1} b)^d computed
# by mpmath at 25 digits and checked against a run at 35 digits, never a hard coded constant); e^x in d = 1 and e^{x+y/2} in d = 2 (exact integrals 2 sinh 1 and 2 sinh 1 * 4 sinh(1/2)), which do not vanish on dQ.
# Rules. d = 1: Voronoi (cells bounded by the midpoints of consecutive nodes and clipped to [-1,1], so the weights are (x_(i+1) - x_(i-1))/2 inside and the end cells are [-1, (x_(1)+x_(2))/2] and [(x_(N-1)+x_(N))/2, 1], which coincide with the trapezoid weights except at the two end nodes)
# and the sorted trapezoid rule with the end points -1 and 1 added as nodes (Lemma bsy:lem:trap). d = 2: Voronoi with the cells clipped to Q and Delaunay (linear interpolation on the
# Delaunay triangulation of the nodes and the four corners). Voronoi cell areas: scipy.spatial.Voronoi of the nodes together with the reflections in the four sides of the nodes within a band b
# of a side, b = min(2, 8 sqrt(log N / N)); the bisector of a node and its reflection is the side itself and a reflected node is never closer to a point of Q than its original, so the cell
# of the node in this diagram is the cell clipped to Q. The area of a cell is the sum over its finite Voronoi ridges of d_ij l_ij / 4 (d_ij the distance of the two nodes, l_ij the ridge length). Checks logged:
# the weights sum to 2^d (Voronoi, Delaunay, trapezoid), and at N = 1000 and 4000 the band areas equal those with all nodes reflected. All three integrands of each dimension are integrated with the same nodes of a replication.
# Rows of the table: bump d=1 Voronoi, e^x d=1 trapezoid, bump d=2 Voronoi, bump d=2 Delaunay, e^{x+y/2} d=2 Delaunay; the other combinations (bump d=1 trapezoid, e^x d=1 Voronoi, e^{x+y/2} d=2 Voronoi) are stored too.
# Replications: 200 in d = 1, 30 in d = 2, seeds SeedSequence([20261025, d, rep]) (the nodes for the successive N of a replication are consecutive draws of that generator), stored.
# Bias checks (logged under 'bias'): for e^x the mean error of the trapezoid rule against the leading term of Lemma bsy:lem:trap, 4 (e - 1/e) / (2 (N+1)(N+2)) ~ 4.7008 / N^2 (ratio measured / predicted and mean error * N^2);
# for e^{x+y/2} the mean error * N (Remark bsy:rem:flat: order N^(-2/d) = N^(-1)); for the bump the mean error against its standard error (case (a): exponentially small).
# Precision: the sums of weights times values are repeated with math.fsum (exact summation) on 3 replications per N (d = 1 all N, d = 2 up to 4000) and the largest absolute difference is logged
# next to the RMSE; the reference integrals are checked at two mpmath precisions. Run times: d = 1 and d = 2 up to N = 16000 take seconds with --procs; N = 256000 and 1024000 in d = 2 (about a minute per replication) go to Hopper.
# Options: --quick (N = 1000, 4000, 4 replications), --out FILE (json; the log has the extension .log), --procs P, --dims 1|2|12, --maxN1 / --maxN2 (largest N per dimension, default 1024000).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Run: python3 -I code/bsy_wmc.py [--quick] [--out data/bsy_wmc.json] [--procs 12]   (Hopper: bash misc/tools/hopper_run.sh -c 16 -m 64G -t 120 -n bsy_wmc -g bsy_wmc.json -g bsy_wmc.log code/bsy_wmc.py -- --procs 12 --out bsy_wmc.json)
import argparse, json, math, os, sys, time
import numpy as np
import multiprocessing as mp
import mpmath as mpm
from scipy.spatial import Voronoi, Delaunay

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
ap.add_argument('--out', default=None)
ap.add_argument('--procs', type=int, default=1)
ap.add_argument('--dims', default='12')
ap.add_argument('--maxN1', type=int, default=1024000)
ap.add_argument('--maxN2', type=int, default=1024000)
ap.add_argument('--reps1', type=int, default=200)
ap.add_argument('--reps2', type=int, default=30)
args = ap.parse_args()
OUT = args.out or ('misc/reviews/nick-b4/quick/BSY-E/bsy_wmc.json' if args.quick else 'data/bsy_wmc.json')
LOGF = os.path.splitext(OUT)[0] + '.log'
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
LOG = open(LOGF, 'w')
SEED0 = 20261025
NFIT = 4000


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True); LOG.write(s + '\n'); LOG.flush()


def Nlist(maxN):
    out, n = [], 1000
    while n <= maxN:
        out.append(n); n *= 4
    return out


if args.quick:
    N1 = N2 = [1000, 4000]
    REPS = {1: 4, 2: 4}
else:
    N1, N2 = Nlist(args.maxN1), Nlist(args.maxN2)
    REPS = {1: args.reps1, 2: args.reps2}
NS = {1: N1, 2: N2}


def bump1(u):
    out = np.zeros_like(u)
    m = np.abs(u) < 1
    out[m] = np.exp(-1.0 / (1.0 - u[m] ** 2))
    return out


def g_bump(P):
    # P: (n,) or (n,2)
    if P.ndim == 1:
        return bump1(P / 0.9)
    return bump1(P[:, 0] / 0.9) * bump1(P[:, 1] / 0.9)


def g_exp(P):
    return np.exp(P) if P.ndim == 1 else np.exp(P[:, 0] + 0.5 * P[:, 1])


def exact(dps):
    mpm.mp.dps = dps
    b = lambda u: mpm.exp(-1 / (1 - u ** 2))
    I1 = mpm.mpf('0.9') * mpm.quad(b, [-1, 0, 1])
    s1 = 2 * mpm.sinh(1)
    s2 = 2 * mpm.sinh(1) * 4 * mpm.sinh(mpm.mpf(1) / 2)
    return {('bump', 1): I1, ('bump', 2): I1 ** 2, ('exp', 1): s1, ('exp', 2): s2}


# ----------------------------------------------------------------------------- weights
def w1_voronoi(x):
    # x sorted, in (-1, 1)
    mid = 0.5 * (x[1:] + x[:-1])
    edges = np.concatenate([[-1.0], mid, [1.0]])
    return np.diff(edges)


def w1_trap(x):
    e = np.concatenate([[-1.0], x, [1.0]])
    w = np.empty(len(e))
    d = np.diff(e)
    w[0] = 0.5 * d[0]; w[-1] = 0.5 * d[-1]
    w[1:-1] = 0.5 * (d[:-1] + d[1:])
    return e, w


def vor_areas(P, full=False):
    N = len(P)
    b = 2.0 if full else min(2.0, 8 * math.sqrt(math.log(N) / N))
    pts = [P]
    for ax in (0, 1):
        for sg in (-1.0, 1.0):
            m = sg * P[:, ax] > 1 - b
            Q = P[m].copy()
            Q[:, ax] = 2 * sg - Q[:, ax]
            pts.append(Q)
    A = np.vstack(pts)
    V = Voronoi(A, qhull_options='Qbb Qc Qz')
    rp = V.ridge_points
    rv = np.array(V.ridge_vertices)
    ok = (rv >= 0).all(1)
    rp = rp[ok]; rv = rv[ok]
    # only ridges that have at least one original node matter
    l = np.linalg.norm(V.vertices[rv[:, 0]] - V.vertices[rv[:, 1]], axis=1)
    d = np.linalg.norm(A[rp[:, 0]] - A[rp[:, 1]], axis=1)
    c = 0.25 * d * l
    area = np.bincount(rp[:, 0], weights=c, minlength=len(A)) + np.bincount(rp[:, 1], weights=c, minlength=len(A))
    return area[:N]


def del_weights(P):
    C = np.array([[-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0]])
    A = np.vstack([P, C])
    T = Delaunay(A)
    s = T.simplices
    a, b, c = A[s[:, 0]], A[s[:, 1]], A[s[:, 2]]
    ar = 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0]))
    w = np.bincount(s.ravel(), weights=np.repeat(ar / 3.0, 3), minlength=len(A))
    return A, w


# ----------------------------------------------------------------------------- one replication
def rows_of(d):
    return {1: ['bump|voronoi', 'bump|trap', 'exp|voronoi', 'exp|trap'], 2: ['bump|voronoi', 'bump|delaunay', 'exp|voronoi', 'exp|delaunay']}[d]


def one_rep(task):
    d, rep = task
    rng = np.random.default_rng(np.random.SeedSequence([SEED0, d, rep]))
    EX = EXACT
    out = {}
    chk = {}
    for N in NS[d]:
        if d == 1:
            x = np.sort(rng.random(N) * 2 - 1)
            wv = w1_voronoi(x)
            e, wt = w1_trap(x)
            vals = {}
            for gn, gf in (('bump', g_bump), ('exp', g_exp)):
                p1 = wv * gf(x)
                p2 = wt * gf(e)
                vals[gn + '|voronoi'] = float(p1.sum()) - EX[(gn, 1)]
                vals[gn + '|trap'] = float(p2.sum()) - EX[(gn, 1)]
                if rep < 3:
                    chk[(gn, N, 'v')] = abs(math.fsum(p1) - float(p1.sum()))
                    chk[(gn, N, 't')] = abs(math.fsum(p2) - float(p2.sum()))
            chk[('wsum', N)] = max(abs(wv.sum() - 2), abs(wt.sum() - 2))
        else:
            P = rng.random((N, 2)) * 2 - 1
            av = vor_areas(P)
            A, wd = del_weights(P)
            vals = {}
            for gn, gf in (('bump', g_bump), ('exp', g_exp)):
                p1 = av * gf(P)
                p2 = wd * gf(A)
                vals[gn + '|voronoi'] = float(p1.sum()) - EX[(gn, 2)]
                vals[gn + '|delaunay'] = float(p2.sum()) - EX[(gn, 2)]
                if rep < 3 and N <= 4000:
                    chk[(gn, N, 'v')] = abs(math.fsum(p1) - float(p1.sum()))
                    chk[(gn, N, 't')] = abs(math.fsum(p2) - float(p2.sum()))
            chk[('wsum', N)] = max(abs(av.sum() - 4), abs(wd.sum() - 4))
            if rep == 0 and N <= 4000:
                chk[('band_vs_full', N)] = float(np.abs(av - vor_areas(P, full=True)).max())
        out[N] = vals
    return d, rep, out, chk


EXACT = None


def fit(Ns, rm):
    Ns = np.array(Ns, float); rm = np.array(rm, float)
    sel = Ns >= NFIT
    X = np.log(Ns[sel]); Y = np.log(rm[sel])
    if len(X) < 2:
        return float('nan'), float('nan')
    A = np.vstack([X, np.ones_like(X)]).T
    coef = np.linalg.lstsq(A, Y, rcond=None)[0]
    n = len(X)
    r = Y - A @ coef
    se = math.sqrt(float(r @ r) / (n - 2) / float(((X - X.mean()) ** 2).sum())) if n > 2 else float('nan')
    return float(coef[0]), se


def main():
    global EXACT
    t0 = time.time()
    log('bsy_wmc.py', 'quick' if args.quick else 'full', 'N1', N1, 'N2', N2, 'reps', REPS, 'procs', args.procs, 'dims', args.dims)
    e25, e35 = exact(25), exact(35)
    EXACT = {k: float(v) for k, v in e25.items()}
    dif = {str(k): float(abs(e25[k] - e35[k])) for k in e25}
    log('reference integrals (mpmath 25 digits):', {str(k): mpm.nstr(v, 20) for k, v in e25.items()})
    log('difference 25 against 35 digits:', dif)
    res = {'quick': args.quick, 'NS': NS, 'reps': REPS, 'seed': [SEED0, 'd', 'rep'], 'exact': {str(k): mpm.nstr(v, 25) for k, v in e25.items()}, 'exact_prec_diff': dif, 'rows': {}}
    dims = [int(c) for c in args.dims]
    tasks = [(d, r) for d in dims for r in range(REPS[d])]
    tasks.sort(key=lambda t: -t[0])  # the slow plane replications first
    if args.procs > 1:
        with mp.get_context('fork').Pool(args.procs) as pool:
            outs = pool.map(one_rep, tasks, chunksize=1)
    else:
        outs = [one_rep(t) for t in tasks]
    log('replications done in %.1f s' % (time.time() - t0))
    chk_all = {}
    for d in dims:
        O = sorted([o for o in outs if o[0] == d], key=lambda o: o[1])
        for o in O:
            for k, v in o[3].items():
                chk_all.setdefault(k, []).append(v)
        for row in rows_of(d):
            gn, rule = row.split('|')
            err = np.array([[o[2][N][row] for N in NS[d]] for o in O])  # reps x N
            rm = np.sqrt((err ** 2).mean(0)); me = err.mean(0); se_me = err.std(0, ddof=1) / math.sqrt(len(O))
            sl, se = fit(NS[d], rm)
            key = '%s|d%d|%s' % (gn, d, rule)
            res['rows'][key] = {'N': NS[d], 'rmse': rm.tolist(), 'mean_error': me.tolist(), 'mean_error_se': se_me.tolist(), 'slope': sl, 'slope_se': se, 'reps': len(O)}
            log('%-24s slope %.3f +- %.3f  RMSE %s' % (key, sl, se, ' '.join('%d:%.3e' % (N, r) for N, r in zip(NS[d], rm))))
    # bias checks
    bias = {}
    if 1 in dims:
        r = res['rows']['exp|d1|trap']
        pred = [2 * (math.e - 1 / math.e) / ((N + 1) * (N + 2)) for N in NS[1]]
        bias['exp_d1_trap'] = {'N': NS[1], 'mean_error_times_N2': [m * N * N for m, N in zip(r['mean_error'], NS[1])], 'ratio_to_lemma': [m / p for m, p in zip(r['mean_error'], pred)], 'lemma_constant': 4 * (math.e - 1 / math.e) / 2}
        log('bias, e^x d=1 trapezoid: mean error * N^2 =', ['%.4f' % v for v in bias['exp_d1_trap']['mean_error_times_N2']], '(lemma constant %.4f); ratio to the leading term' % bias['exp_d1_trap']['lemma_constant'], ['%.4f' % v for v in bias['exp_d1_trap']['ratio_to_lemma']])
        r = res['rows']['bump|d1|trap']
        log('bump d=1 trapezoid: slope %.3f; mean error / standard error:' % r['slope'], ['%.2f' % (m / s) for m, s in zip(r['mean_error'], r['mean_error_se'])])
        r = res['rows']['bump|d1|voronoi']
        log('bump d=1 Voronoi mean error / standard error:', ['%.2f' % (m / s) for m, s in zip(r['mean_error'], r['mean_error_se'])])
    if 2 in dims:
        for key in ('exp|d2|delaunay', 'exp|d2|voronoi'):
            r = res['rows'][key]
            bias[key] = {'N': NS[2], 'mean_error_times_N': [m * N for m, N in zip(r['mean_error'], NS[2])], 'mean_error_se_times_N': [s * N for s, N in zip(r['mean_error_se'], NS[2])]}
            log('bias,', key, ': mean error * N =', ['%.4f' % v for v in bias[key]['mean_error_times_N']], 'standard errors', ['%.4f' % v for v in bias[key]['mean_error_se_times_N']])
        for key in ('bump|d2|voronoi', 'bump|d2|delaunay'):
            r = res['rows'][key]
            log(key, 'mean error / standard error:', ['%.2f' % (m / s) for m, s in zip(r['mean_error'], r['mean_error_se'])])
    res['bias'] = bias
    # precision checks
    fs = {}
    for k, v in chk_all.items():
        fs[str(k)] = float(max(v))
    mx_w = max([v for k, v in fs.items() if 'wsum' in k] or [0.0])
    mx_fs = max([v for k, v in fs.items() if 'wsum' not in k and 'band' not in k] or [0.0])
    mx_band = max([v for k, v in fs.items() if 'band' in k] or [0.0])
    res['precision'] = {'max_weight_sum_error': mx_w, 'max_fsum_difference': mx_fs, 'max_band_vs_full_area_difference': mx_band, 'detail': fs}
    log('precision: weights sum to |Q| within %.2e; fsum against pairwise sum differs by at most %.2e; band against full reflection areas differ by at most %.2e' % (mx_w, mx_fs, mx_band))
    res['seconds'] = time.time() - t0
    with open(OUT, 'w') as f:
        json.dump(res, f)
    log('wrote', OUT, 'in %.1f s' % res['seconds'])


if __name__ == '__main__':
    main()
