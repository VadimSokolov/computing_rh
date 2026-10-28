# fc_extra.py: group FC, additions to items 2 to 4 requested after the first runs. Written from the mathematics (the collaborator's scripts were not read, run or copied).
# --mode q1: for the Farey chain, Z_N(beta, x) = 2^(-beta) (1 + x^N + Q(x)), Q(1) = 2^beta sum_{k=1}^{N-1} c_k(beta) = 2^beta (Z_N(beta, 1) - 2 * 2^(-beta));
#   beta*(N) is the root of Q(1) = 2 (Q(1) is decreasing in beta: each term (t/2)^(-beta) with t >= 3 for the non extreme words, checked on a grid), found by
#   bisection in double precision (all terms positive), and Q(1) at beta = 2 and 2.2. Inputs: the sparse tables of fc_ly.py (--tables, several files).
# --mode wind: the winding number of C_N and Z_N around the circle |beta - (1 + rho_n)| = R (default R = 0.3, 160 points, uncertified) for the given N, split over
#   (type, n) tasks; also the maximum phase step (must be < pi for the count to be reliable).
# --mode means: for N = 10, 12, ..., --nmax, the means over the first twenty 1 + rho_n of |Z_N|, |R_N|, |W_N(2-beta)|/zeta(2) and |zeta(beta-1)/zeta(beta)|.
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fc_lib as L
import fc_fisher as FF
from multiprocessing import Pool
import mpmath as mp

def q1_task(args):
    N, tabfiles = args
    for f in tabfiles:
        z = np.load(f)
        if 'k%d' % N in z.files:
            k = z['k%d' % N].astype(np.int64); t = z['t%d' % N]; n = z['n%d' % N]
            break
    lt = np.log(t.astype(np.float64)); nf = n.astype(np.float64)
    mid = (k > 0) & (k < N)
    def Q1(b):
        return 2.0 ** b * float(np.sum(nf[mid] * np.exp(-b * lt[mid])))
    grid = np.arange(0.5, 6.01, 0.25)
    vals = [Q1(b) for b in grid]
    mono = all(a > b for a, b in zip(vals[:-1], vals[1:]))
    lo, hi = 0.5, 6.0
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if Q1(m) > 2.0:
            lo = m
        else:
            hi = m
    return {'N': N, 'beta_star': 0.5 * (lo + hi), 'Q1_at_2': Q1(2.0), 'Q1_at_2.2': Q1(2.2), 'Q1_at_1.5': Q1(1.5), 'monotone_on_grid': mono}

def wind_task(args):
    N, typ, n, births, zerofile, nmaxg, R, npts = args
    z = np.load(births)
    Mn = L.fib(N + 2)
    c = np.zeros(Mn + 1, dtype=np.int64); c[1] = 1
    for j in range(1, N + 1):
        h = z['H%d' % j].astype(np.int64); c[:len(h)] += h
    P = FF.Prep(c, FF.get_phi(nmaxg))
    g = float(L.zeros_list(zerofile, 20)[n - 1])
    cen = complex(1.5, g)
    ph = np.array([P.F(cen + R * np.exp(2j * np.pi * k / npts), typ)[0] for k in range(npts + 1)])
    dang = np.angle(ph[1:] / ph[:-1])
    return {'N': N, 'typ': typ, 'n': n, 'winding': float(dang.sum() / (2 * np.pi)), 'max_phase_step': float(np.max(np.abs(dang)))}

def means_task(args):
    N, births, zerofile, nmaxg = args
    z = np.load(births)
    Mn = L.fib(N + 2)
    c = np.zeros(Mn + 1, dtype=np.int64); c[1] = 1
    for j in range(1, N + 1):
        h = z['H%d' % j].astype(np.int64); c[:len(h)] += h
    P = FF.Prep(c, FF.get_phi(nmaxg))
    mp.mp.dps = 20
    zs = [float(s) for s in L.zeros_list(zerofile, 20)]
    aZ = []; aR = []; aP = []; aZR = []
    for g in zs:
        b = complex(1.5, g)
        A, B, _, _ = P.AB(b)
        W = B / (2 - b)
        zr = complex(mp.zeta(b - 1) / mp.zeta(b))
        prof = W / P.z2
        aZ.append(abs(A)); aP.append(abs(prof)); aZR.append(abs(zr)); aR.append(abs(A - zr - prof))
    return {'N': N, 'mean_absZ': float(np.mean(aZ)), 'mean_absR': float(np.mean(aR)), 'max_absR': float(np.max(aR)), 'mean_prof': float(np.mean(aP)), 'mean_zratio': float(np.mean(aZR)),
            'min_absZ': float(np.min(aZ)), 'max_absZ': float(np.max(aZ))}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', required=True)
    ap.add_argument('--tables', nargs='*', default=[])
    ap.add_argument('--nlist', default='30,32,34')
    ap.add_argument('--births', default='fc_births.npz')
    ap.add_argument('--zeros', default='zeros_1700.txt')
    ap.add_argument('--nmax', type=int, default=34)
    ap.add_argument('--radius', type=float, default=0.3)
    ap.add_argument('--npts', type=int, default=160)
    ap.add_argument('--nproc', type=int, default=8)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    if a.mode == 'q1':
        Ns = list(range(10, 37, 2))
        with Pool(min(a.nproc, len(Ns))) as p:
            res = p.map(q1_task, [(N, a.tables) for N in Ns])
    elif a.mode == 'wind':
        Ns = [int(x) for x in a.nlist.split(',')]
        tasks = [(N, typ, n, a.births, a.zeros, a.nmax, a.radius, a.npts) for N in sorted(Ns, reverse=True) for typ in ('raw', 'corr') for n in range(1, 21)]
        with Pool(a.nproc) as p:
            res = p.map(wind_task, tasks, chunksize=1)
    else:
        Ns = list(range(10, a.nmax + 1, 2))
        with Pool(a.nproc) as p:
            res = sorted(p.map(means_task, [(N, a.births, a.zeros, a.nmax) for N in sorted(Ns, reverse=True)], chunksize=1), key=lambda r: r['N'])
    with open(a.out, 'w') as f:
        json.dump(res, f, indent=1)
    for r in res[:60]:
        if a.mode != 'wind':
            print(r, flush=True)
    if a.mode == 'wind':
        bad = [r for r in res if abs(r['winding'] - 1) > 1e-6 or r['max_phase_step'] > 2.5]
        print('wind: %d tasks, %d not winding 1 or phase step > 2.5' % (len(res), len(bad)), bad[:5], flush=True)

if __name__ == '__main__':
    main()
