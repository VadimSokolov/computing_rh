# fc_perron.py: group FC, item 2. Observation obs:perron and Conjecture conj:R of the consolidated ferromagnets chapter: the Perron
# decomposition Z_N(beta) = zeta(beta-1)/zeta(beta) + W_N(2-beta)/zeta(2) + R_N(beta) on Re beta = 3/2.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Inputs: the births file of fc_stern.py (exact counts c_N(m)), data/zeros_1700.txt. Z_N(beta) = sum_m c_N(m) m^(-beta); the profile
# w_N = c_N/phi; W_N(z) = (1/z) sum_{m>=2} (w_N(m-1) - w_N(m)) m^z (fc_lib.py); R_N = Z_N - zeta(beta-1)/zeta(beta) - W_N(2-beta)/zeta(2).
# Points: beta = 3/2 + i t for t in {0.3, 1, 2, 5, 10} (the table) and beta = 1 + rho_n = 3/2 + i gamma_n, n = 1..20 (the Riemann heights).
# Precision: double precision numpy for every N from 10 to --nmax (error of the sums below 1e-9, checked against the next line), and
# python-flint ball arithmetic at 128 bits (acb) for N in --arbN (the logged column ball_vs_double is the largest difference between
# the ball midpoint and the double value over A, B at all 25 points, and the ball radius is logged).
# Also the crossing t_N where |W_N(2-beta)|/zeta(2) first drops below |zeta(beta-1)/zeta(beta)| (grid 0.01, N <= --crossmax).
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fc_lib as L
from multiprocessing import Pool
import mpmath as mp

TS = [0.3, 1.0, 2.0, 5.0, 10.0]

def arb_sums(c, phi, betas, prec=128):
    from flint import arb, acb, fmpq, ctx
    ctx.prec = prec
    Mn = len(c) - 1
    lm = [arb(0)] + [arb(m).log() for m in range(1, Mn + 2)]
    def wq(m):
        return fmpq(int(c[m]), int(phi[m])) if 1 <= m <= Mn else fmpq(0)
    cA = [0] * (Mn + 2); cB = [arb(0)] * (Mn + 2)
    prev = wq(1)
    for m in range(1, Mn + 2):
        cA[m] = int(c[m]) if m <= Mn else 0
        if m >= 2:
            cur = wq(m)
            dq = (prev - cur) * (m * m)
            cB[m] = arb(int(dq.p)) / int(dq.q)
            prev = cur
    out = []
    for (bre, bim) in betas:
        b = acb(arb(bre), arb(bim))
        A = acb(0); B = acb(0)
        for m in range(1, Mn + 2):
            e = (-(b * lm[m])).exp()
            if cA[m]:
                A += cA[m] * e
            if m >= 2:
                B += cB[m] * e
        out.append((A, B))
    return out

def task(args):
    N, births, zerofile, arbN, crossmax, quick, nmaxg = args
    t0 = time.time()
    z = np.load(births)
    Mn = L.fib(N + 2)
    c = np.zeros(Mn + 1, dtype=np.int64); c[1] = 1
    for j in range(1, N + 1):
        h = z['H%d' % j].astype(np.int64); c[:len(h)] += h
    phi = get_phi(nmaxg)
    zs = [float(s) for s in L.zeros_list(zerofile, 20)]
    zstr = L.zeros_list(zerofile, 20)
    pts = [('t', t) for t in TS] + [('g', g) for g in zs]
    betas = [(1.5, p[1]) for p in pts]
    bc = [complex(*b) for b in betas]
    vals = L.Z_and_W(c, phi, bc)
    res = {'N': N, 'pts': []}
    for (kind, t), b, (A, B) in zip(pts, bc, vals):
        W = B / (2 - b)
        res['pts'].append({'kind': kind, 't': t, 'Z': [A.real, A.imag], 'W': [W.real, W.imag]})
    if N in arbN:
        bpts = [(1.5, TS[i]) for i in range(5)]
        from flint import arb, acb, ctx
        # exact-decimal t for the table; zeros from the 68 digit strings
        ctx.prec = 128
        arbbetas = []
        for (kind, t), s in zip(pts, [None] * 5 + zstr):
            pass
        arb_in = []
        for i, (kind, t) in enumerate(pts):
            if kind == 't':
                arb_in.append((arb('1.5'), arb(str(t))))
            else:
                arb_in.append((arb('1.5'), arb(zstr[i - 5])))
        # reuse arb_sums with (arb, arb) pairs
        def conv(bb):
            return bb
        from flint import fmpq
        Mnn = len(c) - 1
        res_arb = arb_sums_pairs(c, phi, arb_in)
        diffs = []; rads = []
        zeta2 = arb.pi() ** 2 / 6
        Rs = []
        for i, ((A, B), b) in enumerate(zip(res_arb, arb_in)):
            bb = acb(b[0], b[1])
            W = B / (2 - bb)
            zr = (bb - 1).zeta() / bb.zeta()
            R = A - zr - W / zeta2
            Rs.append(R)
            A0, B0 = vals[i]
            diffs.append(max(abs(float(A.real.mid()) - A0.real), abs(float(A.imag.mid()) - A0.imag), abs(float(B.real.mid()) - B0.real) / max(1, abs(B0)), abs(float(B.imag.mid()) - B0.imag) / max(1, abs(B0))))
            rads.append(float(max(R.real.rad(), R.imag.rad())))
            res['pts'][i]['R_ball'] = [float(R.real.mid()), float(R.imag.mid())]
            res['pts'][i]['zr_ball'] = [float(zr.real.mid()), float(zr.imag.mid())]
            res['pts'][i]['prof_ball'] = [float((W / zeta2).real.mid()), float((W / zeta2).imag.mid())]
        res['ball_vs_double'] = max(diffs)
        res['ball_radius_R'] = max(rads)
    # crossing
    if N <= crossmax:
        ts = np.arange(0.02, 8.0, 0.02)
        bb = [complex(1.5, t) for t in ts]
        v = L.Z_and_W(c, phi, bb)
        mp.mp.dps = 20
        prof = []; zr = []
        for t, (A, B), b in zip(ts, v, bb):
            prof.append(abs(B / (2 - b)) / (np.pi ** 2 / 6))
            zr.append(abs(complex(mp.zeta(b - 1) / mp.zeta(b))))
        prof = np.array(prof); zr = np.array(zr)
        below = np.nonzero(prof < zr)[0]
        res['cross_t'] = float(ts[below[0]]) if len(below) else None
        # last time profile exceeds ratio
        above = np.nonzero(prof >= zr)[0]
        res['last_above_t'] = float(ts[above[-1]]) if len(above) else None
    res['seconds'] = time.time() - t0
    print('N=%d done %.1f s' % (N, res['seconds']), flush=True)
    return res

def arb_sums_pairs(c, phi, arb_in, prec=128):
    from flint import arb, acb, fmpq, ctx
    ctx.prec = prec
    Mn = len(c) - 1
    lm = [arb(0)] + [arb(m).log() for m in range(1, Mn + 2)]
    cA = [0] * (Mn + 2); cB = [arb(0)] * (Mn + 2)
    def wq(m):
        return fmpq(int(c[m]), int(phi[m])) if 1 <= m <= Mn else fmpq(0)
    prev = wq(1)
    for m in range(1, Mn + 2):
        cA[m] = int(c[m]) if m <= Mn else 0
        if m >= 2:
            cur = wq(m)
            dq = (prev - cur) * (m * m)
            cB[m] = arb(int(dq.p)) / int(dq.q)
            prev = cur
    out = []
    for (br, bi) in arb_in:
        b = acb(br, bi)
        A = acb(0); B = acb(0)
        for m in range(1, Mn + 2):
            e = (-(b * lm[m])).exp()
            if cA[m]:
                A += cA[m] * e
            if m >= 2:
                B += cB[m] * e
        out.append((A, B))
    return out

_PHI = {}

def get_phi(nmax):
    if nmax not in _PHI:
        _PHI[nmax] = L.phi_sieve(L.fib(nmax + 2) + 2)
    return _PHI[nmax]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--births', default='fc_births.npz')
    ap.add_argument('--zeros', default='zeros_1700.txt')
    ap.add_argument('--nmin', type=int, default=10)
    ap.add_argument('--nmax', type=int, default=34)
    ap.add_argument('--arbN', default='10,12,16,20,23,24,26')
    ap.add_argument('--crossmax', type=int, default=24)
    ap.add_argument('--nproc', type=int, default=16)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='fc_perron.json')
    args = ap.parse_args()
    arbN = set(int(x) for x in args.arbN.split(',') if x)
    if args.quick:
        args.nmax = 14; args.nproc = 2; arbN = {10}; args.crossmax = 12
    Ns = list(range(args.nmin, args.nmax + 1))
    Ns = sorted(set(Ns))
    Ns = sorted(Ns, key=lambda N: -N)
    tasks = [(N, args.births, args.zeros, arbN, args.crossmax, args.quick, args.nmax) for N in Ns]
    with Pool(args.nproc) as pool:
        results = pool.map(task, tasks, chunksize=1)
    results.sort(key=lambda r: r['N'])
    zeta2 = np.pi ** 2 / 6
    mp.mp.dps = 25
    betas = {}
    out = {'N': {}}
    # zeta ratios (double) at the points
    pts = results[0]['pts']
    zr = []
    for p in pts:
        b = mp.mpf(3) / 2 + 1j * mp.mpf(p['t'])
        zr.append(complex(mp.zeta(b - 1) / mp.zeta(b)))
    for r in results:
        N = r['N']
        rows = []
        for p, z0 in zip(r['pts'], zr):
            Z = complex(*p['Z']); W = complex(*p['W'])
            prof = W / zeta2
            R = Z - z0 - prof
            rows.append({'kind': p['kind'], 't': p['t'], 'absZmzr': abs(Z - z0), 'absprof': abs(prof), 'ratio_act_pred': abs(Z - z0) / abs(prof),
                         'abszr': abs(z0), 'absR': abs(R)})
        hts = [x for x in rows if x['kind'] == 'g']
        entry = {'table': [x for x in rows if x['kind'] == 't'], 'heights': hts,
                 'meanR': float(np.mean([x['absR'] for x in hts])), 'maxR': float(np.max([x['absR'] for x in hts])),
                 'minprof': float(np.min([x['absprof'] for x in hts])), 'maxprof': float(np.max([x['absprof'] for x in hts])),
                 'meanprof': float(np.mean([x['absprof'] for x in hts])),
                 'cross_t': r.get('cross_t'), 'last_above_t': r.get('last_above_t')}
        if 'ball_vs_double' in r:
            entry['ball_vs_double'] = r['ball_vs_double']; entry['ball_radius_R'] = r['ball_radius_R']
            # R from balls
            Rb = [abs(complex(*p['R_ball'])) for p in r['pts'][5:]]
            entry['meanR_ball'] = float(np.mean(Rb)); entry['maxR_ball'] = float(np.max(Rb))
            entry['table_ball'] = [{'t': p['t'], 'absZmzr': abs(complex(*p['Z']) - complex(*p['zr_ball'])), 'absprof': abs(complex(*p['prof_ball']))} for p in r['pts'][:5]]
        out['N'][str(N)] = entry
        print('N=%d meanR=%.5f maxR=%.5f prof %.4f..%.4f cross=%s ball=%s' % (N, entry['meanR'], entry['maxR'], entry['minprof'], entry['maxprof'], entry['cross_t'], entry.get('ball_vs_double')), flush=True)
    # decay fits
    Nl = np.array([int(k) for k in out['N']]); mR = np.array([out['N'][k]['meanR'] for k in out['N']])
    fits = {}
    for lo, hi in [(12, 24), (16, 34), (24, 34), (12, 34), (20, 34)]:
        sel = (Nl >= lo) & (Nl <= hi)
        if sel.sum() >= 3:
            p = np.polyfit(np.log(Nl[sel]), np.log(mR[sel]), 1)
            q = np.polyfit(Nl[sel], np.log(mR[sel]), 1)
            fits['%d_%d' % (lo, hi)] = {'power_exponent': float(p[0]), 'exp_rate_per_N': float(q[0])}
    out['fits_meanR'] = fits
    loc = {}
    for N in Nl:
        if N + 4 in Nl and N >= 12:
            a = out['N'][str(N)]['meanR']; b = out['N'][str(N + 4)]['meanR']
            loc[str(N)] = float(np.log(b / a) / np.log((N + 4) / N))
    out['local_exponent_N_to_N+4'] = loc
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1)

if __name__ == '__main__':
    main()
