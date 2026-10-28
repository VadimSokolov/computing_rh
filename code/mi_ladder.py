# mi_ladder.py: the cosh ladder for zeta in Chapter 22 (Nick Polson's Miyatake chapter), Section "The ladder and the quasi Riemann abscissa":
# the Thorin path gamma_k(b) (first twelve atoms against c = 1/2 + b, the figure of Proposition BTheta), monotonicity in b, the number of real
# zeros of A_c(t) = Re xi(c + i t) below height 100 (29 = N(100)), the strict interlacing of the zeros of A_c and B_c = Im xi(c + i t) on [0, 100]
# (Theorem paired, C:605 to C:609), the phase count theta_c(T) - theta_c(-T) = pi (N(T) - 2 N_{>c}(T)) + O(log T) (Theorem density, C:620), and the
# arithmetic of the last paragraph (C:654, an off line zero above 3e12).
# Conventions. xi(s) = (1/2) s (s - 1) pi^(-s/2) Gamma(s/2) zeta(s), Xi(t) = xi(1/2 + i t), F_b(z) = (Xi(z + i b) + Xi(z - i b))/2, on R F_b(x) = Re xi(c + i x)
# with c = 1/2 + b; gamma_k(b) the positive zeros of F_b in increasing order; A_c(t) = Re xi(c + i t), B_c(t) = Im xi(c + i t) (B_c(0) = 0 counts as the
# first zero of B_c, A_c(0) = xi(c) > 0). theta_c(t) = arg xi(c + i t) continuous with theta_c(0) = 0; theta_c(-T) = -theta_c(T) since xi(conj s) = conj xi(s).
# Method. Zeros: sign changes on a grid (step 0.05, repeated at 0.04) refined by regula falsi; Thorin path: continuation in c from 0.5 (gamma_k(1/2) =
# the zeros of zeta, mp.zetazero) in steps 0.01 with secant steps, checked against the sign change scan at the tabulated c. Phase: sum of arg of ratios of
# consecutive values on a grid of step 0.02 (the largest step is logged; each step is far below pi), N(T) from data/zeros_1700.txt.
# Error control: the tabulated gamma_k are computed at 20 and 30 digits (difference logged); counts at two grid steps.
# Options: --task {path,counts,density,arith,all} --quick --out FILE.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse, os
from multiprocessing import Pool
import mpmath as mp


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def xval(job):
    c, t, dps = job
    mp.mp.dps = dps
    v = xi(mp.mpc(c, t))
    return (float(v.real), float(v.imag), float(mp.arg(v)))


def pmap(jobs, chunk=16):
    with Pool(min(8, os.cpu_count() or 1)) as p:
        return p.map(xval, jobs, chunksize=chunk)


class Log:
    def __init__(self, path):
        self.f = open(path, 'a')

    def __call__(self, *a):
        s = ' '.join(str(x) for x in a)
        print(s)
        self.f.write(s + '\n')
        self.f.flush()


def zeros_of(vals, ts, part):
    out = []
    for i in range(len(ts) - 1):
        a, b = vals[i][part], vals[i + 1][part]
        if a == 0 or a * b < 0:
            out.append((ts[i], ts[i + 1]))
    return out


def refine(c, a, b, part, dps=20):
    mp.mp.dps = dps
    g = (lambda t: mp.re(xi(mp.mpc(c, t)))) if part == 0 else (lambda t: mp.im(xi(mp.mpc(c, t))))
    return mp.findroot(g, (a, b), solver='illinois', tol=mp.mpf(10) ** (-dps + 4), maxsteps=60)


def task_path(log, R, quick):
    cs_tab = [1.0, 0.9, 0.8, 0.7, 0.6, 0.55, 0.5]
    nz = 12
    out = {}
    res = {}
    for dps in (20, 30):
        mp.mp.dps = dps
        g = [mp.zetazero(k).imag for k in range(1, nz + 1)]
        path = {}
        cgrid = [mp.mpf(50 - i) / 100 + mp.mpf(50) / 100 for i in range(0, 51)]  # 1.00 down to 0.50
        cgrid = [mp.mpf(100 - i) / 100 for i in range(0, 51)]
        # continuation from c = 0.5 upward
        cur = list(g)
        up = {}
        up[0] = list(cur)
        for i in range(1, 51):
            c = mp.mpf(50 + i) / 100
            new = []
            for k in range(nz):
                f = lambda t: mp.re(xi(mp.mpc(c, t)))
                # secant from previous c root and a nearby point
                new.append(mp.findroot(f, (cur[k], cur[k] - mp.mpf('0.02')), solver='secant', tol=mp.mpf(10) ** (-dps + 6), maxsteps=40))
            cur = new
            up[i] = list(cur)
        res[dps] = up
    up20, up30 = res[20], res[30]
    # monotone decreasing in c
    cs_sorted = sorted(up30.keys())
    mono = True
    viol = []
    for k in range(nz):
        for a, b in zip(cs_sorted[:-1], cs_sorted[1:]):
            if not (up30[b][k] < up30[a][k]):
                mono = False
                viol.append((k + 1, 0.5 + a / 100, 0.5 + b / 100))
    log('gamma_k(c) strictly decreasing in c on the grid c = 0.50, 0.51, ..., 1.00 for k = 1..12:', mono, viol[:5])
    diff = max(abs(up20[c][k] - up30[c][k]) for c in cs_sorted for k in range(nz))
    log('max |gamma_k(20 digits) - gamma_k(30 digits)| over the grid:', mp.nstr(diff, 3))
    out['monotone'] = mono
    out['dps_diff'] = mp.nstr(diff, 3)
    # gaps: smallest gap between consecutive gamma_k on the grid (distinct atoms)
    mingap = min(up30[c][k + 1] - up30[c][k] for c in cs_sorted for k in range(nz - 1))
    log('smallest gap between consecutive atoms gamma_{k+1} - gamma_k on the grid:', mp.nstr(mingap, 5))
    out['mingap'] = mp.nstr(mingap, 5)
    rows = []
    log('c      ' + ' '.join('g%02d' % (k + 1) for k in range(nz)))
    for c in cs_tab:
        key = int(round((c - 0.5) * 100))
        row = [mp.nstr(v, 8) for v in up30[key]]
        log(c, ' '.join(row))
        rows.append({'c': c, 'gamma': row})
    out['table'] = rows
    # check against sign change scan at the tabulated c (first 12 sign changes on t in (0, 60))
    chk = []
    for c in (1.0, 0.8, 0.6):
        ts = [0.05 * k for k in range(1, 1201)]
        vals = pmap([(c, t, 15) for t in ts])
        zs = zeros_of(vals, ts, 0)[:nz]
        mx = max(abs(float(mp.findroot(lambda t: mp.re(xi(mp.mpc(c, t))), (a, b), solver='illinois')) - float(up30[int(round((c - 0.5) * 100))][k])) for k, (a, b) in enumerate(zs))
        log('c =', c, 'first 12 sign changes of Re xi(c + it) agree with the continuation to', mx)
        chk.append((c, mx, len(zs)))
    out['scan_check'] = chk
    R['path'] = out


def task_counts(log, R, quick):
    cl = [0.5, 0.55, 0.6, 0.7, 0.8, 1.0, 1.25, 1.5]
    if os.environ.get('MI_CS'):
        cl = [float(x) for x in os.environ['MI_CS'].split(',')]
    T = 100.0
    out = []
    for c in cl:
        row = {'c': c}
        for h in (0.05, 0.04):
            ts = [h * k for k in range(1, int(T / h) + 1)]
            vals = pmap([(c, t, 15) for t in ts])
            za = zeros_of(vals, ts, 0)
            row['A_h%g' % h] = len(za)
            if c > 0.5:
                zb = zeros_of(vals, ts, 1)
                row['B_pos_h%g' % h] = len(zb)
                if h == 0.05:
                    A = [float(refine(c, a, b, 0)) for a, b in za]
                    B = [0.0] + [float(refine(c, a, b, 1)) for a, b in zb]
                    ev = sorted([(z, 'A') for z in A] + [(z, 'B') for z in B])
                    seq = ''.join(e[1] for e in ev)
                    alt = all(seq[i] != seq[i + 1] for i in range(len(seq) - 1))
                    row['order_alternates'] = alt
                    row['seq_head'] = seq[:12]
                    row['A_first'] = A[0]
                    row['min_gap_A_B'] = min(abs(x - y) for x in A for y in B)
        log('c =', c, row)
        out.append(row)
    R['counts'] = out


def task_density(log, R, quick):
    zs = [mp.mpf(l.strip()[:30]) for l in open('data/zeros_1700.txt') if l.strip()]
    def Nplus(T):
        return sum(1 for z in zs if z <= T)
    out = []
    h = 0.02
    for c in (0.6, 0.75, 1.0):
        ts = [h * k for k in range(0, int(100 / h) + 1)]
        ts[0] = 1e-6   # xi(c) is real positive; avoid the removable pole of zeta at s = 1 for c = 1
        vals = pmap([(c, t, 15) for t in ts], chunk=64)
        th = 0.0
        phase = {}
        worst = 0.0
        import cmath
        prev = complex(vals[0][0], vals[0][1])
        for i in range(1, len(ts)):
            cur = complex(vals[i][0], vals[i][1])
            d = cmath.phase(cur / prev)
            worst = max(worst, abs(d))
            th += d
            prev = cur
            for T in (20, 40, 60, 80, 100):
                if abs(ts[i] - T) < 1e-9:
                    phase[T] = th
        for T in (20, 40, 60, 80, 100):
            lhs = 2 * phase[T]                    # theta_c(T) - theta_c(-T)
            n = 2 * Nplus(T)                      # N(T), zeros with |gamma| <= T (N_{>c} = 0 for zeta, all known zeros on the line)
            err = lhs - cmath.pi * n
            log('c = %g T = %d: theta_c(T) - theta_c(-T) = %.4f, pi N(T) = %.4f (N = %d), error %.4f, error/log T = %.4f (max phase step %.3f)' % (c, T, lhs, cmath.pi * n, n, err, err / cmath.log(T).real, worst))
            out.append({'c': c, 'T': T, 'lhs': lhs, 'piN': cmath.pi * n, 'N': n, 'err': err, 'err_over_logT': err / cmath.log(T).real, 'max_step': worst})
    R['density'] = out


def task_arith(log, R, quick):
    mp.mp.dps = 30
    x = mp.mpf('3e12')
    L = mp.log(x / (2 * mp.pi))
    log('L = log(3e12/(2 pi)) =', mp.nstr(L, 6), ' 2/L =', mp.nstr(2 / L, 5))
    rows = []
    for eta in ('0.05', '0.0737', '0.1', '0.2', '0.3', '0.5'):
        e = mp.mpf(eta)
        ratio = (e / L) / (e * e / 2)
        log('eta = %s: eta/L = %s, eta^2/2 = %s, (eta/L)/(eta^2/2) = 2/(L eta) = %s; heat time below the local value: %s' % (eta, mp.nstr(e / L, 5), mp.nstr(e * e / 2, 5), mp.nstr(ratio, 5), bool(e / L < e * e / 2)))
        rows.append({'eta': eta, 'eta_over_L': mp.nstr(e / L, 5), 'half_eta2': mp.nstr(e * e / 2, 5), 'ratio': mp.nstr(ratio, 5)})
    R['arith'] = {'L': mp.nstr(L, 6), 'two_over_L': mp.nstr(2 / L, 5), 'rows': rows}


TASKS = ['path', 'counts', 'density', 'arith']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    base = a.out or ('data/mi_ladder' if not a.quick else 'misc/reviews/nick-miyatake/quick/G4/mi_ladder')
    os.makedirs(os.path.dirname(base) or '.', exist_ok=True)
    tasks = TASKS if a.task == 'all' else a.task.split(',')
    for t in tasks:
        log = Log(base + '.log')
        R = {}
        t0 = time.time()
        log('== task', t, 'quick' if a.quick else 'full', time.strftime('%Y-%m-%d %H:%M'))
        globals()['task_' + t](log, R, a.quick)
        log('== task', t, 'done in %.1f s' % (time.time() - t0))
        json.dump(R, open(base + '_' + t + '.json', 'w'), indent=1, default=str)


if __name__ == '__main__':
    main()
