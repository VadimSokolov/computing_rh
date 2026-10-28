# fc_fisher.py: group FC, item 3. Fisher zeros of Knauf's chain near beta = 1 + rho_n, n = 1..20: the zeros of Z_N(beta) = sum c_N(m) m^(-beta)
# ("raw") and of the profile corrected C_N(beta) = Z_N(beta) - W_N(2-beta)/zeta(2) ("corr", Observation obs:corrected), N = 10, 12, ..., --nmax.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Method. The counts c_N(m) are the exact integers of fc_stern.py. C_N = A - B/(zeta(2)(2-beta)) with A = sum c m^-beta, B = sum d(m) m^(2-beta),
# d(m) = w(m-1) - w(m) (fc_lib.py), so C_N' = A' - (B'/(2-beta) + B/(2-beta)^2)/zeta(2). Double precision Newton from 1 + rho_n (68 digit ordinates
# of data/zeros_1700.txt, rounded to double). Enclosure (N <= --certmax): interval Newton in ball arithmetic (python-flint acb, 100 bits,
# exact rational weights d(m) m^2): X = square of half side r = 1e-8 around the double root c; the operator N(X) = c - F(c)/F'(X) with F(c) from
# the exact centre and F'(X) from the ball; if N(X) lies in the interior of X, X contains exactly one zero of F (Moore, Krawczyk). The
# check compares mid + rad of N(X) - c with r in double precision, with a margin of at least a factor 10 (logged as margin). Count check
# (N <= --argmax, uncertified, double): the winding number of F around the circle |beta - (1 + rho_n)| = 0.3 (240 points) must be 1.
# The predicted displacement is |zeta(1+rho_n)|/|zeta'(rho_n)| |W_N(1/2 - i gamma_n)|/zeta(2) (mpmath, 30 digits).
import sys, os, json, time, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fc_lib as L
from multiprocessing import Pool
import mpmath as mp

_PHI = {}
def get_phi(nmax):
    if nmax not in _PHI:
        _PHI[nmax] = L.phi_sieve(L.fib(nmax + 2) + 2)
    return _PHI[nmax]

class Prep:
    def __init__(self, c, phi):
        Mn = len(c) - 1
        d, w = L.profile_terms(c, phi)
        self.m1 = np.arange(1, Mn + 1, dtype=np.float64)
        self.lm1 = np.log(self.m1)
        self.cm = c[1:Mn + 1].astype(np.float64)
        m2 = np.arange(2, Mn + 2, dtype=np.float64)
        self.lm2 = np.log(m2)
        self.d2 = d * m2 * m2
        self.cml = self.cm * self.lm1
        self.d2l = self.d2 * self.lm2
        self.z2 = np.pi ** 2 / 6
    def AB(self, b):
        e1 = np.exp(-b * self.lm1); e2 = np.exp(-b * self.lm2)
        return np.dot(self.cm, e1), np.dot(self.d2, e2), -np.dot(self.cml, e1), -np.dot(self.d2l, e2)
    def F(self, b, typ):
        A, B, Ap, Bp = self.AB(b)
        if typ == 'raw':
            return A, Ap
        z = 2 - b
        return A - B / (self.z2 * z), Ap - (Bp / z + B / z ** 2) / self.z2

def newton(P, typ, b0):
    b = complex(b0)
    for it in range(60):
        f, fp = P.F(b, typ)
        step = f / fp
        b -= step
        if abs(step) < 1e-13:
            break
    f, fp = P.F(b, typ)
    return b, abs(f), it

def winding(P, typ, center, rad=0.3, npts=240):
    ph = [P.F(center + rad * np.exp(2j * np.pi * k / npts), typ)[0] for k in range(npts + 1)]
    ph = np.array(ph)
    dang = np.angle(ph[1:] / ph[:-1])
    return float(dang.sum() / (2 * np.pi)), float(np.max(np.abs(dang)))

def arb_tables(c, phi, prec=100):
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
    return lm, cA, cB

def arb_F(tab, b, typ, deriv):
    from flint import arb, acb
    lm, cA, cB = tab
    Mn = len(lm) - 2
    A = acb(0); B = acb(0); Ap = acb(0); Bp = acb(0)
    for m in range(1, Mn + 2):
        e = (-(b * lm[m])).exp()
        if cA[m]:
            t = cA[m] * e
            A += t
            if deriv:
                Ap -= t * lm[m]
        if m >= 2:
            t = cB[m] * e
            B += t
            if deriv:
                Bp -= t * lm[m]
    z2 = arb.pi() ** 2 / 6
    if typ == 'raw':
        return A, Ap
    z = 2 - b
    return A - B / (z2 * z), Ap - (Bp / z + B / (z * z)) / z2

def certify(tab, typ, c0, r=1e-8):
    from flint import arb, acb
    cc = acb(arb(c0.real), arb(c0.imag))
    X = acb(arb(c0.real, r), arb(c0.imag, r))
    Fc, _ = arb_F(tab, cc, typ, False)
    _, dX = arb_F(tab, X, typ, True)
    Nb = cc - Fc / dX
    x = Nb - cc
    bre = abs(float(x.real.mid())) + float(x.real.rad())
    bim = abs(float(x.imag.mid())) + float(x.imag.rad())
    return (bre < r and bim < r), max(bre, bim) / r

_CACHE = {}

def get_prep(N, births, nmaxg, certmax):
    key = (N, births)
    if key not in _CACHE:
        _CACHE.clear()
        z = np.load(births)
        Mn = L.fib(N + 2)
        c = np.zeros(Mn + 1, dtype=np.int64); c[1] = 1
        for j in range(1, N + 1):
            h = z['H%d' % j].astype(np.int64); c[:len(h)] += h
        phi = get_phi(nmaxg)
        P = Prep(c, phi)
        tab = arb_tables(c, phi) if N <= certmax else None
        _CACHE[key] = (P, tab)
    return _CACHE[key]

def task(args):
    """one zero n of one N: both chains"""
    N, n, births, zerofile, certmax, argmax, nmaxg = args
    t0 = time.time()
    P, tab = get_prep(N, births, nmaxg, certmax)
    zstr = L.zeros_list(zerofile, 20)
    mp.mp.dps = 30
    rho = mp.mpc(0.5, mp.mpf(zstr[n - 1]))
    g = float(mp.mpf(zstr[n - 1]))
    target = complex(1.5, g)
    rec = {'N': N, 'n': n, 'gamma': g}
    for typ in ('raw', 'corr'):
        b, res, it = newton(P, typ, target)
        rec[typ] = {'zero': [b.real, b.imag], 'dist': abs(b - target), 'resid': res, 'iters': it}
        if tab is not None:
            ok, margin = certify(tab, typ, b)
            rec[typ]['certified'] = bool(ok); rec[typ]['margin_ratio'] = margin
        if N <= argmax:
            wnd, mx = winding(P, typ, target)
            rec[typ]['winding_disc03'] = wnd; rec[typ]['max_phase_step'] = mx
    A, B, _, _ = P.AB(target)
    W = B / (2 - target)
    rec['pred'] = float(abs(mp.zeta(1 + rho)) / abs(mp.zeta(rho, derivative=1)) * abs(W) / P.z2)
    rec['seconds'] = time.time() - t0
    print('N=%d n=%d raw %.5f corr %.5f %.0f s' % (N, n, rec['raw']['dist'], rec['corr']['dist'], rec['seconds']), flush=True)
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--births', default='fc_births.npz')
    ap.add_argument('--zeros', default='zeros_1700.txt')
    ap.add_argument('--nmin', type=int, default=10)
    ap.add_argument('--nmax', type=int, default=34)
    ap.add_argument('--certmax', type=int, default=26)
    ap.add_argument('--argmax', type=int, default=28)
    ap.add_argument('--nproc', type=int, default=16)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='fc_fisher.json')
    args = ap.parse_args()
    if args.quick:
        args.nmax = 14; args.nproc = 2
    Ns = sorted(range(args.nmin, args.nmax + 1, 2), key=lambda N: -N)
    tasks = [(N, n, args.births, args.zeros, args.certmax, args.argmax, args.nmax) for N in Ns for n in range(1, 21)]
    with Pool(args.nproc) as pool:
        recs = pool.map(task, tasks, chunksize=1)
    res = []
    for N in sorted(Ns):
        zr = sorted([r for r in recs if r['N'] == N], key=lambda r: r['n'])
        out = {'N': N, 'zeros': zr}
        for typ in ('raw', 'corr'):
            d = np.array([r[typ]['dist'] for r in zr])
            out[typ + '_mean'] = float(d.mean()); out[typ + '_max'] = float(d.max()); out[typ + '_within_0.005'] = int((d <= 0.005).sum())
            out[typ + '_median'] = float(np.median(d))
        out['corr_better'] = int(sum(r['corr']['dist'] < r['raw']['dist'] for r in zr))
        print('N=%d raw mean %.4f max %.4f w005 %d | corr mean %.4f max %.4f w005 %d | corr better %d' % (
            N, out['raw_mean'], out['raw_max'], out['raw_within_0.005'], out['corr_mean'], out['corr_max'], out['corr_within_0.005'], out['corr_better']), flush=True)
        res.append(out)
    with open(args.out, 'w') as f:
        json.dump({'results': res}, f, indent=1)

if __name__ == '__main__':
    main()
