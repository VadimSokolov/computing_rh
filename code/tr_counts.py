# tr_counts.py: the zero counts of Chapter 28, Sections "Transient zeros" and "Symmetrised truncations": the table of left and right zero counts
# of the hard cutoff and mean corrected truncations (argument principle on 0.02 < Re z < 0.2499 and 0.2501 < Re z < 0.48), the symmetrised
# hard cutoff table (zeros in the strip, on the line, left zeros of g_n), the first zeros of F_8, the modulus inequality minimum of g_8, and the
# mean corrected F_64, F_32 below height 40.
# Conventions. lambda_k = 2 pi^2 (k - 1/2)^2, A_k = prod_{j<=n, j!=k} lambda_j/(lambda_j - lambda_k) (an exact rational), M_n^(0)(z) = Gamma(1-z)
# sum_k A_k lambda_k^z, M_n^(1)(z) = sum_k A_k lambda_k^z e^(lambda_k c) Gamma(1-z, lambda_k c), c = zeta(2, n+1/2)/(2 pi^2) so that
# lambda_k c = (k - 1/2)^2 zeta(2, n + 1/2). g_n(z) = (2 pi)^(-z) M_n(z), F_n(z) = g_n(z) + g_n(1/2 - z). "Heights" are Im z (gamma = 2 Im z).
# Method. Zeros in a rectangle by the argument principle (python-flint acb at a precision of PREC bits, adaptive subdivision until the argument
# change of consecutive points is below 1 rad; the sum over the boundary divided by 2 pi must be within 0.05 of an integer, logged);
# zeros on the line by sign changes of F_n(1/4 + i y) (real) at step 0.005. Every count is recomputed at 1.5 times the precision and the two
# counts are compared (logged). The hard cutoff factor Gamma(1-z) has no zeros and is omitted for g_n counts, included for F_n.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Output: data/tr_counts.json, log data/tr_counts.log. Usage: python3 tr_counts.py [--quick] [--out FILE] [--cpus N]; --quick: n = 8, 16 only.
# Full run on Hopper: bash misc/tools/hopper_run.sh -c 16 -t 240 -g tr_counts.json -g tr_counts.log code/tr_counts.py -- --out tr_counts.json --cpus 16
import argparse, json, math, multiprocessing as mpr, sys, time
from flint import acb, arb, ctx, fmpq

def setup(n, kind, prec):
    ctx.prec = prec
    pi = arb.pi()
    lam = [2 * pi * pi * (arb(2 * k - 1) / 2) ** 2 for k in range(1, n + 1)]
    A = []
    for k in range(1, n + 1):
        q = fmpq(1)
        for j in range(1, n + 1):
            if j != k:
                q *= fmpq((2 * j - 1) ** 2, (2 * j - 1) ** 2 - (2 * k - 1) ** 2)
        A.append(arb(q.p) / arb(q.q))
    loglam = [l.log() for l in lam]
    zz = acb(2).zeta(acb(n) + acb(1) / 2) if kind == 'mean' else None
    lc = [arb(((2 * k - 1) ** 2)) / 4 * zz.real for k in range(1, n + 1)] if kind == 'mean' else None
    return A, loglam, lc

def make_M(n, kind, prec):
    A, loglam, lc = setup(n, kind, prec)
    def M(z):
        s = acb(0)
        if kind == 'hard':
            for a, ll in zip(A, loglam):
                s += a * (z * ll).exp()
            return s
        one = acb(1) - z
        for a, ll, x in zip(A, loglam, lc):
            s += a * (z * ll).exp() * x.exp() * acb(x).gamma_upper(one)
        return s
    return M

def winding(f, corners, step):
    tot = 0.0; worst = 0.0
    def seg(p0, f0, p1, f1, d):
        r = (f1 / f0).arg()
        r = float(r)
        if abs(r) < 1.0 or d > 40:
            return r
        pm = (p0 + p1) / 2; fm = f(pm)
        return seg(p0, f0, pm, fm, d + 1) + seg(pm, fm, p1, f1, d + 1)
    for c0, c1 in zip(corners, corners[1:] + corners[:1]):
        L = abs(complex(c1) - complex(c0)); m = max(1, int(L / step))
        pts = [c0 + (c1 - c0) * acb(i) / m for i in range(m + 1)]
        fs = [f(p) for p in pts]
        for i in range(m):
            tot += seg(pts[i], fs[i], pts[i + 1], fs[i + 1], 0)
    return tot / (2 * math.pi)

def rect(f, x1, x2, y1, y2, step=0.05):
    corners = [acb(x1, y1), acb(x2, y1), acb(x2, y2), acb(x1, y2)]
    w = winding(f, corners, step)
    return round(w), abs(w - round(w))

def count_job(job):
    kind, n, y1, y2, side, prec = job
    outs = []
    for p in (prec, prec * 3 // 2):
        M = make_M(n, kind, p)
        x1, x2 = (0.02, 0.2499) if side == 'L' else (0.2501, 0.48)
        c, err = rect(M, x1, x2, y1, y2)
        outs.append((c, err))
    return job, outs

def F_factory(n, kind, prec):
    M = make_M(n, kind, prec)
    pi = arb.pi()
    def g(z):
        v = M(z) * (-(z) * (2 * pi).log()).exp()
        if kind == 'hard': v = v * (acb(1) - z).gamma()
        return v
    def F(z):
        return g(z) + g(acb(1) / 2 - z)
    return g, F

def sym_job(job):
    kind, n, y1, y2, prec = job
    res = {}
    for p in (prec, prec * 3 // 2):
        g, F = F_factory(n, kind, p)
        strip, e1 = rect(F, 0.02, 0.48, y1, y2, 0.05)
        left, e2 = rect(g, 0.02, 0.2499, y1, y2) if kind == 'hard' else (None, 0)
        # sign changes on the line
        N = int((y2 - y1) / 0.005); prev = None; sc = 0; roots = []
        for i in range(N + 1):
            y = y1 + i * 0.005
            v = float(F(acb(0.25, y)).real)
            if prev is not None and prev * v < 0:
                sc += 1
                if len(roots) < 12: roots.append(round(2 * (y - 0.005 * v / (v - prev)), 4))
            prev = v
        res[p] = (strip, left, sc, roots, max(e1, e2))
    return job, res

def modulus_job(prec):
    g, F = F_factory(8, 'hard', prec)
    best = (1e9, None)
    for i in range(0, 25):
        x = 0.01 + 0.01 * i
        for j in range(0, 1200):
            y = 0.5 + 0.05 * j
            z = acb(x, y)
            d = float(g(z).abs_upper().log()) - float(g(acb(0.5) - z).abs_upper().log())
            if d < best[0]: best = (d, (x, y))
    return best

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None)
    ap.add_argument('--cpus', type=int, default=8); ap.add_argument('--modulus-only', action='store_true')
    a = ap.parse_args()
    out = a.out or ('/tmp/tr_counts.json' if a.quick else 'data/tr_counts.json')
    logf = open(out.replace('.json', '.log'), 'w')
    def log(s): print(s); logf.write(s + '\n'); logf.flush()
    t0 = time.time(); pool = mpr.Pool(a.cpus); R = {}
    if a.modulus_only:
        b = modulus_job(400 + 24 * 8 * 2); log(f"modulus inequality g_8: min over x in [0.01,0.25], y in [0.5,60.5] of log|g_8(z)| - log|g_8(1/2-z)| = {b[0]:.4f} at {b[1]}"); return
    P = lambda n: 400 + 12 * n * 2
    if a.quick:
        rows = [('hard', 8, 5, 10), ('mean', 8, 5, 10)]
        cj = [(k, n, y1, y2, s, P(n)) for k, n, y1, y2 in rows for s in 'LR']
        sj = [('hard', 8, 0.5, 15, P(8))]
    else:
        rows = [('hard', 32, 15, 60), ('hard', 64, 15, 60), ('mean', 32, 0.5, 20), ('mean', 16, 20, 40), ('mean', 32, 20, 40),
                ('mean', 64, 20, 40), ('mean', 64, 40, 75), ('mean', 128, 20, 40), ('hard', 32, 7.5, 30), ('hard', 64, 7.5, 30), ('hard', 32, 30, 60)]
        cj = [(k, n, y1, y2, s, P(n)) for k, n, y1, y2 in rows for s in 'LR']
        sj = [('hard', n, 0.5, 60, P(n)) for n in (8, 16, 32, 64)] + [('mean', 64, 0.5, 40, P(64)), ('mean', 32, 0.5, 40, P(32))]
    for job, outs in pool.imap_unordered(count_job, cj):
        k, n, y1, y2, s, p = job
        log(f"count {k} n={n} Im z in ({y1},{y2}) side {s}: {outs[0][0]} (arg residual {outs[0][1]:.3g}); at 1.5x precision {outs[1][0]} ({outs[1][1]:.3g}); {'agree' if outs[0][0] == outs[1][0] else 'DIFFER'}")
        R[f"{k}_{n}_{y1}_{y2}_{s}"] = outs
    for job, res in pool.imap_unordered(sym_job, sj):
        k, n, y1, y2, p = job
        (s1, l1, c1, r1, e1), (s2, l2, c2, r2, e2) = res[p], res[p * 3 // 2]
        log(f"sym {k} n={n} Im z in ({y1},{y2}): strip {s1}, on line {c1}, left zeros of g {l1}; at 1.5x: {s2}, {c2}, {l2}; first line zeros (gamma=2y) {r1[:8]}")
        R[f"sym_{k}_{n}"] = {str(p): res[p] for p in res}
    if not a.quick:
        b = pool.apply(modulus_job, (P(8),)); log(f"modulus inequality g_8: min over x in [0.01,0.25], y in [0.5,60.5] of log|g_8(z)| - log|g_8(1/2-z)| = {b[0]:.4f} at {b[1]}")
        R['modulus'] = b
    json.dump(R, open(out, 'w'), default=str); log(f"time {time.time() - t0:.1f} s")
if __name__ == '__main__':
    main()
