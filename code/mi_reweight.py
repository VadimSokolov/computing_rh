# mi_reweight.py: the reweighting experiment of Chapter 22 (Nick Polson's Miyatake chapter, prefix mi:), part G1: for P = 2, weights w_n on the
# supports U = {1}, {1,2}, {1,2,4} chosen to match the deficit at k heights, the number of real zeros tracked (7, 9, 10 against 6, one spurious each) and
# w_4 - 1 about 1.9e9. Added for the book on 22 Oct 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. Xi_w = Xi + D_w, D_w(t) = (t^2 + 1/4) sum_n (1 - w_n) C_(pi n^2)(t) with w_n = 0 outside U (C_a as in code/mi_friable.py), so
# "matching the deficit at t_i" is D_w(t_i) = 0, i.e. sum_(n in U) w_n C_n(t_i) = sum_(all n) C_n(t_i) (n up to 6, e^(-pi 36) negligible). The chapter does not say
# how the heights are chosen: variant A uses k = |U| heights t_i = 40 i / k (exact solve), variant B uses the same three heights 40/3, 80/3, 40 for
# every U (least squares when |U| < 3). Real zeros of Xi_w: sign changes on [0, 80] at step 0.1 (a pair closer than the step is missed); a zero is
# "tracked" if it lies within 0.5 of a zero of Xi, "spurious" otherwise; the tracked count is the number of distinct zeros of Xi within 0.5 of a zero of Xi_w
# (all of them below the first spurious or missing one are listed). Error control: weights and counts are recomputed at 10 more digits.
# Usage: python3 -I code/mi_reweight.py VARIANT SIZE [--out FILE], VARIANT in A B, SIZE in 1 2 3 (|U|). Each call below 30 s.
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mi_friable as F
from mpmath import mp, mpf, matrix, lu_solve, qr_solve, pi, fabs

def run(variant, size, dps, step=mpf('0.1'), tmax=80):
    mp.dps = dps
    U = [1, 2, 4][:size]
    k = size if variant == 'A' else 3
    heights = [mpf(40) * i / k for i in range(1, k + 1)]
    C = lambda n, t: F.Ca(pi * n * n, t)
    A = matrix(k, size); b = matrix(k, 1)
    for i, t in enumerate(heights):
        for j, n in enumerate(U): A[i, j] = C(n, t)
        b[i] = sum(C(n, t) for n in range(1, 7))
    w = lu_solve(A, b) if k == size else qr_solve(A, b)[0]
    w = [w[j] for j in range(size)]
    def Xiw(t):
        return F.Xi_hardy(t) + (t * t + mpf(1) / 4) * (sum(C(n, t) for n in range(1, 7)) - sum(w[j] * C(U[j], t) for j in range(size)))
    zl = []; prev = Xiw(mpf('0.05'))
    ts = [mpf('0.05') + i * step for i in range(0, int((tmax - 0.05) / step))]
    for t in ts[1:]:
        v = Xiw(t)
        if (v > 0) != (prev > 0): zl.append(t - step / 2)
        prev = v
    G = F.zeros(60)
    near = [min(range(len(G)), key=lambda i: fabs(G[i] - z)) for z in zl]
    trk = [i for z, i in zip(zl, near) if fabs(G[i] - z) < 0.5]
    spur = [z for z, i in zip(zl, near) if fabs(G[i] - z) >= 0.5]
    return {'variant': variant, 'support': U, 'heights': heights, 'weights': w, 'w_last - 1': w[-1] - 1, 'real zeros of Xi_w in (0,%d)' % tmax: len(zl),
            'tracked (distinct zeros of Xi within 0.5)': len(set(trk)), 'indices of tracked': sorted(set(trk)), 'spurious': spur, 'zeros': zl}

if __name__ == '__main__':
    v, size = sys.argv[1], int(sys.argv[2]); t0 = time.time()
    r1 = run(v, size, 30)
    mp.dps = 30
    r2 = run(v, size, 40) if '--nocheck' not in sys.argv else None
    out = {'result': r1, 'seconds': round(time.time() - t0, 1)}
    if r2:
        out['check'] = {'w_last-1 rel diff (30 vs 40 digits)': mp.nstr(fabs(r1['w_last - 1'] - r2['w_last - 1']) / fabs(r2['w_last - 1']), 3), 'zero counts equal': len(r1['zeros']) == len(r2['zeros'])}
    def ts(o):
        if isinstance(o, dict): return {str(k): ts(x) for k, x in o.items()}
        if isinstance(o, (list, tuple)): return [ts(x) for x in o]
        if isinstance(o, mpf): return mp.nstr(o, 10)
        return o
    txt = json.dumps(ts(out), indent=1)
    if '--out' in sys.argv: open(sys.argv[sys.argv.index('--out') + 1], 'w').write(txt + '\n')
    print(txt[:2500])
