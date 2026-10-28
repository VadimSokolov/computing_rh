# ev_prony2.py: Prony fits to the reciprocal Xi hitting time tail (E:516-532, tab:prony2) and the digits gained per added atom or pair (M5 of review A2).
# Claims supported: the four rows of tab:prony2 (K = 6, 8, 10, 10; x0 = 0.02, 0.02, 0.015, 0.02; step 0.003; errors of gamma_1..gamma_6, gamma_8) and
# "each extra pair of atoms buys roughly five more digits" (E:532).
# Conventions. X = sum E_k/gamma_k^2; the tail is taken in two forms, the survival function S(x) = sum c_k exp(-lam_k x) and the density
# f(x) = sum c_k lam_k exp(-lam_k x) (the chapter says "tail of the hitting time" and does not say which), lam_k = gamma_k^2,
# c_k = 2 i xi(1/2)/(gamma_k xi'(1/2 + i gamma_k)). Data y_j = S(x0 + j d), j < 2K, d = 0.003 (E:516).
# Method. The residue series is summed exactly (instead of Talbot inversion) over the first NZ ordinates, each refined by Newton on Xi from the 68
# digit file to the working precision; the series terms |c_k| exp(-lam_k x) are below 10^(-dps-5) beyond NZ. Prony: the linear prediction
# polynomial of degree K from the K x K Hankel system, its roots z_k = exp(-lam_k d), lam_k = -log z_k/d; the K roots are matched to the
# true lam_1..lam_K by sorting on the real part; error of gamma_k = |sqrt(lam_hat_k) - gamma_k|.
# Error control. Every row is run at two precisions (dps and dps + 30 with all zeros refined to the higher one); the digits of the errors that
# agree between the two runs are logged (the error of a row is trustworthy where the two runs agree to the printed digit).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 -I ev_prony2.py [--quick] [--out FILE] [--dps D]; quick = dps 70 (zeros not refined, 68 digit file), rows K = 6, 8 only; full = dps 150, 170.
import argparse, json, sys, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, mpc, exp, log, sqrt, nstr, matrix, lu_solve, polyroots
from ev_evt_lib import *

ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None); ap.add_argument('--dps', type=int, default=150)
args = ap.parse_args()
logf = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.append(s)

ROWS = [(6, '0.02'), (8, '0.02'), (10, '0.015'), (10, '0.02')]
CH = {  # chapter's tab:prony2 errors (gamma_1..6, gamma_8)
    (6, '0.02'): ['9e-14', '1e-9', '3e-7', '6e-4', '0.013', '0.52', None],
    (8, '0.02'): ['5e-25', '1e-19', '3e-16', '3e-11', '4e-9', '2e-5', '0.11'],
    (10, '0.015'): ['8e-34', '1e-27', '1e-23', '7e-18', '4e-15', '2e-10', '5e-5'],
    (10, '0.02'): ['2e-39', '7e-33', '2e-28', '5e-22', '6e-19', '2e-13', '4e-7']}

def prepare(dps, nz, refine):
    mp.dps = dps
    xh = xi(mpf(1) / 2)
    t = [mpf(z) for z in zero_strings(nz)]
    if refine:
        t = [refine_zero(g, 3) for g in t]
    c = [ck_from_zero(g, xh) for g in t]
    return t, c

def prony(y, K, d):
    A = matrix(K, K); b = matrix(K, 1)
    for i in range(K):
        for j in range(K):
            A[i, j] = y[i + j]
        b[i] = -y[i + K]
    p = lu_solve(A, b)
    coeffs = [mpf(1)] + [p[K - 1 - i] for i in range(K)]  # highest degree first: z^K + p_{K-1} z^{K-1} + ... + p_0
    roots = polyroots(coeffs, maxsteps=500, extraprec=4 * mp.prec)
    lam = sorted([(-log(z) / d) for z in roots], key=lambda v: v.real)
    return lam

def run(dps, nz, refine, rows, form):
    t, c = prepare(dps, nz, refine)
    lam = [g * g for g in t]
    d = mpf('0.003'); res = {}
    for K, x0 in rows:
        x0 = mpf(x0)
        if form == 'S':
            y = [sum(c[k] * exp(-lam[k] * (x0 + j * d)) for k in range(nz)) for j in range(2 * K)]
        else:
            y = [sum(c[k] * lam[k] * exp(-lam[k] * (x0 + j * d)) for k in range(nz)) for j in range(2 * K)]
        lh = prony(y, K, d)
        errs = [abs(sqrt(lh[k]) - t[k]) for k in range(K)]
        res[(K, str(x0)[:5])] = errs
    return res

t0 = time.time()
rows = ROWS[:2] if args.quick else ROWS
d1, d2 = (70, 90) if args.quick else (args.dps, args.dps + 30)
nz = 70 if args.quick else 100
out = {}
for form in (('S', 'f') if args.quick else ('S',)):  # the density form only in the quick run (time)
    A = run(d1, nz, not args.quick, rows, form)
    mp.dps = 20
    B = run(d2, nz + (0 if args.quick else 15), True, rows, form)
    mp.dps = 20
    for (K, x0) in rows:
        key = (K, {'0.02': '0.02', '0.015': '0.015'}[x0])
        ea = A[(K, x0)] if (K, x0) in A else A[(K, x0[:5])]
        eb = B[(K, x0)] if (K, x0) in B else B[(K, x0[:5])]
        idx = [0, 1, 2, 3, 4, 5, 7]
        mine = [nstr(ea[i], 2) if i < K else '-' for i in idx]
        agree = [nstr(abs(ea[i] - eb[i]) / ea[i], 2) if i < K else '-' for i in idx]
        P('form', form, 'K', K, 'x0', x0, 'errors g1..g6,g8 (%d digits)' % d1, mine)
        P('   chapter', CH[key])
        P('   relative change %d -> %d digits' % (d1, d2), agree)
        out['%s_K%d_x0%s' % (form, K, x0)] = {'errors': mine, 'chapter': CH[key], 'relchange': agree}
# digits per added atom and pair from the chapter's own table (M5)
import math
ch = {k: [float(x) for x in v if x] for k, v in CH.items()}
for g in range(3):
    e6, e8, e10 = ch[(6, '0.02')][g], ch[(8, '0.02')][g], ch[(10, '0.02')][g]
    P('chapter table, gamma_%d at x0 = 0.02: digits gained K 6->8 = %.1f, 8->10 = %.1f (per atom %.1f, %.1f)' % (g + 1, math.log10(e6 / e8), math.log10(e8 / e10), math.log10(e6 / e8) / 2, math.log10(e8 / e10) / 2))
P('time', round(time.time() - t0, 1))
of = args.out or ('/tmp/ev_prony2.json' if args.quick else os.path.join(HERE, '..', 'data', 'ev_prony2.json'))
json.dump(out, open(of, 'w'), indent=1)
open(of.replace('.json', '.log'), 'w').write('\n'.join(logf) + '\n')
