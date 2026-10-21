# zh_detect.py: the detection experiment of Section zh:sec:horizon of Chapter 18 (Zhu's compact windows, the Landau Widom horizon and the Thorin program), Table zh:tab:det (the smallest window half length a at which Yoshida's form on K(a) detects a hypothetical off line quadruple 1/2 +- delta +- i gamma_0 placed midway between two consecutive zeros) and Figure zh:fig:horizon (drawn here as fig/zh_horizon.pdf), with the sensitivity runs the reviewers asked for. Added for the book in October 2026 and written from the mathematics; the collaborator's scripts were neither read nor run.
# Model. The zero side matrix of Theorem yo:thm:lowner(i) on the basis chi_n = (2a)^(-1/2) e^(i pi n y/a), |n| <= N, nu_n = pi n/a: M_nm = (-1)^(n+m) (2/a) sum_rho sin^2(a gamma_rho)/((gamma_rho + nu_n)(gamma_rho + nu_m)), gamma_rho = (rho - 1/2)/i, over the zeros of both signs. The sign factor is a congruence by diag((-1)^n) and is dropped. Write u_n(z) = sqrt(2/a) sin(a z)/(z + nu_n); then M = sum over zeros of u(z) u(z)^T (bilinear, not sesquilinear), and u_n(-z) = u_(-n)(z), so the mirror symmetry n -> -n commutes with M and M splits exactly into an even sector (modes n = 0..N, orthonormal coordinates (u_n + u_(-n))/sqrt 2) and an odd sector (n = 1..N). The on line zeros give a positive semidefinite part A (the Gram matrix of the vectors of the zeros), positive definite as soon as the number of zero vectors is at least the sector size (a Cauchy matrix has all minors nonzero). The quadruple gamma_0 - i delta, gamma_0 + i delta, -gamma_0 - i delta, -gamma_0 + i delta adds 2 Re(u(z) u(z)^T) + 2 Re(u(-z) u(-z)^T), z = gamma_0 - i delta, which is of rank at most four; in each sector it is 4 p p^T - 4 q q^T with p + i q the sector coordinates of u(z), a rank two perturbation W D W^T, D = diag(4, -4). The chapter's model keeps every on line zero and adds the quadruple, which puts two more zeros in the gap (four with the mirror images) and so changes the zero count N(T); the variant mode replace removes the two on line zeros that flank gamma_0 (the count is then preserved, as in Theorem thm:ch13:perturb).
# Exact sign. With A positive definite and D = diag(d, -d), let G = W^T A^(-1) W (2 by 2, positive semidefinite) and S = D^(-1) + G. By Haynsworth inertia additivity the sector matrix has n_+(S) - 1 negative eigenvalues, and since S_11 = 1/d + G_11 > 0 it is indefinite if and only if det S > 0 (positive definite if and only if det S < 0). The same holds for A + mu I, so the event lambda_min < -mu is decided exactly for every mu >= 0: mu = 0 is criterion (ii) and mu = 1e-11 is the chapter's cut decided in exact arithmetic (key ex11). A is assembled in python-flint arb ball arithmetic from the partial fraction (Loewner) form: with t_j = gamma_j^2, tau_n = nu_n^2, S_j = sin^2(a gamma_j), R(tau) = sum_j S_j/(t_j - tau), the even sector is (4/a) kappa_n kappa_m (R_e(tau_n) - R_e(tau_m))/(tau_n - tau_m) with R_e = sum S_j + tau R, kappa_0 = 1, kappa_n = sqrt 2, and the odd sector is (8/a) nu_n nu_m (R(tau_n) - R(tau_m))/(tau_n - tau_m), the diagonals being the derivatives. The linear solve is certified by arb_mat.solve; the precision starts at 192 bits and doubles (to 4096) while a sign is not decided, so a decision is a certificate for the model with the stored zero ordinates. For the heights 1000 and 2000 the exact test uses the band |nu_n - gamma_0| <= B (default 150, the mirror band comes with the sector); by Cauchy interlacing a negative eigenvalue of the band submatrix is a negative eigenvalue of the full model, so band detection is certified and band absence is not (the runs with B = 75, 300 and 600 test how fast the band threshold falls with B). At the heights 1000 and 2000 the exact sign of the full model is therefore only bounded above, by the smaller of the band test and the certified float64 test f64c.
# Float64 criterion (key f64): the chapter's computation. The sector matrices are built as Gram matrices V V^T in float64 (the weights sin^2(a gamma_j) come from arb, the differences gamma_j - nu_n from a double double representation so that near resonances they keep full relative accuracy), the pair term is added, and the smallest eigenvalue of each sector is computed by LAPACK; the event is lambda_min < -1e-11. A rigorous looking bound is reported next to it: Weyl's inequality with the componentwise rounding bound of the Gram product (gamma_K times the one norm of |V||V|^T) plus the representation error of the pair term and of the sum, and n times the unit roundoff times the one norm of the matrix for the eigensolver (the textbook constant, LAPACK has no tabulated one); a computed value below -1e-11 minus the bound is a certified float64 detection (the key f64c is the event lambda_min < -bound, a rigorous detection of the sign change by the float64 eigenvalue alone). The bound is of order 1e-11 for sectors of size 2000, so near the cut the exact keys decide.
# Scan. For each (gamma_0, delta) the half length a runs over a grid of step 0.05 from a low value until all three keys have fired and two more grid points have been scanned, then each key is refined by bisection to width 0.005 between the last grid point without and the first with detection; the pattern of the three keys along the grid is stored, so a lack of monotonicity is visible. N(a, gamma_0) = ceil(a (gamma_0 + extra)/pi) + 2 with extra = 45 (90 in the sensitivity run), M = 2100 zeros (4000 and 8000 in the sensitivity runs; more zeros can only raise the eigenvalues, so the thresholds can only move up), gamma_0 the midpoint of the gap that rounds to 51.37, 199.64, 1000.57 and 1999.99, delta in 0.001, 0.003, 0.01, 0.05, 0.2, 0.45. The comparison scale x_*/2 of Remark yo:rem:scale is computed by the xstar task from D(x) = 8 sum_j sin^2(gamma_j x/2)/gamma_j^2 + 16 Re[sin^2(z x/2)/z^2], z = gamma_0 - i delta, in two readings: x_env, the first x at which the pair term can push D negative for some phase of its fast oscillation (period 2 pi/gamma_0 in x; the envelope minimum is exact and grid independent), and x_first, the first negative value on a grid finer than that period (the printed 32.9 of the remark for gamma = 1000, delta = 0.3 is an aliasing artefact of a step 0.01 grid; the resolved value is 32.41).
# Checks. selftest compares the sector matrices of the Loewner construction with the direct sums over all zeros and the explicit quadruple (complex terms included), the rank two form of the pair term, the exact decisions with mpmath eigenvalues of the same matrix (shifts mu around the smallest eigenvalue), the float64 eigenvalue with the mpmath one, and the zero side formula of Theorem yo:thm:lowner(i) with the prime side Weil explicit formula (yo:EF) for five entries at a = 1 (agreement up to the zero tail, about 1e-5 once the mean tail is added).
# Output. Every task writes out/task_KKK.json (values as strings; parameters, zero counts, precision of each decision, the criterion that produced each threshold, the scan pattern) and out/task_KKK.log; --gather merges them into data/zh_detect.json and data/zh_detect.log and prints the comparison with the chapter's table; --figure draws fig/zh_horizon.pdf from data/zh_detect.json. The zeros come from data/zeros_8000.txt (code/zh_zeros.py, run on Hopper) or, failing that, data/zeros_1700.txt (then M is at most 1700 and the run says so).
# Usage: python3 zh_detect.py --quick [--delta D] [--out FILE] (gamma_0 = 51.37, a = 0.9 .. 1.4, about 25 s, writes to /tmp unless --out is given); python3 zh_detect.py --task K [--outdir DIR] (one entry of the table printed by --list-tasks, meant for the sbatch array of misc/reviews/nick-zhu/z2-run); --count-tasks; --selftest; --gather [--tasks-dir DIR] [--out FILE]; --figure [--data FILE] [--fig FILE]. python-flint is needed for everything except --list-tasks, --count-tasks, --gather and --figure.
import argparse
import glob
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

try:
    from flint import arb, acb, arb_mat, ctx
except ImportError:
    arb = acb = arb_mat = ctx = None

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
T0 = time.time()
LOGFH = None
U = 2.0 ** -53            # unit roundoff of float64
MU_CUT = 1e-11            # the chapter's cut
KEYS = ('f64', 'f64c', 'ex0', 'ex11')
KEYNAME = {'f64': 'float64 eigenvalue below -1e-11 (the chapter)', 'f64c': 'float64 eigenvalue below minus the Weyl bound (certified lambda_min < 0)',
           'ex0': 'exact, lambda_min < 0', 'ex11': 'exact, lambda_min < -1e-11'}
TARGETS = ['51.37', '199.64', '1000.57', '1999.99']
DELTAS = [0.001, 0.003, 0.01, 0.05, 0.2, 0.45]
CHAPTER = {   # Table zh:tab:det as printed in the chapter
    '51.37': {0.01: 1.20, 0.05: 1.15, 0.2: 1.15, 0.45: 1.10},
    '199.64': {0.01: 1.85, 0.05: 1.80, 0.2: 1.75, 0.45: 1.75},
    '1000.57': {0.01: 2.65, 0.05: 2.65, 0.2: 2.60, 0.45: 2.50},
    '1999.99': {0.01: 3.05, 0.05: 3.00, 0.2: 2.95, 0.45: 2.45},
}
SCAN = {'51.37': (0.50, 2.00), '199.64': (1.00, 2.40), '1000.57': (1.50, 3.30), '1999.99': (1.50, 3.40)}
STEP = 0.05
TOL = 0.005
KAPPA1 = 0.0231049931


def fs(x, n=12):
    return '%.*g' % (n, x)


def log(*a):
    s = '[%7.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(s, flush=True)
    if LOGFH is not None:
        LOGFH.write(s + '\n')
        LOGFH.flush()


def need_flint():
    if arb is None:
        sys.exit('zh_detect.py needs python-flint for this mode (the book venv on Hopper has it; locally: pip install python-flint)')


# ---------------------------------------------------------------------------------------------------------------------------
# task table
# ---------------------------------------------------------------------------------------------------------------------------
def task_table():
    tasks = []
    for gi in range(4):
        for d in DELTAS:
            tasks.append(dict(kind='scan', variant='main', gi=gi, delta=d, M=2100, extra=45, band=(None if gi < 2 else 150.0), mode='keep'))
    for variant in ('M4000', 'M8000', 'x90', 'replace', 'band75', 'band150', 'band300', 'band600'):
        for gi in range(4):
            if variant == 'band150' and gi >= 2:
                continue            # the main run already is the band 150 run at these heights
            if variant == 'band600' and gi < 2:
                continue            # at the two lower heights the band 300 already covers every mode
            for d in (0.01, 0.2):
                t = dict(kind='scan', variant=variant, gi=gi, delta=d, M=2100, extra=45, band=(None if gi < 2 else 150.0), mode='keep')
                if variant == 'M4000':
                    t['M'] = 4000
                elif variant == 'M8000':
                    t['M'] = 8000
                elif variant == 'x90':
                    t['extra'] = 90
                elif variant == 'replace':
                    t['mode'] = 'replace'
                elif variant.startswith('band'):
                    t['band'] = float(variant[4:])
                tasks.append(t)
    tasks.append(dict(kind='xstar', variant='xstar'))
    for t in tasks:
        if t['kind'] == 'scan':
            lo, hi = SCAN[TARGETS[t['gi']]]
            if t['delta'] <= 0.003:
                hi += 0.5
            if t['variant'] == 'replace':
                lo = 0.2 if t['gi'] < 2 else 0.5       # removing the two flanking zeros makes the form negative at much smaller a
            t['scan_lo'], t['scan_hi'] = lo, hi
    return tasks


# ---------------------------------------------------------------------------------------------------------------------------
# zeros and parameters
# ---------------------------------------------------------------------------------------------------------------------------
def find_file(names):
    for n in names:
        for base in (Path.cwd(), HERE, ROOT, ROOT / 'data'):
            p = base / n
            if p.is_file():
                return p
    return None


def load_zero_strings(need, path=None):
    p = Path(path) if path else find_file(['data/zeros_8000.txt', 'zeros_8000.txt'])
    if p is None or not p.is_file():
        p = find_file(['data/zeros_1700.txt', 'zeros_1700.txt'])
    if p is None:
        sys.exit('no zeros file (data/zeros_8000.txt or data/zeros_1700.txt)')
    lines = [s.strip() for s in open(p) if s.strip()]
    if len(lines) < need:
        log('WARNING: %s has %d zeros, %d requested; using %d' % (p, len(lines), need, len(lines)))
    return lines[:need], str(p)


def arb_zeros(P, prec):
    cache = P['_arbz']
    if prec not in cache:
        ctx.prec = prec
        cache[prec] = [arb(s) for s in P['strings']]
    return cache[prec]


def make_params(spec, zstr_all, source):
    need_flint()
    M = min(spec['M'], len(zstr_all))
    strings = zstr_all[:M]
    fl = np.array([float(s) for s in strings])
    mids = 0.5 * (fl[:-1] + fl[1:])
    gi = spec['gi']
    j0 = int(np.argmin(np.abs(mids - float(TARGETS[gi]))))
    if abs(mids[j0] - float(TARGETS[gi])) > 0.006:
        sys.exit('no gap midpoint near %s among the first %d zeros' % (TARGETS[gi], M))
    use = np.ones(M, dtype=bool)
    if spec['mode'] == 'replace':
        use[j0] = False
        use[j0 + 1] = False
    ctx.prec = 300
    lo = np.array([float(arb(s) - arb(float(s))) for s in strings])
    P = dict(spec=spec, M=M, strings=strings, source=source, delta=spec['delta'], extra=spec['extra'], band=spec['band'],
             mode=spec['mode'], j0=j0, label=TARGETS[gi], use=use, use_idx=[j for j in range(M) if use[j]],
             Zhi=fl[use], Zlo=lo[use], g0f=float(0.5 * (fl[j0] + fl[j0 + 1])), nvec=int(use.sum()),
             _cache={}, _arbz={}, _pstart=192)
    return P


def nodes_for(P, a, N):
    if P['band'] is None:
        return list(range(0, N + 1)), list(range(1, N + 1))
    J = [n for n in range(1, N + 1) if abs(math.pi * n / a - P['g0f']) <= P['band']]
    return J, J


# ---------------------------------------------------------------------------------------------------------------------------
# the pair term (arb)
# ---------------------------------------------------------------------------------------------------------------------------
def pair_vectors(P, a, p, nodes_even, nodes_odd):
    """Real and imaginary parts of the sector coordinates of u(z), z = gamma_0 - i delta, as arb lists."""
    ctx.prec = p
    a_arb = arb(a)
    pi = arb.pi()
    g0 = (arb(P['strings'][P['j0']]) + arb(P['strings'][P['j0'] + 1])) / 2
    z = acb(g0, -arb(str(P['delta'])))
    sz = (a_arb * z).sin()
    c = (2 / a_arb).sqrt()
    sq2 = arb(2).sqrt()
    z2 = z * z
    pe, qe, po, qo = [], [], [], []
    for n in nodes_even:
        if n == 0:
            u = c * sz / z
        else:
            nu = pi * n / a_arb
            u = c * sz * sq2 * z / (z2 - nu * nu)
        pe.append(u.real)
        qe.append(u.imag)
    for n in nodes_odd:
        nu = pi * n / a_arb
        u = -(c * sz * sq2 * nu) / (z2 - nu * nu)
        po.append(u.real)
        qo.append(u.imag)
    return dict(even=(pe, qe), odd=(po, qo))


# ---------------------------------------------------------------------------------------------------------------------------
# exact sign test (arb)
# ---------------------------------------------------------------------------------------------------------------------------
def exact_build(P, a, N, prec):
    """Sector matrices of the on line zeros (arb, Loewner form) for the nodes of the scope (full or band)."""
    ctx.prec = prec
    a_arb = arb(a)
    pi = arb.pi()
    Z = arb_zeros(P, prec)
    idx = P['use_idx']
    even_nodes, odd_nodes = nodes_for(P, a, N)
    allnodes = sorted(set(even_nodes) | set(odd_nodes))
    nu = {n: pi * n / a_arb for n in allnodes}
    tau = {n: nu[n] * nu[n] for n in allnodes}
    t = [Z[j] * Z[j] for j in idx]
    S = []
    for j in idx:
        s = (a_arb * Z[j]).sin()
        S.append(s * s)
    SS = arb(0)
    for s in S:
        SS += s
    Ro, Ro1 = {}, {}
    for n in allnodes:
        tn = tau[n]
        r0 = arb(0)
        r1 = arb(0)
        for tj, sj in zip(t, S):
            r = 1 / (tj - tn)
            w = sj * r
            r0 += w
            r1 += w * r
        Ro[n] = r0
        Ro1[n] = r1
    c4 = 4 / a_arb
    c8 = 8 / a_arb
    kap = {n: (arb(1) if n == 0 else arb(2).sqrt()) for n in even_nodes}
    Re = {n: SS + tau[n] * Ro[n] for n in even_nodes}
    Re1 = {n: Ro[n] + tau[n] * Ro1[n] for n in even_nodes}
    ent_e = []
    for n in even_nodes:
        for m in even_nodes:
            if n == m:
                ent_e.append(c4 * kap[n] * kap[n] * Re1[n])
            else:
                ent_e.append(c4 * kap[n] * kap[m] * (Re[n] - Re[m]) / (tau[n] - tau[m]))
    ent_o = []
    for n in odd_nodes:
        for m in odd_nodes:
            if n == m:
                ent_o.append(c8 * tau[n] * Ro1[n])
            else:
                ent_o.append(c8 * nu[n] * nu[m] * (Ro[n] - Ro[m]) / (tau[n] - tau[m]))
    return dict(even=(even_nodes, ent_e), odd=(odd_nodes, ent_o))


def decide(entries, n, W, mu, p):
    """Exact sign test of A + mu I + W D W^T, D = diag(4, -4): ('neg' | 'pos' | 'undecided', det S as text, precision)."""
    ctx.prec = p
    ent = list(entries)
    if mu:
        m = arb(str(mu))
        for i in range(n):
            ent[i * (n + 1)] = ent[i * (n + 1)] + m
    A = arb_mat(n, n, ent)
    try:
        X = A.solve(W)
    except ZeroDivisionError:
        return ('undecided', 'solve failed', p)
    G = W.transpose() * X
    s11 = arb(1) / 4 + G[0, 0]
    s22 = G[1, 1] - arb(1) / 4
    s12 = (G[0, 1] + G[1, 0]) / 2
    det = s11 * s22 - s12 * s12
    if det > 0:
        return ('neg', det.str(8, radius=False), p)
    if det < 0:
        return ('pos', det.str(8, radius=False), p)
    return ('undecided', det.str(6), p)


def make_W(pv, sec):
    pl, ql = pv[sec]
    n = len(pl)
    flat = []
    for i in range(n):
        flat.append(pl[i])
        flat.append(ql[i])
    return arb_mat(n, 2, flat)


def exact_eval(P, a, N, mus=(0.0, MU_CUT)):
    even_nodes, odd_nodes = nodes_for(P, a, N)
    if len(even_nodes) == 0 or len(odd_nodes) == 0:
        return None
    state = {(mu, sec): ('undecided', '', 0) for mu in mus for sec in ('even', 'odd')}

    def done(mu):
        st = [state[(mu, 'even')][0], state[(mu, 'odd')][0]]
        return ('neg' in st) or (st == ['pos', 'pos'])
    p = P['_pstart']
    while True:
        B = exact_build(P, a, N, p)
        pv = pair_vectors(P, a, p, B['even'][0], B['odd'][0])
        Ws = {sec: make_W(pv, sec) for sec in ('even', 'odd')}
        for sec in ('even', 'odd'):
            nodes, ent = B[sec]
            for mu in mus:
                if done(mu) or state[(mu, sec)][0] != 'undecided':
                    continue
                state[(mu, sec)] = decide(ent, len(nodes), Ws[sec], mu, p)
        if all(done(mu) for mu in mus) or p >= 4096:
            break
        p *= 2
    P['_pstart'] = max(128, min(p, 2048))
    out = {}
    for mu in mus:
        st = {sec: state[(mu, sec)] for sec in ('even', 'odd')}
        decs = [st['even'][0], st['odd'][0]]
        dec = 'neg' if 'neg' in decs else ('pos' if decs == ['pos', 'pos'] else 'undecided')
        out[fs(mu, 3)] = dict(dec=dec, prec=max(st['even'][2], st['odd'][2]), even=[st['even'][0], st['even'][1]],
                              odd=[st['odd'][0], st['odd'][1]])
    out['sizes'] = [len(even_nodes), len(odd_nodes)]
    return out


# ---------------------------------------------------------------------------------------------------------------------------
# float64 criterion (the chapter's computation)
# ---------------------------------------------------------------------------------------------------------------------------
def f64_eval(P, a, N):
    from scipy.linalg import eigh
    ctx.prec = 256
    a_arb = arb(a)
    pi = arb.pi()
    nus = [pi * n / a_arb for n in range(N + 1)]
    nu_hi = np.array([float(x) for x in nus])
    nu_lo = np.array([float(x - arb(float(x))) for x in nus])
    Z = arb_zeros(P, 256)
    S = np.empty(P['nvec'])
    for k, j in enumerate(P['use_idx']):
        s = (a_arb * Z[j]).sin()
        S[k] = float(s * s)
    gh, gl = P['Zhi'], P['Zlo']
    t = gh * gh
    d1 = (gh[None, :] - nu_hi[:, None]) + (gl[None, :] - nu_lo[:, None])
    Dn = d1 * (gh[None, :] + nu_hi[:, None])
    del d1
    kap = np.full(N + 1, math.sqrt(2.0))
    kap[0] = 1.0
    Ve = (kap[:, None] * np.sqrt((4.0 / a) * S * t)[None, :]) / Dn
    Vo = (nu_hi[1:, None] * np.sqrt((8.0 / a) * S)[None, :]) / Dn[1:]
    del Dn
    pv = pair_vectors(P, a, 256, list(range(0, N + 1)), list(range(1, N + 1)))
    lam, bnd = {}, {}
    for sec, V in (('even', Ve), ('odd', Vo)):
        pvec = np.array([float(x) for x in pv[sec][0]])
        qvec = np.array([float(x) for x in pv[sec][1]])
        A = V @ V.T
        Pm = 4.0 * np.outer(pvec, pvec) - 4.0 * np.outer(qvec, qvec)
        Mf = A + Pm
        lam[sec] = float(eigh(Mf, eigvals_only=True, subset_by_index=[0, 0])[0])
        absV = np.abs(V)
        N1 = float(np.max(absV @ absV.sum(axis=0)))
        K = V.shape[1]
        gK = (K + 8) * U / (1.0 - (K + 8) * U)
        b_con = (gK + 17 * U) * N1
        b_pair = 16 * U * (float(pvec @ pvec) + float(qvec @ qvec))
        b_add = U * (float(np.abs(A).sum(axis=0).max()) + float(np.abs(Pm).sum(axis=0).max()))
        b_eig = A.shape[0] * U * float(np.abs(Mf).sum(axis=0).max())
        bnd[sec] = b_con + b_pair + b_add + b_eig
        del A, Pm, Mf, absV
    lam_min = min(lam.values())
    bound = bnd['even'] if lam['even'] <= lam['odd'] else bnd['odd']
    return dict(lam_even=fs(lam['even']), lam_odd=fs(lam['odd']), lam_min=fs(lam_min), bound=fs(bound, 4),
                flag=bool(lam_min < -MU_CUT), certified=bool(lam_min < -MU_CUT - bound),
                negative_certified=bool(lam_min < -bound))


# ---------------------------------------------------------------------------------------------------------------------------
# evaluation, scan, bisection
# ---------------------------------------------------------------------------------------------------------------------------
def evaluate(P, a):
    key = round(a, 9)
    if key in P['_cache']:
        return P['_cache'][key]
    t_ev = time.time()
    N = int(math.ceil(a * (P['g0f'] + P['extra']) / math.pi)) + 2
    rec = {'a': fs(a, 10), 'N': N}
    if N + 1 > P['nvec']:
        rec['rank_guard'] = True
    rec['f64'] = f64_eval(P, a, N)
    ex = None if rec.get('rank_guard') else exact_eval(P, a, N)
    rec['exact'] = ex
    flags = {'f64': rec['f64']['flag'], 'f64c': rec['f64']['negative_certified'], 'ex0': False, 'ex11': False}
    und = {'ex0': False, 'ex11': False}
    if ex is not None:
        flags['ex0'] = ex[fs(0.0, 3)]['dec'] == 'neg'
        flags['ex11'] = ex[fs(MU_CUT, 3)]['dec'] == 'neg'
        und['ex0'] = ex[fs(0.0, 3)]['dec'] == 'undecided'
        und['ex11'] = ex[fs(MU_CUT, 3)]['dec'] == 'undecided'
    rec['flags'] = flags
    rec['undecided'] = und
    rec['seconds'] = round(time.time() - t_ev, 2)
    P['_cache'][key] = rec
    log('a=%s N=%d f64 lam=%s bound=%s flag=%d cert=%d | ex0=%s ex11=%s prec=%s | %.1fs' % (
        rec['a'], N, rec['f64']['lam_min'], rec['f64']['bound'], rec['f64']['flag'], rec['f64']['certified'],
        ('guard' if ex is None else ex[fs(0.0, 3)]['dec']), ('guard' if ex is None else ex[fs(MU_CUT, 3)]['dec']),
        ('-' if ex is None else ex[fs(0.0, 3)]['prec']), rec['seconds']))
    return rec


def grid_points(lo, hi):
    k0 = int(round(lo / STEP))
    k1 = int(round(hi / STEP))
    return [k / 20.0 for k in range(k0, k1 + 1)]


def run_scan(P, lo, hi, bisect=True):
    pts = []
    first = {}
    for a in grid_points(lo, hi):
        r = evaluate(P, a)
        pts.append(a)
        for k in KEYS:
            if k not in first and r['flags'][k]:
                first[k] = len(pts) - 1
        if r.get('rank_guard'):
            log('rank guard: the sector is larger than the number of zero vectors, scan stops at a = %s' % r['a'])
            break
        if len(first) == len(KEYS) and len(pts) - 1 >= max(first.values()) + 2:
            break
    res = {}
    for k in KEYS:
        flags = [P['_cache'][round(a, 9)]['flags'][k] for a in pts]
        und = [P['_cache'][round(a, 9)]['undecided'].get(k, False) for a in pts]
        pattern = ''.join('?' if u else ('1' if f else '0') for f, u in zip(flags, und))
        if k not in first:
            res[k] = dict(first_grid=None, value=None, bracket=None, pattern=pattern,
                          note='not detected for a up to %s' % fs(pts[-1], 4))
            continue
        i = first[k]
        if i == 0:
            res[k] = dict(first_grid=fs(pts[0], 4), value=fs(pts[0], 4), bracket=None, pattern=pattern,
                          note='detected at the smallest a scanned')
            continue
        lo_a, hi_a = pts[i - 1], pts[i]
        if bisect:
            while hi_a - lo_a > TOL + 1e-12:
                mid = 0.5 * (lo_a + hi_a)
                if evaluate(P, mid)['flags'][k]:
                    hi_a = mid
                else:
                    lo_a = mid
        res[k] = dict(first_grid=fs(pts[i], 4), value=fs(hi_a, 6), bracket=[fs(lo_a, 6), fs(hi_a, 6)], pattern=pattern,
                      note=('bisected to width %s' % fs(hi_a - lo_a, 3)) if bisect else 'grid only')
        log('threshold %-5s first grid %s refined %s  pattern %s' % (k, fs(pts[i], 4), fs(hi_a, 6), pattern))
    return res, pts


# ---------------------------------------------------------------------------------------------------------------------------
# the comparison scale x_* of Remark yo:rem:scale
# ---------------------------------------------------------------------------------------------------------------------------
def d_on(gam, xs):
    return 8.0 * ((np.sin(0.5 * np.outer(xs, gam)) ** 2) / gam ** 2).sum(axis=1)


def xstar_first(gam, g0, delta, x_max, step_env, step_fine):
    """Two readings of 'D plus the pair term first turns negative' (Remark yo:rem:scale): x_env, the first x on a coarse grid at which the pair term can
    push D negative for some phase (the pair term is 8 Re[(1 - cos(z x))/z^2], z = g0 - i delta, whose minimum over the fast phase g0 x is
    8 (A - sqrt(cosh^2(delta x) A^2 + sinh^2(delta x) B^2)), A + iB = 1/z^2 conjugated; grid independent), and x_first, the first x with D < 0 on a grid
    of step step_fine (a fraction of the period 2 pi/g0 of the pair term) in a window after x_env."""
    z = complex(g0, -delta)
    A = (g0 * g0 - delta * delta) / (g0 * g0 + delta * delta) ** 2
    B = 2 * g0 * delta / (g0 * g0 + delta * delta) ** 2
    x_env = None
    C = 1000
    xs_all = np.arange(0.5, x_max + 0.5 * step_env, step_env)
    for i in range(0, len(xs_all), C):
        xs = xs_all[i:i + C]
        env = 8.0 * (A - np.hypot(np.cosh(delta * xs) * A, np.sinh(delta * xs) * B))
        neg = np.nonzero(d_on(gam, xs) + env < 0)[0]
        if len(neg):
            x_env = float(xs[neg[0]])
            break
    if x_env is None:
        return None, None
    x_first = None
    lo = max(0.5, x_env - 2 * 2 * math.pi / g0)
    hi = lo + min(6.0 / delta, 3e5 * step_fine)
    xs_all = np.arange(lo, hi, step_fine)
    for i in range(0, len(xs_all), C):
        xs = xs_all[i:i + C]
        pair = 16.0 * ((np.sin(0.5 * z * xs) ** 2) / (z * z)).real
        neg = np.nonzero(d_on(gam, xs) + pair < 0)[0]
        if len(neg):
            x_first = float(xs[neg[0]])
            break
    return x_env, x_first


def xstar_cases(zstr, cases):
    gam = np.array([float(s) for s in zstr])
    out = []
    for (label, g0, delta) in cases:
        x_max = max(60.0, 3.0 * (2 * math.log(g0) + math.log(KAPPA1)) / delta)
        # the pair term oscillates with period 2 pi / gamma_0 in x, so the fine grid must be a fraction of that or the first negative x aliases
        step_env = max(0.01, x_max / 2e5)
        step_fine = min(0.01, 2 * math.pi / (8 * g0))
        xe, xf = xstar_first(gam, g0, delta, x_max, step_env, step_fine)
        pred = (2 * math.log(g0) + math.log(KAPPA1)) / delta
        out.append(dict(label=label, gamma0=fs(g0, 10), delta=fs(delta, 4), x_env=(None if xe is None else fs(xe, 6)),
                        half_env=(None if xe is None else fs(xe / 2, 6)), x_first_negative=(None if xf is None else fs(xf, 6)),
                        half_first_negative=(None if xf is None else fs(xf / 2, 6)), predicted_x=fs(pred, 6), x_max=fs(x_max, 6),
                        step_env=fs(step_env, 4), step_fine=fs(step_fine, 4), zeros=len(gam)))
        log('xstar %-9s gamma=%-10s delta=%-6s x_env=%s, first negative on the fine grid %s (predicted %.2f); x_env/2 = %s' % (
            label, fs(g0, 8), fs(delta, 4), out[-1]['x_env'], out[-1]['x_first_negative'], pred, out[-1]['half_env']))
    return out


def midpoints_for_targets(zstr):
    fl = np.array([float(s) for s in zstr])
    mids = 0.5 * (fl[:-1] + fl[1:])
    res = {}
    for g in TARGETS:
        j = int(np.argmin(np.abs(mids - float(g))))
        if abs(mids[j] - float(g)) < 0.006:
            res[g] = float(mids[j])
    return res


# ---------------------------------------------------------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------------------------------------------------------
def versions():
    import flint
    import scipy
    return dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__, flint=flint.__version__)


def run_task(spec, outdir, zeros_path, k, nobisect=False):
    global LOGFH
    os.makedirs(outdir, exist_ok=True)
    name = 'task_%03d' % k
    LOGFH = open(os.path.join(outdir, name + '.log'), 'w')
    log('task', k, json.dumps(spec))
    if spec['kind'] == 'xstar':
        zstr, src = load_zero_strings(8000, zeros_path)
        mids = midpoints_for_targets(zstr)
        cases = [('yoshida', 100.0, 0.3), ('yoshida', 100.0, 0.2), ('yoshida', 1000.0, 0.3)]
        for g in TARGETS:
            if g in mids:
                for d in DELTAS:
                    cases.append((g, mids[g], d))
        res = xstar_cases(zstr, cases)
        out = dict(task=k, spec=spec, source=src, cases=res, versions=versions())
    else:
        need_flint()
        zstr, src = load_zero_strings(spec['M'], zeros_path)
        P = make_params(spec, zstr, src)
        log('zeros: %d from %s; gap %d, gamma_0 = %s (label %s), on line vectors %d, mode %s, band %s, extra %d' % (
            P['M'], src, P['j0'] + 1, fs(P['g0f'], 10), P['label'], P['nvec'], P['mode'], P['band'], P['extra']))
        res, pts = run_scan(P, spec['scan_lo'], spec['scan_hi'], bisect=not nobisect)
        evals = [P['_cache'][key] for key in sorted(P['_cache'])]
        out = dict(task=k, spec=spec, source=src, zeros_used=P['M'], on_line_vectors=P['nvec'], j0=P['j0'] + 1,
                   gamma0=fs(P['g0f'], 12), label=P['label'], delta=spec['delta'], thresholds=res, evals=evals,
                   scope=('full' if spec['band'] is None else 'band %s' % fs(spec['band'], 4)), versions=versions(),
                   chapter=CHAPTER[P['label']].get(spec['delta']), seconds=round(time.time() - T0, 1))
    with open(os.path.join(outdir, name + '.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    log('wrote', os.path.join(outdir, name + '.json'))
    LOGFH.close()
    LOGFH = None
    return out


def run_quick(args):
    need_flint()
    spec = dict(kind='scan', variant='quick', gi=0, delta=args.delta, M=2100, extra=45, band=None, mode='keep',
                scan_lo=0.9, scan_hi=1.4)
    zstr, src = load_zero_strings(2100, args.zeros)
    P = make_params(spec, zstr, src)
    log('quick: %d zeros from %s; gap %d, gamma_0 = %s, delta = %s' % (P['M'], src, P['j0'] + 1, fs(P['g0f'], 10), args.delta))
    res, pts = run_scan(P, spec['scan_lo'], spec['scan_hi'], bisect=True)
    evals = [P['_cache'][key] for key in sorted(P['_cache'])]
    cases = [('yoshida', 100.0, 0.3), ('yoshida', 100.0, 0.2), ('yoshida', 1000.0, 0.3), ('51.37', P['g0f'], args.delta)]
    xs = xstar_cases(zstr, cases)
    out = dict(spec=spec, source=src, zeros_used=P['M'], j0=P['j0'] + 1, gamma0=fs(P['g0f'], 12), label=P['label'],
               delta=args.delta, thresholds=res, evals=evals, xstar=xs, versions=versions(), seconds=round(time.time() - T0, 1))
    dst = args.out or os.path.join('/tmp', 'zh_detect_quick.json')
    with open(dst, 'w') as fh:
        json.dump(out, fh, indent=1)
    log('wrote', dst)
    return out


# ---------------------------------------------------------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------------------------------------------------------
def von_mangoldt(n):
    for p in range(2, int(n ** 0.5) + 2):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
    return math.log(n) if n > 1 else 0.0


def prime_side_entry(n, m, a):
    import mpmath as mp
    mp.mp.dps = 25
    k = n - m
    A = mp.mpf(a)

    def F(y):
        y = mp.mpf(y)
        if abs(y) > 2 * A:
            return mp.mpc(0)
        ph = mp.expj(mp.pi * m * y / A)
        if k == 0:
            return ph * (2 * A - abs(y)) / (2 * A)
        if y >= 0:
            return (-1) ** k * ph * (1 - mp.expj(mp.pi * k * y / A)) / (2j * mp.pi * k)
        return (-1) ** k * ph * (mp.expj(mp.pi * k * y / A) - 1) / (2j * mp.pi * k)
    pole = mp.quad(lambda y: F(y) * 2 * mp.cosh(y / 2), [-2 * A, 0, 2 * A])
    F0 = F(0)
    nmax = int(mp.floor(mp.exp(2 * A)))
    primes = sum(von_mangoldt(q) / mp.sqrt(q) * (F(mp.log(q)) + F(-mp.log(q))) for q in range(2, nmax + 1))
    u1 = min(mp.mpf(1), 2 * A)
    pf = mp.quad(lambda y: (F(y) + F(-y) - 2 * F0) / y, [0, u1])
    if 2 * A > u1:
        pf += mp.quad(lambda y: (F(y) + F(-y)) / y, [u1, 2 * A])
    def kern(y):
        if y < mp.mpf('1e-8'):
            return mp.mpf(1) / 4 - y / 48
        return mp.exp(y / 2) / (2 * mp.sinh(y)) - 1 / (2 * y)
    arch = mp.quad(lambda y: kern(y) * (F(y) + F(-y)), [0, 2 * A])
    V1 = -pf / 2 - (mp.euler + mp.log(2)) * F0 - arch
    return pole - mp.log(mp.pi) * F0 - primes + V1


def zero_side_entry(n, m, a, zstr):
    import mpmath as mp
    mp.mp.dps = 25
    nun, num = mp.pi * n / a, mp.pi * m / a
    s = mp.mpf(0)
    for z in zstr:
        g = mp.mpf(z)
        w = mp.sin(a * g) ** 2
        s += w / ((g + nun) * (g + num)) + w / ((-g + nun) * (-g + num))
    T = mp.mpf(zstr[-1])
    tail = (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)         # sum over gamma > T of gamma^-2 (mean density)
    return (-1) ** (n + m) * (2 / mp.mpf(a)) * (s + tail), (-1) ** (n + m) * (2 / mp.mpf(a)) * s


def selftest():
    need_flint()
    import mpmath as mp
    ok = True
    zstr, src = load_zero_strings(2100)
    # (A) the zero side formula against the explicit formula (yo:EF), a = 1, five entries
    for (n, m) in [(0, 0), (1, 0), (1, 1), (2, -1), (3, 3)]:
        pr = prime_side_entry(n, m, 1)
        zs, zs0 = zero_side_entry(n, m, 1, zstr[:1700])
        d = abs(pr - zs)
        log('(A) entry (%d,%d): prime side %s, zero side %s (with mean tail; without %s), |difference| = %s' % (
            n, m, mp.nstr(pr, 12), mp.nstr(zs, 12), mp.nstr(zs0, 8), mp.nstr(d, 3)))
        ok &= bool(d < 5e-5)
    # (B), (C) sector construction, pair term and exact decisions on a small model (45 zeros, delta = 0.45)
    spec = dict(kind='scan', variant='selftest', gi=0, delta=0.45, M=45, extra=45, band=None, mode='keep')
    P = make_params(spec, zstr, src)
    a = 1.0
    N = int(math.ceil(a * (P['g0f'] + 45) / math.pi)) + 2
    p = 256
    ctx.prec = p
    a_arb = arb(a)
    pi = arb.pi()
    Z = arb_zeros(P, p)
    g0 = (Z[P['j0']] + Z[P['j0'] + 1]) / 2
    dl = arb(str(P['delta']))
    quad = [acb(g0, -dl), acb(g0, dl), acb(-g0, -dl), acb(-g0, dl)]
    nu = {i: pi * i / a_arb for i in range(-N, N + 1)}

    def term(z, i, j):
        s = (a_arb * z).sin()
        return (2 / a_arb) * s * s / ((z + nu[i]) * (z + nu[j]))
    on_cache, pr_cache = {}, {}

    def Mon(i, j):
        if (i, j) not in on_cache:
            s = arb(0)
            for jj in P['use_idx']:
                g = Z[jj]
                s += term(acb(g), i, j).real + term(acb(-g), i, j).real
            on_cache[(i, j)] = on_cache[(j, i)] = s
        return on_cache[(i, j)]

    def Mpair(i, j):
        if (i, j) not in pr_cache:
            s = arb(0)
            for z in quad:
                s += term(z, i, j).real
            pr_cache[(i, j)] = pr_cache[(j, i)] = s
        return pr_cache[(i, j)]
    half = arb(1) / 2
    sq2 = arb(2).sqrt()

    def sector(fun, sec):
        nodes = list(range(0, N + 1)) if sec == 'even' else list(range(1, N + 1))
        rows = []
        for n in nodes:
            for m in nodes:
                if sec == 'even':
                    cn = arb(1) if n == 0 else sq2
                    cm = arb(1) if m == 0 else sq2
                    # <f_n, M f_m> with f_n = (e_n + e_-n)/c_n' (c_n' = sqrt 2 for n >= 1, 1 for n = 0)
                    if n == 0 and m == 0:
                        v = fun(0, 0)
                    elif n == 0:
                        v = (fun(0, m) + fun(0, -m)) / sq2
                    elif m == 0:
                        v = (fun(n, 0) + fun(-n, 0)) / sq2
                    else:
                        v = half * (fun(n, m) + fun(n, -m) + fun(-n, m) + fun(-n, -m))
                else:
                    v = half * (fun(n, m) - fun(n, -m) - fun(-n, m) + fun(-n, -m))
                rows.append(v)
        return nodes, rows
    t0 = time.time()
    res = {}
    for sec in ('even', 'odd'):
        nodes, dir_on = sector(Mon, sec)
        _, dir_pair = sector(Mpair, sec)
        res[sec] = (nodes, dir_on, dir_pair)
    log('(B) direct sector matrices built in %.1fs' % (time.time() - t0))
    B = exact_build(P, a, N, p)
    pv = pair_vectors(P, a, p, B['even'][0], B['odd'][0])
    for sec in ('even', 'odd'):
        nodes, dir_on, dir_pair = res[sec]
        n = len(nodes)
        e1 = max(abs(float(B[sec][1][i] - dir_on[i])) for i in range(n * n))
        pl, ql = pv[sec]
        e2 = max(abs(float(4 * pl[i] * pl[j] - 4 * ql[i] * ql[j] - dir_pair[i * n + j])) for i in range(n) for j in range(n))
        log('(B) %s sector (size %d): max |Loewner - direct| on line = %.2e, max |4pp^T - 4qq^T - direct quadruple| = %.2e' % (sec, n, e1, e2))
        ok &= (e1 < 1e-40 and e2 < 1e-40)
    # (C) exact decisions against mpmath eigenvalues, both sectors, shifts around the smallest eigenvalue
    mp.mp.dps = 40
    for sec in ('even', 'odd'):
        nodes, dir_on, dir_pair = res[sec]
        n = len(nodes)
        Mm = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                Mm[i, j] = mp.mpf(str(dir_on[i * n + j].mid().str(45, radius=False))) + mp.mpf(str(dir_pair[i * n + j].mid().str(45, radius=False)))
        ev = mp.eigsy(Mm, eigvals_only=True)
        lam = min(ev[i] for i in range(n))
        W = make_W(pv, sec)
        shifts = [0.0, MU_CUT]
        if lam < 0:
            shifts += [float(0.9 * -lam), float(1.1 * -lam)]
        row = []
        for mu in shifts:
            dec = decide(B[sec][1], n, W, mu, p)[0]
            want = 'neg' if lam < -mp.mpf(mu) else 'pos'
            row.append('mu=%s %s (expect %s)' % (fs(mu, 3), dec, want))
            ok &= (dec == want)
        log('(C) %s sector: smallest mpmath eigenvalue %s; %s' % (sec, mp.nstr(lam, 10), '; '.join(row)))
        if sec == 'even':
            lam_even = lam
    # (D) float64 sector eigenvalues against mpmath on the same model
    f = f64_eval(P, a, N)
    log('(D) float64 even sector lambda_min %s against mpmath %s' % (f['lam_even'], mp.nstr(lam_even, 10)))
    ok &= abs(float(f['lam_even']) - float(lam_even)) < 1e-9
    # (E) independent inertia check at (199.64, 0.2): unpivoted LDL^T of the even sector matrix in ball arithmetic against the Schur test
    spec2 = dict(kind='scan', variant='selftest2', gi=1, delta=0.2, M=1700, extra=45, band=None, mode='keep')
    P2 = make_params(spec2, zstr, src)
    for a2, want in ((1.45, 0), (1.6, 1)):
        N2 = int(math.ceil(a2 * (P2['g0f'] + 45) / math.pi)) + 2
        p2 = 1536
        B2 = exact_build(P2, a2, N2, p2)
        pv2 = pair_vectors(P2, a2, p2, B2['even'][0], B2['odd'][0])
        nodes2, ent2 = B2['even']
        n2 = len(nodes2)
        pl, ql = pv2['even']
        ctx.prec = p2
        Ms = [ent2[i * n2 + j] + 4 * pl[i] * pl[j] - 4 * ql[i] * ql[j] for i in range(n2) for j in range(n2)]
        neg = ldl_negative_pivots(Ms, n2)
        dec = decide(ent2, n2, make_W(pv2, 'even'), 0.0, p2)[0]
        log('(E) a = %s, even sector of size %d at %d bits: negative pivots of LDL^T = %s, Schur test says %s (expected %d negative eigenvalue)' % (a2, n2, p2, neg, dec, want))
        ok &= (neg == want and dec == ('neg' if want else 'pos'))
    log('selftest', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


def ldl_negative_pivots(entries, n):
    """Number of negative pivots of an unpivoted LDL^T factorisation in ball arithmetic (Sylvester inertia); None if a pivot ball contains zero."""
    M = [[entries[i * n + j] for j in range(n)] for i in range(n)]
    d = []
    L = [[arb(0)] * n for _ in range(n)]
    neg = 0
    for j in range(n):
        s = M[j][j]
        for k in range(j):
            s -= L[j][k] * L[j][k] * d[k]
        if s < 0:
            neg += 1
        elif not (s > 0):
            return None
        d.append(s)
        for i in range(j + 1, n):
            t = M[i][j]
            for k in range(j):
                t -= L[i][k] * L[j][k] * d[k]
            L[i][j] = t / s
    return neg


# ---------------------------------------------------------------------------------------------------------------------------
# gather, figure
# ---------------------------------------------------------------------------------------------------------------------------
def fv(x):
    return None if x is None else float(x)


def gather(args):
    files = sorted(glob.glob(os.path.join(args.tasks_dir, 'task_*.json')))
    if not files:
        sys.exit('no task files in %s' % args.tasks_dir)
    cases, xstar = [], None
    for f in files:
        d = json.load(open(f))
        if d['spec']['kind'] == 'xstar':
            xstar = d
            continue
        c = {k: d[k] for k in ('task', 'spec', 'source', 'zeros_used', 'on_line_vectors', 'j0', 'gamma0', 'label', 'delta',
                               'thresholds', 'scope', 'chapter', 'seconds', 'versions')}
        c['variant'] = d['spec']['variant']
        if c['variant'] == 'main':
            c['evals'] = d['evals']
        cases.append(c)
    logname = os.path.splitext(args.out)[0] + '.log'
    lines = []

    def pr(s=''):
        print(s)
        lines.append(s)
    pr('zh_detect gather: %d scan tasks%s' % (len(cases), ', xstar present' if xstar else ', xstar missing'))
    main = [c for c in cases if c['variant'] == 'main']
    pr('')
    pr('Table: smallest a on the grid of step 0.05 (refined value after bisection to 0.005 in brackets); f64 = chapter criterion in float64,')
    pr('ex0 = exact lambda_min < 0, ex11 = exact lambda_min < -1e-11; scope is the full matrix for gamma_0 <= 200 and the band 150 above (band detection is certified, band absence is not)')
    pr('f64c = float64 eigenvalue below minus the Weyl bound (rigorous detection of lambda_min < 0, not an exact threshold)')
    pr('%-8s %-6s %-8s %-8s %-16s %-16s %-16s %-16s %s' % ('gamma_0', 'delta', 'chapter', 'horizon', 'f64', 'f64c', 'ex0', 'ex11', 'x*/2'))
    xmap = {}
    if xstar:
        for x in xstar['cases']:
            xmap[(x['label'], x['delta'])] = x['half_env']
    for c in sorted(main, key=lambda c: (float(c['label']), c['delta'])):
        hz = 0.5 * math.log(float(c['gamma0']) / (2 * math.pi))
        cells = []
        for k in KEYS:
            t = c['thresholds'][k]
            cells.append('%s (%s)' % (t['first_grid'], t['value']) if t['first_grid'] else 'none')
        pr('%-8s %-6s %-8s %-8.3f %-16s %-16s %-16s %-16s %s' % (c['label'], fs(c['delta'], 4), c['chapter'] if c['chapter'] is not None else '-',
                                                               hz, cells[0], cells[1], cells[2], cells[3], xmap.get((c['label'], fs(c['delta'], 4)), '-')))
    pr('')
    pr('Sensitivity (refined values; difference to the main run in brackets)')
    mainmap = {(c['label'], c['delta']): c for c in main}
    for c in sorted([c for c in cases if c['variant'] != 'main'], key=lambda c: (c['variant'], float(c['label']), c['delta'])):
        m = mainmap.get((c['label'], c['delta']))
        cells = []
        for k in KEYS:
            v = fv(c['thresholds'][k]['value'])
            w = fv(m['thresholds'][k]['value']) if m else None
            cells.append('%s (%s)' % (fs(v, 5), fs(v - w, 3) if (v is not None and w is not None) else '-') if v is not None else 'none')
        pr('%-8s %-8s %-6s %-20s %-20s %-20s %-20s' % (c['variant'], c['label'], fs(c['delta'], 4), cells[0], cells[1], cells[2], cells[3]))
    meta = dict(script='code/zh_detect.py', tasks_dir=args.tasks_dir, scan_tasks=len(cases), keys=KEYNAME, chapter_table=CHAPTER,
                grid_step=STEP, bisection_tolerance=TOL, cut=MU_CUT, note='thresholds are the smallest tested half length a with detection; value is the upper end of the final bisection bracket')
    out = dict(meta=meta, cases=cases, xstar=(xstar['cases'] if xstar else None))
    with open(args.out, 'w') as fh:
        json.dump(out, fh, indent=1)
    with open(logname, 'w') as fh:
        for f in files:
            lf = os.path.splitext(f)[0] + '.log'
            if os.path.isfile(lf):
                fh.write('=== %s ===\n' % os.path.basename(lf))
                fh.write(open(lf).read())
        fh.write('=== summary ===\n' + '\n'.join(lines) + '\n')
    pr('wrote %s and %s' % (args.out, logname))


def make_figure(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 0.9})
    d = json.load(open(args.data))
    if 'cases' in d:
        cases = [c for c in d['cases'] if c['variant'] == 'main']
    else:
        cases = [dict(variant='quick', gamma0=d['gamma0'], delta=d['delta'], thresholds=d['thresholds'])]
    marks = {0.001: ('x', 4.0), 0.003: ('+', 4.5), 0.01: ('o', 3.2), 0.05: ('s', 3.0), 0.2: ('^', 3.4), 0.45: ('D', 2.8)}
    fig, axes = plt.subplots(1, 2, figsize=(5.8, 2.6), sharey=True)
    gg = np.linspace(40, 2600, 300)
    def best(c, key):
        # right panel: the smallest rigorous detection of lambda_min < 0 (exact band or full test, or certified float64)
        ks = ('ex0', 'f64c') if key == 'ex0' else (key,)
        vals = [float(c['thresholds'][k]['value']) for k in ks if c['thresholds'].get(k, {}).get('value') is not None]
        return min(vals) if vals else None
    for ax, key, title in zip(axes, ('f64', 'ex0'), ('float64, cut $-10^{-11}$', 'exact sign (bound for $\\gamma\\geq1000$)')):
        ax.plot(gg, 0.5 * np.log(gg / (2 * np.pi)), color='black', ls='--', lw=0.8, label=r'$\frac{1}{2}\log(\gamma/2\pi)$')
        for delta, (mk, ms) in marks.items():
            xs, ys = [], []
            for c in sorted(cases, key=lambda c: float(c['gamma0'])):
                v = best(c, key)
                if abs(float(c['delta']) - delta) < 1e-12 and v is not None:
                    xs.append(float(c['gamma0']))
                    ys.append(v)
            if xs:
                filled = mk in 'osD^'
                ax.plot(xs, ys, color='0.6', lw=0.5)
                ax.plot(xs, ys, mk, color='black', mfc=('black' if delta in (0.01, 0.2) else 'white') if filled else 'black', ms=ms,
                        mew=0.7, ls='none', label=r'$\delta=%g$' % delta)
        ax.set_xscale('log')
        ax.set_xlim(40, 2600)
        ax.set_xlabel(r'height $\gamma_0$ of the quadruple')
        ax.set_title(title, fontsize=8)
    axes[0].set_ylabel(r'smallest half length $a$ that detects')
    axes[0].legend(loc='upper left', frameon=False, fontsize=6, ncol=1)
    fig.tight_layout()
    fig.savefig(args.fig)
    plt.close(fig)
    log('wrote', args.fig)


def main():
    global LOGFH
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--list-tasks', action='store_true')
    ap.add_argument('--count-tasks', action='store_true')
    ap.add_argument('--task', type=int)
    ap.add_argument('--gather', action='store_true')
    ap.add_argument('--figure', action='store_true')
    ap.add_argument('--delta', type=float, default=0.2)
    ap.add_argument('--scan', type=float, nargs=2, metavar=('LO', 'HI'))
    ap.add_argument('--M', type=int)
    ap.add_argument('--nobisect', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--outdir', default='out')
    ap.add_argument('--zeros')
    ap.add_argument('--tasks-dir', default='out')
    ap.add_argument('--data', default=None)
    ap.add_argument('--fig', default=None)
    args = ap.parse_args()
    tasks = task_table()
    if args.count_tasks:
        print(len(tasks))
        return 0
    if args.list_tasks:
        for k, t in enumerate(tasks):
            print(k, json.dumps(t))
        return 0
    if args.selftest:
        return selftest()
    if args.quick:
        run_quick(args)
        return 0
    if args.task is not None:
        spec = dict(tasks[args.task])
        if args.scan:                      # overrides for local probes; the sbatch array never uses them
            spec['scan_lo'], spec['scan_hi'] = args.scan
        if args.M:
            spec['M'] = args.M
            spec['variant'] = spec['variant'] + '_M%d' % args.M
        run_task(spec, args.outdir, args.zeros, args.task, nobisect=args.nobisect)
        return 0
    if args.gather:
        args.out = args.out or str(ROOT / 'data' / 'zh_detect.json')
        gather(args)
        return 0
    if args.figure:
        args.data = args.data or str(ROOT / 'data' / 'zh_detect.json')
        args.fig = args.fig or str(ROOT / 'fig' / 'zh_horizon.pdf')
        make_figure(args)
        return 0
    ap.print_help()
    return 2


if __name__ == '__main__':
    sys.exit(main())
