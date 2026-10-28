# fc_ly.py: group FC, item 4. The Farey fraction chain of Observation obs:ly: Z_N(beta, x) = sum over words sigma in {A,B}^N of (Tr M_N(sigma))^(-beta) x^|sigma|,
# A = [[1,0],[1,1]], B = A^T, |sigma| = number of letters A; beta_LY(N) = least beta such that all zeros of Z_N(beta, .) lie on |x| = 1 (for all larger beta);
# the merging angle theta_N: the two zeros that leave the circle merge on it at x = exp(+-i theta_N) at beta = beta_LY(N).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Stage A (exact): the table n_N(k, t) = #{words with k letters A and trace t}, by splitting a word into a prefix (N//2 letters) and a suffix and
# Tr(PS) = P00 S00 + P01 S10 + P10 S01 + P11 S11 in int64; the sparse table is written to --tables (npz) and checked: sum n = 2^N, symmetry
# n(k,t) = n(N-k,t) (the flip + reversal bijection), max trace = Lucas(N), n(0,.) = n(N,.) at trace 2.
# Stage B (even N): c_k(beta) = sum_t n(k,t) t^(-beta), palindromic c_k = c_(N-k). With y = x + 1/x, x^(-N/2) Z = q(y) = c_M + sum_(j=1..M) c_(M+j) C_j(y), M = N/2,
# C_0 = 2, C_1 = y, C_(j+1) = y C_j - C_(j-1) (C_j(2 cos phi) = 2 cos j phi). All zeros on |x| = 1 (simple) iff q has M distinct real roots in (-2, 2).
# (i) double precision scan on a beta grid (0.2 to 4, step 0.1, mpmath polyroots at 40 digits, a root is real if |Im| < 1e-12): the number of transitions
# of the predicate (expected 1), then bisection of the transition to 1e-12. (ii) certified: for N <= --certmax, python-flint balls (--prec bits), the
# coefficients c_k(beta) by ball arithmetic from the exact table, and the number of distinct real roots of q in (-2, 2) by a Sturm sequence in arb_poly
# (all signs certified, otherwise the step is reported UNCERTIFIED); the bracket [lo, hi] has "count = M" certified at hi and "count < M" certified at lo,
# then 6 more certified bisection steps. The monotone assumption (Lee Yang for all beta >= hi) is checked on the grid only, not proved.
# theta_N: at hi, the real roots of q (mpmath), the closest pair y1 < y2, theta = arccos((y1+y2)/4), uncertainty from arccos of y1/2 and y2/2.
# Extra: the quartet at N = 20, beta = 1.2 (zeros h = log x with |x| != 1 of Z_20(1.2, x)).
import sys, os, json, time, argparse
import numpy as np
from multiprocessing import Pool
import mpmath as mp

def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def letters():
    return np.array([[1, 0], [1, 1]], dtype=np.int64), np.array([[1, 1], [0, 1]], dtype=np.int64)

def all_products(n):
    """arrays (m00, m01, m10, m11, k) over all 2^n words of length n"""
    A, B = letters()
    m = np.array([[1, 0, 0, 1]], dtype=np.int64)
    k = np.zeros(1, dtype=np.int64)
    for _ in range(n):
        def mul(m, X):
            a, b, c, d = m[:, 0], m[:, 1], m[:, 2], m[:, 3]
            return np.stack([a * X[0, 0] + b * X[1, 0], a * X[0, 1] + b * X[1, 1], c * X[0, 0] + d * X[1, 0], c * X[0, 1] + d * X[1, 1]], axis=1)
        m = np.concatenate([mul(m, A), mul(m, B)])
        k = np.concatenate([k + 1, k])
    return m, k

def slice_task(args):
    N, i0, i1, tmax = args
    n1 = N // 2; n2 = N - n1
    P, kp = all_products(n1)
    S, ks = all_products(n2)
    size = (N + 1) * (tmax + 1)
    tab = np.zeros(size, dtype=np.int64)
    S00, S01, S10, S11 = S[:, 0], S[:, 1], S[:, 2], S[:, 3]
    blk = max(1, (1 << 21) // len(S))
    for lo in range(i0, i1, blk):
        hi = min(i1, lo + blk)
        p = P[lo:hi]
        tr = p[:, 0:1] * S00[None, :] + p[:, 1:2] * S10[None, :] + p[:, 2:3] * S01[None, :] + p[:, 3:4] * S11[None, :]
        key = (kp[lo:hi, None] + ks[None, :]) * (tmax + 1) + tr
        key = key.ravel()
        assert key.max() < size
        u, cnt = np.unique(key, return_counts=True)
        tab[u] += cnt
    nz = np.nonzero(tab)[0]
    return nz, tab[nz]

def stage_a(N, nproc):
    tmax = lucas(N)
    n1 = N // 2
    ntask = min(2 ** n1, 4 * nproc)
    edges = np.linspace(0, 2 ** n1, ntask + 1).astype(int)
    tasks = [(N, int(edges[i]), int(edges[i + 1]), tmax) for i in range(ntask)]
    tab = np.zeros((N + 1) * (tmax + 1), dtype=np.int64)
    with Pool(nproc) as pool:
        for nz, v in pool.imap_unordered(slice_task, tasks):
            tab[nz] += v
    nz = np.nonzero(tab)[0]
    k = nz // (tmax + 1); t = nz % (tmax + 1)
    cnt = tab[nz]
    chk = {'sum_n': int(cnt.sum()), 'two_pow_N': 2 ** N, 'max_trace': int(t.max()), 'lucas_N': tmax, 'nnz': int(len(nz))}
    # symmetry n(k,t) = n(N-k,t)
    tab2 = tab.reshape(N + 1, tmax + 1)
    chk['symmetric'] = bool(np.array_equal(tab2, tab2[::-1]))
    chk['n_k0_trace2'] = int(tab2[0, 2]); chk['n_kN_trace2'] = int(tab2[N, 2])
    return k.astype(np.int16), t.astype(np.int64), cnt, chk

def coeffs_float(k, lt, cnt, N, beta):
    w = cnt * np.exp(-beta * lt)
    return np.bincount(k, weights=w, minlength=N + 1)

def cheb_polys(M):
    C = [[2], [0, 1]]
    for j in range(1, M):
        a = [0] + C[j]
        b = C[j - 1] + [0] * (len(a) - len(C[j - 1]))
        C.append([x - y for x, y in zip(a, b)])
    return C  # monomial coefficients low -> high

def q_coeffs(c, N):
    """monomial coefficients (low->high, mp.mpf) of q(y) from c_0..c_N (palindromic assumed: uses c_M..c_N)"""
    M = N // 2
    C = cheb_polys(M)
    q = [mp.mpf(0)] * (M + 1)
    q[0] += mp.mpf(c[M])
    for j in range(1, M + 1):
        for i, a in enumerate(C[j]):
            q[i] += mp.mpf(c[M + j]) * a
    return q

def real_roots_float(q):
    mp.mp.dps = 50
    roots = mp.polyroots(list(reversed(q)), maxsteps=500, extraprec=300)
    return roots

def ly_holds(c, N):
    q = q_coeffs(c, N)
    roots = real_roots_float(q)
    ok = all(abs(mp.im(r)) < 1e-12 and abs(mp.re(r)) < 2 + 1e-12 for r in roots)
    return ok, roots

# ------- certified stage
def arb_coeffs(k, t, cnt, N, beta_str, prec, cache):
    from flint import arb, ctx
    ctx.prec = prec
    if 'ut' not in cache:
        ut, inv = np.unique(t, return_inverse=True)
        cache['ut'] = ut; cache['inv'] = inv
        cache['lt'] = [arb(int(x)).log() for x in ut]
    beta = arb(beta_str)
    w = [(-(beta * l)).exp() for l in cache['lt']]
    cs = [arb(0)] * (N + 1)
    inv = cache['inv']
    for kk, ii, nn in zip(k.tolist(), inv.tolist(), cnt.tolist()):
        cs[kk] += nn * w[ii]
    return cs

def sturm_count(cs, N, prec):
    from flint import arb, arb_poly, ctx
    ctx.prec = prec
    M = N // 2
    C = cheb_polys(M)
    q = [arb(0)] * (M + 1)
    q[0] = q[0] + cs[M]
    for j in range(1, M + 1):
        for i, a in enumerate(C[j]):
            if a:
                q[i] = q[i] + cs[M + j] * a
    p0 = arb_poly(q)
    seq = [p0, p0.derivative()]
    while seq[-1].degree() > 0:
        qq, r = divmod(seq[-2], seq[-1])
        r = -r
        # strip top coefficients that are certified zero? require lead nonzero
        d = r.degree()
        lead = r[d] if d >= 0 else None
        if d < 0:
            raise RuntimeError('UNCERTIFIED: zero remainder')
        if not (lead > 0 or lead < 0):
            raise RuntimeError('UNCERTIFIED: leading coefficient contains 0 at degree %d' % d)
        if d != seq[-1].degree() - 1:
            raise RuntimeError('UNCERTIFIED/degenerate: degree drop %d -> %d' % (seq[-1].degree(), d))
        seq.append(r)
    def signs(x):
        out = []
        for p in seq:
            v = p(arb(x))
            if v > 0:
                out.append(1)
            elif v < 0:
                out.append(-1)
            else:
                raise RuntimeError('UNCERTIFIED: sign at %s' % x)
        return out
    def var(s):
        return sum(1 for a, b in zip(s[:-1], s[1:]) if a != b)
    return var(signs(-2)) - var(signs(2))

def ly_task(args):
    N, tables, nproc_unused, certmax, prec, nsteps = args
    t0 = time.time()
    z = np.load(tables)
    k = z['k%d' % N].astype(np.int64); t = z['t%d' % N]; cnt = z['n%d' % N]
    lt = np.log(t.astype(np.float64)); cntf = cnt.astype(np.float64)
    M = N // 2
    res = {'N': N}
    # grid scan
    grid = [round(0.2 + 0.1 * i, 3) for i in range(0, 39)]
    pat = []
    for b in grid:
        c = coeffs_float(k, lt, cntf, N, b)
        ok, _ = ly_holds(c, N)
        pat.append(ok)
    res['grid_pattern'] = ''.join('1' if x else '0' for x in pat)
    trans = [i for i in range(len(pat) - 1) if pat[i] != pat[i + 1]]
    res['grid_transitions'] = [grid[i] for i in trans]
    # bisect first False->True transition
    i0 = [i for i in trans if (not pat[i]) and pat[i + 1]]
    if not i0:
        res['error'] = 'no False->True transition on the grid'
        return res
    lo, hi = grid[i0[-1]], grid[i0[-1] + 1]
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        c = coeffs_float(k, lt, cntf, N, mid)
        if ly_holds(c, N)[0]:
            hi = mid
        else:
            lo = mid
    res['float_lo'] = lo; res['float_hi'] = hi
    # theta at hi
    c = coeffs_float(k, lt, cntf, N, hi)
    ok, roots = ly_holds(c, N)
    ys = sorted([mp.re(r) for r in roots])
    gaps = [(ys[i + 1] - ys[i], i) for i in range(len(ys) - 1)]
    g, i = min(gaps)
    ym = (ys[i] + ys[i + 1]) / 2
    res['closest_pair_y'] = [float(ys[i]), float(ys[i + 1])]
    res['theta_float'] = float(mp.acos(ym / 2)); res['theta_halfwidth'] = float(abs(mp.acos(ys[i] / 2) - mp.acos(ys[i + 1] / 2)) / 2)
    res['pair_index_from_top'] = len(ys) - 1 - i
    res['N_theta'] = N * res['theta_float']
    res['theta_over_2pi_N'] = res['theta_float'] / (2 * np.pi / N)
    # extra: quartet at N = 20, beta = 1.2
    if N == 20:
        c = coeffs_float(k, lt, cntf, N, 1.2)
        q = q_coeffs(c, N)
        roots = real_roots_float(q)
        quartet = []
        for r in roots:
            if abs(mp.im(r)) > 1e-9 or abs(mp.re(r)) > 2:
                x = (r + mp.sqrt(r * r - 4)) / 2
                h = mp.log(x)
                quartet.append([float(mp.re(h)), float(mp.im(h))])
                x2 = (r - mp.sqrt(r * r - 4)) / 2
                h2 = mp.log(x2)
                quartet.append([float(mp.re(h2)), float(mp.im(h2))])
        res['quartet_h_N20_beta1.2'] = quartet
    # certified stage
    if N <= certmax:
        cache = {}
        try:
            bhi = '%.12f' % (hi + 1e-9); blo = '%.12f' % (lo - 1e-9)
            c_hi = arb_coeffs(k, t, cnt, N, bhi, prec, cache)
            n_hi = sturm_count(c_hi, N, prec)
            c_lo = arb_coeffs(k, t, cnt, N, blo, prec, cache)
            n_lo = sturm_count(c_lo, N, prec)
            res['cert_hi'] = float(bhi); res['cert_hi_count'] = int(n_hi); res['cert_lo'] = float(blo); res['cert_lo_count'] = int(n_lo); res['M'] = M
            if n_hi == M and n_lo < M:
                a, b = float(blo), float(bhi)
                A, B = blo, bhi
                for _ in range(nsteps):
                    mid = '%.14f' % ((float(A) + float(B)) / 2)
                    cm = arb_coeffs(k, t, cnt, N, mid, prec, cache)
                    n = sturm_count(cm, N, prec)
                    if n == M:
                        B = mid
                    else:
                        A = mid
                res['cert_bracket'] = [float(A), float(B)]
                res['certified'] = True
            else:
                res['certified'] = False
                res['note'] = 'bracket counts not (M, <M)'
        except RuntimeError as e:
            res['certified'] = False; res['note'] = str(e)
    res['seconds'] = time.time() - t0
    print('N=%d float beta_LY in [%.10f, %.10f] theta %.6f (+-%.1e) transitions %s cert %s %s | %.0f s' % (N, lo, hi, res['theta_float'], res['theta_halfwidth'], res['grid_transitions'], res.get('certified'), res.get('cert_bracket'), res['seconds']), flush=True)
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmin', type=int, default=10)
    ap.add_argument('--nmax', type=int, default=30)
    ap.add_argument('--certmax', type=int, default=26)
    ap.add_argument('--prec', type=int, default=400)
    ap.add_argument('--steps', type=int, default=6)
    ap.add_argument('--nproc', type=int, default=16)
    ap.add_argument('--tables', default='fc_ly_tables.npz')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--stage', default='ab')
    ap.add_argument('--out', default='fc_ly.json')
    args = ap.parse_args()
    if args.quick:
        args.nmax = 12; args.certmax = 12; args.nproc = 2; args.steps = 2; args.prec = 200
    Ns = list(range(args.nmin, args.nmax + 1, 2))
    out = {'stageA': {}, 'stageB': []}
    if 'a' in args.stage:
        save = {}
        for N in Ns:
            t0 = time.time()
            k, t, cnt, chk = stage_a(N, args.nproc)
            chk['seconds'] = time.time() - t0
            out['stageA'][str(N)] = chk
            save['k%d' % N] = k; save['t%d' % N] = t; save['n%d' % N] = cnt
            print('A N=%d %s' % (N, chk), flush=True)
        np.savez(args.tables, **save)
    if 'b' in args.stage:
        tasks = [(N, args.tables, 0, args.certmax, args.prec, args.steps) for N in sorted(Ns, key=lambda N: -N)]
        with Pool(min(args.nproc, len(tasks))) as pool:
            out['stageB'] = sorted(pool.map(ly_task, tasks, chunksize=1), key=lambda r: r['N'])
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1, default=float)

if __name__ == '__main__':
    main()
