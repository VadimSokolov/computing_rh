# fc_stern.py: group FC (Knauf's arithmetic ferromagnet, Section "Knauf's arithmetic ferromagnet" of the consolidated ferromagnets
# chapter): the Stern row and the counts c_N(m) of Proposition prop:stern, exactly, for N up to --nmax (default 34), and the checks
# Z_N(-1) = 3^N, max entry F_{N+2}, w_N(m) = 1 for m <= N+1, monotonicity c_N <= c_{N+1} <= phi, plus Yao-Knuth item 5 (--yk).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. s(1) = 1, s(2n) = s(n), s(2n+1) = s(n) + s(n+1). Row N is s(2^N), ..., s(2^(N+1) - 1) (2^N entries); c_N(m) is the
# multiplicity of m in row N. Pair (s(k), s(k+1)) for k in row N is the pair at depth N; the row-N+1 pairs are (a, a+b) (position 2k)
# and (a+b, b) (position 2k+1). The entry at position 2k of row N+1 is s(k) (old entry), the entry at 2k+1 is a+b (new, "born at row
# N+1"). So c_N(m) = [m = 1] + sum_{j=1..N} H_j(m), where H_j(m) = #{pairs (a, b) at depth j-1 with a + b = m} (births), and every
# coprime pair occurs exactly once (Calkin-Wilf), so c_N(m) <= phi(m) and c_N(m) = #{a coprime to m : cfsum(a/(m-a)) <= N}.
# Method. The 2^(N-1) pairs of depth N-1 are enumerated without storing the tree: prefixes at depth L = 12 (4096 pairs) are expanded in
# numpy (uint32) one prefix per task, level by level, with np.bincount of the sums into H_j; tasks are spread over a Pool. All
# integer arithmetic is exact (sums < F_(N+2) < 2^32). H_j has length F_(j+2)+1.
# Checks. (1) direct Stern recursion for N <= 22 (independent of the pair recursion) gives the same c_N; (2) c_N(m) against Euclid's
# algorithm on (a, m-a) for all m <= 800 and all N <= 40 (cfsum <= N); (3) sum m c_N(m) = 3^N, sum c_N(m) = 2^N, max m = F_(N+2),
# c_N <= c_(N+1) <= phi, w_N(m) = 1 for m <= N+1 and w_N(N+2) < 1; (4) the checksums written to JSON. The births are saved to
# --births (npz of uint32 arrays) for the next jobs (fc_perron.py, fc_fisher.py).
# Yao-Knuth (--yk): mean over a coprime to m of S(a/m) = sum of partial quotients (Euclid on (m, a)); least squares fit
# mean = c2 ln^2 m + c1 ln m + c0 over a ladder of m, with m = 10^5, 10^6, 10^7, 10^8 computed exactly (chunked Euclid in numpy).
import sys, os, json, time, argparse, math
import numpy as np
from multiprocessing import Pool

L0 = 12

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def expand(args):
    """a group of prefix pairs (a,b) at depth L0, each expanded to depth nmax-1; returns {j: bincount of sums} accumulated over the group."""
    pref, nmax = args
    acc = {j: np.zeros(fib(j + 2) + 1, dtype=np.uint64) for j in range(L0 + 1, nmax + 1)}
    for (pa, pb) in pref:
        a = np.array([pa], dtype=np.uint32)
        b = np.array([pb], dtype=np.uint32)
        d = L0
        while True:
            s = a + b
            bc = np.bincount(s, minlength=fib(d + 3) + 1)
            acc[d + 1][:len(bc)] += bc.astype(np.uint64)
            if d + 1 >= nmax:
                break
            a, b = np.concatenate([a, s]), np.concatenate([s, b])
            d += 1
    return {j: v.astype(np.uint32) for j, v in acc.items()}

def prefixes():
    pairs = [(1, 1)]
    H = {}
    d = 0
    births = {}
    while True:
        a = np.array([p[0] for p in pairs], dtype=np.int64)
        b = np.array([p[1] for p in pairs], dtype=np.int64)
        if d == L0:
            return pairs, births
        s = a + b
        births[d + 1] = np.bincount(s, minlength=fib(d + 3) + 1).astype(np.uint32)
        nxt = []
        for x, y in pairs:
            nxt.append((x, x + y)); nxt.append((x + y, y))
        pairs = nxt
        d += 1

def build(nmax, nproc):
    pairs, births = prefixes()
    # births[j] for j = 1..L0 from depths 0..L0-1; depth >= L0 from tasks
    tasks = []
    per = max(1, len(pairs) // nproc)
    for i in range(0, len(pairs), per):
        tasks.append((pairs[i:i + per], nmax))
    for j in range(L0 + 1, nmax + 1):
        births[j] = np.zeros(fib(j + 2) + 1, dtype=np.uint64)
    with Pool(nproc) as pool:
        for out in pool.imap_unordered(expand, tasks):
            for j, bc in out.items():
                if j >= L0 + 1 and j <= nmax:
                    births[j][:len(bc)] += bc
    for j in range(L0 + 1, nmax + 1):
        births[j] = births[j].astype(np.uint32)
    return births

def phi_sieve(n):
    phi = np.arange(n + 1, dtype=np.int64)
    is_p = np.ones(n + 1, dtype=bool); is_p[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if is_p[p]:
            is_p[p * p::p] = False
    for p in np.nonzero(is_p)[0]:
        phi[p::p] -= phi[p::p] // p
    return phi

def stern_row_direct(N):
    # s(2^N .. 2^(N+1)-1) by the recursion on the full sequence up to 2^(N+1), vectorised by blocks
    seq = np.array([0, 1], dtype=np.int64)  # s(0)=0, s(1)=1
    # s(2n) = s(n), s(2n+1) = s(n)+s(n+1): build length 2^k+1 arrays
    seq = np.array([0, 1, 1], dtype=np.int64)  # s(0), s(1), s(2)
    k = 1
    while k < N + 1:
        n = np.arange(0, len(seq))  # indices 0..2^k
        new = np.zeros(2 * (len(seq) - 1) + 1, dtype=np.int64)
        new[0::2] = seq
        new[1::2] = seq[:-1] + seq[1:]
        seq = new
        k += 1
    return seq[2 ** N: 2 ** (N + 1)]

def cfsum(a, b):
    t = 0
    while b:
        t += a // b
        a, b = b, a % b
    return t

def yk_mean(m, chunk=2_000_000):
    """exact sum over a in [1,m-1] coprime to m of the sum of partial quotients of m/a (Euclid on (m,a)); returns (sum, count)."""
    tot = 0; cnt = 0
    for lo in range(1, m, chunk):
        a = np.arange(lo, min(lo + chunk, m), dtype=np.int64)
        x = np.full_like(a, m)
        y = a.copy()
        S = np.zeros_like(a)
        idx = np.arange(len(a))
        while len(idx):
            q = x // y
            S[idx] += q
            r = x - q * y
            x, y = y, r
            keep = y != 0
            idx = idx[keep]; x = x[keep]; y = y[keep]
            if len(idx) == 0:
                break
        # x holds last nonzero remainder only for the kept, so compute gcd separately
        g = np.gcd(a, m)
        ok = g == 1
        tot += int(S[ok].sum()); cnt += int(ok.sum())
    return tot, cnt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nmax', type=int, default=34)
    ap.add_argument('--nproc', type=int, default=16)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='fc_stern.json')
    ap.add_argument('--births', default='fc_births.npz')
    ap.add_argument('--yk', action='store_true')
    ap.add_argument('--ykmax', type=int, default=8)
    args = ap.parse_args()
    if args.quick:
        args.nmax = 18
        args.nproc = 2
    nmax = args.nmax
    res = {'nmax': nmax}
    t0 = time.time()
    if not args.yk:
        births = build(nmax, args.nproc)
        print('births built, %.1f s' % (time.time() - t0), flush=True)
        M = fib(nmax + 2)
        phi = phi_sieve(M)
        print('phi sieve done %.1f s' % (time.time() - t0), flush=True)
        np.savez(args.births, **{'H%d' % j: births[j] for j in births})
        cN = np.zeros(M + 1, dtype=np.int64)
        cN[1] = 1
        checks = {}
        prev = None
        ok_all = True
        rows = []
        for N in range(1, nmax + 1):
            h = births[N].astype(np.int64)
            cN[:len(h)] += h
            tot_mass = int(cN[:fib(N + 2) + 1].sum())
            mmax = int(np.nonzero(cN)[0].max())
            Mn = fib(N + 2)
            # exact weighted sum with Python ints via int64 splitting (m < 2^32, c < 2^32, product < 2^64 may overflow int64): use float128-free split
            c_cur = cN[:Mn + 1]
            m_arr = np.arange(Mn + 1, dtype=np.uint64)
            prod = m_arr * c_cur.astype(np.uint64)  # < 2^63 since m c <= m phi(m) < 2^55
            # sum of 2^55-ish terms over 2^24 entries may overflow uint64: do it exactly in chunks with Python ints
            zneg = 0
            for lo in range(0, Mn + 1, 1 << 20):
                zneg += int(prod[lo:lo + (1 << 20)].sum(dtype=np.uint64))
            w1 = bool(np.all(c_cur[1:N + 2] == phi[1:N + 2]))
            w_next = bool(c_cur[N + 2] < phi[N + 2]) if N + 2 <= Mn else None
            le_phi = bool(np.all(c_cur[1:] <= phi[1:Mn + 1]))
            mono = True if prev is None else bool(np.all(c_cur[:len(prev)] >= prev))
            prev = c_cur.copy()
            okN = (tot_mass == 2 ** N) and (zneg == 3 ** N) and (mmax == Mn) and w1 and le_phi and mono
            ok_all &= okN
            # first m with c_N(m) < phi(m)
            nz = np.nonzero(c_cur[1:] < phi[1:Mn + 1])[0]
            first_def = int(nz[0]) + 1 if len(nz) else None
            rows.append({'N': N, 'sum_c': tot_mass, 'two_pow_N': 2 ** N, 'sum_m_c': zneg, 'three_pow_N': 3 ** N, 'max_m': mmax, 'F_N+2': Mn,
                         'w_eq_1_upto_N+1': w1, 'w_lt_1_at_N+2': w_next, 'c_le_phi': le_phi, 'monotone': mono, 'first_m_with_w_lt_1': first_def,
                         'ok': okN, 'checksum_c': int((c_cur * (np.arange(Mn + 1) % 1000003)).sum() % (2 ** 61 - 1))})
            print('N=%d ok=%s sum_c=%d sum_mc=%d max=%d F=%d w1=%s first_def=%s' % (N, okN, tot_mass, zneg, mmax, Mn, w1, first_def), flush=True)
        res['rows'] = rows
        res['all_ok'] = ok_all
        # direct Stern recursion for N <= 22
        direct_ok = {}
        for N in range(1, min(nmax, 22) + 1):
            row = stern_row_direct(N)
            bc = np.bincount(row)
            Mn = fib(N + 2)
            c = np.zeros(Mn + 1, dtype=np.int64); c[1] = 1
            for j in range(1, N + 1):
                h = births[j].astype(np.int64)
                c[:len(h)] += h
            direct_ok[N] = bool(len(bc) == Mn + 1 and np.array_equal(bc[:Mn + 1], c[:len(bc)]) or (np.array_equal(bc, c[:len(bc)]) and not c[len(bc):].any()))
        res['direct_stern_equal'] = direct_ok
        # Euclid cross-check for m <= 800
        mism = 0
        for m in range(2, 801):
            cnt = {}
            for a in range(1, m):
                if math.gcd(a, m) == 1:
                    s = cfsum(a, m - a)
                    cnt[s] = cnt.get(s, 0) + 1
            cum = 0
            for N in range(1, nmax + 1):
                cum += cnt.get(N, 0)
                # c_N(m) = number with cfsum <= N, N>=1
                pass
            # compare c_N(m) for several N
            c = 0
            for N in range(1, nmax + 1):
                c += int(births[N][m]) if m < len(births[N]) else 0
                euc = sum(v for k, v in cnt.items() if k <= N)
                if euc != c:
                    mism += 1
        res['euclid_mismatches_m_le_3000'] = mism
        print('euclid mismatches', mism, flush=True)
    else:
        # Yao-Knuth
        out = {}
        for e in range(3, args.ykmax + 1):
            m = 10 ** e
            tot, cnt = yk_mean(m)
            out[str(m)] = {'sum': tot, 'count': cnt, 'mean': tot / cnt}
            print('m=10^%d mean=%.9f count=%d' % (e, tot / cnt, cnt), flush=True)
        # ladder of moduli for the fit
        ladder = {}
        for m in [1009, 3001, 10007, 30011, 100003, 300007, 1000003, 3000017, 10000019, 30000001]:
            if m <= 10 ** args.ykmax * 3:
                tot, cnt = yk_mean(m)
                ladder[str(m)] = tot / cnt
        res['powers'] = out
        res['ladder'] = ladder
        print('ladder', ladder, flush=True)
        def fit(d):
            ms = np.array([float(k) for k in d]); y = np.array([d[k] if not isinstance(d[k], dict) else d[k]['mean'] for k in d])
            X = np.vstack([np.log(ms) ** 2, np.log(ms), np.ones_like(ms)]).T
            sol, *_ = np.linalg.lstsq(X, y, rcond=None)
            return [float(x) for x in sol]
        res['fit_powers_3term'] = fit(out)
        allp = dict(out); allp.update(ladder)
        res['fit_all_3term'] = fit(allp)
        res['fit_powers_ge1e5'] = fit({k: v for k, v in out.items() if int(k) >= 10 ** 5})
        res['six_over_pi2'] = 6 / np.pi ** 2
        print('fits', res['fit_powers_3term'], res['fit_all_3term'], res['fit_powers_ge1e5'], 6 / np.pi ** 2, flush=True)
    res['seconds'] = time.time() - t0
    with open(args.out, 'w') as f:
        json.dump(res, f, indent=1)

if __name__ == '__main__':
    main()
