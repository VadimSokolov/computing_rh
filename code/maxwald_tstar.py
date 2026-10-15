# maxwald_tstar.py: zeros of the symmetric finite stages F_N^sym on and off the critical line, for Table mw:tab:Tstar of
# Section sec:maxwald (ch/maxwald3.tex), subsection "The symmetric stages without poles, and how far they stay on the line".
# Added for the book in October 2026 from the checks of the authors' chapter, misc/reviews/nick-maxwald/ (review-3 items
# M4 to M7: poles of F_N, 72 terms at N = 8, the counts of N = 16 above t = 900, the N = 64 row). With P_N(q) = prod_{n <= N}
# (1 - q^n) prod_{m <= N} (1 - q^(2m-1)) = sum p_m q^m (exact integers), D_N(w) = sum_{m >= 1} p_m m^(-w) and the gamma factor
# A(s) = -(s - 1)/2 Gamma(1 + s/2) pi^(-s/2), the stage is F_N(s) = A(s) D_N(s/2) and F_N^sym(s) = (F_N(s) + F_N(1 - s))/2
# = A(s) H(s)/2 with H(s) = D_N(s/2) + (A(1 - s)/A(s)) D_N((1 - s)/2). A is analytic and zero free for Im s > 0, and
# F_N^sym(1/2 + it) = Re[A D_N] is real, so for each window (t1, t2] of height 10:
#  line:  the zeros on the critical line are the sign changes of Re[e^(i arg A) D_N(1/4 + it/2)] on a grid aligned with the
#         window edges (step 0.005, or 0.02 in tier B), so a zero just above an edge belongs to the window above it;
#  rect:  the zeros in the rectangle [-1, 2] x [t1, t2] are the winding number of H round its boundary (arguments of the
#         ratios of consecutive values, every step above 1 radian bisected), the number the chapter counts;
#  off line zeros = rect - line, an even number in every window because they come in pairs sigma + it, 1 - sigma + it.
# The same winding is taken round [0, 1] x [t1, t2] (rect_crit), the critical strip, which is where Corollary mw:cor:hurwitz
# (c) needs zero freeness off the line: the chapter's rectangle also holds the zeros with Re s in [-1, 0) and (1, 2].
# Nothing underflows: the factor |A| of size e^(-pi t/4) is never formed, only log A (scipy loggamma), and the windings use H.
# Tier A evaluates D_N in double precision (numpy); its absolute error is about 1e-15 S with S = sum |p_m| m^(-1/4) (2.4e9 at
# N = 64, 7e13 at N = 96, 2e18 at N = 128), so each window reports log10 of max|D_N| over that error and a flag when it is
# below 3; tier B (mpmath, 30 digits beyond S) is used for the low windows of N = 96 and 128, where D_N is a difference of
# huge terms. The first off line window of each N is then searched for its off line pair: local minima of |H| for
# 1/2 < Re s <= 2, polished by the secant method on F_N^sym in mpmath, and the zero of F_N^sym on the line nearest 280.004 at
# N = 48 (a window edge case of the chapter) is polished too.
# Usage (each in its own folder, one task of the table per sbatch array index; the index is mapped explicitly):
#   python3 maxwald_tstar.py --list                         prints the task table
#   python3 maxwald_tstar.py --task K                       runs task K (N, lo, hi, tier), writes maxwald_tstar_N<N>_<lo>_<hi>_<tier>.json
#   python3 maxwald_tstar.py --block N lo hi [tier [hl hr [siglo sighi]]]  one block by hand (tier A or B, rectangle sigma range)
#   python3 maxwald_tstar.py --merge                        reads the shards, locates the pairs, writes maxwald_tstar.json
#   python3 maxwald_tstar.py --locate N lo hi               only the pair search in one window
import glob
import json
import math
import sys
import time

import numpy as np
import mpmath as mp
from scipy.special import loggamma

LOGPI = math.log(math.pi)
TWO_PI = 2 * math.pi
SIG_LO, SIG_HI = -1.0, 2.0
T0 = 0.02
CHAPTER_HEIGHT = {8: 400, 16: 1000, 32: 700, 48: 900, 64: 400}


def pn_coeffs(N):
    """Exact integer coefficients of P_N, degree N(3N+1)/2."""
    deg = N * (3 * N + 1) // 2
    p = [0] * (deg + 1)
    p[0] = 1
    cur = 0
    for a in list(range(1, N + 1)) + list(range(1, 2 * N, 2)):
        for m in range(cur + a, a - 1, -1):
            p[m] -= p[m - a]
        cur += a
    assert cur == deg
    return p


class Stage:
    def __init__(self, N):
        self.N = N
        p = pn_coeffs(N)
        idx = [m for m in range(1, len(p)) if p[m] != 0]
        self.terms = len(idx)
        self.m = np.array(idx, dtype=np.int64)
        self.cint = [p[m] for m in idx]
        self.c = np.array([float(x) for x in self.cint])
        self.logm = np.log(self.m.astype(float))
        self.S14 = float(np.sum(np.abs(self.c) * self.m ** -0.25))
        self.dps = 30 + int(math.ceil(math.log10(max(self.S14, 10.0))))
        self._L = None

    def D_np(self, w):
        w = np.atleast_1d(np.asarray(w, dtype=complex))
        out = np.empty(w.shape, dtype=complex)
        B = max(1, int(2_000_000 // len(self.c)))
        for i in range(0, len(w), B):
            out[i:i + B] = np.exp(-np.multiply.outer(w[i:i + B], self.logm)) @ self.c
        return out

    def D_mp(self, w):
        mp.mp.dps = self.dps
        if self._L is None:
            self._L = [mp.log(int(m)) for m in self.m]
        w = mp.mpmathify(w)
        acc = mp.mpc(0)
        for c, L in zip(self.cint, self._L):
            acc += c * mp.exp(-w * L)
        return acc

    def D(self, w, tier):
        w = np.atleast_1d(w)
        if tier == 'A':
            return self.D_np(w)
        return np.array([complex(self.D_mp(mp.mpc(z.real, z.imag))) for z in w])


def logA(s):
    s = np.asarray(s, dtype=complex)
    return np.log(-(s - 1) / 2) + loggamma(1 + s / 2) - (s / 2) * LOGPI


def make_H(st, tier):
    def H(s):
        s = np.atleast_1d(np.asarray(s, dtype=complex))
        R = np.exp(logA(1 - s) - logA(s))
        return st.D(s / 2, tier) + R * st.D((1 - s) / 2, tier)
    return H


def boundary(sig_lo, sig_hi, t1, t2, h):
    c = [complex(sig_lo, t1), complex(sig_hi, t1), complex(sig_hi, t2), complex(sig_lo, t2)]
    pts = []
    for k in range(4):
        a, b = c[k], c[(k + 1) % 4]
        n = max(2, int(math.ceil(abs(b - a) / h)))
        pts.append(a + (b - a) * np.arange(n) / n)
    pts.append(np.array([c[0]]))
    return np.concatenate(pts)


def winding(Hf, pts, thr=1.0, maxit=12):
    P = pts.copy()
    V = Hf(P)
    for _ in range(maxit):
        inc = np.angle(V[1:] / V[:-1])
        bad = np.nonzero(np.abs(inc) > thr)[0]
        if len(bad) == 0:
            break
        mids = 0.5 * (P[bad] + P[bad + 1])
        P = np.insert(P, bad + 1, mids)
        V = np.insert(V, bad + 1, Hf(mids))
    inc = np.angle(V[1:] / V[:-1])
    tot = inc.sum() / TWO_PI
    return int(round(tot)), float(abs(tot - round(tot))), int(np.sum(np.abs(inc) > thr)), len(P)


def window(st, lo, hi, tier, hl, hr, sig=(SIG_LO, SIG_HI)):
    t1 = max(lo, T0)
    K = int(round((hi - t1) / hl))
    ts = t1 + hl * np.arange(K + 1)
    s = 0.5 + 1j * ts
    D = st.D(s / 2, tier)
    lv = (np.exp(1j * logA(s).imag) * D).real
    sg = np.sign(lv)
    chg = np.nonzero(sg[1:] != sg[:-1])[0]
    on = len(chg)
    zz = ts[chg] + hl * lv[chg] / (lv[chg] - lv[chg + 1])
    edge = [round(float(z), 6) for z in zz if abs(z - lo) < 0.02 or abs(z - hi) < 0.02]
    eps = 2.2e-16 * 4 * st.S14
    lf = np.log10(np.abs(lv) + 1e-300) + logA(s).real / math.log(10)
    Hf = make_H(st, tier)
    w, resid, unres, npts = winding(Hf, boundary(sig[0], sig[1], t1, hi, hr))
    wc = w if tuple(sig) == (0.0, 1.0) else winding(Hf, boundary(0.0, 1.0, t1, hi, hr))[0]
    return {'lo': lo, 'hi': hi, 'rect': w, 'line': on, 'off': w - on, 'rect_crit': wc, 'off_crit': wc - on,
            'winding_resid': resid,
            'unresolved_steps': unres, 'boundary_points': npts, 'edge_line_zeros': edge,
            'log10_margin': float(math.log10(max(np.max(np.abs(D)), 1e-300) / eps)),
            'log10_absF_median': float(np.median(lf)), 'log10_absF_min_grid': float(np.min(lf))}


def run_block(N, lo, hi, tier='A', hl=None, hr=None, sig=(SIG_LO, SIG_HI)):
    hl = hl or (0.005 if tier == 'A' else 0.02)
    hr = hr or (0.02 if tier == 'A' else 0.1)
    t0 = time.time()
    st = Stage(N)
    wins = []
    t = float(lo)
    print(f'== N={N}: {st.terms} terms, S={st.S14:.3g}, tier {tier}, heights {lo} to {hi}', flush=True)
    while t < hi - 1e-9:
        w = window(st, t, min(t + 10.0, hi), tier, hl, hr, sig)
        wins.append(w)
        flag = ''
        if w['off']:
            flag += f"   OFF LINE: {w['off']}"
        flag += f"   (in 0 <= Re s <= 1: {w['rect_crit']}, off line {w['off_crit']})" if w['off_crit'] != w['off'] or w['off'] else ''
        if w['unresolved_steps'] or w['winding_resid'] > 1e-3:
            flag += '   WINDING UNRESOLVED'
        if tier == 'A' and w['log10_margin'] < 3:
            flag += '   PRECISION FLAG'
        print(f"  ({w['lo']:7.1f},{w['hi']:7.1f}]: rectangle {w['rect']:3d}  line {w['line']:3d}{flag}", flush=True)
        t += 10.0
    out = {'N': N, 'lo': lo, 'hi': hi, 'tier': tier, 'terms': st.terms, 'S14': st.S14, 'step_line': hl,
           'step_rect': hr, 'sigma_range': list(sig), 'windows': wins, 'seconds': round(time.time() - t0, 1)}
    json.dump(out, open(f'maxwald_tstar_N{N}_{int(lo)}_{int(hi)}_{tier}.json', 'w'))
    return out


def tasks():
    T = [(8, 0, 400, 'A')]
    T += [(16, a, a + 250, 'A') for a in range(0, 1000, 250)]
    T += [(32, a, a + 200, 'A') for a in range(0, 800, 200)]
    T += [(48, a, a + 200, 'A') for a in range(0, 1000, 200)]
    T += [(64, a, a + 250, 'A') for a in range(0, 1500, 250)]
    T += [(96, a, a + 10, 'B') for a in range(0, 100, 10)] + [(96, a, a + 200, 'A') for a in range(100, 2500, 200)]
    T += [(128, a, a + 10, 'B') for a in range(0, 130, 10)] + [(128, a, a + 200, 'A') for a in range(130, 3130, 200)]
    return T


def Hmp(st, s):
    mp.mp.dps = st.dps
    s = mp.mpmathify(s)

    def lA(z):
        return mp.log(-(z - 1) / 2) + mp.loggamma(1 + z / 2) - (z / 2) * mp.log(mp.pi)
    R = mp.exp(lA(1 - s) - lA(s))
    return st.D_mp(s / 2) + R * st.D_mp((1 - s) / 2)


def locate(N, lo, hi, nmin=12):
    """Off line zeros of F_N^sym in the window (lo, hi] with Re s > 1/2: minima of |H| on a grid, polished in mpmath."""
    st = Stage(N)
    Hf = make_H(st, 'A')
    sig = np.arange(0.5, 2.0 + 1e-9, 0.025)
    ts = np.arange(max(lo, T0), hi + 1e-9, 0.05)
    S, T = np.meshgrid(sig, ts, indexing='ij')
    Z = np.abs(Hf((S + 1j * T).ravel())).reshape(S.shape)
    cands = []
    for i in range(len(sig)):
        for j in range(len(ts)):
            nb = Z[max(i - 1, 0):i + 2, max(j - 1, 0):j + 2]
            if Z[i, j] <= nb.min() and sig[i] > 0.5:
                cands.append((Z[i, j], sig[i], ts[j]))
    cands.sort()
    found = []
    for _, s0, t0 in cands[:nmin]:
        try:
            z = mp.findroot(lambda x: Hmp(st, x), (mp.mpc(s0, t0), mp.mpc(s0 + 1e-3, t0 + 1e-3)), solver='secant',
                            tol=1e-18, maxsteps=40, verify=False)
        except Exception:
            continue
        mp.mp.dps = st.dps
        zr, zi = float(z.real), float(z.imag)
        res = abs(Hmp(st, z)) / (abs(st.D_mp(z / 2)) + 1e-300)
        if res > 1e-12 or not (lo - 1e-9 <= zi <= hi + 1e-9) or abs(zr - 0.5) < 1e-8 or zr > SIG_HI + 1e-9:
            continue
        if any(abs(zr - a) + abs(zi - b) < 1e-6 for a, b, _ in found):
            continue
        mirror = abs(Hmp(st, 1 - mp.conj(z))) / (abs(st.D_mp((1 - mp.conj(z)) / 2)) + 1e-300)
        found.append((zr, zi, float(mirror)))
    found.sort(key=lambda x: x[1])
    return [{'re': a, 'im': b, 're_minus_half': a - 0.5, 'mirror_relative_residual': c} for a, b, c in found]


def line_zero_near(N, t0):
    st = Stage(N)
    mp.mp.dps = st.dps

    def f(t):
        s = mp.mpc(0.5, t)
        lA = mp.log(-(s - 1) / 2) + mp.loggamma(1 + s / 2) - (s / 2) * mp.log(mp.pi)
        return (mp.exp(1j * lA.imag) * st.D_mp(s / 2)).real
    return float(mp.findroot(f, (t0 - 1e-3, t0 + 1e-3), solver='secant', tol=1e-20, maxsteps=40, verify=False))


def merge():
    shards = [json.load(open(f)) for f in sorted(glob.glob('maxwald_tstar_N*_*_*_*.json'))]
    byN = {}
    for sh in shards:
        d = byN.setdefault(sh['N'], {'terms': sh['terms'], 'S14': sh['S14'], 'win': {}})
        for w in sh['windows']:
            if w['lo'] not in d['win'] or sh['tier'] == 'B':
                d['win'][w['lo']] = dict(w, tier=sh['tier'])
    res = {'working_precision': 'double (numpy) in tier A, mpmath at 30 digits beyond log10 S in tier B and for the pair search',
           'window_height': 10.0, 'sigma_range': [SIG_LO, SIG_HI]}
    for N in sorted(byN):
        d = byN[N]
        wl = [d['win'][k] for k in sorted(d['win'])]
        contiguous = all(abs(wl[i]['hi'] - wl[i + 1]['lo']) < 1e-9 for i in range(len(wl) - 1))
        top = wl[-1]['hi']
        e = {'terms_in_F_N': d['terms'], 'S_sum_abs_p_m_m_pow_minus_quarter': d['S14'], 'height_covered': top,
             'windows_contiguous_from_0': bool(contiguous and wl[0]['lo'] == 0),
             'windows': [[w['lo'], w['rect'], w['line']] for w in wl],
             'odd_off_windows': [w['lo'] for w in wl if w['off'] % 2],
             'unresolved_windows': [w['lo'] for w in wl if w['unresolved_steps'] or w['winding_resid'] > 1e-3],
             'precision_flagged_windows': [w['lo'] for w in wl if w['tier'] == 'A' and w['log10_margin'] < 3],
             'min_log10_margin_tier_A': min([w['log10_margin'] for w in wl if w['tier'] == 'A'] or [None])}

        def upto(h):
            sel = [w for w in wl if w['hi'] <= h + 1e-9]
            return {'height': h, 'rect': sum(w['rect'] for w in sel), 'line': sum(w['line'] for w in sel),
                    'complete': bool(sel and sel[-1]['hi'] >= h - 1e-9)}
        e['totals_to_height'] = {str(h): upto(h) for h in sorted({CHAPTER_HEIGHT.get(N, top), 900, top} & set(range(0, int(top) + 1, 10)))}
        off = [w for w in wl if w['off'] > 0]
        if off:
            f = off[0]
            e['first_off_line_window'] = [f['lo'], f['hi']]
            e['first_off_line_window_20_grid'] = [20 * math.floor(f['lo'] / 20), 20 * math.floor(f['lo'] / 20) + 20]
            e['T_star'] = f['lo']
            e['T_star_over_N'] = f['lo'] / N
            e['off_counts_first_windows'] = [[w['lo'], w['off']] for w in off[:6]]
            print(f'N={N}: first off line window ({f["lo"]}, {f["hi"]}], locating the pair', flush=True)
            e['first_pair'] = locate(N, f['lo'], f['hi'])
        else:
            e['first_off_line_window'] = None
            e['T_star'] = f'none below {top}'
        offc = [w for w in wl if w['off_crit'] > 0]
        e['T_star_critical_strip_0_1'] = offc[0]['lo'] if offc else f'none below {top}'
        e['off_line_zeros_in_critical_strip_per_window'] = [[w['lo'], w['off_crit']] for w in offc[:6]]
        e['windows_crit'] = [[w['lo'], w['rect_crit'], w['line']] for w in wl]
        if N == 48:
            e['line_zero_near_280_004'] = line_zero_near(48, 280.004)
        res[f'N{N}'] = e
    json.dump(res, open('maxwald_tstar.json', 'w'), indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else {a: b for a, b in v.items() if a not in ('windows',)}) for k, v in res.items()}, indent=1)[:6000])


if __name__ == '__main__':
    a = sys.argv[1:]
    if not a or a[0] in ('-h', '--help'):
        print(open(__file__).read().split('import glob')[0])
    elif a[0] == '--list':
        for i, t in enumerate(tasks()):
            print(i, *t)
        print(len(tasks()), 'tasks')
    elif a[0] == '--task':
        N, lo, hi, tier = tasks()[int(a[1])]
        run_block(N, lo, hi, tier)
    elif a[0] == '--block':
        run_block(int(a[1]), float(a[2]), float(a[3]), a[4] if len(a) > 4 else 'A',
                  float(a[5]) if len(a) > 5 else None, float(a[6]) if len(a) > 6 else None,
                  (float(a[7]), float(a[8])) if len(a) > 8 else (SIG_LO, SIG_HI))
    elif a[0] == '--merge':
        merge()
    elif a[0] == '--locate':
        print(json.dumps(locate(int(a[1]), float(a[2]), float(a[3])), indent=1))
