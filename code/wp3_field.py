# wp3_field.py: work package 3 (run book rh_compute.tex, Section sec:spin; chapter ch/nc_cauchy.tex, Theorem cm:thm:exact, Proposition cm:prop:gap).
# The field phi of all other zeros at a hypothetical pair of zeros near the ordinate gamma_k, for the first 10^4 ordinates, in ball arithmetic.
# Conventions (those of nc_cauchy_gap_theory.py and nc_cauchy_field_zero.py): F(z) = Xi(z) (z^2-c^2)(z^2-cbar^2)/((z^2-g_k^2)(z^2-g_{k+1}^2)) is Xi with the real zeros +-g_k, +-g_{k+1}
# replaced by the pair c = m + i eta and its mirror images, m = (g_k+g_{k+1})/2, G = F/((z-c)(z-cbar)), and
#   phi_A(k, eta) = (G'/G)(m) = Xi'/Xi(m) - sum_{a in {g_k, g_{k+1}}} 2m/(m^2-a^2) + 4m/(4m^2+eta^2)         (field at the midpoint pair, the book's definition)
#   phi_B(k)      = the regular part of Xi'/Xi at g_k, = lim_{z->g_k} (Xi'/Xi(z) - 1/(z-g_k)) = a2/a1 for Xi(g_k+x) = a1 x + a2 x^2 + ...   (field of all other zeros, mirror image included, at the ordinate itself)
# Xi'/Xi(x) = Re[ i dlogxi(1/2+ix) ], dlogxi(s) = 1/s + 1/(s-1) - (log pi)/2 + psi(s/2)/2 + zeta'/zeta(s), at the exact binary point x.
# The ordinates are the certified balls of acb.zeta_zero(s) (flint: rigorous enumeration and isolation), 128 bit; the removed ordinates enter only through the ball radius (about 1e-35).
# Also: zeros of phi_A(x, eta=0.1) as a function of the pair position x in (g_{k-1}, g_{k+2}) (by sign changes on a grid, then bisection; computed, not certified) with the slope there, and the share of the window where |phi_A| < 1e-2.
# Independent check (method 2): sum over data/zeros_8000.txt plus the integral of the smooth density for k <= 1000; and the ordinates against the 50 digit file.
# Usage: python3 wp3_field.py [--n 10000] [--workers 16] [--quick] [--out wp3_field.json] [zeros_8000.txt]
import argparse, json, math, sys, time
from multiprocessing import Pool
from flint import acb, acb_series, arb, ctx

PREC = 128
ETA = 0.1
FILE = []

def dlogxi(s):
    S = acb_series([s, 1], prec=2)
    co = S.zeta().coeffs()
    return 1 / s + 1 / (s - 1) - arb.pi().log() / 2 + (s / 2).digamma() / 2 + co[1] / co[0]

def xi_logder(x):
    """Xi'/Xi at the real ball x, as an acb whose imaginary part must contain 0"""
    return acb(0, 1) * dlogxi(acb(arb(0.5), x))

def phiA(x, ga, gb, eta=ETA):
    v = xi_logder(x)
    e = arb(eta)
    return v.real - 2 * x / (x * x - ga * ga) - 2 * x / (x * x - gb * gb) + 4 * x / (4 * x * x + e * e), v.imag

def reg_at(g):
    """regular part of Xi'/Xi at the zero ball g: a2/a1"""
    s = acb(arb(0.5), g)
    S = acb_series([s, acb(0, 1)], prec=3)
    two = acb(2)
    v = S * (S - 1) / 2 * (-(S / 2) * arb.pi().log()).exp() * (S / 2).gamma() * S.zeta()
    co = v.coeffs()
    return (co[2] / co[1]).real

def fl(a): return float(a.mid())

def work(args):
    k0, k1, nmax = args
    ctx.prec = PREC
    lo = max(1, k0 - 1); cnt = k1 + 2 - lo + 1
    zs = acb.zeta_zeros(lo, cnt)
    g = {lo + i: z.imag for i, z in enumerate(zs)}
    rows = []
    for k in range(k0, k1 + 1):
        gk, gk1 = g[k], g[k + 1]
        m = arb((gk + gk1).mid() / 2)
        pa, imag = phiA(m, gk, gk1)
        pa0 = pa - 4 * m / (4 * m * m + arb(ETA) ** 2) + 1 / m
        pb = reg_at(gk)
        row = dict(k=k, gamma=fl(gk), grad=float(gk.rad()), m=fl(m), phiA=fl(pa), phiA_rad=float(pa.rad()), phiA_eta0=fl(pa0), phiB=fl(pb), phiB_rad=float(pb.rad()),
                   imag_contains0=bool(imag.contains(0)))
        # ordinate against the 50 digit file
        if k <= len(FILE):
            ctx.prec = 300; d = arb(FILE[k - 1]) - gk; ctx.prec = PREC
            row['file_diff'] = float(abs(d.mid()))
            row['file_ok'] = row['file_diff'] <= float(gk.rad()) * (1 + 1e-9) + 1e-45
        # zeros of the field in the pair position x, window (g_{k-1}, g_{k+2})
        a = g.get(k - 1, arb(0)) if k >= 2 else arb(0.5)
        b = g[k + 2]
        a = float(a.mid()) + 0.02 if k >= 2 else 1.0
        b = float(b.mid()) - 0.02
        ng = 24
        xs = [a + (b - a) * i / ng for i in range(ng + 1)]
        def f(x): return fl(phiA(arb(x), gk, gk1)[0])
        vals = [f(x) for x in xs]
        roots = []
        for i in range(ng):
            if vals[i] * vals[i + 1] < 0:
                lo_, hi_, flo = xs[i], xs[i + 1], vals[i]
                for _ in range(45):
                    mid = (lo_ + hi_) / 2; fm = f(mid)
                    if fm * flo > 0: lo_, flo = mid, fm
                    else: hi_ = mid
                r = (lo_ + hi_) / 2; h = 1e-5
                sl = (f(r + h) - f(r - h)) / (2 * h)
                roots.append((r, sl))
        width = sum(2e-2 / abs(s) for r, s in roots)
        row.update(nroots=len(roots), roots=[r for r, s in roots], slopes=[s for r, s in roots], window=b - a, hide_width=width,
                   gap=fl(gk1 - gk))
        rows.append(row)
    return rows

def method2(rows, zfile):
    """phi_A by the sum over the 8000 ordinates of the file plus the smooth tail, for k <= 1000"""
    import mpmath as mp
    mp.mp.dps = 30
    G = [mp.mpf(l.split()[-1]) for l in open(zfile) if l.strip()]
    T = G[-1] + mp.pi / mp.log(G[-1] / (2 * mp.pi))
    out = []
    byk = {r['k']: r for r in rows}
    for k in [1, 2, 3, 5, 8, 12, 15, 30, 50, 100, 200, 500, 1000]:
        if k not in byk: continue
        m = (G[k - 1] + G[k]) / 2
        s = mp.fsum(2 * m / (m * m - G[j] ** 2) for j in range(len(G)) if j not in (k - 1, k))
        tail = mp.quad(lambda t: mp.log(t / (2 * mp.pi)) / (2 * mp.pi) * 2 * m / (m * m - t * t), [T, 10 * T, mp.inf])
        # Xi'/Xi total = sum over all zeros; the sum above is over zeros 1..8000 except k,k+1 (+ tail); add mirror term
        val = s + tail + 4 * m / (4 * m * m + mp.mpf(ETA) ** 2)
        out.append(dict(k=k, method2=float(val), method1=byk[k]['phiA'], diff=float(val - byk[k]['phiA'])))
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=10000)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='wp3_field.json')
    ap.add_argument('zeros', nargs='?', default='zeros_8000.txt')
    a = ap.parse_args()
    FILE.extend(l.split()[-1] for l in open(a.zeros) if l.strip())
    n = 40 if a.quick else a.n
    t0 = time.time()
    chunk = max(1, n // (a.workers * 4))
    jobs = [(s, min(n, s + chunk - 1), 0) for s in range(1, n + 1, chunk)]
    with Pool(a.workers) as p:
        res = p.map(work, jobs, chunksize=1)
    rows = sorted([r for rr in res for r in rr], key=lambda r: r['k'])
    print(f"{len(rows)} rows in {time.time()-t0:.0f} s", flush=True)
    import statistics as st
    edges_h = [0, 100, 300, 1000, 3000, 6000, 10000]
    edges_p = [0, 1e-3, 1e-2, 1e-1, 1, 3, 10, 1e9]
    summ = []
    def stats(rs, label):
        v = [abs(r['phiA']) for r in rs]; vb = [abs(r['phiB']) for r in rs]
        hist = [sum(1 for x in v if edges_p[i] <= x < edges_p[i + 1]) for i in range(len(edges_p) - 1)]
        histb = [sum(1 for x in vb if edges_p[i] <= x < edges_p[i + 1]) for i in range(len(edges_p) - 1)]
        hw = [r['hide_width'] / r['window'] for r in rs]
        d = dict(label=label, count=len(rs), hist_phiA=hist, hist_phiB=histb, median_A=st.median(v), mean_A=st.mean(v), max_A=max(v), min_A=min(v),
                 frac_A_lt_1e2=sum(1 for x in v if x < 1e-2) / len(v), frac_A_lt_1e1=sum(1 for x in v if x < 1e-1) / len(v),
                 median_B=st.median(vb), frac_B_lt_1e2=sum(1 for x in vb if x < 1e-2) / len(vb),
                 mean_hide_share=st.mean(hw), median_nroots=st.median([r['nroots'] for r in rs]), min_nroots=min(r['nroots'] for r in rs),
                 mean_logT=st.mean(math.log(r['gamma'] / (2 * math.pi)) for r in rs),
                 median_A_over_logT=st.median(abs(r['phiA']) / math.log(r['gamma'] / (2 * math.pi)) for r in rs))
        summ.append(d)
    stats(rows, 'all')
    for i in range(len(edges_h) - 1):
        rs = [r for r in rows if edges_h[i] <= r['gamma'] < edges_h[i + 1]]
        if rs: stats(rs, f"gamma in [{edges_h[i]},{edges_h[i+1]})")
    chk = dict(all_imag_contain0=all(r['imag_contains0'] for r in rows), max_phiA_rad=max(r['phiA_rad'] for r in rows), max_phiB_rad=max(r['phiB_rad'] for r in rows),
               max_gamma_rad=max(r['grad'] for r in rows), file_checked=sum(1 for r in rows if 'file_ok' in r), file_bad=sum(1 for r in rows if 'file_ok' in r and not r['file_ok']),
               file_maxdiff=max((r['file_diff'] for r in rows if 'file_diff' in r), default=None),
               gamma_last=rows[-1]['gamma'], windows_without_root=[r['k'] for r in rows if r['nroots'] == 0][:20], n_windows_without_root=sum(1 for r in rows if r['nroots'] == 0))
    m2 = method2(rows, a.zeros)
    json.dump(dict(summary=summ, checks=chk, method2=m2, rows=rows), open(a.out, 'w'))
    for s in summ: print(json.dumps(s))
    print(json.dumps(chk)); print(json.dumps(m2))
