# px_dyck.py: PX.4 of the reconciliation of the notes on the pole cancellation and the excursion (Chapter 25 on extreme values, Section
# ev:sec:decomp, the Dyck chain): the zeros of the Dirichlet polynomial F_n(s) = E[Y_n^s] = sum_h P(H_n = h) (sqrt(2/pi) h / sqrt(2n))^s,
# H_n the height (maximal level) of a uniform Dyck path of length 2n, for n = 800 ... 25600 (n = 1600, 6400, 25600 are the required ones).
# Exact law. The number of Dyck paths of length 2n and height at most h is N(h) = sum_k [C(2n, n + k(h+2)) - C(2n, n - 1 + k(h+2))] (reflection principle,
# barriers at -1 and h + 1), out of Cat(n) = C(2n, n)/(n+1); P(H_n = h) = (N(h) - N(h-1))/Cat(n) is computed in Python integers (checked against a
# direct dynamic programme for n <= 12 with --selftest), and converted to balls at the working precision; the heights run up to the first h with
# 1 - N(h)/Cat(n) < 1e-62 and those with P(H_n = h) < 1e-62 are dropped (their contribution to F_n on the contours is below 1e-50).
# Census. The argument principle on the boundary of -4 < Re s < 6, 0 < Im s < 60 and on the sub-boxes (heights 0-20, 20-40, 40-60; Re s < 1/2 and
# Re s > 1/2), by px_lib.sampled_winding at 256 bits (77 digits) and again at 192 bits (computed: adaptive sampling, not certified; the minimum of
# |F_n| on the contour is logged). F_n is real on the real axis and positive there (E Y^s > 0), so the bottom edge Im s = 0 carries no zero.
# Zeros near 1/2 + i gamma_k (k <= 13). Newton at 500 bits (central difference, h = 1e-30) from 1/2 + i gamma_k (gamma_k from data/zeros_1700.txt); if that does not
# converge to a zero within 1.5 of the start or two starts converge to the same zero, continuation from the next larger n is used and logged;
# each zero is confirmed alone by the argument principle on the square of half side 0.4 about it (count 1), and the real part is compared
# with 1/2. Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 px_dyck.py [--quick] [--selftest] [--out FILE] [zeros_file]; --quick: n = 400 and 1600 only, census on one box.
import argparse, json, sys, time, math
from concurrent.futures import ProcessPoolExecutor
from flint import acb, arb, ctx, fmpz
sys.path.insert(0, '.')
import px_lib as L


def dyck_counts(n, tail=10 ** 62):
    two = 2 * n
    row = [1] * (n + 1)
    for j in range(n):
        row[j + 1] = row[j] * (two - j) // (j + 1)

    def C(j):
        if j < 0 or j > two:
            return 0
        return row[j] if j <= n else row[two - j]
    cat = row[n] // (n + 1)

    def Nle(h):
        m = h + 2
        s = 0
        kmax = n // m + 1
        for k in range(-kmax, kmax + 1):
            s += C(n + k * m) - C(n - 1 + k * m)
        return s
    D = {}
    prev = Nle(0)
    h = 0
    while True:
        h += 1
        cur = Nle(h)
        D[h] = cur - prev
        prev = cur
        if (cat - cur) * tail < cat or h > n:
            break
    assert cur == cat or (cat - cur) * tail < cat
    return cat, D


def brute(n):
    # DP over paths: count by maximal height
    from collections import defaultdict
    st = {(0, 0): 1}
    for _ in range(2 * n):
        nx = defaultdict(int)
        for (x, m), c in st.items():
            nx[(x + 1, max(m, x + 1))] += c
            if x > 0:
                nx[(x - 1, m)] += c
        st = nx
    res = defaultdict(int)
    for (x, m), c in st.items():
        if x == 0:
            res[m] += c
    return res


class Dyck:
    def __init__(self, n, prec, cat=None, D=None):
        ctx.prec = prec
        self.n, self.prec = n, prec
        if D is None:
            cat, D = dyck_counts(n)
        self.cat = cat
        catb = arb(fmpz(cat))
        c0 = (arb(2) / arb.pi()).log() / 2 - (arb(2 * n)).log() / 2
        self.P, self.Lg = [], []
        self.nterms = 0
        for h, d in sorted(D.items()):
            if d * 10 ** 62 < cat:
                continue
            self.P.append(arb(fmpz(d)) / catb)
            self.Lg.append(c0 + arb(h).log())
        self.nterms = len(self.P)

    def F(self, s):
        ctx.prec = self.prec
        tot = acb(0)
        for p, l in zip(self.P, self.Lg):
            tot += p * (s * l).exp()
        return tot


def census(n, prec, boxes):
    d = Dyck(n, prec)
    res = {}
    for name, box in boxes.items():
        ctx.prec = prec
        t = time.time()
        c, tot, mn, ne = L.sampled_winding(d.F, box)
        res[name] = dict(box=box, count=c, total=tot, minF=mn, evals=ne, sec=time.time() - t)
    return n, prec, d.nterms, res


def zeros_for(n, gam, starts=None, prec=500):
    d = Dyck(n, prec)
    out = []
    for k in range(1, 14):
        z0 = complex(0.5, gam[k - 1]) if starts is None else starts[k - 1]
        z, r, it = L.newton(d.F, z0, prec, maxstep=0.5)
        ctx.prec = 256
        d2 = Dyck(n, 256) if False else d
        ctx.prec = prec
        box = (float(z.real) - 0.4, float(z.real) + 0.4, float(z.imag) - 0.4, float(z.imag) + 0.4)
        ctx.prec = 256
        d256 = Dyck(n, 256, d.cat, None) if False else None
        out.append(dict(k=k, re=z.real.str(32, radius=False), im=z.imag.str(32, radius=False), resid=r, steps=it,
                        dist=float(abs(z - acb(arb(0.5), arb(gam[k - 1])))), left=bool(z.real < arb(1) / 2), x=float(z.real), y=float(z.imag)))
    return n, out, d


def square_check(n, zs):
    d = Dyck(n, 256)
    res = []
    for z in zs:
        box = (z['x'] - 0.4, z['x'] + 0.4, z['y'] - 0.4, z['y'] + 0.4)
        ctx.prec = 256
        c, tot, mn, ne = L.sampled_winding(d.F, box)
        res.append((z['k'], c, mn))
    return n, res


def run_n(n, gam, starts=None):
    t = time.time()
    n_, out, d = zeros_for(n, gam, starts)
    # confirm alone
    n_, sq = square_check(n, out)
    for z, (k, c, mn) in zip(out, sq):
        z['square_count'] = c
        z['square_min'] = mn
    return n, out, time.time() - t, d.nterms


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--out', default='px_dyck.json')
    ap.add_argument('zeros', nargs='?', default='zeros_1700.txt')
    args = ap.parse_args()
    if args.selftest or args.quick:
        for n in (1, 2, 3, 5, 8, 12):
            cat, D = dyck_counts(n, tail=10 ** 62)
            b = brute(n)
            ok = all(D.get(h, 0) == b.get(h, 0) for h in range(1, 2 * n + 2)) and sum(D.values()) == cat
            print('selftest n=%d Cat=%d dyck_counts vs DP: %s' % (n, cat, ok), flush=True)
            assert ok
        if args.selftest:
            sys.exit(0)
    gam = [float(l.split()[-1]) for l in open(args.zeros) if l.strip()]
    out = {}
    ns = [400, 1600] if args.quick else [800, 1600, 3200, 6400, 12800, 25600]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=8) as ex:
        # zeros by Newton from 1/2 + i gamma_k at every n
        fz = {n: ex.submit(run_n, n, gam) for n in ns}
        # census at the required n
        req = [n for n in ns if n in (1600, 6400, 25600)] if not args.quick else [1600]
        boxfull = {'full': (-4.0, 6.0, 0.0, 60.0), 'h0_20': (-4.0, 6.0, 0.0, 20.0), 'h20_40': (-4.0, 6.0, 20.0, 40.0), 'h40_60': (-4.0, 6.0, 40.0, 60.0),
                   'left_of_half': (-4.0, 0.5, 0.0, 60.0), 'right_of_half': (0.5, 6.0, 0.0, 60.0)}
        if args.quick:
            boxfull = {'h0_20': (-4.0, 6.0, 0.0, 20.0), 'right_of_half': (0.5, 6.0, 0.0, 20.0)}
        fc = [ex.submit(census, n, p, boxfull) for n in req for p in (256, 192)]
        census_res = []
        for f in fc:
            n, prec, nt, res = f.result()
            census_res.append(dict(n=n, prec=prec, nterms=nt, res=res))
            for name, r in res.items():
                print('CENSUS n=%d prec=%d %s count=%d total=%.6f minF=%.3g evals=%d sec=%.0f' % (n, prec, name, r['count'], r['total'], r['minF'], r['evals'], r['sec']), flush=True)
        zres = {}
        for n in ns:
            n_, zs, sec, nt = fz[n].result()
            zres[n] = zs
            print('DIRECT ZEROS n=%d (%d terms, %.0f s)' % (n, nt, sec), flush=True)
            for z in zs:
                print('   k=%2d %s %s dist=%.4f left=%s resid=%.2g steps=%d square=%d (min %.2g)' % (z['k'], z['re'][:12], z['im'][:12], z['dist'], z['left'], z['resid'], z['steps'], z['square_count'], z['square_min']), flush=True)
    # continuation in n from the largest n: the start of the zero k at n is the zero k found at the next larger n
    chain = {}
    prev = None
    for n in sorted(ns, reverse=True):
        starts = None if prev is None else [complex(z['x'], z['y']) for z in prev]
        n_, zs, sec, nt = run_n(n, gam, starts)
        chain[n] = zs
        prev = zs
        distinct = len(set((z['re'][:20], z['im'][:20]) for z in zs))
        print('CHAIN ZEROS n=%d (%d terms, %.0f s): distinct %d of 13, all squares count 1: %s, all left: %s' % (n, nt, sec, distinct, all(z['square_count'] == 1 for z in zs), all(z['left'] for z in zs)), flush=True)
        for z in zs:
            print('   k=%2d %s %s dist=%.4f left=%s resid=%.2g steps=%d square=%d (min %.2g)' % (z['k'], z['re'][:12], z['im'][:12], z['dist'], z['left'], z['resid'], z['steps'], z['square_count'], z['square_min']), flush=True)
    out['census'] = census_res
    out['zeros_direct'] = {str(n): v for n, v in zres.items()}
    out['zeros_chain'] = {str(n): v for n, v in chain.items()}
    json.dump(out, open(args.out, 'w'), indent=1, default=str)
    print('total sec', time.time() - t0, flush=True)
