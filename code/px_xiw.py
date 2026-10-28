# px_xiw.py: PX.3 of the reconciliation of the notes on the pole cancellation and the excursion (Chapter 22, the paragraph on the truncations
# after Theorem mi:thm:noproxy): the truncated Riemann functions Xi_w(t) = 1/2 - (1/2)(t^2 + 1/4) sum_{n<=N} (F_n(t) + F_n(-t)),
# F_n(t) = (pi n^2)^(-1/4 - i t/2) Gamma(1/4 + i t/2, pi n^2), for N = 1 and 2 (the weights w_n = 1 for n <= N; for real t this is
# 1/2 - (t^2 + 1/4) sum Re F_n(t), the function of eq. mi:eq:riemannt with psi replaced by sum_{n<=N} exp(-pi n^2 x)).
# PX.3a c_E = sum_{m>N} (4 pi m^2 - 1) exp(-pi m^2) (Theorem mi:thm:kink and the remark after it), in balls.
# PX.3b the real zeros in (0, 60]: sign changes on the grid t = j/20 (every grid value has its sign certified by a ball evaluation, 256 bits),
# each zero refined by Newton (500 bits, central difference, checked at 700 bits); the count is completed by the argument principle on the thin box
# [0.0125, 60.0125] x [-0.01, 0.01] (sampled, see below), which counts real zeros and nonreal zeros with |Im t| < 0.01 together; the number of
# zeros of Xi (the ordinates gamma_k <= 60) is read from data/zeros_1700.txt.
# PX.3c the nonreal zeros with 0 < Re t <= 60 and 0.01 <= Im t <= 12: argument principle on the boxes [10j, 10j + 10] x [0.01, 12], j = 0..5, at
# 256 bits and again at 384 bits (counts must agree), each zero isolated by a quadtree of boxes with count 1 and refined by Newton as above.
# The argument principle here is the SAMPLED one (px_lib.sampled_winding: adaptive steps, each half step turning by less than 0.35 rad and the
# half turns adding up to the whole turn to 1e-4, winding checked to be an integer to 1e-3): a ball enclosure of whole segments is not
# feasible, because |Xi_w| falls to 1e-11 on the contour (N = 2) while the ball radius of the evaluation is of the order of the segment
# length times 1 (the pieces F_n are O(1) before cancelling). Status of the counts: computed; sign changes on the grid and the refined zeros: certified
# signs, computed locations.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 px_xiw.py [--quick] [--out FILE] [zeros_file] ; --quick: N = 1, window j = 0 and 2 and the grid up to t = 25.
import argparse, json, sys, time, math
from concurrent.futures import ProcessPoolExecutor
from flint import acb, arb, ctx
sys.path.insert(0, '.')
import px_lib as L


def xiw(N):
    def f(t):
        tot = acb(0)
        for n in range(1, N + 1):
            a = arb.pi() * (n * n)
            la = a.log()
            ac = acb(a)
            for sg in (1, -1):
                s = acb(arb(1) / 4) + (acb(0, 1) * t * sg) / 2
                tot += (-s * la).exp() * ac.gamma_upper(s)
        return acb(arb(1) / 2) - acb(arb(1) / 2) * (t * t + acb(arb(1) / 4)) * tot
    return f


def cE(N):
    ctx.prec = 300
    tot = arb(0)
    for m in range(N + 1, 40):
        a = arb.pi() * (m * m)
        tot += (4 * a - 1) * (-a).exp()
    return tot


def task_real(N, tmax, prec=256):
    ctx.prec = prec
    f = xiw(N)
    j = 0
    vals = []
    for j in range(0, int(tmax * 20) + 1):
        v = f(acb(arb(j) / 20))
        r = v.real
        if r.contains(0):
            raise L.Fail('sign not certified at t=%g' % (j / 20))
        vals.append((j, r.mid() > 0))
    zeros = []
    for (j0, s0), (j1, s1) in zip(vals, vals[1:]):
        if s0 != s1:
            z, res, it = L.newton(f, complex((j0 + j1) / 40, 0.0), 500)
            z2, res2, it2 = L.newton(f, complex(float(z.real), 0.0), 700)
            ok = j0 / 20 < float(z.real) < j1 / 20
            zeros.append(dict(t=z2.real.str(30, radius=False), interval=(j0 / 20, j1 / 20), inside=ok, resid=res2, diff=float(abs(z2 - z)), im=float(abs(z2.imag))))
    ctx.prec = prec
    box = (0.0125, tmax + 0.0125, -0.01, 0.01)
    c, tot, mn, n = L.sampled_winding(f, box)
    return dict(N=N, kind='real', grid_sign_changes=len(zeros), zeros=zeros, thinbox=box, thinbox_count=c, thinbox_total=tot, thinbox_min=mn, evals=n,
                sign_at_0=bool(vals[0][1]), sign_at_end=bool(vals[-1][1]))


def task_window(N, j, ymax=12.0, ylow=0.01):
    f = xiw(N)
    box = (10.0 * j, 10.0 * j + 10.0, ylow, ymax)
    res = dict(N=N, kind='window', j=j, box=box)
    counts = {}
    for prec in (256, 384):
        ctx.prec = prec
        c, tot, mn, n = L.sampled_winding(f, box)
        counts[prec] = dict(count=c, total=tot, min=mn, evals=n)
    res['counts'] = counts
    ctx.prec = 256
    log = []
    tiles = L.isolate(lambda b: L.sampled_winding(f, b)[0], box, 2.5, log=log.append)
    zs = []
    for b, c in tiles:
        x0, x1, y0, y1 = b
        z, r, it = L.newton(f, complex((x0 + x1) / 2, (y0 + y1) / 2), 500)
        z2, r2, it2 = L.newton(f, complex(float(z.real), float(z.imag)), 700)
        zs.append(dict(tile=b, tilecount=c, re=z2.real.str(30, radius=False), im=z2.imag.str(30, radius=False), resid=r2, diff=float(abs(z2 - z)),
                       inside=(x0 < float(z.real) < x1 and y0 < float(z.imag) < y1)))
    res['zeros'] = zs
    res['log'] = log
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='px_xiw.json')
    ap.add_argument('zeros', nargs='?', default='zeros_1700.txt')
    args = ap.parse_args()
    out = {}
    Ns = [1] if args.quick else [1, 2]
    for N in (1, 2):
        c = cE(N)
        out['cE_%d' % N] = c.str(30, radius=False)
        print('PX.3a N=%d c_E =' % N, c.str(15, radius=False), flush=True)
    gam = [float(l.split()[-1]) for l in open(args.zeros) if l.strip()]
    nxi = sum(1 for g in gam if g <= 60)
    out['n_gamma_le_60'] = nxi
    out['gamma_13'] = gam[12]
    print('zeros of Xi with gamma <= 60:', nxi, 'gamma_13 =', gam[12], 'gamma_14 =', gam[13], flush=True)
    tmax = 25.0 if args.quick else 60.0
    wins = [0, 2] if args.quick else list(range(6))
    with ProcessPoolExecutor(max_workers=8) as ex:
        fr = [ex.submit(task_real, N, tmax) for N in Ns]
        fw = [ex.submit(task_window, N, j) for N in Ns for j in wins]
        res = []
        for f in fr + fw:
            r = f.result()
            res.append(r)
            if r['kind'] == 'real':
                print('REAL N=%d: sign changes %d, thin box count %d (total %.6f, min|f| %.3g)' % (r['N'], r['grid_sign_changes'], r['thinbox_count'], r['thinbox_total'], r['thinbox_min']), flush=True)
                for z in r['zeros']:
                    print('   real zero', z, flush=True)
            else:
                print('WINDOW N=%d j=%d counts %r' % (r['N'], r['j'], r['counts']), flush=True)
                for l in r['log']:
                    print('  ', l, flush=True)
                for z in r['zeros']:
                    print('   zero', z['re'][:14], z['im'][:12], 'resid %.2g diff %.2g inside %s' % (z['resid'], z['diff'], z['inside']), flush=True)
    out['results'] = res
    json.dump(out, open(args.out, 'w'), indent=1, default=str)
