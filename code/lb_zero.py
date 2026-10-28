# lb_zero.py: LB.2 and LB.3 of the checks of the note on the lattice bias (review-LB.md M1 to M4 and section 7; Chapter 25 on extreme values, Section
# ev:sec:decomp, the Dyck chain): zeros of the Dirichlet polynomial F_n(s) = E[Y_n^s] = sum_h P(H_n = h) (sqrt(2/pi) h / sqrt(2n))^s of the Dyck chain (H_n the height
# of a uniform Dyck path of length 2n; ch/ev1.tex, the normalisation of px_dyck.py: F_n -> 2 xi(s), F_n(1) -> 1) near the zeros rho_k = 1/2 + i gamma_k of xi, against the
# first order law rho_k(n) = rho_k + delta_n c(rho_k), delta_n = 3 / (2 sqrt(2n)), c(rho) = rho sqrt(2/pi) xi(rho-1)/xi'(rho) (lb_c.py). The zeros of F_n do not depend on
# the constant sqrt(2/pi)^s, so the normalisation of the note (zeros of E (H_n/sqrt(2n))^s) gives the same zeros and the same delta_n (shift of the scale of M, not of Y).
# Exact law at large n. By the reflection principle P(H_n <= h) = G(h) = (n+1) sum_{k in Z} [r(k m) - r(k m - 1)], m = h + 2, r(j) = C(2n, n+j)/C(2n, n) = r(-j),
# r(j+1) = r(j) (n - j)/(n + j + 1) (a ball recurrence); the k = 0 term is 1/(n+1), so G(h) = 1 - (n+1) T(h), T(h) = sum_{k>=1} [r(km-1) - 2 r(km) + r(km+1)] and
# P(H_n = h) = (n+1)(T(h-1) - T(h)): no cancellation of the big terms. Terms with P < 10^-cut are dropped (law A: cut 110 at 512 bits; law B: cut 90 at 320 bits), the first h is
# chosen so that P(H_n = h) < 10^-(cut+30) below it, the last when 1 - G < 10^-cut. --selftest compares the law with the exact integers of px_dyck.dyck_counts at n = 1600
# and logs sum P and F_n(1), F_n(2) against 1 and pi/3.
# Zero finding. Newton with the exact derivative F_n'(s) = sum P L e^{sL} from the first order prediction (and, with --cont, also from the zero found at the next larger n of the
# list, which is the zero used: the branch followed in n); the two are compared. The zero is then polished with law B (different precision and cut): the difference of the
# two zeros is logged (the error control). Uniqueness: the argument principle by px_lib.sampled_winding with law B on the square of half side --box about the zero (computed:
# adaptive sampling, not certified; the ball enclosure of segments fails at once, |F_n| ~ 1e-58 at height 170 against terms of size 1). The minimum of |F_n| on the contour is logged.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 lb_zero.py [--selftest] [--quick] --ns 100000,1000000 --ks 63,64,65 [--cont] [--box 0.25] [--out FILE] [--workers 8] [zeros_file]
import argparse, json, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor
from flint import acb, arb, ctx, fmpz
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '.')
import lb_c
import px_lib as L


def dyck_law(n, prec, cut):
    ctx.prec = prec
    thr = 10.0 ** (-cut)
    jmax = int(math.sqrt(n * (cut + 25) * math.log(10))) + 10
    r = [arb(1)]
    nb = arb(n)
    for j in range(jmax + 2):
        r.append(r[-1] * (nb - j) / (nb + j + 1))
    xlo = math.pi / math.sqrt(2 * (cut + 30) * math.log(10))
    hlo = max(2, int(0.9 * xlo * math.sqrt(2 * n)))
    n1 = arb(n + 1)

    def T(h):
        m = h + 2
        s = arb(0)
        k = 1
        while k * m - 1 <= jmax:
            j = k * m
            s += r[j - 1] - 2 * r[j] + (r[j + 1] if j + 1 <= jmax + 1 else arb(0))
            k += 1
        return s
    c0 = (arb(2) / arb.pi()).sqrt() / arb(2 * n).sqrt()
    P, Lg, H = [], [], []
    Tprev = T(hlo - 1)
    h = hlo
    h_end = None
    while True:
        Tc = T(h)
        p = n1 * (Tprev - Tc)
        if float(p.mid()) > thr:
            P.append(p)
            Lg.append((c0 * h).log())
            H.append(h)
        Tprev = Tc
        if h > math.sqrt(2 * n) and float((n1 * Tc).mid()) < thr:
            h_end = h
            break
        h += 1
        if h + 2 > jmax:
            raise RuntimeError('jmax too small')
    return dict(n=n, prec=prec, cut=cut, P=P, L=Lg, H=H, nterms=len(P), hlo=hlo, hend=h_end)


def FdF(law, s, deriv=True):
    ctx.prec = law['prec']
    tot = acb(0)
    dt = acb(0)
    if deriv:
        for p, l in zip(law['P'], law['L']):
            e = p * (s * l).exp()
            tot += e
            dt += e * l
        return tot, dt
    for p, l in zip(law['P'], law['L']):
        tot += p * (s * l).exp()
    return tot


def F(law):
    return lambda s: FdF(law, s, False)


def newton(law, z0, iters=60, tol=1e-40, maxstep=0.5):
    ctx.prec = law['prec']
    z = acb(arb(z0.real), arb(z0.imag)) if not isinstance(z0, acb) else acb(z0.real.mid(), z0.imag.mid())
    steps = 0
    for it in range(iters):
        f, d = FdF(law, z)
        dz = f / d
        nz = float(abs(dz).mid())
        if nz > maxstep:
            dz = dz * (maxstep / nz)
        z = z - dz
        z = acb(z.real.mid(), z.imag.mid())
        steps = it + 1
        if nz < tol:
            break
    f, d = FdF(law, z)
    return z, float(abs(f).mid()), steps


def cplx(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))


def fmt(z, d=30):
    return z.real.str(d, radius=False) + (' + ' if float(z.imag.mid()) >= 0 else ' - ') + abs(z.imag).str(d, radius=False) + 'i'


def task(k, ns, gam, cont, box, quick):
    t00 = time.time()
    ctx.prec = 600
    gb = arb(gam[k - 1])
    rho = acb(arb(0.5), gb)
    ctx.prec = 256
    c = lb_c.c_of_rho(acb(arb(0.5) + arb(0, 1e-60), gb + arb(0, 1e-60)))[0]
    c = acb(c.real.mid(), c.imag.mid())
    out = []
    prev = None
    for n in sorted(ns, reverse=True):
        t0 = time.time()
        ctx.prec = 512
        delta = 3 / (2 * math.sqrt(2 * n))
        dl = arb(3) / (arb(2) * arb(2 * n).sqrt())
        pred = rho + dl * c
        A = dyck_law(n, 512 if not quick else 256, 110 if not quick else 80)
        tb = time.time() - t0
        zp, rp, itp = newton(A, pred)
        ctx.prec = 512
        rec = dict(k=k, n=n, delta=delta, nterms=A['nterms'], hlo=A['hlo'], hend=A['hend'], sec_law=tb,
                   c_re=c.real.str(15, radius=False), c_im=c.imag.str(15, radius=False),
                   pred=fmt(pred, 20), pred_re_m_half=(pred.real - arb(0.5)).str(15, radius=False),
                   newton_pred=fmt(zp), resid_pred=rp, steps_pred=itp)
        z = zp
        rec['start'] = 'prediction'
        if cont and prev is not None:
            zc, rc, itc = newton(A, prev)
            rec['newton_cont'] = fmt(zc)
            rec['resid_cont'] = rc
            rec['steps_cont'] = itc
            ctx.prec = 512
            rec['cont_equals_pred'] = bool(abs(zc - zp) < 1e-20)
            rec['dist_cont_pred'] = float(abs(zc - zp).mid())
            z = zc
            rec['start'] = 'continuation'
        ctx.prec = 512
        rec['re'] = z.real.str(40, radius=False)
        rec['im'] = z.imag.str(40, radius=False)
        rec['re_m_half'] = (z.real - arb(0.5)).str(25, radius=False)
        rec['pred_ratio'] = float(((z.real - arb(0.5)) / (pred.real - arb(0.5))).mid())
        rec['resid_vs_law'] = float(abs(z - pred).mid())
        rec['resid_vs_law_re'] = float((z.real - pred.real).mid())
        rec['resid_vs_law_im'] = float((z.imag - pred.imag).mid())
        rec['dist_from_rho'] = float(abs(z - rho).mid())
        rec['obs_minus_rho'] = fmt(z - rho, 20)
        # second precision and cut
        B = dyck_law(n, 320, 90 if not quick else 70)
        zb, rb, itb = newton(B, z, iters=3, tol=1e-60)
        ctx.prec = 512
        rec['polish_B_diff'] = float(abs(zb - z).mid())
        rec['nterms_B'] = B['nterms']
        # uniqueness
        ctx.prec = 320
        zf = cplx(zb)
        bx = (zf.real - box, zf.real + box, zf.imag - box, zf.imag + box)
        try:
            cnt, tot, mn, ne = L.sampled_winding(F(B), bx)
            rec.update(square_count=cnt, square_total=tot, square_min=mn, square_evals=ne)
        except L.Fail as e:
            rec.update(square_count=None, square_fail=str(e))
        rec['box'] = box
        rec['sec'] = time.time() - t0
        print('k=%d n=%d terms=%d (%.0fs law) re-1/2=%s pred=%s ratio=%.5f |z-pred|=%.3g polishdiff=%.2g steps=%d square=%s sec=%.0f' % (
            k, n, A['nterms'], tb, rec['re_m_half'][:14], rec['pred_re_m_half'][:12], rec['pred_ratio'], rec['resid_vs_law'], rec['polish_B_diff'], itp,
            rec.get('square_count'), rec['sec']), flush=True)
        out.append(rec)
        prev = z
    return k, out


def selftest():
    import px_dyck
    n = 1600
    for cut in (62, 90):
        law = dyck_law(n, 384, cut)
        cat, D = px_dyck.dyck_counts(n, tail=10 ** cut)
        ctx.prec = 384
        catb = arb(fmpz(cat))
        md = 0.0
        for h, p in zip(law['H'], law['P']):
            d = arb(fmpz(D.get(h, 0))) / catb
            md = max(md, float(abs(p - d).mid()))
        sp = sum((p for p in law['P']), arb(0))
        f1 = FdF(law, acb(1), False)
        f2 = FdF(law, acb(2), False)
        print('SELFTEST n=%d cut=%d terms=%d max|P - exact integers|=%.2e  1-sum P=%.2e  F(1)=%s  F(2)=%s (pi/3=%.8f)' % (
            n, cut, law['nterms'], md, float((1 - sp).mid()), f1.real.str(10, radius=False), f2.real.str(10, radius=False), math.pi / 3), flush=True)
    for n in (25600, 10 ** 6):
        law = dyck_law(n, 384, 80)
        f1 = FdF(law, acb(1), False)
        f2 = FdF(law, acb(2), False)
        print('SELFTEST n=%d terms=%d F(1)=%s F(2)=%s (limits 1 and pi/3 = 1.04719755)' % (n, law['nterms'], f1.real.str(10, radius=False), f2.real.str(10, radius=False)), flush=True)


def parse_ks(s):
    out = []
    for part in s.split(','):
        if '-' in part:
            a, b = part.split('-')
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--ns', default='100000')
    ap.add_argument('--ks', default='64')
    ap.add_argument('--cont', action='store_true')
    ap.add_argument('--box', type=float, default=0.25)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--out', default='lb_zero.json')
    ap.add_argument('zeros', nargs='?', default='zeros_1700.txt')
    a = ap.parse_args()
    if a.selftest:
        selftest()
        sys.exit(0)
    gam = [l.split()[-1] for l in open(a.zeros) if l.strip()]
    ns = [int(x) for x in a.ns.split(',')]
    ks = parse_ks(a.ks)
    t0 = time.time()
    res = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        fs = [ex.submit(task, k, ns, gam, a.cont, a.box, a.quick) for k in ks]
        for f in fs:
            k, out = f.result()
            res[k] = out
    # residual ratios between successive n for each k
    for k in ks:
        rows = sorted(res[k], key=lambda r: r['n'])
        for r1, r2 in zip(rows, rows[1:]):
            print('RATIO k=%d n %d -> %d: |z-pred| %.4g -> %.4g ratio %.3f (n ratio %.2f)' % (k, r1['n'], r2['n'], r1['resid_vs_law'], r2['resid_vs_law'],
                  r1['resid_vs_law'] / r2['resid_vs_law'] if r2['resid_vs_law'] else float('nan'), r2['n'] / r1['n']), flush=True)
    json.dump({str(k): v for k, v in res.items()}, open(a.out, 'w'), indent=1, default=str)
    print('total sec %.0f' % (time.time() - t0))
