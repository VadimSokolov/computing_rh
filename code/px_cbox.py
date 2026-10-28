# px_cbox.py: PX.2 of the reconciliation of the notes on the pole cancellation and the excursion (Proposition fm:prop:nocomponent, Chapter
# on ferromagnets): the zeros of the characteristic function of the first component of Riemann's theta series. Up to a constant,
# C(theta) = F(theta) + F(-theta), F(theta) = pi^(-1/4 - i theta/2) Gamma(1/4 + i theta/2, pi), the Fourier transform of the density
# proportional to exp(|u|/2) exp(-pi exp(2|u|)) (u = log x / 2, x = e^(2u)). C is even and real on the real and imaginary axes, positive there
# (a characteristic function of a symmetric law is positive on the imaginary axis too: E cosh(yX) > 0), so the zeros in the closed first quadrant
# are interior to it. This script counts them in the box [0,40] x [0,30] and in [0,25] x [0,20] by the CERTIFIED argument principle (each
# boundary segment is enclosed by a ball evaluation of C, python-flint acb gamma_upper at 256 bits; the turn along a segment whose image ball
# excludes 0 is the principal argument of C(b)/C(a)), isolates the zeros by a quadtree of boxes with certified count 1, refines each by Newton at
# 500 bits (central difference derivative) and confirms by a second Newton at 700 bits. Written from the mathematics (the collaborator's scripts
# were not read, run or copied).
# Usage: python3 px_cbox.py [--quick] [--out FILE] ; --quick uses the box [0,25] x [0,20] only. Run: see checks-PX.md.
import argparse, json, sys, time
from concurrent.futures import ProcessPoolExecutor
from flint import acb, arb, ctx
sys.path.insert(0, '.')
import px_lib as L


def C(z):
    tot = acb(0)
    a = acb(arb.pi())
    la = arb.pi().log()
    for sg in (1, -1):
        s = acb(arb(1) / 4) + (acb(0, 1) * z * sg) / 2
        tot += (-s * la).exp() * a.gamma_upper(s)
    return tot


def count_box(box, prec=256):
    ctx.prec = prec
    t = time.time()
    c, tot, lb, n = L.cert_winding(C, box)
    return dict(box=box, count=c, total=tot, lb=lb, segs=n, sec=time.time() - t)


def task_locate(box):
    ctx.prec = 256
    log = []
    tiles = L.isolate(lambda b: L.cert_winding(C, b)[0], box, 3.0, log=log.append)
    res = []
    for b, c in tiles:
        x0, x1, y0, y1 = b
        z0 = complex((x0 + x1) / 2, (y0 + y1) / 2)
        z, r, it = L.newton(C, z0, 500)
        z2, r2, it2 = L.newton(C, complex(float(z.real), float(z.imag)), 700)
        inside = x0 < float(z.real) < x1 and y0 < float(z.imag) < y1
        res.append(dict(tile=b, tilecount=c, re=z2.real.str(32, radius=False), im=z2.imag.str(32, radius=False), resid500=r, resid700=r2,
                        diff=float(abs(z2 - z)), inside=inside))
    return res, log


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='px_cbox.json')
    args = ap.parse_args()
    out = {}
    # C is even, real positive on both axes and singular for a ball evaluation of Gamma(s, pi) when s is a ball around 0, -1, ... (the
    # points theta = i(2k + 1/2) of the imaginary axis), so no contour runs along the imaginary axis: the symmetric boxes [-X, X] x [0, Y]
    # contain twice the zeros of the quadrant (no zero lies on the axes), and the zeros are located in [2, X] x [0, Y].
    boxes = [(-25.0, 25.0, 0.0, 20.0), (2.0, 25.0, 0.0, 20.0), (-2.0, 2.0, 0.0, 20.0)]
    if not args.quick:
        boxes += [(-40.0, 40.0, 0.0, 30.0), (2.0, 40.0, 0.0, 30.0), (-2.0, 2.0, 0.0, 30.0)]
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(count_box, b) for b in boxes]
        futs2 = [ex.submit(count_box, b, 384) for b in boxes[:2] + boxes[3:5]]
        for f in futs + futs2:
            r = f.result()
            print('count', r, flush=True)
            out.setdefault('counts', []).append(r)
        big = (2.0, 25.0, 0.0, 20.0) if args.quick else (2.0, 40.0, 0.0, 30.0)
        x0, x1, y0, y1 = big
        sub = [(x0, 20.0, y0, 15.0), (20.0, x1, y0, 15.0), (x0, 20.0, 15.0, y1), (20.0, x1, 15.0, y1)] if not args.quick else [big]
        futs3 = [ex.submit(task_locate, b) for b in sub]
        zeros = []
        for b, f in zip(sub, futs3):
            res, log = f.result()
            for l in log:
                print(l, flush=True)
            for z in res:
                print('zero', z, flush=True)
            zeros += res
    zeros.sort(key=lambda z: (z['re'], z['im']))
    out['zeros'] = zeros
    out['n_located'] = len(zeros)
    json.dump(out, open(args.out, 'w'), indent=1, default=str)
    print('located', len(zeros), flush=True)
