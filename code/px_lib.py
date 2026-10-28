# px_lib.py: shared helpers of the computations PX.2 to PX.4 (Chapters on Miyatake, ferromagnets and extreme values): argument principle
# on rectangles (a sampled version with adaptive step control, and a certified version with ball enclosures of the segments), Newton
# with a central difference derivative, and quadtree isolation of the zeros. Written from the mathematics (the collaborator's scripts
# were not read, run or copied). The caller sets flint.ctx.prec.
import math
from flint import acb, arb, ctx


class Fail(Exception):
    pass


def pt(z):
    return acb(arb(z.real), arb(z.imag))


def dphase(fa, fb):
    return float((fb / fa).arg())


def sampled_winding(f, box, tol=0.35, maxdepth=44):
    """Winding number of f(z) along the boundary of box=(x0,x1,y0,y1), counterclockwise, by adaptive sampling: a piece [a,b] is accepted
    when both half pieces turn by less than tol and the half turns add up to the whole turn (to 1e-4). Returns (count, total/2pi,
    min |f| over the sampled points, evaluations). The result is exact only if f has no zero between samples (computed, not certified)."""
    x0, x1, y0, y1 = box
    corners = [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0)]
    st = {'n': 0, 'min': float('inf')}
    cache = {}

    def ev(z):
        k = (z.real, z.imag)
        if k not in cache:
            v = f(pt(z))
            st['n'] += 1
            a = float(abs(v))
            if a < st['min']:
                st['min'] = a
            cache[k] = v
        return cache[k]

    def piece(a, b, fa, fb, depth):
        m = (a + b) / 2
        fm = ev(m)
        d1 = dphase(fa, fm)
        d2 = dphase(fm, fb)
        d = dphase(fa, fb)
        if abs(d1) < tol and abs(d2) < tol and abs(d1 + d2 - d) < 1e-4:
            return d1 + d2
        if depth >= maxdepth:
            raise Fail('sampling depth exceeded near %s' % (m,))
        return piece(a, m, fa, fm, depth + 1) + piece(m, b, fm, fb, depth + 1)

    tot = 0.0
    for i in range(4):
        a, b = corners[i], corners[i + 1]
        # start with 8 equal pieces
        pts = [a + (b - a) * j / 8 for j in range(9)]
        for j in range(8):
            tot += piece(pts[j], pts[j + 1], ev(pts[j]), ev(pts[j + 1]), 0)
    c = tot / (2 * math.pi)
    if abs(c - round(c)) > 1e-3:
        raise Fail('winding %r is not an integer' % c)
    return int(round(c)), c, st['min'], st['n']


def cert_winding(f, box, minlen=1e-12, maxseg=3_000_000):
    """Certified winding number of f along the boundary of box: every segment is enclosed by a ball evaluation of f (the argument
    is a complex ball covering the segment); when the image ball excludes 0 the turn along the segment is the principal argument of
    f(b)/f(a). Raises Fail if a segment of length minlen still meets 0 or maxseg segments are used. Returns (count, total/2pi,
    certified lower bound of |f| on the contour, segments)."""
    x0, x1, y0, y1 = box
    corners = [complex(x0, y0), complex(x1, y0), complex(x1, y1), complex(x0, y1), complex(x0, y0)]
    st = {'n': 0, 'lb': float('inf')}
    cache = {}

    def pv(z):
        k = (z.real, z.imag)
        if k not in cache:
            cache[k] = f(pt(z))
        return cache[k]

    def seg(a, b):
        st['n'] += 1
        if st['n'] > maxseg:
            raise Fail('too many segments')
        m = (a + b) / 2
        rx = abs(b.real - a.real) / 2
        ry = abs(b.imag - a.imag) / 2
        v = f(acb(arb(m.real, rx), arb(m.imag, ry)))
        mid = abs(v.mid())
        r = float(v.real.rad()) + float(v.imag.rad())
        mabs = float(mid)
        if mabs > r / 0.7 and r == r:
            lb = mabs - r
            if lb < st['lb']:
                st['lb'] = lb
            return dphase(pv(a), pv(b))
        if abs(b - a) < minlen:
            raise Fail('zero not excluded on a segment of length %g near %s' % (abs(b - a), m))
        return seg(a, m) + seg(m, b)

    tot = 0.0
    for i in range(4):
        tot += seg(corners[i], corners[i + 1])
    c = tot / (2 * math.pi)
    if abs(c - round(c)) > 1e-6:
        raise Fail('certified winding %r is not an integer' % c)
    return int(round(c)), c, st['lb'], st['n']


def newton(f, z0, prec, h=1e-30, iters=60, tol=1e-45, maxstep=0.5):
    """Newton for f(z)=0 with a central difference derivative (h=1e-30, prec 500 gives an error below 1e-58). Returns (z, |f(z)|, steps)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        z = acb(arb(z0.real), arb(z0.imag))
        hh = acb(arb(h))
        for it in range(iters):
            fz = f(z)
            d = (f(z + hh) - f(z - hh)) / (2 * hh)
            dz = fz / d
            nz = float(abs(dz))
            if nz > maxstep:
                dz = dz * (maxstep / nz)
            z = z - dz
            z = acb(z.real.mid(), z.imag.mid())
            if nz < tol:
                break
        return z, float(abs(f(z))), it + 1
    finally:
        ctx.prec = old


def isolate(count, box, maxsize, depth=0, log=None):
    """Quadtree: returns a list of (box, 1) with certified/computed count 1 and size <= maxsize; boxes with count >1 at depth 40 are
    returned with their count. count(box) -> int may raise Fail when a zero is too near the boundary: other split fractions are tried."""
    c = count(box)
    if c == 0:
        return []
    x0, x1, y0, y1 = box
    if c == 1 and max(x1 - x0, y1 - y0) <= maxsize:
        return [(box, 1)]
    if depth >= 40:
        return [(box, c)]
    for fx, fy in ((0.4137, 0.4137), (0.3011, 0.5821), (0.5821, 0.3011), (0.6563, 0.2671)):
        xm = x0 + (x1 - x0) * fx
        ym = y0 + (y1 - y0) * fy
        kids = [(x0, xm, y0, ym), (xm, x1, y0, ym), (x0, xm, ym, y1), (xm, x1, ym, y1)]
        try:
            counts = [count(k) for k in kids]
        except Fail:
            continue
        if sum(counts) != c and log:
            log('WARNING: children counts %r do not add up to %d for box %r' % (counts, c, box))
        out = []
        for k, ck in zip(kids, counts):
            if ck:
                out += isolate(count, k, maxsize, depth + 1, log)
        return out
    raise Fail('no split fraction works for box %r' % (box,))
