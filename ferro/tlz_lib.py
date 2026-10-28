# tlz_lib.py: helper for the TLZ computations (plan reconcile-TLZ.md): a rigorous argument principle count in ball arithmetic.
# winding(F, corners): number of zeros minus poles of an analytic F inside the rectangle with the given exact corners (a closed counterclockwise contour).
# Method: each side is split recursively; on a piece the ball F(segment box) must exclude 0 and have radius < 0.7 |midpoint|, so that every value on the piece
# lies in a cone of half angle < 45 degrees about the midpoint; then the change of the argument over the piece is the principal argument of F(end)/F(start)
# (|change| < 90 degrees), computed in ball arithmetic. The sum over the contour is a ball that must contain one integer (checked).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
from flint import acb, arb, ctx

def _seg_box(a, b):
    ar, ai, br, bi = a.real, a.imag, b.real, b.imag
    re = arb((ar + br) / 2, abs(br - ar) / 2)
    im = arb((ai + bi) / 2, abs(bi - ai) / 2)
    return acb(re, im)

def _arg_ratio(Fa, Fb):
    return (Fb / Fa).log().imag

class Counter:
    def __init__(self, F, maxdepth=40):
        self.F, self.maxdepth = F, maxdepth
        self.evals = 0
        self.minabs = None      # lower bound of min |F| over the contour pieces (ball lower bounds)
        self.maxdepth_used = 0

    def seg(self, a, b, Fa, Fb, depth=0):
        self.maxdepth_used = max(self.maxdepth_used, depth)
        B = self.F(_seg_box(a, b)); self.evals += 1
        ok = False
        if B.is_finite() and B.abs_lower() > 0:
            m = B.mid()
            ok = bool(B.rad() < arb('0.7') * abs(m)) if hasattr(B, 'rad') else False
        if ok:
            lo = B.abs_lower()
            if self.minabs is None or lo < self.minabs: self.minabs = lo
            return _arg_ratio(Fa, Fb)
        if depth >= self.maxdepth:
            raise RuntimeError('argument principle: depth exceeded (zero on the contour?) near %s' % a)
        c = (a + b) / 2
        Fc = self.F(c); self.evals += 1
        return self.seg(a, c, Fa, Fc, depth + 1) + self.seg(c, b, Fc, Fb, depth + 1)

    def winding(self, x0, x1, y0, y1):
        """x0,x1,y0,y1 exact arb numbers; counterclockwise"""
        P = [acb(x0, y0), acb(x1, y0), acb(x1, y1), acb(x0, y1), acb(x0, y0)]
        FP = [self.F(p) for p in P]
        tot = arb(0)
        for i in range(4):
            tot += self.seg(P[i], P[i + 1], FP[i], FP[i + 1])
        w = tot / (2 * arb.pi())
        n = int(w.mid().floor() if False else round(float(w.mid())))
        if not (w.overlaps(arb(n)) and w.rad() < arb('0.01')):
            raise RuntimeError('winding not an integer: %s' % w)
        return n, w

def grid(x, bits=60):
    """exact dyadic arb near the float/decimal x"""
    from fractions import Fraction
    q = round(Fraction(str(x)) * 2 ** bits)
    return arb(q) / arb(2) ** bits
