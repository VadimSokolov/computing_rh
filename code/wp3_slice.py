# wp3_slice.py: work package 3 (run book Section sec:spin): the slice threshold w* of a planted pair at heights up to the 10^4-th ordinate, by the exact criterion
# of nc_cauchy_exact_criterion.py (Theorem cm:thm:exact: the tilt F(z+iw)+F(z-iw) has real zeros near m iff Phi_w(u) = atan(2wu/(u^2+eta^2-w^2)) + arg G(m+u+iw) takes a value in pi/2+pi Z),
# with flint (160 bit balls, ordinates from the certified acb.zeta_zeros, the argument taken in double precision) in place of mpmath. Compared with the gap law (1/2) phi^2 eta^3 and the
# fixed point of Proposition cm:prop:gap (as in nc_cauchy_gap_theory.py). Two thresholds are reported: "grid" (the 57 point scan of the original script, bisection 45 steps) and "refined"
# (the same with a golden section search of the maximum of the scan near its best point). Floating point: computed, not certified. Also: the pair planted at the zero of the field phi.
# Usage: python3 wp3_slice.py [--quick] [--workers 16] [--out wp3_slice.json]   (needs wp3_field.py beside it)
import argparse, json, math, sys, time
from multiprocessing import Pool
from flint import acb, arb, ctx
from wp3_field import dlogxi
ctx.prec = 160
PI = math.pi
LOGPI = None
CNT = [0]
def Xi(z):
    CNT[0] += 1
    s = acb(arb(0.5)) + acb(0, 1) * z
    return s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
class Setup:
    def __init__(self, k, eta, m=None):
        gs = acb.zeta_zeros(k, 2)
        self.gk, self.gk1 = gs[0].imag, gs[1].imag
        self.eta = eta
        self.m = float(((self.gk + self.gk1) / 2).mid()) if m is None else m
        mm = arb(self.m)
        self.c = acb(mm, arb(eta)); self.cb = acb(mm, -arb(eta))
    def G(self, x, w):
        z = acb(arb(x), arb(w))
        return Xi(z) * (z + self.c) * (z + self.cb) / ((z * z - self.gk ** 2) * (z * z - self.gk1 ** 2))
    def argG(self, x, w): return float(self.G(x, w).arg().mid())
    def logGd(self, z):
        gk, gk1 = self.gk, self.gk1
        return acb(0, 1) * dlogxi(acb(arb(0.5)) + acb(0, 1) * z) - (2 * z / (z * z - gk * gk) + 2 * z / (z * z - gk1 * gk1)) + 1 / (z + self.c) + 1 / (z + self.cb)
    def phi(self, x=None):
        z = acb(arb(self.m if x is None else x))
        return float(self.logGd(z).real.mid())
    def level_gap(self, w, refine):
        eta, m = self.eta, self.m
        psi0 = self.argG(m, w); psi0 -= PI * round(psi0 / PI); sgn = 1 if psi0 > 0 else -1
        u0 = math.sqrt(max(eta * eta - w * w, 1e-30))
        def val(t):
            u = sgn * u0 * t
            ps = self.argG(m + u, w); ps -= PI * round((ps - psi0) / PI)
            return sgn * (math.atan(2 * w * u / (u * u + eta * eta - w * w)) + ps) - PI / 2
        ts = [0.2 + 2.8 * i / 56 for i in range(57)]
        vs = [val(t) for t in ts]
        best = max(vs)
        if refine and best < 0:
            i = vs.index(best); a, b = ts[max(0, i - 1)], ts[min(56, i + 1)]
            g = (math.sqrt(5) - 1) / 2
            c, d = b - g * (b - a), a + g * (b - a); fc, fd = val(c), val(d)
            for _ in range(28):
                if fc > fd: b, d, fd = d, c, fc; c = b - g * (b - a); fc = val(c)
                else: a, c, fc = c, d, fd; d = a + g * (b - a); fd = val(d)
                best = max(best, fc, fd)
        return best
    def threshold(self, refine, steps=45):
        eta = self.eta
        lo, hi = eta * 0.9, eta * (1 - 1e-12)
        for _ in range(steps):
            mid = (lo + hi) / 2
            if self.level_gap(mid, refine) >= 0: hi = mid
            else: lo = mid
        return (lo + hi) / 2
    def fixed_point(self):
        eta, m = self.eta, self.m
        eps = 0.0
        for _ in range(30):
            w = eta - eps; z = acb(arb(m), arb(w))
            B = float(self.G(m, w).arg().mid()); d = abs(B - PI * round(B / PI))
            b = -float(self.logGd(z).imag.mid()); c = b + 1 / (2 * eta - eps)
            eps = d * d / (4 * c)
        return eps

def job(a):
    k, eta, place = a
    t0 = time.time(); CNT[0] = 0
    S = Setup(k, eta)
    row = dict(k=k, eta=eta, place=place, gamma=float(S.gk.mid()), m=S.m)
    if place == 'fieldzero':
        # pair position x with phi(x) = 0, near the midpoint, by bisection on the field of the planted configuration (pair moves with x)
        def field(x):
            T = Setup.__new__(Setup); T.__dict__.update(S.__dict__)
            mm = arb(x); T.m = x; T.c = acb(mm, arb(eta)); T.cb = acb(mm, -arb(eta))
            return T.phi(x)
        zz4 = acb.zeta_zeros(max(1, k - 1), 4)
        lo = float(zz4[0].imag.mid()) + 0.05 if k >= 2 else 1.0; hi = float(zz4[3].imag.mid()) - 0.05
        xs = [lo + (hi - lo) * i / 30 for i in range(31)]; vs = [field(x) for x in xs]
        roots = []
        for i in range(30):
            if vs[i] * vs[i + 1] < 0:
                a, b, fa = xs[i], xs[i + 1], vs[i]
                for _ in range(50):
                    mid = (a + b) / 2; fm = field(mid)
                    if fm * fa > 0: a, fa = mid, fm
                    else: b = mid
                roots.append((a + b) / 2)
        if not roots: row['error'] = 'no root of phi between the ordinates'; return row
        x0 = min(roots, key=lambda r: abs(r - S.m)); row['x_root'] = x0
        S = Setup(k, eta, m=x0)
        row['m'] = x0
    row['phi'] = S.phi(); row['lead'] = 0.5 * row['phi'] ** 2 * eta ** 3
    if place == 'midpoint': row['fixed'] = S.fixed_point()
    w1 = S.threshold(False); row['w_grid'] = w1; row['gap_grid'] = eta - w1
    w2 = S.threshold(True); row['w_refined'] = w2; row['gap_refined'] = eta - w2
    row['evals'] = CNT[0]; row['seconds'] = time.time() - t0
    print(json.dumps(row), flush=True)
    return row

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--workers', type=int, default=16); ap.add_argument('--out', default='wp3_slice.json'); ap.add_argument('--fz-only', dest='fz_only', action='store_true')
    a = ap.parse_args()
    ks = [1, 3] if a.quick else [1, 3, 8, 15, 30, 100, 300, 1000, 3000, 6000, 9999]
    tasks = [(k, e, 'midpoint') for k in ks for e in (0.1, 0.3)]
    if not a.quick: tasks += [(k, e, 'fieldzero') for k in (8, 100, 1000, 3000, 9999) for e in (0.1, 0.3)]
    if a.fz_only: tasks = [(k, e, 'fieldzero') for k in (8, 100, 1000, 3000, 9999) for e in (0.1, 0.3)]
    with Pool(a.workers) as p: rows = p.map(job, tasks, chunksize=1)
    json.dump(rows, open(a.out, 'w'), indent=0)
