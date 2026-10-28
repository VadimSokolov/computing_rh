# Note W (Chapter 13, Section fm:sec:bessel, Table fm:tab:bessel and Proposition fm:prop:tvcert), computation N2:
# certified enclosures of TV(q, p) = int_0^infty |q - p| for the members of the family of Theorem fm:thm:bessel
# fitted by wf_fit.py, at the exact seven digit decimals that wf_fit.py wrote to wf_fit.json, in the book's units
# (t, a, l, lam, mu): q = m(t) exp(-a cosh(l t)) / C, m(t) = exp(mu t^2) prod_k cosh(lam_k t),
# C = int_R m exp(-a cosh l t) dt = (2/l) E K_{Y/l}(a) (Theorem fm:thm:bessel(ii)), p = Phi/Xi(0).
# Method: that of fm_tvcert.py (which is imported: its enclosure of p, of its tail and of d''' for the part of p, the
# interval cover, and its exact integral of |alpha + beta h|). On [0, 3], cut into intervals [c - r, c + r], r = 2^-12,
# with d = q - p: | int |d| - int_{-r}^{r} |d(c) + d'(c) h| dh | <= |d''(c)| r^3/3 + M3 r^4/12, M3 >= sup |d'''| from
# interval evaluation of the Taylor series of q on a cover of [0, 3] (768 balls) plus that of p. Tail beyond 3: the
# integrand of q is at most g(t) = exp(mu t^2 + Lam t - (a/2) e^{l t})/C with Lam = sum lam_k, and g decreases at a
# rate at least h(3) = (a l/2) e^{3l} - 6 mu - Lam (h is increasing for t >= 3); so int_3^inf q <= g(3)/h(3). C is
# enclosed by the Bessel formula (Arb K_nu(a) at real order, mu = 0, an average over the 2^N values of Y) and
# independently by Arb's quadrature of 2 int_0^5 plus the same kind of tail bound; both enclosures must overlap
# (the Bessel one is used when mu = 0, the quadrature when mu > 0). Working precision 128 bits (ball arithmetic: every
# printed bound is an upper bound of the distance of the law with exactly these parameters).
# Check of the machinery: the same code on P\'olya's Phi* (a = 2 pi, l = 2, m = cosh(9t/2)) must reproduce the
# enclosure of fm_tvcert.json (name polya_star).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 wf_tvcert.py [--quick]  reads wf_fit.json and fm_tvcert.json; writes wf_tvcert.json
#   --quick uses r = 2^-6 and writes /tmp/wf_tvcert_quick.json (a seconds long run of the same code)
import json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fm_tvcert as T
from flint import arb, acb, arb_series, ctx

QUICK = '--quick' in sys.argv
PREC = 128


def parse(x):
    return arb.pi()*2 if x == '2pi' else arb(x)


class Member:
    """a, l, mu arb; lams list of arb; the law q = m exp(-a cosh l t)/C."""
    def __init__(self, name, a, l, lams, mu):
        self.name, self.a, self.l, self.lams, self.mu = name, a, l, lams, mu
        self.Lam = sum(lams, arb(0))

    @classmethod
    def from_row(cls, name, e):
        return cls(name, parse(e['a']), parse(e['l']), [arb(s) for s in e['lams']], arb(e['mu']))

    def ys(self):
        """The 2^N values of sum eps_k lam_k (mu = 0)."""
        ys = [arb(0)]
        for lam in self.lams:
            ys = [y + lam for y in ys] + [y - lam for y in ys]
        return ys

    def mser(self, X):
        """m(t) and exp(-a cosh l t) as one series in X (arb_series)."""
        s = (X*X*self.mu).exp()
        for lam in self.lams:
            s = s*T.cosh_ser(X, lam)
        return s*(T.cosh_ser(X, self.l)*(-self.a)).exp()

    def coef(self, x0):
        return T.coeffs4(self.mser(T.ser(x0)))[:4]

    def f(self, z):
        """The unnormalised density at an acb z."""
        s = (z*z*self.mu).exp()
        for lam in self.lams:
            s = s*((z*lam).cosh())
        return s*((-(z*self.l).cosh()*self.a)).exp()

    def tail(self, T0, k=0):
        """Bound of int_{T0}^inf t^k m(t) exp(-a cosh l t) dt (unnormalised), T0 >= 1."""
        ex = self.Lam + k
        g = (self.mu*T0*T0 + ex*T0 - self.a/2*(self.l*T0).exp()).exp()
        h = self.a*self.l/2*(self.l*T0).exp() - 2*self.mu*T0 - ex
        assert h > 0 and self.a*self.l*self.l/2*(self.l*T0).exp() >= 2*self.mu
        return g/h

    def C_bessel(self):
        assert self.mu == 0
        tot = arb(0); ys = self.ys()
        for y in ys:
            if y.contains(0):             # equal fields cancel exactly (lam - lam = 0): order exactly 0, as K_nu is
                assert y.rad() < 1e-30    # singular as a ball containing 0 in Arb's formula
                y = arb(0)
            tot += acb(self.a).bessel_k(acb(abs(y)/self.l)).real
        return 2/self.l*tot/len(ys)

    def C_quad(self, T0=5):
        I = acb.integral(lambda z, _: self.f(z), 0, T0, rel_tol=2**-(PREC - 20))
        return 2*(I.real + arb(0, self.tail(T0).upper()))


def overlap(x, y):
    return x.overlaps(y)


def enclosure(x):
    return [T.lo(x), T.hi(x)]


def setup():
    ctx.prec = PREC
    T.TB = (4*(2*arb.pi()**2*9**4 + 3*arb.pi()*81)*(2*arb.pi()*81 + arb('6.5'))**3*(-81*arb.pi()).exp()).upper()
    T.XI0 = T.xi_half()
    T.Pc = []
    for i in range(T.NI):
        c = T.Phi_ser(T.centre(i))
        T.Pc.append([c[0]/T.XI0, c[1]/T.XI0, c[2]/T.XI0])
    T.M3P = T.M3_cover(lambda x0: [ci/T.XI0 for ci in T.Phi_ser(x0)])
    T.TP = T.tail_Phi()/T.XI0


def certify(mem):
    t0 = time.time()
    ctx.prec = PREC
    Cq = mem.C_quad()
    if mem.mu == 0:
        Cb = mem.C_bessel()
        assert overlap(Cb, Cq), (Cb, Cq)
        C = Cb
        cinfo = dict(C_bessel=[T.lo(Cb), T.hi(Cb)], C_quad=[T.lo(Cq), T.hi(Cq)], C_used='bessel')
    else:
        C = Cq
        cinfo = dict(C_quad=[T.lo(Cq), T.hi(Cq)], C_used='quad')
    M3q = T.M3_cover(mem.coef)/C
    M3 = arb(M3q.upper()) + T.M3P
    cs = [mem.coef(T.centre(i)) for i in range(T.NI)]
    def cd(i):
        a = cs[i]; pp = T.Pc[i]
        return [a[0]/C - pp[0], a[1]/C - pp[1], a[2]/C - pp[2]]
    s = T.tv_sum(cd, M3)
    tl = mem.tail(T.A)/C + T.TP
    tv = s + arb(0, tl.upper())
    out = dict(name=mem.name, TV=[T.lo(tv), T.hi(tv)], M3=T.hi(M3), tails=T.hi(tl), secs=round(time.time() - t0, 1), **cinfo)
    out['TV_ub_2sig'] = ceil_sig(out['TV'][1], 2)
    out['TV_ub_4sig'] = ceil_sig(out['TV'][1], 4)
    print(json.dumps(out), flush=True)
    return out


def ceil_sig(x, n):
    e = math.floor(math.log10(x)) - n + 1
    return float('%.*e' % (n - 1, math.ceil(x/10**e)*10.0**e))


if __name__ == '__main__':
    t0 = time.time()
    ctx.prec = PREC                   # before any decimal is parsed: the parameters are balls at this precision
    if QUICK:
        T.R = arb(1)/64; T.NI = 96
    fit = json.load(open('wf_fit.json'))
    members = []
    for k in ('row1', 'row2', 'row3', 'row4', 'row5'):
        members.append(Member.from_row(k, fit[k]['exact']))
    r2 = fit['row2']['tv_at_decimals']
    for k in ('two_fields', 'tilt_free_kernel'):
        if k in fit and fit[k]['tv_at_decimals'] < 0.99*r2:
            members.append(Member.from_row(k, fit[k]['exact']))
            print('exploratory member %s beats row 2 (%.4g < %.4g): certified too' % (k, fit[k]['tv_at_decimals'], r2))
    members.append(Member('check_polya_star', arb.pi()*2, arb(2), [arb('4.5')], arb(0)))
    setup()
    print('Xi(0) = %s; sup|p\'\'\'| <= %.4g; tail of p beyond %d <= %.3g; r = %s; %d intervals; setup %.0f s'
          % (T.XI0.str(15), T.hi(T.M3P), T.A, T.hi(T.TP), T.R.str(5), T.NI, time.time() - t0), flush=True)
    res = [certify(m) for m in members]
    chk = [r for r in res if r['name'] == 'check_polya_star'][0]
    ref = [r for r in json.load(open('fm_tvcert.json'))['results'] if r['name'] == 'polya_star'][0]
    print('check: wf_tvcert on Phi* gives TV in %s; fm_tvcert.json gives %s' % (chk['TV'], ref['TV']))
    out = dict(prec=PREC, r=T.R.str(8), NI=T.NI, Xi0=T.XI0.str(20), results=res, check_ref=ref['TV'],
               secs=round(time.time() - t0, 1))
    json.dump(out, open('/tmp/wf_tvcert_quick.json' if QUICK else 'wf_tvcert.json', 'w'), indent=1)
    print('done in %.0f s' % (time.time() - t0))
