"""Track E: enclosures of (log xi)'''(1/2+h) and (log Lambda_DH)'''(1/2+h) at sample h, by acb_series arithmetic
(python-flint 0.9.0 at /tmp/d1lib; Arb trusted). Each printed value is a ball (midpoint +- radius)."""
import sys; sys.path.insert(0, '/tmp/d1lib')
from flint import acb_series, acb, arb, ctx
ctx.prec = 256; ctx.cap = 4
half = acb(1)/2
pi = acb.pi()
kappa = ((arb(10) - 2*arb(5).sqrt()).sqrt() - 2)/(arb(5).sqrt() - 1)
coef = {1: arb(1), 2: kappa, 3: -kappa, 4: arb(-1)}
def third(F):
    return F.coeffs()[3]*6
for h in [1, 2, 8, 24, 40]:
    s = acb_series([half + h, 1])
    try:
        xi = s*(s-1)/2 * (-(s/2)*pi.log()).exp() * (s/2).gamma() * s.zeta()
        gP = third(xi.log()).real
        f = (-(s*acb(5).log())).exp() * sum((acb_series([acb(c)])*s.zeta(acb(r)/5) for r, c in coef.items()), acb_series([0]))
        Lam = ((s/2)*(acb(5)/pi).log()).exp() * ((s+1)/2).gamma() * f
        gD = third(Lam.log()).real
        print("h=%3d  Polya (log xi)''' in %s   DH (log Lambda)''' in %s" % (h, gP.str(12, radius=True), gD.str(12, radius=True)))
    except Exception as e:
        print("h=%d failed: %r" % (h, e))
