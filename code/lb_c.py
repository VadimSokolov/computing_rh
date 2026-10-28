# lb_c.py: LB.1 of the checks of the note on the lattice bias (review-LB.md, M1 to M4 and section 7; Chapter 25 on extreme values, Section
# ev:sec:decomp, the Dyck chain and the first order law rho_n(N) = rho_n + delta_N c(rho_n)). It computes, in ball arithmetic (python-flint acb),
# c(rho) = rho * sqrt(2/pi) * xi(rho - 1) / xi'(rho),   xi(s) = s (s - 1) pi^(-s/2) Gamma(s/2) zeta(s) / 2   (xi(0) = xi(1) = 1/2, E Y^s = 2 xi(s)),
# for the first 200 zeros rho_n = 1/2 + i gamma_n, and certifies the sign of Re c(rho_n) for every n. The Dyck normalisation: the zero of
# s -> E[(M - delta)^s] near rho is rho + delta rho E M^(rho-1) / (E M^s)'(rho), and with M = sqrt(pi/2) Y this ratio is c(rho) (the factors (pi/2)^(s/2) cancel).
# Conventions. gamma_n: the 70 digit ordinates of data/zeros_1700.txt, used as balls of radius 1e-60 in both coordinates (the real part 1/2 is exact for
# these zeros, the ball only has to contain the zero). xi'(rho) is the coefficient of x in the power series of xi(rho + x) (acb_series; the derivative is enclosed on the ball without using zeta(rho) = 0);
# xi(rho - 1) is evaluated twice: directly at the ball rho - 1 (real part -1/2) and through xi(rho - 1) = xi(2 - rho) (real part 3/2); the two balls must overlap.
# Error control. Every number is a ball; the sign of Re c is certified when the ball of Re c excludes 0; ball radius of c is logged; the run is repeated at 128 and
# 256 bits and the midpoints are compared. simple zero: the ball of xi'(rho) excludes 0. The ball contains a zero: xi(ball) is logged (modulus below 1e-60).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 lb_c.py [--quick] [--out FILE] [--count 200] [zeros_file]
import argparse, json, sys, time
from flint import acb, acb_series, arb, ctx


def xi(s):
    return s * (s - 1) * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta() / 2


def xi_series(rho):
    """first two Taylor coefficients (xi(rho), xi'(rho)) of xi at the ball rho, by power series arithmetic in acb_series (no use of zeta(rho) = 0)"""
    S = acb_series([rho, 1], prec=2)
    lp = acb(arb.pi().log())
    v = S * (S - 1) / 2 * (-(S / 2) * lp).exp() * (S / 2).gamma() * S.zeta()
    co = v.coeffs()
    return co[0], co[1]


def c_of_rho(rho):
    """returns (c, xi(rho-1) direct, xi(rho-1) through 2-rho, xi', zeta(mid))"""
    x1 = xi(rho - 1)
    x2 = xi(2 - rho)
    x0, xp = xi_series(rho)
    c = rho * (arb(2) / arb.pi()).sqrt() * x2 / xp
    return c, x1, x2, xp, x0



if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='lb_c.json')
    ap.add_argument('--count', type=int, default=200)
    ap.add_argument('zeros', nargs='?', default='zeros_1700.txt')
    a = ap.parse_args()
    gam = [l.split()[-1] for l in open(a.zeros) if l.strip()]
    cnt = 3 if a.quick else a.count
    gam = gam[:cnt]
    t0 = time.time()
    res = {}
    # arb(str, rad) is not supported by flint 0.9, so build the balls from strings explicitly
    def run2(prec):
        ctx.prec = prec
        rows = []
        for n, g in enumerate(gam, 1):
            ctx.prec = 400
            gb = arb(g)
            ctx.prec = prec
            rho = acb(arb(0.5) + arb(0, 1e-60), gb + arb(0, 1e-60))
            c, x1, x2, xp, x0 = c_of_rho(rho)
            zmid = abs(x0)
            rows.append(dict(n=n, c=c, overlap=bool(x1.overlaps(x2)), simple=bool(not xp.contains(0)), zmid=float(zmid.mid()),
                             rad=float(abs(c).rad()) if hasattr(abs(c), 'rad') else 0.0))
        return rows
    r256 = run2(256)
    r128 = run2(128)
    pos = []
    summary = []
    maxdiff = 0.0
    for a1, a2 in zip(r256, r128):
        c = a1['c']
        cm = complex(float(c.real.mid()), float(c.imag.mid()))
        d = float(abs(a1['c'] - a2['c']).mid())
        # difference of midpoints, to be compared with the radii
        maxdiff = max(maxdiff, float(abs(acb(a1['c'].real.mid(), a1['c'].imag.mid()) - acb(a2['c'].real.mid(), a2['c'].imag.mid())).mid()))
        neg = bool(c.real < 0)
        p = bool(c.real > 0)
        if p:
            pos.append(a1['n'])
        summary.append(dict(n=a1['n'], re=c.real.str(15, radius=False), im=c.imag.str(15, radius=False), re_rad=float(c.real.rad()), neg=neg, pos=p,
                            sign_certified=bool(neg or p), sign128=('neg' if a2['c'].real < 0 else 'pos' if a2['c'].real > 0 else '?'),
                            overlap=a1['overlap'], simple=a1['simple'], zeta_mid=a1['zmid']))
    for s in summary:
        print('n=%3d Re c=%s Im c=%s rad(Re)=%.1e sign=%s (128 bit: %s) xi(rho-1) balls overlap=%s simple=%s |xi(ball)|=%.1e' % (
            s['n'], s['re'], s['im'], s['re_rad'], 'neg' if s['neg'] else 'POS' if s['pos'] else 'UNDECIDED', s['sign128'], s['overlap'], s['simple'], s['zeta_mid']), flush=True)
    print('N ZEROS', len(summary), ' signs certified:', sum(s['sign_certified'] for s in summary), ' all overlap:', all(s['overlap'] for s in summary),
          ' all simple:', all(s['simple'] for s in summary), ' max |xi(ball)|: %.2e' % max(s['zeta_mid'] for s in summary), flush=True)
    print('Re c > 0 for n =', pos, flush=True)
    print('max difference of the 256 and 128 bit midpoints of c: %.2e' % maxdiff, flush=True)
    print('first n with Re c > 0:', pos[0] if pos else None, ' Re c < 0 for all n below it:', all(s['neg'] for s in summary if pos and s['n'] < pos[0]), flush=True)
    json.dump(dict(rows=summary, positive=pos, maxdiff=maxdiff, sec=time.time() - t0), open(a.out, 'w'), indent=1)
    print('sec %.1f' % (time.time() - t0))
