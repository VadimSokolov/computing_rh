"""Reproduce D10 examples. Run: timeout 30 python3 -I verify_gap.py."""

from decimal import Decimal as D, getcontext
from fractions import Fraction as Q
import argparse


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--precision', type=int, default=80)
args = parser.parse_args()
getcontext().prec = args.precision


def sinh(x):
    return (x.exp() - (-x).exp()) / 2


def cosh(x):
    return (x.exp() + (-x).exp()) / 2


def asinh(x):
    return (x + (1 + x*x).sqrt()).ln()


def atan(x):
    # Half-angle reduction followed by the alternating power series.
    y = x / (1 + (1+x*x).sqrt())
    term, total = y, y
    for k in range(1, 300):
        term *= -y*y
        new = total + term / (2*k+1)
        if new == total:
            break
        total = new
    return 2*total


def bisect(f, lo, hi, steps=240):
    lo, hi = D(lo), D(hi)
    flo = f(lo)
    assert flo * f(hi) < 0
    for _ in range(steps):
        mid = (lo + hi) / 2
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2


def quartic_rhs(w):
    phase = atan(4 * w / (5 - w * w))
    b = 2 * w * (3 + w * w) / ((5 - w * w)**2 + 16 * w * w)
    return phase**2 / (4 * (b + 1 / (1 + w)))


def law_data(eta):
    eta = D(eta)
    d = atan(sinh(eta))
    B = sinh(eta)**2 / cosh(eta)
    t = (B * B + B * (B * B + 2).sqrt()) / 2
    exact = asinh(t.sqrt())
    h = 1/cosh(eta)
    v = sinh(eta)/cosh(eta)

    def rhs(w):
        C, S = cosh(w), sinh(w)
        a, c = 1-h*C, v*S
        real, imag = a*a-c*c+B*B, -2*a*c
        preal = -2*(a*v*C-c*h*S)
        pimag = 2*(a*h*S+c*v*C)
        assert real > 0
        phase = atan(-imag/real)
        b = -(pimag*real-preal*imag)/(real*real+imag*imag) + 2*w/(eta*eta-w*w)
        return phase**2 / (4 * (b + 1 / (eta + w)))

    fp = bisect(lambda w: eta - w - rhs(w), eta * D('0.75'), eta * D('0.95'))
    return d, B, exact, fp, rhs(exact)


if __name__ == '__main__':
    exact = 1/D(2).sqrt().sqrt()
    fp = bisect(lambda w: 1 - w - quartic_rhs(w), '0.8', '0.9')
    print('Quartic: exact w, exact gap, fixed point w, fixed point gap')
    print(*(f'{v:.30f}' for v in (exact, 1-exact, fp, 1-fp)))
    for w in ('0.84', '0.85', '0.86'):
        print('quartic residual at w', w, f'{1-D(w)-quartic_rhs(D(w)):.30f}')
    print('Law: eta, d, B, exact w/eta, actual gap/eta, fixed point gap/eta, RHS at exact w / eta')
    for eta in ('0.01', '0.1', '0.3', D(2).ln()):
        eta = D(eta)
        d, B, exact, fp, rhs_exact = law_data(eta)
        print(*(f'{v:.24f}' for v in (eta, d, B, exact/eta, 1-exact/eta, 1-fp/eta, rhs_exact/eta)))
    w = Q(17, 20)
    x = 4*w/(5-w*w)
    b = 2*w*(3+w*w)/((5-w*w)**2+16*w*w)
    # sqrt(1711^2+1360^2) > 2185 gives tan(atan(x)/2) < 170/487.
    y = Q(170, 487)
    assert 1711**2 + 1360**2 > 2185**2
    upper_atan = 2*(y-y**3/3+y**5/5)
    upper_rhs = upper_atan**2/(4*(b+1/(1+w)))
    assert upper_atan < Q(6719, 10000)
    assert Q(6719, 10000)**2 < Q(3, 5)*(b+1/(1+w))
    assert 2*w**4 > 1
    assert upper_rhs < Q(3, 20)
    print('Exact rational inequalities: 2*(17/20)^4 > 1 and H(17/20) < 3/20: PASS')
    print('Rational upper bound for H(17/20):', f'{D(upper_rhs.numerator)/D(upper_rhs.denominator):.24f}')
    print('Sharpness among simple zeros: t, eta, wstar/eta for A=sqrt(1+t^2), B=t^3')
    for text in ('0.1', '0.01', '0.001'):
        t = D(text)
        a2, b2 = 1+t*t, t**6
        c = 1-a2-b2
        u = (-c+(c*c+4*b2).sqrt())/2
        k = 1-a2-2*b2
        v = (-k+(k*k+8*b2).sqrt())/4
        eta = asinh(u.sqrt())
        ratio = asinh(v.sqrt())/eta
        print(*(f'{number:.24f}' for number in (t, eta, ratio)))
