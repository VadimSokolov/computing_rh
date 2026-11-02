"""Original Track CS ball computations. Run: timeout 600 python3 -I cs_certify.py."""
import sys
import time
import math
from pathlib import Path
sys.path.insert(0, '/tmp/cs_lib')
import flint
from flint import arb, acb, arb_poly, ctx

ctx.prec = 160
out = []
def record(*items):
    line = ' '.join(str(x) for x in items)
    out.append(line)
    print(line, flush=True)

def ghs_coefficients():
    # A deliberately coarse enclosure of pi makes the certificate easy to audit.
    # The standard rational bounds 3.14159 < pi < 3.14160 suffice.
    pi = arb('3.141595', '0.000005')
    assert pi.contains(arb.pi())
    r = 3 / (2*pi)
    t = arb_poly([0, 1])
    f = 1-r*t**2-r*t**7+t**9
    def d(p):
        return t*p.derivative()
    h = d(d(d(f)))*f**2-3*d(d(f))*d(f)*f+2*d(f)**3
    p = 8*pi*(t**4-1)*f**3-t**2*h
    q = p(arb_poly([1, 1]))
    record('GHS degree', q.degree(), 'constant is exactly zero by evenness')
    assert q.degree() == 31
    assert q[0].contains(0)
    for j in range(1, 32):
        assert q[j] > 0, (j, q[j])
        record('GHS coefficient', j, q[j].str(15))
    bounds = [math.floor(float(q[j].lower()))-1 for j in range(1, 32)]
    assert all(q[j] > bounds[j-1] > 0 for j in range(1, 32))
    record('GHS exact integer lower bounds', bounds)
    record('GHS all 31 nonconstant coefficients strictly positive: PASS')
    pi = arb.pi()
    ratio = 3/(2*pi)
    u2 = 32*pi-(81-25*ratio)/(1-ratio)
    record('refined U_second_at_zero', u2.str(25))

def moments(name, a, b, refined=False):
    a, b = arb(a), arb(b)
    cutoff = arb(2 if refined else 8)
    p = arb.pi()
    def weight(z):
        if refined:
            return (2*p*p*(9*z).cosh()-3*p*(5*z).cosh()) * (-2*p*(4*z).cosh()).exp()
        # Removing exp(-a) improves conditioning without changing normalized moments.
        return (-a*(z.cosh()-1)-b*z*z).exp()

    raw = []
    for k in (0, 2, 4, 6):
        count = [0]
        def fun(z, analytic):
            # Entire integrand, so analytic=True adds no restriction.
            count[0] += 1
            return z**k * weight(z)
        integral = acb.integral(fun, 0, cutoff, rel_tol=arb(2)**-125,
                                abs_tol=arb(2)**-135, eval_limit=200000,
                                depth_limit=40)
        assert integral.is_finite() and integral.imag.contains(0)
        # For x>=L, x^k <= exp(k*x). The exponents below are concave and
        # bounded above by their tangent at L, of negative slope.
        if refined:
            # 0 < w(x) <= 2*pi^2 exp(9x-pi exp(4x)).
            slope = 4*p*(4*cutoff).exp()-(k+9)
            tail = 2*p*p*((k+9)*cutoff-p*(4*cutoff).exp()).exp()/slope
        else:
            # 0 < w(x) <= exp(a-a exp(x)/2); b>=0.
            assert b >= 0
            slope = a*cutoff.exp()/2-k
            tail = (a+k*cutoff-a*cutoff.exp()/2).exp()/slope
        assert slope > 0 and tail > 0
        enclosure = 2*integral.real + arb(tail, tail)
        raw.append(enclosure)
        record(name, 'raw_moment', k, enclosure.str(30), 'one_sided_tail<=', tail.str(8), 'evals', count[0])
    m2, m4, m6 = [x/raw[0] for x in raw[1:]]
    k4 = m4-3*m2*m2
    k6 = m6-15*m4*m2+30*m2**3
    hankel = 3*m2*k6-10*k4*k4
    record(name, 'variance', m2.str(30))
    record(name, 'kappa4', k4.str(30))
    record(name, 'kappa6', k6.str(30))
    record(name, '3_variance_kappa6_minus_10_kappa4_squared', hankel.str(30))
    if name == 'cosh_a1_b10':
        assert k6 < 0
        record(name, 'GS exclusion by sixth cumulant: PASS')
    else:
        assert k4 < 0 and k6 > 0 and hankel > 0
        record(name, 'necessary moment inequalities: PASS; membership unresolved')

def main():
    started = time.monotonic()
    record('Python', sys.version.split()[0], 'python-flint', flint.__version__, 'bits', ctx.prec)
    ghs_coefficients()
    moments('cosh_a1_b0', 1, 0)
    moments('cosh_a2pi_b0', 2*arb.pi(), 0)
    moments('refined', 1, 0, True)
    moments('cosh_a1_b10', 1, 10)
    record('elapsed_seconds', round(time.monotonic()-started, 3))
    directory = Path(__file__).resolve().parent/'results'
    directory.mkdir(exist_ok=True)
    (directory/'numbers.txt').write_text('\n'.join(out)+'\n')

if __name__ == '__main__':
    main()
