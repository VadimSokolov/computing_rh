"""Independent calculations for report.md. No manuscript code is imported.

Run with a Python containing mpmath and scipy, and python-flint in ./pylib.
Ball outputs enclose exact evaluations. Quadrature outputs are computed only.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "pylib"))
import flint
from flint import acb, arb, acb_series, ctx
import mpmath as mp
import scipy
from scipy.integrate import quad

ctx.prec = 200
mp.mp.dps = 60
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
results = {}


def save(key, value):
    results[key] = value
    (OUT / "numbers.json").write_text(json.dumps(results, indent=2) + "\n")


def bs(x):
    return x.str(45)


def ms(x):
    return mp.nstr(x, 45)


def xi_mp(s):
    if s == 0 or s == 1:
        return mp.mpf(1) / 2
    return s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s) / 2


def logder(s):
    """A ball enclosure of xi'/xi, with removable values handled explicitly."""
    if s == acb(1):
        return acb(1 + arb.const_euler() / 2 - (4 * arb.pi()).log() / 2)
    if s == acb(0):
        return -logder(acb(1))
    z = acb_series([s, 1], prec=2).zeta()
    return 1 / s + 1 / (s - 1) - arb.pi().log() / 2 + (s / 2).digamma() / 2 + z[1] / z[0]


def density(alpha, w):
    return logder(acb(alpha, w)).real / arb.pi()


def dh_series(s, order=3):
    v = acb_series([s, 1], prec=order)
    k = ((10 - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
    h = v.zeta(arb(1) / 5, deflate=True) + k * v.zeta(arb(2) / 5, deflate=True)
    h -= k * v.zeta(arb(3) / 5, deflate=True) + v.zeta(arb(4) / 5, deflate=True)
    return (-v * arb(5).log()).exp() * h


def dh_logder(s):
    v = dh_series(s, 2)
    return (arb(5) / arb.pi()).log() / 2 + ((s + 1) / 2).digamma() / 2 + v[1] / v[0]


def dh_completed(s):
    return ((arb(5) / arb.pi()).log() * s / 2).exp() * ((s + 1) / 2).gamma() * dh_series(s, 1)[0]


def dh_mp(s):
    k = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
    coeff = [1, k, -k, -1]
    if s == 1:
        return -sum(c * mp.digamma(mp.mpf(j) / 5) for j, c in enumerate(coeff, 1)) / 5
    return sum(c * mp.zeta(s, mp.mpf(j) / 5) for j, c in enumerate(coeff, 1)) / mp.power(5, s)


def basic_checks():
    vals = []
    for a, w in [(1, "0"), (1, "14.1347"), (1, "30.4249"), (1, "45"),
                 (2, "0"), (2, "14.1347"), (2, "25"), (1, "40")]:
        u = density(arb(a), arb(w))
        vals.append({"alpha": a, "w": w, "u_ball": bs(u)})
    a = acb_series([50, 1], prec=3)
    zz = a.zeta().log()
    lg = (a - 1).log() - a * arb.pi().log() / 2 + (a / 2 + 1).lgamma() + zz
    save("xi_values", vals)
    save("moments", {"mean_1_ball": bs(logder(acb(1))), "mean_50_ball": bs(lg[1]),
                     "variance_50_ball": bs(2 * lg[2]),
                     "saddle_mean_50_ball": bs((arb(54) / (2 * arb.pi())).log() / 2),
                     "saddle_variance_50_ball": bs(arb(1) / 108)})

    def phi(u):
        return 2 * mp.fsum((2 * mp.pi ** 2 * n ** 4 * mp.exp(mp.mpf('4.5') * u)
                           - 3 * mp.pi * n ** 2 * mp.exp(mp.mpf('2.5') * u))
                          * mp.exp(-mp.pi * n * n * mp.exp(2 * u)) for n in range(1, 13))
    rows = []
    for alpha in [mp.mpf('.5'), mp.mpf('1'), mp.mpf('2')]:
        integ = 2 * mp.quad(lambda u: phi(u) * mp.cosh((alpha - mp.mpf('.5')) * u), [0, .5, 1, 2, 3])
        rows.append({"alpha": str(alpha), "integral_computed": ms(integ), "xi_computed": ms(xi_mp(alpha)),
                     "difference_computed": ms(integ - xi_mp(alpha))})
    save("kernel_normalization", rows)
    c, w = mp.mpf(2), mp.mpf('1.3')
    integral = mp.quad(lambda x: 2 * mp.sin(w*x/2)**2 * mp.exp(-c*x)/x, [0,1,mp.inf])
    q, aa, gg = mp.mpf('1.3'), mp.mpf('.7'), mp.mpf('2.1')
    upair = lambda v: (aa/(aa*aa+(v-gg)**2)+aa/(aa*aa+(v+gg)**2))/mp.pi
    pairint = mp.quad(lambda v: upair(v)*mp.log(1+q*q/(v*v)), [0,1,gg,4,10,mp.inf])
    save("elementary_integrals", {"one_exponential_computed": ms(integral),
         "one_exponential_formula_computed": ms(mp.log(1+w*w/(c*c))/2),
         "pair_Thorin_integral_computed": ms(pairint),
         "pair_formula_computed": ms(mp.log(((q+aa)**2+gg**2)/(aa**2+gg**2)))})


def zero_and_prime_checks():
    zeros = acb.zeta_zeros(1, 300)
    save("zero_heights", {"first_ball": bs(zeros[0]), "number_300_ball": bs(zeros[-1])})
    T = mp.mpf(str(float(zeros[-1].imag.mid())))
    rows = []
    for w in ["0", "14.1347", "30.4249", "45"]:
        ww = arb(w)
        total = arb(0)
        for z in zeros:
            for sign in [1, -1]:
                total += (1 / (acb(1, ww) - acb(z.real, sign * z.imag))).real / arb.pi()
        v = mp.mpf(w)
        tail = mp.quad(lambda t: mp.log(t/(2*mp.pi))/(2*mp.pi**2) *
                       (mp.mpf('.5')/(mp.mpf('.25')+(v-t)**2) +
                        mp.mpf('.5')/(mp.mpf('.25')+(v+t)**2)), [T,2*T,mp.inf])
        rows.append({"w":w, "first_300_sum_ball":bs(total), "direct_minus_first_300_ball":bs(density(arb(1),ww)-total),
                     "smooth_counting_tail_computed":ms(tail),
                     "sum_plus_smooth_tail_computed":ms(mp.mpf(str(float(total.mid())))+tail)})
    save("zero_sum_checks", rows)
    bound = 200000
    sieve = bytearray(b'\x01') * (bound + 1)
    sieve[0:2] = b'\x00\x00'
    for p in range(2, int(bound**.5)+1):
        if sieve[p]:
            sieve[p*p:bound+1:p] = b'\x00' * ((bound-p*p)//p+1)
    powers = []
    for p in range(2, bound+1):
        if sieve[p]:
            n = p
            while n <= bound:
                powers.append((n, arb(p).log()))
                n *= p
    rows = []
    for w in ["0", "14.1347", "25"]:
        ww = arb(w)
        s = acb(2,ww)
        term = (1/s+1/(s-1)-arb.pi().log()/2+(s/2).digamma()/2).real/arb.pi()
        term -= sum((lp * (ww*arb(n).log()).cos()/(n*n) for n,lp in powers),arb(0))/arb.pi()
        rows.append({"w":w, "truncated_prime_formula_ball":bs(term),
                     "truncated_minus_exact_ball":bs(term-density(arb(2),ww))})
    save("prime_checks", {"cutoff":bound, "rows":rows,
         "uniform_tail_bound_ball":bs((arb(bound).log()+1)/(arb.pi()*bound))})
    return zeros


def convolution_and_grid(zeros):
    old = ctx.prec
    ctx.prec = 100
    heights = [float(z.imag.mid()) for z in zeros if float(z.imag.mid()) < 250]
    rows = []
    for w in [0.0,14.1347,40.0]:
        h = .25
        count = [0]
        def integrand(y):
            count[0] += 1
            val = float(density(arb('.75'),arb(y)).mid())
            return val * h/mp.pi.__float__() * (1/(h*h+(w-y)**2)+1/(h*h+(w+y)**2))
        pts = sorted(set([0.0,250.0,w]+heights+[max(0.0,w-.25),min(250.0,w+.25)]))
        val, err = quad(integrand,0,250,epsabs=2e-10,epsrel=2e-10,points=pts,limit=2500)
        exact = float(density(arb(1),arb(w)).mid())
        rows.append({"w":w, "window_250_computed":val, "scipy_error_estimate":err,
                     "exact_minus_window_computed":exact-val,"evaluations":count[0]})
    save("finite_window_convolution", rows)
    rows = []
    for a in ['.55','.75','1','1.5','2']:
        best = None
        for j in range(601):
            w = arb(j)/10
            u = density(arb(a),w)
            if best is None or u < best[0]:
                best=(u,j)
            assert u > 0
        rows.append({"alpha":a, "grid":"j/10, 0 <= j <= 600", "minimum_w":best[1]/10,
                     "minimum_ball":bs(best[0])})
    save("positive_grid_only",rows)
    ctx.prec=old


def dh_checks():
    root = mp.findroot(dh_mp,mp.mpc('.8085','85.6993'))
    # The exact decimal center is enclosed before bounding its surrounding box.
    cr, ci = ms(root.real), ms(root.imag)
    center = acb(cr,ci)
    r = arb('1e-20')
    box = acb(arb(center.real.mid(),r+center.real.rad()),arb(center.imag.mid(),r+center.imag.rad()))
    f = dh_series(center,3)
    fbox=dh_series(box,3)
    lhs = f[0].abs_upper()+(2*fbox[2]).abs_upper()*r*r/2
    rhs = f[1].abs_lower()*r
    assert lhs < rhs
    assert box.real > arb('.8')
    save("dh_root", {"center_real":cr,"center_imag":ci,"radius":"1e-20",
        "D_center_ball":bs(f[0]),"Dprime_abs_lower":bs(f[1].abs_lower()),
        "Dsecond_box_abs_upper":bs((2*fbox[2]).abs_upper()),
        "Rouche_left_ball":bs(lhs),"Rouche_right_ball":bs(rhs),
        "unique_zero_in_disk":True,"multiplicity":1})
    rows=[]
    for a in ['.5','.5001','.6','.75','1','2']:
        s=acb(a,ci)
        rows.append({"alpha":a,"ordinate":ci,"u_ball":bs(dh_logder(s).real/arb.pi())})
    save("dh_boundary_values", rows)
    s=acb('.7','10')
    save("dh_functional_equation", {"point":"0.7+10i", "H_s_ball":bs(dh_completed(s)),
         "H_1minus_s_ball":bs(dh_completed(1-s)),"difference_ball":bs(dh_completed(s)-dh_completed(1-s)),
         "H_half_ball":bs(dh_completed(acb('.5'))),
         "H_1_ball":bs(dh_completed(acb(1))),
         "positive_kernel_lower_factor_ball":bs(1-arb(5)/(2*arb.pi()))})
    k=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
    coeff=[0,1,k,-k,-1]
    theta=lambda t:mp.fsum(n*coeff[n%5]*mp.exp(-mp.pi*n*n*t/5) for n in range(1,31))
    ker=lambda u:2*mp.sqrt(mp.pi/5)*mp.exp(mp.mpf('1.5')*u)*theta(mp.exp(2*u))
    rows=[]
    for a in [mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]:
        integral=2*mp.quad(lambda u:ker(u)*mp.cosh((a-mp.mpf('.5'))*u),[0,.5,1,2,3])
        exact=mp.power(5/mp.pi,a/2)*mp.gamma((a+1)/2)*dh_mp(a)
        rows.append({"alpha":str(a),"kernel_integral_computed":ms(integral),"H_computed":ms(exact)})
    save("dh_kernel_integrals",rows)
    # Nonmultiplicativity: 2 and 3 are coprime, but a_6 != a_2 a_3.
    kk=((10-2*arb(5).sqrt()).sqrt()-2)/(arb(5).sqrt()-1)
    save("dh_no_euler_product",{"a_6":"1","a_2_times_a_3_ball":bs(-kk*kk)})


def main():
    start=time.monotonic()
    save("environment",{"python":sys.version,"platform":platform.platform(),"mpmath":mp.__version__,
         "python_flint":flint.__version__,"scipy":scipy.__version__,"ball_precision_bits":ctx.prec,
         "mpmath_decimal_digits":mp.mp.dps,"manuscript_sha256":hashlib.sha256((ROOT/'note.tex').read_bytes()).hexdigest()})
    basic_checks()
    print('Basic identities computed.',flush=True)
    zeros=zero_and_prime_checks()
    print('Zero and prime sums computed.',flush=True)
    convolution_and_grid(zeros)
    print('Convolution and finite grids computed.',flush=True)
    dh_checks()
    print('Davenport and Heilbronn ball checks computed.',flush=True)
    save("elapsed_seconds",time.monotonic()-start)
    with (OUT/'numbers.txt').open('w') as f:
        for key,value in results.items():
            f.write(key+'\n'+json.dumps(value,indent=2)+'\n\n')
    print('Wrote results/numbers.json and results/numbers.txt.',flush=True)


if __name__ == '__main__':
    main()
