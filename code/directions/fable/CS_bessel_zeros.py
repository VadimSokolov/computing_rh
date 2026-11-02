# Track CS. First Lee and Yang zero alpha_1 of rho_a (zeros of z -> K_{iz}(a) = int_0^oo e^{-a cosh t} cos(zt) dt)
# and of rho*, and certified non vanishing at 2 alpha_1 and 3 alpha_1 (the uniform and Ising progressions).
# python-flint balls; the transforms are evaluated by rigorous quadrature (acb.integral) on [0, T] with a
# written tail bound, because Arb's bessel_k with a ball order near the imaginary axis loses all precision.
import sys
sys.path.insert(0, '/tmp/cs_lib')
from flint import arb, acb, ctx
ctx.dps = 25
PI = arb.pi()

def transform(weight, T, tail):
    """z -> 2 int_0^T weight(t) cos(z t) dt + [+/- tail], for a real ball z."""
    def f(z):
        zc = acb(z)
        val = acb.integral(lambda t, a: weight(t) * (zc * t).cos(), acb(0), acb(T))
        return val.real * 2 + arb(0, tail)
    return f

def rho_a_weight(a):
    return lambda t: (-a * t.cosh()).exp()

# tails: for t >= T, e^{-a cosh t} <= e^{-(a/2) e^t}, and int_T^oo e^{-(a/2)e^t} dt <= (2/(a e^T)) e^{-(a/2) e^T}
def tail_a(a, T):
    eT = arb(T).exp()
    return float(((2 / (a * eT)) * (-(a / 2) * eT).exp()).str(5, radius=False)) * 2 + 1e-200

f_1 = transform(rho_a_weight(arb(1)), 6, tail_a(arb(1), 6))
f_2pi = transform(rho_a_weight(2 * PI), 4, tail_a(2 * PI, 4))
# rho*: weight (2 pi^2 cosh 9u - 3 pi cosh 5u) e^{-2 pi cosh 4u}; for u >= 1.5 it is below 2 pi^2 e^{9u - pi e^{4u}} < 1e-500
def star_weight(u):
    return (2 * PI**2 * (9 * u).cosh() - 3 * PI * (5 * u).cosh()) * (-2 * PI * (4 * u).cosh()).exp()
f_star = transform(star_weight, 1.5, 1e-200)

# cross check against Arb's Bessel K at an exact point: 2 K_{i z}(a) = int_R e^{-a cosh t} e^{izt} dt
chk = acb(0, 1).bessel_k(acb(1)); chk2 = acb(1).bessel_k(acb(0, 1))
print("cross check at z=1, a=1: quadrature", f_1(arb(1)).str(12), " Arb K_i(1) candidates", chk.str(12), chk2.str(12))

def lower(x): return x.mid() - x.rad()
def upper(x): return x.mid() + x.rad()
def pos(x): return (lower(x).mid() - lower(x).rad()) > 0
def neg(x): return (upper(x).mid() + upper(x).rad()) < 0
def nonzero(x): return pos(x) or neg(x)
def ball(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return arb((lo + hi) / 2, (hi - lo) / 2)

def certify_positive(f, lo, hi, width=0.05, depth=0):
    n = max(1, int((hi - lo) / width))
    for i in range(n):
        a_, b_ = lo + (hi - lo) * i / n, lo + (hi - lo) * (i + 1) / n
        if not pos(f(ball(a_, b_))):
            assert depth < 4, ("not certified positive near", a_, b_, f(ball(a_, b_)))
            certify_positive(f, a_, b_, width / 8, depth + 1)

def first_zero(f, step=0.05, zmax=60.0):
    z = 0.0
    assert pos(f(arb(0)))
    while z < zmax:
        v = f(arb(z + step))
        if neg(v):
            lo, hi = z, z + step
            break
        assert pos(v), (z + step, v)
        z += step
    else:
        raise SystemExit("no sign change below zmax")
    certify_positive(f, 0.0, lo)
    for _ in range(50):
        m = (lo + hi) / 2
        v = f(arb(m))
        if pos(v): lo = m
        elif neg(v): hi = m
        else: break
    return lo, hi

for name, f in [("rho_1 (a=1)", f_1), ("rho_2pi (a=2pi)", f_2pi), ("rho_star", f_star)]:
    lo, hi = first_zero(f)
    print("\n%s: f > 0 on [0, %.4f] certified; first zero alpha_1 in [%.12f, %.12f]" % (name, lo, lo, hi))
    for k in (2, 3):
        v = f(ball(k * lo, k * hi))
        print("  f(%d alpha_1) = %s   certified nonzero: %s" % (k, v.str(8), nonzero(v)))
        if not nonzero(v) and name.startswith("rho_"):
            # the value is below the quadrature tolerance: use Arb's Bessel K at the bracket, 2 K_{iz}(a),
            # at higher precision (the order i z is far from the integers here)
            a_val = arb(1) if "a=1" in name else 2 * PI
            old = ctx.dps; ctx.dps = 80
            w = 2 * acb(a_val).bessel_k(acb(0, ball(k * lo, k * hi))).real
            ctx.dps = old
            print("  via Arb K at 80 digits: 2 K_{i %d alpha_1}(a) = %s   certified nonzero: %s" % (k, w.str(8), nonzero(w)))
    # further sign changes (computed, exact points, for information)
    zs = []; zp = lo + 1e-9; prev = -1.0; g = 0.01
    while zp < 40 and len(zs) < 8:
        zn = zp + g
        cur = 1.0 if pos(f(arb(zn))) else -1.0
        if prev * cur < 0: zs.append(zn - g / 2)
        prev, zp = cur, zn
    print("  next sign changes (computed): %s" % ["%.3f" % z for z in zs])
