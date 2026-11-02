# Track CS. Certify U*''' > 0 on (0, 2.5] for U*(u) = 2 pi cosh 4u - log(2 pi^2 cosh 9u - 3 pi cosh 5u)
# (Ellis, Monroe, Newman condition for rho*), with python-flint ball arithmetic and power series.
# Strategy: U* is even so U*'''(0) = 0 exactly; certify U*'''' > 0 on [0, d0] (then U*''' > 0 on (0, d0]),
# and U*''' > 0 on [d0, 2.5] by adaptive bisection. Beyond u = 2 an analytic bound is in CS.md.
import sys
sys.path.insert(0, '/tmp/cs_lib')
from flint import arb, arb_series, ctx
ctx.dps = 40
ctx.cap = 6
PI = arb.pi()

def coeffs(s, n):
    try:
        return [s[i] for i in range(n)]
    except Exception:
        c = s.coeffs()
        return list(c) + [arb(0)] * (n - len(c))

def cosh_s(s):
    return (s.exp() + (-s).exp()) / 2

def taylor(u):
    """Taylor coefficients c_0..c_5 of U* at the ball u (valid for every point of the ball)."""
    s = arb_series([u, 1])
    h = 2 * PI**2 * cosh_s(9 * s) - 3 * PI * cosh_s(5 * s)
    U = 2 * PI * cosh_s(4 * s) - h.log()
    return coeffs(U, 6)

def lower(x):
    return x.mid() - x.rad()

def positive(x):
    lb = lower(x)
    return (lb.mid() - lb.rad()) > 0

def ball(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return arb((lo + hi) / 2, (hi - lo) / 2 + arb(0).rad())

# derivative values at 0 (exact point)
c = taylor(arb(0))
print("U*''(0)  =", 2 * c[2])
print("U*'''(0) =", 6 * c[3])
print("U*''''(0)=", 24 * c[4])

# Step 1: U*'''' > 0 on [0, d0], by adaptive bisection on the fourth coefficient
d0 = arb(1) / 64
stack4 = [(0.0, float(d0))]
ok4 = True; n4 = 0; min4 = None
while stack4:
    lo, hi = stack4.pop()
    v = 24 * taylor(ball(lo, hi))[4]
    if positive(v):
        n4 += 1
        lb = float(lower(v).mid().str(10, radius=False))
        min4 = lb if min4 is None else min(min4, lb)
        continue
    if hi - lo < 1e-7:
        ok4 = False; print("FAILED (4th derivative) on", lo, hi, v); break
    m = (lo + hi) / 2
    stack4.append((lo, m)); stack4.append((m, hi))
print("U*'''' > 0 on [0, 1/64]: certified %s with %d intervals, least lower bound %s" % (ok4, n4, min4))

# Step 2: U*''' > 0 on [d0, 2.5] by bisection
stack = [(float(d0), 2.5)]
n_int = 0
minw = 1e-6
worst = None
failed = False
while stack:
    lo, hi = stack.pop()
    c = taylor(ball(lo, hi))
    v = 6 * c[3]
    if positive(v):
        n_int += 1
        lb = float(lower(v).mid().str(10, radius=False))
        if worst is None or lb < worst[0]:
            worst = (lb, lo, hi)
        continue
    if hi - lo < minw:
        failed = True
        print("FAILED on", lo, hi, v)
        break
    m = (lo + hi) / 2
    stack.append((lo, m)); stack.append((m, hi))
print("intervals certified:", n_int, "failed:", failed)
print("smallest certified lower bound of U*''' and where:", worst)
print("RESULT:", "U*''' > 0 on (0, 2.5] CERTIFIED" if (ok4 and not failed) else "NOT certified")
# grid values for the record
for u in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]:
    c = taylor(arb(u))
    print("u=%.2f  U*'''=%s" % (u, (6 * c[3]).str(8, radius=False)))
