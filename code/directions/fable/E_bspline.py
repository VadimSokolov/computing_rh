"""Track E, splines. Symmetric four knot B spline with knots -a,-b,b,a (a=1, r=b/a):
transform L(z) = 6[sinh(z) - sinh(rz)/r]/((1-r^2) z^3); on the imaginary axis the zeros are
those of f_r(theta) = sin(theta) - sin(r theta)/r away from the triple zero at 0.
(a) identities at r = 1/3, 1/2, 2/3; (b) scan r in (0,1): zeros in [0.5,60]x[-8,8] by the argument
principle against zeros in the thin box [0.5,60]x[-0.05,0.05]; (c) the limit r -> 0 and an asymmetric
hat function (knots 0,1,3) as checks of the obstruction by asymmetry. numpy only."""
import numpy as np
from fractions import Fraction

def f(th, r):
    return np.sin(th) - np.sin(r*th)/r

def winding(fun, x0, x1, y0, y1, m=40000):
    # number of zeros inside the rectangle, by the change of argument along its boundary
    xs = np.linspace(x0, x1, m); ys = np.linspace(y0, y1, m//4)
    path = np.concatenate([xs + 1j*y0, x1 + 1j*ys, xs[::-1] + 1j*y1, x0 + 1j*ys[::-1]])
    vals = fun(path)
    d = np.angle(vals[1:]/vals[:-1])
    if np.max(np.abs(d)) > 2.5:
        return -1   # sampling too coarse or a zero on the contour
    return int(round(d.sum()/(2*np.pi)))

print("(a) identities, a = 1, z random complex")
rng = np.random.default_rng(1)
for _ in range(3):
    z = complex(rng.normal(), rng.normal())*2
    L = lambda r: 6*(np.sinh(z) - np.sinh(r*z)/r)/((1-r*r)*z**3)
    sc = lambda c: np.sinh(c*z)/(c*z)
    print("  r=1/3: L=%.12g%+.12gi  U^3=%.12g%+.12gi" % (L(1/3).real, L(1/3).imag, (sc(1/3)**3).real, (sc(1/3)**3).imag))
    print("  r=1/2: L=%.12g%+.12gi  U_1/2 U_1/4^2=%.12g%+.12gi" % (L(1/2).real, L(1/2).imag, (sc(1/2)*sc(1/4)**2).real, (sc(1/2)*sc(1/4)**2).imag))
    spin = (4*np.cosh(z/3) + 1)/5     # S in {-1/3,0,1/3} with weights 2/5,1/5,2/5
    v = sc(1/3)*sc(1/6)**2*spin
    print("  r=2/3: L=%.12g%+.12gi  U_1/3 U_1/6^2 S=%.12g%+.12gi" % (L(2/3).real, L(2/3).imag, v.real, v.imag))

print("(b) scan: r, zeros in the full box, zeros in the thin box, all real?")
rs = sorted(set(list(np.round(np.linspace(0.01, 0.99, 197), 6)) + [float(Fraction(p, q)) for q in range(2, 9) for p in range(1, q)]))
allreal = []
for r in rs:
    fun = lambda w: f(w, r)
    nfull = winding(fun, 0.5, 60.0, -8.0, 8.0)
    nthin = winding(fun, 0.5, 60.0, -0.05, 0.05)
    ok = (nfull == nthin and nfull >= 0)
    if ok:
        allreal.append((r, nfull))
    if abs(r - round(r*12)/12) < 1e-9 or ok:
        print("  r=%.6f  full=%3d  thin=%3d  %s" % (r, nfull, nthin, "ALL REAL" if ok else ""))
print("  ratios with all zeros real (to 0.05 in the imaginary part), among %d tested:" % len(rs))
print("  ", [(("%.6f" % r), n) for r, n in allreal])
# refine the found ratios with a thinner box and a longer range
print("  refinement, thin box 0.002, range [0.5,120]:")
for r, _ in allreal:
    fun = lambda w: f(w, r)
    print("   r=%.6f full=%d thin=%d" % (r, winding(fun, 0.5, 120.0, -8.0, 8.0, 80000), winding(fun, 0.5, 120.0, -0.002, 0.002, 80000)))

print("(c) limit r -> 0: zeros of sin(theta) - theta in [0.5,60]x[-8,8] vs thin box")
g = lambda w: np.sin(w) - w
print("  full=%d thin=%d" % (winding(g, 0.5, 60.0, -8.0, 8.0), winding(g, 0.5, 60.0, -0.05, 0.05)))
print("  asymmetric hat, knots 0,1,3, centred: phi(theta) = e^{-4i theta/3}(2 - 3e^{i theta} + e^{3i theta})/(3 (i theta)^2)")
h = lambda w: np.exp(-4j*w/3)*(2 - 3*np.exp(1j*w) + np.exp(3j*w))/(3*(1j*w)**2)
th = np.linspace(0.3, 60, 200000)
print("  min |phi| on [0.3,60] = %.3e (no real zero if not small); zeros in [0.3,60]x[-8,8] = %d" % (np.min(np.abs(h(th))), winding(h, 0.3, 60.0, -8.0, 8.0)))
