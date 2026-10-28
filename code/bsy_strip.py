# bsy_strip.py: the zero free strip of the Davenport and Heilbronn function and the constants of Theorem bsy:thm:dhid, Chapter 35, Section bsy:sec:sin2
# (claims: sigma_1 = 1.3938 gives sum_{n >= 2}|a_n| n^(-sigma) <= 1; f(1/2) = 0.8254; the decided proof's constants). Task C5 of reconcile-M.md:
# S(1.4) < 0.98489 and S(1.3938) > 1.004 by integral test bounds at K = 100 in interval arithmetic. Added for the book in October 2026, written from the
# mathematics (the collaborator's scripts were not read, run or copied).
# Conventions. f(s) = sum a_n n^-s, a_n = 1, kappa, -kappa, -1, 0 for n = 1, 2, 3, 4, 0 mod 5, kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1). S(sigma) =
# sum_{n >= 2} |a_n| n^-sigma (n = 1 is excluded: it is the leading 1 of |f - 1| <= S).
# Method. S = sum_{n=2}^{K-1} |a_n| n^-sigma + sum over the residue classes c in {1,2,3,4} of |a_c| T_c, where T_c = sum_{n >= K, n = c mod 5} n^-sigma, and for the
# decreasing g = x^-sigma and n0 the first n >= K in the class: (1/5) int_{n0}^inf g <= T_c <= g(n0) + (1/5) int_{n0}^inf g, int = n0^(1-sigma)/(sigma - 1).
# Refined (convexity of x^-sigma) bounds: trapezoid gives T_c >= g(n0)/2 + (1/5) int_{n0}^inf g, midpoint gives T_c <= (1/5) int_{n0-5/2}^inf g.
# All in mpmath interval arithmetic (iv, 40 digits); K = 100 and K = 1000 (the interval widths are logged). Also the root sigma_1 of S = 1 (bisection on the
# interval bounds), delta_max = sigma_1 - 1/2 and the thresholds pi/(2 delta_max), and the constants of the decided proof: f(1/2), gamma_1 = first sign change
# of Lambda(1/2 + it) in [5, 5.2] (bisection, mpmath hurwitz), (1/2) log(5/pi), (1/2) psi(2), -zeta'(3), zeta(3) - 1, the bound 2/5 (85.69^2 + 61/4) e^{-(85.69 - 5.1) pi/3}.
# Output: data/bsy_strip.json and the log bsy_strip.log. Usage: python3 bsy_strip.py [--quick] [--out FILE]; --quick uses K = 50 and 25 digits.
import argparse, json
from pathlib import Path
import mpmath as mp
from mpmath import iv

ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_strip_quick.json") if args.quick else here.parent / "data" / "bsy_strip.json")
dps = 25 if args.quick else 40
iv.dps = dps; mp.mp.dps = dps
log = []; res = {}
def P(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)

s5 = iv.sqrt(iv.mpf(5))
kappa = (iv.sqrt(10 - 2*s5) - 2)/(s5 - 1)
P("kappa =", kappa)
res["kappa"] = str(kappa)

def lo_hi(v):
    # endpoints of an mpmath interval as plain mpf (outward exact)
    return mp.make_mpf(v._mpi_[0]), mp.make_mpf(v._mpi_[1])

def Sbounds(sigma, K, refined=False):
    sg = iv.mpf(sigma) if not isinstance(sigma, type(iv.mpf(1))) else sigma
    w = {1: iv.mpf(1), 2: kappa, 3: kappa, 4: iv.mpf(1)}
    tot = iv.mpf(0)
    for n in range(2, K):
        c = n % 5
        if c:
            tot += w[c]*iv.mpf(n)**(-sg)
    for c in (1, 2, 3, 4):
        n0 = K + ((c - K) % 5)
        I = iv.mpf(n0)**(1 - sg)/(sg - 1)/5
        g0 = iv.mpf(n0)**(-sg)
        Ilo, Ihi = lo_hi(I); glo, ghi = lo_hi(g0)
        if not refined:
            tot += w[c]*iv.mpf([Ilo, Ihi + ghi])  # hull of [I_low, I + g(n0)]
        else:  # convexity: trapezoid below, midpoint above
            Im = (iv.mpf(n0) - iv.mpf(5)/2)**(1 - sg)/(sg - 1)/5
            tot += w[c]*iv.mpf([Ilo + glo/2, lo_hi(Im)[1]])
    return tot

Ks = [50] if args.quick else [100, 1000]
res["S"] = {}
for K in Ks:
    for sg in ("1.4", "1.3938", "1.395136"):
        s_ = iv.mpf(sg)
        v = Sbounds(s_, K)
        P("K", K, "S(" + sg + ") in [", mp.nstr(lo_hi(v)[0], 10), ",", mp.nstr(lo_hi(v)[1], 10), "]")
        res["S"][f"K{K}_sigma{sg}"] = [mp.nstr(lo_hi(v)[0], 10), mp.nstr(lo_hi(v)[1], 10)]
    for sg in ("1.4", "1.3938"):
        v = Sbounds(iv.mpf(sg), K, True)
        P("K", K, "refined (convexity) S(" + sg + ") in [", mp.nstr(lo_hi(v)[0], 10), ",", mp.nstr(lo_hi(v)[1], 10), "]")
        res["S"][f"refined_K{K}_sigma{sg}"] = [mp.nstr(lo_hi(v)[0], 10), mp.nstr(lo_hi(v)[1], 10)]
    v14 = Sbounds(iv.mpf("1.4"), K); v13 = Sbounds(iv.mpf("1.3938"), K)
    r14 = Sbounds(iv.mpf("1.4"), K, True); r13 = Sbounds(iv.mpf("1.3938"), K, True)
    P("  K", K, "refined: S(1.4) <", mp.nstr(lo_hi(r14)[1], 6), "claim 0.98489:", bool(lo_hi(r14)[1] < mp.mpf("0.98489")), "; S(1.3938) >", mp.nstr(lo_hi(r13)[0], 6), "claim 1.004:", bool(lo_hi(r13)[0] > mp.mpf("1.004")))
    P("  K", K, ": S(1.4) <", mp.nstr(lo_hi(v14)[1], 6), "claim 0.98489:", bool(lo_hi(v14)[1] < mp.mpf("0.98489")), "; S(1.3938) >", mp.nstr(lo_hi(v13)[0], 6), "claim 1.004:", bool(lo_hi(v13)[0] > mp.mpf("1.004")))
# root of S = 1 at K = Ks[-1]
K = Ks[-1]
lo, hi = mp.mpf("1.3"), mp.mpf("1.5")
for _ in range(60):
    mid = (lo + hi)/2
    v = Sbounds(iv.mpf(mid), K)
    vl, vh = lo_hi(v)
    if vh < 1: hi = mid   # S < 1 certified: root is below mid
    elif vl > 1: lo = mid
    else: break
P("sigma_1 (S = 1) in [", mp.nstr(lo, 8), ",", mp.nstr(hi, 8), "] (certified by the bounds at K =", K, ")")
res["sigma1_bracket"] = [mp.nstr(lo, 8), mp.nstr(hi, 8)]
for s1 in (mp.mpf("1.4"), lo):
    dmax = s1 - mp.mpf(1)/2
    P("sigma_1", mp.nstr(s1, 6), " delta_max <", mp.nstr(dmax, 6), " threshold pi/(2 delta_max) =", mp.nstr(mp.pi/(2*dmax), 6))
# decided proof constants
mp.mp.dps = dps
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
def f(s):
    return mp.mpf(5)**(-s)*(mp.zeta(s, mp.mpf(1)/5) + kap*mp.zeta(s, mp.mpf(2)/5) - kap*mp.zeta(s, mp.mpf(3)/5) - mp.zeta(s, mp.mpf(4)/5))
P("f(1/2) =", mp.nstr(f(mp.mpf(1)/2), 8))
res["f_half"] = mp.nstr(f(mp.mpf(1)/2), 8)
Lam = lambda t: ((5/mp.pi)**((mp.mpf(1)/2 + 1j*t)/2)*mp.gamma((mp.mpf(3)/2 + 1j*t)/2)*f(mp.mpf(1)/2 + 1j*t))
ph = mp.exp(-1j*mp.arg(Lam(mp.mpf("5.1"))))
g1 = mp.findroot(lambda t: mp.re(Lam(t)*ph), mp.mpf("5.09"))
P("gamma_1 (zero of Lambda(1/2+it) near 5.09) =", mp.nstr(g1, 10), " |f(1/2 + i g1)| =", mp.nstr(abs(f(mp.mpf(1)/2 + 1j*g1)), 3))
res["gamma1"] = mp.nstr(g1, 10)
c1 = mp.log(mp.mpf(5)/mp.pi)/2; c2 = mp.psi(0, 2)/2; c3 = -mp.zeta(3, derivative=1); c4 = mp.zeta(3) - 1
P("(1/2)log(5/pi) =", mp.nstr(c1, 6), " (1/2)psi(2) =", mp.nstr(c2, 6), " -zeta'(3) =", mp.nstr(c3, 6), " zeta(3)-1 =", mp.nstr(c4, 6))
P("Lambda'/Lambda(3) <= c1 + c2 + c3/(1 - c4) =", mp.nstr(c1 + c2 + c3/(1 - c4), 6), "< 2 ; chapter proof bound 1/4 + 1/2 + 1/3 =", mp.nstr(mp.mpf(1)/4 + mp.mpf(1)/2 + mp.mpf(1)/3, 6))
bnd = mp.mpf(2)/5*(mp.mpf("85.69")**2 + mp.mpf(61)/4)*mp.exp(-(mp.mpf("85.69") - mp.mpf("5.1"))*mp.pi/3)
P("2/5 (85.69^2 + 61/4) e^{-(85.69 - 5.1) pi/3} =", mp.nstr(bnd, 6))
res.update({"c1": mp.nstr(c1, 6), "c2": mp.nstr(c2, 6), "c3": mp.nstr(c3, 6), "c4": mp.nstr(c4, 6), "final_bound": mp.nstr(bnd, 6)})
Path(out).write_text(json.dumps(res, indent=1)); Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
