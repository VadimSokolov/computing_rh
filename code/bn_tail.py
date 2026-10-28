"""Saddlepoint tails of H and the two explicit formulas for its law.
(a) Lugannani and Rice tail probabilities against exact Talbot inversion of (1 - phi)/s and phi/s.
(b) zero side: P(H > y) = sum_j (a_j/gamma_j^2) exp(-gamma_j^2 y), a_j = -2 xi(1/2) gamma_j/eta_j; convergence in J.
(c) prime side: f_H(y) = sum_n mu(n) n^{-1/2} g_n(y), g_n = L^{-1}[phi_G(s) e^{-sqrt s log n}] (Bromwich line right of the pole at 1/4),
    phi_G = xi(1/2)/G(1/2 + sqrt s), G(s) = s(s-1) pi^{-s/2} Gamma(s/2)/2.
(d) where the saddlepoint sits: y < E H at 1/2 + w on the real axis (w = sqrt(-theta)); y > E H at 1/2 + i sqrt(theta).
Usage: python3 bn_tail.py PART [DATADIR]  (DATADIR holds zeros500.json)"""
import mpmath as mp, json, sys, os
part = sys.argv[1]
if len(sys.argv) > 2: os.chdir(sys.argv[2])
def xi(s): return 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def G(s): return 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)
mp.mp.dps = 30
h = mp.mpf(1)/2; x12 = xi(h)
phi = lambda s: x12/xi(h + mp.sqrt(s))
K = lambda th: -mp.log(mp.re(xi(h + mp.sqrt(-mp.re(th))))/x12)
Z = [mp.mpf(v) for v in json.load(open("zeros500.json"))]
def saddle(y, lo=-mp.mpf(10)**6, hi=mp.mpf("199.78")):
    """K' is increasing on (-inf, gamma_1^2): bisection for K'(theta) = y"""
    for _ in range(90):
        mid = (lo + hi)/2
        if mp.diff(K, mid) < y: lo = mid
        else: hi = mid
    return (lo + hi)/2
ys = [0.008, 0.012, 0.016, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08]
res = {}
if part == "a":
    rows = []; th = mp.mpf(0)
    for y in ys:
        y = mp.mpf(y)
        th = saddle(y)
        K2 = mp.diff(K, th, 2)
        r = mp.sign(th)*mp.sqrt(2*(th*y - K(th))); u = th*mp.sqrt(K2)
        lr = mp.ncdf(-r) + mp.npdf(r)*(1/u - 1/r)          # P(H > y)
        up = mp.invertlaplace(lambda s: (1 - phi(s))/s, y, method="talbot")
        lo = mp.invertlaplace(lambda s: phi(s)/s, y, method="talbot")
        rows.append([float(y), float(th), float(up), float(lr), float(lo), float(1 - lr)])
        print(rows[-1], "upper ratio", float(lr/up), "lower ratio", float((1 - lr)/lo), flush=True)
    res["rows"] = rows
if part == "b":
    Xi = lambda t: mp.re(xi(h + 1j*t))
    J = 120
    c = []
    for j in range(J):
        eta = mp.diff(Xi, Z[j]); a = -2*x12*Z[j]/eta
        c.append(a/Z[j]**2)
    res["a_over_g2"] = [mp.nstr(v, 10) for v in c[:12]]
    rows = []
    for y in ys:
        y = mp.mpf(y)
        ex = mp.invertlaplace(lambda s: (1 - phi(s))/s, y, method="talbot")
        part_sums = {}
        for Jn in [1, 2, 5, 10, 20, 50, 120]:
            part_sums[Jn] = float(mp.fsum(c[j]*mp.e**(-Z[j]**2*y) for j in range(Jn)))
        need = next((Jn for Jn in [1, 2, 5, 10, 20, 50, 120] if abs(part_sums[Jn] - ex) < 1e-8*abs(ex)), None)
        rows.append({"y": float(y), "exact": float(ex), "partial": part_sums, "J_for_1e-8": need})
        print(rows[-1], flush=True)
    res["rows"] = rows
if part == "c":
    phiG = lambda s: x12/G(h + mp.sqrt(s))
    def hn(n, y):
        a = mp.log(n)
        # g_n = L^{-1}[phi_G(s) e^{-sqrt s log n}]; phi_G has a pole at s = 1/4, so shift the contour by 1
        F = lambda s: phiG(s + 1)*mp.e**(-mp.sqrt(s + 1)*a)
        return mp.e**y*mp.invertlaplace(F, y, method="talbot")
    sqf = [n for n in range(1, 400) if mp.moebius(n) != 0] if hasattr(mp, "moebius") else None
    def mobius(n):
        m, p, k = n, 2, 1
        while p*p <= m:
            if m % p == 0:
                m //= p
                if m % p == 0: return 0
                k = -k
            p += 1
        return -k if m > 1 else k
    sqf = [n for n in range(1, 400) if mobius(n) != 0]
    rows = []
    for y in [0.008, 0.012, 0.016, 0.02, 0.03, 0.04]:
        y = mp.mpf(y)
        ex = mp.invertlaplace(phi, y, method="talbot")
        acc = mp.mpf(0); trace = []
        for n in sqf:
            t = mobius(n)*hn(n, y)/mp.sqrt(n)
            acc += t
            trace.append([n, float(t), float(acc)])
            if n > 3 and abs(t) < 1e-14*abs(ex): break
        rows.append({"y": float(y), "exact": float(ex), "sum": float(acc), "relerr": float((acc - ex)/ex), "terms": trace[:8], "last_n": trace[-1][0]})
        print(rows[-1]["y"], rows[-1]["exact"], rows[-1]["sum"], rows[-1]["relerr"], rows[-1]["last_n"], flush=True)
    res["rows"] = rows
if part == "d":
    rows = []; th = mp.mpf(-20000)
    for y in [0.004, 0.006, 0.01, 0.015, 0.02, 0.0231, 0.03, 0.05, 0.08, 0.12, 0.2]:
        y = mp.mpf(y)
        th = saddle(y)
        if th < 0:
            w = mp.sqrt(-th); where = "real axis at 1/2 + %.3f" % float(w)
            primes = [float(mp.log(p)/mp.log(p)*p**(-h - w)) for p in [2, 3, 5]]
        else:
            w = mp.sqrt(th); where = "critical line at height %.4f" % float(w); primes = None
        rows.append({"y": float(y), "theta": float(th), "where": where, "p^(-1/2-w) for p=2,3,5": primes})
        print(rows[-1], flush=True)
    res["rows"] = rows
json.dump(res, open("bn_tail_%s.json" % part, "w"), indent=1, default=str)
