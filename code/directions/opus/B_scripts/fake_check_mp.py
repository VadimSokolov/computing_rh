#!/usr/bin/env python3
"""High precision checks of Theorem 1 of B.md with the closed form vector. mpmath, floating point, NOT certified.

(i) self duality of the closed form, residual in 60 digits; ratio to the shifted Gaussian (-1)^n exp(-pi (n-q/2)^2/q).
(ii) functional equation M(s) = M(1-s) at sample points. (iii) the bound of Theorem 1(e) at sample points.
(iv) refinement of the zeros with Re s > 1 from fake_exact_N*.json, with a numerical winding number.
Usage: python3 -I fake_check_mp.py
"""
import json, os, time
import mpmath as mp
mp.mp.dps = 60
here = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()

def vector(N):
    r = 4 * N; q = r * r; K = N * r; half = q // 2
    y = [mp.cos(2 * mp.pi * n / q) for n in range(K, half + 1)]
    w = []
    for i, n in enumerate(range(K, half + 1)):
        p = mp.mpf(1)
        for j in range(len(y)):
            if j != i: p *= y[i] - y[j]
        w.append(1 / ((1 if n == half else 2) * p))
    c = [mp.mpf(0)] * q
    for i, n in enumerate(range(K, half + 1)):
        c[n] = w[i] / w[-1]; c[(q - n) % q] = c[n]
    return r, q, K, c

def rL(s, r, q, c):           # r^s L_c(s) = (r/q)^s sum_b c_b zeta(s, b/q)
    tot = mp.mpf(0)
    for b in range(1, q + 1):
        v = c[b % q]
        if v != 0: tot += v * mp.zeta(s, mp.mpf(b) / q)
    return (mp.mpf(r) / q) ** s * tot
def M(s, r, q, c):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * rL(s, r, q, c)

print("(i) self duality and shape")
for N in (2, 3, 4, 5, 6):
    r, q, K, c = vector(N)
    cs = [mp.cos(2 * mp.pi * k / q) for k in range(q)]
    res = mp.mpf(0)
    for m in range(q // 2 + 1):
        v = sum(c[b] * cs[(m * b) % q] for b in range(q) if c[b] != 0) / mp.sqrt(q) - c[m]
        res = max(res, abs(v))
    half = q // 2
    alt = all(mp.sign(c[n]) == (-1) ** n for n in range(K, half + 1))
    rat = {j: float(abs(c[half - j]) / mp.exp(-mp.pi * j * j / q)) for j in (0, r // 4, r // 2, r, 3 * r // 2, half - K)}
    print(f" N={N} q={q}: max|Fc-c|={mp.nstr(res,3)}, c_n=(-1)^n|c_n| on [K,q/2]: {alt}, max|c|={mp.nstr(max(abs(x) for x in c),6)} at n={max(range(q),key=lambda n:abs(c[n]))}, |c_K|={mp.nstr(abs(c[K]),4)}, exp(-pi N^2)={mp.nstr(mp.exp(-mp.pi*N*N),4)}")
    print("    |c_{q/2-j}|/exp(-pi j^2/q) at j:", {k: round(v, 4) for k, v in rat.items()}, f"{time.time()-T0:.0f}s")

print("(ii) functional equation, relative difference of M(s), M(1-s), 60 digits")
for N, pts in ((2, [mp.mpc(1.3, 88), mp.mpc(0.7, 20)]), (3, [mp.mpc(1.2, 180)])):
    r, q, K, c = vector(N)
    for s in pts:
        a, b = M(s, r, q, c), M(1 - s, r, q, c)
        print(f" N={N} s={mp.nstr(s,5)}: |M(s)|={mp.nstr(abs(a),4)}, rel diff={mp.nstr(abs(a-b)/abs(a),3)}", f"{time.time()-T0:.0f}s")

print("(iii) |M| over the bound of Theorem 1(e)")
def bound(s, N):
    t = s.imag; nu = t / (2 * mp.pi * N * N); om = mp.acos(nu); so = mp.sin(om)
    g = mp.sqrt(1 - nu * nu) - nu * om
    return abs(s * (s - 1)) * (1 + 2 / (mp.pi * so)) / (mp.pi * N * N * so) * mp.exp(-mp.pi * t / 4) * mp.exp(-mp.pi * N * N * g)
for N in (2, 3):
    r, q, K, c = vector(N)
    for nu in (0.1, 0.5, 0.9):
        t = nu * 2 * mp.pi * N * N
        row = []
        for sg in (-1, 0.5, 2):
            s = mp.mpc(sg, t)
            row.append(mp.nstr(abs(M(s, r, q, c)) / bound(s, N), 3))
        print(f" N={N} nu={nu} t={float(t):.2f}: sigma=-1,1/2,2 ratios", row, f"{time.time()-T0:.0f}s")

print("(iv) zeros with Re s > 1 from the double precision scan, refined at 40 digits, winding on a square of half width 1e-6")
mp.mp.dps = 40
for N in (2, 3, 4):
    d = json.load(open(os.path.join(here, f"fake_exact_N{N}.json")))
    r, q, K, c = vector(N)
    for key in ("1.0", "-1.0"):
        zs = [z for z in d["offline_scan"].get(key, []) if z[0] > 1.0]
        if not zs: 
            print(f" N={N} eps={key}: none in the scan"); continue
        z = zs[0]; eps = mp.mpf(key)
        f = lambda s: mp.zeta(s) + eps * rL(s, r, q, c)
        root = mp.findroot(f, (mp.mpc(z[0], z[1]), mp.mpc(z[0] + 1e-6, z[1] + 1e-6)), solver="secant", tol=1e-30, verify=False, maxsteps=12)
        h = mp.mpf("1e-6"); cor = [root + h * mp.mpc(a, b) for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
        total = mp.mpf(0); prev = f(cor[0])
        for k in range(1, 5):
            cur = f(cor[k % 4]); total += mp.arg(cur / prev); prev = cur
        print(f" N={N} eps={key}: start {z}, root {mp.nstr(root,16)}, |f|={mp.nstr(abs(f(root)),3)}, winding={mp.nstr(total/(2*mp.pi),6)}", f"{time.time()-T0:.0f}s")
print(f"total {time.time()-T0:.0f}s")
