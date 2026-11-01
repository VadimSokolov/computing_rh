# A planted pair of multiplicity n: remove gamma_k, gamma_{k+1} from Xi, plant (z^2-c^2)^n (z^2-cbar^2)^n, c = m + i eta,
# m the midpoint. The tilt F(x+iw)+F(x-iw) = 2|F(x+iw)| cos arg F(x+iw) on the real axis; count its real zeros in
# (m-r, m+r) from the argument (each crossing of pi/2 + pi Z), and find the least w with all 2n local zeros real, and the
# least w with at least two real. Compare with eta sin((n-1)pi/(2n) + |psi|/n) (n even) and the eta^3 law (n odd).
import mpmath as mp, json, math
mp.mp.dps = 20
G = [mp.im(mp.zetazero(j)) for j in range(1, 20)]
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
k = 8; gk, gk1 = G[k-1], G[k]; m = (gk+gk1)/2; r = mp.mpf('0.8')
def Gfun(z): return xi(mp.mpf(1)/2+1j*z)/((z*z-gk*gk)*(z*z-gk1*gk1))
xs = [m - r + 2*r*i/800 for i in range(801)]
Gcache = {}
def count(n, eta, w):
    c = m + 1j*eta
    if w not in Gcache: Gcache[w] = [Gfun(x + 1j*w) for x in xs]
    tot = 0; prev = None
    for x, g in zip(xs, Gcache[w]):
        z = x + 1j*w
        v = g*((z*z-c*c)*(z*z-mp.conj(c)**2))**n
        s = mp.re(v)   # sign of cos arg F
        if prev is not None and s*prev < 0: tot += 1
        prev = s
    return tot
def least(n, eta, target):
    lo, hi = mp.mpf(0), mp.mpf(eta)*mp.mpf('1.05')
    if count(n, eta, hi) < target: return None
    for _ in range(20):
        mid = (lo+hi)/2
        if count(n, eta, mid) >= target: hi = mid
        else: lo = mid
    return float(hi)
out = []
import sys
ETAS = [float(a) for a in sys.argv[1].split(',')] if len(sys.argv) > 1 else [0.1, 0.2]
NS = [int(a) for a in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1, 2, 4]
for eta in ETAS:
    w0 = mp.mpf(eta)*0.995
    psi = mp.arg(Gfun(m + 1j*w0)); psi = float(abs(psi - mp.pi*mp.nint(psi/mp.pi)))
    for n in NS:
        all_real = least(n, eta, 2*n); some = least(n, eta, 2)
        pred = eta*math.sin((n-1)*math.pi/(2*n) + psi/n) if n % 2 == 0 else eta*(1 - 0.5*(psi/n)**2)
        out.append(dict(n=n, eta=eta, all_real=all_real, first_two=some, pred=pred, psi=psi))
        if all_real is None or some is None: print(f"eta={eta} n={n}: grid count failed", flush=True); continue
        print(f"eta={eta} n={n}: all 2n real from w={all_real:.5f} (w/eta={all_real/eta:.4f}), first real pair at {some:.5f};"
              f" prediction {pred:.5f}", flush=True)
json.dump(out, open(sys.argv[3] if len(sys.argv) > 3 else "multiplicity.json", "w"), indent=0)
