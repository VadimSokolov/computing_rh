# Plant the pair where the full field of the other zeros (including the mirror pair) vanishes; the theory predicts a gap
# of higher order than eta^3. Compare with the midpoint placement. Floating point, mpmath.
import mpmath as mp, json
mp.mp.dps = 20
G = [mp.im(mp.zetazero(j)) for j in range(1, 20)]
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def dlogxi(s): return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
def threshold(F, c, eta):
    A = lambda z, w: F(z+1j*w) + F(z-1j*w)
    z = mp.mpc(c); w = mp.mpf(0); h = mp.mpf(eta)/20
    while h > mp.mpf(eta)*1e-8 and w < 1.05*eta:
        wn = w + h
        try: zn = mp.findroot(lambda u: A(u, wn), z)
        except Exception: h /= 2; continue
        if mp.im(zn) > 1e-11 and abs(zn - z) < 0.2: w, z = wn, zn
        else: h /= 2
    return w
out = []
for k in [1, 3, 5]:
    gk, gk1 = G[k-1], G[k]
    lo = (G[k-2] if k >= 2 else mp.mpf(1)) + mp.mpf('0.3'); hi = G[k+1] - mp.mpf('0.3')
    for eta in [0.1, 0.2, 0.3, 0.4]:
        def field(x, eta=eta):
            c0 = x + 1j*eta   # pair moves with x
            return mp.re(1j*dlogxi(mp.mpf(1)/2+1j*x) - (2*x/(x*x-gk*gk) + 2*x/(x*x-gk1*gk1)) + 1/(x+c0) + 1/(x+mp.conj(c0)))
        xs = [lo + (hi-lo)*i/400 for i in range(401)]; vals = [field(x) for x in xs]
        roots = [mp.findroot(field, (xs[i], xs[i+1]), solver='anderson') for i in range(400) if vals[i]*vals[i+1] < 0]
        for label, m in [("midpoint", (gk+gk1)/2)] + [("field zero", r) for r in roots]:
            c = m + 1j*eta; cb = m - 1j*eta
            F = lambda z: xi(mp.mpf(1)/2+1j*z)*(z*z-c*c)*(z*z-cb*cb)/((z*z-gk*gk)*(z*z-gk1*gk1))
            w = threshold(F, c, eta)
            out.append(dict(k=k, place=label, m=float(m), eta=eta, gap=float(eta-w), field=float(field(m))))
            print(f"k={k} {label:10s} m={float(m):8.4f} field={float(field(m)):+.2e} eta={eta}: gap {float(eta-w):.3e}", flush=True)
json.dump(out, open("extremum2.json", "w"), indent=0)
