# Threshold of the slice by continuation: track the nonreal zero born from the planted c = m + i eta as w grows,
# and locate the w at which its imaginary part reaches 0 (square root law). Floating point, mpmath 15 digits.
import mpmath as mp, json, sys
mp.mp.dps = 15
G = [mp.im(mp.zetazero(j)) for j in range(1, 40)]


def Xi(z):
    s = mp.mpf(1) / 2 + 1j * z
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


res = []
for k in [int(x) for x in sys.argv[1].split(",")]:
    for eta in [0.1, 0.2, 0.3, 0.4]:
        gk, gk1 = G[k - 1], G[k]
        m = (gk + gk1) / 2
        c = m + 1j * eta
        cb = m - 1j * eta
        F = lambda z: Xi(z) * (z * z - c * c) * (z * z - cb * cb) / ((z * z - gk * gk) * (z * z - gk1 * gk1))
        A = lambda z, w: F(z + 1j * w) + F(z - 1j * w)
        z = mp.mpc(c)
        w = mp.mpf(0)
        h = mp.mpf(eta) / 20
        while h > mp.mpf(eta) * 1e-6 and w < 1.05 * eta:
            wn = w + h
            try:
                zn = mp.findroot(lambda u: A(u, wn), z)
            except Exception:
                h /= 2
                continue
            if mp.im(zn) > 1e-9 and abs(zn - z) < 0.2:
                w, z = wn, zn
            else:
                h /= 2
        res.append(dict(k=k, height=float(m), eta=eta, threshold=float(w), rel=float((eta - w) / eta)))
        print(f"k={k} m={float(m):.3f} eta={eta}: w*={float(w):.5f} (eta-w*)/eta={float((eta - w) / eta):+.4f}", flush=True)
json.dump(res, open(f"survey2_{sys.argv[1].replace(',', '_')}.json", "w"), indent=0)
