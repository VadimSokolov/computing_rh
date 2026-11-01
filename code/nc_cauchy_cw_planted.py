# Strip Lee-Yang for spins of Xi type, PSD couplings. Two copies of the planted spin (strip 0.3) in a Curie Weiss coupling
# with matrix a*[[1,1],[1,1]] (positive semidefinite): Z(h) = E_Z[ phi(h + sqrt(a) Z)^2 ], Z standard normal,
# phi(s) = F(-i s)/F(0) the moment generating function of the planted law. Track the zero of Z near the planted zero
# s = 0.3 + i m and record its real part as a grows: the theorem says it never exceeds 0.3. Floating point, mpmath.
import mpmath as mp, numpy as np, json
mp.mp.dps = 18
g1, g2 = mp.im(mp.zetazero(1)), mp.im(mp.zetazero(2)); m = (g1 + g2) / 2; eta = mp.mpf('0.3')
c = m + 1j * eta; cb = m - 1j * eta


def Xi(z):
    s = mp.mpf(1) / 2 + 1j * z
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def F(z):
    return Xi(z) * (z * z - c * c) * (z * z - cb * cb) / ((z * z - g1 * g1) * (z * z - g2 * g2))


F0 = F(mp.mpf(0))
phi = lambda s: F(-1j * s) / F0
xs, ws = np.polynomial.hermite_e.hermegauss(40)          # nodes for E f(Z), Z ~ N(0,1)
ws = ws / ws.sum()


def Zfun(h, a):
    if a == 0:
        return phi(h) ** 2
    return sum(mp.mpf(w) * phi(h + mp.sqrt(a) * mp.mpf(x)) ** 2 for x, w in zip(xs, ws))


out = []
z = mp.mpc(eta, m)       # zero of phi at s = eta + i m (double for the product)
for a in [0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.07, 0.1]:
    try:
        z = mp.findroot(lambda h: Zfun(h, mp.mpf(a)), z + mp.mpc(0, 1e-3))
    except Exception as e:
        print("a", a, "failed", e); continue
    out.append((a, float(mp.re(z)), float(mp.im(z))))
    print(f"a={a}: zero at Re h = {float(mp.re(z)):.5f}, Im h = {float(mp.im(z)):.4f}", flush=True)
json.dump(out, open("cw_planted.json", "w"))
