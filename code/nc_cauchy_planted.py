# Problem 2 on Xi with a planted off line pair: F(z) = Xi(z) (z^2-c^2)(z^2-cb^2)/((z^2-g1^2)(z^2-g2^2)),
# c = m + i eta, m = (g1+g2)/2. A_w(z) = F(z+iw)+F(z-iw) (the cosine tilt). Count nonreal zeros of A_w in
# [1,45.3]x[-2.2,2.2] by the argument principle minus real sign changes. Floating point (mpmath, 15 digits).
import mpmath as mp, numpy as np, json, sys
mp.mp.dps = 15
g1, g2 = float(mp.im(mp.zetazero(1))), float(mp.im(mp.zetazero(2)))
m = (g1+g2)/2
def Xi(z): 
    s = mp.mpf(1)/2 + 1j*z
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def F(z, eta):
    c = m + 1j*eta; cb = m - 1j*eta
    return Xi(z)*(z*z-c*c)*(z*z-cb*cb)/((z*z-g1*g1)*(z*z-g2*g2))
def A(z, eta, w): return F(z+1j*w, eta) + F(z-1j*w, eta)
def count(eta, w, X0=1.0, X1=45.3, Y=2.2, n=1500):
    t = np.linspace(0, 1, n, endpoint=False)
    pts = list(X0+(X1-X0)*t - 1j*Y) + list(X1+1j*(-Y+2*Y*t[::5])) + list(X1-(X1-X0)*t + 1j*Y) + list(X0+1j*(Y-2*Y*t[::5]))
    v = np.array([complex(A(mp.mpc(p), eta, w)) for p in pts])
    a = np.angle(v); d = np.diff(np.concatenate([a, a[:1]])); d = (d+np.pi) % (2*np.pi) - np.pi
    if np.max(np.abs(d)) > 2.5: print("warning: coarse step", eta, w, np.max(np.abs(d)), file=sys.stderr)
    tot = int(round(d.sum()/(2*np.pi)))
    xs = np.linspace(X0, X1, 2200); f = np.array([float(mp.re(A(mp.mpf(x), eta, w))) for x in xs])
    real = int(np.sum(np.sign(f[1:]) != np.sign(f[:-1])))
    return tot, real
eta = float(sys.argv[1]); out = []
for w in [float(x) for x in sys.argv[2].split(",")]:
    tot, real = count(eta, w); out.append((w, tot, real, tot-real)); print(eta, w, tot, real, tot-real, flush=True)
json.dump(out, open(f"planted_{eta}_{sys.argv[3]}.json", "w"))
