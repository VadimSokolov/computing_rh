# Phase criterion: real zeros of A_w(x) = 2 Re F(x+iw) are where Theta_w(x) = arg F(x+iw) crosses pi/2 + k pi.
import mpmath as mp, numpy as np, json
mp.mp.dps = 15
g1, g2 = mp.im(mp.zetazero(1)), mp.im(mp.zetazero(2)); m = (g1+g2)/2; eta = mp.mpf('0.3')
def Xi(z):
    s = mp.mpf(1)/2 + 1j*z
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
c = m + 1j*eta; cb = m - 1j*eta
def F(z): return Xi(z)*(z*z-c*c)*(z*z-cb*cb)/((z*z-g1*g1)*(z*z-g2*g2))
xs = np.linspace(float(m)-2.5, float(m)+2.5, 1201)
out = {}
for w in [0.2, 0.28, 0.2951, 0.31]:
    th = np.unwrap(np.array([float(mp.arg(F(mp.mpf(x)+1j*w))) for x in xs]))
    out[str(w)] = th.tolist()
json.dump(dict(x=xs.tolist(), th=out, m=float(m)), open('phase.json', 'w'))
