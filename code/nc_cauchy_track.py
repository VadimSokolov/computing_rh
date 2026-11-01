# Track the nonreal zero of A_w(z) = F(z+iw)+F(z-iw) born from the planted zero c = m + 0.3 i, as w grows.
import mpmath as mp, json
mp.mp.dps = 20
g1, g2 = mp.im(mp.zetazero(1)), mp.im(mp.zetazero(2)); m = (g1+g2)/2; eta = mp.mpf('0.3')
def Xi(z):
    s = mp.mpf(1)/2 + 1j*z
    return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
c = m + 1j*eta; cb = m - 1j*eta
def F(z): return Xi(z)*(z*z-c*c)*(z*z-cb*cb)/((z*z-g1*g1)*(z*z-g2*g2))
def A(z, w): return F(z+1j*w) + F(z-1j*w)
z = c; out = []
w = mp.mpf(0)
while w <= mp.mpf('0.2965'):
    z = mp.findroot(lambda u: A(u, w), z)
    out.append((float(w), float(mp.re(z)), float(mp.im(z))))
    w += mp.mpf('0.0025') if w < mp.mpf('0.28') else mp.mpf('0.0005')
print(out[0], out[len(out)//2], out[-1])
json.dump(out, open('track.json', 'w'))
