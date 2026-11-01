"""Checks of the half line Thorin density m for ch/nc_riemanncont.tex (mpmath, 30 digits, computed)."""
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
g = mp.im(mp.zetazero(1))
k = lambda th: th**2/(pi*g*(g*g+th**2))
for u in [mp.mpf(3), mp.mpf(150), mp.mpf(500)]:
    r = mp.sqrt(u)
    I = -2*mp.quad(lambda th: mp.log(abs(1-u/th**2))*k(th), [0, r/2, r, 2*r, g, mp.inf])
    print("u=", u, " -2 int log|1-u/th^2| m =", mp.nstr(I, 20), "  log(1+u/g^2) =", mp.nstr(mp.log(1+u/g**2), 20))
for s in [mp.mpf(1), mp.mpf(10)]:
    I = mp.quad(lambda th: mp.log(1+s*s/th**2)*k(th), [0, s, g, mp.inf])
    P = (1+s/g)*mp.exp(-s/g)
    print("s=", s, " int log(1+s^2/th^2) m =", mp.nstr(I, 20), "  -log P(s) =", mp.nstr(-mp.log(P), 20))
# off line pair: m(sqrt w)/w has poles at w=-a, a=(gamma - i delta)^2
c = mp.mpc(100, -0.45)
f = lambda w: 2*mp.re(1/(c*(c*c+w)))  # real for real w
for w in [mp.mpf(1), mp.mpf(10**4)]:
    print("pair contribution to pi m(sqrt w)/w at w=", w, ":", mp.nstr(f(w), 12))
print("pole location -a =", mp.nstr(-c*c, 12))
