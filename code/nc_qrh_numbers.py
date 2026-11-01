# nc_qrh_numbers.py: numbers printed in ch/nc_qrh.tex (the October 2026 half plane claims and the Thorin programme).
# Exact rational arithmetic (fractions) for the exponent bookkeeping; mpmath at 40 digits for xi'/xi and the root of
# Liu's cubic. Floating point values are computed, not certified.
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 40

out = {}
# 1. Exponent bookkeeping of OAI26b, (2.5): |A_1|^2 << D^{1+eps} H/Y + D^2/Y^2, Y = H^{1/6}, H = D^h.
def exps(h):
    return (F(1) + F(5, 6) * h) / 2, 1 - h / 6
for h in [F(1), F(11, 10), F(6, 7)]:
    a, b = exps(h)
    out['h=%s' % h] = (str(a), str(b), str(max(a, b)))
# 2. Part II of OAI26a: C(s) = s - 11/16, C(7/8) = 3/16; Part I C(s) = s - 2/3, C(11/12) = 1/4.
out['C_II(7/8)'] = str(F(7, 8) - F(11, 16)); out['C_I(11/12)'] = str(F(11, 12) - F(2, 3))
# Liu's comparison B(l) = 11/12 - l/4, l = 1/6 gives 7/8; l < 1 forces B > 2/3; B = 1/2 needs l = 5/3.
out['B(1/6)'] = str(F(11, 12) - F(1, 6) / 4); out['B(1)'] = str(F(11, 12) - F(1, 4)); out['l(1/2)'] = str(4 * (F(11, 12) - F(1, 2)))
# 3. Liu's constants.
f = lambda x: 1188 * x**3 - 2160 * x**2 + 171 * x + 190
xi = mp.findroot(f, mp.mpf('0.38858'))
Bnew = (18 * xi + 7) / (18 * xi + 9)
xi2 = mp.mpf(20) / 33 + mp.sqrt(1391) / 33 * mp.cos(mp.acos(mp.mpf(28470) / mp.mpf(1391)**1.5) / 3 - 2 * mp.pi / 3)
Balg = (1507 - 2 * mp.sqrt(921)) / 1653
out['xi_new'] = mp.nstr(xi, 15); out['xi_new closed form'] = mp.nstr(xi2, 15)
out['f(0.38858), f(0.38859)'] = (mp.nstr(f(mp.mpf('0.38858')), 6), mp.nstr(f(mp.mpf('0.38859')), 6))
out['B_new'] = mp.nstr(Bnew, 18); out['B_alg'] = mp.nstr(Balg, 15); out['B_r'] = mp.nstr(mp.mpf(34999) / 40000, 10)
out['B_alg - B_new'] = mp.nstr(Balg - Bnew, 4)
# 4. Conditional consequences: de Bruijn bound 2 (theta - 1/2)^2, and b_theta = xi'/xi(theta), u_theta(0) = b/pi.
def logxi(s):
    return mp.log(s * (s - 1) / 2) - s / 2 * mp.log(mp.pi) + mp.loggamma(s / 2) + mp.log(mp.zeta(s))
thetas = {'11/12': mp.mpf(11) / 12, '7/8': mp.mpf(7) / 8, 'B_alg': Balg, 'B_new': Bnew, '1': mp.mpf(1)}
for k, t in thetas.items():
    b = mp.diff(logxi, t) if t != 1 else None
    if t == 1:
        b = 1 + mp.euler / 2 - mp.log(4 * mp.pi) / 2
    out['theta=%s' % k] = dict(dBN=mp.nstr(2 * (t - mp.mpf(1) / 2)**2, 8), b=mp.nstr(b, 10), u0=mp.nstr(b / mp.pi, 10),
                               ratio_b_over_width=mp.nstr(b / (t - mp.mpf(1) / 2), 8))
# 5. Condition number of the Cauchy deconvolution of width w at grid step h: exp(pi w / h) (Section cx:sec:problem).
for w, name in [(mp.mpf(1) / 8, '1 to 7/8'), (mp.mpf(3) / 8, '7/8 to 1/2'), (mp.mpf(1) / 2, '1 to 1/2')]:
    out['log cond %s, h=0.05' % name] = mp.nstr(mp.pi * w / mp.mpf('0.05'), 6)
for k, v in out.items():
    print(k, ':', v)
# 6. Sign of the canonical measure between prime powers: density e^{-(alpha-1)x} (e^{-3x}/(1-e^{-2x}) - 1);
# the bracket vanishes where e^{-3x} + e^{-2x} = 1.
x0 = mp.findroot(lambda x: mp.e**(-3 * x) + mp.e**(-2 * x) - 1, 0.3)
print('x0 :', mp.nstr(x0, 12), ' log 2 :', mp.nstr(mp.log(2), 12))
for x in [mp.mpf('0.5'), mp.log(2) - mp.mpf('1e-9')]:
    print('bracket at', mp.nstr(x, 8), ':', mp.nstr(mp.e**(-3 * x) / (1 - mp.e**(-2 * x)) - 1, 10))
# 7. Comparison threshold of Lemma lem:dbeta:comp(b) under Q_theta: alpha >= 1/2 + sqrt2 (theta - 1/2).
print('alpha threshold at 7/8 :', mp.nstr(mp.mpf(1) / 2 + mp.sqrt(2) * mp.mpf(3) / 8, 8))
