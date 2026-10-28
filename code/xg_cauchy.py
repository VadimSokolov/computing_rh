# code/xg_cauchy.py (October 2026, from misc/reviews/nick-xiggc/: review-A2.md, the rows thm:cauchy and tab:cauchy, item (d)9 and errors E3 and E15):
# the numbers of Theorem xg:thm:cauchy (Section sec:xg:cauchy, the end of Section 13.4). The program is written from the mathematics, not from the
# verify scripts of the authors' chapter. Notation of Section ch:dbn: Phi_dBN(u) = sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) exp(-pi n^2 e^{4u})
# (even in u by the functional equation of theta, so the series is valid for every real u), H_0(z) = int_0^oo Phi_dBN(u) cos(zu) du = Xi(z/2)/8, and for
# a real mu the tilted kernel H^(mu)(z) = int_0^oo e^{-mu u} Phi_dBN(u) cos(zu) du (mu > 0 is the cusped weight e^{-mu|u|} of the theorem, mu < 0 the weight
# e^{+|mu||u|}). Phi_dBN is negligible beyond u = 1.6 (the first term is exp(-pi e^{6.4}) = e^{-1888}).
# (1) The rule. [0, 1.6] is cut into 160 pieces of length 0.01 and each piece gets the 24 node Gauss Legendre rule of mpmath (40 digits). For
# |z| <= 400 the error of the rule is below e^{-100} (the integrand is analytic in |Im u| < pi/8 and the Bernstein ellipse of the piece reaches 0.1 into the plane,
# where e^{400 * 0.1} e^{-48 log 40} = e^{-137}). The rule is checked three ways: H_0(z) against Xi(z/2)/8 computed from mpmath's zeta at five heights
# (Xi(t) = xi(1/2 + it), xi(s) = s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)/2), H_0(2 gamma_j) = 0 at the six ordinates gamma_j <= 40 of data/zeros_8000.txt, and
# the evenness of the series.
# (2) The tail law. Integrating by parts twice, int_0^oo f(u) cos(zu) du = -f'(0)/z^2 + f'''(0)/z^4 - ... for f = e^{-mu u} Phi_dBN, with f'(0) = -mu Phi_dBN(0) and
# f'''(0) = -3 mu Phi_dBN''(0) - mu^3 Phi_dBN(0) (Phi_dBN is even and analytic), so z^2 H^(mu)(z) / (mu Phi_dBN(0)) = 1 + (3 |Phi_dBN''(0)| / Phi_dBN(0) - mu^2)/z^2 + O(z^-4).
# Phi_dBN(0) and Phi_dBN''(0) are summed from the series (the second derivative termwise, 32 c^4 - 224 c^3 + 330 c^2 - 75 c times e^{-c} with c = pi n^2, and again by
# mpmath's numerical derivative); the four rows of the chapter's table and the Laguerre quantity (H')^2 - H H'' at z = 33.5, mu = 1, are evaluated at 40 digits.
# (3) The scans, in double precision with the same nodes (the values are at least 1e-9 against rounding errors of 1e-14; the points that decide a claim are
# redone at 40 digits): mu = 1 on [0, 80], step 0.05 (sign changes of H^(1), sign of the Laguerre quantity), mu = 0.01 on [0, 400], step 0.05 (sign changes
# and the last one), H_0 on [0, 80] (six sign changes), and the weight e^{+u} (mu = -1) on [0, 80] (Laguerre quantity). The Cauchy representation
# H^(mu)(z) = int H_0(y) mu / (pi (mu^2 + (z-y)^2)) dy is checked in double precision with Gauss Legendre pieces of length 5 on [-150, 150], which also checks
# int H_0 = pi Phi_dBN(0).
# Run on Hopper from a folder that holds zeros_8000.txt: python3 xg_cauchy.py (about a minute, one core); writes xg_cauchy.json there. The argument quick
# coarsens the scans (step 0.5) for a test on the laptop.
import json
import os
import sys
import time

import mpmath as mp
import numpy as np
from mpmath.calculus.quadrature import GaussLegendre

QUICK = 'quick' in sys.argv[1:]
mp.mp.dps = 40
t0 = time.time()
PI = mp.pi
STEP = 0.5 if QUICK else 0.05


def say(*a):
    print(*a, flush=True)


def sci(x, d=4):
    """x as d significant digits in the form the book prints, 7.237\\times10^{-5}."""
    m, e = ('%.*e' % (d - 1, float(x))).split('e')
    e = int(e)
    return '%s\\times10^{%d}' % (m, e) if e else m


def phi(u):
    """Phi_dBN(u): the series of Section ch:dbn, summed until exp(-pi n^2 e^{4u}) < e^{-200}."""
    u = mp.mpf(u)
    e4, e5, e9 = mp.exp(4 * u), mp.exp(5 * u), mp.exp(9 * u)
    s, n = mp.mpf(0), 1
    while n == 1 or PI * n * n * e4 <= 200:
        c = PI * n * n
        s += (2 * c * c * e9 - 3 * c * e5) * mp.exp(-c * e4)
        n += 1
    return s


# ---------------------------------------------------------------- 0. the kernel at the origin
phi0 = phi(0)
d2_series = mp.mpf(0)
for n in range(1, 40):
    c = PI * n * n
    d2_series += mp.exp(-c) * (32 * c ** 4 - 224 * c ** 3 + 330 * c ** 2 - 75 * c)
d2_numeric = mp.diff(phi, 0, 2)
assert abs(d2_series - d2_numeric) < mp.mpf(10) ** -25, (d2_series, d2_numeric)
assert d2_series < 0
assert abs(phi(mp.mpf('0.1')) - phi(mp.mpf('-0.1'))) < mp.mpf(10) ** -30       # even
assert abs(phi(mp.mpf('0.37')) - phi(mp.mpf('-0.37'))) < mp.mpf(10) ** -30
C3 = 3 * abs(d2_series) / phi0
say('Phi_dBN(0) = %s ; Phi_dBN\'\'(0) = %s (series and numerical derivative agree) ; 3|Phi\'\'(0)|/Phi(0) = %s ; Phi_pol(0) = 2 Phi_dBN(0) = %s'
    % (mp.nstr(phi0, 15), mp.nstr(d2_series, 15), mp.nstr(C3, 12), mp.nstr(2 * phi0, 12)))

# ---------------------------------------------------------------- 1. the rule
rule = GaussLegendre(mp.mp).calc_nodes(4, mp.mp.prec)       # 24 nodes on [-1, 1]
assert len(rule) == 24 and abs(sum(w for x, w in rule) - 2) < mp.mpf(10) ** -35
PIECE = mp.mpf(1) / 100
U, W = [], []
for k in range(160):
    mid, half = k * PIECE + PIECE / 2, PIECE / 2
    for x, w in rule:
        u = mid + half * x
        U.append(u)
        W.append(half * w * phi(u))                           # weight times kernel
Uf = np.array([float(u) for u in U])
Wf = np.array([float(w) for w in W])
say('rule: %d nodes on [0, 1.6] (%.1f s)' % (len(U), time.time() - t0))
EXPS = {}


def Hmp(mu, z, derivs=False):
    """H^(mu)(z) (and H', H'') at 40 digits from the rule."""
    mu, z = mp.mpf(mu), mp.mpf(z)
    key = str(mu)
    if key not in EXPS:
        EXPS[key] = [wp * mp.exp(-mu * u) for u, wp in zip(U, W)]
    h0 = h1 = h2 = mp.mpf(0)
    for u, c in zip(U, EXPS[key]):
        co, si = mp.cos_sin(z * u)
        h0 += c * co
        if derivs:
            h1 -= c * u * si
            h2 -= c * u * u * co
    return (h0, h1, h2) if derivs else h0


def scan(mu, zs):
    """H, H', H'' for the array zs in double precision."""
    cv = Wf * np.exp(-mu * Uf)
    c1, c2 = cv * Uf, cv * Uf ** 2
    H, H1, H2 = (np.empty(len(zs)) for _ in range(3))
    for k in range(0, len(zs), 400):
        arg = zs[k:k + 400, None] * Uf[None, :]
        C, S = np.cos(arg), np.sin(arg)
        H[k:k + 400], H1[k:k + 400], H2[k:k + 400] = C @ cv, -(S @ c1), -(C @ c2)
    return H, H1, H2


def sign_changes(zs, H):
    idx = np.nonzero(H[:-1] * H[1:] < 0)[0]
    return [float(zs[i] - H[i] * (zs[i + 1] - zs[i]) / (H[i + 1] - H[i])) for i in idx]


def xi_cap(t):
    """The book's Xi(t) = xi(1/2 + it), real."""
    s = mp.mpf(1) / 2 + 1j * mp.mpf(t)
    return (s * (s - 1) / 2 * PI ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)).real


out = dict(script='xg_cauchy.py', dps=40, nodes=len(U), quick=QUICK,
           Phi_dBN_0=mp.nstr(phi0, 25), Phi_dBN_second_derivative_0=mp.nstr(d2_series, 25), three_abs_Phi2_over_Phi0=mp.nstr(C3, 20))

# ---------------------------------------------------------------- 2. the rule against Xi and the zeros
errs = []
for z in ['0', '10', '28.5', '50', '80']:
    errs.append(abs(Hmp(0, z) - xi_cap(mp.mpf(z) / 2) / 8))
say('H_0(z) from the rule against Xi(z/2)/8 from zeta at z = 0, 10, 28.5, 50, 80: largest difference %s' % mp.nstr(max(errs), 3))
assert max(errs) < mp.mpf(10) ** -30
path = next(p for p in ('zeros_8000.txt', 'data/zeros_8000.txt') if os.path.exists(p))
gam = [mp.mpf(line.split()[0]) for line in open(path) if line.strip()]
below = [g for g in gam if 2 * g <= 80]
say('ordinates gamma_j with 2 gamma_j <= 80: %d (%s) ; the next one: 2 gamma = %s'
    % (len(below), ', '.join(mp.nstr(g, 7) for g in below), mp.nstr(2 * gam[len(below)], 7)))
zerr = max(abs(Hmp(0, 2 * g)) for g in below)
say('largest |H_0(2 gamma_j)| over these ordinates: %s' % mp.nstr(zerr, 3))
assert len(below) == 6 and zerr < mp.mpf(10) ** -30
out.update(rule_vs_Xi_max_difference=mp.nstr(max(errs), 5), ordinates_with_2gamma_le_80=len(below), max_abs_H0_at_2gamma=mp.nstr(zerr, 5))

# ---------------------------------------------------------------- 3. the four rows of the chapter's table
chapter = {('1', 80): ('7.237\\times10^{-5}', '6.980\\times10^{-5}'), ('1', 400): ('2.796\\times10^{-6}', '2.792\\times10^{-6}'),
           ('3', 400): ('8.387\\times10^{-6}', '8.376\\times10^{-6}'), ('0.01', 400): ('2.796\\times10^{-8}', '2.792\\times10^{-8}')}
rows = []
say('rows of the table (mu, z, H^(mu)(z), mu Phi(0)/z^2), as the book would print them, and the next term of the expansion:')
for (mu, z), (hc, ac) in chapter.items():
    h = Hmp(mu, z)
    asym = mp.mpf(mu) * phi0 / mp.mpf(z) ** 2
    ratio_obs = h / asym - 1
    ratio_pred = (C3 - mp.mpf(mu) ** 2) / mp.mpf(z) ** 2
    row = '%s & %d & $%s$ & $%s$\\\\' % (mu, z, sci(h), sci(asym))
    same = sci(h) == hc and sci(asym) == ac
    say('  ' + row + '   ratio - 1 = %s against the next term %s (quotient %s) ; chapter: %s %s' % (
        mp.nstr(ratio_obs, 5), mp.nstr(ratio_pred, 5), mp.nstr(ratio_obs / ratio_pred, 5), (hc, ac), 'AGREES' if same else 'DIFFERS'))
    rows.append(dict(mu=mu, z=z, H=mp.nstr(h, 12), asymptotic=mp.nstr(asym, 12), ratio_minus_1=mp.nstr(ratio_obs, 8),
                     next_term=mp.nstr(ratio_pred, 8), latex=row, chapter=dict(H=hc, asymptotic=ac), agrees_with_chapter=bool(same)))
out['table'] = rows

# ---------------------------------------------------------------- 4. the Laguerre quantity at z = 33.5, mu = 1
h0, h1, h2 = Hmp(1, '33.5', True)
lag = h1 ** 2 - h0 * h2
say('mu = 1, z = 33.5: H = %s, H\' = %s, H\'\' = %s, (H\')^2 - H H\'\' = %s (reviewer: -2.9e-9) ; -2 a^2/z^6 with a = Phi(0): %s'
    % (mp.nstr(h0, 8), mp.nstr(h1, 8), mp.nstr(h2, 8), mp.nstr(lag, 8), mp.nstr(-2 * phi0 ** 2 / mp.mpf('33.5') ** 6, 5)))
out['laguerre_mu1_z33_5'] = dict(H=mp.nstr(h0, 12), H1=mp.nstr(h1, 12), H2=mp.nstr(h2, 12), laguerre=mp.nstr(lag, 12), reviewer_value=-2.9e-9)

# ---------------------------------------------------------------- 5. the scans
zsA = np.round(np.arange(0, int(round(80 / STEP)) + 1) * STEP, 10)
H, H1, H2 = scan(1.0, zsA)
L = H1 ** 2 - H * H2
neg = zsA[L < 0]
runs = []
for z in neg:
    if runs and abs(z - runs[-1][1] - STEP) < 1e-9:
        runs[-1][1] = float(z)
    else:
        runs.append([float(z), float(z)])
say('mu = 1 on [0, 80], step %g: sign changes of H^(1): %s ; min H = %.4e at z = %g ; Laguerre quantity negative on %d grid points, first at z = %g, in runs %s'
    % (STEP, sign_changes(zsA, H), H.min(), zsA[H.argmin()], len(neg), neg[0] if len(neg) else None, runs[:6]))
Hz0, _, _ = scan(0.0, zsA)
sc0 = sign_changes(zsA, Hz0)
say('mu = 0 (H_0) on [0, 80]: %d sign changes at %s ; 2 gamma_j: %s' % (len(sc0), [round(x, 3) for x in sc0], [mp.nstr(2 * g, 6) for g in below]))
assert len(sc0) == 6
claimA = bool(sign_changes(zsA, H) == [] and H.min() > 0)
say('CHAPTER CLAIM (no sign change of H^(1) on [0, 80], where H_0 has six): %s' % ('HOLDS' if claimA else 'FAILS'))
out['scan_mu1'] = dict(claim_holds=claimA, step=STEP, sign_changes=sign_changes(zsA, H), min_H=float(H.min()), argmin_z=float(zsA[H.argmin()]),
                       laguerre_negative_points=int(len(neg)), laguerre_negative_runs=runs, H0_sign_changes=sc0)

zsB = np.round(np.arange(0, int(round(400 / STEP)) + 1) * STEP, 10)
HB, _, _ = scan(0.01, zsB)
scB = sign_changes(zsB, HB)
i80 = int(round(80 / STEP))
at40 = Hmp('0.01', 40)
say('mu = 0.01 on [0, 400], step %g: %d sign changes, the last at z = %.3f (all at %s) ; H(40) = %s (mp) ; H > 0 on every grid point z >= 80: %s (min %.4e)'
    % (STEP, len(scB), scB[-1], [round(x, 2) for x in scB], mp.nstr(at40, 6), bool((HB[i80:] > 0).all()), HB[i80:].min()))
claimB = bool(at40 < 0 and (HB[i80:] > 0).all())
say('CHAPTER CLAIM (H^(0.01) negative at 40 and positive from 80 on, up to 400 on the grid): %s' % ('HOLDS' if claimB else 'FAILS'))
out['scan_mu0_01'] = dict(claim_holds=claimB, step=STEP, sign_changes=scB, H_at_40=mp.nstr(at40, 12), positive_from_80_to_400=bool((HB[i80:] > 0).all()), min_H_on_80_400=float(HB[i80:].min()))

roots = [mp.findroot(lambda x: Hmp('0.01', x), (z0 - STEP, z0 + STEP), solver='anderson', tol=1e-30, verify=False) for z0 in scB]
say('real zeros of H^(0.01) in [0, 400], refined at 40 digits: %s (residuals %s)' % ([mp.nstr(r, 10) for r in roots], [mp.nstr(abs(Hmp('0.01', r)), 2) for r in roots]))
out['scan_mu0_01']['refined_zeros'] = [mp.nstr(r, 15) for r in roots]

HP, HP1, HP2 = scan(-1.0, zsA)
LP = HP1 ** 2 - HP * HP2
negP = zsA[LP < 0]
say('weight e^{+u} (mu = -1) on [0, 80]: sign changes of H %s, Laguerre quantity negative on %d grid points, first at z = %s'
    % (sign_changes(zsA, HP), len(negP), negP[0] if len(negP) else None))
zP = float(negP[0]) if len(negP) else None
if zP is not None:
    a, b, c = Hmp(-1, repr(zP), True)
    say('  redone at 40 digits at z = %g: (H\')^2 - H H\'\' = %s' % (zP, mp.nstr(b ** 2 - a * c, 8)))
    out['scan_plus_weight'] = dict(first_negative_z=zP, laguerre_at_first_negative=mp.nstr(b ** 2 - a * c, 12), sign_changes_H=sign_changes(zsA, HP))

# ---------------------------------------------------------------- 6. double precision against 40 digits, and the Cauchy representation
worst, worst_d = 0.0, 0.0
for mu, zlist in [(1, [5, 20, 33.5, 60, 79.95]), (0.01, [10, 40, 120, 399.95])]:
    for z in (zlist[:2] if QUICK else zlist):
        a, b, c = Hmp(mu, repr(z), True)
        f = scan(float(mu), np.array([float(z)]))
        worst = max(worst, float(abs(f[0][0] / a - 1)))
        if z <= 80:
            worst_d = max(worst_d, float(abs(f[1][0] / b - 1)), float(abs(f[2][0] / c - 1)))
say('double precision against 40 digits (H at 9 points, H\', H\'\' at 5): largest relative differences %.2e and %.2e' % (worst, worst_d))
assert worst < 1e-5 and worst_d < 1e-5
xg, wg = np.polynomial.legendre.leggauss(24)
edges = np.arange(-150.0, 150.0 + 2.5, 5.0)
yk = np.concatenate([(a + b) / 2 + (b - a) / 2 * xg for a, b in zip(edges[:-1], edges[1:])])
wk = np.concatenate([(b - a) / 2 * wg for a, b in zip(edges[:-1], edges[1:])])
H0y = scan(0.0, np.abs(yk))[0]
mass = float((wk * H0y).sum())
say('int H_0(y) dy = %.12f against pi Phi_dBN(0) = %.12f' % (mass, float(PI * phi0)))
cau = []
for mu, z in [(1, 80), (1, 400), (3, 400), (0.01, 400)]:
    v = float((wk * H0y * mu / (np.pi * (mu * mu + (z - yk) ** 2))).sum())
    d = float(Hmp(mu, z))
    cau.append(dict(mu=mu, z=z, cauchy_average=v, direct=d, relative_difference=abs(v / d - 1)))
    say('E H_0(z + mu C) at (mu, z) = (%g, %g): %.10e against the direct integral %.10e (relative difference %.1e)' % (mu, z, v, d, abs(v / d - 1)))
assert abs(mass / float(PI * phi0) - 1) < 1e-9 and all(r['relative_difference'] < 1e-6 for r in cau)
out['cauchy_representation'] = cau
out['int_H0'] = dict(numerical=mass, pi_Phi0=float(PI * phi0))
say('done in %.1f s' % (time.time() - t0))
out['seconds'] = round(time.time() - t0, 1)
with open('xg_cauchy.json', 'w') as f:
    json.dump(out, f, indent=1)
