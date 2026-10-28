# ht_dhsharp.py: the angle theorem numbers of Nick Polson's residue note (R:440 to 535): the Davenport Heilbronn zero and Taylor coefficients
# (R:484), tab:dh (R:492 to 493, smallest minors for xi and DH, orders 2 to 6 and 10, 1800 random minors), theta = 0.0072 and k <= 435 (R:500),
# the prediction near order 4000 at position 44 (R:527), tab:sharp (R:518 to 520, first failing orders when the m-th zero of the first 30 is
# replaced by a pair at angle theta), and the lemma on three term sequences (the bare thresholds 4, 6, 10).
# Conventions. F(s) = xi(1/2 + sqrt s)/xi(1/2), c_n its Taylor coefficients. DH: f_DH(s) = sum a_n n^(-s), a_n period 5 with values
# (1, kappa, -kappa, -1, 0), kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1), so f_DH(s) = 5^(-s) [zeta(s, 1/5) + kappa zeta(s, 2/5) - kappa zeta(s, 3/5)
# - zeta(s, 4/5)]; xi_DH(s) = (5/pi)^(s/2) Gamma((s+1)/2) f_DH(s), xi_DH(s) = xi_DH(1 - s). F_DH(s) = xi_DH(1/2 + sqrt s)/xi_DH(1/2), entire in s.
# Coefficients by the Cauchy integral on |s| = r with N nodes (real coefficients, so half the nodes), r = 1000, N = 1024, 220 digits (the chapter's
# r = 45 at 140 digits cannot resolve c_104 ~ 1e-300 r^104 ~ 1e-164 against |F| ~ 3: a large radius is needed, F(1000) is only ~ 1e7); the check
# is the same computation for xi against the Polya integral (trapezoid step 1/200 on [-4, 4]) and the second radius r = 700 for both (the maximum relative
# difference of c_n for n <= 104 is logged and written).
# Minors. Contiguous Toeplitz minors det[c_(n+i-j)]_(i,j<k) for n from k-1 to 104-k-off (off = 0 and 1 are both reported, the chapter does not give
# the range), normalised by c_n^k (diag) as in the first table of R. Random minors: random row and column sets of size k in 0..104 (index
# differences outside 0..109 give zero), normalised by the product of the largest entry of each row.
# tab:sharp: polynomial prod_(n<=30, n != m)(1 + s/lambda_n) (1 + 2 cos(theta) s/lambda_m + s^2/lambda_m^2), lambda_n = gamma_n^2, zeros from
# acb.zeta_zeros at 800 bits; first order k <= 40 with a contiguous Toeplitz minor det[c_(n+i-j)] (c_j = 0 outside 0..degree) whose ball is
# certainly negative, n from 0 to degree + 1; the n of the first failure is reported.
# Options: --quick, --out FILE. Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, time, argparse, random
from mpmath import mp, mpf, mpc, zeta, gamma, pi, exp, log, sqrt, findroot, matrix, det, nstr, atan, floor
ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/ht_dhsharp.json')
a = ap.parse_args()
LOG = open(a.out.replace('.json', '.log'), 'w')
def lg(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
R = {}; t0 = time.time()
DPS = 80 if a.quick else 220
mp.dps = DPS
NODES = 128 if a.quick else 1024
RAD = 100 if a.quick else 1000
NC = 40 if a.quick else 109
kap = (sqrt(10 - 2 * sqrt(5)) - 2) / (sqrt(5) - 1)
def fDH(s):
    return mpf(5) ** (-s) * (zeta(s, mpf(1) / 5) + kap * zeta(s, mpf(2) / 5) - kap * zeta(s, mpf(3) / 5) - zeta(s, mpf(4) / 5))
def xiDH(s):
    return (mpf(5) / pi) ** (s / 2) * gamma((s + 1) / 2) * fDH(s)
def xi(s):
    return s * (s - 1) / 2 * pi ** (-s / 2) * gamma(s / 2) * zeta(s)
mp.dps = 60
s0 = mpc('0.3', '2.1')
lg('functional equation: |xi_DH(s) - xi_DH(1-s)|/|xi_DH(s)| at s = 0.3 + 2.1 i (60 digits):', nstr(abs(xiDH(s0) - xiDH(1 - s0)) / abs(xiDH(s0)), 3), ' kappa =', nstr(kap, 15))
rho = findroot(xiDH, mpc('0.8085', '85.699'), tol=mpf(10) ** -50)
lg('N22 DH zero: %s ; |xi_DH| there %s' % (nstr(rho, 16), nstr(abs(xiDH(rho)), 3)))
R['dh_zero'] = nstr(rho, 16)
dl = rho.real - mpf(1) / 2; gm = rho.imag
th = 2 * atan(abs(dl) / gm)
th2 = atan(2 * abs(dl) * gm / (gm ** 2 - dl ** 2))
lg('N24 theta = 2 arctan(delta/gamma) = %s ; arctan(2 delta gamma/(gamma^2-delta^2)) = %s ; pi/theta - 1 = %s ; (1/gamma = %s)' % (nstr(th, 8), nstr(th2, 8), nstr(pi / th - 1, 7), nstr(1 / gm, 6)))
R['theta'] = nstr(th, 8); R['k_bound'] = nstr(pi / th - 1, 7)
Nt = gm / (2 * pi) * log(5 * gm / (2 * pi * mp.e)) + mpf(7) / 8
lg('N27 pi/theta * sqrt(44) = %s ; 1.2 and 1.7 times: %s, %s ; smooth count (T/2pi) log(5T/(2 pi e)) + 7/8 at 85.7 = %s (no 1/4 phase term; rough)' % (nstr(pi / th * sqrt(44), 6), nstr(1.2 * pi / th * sqrt(44), 5), nstr(1.7 * pi / th * sqrt(44), 5), nstr(Nt, 5)))

mp.dps = DPS
def cauchy(fun, r, N, nmax):
    # c_n = (1/N r^n) sum_k fun(r e^{i th_k}) e^{-i n th_k}, th_k = 2 pi k/N; fun real on the real axis, so use k = 0..N/2 and conjugates
    vals = []
    for k in range(N // 2 + 1):
        th_ = 2 * pi * k / N
        s = r * mpc(mp.cos(th_), mp.sin(th_))
        vals.append(fun(s))
    cs = []
    for n in range(nmax + 1):
        tot = vals[0].real + (vals[N // 2].real * (-1) ** n)
        for k in range(1, N // 2):
            th_ = 2 * pi * k / N
            tot += 2 * (vals[k] * mpc(mp.cos(n * th_), -mp.sin(n * th_))).real
        cs.append(tot / (N * mpf(r) ** n))
    return cs
xih = xi(mpf(1) / 2); xidh = xiDH(mpf(1) / 2)
Fxi = lambda s: xi(mpf(1) / 2 + sqrt(s)) / xih
Fdh = lambda s: xiDH(mpf(1) / 2 + sqrt(s)) / xidh
cx = cauchy(Fxi, RAD, NODES, NC); cd = cauchy(Fdh, RAD, NODES, NC)
cx2 = cauchy(Fxi, mpf(700 if not a.quick else 70), NODES, NC); cd2 = cauchy(Fdh, mpf(700 if not a.quick else 70), NODES, NC)
lg('Cauchy coefficients done, time %.0f s' % (time.time() - t0))
def polya(dps, hinv, nmax):
    mp.dps = dps; h = mpf(1) / hinv; J = 4 * hinv; Phi = []
    for j in range(J + 1):
        u = j * h; e2 = exp(2 * u); s = mpf(0); m = 1
        while True:
            t = (2 * pi ** 2 * m ** 4 * exp(9 * u / 2) - 3 * pi * m ** 2 * exp(5 * u / 2)) * exp(-pi * m * m * e2)
            s += t
            if m > 2 and abs(t) < mpf(10) ** (-dps - 30): break
            m += 1
        Phi.append(2 * s)
    mom = lambda k: h * ((Phi[0] if k == 0 else mpf(0)) + 2 * sum(Phi[j] * (j * h) ** k for j in range(1, J + 1)))
    m0 = mom(0)
    return [mom(2 * n) / (m0 * mp.factorial(2 * n)) for n in range(nmax + 1)]
cp = polya(DPS, 200, NC)
mp.dps = DPS
errx = max(abs(cx[n] / cp[n] - 1) for n in range(NC + 1)); errx2 = max(abs(cx[n] / cx2[n] - 1) for n in range(NC + 1)); errd2 = max(abs(cd[n] / cd2[n] - 1) for n in range(NC + 1))
lg('N22 validation: xi Cauchy (r=%d) vs Polya, max rel diff of c_n, n <= %d: %s ; r=%d vs r=700: xi %s, DH %s' % (RAD, NC, nstr(errx, 3), RAD, nstr(errx2, 3), nstr(errd2, 3)))
R['val_xi_polya'] = nstr(errx, 3); R['val_radii'] = [nstr(errx2, 3), nstr(errd2, 3)]
lg('c_1: xi %s, DH %s ; c_2: xi %s, DH %s' % (nstr(cp[1], 10), nstr(cd[1], 10), nstr(cp[2], 10), nstr(cd[2], 10)))
N = NC - 5 if True else NC
NM = 104 if not a.quick else 35
def minors_tab(cs, ks, off):
    out = {}
    for k in ks:
        best, nb = None, None
        for n in range(k - 1, NM - k + 1 - off):
            M = matrix(k, k)
            for i in range(k):
                for j in range(k):
                    M[i, j] = cs[n + i - j]
            v = det(M) / cs[n] ** k
            if best is None or v < best: best, nb = v, n
        out[k] = (best, nb)
    return out
ks = [2, 3, 4, 5, 6, 10]
for off in (0, 1):
    mx = minors_tab(cx, ks, off); md = minors_tab(cd, ks, off)
    lg('N23 tab:dh with n from k-1 to %d - k - off... off=%d (n up to %d - k): xi: %s' % (NM, off, NM - off, ['k=%d %s (n=%d)' % (k, nstr(v[0], 3), v[1]) for k, v in mx.items()]))
    lg('                                                              DH: %s' % (['k=%d %s (n=%d)' % (k, nstr(v[0], 3), v[1]) for k, v in md.items()]))
    R['tabdh_off%d' % off] = {'xi': {k: nstr(v[0], 3) for k, v in mx.items()}, 'dh': {k: nstr(v[0], 3) for k, v in md.items()}}
random.seed(7)
for name, cs in (('xi', cx), ('DH', cd)):
    negs = 0; worst = mpf(1); tot = 0
    for k in (3, 4, 6):
        for _ in range(60 if a.quick else 600):
            rows = sorted(random.sample(range(0, NM + 1), k)); cols = sorted(random.sample(range(0, NM + 1), k))
            M = matrix(k, k)
            for i in range(k):
                for j in range(k):
                    d_ = rows[i] - cols[j]
                    M[i, j] = cs[d_] if 0 <= d_ <= NC else mpf(0)
            rm = mpf(1)
            for i in range(k): rm *= max(abs(M[i, j]) for j in range(k))
            v = det(M) / rm if rm != 0 else mpf(0)
            tot += 1; worst = min(worst, v)
            if v < -mpf(10) ** -80: negs += 1
    lg('N23 random minors (%s): %d of orders 3,4,6; negative below -1e-80: %d; smallest normalised %s' % (name, tot, negs, nstr(worst, 3)))
    R['random_' + name] = [tot, negs, nstr(worst, 3)]
lg('DH part done, time %.0f s' % (time.time() - t0))

# ---------------------------------------------------------------- tab:sharp (flint)
from flint import arb, acb, ctx, arb_mat
ctx.prec = 400 if a.quick else 800
Z = acb.zeta_zeros(1, 30)
lam = [(z.imag) ** 2 for z in Z]
pi_a = arb.pi()
def polymul(p, q):
    out = [arb(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out
def first_fail(coefs, kmax):
    deg = len(coefs) - 1
    for k in range(2, kmax + 1):
        for n in range(0, deg + 2):
            M = arb_mat(k, k)
            for i in range(k):
                for j in range(k):
                    idx = n + i - j
                    M[i, j] = coefs[idx] if 0 <= idx <= deg else arb(0)
            d = M.det()
            if d < 0:
                return k, n
    return None, None
ms = [1, 3] if a.quick else [1, 3, 6, 11, 16]
tabs = {}
for tname, q in (('pi/4', 4), ('pi/6', 6), ('pi/10', 10)):
    if a.quick and q != 4: continue
    cth = (pi_a / q).cos()
    bare, nb = first_fail([arb(1), 2 * cth, arb(1)], 40)
    row = []
    for m in ms:
        P = [arb(1)]
        for n in range(30):
            if n == m - 1: continue
            P = polymul(P, [arb(1), 1 / lam[n]])
        P = polymul(P, [arb(1), 2 * cth / lam[m - 1], 1 / lam[m - 1] ** 2])
        k, n = first_fail(P, 40)
        row.append((m, k, n))
    lg('N25 tab:sharp theta=%s: bare threshold %s (first failing order of the quadratic factor alone); (m, first failing order k, index n): %s' % (tname, bare, row))
    tabs[tname] = {'bare': bare, 'row': row}
R['sharp'] = tabs
lg('time %.0f s' % (time.time() - t0))
json.dump(R, open(a.out, 'w'), indent=1, default=str)
