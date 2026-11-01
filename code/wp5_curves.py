# wp5_curves.py: Work package 5 of the run book (rh_compute.tex Section sec:curves): the Weil ladder for zeta functions of curves over F_q.
# Written from the mathematics (the collaborator's code and nc_route_genus2.py were not read).
#
# Curve: y^2 + h(x) y = f(x), deg f = 2g+1, deg h <= g (h = 0 in odd characteristic), one point at infinity.  Zeta Z(T) = P(T)/((1-T)(1-qT)),
# P(T) = prod (1 - alpha_j T) = 1 + a_1 T + ... + q^g T^{2g}, N_m = q^m + 1 - S_m, S_m = sum alpha_j^m.   zeta_C(s) = Z(q^{-s}), ell = log q.
# Zeros: s = log(alpha_j)/ell + 2 pi i k/ell; r = -i (s - 1/2) is real for all of them iff RH.  Poles r = 2 pi k/ell +- i/2.
# Weil form for even real f supported in [-L, L], g(u) = int f(v) f(v-u) dv, ghat(r) = F(r) F(-r):
#   W(f) = sum_zeros ghat(r)  =  sum_k [ghat(2 pi k/ell + i/2) + ghat(2 pi k/ell - i/2)]  +  (2g-2) ell g(0)  -  2 ell sum_{d,k>=1} d P_d q^{-kd/2} g(k d ell)
# (closed points of degree d, P_d of them: N_m = sum_{d|m} d P_d; the middle term is the conductor q^{(g-1)s}).  Poisson summation turns the pole sum into a finite sum:
#   W(f) = ell [ 2g g(0) + 2 sum_{m>=1} ( (q^{m/2} + q^{-m/2}) - N_m q^{-m/2} ) g(m ell) ] = ell [ 2g g(0) + 2 sum_m s_m g(m ell) ],  s_m = S_m q^{-m/2}.
# Ladder: modes phi_0 = (2L)^{-1/2}, phi_n = L^{-1/2} cos(pi n u / L), 2L = tau ell (lambda^2 = q^tau, closed points of norm up to lambda^2), matrix W_jk with g -> g_jk.
#
# Subcommands (all from data/ on Hopper; numpy and mpmath always, python-flint for the certified brackets):
#   counts   --out F          point counts of curves of genus 1,2,3 over every prime power q <= 49 (own GF(p^n) arithmetic), numerators, checks
#   validate --out F          explicit formula against the zeros for the bump function, matrix entries against quadrature
#   toeplitz --out F          the Toeplitz test of the route chapter (y^2 = x^5 + x + 1 and fake numerators)
#   ladder   --out F --spec S the Weil ladder for a list of curves / fakes (see main)
import sys, json, math, random, time, argparse
import numpy as np
import mpmath as mp


# ---------------------------------------------------------------- finite fields
def factorint(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return sorted(set(out))


class GF:
    """F_{p^n} as F_p[x]/(x^n + c_{n-1} x^{n-1} + ... + c_0), x primitive; elements are ints (base p digits); exp/log tables."""
    _cache = {}

    @classmethod
    def get(cls, p, n):
        if (p, n) not in cls._cache:
            cls._cache[(p, n)] = GF(p, n)
        return cls._cache[(p, n)]

    def __init__(s, p, n):
        s.p, s.n, s.Q = p, n, p ** n
        s.pw = np.array([p ** i for i in range(n)], dtype=np.int64)
        rng = random.Random(1000 * p + n)
        pr = factorint(s.Q - 1)
        while True:
            c = [rng.randrange(p) for _ in range(n)]
            if c[0] == 0:
                continue
            s.red = [(-x) % p for x in c]  # x^n = sum red_i x^i
            if s._powx(s.Q - 1) != s._one() and s.Q > 2:
                continue
            if s.Q > 2 and any(s._powx((s.Q - 1) // r) == s._one() for r in pr):
                continue
            break
        # tables by doubling
        D = np.zeros((s.Q - 1, n), dtype=np.int64)
        D[0, 0] = 1
        L = 1
        while L < s.Q - 1:
            g = np.array(s._powx(L), dtype=np.int64)
            m = min(L, s.Q - 1 - L)
            D[L:L + m] = s._mulconst(D[:m], g)
            L += m
        s.exp = D @ s.pw
        s.log = -np.ones(s.Q, dtype=np.int64)
        s.log[s.exp] = np.arange(s.Q - 1)

    def _one(s):
        return [1] + [0] * (s.n - 1)

    def _mulx(s, A):
        top = A[..., -1]
        new = np.zeros_like(A)
        new[..., 1:] = A[..., :-1]
        return (new + top[..., None] * np.array(s.red, dtype=np.int64)) % s.p

    def _mulconst(s, A, g):
        res = np.zeros_like(A)
        cur = A.copy()
        for i in range(s.n):
            res = (res + g[i] * cur) % s.p
            cur = s._mulx(cur)
        return res

    def _polymul(s, a, b):
        A = np.array(a, dtype=np.int64)[None, :]
        return list(map(int, s._mulconst(A, np.array(b, dtype=np.int64))[0]))

    def _powx(s, e):
        res = s._one()
        base = [0, 1] + [0] * (s.n - 2) if s.n > 1 else list(s.red)
        while e:
            if e & 1:
                res = s._polymul(res, base)
            base = s._polymul(base, base)
            e >>= 1
        return res

    def mul(s, a, b):
        a, b = np.asarray(a), np.asarray(b)
        a, b = np.broadcast_arrays(a, b)
        res = np.zeros(a.shape, dtype=np.int64)
        m = (a != 0) & (b != 0)
        res[m] = s.exp[(s.log[a[m]] + s.log[b[m]]) % (s.Q - 1)]
        return res

    def add(s, a, b):
        a, b = np.asarray(a), np.asarray(b)
        if s.p == 2:
            return a ^ b
        res = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        pa = 1
        for _ in range(s.n):
            res = res + (((a // pa) % s.p + (b // pa) % s.p) % s.p) * pa
            pa *= s.p
        return res

    def chi(s, a):
        a = np.asarray(a)
        r = np.zeros(a.shape, dtype=np.int64)
        m = a != 0
        r[m] = 1 - 2 * (s.log[a[m]] & 1)
        return r

    def trace2(s, a):  # trace to F_2 (p = 2)
        t, tot = np.asarray(a), np.asarray(a)
        for _ in range(s.n - 1):
            t = s.mul(t, t)
            tot = tot ^ t
        return tot

    def horner(s, coefs, x):
        r = np.full(x.shape, coefs[-1], dtype=np.int64)
        for c in coefs[-2::-1]:
            r = s.add(s.mul(r, x), np.int64(c))
        return r


def embed_map(p, k, m):
    """map from F_q = GF(p,k) (ints) to GF(p, k m): sends the abstract element sum d_i x^i to sum d_i theta^i, theta a root of the defining polynomial of GF(p,k)."""
    F, G = GF.get(p, k), GF.get(p, k * m)
    if k == 1:
        return lambda c: int(c)
    sub = np.concatenate([[0], G.exp[np.arange(0, G.Q - 1, (G.Q - 1) // (F.Q - 1))]]) if F.Q > 1 else np.array([0])
    # evaluate x^k + sum_i c_i x^i (c = -red) at every element of the subfield
    val = np.zeros(sub.shape, dtype=np.int64)
    val = G.add(val, G.horner([(-r) % p for r in F.red] + [1], sub))
    theta = int(sub[np.where(val == 0)[0][0]])
    pw = [1]
    for _ in range(k - 1):
        pw.append(int(G.mul(np.int64(pw[-1]), np.int64(theta))))

    def f(c):
        d = [(int(c) // p ** i) % p for i in range(k)]
        r = np.int64(0)
        for di, ti in zip(d, pw):
            r = G.add(r, G.mul(np.int64(di), np.int64(ti)))
        return int(r)
    return f


def count_points(p, k, m, fco, hco):
    """#C(F_{q^m}), q = p^k, y^2 + h y = f, coefficients abstract ints in GF(p,k) (low degree first); point at infinity included (deg f odd)."""
    G = GF.get(p, k * m)
    em = embed_map(p, k, m)
    f = [em(c) for c in fco]
    x = np.arange(G.Q, dtype=np.int64)
    fx = G.horner(f, x)
    if p != 2:
        return 1 + int(G.Q + G.chi(fx).sum())
    h = [em(c) for c in hco]
    hx = G.horner(h, x)
    z0 = hx == 0
    cnt = int(z0.sum())
    nz = ~z0
    lf = np.zeros(int(nz.sum()), dtype=np.int64)
    fn, hn = fx[nz], hx[nz]
    z = np.zeros(fn.shape, dtype=np.int64)
    mm = fn != 0
    z[mm] = G.exp[(G.log[fn[mm]] - 2 * G.log[hn[mm]]) % (G.Q - 1)]
    tr = G.trace2(z)
    cnt += 2 * int((tr == 0).sum())
    return 1 + cnt


# ---------------------------------------------------------------- numerators
def numerator_from_counts(q, g, Ns):
    """a_1..a_{2g} (ints) of P(T) from N_1..N_g by Newton's identities and the functional equation a_{2g-k} = q^{g-k} a_k."""
    S = [None] + [q ** m + 1 - Ns[m - 1] for m in range(1, g + 1)]
    a = [1]
    for k in range(1, g + 1):
        t = -S[k] - sum(a[i] * S[k - i] for i in range(1, k))
        assert t % k == 0
        a.append(t // k)
    full = a + [0] * g
    for k in range(g):
        full[2 * g - k] = q ** (g - k) * full[k]
    return full  # list of length 2g+1, a_0 = 1


def power_sums(q, g, a, M):
    a = [ctx.mpf(x) for x in a]
    S = [None]
    for k in range(1, M + 1):
        t = -(k * a[k] if k <= 2 * g else 0) - sum(a[i] * S[k - i] for i in range(1, min(k - 1, 2 * g) + 1))
        S.append(t)
    return S


def mobius(n):
    r, d, m = 1, 2, n
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:
                return 0
            r = -r
        d += 1
    return -r if m > 1 else r


def closed_points(q, S, M):
    """P_d, d = 1..M, from N_m = q^m + 1 - S_m; asserts integrality."""
    N = [None] + [q ** m + 1 - S[m] for m in range(1, M + 1)]
    P = [None]
    for d in range(1, M + 1):
        t = sum(mobius(d // e) * N[e] for e in range(1, d + 1) if d % e == 0)
        assert t % d == 0 and t >= 0, (d, t)
        P.append(t // d)
    return N, P


def roots_of_numerator(q, g, a, dps=40):
    mp.mp.dps = dps
    # P(T) = prod(1 - alpha T): alpha are the roots of x^{2g} + a_1 x^{2g-1} + ... + a_{2g}
    return mp.polyroots([mp.mpf(c) for c in a], maxsteps=500, extraprec=200)


def rh_holds(q, g, a):
    al = roots_of_numerator(q, g, a)
    dev = max(abs(abs(x) - mp.sqrt(q)) / mp.sqrt(q) for x in al)
    return float(dev), al


# ---------------------------------------------------------------- curves
PRIME_POWERS = [(2, 1), (3, 1), (2, 2), (5, 1), (7, 1), (2, 3), (3, 2), (11, 1), (13, 1), (2, 4), (17, 1), (19, 1), (23, 1), (5, 2), (3, 3), (29, 1), (31, 1),
                (2, 5), (37, 1), (41, 1), (43, 1), (47, 1), (7, 2)]
QCAP = 7_000_000


def random_curve(p, k, g, rng):
    q = p ** k
    while True:
        if p != 2:
            f = [rng.randrange(q) for _ in range(2 * g + 1)] + [1]
            if f[-2] == 0 and False:
                continue
            return f, [0]
        f = [rng.randrange(q) for _ in range(2 * g + 1)] + [1]
        h = [rng.randrange(q) for _ in range(g)] + [rng.randrange(1, q)]
        if h[g] == 0:
            continue
        return f, h


def analyse_curve(p, k, g, fco, hco):
    q = p ** k
    Ns = [count_points(p, k, m, fco, hco) for m in range(1, g + 1)]
    a = numerator_from_counts(q, g, Ns)
    rec = dict(q=q, p=p, k=k, g=g, f=fco, h=hco, N=Ns, a=a)
    chk = None
    if q ** (g + 1) <= QCAP:
        Ng1 = count_points(p, k, g + 1, fco, hco)
        pred = q ** (g + 1) + 1 - power_sums(q, g, a, g + 1)[g + 1]
        chk = (Ng1 == pred)
        rec['N_g1_counted'], rec['N_g1_predicted'] = Ng1, pred
    rec['check_Ng1'] = chk
    dev, _ = rh_holds(q, g, a)
    rec['rh_dev'] = dev
    rec['hasse_weil_ok'] = bool(dev < 1e-12)
    return rec


def cmd_counts(args):
    out = []
    t0 = time.time()
    for g in (1, 2, 3):
        for (p, k) in PRIME_POWERS:
            q = p ** k
            if q ** (g + 1) > QCAP:
                pass
            rng = random.Random(7 * q + 31 * g)
            tries = 0
            while True:
                tries += 1
                fco, hco = random_curve(p, k, g, rng)
                rec = analyse_curve(p, k, g, fco, hco)
                if rec['check_Ng1'] is not False and rec['hasse_weil_ok']:
                    break
                rec.update(rejected=True)
                out.append(rec)
                print('rejected', g, q, rec['N'], rec['check_Ng1'], rec['rh_dev'], flush=True)
                if tries > 12:
                    break
            rec['tries'] = tries
            out.append(rec)
            print(g, q, rec['N'], rec['a'], rec['check_Ng1'], '%.1e' % rec['rh_dev'], '%.0fs' % (time.time() - t0), flush=True)
    # the named curve y^2 = x^5 + x + 1 over F_p, p = 3..47 (also the bad primes: discriminant 5^5 + 4^4 = 3381 = 3 7^2 23)
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        rec = analyse_curve(p, 1, 2, [1, 1, 0, 0, 0, 1], [0])
        rec['named'] = 'y^2=x^5+x+1'
        out.append(rec)
        print('named', p, rec['N'], rec['a'], rec['check_Ng1'], '%.1e' % rec['rh_dev'], flush=True)
    json.dump(out, open(args.out, 'w'), indent=1)


# ---------------------------------------------------------------- the form (generic over mp.mp, floating, and mp.iv, rigorous intervals)
def gjk_fun(ctx, L, N):
    """G(j, k, u) = int phi_j(v) phi_k(v-u) dv, even modes on [-L, L], 0 <= u <= 2L (derived with the product formula; checked against quadrature in validate)."""
    a = [ctx.pi * j / L for j in range(N + 1)]
    c = [1 / ctx.sqrt(2 * L)] + [1 / ctx.sqrt(L)] * N

    def G(j, k, u):
        sj, sk = ctx.sin(a[j] * u), ctx.sin(a[k] * u)
        if j == k == 0:
            return c[0] ** 2 * (2 * L - u)
        if j == k:
            return c[k] ** 2 / 2 * (-sj / a[k] + (2 * L - u) * ctx.cos(a[k] * u))
        sg = -1 if (j + k) % 2 else 1
        return c[j] * c[k] / 2 * sg * (-(sj + sk) / (a[j] + a[k]) + (sk - sj) / (a[j] - a[k]))
    return G


def s_coeffs(ctx, q, g, a, M):
    """s_m = S_m q^{-m/2}, m = 1..M, from the numerator coefficients a (list a_0..a_{2g}, ints or mp numbers) by Newton's identities."""
    a = [ctx.mpf(x) for x in a]
    S = [None]
    for k in range(1, M + 1):
        t = -(k * a[k] if k <= 2 * g else 0) - sum(a[i] * S[k - i] for i in range(1, min(k - 1, 2 * g) + 1))
        S.append(t)
    rq = ctx.sqrt(ctx.mpf(q))
    return [None] + [ctx.mpf(S[m]) / rq ** m if not isinstance(S[m], int) else ctx.mpf(S[m]) / rq ** m for m in range(1, M + 1)]


def build_matrix(ctx, q, g, a, tau, N):
    ell = ctx.log(ctx.mpf(q))
    L = ell * tau / 2
    M = int(math.floor(tau - 1e-12))
    s = s_coeffs(ctx, q, g, a, max(M, 1))
    G = gjk_fun(ctx, L, N)
    W = [[None] * (N + 1) for _ in range(N + 1)]
    us = [ell * m for m in range(1, M + 1)]
    for j in range(N + 1):
        for k in range(j, N + 1):
            t = ctx.mpf(0)
            for m in range(1, M + 1):
                t = t + s[m] * G(j, k, us[m - 1])
            w = ell * (2 * t + (2 * g if j == k else 0))
            W[j][k] = W[k][j] = w
    return W


def ldl_pd(A, shift, ctx):
    """True if A - shift I has an LDL^T factorisation with all pivot intervals positive (interval arithmetic: then every symmetric matrix inside A - shift I is positive definite)."""
    n = len(A)
    Lm = [[ctx.mpf(0)] * n for _ in range(n)]
    D = [None] * n
    for i in range(n):
        d = A[i][i] - shift
        for k in range(i):
            d = d - Lm[i][k] ** 2 * D[k]
        if not (d.a > 0):
            return False
        D[i] = d
        for j in range(i + 1, n):
            t = A[j][i]
            for k in range(i):
                t = t - Lm[j][k] * Lm[i][k] * D[k]
            Lm[j][i] = t / d
    return True


def rayleigh(A, v, ctx):
    n = len(A)
    vi = [ctx.mpf(x) for x in v]
    num = ctx.mpf(0)
    for i in range(n):
        row = ctx.mpf(0)
        for j in range(n):
            row = row + A[i][j] * vi[j]
        num = num + vi[i] * row
    den = ctx.mpf(0)
    for x in vi:
        den = den + x * x
    return num / den


def rung(q, g, a, tau, N, dps, certify=True):
    """one rung: computed least eigenvalue (mp, dps digits), and a certified bracket from mp.iv at the same precision."""
    mp.mp.dps = dps
    A = build_matrix(mp.mp, q, g, a, tau, N)
    E, Q = mp.eigsy(mp.matrix(A))
    ev = [E[i] for i in range(N + 1)]
    i0 = min(range(N + 1), key=lambda i: ev[i])
    eps = ev[i0]
    out = dict(q=q, g=g, tau=tau, N=N, dps=dps, eps=mp.nstr(eps, 12), log10eps=(float(mp.log10(abs(eps))) if eps != 0 else None), neg=bool(eps < 0))
    out['eps2'] = mp.nstr(sorted(ev)[1], 8)
    if certify:
        mp.iv.dps = dps
        Ai = build_matrix(mp.iv, q, g, a, tau, N)
        v = [Q[j, i0] for j in range(N + 1)]
        rq = rayleigh(Ai, v, mp.iv)
        out['hi'] = mp.nstr(rq.b, 8)  # certified: least eigenvalue of the exact matrix <= hi
        lo = None
        if eps > 0:
            for f in (0.999, 0.99, 0.9, 0.5, 0.0):
                if ldl_pd(Ai, mp.iv.mpf(eps * f), mp.iv):
                    lo = f
                    break
            out['lo_frac'] = lo  # certified: least eigenvalue >= lo_frac * eps (None: not certified at this precision)
        out['cert_neg'] = bool(rq.b < 0)
    return out


# ---------------------------------------------------------------- validation
def bump_funcs(L):
    f = lambda u: mp.cos(mp.pi * u / (2 * L)) ** 2 if abs(u) <= L else mp.mpf(0)
    c = mp.pi ** 2 / L ** 2

    def F(r):  # int f e^{iru}, f even
        r = mp.mpmathify(r)
        if abs(r) < mp.mpf(10) ** -20:
            r = r + mp.mpf(10) ** -20
        return mp.sin(r * L) * (1 / r - r / (r * r - c))
    return f, F


def closed_point_counts(q, S, M):
    N, P = closed_points(q, S, M)
    N2 = [None] + [sum(d * P[d] for d in range(1, m + 1) if m % d == 0) for m in range(1, M + 1)]
    assert N == N2
    return N, P


def validate_one(q, g, a, tau, K=3000):
    mp.mp.dps = 30
    ell = mp.log(q)
    L = ell * tau / 2
    M = int(math.floor(tau))
    f, F = bump_funcs(L)
    S = power_sums(q, g, a, M)
    N, P = closed_point_counts(q, S, M)
    gq = lambda u: mp.quad(lambda v: f(v) * f(v - u), [u - L, L]) if u < 2 * L else mp.mpf(0)
    gm = [gq(ell * m) for m in range(0, M + 1)]
    # zero side
    al = roots_of_numerator(q, g, a, 30)
    Zs = mp.mpf(0)
    ks = range(-K, K + 1)
    for x in al:
        r0 = -1j * (mp.log(x) / ell - mp.mpf(1) / 2)
        r0 = complex(r0)
        rr = r0 + 2 * np.pi * np.arange(-K, K + 1) / float(ell)
        Zs += mp.mpf(float(np.sum(_Fnp(rr, float(L)) ** 2).real))
    # explicit formula side: poles + conductor + primes (closed points)
    rk = 2 * np.pi * np.arange(-K, K + 1) / float(ell)
    poles = np.sum(_Fnp(rk + 0.5j, float(L)) ** 2 + _Fnp(rk - 0.5j, float(L)) ** 2).real
    prime = 2 * ell * sum(d * P[d] * mp.mpf(q) ** (-mp.mpf(k * d) / 2) * gm[k * d] for d in range(1, M + 1) for k in range(1, M // d + 1))
    E = mp.mpf(float(poles)) + (2 * g - 2) * ell * gm[0] - prime
    # Poisson finite form
    sm = s_coeffs(mp.mp, q, g, a, M)
    Pf = ell * (2 * g * gm[0] + 2 * sum(sm[m] * gm[m] for m in range(1, M + 1)))
    return dict(q=q, g=g, a=a, tau=tau, zeros=mp.nstr(Zs, 15), explicit=mp.nstr(E, 15), poisson=mp.nstr(Pf, 15), rel_ze=mp.nstr(abs(Zs - E) / abs(Zs), 3),
                rel_zp=mp.nstr(abs(Zs - Pf) / abs(Zs), 3))


def _Fnp(r, L):
    r = np.asarray(r, dtype=complex)
    c = np.pi ** 2 / L ** 2
    r = np.where(np.abs(r) < 1e-12, 1e-12, r)
    return np.sin(r * L) * (1 / r - r / (r * r - c))


def validate_matrix(q, g, a, tau, N=3):
    mp.mp.dps = 30
    ell = mp.log(q)
    L = ell * tau / 2
    M = int(math.floor(tau))
    sm = s_coeffs(mp.mp, q, g, a, M)
    W = build_matrix(mp.mp, q, g, a, tau, N)
    rng = random.Random(3)
    xi = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N + 1)]
    phi = lambda j, u: (1 / mp.sqrt(2 * L) if j == 0 else mp.cos(mp.pi * j * u / L) / mp.sqrt(L)) if abs(u) <= L else mp.mpf(0)
    fx = lambda u: sum(xi[j] * phi(j, u) for j in range(N + 1))
    gq = lambda u: mp.quad(lambda v: fx(v) * fx(v - u), [u - L, L]) if u < 2 * L else mp.mpf(0)
    val = ell * (2 * g * gq(0) + 2 * sum(sm[m] * gq(ell * m) for m in range(1, M + 1)))
    qf = sum(xi[i] * W[i][j] * xi[j] for i in range(N + 1) for j in range(N + 1))
    return dict(q=q, g=g, tau=tau, N=N, quad=mp.nstr(val, 15), matrix=mp.nstr(qf, 15), rel=mp.nstr(abs(val - qf) / abs(val), 3))


def cmd_validate(args):
    out = dict(forms=[], matrices=[])
    for (q, g, a) in [(13, 2, [1, 1, 4, 13, 169]), (9, 1, [1, 1, 9]), (49, 3, [1, -11, 67, -370, 3283, -26411, 117649])]:
        for tau in (2.7, 5.5, 9.3):
            r = validate_one(q, g, a, tau)
            out['forms'].append(r)
            print(r, flush=True)
            r = validate_matrix(q, g, a, tau)
            out['matrices'].append(r)
            print(r, flush=True)
    json.dump(out, open(args.out, 'w'), indent=1)


# ---------------------------------------------------------------- Toeplitz test of the route chapter
def toeplitz_min(a, q, g, sizes, dps=60):
    mp.mp.dps = dps
    S = power_sums(q, g, a, max(sizes))
    s = [mp.mpf(2 * g)] + [mp.mpf(S[m]) / mp.sqrt(q) ** m for m in range(1, max(sizes))]
    out = {}
    for n in sizes:
        T = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                T[i, j] = s[abs(i - j)]
        E = mp.eigsy(T, eigvals_only=True)
        ev = sorted(E[i] for i in range(n))
        out[n] = (mp.nstr(ev[0], 6), sum(1 for e in ev if abs(e) > mp.mpf(10) ** -(dps // 2)))
    return out


def cmd_toeplitz(args):
    res = {}
    sizes = [2, 3, 4, 5, 6, 8, 12, 20, 30, 40]
    named = []
    for p in (7, 11, 13, 17, 19, 23, 29, 31):
        rec = analyse_curve(p, 1, 2, [1, 1, 0, 0, 0, 1], [0])
        tm = toeplitz_min(rec['a'], p, 2, sizes)
        named.append(dict(p=p, a=rec['a'], N=rec['N'], check_Ng1=rec['check_Ng1'], rh_dev=rec['rh_dev'], toeplitz=tm))
        print(p, rec['a'], rec['check_Ng1'], '%.1e' % rec['rh_dev'], {n: tm[n][0] for n in (2, 3, 5, 40)}, 'rank at 40:', tm[40][1], flush=True)
    res['named'] = named
    # fakes at p = 13: same c1, c2 = 29 and -26
    a13 = [r for r in named if r['p'] == 13][0]['a']
    c1 = a13[1]
    fakes = []
    for c2 in (29, -26, a13[2] + 5, 40):
        a = [1, c1, c2, 13 * c1, 169]
        dev, al = rh_holds(13, 2, a)
        mods = sorted(float(abs(x) / mp.sqrt(13)) for x in al)
        tm = toeplitz_min(a, 13, 2, list(range(2, 41)))
        first = next((n for n in range(2, 41) if mp.mpf(tm[n][0]) < -1e-9), None)
        fakes.append(dict(c1=c1, c2=c2, moduli_over_sqrtp=[mods[0], mods[-1]], first_fail=first, min3=tm[3][0], min40=tm[40][0], rh_dev=dev))
        print('fake', c1, c2, mods[0], mods[-1], 'first failing size', first, tm[3][0], tm[40][0], flush=True)
    res['fakes'] = fakes
    json.dump(res, open(args.out, 'w'), indent=1)


# ---------------------------------------------------------------- the ladder
def fake_genus1(q, a1):
    return [1, a1, q]


def fake_genus2(q, c1, c2):
    return [1, c1, c2, q * c1, q * q]


def fake_genus3(q, b, c1, c2):
    A, B = [1, b, q], fake_genus2(q, c1, c2)
    out = [0] * 7
    for i, x in enumerate(A):
        for j, y in enumerate(B):
            out[i + j] += x * y
    return out


def worker(task):
    label, q, g, a, tau, N = task
    a = [int(x) if isinstance(x, (int,)) or (isinstance(x, str) and x.lstrip('-').isdigit()) else x for x in a]
    d0 = 50 + int(0.16 * N * tau * math.log(q))
    t0 = time.time()
    r1 = rung(q, g, a, tau, N, d0, certify=False)
    dps = d0
    for _ in range(4):
        dps += 40
        r2 = rung(q, g, a, tau, N, dps, certify=True)
        e1, e2 = mp.mpf(r1['eps']), mp.mpf(r2['eps'])
        agree = abs(e1 - e2) <= mp.mpf(1e-6) * abs(e2)
        r1 = r2
        if agree:
            break
    r2.update(label=label, a=[str(x) for x in a], stable=bool(agree), secs=round(time.time() - t0, 1))
    return r2


def cmd_ladder(args):
    from multiprocessing import Pool
    counts = json.load(open(args.counts))
    pick = lambda g, q: [r for r in counts if r['g'] == g and r['q'] == q and not r.get('rejected') and 'named' not in r][0]
    curves = []
    for g in (1, 2, 3):
        for q in (3, 5, 13, 49):
            r = pick(g, q)
            curves.append(('genuine g%d q%d' % (g, q), q, g, r['a']))
    named = [r for r in counts if r.get('named') and r['p'] == 13][0]
    curves.append(('genuine y2=x5+x+1 q13', 13, 2, named['a']))
    c1 = named['a'][1]
    for a1 in (8, 9, 12):
        curves.append(('fake g1 q13 a1=%d' % a1, 13, 1, fake_genus1(13, a1)))
    curves.append(('fake g1 q5 a1=5', 5, 1, fake_genus1(5, 5)))
    for c2 in (29, -26, 27, 26.5, 26.3, 26.26):
        a = fake_genus2(13, c1, c2) if isinstance(c2, int) else [1, c1, str(c2), str(13 * c1), 169]
        curves.append(('fake g2 q13 c1=%d c2=%s' % (c1, c2), 13, 2, a))
    b = pick(1, 13)['a'][1]
    curves.append(('fake g3 q13 (g1 factor b=%d, fake g2 c2=29)' % b, 13, 3, fake_genus3(13, b, c1, 29)))
    taus = [float(x) for x in args.taus.split(',')]
    Ns = [int(x) for x in args.Ns.split(',')]
    tasks = [(lab, q, g, [str(x) if not isinstance(x, int) else x for x in a], tau, N) for (lab, q, g, a) in curves for tau in taus for N in Ns]
    print(len(tasks), 'tasks', flush=True)
    res = []
    with Pool(args.procs) as pool:
        for i, r in enumerate(pool.imap_unordered(worker, tasks, chunksize=1)):
            res.append(r)
            if i % 50 == 0:
                print(i, r['label'], r['tau'], r['N'], r['eps'], r['secs'], flush=True)
    json.dump(res, open(args.out, 'w'), indent=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('--out', default='wp5_out.json')
    ap.add_argument('--counts', default='wp5_counts.json')
    ap.add_argument('--taus', default='2,3,4,6,8,12,16,24')
    ap.add_argument('--Ns', default='2,4,8,12,16,24')
    ap.add_argument('--procs', type=int, default=8)
    args = ap.parse_args()
    dict(counts=cmd_counts, validate=cmd_validate, toeplitz=cmd_toeplitz, ladder=cmd_ladder)[args.cmd](args)


if __name__ == '__main__':
    main()
