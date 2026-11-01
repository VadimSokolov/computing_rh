# wl_weil.py: Work package 1 of the run book (the Weil ladder, rh_compute.tex Section sec:weil): the least eigenvalue eps_N(lambda) of the Weil form on the
# first N+1 even modes of [1/lambda, lambda], X = lambda^2, L = log(lambda) = log(X)/2. Written from the mathematics of the note (the collaborator's
# weil_ladder.py and ladder_report.py were not read, run or copied; his pilot/*.json were read only as numbers to compare against).
#
# The form (u = log x, f supported in [-L, L], g(u) = int f(v) f(v-u) dv, F(r) = int f e^{iru}):
#   W(f) = 2 F(i/2)^2 - g(0) (log pi + gamma_E + log(1 - e^{-4L})) + int_0^{2L} 2 (g(0) e^{-2u} - e^{-u/2} g(u)) / (1 - e^{-2u}) du - 2 sum_{n <= X} Lambda(n) n^{-1/2} g(log n).
# Modes phi_0 = (2L)^{-1/2}, phi_n = L^{-1/2} cos(a_n u), a_n = pi n / L.  c_0 = (2L)^{-1/2}, c_n = L^{-1/2}.  s_j(u) = sin(a_j u).
# Cross correlation g_jk(u) = int phi_j(v) phi_k(v-u) dv over [u-L, L] (derived here by the product formula, only s_j(u) enters, the endpoints are at multiples of pi):
#   j != k:  g_jk(u) = (c_j c_k / 2) (-1)^(j+k) [ -(s_j + s_k)/(a_j + a_k) + (s_k - s_j)/(a_j - a_k) ]
#   j = k >= 1:  g_kk(u) = (c_k^2 / 2) [ -s_k / a_k + (2L - u) cos(a_k u) ],   g_00(u) = c_0^2 (2L - u),   g_jk(0) = delta_jk.
# Because g_jk is linear in the s_j (off the diagonal), the whole off diagonal of the matrix needs only the N functionals Z_j = int_0^{2L} e^{-u/2} s_j/(1 - e^{-2u}) du + sum_n w_n s_j(log n),
# w_n = Lambda(n)/sqrt(n) (the primes up to X), and the diagonal needs N+1 further integrals D_k (the singular part cancels inside D_k) and the sums Q_k = sum_n w_n (2L - log n) cos(a_k log n).
# F_n(r) = c_n ( sin((r-a_n)L)/(r-a_n) + sin((r+a_n)L)/(r+a_n) ),  F_n(i/2) = c_n (-1)^n sinh(L/2)/(a_n^2 + 1/4).
#
# Rigorous quadrature. The 2N+1 integrands are analytic in u except for the poles of 1/(1 - e^{-2u}) at u = i pi m, m != 0 (u = 0 is removable: the numerator of D_k vanishes there
# because g_kk(0) = 1, and s_j(0) = 0). With u = L (1 + x), x in [-1, 1], and E_rho the Bernstein ellipse x = (rho e^{it} + rho^{-1} e^{-it})/2 with rho below the pole value rho_pole = |z + sqrt(z^2 - 1)|,
# z = -1 + i pi/L, the maximum modulus principle gives |f| <= M inside E_rho, M = max over the boundary, enclosed by evaluating the integrand on K complex balls covering the boundary.
# Chebyshev coefficients satisfy |a_k| <= 2 M rho^{-k}; n point Gauss-Legendre is exact up to degree 2n-1 and |int T_k - sum w_i T_k(x_i)| <= 2 + 2 = 4 for the rest, hence
#   |I - G_n| <= L * 8 M rho^{-2n} * rho/(rho - 1)     (all quantities balls; n is raised until this is below 10^-D; the bound is added to every integral as a ball radius).
# The prime sums are finite sums in ball arithmetic. So the matrix is a matrix of balls containing the exact Weil matrix.
# Positivity: an LDL^T factorisation in ball arithmetic; if all pivot balls are positive, every symmetric matrix in the box is positive definite (in particular the true one).
# Eigenvalue bracket: eps_up = Rayleigh quotient (ball) of the inverse iteration vector (exact binary rationals, upper bound for the true eps), eps_lo = largest t tried with A - t I certified positive definite.
#
# Usage (from data/ on Hopper, scripts and zeros_1700.txt beside):
#   python3 wl_weil.py validate zeros_1700.txt --out wl_validate.json
#   python3 wl_weil.py float --tasks 13:8,13:12 --procs 8 --out wl_float.json [--zeros zeros_1700.txt]
#   python3 wl_weil.py cert  --tasks 13:8,13:12 --procs 8 --out wl_cert.json
import argparse, json, math, random, sys, time
from multiprocessing import Pool
from flint import arb, acb, ctx


def tof(x):
    return float(x.mid().str(17, radius=False))


def bits_for(D):
    return int(D * 3.3219281) + 40


def prime_powers(X):
    """list of (n, p): n = p^m <= X"""
    sv = bytearray([1]) * (X + 1)
    out = []
    for p in range(2, X + 1):
        if sv[p]:
            for q in range(p * p, X + 1, p):
                sv[q] = 0
            q = p
            while q <= X:
                out.append((q, p))
                q *= p
    return out


class Geo:
    def __init__(s, X, N):
        s.X, s.N = X, N
        s.L = arb(X).log() / 2
        s.a = [arb.pi() * j / s.L for j in range(N + 1)]
        s.c = [1 / (2 * s.L).sqrt()] + [1 / s.L.sqrt()] * N
        s.K0 = arb.pi().log() + arb.const_euler() + (1 - (-4 * s.L).exp()).log()
        s.cfloat_L = tof(s.L)


def fam(g, u):
    """the 2N+1 integrands at u (arb or acb): Sv[j] (j = 1..N) = e^{-u/2} s_j/(1-e^{-2u}), Dv[k] (k = 0..N) = 2 (e^{-2u} - e^{-u/2} g_kk(u))/(1-e^{-2u})"""
    L, N, a, c = g.L, g.N, g.a, g.c
    e1 = (-u / 2).exp()
    e2 = (-2 * u).exp()
    R = 1 / (1 - e2)
    Sv = [None] * (N + 1)
    Dv = [None] * (N + 1)
    t = 2 * L - u
    Dv[0] = 2 * (e2 - e1 * t / (2 * L)) * R
    for k in range(1, N + 1):
        s, co = (a[k] * u).sin_cos()
        Sv[k] = e1 * s * R
        gk = c[k] ** 2 / 2 * (-s / a[k] + t * co)
        Dv[k] = 2 * (e2 - e1 * gk) * R
    return Sv, Dv


def boundary_M(g, rho, K):
    """upper bound (arb) for the max modulus of all integrands on the ellipse boundary, by K balls over t in [0, pi] (symmetry in conjugation)"""
    rho = arb(rho)
    A_ = (rho + 1 / rho) / 2
    B_ = (rho - 1 / rho) / 2
    best = arb(0)
    for i in range(K):
        th = arb((i + 0.5) * math.pi / K, math.pi / (2 * K))
        x = acb(A_ * th.cos(), B_ * th.sin())
        u = g.L * (1 + x)
        Sv, Dv = fam(g, u)
        for v in Sv[1:] + Dv:
            m = abs(v).upper()
            if not m.is_finite():
                raise ValueError('boundary ball not finite')
            if m > best:
                best = m
    return best


def choose_rho(g, D, K):
    """pick rho (below the pole) and the node count n minimising n"""
    Lf = g.cfloat_L
    z = complex(-1, math.pi / Lf)
    r1 = abs(z + (z * z - 1) ** 0.5)
    r2 = abs(z - (z * z - 1) ** 0.5)
    rp = max(r1, r2)
    best = None
    old = ctx.prec
    ctx.prec = 120
    for f in (0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9):
        rho = 1 + f * (rp - 1)
        M = boundary_M(g, rho, K)
        logM = tof(M.log())
        n = int(((math.log(8 * Lf * rho / (rho - 1)) + logM) + D * math.log(10)) / (2 * math.log(rho))) + 2
        if best is None or n < best[0]:
            best = (n, rho, M)
    ctx.prec = old
    return best, rp


def arch_integrals(g, D, bits):
    N = g.N
    K = 64 + 3 * N
    (n, rho, M), rp = choose_rho(g, D, K)
    ctx.prec = bits
    err = 8 * M * g.L * arb(rho) / (arb(rho) - 1) * arb(rho) ** (-2 * n)
    assert err < arb(10) ** (-D), 'quadrature bound'
    S = [arb(0)] * (N + 1)
    Dd = [arb(0)] * (N + 1)
    for k in range(n):
        x, w = arb.legendre_p_root(n, k, weight=True)
        u = g.L * (1 + x)
        w = w * g.L
        Sv, Dv = fam(g, u)
        for j in range(1, N + 1):
            S[j] += w * Sv[j]
        for j in range(N + 1):
            Dd[j] += w * Dv[j]
    z = err * arb(0, 1)
    S = [x + z for x in S]
    Dd = [x + z for x in Dd]
    info = dict(nodes=n, rho=rho, rho_pole=rp, K=K, M_log10=tof(M.log() / arb(10).log()), quad_err_log10=tof(err.log() / arb(10).log()))
    return S, Dd, info


def prime_sums(g):
    N = g.N
    P = [arb(0)] * (N + 1)
    Q = [arb(0)] * (N + 1)
    for n, p in prime_powers(g.X):
        u = arb(n).log()
        w = arb(p).log() / arb(n).sqrt()
        t = 2 * g.L - u
        Q[0] += w * t
        for k in range(1, N + 1):
            s, co = (g.a[k] * u).sin_cos()
            P[k] += w * s
            Q[k] += w * t * co
    return P, Q


def build(X, N, D):
    """matrix of balls (list of lists of arb, symmetric) and info"""
    bits = bits_for(D)
    ctx.prec = bits
    g = Geo(X, N)
    S, Dd, info = arch_integrals(g, D, bits)
    ctx.prec = bits
    P, Q = prime_sums(g)
    L, a, c = g.L, g.a, g.c
    sh = (L / 2).sinh()
    F = [c[j] * (-1) ** j * sh / (a[j] ** 2 + arb(1) / 4) for j in range(N + 1)]
    Z = [arb(0)] + [S[j] + P[j] for j in range(1, N + 1)]
    A = [[None] * (N + 1) for _ in range(N + 1)]
    for j in range(N + 1):
        for k in range(j, N + 1):
            if j != k:
                sg = (-1) ** (j + k)
                v = 2 * F[j] * F[k] - c[j] * c[k] * sg * (-(Z[j] + Z[k]) / (a[j] + a[k]) + (Z[k] - Z[j]) / (a[j] - a[k]))
            else:
                pr = 2 * c[0] ** 2 * Q[0] if j == 0 else c[j] ** 2 * (-P[j] / a[j] + Q[j])
                v = 2 * F[j] ** 2 - g.K0 + Dd[j] - pr
            A[j][k] = v
            A[k][j] = v
    return A, g, info


def Fn(g, n, r):
    """F_n(r), r real arb"""
    L, a, c = g.L, g.a, g.c
    return c[n] * ((((r - a[n]) * L).sin()) / (r - a[n]) + (((r + a[n]) * L).sin()) / (r + a[n]))


def ldl(A, ball):
    """LDL^T of a symmetric list-of-lists matrix; ball=False strips radii (floating point), ball=True keeps them. Returns (D, Lrows) or None when a pivot is not certified positive (ball)"""
    n = len(A)
    D = []
    Lr = []
    M = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            s = A[i][j]
            if i == j:
                if j:
                    s = s - sum((row[k] * row[k] * D[k] for k in range(j)), arb(0))
            elif j:
                mj = M[j]
                s = s - sum((row[k] * mj[k] for k in range(j)), arb(0))
            if not ball:
                s = s.mid()
            if i == j:
                if not (s > 0):
                    return None
                D.append(s)
                row.append(arb(1))
            else:
                v = s / D[j]
                row.append(v if ball else v.mid())
        Lr.append(row)
        M.append([row[k] * D[k] for k in range(i + 1)])
    return D, Lr


def ldl_solve(D, Lr, b):
    n = len(D)
    y = [None] * n
    for i in range(n):
        s = b[i]
        if i:
            s = s - sum((Lr[i][k] * y[k] for k in range(i)), arb(0))
        y[i] = s.mid()
    z = [(y[i] / D[i]).mid() for i in range(n)]
    x = [None] * n
    for i in range(n - 1, -1, -1):
        s = z[i]
        if i < n - 1:
            s = s - sum((Lr[k][i] * x[k] for k in range(i + 1, n)), arb(0))
        x[i] = s.mid()
    return x


def eps_upper(A, iters=8):
    """inverse iteration in floating arb, Rayleigh quotient with the ball matrix: returns (x, ball of the quotient)"""
    n = len(A)
    Amid = [[x.mid() for x in row] for row in A]
    f = ldl(Amid, False)
    if f is None:
        return None, None
    D, Lr = f
    x = [arb(1) / (1 + j * j) for j in range(n)]
    for it in range(iters):
        y = ldl_solve(D, Lr, x)
        nr = sum((v * v for v in y), arb(0)).sqrt()
        x = [(v / nr).mid() for v in y]
    num = sum((x[j] * A[j][k] * x[k] for j in range(n) for k in range(n)), arb(0))
    den = sum((v * v for v in x), arb(0))
    return x, num / den


def shifted(A, t):
    B = [row[:] for row in A]
    for i in range(len(B)):
        B[i][i] = B[i][i] - t
    return B


def run_cert(task):
    X, N = task
    t0 = time.time()
    out = dict(X=X, N=N)
    for extra in (60, 120, 240, 480, 960, 1920, 3840):
        D = int(3.5 * N + extra)
        A, g, info = build(X, N, D)
        r = ldl(A, True)
        out.update(D=D, **info)
        if r is not None:
            break
    out['t_build_s'] = round(time.time() - t0, 1)
    if r is None:
        out.update(certified_positive=False, note='LDL pivots not certified at D=%d' % D)
        return out
    out['certified_positive'] = True
    out['min_pivot_log10'] = tof(min(d.lower() for d in r[0]).log() / arb(10).log())
    out['max_entry_radius_log10'] = tof(max(A[i][j].rad() for i in range(N + 1) for j in range(N + 1)).log() / arb(10).log()) if any(A[i][j].rad() > 0 for i in range(N + 1) for j in range(N + 1)) else None
    x, q = eps_upper(A)
    up = q.upper()
    out['eps_up'] = up.str(12, radius=False)
    out['eps_up_log10'] = tof(up.log() / arb(10).log())
    out['rayleigh_radius_rel'] = tof(q.rad() / up) if q.rad() > 0 else 0.0
    lo = None
    for m in (0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
        t = up * (1 - arb(m))
        if ldl(shifted(A, t), True) is None:
            break
        lo = (t, m)
    if lo:
        out['eps_lo'] = lo[0].mid().str(12, radius=False)
        out['rel_gap'] = lo[1]
    else:
        out['eps_lo'] = None
    out['seconds'] = round(time.time() - t0, 1)
    return out


def load_zeros(path):
    return [l.split()[-1] for l in open(path) if l.strip()]


def run_float(task, zeros=None):
    import mpmath as mp
    X, N, *rest = task
    D = int(3.5 * N + 30)
    t0 = time.time()
    A, g, info = build(X, N, D)
    mp.mp.dps = D + 10
    n = N + 1
    Am = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Am[i, j] = mp.mpf(A[i][j].mid().str(D + 10, radius=False))
    out = dict(X=X, N=N, dps=D, L=tof(g.L), **info)
    x, q = eps_upper(A)
    out['eps_inv_iter'] = q.mid().str(12, radius=False)
    if N <= 40:
        E, Q = mp.eigsy(Am)
        order = sorted(range(n), key=lambda i: E[i])
        eps = E[order[0]]
        xi = [Q[j, order[0]] for j in range(n)]
        out['eps'] = mp.nstr(eps, 12)
        out['second'] = mp.nstr(E[order[1]], 8)
        out['all_eigs_positive'] = bool(E[order[0]] > 0)
    else:
        xi = [mp.mpf(v.str(D + 10, radius=False)) for v in x]
        out['eps'] = out['eps_inv_iter']
    L = mp.mpf(g.L.str(D + 10, radius=False))
    cc = [1 / mp.sqrt(2 * L)] + [1 / mp.sqrt(L)] * N
    aa = [mp.pi * j / L for j in range(n)]

    def xh(r):
        return sum(xi[j] * cc[j] * (L * mp.sinc((r - aa[j]) * L) + L * mp.sinc((r + aa[j]) * L)) for j in range(n))
    top = mp.pi * N / L
    roots = []
    step = mp.mpf(1) / 40
    r0 = mp.mpf(step)
    f0 = xh(r0)
    r = r0
    mp.mp.dps = D + 10
    while r < top:
        r1 = r + step
        f1 = xh(r1)
        if f0 * f1 < 0:
            rt = mp.findroot(xh, (r, r1), solver='anderson', tol=mp.mpf(10) ** (-(D - 5)), maxsteps=200)
            roots.append(rt)
        r, f0 = r1, f1
    out['nroots'] = len(roots)
    out['roots'] = [mp.nstr(z, 40) for z in roots[:6]]
    if zeros is not None:
        gam = [mp.mpf(z) for z in zeros]
        def near(gm):
            return min((abs(z - gm) for z in roots), default=None)
        for k in (1, 5):
            gk = gam[k - 1] if D < 60 else mp.zetazero(k).imag
            d = near(gk)
            out['err_gamma%d' % k] = mp.nstr(d, 4) if d is not None else None
        m = 0
        cnt = 0
        for gm in gam:
            if gm > top:
                break
            cnt += 1
            d = near(gm)
            if d is not None and d < mp.mpf(10) ** -3:
                m += 1
        out['matched_1e-3'] = m
        out['ordinates_below_top'] = cnt
    out['seconds'] = round(time.time() - t0, 1)
    return out


def run_validate(zpath):
    import mpmath as mp
    zer = load_zeros(zpath)
    res = []
    random.seed(1)
    for X, N in ((5, 8), (13, 8), (50, 12), (200, 12)):
        D = 60
        A, g, info = build(X, N, D)
        ctx.prec = bits_for(D)
        L = g.L
        # test vector with f(L) = sum xi_n c_n (-1)^n = 0: then f is C^1 and F ~ r^-3
        xi = [arb(random.uniform(-1, 1)) / (1 + n) ** 2 for n in range(N + 1)]
        s = sum((xi[n] * g.c[n] * (-1) ** n for n in range(N + 1)), arb(0))
        xi[0] = xi[0] - s / g.c[0]
        W = sum((xi[j] * A[j][k] * xi[k] for j in range(N + 1) for k in range(N + 1)), arb(0))

        def Fhat(r):
            return sum((xi[n] * Fn(g, n, r) for n in range(N + 1)), arb(0))
        zs = 0
        part = {}
        for i, z in enumerate(zer, 1):
            r = arb(z)
            zs += 2 * Fhat(r) ** 2
            if i in (400, 1700):
                part[i] = zs
        row = dict(X=X, N=N, W_matrix=W.str(14, radius=False), W_matrix_rad=W.rad().str(3), zero_sum_400=part[400].str(14, radius=False), zero_sum_1700=part[1700].str(14, radius=False),
                   rel_diff_400=tof(abs(W - part[400]) / abs(W)), rel_diff_1700=tof(abs(W - part[1700]) / abs(W)),
                   last_term=tof(2 * Fhat(arb(zer[-1])) ** 2))
        # the bump f = cos^2(pi u/2L) = sqrt(2L)/2 phi_0 + sqrt(L)/2 phi_1 (the run book's check)
        b = [(2 * L).sqrt() / 2, L.sqrt() / 2] + [arb(0)] * (N - 1)
        Wb = sum((b[j] * A[j][k] * b[k] for j in range(N + 1) for k in range(N + 1)), arb(0))
        zb = 0
        pb = {}
        for i, z in enumerate(zer, 1):
            r = arb(z)
            zb += 2 * (b[0] * Fn(g, 0, r) + b[1] * Fn(g, 1, r)) ** 2
            if i in (400, 1700):
                pb[i] = zb
        row.update(bump_W=Wb.str(14, radius=False), bump_zero_sum_400=pb[400].str(14, radius=False), bump_zero_sum_1700=pb[1700].str(14, radius=False))
        # F_n(i/2) closed form against the integral
        mp.mp.dps = 30
        n1 = 1
        val = mp.quad(lambda u: mp.mpf(g.c[n1].str(30, radius=False)) * mp.cos(mp.mpf(g.a[n1].str(30, radius=False)) * u) * mp.cosh(u / 2), [-mp.mpf(L.str(30, radius=False)), mp.mpf(L.str(30, radius=False))])
        row['Fn_half_check'] = [mp.nstr(val, 15), (g.c[n1] * (-1) ** n1 * (L / 2).sinh() / (g.a[n1] ** 2 + arb(1) / 4)).str(15, radius=False)]
        res.append(row)
        print(json.dumps(row), flush=True)
    return res


def run_selftest():
    """entries built at two precisions (different node counts and ellipses) must have overlapping balls; reports the actual difference of the midpoints against the claimed radius"""
    res = []
    for X, N in ((2, 16), (13, 20), (200, 20), (1000, 24), (10000, 20)):
        D1 = int(3.5 * N + 30)
        A1, g1, i1 = build(X, N, D1)
        A2, g2, i2 = build(X, N, D1 + 100)
        ok = True
        worst = 0.0
        ctx.prec = bits_for(D1 + 100)
        for j in range(N + 1):
            for k in range(N + 1):
                x, y = A1[j][k], A2[j][k]
                if not x.overlaps(y):
                    ok = False
                d = abs(x.mid() - y.mid())
                r = x.rad() + y.rad()
                worst = max(worst, tof(d / r))
        row = dict(X=X, N=N, D1=D1, D2=D1 + 100, nodes1=i1['nodes'], nodes2=i2['nodes'], all_overlap=ok, max_middiff_over_radius=worst)
        print(json.dumps(row), flush=True)
        res.append(row)
    return res


def _cert(t):
    return run_cert(t)


def _float(args):
    t, zeros = args
    return run_float(t, zeros)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['validate', 'float', 'cert', 'selftest'])
    ap.add_argument('zeros', nargs='?', default='zeros_1700.txt')
    ap.add_argument('--tasks', default='')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    tasks = [tuple(int(v) for v in t.split(':')) for t in a.tasks.split(',') if t]
    if a.mode == 'validate':
        res = run_validate(a.zeros)
    elif a.mode == 'selftest':
        res = run_selftest()
    elif a.mode == 'cert':
        with Pool(a.procs) as p:
            res = []
            for r in p.imap_unordered(_cert, tasks):
                print(json.dumps(r), flush=True)
                res.append(r)
    else:
        zs = load_zeros(a.zeros)
        with Pool(a.procs) as p:
            res = []
            for r in p.imap_unordered(_float, [(t, zs) for t in tasks]):
                print(json.dumps(r), flush=True)
                res.append(r)
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)
