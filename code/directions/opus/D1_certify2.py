#!/usr/bin/env python3
"""D1: second certificate, with the prime term kept as an operator (D1.md, Section 3.6).

Lower bound in ball arithmetic. With m_s just below m(T) (m is increasing in |t|),
  W(f) >= m_s |f|^2 + K(f) - kappa T(f),
  K(f) = Pole(f) + (1/2 pi) int_{|t|<=T} |F|^2 (m - m_s) dt,
where T(f) = sum_j w_j g(c_j) is the exact prime form. By the Schur test T(v) <= tau |v|^2 with
tau = sup_x V(x), V(x) = (1/2)[A(L+|x|) + A(L-|x|)], A(y) = sum_{c_j <= y} w_j. In the normalised
Legendre basis with V = orders <= Nmax of one parity, Q the projection on the complement,
R = (T^2)_V - (T_V)^2 (so |QTu|^2 = u^t R u), and D0 = m_s - epsD - kappa tau - c, the bound
W >= c |f|^2 on the sector follows from positivity of the finite matrix
  m_s I + K_V - kappa T_V - c I - (1+s) kappa^2 R/D0 - (1+1/s)(epsB^2/D0) I.
T_V and (T^2)_V are exact integrals of piecewise polynomials (Gauss Legendre on each piece).
Positivity is verified by a midpoint Cholesky factor Lt with exact dyadic entries, and the ball
test |A - Lt Lt^t|_F |Lt^{-1}|_F^2 < 1, with |Lt^{-1}|_F bounded through an approximate inverse.

Usage:
  timeout 600 python3 -I D1_certify2.py rung X T Nmax ng prec
  timeout 600 python3 -I D1_certify2.py window L T Nmax ng prec
  timeout 600 python3 -I D1_certify2.py arch L T Nmax ng prec        kappa = 0
  timeout 600 python3 -I D1_certify2.py kappa X T Nmax ng prec k1,k2,... each kappa tested at c = 0
  timeout 600 python3 -I D1_certify2.py kedge X T Nmax ng prec       edges of the weight of the prime 2 (X = 3)
  timeout 600 python3 -I D1_certify2.py upper2 X_or_L N 0 0 prec     upper bound, ball Rayleigh quotient in N cosine modes
  timeout 600 python3 -I D1_certify2.py hyper X T Nmax ng prec rho   M - kappa T >= c on a fixed hyperplane, kappa = +1, -1
  timeout 600 python3 -I D1_certify2.py rungw X T Nmax ng prec [e|o] as rung, tau from the weighted Schur test; T = 0 for automatic T, Nmax; ng = 0 prints tau_w only
T must be a multiple of 8, at least 16.
"""
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import D1_certify as C
from D1_certify import arb, acb, arb_mat, ctx, up, lo, m_of, gl_nodes, bessel_table, spherical_i, record


def panels(T):
    """(a, b, rho, y, Phi_psi) for the graded panels. See D1.md, Section 3.4."""
    golden = (1 + arb(5).sqrt()) / 2
    out = []
    for j in range(8):
        out.append((j, j + 1, golden, arb(1) / 4, (arb(8) / 2 + arb(14) / 10).log() + arb("10.09")))
    Phi = (arb(T) / 2 + 6).log() + arb("2.4")
    for a in (8, 12):
        out.append((a, a + 4, arb(2), arb(3) / 2, Phi))
    a = 16
    while a < T:
        out.append((a, a + 8, arb(2), arb(3), Phi))
        a += 8
    return out


def assemble_graded(L, T, Nmax, ng):
    t0 = time.time()
    near = gl_nodes(ng)
    ngf = max(40, int(0.7 * ng))
    far = gl_nodes(ngf)
    half = arb(1) / 2
    nodes, weights = [], []
    for (a, b, rho, y, Phi) in panels(T):
        xs, ws = near if b - a == 1 else far
        w = arb(b - a)
        c = arb(a + b) / 2
        for x, wg in zip(xs, ws):
            nodes.append(c + half * w * x)
            weights.append(half * w * wg)
    pref = [(2 * L * (2 * n + 1)).sqrt() for n in range(Nmax + 1)]
    table = []
    for t in nodes:
        tab = bessel_table(L * t, Nmax)
        table.append([pref[n] * tab[n] for n in range(Nmax + 1)])
    mt = [m_of(t) for t in nodes]
    return dict(nodes=nodes, weights=weights, table=table, mt=mt, cosm=[], nn=len(nodes), t_assemble=time.time() - t0, ngf=ngf)


def quad_error_graded(L, T, Nmax, ng, ngf, ms):
    tot = arb(0)
    lp = C.PI.log()
    for (a, b, rho, y, Phi) in panels(T):
        n = ng if b - a == 1 else ngf
        Mh = 2 * L * (2 * Nmax + 1) * (2 * L * y).exp() * (Phi + lp + abs(ms))
        tot += (arb(b - a) / 2) * 8 * Mh * rho ** (-2 * n) / (1 - 1 / rho) / C.PI
    return tot


def pieces(L, comb):
    """Pieces of [0,L] on which all the indicators are constant, with the rigorous ends."""
    Lf = float(L)
    pts = []
    for (_, c, _) in comb:
        f = abs(Lf - float(c))
        if 1e-13 < f < Lf - 1e-13:
            pts.append((f, abs(L - c)))
    pts.sort(key=lambda z: z[0])
    merged = []
    for f, b in pts:
        if merged and abs(merged[-1][0] - f) < 1e-11:
            continue
        merged.append((f, b))
    edges = [(0.0, arb(0))] + merged + [(Lf, L)]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b[0] - a[0] < 1e-11:
            continue
        cen = 0.5 * (a[0] + b[0])
        act = []
        for (_, c, w) in comb:
            cf = float(c)
            act.append((c, w, cen >= cf - Lf, cen <= Lf - cf))
        out.append((a[1], b[1], act))
    return out


def tau_sup(L, comb):
    """sup_x V(x), the Schur bound for T, over the pieces (V is constant on each piece)."""
    best = arb(0)
    for (_, _, act) in pieces(L, comb):
        v = arb(0)
        for (c, w, im, ip) in act:
            if im:
                v += w / 2
            if ip:
                v += w / 2
        if up(v) > up(best):
            best = v
    return best


def tau_weighted(L, comb, deg=6, M=1500):
    """Weighted Schur test: |T(u)| <= sup_x (T h)(x)/h(x) |u|^2 for any positive h. D1.md, Section 3.6.
    h(x) = exp(sum_k b_k (x/L)^(2k)) with b chosen in double precision; the supremum is then bounded in
    ball arithmetic on M cells of [0,L], each ball covering its cell, so the bound holds for the chosen h."""
    import numpy as np
    from scipy.optimize import minimize
    Lf = float(L)
    P = pieces(L, comb)
    grids = []
    for (a, b, act) in P:
        af, bf = float(a), float(b)
        xs = np.linspace(af, bf, 80, endpoint=False) + (bf - af) / 160
        grids.append((xs, [(float(c), float(w), im, ip) for (c, w, im, ip) in act]))

    def logh(x, bb):
        y = (x / Lf) ** 2
        out = 0 * y
        q = y.copy()
        for b_ in bb:
            out = out + b_ * q
            q = q * y
        return out

    def ratio(bb):
        best = 0.0
        for xs, act in grids:
            den = logh(xs, bb)
            num = 0 * xs
            for c, w, im, ip in act:
                if im:
                    num = num + 0.5 * w * np.exp(logh(xs - c, bb) - den)
                if ip:
                    num = num + 0.5 * w * np.exp(logh(xs + c, bb) - den)
            best = max(best, float(num.max()))
        return best
    r = minimize(ratio, np.zeros(deg), method="Nelder-Mead", options=dict(maxiter=6000, xatol=1e-3, fatol=1e-4))
    r = minimize(ratio, r.x, method="Powell", options=dict(maxiter=3000))
    bb = [arb(float(x)) for x in r.x]
    keep = ctx.prec
    ctx.prec = 128

    def hlog(x):
        y = (x / L) * (x / L)
        out = arb(0)
        q = y
        for b_ in bb:
            out = out + b_ * q
            q = q * y
        return out
    best = 0.0
    for (a, b, act) in P:
        af, bf = float(a), float(b)
        m = max(8, int(M * (bf - af) / Lf))
        wd = (b - a) / m
        rad = up(wd) * 0.5 * (1 + 1e-9)
        for k in range(m):
            x = a + (k + arb(1) / 2) * wd + arb(0, rad)
            den = hlog(x)
            num = arb(0)
            for (c, w, im, ip) in act:
                if im:
                    num = num + (w / 2) * (hlog(x - c) - den).exp()
                if ip:
                    num = num + (w / 2) * (hlog(x + c) - den).exp()
            best = max(best, up(num))
    ctx.prec = keep
    return arb(best) * (1 + arb(1) / 10 ** 12), [float(x) for x in r.x], float(r.fun)  # the factor covers the rounding of up()


def prime_tables(L, comb, Nmax):
    """Values at the Gauss nodes of [0,L]: (node weight, p_n(x), (T p_n)(x)), n = 0..Nmax."""
    t0 = time.time()
    nx = Nmax + 1
    xs, ws = gl_nodes(nx)
    half = arb(1) / 2
    nrm = [((2 * n + 1) / (2 * L)).sqrt() for n in range(Nmax + 1)]
    A1 = [arb(2 * n + 1) / (n + 1) for n in range(Nmax + 1)]
    A2 = [arb(n) / (n + 1) for n in range(Nmax + 1)]
    iL = 1 / L

    def leg(u):
        out = [arb(1)] * (Nmax + 1)
        if Nmax >= 1:
            out[1] = u
        pm, pc = arb(1), u
        for n in range(1, Nmax):
            pn = A1[n] * u * pc - A2[n] * pm
            out[n + 1] = pn
            pm, pc = pc, pn
        return out

    rows0, rowsT, wts = [], [], []
    for (a, b, act) in pieces(L, comb):
        for x, wg in zip(xs, ws):
            xt = half * (b + a) + half * (b - a) * x
            base = leg(xt * iL)
            row = [arb(0)] * (Nmax + 1)
            for (c, w, im, ip) in act:
                wh = w * half
                if im:
                    v = leg((xt - c) * iL)
                    row = [r + wh * q for r, q in zip(row, v)]
                if ip:
                    v = leg((xt + c) * iL)
                    row = [r + wh * q for r, q in zip(row, v)]
            rows0.append([nrm[n] * base[n] for n in range(Nmax + 1)])
            rowsT.append([nrm[n] * row[n] for n in range(Nmax + 1)])
            wts.append(2 * half * (b - a) * wg)
    return rows0, rowsT, wts, time.time() - t0


def prime_blocks(tables, idx):
    rows0, rowsT, wts, _ = tables
    nn = len(wts)
    nm = len(idx)
    P0 = arb_mat(nn, nm, [rows0[i][j] for i in range(nn) for j in idx])
    E = arb_mat(nn, nm, [rowsT[i][j] for i in range(nn) for j in idx])
    Ew = arb_mat(nn, nm, [wts[i] * rowsT[i][j] for i in range(nn) for j in idx])
    return P0.transpose() * Ew, E.transpose() * Ew, nn


def chol_mid(A):
    """Cholesky of the midpoints of a lower triangular list of balls; dyadic factor or None."""
    n = len(A)
    Lm = [[None] * (i + 1) for i in range(n)]
    for i in range(n):
        Li = Lm[i]
        Ai = A[i]
        for j in range(i):
            Lj = Lm[j]
            x = Ai[j].mid()
            for k in range(j):
                x -= Li[k] * Lj[k]
            Li[j] = (x / Lj[j]).mid()
        x = Ai[i].mid()
        for k in range(i):
            x -= Li[k] * Li[k]
        if not x > 0:
            return None
        Li[i] = x.sqrt().mid()
    return Lm


def verify_pd(A):
    """True only if every real symmetric matrix in the balls of A is positive definite."""
    Lm = chol_mid(A)
    if Lm is None:
        return False, None
    n = len(A)
    flat = []
    for i in range(n):
        flat.extend(Lm[i])
        flat.extend([arb(0)] * (n - i - 1))
    LL = arb_mat(n, n, flat)
    B = LL * LL.transpose()
    F2 = arb(0)
    for i in range(n):
        for j in range(i + 1):
            e = A[i][j] - B[i, j]
            eu = abs(e.mid()) + e.rad()  # upper bound for |e|, so F2 is a sum of positive balls
            F2 += (eu * eu) * (1 if i == j else 2)
    # approximate inverse X with exact dyadic entries; no growth of radii, unlike a ball inverse
    X = [[arb(0)] * n for _ in range(n)]
    for col in range(n):
        for i in range(col, n):
            t = arb(1) if i == col else arb(0)
            Li_ = Lm[i]
            for k in range(col, i):
                t -= Li_[k] * X[k][col]
            X[i][col] = (t / Li_[i]).mid()
    XM = arb_mat(n, n, [X[i][j] for i in range(n) for j in range(n)])
    Sm = XM * LL
    S2 = arb(0)
    X2 = arb(0)
    for i in range(n):
        for j in range(n):
            e = Sm[i, j] - (1 if i == j else 0)
            eu = abs(e.mid()) + e.rad()
            S2 += eu * eu
            if j <= i:
                X2 += X[i][j] * X[i][j]
    sn = S2.sqrt()
    if not sn < 1:
        return False, Lm
    G2 = X2 / ((1 - sn) * (1 - sn))  # |Lt^{-1}|_F <= |X|_F / (1 - |I - X Lt|_F)
    return bool(F2.sqrt() * G2 < 1), Lm


def lam_estimate(A, Lm, iters=4):
    """Least eigenvalue of the midpoints by inverse iteration with the Cholesky factor."""
    n = len(A)
    x = [arb(1) / (i + 1) ** 2 for i in range(n)]
    def solve(b):
        y = [None] * n
        for i in range(n):
            t = b[i]
            Li = Lm[i]
            for k in range(i):
                t -= Li[k] * y[k]
            y[i] = (t / Li[i]).mid()
        z = [None] * n
        for i in range(n - 1, -1, -1):
            t = y[i]
            for k in range(i + 1, n):
                t -= Lm[k][i] * z[k]
            z[i] = (t / Lm[i][i]).mid()
        return z
    lam = None
    for _ in range(iters):
        z = solve(x)
        nz = sum((a * a for a in z), arb(0)).sqrt()
        nx = sum((a * a for a in x), arb(0)).sqrt()
        lam = float((nx / nz).mid())
        x = [(a / nz).mid() for a in z]
    return lam


class Sector:
    def __init__(self, L, T, Nmax, nd, tables, sector, ms, eq):
        self.sign = 1 if sector == "e" else -1
        idx, K0, _, a = C.sector_matrices(L, T, Nmax, nd, [], sector, ms)
        self.idx = idx
        self.n = len(idx)
        n = self.n
        wide = arb(0, up(eq) * 1.01)
        self.base = []
        for i in range(n):
            row = []
            for j in range(i + 1):
                x = K0[i, j] + self.sign * 2 * a[i] * a[j] + wide
                if i == j:
                    x = x + ms
                row.append(x)
            self.base.append(row)
        if tables is not None:
            TV, T2, _ = prime_blocks(tables, idx)
            R = T2 - TV * TV
            self.TV = [[TV[i, j] for j in range(i + 1)] for i in range(n)]
            self.R = [[R[i, j] for j in range(i + 1)] for i in range(n)]
        else:
            self.TV = None
            self.R = None

    def matrix(self, kappa, c, ms, tau, epsB, epsD, s_split=1e-6):
        D0 = ms - epsD - abs(kappa) * tau - arb(c)
        if not D0 > 0:
            raise SystemExit("D0 not positive: T too small")
        f1 = (1 + arb(s_split)) / D0
        shift = arb(c) + (1 + 1 / arb(s_split)) * epsB * epsB / D0
        A = []
        for i in range(self.n):
            row = []
            bi = self.base[i]
            for j in range(i + 1):
                x = bi[j]
                if self.TV is not None and kappa != 0:
                    x = x - kappa * self.TV[i][j] - f1 * kappa * kappa * self.R[i][j]
                if i == j:
                    x = x - shift
                row.append(x)
            A.append(row)
        return A

    def certify(self, kappa, ms, tau, epsB, epsD):
        A0 = self.matrix(kappa, 0, ms, tau, epsB, epsD)
        Lm = chol_mid(A0)
        if Lm is None:
            return None
        lam = lam_estimate(A0, Lm)
        c = 0.9995 * lam
        for _ in range(12):
            ok, _ = verify_pd(self.matrix(kappa, c, ms, tau, epsB, epsD))
            if ok:
                return c
            c *= 0.8
        return None


def setup_window(L, comb, T, Nmax, ng, tag, label, kappa_max=1, tau_value=None):
    t0 = time.time()
    tau = tau_value if tau_value is not None else (tau_sup(L, comb) if comb else arb(0))
    mT = m_of(arb(T))
    m0 = acb(arb(1) / 4, arb(0)).digamma().real - C.PI.log()
    ms = arb(float(mT.lower()) - 1e-12)
    epsB, epsD = C.tail_constants(L, T, Nmax, ms, arb(0), mT, m0)
    nd = assemble_graded(L, T, Nmax, ng)
    eq = quad_error_graded(L, T, Nmax, ng, nd["ngf"], ms)
    Atot = sum((w for (_, _, w) in comb), arb(0))
    print("%s: L=%.6f T=%d Nmax=%d ng=%d/%d prec=%d m_s=%.5f tau=%.5f A_L=%.5f D0=%.5f epsB=%.1e eq=%.1e primes=%s; nodes %d in %.1fs" % (label, float(L), T, Nmax, ng, nd["ngf"], ctx.prec, float(ms), float(tau), float(Atot), float(ms - kappa_max * tau), up(epsB), up(eq), [n for (n, _, _) in comb], nd["nn"], nd["t_assemble"]), flush=True)
    if not (ms - kappa_max * tau) > 0:
        raise SystemExit("need m(T) > tau")
    tables = None
    if comb:
        tables = prime_tables(L, comb, Nmax)
        print("  prime tables: %d space nodes in %.1fs" % (len(tables[2]), tables[3]), flush=True)
    return dict(tau=tau, ms=ms, epsB=epsB, epsD=epsD, nd=nd, eq=eq, tables=tables, A=Atot, t0=t0)


def run(L, comb, T, Nmax, ng, label, tag="cert2", kappa=1, tau_value=None, sectors="eo"):
    S = setup_window(L, comb, T, Nmax, ng, tag, label, tau_value=tau_value)
    res = dict(tag=tag, label=label, L=float(L), T=T, Nmax=Nmax, ng=ng, prec=ctx.prec, m_s=float(S["ms"]), tau=float(S["tau"]), A_L=float(S["A"]), epsB=up(S["epsB"]), eq=up(S["eq"]), primes=[n for (n, _, _) in comb])
    for sector in sectors:
        ts = time.time()
        sec = Sector(L, T, Nmax, S["nd"], S["tables"], sector, S["ms"], S["eq"])
        c = sec.certify(kappa if comb else 0, S["ms"], S["tau"], S["epsB"], S["epsD"])
        print("  sector %s modes=%d: certified W >= c |f|^2 with c=%s (%.1fs)" % (sector, sec.n, ("%.5e" % c) if c else "NONE", time.time() - ts), flush=True)
        res["c_" + sector] = c
        res["modes_" + sector] = sec.n
    res["seconds"] = time.time() - S["t0"]
    record(res)
    print("  total %.1fs" % res["seconds"], flush=True)
    return res


def main():
    mode = sys.argv[1]
    T = int(sys.argv[3]); Nmax = int(sys.argv[4]); ng = int(sys.argv[5])
    prec = int(sys.argv[6]) if len(sys.argv) > 6 else 256
    C.setup(prec)
    if mode == "rung":
        X = int(sys.argv[2])
        L = arb(X).log() / 2
        run(L, C.comb_terms(None, strict_int=X), T, Nmax, ng, "rung X=%d" % X)
    elif mode == "window":
        Ls = sys.argv[2]
        L = arb(Ls)
        run(L, C.comb_terms(2 * float(Ls)), T, Nmax, ng, "window L=%s" % Ls)
    elif mode == "arch":
        Ls = sys.argv[2]
        L = arb(Ls)
        run(L, [], T, Nmax, ng, "arch only L=%s" % Ls, tag="arch2")
    elif mode == "kappa":
        X = int(sys.argv[2])
        L = arb(X).log() / 2
        comb = C.comb_terms(None, strict_int=X)
        kaps = sys.argv[7].split(",")
        S = setup_window(L, comb, T, Nmax, ng, "kappa", "kappa X=%d" % X, kappa_max=1.1)
        out = {}
        for sector in "eo":
            sec = Sector(L, T, Nmax, S["nd"], S["tables"], sector, S["ms"], S["eq"])
            for ks in kaps:
                A = sec.matrix(arb(ks), 0, S["ms"], S["tau"], S["epsB"], S["epsD"])
                ok, _ = verify_pd(A)
                print("  sector %s kappa %s: W_kappa positive definite certified: %s" % (sector, ks, ok), flush=True)
                out["%s_%s" % (sector, ks)] = bool(ok)
        record(dict(tag="kappa2", X=X, T=T, Nmax=Nmax, ng=ng, prec=prec, res=out))
    elif mode == "upper2":
        # sys.argv[2]: X integer for L = (1/2) log X, or decimal L; T slot = cosine modes N (argv 3), argv 4,5 ignored
        arg = sys.argv[2]
        N = T
        if "." in arg:
            L = arb(arg)
            comb = C.comb_terms(2 * float(arg))
        else:
            L = arb(int(arg)).log() / 2
            comb = C.comb_terms(None, strict_int=int(arg))
        cm = C.cos_matrices(L, N, comb)
        res = dict(tag="upper2", L=float(L), N=N, prec=prec, primes=[n for (n, _, _) in comb])
        for sector in "eo":
            M = C.form_matrix(cm, sector, [1] * len(comb))
            n = len(M)
            low = [[M[i][j] for j in range(i + 1)] for i in range(n)]
            Lm = chol_mid(low)
            if Lm is None:
                print("  sector %s: finite section not positive at N=%d, N too large for prec or form negative" % (sector, N), flush=True)
                continue
            x = [arb(1) / (i + 1) ** 2 for i in range(n)]
            for _ in range(6):
                y = [None] * n
                for i in range(n):
                    t = x[i]
                    for k in range(i):
                        t -= Lm[i][k] * y[k]
                    y[i] = (t / Lm[i][i]).mid()
                z = [None] * n
                for i in range(n - 1, -1, -1):
                    t = y[i]
                    for k in range(i + 1, n):
                        t -= Lm[k][i] * z[k]
                    z[i] = (t / Lm[i][i]).mid()
                nz = sum((a * a for a in z), arb(0)).sqrt()
                x = [(a / nz).mid() for a in z]
            num = arb(0)
            for i in range(n):
                t = arb(0)
                for j in range(n):
                    t += M[i][j] * x[j]
                num += x[i] * t
            q = num / sum((a * a for a in x), arb(0))
            print("  L=%.6f N=%d sector %s: ball Rayleigh quotient, inf W/|f|^2 <= %.6e" % (float(L), N, sector, up(q)), flush=True)
            res["U_" + sector] = up(q)
        record(res)
    elif mode == "kedge":
        # edges of the weight of the first comb term, all other weights fixed at their true value 1 (one term here)
        X = int(sys.argv[2])
        L = arb(X).log() / 2
        comb = C.comb_terms(None, strict_int=X)
        S = setup_window(L, comb, T, Nmax, ng, "kedge", "kedge X=%d" % X, kappa_max=1.3)
        out = {}
        for sector in "eo":
            sec = Sector(L, T, Nmax, S["nd"], S["tables"], sector, S["ms"], S["eq"])
            def pd(k):
                return chol_mid(sec.matrix(arb(k), 0, S["ms"], S["tau"], S["epsB"], S["epsD"])) is not None
            for name, far in (("up", 1.3), ("down", 0.2)):
                good, bad = 1.0, far
                if pd(far):
                    good = far
                else:
                    for _ in range(34):
                        mid = 0.5 * (good + bad)
                        if pd(mid):
                            good = mid
                        else:
                            bad = mid
                safe = 1 + (good - 1) * (1 - 1e-6)
                ok, _ = verify_pd(sec.matrix(arb(safe), 0, S["ms"], S["tau"], S["epsB"], S["epsD"]))
                print("  sector %s edge %s: W_kappa >= 0 certified at kappa=%.12f: %s" % (sector, name, safe, ok), flush=True)
                out["%s_%s" % (sector, name)] = (safe, bool(ok))
        record(dict(tag="kedge", X=X, T=T, Nmax=Nmax, ng=ng, prec=prec, res=out))
    elif mode == "rungw":
        # as rung, with tau from the weighted Schur test; T = 0 chooses T and Nmax from tau_w. Optional argv[7]: sectors
        X = int(sys.argv[2])
        L = arb(X).log() / 2
        comb = C.comb_terms(None, strict_int=X)
        tw, bb, tf = tau_weighted(L, comb)
        print("  X=%d L=%.5f: weighted Schur tau_w=%.5f (double %.5f), unweighted sup V=%.5f; h coefficients %s" % (X, float(L), up(tw), tf, up(tau_sup(L, comb)), ["%.4f" % x for x in bb]), flush=True)
        if T == 0:
            T = 8 * math.ceil(2 * math.pi * math.exp(up(tw) + 0.75) / 8)
            Nmax = int(1.25 * math.e * float(L) * T / 2) + 10
        if T > 8 and ng > 0:
            run(L, comb, T, Nmax, ng, "rungw X=%d" % X, tag="cert2w", tau_value=tw, sectors=sys.argv[7] if len(sys.argv) > 7 else "eo")
    elif mode == "hyper":
        # first lemma on a hyperplane: M - kappa T + rho v v^t >= c for kappa = +1 and -1, with one
        # fixed polynomial v in each sector, the ground state of M in the Legendre basis (double
        # precision eigenvector, exact dyadic coefficients). Since v lies in V, <f,v> = u.v.
        import json
        import numpy as np
        X = int(sys.argv[2])
        L = arb(X).log() / 2
        comb = C.comb_terms(None, strict_int=X)
        rho = arb(sys.argv[7]) if len(sys.argv) > 7 else arb(100)
        S = setup_window(L, comb, T, Nmax, ng, "hyper", "hyper X=%d" % X)
        out = {}
        vecs = {}
        for sector in "eo":
            sec = Sector(L, T, Nmax, S["nd"], S["tables"], sector, S["ms"], S["eq"])
            n = sec.n
            def tofloat(A):
                F = np.zeros((n, n))
                for i in range(n):
                    for j in range(i + 1):
                        F[i, j] = F[j, i] = float(A[i][j].mid())
                return F
            F0 = tofloat(sec.matrix(arb(0), 0, S["ms"], S["tau"], S["epsB"], S["epsD"]))
            ev, U = np.linalg.eigh(F0)
            vf = [float(x) for x in U[:, 0]]
            v = [arb(x) for x in vf]
            vecs[sector] = dict(orders=[int(k) for k in sec.idx], coeffs=vf)
            print("  sector %s: M in the Legendre basis (lower form): lowest %.6e, next %.6e" % (sector, ev[0], ev[1]), flush=True)
            for kap in (1, -1):
                def mat(c):
                    A = sec.matrix(arb(kap), c, S["ms"], S["tau"], S["epsB"], S["epsD"])
                    return [[A[i][j] + rho * v[i] * v[j] for j in range(i + 1)] for i in range(n)]
                lam = float(np.linalg.eigvalsh(tofloat(mat(0)))[0])
                c = 0.999 * lam
                done = None
                for _ in range(10):
                    ok, _ = verify_pd(mat(c))
                    if ok:
                        done = c
                        break
                    c *= 0.9
                print("  sector %s kappa %+d: on v^perp certified M - kappa T >= c |f|^2 with c=%s (float estimate %.5e)" % (sector, kap, ("%.5e" % done) if done else "NONE", lam), flush=True)
                out["%s_%+d" % (sector, kap)] = done
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "D1_hyper_vectors.json"), "w") as fh:
            json.dump(dict(X=X, T=T, Nmax=Nmax, ng=ng, prec=prec, rho=sys.argv[7] if len(sys.argv) > 7 else "100", basis="normalised Legendre p_n(x)=sqrt((2n+1)/(2L)) P_n(x/L)", vectors=vecs), fh)
        record(dict(tag="hyper", X=X, T=T, Nmax=Nmax, ng=ng, prec=prec, res=out))
    else:
        raise SystemExit("unknown mode")


if __name__ == "__main__":
    main()
