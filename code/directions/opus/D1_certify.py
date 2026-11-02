#!/usr/bin/env python3
"""D1: ball arithmetic certificate for the Weil form on a window [-L,L].

Everything here is rigorous ball arithmetic (Arb through python-flint), except the
choice of test vectors and of the constants to test, which come from double precision
and need no justification: a pass of the ball LDL^T or a ball Rayleigh quotient is what
counts. The mathematics is written out in D1.md, Sections 2 and 3.

Lower bound (mode "cert"): for |t| >= T the symbol Psi(t) = m(t) - sum_j kappa_j w_j cos(t c_j)
is at least beta = m(T) - sum_j |kappa_j| w_j, because m is increasing in |t|. So
   W(f) >= beta |f|^2 + Pole(f) + (1/2 pi) int_{|t|<=T} |F|^2 (Psi - beta) dt,
and this bounded form is assembled in the normalised Legendre basis, with Gauss Legendre
on unit panels (explicit remainder from the Bernstein ellipse), an explicit tail beyond
order Nmax, and a ball LDL^T of M - (c + delta) I.

Upper bound (mode "upper"): the closed forms in the periodic cosine and sine modes
(D1.md, Section 2.3), and a ball Rayleigh quotient at an explicit vector.

Install (once): python3 -m pip install --only-binary=:all: --no-deps --target /tmp/d1lib python-flint
Usage:
  timeout 600 python3 -I D1_certify.py rung X T Nmax ng [prec]   window L = (1/2) log X, primes n < X
  timeout 600 python3 -I D1_certify.py window L T Nmax ng [prec] window L decimal, primes log n < 2L
  timeout 600 python3 -I D1_certify.py kappa T Nmax ng           weight interval of the prime 2 at L*=(1/2) log 3
  timeout 600 python3 -I D1_certify.py arch L T Nmax ng          arch form A+P alone (kappa = 0) at decimal L
  timeout 600 python3 -I D1_certify.py upper L_or_logX N         L decimal, or logX:3 for (1/2) log 3
"""
import sys, os, json, time, math
if os.path.isdir("/tmp/d1lib"):
    sys.path.insert(0, "/tmp/d1lib")
from flint import arb, acb, arb_mat, ctx
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "D1_certificates.jsonl")


def record(d):
    with open(OUT, "a") as fh:
        fh.write(json.dumps(d) + "\n")


def s(x, d=12):
    return x.str(d, radius=True)


def up(x):
    return float(x.upper())


def lo(x):
    return float(x.lower())


PI = None


def setup(prec):
    global PI
    ctx.prec = prec
    PI = arb.pi()


def von_mangoldt(n):
    """log p if n = p^k, else None (returns p)."""
    for p in range(2, n + 1):
        if n % p == 0:
            q = n
            while q % p == 0:
                q //= p
            return p if q == 1 else None
    return None


def comb_terms(limit_float, strict_int=None):
    """List of (arb c = log n, arb w = 2 Lambda(n)/sqrt n) for n with log n < 2L."""
    out = []
    n = 2
    while True:
        if strict_int is not None:
            if n >= strict_int:
                break
        elif math.log(n) >= limit_float:
            break
        p = von_mangoldt(n)
        if p is not None:
            out.append((n, arb(n).log(), 2 * arb(p).log() / arb(n).sqrt()))
        n += 1
    return out


def m_of(t):
    """m(t) = Re psi(1/4 + i t/2) - log pi at a real ball t."""
    return acb(arb(1) / 4, t / 2).digamma().real - PI.log()


def gl_nodes(ng):
    xs, ws = [], []
    for k in range(ng):
        x, w = arb.legendre_p_root(ng, k, weight=True)
        xs.append(x)
        ws.append(w)
    return xs, ws


def bessel_table(z, Nmax):
    """j_0(z),...,j_Nmax(z) at a positive ball z, by the exact downward recurrence
    from two rigorous starting values."""
    pref = (PI / (2 * z)).sqrt()
    half = arb(1) / 2
    jp = pref * z.bessel_j(arb(Nmax + 1) + half)
    jc = pref * z.bessel_j(arb(Nmax) + half)
    out = [None] * (Nmax + 1)
    out[Nmax] = jc
    iz = 1 / z
    for n in range(Nmax, 0, -1):
        jm = (2 * n + 1) * iz * jc - jp
        out[n - 1] = jm
        jp, jc = jc, jm
    return out


def spherical_i(n, z):
    return (PI / (2 * z)).sqrt() * z.bessel_i(arb(n) + arb(1) / 2)


def tail_sum(x, Nmax, L):
    """sum_{n>Nmax} 2L (2n+1) (x^n/(2n+1)!!)^2, rigorous with a geometric remainder."""
    u = arb(1)
    for n in range(1, Nmax + 2):
        u = u * x / (2 * n + 1)
    n = Nmax + 1
    term = 2 * L * (2 * n + 1) * u * u
    total = arb(0)
    x2 = x * x
    while True:
        total += term
        r = x2 / ((2 * n + 1) * (2 * n + 3))
        if up(r) < 0.5:
            total += term * r / (1 - r)
            return total
        term = term * r
        n += 1
        if n > Nmax + 5000:
            raise RuntimeError("tail not converging")


def assemble(L, T, Nmax, ng, comb):
    """Node data shared by both sectors. comb: list of (n, c, w)."""
    t0 = time.time()
    xs, ws = gl_nodes(ng)
    nodes = []
    weights = []
    half = arb(1) / 2
    for j in range(T):
        for x, w in zip(xs, ws):
            nodes.append(arb(j) + half + half * x)
            weights.append(half * w)
    nn = len(nodes)
    # q_n(t_k) = sqrt(2L(2n+1)) j_n(L t_k)
    pref = [(2 * L * (2 * n + 1)).sqrt() for n in range(Nmax + 1)]
    table = []
    for t in nodes:
        tab = bessel_table(L * t, Nmax)
        table.append([pref[n] * tab[n] for n in range(Nmax + 1)])
    mt = [m_of(t) for t in nodes]
    cosm = [[(t * c).cos() for t in nodes] for (_, c, _) in comb]
    return dict(nodes=nodes, weights=weights, table=table, mt=mt, cosm=cosm, nn=nn, t_assemble=time.time() - t0)


def sector_matrices(L, T, Nmax, nd, comb, sector, beta0):
    """Return K0 = Q diag(W (m-beta0)) Q^T, K1 list = Q diag(W w_j cos) Q^T, pole vector."""
    idx = list(range(0 if sector == "e" else 1, Nmax + 1, 2))
    sgn = [(-1) ** (n // 2) for n in idx]
    nm = len(idx)
    nn = nd["nn"]
    Wk = [w / PI for w in nd["weights"]]
    Qflat = []
    for i, n in enumerate(idx):
        sg = sgn[i]
        for k in range(nn):
            Qflat.append(sg * nd["table"][k][n])
    Q = arb_mat(nm, nn, Qflat)
    Qt_scaled = []
    # columns of Q^T scaled: build (nn x nm) with entry W_k phi_k Q_{i,k}
    def scaled(phi):
        ent = []
        for k in range(nn):
            f = Wk[k] * phi[k]
            for i, n in enumerate(idx):
                ent.append(f * (sgn[i] * nd["table"][k][n]))
        return arb_mat(nn, nm, ent)
    phi0 = [m - beta0 for m in nd["mt"]]
    K0 = Q * scaled(phi0)
    K1 = []
    for (j, (_, c, w)) in enumerate(comb):
        phi = [w * cs for cs in nd["cosm"][j]]
        K1.append(Q * scaled(phi))
    a = [(2 * L * (2 * n + 1)).sqrt() * spherical_i(n, L / 2) for n in idx]
    return idx, K0, K1, a


def quad_error(L, T, Nmax, ng, beta, wtot, cmax=None):
    """Entrywise bound for the quadrature error of (1/pi) int_0^T Q_m Q_n (Psi-beta) dt.
    Unit panels, Bernstein parameter rho = golden ratio, |Im t| <= 1/4. See D1.md."""
    rho = (1 + arb(5).sqrt()) / 2
    cmax = 2 * L if cmax is None else cmax  # every comb frequency c_j is below 2L
    Phi = (arb(T) / 2 + arb(14) / 10).log() + arb("10.09") + PI.log()
    Phi += wtot * (cmax / 4).cosh() + abs(beta)
    Mh = 2 * L * (2 * Nmax + 1) * (L / 2).exp() * Phi
    return (arb(T) / PI) * arb(1) / 2 * 8 * Mh * rho ** (-2 * ng) / (1 - 1 / rho)


def tail_constants(L, T, Nmax, beta, wtot, mT, m0):
    SigmaT = tail_sum(L * T, Nmax, L)
    sig_a2 = tail_sum(L / 2, Nmax, L) * (L * L / 4).exp()
    CPsi = max(up(beta - m0 + wtot), up(mT + wtot - beta))
    CPsi = arb(CPsi) * (1 + arb(1) / 10 ** 6)
    Acosh = (L + L.sinh()).sqrt()
    epsB = CPsi * (T * SigmaT / PI).sqrt() + 2 * Acosh * sig_a2.sqrt()
    epsD = CPsi * (T / PI) * SigmaT + 2 * sig_a2
    return epsB, epsD


def ldlt_positive(A):
    """Ball LDL^T of a symmetric list-of-lists (lower triangle used). True if every pivot is
    certainly positive, in which case every real symmetric matrix in the balls is positive definite."""
    n = len(A)
    Lm = [[None] * n for _ in range(n)]  # L D part stored as l_ij d_j
    D = [None] * n
    Lu = [[None] * n for _ in range(n)]
    for i in range(n):
        Ai = A[i]
        Lui = Lu[i]
        Lmi = Lm[i]
        for j in range(i):
            x = Ai[j]
            Luj = Lu[j]
            for k in range(j):
                x -= Lmi[k] * Luj[k]
            Lmi[j] = x
            Lui[j] = x / D[j]
        x = Ai[i]
        for k in range(i):
            x -= Lmi[k] * Lui[k]
        if not x > 0:
            return False, i, x
        D[i] = x
    return True, n, min(D, key=lambda z: up(z))


def build_test(K0, K1list, kappas, a, poleSign, beta0, c_total, eq):
    n = K0.nrows()
    wide = arb(0, eq) if isinstance(eq, float) else arb(0, up(eq) * 1.01)
    A = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            x = K0[i, j]
            for K1, kap in zip(K1list, kappas):
                x = x - kap * K1[i, j]
            x = x + poleSign * 2 * a[i] * a[j] + wide
            if i == j:
                x = x + beta0 - c_total
            row.append(x)
        A.append(row)
    return A


def search_c(K0, K1list, kappas, a, poleSign, beta0, eq, delta_of_c, c0=1e-2, steps=7):
    """Largest c (to a factor) with M - (c+delta) I certainly positive definite."""
    def test(c):
        A = build_test(K0, K1list, kappas, a, poleSign, beta0, arb(c) + delta_of_c(c), eq)
        return ldlt_positive(A)
    c = c0
    ok = test(c)[0]
    while not ok:
        c /= 10
        if c < 1e-60:
            return None
        ok = test(c)[0]
    hi = c * 10 if c < c0 else c
    lo_ = c
    if c < c0:
        for _ in range(steps):
            mid = math.sqrt(lo_ * hi)
            if test(mid)[0]:
                lo_ = mid
            else:
                hi = mid
    return lo_


def run_window(L, label, comb, T, Nmax, ng, kappas=None, tag="rung", c0=1e-2):
    t0 = time.time()
    if kappas is None:
        kappas = [1] * len(comb)
    wtot = sum((abs(kap) * w for kap, (_, _, w) in zip(kappas, comb)), arb(0))
    tT = arb(T)
    mT = m_of(tT)
    m0 = acb(arb(1) / 4, arb(0)).digamma().real - PI.log()
    beta = arb(lo(mT - wtot) - 1e-12)  # exact dyadic below m(T) - wtot, so Psi >= beta for |t| >= T
    epsB, epsD = tail_constants(L, T, Nmax, beta, wtot, mT, m0)
    eq = quad_error(L, T, Nmax, ng, beta, wtot)
    print("%s: L=%s T=%d Nmax=%d ng=%d prec=%d beta=%.6f epsB=%.2e epsD=%.2e eq=%.2e comb=%s" % (label, s(L, 10), T, Nmax, ng, ctx.prec, float(beta), up(epsB), up(epsD), up(eq), [n for (n, _, _) in comb]), flush=True)
    nd = assemble(L, T, Nmax, ng, comb)
    print("  assembled nodes=%d in %.1fs" % (nd["nn"], nd["t_assemble"]), flush=True)
    res = dict(tag=tag, label=label, L=float(L), T=T, Nmax=Nmax, ng=ng, prec=ctx.prec, beta=float(beta), epsB=up(epsB), epsD=up(epsD), eq=up(eq), kappas=[float(k) for k in kappas], primes=[n for (n, _, _) in comb])
    data = {}
    for sector, sign in (("e", 1), ("o", -1)):
        ts = time.time()
        idx, K0, K1, a = sector_matrices(L, T, Nmax, nd, comb, sector, beta)
        def delta(c):
            return epsB * epsB / (beta - epsD - arb(c))
        c = search_c(K0, K1, kappas, a, sign, beta, eq, delta, c0=c0)
        data[sector] = (idx, K0, K1, a, sign)
        rmax = max(float(K0[i, j].rad()) for i in range(K0.nrows()) for j in range(i + 1))
        print("  sector %s modes=%d: certified W >= c |f|^2 with c=%.6e (max ball radius %.1e, %.1fs)" % (sector, len(idx), c if c else float("nan"), rmax, time.time() - ts), flush=True)
        res["c_" + sector] = c
        res["modes_" + sector] = len(idx)
    res["seconds"] = time.time() - t0
    record(res)
    return res, data, (beta, epsB, epsD, eq)


# ---------- cosine and sine modes: closed forms in balls ----------

def cos_matrices(L, N, comb, ngx=600):
    """Ball matrices of A, P and the comb terms in the periodic modes. Section 2.3 of D1.md."""
    assert float(L) <= 2
    ngx = max(ngx, 3 * N)  # the node count must grow with the top frequency pi N / L
    # sup of |h| on the ellipse, rho = 2, foci 0 and 2L, by balls on its boundary and the maximum
    # principle: h is analytic inside since the semi-minor axis 0.75 L is below pi and z = 0 is removable
    keep = ctx.prec
    ctx.prec = 64
    Hsup = 0.0
    MB = 720
    for k in range(MB):
        th = (2 * PI * (k + arb(1) / 2) / MB) + arb(0, 3.1416 / MB)  # ball covering the k-th arc
        co, sn_ = th.cos(), th.sin()
        z = acb(L + L * arb("1.25") * co, L * arb("0.75") * sn_)
        val = (z / 2).exp() / z.sinh() - 1 / z
        Hsup = max(Hsup, up(abs(val)))
    ctx.prec = keep
    xs, ws = gl_nodes(ngx)
    nodes = [L * (1 + x) for x in xs]
    wts = [L * w for w in ws]
    hx = [(z / 2).exp() / z.sinh() - 1 / z for z in nodes]
    tri = [1 - z / (2 * L) for z in nodes]
    # remainder: rho = 2 on [0,2L], |h| <= Hsup, |1 - z/2L| <= 2.125 <= 2.4, trig <= cosh(om b), b = 0.75 L
    rho = arb(2)
    b = L * (rho - 1 / rho) / 2
    def eps(n):
        om = PI * n / L
        return L * 8 * arb(Hsup) * arb("2.4") * (om * b).cosh() * rho ** (-2 * ngx) / (1 - 1 / rho)
    In = [arb(0)] * (N + 1)
    Cn = [arb(0)] * (N + 1)
    for n in range(N + 1):
        om = PI * n / L
        cs = arb(0)
        sn = arb(0)
        for z, w, h, tr in zip(nodes, wts, hx, tri):
            si, co = (om * z).sin_cos()
            cs += w * h * tr * co
            sn += w * h * si
        e = arb(0, up(eps(n)))
        Cn[n] = cs + e
        In[n] = (sn + e + (2 * PI * n).si()) if n > 0 else arb(0)
    gam = arb.const_euler() if hasattr(arb, "const_euler") else arb("0.57721566490153286060651209008240243104215933593992")
    n1 = list(range(1, N + 1))
    Ae = [[arb(0)] * (N + 1) for _ in range(N + 1)]
    Ao = [[arb(0)] * N for _ in range(N)]
    logpi = PI.log()
    Ae[0][0] = -gam - logpi - (4 * L).log() + 1 - Cn[0]
    for n in n1:
        base = (arb(n) / (2 * L)).log() - (2 * PI * n).ci() - Cn[n]
        fr = In[n] / (2 * PI * n)
        Ae[n][n] = base + fr
        Ao[n - 1][n - 1] = base - fr
        v = (-1) ** n * In[n] / (arb(2).sqrt() * PI * n)
        Ae[0][n] = v
        Ae[n][0] = v
    for m in n1:
        for n in range(1, m):
            sg = (-1) ** (m + n)
            den = PI * (m * m - n * n)
            xe = -sg * (n * In[n] - m * In[m]) / den
            xo = -sg * (m * In[n] - n * In[m]) / den
            Ae[m][n] = xe
            Ae[n][m] = xe
            Ao[m - 1][n - 1] = xo
            Ao[n - 1][m - 1] = xo
    sh = (L / 2).sinh()
    ap = [(2 * L).sqrt() ** -1 * sh * 4] + [L.sqrt() ** -1 * (-1) ** n * sh / ((PI * n / L) ** 2 + arb(1) / 4) for n in n1]
    am = [-(L.sqrt() ** -1) * (-1) ** n * 2 * (PI * n / L) * sh / ((PI * n / L) ** 2 + arb(1) / 4) for n in n1]
    Te = []
    To = []
    for (_, c, w) in comb:
        Sn = [(PI * n * c / L).sin() for n in range(N + 1)]
        Co = [(PI * n * c / L).cos() for n in range(N + 1)]
        tric = 1 - c / (2 * L)
        Ge = [[arb(0)] * (N + 1) for _ in range(N + 1)]
        Go = [[arb(0)] * N for _ in range(N)]
        Ge[0][0] = w * tric
        for n in n1:
            sn = Sn[n] / (2 * PI * n)
            Ge[n][n] = w * (tric * Co[n] - sn)
            Go[n - 1][n - 1] = w * (tric * Co[n] + sn)
            v = -w * (-1) ** n * Sn[n] / (arb(2).sqrt() * PI * n)
            Ge[0][n] = v
            Ge[n][0] = v
        for m in n1:
            for n in range(1, m):
                sg = (-1) ** (m + n)
                den = PI * (m * m - n * n)
                xe = w * sg * (n * Sn[n] - m * Sn[m]) / den
                xo = w * sg * (m * Sn[n] - n * Sn[m]) / den
                Ge[m][n] = xe
                Ge[n][m] = xe
                Go[m - 1][n - 1] = xo
                Go[n - 1][m - 1] = xo
        Te.append(Ge)
        To.append(Go)
    return dict(Ae=Ae, Ao=Ao, ap=ap, am=am, Te=Te, To=To)


def tofloat(M):
    return np.array([[float(x.mid()) for x in row] for row in M])


def rayleigh(M, v):
    n = len(v)
    vv = [arb(float(x)) for x in v]
    num = arb(0)
    for i in range(n):
        Mi = M[i]
        t = arb(0)
        for j in range(n):
            t += Mi[j] * vv[j]
        num += vv[i] * t
    den = sum((x * x for x in vv), arb(0))
    return num / den


def form_matrix(cm, sector, kappas):
    """Python list ball matrix of A + P - sum kappa_j T_j."""
    if sector == "e":
        A, T, pv, sg = cm["Ae"], cm["Te"], cm["ap"], 2
    else:
        A, T, pv, sg = cm["Ao"], cm["To"], cm["am"], -2
    n = len(A)
    M = [[A[i][j] + sg * pv[i] * pv[j] for j in range(n)] for i in range(n)]
    for Tj, kap in zip(T, kappas):
        if kap != 0:
            for i in range(n):
                for j in range(n):
                    M[i][j] = M[i][j] - kap * Tj[i][j]
    return M


def upper_bound(cm, sector, kappas):
    M = form_matrix(cm, sector, kappas)
    F = tofloat(M)
    w, V = np.linalg.eigh(0.5 * (F + F.T))
    return rayleigh(M, V[:, 0]), w[0]


def main():
    mode = sys.argv[1]
    if mode == "rung":
        X = int(sys.argv[2]); T = int(sys.argv[3]); Nmax = int(sys.argv[4]); ng = int(sys.argv[5])
        prec = int(sys.argv[6]) if len(sys.argv) > 6 else 192
        setup(prec)
        L = arb(X).log() / 2
        comb = comb_terms(None, strict_int=X)
        run_window(L, "rung X=%d" % X, comb, T, Nmax, ng)
    elif mode == "window":
        Ls = sys.argv[2]; T = int(sys.argv[3]); Nmax = int(sys.argv[4]); ng = int(sys.argv[5])
        prec = int(sys.argv[6]) if len(sys.argv) > 6 else 192
        setup(prec)
        L = arb(Ls)
        comb = comb_terms(2 * float(Ls))
        run_window(L, "window L=%s" % Ls, comb, T, Nmax, ng)
    elif mode == "arch":
        Ls = sys.argv[2]; T = int(sys.argv[3]); Nmax = int(sys.argv[4]); ng = int(sys.argv[5])
        setup(192)
        L = arb(Ls)
        run_window(L, "arch only L=%s" % Ls, [], T, Nmax, ng, tag="arch", c0=1e-2)
    elif mode == "kappa":
        T = int(sys.argv[2]); Nmax = int(sys.argv[3]); ng = int(sys.argv[4])
        setup(192)
        L = arb(3).log() / 2
        comb = comb_terms(None, strict_int=3)
        (n2, c2, w2) = comb[0]
        kapmax = arb("1.1")
        wtot = kapmax * w2
        mT = m_of(arb(T))
        m0 = acb(arb(1) / 4, arb(0)).digamma().real - PI.log()
        beta = arb(lo(mT - wtot) - 1e-12)
        epsB, epsD = tail_constants(L, T, Nmax, beta, wtot, mT, m0)
        eq = quad_error(L, T, Nmax, ng, beta, wtot)
        nd = assemble(L, T, Nmax, ng, comb)
        delta = epsB * epsB / (beta - epsD)
        print("beta", float(beta), "epsB", up(epsB), "eq", up(eq), "assembled", nd["t_assemble"], flush=True)
        out = {}
        for sector, sign in (("e", 1), ("o", -1)):
            idx, K0, K1, a = sector_matrices(L, T, Nmax, nd, comb, sector, beta)
            for kap_s in sys.argv[5].split(","):
                kap = arb(kap_s)
                A = build_test(K0, K1, [kap], a, sign, beta, delta, eq)
                ok, i, x = ldlt_positive(A)
                print("sector %s kappa %s: positive definite %s" % (sector, kap_s, ok), flush=True)
                out["%s_%s" % (sector, kap_s)] = ok
        record(dict(tag="kappa", T=T, Nmax=Nmax, ng=ng, beta=float(beta), res=out))
    elif mode == "upper":
        arg = sys.argv[2]
        N = int(sys.argv[3])
        setup(192)
        if arg.startswith("logX:"):
            X = int(arg[5:])
            L = arb(X).log() / 2
            comb = comb_terms(None, strict_int=X)
        else:
            L = arb(arg)
            comb = comb_terms(2 * float(arg))
        t0 = time.time()
        cm = cos_matrices(L, N, comb)
        res = dict(tag="upper", L=float(L), N=N, primes=[n for (n, _, _) in comb])
        for sector in "eo":
            q, w = upper_bound(cm, sector, [1] * len(comb))
            print("L=%.6f N=%d sector %s: W form Rayleigh <= %.6e (float eig %.6e)" % (float(L), N, sector, up(q), w), flush=True)
            res["W_" + sector] = up(q)
            q, w = upper_bound(cm, sector, [0] * len(comb))
            print("L=%.6f N=%d sector %s: A+P Rayleigh <= %.6e (float eig %.6e)" % (float(L), N, sector, up(q), w), flush=True)
            res["AP_" + sector] = up(q)
        # optional weights beyond the edges
        if len(sys.argv) > 4:
            for spec in sys.argv[4].split(","):
                sector, kap_s = spec.split(":")
                q, w = upper_bound(cm, sector, [arb(kap_s)] + [1] * (len(comb) - 1))
                print("L=%.6f N=%d sector %s kappa %s: Rayleigh <= %.6e (float eig %.6e)" % (float(L), N, sector, kap_s, up(q), w), flush=True)
                res["W_%s_%s" % (sector, kap_s)] = up(q)
        res["seconds"] = time.time() - t0
        record(res)
    else:
        raise SystemExit("unknown mode")


if __name__ == "__main__":
    main()
