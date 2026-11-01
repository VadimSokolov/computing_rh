# wl_dh.py: Work package 1, the part "prime dependence" of rh_compute.tex Section "What to fit": the Weil form of the Davenport and Heilbronn function f, its odd modes, and the
# general cell routine (float eigenvalues, ball certification of the sign, zeros of xi_hat and their distances to the ordinates) used by wl_primes.py and wl_accuracy.py.
# Extends code/wl_weil.py by importing it (monkey patching its integrand family fam for the exponent and parity); nothing of Nick's misc/reviews code is used.
#
# The Davenport and Heilbronn function. kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1), a_n periodic mod 5 with (a_1..a_5) = (1, kappa, -kappa, -1, 0), f(s) = sum a_n n^-s
#   = 5^-s [zeta(s,1/5) + kappa zeta(s,2/5) - kappa zeta(s,3/5) - zeta(s,4/5)] = ((1 - i kappa)/2) L(s, chi) + ((1 + i kappa)/2) L(s, conj chi), chi the character mod 5 with chi(2) = i.
# chi is odd (chi(-1) = chi(4) = -1), conductor 5, root number eps(chi) = tau(chi)/(i sqrt 5), and the combination is chosen so that b/a = (1 + i kappa)/(1 - i kappa) = eps(chi): then
#   Lambda_f(s) = (5/pi)^((s+1)/2) Gamma((s+1)/2) f(s) = Lambda_f(1 - s)    (entire, real on Re s = 1/2, Lambda_f(conj s) = conj Lambda_f(s)),
# checked numerically in the `zeros` mode. f has no Euler product and no pole.
# Explicit formula (derivation as for zeta, rh_compute.tex eq:W, with the three changes): put rho = 1/2 + i gamma, g(u) = int f(v) f(v-u) dv, ghat(z) = int g e^{izu} = F(z) F(-z) (entire).
#   sum over the zeros rho of Lambda_f of ghat(gamma_rho) = (1/2 pi i) [int_{Re s = c} - int_{Re s = 1-c}] ghat(-i(s - 1/2)) Lambda_f'/Lambda_f(s) ds, c large; the two lines contribute equally
#   (Lambda_f(s) = Lambda_f(1-s), ghat even), and on Re s = c: Lambda_f'/Lambda_f = (1/2) log(5/pi) + (1/2) psi((s+1)/2) - sum Lambda_f(n) n^-s, where -f'/f = sum Lambda_f(n) n^-s (n >= 2; Dirichlet
#   coefficients of the logarithmic derivative, not supported on primes: Lambda_f(n) = a_n log n - sum_{d | n, d < n} Lambda_f(d) a_{n/d}). Since (1/2 pi) int ghat(r - i sigma) n^(-1/2 - sigma - i r) dr = n^(-1/2) g(log n),
#   W_DH(f) = sum_rho ghat(gamma_rho) = g(0) (log(5/pi) - gamma_E - log(1 - e^{-4L})) + int_0^{2L} 2 (g(0) e^{-2u} - e^{-3u/2} g(u))/(1 - e^{-2u}) du - 2 sum_{n <= X} Lambda_f(n) n^{-1/2} g(log n).
#   Changes against zeta: no pole term 2F(i/2)^2 (f is entire), log pi -> log pi - log 5 (conductor), psi(s/2) -> psi((s+1)/2) so e^{-u/2} -> e^{-3u/2}, and Lambda(n) -> Lambda_f(n) (all n, real).
#   The zeros rho run over all zeros of Lambda_f: on the line (gamma real), and the quadruples rho, conj rho, 1 - rho, 1 - conj rho (gamma = +-t +- i delta, delta = Re rho - 1/2) when Re rho != 1/2.
# Odd modes psi_n = L^{-1/2} sin(a_n u), n = 1..N (they vanish at +-L). Product formula: g_jk(u) = (c_j c_k/2)(-1)^(j+k) [ (s_j + s_k)/(a_j + a_k) + (s_k - s_j)/(a_j - a_k) ], g_kk = (c_k^2/2)[ s_k/a_k + (2L-u) cos(a_k u) ]
#   (even case: the sign of the (s_j+s_k) term and of s_k/a_k flipped). F_n(z) = i c [sin((z-a)L)/(z-a) - sin((z+a)L)/(z+a)], F_n(i/2) = c 2 a_n (-1)^n sinh(L/2)/(a_n^2 + 1/4) (real), and F is odd so ghat(i/2) = F(i/2)F(-i/2) = -F(i/2)^2:
#   the pole term of zeta on odd functions is -2 F_j(i/2) F_k(i/2) (negative semidefinite), 2 F_j F_k on even ones.
#
# Usage (from data/ on Hopper, scripts and zeros_1700.txt beside):
#   python3 wl_dh.py zeros --T 200 --procs 40 --out wl_dh_zeros.json        (zeros of f: line zeros by sign changes, off line zeros by a complex secant search; counts against cf_dhzeros.json if present)
#   python3 wl_dh.py validate --zeros wl_dh_zeros.json --out wl_dh_validate.json
#   python3 wl_dh.py scan --tasks dh:even:13:16,... --procs 40 --out wl_dh_scan.json     (task = kind:parity:X:N[:Y], kind zeta or dh, Y = prime truncation for zeta)
import argparse, json, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import mpmath as mp
try:
    import wl_weil as W
    from flint import arb, acb, ctx
except ImportError:  # laptop without python-flint: only the zeros mode works
    W = None

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------- coefficients, weights
def kappa_arb():
    s5 = arb(5).sqrt()
    return ((10 - 2 * s5).sqrt() - 2) / (s5 - 1)


def lamf_list(X, kap):
    """Lambda_f(n), n <= X, as arb (index 0 unused)"""
    tab = [arb(0), arb(1), kap, -kap, arb(-1)]
    lam = [arb(0)] * (X + 1)
    acc = [arb(0)] * (X + 1)
    for d in range(1, X + 1):
        lam[d] = tab[d % 5] * arb(d).log() - acc[d]
        for q in range(2, X // d + 1):
            cq = tab[q % 5]
            if q % 5:
                acc[d * q] += lam[d] * cq
    return lam


def make_weights(spec, X):
    """list of (n, w_n = Lambda(n)/sqrt n) as arb at the current precision"""
    if spec['kind'] == 'zeta':
        Y = spec.get('Y') or X
        return [(n, arb(p).log() / arb(n).sqrt()) for n, p in W.prime_powers(min(Y, X))]
    lam = lamf_list(X, kappa_arb())
    return [(n, lam[n] / arb(n).sqrt()) for n in range(2, X + 1) if q_nonzero(lam[n])]


def q_nonzero(x):
    return not (x == 0)


def prime_sums_w(g, wts, odd):
    N = g.N
    P = [arb(0)] * (N + 1)
    Q = [arb(0)] * (N + 1)
    for n, w in wts:
        u = arb(n).log()
        t = 2 * g.L - u
        Q[0] += w * t
        for k in range(1, N + 1):
            s, co = (g.a[k] * u).sin_cos()
            P[k] += w * s
            Q[k] += w * t * co
    return P, Q


def make_fam(ex, odd):
    sgn = 1 if odd else -1

    def fam(g, u):
        L, N, a, c = g.L, g.N, g.a, g.c
        e1 = (-ex * u).exp()
        e2 = (-2 * u).exp()
        R = 1 / (1 - e2)
        Sv = [None] * (N + 1)
        Dv = [None] * (N + 1)
        t = 2 * L - u
        Dv[0] = 2 * (e2 - e1 * t / (2 * L)) * R
        for k in range(1, N + 1):
            s, co = (a[k] * u).sin_cos()
            Sv[k] = e1 * s * R
            gk = c[k] ** 2 / 2 * (sgn * s / a[k] + t * co)
            Dv[k] = 2 * (e2 - e1 * gk) * R
        return Sv, Dv
    return fam


def build_gen(X, N, D, spec):
    """matrix of balls of the form for spec = dict(kind='zeta'|'dh', parity='even'|'odd', Y=prime truncation); even: indices 0..N, odd: modes 1..N (size N)"""
    odd = spec['parity'] == 'odd'
    ex = 0.5 if spec['kind'] == 'zeta' else 1.5
    bits = W.bits_for(D)
    ctx.prec = bits
    g = W.Geo(X, N)
    old = W.fam
    W.fam = make_fam(arb(ex), odd)
    try:
        S, Dd, info = W.arch_integrals(g, D, bits)
    finally:
        W.fam = old
    ctx.prec = bits
    wts = make_weights(spec, X)
    P, Q = prime_sums_w(g, wts, odd)
    L, a, c = g.L, g.a, g.c
    sh = (L / 2).sinh()
    K0 = g.K0 if spec['kind'] == 'zeta' else g.K0 - arb(5).log()
    if spec['kind'] == 'zeta':
        if odd:
            F = [arb(0)] + [c[j] * 2 * a[j] * (-1) ** j * sh / (a[j] ** 2 + arb(1) / 4) for j in range(1, N + 1)]
            ps = -2
        else:
            F = [c[j] * (-1) ** j * sh / (a[j] ** 2 + arb(1) / 4) for j in range(N + 1)]
            ps = 2
    else:
        F = [arb(0)] * (N + 1)
        ps = 0
    Z = [arb(0)] + [S[j] + P[j] for j in range(1, N + 1)]
    lo = 1 if odd else 0
    A = [[None] * (N + 1) for _ in range(N + 1)]
    for j in range(lo, N + 1):
        for k in range(j, N + 1):
            if j != k:
                sg = (-1) ** (j + k)
                if odd:
                    v = ps * F[j] * F[k] - c[j] * c[k] * sg * ((Z[j] + Z[k]) / (a[j] + a[k]) + (Z[k] - Z[j]) / (a[j] - a[k]))
                else:
                    v = ps * F[j] * F[k] - c[j] * c[k] * sg * (-(Z[j] + Z[k]) / (a[j] + a[k]) + (Z[k] - Z[j]) / (a[j] - a[k]))
            else:
                if j == 0:
                    pr = 2 * c[0] ** 2 * Q[0]
                else:
                    pr = c[j] ** 2 * ((P[j] / a[j] if odd else -P[j] / a[j]) + Q[j])
                v = ps * F[j] ** 2 - K0 + Dd[j] - pr
            A[j][k] = v
            A[k][j] = v
    if odd:
        A = [row[1:] for row in A[1:]]
    info['nweights'] = len(wts)
    return A, g, info


# ----------------------------------------------------------------------------- transforms and zero sums
def Fz(g, b, z, odd):
    """F(z) = sum b_n F_n(z) (even) or G(z) = sum b_n G_n(z) with F = i G (odd; b indexed 1..N as list positions 0..N-1), z acb or arb"""
    L, a, c = g.L, g.a, g.c
    tot = acb(0)
    for pos, bn in enumerate(b):
        n = pos + 1 if odd else pos
        if odd:
            tot += bn * c[n] * (((z - a[n]) * L).sin() / (z - a[n]) - ((z + a[n]) * L).sin() / (z + a[n]))
        else:
            tot += bn * c[n] * (((z - a[n]) * L).sin() / (z - a[n]) + ((z + a[n]) * L).sin() / (z + a[n]))
    return tot


def ghat(g, b, z, odd):
    """ghat(z) = F(z) F(-z); odd: F = iG with G odd, so ghat = i G(z) i G(-z) = G(z)^2"""
    if odd:
        v = Fz(g, b, z, True)
        return v * v
    return Fz(g, b, z, False) * Fz(g, b, -z, False)


def zero_terms(g, b, odd, on, off):
    """returns (sum over the line zeros of 2 ghat(t), list of the quadruple sums for the off line zeros)"""
    s_on = acb(0)
    for t in on:
        s_on += 2 * ghat(g, b, acb(t), odd)
    qs = []
    for a_, t in off:
        d = a_ - arb(1) / 2
        tot = acb(0)
        for sr in (1, -1):
            for si in (1, -1):
                tot += ghat(g, b, acb(sr * t, si * d), odd)
        qs.append(tot)
    return s_on, qs


# ----------------------------------------------------------------------------- zeros of f (mpmath)
_kap = {}


def kap_mp():
    d = mp.mp.dps
    if d not in _kap:
        s5 = mp.sqrt(5)
        _kap[d] = (mp.sqrt(10 - 2 * s5) - 2) / (s5 - 1)
    return _kap[d]


def fval(s):
    k = kap_mp()
    q = mp.mpf(1) / 5
    return mp.mpf(5) ** (-s) * (mp.zeta(s, q) + k * mp.zeta(s, 2 * q) - k * mp.zeta(s, 3 * q) - mp.zeta(s, 4 * q))


def Lam(s):
    return (5 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * fval(s)


def _scan_seg(args):
    t0, t1, h = args
    mp.mp.dps = 20
    out = []
    n = int(round((t1 - t0) / h))
    prev = Lam(mp.mpf(1) / 2 + 1j * mp.mpf(t0)).real
    tp = mp.mpf(t0)
    for i in range(1, n + 1):
        t = mp.mpf(t0) + i * mp.mpf(h)
        cur = Lam(mp.mpf(1) / 2 + 1j * t).real
        if prev * cur < 0:
            out.append((float(tp), float(t)))
        prev, tp = cur, t
    return out


def _refine_line(br):
    mp.mp.dps = 70
    a, b = mp.mpf(br[0]), mp.mpf(br[1])
    zf = lambda t: Lam(mp.mpf(1) / 2 + 1j * t).real
    try:
        r = mp.findroot(zf, (a, b), solver='anderson', tol=mp.mpf(10) ** -62, maxsteps=80)
    except Exception:
        r = mp.findroot(zf, (a, b), solver='illinois', tol=mp.mpf(10) ** -60, maxsteps=200, verify=False)
    return mp.nstr(r, 62)


def _newton_off(st):
    mp.mp.dps = 30

    def fb(s):
        if not (0.2 < s.real < 2.2) or not (st[1] - 6 < s.imag < st[1] + 6):
            raise ValueError('left the box')
        return fval(s)
    try:
        r = mp.findroot(fb, mp.mpc(st[0], st[1]), tol=mp.mpf(10) ** -24, maxsteps=30, verify=False)
    except Exception:
        return None
    if abs(fval(r)) > mp.mpf(10) ** -15 or not (mp.mpf(1) / 2 + mp.mpf(10) ** -6 < r.real < 2) or r.imag <= 0:
        return None
    return (float(r.real), float(r.imag))


def _refine_off(z):
    mp.mp.dps = 70
    r = mp.findroot(fval, mp.mpc(z[0], z[1]), tol=mp.mpf(10) ** -62, maxsteps=60, verify=False)
    return (mp.nstr(r.real, 62), mp.nstr(r.imag, 62))


def run_zeros(T, procs):
    t0 = time.time()
    mp.mp.dps = 30
    chk = dict(kappa=mp.nstr(kap_mp(), 25))
    mx = 0
    for s in (mp.mpc(0.3, 2), mp.mpc(0.7, 31), mp.mpc(1.4, 90)):
        mx = max(mx, abs(Lam(s) - Lam(1 - s)) / abs(Lam(s)))
    chk['functional_equation_max_rel'] = mp.nstr(mx, 5)
    chk['Lambda_imag_on_line_max_rel'] = mp.nstr(max(abs(Lam(mp.mpf(1) / 2 + 1j * t).imag / Lam(mp.mpf(1) / 2 + 1j * t).real) for t in (5, 40, 100)), 5)
    h = 0.02
    segs = [(i, min(i + 5, T), h) for i in range(0, int(T), 5)]
    with Pool(procs) as p:
        brs = [b for part in p.imap(_scan_seg, segs) for b in part]
        line = p.map(_refine_line, brs)
        starts = [(sg, t) for t in [0.5 + 0.8 * i for i in range(int(T / 0.8) + 2)] for sg in (0.55, 0.65, 0.75, 0.85, 0.95)]
        cand = [r for r in p.map(_newton_off, starts) if r is not None]
        uniq = []
        for r in sorted(cand, key=lambda z: z[1]):
            if r[1] <= T and not any(abs(r[0] - u[0]) + abs(r[1] - u[1]) < 1e-6 for u in uniq):
                uniq.append(r)
        off = p.map(_refine_off, uniq)
    out = dict(T=T, step=h, kappa=chk['kappa'], checks=chk, online=line, offline=off, n_online=len(line), n_offline=len(off), seconds=round(time.time() - t0, 1))
    cf = os.path.join(HERE, '..', 'data', 'cf_dhzeros.json')
    if not os.path.exists(cf):
        cf = os.path.join(HERE, 'cf_dhzeros.json')
    if os.path.exists(cf) and T == 200:
        d = json.load(open(cf))['T2']
        out['cf_check'] = dict(cf_sign_changes=d['sign_changes_line'], cf_offline=[(z['center_re'][:20], z['center_im'][:20]) for z in d['zeros']],
                               mine_offline=[(a[:20], b[:20]) for a, b in off])
    return out


# ----------------------------------------------------------------------------- validation
def load_dh_zeros(path):
    d = json.load(open(path))
    return d['online'], d['offline']


def run_validate(zpath, zeta_zeros):
    res = []
    on, off = load_dh_zeros(zpath)
    # (0) the generalised builder against wl_weil.build for zeta even (same entries)
    for X, N in ((13, 8), (200, 12)):
        D = 60
        A1, g1, _ = W.build(X, N, D)
        A2, g2, _ = build_gen(X, N, D, dict(kind='zeta', parity='even'))
        ok = all(A1[j][k].overlaps(A2[j][k]) for j in range(N + 1) for k in range(N + 1))
        dm = max(abs(A1[j][k].mid() - A2[j][k].mid()) for j in range(N + 1) for k in range(N + 1))
        row = dict(test='builder_vs_wl_weil', X=X, N=N, all_overlap=bool(ok), max_mid_diff_log10=W.tof(dm.log() / arb(10).log()) if dm > 0 else None)
        print(json.dumps(row), flush=True)
        res.append(row)
    zz = [l.split()[-1] for l in open(zeta_zeros) if l.strip()]
    for kind, parity in (('dh', 'even'), ('dh', 'odd'), ('zeta', 'odd')):
        for X in (13, 50, 200):
            N = 12
            D = 70
            spec = dict(kind=kind, parity=parity)
            A, g, info = build_gen(X, N, D, spec)
            ctx.prec = W.bits_for(D)
            L = g.L
            odd = parity == 'odd'
            if odd:  # f = sin x cos^3 x = sin(2x)/4 + sin(4x)/8, x = pi u/2L, modes psi_1, psi_2
                b = [L.sqrt() / 4, L.sqrt() / 8] + [arb(0)] * (N - 2)
            else:    # f = cos^4 x = 3/8 + cos(2x)/2 + cos(4x)/8
                b = [3 * (2 * L).sqrt() / 8, L.sqrt() / 2, L.sqrt() / 8] + [arb(0)] * (N - 2)
            n = len(A)
            Wm = sum((b[j] * A[j][k] * b[k] for j in range(n) for k in range(n)), arb(0))
            if kind == 'dh':
                s_on, qs = zero_terms(g, b, odd, [arb(t) for t in on], [(arb(a_), arb(t)) for a_, t in off])
                tot = s_on + sum(qs, acb(0))
                row = dict(test='explicit_formula', kind=kind, parity=parity, X=X, N=N, D=D, W_matrix=Wm.str(16, radius=False),
                           zero_sum_all=tot.real.str(16, radius=False), zero_sum_imag=tot.imag.str(3),
                           rel_diff=W.tof(abs(Wm - tot.real) / abs(Wm)), line_zeros_sum=s_on.real.str(16, radius=False),
                           offline_quadruples=[q.real.str(8, radius=False) for q in qs],
                           rel_diff_without_offline=W.tof(abs(Wm - s_on.real) / abs(Wm)), n_line=len(on), n_off=len(off), T=200)
            else:
                s_on, _ = zero_terms(g, b, odd, [arb(z) for z in zz[:1700]], [])
                row = dict(test='explicit_formula', kind=kind, parity=parity, X=X, N=N, D=D, W_matrix=Wm.str(16, radius=False),
                           zero_sum_all=s_on.real.str(16, radius=False), rel_diff=W.tof(abs(Wm - s_on.real) / abs(Wm)), n_line=1700)
            print(json.dumps(row), flush=True)
            res.append(row)
    return res


# ----------------------------------------------------------------------------- the cell routine
_ord_cache = {}
_dh_on = None


def ordinates(kind, D, top, zfile):
    """ordinates below top at D+10 digits: zeta by mp.zetazero; dh from the zero file (62 digits)"""
    global _dh_on
    if kind == 'dh':
        if _dh_on is None:
            _dh_on = json.load(open(zfile))['online']
        return [mp.mpf(s) for s in _dh_on if mp.mpf(s) < top]
    mp.mp.dps = D + 10
    out = []
    k = 1
    while True:
        c = _ord_cache.get(k)
        if c is None or c[0] < D + 10:
            c = (D + 10, mp.zetazero(k).imag)
            _ord_cache[k] = c
        if c[1] > top:
            break
        out.append(c[1])
        k += 1
    return out


def xi_roots(g, xi, D, odd, N):
    mp.mp.dps = D + 10
    L = mp.mpf(g.L.str(D + 10, radius=False))
    lo = 1 if odd else 0
    cc = {j: (1 / mp.sqrt(2 * L) if j == 0 else 1 / mp.sqrt(L)) for j in range(lo, N + 1)}
    aa = {j: mp.pi * j / L for j in range(lo, N + 1)}
    sg = -1 if odd else 1

    def xh(r):
        return sum(xi[j - lo] * cc[j] * (L * mp.sinc((r - aa[j]) * L) + sg * L * mp.sinc((r + aa[j]) * L)) for j in range(lo, N + 1))
    top = mp.pi * N / L
    roots = []
    step = mp.mpf(1) / 40
    r = step
    f0 = xh(r)
    while r < top:
        r1 = r + step
        f1 = xh(r1)
        if f0 * f1 < 0:
            roots.append(mp.findroot(xh, (r, r1), solver='anderson', tol=mp.mpf(10) ** (-(D - 5)), maxsteps=200))
        r, f0 = r1, f1
    return roots, top


def neg_vector(Am):
    """float LDL^T without abort: x with x^T A x = D_k <= 0 for the first nonpositive pivot k, or None"""
    n = len(Am)
    Dg, Lr = [], []
    for i in range(n):
        row = []
        for j in range(i + 1):
            s = Am[i][j]
            if j:
                lj = row if i == j else Lr[j]
                s = s - sum((row[k] * lj[k] * Dg[k] for k in range(j)), arb(0))
            if i == j:
                if not (s > 0):
                    x = [arb(0)] * n
                    x[i] = arb(1)
                    for m in range(i - 1, -1, -1):
                        x[m] = -sum((Lr[q][m] * x[q] for q in range(m + 1, i)), arb(0)) - row[m] * x[i]
                        x[m] = x[m].mid()
                    return x
                Dg.append(s)
                row.append(arb(1))
            else:
                row.append((s / Dg[j]).mid())
        Lr.append(row)
    return None


def rayleigh(A, x):
    n = len(A)
    num = sum((x[j] * A[j][k] * x[k] for j in range(n) for k in range(n)), arb(0))
    den = sum((v * v for v in x), arb(0))
    return num / den


def cell(task):
    """task: dict(kind, parity, X, N, Y, cert, roots, neig, zfile, zeta_file). Returns a json dict."""
    kind, parity, X, N = task['kind'], task['parity'], task['X'], task['N']
    spec = dict(kind=kind, parity=parity, Y=task.get('Y'))
    odd = parity == 'odd'
    neig = task.get('neig', 60)
    t0 = time.time()
    out = dict(kind=kind, parity=parity, X=X, N=N, Y=task.get('Y'))
    for extra in task.get('extras', (60, 120, 240, 480, 960)):
        D = int(3.5 * N + extra)
        A, g, info = build_gen(X, N, D, spec)
        n = len(A)
        out.update(D=D, **info)
        Amid = [[x.mid() for x in row] for row in A]
        xi = None
        mineig = None
        if n <= neig:
            mp.mp.dps = D + 10
            Am = mp.matrix(n, n)
            for i in range(n):
                for j in range(n):
                    Am[i, j] = mp.mpf(Amid[i][j].str(D + 10, radius=False))
            E, Q = mp.eigsy(Am)
            order = sorted(range(n), key=lambda i: E[i])
            mineig = E[order[0]]
            xi = [Q[j, order[0]] for j in range(n)]
            out['eps'] = mp.nstr(mineig, 12)
            out['second'] = mp.nstr(E[order[1]], 8)
            out['n_negative'] = sum(1 for i in range(n) if E[i] < 0)
            out['eigs_lowest'] = [mp.nstr(E[order[i]], 8) for i in range(min(3, n))]
            positive = mineig > 0
            vec = [arb(mp.nstr(v, D + 10)) for v in xi]
        else:
            f = W.ldl(Amid, False)
            if f is not None:
                x, q = W.eps_upper(A)
                xi = [mp.mpf(v.str(D + 10, radius=False)) for v in x]
                out['eps'] = q.mid().str(12, radius=False)
                positive = True
                vec = None
            else:
                vec = neg_vector(Amid)
                positive = False
                q = rayleigh(A, vec)
                out['eps'] = None
                out['rayleigh_float'] = q.mid().str(8, radius=False)
        done = False
        if positive:
            r = W.ldl(A, True)
            if r is not None:
                out['certified_positive'] = True
                out['certified_negative'] = False
                x, q = W.eps_upper(A)
                up = q.upper()
                out['eps_up'] = up.str(12, radius=False)
                lo = None
                for m in (0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
                    t = up * (1 - arb(m))
                    if W.ldl(W.shifted(A, t), True) is None:
                        break
                    lo = (t, m)
                out['eps_lo'] = lo[0].mid().str(12, radius=False) if lo else None
                done = True
        else:
            q = rayleigh(A, vec)
            up = q.upper()
            if up < 0:
                out['certified_negative'] = True
                out['certified_positive'] = False
                out['rayleigh_upper'] = up.str(10, radius=False)
                out['rayleigh_log10_abs'] = W.tof((-up).log() / arb(10).log())
                done = True
        if done or not task.get('cert', True):
            break
    if not out.get('certified_positive') and not out.get('certified_negative'):
        out['certified_positive'] = False
        out['certified_negative'] = False
        out['note'] = 'sign not certified up to D=%d' % D
    if task.get('roots', True) and xi is not None and n <= 100:
        roots, top = xi_roots(g, xi, D, odd, N)
        ords = ordinates(task.get('ords', kind), D, top, task.get('zfile', 'wl_dh_zeros.json'))
        mp.mp.dps = D + 10
        errs = []
        for gm in ords:
            d = min((abs(z - gm) for z in roots), default=None)
            errs.append(float(mp.log10(d)) if d is not None and d > 0 else None)
        out.update(nroots=len(roots), top=float(top), roots=[mp.nstr(z, 30) for z in roots[:8]], n_ordinates=len(ords),
                   err_log10=errs, matched_1e_3=sum(1 for e in errs if e is not None and e < -3))
    if task.get('decomp') and xi is not None and not odd and kind == 'dh':
        out['decomp'] = decomposition(g, xi, D, task['zfile'])
    out['seconds'] = round(time.time() - t0, 1)
    return out


def decomposition(g, xi, D, zfile):
    """W(f) for the least eigenvector (norm 1) = sum over the zeros of f of ghat: the line zeros, each off line quadruple, and the remainder (zeros above the file's height)"""
    on, off = load_dh_zeros(zfile)
    b = [arb(mp.nstr(v, D + 10)) for v in xi]
    s_on, qs = zero_terms(g, b, False, [arb(t) for t in on], [(arb(a_), arb(t)) for a_, t in off])
    return dict(line_sum=s_on.real.str(10, radius=False), offline=[dict(re=a_[:12], im=t[:12], quad=q.real.str(8, radius=False)) for (a_, t), q in zip(off, qs)])


def _cell(t):
    r = cell(t)
    print(json.dumps(r), flush=True)
    return r


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['zeros', 'validate', 'scan'])
    ap.add_argument('--T', type=int, default=200)
    ap.add_argument('--tasks', default='')
    ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--zeros', default='wl_dh_zeros.json')
    ap.add_argument('--zeta-zeros', default='zeros_1700.txt')
    ap.add_argument('--out', default=None)
    ap.add_argument('--nocert', action='store_true')
    ap.add_argument('--decomp', action='store_true')
    ap.add_argument('--noroots', action='store_true')
    a = ap.parse_args()
    if a.mode == 'zeros':
        res = run_zeros(a.T, a.procs)
        print(json.dumps({k: v for k, v in res.items() if k not in ('online', 'offline')}), flush=True)
    elif a.mode == 'validate':
        res = run_validate(a.zeros, a.zeta_zeros)
    else:
        tasks = []
        for t in a.tasks.split(','):
            if not t:
                continue
            p = t.split(':')
            d = dict(kind=p[0], parity=p[1], X=int(p[2]), N=int(p[3]), zfile=a.zeros, cert=not a.nocert, decomp=a.decomp, roots=not a.noroots)
            if len(p) > 4:
                d['Y'] = int(p[4])
            tasks.append(d)
        with Pool(a.procs) as pool:
            res = list(pool.imap_unordered(_cell, tasks))
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)
