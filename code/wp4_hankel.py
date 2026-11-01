# wp4_hankel.py: work package 4 (run book Section sec:hankel; chapter ch/nc_route.tex, Proposition ev:prop:hankel): the Hankel detection test of nc_route_hankel_order_hi.py pushed to N = 2000
# ordinates and about 700 digits, in ball arithmetic. Stieltjes moments M_j = sum_{k<=N} g_k^(-2j) + tail_j, tail_j = int_T^inf t^(-2j) (1/2pi) log(t/2pi) dt (closed form, T = g_N: a positive density model for the zeros
# above g_N); RH iff det[M_{i+j+1}] > 0 and det[M_{i+j+2}] > 0 for every order n. Planted pair: c = (g_k+g_{k+1})/2 + i eta replaces g_k, g_{k+1}:
# M'_j = M_j - g_k^(-2j) - g_{k+1}^(-2j) + 2 Re c^(-2j). The least order n at which a determinant is not positive is found by the LDL^T factorisation without pivoting of the two Hankel matrices in ball arithmetic:
# the leading minors are positive iff all pivots are, and the order is certified when the pivot of that order is a ball below 0 and all earlier pivots are balls above 0. A pivot ball containing 0 is "uncertain": the
# working precision is doubled and the case rerun. Ordinates: certified balls of acb.zeta_zeros at 2400 bits (about 720 digits), their midpoints used as exact numbers (the measure is then an atomic measure with these atoms).
# The genuine measure is checked to be positive definite (all pivots certified > 0) up to order NB. What is certified: the order of the first nonpositive determinant of the model measure with the midpoint ordinates and the tail
# density; not certified: that this equals the order for the true zeros above g_N (tail model), and the ordinates enter as the midpoints of balls of radius about 1e-720.
# Usage: python3 wp4_hankel.py [--quick] [--workers 16] [--trunc DIGITS] [--out wp4_hankel.json]
import argparse, json, math, sys, time
from multiprocessing import Pool
from flint import acb, arb, ctx
ORDP = 2400
G = []           # ordinates (arb, radius 0)
GRAD = [-1e9]  # log10 of the largest ball radius
BASE = {}
def get_ordinates(N, trunc=None):
    ctx.prec = ORDP
    out = []
    for s in range(1, N + 1, 100):
        zs = acb.zeta_zeros(s, min(100, N - s + 1))
        for z in zs:
            GRAD[0] = max(GRAD[0], float(z.imag.rad().log() / arb(10).log()))
            out.append(z.imag.mid())
    if trunc:
        out = [arb(g.str(trunc, radius=False)) for g in out]
    return out
def tail(j, T, L):
    return T ** (1 - 2 * j) * (L / (2 * j - 1) + arb(1) / (2 * j - 1) ** 2) / (2 * arb.pi())
def base_moments(N, JM):
    xs = [1 / (g * g) for g in G[:N]]
    pw = list(xs); M = [None]
    T = G[N - 1]; L = (T / (2 * arb.pi())).log()
    for j in range(1, JM + 1):
        M.append(sum(pw, arb(0)) + tail(j, T, L))
        pw = [p * x for p, x in zip(pw, xs)]
    return M
def planted_moments(M, k, eta, JM):
    gk, gk1 = G[k - 1], G[k]
    c = acb((gk + gk1) / 2, arb(eta)); y = 1 / (c * c)
    xa, xb = 1 / (gk * gk), 1 / (gk1 * gk1)
    pa, pb, py = xa, xb, y
    out = [None]
    for j in range(1, JM + 1):
        out.append(M[j] - pa - pb + 2 * py.real)
        pa, pb, py = pa * xa, pb * xb, py * y
    return out
def first_fail(M, off, nmax):
    """least n with a nonpositive leading minor of [M_{i+j+off}], by LDL^T; returns (n or None, status, pivot)"""
    L = []; d = []
    for i in range(nmax):
        row = []
        for j in range(i):
            s = M[i + j + off]
            for k in range(j): s -= row[k] * L[j][k] * d[k]
            row.append(s / d[j])
        di = M[2 * i + off]
        for k in range(i): di -= row[k] * row[k] * d[k]
        if di > 0: pass
        elif di < 0: return i + 1, 'neg', di
        else: return i + 1, 'unc', di
        L.append(row); d.append(di)
    return None, 'pos', None
def run_case(a):
    k, eta, N, nmax, P = a
    ctx.prec = P
    t0 = time.time()
    M = BASE[(N, P)]
    JM = 2 * nmax + 2
    Mp = planted_moments(M, k, eta, JM)
    n1, s1, d1 = first_fail(Mp, 1, nmax)
    n2, s2, d2 = first_fail(Mp, 2, nmax if n1 is None else n1)
    cands = [x for x in ((n1, s1, d1, 'H1'), (n2, s2, d2, 'H2')) if x[0] is not None]
    if not cands:
        return dict(k=k, eta=eta, N=N, P=P, nmax=nmax, order=None, status='pos', seconds=time.time() - t0)
    n, s, d, which = min(cands, key=lambda x: x[0])
    # the other matrix must not be uncertain below n
    unc = (s1 == 'unc' and n1 is not None and n1 <= n) or (s2 == 'unc' and n2 is not None and n2 <= n)
    return dict(k=k, eta=eta, N=N, P=P, nmax=nmax, order=n, status='unc' if (s == 'unc' or unc) else 'neg', which=which, orders=(n1, n2), pivot=d.str(8), seconds=time.time() - t0)
def init_base(P, Ns, JM):
    ctx.prec = P
    for N in Ns:
        BASE[(N, P)] = base_moments(N, JM)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--trunc', type=int, default=None); ap.add_argument('--out', default='wp4_hankel.json')
    ap.add_argument('--ks', default='1,2,3,5,8,12,17,24,32,48,64,96,128'); ap.add_argument('--etas', default='0.3,0.1')
    ap.add_argument('--Ns', default='2000,220'); ap.add_argument('--P0', type=int, default=3400); ap.add_argument('--nb', type=int, default=250)
    a = ap.parse_args()
    t0 = time.time()
    Nmax = 2000 if not a.quick else 300
    G.extend(get_ordinates(Nmax, a.trunc))
    print(f"{len(G)} ordinates in {time.time()-t0:.0f} s, max ordinate ball radius 10^{GRAD[0]:.1f}, g_N = {float(G[-1]):.6f}", flush=True)
    ks = [int(x) for x in a.ks.split(',')]; etas = [float(x) for x in a.etas.split(',')]; Ns = [int(x) for x in a.Ns.split(',')]
    if a.quick: ks, Ns, a.nb = [1, 2, 3, 5, 8], [300], 40
    nmaxof = lambda k: int(2.6 * k + 16)
    JM = 2 * max(nmaxof(k) for k in ks) + 2
    JM = max(JM, 2 * a.nb + 2)
    tasks = [(k, eta, N) for N in Ns for eta in etas for k in ks if not (N != Ns[0] and (eta != etas[0] or k > 24))]
    results = {}; base_check = {}
    P = a.P0
    pending = tasks
    while pending and P <= 30000:
        init_base(P, Ns, JM)
        print(f"P = {P} bits: base moments {time.time()-t0:.0f} s", flush=True)
        ctx.prec = P
        # baseline positivity of the genuine measure
        for N in Ns:
            b1 = first_fail(BASE[(N, P)], 1, a.nb); b2 = first_fail(BASE[(N, P)], 2, a.nb)
            base_check[(N, P)] = dict(H1=(b1[0], b1[1]), H2=(b2[0], b2[1]))
            print(f"  baseline N={N} P={P}: order<= {a.nb}: H1 {b1[:2]} H2 {b2[:2]}", flush=True)
        with Pool(a.workers) as p:
            res = p.map(run_case, [(k, eta, N, nmaxof(k), P) for (k, eta, N) in pending], chunksize=1)
        nxt = []
        for t, r in zip(pending, res):
            print(json.dumps(r), flush=True)
            if r['status'] == 'unc': nxt.append(t)
            else: results[t] = r
        pending = nxt; P *= 2
        BASE.clear()
    rows = [dict(r) for r in results.values()]
    fits = {}
    try:
        import numpy as np
        for N in Ns:
            for eta in etas:
                pts = sorted((r['k'], r['order']) for r in rows if r['N'] == N and r['eta'] == eta and r['order'])
                if len(pts) >= 3:
                    kk = np.array([p[0] for p in pts], float); nn = np.array([p[1] for p in pts], float)
                    A = np.vstack([kk, np.ones_like(kk)]).T; sl, ic = np.linalg.lstsq(A, nn, rcond=None)[0]
                    pw = np.polyfit(np.log(kk), np.log(nn), 1)
                    A2 = np.vstack([kk * kk, kk, np.ones_like(kk)]).T; q = np.linalg.lstsq(A2, nn, rcond=None)[0]
                    fits[f"N={N} eta={eta}"] = dict(points=pts, linear=[float(sl), float(ic)], rms_linear=float(np.sqrt(np.mean((A @ [sl, ic] - nn) ** 2))), powerlaw=[float(pw[0]), float(math.exp(pw[1]))], quadratic=[float(x) for x in q],
                                                   ratio_n_over_k=[float(n / k) for k, n in pts])
    except Exception as e:
        fits['error'] = str(e)
    json.dump(dict(max_ordinate_ball_radius=GRAD[0], trunc=a.trunc, baseline={f"{k[0]},{k[1]}": v for k, v in base_check.items()}, rows=rows, fits=fits, unresolved=[list(t) for t in pending]), open(a.out, 'w'), indent=0)
    print(json.dumps(fits)); print(f"done {time.time()-t0:.0f} s")
