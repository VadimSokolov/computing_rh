#!/usr/bin/env python3
"""D1: zero side certificate, Theorem Z of D1.md, Section 3.8.

Proves W_L(f) >= c |f|^2 on L^2[-L,L] from the zeros of zeta alone, in ball arithmetic:
  W_L(f) >= omega sum_{|gamma|<=T3} |F(gamma)|^2 + beta (1/2 pi) int_{|t|>T3} |F|^2 dt - tau |f|^2.
Input: Arb's rigorous zeros (acb.zeta_zeros, on the line, complete up to the last one used), the
gap condition (no prime power in [exp(2L), exp(2L+2 eta))), and the single number A_L.
No prime term is evaluated. This reaches only L <= (1/2) log X - eta, not the closed rung.

Usage:
  timeout 600 python3 -I D1_zero_cert.py auto L eta eps betamin varpi prec ng
eta: width of the multiplier a; eps: target size of the tails; betamin: floor for beta_0;
varpi: share of the zeros above t_1 given to the high frequency bound.
Importing D1_certify and D1_certify2 only for the common ball routines.
"""
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import D1_certify as C
import D1_certify2 as Z
from D1_certify import arb, acb, arb_mat, ctx, up, gl_nodes, bessel_table, record

ZEROS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "data", "zeros_8000.txt")


def log_tail_f(eta, r, R):
    Cr = eta * math.sqrt(3.0) / (8 * math.sqrt(r))
    return math.log(2 * Cr * R / (2 * r - 1)) + 2 * r * math.log(2 * r / (eta * R))


def best_r(eta, R):
    r0 = max(3, int(eta * R / (2 * math.e)))
    return min(range(max(3, r0 - 2), r0 + 3), key=lambda r: log_tail_f(eta, r, R))


def radius_for(eta, eps):
    lo_, hi = 1.0, 1e7
    for _ in range(80):
        mid = math.sqrt(lo_ * hi)
        if log_tail_f(eta, best_r(eta, mid), mid) <= math.log(eps):
            hi = mid
        else:
            lo_ = mid
    return hi, best_r(eta, hi)


def tail(eta, r, R):
    """Upper bound for int_{|u|>=R} Phi_{eta,r}(u) du, Lemma in Section 3.8."""
    Cr = eta * arb(3).sqrt() / (8 * arb(r).sqrt())
    return 2 * Cr * R / (2 * r - 1) * (2 * r / (eta * R)) ** (2 * r)


def is_prime_power(n):
    for p in range(2, n + 1):
        if n % p == 0:
            while n % p == 0:
                n //= p
            return n == 1
    return False


def plan(Lf, etaf, eps, betamin, varpi):
    gam = [float(x) for x in open(ZEROS).read().split()]
    gam = [g for g in gam if g > 1]
    Ra, ra = radius_for(etaf, eps)
    best = None
    K = len(gam)
    # largest gap above each index, up to the top zero that will be used; refined after the choice
    topidx = K - 1
    maxgap = [0.0] * K
    m = 0.0
    for k in range(topidx - 1, -1, -1):
        m = max(m, gam[k + 1] - gam[k])
        maxgap[k] = m
    for k in range(3, 800):
        d = maxgap[k]
        for i in range(2, 200):
            e2 = i * 0.01
            Lb = Lf + etaf + e2
            kap = d * Lb / math.pi
            if kap >= 1:
                break
            b0 = 2 * math.pi * varpi * (1 - kap) ** 2 / d
            if b0 < betamin:
                continue
            Rc, rc = radius_for(e2, eps)
            T3 = 8 * math.ceil((gam[k] + 2 * Rc) / 8)
            if best is None or T3 < best[0]:
                best = (T3, k, e2, Rc, rc, kap, b0, d)
    return gam, best, Ra, ra


def main():
    if sys.argv[1] != "auto":
        raise SystemExit("usage: auto L eta eps betamin varpi prec ng")
    Ls, etas = sys.argv[2], sys.argv[3]
    eps = float(sys.argv[4])
    betamin = float(sys.argv[5])
    varpis = sys.argv[6]
    prec = int(sys.argv[7])
    ng = int(sys.argv[8])
    Lf, etaf, varpif = float(Ls), float(etas), float(varpis)
    for n in range(2, 5000):
        if 2 * Lf - 1e-9 <= math.log(n) < 2 * Lf + 2 * etaf + 1e-9 and is_prime_power(n):
            raise SystemExit("gap condition fails at n=%d" % n)
    t0 = time.time()
    gam, best, Raf, ra = plan(Lf, etaf, eps, betamin, varpif)
    T3, k1, e2f, Rcf, rc, kapf, b0f, df = best
    Yaf = T3 + Raf
    top = Yaf + Raf
    ktop = next(i for i, g in enumerate(gam) if g >= top) + 1
    kz = next(i for i, g in enumerate(gam) if g > T3 + 1)
    print("plan: T3=%d t1=gamma_%d=%.3f eta2=%.2f Rc=%.1f rc=%d float kappa=%.4f beta0=%.4f delta=%.3f; Ra=%.1f ra=%d; zeros to index %d (height %.1f); %d zeros below T3" % (
        T3, k1 + 1, gam[k1], e2f, Rcf, rc, kapf, b0f, df, Raf, ra, ktop, gam[ktop - 1], kz), flush=True)
    # rigorous zeros
    ctx.prec = 200
    low = acb.zeta_zeros(1, kz)
    ctx.prec = 64
    high = []
    cur = kz + 1
    while cur <= ktop:
        num = min(500, ktop - cur + 1)
        high.extend(acb.zeta_zeros(cur, num))
        cur += num
    print("zeros by Arb: %d in %.1fs; last %s" % (kz + len(high), time.time() - t0, high[-1].imag.str(12)), flush=True)
    C.setup(prec)
    L = arb(Ls)
    eta = arb(etas)
    eta2 = arb("%.2f" % e2f)
    varpi = arb(varpis)
    gamz = [z.imag for z in low] + [z.imag for z in high]
    for z in low + high:
        assert (z.real - arb(1) / 2).contains(0)
    # sampling constants from the rigorous list
    d_up = 0.0
    for k in range(k1, len(gamz) - 1):
        d_up = max(d_up, up(gamz[k + 1] - gamz[k]))
    delta = arb(d_up) * (1 + arb(1) / 10 ** 9)
    Lb = L + eta + eta2
    kappa = delta * Lb / C.PI
    if not kappa < 1:
        raise SystemExit("kappa >= 1")
    t1 = arb(up(gamz[k1])) + arb(1) / 10 ** 9
    tK = arb(float(gamz[-1].mid()) - 1e-6)
    Rc = arb(Rcf)
    Yc = t1 + Rc
    if not Yc + Rc < T3:
        Rc = (arb(T3) - t1) / 2 * (1 - arb(1) / 10 ** 9)
        Yc = t1 + Rc
    eps1 = tail(eta2, rc, Rc)
    eps4 = eps1
    Ra = arb(Raf)
    Ya = arb(T3) + Ra
    eps3 = 2 * tail(eta, ra, Ra)
    s0 = tK - Ya
    eps2 = tail(eta, ra, s0)
    epsi = (eps1 * eps1 + eps2 * eps2).sqrt()
    eps5 = (eta * eta / (16 * ra)).exp() * tail(eta, ra, Ya)
    # tail of the zero sum beyond H = tK
    Cr = eta * arb(3).sqrt() / (8 * arb(ra).sqrt())
    D0 = (2 * Cr / (2 * ra - 1)) ** 2 * (2 * ra / eta) ** (4 * ra)
    p = 4 * ra - 2
    H = tK
    lH = H.log()
    SE = 2 * (eta * eta / (8 * ra)).exp() * (H * lH * D0 * s0 ** (-p) + D0 * ((1 + lH) * s0 ** (1 - p) / (p - 1) + s0 ** (2 - p) / (H * (p - 1) * (p - 2))))
    comb = C.comb_terms(2 * Lf)
    AL = sum((w for (_, _, w) in comb), arb(0))
    m0 = acb(arb(1) / 4, arb(0)).digamma().real - C.PI.log()
    psi0 = C.m_of(arb(T3)) - AL
    beta0 = 2 * C.PI * varpi * (1 - kappa) ** 2 * (1 - 2 * eps4) / delta
    beta = arb(up(beta0) * (1 - 1e-9))
    if not beta < psi0:
        beta = psi0 * (1 - arb(1) / 10 ** 9)
    if not beta <= beta0:
        raise SystemExit("beta above beta0")
    tau = 4 * C.PI * varpi * (1 - kappa) * epsi / delta + eps3 * (abs(m0) + AL) + 2 * L * L.exp() * SE + eps5 * (2 + eps5) * 2 * (L + L.sinh())
    omega = (1 - varpi) * (1 - eps3)
    print("kappa=%.5f delta=%.5f L_b=%.4f beta=%.5e psi0=%.4f omega=%.5f eps1=%.2e eps2=%.2e eps3=%.2e eps5=%.2e SE=%.2e tau=%.4e" % (
        float(kappa), float(delta), float(Lb), float(beta), float(psi0), float(omega), up(eps1), up(eps2), up(eps3), up(eps5), up(SE), up(tau)), flush=True)
    Nmax = int(1.5 * Lf * T3) + 40
    SigT = C.tail_sum(L * T3, Nmax, L)
    Kz = sum(1 for g in gamz[:kz] if g < arb(T3))
    epsB = omega * 2 * Kz * (2 * L * SigT).sqrt() + beta * (T3 * SigT / C.PI).sqrt()
    epsD = beta * (T3 / C.PI) * SigT
    print("Nmax=%d prec=%d Kz=%d SigT=%.2e epsB=%.2e epsD=%.2e" % (Nmax, prec, Kz, up(SigT), up(epsB), up(epsD)), flush=True)
    # B = (1/2 pi) int_{|t|<=T3} q q^t, by the graded panels, entire integrand
    near = gl_nodes(ng)
    ngf = max(40, int(0.7 * ng))
    far = gl_nodes(ngf)
    half = arb(1) / 2
    nodes, weights = [], []
    eqB = arb(0)
    for (a, b, rho, y, _) in Z.panels(T3):
        xs, ws = near if b - a == 1 else far
        w = arb(b - a)
        cen = arb(a + b) / 2
        for x, wg in zip(xs, ws):
            nodes.append(cen + half * w * x)
            weights.append(half * w * wg)
        n = ng if b - a == 1 else ngf
        eqB += (w / 2) * 8 * (2 * L * (2 * Nmax + 1) * (2 * L * y).exp()) * rho ** (-2 * n) / (1 - 1 / rho) / C.PI
    pref = [(2 * L * (2 * n + 1)).sqrt() for n in range(Nmax + 1)]
    table = []
    for t in nodes:
        tab = bessel_table(L * t, Nmax)
        table.append([pref[n] * tab[n] for n in range(Nmax + 1)])
    nd = dict(nodes=nodes, weights=weights, table=table, mt=[arb(1)] * len(nodes), cosm=[], nn=len(nodes))
    ztab = []
    for g in gamz[:Kz]:
        # evaluate at the exact midpoint and add the radius by hand: |j_n'| <= 1 on the real line.
        # Passing the ball itself would multiply its radius by the growth of the ball recurrence.
        g0 = g.mid()
        wid = arb(0, up(L * (abs(g - g0) + g.rad())))
        tab = bessel_table(L * g0, Nmax)
        ztab.append([pref[n] * (tab[n] + wid) for n in range(Nmax + 1)])
    print("tables in %.1fs, %d nodes, eqB=%.2e" % (time.time() - t0, len(nodes), up(eqB)), flush=True)
    res = dict(tag="zero", L=Lf, eta=etaf, T3=T3, Nmax=Nmax, prec=prec, beta=float(beta), tau=up(tau), omega=float(omega), kappa=float(kappa), delta=float(delta), ktop=ktop)
    for sector in "eo":
        ts = time.time()
        idx, B, _, _ = C.sector_matrices(L, T3, Nmax, nd, [], sector, arb(0))
        sg = [(-1) ** (n // 2) for n in idx]
        nm = len(idx)
        QZ = arb_mat(nm, Kz, [sg[i] * ztab[k][idx[i]] for i in range(nm) for k in range(Kz)])
        S = 2 * (QZ * QZ.transpose())
        wide = arb(0, up(beta * eqB) * 1.01)
        base = [[omega * S[i, j] - beta * B[i, j] + wide for j in range(i + 1)] for i in range(nm)]

        def mat(c):
            gap_ = beta - tau - arb(c)
            D = gap_ - epsD
            if not D > 0:
                raise SystemExit("D not positive")
            sh = gap_ - epsB * epsB / D
            return [[(base[i][j] + sh) if i == j else base[i][j] for j in range(i + 1)] for i in range(nm)]
        A0 = mat(0)
        Lm = Z.chol_mid(A0)
        if Lm is None:
            print("  sector %s modes=%d: not positive" % (sector, nm), flush=True)
            res["c_" + sector] = None
            continue
        lam = Z.lam_estimate(A0, Lm)
        c = 0.999 * lam
        done = None
        for _ in range(12):
            ok, _ = Z.verify_pd(mat(c))
            if ok:
                done = c
                break
            c *= 0.8
        print("  sector %s modes=%d: float lowest %.5e; certified W >= c |f|^2 with c=%s (%.1fs)" % (sector, nm, lam, ("%.5e" % done) if done else "NONE", time.time() - ts), flush=True)
        res["c_" + sector] = done
    res["seconds"] = time.time() - t0
    record(res)
    print("  total %.1fs" % res["seconds"], flush=True)


if __name__ == "__main__":
    main()
