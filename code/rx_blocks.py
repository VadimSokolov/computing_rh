# rx_blocks.py: the numerical claims of the one prime blocks (Remark ht:rem:blocks, R1) and of the pole needing the gamma factor (Remark pol:rem:pole, R3); Sections 3.1 and 3.3 of
# misc/reviews/nick-b4/reconcile-N.md, task E5; source: the reciprocal Xi note, cited N:n. Added for the book in October 2026, written from the mathematics (the collaborator's scripts were
# not read, run or copied).
# Claims. For a prime p, L = log p, q = p^(-alpha), vartheta = tau L, w = sqrt(2u): (1) the Mobius density (1/pi) sum_k L q^k (1 - cos k vartheta) = (1/pi) L q (1+q)(1 - cos vartheta)/((1-q)(1 - 2q cos vartheta + q^2))
# = (L/2pi) coth(alpha L/2) - sum_k alpha/(pi (alpha^2 + (tau - 2 pi k/L)^2)); the Liouville density (1/pi) sum_k (-1)^(k+1) L q^k (1 - cos k vartheta) = (1/pi) L q (1-q)(1 - cos vartheta)/((1+q)(1 + 2q cos vartheta + q^2)) >= 0,
# the Thorin density of B+_{p,alpha}(u) = (1+q)/(1 + q p^(-w)) exp(-w L q/(1+q)) (so -log B+ = int_0^inf log(1 + w^2/tau^2) rho+(tau) d tau); the mass of one Liouville bump; B+ -> 1/cosh(L sqrt(u/2)) as alpha -> 0.
# (2) Phi^lambda_2(1/2) = zeta(2 alpha + 2w) zeta(alpha)/(zeta(alpha + w) zeta(2 alpha)) exp(+w G_lambda(alpha)) and Phi^mu_2(1/2) = zeta(alpha)/zeta(alpha + w) exp(+w zeta'/zeta(alpha)) at alpha = 2, u = 1/2 (w = 1), each the
# product over primes of B+ resp. B-, with B-_{p,alpha}(u) = (1 - q p^(-w))/(1 - q) exp(-w L q/(1-q)) (note: the Mobius factor has the zeta factor at alpha + w in the numerator, since sum mu(n) n^(-s) = 1/zeta(s)); also with the sign of the exponent as printed in the source (-w G). (3) For f(s) = (s-1) zeta(s): (log f)'(alpha) = zeta'/zeta(alpha) + 1/(alpha-1) > 1/(2 f(alpha)) (the proof of R3),
# f(alpha) = E h(G_(alpha-1)) with h(x) = x/(1 - e^-x), and the function D^f_alpha(r) = Re f'/f(alpha + i r) - f'/f(alpha) is negative somewhere with mean over [0, R] tending to -f'/f(alpha).
# Conventions. zeta by mpmath (zeta and its derivative); p^(-w) with w = 1; G_lambda(s) = zeta'/zeta(s) - 2 zeta'/zeta(2s); G_mu = zeta'/zeta; the Euler product over primes p <= 1e7 in double precision (math.fsum),
# its tail beyond X estimated as -w/X (the sum of the blocks' leading term -w L p^-2 over primes > X, by the prime number theorem, error O(1/(X log X)), logged as the expected difference).
# Method and error control. Closed forms against their series (mp.nsum, 30 digits) and against the independent forms (L/pi)[q/(1 -+ q) - Re z/(1 -+ z)], z = q e^(-i vartheta) (from g'/g); the Thorin integral against -log B+ by quadrature over
# 300 periods plus the mean-density tail (error about 1e-7, logged); every mp quantity at two precisions (30 and 40 digits) with the largest difference logged.
# Output: data/rx_blocks.json (numbers as strings) and the log rx_blocks.log next to it.
# Usage: python3 rx_blocks.py [--quick] [--out FILE]; --quick uses 20 and 26 digits, primes to 1e6, 60 grid points and writes to /tmp unless --out is given; the full case takes about 20 seconds; on the cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g rx_blocks.json -g rx_blocks.log code/rx_blocks.py -- --out rx_blocks.json
import argparse
import json
import math
import os
import time
from multiprocessing import Pool

import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def dens_mob_closed(L, q, th):
    return L * q * (1 + q) * (1 - mp.cos(th)) / (mp.pi * (1 - q) * (1 - 2 * q * mp.cos(th) + q * q))


def dens_liou_closed(L, q, th):
    return L * q * (1 - q) * (1 - mp.cos(th)) / (mp.pi * (1 + q) * (1 + 2 * q * mp.cos(th) + q * q))


def dens_series(L, q, th, sgn):
    return mp.nsum(lambda k: (sgn ** (int(k) + 1)) * L * q ** k * (1 - mp.cos(k * th)), [1, mp.inf]) / mp.pi


def dens_re(L, q, th, sgn):
    z = q * mp.expjpi(-th / mp.pi)
    return L / mp.pi * (q / (1 + sgn * q) - mp.re(z / (1 + sgn * z)))


def block_checks(dps):
    mp.mp.dps = dps
    out = {}
    for (p, al, tau) in ((2, "1.5", "3"), (3, "0.05", "1")):
        L = mp.log(p)
        q = mp.mpf(p) ** (-mp.mpf(al))
        th = mp.mpf(tau) * L
        # Mobius density: series needs sign +1 in (-1)^(k+1) -> use direct sum
        smob = mp.nsum(lambda k: L * q ** k * (1 - mp.cos(k * th)), [1, mp.inf]) / mp.pi
        sli = dens_series(L, q, th, -1)
        cm, cl = dens_mob_closed(L, q, th), dens_liou_closed(L, q, th)
        rm = dens_re(L, q, th, -1) if False else L / mp.pi * (q / (1 - q) - mp.re((q * mp.expjpi(-th / mp.pi)) / (1 - q * mp.expjpi(-th / mp.pi))))
        rl = dens_re(L, q, th, +1)
        # Poisson form of the Mobius density: constant minus a row of Cauchy kernels
        a = mp.mpf(al)
        cauchy = mp.nsum(lambda k: a / (mp.pi * (a * a + (mp.mpf(tau) - 2 * mp.pi * k / L) ** 2)), [-mp.inf, mp.inf])
        pois = L / (2 * mp.pi) * mp.coth(a * L / 2) - cauchy
        out["%d_%s" % (p, al)] = dict(mob_closed=cm, mob_series=smob, mob_re=rm, mob_poisson=pois, liou_closed=cl, liou_series=sli, liou_re=rl)
    # minimum of the Liouville density over a grid of vartheta in [0, 2 pi] at p = 3, alpha = 0.05 (nonnegativity)
    L, q = mp.log(3), mp.mpf(3) ** (-mp.mpf("0.05"))
    out["liou_min"] = min(dens_liou_closed(L, q, 2 * mp.pi * j / 400) for j in range(401))
    # -log B+ by series and by the Thorin integral, p = 2, alpha = 3/2, u = 1/2 (w = 1)
    L = mp.log(2)
    al = mp.mpf("1.5")
    q = mp.mpf(2) ** (-al)
    w = mp.mpf(1)
    mlog_closed = -(mp.log(1 + q) - mp.log(1 + q * mp.mpf(2) ** (-w)) - w * L * q / (1 + q))
    mlog_series = mp.nsum(lambda k: (-1) ** (int(k) + 1) * q ** k / k * (mp.exp(-w * k * L) - 1 + w * k * L), [1, mp.inf])
    rho = lambda tau: dens_liou_closed(L, q, tau * L)
    per = 2 * mp.pi / L
    Np = 300
    integ = mp.fsum(mp.quad(lambda tau: mp.log(1 + w * w / tau ** 2) * rho(tau), [per * (j + mp.mpf(1) / 2) - per / 2 + (0 if j else 0), per * (j + 1) - per / 2 + per / 2 * 0] if False else [per * j, per * (j + mp.mpf(1) / 2), per * (j + 1)]) for j in range(Np))
    mean = L / (2 * mp.pi) * 2 * q / (1 + q)
    tailint = mean * mp.quad(lambda tau: mp.log(1 + w * w / tau ** 2), [per * Np, mp.inf])
    out["thorin"] = dict(minus_log_series=mlog_series, minus_log_closed=-mlog_closed * -1 if False else -(mp.log(1 + q) - mp.log(1 + q * mp.mpf(2) ** (-w)) - w * L * q / (1 + q)) * -1, thorin_integral=integ + tailint)
    out["thorin"]["minus_log_closed"] = -(mp.log(1 + q) - mp.log(1 + q * mp.mpf(2) ** (-w)) - w * L * q / (1 + q))
    # alpha -> 0 limit of B+ at p = 2, u = 1/2
    al0 = mp.mpf("1e-12")
    q0 = mp.mpf(2) ** (-al0)
    B0 = (1 + q0) / (1 + q0 * mp.mpf(2) ** (-w)) * mp.exp(-w * L * q0 / (1 + q0))
    out["B0"] = (B0, 1 / mp.cosh(L * mp.sqrt(mp.mpf(1) / 4)))
    # masses of one Liouville bump (p = 3), full period and windows around tau_k = pi/L
    L = mp.log(3)
    masses = {}
    for al in ("0.05", "0.01"):
        q = mp.mpf(3) ** (-mp.mpf(al))
        f = lambda th: q * (1 - q) * (1 - mp.cos(th)) / (mp.pi * (1 + q) * (1 + 2 * q * mp.cos(th) + q * q))
        full = mp.quad(f, [0, mp.pi / 2, mp.pi, 3 * mp.pi / 2, 2 * mp.pi])
        wins = {}
        for name, h in (("pi/2", mp.pi / 2), ("pi/4", mp.pi / 4), ("pi/8", mp.pi / 8)):
            wins[name] = mp.quad(f, [mp.pi - h, mp.pi - h / 2, mp.pi, mp.pi + h / 2, mp.pi + h])
        masses[al] = dict(full=full, q2=2 * q / (1 + q), wins=wins)
    out["masses"] = masses
    # zeta values at alpha = 2, u = 1/2
    a, w = mp.mpf(2), mp.mpf(1)
    r = lambda s: mp.zeta(s, derivative=1) / mp.zeta(s)
    G_l = r(a) - 2 * r(2 * a)
    Pl = mp.zeta(2 * a + 2 * w) * mp.zeta(a) / (mp.zeta(a + w) * mp.zeta(2 * a)) * mp.exp(w * G_l)
    Pm = mp.zeta(a) / mp.zeta(a + w) * mp.exp(w * r(a))
    Pl_pr = mp.zeta(2 * a + 2 * w) * mp.zeta(a) / (mp.zeta(a + w) * mp.zeta(2 * a)) * mp.exp(-w * G_l)
    Pm_pr = mp.zeta(a) / mp.zeta(a + w) * mp.exp(-w * r(a))
    out["phi"] = dict(lam=Pl, mu=Pm, lam_printed_sign=Pl_pr, mu_printed_sign=Pm_pr, G_lambda=G_l, zetap_over_zeta=r(a))
    # (log f)' and 1/(2 f), f(alpha) = E h(G_{alpha-1})
    lf = {}
    for al in ("1.01", "1.2", "2", "3", "6", "20"):
        x = mp.mpf(al)
        d = r(x) + 1 / (x - 1)
        f = (x - 1) * mp.zeta(x)
        hh = lambda t: t / (-mp.expm1(-t))
        if al == "1.01":
            Eh = None  # shape 0.01: the integrable singularity t^(-0.99) defeats the quadrature
        elif x < 3:
            # v = t^(x-1) removes the singularity t^(x-2) dt = dv/(x-1)
            k = 1 / (x - 1)
            Eh = mp.quad(lambda v: hh(v ** k) * mp.exp(-(v ** k)), [0, 1, 4, 16, 100, mp.inf]) / ((x - 1) * mp.gamma(x - 1))
        else:
            Eh = mp.quad(lambda t: hh(t) * t ** (x - 2) * mp.exp(-t), [0, 1, 10, 100, mp.inf]) / mp.gamma(x - 1)
        lf[al] = dict(dlogf=d, half_over_f=1 / (2 * f), f=f, Eh=Eh)
    out["logf"] = lf
    return out


def prime_euler(X, w=1.0, alpha=2.0):
    isp = np.ones(X + 1, dtype=bool)
    isp[:2] = False
    for k in range(2, int(X ** 0.5) + 1):
        if isp[k]:
            isp[k * k::k] = False
    p = np.nonzero(isp)[0].astype(float)
    L = np.log(p)
    q = p ** (-alpha)
    pw = p ** (-w)
    lp = np.log1p(q) - np.log1p(q * pw) - w * L * q / (1 + q)
    lm = np.log1p(-q * pw) - np.log1p(-q) - w * L * q / (1 - q)
    return math.fsum(lp.tolist()), math.fsum(lm.tolist())


def Dgrid(args):
    lo, hi, step, dps = args
    mp.mp.dps = dps
    a = mp.mpf(2)
    r = lambda s: mp.zeta(s, derivative=1) / mp.zeta(s)
    b = r(a) + 1 / (a - 1)
    vals = []
    n = int(round((hi - lo) / step))
    for i in range(n + 1):
        x = lo + i * step
        s = a + 1j * x
        vals.append(float(mp.re(r(s) + 1 / (s - 1) - b)))
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    d1, d2, X = (20, 26, 10 ** 6) if a.quick else (30, 40, 10 ** 7)
    out = a.out or ("/tmp/rx_blocks.json" if a.quick else os.path.join(HERE, "..", "data", "rx_blocks.json"))
    log = open(os.path.splitext(out)[0] + ".log", "w")

    def P(*v):
        line = " ".join(str(x) for x in v)
        print(line)
        log.write(line + "\n")
        log.flush()

    t0 = time.time()
    R = 100 if not a.quick else 40
    nchunk = 8
    chunks = [(R * i / nchunk, R * (i + 1) / nchunk - 0.05, 0.05, 15) for i in range(nchunk)]
    with Pool(8) as pool:
        ra = pool.apply_async(block_checks, (d1,))
        rb = pool.apply_async(block_checks, (d2,))
        grid = pool.map(Dgrid, chunks)
        r1, r2 = ra.get(), rb.get()
    mp.mp.dps = d2
    res = {}
    P("rx_blocks.py: R1 and R3; dps", d1, d2, "primes to", X, "quick" if a.quick else "full")
    for key in ("2_1.5", "3_0.05"):
        o, o1 = r2[key], r1[key]
        p, al = key.split("_")
        P("one prime p=%s alpha=%s (tau = 3 resp. 1):" % (p, al))
        P("  Mobius: closed %s series %s Re-form %s Poisson form %s" % tuple(mp.nstr(o[k], 18) for k in ("mob_closed", "mob_series", "mob_re", "mob_poisson")))
        P("    max |diff| to closed: %s  precision diff %s" % (mp.nstr(max(abs(o[k] - o["mob_closed"]) for k in ("mob_series", "mob_re", "mob_poisson")), 3), mp.nstr(abs(o["mob_closed"] - o1["mob_closed"]), 3)))
        P("  Liouville: closed %s series %s Re-form %s" % tuple(mp.nstr(o[k], 18) for k in ("liou_closed", "liou_series", "liou_re")))
        P("    max |diff| to closed: %s" % mp.nstr(max(abs(o[k] - o["liou_closed"]) for k in ("liou_series", "liou_re")), 3))
        res[key] = {k: mp.nstr(o[k], 20) for k in o}
    P("min of the Liouville density over a grid of 401 points, p=3 alpha=0.05:", mp.nstr(r2["liou_min"], 6), "(must be >= 0)")
    th = r2["thorin"]
    P("-log B+_{2,3/2}(1/2): closed %s series %s Thorin integral (300 periods + mean tail) %s |diff| %s" % (mp.nstr(th["minus_log_closed"], 15), mp.nstr(th["minus_log_series"], 15), mp.nstr(th["thorin_integral"], 15), mp.nstr(abs(th["thorin_integral"] - th["minus_log_closed"]), 3)))
    P("B+_{2,alpha}(1/2) at alpha = 1e-12:", mp.nstr(r2["B0"][0], 15), "vs 1/cosh(L/2) =", mp.nstr(r2["B0"][1], 15))
    for al in ("0.05", "0.01"):
        m = r2["masses"][al]
        P("mass of one Liouville bump p=3 alpha=%s: full period %s (= 2q/(1+q) = %s); windows around tau_k half-width pi/(2L), pi/(4L), pi/(8L): %s" % (al, mp.nstr(m["full"], 8), mp.nstr(m["q2"], 8), ", ".join(mp.nstr(m["wins"][k], 6) for k in ("pi/2", "pi/4", "pi/8"))))
    ph = r2["phi"]
    P("Phi^lambda_2(1/2) = %s (-log = %s); Phi^mu_2(1/2) = %s; printed sign (-wG): %s and %s; precision diff %s" % (mp.nstr(ph["lam"], 12), mp.nstr(-mp.log(ph["lam"]), 12), mp.nstr(ph["mu"], 12), mp.nstr(ph["lam_printed_sign"], 9), mp.nstr(ph["mu_printed_sign"], 9), mp.nstr(abs(ph["lam"] - r1["phi"]["lam"]), 3)))
    lp, lm = prime_euler(X)
    tail = -1.0 / X
    P("Euler products over primes <= %d: Phi^lambda = %.12f, Phi^mu = %.12f; with the tail correction -w/X in the logarithm: %.12f, %.12f; closed forms %s, %s" % (X, math.exp(lp), math.exp(lm), math.exp(lp + tail), math.exp(lm + tail), mp.nstr(ph["lam"], 12), mp.nstr(ph["mu"], 12)))
    P("   differences of the log (with tail correction): lambda %.3e, mu %.3e (expected size 1/(X log X) = %.1e)" % (lp + tail - float(mp.log(ph["lam"])), lm + tail - float(mp.log(ph["mu"])), 1 / (X * math.log(X))))
    res["phi"] = {k: mp.nstr(ph[k], 15) for k in ph}
    res["euler_products"] = dict(X=X, lam=lp, mu=lm)
    P("(log f)'(alpha), 1/(2 f(alpha)), f(alpha), E h(G_(alpha-1)):")
    res["logf"] = {}
    for al, v in r2["logf"].items():
        P("  alpha=%s: %s  %s  %s  %s (|f - E h| = %s)  above: %s" % (al, mp.nstr(v["dlogf"], 8), mp.nstr(v["half_over_f"], 8), mp.nstr(v["f"], 10), mp.nstr(v["Eh"], 10) if v["Eh"] is not None else "n/a", mp.nstr(abs(v["f"] - v["Eh"]), 3) if v["Eh"] is not None else "n/a", v["dlogf"] > v["half_over_f"]))
        res["logf"][al] = dict(dlogf=mp.nstr(v["dlogf"], 10), half=mp.nstr(v["half_over_f"], 10))
    D = np.concatenate([np.array(g) for g in grid])
    n = len(D)
    step = 0.05
    xs = np.concatenate([np.arange(len(g)) * step + R * i / nchunk for i, g in enumerate(grid)])
    b2 = float(r2["logf"]["2"]["dlogf"])
    i0 = int(np.argmin(D))
    mean = float(np.trapezoid(D, xs) / (xs[-1] - xs[0])) if hasattr(np, "trapezoid") else float(np.trapz(D, xs) / (xs[-1] - xs[0]))
    P("D^f_2(r) on r in [0, %g] (step %g): min %.6f at r = %.2f; negative on %d of %d grid points; mean over [0, %g] = %.5f vs -f'/f(2) = %.5f (the mean converges like 1/R)" % (R, step, D[i0], xs[i0], int((D < 0).sum()), n, R, mean, -b2))
    res["Df"] = dict(min=float(D[i0]), argmin=float(xs[i0]), mean=mean, target=-b2, R=R)
    P("total seconds", round(time.time() - t0, 1))
    json.dump(res, open(out, "w"), indent=1)
    log.close()


if __name__ == "__main__":
    main()
