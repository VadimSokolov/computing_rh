# rx_tail.py: the numerical claims of the exact heavy tail of the clock T_alpha (Section 3.5 of misc/reviews/nick-b4/reconcile-N.md, task E3; source: the
# reciprocal Xi note, cited N:n; chapter: ch/basepoint.tex after line 341). Added for the book in October 2026, written from the mathematics (the
# collaborator's scripts were not read, run or copied).
# Claim. With v_a(theta) = Re (xi'/xi)(a + i theta), b_a = v_a(0) and c(a) = -(1/2)(log xi)'''(a) (so v_a(theta) - b_a = c(a) theta^2 + O(theta^4)),
# t times the Levy density of T_a is pi^(-1) int_0^inf e^(-t theta^2)(v_a(theta) - b_a) d theta ~ c(a) t^(-3/2)/(4 sqrt(pi)) as t -> inf (Watson's lemma, v_a
# grows logarithmically). The check computes c(2) and the ratio of the left side to c(2)/(4 sqrt(pi) t^(3/2)) at t = 25, 100, 400 (expected 1 + O(1/t)).
# Conventions. xi(s) = s (s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2, (xi'/xi)(s) = 1/s + 1/(s-1) - (1/2) log pi + (1/2) psi(s/2) + zeta'/zeta(s).
# Method. c(a) from the closed form (log xi)'''(a) = 2/a^3 + 2/(a-1)^3 + psi''(a/2)/8 + (zeta'/zeta)''(a), (zeta'/zeta)'' = zeta'''/zeta - 3 r zeta''/zeta + 2 r^3,
# r = zeta'/zeta (mp.zeta with derivative = 1, 2, 3), checked against mp.diff of xi'/xi and against (v_a(0.01) - b_a)/10^-4 (which is c(a) (1 + O(10^-4))).
# The ratio by mp.quad of e^(-t theta^2)(v_2(theta) - b_2) on [0, Theta] with e^(-t Theta^2) = 10^-45 (v_2 grows like log theta, so the omitted part is below 1e-43).
# Error control. Everything is computed at two working precisions (20 and 30 digits) and the largest difference is logged; mp.quad reports its error estimate.
# Output: data/rx_tail.json (numbers as strings) and the log rx_tail.log next to it.
# Usage: python3 rx_tail.py [--quick] [--out FILE]; --quick uses 12 and 18 digits and writes to /tmp unless --out is given; the full case takes about 20 seconds; on the cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 30 -g rx_tail.json -g rx_tail.log code/rx_tail.py -- --out rx_tail.json
import argparse
import json
import os
import time
from multiprocessing import Pool

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def L(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2 + mp.zeta(s, derivative=1) / mp.zeta(s)


def c_closed(a):
    a = mp.mpf(a)
    z, z1, z2, z3 = mp.zeta(a), mp.zeta(a, derivative=1), mp.zeta(a, derivative=2), mp.zeta(a, derivative=3)
    r = z1 / z
    rpp = z3 / z - 3 * r * z2 / z + 2 * r ** 3
    L2 = 2 / a ** 3 + 2 / (a - 1) ** 3 + mp.psi(2, a / 2) / 8 + rpp
    return -L2 / 2


def job(args):
    dps, quick = args
    mp.mp.dps = dps
    out = {}
    for a in ("1.5", "2", "3"):
        am = mp.mpf(a)
        cc = c_closed(am)
        cd = -mp.diff(L, am, 2) / 2
        b = mp.re(L(am))
        th = mp.mpf("0.01")
        cs = (mp.re(L(am + 1j * th)) - b) / th ** 2
        out["c" + a] = (cc, cd, cs)
    b2 = mp.re(L(mp.mpf(2)))
    c2 = out["c2"][0]
    for t in (25, 100, 400):
        t = mp.mpf(t)
        Th = mp.sqrt(mp.mpf(45) * mp.log(10) / t)
        f = lambda th: mp.exp(-t * th * th) * (mp.re(L(2 + 1j * th)) - b2)
        pts = [Th * k / 8 for k in range(0, 9)]
        I, e = mp.quad(f, pts, error=True)
        lev = I / mp.pi
        ratio = lev / (c2 / (4 * mp.sqrt(mp.pi) * t ** mp.mpf("1.5")))
        out["t%d" % int(t)] = (ratio, e, lev)
    out["b2"] = b2
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    d1, d2 = (12, 18) if a.quick else (20, 30)
    out = a.out or ("/tmp/rx_tail.json" if a.quick else os.path.join(HERE, "..", "data", "rx_tail.json"))
    log = open(os.path.splitext(out)[0] + ".log", "w")

    def P(*v):
        line = " ".join(str(x) for x in v)
        print(line)
        log.write(line + "\n")
        log.flush()

    t0 = time.time()
    with Pool(2) as pool:
        r1, r2 = pool.map(job, [(d1, a.quick), (d2, a.quick)])
    mp.mp.dps = d2
    P("rx_tail.py: exact heavy tail, alpha = 2; dps", d1, d2, "quick" if a.quick else "full")
    P("b_2 = (xi'/xi)(2) =", mp.nstr(r2["b2"], 20))
    res = {}
    for al in ("1.5", "2", "3"):
        cc, cd, cs = r2["c" + al]
        P("c(%s) = -(1/2)(log xi)'''(%s): closed form" % (al, al), mp.nstr(cc, 15), " mp.diff", mp.nstr(cd, 15), " (v(0.01)-b)/1e-4", mp.nstr(cs, 15),
          " |closed - diff|", mp.nstr(abs(cc - cd), 3), " | precision diff", mp.nstr(abs(cc - r1["c" + al][0]), 3))
        res["c" + al] = mp.nstr(cc, 15)
    for tt in (25, 100, 400):
        ratio, e, lev = r2["t%d" % tt]
        ra1 = r1["t%d" % tt][0]
        P("t = %d: t * Levy density = %s, ratio to c(2)/(4 sqrt(pi) t^(3/2)) = %s (quad error %s, precision diff %s)" % (tt, mp.nstr(lev, 12), mp.nstr(ratio, 12), mp.nstr(e, 3), mp.nstr(abs(ratio - ra1), 3)))
        res["ratio%d" % tt] = mp.nstr(ratio, 12)
        res["lev%d" % tt] = mp.nstr(lev, 12)
    P("total seconds", round(time.time() - t0, 1))
    json.dump(res, open(out, "w"), indent=1)
    log.close()


if __name__ == "__main__":
    main()
