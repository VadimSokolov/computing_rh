# pg10_pairs.py: the numerical claims N10.2 of the interlacing across rungs in Chapter 22 (Section mi:sec:interlace, the example after
# Proposition mi:prop:pairs): the real zeros of A_2(t) = Re xi(2 + it) and B_{3.9}(t) = Im xi(3.9 + it) on (0, 1000], their counts on
# (0, 100], the height up to which they alternate and the first failure of alternation, and all failures to 1000; as checks (not
# printed) the pairs (A_5, A_11), (A_5, B_11) and (A_2, B_2). Added for the book in October 2026, written from the mathematics (the
# collaborator's scripts were not read, run or copied).
# Conventions. xi(s) = s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)/2; A_a(t) = Re xi(a + it), B_a(t) = Im xi(a + it) for real t (the book's A_alpha,
# B_alpha, mi:sec:paired); the zero t = 0 of the odd function B_a is not counted. Two zero sets alternate on an interval when, in the merged
# list ordered by height, no two consecutive zeros belong to the same function; a failure is a pair of consecutive zeros of the same function.
# Method. For a > 1 write xi(a + it) = |xi| e^(i phi) with phi = Im log Gamma((a+it)/2) - (t/2) log pi + arg(s(s-1)/2) + arg zeta(s), so that
# A_a and B_a vanish where cos phi and sin phi change sign. zeta(s) is evaluated in double precision by the Euler Maclaurin formula
# (N terms, K Bernoulli terms, vectorised, the table of n^(-it) shared by all abscissas), phi on the grid t_j = j h (h = 0.005); a sign change
# of cos phi or sin phi between neighbours brackets a zero (phi' >= 0.5 log(t/2 pi) is below 3, so h phi' < 0.02 and no pair of zeros of one
# function hides between grid points); every bracket is refined by bisection to 1.5e-10 (25 steps). The zeros of the pair (A_2, B_3.9) are
# then polished by Newton's method in ball arithmetic (python-flint, 128 bits; derivative from acb_series) to 25 digits, the last Newton
# step is logged. Not certified: the scan is a grid in double precision, so a closely spaced pair of zeros of one function could be
# missed. The counts of A_a and B_a are compared (they differ by at most 1 when the phase is monotone, as for a > 1) and logged.
# Error control. (1) The double precision phi is compared with the ball arithmetic phi at 12 heights (maximum difference of phi mod 2 pi logged).
# (2) On (0, 100] the whole scan is repeated with h/2 and N + 100 Euler Maclaurin terms and the zero lists must agree (maximum difference logged).
# (3) Polished against bisection zeros: maximum difference logged.
# Usage: python3 pg10_pairs.py [--tmax T] [--quick] [--out FILE] [--polish-max T] [--N N]; --quick is --tmax 100 (5 seconds, laptop); the full
# case (T = 1000, N = 400) takes 7 seconds; on the cluster: bash misc/tools/hopper_run.sh -c 2 -t 30 -g pg10_pairs.json -g pg10_pairs.log code/pg10_pairs.py -- --out pg10_pairs.json
# (T = 6000 needs --N 2000 and --polish-max 0, about 2 minutes: it shows where the failures of (A_5, B_11) begin).
import argparse
import json
import math
import sys
import time

import numpy as np
from scipy.special import bernoulli, loggamma

LOG = None


def zeta_np(sigs, t, N=400, K=10):
    """zeta(sig + i t) for each sig in sigs, by Euler Maclaurin; t a 1d array. Returns dict sig -> complex array."""
    t = np.asarray(t, dtype=float)
    n = np.arange(1, N, dtype=float)
    E = np.exp(-1j * np.outer(np.log(n), t))
    B = bernoulli(2 * K)
    out = {}
    for sig in sigs:
        w = n ** (-float(sig))
        S = w @ E
        s = sig + 1j * t
        NmS = float(N) ** (-sig) * np.exp(-1j * t * math.log(N))
        z = S + N * NmS / (s - 1) + 0.5 * NmS
        rise = s.copy()
        pw = NmS / N
        fact = 2.0
        for k in range(1, K + 1):
            z = z + B[2 * k] / fact * rise * pw
            rise = rise * (s + 2 * k - 1) * (s + 2 * k)
            pw = pw / (float(N) ** 2)
            fact = fact * (2 * k + 1) * (2 * k + 2)
        out[sig] = z
    return out


def phase_np(sig, t, N=400, K=10, chunk=3000):
    t = np.asarray(t, dtype=float)
    res = np.empty(len(t))
    for a in range(0, len(t), chunk):
        tt = t[a:a + chunk]
        z = zeta_np([sig], tt, N, K)[sig]
        s = sig + 1j * tt
        res[a:a + chunk] = np.imag(loggamma(s / 2)) - 0.5 * tt * math.log(math.pi) + np.angle(s * (s - 1) / 2) + np.angle(z)
    return res


def grid_phases(sigs, h, M, N=400, K=10, chunk=2500):
    """phi for each sig on the grid t_j = h j, j = 1..M, reusing the table n^(-i t) (built by a running product, re-anchored every chunk)."""
    n = np.arange(1, N, dtype=float)
    ln = np.log(n)
    fac = np.exp(-1j * h * ln)
    B = bernoulli(2 * K)
    res = {sig: np.empty(M) for sig in sigs}
    for a in range(0, M, chunk):
        cnt = min(chunk, M - a)
        t = h * (a + 1 + np.arange(cnt))
        E0 = np.exp(-1j * (h * (a + 1)) * ln)
        E = np.cumprod(np.concatenate([E0[:, None], np.tile(fac[:, None], (1, cnt - 1))], axis=1), axis=1)
        for sig in sigs:
            z = em_sum(sig, t, E, n, N, K, B)
            s_ = sig + 1j * t
            res[sig][a:a + cnt] = np.imag(loggamma(s_ / 2)) - 0.5 * t * math.log(math.pi) + np.angle(s_ * (s_ - 1) / 2) + np.angle(z)
    return res


def em_sum(sig, t, E, n, N, K, B):
    w = n ** (-float(sig))
    S = w @ E
    s = sig + 1j * t
    NmS = float(N) ** (-sig) * np.exp(-1j * t * math.log(N))
    z = S + N * NmS / (s - 1) + 0.5 * NmS
    rise = s.copy()
    pw = NmS / N
    fact = 2.0
    for k in range(1, K + 1):
        z = z + B[2 * k] / fact * rise * pw
        rise = rise * (s + 2 * k - 1) * (s + 2 * k)
        pw = pw / (float(N) ** 2)
        fact = fact * (2 * k + 1) * (2 * k + 2)
    return z


def sign_changes(sg):
    return np.nonzero(sg[:-1] * sg[1:] < 0)[0]


def zeros_for(sig, tgrid, kind, N, K, ph=None):
    """Zeros of A (kind 'A', cos phi) or B (kind 'B', sin phi) at abscissa sig on the grid, bisected to 25 steps."""
    if ph is None:
        ph = phase_np(sig, tgrid, N, K)
    f = np.cos(ph) if kind == "A" else np.sin(ph)
    idx = sign_changes(np.sign(f))
    lo, hi = tgrid[idx].copy(), tgrid[idx + 1].copy()
    slo = np.sign(f[idx])
    for _ in range(25):
        mid = (lo + hi) / 2
        pm = phase_np(sig, mid, N, K)
        fm = np.cos(pm) if kind == "A" else np.sin(pm)
        left = np.sign(fm) == slo
        lo = np.where(left, mid, lo)
        hi = np.where(left, hi, mid)
    return (lo + hi) / 2, ph


def merge_failures(za, zb, na="1", nb="2"):
    """na, nb: names of the two functions; a failure is two consecutive zeros of the same function."""
    ev = sorted([(float(t), na) for t in za] + [(float(t), nb) for t in zb])
    fails = []
    for i in range(len(ev) - 1):
        if ev[i][1] == ev[i + 1][1]:
            fails.append((ev[i][0], ev[i + 1][0], ev[i][1]))
    return ev, fails


def polish(sig, t0, kind, bits=128, steps=4):
    from flint import acb, acb_series, arb, ctx
    ctx.prec = bits
    t = arb(repr(float(t0)))
    last = None
    for _ in range(steps):
        s0 = acb(arb(repr(float(sig))), t)
        s = acb_series([s0, 1j], prec=2)
        xi = s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
        c = xi.coeffs()
        if kind == "A":
            f, d = c[0].real, c[1].real
        else:
            f, d = c[0].imag, c[1].imag
        step = (f / d).mid()
        t = (t - step).mid()
        last = abs(step)
    return t, last


def xi_phase_flint(sig, t, bits=128):
    from flint import acb, arb, ctx
    ctx.prec = bits
    s = acb(arb(repr(float(sig))), arb(repr(float(t))))
    xi = s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
    x = xi.mid()
    m = max(abs(x.real), abs(x.imag))   # rescale first: |xi| underflows a double for t > 900
    return math.atan2(float(x.imag / m), float(x.real / m))


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmax", type=float, default=1000.0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    ap.add_argument("--polish-max", type=float, default=1000.0)
    ap.add_argument("--step", type=float, default=0.005)
    ap.add_argument("--N", type=int, default=400, help="Euler Maclaurin terms (needs N >= tmax/pi)")
    a = ap.parse_args()
    if a.quick:
        a.tmax = 100.0
    out = a.out or ("/tmp/pg10_pairs.json" if a.quick else "data/pg10_pairs.json")
    logf = open(out.rsplit(".", 1)[0] + ".log", "w")
    def log(s):
        print(s)
        logf.write(s + "\n")
        logf.flush()
    t0 = time.time()
    h = a.step
    M = int(round(a.tmax / h))
    tgrid = h * np.arange(1, M + 1)
    N, K = a.N, 10
    assert N * math.pi >= a.tmax, "increase --N (the Euler Maclaurin formula needs |s| < 2 pi N / 2)"
    log("pg10_pairs: T = %g, step %g, Euler Maclaurin N = %d, K = %d" % (a.tmax, h, N, K))
    res = {"tmax": a.tmax, "step": h}
    # error control (1): double against ball arithmetic phase
    worst = 0.0
    for tt in np.linspace(3.1, a.tmax, 12):
        for sig in (2.0, 3.9, 5.0, 11.0):
            p1 = phase_np(sig, np.array([tt]))[0]
            p2 = xi_phase_flint(sig, tt)
            d = abs(((p1 - p2 + math.pi) % (2 * math.pi)) - math.pi)
            worst = max(worst, d)
    log("error control (1): double precision phi against ball arithmetic at 48 points (4 abscissas, 12 heights): max difference mod 2 pi %.2e" % worst)
    res["phi_error"] = worst
    sigs = {"2": 2.0, "3.9": 3.9, "5": 5.0, "11": 11.0}
    Z = {}
    PH = grid_phases(list(sigs.values()), h, M, N, K)
    for name, sig in sigs.items():
        for kind in ("A", "B"):
            z, ph = zeros_for(sig, tgrid, kind, N, K, PH[sig])
            Z[(name, kind)] = z
            res["count_%s_%s" % (kind, name)] = int(len(z))
            log("zeros of %s_%s on (0, %g]: %d (%d in (0, 100])" % (kind, name, a.tmax, len(z), int(np.sum(z <= 100))))
        nA_, nB_ = len(Z[(name, "A")]), len(Z[(name, "B")])
        log("   count check at abscissa %s: #A - #B = %d (|#A - #B| <= 1 expected; the zeros of A and B alternate where the phase is monotone)" % (name, nA_ - nB_))
    # error control (2): half step and N + 100 on (0, min(tmax, 100)]
    T2 = min(a.tmax, 100.0)
    M2 = int(round(T2 / (h / 2)))
    tg2 = (h / 2) * np.arange(1, M2 + 1)
    worst2 = 0.0
    same = True
    for name, sig in (("2", 2.0), ("3.9", 3.9)):
        for kind in ("A", "B"):
            z2, _ = zeros_for(sig, tg2, kind, N + 100, K + 2)
            z1 = Z[(name, kind)]
            z1 = z1[z1 <= T2]
            if len(z1) != len(z2):
                same = False
                log("error control (2): COUNT DIFFERS for %s_%s: %d against %d" % (kind, name, len(z1), len(z2)))
            else:
                worst2 = max(worst2, float(np.max(np.abs(z1 - z2))))
    log("error control (2): on (0, %g] with step %g and N = %d: same counts: %s, max zero difference %.2e" % (T2, h / 2, N + 100, same, worst2))
    res["control2_same_counts"] = same
    res["control2_max_diff"] = worst2
    # the pair (A_2, B_3.9)
    ev, fails = merge_failures(Z[("2", "A")], Z[("3.9", "B")], "A_2", "B_3.9")
    nA = int(np.sum(Z[("2", "A")] <= 100))
    nB = int(np.sum(Z[("3.9", "B")] <= 100))
    log("N10.2a  zeros of A_2 in (0, 100]: %d" % nA)
    log("N10.2b  zeros of B_3.9 in (0, 100]: %d" % nB)
    log("merged list starts: %s" % " ".join(lab for _, lab in ev[:20]))
    if fails:
        f0 = fails[0]
        k = [i for i, e in enumerate(ev) if e[0] == f0[0]][0]
        log("N10.2c  alternation holds on the merged list up to the zero at height %.10f (zero number %d of the merged list); the first failure is the pair of consecutive zeros of %s at %.10f and %.10f" % (f0[0], k + 1, f0[2], f0[0], f0[1]))
        log("N10.2d  first failure at heights %.10f, %.10f" % (f0[0], f0[1]))
    log("failures of (A_2, B_3.9) on (0, %g]: %d" % (a.tmax, len(fails)))
    for f in fails:
        log("   failure: two consecutive zeros of %s at %.10f and %.10f" % (f[2], f[0], f[1]))
    res["pair_A2_B39"] = {"nA100": nA, "nB100": nB, "failures": [list(f) for f in fails], "first_failure": list(fails[0]) if fails else None}
    # polishing
    pm = a.polish_max
    pol = {}
    worst3 = 0.0
    for name, kind in (("2", "A"), ("3.9", "B")):
        zz = Z[(name, kind)]
        zz = zz[zz <= pm]
        vals = []
        lastmax = 0.0
        for z0 in zz:
            tp, last = polish(sigs[name], z0, kind)
            worst3 = max(worst3, abs(float(tp) - float(z0)))
            lastmax = max(lastmax, float(last))
            vals.append(tp.str(26, radius=False))
        pol["%s_%s" % (kind, name)] = vals
        log("polished %d zeros of %s_%s up to %g: largest last Newton step %.1e" % (len(zz), kind, name, pm, lastmax))
    log("error control (3): polished against bisection zeros: max difference %.2e" % worst3)
    res["polished_zeros"] = pol
    res["polish_max_diff"] = worst3
    # first failure polished
    if fails:
        f0 = fails[0]
        nm = "2" if f0[2] == "A_2" else "3.9"
        tp1, _ = polish(sigs[nm], f0[0], f0[2][0])
        tp2, _ = polish(sigs[nm], f0[1], f0[2][0])
        log("first failure, polished: %s %s and %s" % (f0[2], tp1.str(20, radius=False), tp2.str(20, radius=False)))
        res["first_failure_polished"] = [tp1.str(26, radius=False), tp2.str(26, radius=False)]
    # checks
    for (n1, k1, n2, k2, label) in (("5", "A", "11", "A", "A_5 with A_11"), ("5", "A", "11", "B", "A_5 with B_11"), ("2", "A", "2", "B", "A_2 with B_2")):
        ev2, f2 = merge_failures(Z[(n1, k1)], Z[(n2, k2)], "%s_%s" % (k1, n1), "%s_%s" % (k2, n2))
        last = f2[-1][1] if f2 else None
        first = f2[0][0] if f2 else None
        log("check %s: %d failures on (0, %g]; first at %s, last at %s; merged list starts %s" % (label, len(f2), a.tmax, "%.6f" % first if first else "-", "%.6f" % last if last else "-", " ".join(l for _, l in ev2[:12])))
        res["check_" + label.replace(" ", "_")] = {"failures": len(f2), "first": first, "last": last, "tail_failures_above_500": int(sum(1 for f in f2 if f[0] > 500)), "list": [list(f) for f in f2[:200]]}
    # drift diagnostic for the pair (5, 11): D - 3 pi/2 = (phi_11 - phi_5) - 3 pi/2 (mod 2 pi, taken in (-pi, pi]); (iii) needs it to change sign
    dev = np.angle(np.exp(1j * (PH[11.0] - PH[5.0] - 1.5 * math.pi)))
    for lo_, hi_ in ((100.0, a.tmax), (a.tmax / 2, a.tmax), (a.tmax * 0.9, a.tmax)):
        m_ = (tgrid >= lo_) & (tgrid <= hi_)
        if m_.any():
            log("diagnostic (5, 11): D - 3pi/2 on [%g, %g]: min %.5f, max %.5f, t*(D - 3pi/2) between %.3f and %.3f; J(5,11) = sum_p arctan p^-5 - arctan p^-11 is about 0.035 (p = 2 alone 0.0308)" % (lo_, hi_, dev[m_].min(), dev[m_].max(), (tgrid[m_] * dev[m_]).min(), (tgrid[m_] * dev[m_]).max()))
    sc = int(np.sum(np.sign(dev[tgrid >= 50])[:-1] * np.sign(dev[tgrid >= 50])[1:] < 0))
    log("diagnostic (5, 11): sign changes of D - 3pi/2 on [50, %g]: %d" % (a.tmax, sc))
    res["dev_5_11_sign_changes"] = sc
    json.dump(res, open(out, "w"), indent=1)
    log("wrote %s (%.1f s)" % (out, time.time() - t0))


if __name__ == "__main__":
    main()
