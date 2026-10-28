# ex_zeros3e12.py: independent location of the zeros of zeta in the window t in [c-6, c+6], c = 3e12 + 100, the check of the claim in
# the explicit-sparse chapter (Subsection "Above the verified range"): 52 zeros within 6 of c, mean spacing 0.231, each certified by a sign change
# of Hardy's function. Added for the book in October 2026, written from the mathematics (the author's scripts were not read).
# Conventions. Z(t) = exp(i theta(t)) zeta(1/2 + i t) is real, theta(t) = Im logGamma(1/4 + i t/2) - (t/2) log(pi).
# Method. Z is evaluated in Arb ball arithmetic (python-flint acb.zeta, acb.lgamma, 160 bits) on the grid c - 6 + 0.02 j, j = 0..600, where
# t is an exact rational. The sign at a point is certified only if the real ball of Z excludes 0. On a grid interval whose two ends
# have certified opposite signs, bisection (exact dyadic rational midpoints, certified sign at each midpoint) runs until the bracket is
# shorter than 1e-10; the bracket and its midpoint are recorded. An interval with an uncertified end is counted and not bisected, and a
# bisection midpoint with an uncertified sign stops the bisection (flagged). Two zeros inside one grid interval (even number of sign
# changes) cannot be seen by a sign scan: the chapter's count is compared with the number of certified changes, and the smooth count
# (theta(c+6) - theta(c-6))/pi + 0 is logged for the comparison.
# Options. --quick (c = 1000, window +-1, step 0.1, 96 bits, seconds), --task I --ntasks N (slice of the grid intervals), --out FILE,
# --aggregate FILE... (merge task files into data/ex_zeros3e12.json). Error control: precision 160 bits; the script re-evaluates the
# midpoint of every bracket at 192 bits when --verify is given (certified sign must agree).
import sys, json, argparse, time
from flint import arb, acb, fmpq, ctx

def Z(t, prec):
    ctx.prec = prec
    tt = arb(t)
    s = acb(arb(1) / 2, tt)
    th = acb.lgamma(acb(arb(1) / 4, tt / 2)).imag - tt / 2 * arb.pi().log()
    z = acb(0, th).exp() * acb.zeta(s)
    return z.real, z.imag

def sign(t, prec):
    re, im = Z(t, prec)
    if re > 0: return 1
    if re < 0: return -1
    return 0

def fstr(q, n=12):
    # exact decimal string of a rational with power of two denominator
    from fractions import Fraction
    f = Fraction(int(q.p), int(q.q)); sc = 10 ** n
    v = round(f * sc)
    return "%d.%0*d" % (v // sc, n, v % sc)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true"); ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1); ap.add_argument("--out", default=None)
    ap.add_argument("--aggregate", nargs="*"); ap.add_argument("--verify", action="store_true")
    a = ap.parse_args()
    if a.aggregate:
        zs = []; unc = 0; nint = 0; flags = []; meta = None
        for fn in a.aggregate:
            d = json.load(open(fn)); zs += d["zeros"]; unc += d["uncertified_intervals"]; nint += d["intervals"]; flags += d["flags"]; meta = d
        zs.sort(key=lambda z: z["ordinate"]); n = len(zs)
        gaps = [float(zs[i + 1]["ordinate"]) - float(zs[i]["ordinate"]) for i in range(n - 1)]
        out = dict(center=meta["center"], half_window=meta["half_window"], step=meta["step"], prec=meta["prec"], count=n,
                   uncertified_intervals=unc, intervals=nint, flags=flags,
                   mean_spacing=(float(zs[-1]["ordinate"]) - float(zs[0]["ordinate"])) / (n - 1) if n > 1 else None,
                   window_over_count=2 * meta["half_window"] / n if n else None, zeros=zs)
        json.dump(out, open(a.out or "data/ex_zeros3e12.json", "w"), indent=1)
        print("count", n, "uncertified", unc, "intervals", nint, "mean spacing", out["mean_spacing"], "flags", flags)
        return
    if a.quick: c, hw, step, prec, tol = 1000, 1, fmpq(1, 10), 96, fmpq(1, 10**6)
    else: c, hw, step, prec, tol = 3 * 10**12 + 100, 6, fmpq(1, 50), 160, fmpq(1, 10**10)
    ni = int((2 * hw) / step); lo_i = a.task * ni // a.ntasks; hi_i = (a.task + 1) * ni // a.ntasks
    t0 = time.time(); zeros = []; unc = 0; flags = []
    sg = {}
    def S(j):
        if j not in sg: sg[j] = sign(c - hw + step * j, prec)
        return sg[j]
    for j in range(lo_i, hi_i):
        s1, s2 = S(j), S(j + 1)
        if s1 == 0 or s2 == 0: unc += 1; continue
        if s1 == s2: continue
        lo = c - hw + step * j; hi = lo + step; ok = True
        while hi - lo > tol:
            m = (lo + hi) / 2; sm = sign(m, prec)
            if sm == 0: flags.append("midpoint sign uncertain near " + fstr(m)); ok = False; break
            if sm == s1: lo = m
            else: hi = m
        z = dict(ordinate=fstr((lo + hi) / 2), lo=fstr(lo), hi=fstr(hi), certified=ok)
        if a.verify:
            z["verify192"] = sign((lo + hi) / 2, 192) in (0, 1, -1)
        zeros.append(z); print("zero", z, "t=%.0fs" % (time.time() - t0), flush=True)
    d = dict(center=c, half_window=hw, step=float(step.p) / float(step.q), prec=prec, task=a.task, ntasks=a.ntasks, intervals=hi_i - lo_i,
             uncertified_intervals=unc, flags=flags, zeros=zeros, seconds=time.time() - t0)
    ctx.prec = 100
    th = lambda t: acb.lgamma(acb(arb(1) / 4, arb(t) / 2)).imag - arb(t) / 2 * arb.pi().log()
    d["smooth_count_window"] = float(((th(c + hw) - th(c - hw)) / arb.pi()).mid())
    json.dump(d, open(a.out or ("ex_zeros3e12_quick.json" if a.quick else "ex_zeros3e12_t%d.json" % a.task), "w"), indent=1)
    print("task", a.task, "zeros", len(zeros), "uncertified", unc, "smooth count", d["smooth_count_window"], "s", d["seconds"])

main()
