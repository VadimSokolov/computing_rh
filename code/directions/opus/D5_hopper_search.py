#!/usr/bin/env python3
# Proposal for Hopper, not run there: adversarial search against Conjecture 14 of D5.md, delta=1.
# Each task fixes (n,J) from its id, hill climbs Y^2-f(J) in floating point over heights t_j and c_j<=1 with c_0=1,
# and rechecks the best point of each restart. With python-flint the recheck is ball arithmetic: excess_lo>0 is a certified counterexample.
# Without python-flint it falls back to mpmath and marks the line "certified": false.
# Usage: python3 -I D5_hopper_search.py --task 17 --seconds 14400 --out results/   (SLURM array in D5.md, Section 6)
import argparse, json, os, time
import numpy as np
try:
    from flint import acb, arb, acb_poly, ctx
    HAVE_FLINT = True
except Exception:
    import mpmath as mp
    HAVE_FLINT = False
NS = [2, 3, 4, 5, 6, 8, 10, 12]
JS = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 1.95]
def f(J): return 0.5 * (1 - 2 * J + np.sqrt(1 + 2 * J * J))
def qpoly(ts, cs):
    Q = np.array([1.0 + 0j])
    for t, c in zip(ts, cs): Q = np.convolve(Q, np.array([1, 2j * t, -t * t - c], dtype=complex))
    return Q
def zpoly(Q, J):
    d = Q.copy(); s = np.convolve(d, d); term = 1.0; k = 0
    while len(d) > 1:
        k += 1; d = np.polyder(d); term *= J / k; s = np.polyadd(s, term * np.convolve(d, d))
    return s
def y2(ts, cs, J):
    r = np.roots(zpoly(qpoly(ts, cs), J)); y = r.real.max()
    return float(y * y) if y > 0 else 0.0
def recheck(ts, cs, J, bits=256):
    """Lower bound of Y^2-f(J), from coefficients built exactly from the double precision parameters."""
    if HAVE_FLINT:
        ctx.prec = bits
        Q = acb_poly([acb(1)])
        for t, c in zip(ts, cs): Q = Q * acb_poly([-(acb(float(t)) ** 2) - acb(float(c)), acb(0, 2 * arb(float(t))), acb(1)])
        d = Q; S = d * d; term = arb(1); k = 0
        while d.degree() > 0:
            k += 1; d = d.derivative(); term = term * arb(float(J)) / k; S = S + d * d * acb(term)
        try: roots = S.roots()
        except ValueError: return None, True  # cluster not separated at this precision: inspect by hand
        y_lo = max((z.real.lower() for z in roots), key=float)
        fJ = (1 - 2 * arb(float(J)) + (1 + 2 * arb(float(J)) ** 2).sqrt()) / 2
        ex = y_lo * y_lo - fJ.upper()
        return (float(ex) if bool(y_lo > 0) else -1.0), True  # a counterexample requires ex > 0 as a ball: rerun with larger bits before believing it
    mp.mp.dps = 80
    Q = [mp.mpc(1)]
    for t, c in zip(ts, cs):
        q = [mp.mpc(1), mp.mpc(0, 2 * mp.mpf(float(t))), -(mp.mpf(float(t)) ** 2) - mp.mpf(float(c))]
        R = [mp.mpc(0)] * (len(Q) + 2)
        for i, x in enumerate(Q):
            for j, u in enumerate(q): R[i + j] += x * u
        Q = R
    d = Q; term = mp.mpf(1); k = 0; S = [mp.mpc(0)] * (2 * len(Q) - 1)
    while True:
        D = [mp.mpc(0)] * (2 * len(d) - 1)
        for i, x in enumerate(d):
            for j, u in enumerate(d): D[i + j] += x * u
        off = len(S) - len(D)
        for i, x in enumerate(D): S[off + i] += term * x
        if len(d) == 1: break
        m = len(d) - 1; d = [d[i] * (m - i) for i in range(m)]; k += 1; term = term * mp.mpf(float(J)) / k
    r = mp.polyroots(S, maxsteps=5000, extraprec=1500); y = max(mp.re(z) for z in r)
    return float(y * y - (1 - 2 * mp.mpf(float(J)) + mp.sqrt(1 + 2 * mp.mpf(float(J)) ** 2)) / 2) if y > 0 else -1.0, False
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--task", type=int, required=True); ap.add_argument("--seconds", type=float, default=60); ap.add_argument("--steps", type=int, default=400); ap.add_argument("--out", default=".")
    a = ap.parse_args(); n = NS[a.task % len(NS)]; J = JS[(a.task // len(NS)) % len(JS)]; rng = np.random.default_rng(a.task); fJ = f(J); t0 = time.time(); best = -9.0; rs = 0
    os.makedirs(a.out, exist_ok=True); path = os.path.join(a.out, f"d5_search_{a.task}.jsonl")
    with open(path, "a") as out:
        while time.time() - t0 < a.seconds:
            rs += 1; scale = rng.choice([0.03, 0.1, 0.3, 1, 3, 10])
            ts = rng.normal(0, scale, n); cs = np.minimum(1, rng.uniform(-3, 1, n)); cs[0] = 1; cur = y2(ts, cs, J) - fJ; sig = 0.5
            for _ in range(a.steps):
                t2 = ts + sig * rng.normal(0, 1, n) * rng.choice([0.03, 0.3, 1]); c2 = np.minimum(1, cs + sig * rng.normal(0, 1, n) * rng.choice([0.03, 0.3, 1])); c2[0] = 1
                v = y2(t2, c2, J) - fJ
                if v > cur: ts, cs, cur = t2, c2, v; sig = min(sig * 1.2, 1)
                else: sig = max(sig * 0.985, 1e-3)
            if cur > best - 1e-3:
                ex, cert = recheck(ts, cs, J); best = max(best, cur)
                out.write(json.dumps({"n": n, "J": J, "float_excess": cur, "excess_lo": ex, "certified": cert, "t": list(map(float, ts)), "c": list(map(float, cs))}) + "\n"); out.flush()
    print(f"task {a.task} n={n} J={J} restarts {rs} best float excess {best:.3e} flint {HAVE_FLINT} file {path}")
if __name__ == "__main__": main()
