# pg10_window.py: the numerical claims N10.1 of the interlacing across rungs in Chapter 22 (Section mi:sec:interlace, Proposition
# mi:prop:pairs and its examples): I(5, 11), arctan(1/4) - arctan(2^-3.9), the endpoints in h of the cells for alpha = 2, and the abscissa
# sigma_1 with M(sigma_1) = pi/4. Added for the book in October 2026, written from the mathematics (the collaborator's scripts were not
# read, run or copied).
# Conventions. M(s) = -zeta'(s)/zeta(s) = sum_n Lambda(n) n^-s for s > 1 (decreasing); I(a, a') = int_a^a' M = log(zeta(a)/zeta(a')); h = a' - a;
# J(a, a') = sum_p (arctan p^-a - arctan p^-a'), so J <= I. With u(h) = h pi/4 + I(2, 2+h) and l(h) = h pi/4 - I(2, 2+h) the conditions of
# Proposition mi:prop:pairs read: (i)_k  l(h) > k pi and u(h) < (k+1) pi;  (ii)_k  l(h) > k pi - pi/2 and u(h) < k pi + pi/2.
# Method. All values in ball arithmetic (python-flint, arb/acb_series for zeta and zeta'). The endpoints are roots of u(h) = c pi/2 or
# l(h) = c pi/2; u and l are increasing for alpha = 2 because M(2+h) <= M(2) < pi/4 (the script certifies M(2) < pi/4), so every root is
# unique and is enclosed by interval bisection with the sign of each midpoint certified by the ball (an undecided midpoint stops the run).
# For k = 0 the lower condition l(h) > 0 (resp. l(h) > -pi/2) holds for all h > 0 because l(0) = 0 and l is increasing. sigma_1 is the unique
# root of M(s) = pi/4 on (1, oo), M being decreasing (bisection with certified signs). Not printed in the chapter: M(2), M(3) and
# J(2, 3.9) summed over the primes p < 2e5 (prime powers n <= 2e5 for M) with the tail bound sum_{n > N} <= N^-1 for J (terms <= p^-2)
# and sum_{n > N} log n n^-s <= N^(1-s)(log N/(s-1) + 1/(s-1)^2) for M, against the ball value.
# Status: the roots, I(5, 11), the arctan difference and sigma_1 are CERTIFIED (the printed digits lie inside the enclosing balls);
# the prime sums are a consistency check.
# Error control: the whole run is repeated at a second precision (the enclosures must overlap; the difference of the midpoints is logged).
# Usage: python3 pg10_window.py [--quick] [--out FILE]; --quick uses 100 and 130 bits and stops the bisection at 1e-15 (about 3 seconds);
# the full run (200 and 260 bits, 1e-35) takes about 20 seconds; on the cluster:
# bash misc/tools/hopper_run.sh -c 1 -t 20 -g pg10_window.json -g pg10_window.log code/pg10_window.py -- --out pg10_window.json
import argparse
import json
import math
import sys
import time

from flint import acb, acb_series, arb, ctx


def zeta_and_deriv(sig):
    """zeta(sig), zeta'(sig) as arb balls (sig an arb)."""
    z = acb_series([acb(sig), 1], prec=2).zeta()
    c = z.coeffs()
    return c[0].real, c[1].real


def Mfun(sig):
    z, d = zeta_and_deriv(sig)
    return -d / z


def Ifun(a, ap):
    return (arb(a).zeta() / arb(ap).zeta()).log()


def sign_of(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def bisect(f, lo, hi, width):
    """f increasing (or decreasing), certified signs at lo and hi; returns (lo, hi) with f(lo), f(hi) of opposite certified signs."""
    slo, shi = sign_of(f(lo)), sign_of(f(hi))
    if slo == 0 or shi == 0 or slo == shi:
        raise RuntimeError("bracket not certified")
    while hi - lo > width:
        mid = (lo + hi) / 2
        s = sign_of(f(mid))
        if s == 0:
            raise RuntimeError("midpoint undecided at width %s" % (hi - lo))
        if s == slo:
            lo = mid
        else:
            hi = mid
    return lo, hi


def run(bits, width, pmax, log):
    ctx.prec = bits
    pi = arb.pi()
    out = {"bits": bits}
    # certificates for monotonicity
    M2 = Mfun(arb(2))
    assert M2 < pi / 4
    out["M2_ball"] = M2.str(30, radius=False)
    log("M(2) = -zeta'(2)/zeta(2) = %s < pi/4 = %s (certified)" % (M2.str(25, radius=False), (pi / 4).str(25, radius=False)))
    # N10.1a
    I511 = Ifun(5, 11)
    assert I511 < pi / 2
    out["I_5_11"] = I511
    log("N10.1a  I(5,11) = log(zeta(5)/zeta(11)) = %s  (< pi/2: %s)" % (I511.str(30, radius=True), I511 < pi / 2))
    # N10.1b
    D = arb(1) / 4
    D = D.atan() - (arb(2) ** (arb(-39) / 10)).atan()
    out["arctan_diff"] = D
    log("N10.1b  arctan(1/4) - arctan(2^-3.9) = %s ; pi/40 = %s ; difference exceeds pi/40: %s" % (D.str(30, radius=True), (pi / 40).str(20, radius=False), D > pi / 40))
    # N10.1c..h and further cells
    def u(h): return h * pi / 4 + Ifun(2, 2 + h)
    def l(h): return h * pi / 4 - Ifun(2, 2 + h)
    def root(kind, c2):
        """root h of kind(h) = c2 * pi/2, kind in {u, l}"""
        f = (lambda h: u(h) - c2 * pi / 2) if kind == "u" else (lambda h: l(h) - c2 * pi / 2)
        return bisect(f, arb(0), arb(40), width)
    cells = {}
    for k in range(0, 4):
        # (i)_k: l > k pi, u < (k+1) pi ; (ii)_k: l > k pi - pi/2, u < k pi + pi/2
        lo_i = None if k == 0 else root("l", 2 * k)
        up_i = root("u", 2 * (k + 1))
        lo_ii = None if k == 0 else root("l", 2 * k - 1)
        up_ii = root("u", 2 * k + 1)
        cells["i_%d" % k] = (lo_i, up_i)
        cells["ii_%d" % k] = (lo_ii, up_ii)
    names = {("i_0", 1): "N10.1c", ("i_1", 0): "N10.1d", ("i_1", 1): "N10.1e", ("ii_0", 1): "N10.1f", ("ii_1", 0): "N10.1g", ("ii_1", 1): "N10.1h"}
    out["cells"] = {}
    for key, (a, b) in cells.items():
        ent = {}
        for j, e in ((0, a), (1, b)):
            if e is None:
                ent[["lower", "upper"][j]] = None
                continue
            mid = (e[0] + e[1]) / 2
            ent[["lower", "upper"][j]] = [e[0], e[1]]
            lab = names.get((key, j), "")
            log("%-7s cell %-5s %s endpoint h in [%s, %s]  midpoint %s" % (lab, key, ["lower", "upper"][j], e[0].str(25, radius=False), e[1].str(25, radius=False), mid.str(20, radius=False)))
        out["cells"][key] = ent
    # sigma_1
    s1 = bisect(lambda s: Mfun(s) - pi / 4, arb(12) / 10, arb(3), width)
    out["sigma_1"] = [s1[0], s1[1]]
    log("N10.1i  sigma_1 with M(sigma_1) = pi/4: in [%s, %s]" % (s1[0].str(25, radius=False), s1[1].str(25, radius=False)))
    # prime sums (consistency)
    N = pmax
    sieve = bytearray([1]) * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    primes = [i for i in range(N + 1) if sieve[i]]
    for sig in (2, 3):
        S = arb(0)
        for p in primes:
            lp = arb(p).log()
            q = p
            while q <= N:
                S += lp * arb(q) ** (-sig)
                q *= p
        tail = arb(N) ** (1 - sig) * (arb(N).log() / (sig - 1) + arb(1) / (sig - 1) ** 2)
        exact = Mfun(arb(sig))
        ok = (exact >= S) and (exact <= S + tail)
        out["M%d_sum" % sig] = S
        log("check   M(%d): sum over n <= %d = %s, tail bound %s, ball value %s, inside [sum, sum + tail]: %s" % (sig, N, S.str(15, radius=False), tail.str(3, radius=False), exact.str(15, radius=False), ok))
        assert ok
    Jp = arb(0)
    for p in primes:
        Jp += (arb(p) ** (-2)).atan() - (arb(p) ** (arb(-39) / 10)).atan()
    tail = arb(1) / (N - 1)
    log("check   J(2,3.9): sum over p < %d = %s, tail bound %s, so J in [%s, %s]; I(2,3.9) = %s; J <= I: %s; J > pi/40: %s" % (N, Jp.str(15, radius=False), tail.str(3, radius=False), Jp.str(10, radius=False), (Jp + tail).str(10, radius=False), Ifun(2, arb(39) / 10).str(15, radius=False), (Jp + tail) <= Ifun(2, arb(39) / 10) or "not decided", Jp > pi / 40))
    out["J_2_3.9_partial"] = Jp
    return out


def ser(x):
    if isinstance(x, arb):
        return x.str(40, radius=True)
    if isinstance(x, (list, tuple)):
        return [ser(y) for y in x]
    if isinstance(x, dict):
        return {k: ser(v) for k, v in x.items()}
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or ("/tmp/pg10_window.json" if a.quick else "data/pg10_window.json")
    logf = open(out.rsplit(".", 1)[0] + ".log", "w")
    def log(s):
        print(s)
        logf.write(s + "\n")
        logf.flush()
    t0 = time.time()
    bits1, bits2, width, pmax = (100, 130, arb(10) ** -15, 20000) if a.quick else (200, 260, arb(10) ** -35, 200000)
    log("pg10_window: precision %d bits (first run), %d bits (second run); bisection width %s; primes to %d" % (bits1, bits2, "1e-15" if a.quick else "1e-35", pmax))
    r1 = run(bits1, width, pmax, log)
    log("--- second run at %d bits ---" % bits2)
    ctx.prec = bits2
    r2 = run(bits2, width, pmax, lambda s: None)
    # compare
    ctx.prec = bits2
    worst = arb(0)
    for key, ent in r1["cells"].items():
        for kk in ("lower", "upper"):
            e1, e2 = ent[kk], r2["cells"][key][kk]
            if e1 is None:
                continue
            ok = e1[0] <= e2[1] and e2[0] <= e1[1]
            d = abs((e1[0] + e1[1]) / 2 - (e2[0] + e2[1]) / 2)
            worst = d if d > worst else worst
            assert ok, (key, kk)
    s1, s2 = r1["sigma_1"], r2["sigma_1"]
    assert s1[0] <= s2[1] and s2[0] <= s1[1]
    d = abs((s1[0] + s1[1]) / 2 - (s2[0] + s2[1]) / 2)
    worst = d if d > worst else worst
    log("second run: all enclosures overlap; largest difference of midpoints (cells and sigma_1) %s" % worst.str(3, radius=False))
    r1["largest_midpoint_difference_between_precisions"] = worst
    r1["quick"] = a.quick
    json.dump(ser(r1), open(out, "w"), indent=1)
    log("wrote %s (%.1f s)" % (out, time.time() - t0))


if __name__ == "__main__":
    main()
