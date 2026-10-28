# pg8_E.py: the zeros of Xi^E = 8 pi^2 r_0 calE_{2 pi/9} (a_0 = 2 pi) on (0, 750] and in the box [0.25, 750] x [-3, 3] (Chapter 22, Section
# mi:sec:mixture, the first of the three computations of the last paragraph; placeholders [N8.4a] to [N8.4d] of reconcile-X8.md). Added for the
# book in October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Object. calE_theta(z) = int cosh(9u/2) e^{-a_0 (cosh 2u - cos theta)} / (cosh 2u - cos theta) e^{izu} du with theta = 2 pi/9 (Proposition
# mi:prop:rigid), so cos theta = c_0 and the apparent pole of the kernel at u = i pi/9 is cancelled by the zero of cosh(9u/2) (it is removable:
# by mi:eq:nine the kernel is U(u) e^{-2 pi (cosh 2u - c_0)}/(cosh u + c_1) with U entire and the first pole of 1/(cosh u + c_1) at 8 i pi/9). The
# normalising constant 8 pi^2 r_0 does not move the zeros and is not applied; r_0 itself is [N8.1] (pg8_constants.py).
# Method. As in pg8_meixner.py (contour Im u = 0.77, e^{0.77 z} scaling, H = 1000 nodes per unit, v <= 4.2, 160 bits; real zeros by sign changes
# on a grid of step 0.01 refined and certified as sign changes of the discretised function; the argument principle on [0.25, 750] x [-3, 3] with
# the Taylor/Lipschitz test of track(); second quadrature H = 1250, v <= 4.6). Least gap: the smallest difference of consecutive real zeros
# and the lower zero of that pair. Main term: M(T) = (T/2 pi) log(T/(2 pi e)), reported as M(750) - M(50), against the number of zeros in (50, 750].
# Status: computed (ball arithmetic for the rounding error, quadrature error estimated by the second run); the winding number is exact for the
# discretised function.
# Output: data/pg8_E.json and data/pg8_E.log. Usage: python3 pg8_E.py [--quick] [--out FILE] [--cpus N] [--x1 X]; --quick (H = 250, 64 bits)
# writes to /tmp unless --out is given. Hopper: bash misc/tools/hopper_run.sh -c 8 -t 240 -n pg8_E -g pg8_E.json -g pg8_E.log code/pg8_E.py code/pg8_constants.py -- --out pg8_E.json
import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
import pg8_constants as C
from flint import acb, arb, ctx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--H", type=int, default=None, help="nodes per unit (default 1000; for tests only)")
    ap.add_argument("--V", type=float, default=None, help="cut in v (default 4.2)")
    ap.add_argument("--prec", type=int, default=None, help="bits (default 160)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("--tmax", type=float, default=750.0)
    ap.add_argument("--x1", type=float, default=None)
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or (os.path.join(here, "..", "data", "pg8_E.json") if not a.quick else "/tmp/pg8_E.json")
    log = C.Log(os.path.splitext(out)[0] + ".log")
    if a.quick:
        H, V, H2, V2, prec, x1 = 250, 3.6, 300, 3.9, 64, (a.x1 or 60.0)
    else:
        H, V, prec = (a.H or 1000), (a.V or 4.2), (a.prec or 160)
        H2, V2 = int(H * 1.25), V + 0.4
        x1 = (a.x1 or a.tmax)
    ctx.prec = prec
    t0 = time.time()
    qa = C.Quad(C.kern_E(), H, V)
    qb = None if a.quick else C.Quad(C.kern_E(), H2, V2)
    log("nodes built (%.1f s); H = %d, v <= %g, %d bits" % (time.time() - t0, H, V, prec))
    r = C.analyse("Xi^E", qa, qb, a.tmax, (0.25, x1, 3.0), a.cpus, log, quick=a.quick)
    rz = [float(x) for x in r["roots"]]
    r["real_zeros_0_tmax"] = sum(1 for x in rz if x <= a.tmax)
    n50 = sum(1 for x in rz if 50 < x <= a.tmax)
    M = lambda T: (T / (2 * math.pi)) * math.log(T / (2 * math.pi * math.e))
    r["real_zeros_50_tmax"] = n50
    r["main_term_50_tmax"] = "%.4f" % (M(a.tmax) - M(50.0))
    r["parameters"] = {"H": H, "V": V, "H2": H2, "V2": V2, "prec": prec, "tmax": a.tmax, "tau": 0.77, "quick": a.quick}
    log("real zeros in (0,%g]: %d; in (50,%g]: %d; main term M(%g) - M(50) = %.4f" % (a.tmax, r["real_zeros_0_tmax"], a.tmax, n50, a.tmax, M(a.tmax) - M(50.0)))
    if "least_gap" in r:
        log("least gap %s at height %s (lower zero), %s (upper zero)" % (r["least_gap"]["gap"], r["least_gap"]["lower_zero"][:12], r["least_gap"]["upper_zero"][:12]))
    with open(out, "w") as fh:
        json.dump(r, fh, indent=1)
    log("total time %.1f s; wrote %s" % (time.time() - t0, out))


if __name__ == "__main__":
    main()
