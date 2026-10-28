# pg8_main.py: the zeros of xi^(inf) + Xi^E, and of xi^(inf) alone, in (50, 500] and in the box [50, 500] x [-3, 3] (Chapter 22, Section
# mi:sec:mixture, the third computation of the last paragraph; placeholders [N8.5a] to [N8.5c] of reconcile-X8.md). Added for the book in
# October 2026, written from the mathematics (the collaborator's scripts were not read, run or copied).
# Object. By Proposition mi:prop:polar (psi_E) and the definition of xi^(inf) (the atoms of mi:eq:psiE with every n), Xi_main = xi^(inf) + Xi^E is the
# transform of 8 pi^2 cosh(9u/2) [ sum_{n>=1} n^4 e^{-2 pi n^2 cosh 2u} + r_0 e^{-2 pi (cosh 2u - c_0)}/(cosh 2u - c_0) ], that is, with the relative
# normalisation of mi:eq:psiE, the sum of 8 pi^2 n^4 (1/2){K_{9/4+iz/2}(2 pi n^2) + K_{9/4-iz/2}(2 pi n^2)} (n >= 1) and 8 pi^2 r_0 calE_{2 pi/9}(z); r_0 is
# recomputed here (certified, as in pg8_constants.py). Truncation: the terms n <= 8 are used; the second quadrature uses n <= 10 as well as another
# step and cut, so the difference logged includes the tail n >= 9 (the n-th atom at height t is below e^{-2 pi n^2} times a factor of size at most e^{t pi/4},
# i.e. 1e-220 relative to the data for n = 9 at t = 500, so the tail is far below the working precision).
# Method. As in pg8_meixner.py (contour Im u = 0.77, scaling e^{0.77 z}, H = 1000, v <= 4.2, 160 bits; sign changes, ball refinement and certification
# of real zeros; argument principle with the Taylor/Lipschitz test on [50, 500] x [-3, 3]; nonreal zeros from grid minima plus Newton; the first
# nonreal zero enclosed by interval Newton with an allowance for the quadrature error). Both functions are analysed: Xi_main, and xi^(inf) alone
# (the r_0 term omitted), the latter for the editor's information.
# Status: computed (ball arithmetic for the rounding error, quadrature error estimated by the second run); the winding number and the enclosure are
# exact for the discretised function.
# Output: data/pg8_main.json and data/pg8_main.log. Usage: python3 pg8_main.py [--quick] [--out FILE] [--cpus N] [--x1 X]; --quick (H = 250, 64 bits).
# Hopper: bash misc/tools/hopper_run.sh -c 8 -t 240 -n pg8_main -g pg8_main.json -g pg8_main.log code/pg8_main.py code/pg8_constants.py -- --out pg8_main.json
import argparse
import json
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
    ap.add_argument("--tmax", type=float, default=500.0)
    ap.add_argument("--x1", type=float, default=None)
    ap.add_argument("--only", default=None, help="xi_inf+Xi^E or xi_inf")
    ap.add_argument("--noap", action="store_true", help="skip the argument principle (quick tests only)")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    out = a.out or (os.path.join(here, "..", "data", "pg8_main.json") if not a.quick else "/tmp/pg8_main.json")
    log = C.Log(os.path.splitext(out)[0] + ".log")
    if a.quick:
        H, V, H2, V2, prec, x1 = 250, 3.6, 300, 3.9, 64, (a.x1 or 120.0)
    else:
        H, V, prec = (a.H or 1000), (a.V or 4.2), (a.prec or 160)
        H2, V2 = int(H * 1.25), V + 0.4
        x1 = (a.x1 or a.tmax)
    t0 = time.time()
    ctx.prec = 260
    ph, im, tail = C.phi_imag_axis(arb.pi() / 9, 14, 260)
    r0 = (2 * arb.pi() / 9).sin() * ph / (18 * arb.pi() ** 2)
    log("r_0 = %s (certified ball radius %s)" % (C.S(r0, 40), r0.rad().str(3)))
    ctx.prec = prec
    results = {"parameters": {"H": H, "V": V, "H2": H2, "V2": V2, "prec": prec, "tmax": a.tmax, "tau": 0.77, "quick": a.quick, "atoms": 8, "atoms_second": 10},
               "r0": C.S(r0, 40)}
    for name, withE in [("xi_inf+Xi^E", True), ("xi_inf", False)]:
        if a.only and name != a.only:
            continue
        qa = C.Quad(C.kern_main(r0, 8, withE), H, V)
        qb = None if a.quick else C.Quad(C.kern_main(r0, 10, withE), H2, V2)
        log("== %s: nodes built (%.1f s)" % (name, time.time() - t0))
        r = C.analyse(name, qa, qb, a.tmax, (50.0, x1, 3.0), a.cpus, log, quick=a.quick, ap=not a.noap)
        rz = [float(x) for x in r["roots"]]
        r["real_zeros_50_tmax"] = sum(1 for x in rz if 50 < x <= a.tmax)
        r["real_zeros_0_tmax"] = sum(1 for x in rz if x <= a.tmax)
        log("[%s] real zeros in (50,%g]: %d (in (0,%g]: %d)" % (name, a.tmax, r["real_zeros_50_tmax"], a.tmax, r["real_zeros_0_tmax"]))
        if r.get("nonreal_zeros_upper") and not a.quick:
            zs = r["nonreal_zeros_upper"]
            z1 = C.newton_c(qa, acb(arb(zs[0][0]), arb(zs[0][1])), tol=1e-34)
            qerr = float(abs(qa.eval(z1) - qb.eval(z1)).mid()) * 10 + 1e-40
            ok, rad, Nx = C.enclose(qa, z1, qerr)
            r["first_nonreal_zero"] = {"re": C.S(z1.real, 30), "im": C.S(z1.imag, 30), "enclosure_radius": "%g" % rad, "krawczyk_contains": ok, "quadrature_allowance": "%.3g" % qerr}
            log("[%s] first nonreal zero z = %s + %s i; interval Newton in the box +-%g (allowance %.3g for the quadrature error): %s"
                % (name, C.S(z1.real, 30), C.S(z1.imag, 30), rad, qerr, "encloses exactly one simple zero" if ok else "FAILED"))
        results[name] = r
        with open(out, "w") as fh:
            json.dump(results, fh, indent=1)
    log("total time %.1f s; wrote %s" % (time.time() - t0, out))


if __name__ == "__main__":
    main()
