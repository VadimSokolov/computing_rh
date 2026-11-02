#!/usr/bin/env python3
"""Track C, script 3b. How close to the floor can a sum of uniforms get? Starts from the coefficient
pattern pi/gamma with a double zero at 2 gamma_1, then free search. The final law is an explicit
sum of uniforms: the numbers gamma only serve as a starting point of the search, and no claim about
zeros rests on them. Run: timeout 500 python3 -I C_uniform_fit2.py (computed, double precision)
"""
import time
import numpy as np
from scipy import optimize
import importlib.util
import os
import sys
here = os.path.dirname(os.path.abspath(sys.argv[0]))
spec = importlib.util.spec_from_file_location("C_uniform_fit", os.path.join(here, "C_uniform_fit.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

G = np.array([14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073, 48.005151, 49.773832])
EG = m.EG


def fit(start, J, label, tmax):
    def obj(v):
        cs = np.abs(v[:J])
        b = v[J]
        if b <= -J + 0.05:
            return 1.0
        return m.tv_of_phi(m.sinc_prod(cs) * m.tail_phi(0.5, b, J + 1))
    t0 = time.time()
    v = np.array(start, dtype=float)
    best = obj(v)
    print(f"{label}: start TV {best:.4e}", flush=True)
    while time.time() - t0 < tmax:
        res = optimize.minimize(obj, v, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 1200, "adaptive": True})
        if res.fun < best * (1 - 1e-4):
            best, v = res.fun, res.x
        else:
            res = optimize.minimize(obj, v, method="Powell", options={"xtol": 1e-7, "ftol": 1e-10, "maxiter": 30})
            if res.fun < best * (1 - 1e-4):
                best, v = res.fun, res.x
            else:
                break
    cs = np.abs(v[:J])
    k2, k4, k6 = m.cum_of(cs, 0.5, v[J], J + 1)
    print(f"{label}: TV={best:.4e} b={v[J]:.5f} Var rel {k2/m.VARX-1:+.1e} k4 rel {k4/m.KAPPA4-1:+.1e}")
    print("   pi/c_k =", np.array2string(np.pi / cs, precision=4, max_line_width=180), flush=True)
    return best, v


def main():
    from scipy import special
    # restart from the first run (C_uniform_fit2_first_run.log, TV 5.1347e-07), whose start was the pattern
    # pi/gamma_1, pi/gamma_2, pi/gamma_3, pi/(2 gamma_1), pi/gamma_4, pi/gamma_5, pi/gamma_6 with harmonic tail
    th = np.array([14.1348, 36.249, 38.4304, 21.0234, 25.3371, 43.2102, 59.9942])
    cs = np.pi / th
    best, v = fit(list(cs) + [1.02134], 7, "restart (J=7)", 90)
    for Jn in (9, 12):
        Jold = len(v) - 1
        b = v[-1]
        extra = [0.5 / (k + b) for k in range(Jold + 1, Jn + 1)]
        best, v = fit(list(np.abs(v[:-1])) + extra + [b], Jn, f"extended (J={Jn})", 170)
        print("coefficients c_k:", np.array2string(np.abs(v[:-1]), precision=10, max_line_width=180), "b =", v[-1], flush=True)


if __name__ == "__main__":
    main()
