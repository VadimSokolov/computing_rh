# wp2_check.py: independent check of code/wp2_u_certify.py by mpmath (floating point, 30 digits, not certified): u_alpha(v) = pi^-1 Re(xi'/xi)(alpha + i v) from
# xi'/xi = 1/s + 1/(s-1) - log(pi)/2 + psi(s/2)/2 + zeta'(s)/zeta(s), at v = 0 and at the argmin points of the aggregate file (wp2_agg.json), compared with the
# stored midpoint value of the ball run (agreement to the logged float digits expected). Usage: python3 wp2_check.py wp2_agg.json OUT.json
import json, sys
from mpmath import mp, mpf, mpc, zeta, psi, log, pi
mp.dps = 30
def u(alpha, v):
    s = mpc(alpha, v)
    L = 1 / s + 1 / (s - 1) - log(pi) / 2 + psi(0, s / 2) / 2 + zeta(s, derivative=1) / zeta(s)
    return L.real / pi
agg = json.load(open(sys.argv[1]))
res = []
for k, d in agg.items():
    al = d["alpha"]
    a = mpf(eval(al + ".0")) if "/" in al else mpf(al)
    v = mpf(d["least_mid_at"])
    t = u(a, v)
    # the stored point is a cell midpoint, so the two floats should agree to the cell-midpoint value (the ball run stores u at that midpoint)
    r = {"key": k, "alpha": al, "v": float(v), "u_mpmath": float(t), "u_ball_least_mid": d["least_mid"], "rel_diff": float(abs(t - d["least_mid"]) / abs(t)),
         "u0_mpmath": float(u(a, 0)), "u0_ball": d["u0"]}
    res.append(r)
    print(r, flush=True)
json.dump(res, open(sys.argv[2], "w"), indent=1)
