"""Scan of F(x)=2*int_0^inf e^{-u} Phi(u) cos(xu) du on [0,80] for ch/nc_riemanncont.tex.
Floating point (mpmath), not certified. Writes data/nc_riemanncont_scan.json (paths relative to this script).
|F''(x)| <= D2 = 2 int_0^inf u^2 e^{-u} Phi(u) du, so linear interpolation at step h errs by at most h^2*D2/8."""
import json, os
import mpmath as mp
mp.mp.dps = 25
pi = mp.pi
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "nc_riemanncont_scan.json")

def Phi(u):
    u = abs(u)
    return 2*mp.fsum((2*pi**2*n**4*mp.exp(4.5*u)-3*pi*n*n*mp.exp(2.5*u))*mp.exp(-pi*n*n*mp.exp(2*u)) for n in range(1, 12))

U = mp.mpf(3)          # Phi(3) is below 1e-100000, so [0,3] suffices
panels = 120
nodes, weights = [], []
# Gauss-Legendre of degree 16 on [-1,1]
gl = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(4, mp.mp.prec)  # 3*2^3 = 24 nodes per panel
h = U/panels
for p in range(panels):
    a = p*h
    for (t, w) in gl:
        nodes.append(a + h*(t+1)/2); weights.append(w*h/2)
vals = [mp.exp(-u)*Phi(u) for u in nodes]
def F(x):
    return 2*mp.fsum(w*v*mp.cos(x*u) for u, w, v in zip(nodes, weights, vals))
D2 = 2*mp.fsum(w*v*u*u for u, w, v in zip(nodes, weights, vals))
Phi0 = Phi(0)
# validation of the fixed quadrature against adaptive quad
chk = {}
for x in [mp.mpf(1), mp.mpf(37), mp.mpf(80)]:
    q = 2*mp.quad(lambda u: mp.exp(-u)*Phi(u)*mp.cos(x*u), mp.linspace(0, 3, 61))
    chk[str(x)] = mp.nstr(abs(q-F(x)), 5)
step = mp.mpf(1)/4
xs = [k*step for k in range(0, 321)]
Fv = [F(x) for x in xs]
imin = min(range(len(Fv)), key=lambda i: Fv[i])
signs = [mp.nstr(xs[i], 6) for i in range(len(xs)-1) if Fv[i]*Fv[i+1] < 0]
bound = step**2*D2/8
concl = (not signs) and bound < Fv[imin] and Fv[imin] > 0
res = {"D2": mp.nstr(D2, 12), "Phi0": mp.nstr(Phi0, 18), "step": "1/4", "interval": [0, 80],
       "min_F_grid": mp.nstr(Fv[imin], 12), "argmin_x": mp.nstr(xs[imin], 8),
       "interp_error_bound": mp.nstr(bound, 12), "sign_changes_on_grid": signs,
       "quad_vs_adaptive_abs_diff": chk,
       "F_grid": [mp.nstr(v, 10) for v in Fv],
       "conclusion": "no real zero of F on [0,80] (floating point, not certified)" if concl else "conclusion does not hold",
       "conclusion_holds": bool(concl)}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(res, f, indent=1)
for k in ["D2", "min_F_grid", "argmin_x", "interp_error_bound", "sign_changes_on_grid", "quad_vs_adaptive_abs_diff", "conclusion"]:
    print(k, res[k])
