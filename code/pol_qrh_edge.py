"""pol_qrh_edge.pdf: the Thorin density g(alpha, tau) = Re xi'/xi(alpha + i tau) of xi(alpha)/xi(alpha + sqrt u)
at alpha = 7/8 (the edge of the half plane of Q_{7/8}), 1 and 3/2, for 0 <= tau <= 60, computed from xi with no zeros.
Black and white (line styles only). Writes pol_qrh_edge.pdf and pol_qrh_edge.json (minimum, its location, peak values).
Based on qrh_fig.py of N. Polson; rewritten for the print edition.  Usage: python3 pol_qrh_edge.py [OUTDIR]"""
import sys, os, json
import mpmath as mp, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out = sys.argv[1] if len(sys.argv) > 1 else "."
mp.mp.dps = 20
def g(a, t):
    s = mp.mpc(a, t)
    return float(mp.re(1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.psi(0, s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)))
ts = np.arange(0.0, 60.0001, 0.05); ts[0] = 1e-9
styles = [(0.875, "-", r"$\alpha=7/8$"), (1.0, "--", r"$\alpha=1$"), (1.5, ":", r"$\alpha=3/2$")]
res = {}
plt.rcParams.update({"font.size": 8})
fig, ax = plt.subplots(figsize=(5.3, 2.2))
for a, ls, lab in styles:
    v = np.array([g(a, t) for t in ts])
    pk = [float(v[i]) for i in range(1, len(v)-1) if v[i] > v[i-1] and v[i] > v[i+1]]
    res[str(a)] = {"min": float(v.min()), "argmin": float(ts[v.argmin()]), "min_peak": min(pk), "max_peak": max(pk)}
    ax.plot(ts, v, ls, color="k", lw=0.8, label=lab)
ax.set_xlim(0, 60); ax.set_ylim(0, 3.2)
ax.set_xlabel(r"$\tau$"); ax.set_ylabel(r"$\mathrm{Re}\,\xi'/\xi(\alpha+i\tau)$")
ax.legend(frameon=False, ncol=3, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(out, "pol_qrh_edge.pdf"))
json.dump(res, open(os.path.join(out, "pol_qrh_edge.json"), "w"), indent=1)
print(res)
