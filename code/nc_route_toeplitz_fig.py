import numpy as np, matplotlib, json
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"serif","font.size":9})
def s_from(p, c1, c2, nmax=40):
    alphas = 1/np.roots([p*p, p*c1, c2, c1, 1])
    return [4.0] + [float(np.real(np.sum(alphas**k))/p**(k/2)) for k in range(1, nmax)], np.abs(alphas)/np.sqrt(p)
def mins(s, N=25):
    return [np.linalg.eigvalsh(np.array([[s[abs(k-l)] for l in range(n)] for k in range(n)])).min() for n in range(1, N+1)]
fig, ax = plt.subplots(figsize=(6.6, 2.6))
for (lab, p, c1, c2, col, ls) in [("curve, $p=13$", 13, 1, 4, "#1f4e79", "-"), ("curve, $p=29$", 29, 5, 40, "#548235", "-"),
                                   ("fake, $p=13$, $c_2=29$", 13, 1, 29, "#c55a11", "--"), ("fake, $p=13$, $c_2=-26$", 13, 1, -26, "#7f3f98", "--")]:
    s, m = s_from(p, c1, c2); v = mins(s)
    y = [np.sign(x)*np.log10(1+abs(x)) for x in v]
    ax.plot(range(1, len(v)+1), y, ls, marker="o", ms=2.5, lw=0.9, color=col, label=lab)
ax.axhline(0, color="k", lw=0.4)
ax.set_xlabel("size of the Toeplitz section"); ax.set_ylabel("sign$\\cdot\\log_{10}(1+|\\lambda_{\\min}|)$")
ax.set_title("$y^2=x^5+x+1$: Toeplitz positivity of normalised point counts", fontsize=8)
ax.legend(frameon=False, fontsize=7); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("fig_toeplitz.pdf")
