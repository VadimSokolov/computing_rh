# Figure cm_pairs: (left) two planted pairs of equal width, gap of each against the separation d, measured and from the
# fixed point of the gap proposition; (centre) the same gap divided by eta as the pairs coalesce, against d/eta, with
# the double pair limit 1 - 1/sqrt 2; (right) slice threshold of a pair of multiplicity n over eta, with cos(pi/(2n)).
import json, math, os, sys
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "code")
far = json.load(open(f"{d}/two_pairs_k8.json")); close = json.load(open(f"{d}/two_pairs_close_k8.json"))
mult = json.load(open(f"{d}/multiplicity.json"))
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for eta, col in [(0.1, "C0"), (0.2, "C1")]:
    for j, mk in [(0, "o"), (1, "s")]:
        rows = sorted([r for r in far if r["eta"] == eta and r["other"] == eta and r["j"] == j], key=lambda r: r["d"])
        ax[0].semilogy([r["d"] for r in rows], [r["gap"] for r in rows], mk, color=col, label=f"$\\eta={eta}$, pair {j+1}")
        ax[0].semilogy([r["d"] for r in rows], [r["fixed"] for r in rows], "-", color=col, lw=0.8)
ax[0].set_xlabel("separation $d$"); ax[0].set_ylabel(r"gap $\eta-w^*$"); ax[0].legend(fontsize=7)
ax[0].set_title("two pairs: measured (marks), fixed point (lines)", fontsize=9)
eta = 0.2
rows = [r for r in far + close if r["eta"] == eta and r["other"] == eta]
for j, mk in [(0, "o"), (1, "s")]:
    rr = sorted([r for r in rows if r["j"] == j], key=lambda r: r["d"])
    ax[1].semilogx([r["d"]/eta for r in rr], [r["gap"]/eta for r in rr], mk + "-", ms=4, label=f"pair {j+1}")
ax[1].axhline(1 - 0.14931/0.2, color="C1", ls="--", lw=0.8, label="double pair, all four real")
ax[1].axhline(1 - 0.13449/0.2, color="C0", ls="--", lw=0.8, label="double pair, first two real")
ax[1].axhline(1 - 1/math.sqrt(2), color="k", ls=":", lw=0.8, label=r"$1-1/\sqrt{2}$")
ax[1].set_xlabel(r"$d/\eta$"); ax[1].set_ylabel(r"$(\eta-w^*)/\eta$"); ax[1].legend(fontsize=7)
ax[1].set_title(r"coalescing pairs, $\eta=0.2$", fontsize=9)
for eta, col in [(0.1, "C0"), (0.2, "C1")]:
    rr = sorted([r for r in mult if r["eta"] == eta and r["all_real"] is not None], key=lambda r: r["n"])
    ax[2].plot([r["n"] for r in rr], [r["all_real"]/eta for r in rr], "o", color=col, label=f"$\\eta={eta}$")
ns = [1, 2, 3, 4]
ax[2].plot([2, 4], [math.cos(math.pi/4), math.cos(math.pi/8)], "k_", ms=18, label=r"$\cos\frac{\pi}{2n}$, $n$ even")
ax[2].set_xticks(ns); ax[2].set_xlabel("multiplicity $n$"); ax[2].set_ylabel(r"$w^*/\eta$"); ax[2].legend(fontsize=7)
ax[2].set_title("a pair of multiplicity $n$", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(root, "fig", "cm_pairs.pdf"))
