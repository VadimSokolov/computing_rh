# wp2_aggregate.py: combine the data/wp2_u_*.json task files of code/wp2_u_certify.py into the table per level and height decade
# (certified least cell bound, least midpoint value, comparison with alpha - 1/2, log T, b_alpha = pi u_alpha(0) and V_xi = 0.04620998623).
# Usage: python3 -I wp2_aggregate.py OUT.json file1.json file2.json ...   (the table is printed as markdown)
import json, math, sys
from collections import defaultdict
V = 0.04620998623
out = sys.argv[1]
rows = defaultdict(list)
u0 = {}
for fn in sys.argv[2:]:
    d = json.load(open(fn))
    u0[d["alpha"]] = d["u0_float"]
    for b in d["bins"]:
        rows[d["alpha"]].append(b)
decs = [(0, 10), (10, 100), (100, 1000), (1000, 10000), (10000, 100000)]
agg = {}
print("| alpha | height range | cells (maxj, fails) | least certified cell bound of v = pi u | at v | bound/(alpha-1/2) | least midpoint v | at v | pi u(0) = b_alpha | b_alpha/(alpha-1/2) | log T_hi | last cell end with LB <= u(0) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for al in sorted(rows, key=lambda a: float(eval(a + ".0") if "/" in a else a)):
    a = float(eval(al + ".0")) if "/" in al else float(al)
    for lo, hi in decs:
        bs = [b for b in rows[al] if b["lo"] >= lo and b["hi"] <= hi]
        if not bs:
            continue
        cov = sum(b["hi"] - b["lo"] for b in bs)
        L = min(bs, key=lambda b: b["least_lb"])
        M = min(bs, key=lambda b: b["least_mid"])
        ends = [b["last_end_lb_le_u0"] for b in bs if b["last_end_lb_le_u0"] is not None]
        cells = sum(b["cells"] for b in bs)
        fails = sum(b["fails"] for b in bs)
        maxj = max(b["maxj"] for b in bs)
        pu = math.pi
        key = "%s_%d" % (al, hi)
        agg[key] = {"alpha": al, "lo": lo, "hi": hi, "covered": cov, "least_lb": L["least_lb"], "least_lb_at": L["least_lb_at"], "least_mid": M["least_mid"],
                    "least_mid_at": M["least_mid_at"], "u0": u0[al], "fails": fails, "cells": cells, "maxj": maxj, "end": max(ends) if ends else None}
        print("| %s | %d to %d (covered %d) | %d (%d, %d) | %.6g | %.5f | %.5f | %.6g | %.5f | %.6g | %.5f | %.2f | %s |" % (
            al, lo, hi, cov, cells, maxj, fails, pu * L["least_lb"], L["least_lb_at"], pu * L["least_lb"] / (a - .5), pu * M["least_mid"], M["least_mid_at"],
            pu * u0[al], pu * u0[al] / (a - .5), math.log(hi), ("%.4f" % max(ends)) if ends else "none"))
json.dump(agg, open(out, "w"), indent=1)
