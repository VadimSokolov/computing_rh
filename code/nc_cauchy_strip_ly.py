# Strip Lee-Yang test. Single spin: s in {-1,0,1}, P(s=+-1)=q, P(0)=1-2q, q<1/4: zeros of L(h) at Re h = +-a,
# cosh a = (1-2q)/(2q) (Cauchy defect a). n spins, ferromagnetic couplings, common field h.
# Curie Weiss: Z(h) = sum w(s) exp((J/n) sum_{i<j} s_i s_j + h sum s_i). Chain: transfer matrix, free ends.
# delta(J) = max |Re h| over zeros (Re h = log|z|, z = e^h). Floating point (numpy roots, extended precision via mpmath polyroots).
import numpy as np, mpmath as mp, json
from math import comb, exp, log, acosh
mp.mp.dps = 50
def cw_poly(n, q, J):
    c = {}
    for npl in range(n+1):
        for nm in range(n+1-npl):
            n0 = n-npl-nm; M = npl-nm
            wgt = mp.mpf(comb(n, npl)*comb(n-npl, nm)) * mp.mpf(q)**(npl+nm) * mp.mpf(1-2*q)**n0 * mp.e**(mp.mpf(J)/n*(M*M-(npl+nm))/2)
            c[M] = c.get(M, 0) + wgt
    return [c.get(M, 0) for M in range(n, -n-1, -1)]   # coefficients of z^{M+n}, highest first
def chain_poly(n, q, J):
    # Z as Laurent polynomial in z via dynamic programming over magnetisation
    w = {-1: mp.mpf(q), 0: mp.mpf(1-2*q), 1: mp.mpf(q)}
    st = {(s, s): w[s] for s in (-1, 0, 1)}          # (last spin, M) -> weight
    for _ in range(n-1):
        nx = {}
        for (s, M), v in st.items():
            for t in (-1, 0, 1):
                key = (t, M+t); nx[key] = nx.get(key, 0) + v*w[t]*mp.e**(mp.mpf(J)*s*t)
        st = nx
    c = {}
    for (s, M), v in st.items(): c[M] = c.get(M, 0) + v
    return [c.get(M, 0) for M in range(n, -n-1, -1)]
def defect(coeffs):
    # strip off leading/trailing zeros of the polynomial
    while coeffs and coeffs[-1] == 0: coeffs = coeffs[:-1]
    r = mp.polyroots(coeffs, maxsteps=400, extraprec=400)
    return max(abs(float(mp.log(abs(z)))) for z in r)
out = {}
for q in [0.2, 0.1]:
    a = acosh((1-2*q)/(2*q))
    for kind, fn in [("CW", cw_poly), ("chain", chain_poly)]:
        for n in [2, 4, 8]:
            row = []
            for J in [0, 0.1, 0.3, 0.6, 1.0, 2.0, 4.0]:
                row.append((J, defect(fn(n, q, J))))
            out[f"{kind} q={q} n={n}"] = dict(a=a, row=row)
            print(f"{kind:5s} q={q} a={a:.4f} n={n}: " + "  ".join(f"J={J}:{d:.4f}" for J, d in row), flush=True)
json.dump(out, open('strip_ly.json', 'w'), indent=0)
