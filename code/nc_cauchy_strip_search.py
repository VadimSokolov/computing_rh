# Random search for a counterexample to "ferromagnetic coupling does not widen the zero strip".
# Single spin on {-K..K}, symmetric random weights. n spins, Curie Weiss or chain or two spins, couplings J > 0.
import numpy as np, mpmath as mp, json, itertools, random
mp.mp.dps = 40
random.seed(1)
def single_defect(w):  # w dict value->weight, symmetric
    K = max(w); co = [mp.mpf(w.get(v, 0)) for v in range(K, -K-1, -1)]
    r = mp.polyroots(co, maxsteps=300, extraprec=300)
    return max(abs(float(mp.log(abs(z)))) for z in r)
def system_defect(w, n, J, kind):
    vals = sorted(w); c = {}
    if kind == "CW":
        for conf in itertools.combinations_with_replacement(vals, n):
            # multinomial count
            from collections import Counter
            cnt = Counter(conf); mult = mp.factorial(n)
            for v in cnt.values(): mult /= mp.factorial(v)
            wt = mult
            for v in conf: wt *= w[v]
            M = sum(conf); S2 = sum(x*x for x in conf)
            wt *= mp.e**(mp.mpf(J)/n*(M*M - S2)/2)
            c[M] = c.get(M, 0) + wt
    else:  # chain
        st = {(s, s): mp.mpf(w[s]) for s in vals}
        for _ in range(n-1):
            nx = {}
            for (s, M), v in st.items():
                for t in vals:
                    k = (t, M+t); nx[k] = nx.get(k, 0) + v*w[t]*mp.e**(mp.mpf(J)*s*t)
            st = nx
        for (s, M), v in st.items(): c[M] = c.get(M, 0) + v
    Mmax = max(abs(M) for M in c)
    co = [c.get(M, 0) for M in range(Mmax, -Mmax-1, -1)]
    r = mp.polyroots(co, maxsteps=400, extraprec=400)
    return max(abs(float(mp.log(abs(z)))) for z in r)
worst = []
for trial in range(120):
    K = random.choice([2, 3])
    w = {0: random.random()}
    for v in range(1, K+1):
        w[v] = w[-v] = random.random()*random.choice([1, 0.1, 0.01])
    d0 = single_defect(w)
    if d0 < 1e-6: continue   # already ferromagnetic
    for kind, n in [("CW", 2), ("CW", 3), ("chain", 3), ("chain", 4)]:
        for J in [0.02, 0.1, 0.3]:
            d = system_defect(w, n, J, kind)
            worst.append((d - d0, trial, kind, n, J, d0, d))
worst.sort(reverse=True)
print("largest widening found (d(J) - d0):")
for r in worst[:8]: print(r)
print("cases:", len(worst), " widened (>1e-9):", sum(1 for r in worst if r[0] > 1e-9))
json.dump(worst[:50], open('strip_search.json', 'w'))
