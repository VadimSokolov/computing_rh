# Multivariable strip test (star graph): X a spin one law with strip a, coupled to m Ising spins with couplings w_j > 0,
# fields h0 on X and h_j on the Ising spins. Z = sum_x p(x) e^{h0 x} prod_j 2 cosh(h_j + w_j x).
# Statement tested: Z != 0 whenever Re h0 > a and Re h_j > 0. For random h_j (Re > 0) we compute the zeros of Z in
# z = e^{h0} (a Laurent polynomial of degree 1 in z) ... X takes values -1,0,1 so Z = c_{-1} z^{-1} + c_0 + c_1 z.
import numpy as np, json, random, math
random.seed(3)
worst = []
for trial in range(20000):
    q = random.uniform(0.01, 0.24)
    a = math.acosh((1 - 2 * q) / (2 * q))
    m = random.choice([1, 2, 3, 5])
    w = [random.uniform(0, 3) for _ in range(m)]
    h = [complex(random.uniform(0, 2) ** 2, random.uniform(-6, 6)) for _ in range(m)]
    c = {}
    for x, px in [(-1, q), (0, 1 - 2 * q), (1, q)]:
        val = px
        for wj, hj in zip(w, h):
            val *= 2 * np.cosh(hj + wj * x)
        c[x] = val
    roots = np.roots([c[1], c[0], c[-1]])          # c1 z^2 + c0 z + c_{-1} = 0
    rmax = max(np.log(abs(r)) for r in roots)       # Re h0 of the zeros
    worst.append((rmax - a, q, m, w, [str(x) for x in h]))
worst.sort(key=lambda t: -t[0])
print("largest Re h0 of a zero minus the single spin strip a:", worst[0][0])
print("violations (> 1e-9):", sum(1 for t in worst if t[0] > 1e-9), "of", len(worst))
print(worst[0])
json.dump([(t[0], t[1], t[2]) for t in worst[:20]], open("multifield.json", "w"))
