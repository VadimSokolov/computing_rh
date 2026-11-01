# Genus 2 benchmark. C: y^2 = x^5 + x + 1 over F_p. Count N_1 over F_p and N_2 over F_{p^2}; the numerator of Z is
# P(T) = 1 + c1 T + c2 T^2 + p c1 T^3 + p^2 T^4, c1 = N1 - p - 1, c2 = (N2 - p^2 - 1 + c1^2)/2.
# Normalised power sums s_n = (p^n + 1 - N_n)/p^{n/2} = sum_j 2 cos(n theta_j) under RH. RH for C holds iff the Toeplitz
# matrices [s_{|k-l|}] (s_0 = 4) are positive semidefinite for every size (Caratheodory and Toeplitz). A fake numerator
# with a pair of roots off the circle |T| = p^{-1/2} gives Toeplitz matrices that fail at some size.
import numpy as np, json


def legendre_table(p):
    t = [0] * p
    for y in range(1, p):
        t[y * y % p] = 1
    return [0] + [1 if t[a] else -1 for a in range(1, p)]


def count_Fp(p):
    chi = legendre_table(p)
    return p + 1 + sum(chi[(x ** 5 + x + 1) % p] for x in range(p))


def count_Fp2(p):
    # F_{p^2} = F_p[i]/(i^2 - n) with n a nonresidue; chi on F_{p^2}: a is a square iff Norm(a) is a square in F_p
    chi = legendre_table(p)
    n = next(a for a in range(2, p) if chi[a] == -1)
    def mul(a, b):
        return ((a[0] * b[0] + n * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
    total = 0
    for x0 in range(p):
        for x1 in range(p):
            x = (x0, x1); x2 = mul(x, x); x4 = mul(x2, x2); x5 = mul(x4, x)
            f = ((x5[0] + x0 + 1) % p, (x5[1] + x1) % p)
            if f == (0, 0):
                continue
            norm = (f[0] * f[0] - n * f[1] * f[1]) % p
            total += chi[norm]
    return p * p + 1 + total


def smooth_mod_p(p, f=(1, 1, 0, 0, 0, 1)):
    # f = coefficients of x^5 + x + 1, lowest degree first. The model y^2 = f(x) is a smooth genus 2 curve over F_p
    # only for odd p with gcd(f, f') = 1 mod p. For p = 2 the model y^2 = f(x) is never smooth, and where gcd(f, f')
    # is nontrivial (here p = 3, 7, 23, the primes dividing the discriminant 3381) the curve has a node, the numerator
    # of Z is not a genus 2 quartic and Proposition route:prop:toeplitz does not apply; such primes are skipped.
    if p == 2:
        return False
    a = [c % p for c in f]
    b = [(i * c) % p for i, c in enumerate(f)][1:]
    def trim(u):
        while u and u[-1] == 0:
            u.pop()
        return u
    a, b = trim(a), trim(b)
    while b:
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b):
            q = a[-1] * inv % p; d = len(a) - len(b)
            for i, c in enumerate(b):
                a[i + d] = (a[i + d] - q * c) % p
            a = trim(a)
            if not a:
                break
        a, b = b, a
    return len(a) == 1


def toeplitz_min(s, size):
    T = np.array([[s[abs(k - l)] for l in range(size)] for k in range(size)], float)
    return np.linalg.eigvalsh(T).min()


out = []
for p in [7, 11, 13, 17, 19, 23, 29, 31]:
    if not smooth_mod_p(p):
        print(f"p={p}: singular model, skipped"); continue
    N1, N2 = count_Fp(p), count_Fp2(p)
    c1 = N1 - p - 1; c2 = (N2 - p * p - 1 + c1 * c1) // 2
    roots = np.roots([p * p, p * c1, c2, c1, 1])          # roots of P(T) as polynomial in T (reversed)
    alphas = 1 / roots                                     # Frobenius eigenvalues
    mods = np.abs(alphas) / np.sqrt(p)
    s = [4.0] + [float(np.real(np.sum(alphas ** k)) / p ** (k / 2)) for k in range(1, 40)]
    mins = [toeplitz_min(s, n) for n in (2, 4, 6, 10, 20, 40)]
    out.append(dict(p=p, N1=N1, N2=N2, c1=c1, c2=c2, mod_range=[float(mods.min()), float(mods.max())], toeplitz_min=mins))
    print(f"p={p}: N1={N1} N2={N2} c1={c1} c2={c2} |alpha|/sqrt p in [{mods.min():.6f},{mods.max():.6f}]  "
          f"min eig sizes 2..40: {[round(m, 6) for m in mins]}", flush=True)
# fake numerator at p = 13: keep c1, push c2 outside the region allowed by RH
p = 13; rec = [r for r in out if r["p"] == p][0]
for c2f in [rec["c2"] + 12, rec["c2"] + 25, rec["c2"] - 30]:
    roots = np.roots([p * p, p * rec["c1"], c2f, rec["c1"], 1]); alphas = 1 / roots
    mods = np.abs(alphas) / np.sqrt(p)
    s = [4.0] + [float(np.real(np.sum(alphas ** k)) / p ** (k / 2)) for k in range(1, 40)]
    mins = [(n, toeplitz_min(s, n)) for n in range(2, 41)]
    first_bad = next((n for n, m in mins if m < -1e-9), None)
    out.append(dict(p=p, fake_c2=c2f, mod_range=[float(mods.min()), float(mods.max())], first_negative_size=first_bad,
                    min_at_40=mins[-1][1]))
    print(f"fake p={p} c2={c2f}: |alpha|/sqrt p in [{mods.min():.4f},{mods.max():.4f}], first Toeplitz size with a "
          f"negative eigenvalue: {first_bad}, min eig at size 40: {mins[-1][1]:.3e}")
json.dump(out, open("genus2.json", "w"), indent=0)
