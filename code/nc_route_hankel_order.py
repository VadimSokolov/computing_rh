# How large a finite section does zeta need? Hankel test of the Stieltjes moments M_j = sum_k gamma_k^{-2j} (one atom
# per pair), Proposition ev:prop:hankel: RH iff det[M_{i+j+1}] > 0 and det[M_{i+j+2}] > 0 for every order n.
# Plant a pair c = (g_k+g_{k+1})/2 + i eta in place of g_k, g_{k+1} and find the least order n at which a determinant is
# not positive. Zeros 1..N from mpmath; the tail above g_N as the positive density (1/2pi) log(t/2pi) dt. Floating point.
import mpmath as mp, json, sys
mp.mp.dps = 260
N = 220
G = [mp.im(mp.zetazero(j)) for j in range(1, N+1)]
T = G[-1]
def tail(j):   # integral_T^inf t^{-2j} (1/2pi) log(t/2pi) dt
    return mp.quad(lambda t: t**(-2*j)*mp.log(t/(2*mp.pi))/(2*mp.pi), [T, mp.inf])
JMAX = 90
base = [None] + [mp.fsum(g**(-2*j) for g in G) + tail(j) for j in range(1, JMAX+1)]
def first_failure(M, nmax):
    for n in range(1, nmax+1):
        if 2*n+1 > JMAX: return None
        D1 = mp.det(mp.matrix([[M[i+j+1] for j in range(n)] for i in range(n)]))
        D2 = mp.det(mp.matrix([[M[i+j+2] for j in range(n)] for i in range(n)]))
        if mp.re(D1) <= 0 or mp.re(D2) <= 0:
            return n
    return None
print("genuine zeros: first failure up to order 40:", first_failure(base, 40), flush=True)
out = []
for eta in [mp.mpf('0.3'), mp.mpf('0.05')]:
    for k in [1, 2, 3, 5, 8, 12, 17, 24, 32]:
        gk, gk1 = G[k-1], G[k]; c = (gk+gk1)/2 + 1j*eta
        M = [None] + [base[j] - gk**(-2*j) - gk1**(-2*j) + 2*mp.re(c**(-2*j)) for j in range(1, JMAX+1)]
        n = first_failure(M, 44)
        out.append(dict(k=k, height=float(mp.re(c)), eta=float(eta), first_order=n))
        print(f"eta={float(eta)} k={k:2d} height={float(mp.re(c)):7.2f}: first Hankel order that fails = {n}", flush=True)
json.dump(out, open("hankel_order.json", "w"), indent=0)
