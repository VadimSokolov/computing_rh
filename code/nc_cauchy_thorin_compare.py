# Compare the Thorin density (level y in the variable of the Polya law) of p with that of the
# zero completed ferromagnetic laws f_K of Volume II, Theorem 19.x (keep K zeros, Gaussian drift).
# u_p(y,v) = (1/pi) Re xi'/xi(1/2+y+iv);  u_K(y,v) = (2 beta_K/pi) y + sum_{k<=K} [P_y(v-g_k)+P_y(v+g_k)]
import mpmath as mp, json, sys
mp.mp.dps = 30
VarX = mp.mpf('0.0462099862')
K_list = [16, 150]
Kmax = max(K_list)
g = [mp.im(mp.zetazero(k)) for k in range(1, Kmax+1)]
def xilogd(s):
    return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
def up(y, v):
    if v == 0: v = mp.mpf("1e-10")
    return mp.re(xilogd(mp.mpf(1)/2 + y + 1j*v))/mp.pi
def P(a, x):
    return a/(mp.pi*(a*a + x*x))
def uK(K, y, v):
    sK = sum(1/gk**2 for gk in g[:K])
    beta = VarX/2 - sK
    return 2*beta*y/mp.pi + sum(P(y, v-gk) + P(y, v+gk) for gk in g[:K]), beta
out = {}
for y in [mp.mpf('0.5'), mp.mpf('0.1')]:
    for K in K_list:
        vmax = float(g[K-1]) + 10
        n = 400
        rows = []
        worst_in = 0; worst_rel_in = 0
        for i in range(n+1):
            v = mp.mpf(vmax)*i/n
            a = up(y, v); b, beta = uK(K, y, v)
            rows.append((float(v), float(a), float(b)))
            if v <= g[K-1]/2:
                worst_in = max(worst_in, abs(a-b)); worst_rel_in = max(worst_rel_in, abs(a-b)/abs(a))
        a0 = up(y, 0); b0, beta = uK(K, y, 0)
        out[f"y={float(y)},K={K}"] = dict(beta=float(beta), u_p0=float(a0), u_K0=float(b0),
            max_abs_err_v_le_half_gK=float(worst_in), max_rel_err_v_le_half_gK=float(worst_rel_in),
            gK=float(g[K-1]))
        print(f"y={float(y)} K={K}: beta_K={float(beta):.6e} u_p(0)={float(a0):.8f} u_K(0)={float(b0):.8f} "
              f"max|diff| on [0,gK/2]={float(worst_in):.3e} rel={float(worst_rel_in):.3e} gK={float(g[K-1]):.4f}")
        json.dump(rows, open(f"rows_y{float(y)}_K{K}.json", "w"))
json.dump(out, open("thorin_compare.json", "w"), indent=1)
