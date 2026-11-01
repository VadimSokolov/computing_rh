# Two planted pairs. Remove gamma_k, gamma_{k+1} from Xi and plant c_j = m_j + i eta_j, j = 1, 2, at m_j = M -+ d/2,
# M the midpoint. For each pair track the zero of the tilt F(z+iw)+F(z-iw) born from c_j until it reaches the real axis
# (local slice threshold w*_j) and compare eta_j - w*_j with the one pair law (1/2) phi_j^2 eta_j^3 and with the fixed
# point of the gap proposition, where G_j = F / ((z-c_j)(z-cbar_j)) keeps the other planted pair. Floating point, mpmath.
import mpmath as mp, json, sys
mp.mp.dps = 20
G = [mp.im(mp.zetazero(j)) for j in range(1, 20)]
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)

def make(k, ms, etas):
    gk, gk1 = G[k-1], G[k]
    cs = [m + 1j*e for m, e in zip(ms, etas)]
    def F(z):
        v = xi(mp.mpf(1)/2 + 1j*z)/((z*z-gk*gk)*(z*z-gk1*gk1))
        for c in cs: v *= (z*z-c*c)*(z*z-mp.conj(c)**2)
        return v
    return F, cs

def threshold(F, c, eta):
    A = lambda z, w: F(z+1j*w) + F(z-1j*w)
    z = mp.mpc(c); w = mp.mpf(0); h = mp.mpf(eta)/20
    while h > mp.mpf(eta)*1e-9 and w < 1.05*eta:
        wn = w + h
        try: zn = mp.findroot(lambda u: A(u, wn), z)
        except Exception: h /= 2; continue
        if mp.im(zn) > 1e-12 and abs(zn - z) < 0.1: w, z = wn, zn
        else: h /= 2
    return w

def predictions(F, c, eta):
    m = mp.re(c)
    Gj = lambda z: F(z)/((z-c)*(z-mp.conj(c)))
    phi = mp.re(mp.diff(Gj, m)/Gj(m))
    lead = phi**2*eta**3/2
    # fixed point of the gap proposition
    eps = lead
    for _ in range(60):
        w = eta - eps
        g = Gj(m + 1j*w); a = mp.arg(g); dist = abs(a - mp.pi*mp.nint(a/mp.pi))
        b = -mp.im(mp.diff(Gj, m + 1j*w)/g)
        eps = dist**2/(4*(b + 1/(2*eta - eps)))
    return float(phi), float(lead), float(eps)

out = []
k = int(sys.argv[1]) if len(sys.argv) > 1 else 8
M = (G[k-1] + G[k])/2
for eta in [0.1, 0.2]:
    for d in [3.0, 2.0, 1.2, 0.8, 0.5, 0.3]:
        ms = [M - d/2, M + d/2]
        for etas in [(eta, eta), (eta, eta/2)]:
            F, cs = make(k, ms, [mp.mpf(e) for e in etas])
            for j in range(2):
                w = threshold(F, cs[j], mp.mpf(etas[j]))
                phi, lead, fp = predictions(F, cs[j], mp.mpf(etas[j]))
                gap = float(etas[j] - w)
                row = dict(k=k, d=d, eta=etas[j], other=etas[1-j], j=j, gap=gap, phi=phi, lead=lead, fixed=fp)
                out.append(row)
                print(f"d={d:4.1f} eta_j={etas[j]:.2f} other={etas[1-j]:.2f} j={j}: gap {gap:.3e} "
                      f"lead {lead:.3e} fixed {fp:.3e} phi {phi:+.3f}", flush=True)
json.dump(out, open(f"two_pairs_k{k}.json", "w"), indent=0)
