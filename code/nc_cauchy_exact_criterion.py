# Exact local criterion: real zeros of A_w near m <=> exists u with
#   Phi(u) = atan(2 w u / (u^2 + eta^2 - w^2)) + arg G(m+u+iw)  in  pi/2 + pi Z,   G = F/((z-c)(z-cbar)).
# Threshold w*: the least w at which max_u |Phi - nearest level| reaches 0. Compare with continuation (survey).
import mpmath as mp, json, re, os, subprocess, sys
mp.mp.dps = 25
Gm = [mp.im(mp.zetazero(j)) for j in range(1, 20)]
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
surv = {}
# Continuation thresholds: run nc_cauchy_survey.py (same directory) for k = 1,3,8,15 and parse its "k=.. m=.. eta=..: w*=.." lines.
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, '-I', os.path.join(here, 'nc_cauchy_survey.py'), '1,3,8,15'], capture_output=True, text=True, check=True).stdout
for line in out.splitlines():
    mt = re.match(r"k=(\d+) m=([\d.]+) eta=([\d.]+): w\*=([\d.]+)", line)
    if mt:
        k, mm, e, w = mt.groups()
        surv[(int(k), float(e))] = float(w)
rows = []
for k in [1, 3, 8, 15]:
    gk, gk1 = Gm[k-1], Gm[k]; m = (gk+gk1)/2
    for eta in [0.1, 0.3]:
        c = m + 1j*mp.mpf(eta)
        def argG(x, w):
            z = x + 1j*w
            return mp.arg(xi(mp.mpf(1)/2 + 1j*z)/((z*z-gk*gk)*(z*z-gk1*gk1))*(z+c)*(z+mp.conj(c)))
        def gap_to_level(w):
            # psi(m) sets the side; scan u on a fine grid near the expected extremum and return max over u of
            # (Phi(u) - (pi/2)) on the side where psi(m) > 0, or (-pi/2 - Phi(u)) where psi(m) < 0; >= 0 means real zeros.
            psi0 = argG(m, w); psi0 = psi0 - mp.pi*mp.nint(psi0/mp.pi); sgn = 1 if psi0 > 0 else -1
            u0 = mp.sqrt(max(eta**2 - w**2, mp.mpf(10)**-30))
            best = -mp.inf
            for t in mp.linspace(0.2, 3, 57):
                u = sgn*u0*t
                ps = argG(m+u, w); ps = ps - mp.pi*mp.nint((ps - psi0)/mp.pi)
                Phi = mp.atan(2*w*u/(u*u + eta**2 - w**2)) + ps
                best = max(best, sgn*Phi - mp.pi/2)
            return best
        lo, hi = mp.mpf(eta)*0.9, mp.mpf(eta)*(1-1e-12)
        for _ in range(45):
            mid = (lo+hi)/2
            if gap_to_level(mid) >= 0: hi = mid
            else: lo = mid
        wstar = (lo+hi)/2
        rows.append(dict(k=k, eta=eta, w_exact=float(wstar), w_cont=surv[(k, eta)]))
        print(f"k={k} eta={eta}: threshold from exact criterion {float(wstar):.6f}   by continuation {surv[(k, eta)]:.5f}", flush=True)
json.dump(rows, open('exact_criterion.json', 'w'))
