# Local phase theory for the slice threshold. G = Xi with gamma_k, gamma_{k+1} removed (real zeros only).
# At w near eta: Theta(u) - B = atan(u/eps) - atan(u/(2eta-eps)) - b u, B = arg G(m+iw), b = -Im(G'/G)(m+iw).
# The rise covers a level pi/2 + k pi iff dist(B, pi Z) >= 2 sqrt(eps c) with c = b + 1/(2 eta); so
# eps* = dist(B(m+i w*), pi Z)^2 / (4 c).  Leading order: eps* ~ (G'/G(m))^2 eta^3 / 2.
import mpmath as mp, json, re, os, subprocess, sys
mp.mp.dps = 20
Gm = [mp.im(mp.zetazero(j)) for j in range(1, 40)]
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
def dlogxi(s): return 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
surv = {}
# Continuation thresholds: run nc_cauchy_survey.py (same directory) for k = 1,3,5,8,11,15 and parse its "k=.. m=.. eta=..: w*=.." lines.
here = os.path.dirname(os.path.abspath(__file__))
out = subprocess.run([sys.executable, '-I', os.path.join(here, 'nc_cauchy_survey.py'), '1,3,5,8,11,15'], capture_output=True, text=True, check=True).stdout
for line in out.splitlines():
    mt = re.match(r"k=(\d+) m=([\d.]+) eta=([\d.]+): w\*=([\d.]+)", line)
    if mt:
        k, m, e, w = mt.groups()
        surv[(int(k), float(e))] = float(w)
rows = []
for k in [1, 3, 5, 8, 11, 15]:
    gk, gk1 = Gm[k-1], Gm[k]; m = (gk+gk1)/2
    for eta in [0.1, 0.2, 0.3, 0.4]:
      c0 = m + 1j*eta
      # G now includes the mirror pair -c, -conj(c) of the planted zeros (all factors of F except the pair at m)
      def logG_d(z, c0=c0):
        return 1j*dlogxi(mp.mpf(1)/2 + 1j*z) - (2*z/(z*z-gk*gk) + 2*z/(z*z-gk1*gk1)) + 1/(z+c0) + 1/(z+mp.conj(c0))
      def argG(z, c0=c0):
        return mp.arg(xi(mp.mpf(1)/2 + 1j*z)/((z*z-gk*gk)*(z*z-gk1*gk1))*(z+c0)*(z+mp.conj(c0)))
      a1 = mp.re(logG_d(m))
      if True:
        eps = mp.mpf(0)
        for _ in range(30):   # fixed point for eps
            w = eta - eps; z = m + 1j*w
            B = argG(z); d = abs(B - mp.pi*mp.nint(B/mp.pi))
            b = -mp.im(logG_d(z)); c = b + 1/(2*eta - eps)
            eps = d*d/(4*c)
        lead = a1**2*eta**3/2
        obs = eta - surv[(k, eta)]
        rows.append((k, float(m), eta, float(obs), float(eps), float(lead)))
        print(f"k={k:2d} eta={eta}: observed gap {obs:.3e}  phase theory {float(eps):.3e}  leading (G'/G)^2 eta^3/2 = {float(lead):.3e}")
json.dump(rows, open('gap_theory2.json', 'w'))
