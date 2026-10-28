"""Checks of the closed form claims of the Ki chapter (G5): max x^2/(2cosh x), threshold S, saddle point formulas at n=1e4,
leading orders, 0.02 percent claim, rigidity claims, normalisation of Phi (Xi(0), first zero) and de Bruijn scaling.
Reads data/ki_tilt_chapter.json (from code/ki_tilt.py). Written from the mathematics (the collaborator's scripts were not read, run or copied)."""
import json, sys
from mpmath import mp, mpf, pi, exp, log, cosh, sinh, tanh, findroot, nsum, quad, cos, linspace, gamma, zeta, sqrt
mp.dps = 30
f = lambda x: x**2/(2*cosh(x))
xm = findroot(lambda x: 2*x/(2*cosh(x)) - x**2*sinh(x)/(2*cosh(x)**2), 2)
print('max x^2/(2cosh x)=', f(xm), 'at x=', xm, ' threshold S=xm^2/2=', xm**2/2)
D = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'data/ki_tilt_chapter.json'))
for r in D:
    n = r['n']
    g = lambda u: 2*n/u + 9 - 4*pi*exp(4*u)
    un = findroot(g, log(n)/4 + 0.5)
    s2 = 1/(8*n/un + 2*n/un**2)
    l3 = 4*n/un**3 - 64*pi*exp(4*un)
    ka = l3*s2**3
    print(f"n={n:g} u_n={float(un):.5f} (1/4)log n={float(log(n)/4):.5f} sigma2 quad={r['sigma']**2:.5e} saddle={float(s2):.5e} lead u/8n={float(un/(8*n)):.5e}; kappa quad={r['kappa']:.4e} saddle={float(ka):.4e} lead -u^2/16n^2={float(-un**2/(16*n**2)):.4e}")
    rel = max(abs(r['shift'][k]-r['edgeworth'][k])/abs(r['shift'][k]) for k in range(6))
    print(f"   max_k |shift-edgeworth|/|shift| = {rel*100:.4f} percent; max |shift|={max(abs(x) for x in r['shift']):.3e}; ratio={r['ratio']:.6f}")
# normalisation of the printed Phi
def Pch(u): u = abs(u); return nsum(lambda m: (2*pi**2*m**4*exp(9*u) - 3*pi*m**2*exp(5*u))*exp(-pi*m**2*exp(4*u)), [1, 12])
xi = lambda s: s*(s-1)/2*pi**(-s/2)*gamma(s/2)*zeta(s)
mp.dps = 15
for z in [0, 14.134725141734693, 28.269450283469386]:
    print('chapter Xi(z)=int Phi_ch e^{izu} at z=', z, 2*quad(lambda u: Pch(u)*cos(z*u), linspace(0, 3, 13)))
print('xi(1/2)/4 =', xi(mpf(1)/2)/4)
