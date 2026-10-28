"""Moment tilt of the Polya density by u^(2n): cumulants, zeros of the tilted transform, Laplace error (Ki chapter, G5).

Computes, for the density proportional to u^(2n) Phi(u) on R (symmetric, U), with mu=E|U|, sigma^2=Var|U|, kappa=E(|U|-mu)^3:
  the first six positive zeros z_k of F(z)=int u^(2n) Phi(u) cos(zu) du, the shifts z_k mu - c_k (c_k=(2k-1)pi/2),
  the Edgeworth prediction kappa c_k^3/(6 mu^3), the sup over 0<=s<=20 of |psi_n(s)-1/cosh sqrt(2s)|,
  psi_n(s)=1/E cosh(sqrt(2s)|U|/mu), its ratio to sigma^2/mu^2, the saddle point formulas, max x^2/(2 cosh x).
Supports the tables and numbers of chapter_ki.tex section 'Numerical evidence' and the Conjecture on the Edgeworth rate.
Conventions: --norm book uses Phi of ch/ch02.tex (e^(9u/2), e^(2u), factor 2, Xi(theta)=xi(1/2+i theta)); --norm chapter uses
the chapter's printed Phi (e^(9u), e^(4u), no factor 2, which equals Phi_book(2u)/2). z mu, sigma^2/mu^2, kappa/mu^3 and
the Laplace error are the same in both, mu, sigma, kappa, z scale by 2.
Method: trapezoid rule on [lo,hi] around the mode (spectrally accurate, integrand analytic and negligible at the ends; at lo=0
it vanishes to order 2n), mpmath at dps digits, Phi(u>=0) from the first 12 terms, weights from log Phi to avoid underflow,
zeros by secant iteration from c_k/mu. Error control: every n is run twice (dps, window, step h differ); the max relative
difference of all outputs is logged.
Options: --quick (n=10,100 only), --out FILE (json), --norm book|chapter.
Written from the mathematics (the collaborator's scripts were not read, run or copied).
"""
import sys, json, argparse
from mpmath import mp, mpf, pi, exp, log, cos, cosh, sqrt, findroot, mpf as F, sin

def run(n, norm, dps, W, hdiv):
    mp.dps = dps
    a = 4 if norm == 'chapter' else 2
    fac = 1 if norm == 'chapter' else 2
    def logphi(u):
        t = 0
        for m in range(1, 13):
            t += fac*(2*pi**2*m**4*exp(mpf(9)*a/4*u) - 3*pi*m**2*exp(mpf(5)*a/4*u))*exp(-pi*m**2*exp(a*u))
        return log(t)
    # mode of 2n log u + leading term
    g = lambda u: 2*n/u + mpf(9)*a/4 - a*pi*exp(a*u)
    u0 = findroot(g, mpf(a and 1)*(log(n)/a if n > 1 else 1) if n > 1 else 1)
    s0 = 1/sqrt(2*n/u0**2 + a*a*pi*exp(a*u0))
    lo = max(mpf(0), u0 - W*s0); hi = u0 + W*s0
    N = int((hi-lo)/(s0/hdiv)) + 1
    h = (hi-lo)/N
    us = [lo + i*h for i in range(N+1)]
    lw = [2*n*log(u) + logphi(u) if u > 0 else None for u in us]
    mx = max(x for x in lw if x is not None)
    w = [exp(x-mx) if x is not None else mpf(0) for x in lw]
    tot = sum(w)
    mu = sum(wi*u for wi, u in zip(w, us))/tot
    var = sum(wi*(u-mu)**2 for wi, u in zip(w, us))/tot
    kap = sum(wi*(u-mu)**3 for wi, u in zip(w, us))/tot
    Fz = lambda z: sum(wi*cos(z*u) for wi, u in zip(w, us))
    res = dict(n=n, mu=mu, sigma=sqrt(var), kappa=kap, mode=u0, N=N)
    sh = []; ed = []; zs = []
    for k in range(1, 7):
        c = (2*k-1)*pi/2
        pred = kap*c**3/(6*mu**3)
        z = findroot(Fz, (c/mu*(1+pred/c)*(1+mpf(10)**-8), c/mu*(1-mpf(10)**-8)), solver='secant', tol=mpf(10)**(-dps+6), maxsteps=100)
        zs.append(z); sh.append(z*mu-c); ed.append(pred)
    res['shift'] = sh; res['edgeworth'] = ed; res['atoms'] = [(z*mu)**2/2 for z in zs]
    # Laplace error on [0,20]
    best = mpf(0); arg = None
    M = 400
    for j in range(1, M+1):
        s = mpf(20)*j/M; x = sqrt(2*s)
        e = abs(1/(sum(wi*cosh(x*u/mu) for wi, u in zip(w, us))/tot) - 1/cosh(x))
        if e > best: best, arg = e, s
    res['lap_sup'] = best; res['lap_arg_s'] = arg; res['ratio'] = best/(var/mu**2)
    return res

def flat(r):
    out = {}
    for k, v in r.items():
        if isinstance(v, list):
            for i, x in enumerate(v): out[f'{k}{i+1}'] = x
        else: out[k] = v
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None)
    ap.add_argument('--norm', default='book'); ap.add_argument('--n', default=None)
    A = ap.parse_args()
    ns = [10, 100] if A.quick else [5, 10, 100, 1000, 10000, 100000, 1000000]
    if A.n: ns = [int(x) for x in A.n.split(',')]
    out = []
    for n in ns:
        r1 = run(n, A.norm, 40, 14, 4)
        r2 = run(n, A.norm, 50, 18, 6)
        f1, f2 = flat(r1), flat(r2)
        d = max(abs(f1[k]-f2[k])/abs(f1[k]) for k in f1 if k not in ('N', 'lap_arg_s') and f1[k] != 0)
        print(f'n={n} norm={A.norm} mu={float(r1["mu"]):.6g} sigma={float(r1["sigma"]):.6g} kappa={float(r1["kappa"]):.6g} ratio={float(r1["ratio"]):.6g} maxreldiff_two_runs={float(d):.2e}', flush=True)
        print('  shift ', [f'{float(x):.4e}' for x in r1['shift']]); print('  edgew ', [f'{float(x):.4e}' for x in r1['edgeworth']])
        print('  atoms ', [f'{float(x):.6f}' for x in r1['atoms']])
        r1['reldiff'] = d
        out.append({k: ([float(y) for y in v] if isinstance(v, list) else float(v)) for k, v in r1.items()})
    if A.out: json.dump(out, open(A.out, 'w'), indent=1)
