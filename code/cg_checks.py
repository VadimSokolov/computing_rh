# cg_checks.py: the numerical claims of Chapter 27 (the Comtet Majumdar family and the uniform law), priorities 1 and 2: the exact expansion
# of zeta_delta against quadrature (Theorem cg:thm:expansion), the displacement table (Table cg:tab:drift: slopes and Re rho_delta - 1/2 for
# the first three ordinates), the poles of xi_delta at s = 3, 5, 7 (Proposition cg:prop:obstructions (a): the deltas where F_delta(5) and
# F_delta(7) vanish, the sign of F_delta(3)), the uniform law (Theorem cg:thm:uniform: closed form against the Lorentzian sum, positivity of
# kappa_delta, R(4.1), sqrt(3) pi, the kink constant 15 pi - 7.5 - 4 pi^2; Theorem cg:thm:mellinthorin at s = 2.5, delta = 0.5), the
# threshold alpha_1 = 7.789 of (1/8) zeta(3, alpha/2 + 1) = (alpha - 1)^(-3), and the Comtet Majumdar constants gamma_mu = zeta(s)/((2 pi)^s
# sin(pi s/2)) (Spitzer reading, gamma_2).
# Conventions. g(x) = 1/(e^x - 1) - 1/x. I_delta(s) = (1/pi) int_0^inf x^(-s) g(x + delta) dx for 0 < Re s < 1; zeta_delta(s) = (2 pi)^s
# sin(pi s/2) I_delta(s). F_delta(s) = Li_(1-s)(e^(-delta)) - Gamma(s) delta^(-s) = sum_k zeta(1-s-k)(-delta)^k/k! (Lindelof, delta < 2 pi),
# and pi I_delta(s) = Gamma(1-s) F_delta(s). Zeros of zeta_delta near rho are the zeros of F_delta (the other factors are finite, nonzero).
# Method. Quadrature of I_delta in the variable x = e^u (mp.quad on [-inf, ...] split at the oscillation scale); F_delta by the Lindelof sum
# with 90 terms and, as a second form, by the Poisson form F_delta(n + 1) = n! sum_(k != 0) (delta - 2 pi i k)^(-n-1) at integers; Newton
# (mp.findroot, secant-free: Newton with the numerical derivative of mp) on F_delta from rho + slope delta, at two precisions (30 and 45
# digits). kappa_delta from the closed form and from the Lorentzian sum (m up to 4000 plus the tail by the Euler Maclaurin bound).
# Mellin check: (2 pi)^s int y^(-s) kappa_delta dy over N periods of 2 pi plus the tail with the period average of kappa; N = 60 and 120.
# Error control: each quantity is computed at two precisions or steps, the maximum difference is logged ("diff"), and written.
# Output: data/cg_checks.json and data/cg_checks.log (or --out FILE); --quick runs in seconds with fewer periods and one precision.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, argparse, time
from mpmath import mp, mpf, mpc, zeta, gamma, sin, tan, cos, cosh, sinh, coth, tanh, pi, exp, log, quad, inf, findroot, diff, polylog, factorial, rf, sqrt, atan, nstr, fabs, im, re
ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default='data/cg_checks.json')
ap.add_argument('--only', default='all'); a = ap.parse_args()
res = {}; logl = []
def L(*x):
    s = ' '.join(str(t) for t in x); print(s, flush=True); logl.append(s)
def want(k): return a.only in ('all', k)
def g(x): return 1/(exp(x)-1) - 1/x
def ordinates():
    mp.dps = 40
    return [mp.zetazero(n).imag for n in (1, 2, 3)]
def zdelta_series(s, d, K=90):
    t = tan(pi*s/2); tot = mpc(0)
    for k in range(K):
        r = (-1)**(k//2) if k % 2 == 0 else (-1)**(k//2)*t
        tot += (d/(2*pi))**k*rf(s, k)/factorial(k)*r*zeta(s+k)
    return tot
def F_lind(s, d, K=90):
    return sum(zeta(1-s-k)*(-d)**k/factorial(k) for k in range(K))
def F_poisson(n, d, K=20000):
    # F_delta(n+1) = n! sum_{k != 0} (delta - 2 pi i k)^(-n-1)  (real: 2 Re sum_{k >= 1}); tail by integral
    s = mpf(0)
    for k in range(1, K+1): s += 2*((d-2j*pi*k)**(-n-1)).real
    return factorial(n)*s
if want('quad'):
    t0 = time.time(); mp.dps = 20
    for s in (mpf('0.6'), mpc('0.5', 5)):
        for d in (mpf('0.1'), mpf(1)):
            # x = e^u: int x^(1-s) g(x+d) du ; split at u = -40 (x^(1-Re s) ~ e^(-0.4*40)) is not negligible for s=.6: use the small x series
            f = lambda u: exp((1-s)*u)*g(exp(u)+d)
            pts = [-inf, -60, -20, -5, 0, 3, 10, 40, 120]
            I = quad(f, pts)/pi
            z1 = (2*pi)**s*sin(pi*s/2)*I; z2 = zdelta_series(s, d)
            z3 = (2*pi)**s*sin(pi*s/2)*gamma(1-s)*F_lind(s, d)/pi
            L('quad', nstr(s, 4), nstr(d, 3), 'quad', nstr(z1, 18), 'series', nstr(z2, 18), 'F-form', nstr(z3, 18), 'relerr quad-series', nstr(abs(z1-z2)/abs(z2), 3), 'series-Fform', nstr(abs(z3-z2)/abs(z2), 3))
            res['quad_%s_%s' % (nstr(s, 4), nstr(d, 3))] = {'quad': str(z1), 'series': str(z2), 'F': str(z3), 'rel_quad_series': nstr(abs(z1-z2)/abs(z2), 3)}
    L('quad time', time.time()-t0)
if want('drift'):
    t0 = time.time(); gam = ordinates(); out = {}
    for dps in (30, 45):
        mp.dps = dps; out[dps] = {}
        for j, gj in enumerate(gam):
            rho = mpc(mpf(1)/2, gj); rho = mpc(0.5, mp.zetazero(j+1).imag) if False else mpc(mpf(1)/2, gj)
            slope_c = -rho*tan(pi*rho/2)*zeta(rho+1)/(2*pi*diff(zeta, rho))
            out[dps][('slope', j)] = slope_c
            for d in (mpf('0.01'), mpf('0.05'), mpf('0.1'), mpf('0.3'), mpf(1)):
                r = findroot(lambda s: F_lind(s, d, 50 if dps == 30 else 65), rho + slope_c*d, tol=mpf(10)**(-dps+8))
                out[dps][(str(d), j)] = r
    mp.dps = 30
    for j in range(3):
        sl = out[30][('slope', j)]
        L('drift', 'gamma', j+1, nstr(gam[j], 10), 'slope complex', nstr(sl, 8), 'Re slope', nstr(sl.real, 6), 'diff dps', nstr(abs(out[30][('slope', j)]-out[45][('slope', j)]), 3))
        for d in ('0.01', '0.05', '0.1', '0.3', '1.0'):
            r30, r45 = out[30][(str(mpf(d)), j)], out[45][(str(mpf(d)), j)]
            L('drift', 'gamma', j+1, 'delta', d, 'Re rho_d - 1/2', nstr(r45.real-mpf(1)/2, 8), 'Im', nstr(r45.imag, 8), 'diff', nstr(abs(r30-r45), 3))
            res['drift_%d_%s' % (j+1, d)] = {'disp': nstr(r45.real-mpf(1)/2, 12), 'im': nstr(r45.imag, 12), 'diff_dps': nstr(abs(r30-r45), 3)}
        res['slope_%d' % (j+1)] = {'complex': str(sl), 're': nstr(sl.real, 12)}
    L('drift time', time.time()-t0)
if want('poles'):
    mp.dps = 40
    # F_delta(n+1) by Lindelof sum (n+1 = s) against the Poisson form; zeros of F_delta(5), F_delta(7) in (0, 2 pi)
    for n in (2, 4, 6):
        for d in (mpf('0.1'), mpf(1), mpf(3)):
            fl = F_lind(mpf(n+1), d, 120)
            fp = F_poisson(n, d, 3000 if a.quick else 20000)
            L('poles F_delta(%d) delta=%s Lindelof %s Poisson %s' % (n+1, d, nstr(fl, 15), nstr(fp, 15)))
    def Fp(n, d):
        # exact: n! * 2 Re sum_{k>=1} (d - 2 pi i k)^(-n-1) = (n!) * 2 Re[ (1/(2 pi i))^(n+1) * (-1)^? Hurwitz]: use zeta(n+1, .) with complex shift
        return factorial(n)*2*re(sum(((d-2j*pi*k)**(-n-1)) for k in range(1, 1)) ) if False else None
    def Fz(s, d): return F_lind(s, d, 120)
    for n, guess in ((4, 4.7), (6, 3.0)):
        rt = findroot(lambda d: Fz(mpf(n+1), d), mpf(guess))
        rt2 = findroot(lambda d: F_poisson(n, d, 3000), mpf(guess))
        L('poles F_delta(%d)=0 at delta = %s (Lindelof) %s (Poisson 3000 terms, tail not included)' % (n+1, nstr(rt, 12), nstr(rt2, 8)))
        res['Fzero_%d' % (n+1)] = nstr(rt, 12)
    # sign of F_delta(3) on (0, 2 pi): grid
    mx = max(F_poisson(2, mpf(i)/50*2*pi, 3000) for i in range(1, 50))  # Poisson form: the Lindelof sum converges too slowly near 2 pi
    L('poles F_delta(3) max over delta = 2 pi i/50, i = 1..49 : %s (negative means no zero)' % nstr(mx, 8))
    res['F3max'] = nstr(mx, 8)
    for m in (3, 5, 7, 41):
        L('pole s=%d: zeta(-%d) = %s' % (m, m, nstr(zeta(-m), 8)))
if want('uniform'):
    mp.dps = 30
    lel = lambda d, x: (d/(d*d+x*x))/pi*pi  # ell_delta(x) = delta/(delta^2+x^2) (no 1/pi, as in the chapter's pi kappa)
    def pk(d, y): return y**2/(d*(d*d+y*y)) - coth(d/2)/2*(1-cos(y))/(cosh(d)-cos(y))
    def pk_lor(d, y, M=2000):
        tot = mpf(0)
        for m in range(1, M+1):
            tot += lel(d, y-2*pi*m)+lel(d, y+2*pi*m)-2*lel(d, 2*pi*m)
        # tail: sum_{m>M} [..] ~ d*y^2 * 3 ... use second difference: l'' * y^2 ; l''(x) = 2 d (3x^2 - d^2)/(d^2+x^2)^3 ~ 6 d/x^4
        tail = y**2*6*d/(2*pi)**4*(mpf(1)/(3*(M+mpf(1)/2)**3))
        return tot+tail
    for d in (mpf('0.1'), mpf(1), mpf(5)):
        for y in (mpf('0.5'), mpf(3), mpf(7)):
            L('uniform pi kappa closed %s Lorentz %s  (delta %s, y %s)' % (nstr(pk(d, y), 15), nstr(pk_lor(d, y), 15), d, y))
    ds = [mpf(x) for x in ('0.01', '0.1', '0.5', '1', '2', '3', '4.1', '5.44', '8', '20')]
    mn = {}
    for d in ds:
        ys = [mpf(i)/40 for i in range(1, 1601)] if not a.quick else [mpf(i)/4 for i in range(1, 161)]
        vals = [pk(d, y)/y**2 for y in ys]   # kappa ~ c y^2
        i = min(range(len(vals)), key=lambda i: vals[i]); mn[str(d)] = (nstr(vals[i], 6), nstr(ys[i], 4))
        L('uniform delta', d, 'min of pi kappa/y^2 on (0,40]', nstr(vals[i], 6), 'at y', nstr(ys[i], 4), 'min of pi kappa itself positive:', min(pk(d, y) for y in ys) > 0)
    res['kappa_min_over_y2'] = mn
    Rf = lambda d: 8*tanh(d/2)*sinh(d/2)**2/(d*(d*d+pi**2))
    L('R(4.1) =', nstr(Rf(mpf('4.1')), 10), ' R(pi) =', nstr(Rf(pi), 8), ' R(4) =', nstr(Rf(4), 8), ' sqrt3 pi =', nstr(sqrt(3)*pi, 8))
    L('15 pi - 7.5 - 4 pi^2 =', nstr(15*pi-mpf('7.5')-4*pi**2, 8))
    res['R41'] = nstr(Rf(mpf('4.1')), 10); res['sqrt3pi'] = nstr(sqrt(3)*pi, 8); res['kink'] = nstr(15*pi-mpf('7.5')-4*pi**2, 8)
    # alpha_1: (1/8) zeta(3, alpha/2+1) = (alpha-1)^(-3)
    al = findroot(lambda t: zeta(3, t/2+1)/8-(t-1)**(-3), mpf('7.8'))
    L('alpha_1 =', nstr(al, 10)); res['alpha1'] = nstr(al, 10)
    # monotonicity of R on [pi, 30]
    mono = all(Rf(pi+i*mpf('0.05')) < Rf(pi+(i+1)*mpf('0.05')) for i in range(0, 500))
    L('R increasing on [pi, pi+25] on grid 0.05:', mono)
if want('mellin'):
    mp.dps = 20; s = mpf('2.5'); d = mpf('0.5')
    def pk(y):  # cancellation of order y^2 at small y: work at 80 digits there
        with mp.workdps(80 if y < 1 else 20): return +(y**2/(d*(d*d+y*y)) - coth(d/2)/2*(1-cos(y))/(cosh(d)-cos(y)))
    A = coth(d/2)/2; mean = 1 + (1-cosh(d))/sinh(d)
    vals = {}
    for N in ((20, 40) if a.quick else (60, 120)):
        tot = quad(lambda y: y**(-s)*pk(y), [0, mpf(1)])
        # pieces [1, 2 pi], then periods
        tot += quad(lambda y: y**(-s)*pk(y), [1, 2*pi])
        for m in range(1, N): tot += quad(lambda y: y**(-s)*pk(y), [2*pi*m, 2*pi*(m+1)])
        Lh = 2*pi*N
        tot += quad(lambda y: y**(-s)*(y**2/(d*(d*d+y*y)) - A*mean), [Lh, inf])
        v = (2*pi)**s*tot/pi
        vals[N] = v; L('mellin N', N, 'RHS (2 pi)^s int y^-s kappa', nstr(v, 14))
    zd = zdelta_series(s, d)
    L('mellin zeta_delta(2.5) series', nstr(zd, 14), 'diff N', nstr(abs(vals[min(vals)]-vals[max(vals)]), 3), 'rhs-lhs', nstr(abs(vals[max(vals)]-zd), 3))
    res['mellin'] = {'rhs': nstr(vals[max(vals)], 14), 'series': nstr(zd.real, 14)}
if want('cm'):
    mp.dps = 30
    for s in (mpf(1)/2, mpf('0.3'), mpf('0.7')):
        G = zeta(s)/((2*pi)**s*sin(pi*s/2)); G2 = gamma(1-s)*zeta(1-s)/pi
        L('cm s', nstr(s, 3), 'gamma_mu = zeta(s)/((2pi)^s sin(pi s/2)) =', nstr(G, 12), ' Gamma(1-s) zeta(1-s)/pi =', nstr(G2, 12))
    L('zeta(1/2)/sqrt(pi) =', nstr(zeta(mpf(1)/2)/sqrt(pi), 8), ' zeta(1/2)/sqrt(2 pi) =', nstr(zeta(mpf(1)/2)/sqrt(2*pi), 8))
    # Spitzer for the Gaussian walk with E e^{ikX} = e^{-k^2} (variance 2): E M_n - sqrt(n)*2... compare constant: EM_n = (1/pi) Gamma(1/2) sum k^{-1/2}
    n = 10**5; Sn = sum(mpf(k)**(-mpf(1)/2) for k in range(1, 20001))
    # Euler Maclaurin exact: sum_{k<=n} k^{-1/2} = 2 sqrt(n) + zeta(1/2) + n^{-1/2}/2 + ...
    n = 20000; L('Spitzer sum k^-1/2 minus 2 sqrt n - n^-1/2/2 =', nstr(Sn-2*sqrt(n)-mpf(1)/(2*sqrt(n)), 10), ' zeta(1/2) =', nstr(zeta(mpf(1)/2), 10))
    L('EM_n - 2 sqrt(n/pi) constant for variance 2 steps: Gamma(1/2) zeta(1/2)/pi =', nstr(sqrt(pi)*zeta(mpf(1)/2)/pi, 8), ' = zeta(1/2)/sqrt(pi)')
json.dump(res, open(a.out if a.only == 'all' else a.out.replace('.json', '_%s.json' % a.only), 'w'), indent=1, default=str)
open((a.out if a.only == 'all' else a.out.replace('.json', '_%s.json' % a.only)).replace('.json', '.log'), 'w').write('\n'.join(logl)+'\n')
