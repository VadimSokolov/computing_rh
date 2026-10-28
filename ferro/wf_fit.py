# Note W (Chapter 13, Section fm:sec:bessel, Table fm:tab:bessel), computation N1: the members of the family of
# Theorem fm:thm:bessel closest to the Polya law p = Phi/Xi(0), found by numerical minimisation of
#   TV(q, p) = int_0^infty |q - p|        (both densities are even; this is the total variation distance),
# in the book's units (t, a, l, lam, mu): q(t) = m(t) exp(-a cosh(l t)) / C, m(t) = exp(mu t^2) prod_k cosh(lam_k t),
# C = (2/l) E K_{Y/l}(a), Y = sqrt(2 mu) Z + sum_k lam_k eps_k.  Phi is the series (2.2.3) of the book,
#   Phi(x) = 2 sum_n (2 pi^2 n^4 e^{9x/2} - 3 pi n^2 e^{5x/2}) exp(-pi n^2 e^{2x}),  n <= 6, evaluated at |x|.
# Rows (the rows of Table fm:tab:bessel): 1 m = 1, free (a, l); 2 one field, free (a, l, lam); 3 a = 2 pi, l = 2, one
# field, free lam; 4 a = 2 pi, l = 2, cosh^2(lam t), free lam; 5 a = 2 pi, l = 2, exp(mu t^2) cosh(lam t), free (lam, mu).
# Exploratory (reported, not tabulated): two free fields with free kernel; one field with a tilt and a free kernel.
# Method: Nelder-Mead (scipy) on the double precision objective, Simpson's rule on [0, 3] with step 1e-4 (the tail
# beyond 3 of q and p is below 1e-300), several starts for the exploratory searches. Every optimum is then rounded to
# seven significant digits (exact decimals) and the distance AT THE DECIMALS is recomputed; it must be within 1% of the
# optimum (otherwise more digits are used). Error control: the normaliser is computed twice (Simpson on the same grid
# and scipy kv for mu = 0, scipy quad for mu > 0) and the objective once more with step 2e-4; the differences are logged.
# NOT certified: the certified distances of these exact decimals are in wf_tvcert.py.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 -I wf_fit.py [--quick] ; writes wf_fit.json (and prints the log; redirect to wf_fit.log)
import json, sys, time
import numpy as np, mpmath as mp
from scipy.optimize import minimize
from scipy.special import kv
from scipy.integrate import quad

QUICK = '--quick' in sys.argv
H = 1e-4 if not QUICK else 1e-3
mp.mp.dps = 30
XI0 = float(mp.mpf(1)/2*(mp.mpf(1)/2 - 1)/2*mp.pi**(-mp.mpf(1)/4)*mp.gamma(mp.mpf(1)/4)*mp.zeta(mp.mpf(1)/2))


def phi_series(t):
    t = np.abs(np.asarray(t, float)); s = 0.0
    for n in range(1, 7):
        s = s + 2*(2*np.pi**2*n**4*np.exp(4.5*t) - 3*np.pi*n**2*np.exp(2.5*t))*np.exp(-np.pi*n**2*np.exp(2*t))
    return s


def grid(h):
    n = int(round(3/h)); assert n % 2 == 0
    t = np.linspace(0, 3, n + 1); w = np.ones(n + 1); w[1:-1:2] = 4; w[2:-1:2] = 2
    return t, w*h/3


def logm(t, lams, mu):
    s = mu*t**2
    for l in lams:
        s = s + np.logaddexp(l*t, -l*t) - np.log(2)
    return s


class Obj:
    def __init__(self, h):
        self.t, self.w = grid(h); self.p = phi_series(self.t)/XI0

    def density(self, a, l, lams, mu):
        lq = logm(self.t, lams, mu) - a*np.cosh(l*self.t)
        f = np.exp(lq)
        C = 2*np.dot(self.w, f)                  # integral over the line
        return f/C, C

    def tv(self, a, l, lams, mu):
        q, _ = self.density(a, l, lams, mu)
        return float(np.dot(self.w, np.abs(q - self.p)))


def norm_check(a, l, lams, mu):
    """Relative difference of the normaliser from Simpson and from kv (mu = 0) or quad (mu > 0)."""
    _, C = OB.density(a, l, lams, mu)
    if mu == 0:
        ys = np.array([0.0])
        for lam in lams:
            ys = np.concatenate([ys + lam, ys - lam])
        C2 = 2/l*np.mean(kv(ys/l, a))
    else:
        f = lambda t: np.exp(logm(t, lams, mu) - a*np.cosh(l*t))
        C2 = 2*quad(f, 0, 4, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
    return abs(C/C2 - 1)


def make(kind):
    """returns (unpack: x -> (a, l, lams, mu), start list)"""
    two_pi = 2*np.pi
    if kind == 'row1':
        return (lambda x: (abs(x[0]), abs(x[1]), [], 0.0)), [[2.839, 2.5745]]
    if kind == 'row2':
        return (lambda x: (abs(x[0]), abs(x[1]), [abs(x[2])], 0.0)), [[5.483990, 2.062861, 2.146192]]
    if kind == 'row3':
        return (lambda x: (two_pi, 2.0, [abs(x[0])], 0.0)), [[2.6185]]
    if kind == 'row4':
        return (lambda x: (two_pi, 2.0, [abs(x[0])]*2, 0.0)), [[1.7925]]
    if kind == 'row5':
        return (lambda x: (two_pi, 2.0, [abs(x[0])], x[1]**2)), [[2.1255, np.sqrt(0.94175)]]
    if kind == 'two_fields':
        return (lambda x: (abs(x[0]), abs(x[1]), [abs(x[2]), abs(x[3])], 0.0)), \
            [[5.483990, 2.062861, 2.146192, l2] for l2 in (0.5, 0.2, 1.0, 1.5)]
    if kind == 'tilt_free_kernel':
        return (lambda x: (abs(x[0]), abs(x[1]), [abs(x[2])], x[3]**2)), \
            [[5.483990, 2.062861, 2.146192, np.sqrt(m)] for m in (0.1, 0.02, 0.3, 0.6)]


def run(kind):
    unpack, starts = make(kind)
    best = None
    for x0 in starts:
        f = lambda x: OB.tv(*unpack(x))
        x = np.array(x0, float)
        for _ in range(4):                       # restarts of Nelder-Mead from the last point
            r = minimize(f, x, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-14, maxiter=6000, maxfev=6000,
                                                                   initial_simplex=None))
            x = r.x
        print('  %s start %s -> %.8g at %s' % (kind, np.round(x0, 4).tolist(), r.fun, np.round(x, 6).tolist()), flush=True)
        if best is None or r.fun < best[0]:
            best = (r.fun, x)
    fopt, x = best
    a0, l0, lams0, mu0 = unpack(x)
    fixed_kernel = kind in ('row3', 'row4', 'row5')
    r7 = lambda v: float('%.7g' % v)
    a, l = (a0, l0) if fixed_kernel else (r7(a0), r7(l0))
    lams = [r7(v) for v in lams0]; mu = r7(mu0) if mu0 else 0.0
    exact = dict(a='2pi' if fixed_kernel else '%.7g' % a0, l='2' if fixed_kernel else '%.7g' % l0,
                 lams=['%.7g' % v for v in lams0], mu='%.7g' % mu0 if mu0 else '0')
    tvr = OB.tv(a, l, lams, mu)
    tv2 = OB2.tv(a, l, lams, mu)
    nc = norm_check(a, l, lams, mu)
    out = dict(kind=kind, tv_opt=fopt, params_opt=x.tolist(), exact=exact,
               a=a, l=l, lams=lams, mu=mu, tv_at_decimals=tvr, ratio=tvr/fopt, within_1pct=bool(tvr <= 1.01*fopt),
               tv_step_2h=tv2, tv_step_diff=abs(tv2 - tvr), normaliser_rel_diff=nc)
    print('%s: optimum %.8g; at the decimals %s: %.8g (ratio %.5f); step doubled %.3g apart; normaliser %.2e'
          % (kind, fopt, exact, tvr, tvr/fopt, abs(tv2 - tvr), nc), flush=True)
    return out


if __name__ == '__main__':
    t0 = time.time()
    OB = Obj(H); OB2 = Obj(2*H)
    print('Xi(0) = %.15g (mpmath, 30 digits); Simpson step %g' % (XI0, H))
    which = sys.argv[1:] and [a for a in sys.argv[1:] if not a.startswith('--')] or \
        ['row1', 'row2', 'row3', 'row4', 'row5', 'two_fields', 'tilt_free_kernel']
    R = {}
    for k in which:
        R[k] = run(k)
        if not QUICK:
            json.dump(R, open('wf_fit.json', 'w'), indent=1)
        print('   (%.0f s)' % (time.time() - t0), flush=True)
    if QUICK:
        json.dump(R, open('/tmp/wf_fit_quick.json', 'w'), indent=1)
