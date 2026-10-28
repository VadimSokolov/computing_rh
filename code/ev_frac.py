# ev_frac.py: the cotangent law for fractional shifts (E:452-478, Problem 2): zeros nu_k of Re K_{mu + i nu}(a) at a = 2 pi, the Robin ratios
# beta_k = a K'_{i nu_k}(a)/K_{i nu_k}(a), their slope against sqrt(lambda_k) = nu_k over the last zeros below nu = 45 for mu = 1/2, 3/4, 5/4, 1/4 (chapter: 1.0007,
# .4101, -.4170, 3.04 against cot(pi mu/2) = 1, .4142, -.4142, 2.41), and the three parameter fit beta_k = .37 - .417 sqrt(lambda_k) + 2.43/sqrt(lambda_k) at mu = 5/4
# over the first 24 zeros (relative accuracy 5e-4).
# Conventions. K'_x(a) = d/da K_x(a) = -(K_{x-1}(a) + K_{x+1}(a))/2 with x = i nu; zeros are real roots in nu of the real function Re K_{mu + i nu}(a).
# The chapter does not say how many "last zeros" enter the slope, nor whether the fit is least squares in beta or in relative error: slopes are given by least squares
# over the last m zeros below 45 for m = 6, 8, 10, 12 and the fit for least squares in beta and relative least squares, both reported.
# Method. mpmath besselk at 50 digits (K_{i nu} is of size exp(-pi nu/2), about 10^-31 at nu = 45, so 25 digits are not enough), sign scan in nu with step 0.02 and secant polish. Error control: zeros recomputed at 70 digits (only for mu = 5/4), max difference logged.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Usage: python3 -I ev_frac.py [--out FILE] [--only 1/2|3/4|5/4|1/4] [--dps D]  (one shift about 15 s at 50 digits)
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, mpc, besselk, findroot, pi, nstr, cot, matrix, lu_solve, sqrt
import argparse
ap = argparse.ArgumentParser(); ap.add_argument('--out', default=None); ap.add_argument('--only', default=None); ap.add_argument('--dps', type=int, default=50); args = ap.parse_args()
logf = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); logf.append(s)
t0 = time.time(); a = 2 * pi
def zeros_nu(mu, dps, nmax=45):
    mp.dps = dps
    f = lambda nu: besselk(mu + 1j * nu, a).real
    grid = [mpf(j) / 50 for j in range(1, 50 * nmax)]
    fv = [f(x) for x in grid]; z = []
    for j in range(len(grid) - 1):
        if fv[j] * fv[j + 1] < 0:
            z.append(findroot(f, (grid[j], grid[j + 1]), solver='anderson', tol=mpf(10) ** (-(dps - 6))))
    return z
def beta(nu):
    x = 1j * nu
    return (a * (-(besselk(x - 1, a) + besselk(x + 1, a)) / 2) / besselk(x, a)).real
def slope(nu, b, m):
    nu, b = nu[-m:], b[-m:]; n = len(nu); mx = sum(nu) / n; my = sum(b) / n
    return sum((x - mx) * (y - my) for x, y in zip(nu, b)) / sum((x - mx) ** 2 for x in nu)
res = {}
for mu, name in ((mpf(1) / 2, '1/2'), (mpf(3) / 4, '3/4'), (mpf(5) / 4, '5/4'), (mpf(1) / 4, '1/4')):
    if args.only and args.only != name: continue
    z = zeros_nu(mu, args.dps); b = [beta(x) for x in z]
    P('mu =', name, 'zeros below 45:', len(z), 'cot(pi mu/2) =', nstr(cot(pi * mu / 2), 5), 'slopes over the last m zeros, m = 6, 8, 10, 12:', [nstr(slope(z, b, m), 5) for m in (6, 8, 10, 12)])
    res[name] = [nstr(slope(z, b, m), 5) for m in (6, 8, 10, 12)]
    if name == '5/4':
        mp.dps = args.dps
        z24 = z[:24]; b24 = b[:24]
        for tag, wt in (('plain LS', lambda bb: 1), ('relative LS', lambda bb: 1 / bb)):
            A = matrix(24, 3); y = matrix(24, 1)
            for i in range(24):
                w = wt(b24[i]); A[i, 0] = w; A[i, 1] = w * sqrt(z24[i] ** 2); A[i, 2] = w / z24[i]; y[i] = w * b24[i]
            sol = lu_solve(A.T * A, A.T * y)
            err = max(abs(sol[0] + sol[1] * z24[i] + sol[2] / z24[i] - b24[i]) / abs(b24[i]) for i in range(24))
            P('   fit', tag, ': beta = %s + %s nu + %s/nu, max relative error over 24 zeros %s' % (nstr(sol[0], 4), nstr(sol[1], 4), nstr(sol[2], 4), nstr(err, 3)))
            res['fit_' + tag] = [nstr(sol[i], 5) for i in range(3)] + [nstr(err, 3)]
        for tag, cf in (('chapter coefficients', (mpf('0.37'), mpf('-0.417'), mpf('2.43'))),):
            e = max(abs(cf[0] + cf[1] * z24[i] + cf[2] / z24[i] - b24[i]) / abs(b24[i]) for i in range(24))
            P('   max relative error of the chapter coefficients over 24 zeros:', nstr(e, 3))
        zh = zeros_nu(mu, args.dps + 20)
        P('   zeros at two precisions, max difference:', nstr(max(abs(x - y) for x, y in zip(z, zh)), 3))
P('time', round(time.time() - t0, 1))
of = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'ev_frac.json')
json.dump(res, open(of, 'w'), indent=1); open(of.replace('.json', '.log'), 'w').write('\n'.join(logf) + '\n')
