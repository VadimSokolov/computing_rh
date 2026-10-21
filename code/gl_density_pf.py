# gl_density_pf.py: the numerical claims of Chapter 21 (Thorin Representations, Dirichlet L Functions and the Linnik Constant) about the
# density of Y = sqrt(2H) Z and about the real axis values of 1/xi: Section gl:sec:pfnum (the density f, its value at 0 and at 2, the log
# concavity grid, the minors of order 2 to 5), the Stirling numbers after Lemma gl:lem:stirling, and the smooth model and the elasticities
# of E H^(-1/2) in Problem gl:op:negmoment. Added for the book in October 2026, written from the mathematics (the collaborator's scripts
# were not read, run or copied).
# Conventions. F(t) = xi(1/2)/xi(1/2 + t) = E exp(-t^2 H) is the characteristic function of Y, and f(y) = (1/pi) int_0^infty cos(t y) F(t) dt.
# Method for f. F is even and analytic in the strip |Im t| < gamma_1 = 14.13 (its poles are at t = +-i gamma), so the trapezoid rule of step
# D on the real line converges geometrically: by Poisson summation the rule returns f(y) + f(2 pi/D - y) + f(2 pi/D + y) + ..., and
# f(y) is of order 400 exp(-gamma_1 y), so the aliasing error at y <= 4 is below 1e-50 for D = 1/2 (1e-62 at y <= 2) and below 1e-127 for D = 1/4. Hence
# f(y) = (D/pi)(1/2 + sum_{n >= 1} F(n D) cos(n D y)), with the sum stopped when F(n D) < 10^(-dps - 10) (F decreases for t > 0); f' and
# f'' are the termwise derivatives. Every quantity is computed for two steps D (D = 1/2 and 1/4 in the quick case, 1/4 and 1/8 in the full
# case) and the maximum relative difference of every family (f on the grid, the log concavity numbers, the minors) is logged and written:
# that difference is the error of every grid. The step of every grid in the chapter's tests is exact (decimals j/20 as rationals).
# Further checks of f. mp.quad of the integral at y = 0 and y = 2 (30 digits); 2 sqrt(pi) f(0) against E H^(-1/2) as a Mellin integral;
# f(y) against the residue series f(y) = sum_k r_k exp(-gamma_k y), r_k = -xi(1/2)/Xi'(gamma_k) with Xi(t) = xi(1/2 + i t), which is
# valid (absolutely convergent) for y > pi/4 because |r_k| grows like exp(pi gamma_k/4); at y = 2 this is the chapter's remark that f(2)
# is the first residue term to within the next one, and the residues use the first K ordinates of data/zeros_1700.txt (68 digits).
# Minors. det[f(x_r - y_c)] of order 3 for the 121 pairs (i, j), the Toeplitz minors det[f((r - c) h)] of order 3 to 8 for h = 0.05, 0.2 and
# 0.4 (the chapter prints orders 3 to 5), all positive or not as computed; minors of size 1e-16 are computed at 40 to 70 digits.
# The smooth model of E H^(-1/2) (Problem gl:op:negmoment). E_m replaces, in the product over the zeros, all zeros above a cut T by the
# density theta'(t)/pi, with the first m zeros exact. The text does not say where the cut is for m >= 1, so three versions are computed:
# gram, the cut at the root of theta(T) = (m - 1) pi (T_0 = 9.6669081 for m = 0, the point where N_0 + 1 = 0, and in general the point
# where the smooth count N_0 + 1 equals m); mid, the cut at the midpoint of gamma_m and gamma_(m+1); anch, the midpoint cut with the
# atom of mass N_0(T) + 1 - m at T that the anchored tail of equation gl:eq:anchored adds. Double precision (scipy), checked by the
# error of the same quadrature on the exact integrand. Elasticities: (gamma_k/E) dE/dgamma_k = 2 [1 - int F(t)/(1 + t^2/gamma_k^2) dt/int F(t) dt]
# with E = (2/sqrt(pi)) int_0^infty F(t) dt, which is the derivative of E H^(-1/2) = Gamma(1/2)^(-1) int u^(-1/2) prod (1 + u/gamma^2)^(-1) du
# with respect to the single ordinate gamma_k at fixed other ordinates; the sum over all k is 1 (scaling), checked over the 1700 ordinates.
# Stirling numbers. |Gamma(s) xi(1/2)/xi(1/2 + s)| at s = 0.4 + i tau for several tau, with |zeta(0.9 + i tau)| and the ratio to the
# bound exp(-pi tau/4) tau^(sigma/2 - 9/4) of Lemma gl:lem:stirling.
# Output: data/gl_density_pf.json (numbers as strings, parameters) and the log gl_density_pf.log next to it.
# Usage: python3 gl_density_pf.py [--quick] [--out FILE] [--zeros FILE] [--dps D]; --quick runs everything at 40 digits (steps 1/2 and 1/4) in
# about 4 seconds and writes to /tmp unless --out is given; the full case (70 digits, steps 1/4 and 1/8, data/gl_density_pf.json) takes about
# 8 seconds in one process, so it fits the local limit too; on the cluster: bash misc/tools/hopper_run.sh -c 1 -t 30 -g gl_density_pf.json
# -g gl_density_pf.log code/gl_density_pf.py data/zeros_1700.txt -- --out gl_density_pf.json
import argparse
import json
import math
import sys
import tempfile
import time
import warnings
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import integrate, special

HERE = Path(__file__).resolve().parent
T_START = time.time()
LOGFILE = None
HALF = mp.mpf(1) / 2
XIH = None
DPS = 40


def log(*a):
    line = '[%6.1fs] %s' % (time.time() - T_START, ' '.join(str(x) for x in a))
    print(line, flush=True)
    if LOGFILE is not None:
        LOGFILE.write(line + '\n')
        LOGFILE.flush()


def fs(x, d=25):
    return mp.nstr(x, d)


def rel(a, b):
    return abs(a - b) / abs(b)


def xi(s):
    s = mp.mpmathify(s)
    if s == 1 or s == 0:
        return HALF
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def Fval(t):
    return XIH / xi(HALF + t)


def theta(t):
    return mp.im(mp.loggamma(mp.mpc(mp.mpf(1) / 4, t / 2))) - t / 2 * mp.log(mp.pi)


def F_grid(delta):
    """F(n delta), n = 0, 1, ..., until F < 10^(-dps - 10)."""
    vals = [mp.mpf(1)]
    tol = mp.mpf(10) ** (-(DPS + 10))
    n = 1
    while True:
        v = Fval(n * delta)
        vals.append(v)
        if v < tol:
            return vals
        n += 1


def dens3(y, F, delta):
    """f(y), f'(y), f''(y) from the cosine series with the F grid."""
    a0 = a1 = a2 = mp.mpf(0)
    for n, Fn in enumerate(F):
        w = n * delta
        wt = Fn / 2 if n == 0 else Fn
        c, s = mp.cos_sin(w * y)
        a0 += wt * c
        a1 -= wt * w * s
        a2 -= wt * w * w * c
    k = delta / mp.pi
    return a0 * k, a1 * k, a2 * k


def det_toeplitz(fv, n, h_idx):
    """det[f(|r - c| h)] with f given on the grid j/20 and h = h_idx/20."""
    return mp.det(mp.matrix([[fv[abs(r - c) * h_idx] for c in range(n)] for r in range(n)]))


def det3_shift(fv, d_idx):
    """det[f(x_r - y_c)] with x_r - y_c = (d + r - c)/10, d = i - j: the argument is |d + r - c|/10 = 2 |d + r - c|/20."""
    return mp.det(mp.matrix([[fv[2 * abs(d_idx + r - c)] for c in range(3)] for r in range(3)]))


def evaluate_grid(delta):
    F = F_grid(delta)
    fv, f1, f2 = [], [], []
    for j in range(81):
        a, b, c = dens3(mp.mpf(j) / 20, F, delta)
        fv.append(a)
        f1.append(b)
        f2.append(c)
    phi = [fv[j] * f2[j] - f1[j] ** 2 for j in range(81)]
    minors3 = {d: det3_shift(fv, d) for d in range(-10, 11)}
    toep = {}
    for h_idx, hname in ((1, '0.05'), (4, '0.2'), (8, '0.4')):
        for n in range(3, 9):
            if (n - 1) * h_idx <= 80:
                toep['h=%s,n=%d' % (hname, n)] = det_toeplitz(fv, n, h_idx)
    return dict(F=F, f=fv, f1=f1, f2=f2, phi=phi, minors3=minors3, toep=toep)


def gl_panels(edges, n):
    x, w = np.polynomial.legendre.leggauss(n)
    pts, wts = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        pts.append(0.5 * (b - a) * x + 0.5 * (b + a))
        wts.append(0.5 * (b - a) * w)
    return np.concatenate(pts), np.concatenate(wts)


def main():
    global LOGFILE, XIH, DPS
    ap = argparse.ArgumentParser(description='density of Y, Polya frequency minors, Stirling numbers and the smooth model of E H^(-1/2) (Chapter 21)')
    ap.add_argument('--quick', action='store_true', help='40 digits, steps 1/2 and 1/4, writes to /tmp unless --out is given')
    ap.add_argument('--out', default=None, help='JSON path (default data/gl_density_pf.json; quick: a /tmp file)')
    ap.add_argument('--zeros', default=None, help='path of zeros_1700.txt')
    ap.add_argument('--dps', type=int, default=None)
    a = ap.parse_args()
    quick = a.quick
    DPS = a.dps or (40 if quick else 70)
    out = Path(a.out) if a.out else (Path(tempfile.gettempdir()) / 'gl_density_pf_quick.json' if quick else HERE.parent / 'data' / 'gl_density_pf.json')
    LOGFILE = open(out.with_suffix('.log'), 'w')
    warnings.simplefilter('ignore')
    mp.mp.dps = DPS
    cands = [Path(a.zeros)] if a.zeros else []
    cands += [HERE.parent / 'data' / 'zeros_1700.txt', Path('data/zeros_1700.txt'), Path('zeros_1700.txt'), HERE / 'zeros_1700.txt']
    zf = next((c for c in cands if c.exists()), None)
    if zf is None:
        sys.exit('zeros_1700.txt not found, pass --zeros')
    gam = [mp.mpf(x) for x in zf.read_text().split()]
    XIH = xi(HALF)
    R = dict(script='gl_density_pf.py', mode='quick' if quick else 'full')
    deltas = (mp.mpf(1) / 2, mp.mpf(1) / 4) if quick else (mp.mpf(1) / 4, mp.mpf(1) / 8)
    log('dps', DPS, 'steps', [fs(d, 4) for d in deltas], 'zeros', len(gam))

    # ---- the density on the grid j/20, j = 0..80, for two steps
    G = [evaluate_grid(d) for d in deltas]
    A, B = G
    log('F grid sizes', [len(g['F']) for g in G], 'last values', [fs(g['F'][-1], 3) for g in G])
    dev_f2 = max(rel(A['f'][j], B['f'][j]) for j in range(41))
    dev_f4 = max(rel(A['f'][j], B['f'][j]) for j in range(81))
    dev_phi2 = max(rel(A['phi'][j], B['phi'][j]) for j in range(41))
    dev_phi4 = max(rel(A['phi'][j], B['phi'][j]) for j in range(81))
    dev_m3 = max(rel(A['minors3'][d], B['minors3'][d]) for d in A['minors3'])
    dev_toep = max(rel(A['toep'][k], B['toep'][k]) for k in A['toep'])
    log('two steps agree: f on the chapter grid y <= 2 to relative %s (y <= 4: %s), phi = f f\'\' - f\'^2 for y <= 2 to %s (y <= 4: %s), the order 3 minors to %s, the Toeplitz minors to %s'
        % (fs(dev_f2, 3), fs(dev_f4, 3), fs(dev_phi2, 3), fs(dev_phi4, 3), fs(dev_m3, 3), fs(dev_toep, 3)))
    f = B['f']
    R['grid_error'] = dict(steps=[fs(d, 4) for d in deltas], f_rel_y_le_2=fs(dev_f2, 3), f_rel_y_le_4=fs(dev_f4, 3), phi_rel_y_le_2=fs(dev_phi2, 3), phi_rel_y_le_4=fs(dev_phi4, 3),
                           minors3_rel=fs(dev_m3, 3), toeplitz_rel=fs(dev_toep, 3), aliasing_estimate_at_y_4=fs(400 * mp.exp(-gam[0] * (2 * mp.pi / deltas[0] - 4)), 3))

    # ---- f(0), f(2), quadrature check, relation to E H^(-1/2)
    cache = {}

    def Fc(t):
        v = cache.get(t)
        if v is None:
            v = cache[t] = Fval(t)
        return v

    with mp.workdps(30):
        tmax = 2
        while Fval(mp.mpf(tmax)) > mp.mpf(10) ** (-42):
            tmax += 1
        pts = [mp.mpf(i) for i in range(tmax + 1)]
        f0q = mp.quad(Fc, pts) / mp.pi
        f2q = mp.quad(lambda t: mp.cos(2 * t) * Fc(t), pts) / mp.pi
    two_sqrt_pi_f0 = 2 * mp.sqrt(mp.pi) * f[0]
    V = 2
    while Fval(mp.mpf(V) ** 2) > mp.mpf(10) ** (-(DPS + 12)):
        V += 1
    ehm12 = 4 * mp.quad(lambda v: v * Fval(v * v), [mp.mpf(i) for i in range(V + 1)]) / mp.gamma(HALF)   # E H^(-1/2), u = v^4
    log('f(0) =', fs(f[0], 20), '(quad %s, relative %s); f(2) = %s (quad %s, relative %s)' % (fs(f0q, 15), fs(rel(f[0], f0q), 3), fs(f[40], 15), fs(f2q, 12), fs(rel(f[40], f2q), 3)))
    log('2 sqrt(pi) f(0) =', fs(two_sqrt_pi_f0, 20), ' E H^(-1/2) (Mellin, u = v^4) =', fs(ehm12, 20), ' relative difference', fs(rel(two_sqrt_pi_f0, ehm12), 3))
    R['f0_f2'] = dict(f0=fs(f[0], 30), f0_quad_30_digits=fs(f0q, 25), f2=fs(f[40], 30), f2_quad_30_digits=fs(f2q, 25), two_sqrt_pi_f0=fs(two_sqrt_pi_f0, 30), E_H_minus_half=fs(ehm12, 30),
                      rel_f0_quad=fs(rel(f[0], f0q), 3), rel_f2_quad=fs(rel(f[40], f2q), 3), rel_f0_relation=fs(rel(two_sqrt_pi_f0, ehm12), 3))

    # ---- residue series
    K = 30 if quick else 60
    resid = []
    for k in range(K):
        rho = mp.mpc(HALF, gam[k])
        xip = rho * (rho - 1) / 2 * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1)
        rk = 1j * XIH / xip
        resid.append((mp.re(rk), abs(mp.im(rk)) / abs(rk)))
    r1 = resid[0][0]
    log('residue weights r_k = -xi(1/2)/Xi\'(gamma_k): r_1 =', fs(r1, 12), 'r_2 =', fs(resid[1][0], 12), 'r_3 =', fs(resid[2][0], 12), '; imaginary parts relative below', fs(max(x[1] for x in resid), 3))
    res_rows = {}
    for y in (2, 1.5, 1, mp.mpf('0.9')):
        y = mp.mpf(y)
        j = int(round(float(y * 20)))
        fy = f[j]
        terms = [resid[k][0] * mp.exp(-gam[k] * y) for k in range(K)]
        one = terms[0]
        full = sum(terms)
        res_rows[fs(y, 3)] = dict(f=fs(fy, 20), first_residue=fs(one, 12), first_residue_rel=fs(rel(one, fy), 3), K_terms=K, residue_sum=fs(full, 20), residue_sum_rel=fs(rel(full, fy), 3))
        log('y = %s: f = %s; first residue term %s (relative difference %s); sum of %d residue terms %s (relative difference %s)'
            % (fs(y, 3), fs(fy, 12), fs(one, 10), fs(rel(one, fy), 3), K, fs(full, 12), fs(rel(full, fy), 3)))
    R['residues'] = dict(r=[fs(resid[k][0], 12) for k in range(6)], rows=res_rows, imag_rel_max=fs(max(x[1] for x in resid), 3))

    # ---- log concavity grid (C:436)
    phi = B['phi']
    jmax = max(range(41), key=lambda j: phi[j])
    jmax4 = max(range(81), key=lambda j: phi[j])
    neg2 = all(phi[j] < 0 for j in range(41))
    neg4 = all(phi[j] < 0 for j in range(81))
    log('phi = f f\'\' - f\'^2 on [0, 2] step 0.05: all negative:', neg2, '; largest value', fs(phi[jmax], 6), 'at y =', fs(mp.mpf(jmax) / 20, 4), '; f at that point', fs(f[jmax], 6))
    log('phi on [0, 4] step 0.05: all negative:', neg4, '; largest', fs(phi[jmax4], 6), 'at y =', fs(mp.mpf(jmax4) / 20, 4), '; phi(0) =', fs(phi[0], 8), 'phi(1) =', fs(phi[20], 8), 'phi(2) =', fs(phi[40], 8))
    R['log_concavity'] = dict(grid='y = j/20, j = 0..40', all_negative=bool(neg2), max_phi=fs(phi[jmax], 10), argmax=fs(mp.mpf(jmax) / 20, 4), f_at_argmax=fs(f[jmax], 10),
                              extended_grid_0_4_all_negative=bool(neg4), extended_max=fs(phi[jmax4], 10), extended_argmax=fs(mp.mpf(jmax4) / 20, 4),
                              phi_values=[fs(phi[j], 10) for j in range(41)])

    # ---- minors of order 3 (121 pairs)
    m3 = B['minors3']
    allm = [(m3[i - j], i, j) for i in range(11) for j in range(11)]
    mn = min(allm)
    npos = sum(1 for x in allm if x[0] > 0)
    log('order 3 minors det[f(x_r - y_c)], 121 pairs (i, j): positive', npos, 'of', len(allm), '; minimum', fs(mn[0], 8), 'at (i, j) =', (mn[1], mn[2]), '; at d = i - j = 0:', fs(m3[0], 8), '; d = 5:', fs(m3[5], 8), '; d = 10:', fs(m3[10], 8))
    R['minors3'] = dict(count=len(allm), positive=npos, minimum=fs(mn[0], 10), argmin=[mn[1], mn[2]], by_difference={str(d): fs(m3[d], 10) for d in range(0, 11)})

    # ---- Toeplitz minors
    tp = B['toep']
    R['toeplitz'] = {k: fs(v, 10) for k, v in tp.items()}
    for hname in ('0.05', '0.2', '0.4'):
        log('Toeplitz minors det[f((r - c) h)] at h = %s, n = 3..8: %s' % (hname, [fs(tp['h=%s,n=%d' % (hname, n)], 8) for n in range(3, 9) if 'h=%s,n=%d' % (hname, n) in tp]))
    log('smallest Toeplitz minor over the family:', fs(min(tp.values()), 6), '; all positive:', all(v > 0 for v in tp.values()))

    # ---- Stirling numbers (C:464)
    st = {}
    sig = mp.mpf('0.4')
    for tau in (20, 40, 80, 120, 160, 200):
        s = mp.mpc(sig, tau)
        val = abs(mp.gamma(s) * XIH / xi(HALF + s))
        zt = abs(mp.zeta(HALF + s))
        bound = mp.exp(-mp.pi * tau / 4) * mp.mpf(tau) ** (sig / 2 - mp.mpf(9) / 4)
        st['tau=%d' % tau] = dict(value=fs(val, 8), abs_zeta=fs(zt, 8), value_times_abs_zeta=fs(val * zt, 8), over_bound=fs(val / bound, 6))
        log('|Gamma(s) xi(1/2)/xi(1/2 + s)| at s = 0.4 + %di: %s; |zeta(0.9 + %di)| = %s; value times |zeta| = %s (the Stirling factor); ratio to exp(-pi tau/4) tau^(-2.05) = %s'
            % (tau, fs(val, 6), tau, fs(zt, 6), fs(val * zt, 6), fs(val / bound, 6)))
    R['stirling'] = st

    # ---- smooth model and elasticities of E H^(-1/2) (C:802), double precision
    edges = [0, 2, 5, 10, 20, 40, 80]
    tn, wn = gl_panels(edges, 40)
    Fn = np.array([float(Fval(mp.mpf(float(t)))) for t in tn])
    E_exact = 2 / math.sqrt(math.pi) * float(np.sum(wn * Fn))
    log('E H^(-1/2) by Gauss Legendre on %d nodes: %.15f; Mellin value %s' % (len(tn), E_exact, fs(ehm12, 16)))
    gamf = np.array([float(g) for g in gam[:60]])

    def dth(x):
        return (0.5 * special.digamma(0.25 + 0.5j * x).real - 0.5 * np.log(np.pi)) / np.pi

    def model_E(m, T, atom):
        tot = 0.0
        for t, w in zip(tn, wn):
            fin = float(np.sum(np.log1p(t * t / gamf[:m] ** 2))) if m > 0 else 0.0
            inn = integrate.quad(lambda x: math.log1p(t * t / (x * x)) * dth(x), T, np.inf, limit=200, epsabs=0, epsrel=1e-12)[0]
            tot += w * math.exp(-(fin + inn + atom * math.log1p(t * t / (T * T))))
        return 2 / math.sqrt(math.pi) * tot

    # self test of the model machinery: the inner integral at a large cut against the exact tail of the product
    sm = {}
    claimed = {0: '-0.00696', 1: '-0.002491', 2: '9.974e-5', 5: '-0.00039', 20: '1.369e-5'}
    for m in (0, 1, 2, 5, 20):
        if m == 0:
            cuts = [('gram', mp.findroot(lambda T: theta(T) + mp.pi, 9.7), 0.0)]
        else:
            mid = (gam[m - 1] + gam[m]) / 2
            Tg = mp.findroot(lambda T: theta(T) - (m - 1) * mp.pi, mid)
            atom = float(theta(mid) / mp.pi + 1 - m)
            cuts = [('gram', Tg, 0.0), ('mid', mid, 0.0), ('anch', mid, atom)]
        row = {}
        for name, T, atom in cuts:
            Em = model_E(m, float(T), atom)
            row[name] = dict(T=fs(T, 10), atom=fs(mp.mpf(atom), 5), E=fs(mp.mpf(Em), 12), rel_error=fs(mp.mpf(Em / E_exact - 1), 8))
        sm['m=%d' % m] = row
        log('smooth model, %d exact zeros (chapter %s): %s' % (m, claimed[m], {k: (v['T'], v['atom'], v['rel_error']) for k, v in row.items()}))
    R['smooth_model'] = dict(E_exact=fs(mp.mpf(E_exact), 14), models=sm, claimed=claimed)
    el = []
    for k in range(10):
        g2 = float(gam[k]) ** 2
        el.append(2 * (1 - float(np.sum(wn * Fn / (1 + tn * tn / g2))) / float(np.sum(wn * Fn))))
    allel = []
    for k in range(len(gam)):
        g2 = float(gam[k]) ** 2
        allel.append(2 * (1 - float(np.sum(wn * Fn / (1 + tn * tn / g2))) / float(np.sum(wn * Fn))))
    log('elasticities (gamma_k/E) dE/dgamma_k, k = 1..10:', ['%.6f' % x for x in el], '; sum over the 1700 ordinates %.6f (the sum over all is 1)' % sum(allel))
    R['elasticities'] = dict(first_ten=['%.8f' % x for x in el], sum_1700='%.8f' % sum(allel), claimed_first='0.179943', claimed_sixth='0.031864')

    R['parameters'] = dict(dps=DPS, steps=[fs(d, 4) for d in deltas], grid_points=len(B['F']), zeros_used=K, quick=quick, seconds=fs(time.time() - T_START, 4))
    json.dump(R, open(out, 'w'), indent=1)
    log('wrote', out)


if __name__ == '__main__':
    main()
