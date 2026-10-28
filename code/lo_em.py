"""lo_em.py: Chapter 30 (H = chapter_level_one.tex), Section 4: the EM inversion table (H:240 to 247) and the redesigned injection test (reconcile-H section 2.11, Q3).

What it computes. Job table: EM (Shepp and Vardi) on the real u_1, data on [100, 300], unknown on [80, 320], step 0.05, flat start (total mass matched), rows
(theta', steps) = (0.75, 300), (0.5, 300), (0.5, 3000); for each the KL residual f = sum_i [y_i log(y_i/(Ax)_i) - y_i + (Ax)_i] (plain sum over the data grid; the printed
value is compared in four normalisations: f, f/m, f*h, f/(m h)), the number of peaks of x on (100, 300), the largest distance of a zero to the nearest peak, min x (cell density and cell mass h x).
Job inject: the quartet test. Two consecutive line zeros gamma_a < gamma_b nearest the window centre are removed and the quartet rho, 1-conj(rho), conj(rho), 1-rho, with
rho = sigma + i gamma_m, gamma_m = (gamma_a + gamma_b)/2, takes their place (positive ordinates: widths 1-sigma and sigma at gamma_m, mirror images at -gamma_m): the total mass of u_1 is unchanged.
The control replaces the same two zeros by two line zeros at gamma_m (width 1/2 each). Data: y = u_1(real) + exact Cauchy difference. (theta', sigma) = (0.5, 0.75) and (0.75, 0.85) (sigma > theta'),
plus the old boundary case (0.75, 0.75). Windows: data [100, 300] with unknown [80, 320], and data [150, 250] with unknown [130, 270]. EM steps 300, 1000, 3000, 10000;
the total residual f is reported at every count for the real data, the control and the quartet.
Claims supported: H:240 to 247, reconcile-H rows 35, 36 and Q3, Q4. Conventions of cx_em (the book's code, imported): unknowns are cell densities, A_ij = h P_{1-theta'}(v_i - w_j).
Error control: u_1 by flint balls at 128 bits rechecked at 192 bits (logged by cx_em.u1_half); the FFT product of cx_em equals the dense matrix to rounding (checked on a block here, logged).
Written from the mathematics (the collaborator's scripts were not read, run or copied).
Usage: python3 lo_em.py JOB [--quick] [--procs N] [--out FILE], JOB in table, inject
"""
import sys, os, json, math, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.signal import find_peaks
import lo_common as C
import cx_em

ap = argparse.ArgumentParser()
ap.add_argument('job'); ap.add_argument('--quick', action='store_true'); ap.add_argument('--procs', type=int, default=1)
ap.add_argument('--out', default=None)
a = ap.parse_args()
out_path = a.out or ('lo_em_%s.json' % a.job)
logf = open(out_path.replace('.json', '.log'), 'w')
def log(*s):
    m = ' '.join(str(x) for x in s); print(m, flush=True); logf.write(m + '\n'); logf.flush()

H = 0.05
Z = C.zeros()
VMAX = 120 if a.quick else 320
u = cx_em.u1_half(VMAX, H, a.procs, log)
v = np.arange(len(u)) * H

def window(lo, hi, ulo, uhi):
    return int(round(lo/H)), int(round(hi/H)), int(round(ulo/H)), int(round(uhi/H))

def peaks_stats(P, x, lo, hi):
    w = P.w; sel = (w > lo) & (w < hi)
    pk, _ = find_peaks(x[sel])
    zs = Z[(Z > lo) & (Z < hi)]
    md, nu = C.match(w[sel][pk], zs)
    pk2, _ = find_peaks(x[sel], prominence=0.01*x[sel].max())
    md2, nu2 = C.match(w[sel][pk2], zs)
    return dict(peaks=int(len(pk)), zeros=int(len(zs)), max_err=md, zeros_with_own_peak=nu, peaks_prominence_1pct=int(len(pk2)), max_err_those=md2, zeros_with_own_peak_those=nu2)

res = {}
if a.job == 'table':
    lo, hi, ulo, uhi = (20, 100, 10, 110) if a.quick else (100, 300, 80, 320)
    a0, a1, b0, b1 = window(lo, hi, ulo, uhi)
    # dense versus FFT check on a block
    th = 0.5
    Pt = cx_em.Prob(a0, a0+200, b0, b0+300, H, th, u[a0:a0+201])
    x = np.abs(np.random.RandomState(1).randn(Pt.n)) + 0.1
    Ad = H * C.P(1-th, (np.arange(Pt.m)[:, None] + a0)*H - (np.arange(Pt.n)[None, :] + b0)*H)
    log('FFT against dense product, max relative difference: %.2e' % (np.abs(Ad @ x - Pt.A(x)).max() / np.abs(Ad @ x).max()))
    for th, steps in [(0.75, 300), (0.5, 300), (0.5, 3000)]:
        P = cx_em.Prob(a0, a1, b0, b1, H, th, u[a0:a1+1])
        hist, x = cx_em.run_em(P, steps)
        f = hist[-1]; m = P.m
        d = dict(theta=th, steps=steps, f=float(f), f_per_point=float(f/m), f_times_h=float(f*H), f_per_length=float(f/(hi-lo)),
                 min_x_density=float(x.min()), min_x_mass=float(x.min()*H), min_x_density_in_data_window=float(x[(P.w >= lo) & (P.w <= hi)].min()),
                 f_at=[float(hist[k]) for k in (100, 300) if k <= steps])
        d.update(peaks_stats(P, x, lo, hi))
        # conjugate check by a second run at doubled step on the same window
        res['%g_%d' % (th, steps)] = d
        log(d)
    # grid check: step 0.1 (data every second sample), residual per unit length
    for th, steps in [(0.75, 300)]:
        Hh = 0.1
        a0h, a1h, b0h, b1h = int(round(lo/Hh)), int(round(hi/Hh)), int(round(ulo/Hh)), int(round(uhi/Hh))
        uh = u[::2]
        Ph = cx_em.Prob(a0h, a1h, b0h, b1h, Hh, th, uh[a0h:a1h+1])
        hist, x = cx_em.run_em(Ph, steps)
        res['grid0.1_%g_%d' % (th, steps)] = dict(f=float(hist[-1]), f_per_length=float(hist[-1]/(hi-lo)))
        log('step 0.1 grid, theta', th, 'steps', steps, 'f', hist[-1], 'f per unit length', hist[-1]/(hi-lo))
elif a.job == 'inject':
    cases = [(0.5, 0.75), (0.75, 0.85), (0.75, 0.75)]
    wins = [((100, 300, 80, 320) if not a.quick else (20, 100, 10, 110)), ((150, 250, 130, 270) if not a.quick else (40, 80, 30, 90))]
    STEPS = [300, 1000, 3000, 10000] if not a.quick else [100, 300]
    for (lo, hi, ulo, uhi) in wins:
        a0, a1, b0, b1 = window(lo, hi, ulo, uhi)
        vd = v[a0:a1+1]
        c = 0.5*(lo+hi)
        ia = int(np.argmin(np.abs(Z - c)))
        # two consecutive zeros nearest c: pair (ia, ia+1) or (ia-1, ia), the closer pair midpoint
        j = ia if abs(0.5*(Z[ia]+Z[ia+1]) - c) < abs(0.5*(Z[ia-1]+Z[ia]) - c) else ia-1
        ga, gb = Z[j], Z[j+1]; gm = 0.5*(ga+gb)
        removed = sum(C.P(0.5, vd - s*g) for g in (ga, gb) for s in (1, -1))
        for (thp, sig) in cases:
            quartet = sum(C.P(w_, vd - s*gm) for w_ in (1-sig, sig) for s in (1, -1))
            control = sum(C.P(0.5, vd - s*gm) for _ in range(2) for s in (1, -1))
            datasets = dict(real=u[a0:a1+1], control=u[a0:a1+1] - removed + control, quartet=u[a0:a1+1] - removed + quartet)
            key = 'win%d_%d|theta=%g|sigma=%g' % (lo, hi, thp, sig)
            r = dict(gamma_a=float(ga), gamma_b=float(gb), gamma_m=float(gm), window=[lo, hi, ulo, uhi], mass_diff=float((datasets['quartet']-datasets['real']).sum()*H))
            for name, y in datasets.items():
                P = cx_em.Prob(a0, a1, b0, b1, H, thp, y)
                hist, x = cx_em.run_em(P, STEPS[-1])
                r[name] = {str(k): float(hist[k]) for k in STEPS}
                # local misfit: sum of the residual terms in a window of +-3 around gamma_m
                Ax = P.A(x); loc = np.abs(vd - gm) < 3
                r[name]['local_misfit_pm3'] = float(np.sum(y[loc]*np.log(y[loc]/Ax[loc]) - y[loc] + Ax[loc]))
            r['ratio_quartet_over_control'] = {k: r['quartet'][k]/r['control'][k] for k in r['quartet']}
            res[key] = r
            log(key, 'zeros %.4f %.4f midpoint %.4f' % (ga, gb, gm))
            for name in ('real', 'control', 'quartet'): log('   ', name, {k: ('%.3e' % t) for k, t in r[name].items()})
            log('    ratio quartet/control', {k: '%.2f' % t for k, t in r['ratio_quartet_over_control'].items()})
json.dump(res, open(out_path, 'w'), indent=1)
log('done')
