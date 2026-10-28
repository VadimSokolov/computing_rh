"""lo_flow.py: Chapter 30 (H = chapter_level_one.tex), Sections 1 to 3: the level one density, its masses, the spectral peaks, the band limited backward flow.

What it computes. Job u0: u_1(0) = (2 + gamma_E - log 4 pi)/(2 pi) against u_1(v) at small v by mpmath (30 digits) and u_1(0.05) - u_1(0). Job mass: the mass of u_1 on (0, V]
(V = 400, 300, 200, 100) three ways: Simpson on the grid of step 0.05, the arctan formula sum_{gamma>0} [atan((V-gamma)/(1-sigma)) + atan((V+gamma)/(1-sigma))]/pi
over the 1700 ordinates plus the smooth tail (Riemann von Mangoldt density), and counts N(V). Job spec: tab:spectrum, u_1 on [0, 400] step 0.05 extended by symmetry,
Hann window, FFT zero padded to step 1e-4 in omega, local maxima of |FT| near log n, heights relative to n = 3, against Lambda(n)/n. Job fejer: the table at H:185,
pipeline: 16001 samples v = -400..400, Hann window w, FFT, multiplier exp((1-theta')|omega|) for |omega| <= Omega (sharp cutoff; variant: times the Fejer factor 1 - |omega|/Omega),
inverse FFT, division by w, peaks and minimum on 100 < v < 300; rows with the analytic u_1(0) (the corrected run) and with u_1(0.05) in place of the sample at v = 0
(the variant that reproduces the printed -430 per reconcile-H section 0).
Claims supported: H:85 (mass 201.67, 202 zeros), H:119 to 132 (tab:spectrum), H:183 to 185 (backward flow table), H:190.
Error control: u_1 by flint balls at 128 bits (stencil 1e-5), rechecked at 192 bits and step 2e-5 on a sample (logged by cx_em); the mass is also checked by two quadratures.
Written from the mathematics (the collaborator's scripts were not read, run or copied).
Usage: python3 lo_flow.py [--quick] [--procs N] [--out FILE]
"""
import sys, os, json, math, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.signal import find_peaks
from scipy.integrate import simpson, quad, trapezoid
import lo_common as C

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true'); ap.add_argument('--procs', type=int, default=1)
ap.add_argument('--out', default='lo_flow.json')
a = ap.parse_args()
out = {}
logf = open(a.out.replace('.json', '.log'), 'w')
def log(*s):
    m = ' '.join(str(x) for x in s); print(m, flush=True); logf.write(m + '\n'); logf.flush()

VMAX = 40 if a.quick else 400
H = 0.05
Z = C.zeros()

# ---- u0
import mpmath as mp
mp.mp.dps = 30
def u1mp(v):
    s = mp.mpc(1, v)
    return float((1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)).real / mp.pi)
def u1mp_hi(v):
    mp.mp.dps = 50
    s = mp.mpc(1, v)
    r = (1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)).real / mp.pi
    mp.mp.dps = 30
    return r
exact0 = (2 + mp.euler - mp.log(4*mp.pi)) / (2*mp.pi)
lim = {}
for v in (1e-3, 1e-2, 0.05):
    lim[v] = u1mp_hi(v)
out['u0'] = {'analytic': str(exact0), 'u1_1e-3': str(lim[1e-3]), 'u1_1e-2': str(lim[1e-2]), 'u1_0p05': str(lim[0.05]),
                 'diff_0p05': str(lim[0.05]-exact0), 'diff_1e-3': str(lim[1e-3]-exact0)}
log('u_1(0) analytic', exact0)
log('u_1(1e-3) - u_1(0) =', lim[1e-3]-exact0, ' u_1(0.01)-u_1(0) =', lim[1e-2]-exact0, ' u_1(0.05)-u_1(0) =', lim[0.05]-exact0)

# ---- samples
u = C.u1_grid(VMAX, H, a.procs, log)
v = np.arange(len(u)) * H
out['u1_sample'] = dict(n=len(u), u0=float(u[0]), u_0p05=float(u[1]), u_100=float(u[int(round(min(100, VMAX)/H))]))
# ---- masses
def mass_zero_sum(V, nz=1700):
    g = Z[:nz]
    s = np.sum(np.arctan((V-g)/0.5) + np.arctan((V+g)/0.5)) / math.pi
    # tail: smooth density of zeros above the 1700th; N(T) = T/2pi log(T/2pi e) + 7/8, density log(T/2pi)/2pi, lower limit where N = nz + 1/2
    from scipy.optimize import brentq
    N = lambda T: T/(2*math.pi)*math.log(T/(2*math.pi*math.e)) + 7/8
    T0 = brentq(lambda T: N(T) - (nz + 0.5), g[-1]-5, g[-1]+5)
    f = lambda w: math.log(w/(2*math.pi))/(2*math.pi) * (math.atan((V-w)/0.5) + math.atan((V+w)/0.5)) / math.pi
    t = quad(f, T0, 1e7, limit=500)[0]
    # beyond 1e7 the integrand is ~ V log(w)/(2 pi^2 w^2)
    t += V*(math.log(1e7/(2*math.pi))+1)/(2*math.pi**2*1e7)
    return s, t
masses = {}
for V in ([40, 30, 20, 10] if a.quick else [400, 300, 200, 100]):
    k = int(round(V/H))
    simp = float(simpson(u[:k+1], dx=H))
    trap = float(trapezoid(u[:k+1], dx=H))
    s, t = mass_zero_sum(V)
    N = int((Z < V).sum())
    masses[V] = dict(simpson=simp, trapezoid=trap, arctan_sum=s, tail=t, arctan_total=s+t, N=N)
    log('V=%g: Simpson %.6f trapezoid %.6f arctan sum %.6f + tail %.6f = %.6f; N(V)=%d' % (V, simp, trap, s, t, s+t, N))
out['mass'] = masses
cnt = int(((Z > 100) & (Z < 300)).sum()); out['zeros_100_300'] = cnt; out['N100'] = int((Z<100).sum()); out['N300'] = int((Z<300).sum())
log('zeros in (100,300):', cnt, ' N(100)=', out['N100'], ' N(300)=', out['N300'])
# density of zeros quoted at H:261
for V in (150, 350): log('(2pi)^-1 log(v/2pi) at', V, '=', math.log(V/(2*math.pi))/(2*math.pi))

# ---- spectrum
def vangoldt(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
uu = np.concatenate([u[:0:-1], u])
vv = np.concatenate([-v[:0:-1], v])
w = 0.5*(1 + np.cos(np.pi*vv/VMAX))
NFFT = 2**18 if a.quick else 2**22
F = np.abs(np.fft.rfft(uu*w, NFFT))
om = 2*np.pi*np.arange(len(F))/(NFFT*H)
resolution = om[1]
pk, _ = find_peaks(F, prominence=0.02*F[(om > 1.0) & (om < 1.2)].max())
rows = []
ns = [2,3,4,5,7,8,9,11,13,17,19,23,29,31]
ref = None
for n in ns:
    j = pk[np.argmin(np.abs(om[pk] - math.log(n)))]
    rows.append(dict(n=n, peak_omega=float(om[j]), e_omega=float(math.exp(om[j])), log_n=math.log(n), height=float(F[j]), lam_over_n=vangoldt(n)/n))
h3 = [r for r in rows if r['n'] == 3][0]
for r in rows:
    r['rel_height'] = r['height']/h3['height']; r['rel_lam'] = r['lam_over_n']/h3['lam_over_n']
    log('n=%2d peak omega %.4f (log n %.4f) e^omega %.3f  rel height %.3f  Lambda(n)/n rel %.3f' % (r['n'], r['peak_omega'], r['log_n'], r['e_omega'], r['rel_height'], r['rel_lam']))
out['spectrum'] = dict(rows=rows, omega_step=float(resolution), n_peaks_found=int(len(pk)))
log('omega step', resolution)
# ---- fejer / sharp cutoff table
vi = vv; N = len(uu)
zin = Z[(Z > 100) & (Z < 300)]
sel = (vv > 100) & (vv < 300)
def flow(sample, thp, Om, fejer=False):
    X = np.fft.fft(sample*w)
    wm = 2*np.pi*np.fft.fftfreq(N, d=H)
    m = np.where(np.abs(wm) <= Om, np.exp((1-thp)*np.abs(wm)), 0.0)
    if fejer: m = m*np.maximum(1 - np.abs(wm)/Om, 0)
    y = np.fft.ifft(X*m).real
    return y/w
def stats(y):
    ys = y[sel]; vs = vv[sel]
    pk, _ = find_peaks(ys)
    md, nuniq = C.match(vs[pk], zin)
    pk2, _ = find_peaks(ys, height=0.25*ys.max())
    md2, nu2 = C.match(vs[pk2], zin)
    return dict(peaks=int(len(pk)), zeros=int(len(zin)), max_err=md, zeros_with_own_peak=nuniq, min_value=float(ys.min()), max_value=float(ys.max()),
                peaks_height_ge_quarter_max=int(len(pk2)), max_err_those=md2, zeros_with_own_peak_those=nu2)
bad = uu.copy(); bad[len(u)-1] = u[1]
fl = {}
if not a.quick or True:
    for (thp, Om) in [(0.75, 40), (0.5, 40), (0.5, 60)]:
        for name, smp in (('analytic_u0', uu), ('u0_copied_from_0.05', bad)):
            for fe in (False, True):
                if a.quick: continue
                s = stats(flow(smp, thp, Om, fe)); key = '%s|theta=%g|Omega=%g|%s' % (name, thp, Om, 'fejer' if fe else 'sharp')
                fl[key] = s; log(key, s)
out['flow'] = fl
json.dump(out, open(a.out, 'w'), indent=1)
log('done')
