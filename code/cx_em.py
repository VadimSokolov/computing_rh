"""cx_em.py: EM, log domain Nesterov and SQUAREM for the Cauchy deconvolution of the level one Thorin density.

What it computes: the EM tables of Chapter "Complexity and Acceleration" (window table, full range table, rate products,
grid independence, three point exponents, restart variants, thousand zero run), see jobs below.
Claims supported: I:228 to 275, I:785 to 793, I:811 to 814 of the untrusted chapter (checks file misc/reviews/nick-b4/checks-CX3.md).
Conventions: xi(s)=s(s-1)pi^{-s/2}Gamma(s/2)zeta(s)/2; data y_i=u_1(v_i)=pi^{-1} Re (xi'/xi)(1+i v_i), even in v;
A_ij=h P_{1-theta'}(v_i-w_j) with P_a(v)=(a/pi)/(a^2+v^2), unknowns x_j are cell DENSITIES (mass of cell j is h x_j);
KL residual f(x)=h sum_i [y_i log(y_i/(Ax)_i) - y_i + (Ax)_i] (the integral over the data grid; the unscaled sum gave residuals 1/h times the chapter's, so the factor h is what matches the printed 4.1e-3 at 3000 EM steps);
EM x+ = x/c * A^T(y/Ax), c=A^T 1; flat start x_j = sum(y)/sum(A) (total mass matched).
Method: Ax and A^T r by FFT convolution with the Cauchy kernel (Toeplitz), results identical to the dense matrix to rounding.
Nesterov: z=log x, ext=z+beta(z-z_prev), beta=(t-1)/t+, t+=(1+sqrt(1+4t^2))/2, x+=EM(exp ext);
function restart: if f(x+)>f(x) take a plain EM step from x and reset t=1 (one extra EM evaluation, counted).
SQUAREM: three EM evaluations per iteration, alpha=-|r|/|q|, clip at 1e-30, fall back to x2 if f(EM(x'))>=f(x2).
Steps are counted in EM evaluations (one A and one A^T product each).
Error control: u_1 by python-flint balls at 128 bits (derivative by a five point stencil, step 1e-5) rechecked at 192 bits and step 2e-5 on a sample (logged);
grid dependence is itself a check (jobs gridind, expo).
Written from the mathematics (the collaborator's scripts were not read, run or copied).
Usage: python3 cx_em.py JOB [--quick] [--out FILE]; JOB in window, full, gridind, expo, restart, thousand.
"""
import sys, json, time, math, argparse
import numpy as np
from scipy import fft as sfft
from scipy.signal import find_peaks

ZEROS = None
REFINE = True
def zeros(n=1700):
    global ZEROS
    if ZEROS is None:
        import os
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'zeros_1700.txt')
        if not os.path.exists(p): p = 'zeros_1700.txt'
        ZEROS = np.array([float(l.split()[-1]) for l in open(p) if l.strip()])
    return ZEROS[:n]

# ---------- data layer ----------
def _u1_chunk(args):
    vs, prec, dstep = args
    try:
        from flint import acb, ctx
        ctx.prec = prec
        out = []
        for v in vs:
            s = acb(1, v)
            def lz(t):
                return t.zeta()
            d = dstep
            # five point stencil for zeta'
            zp = (-lz(s + 2*d) + 8*lz(s + d) - 8*lz(s - d) + lz(s - 2*d)) / (12*d)
            zz = lz(s)
            xi = 1/s + 1/(s - 1) - acb(math.pi).log()/2 + (s/2).digamma()/2 + zp/zz
            out.append(float(xi.real.mid()) / math.pi)
        return out
    except ImportError:
        import mpmath as mp
        mp.mp.dps = 20
        out = []
        for v in vs:
            s = mp.mpc(1, v)
            xi = 1/s + 1/(s-1) - mp.log(mp.pi)/2 + mp.digamma(s/2)/2 + mp.zeta(s, derivative=1)/mp.zeta(s)
            out.append(float(xi.real)/float(mp.pi))
        return out

def u1_half(vmax, h, procs=1, log=print):
    """u_1 on v=0,h,..,vmax (even function). The grid point v=0 uses v=1e-4 (error O(1e-8))."""
    n = int(round(vmax/h))
    vs = np.arange(n+1)*h
    vs[0] = 1e-4
    t0 = time.time()
    chunks = np.array_split(vs, max(1, procs*8))
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as p:
            res = p.map(_u1_chunk, [(c, 128, 1e-5) for c in chunks])
    else:
        res = [_u1_chunk((c, 128, 1e-5)) for c in chunks]
    u = np.concatenate([np.array(r) for r in res])
    log('u1 on [0,%g] step %g: %d points in %.1f s' % (vmax, h, n+1, time.time()-t0))
    idx = np.unique(np.linspace(0, n, 40).astype(int)); idx = idx[idx > 0]
    u2 = np.array(_u1_chunk((vs[idx], 192, 2e-5)))
    log('u1 precision check (192 bits, step 2e-5) on %d points: max abs diff %.2e (max u %.3f)' % (len(idx), np.abs(u2-u[idx]).max(), u.max()))
    return u

# ---------- operator ----------
class Prob:
    def __init__(self, a0, a1, b0, b1, h, theta, y):
        """data index a0..a1, unknown index b0..b1 (grid v=h*index), y data array."""
        self.h, self.theta = h, theta
        self.m, self.n = a1-a0+1, b1-b0+1
        self.d0 = a0 - b0
        self.y = np.asarray(y, float); assert len(self.y) == self.m
        a = 1 - theta
        self.kf = lambda mm: h*(a/np.pi)/(a*a + (mm*h)**2)
        self.b0, self.a0 = b0, a0
        self.fa = self._prep(self.n, self.m, self.d0)
        self.ft = self._prep(self.m, self.n, -self.d0)
        self.c = self.At(np.ones(self.m))
        self.x0 = np.full(self.n, self.y.sum()/self.A(np.ones(self.n)).sum())
        self.w = (np.arange(self.n) + b0)*h
    def _prep(self, nv, nout, shift):
        u = np.arange(nv+nout-1)
        kk = self.kf(u - (nv-1) + shift)
        L = sfft.next_fast_len(nv + 2*nout + nv)
        return (kk, nv, nout, L, sfft.rfft(kk, L))
    def _conv(self, v, pr):
        kk, nv, nout, L, fk = pr
        c = sfft.irfft(sfft.rfft(v, L)*fk, L)
        return c[nv-1:nv-1+nout]
    def A(self, x): return self._conv(x, self.fa)
    def At(self, r): return self._conv(r, self.ft)
    def f(self, Ax):
        y = self.y
        with np.errstate(divide='ignore', invalid='ignore'):
            t = np.where(y > 0, y*np.log(y/Ax), 0.0)
        return self.h*float((t - y + Ax).sum())
    def em(self, x, Ax):
        return x*self.At(self.y/Ax)/self.c

# ---------- methods; return (f history by evaluations, final x). hist[k]=f after k evaluations (SQUAREM: at iteration ends, NaN between) ----------
def run_em(P, steps):
    x = P.x0.copy(); Ax = P.A(x); hist = np.empty(steps+1); hist[0] = P.f(Ax)
    for k in range(1, steps+1):
        x = P.em(x, Ax); Ax = P.A(x); hist[k] = P.f(Ax)
    return hist, x

def run_nest(P, steps, restart='function', cap=None, fixed=500):
    x = P.x0.copy(); Ax = P.A(x); fx = P.f(Ax)
    hist = np.empty(steps+1); hist[0] = fx
    z = np.log(x); zp = z.copy(); t = 1.0; k = 0; since = 0
    nrest = 0
    while k < steps:
        tn = (1 + math.sqrt(1 + 4*t*t))/2
        beta = (t-1)/tn
        if cap is not None: beta = min(beta, cap)
        ze = z + beta*(z - zp)
        xe = np.exp(ze); Axe = P.A(xe)
        xn = P.em(xe, Axe); Axn = P.A(xn); fn = P.f(Axn)
        k += 1; since += 1
        do_restart = False
        if restart == 'function' and fn > fx: do_restart = True
        elif restart == 'grad':
            zn = np.log(np.maximum(xn, 1e-300))
            if np.dot(zn - ze, zn - z) < 0: do_restart = True   # step against momentum
        elif restart == 'fixed' and since >= fixed: do_restart = True
        if restart == 'function' and do_restart:
            hist[k] = fn  # the rejected evaluation is counted, its value recorded
            if k < steps:
                xn = P.em(x, Ax); Axn = P.A(xn); fn = P.f(Axn); k += 1
                hist[k] = fn
            zp = np.log(x); z = np.log(xn); x, Ax, fx = xn, Axn, fn; t = 1.0; nrest += 1; continue
        hist[k] = fn
        zp = z; z = np.log(np.maximum(xn, 1e-300)); x, Ax, fx = xn, Axn, fn
        if do_restart: t = 1.0; since = 0; nrest += 1
        else: t = tn
    return hist, x, nrest

def run_squarem(P, evals):
    x0 = P.x0.copy(); Ax0 = P.A(x0)
    hist = np.full(evals+1, np.nan); hist[0] = P.f(Ax0); k = 0
    while k + 3 <= evals:
        x1 = P.em(x0, Ax0); Ax1 = P.A(x1)
        x2 = P.em(x1, Ax1); Ax2 = P.A(x2); f2 = P.f(Ax2)
        r = x1 - x0; q = (x2 - x1) - r
        nq = np.linalg.norm(q)
        if nq > 0:
            al = -np.linalg.norm(r)/nq
            xs = np.maximum(x0 - 2*al*r + al*al*q, 1e-30)
            Axs = P.A(xs)
            x3 = P.em(xs, Axs); Ax3 = P.A(x3); f3 = P.f(Ax3)
        else: f3 = np.inf
        if f3 < f2: x0, Ax0, fx = x3, Ax3, f3
        else: x0, Ax0, fx = x2, Ax2, f2
        k += 3
        hist[k-2] = hist[k-1] = np.nan; hist[k] = fx
    return hist, x0

def steps_to(hist, level):
    ok = np.where(hist <= level)[0]
    return int(ok[0]) if len(ok) else None

# ---------- zeros and peaks ----------
def peaks(P, x, lo, hi, prom=0.05):
    """parabolic refined local maxima of the density x with prominence >= prom, located in [lo,hi]."""
    pk, _ = find_peaks(x, prominence=prom)
    out = []
    for j in pk:
        if 0 < j < len(x)-1:
            d = x[j-1] - 2*x[j] + x[j+1]
            off = 0.5*(x[j-1] - x[j+1])/d if d != 0 else 0
            pos = (P.w[j] + (off*P.h if REFINE else 0.0))
            if lo <= pos <= hi: out.append(pos)
    return np.array(out)

def match(pk, zs):
    """distance of each zero to the nearest peak"""
    if len(pk) == 0: return np.full(len(zs), np.inf)
    i = np.searchsorted(pk, zs); i0 = np.clip(i-1, 0, len(pk)-1); i1 = np.clip(i, 0, len(pk)-1)
    return np.minimum(abs(pk[i0]-zs), abs(pk[i1]-zs))

def basin_masses(P, x, zs, left, right):
    """h*sum x between midpoints of consecutive zeros; first basin starts at `left`, last ends at `right`."""
    mid = np.concatenate([[left], (zs[1:]+zs[:-1])/2, [right]])
    cs = np.concatenate([[0.0], np.cumsum(x)])*P.h
    idx = np.clip(np.round(mid/P.h - P.b0 + 0.5).astype(int), 0, P.n)   # cell j covers [w_j-h/2, w_j+h/2)
    return cs[idx[1:]] - cs[idx[:-1]]

def p3(hist, k):
    r = (hist[k]-hist[2*k])/(hist[2*k]-hist[4*k])
    return math.log2(r) if r > 0 else float('nan')   # nan when the objective is not monotone on the three points

# ---------- problems ----------
def full_problem(h, theta, vmax=400.0, wmax=440.0, procs=1, log=print, u=None):
    if u is None: u = u1_half(vmax, h, procs, log)
    ydata = np.concatenate([u[:0:-1], u])
    M = int(round(vmax/h)); N = int(round(wmax/h))
    return Prob(-M, M, -N, N, h, theta, ydata), u

def window_problem(h, theta, procs=1, log=print):
    u = u1_half(300.0, h, procs, log)
    a0, a1 = int(round(100/h)), int(round(300/h))
    return Prob(a0, a1, int(round(80/h)), int(round(320/h)), h, theta, u[a0:a1+1])

def jr(o, path):
    json.dump(o, open(path, 'w'), indent=1, default=lambda v: None if isinstance(v, float) and math.isnan(v) else (v.tolist() if hasattr(v, 'tolist') else str(v)))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('job'); ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=None); ap.add_argument('--procs', type=int, default=1)
    ap.add_argument('--unrefined', action='store_true'); ap.add_argument('--h', type=float, default=None); ap.add_argument('--steps', type=int, default=None)
    a = ap.parse_args()
    out = a.out or 'cx_em_%s%s.json' % (a.job, '_%g' % a.h if a.h else '')
    lf = open(out.replace('.json', '.log'), 'w')
    def log(*s):
        m = ' '.join(str(t) for t in s); print(m, flush=True); lf.write(m+'\n'); lf.flush()
    R = {}; zs_all = zeros()
    q = a.quick
    global REFINE
    REFINE = not a.unrefined
    if a.job == 'window':
        h = 0.2 if q else 0.05; S = 60 if q else 3000
        levels = [1e-2, 4e-3, 1e-3, 3e-4]
        zs = zs_all[(zs_all > 100) & (zs_all < 300)]; log('zeros in (100,300):', len(zs), ' h', h)
        for th in (0.75, 0.5):
            P = window_problem(h, th, a.procs, log)
            for name, fn in (('EM', lambda: run_em(P, S)), ('Nesterov', lambda: run_nest(P, S)), ('SQUAREM', lambda: run_squarem(P, S))):
                t0 = time.time(); r = fn(); hist, x = r[0], r[1]; el = time.time()-t0
                fin = hist[~np.isnan(hist)][-1]
                pk = peaks(P, x, 100, 300); d = match(pk, zs)
                R['%s_%g' % (name, th)] = dict(steps_to=[steps_to(hist, L) for L in levels], f_final=fin, evals=S, time=el,
                    npeaks=len(pk), maxerr=float(d.max()), min_x=float(x.min()), restarts=(r[2] if len(r) > 2 else None))
                log(th, name, R['%s_%g' % (name, th)])
    elif a.job in ('full', 'gridind'):
        hs = [a.h] if a.h else ([0.2] if q else [0.02])
        S5 = 100 if q else 3000
        for h in hs:
            P3, u = full_problem(h, 0.75, procs=a.procs, log=log)
            zs = zs_all[zs_all < 400]; log('zeros below 400:', len(zs), 'h', h)
            runs = {}
            for th, S in ((0.75, 40 if q else 400), (0.5, S5)):
                P = P3 if th == 0.75 else Prob(P3.a0, P3.a0+P3.m-1, P3.b0, P3.b0+P3.n-1, h, th, P3.y)
                for name in ('EM', 'Nest'):
                    t0 = time.time()
                    r = run_em(P, S) if name == 'EM' else run_nest(P, S); el = time.time()-t0
                    hist, x = r[0], r[1]; runs[(th, name)] = hist
                    pk = peaks(P, x, 0, 400); d = match(pk, zs)
                    m = basin_masses(P, x, zs, 0.0, 400.0)
                    R['%g_%s_h%g' % (th, name, h)] = dict(steps=S, f=hist[-1], npeaks=len(pk), maxerr=float(d.max()), meanerr=float(d.mean()),
                        mass_mean=float(m.mean()), mass_sd=float(m.std()), mass_first=float(m[0]), mass_last=float(m[-1]), time=el,
                        mass_sd_interior=float(m[1:-1].std()), nrest=(r[2] if len(r) > 2 else None))
                    log(th, name, R['%g_%s_h%g' % (th, name, h)])
                fstar = min(runs[(th, 'EM')].min(), runs[(th, 'Nest')].min())
                for fs_name, fs in (('min', fstar), ('0.999min', 0.999*fstar)):
                    em, ne = runs[(th, 'EM')], runs[(th, 'Nest')]
                    ks = [50, 100, 200, 400] if th == 0.75 else [500, 1000, 2000, 3000]
                    R['prod_em_%g_h%g_%s' % (th, h, fs_name)] = {k: k*(em[k]-fs) for k in ks if k < len(em)}
                    ks2 = [500, 1000, 2000]
                    R['prod2_nest_%g_h%g_%s' % (th, h, fs_name)] = {k: k*k*(ne[k]-fs) for k in ks2 if k < len(ne)}
                if th == 0.5:
                    em, ne = runs[(th, 'EM')], runs[(th, 'Nest')]
                    st = steps_to(ne, em[-1]); R['ratio_h%g' % h] = dict(em_steps=len(em)-1, nest_steps=st, ratio=(len(em)-1)/st if st else None)
                    log('ratio', R['ratio_h%g' % h])
                    R['fstar_h%g' % h] = fstar
                log('prods', {k: v for k, v in R.items() if k.startswith('prod')})
    elif a.job == 'expo':
        h = a.h or (0.4 if q else 0.08); S = a.steps or (400 if q else 200000)
        P, u = full_problem(h, 0.5, procs=a.procs, log=log)
        t0 = time.time(); hist, x = run_em(P, S); log('EM %d steps h=%g: %.1f s' % (S, h, time.time()-t0))
        ks = [500, 1000, 2000, 4000, 8000, 16000, 32000, 50000]
        R['p3'] = {k: p3(hist, k) for k in ks if 4*k <= S}
        R['f'] = {k: hist[k] for k in [500, 1000, 2000, 4000, 16000, 50000, 100000, 200000] if k <= S}
        log(R['p3'])
        # geometric tail: local slope of log(f_k - f_inf) estimated by f_k - f_2k ratio: K = k / log((f_k-f_2k)/(f_2k-f_3k))
        Kd = {}
        for k in [50000, 100000]:
            if 3*k <= S: Kd[k] = k/math.log((hist[k]-hist[2*k])/(hist[2*k]-hist[3*k]))
        R['K_geometric_estimate'] = Kd; log('K', Kd)
    elif a.job == 'restart':
        h = a.h or (0.4 if q else 0.04); S = a.steps or (400 if q else 20000)
        P, u = full_problem(h, 0.5, procs=a.procs, log=log)
        hem, x = run_em(P, S); R['f_em'] = {k: hem[k] for k in (1000, 2000, 5000, 10000, 20000) if k <= S}
        variants = [('function', dict(restart='function')), ('none', dict(restart='none')), ('gradient', dict(restart='grad')),
                    ('fixed500', dict(restart='fixed', fixed=500)), ('cap0.9', dict(restart='none', cap=0.9))]
        tgt = hem[S]
        for name, kw in variants:
            t0 = time.time(); r = run_nest(P, S, **kw); hn = r[0]
            st = steps_to(hn, tgt)
            d = dict(steps_to_EM_final=st, f_final=hn[-1], nrest=r[2], time=time.time()-t0)
            ks = [500, 1000, 2000, 4000, 5000]
            d['p3'] = {k: p3(hn, k) for k in ks if 4*k <= S}
            if name == 'function':
                d['speedup'] = {k: (k/steps_to(hn, hem[k]) if steps_to(hn, hem[k]) else None) for k in (1000, 2000, 5000, 10000, 20000) if k <= S}
                d['p_running'] = {(k): p3(hn, k) for k in (2000, 3000, 5000) if 4*k <= S}
            R[name] = d; log(name, d)
    elif a.job == 'thousand':
        h = 0.2 if q else 0.025; vmax = 100.0 if q else 1450.0; S = 100 if q else 4000
        wmax = vmax + 40
        tu = time.time(); u = u1_half(vmax, h, a.procs, log); t_data = time.time()-tu
        P, _ = full_problem(h, 0.5, vmax=vmax, wmax=wmax, log=log, u=u)
        t0 = time.time(); hist, x, nr = run_nest(P, S); t_inv = time.time()-t0
        top = 90.0 if q else 1415.0
        zs = zs_all[zs_all < top]
        pk = peaks(P, x, 0, top); d = match(pk, zs)
        m = basin_masses(P, x, zs, 0.0, top)
        gaps = np.diff(zs); order = np.argsort(gaps)[:4]
        R = dict(h=h, ncells=P.n, nzeros=len(zs), npeaks=len(pk), f_final=hist[-1], t_data=t_data, t_inversion=t_inv, nrest=nr,
                 d_mean=float(d.mean()), d_median=float(np.median(d)), d_p99=float(np.percentile(d, 99)), n_over_003=int((d > 0.03).sum()),
                 mass_mean=float(m.mean()), mass_sd=float(m.std()), mass_max=float(m.max()),
                 closest_pairs=[(float(zs[i]), float(zs[i+1]), float(gaps[i])) for i in order],
                 dist_closest=[float(max(d[i], d[i+1])) for i in order], mass_closest=[float(m[i]+m[i+1]) for i in order], zeros_over_003=[float(z) for z in zs[d > 0.03]])
        log(R)
    jr(R, out); log('wrote', out)

if __name__ == '__main__':
    main()
