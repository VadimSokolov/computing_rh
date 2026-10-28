"""lo_gbc.py: Chapter 30 (H = chapter_level_one.tex), Section 5: generative Bayesian computation for the backward flow, redesigned as in reconcile-H section 2.12 (Q2).

What it computes. A quantile network (numpy, two hidden layers of 256 rectified units, 24 outputs: the counts of zeros in eight unit bins centred in a window of length 16 at
tau = 0.1, 0.5, 0.9, pinball loss, Adam, 25 epochs, 24000 training windows) trained on simulations of the corrected model and tested on u_1 itself.
Model. Zeros: GUE, H = (A + A^*)/2, A = X + iY with X, Y independent standard normal M x M, M = 400, whose off diagonal entries have variance 1, so the spectral edge is
R = 2 sqrt(M) (the unit variance of the real and the imaginary part: E|A_ij|^2 = 2, so the off diagonal variance of H is (2 + 2)/4 = 1). The eigenvalues are unfolded by the semicircle law
of radius R, e_k = M F(x_k/R) with F(s) = 1/2 + (s sqrt(1-s^2) + arcsin s)/pi, to mean spacing 1 (checked below by the measured mean spacing in the window), then divided by the density d ~ U[0.45, 0.8].
Forward map on the grid t = -8, -7.9, .., 8 (spacing 0.1): y(t) = sum over the simulated zeros z within 150/d of the window centre of P_{1/2}(t - z) plus, for the zeros beyond, their mean density d:
(d/pi) [(pi/2 - arctan((W - t)/a)) + (pi/2 - arctan((W + t)/a))], W = 150/d, a = 1/2. There is no smooth part g_1 (reconcile-H section 2.12: the zero sum is all of u_1; the old model counted the mean density twice).
Test data: u_1 itself on windows centred at 150, 158, .., 342 (25 windows, 200 bins), true counts from the zeros of data/zeros_1700.txt; the model error of the training forward map on these windows
(the zero sum with the mean density tail, against the flint u_1) is logged. Misfit statistic (our definition): the root mean square over |t| <= 4 of y minus the forward map of the reconstruction that places
round(median count) zeros at the centres of the bins, with mean density sum(counts)/8 outside the bins and beyond 150/d. Injection: the two zeros nearest the window centre are replaced by the quartet at
sigma = 0.75 (widths 0.25 and 0.75 at the midpoint, with mirror images at negative ordinate, total mass unchanged), the control by two line zeros at the midpoint.
Calibration: on 2000 held out simulated windows, the fraction with N <= q_tau for tau = 0.1, 0.5, 0.9 and the coverage of [q_0.1, q_0.9].
Claims supported: H:257 to 271, reconcile-H rows 37, 38, 40, 46. Written from the mathematics (the collaborator's scripts were not read, run or copied).
Error control: the unfolding is checked by the measured density; u_1 data by cx_em (flint, 128 and 192 bits); two training seeds are run and their test numbers logged.
Usage: python3 lo_gbc.py [--quick] [--procs N] [--out FILE]
"""
import sys, os, json, math, time, argparse
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[_k] = '1'   # one thread per process: eigvalsh in 8 processes
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lo_common as C
import cx_em

ap = argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true'); ap.add_argument('--procs', type=int, default=1); ap.add_argument('--out', default='lo_gbc.json')
a = ap.parse_args()
logf = open(a.out.replace('.json', '.log'), 'w')
def log(*s):
    m = ' '.join(str(x) for x in s); print(m, flush=True); logf.write(m + '\n'); logf.flush()

M = 400; A_W = 0.5; TG = np.round(np.arange(-80, 81) * 0.1, 10)
NTRAIN = 200 if a.quick else 24000; NTEST = 100 if a.quick else 2000; EPOCHS = 3 if a.quick else 25
TAUS = np.array([0.1, 0.5, 0.9])

def semicircle_cdf(s):
    s = np.clip(s, -1, 1)
    return 0.5 + (s*np.sqrt(1-s*s) + np.arcsin(s))/np.pi

def tail_map(d, t=TG, W=None):
    W = 150.0/d if W is None else W
    return (d/np.pi)*((np.pi/2 - np.arctan((W - t)/A_W)) + (np.pi/2 - np.arctan((W + t)/A_W)))

def sim_chunk(args):
    seed, n, scale_mode = args
    rs = np.random.RandomState(seed)
    Y = np.empty((n, len(TG)), np.float32); Nn = np.empty((n, 8), np.float32); D = np.empty(n)
    spac = []
    R = 2*math.sqrt(M)
    for i in range(n):
        d = rs.uniform(0.45, 0.8)
        X = rs.randn(M, M) + 1j*rs.randn(M, M)
        Hm = 0.5*(X + X.conj().T)
        x = np.linalg.eigvalsh(Hm)
        e = M*semicircle_cdf(x/R)
        ec = M/2 + rs.uniform(-40, 40)
        z = (e - ec)/d
        sel = np.abs(e - ec) < 150
        zs = z[sel]
        y = (C.P(A_W, TG[:, None] - zs[None, :])).sum(1) + tail_map(d)
        Y[i] = y; D[i] = d
        Nn[i] = [np.sum((zs >= -4 + k) & (zs < -3 + k)) for k in range(8)]
        if i < 20: spac.append(np.sum(np.abs(zs) < 20)/40.0)
    return Y, Nn, D, spac

def simulate(n, seed0, procs):
    k = max(1, procs*4)
    sizes = [n//k + (1 if i < n % k else 0) for i in range(k)]
    jobs = [(seed0 + i, s, 0) for i, s in enumerate(sizes)]
    if procs > 1:
        from multiprocessing import Pool
        with Pool(procs) as p: r = p.map(sim_chunk, jobs)
    else:
        r = [sim_chunk(j) for j in jobs]
    return (np.concatenate([x[0] for x in r]), np.concatenate([x[1] for x in r]), np.concatenate([x[2] for x in r]), sum([x[3] for x in r], []))

# ---------------- network ----------------
class Net:
    def __init__(self, nin, nh, nout, seed):
        rs = np.random.RandomState(seed)
        self.p = [rs.randn(nin, nh).astype(np.float32)*math.sqrt(2/nin), np.zeros(nh, np.float32),
                  rs.randn(nh, nh).astype(np.float32)*math.sqrt(2/nh), np.zeros(nh, np.float32),
                  rs.randn(nh, nout).astype(np.float32)*math.sqrt(1/nh), np.zeros(nout, np.float32)]
        self.m = [np.zeros_like(q) for q in self.p]; self.v = [np.zeros_like(q) for q in self.p]; self.t = 0
    def fwd(self, x):
        W1, b1, W2, b2, W3, b3 = self.p
        h1 = np.maximum(x@W1 + b1, 0); h2 = np.maximum(h1@W2 + b2, 0)
        return h2@W3 + b3, (x, h1, h2)
    def step(self, x, tgt, taus_out, lr):
        q, (x, h1, h2) = self.fwd(x)
        r = tgt - q
        loss = np.mean(np.maximum(taus_out*r, (taus_out - 1)*r))
        g = np.where(r > 0, -taus_out, 1 - taus_out)/r.size
        W1, b1, W2, b2, W3, b3 = self.p
        gW3 = h2.T@g; gb3 = g.sum(0)
        gh2 = (g@W3.T)*(h2 > 0); gW2 = h1.T@gh2; gb2 = gh2.sum(0)
        gh1 = (gh2@W2.T)*(h1 > 0); gW1 = x.T@gh1; gb1 = gh1.sum(0)
        grads = [gW1, gb1, gW2, gb2, gW3, gb3]
        self.t += 1
        for i, gr in enumerate(grads):
            self.m[i] = 0.9*self.m[i] + 0.1*gr; self.v[i] = 0.999*self.v[i] + 0.001*gr*gr
            mh = self.m[i]/(1 - 0.9**self.t); vh = self.v[i]/(1 - 0.999**self.t)
            self.p[i] = self.p[i] - (lr*mh/(np.sqrt(vh) + 1e-8)).astype(np.float32)
        return loss

def train(Xtr, Ntr, seed):
    net = Net(Xtr.shape[1], 256, 24, seed)
    tgt = np.repeat(Ntr[:, None, :], 3, axis=1).reshape(len(Ntr), 24)  # layout: tau-major below
    taus_out = np.repeat(TAUS, 8)[None, :]
    tgt = np.concatenate([Ntr, Ntr, Ntr], 1)
    rs = np.random.RandomState(seed + 100); B = 128
    for ep in range(EPOCHS):
        perm = rs.permutation(len(Xtr)); tot = 0; nb = 0
        for s in range(0, len(Xtr) - B + 1, B):
            idx = perm[s:s+B]
            tot += net.step(Xtr[idx], tgt[idx], taus_out, 1e-3 if ep < 15 else 3e-4); nb += 1
        if ep % 5 == 4 or ep == EPOCHS - 1: log('epoch %d train pinball loss %.4f' % (ep+1, tot/nb))
    return net

def predict(net, X):
    q, _ = net.fwd(X)
    return np.sort(q.reshape(len(X), 3, 8), axis=1)   # rearranged so q10 <= q50 <= q90

def pinball(q, N):
    r = N[:, None, :] - q
    t = TAUS[None, :, None]
    return float(np.mean(np.maximum(t*r, (t-1)*r)))

# ---------------- real data ----------------
Z = C.zeros()
centres = np.arange(150, 343, 8) if not a.quick else np.arange(150, 200, 8)
vmax = int(centres.max()) + 9
u = cx_em.u1_half(vmax, 0.1, a.procs, log)
def rho(c): return math.log(c/(2*math.pi))/(2*math.pi)
def real_window(c):
    i0 = int(round((c-8)/0.1)); return u[i0:i0+161]
def true_counts(c):
    return np.array([np.sum((Z >= c-4+k) & (Z < c-3+k)) for k in range(8)], float)
def zero_model(c, zs_over=None):
    d = rho(c); W = 150.0/d
    g = Z[np.abs(Z - c) < W]
    return C.P(A_W, TG[:, None] - (g - c)[None, :]).sum(1) + tail_map(d, W=W)
def forward_recon(counts, d_est=None):
    cen = -3.5 + np.arange(8)
    d = counts.sum()/8.0 if d_est is None else d_est
    y = sum(n*C.P(A_W, TG - c0) for n, c0 in zip(counts, cen)) if counts.sum() > 0 else np.zeros_like(TG)
    # outside the bins: mean density d on |t| > 4 up to W, as a continuous density, plus the tail beyond W
    W = 150.0/max(d, 1e-3)
    outside = (d/np.pi)*((np.arctan((-4 - TG)/A_W) - np.arctan((-W - TG)/A_W)) + (np.arctan((W - TG)/A_W) - np.arctan((4 - TG)/A_W)))
    return y + outside + tail_map(max(d, 1e-3), W=W)
def misfit(y, counts):
    r = y - forward_recon(counts)
    c = np.abs(TG) <= 4
    return float(np.sqrt(np.mean(r[c]**2)))

res = {}
t0 = time.time()
Xs, Ns, Ds, spac = simulate(NTRAIN, 1000, a.procs)
log('simulated %d training windows in %.0f s; mean spacing measure (zeros per unit position times 1/d, should be 1/d... ): see density' % (NTRAIN, time.time()-t0))
# density check: mean number of zeros per unit length divided by d, through the bin counts
dens_meas = Ns.sum(1)/8.0
res['density'] = dict(d_min=float(Ds.min()), d_max=float(Ds.max()), mean_ratio_measured_over_d=float(np.mean(dens_meas/Ds)),
                      mean_zeros_per_unit_d_min_half=float(np.mean(dens_meas[Ds < 0.55]/Ds[Ds < 0.55])))
log('density check: mean (zeros per unit length in the 8 bins)/d = %.4f' % res['density']['mean_ratio_measured_over_d'])
res['density']['wrong_edge_density_range'] = [0.45/math.sqrt(2), 0.8/math.sqrt(2)]
mu, sd = float(Xs.mean()), float(Xs.std())
Xtest_sim, Ntest_sim, Dtest_sim, _ = simulate(NTEST, 90000, a.procs)
log('simulation done at %.0f s' % (time.time()-t0))
# real test data
realX = np.array([real_window(c) for c in centres]); realN = np.array([true_counts(c) for c in centres])
mod = np.array([zero_model(c) for c in centres])
res['model_error_real'] = dict(max_abs=float(np.abs(realX - mod).max()), mean_abs=float(np.abs(realX - mod).mean()), max_u=float(realX.max()))
log('training forward map (zero sum + mean density tail) against u_1 on the 25 windows: max abs diff %.2e, mean %.2e (max u %.3f)' % (res['model_error_real']['max_abs'], res['model_error_real']['mean_abs'], realX.max()))
# injected and control windows
injX = []; ctlX = []
for c, y in zip(centres, realX):
    ia = int(np.argmin(np.abs(Z - c)))
    j = ia if abs(0.5*(Z[ia]+Z[ia+1]) - c) < abs(0.5*(Z[ia-1]+Z[ia]) - c) else ia-1
    ga, gb = Z[j] - c, Z[j+1] - c; gm = 0.5*(ga+gb)
    removed = sum(C.P(A_W, TG - s*g) for g in (ga, gb) for s in (1,))  # mirror zeros are far (distance ~ 2c), negligible and identical in all three
    quartet = C.P(0.25, TG - gm) + C.P(0.75, TG - gm)
    ctl = 2*C.P(A_W, TG - gm)
    injX.append(y - removed + quartet); ctlX.append(y - removed + ctl)
injX = np.array(injX); ctlX = np.array(ctlX)

for seed in (1, 2):
    net = train((Xs - mu)/sd, Ns, seed)
    qs = predict(net, (Xtest_sim - mu)/sd)
    cal = {str(t): float(np.mean(Ntest_sim <= qs[:, i, :])) for i, t in enumerate(TAUS)}
    cov = float(np.mean((Ntest_sim >= qs[:, 0, :]) & (Ntest_sim <= qs[:, 2, :])))
    med = np.round(qs[:, 1, :])
    sim = dict(pinball=pinball(qs, Ntest_sim), cal=cal, coverage_10_90=cov, exact=float(np.mean(med == Ntest_sim)), within1=float(np.mean(np.abs(med - Ntest_sim) <= 1)))
    qr = predict(net, (realX - mu)/sd)
    medr = np.round(qr[:, 1, :])
    err = medr - realN
    two = realN >= 2
    real = dict(pinball=pinball(qr, realN), exact=float(np.mean(err == 0)), within1=float(np.mean(np.abs(err) <= 1)),
                n_bins=int(err.size), n_exact=int(np.sum(err == 0)), n_within1=int(np.sum(np.abs(err) <= 1)),
                undercounts=int(np.sum(err < 0)), overcounts=int(np.sum(err > 0)), errors_in_bins_with_ge2=int(np.sum((err != 0) & two)),
                bins_with_ge2=int(two.sum()),
                cal={str(t): float(np.mean(realN <= qr[:, i, :])) for i, t in enumerate(TAUS)},
                coverage_10_90=float(np.mean((realN >= qr[:, 0, :]) & (realN <= qr[:, 2, :]))))
    def mf(X):
        q = predict(net, (X - mu)/sd); mm = np.round(q[:, 1, :])
        return np.array([misfit(y, m) for y, m in zip(X, mm)])
    m_real, m_ctl, m_inj = mf(realX), mf(ctlX), mf(injX)
    real['misfit'] = dict(real_mean=float(m_real.mean()), real_sd=float(m_real.std()), control_mean=float(m_ctl.mean()), injected_mean=float(m_inj.mean()), injected_sd=float(m_inj.std()),
                          n_injected_gt_real=int(np.sum(m_inj > m_real)), n_injected_gt_control=int(np.sum(m_inj > m_ctl)), n=int(len(m_real)))
    # same statistic with the TRUE counts (what a perfect network would give)
    mt = np.array([misfit(y, n) for y, n in zip(realX, realN)])
    mti = np.array([misfit(y, n) for y, n in zip(injX, realN)])
    real['misfit_true_counts'] = dict(real_mean=float(mt.mean()), injected_mean=float(mti.mean()))
    res['seed%d' % seed] = dict(simulated=sim, real=real)
    log('seed', seed, 'simulated:', sim)
    log('seed', seed, 'real:', real)
json.dump(res, open(a.out, 'w'), indent=1)
log('done in %.0f s' % (time.time()-t0))
