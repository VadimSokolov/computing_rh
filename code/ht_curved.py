# ht_curved.py: the numerical claims of Section "Curved boundaries: the Roberts and Shortland view" of Nick Polson's hitting times chapter
# (chapter lines H:404 to H:486) for the inverse first passage boundary b of the reciprocal Xi law mu, E exp(-s tau) = xi(1/2)/xi(1/2 + sqrt s).
# Claims supported: H:404 (b(t)/log(1/t) at t = 1e-4 to 1e-16), H:411 (Talbot against residues), H:413 (linear boundary recovered to 1e-6),
# H:414 (step halving), H:415 (enclosure of Proposition smallb, widths 0.004 and 0.015), H:417 (mass outside [0.0045, 0.12], mean 0.0231),
# H:420 to H:425 (slope -19.99, 3/(2c) = 0.0750, formula values, the exponential law and mu to t = 2 by the complementary form),
# Table tab:ifpt (H:438 to H:443), Table tab:resp and K (H:463), |K| <= 215 (H:486).
# Conventions. xi(s) = s(s-1)/2 pi^(-s/2) Gamma(s/2) zeta(s); LT(s) = xi(1/2)/xi(1/2 + sqrt s) = prod (1 + s/gamma_n^2)^(-1) (every zero rho = 1/2 +- i gamma_n
# gives the factor 1 + w^2/gamma_n^2 in w = sqrt s). F(t) = mu(0, t], S = 1 - F. The residue of LT at s = -gamma_k^2 is
# R_k = 2 i gamma_k xi(1/2)/xi'(rho_k) (ds = 2 w dw), so the density is sum R_k exp(-gamma_k^2 t) and S(t) = sum (R_k/gamma_k^2) exp(-gamma_k^2 t).
# Toy law (H:430): LT_toy = LT * (1 + s/g^2)/((1 + s/lam)(1 + s/conj lam)), -lam = (delta + i g)^2; control: lam replaced by the real pair (g +- 0.001)^2.
# Methods. F(t) for small t: the Bromwich integral on the vertical line through the saddle point sigma (minimum of sigma t + log LT(sigma) - log sigma),
# trapezoid rule of step dy = 2 pi/P in y, whose aliasing error is below exp(-sigma P) F(t + P) (P = 0.04); run again with P = 0.06 and at a second
# precision, the largest relative difference is logged. F(t) for larger t: the residue series with 120 zeros (zeros of data/zeros_1700.txt polished by
# a secant iteration to the working precision of the series, which is chosen from log10 F at the split so that the cancellation e^(95) leaves 40 digits).
# Boundary: the first passage equation Phibar(b(t)/sqrt t) = int_0^t Phibar((b(t) - b(s))/sqrt(t - s)) mu(ds), solved forward in the log scale (scipy log_ndtr) with
# the midpoint rule on cells (t_(j-1), t_j] of mass F_j - F_(j-1) (exact differences in multiprecision), b(s) at the cell midpoint = mean of the two nodal values, the
# current cell implicit, the root bracketed next to the extrapolated previous value and refined by brentq. The mass below the first node t0 is one atom at t0
# with b(t0) = sqrt(t0) Phibar^(-1)(F(t0)) (the lower bound of Proposition smallb); its effect is tested by changing t0 and by adding 0.01 to b(t0).
# Once b < 0 the equivalent complementary form Phi(b/sqrt t) = S(t) + int Phi((b - b(s))/sqrt(t - s)) mu(ds) (all terms positive) is used.
# Error control: step halving (h and h/2), change of t0, of b(t0), two P, two precisions, all logged.
# Options: --task {cdf,smallt,formula,lin,solve,resp,tail,all}  --quick (seconds)  --out BASE (writes BASE.json and BASE.log; default data/ht_curved).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import argparse, json, math, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.special import log_ndtr
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
T0 = time.time()
LOG = None
RES = {}
ZFILE = None
ROWS = None
DOK = True


def log(*a):
    line = '[%7.1fs] %s' % (time.time() - T0, ' '.join(str(x) for x in a))
    print(line, flush=True)
    if LOG:
        LOG.write(line + '\n'); LOG.flush()


def xi(s):
    s = mp.mpmathify(s)
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def logxi(s):
    s = mp.mpmathify(s)
    return mp.log(s * (s - 1) / 2) - s / 2 * mp.log(mp.pi) + mp.loggamma(s / 2) + mp.log(mp.zeta(s))


def LTmu(s):
    return mp.exp(logxi(mp.mpf(1) / 2) - logxi(mp.mpf(1) / 2 + mp.sqrt(s)))


def zeros_file(n):
    return [mp.mpf(l.strip()) for l in open(ZFILE).read().split()[:n]]


def polish(z, dps):
    f = lambda t: mp.re(xi(mp.mpc(mp.mpf(1) / 2, t)))
    return mp.findroot(f, (z, z + mp.mpf(10) ** (-30)), solver='secant', tol=mp.mpf(10) ** (-dps + 8), maxsteps=60)


def poles_base(nz, dps):
    """zeros polished to dps and the residues R_k, as lists of mpf."""
    with mp.workdps(dps):
        zs = zeros_file(nz)
        if dps > 60:
            zs = [polish(z, dps) for z in zs]
        x0 = xi(mp.mpf(1) / 2)
        out = []
        for g in zs:
            rho = mp.mpc(mp.mpf(1) / 2, g)
            dxi = rho * (rho - 1) / 2 * mp.pi ** (-rho / 2) * mp.gamma(rho / 2) * mp.zeta(rho, derivative=1)
            R = mp.re(x0 * 2j * g / dxi)
            out.append((g * g, R))
        return zs, out


class Law:
    def __init__(self, poles, M=None):
        self.poles = poles  # list of (p, c): density sum c exp(-p t), p mpf or mpc
        self.M = M

    def lt(self, s):
        v = LTmu(s)
        return v * self.M(s) if self.M else v

    def S(self, t):
        return mp.re(mp.fsum(c / p * mp.exp(-p * t) for p, c in self.poles))

    def f(self, t):
        return mp.re(mp.fsum(c * mp.exp(-p * t) for p, c in self.poles))


def toy_law(base, zs, m, delta, kind):
    g = zs[m - 1]; g2 = g * g
    if kind == 'toy':
        lam = -(mp.mpf(delta) + 1j * g) ** 2
        a = [lam, mp.conj(lam)]
    else:
        a = [(g + mp.mpf('0.001')) ** 2, (g - mp.mpf('0.001')) ** 2]
    M = lambda s: (1 + s / g2) / ((1 + s / a[0]) * (1 + s / a[1]))
    poles = [(p, c * M(-p)) for k, (p, c) in enumerate(base) if k != m - 1]
    for i in (0, 1):
        ci = a[i] * LTmu(-a[i]) * (1 - a[i] / g2) / (1 - a[i] / a[1 - i])
        poles.append((a[i], ci))
    return Law(poles, M)


def find_sigma(lt, t, dens, hint=None):
    def g(x):
        sg = mp.exp(mp.mpf(x))
        v = sg * t + mp.re(mp.log(lt(sg))) - (0 if dens else mp.log(sg))
        return float(v)
    lo, hi = (math.log(0.5), math.log(1e12)) if hint is None else (math.log(hint) - 0.7, math.log(hint) + 0.7)
    gr = (math.sqrt(5) - 1) / 2
    c = hi - gr * (hi - lo); d = lo + gr * (hi - lo)
    fc, fd = g(c), g(d)
    while hi - lo > 0.01:
        if fc < fd:
            hi, d, fd = d, c, fc; c = hi - gr * (hi - lo); fc = g(c)
        else:
            lo, c, fc = c, d, fd; d = lo + gr * (hi - lo); fd = g(d)
    xm = (lo + hi) / 2
    return math.exp(xm), g(xm)


def inv_vertical(lt, t, P=0.04, dens=False, dps=40, hint=None):
    with mp.workdps(dps):
        t = mp.mpf(t)
        sg0, gmin = find_sigma(lt, t, dens, hint)
        sig = mp.mpf(sg0)
        P = max(P, (45.0 - gmin) / sg0, 1.5 * float(t))
        dy = 2 * mp.pi / P
        tot = mp.mpf(0); peak = mp.mpf(0); j = 0; small = 0
        while True:
            y = j * dy
            s = mp.mpc(sig, y)
            v = mp.re(mp.exp(1j * y * t) * lt(s) / (1 if dens else s))
            tot += v / 2 if j == 0 else v
            peak = max(peak, abs(v))
            if j > 3 and abs(v) < peak * mp.mpf(10) ** (-(dps - 8)):
                small += 1
                if small >= 3:
                    break
            j += 1
            if j > 6000:
                raise RuntimeError('no decay')
        return mp.exp(sig * t) / mp.pi * dy * tot, float(sig)


def Fvals(law, tg, split, dps_res, P=0.04, shift=0.0):
    """mp arrays F, S at nodes tg + shift (Talbot below split, residues above); returns (F, S, info)."""
    Fv, Sv = [], []
    hint = None
    for t in tg:
        tt = t + shift
        if tt <= split:
            v, hint = inv_vertical(law.lt, tt, P=P, hint=hint)
            with mp.workdps(dps_res):
                Fv.append(mp.mpf(v)); Sv.append(1 - mp.mpf(v))
        else:
            with mp.workdps(dps_res):
                s = law.S(mp.mpf(tt))
                Sv.append(s); Fv.append(1 - s)
    return Fv, Sv


def weights(Fv, Sv, dps_res):
    """log of the cell masses, from exact multiprecision differences."""
    lw = [None]
    with mp.workdps(dps_res):
        for j in range(1, len(Fv)):
            d1 = Fv[j] - Fv[j - 1]
            d2 = Sv[j - 1] - Sv[j]
            d = d2 if abs(d2) > 0 and Sv[j] < mp.mpf(1) / 4 else d1
            lw.append(float(mp.log(d)))
    return np.array([0.0] + lw[1:])


def phibar_inv_log(lp):
    """x with log Phibar(x) = lp (lp < 0)."""
    return brentq(lambda x: log_ndtr(-x) - lp, -40.0, 1e10, xtol=1e-14, rtol=1e-14)


def solve_boundary(tg, Fv, Sv, dps_res, b_off=0.0, switch=True):
    """tg nodes (floats), returns b at nodes. Atom at tg[0] of mass F(tg[0])."""
    n = len(tg)
    with mp.workdps(dps_res):
        lF0 = float(mp.log(Fv[0]))
        lS = np.array([float(mp.log(s)) for s in Sv])
    lw = weights(Fv, Sv, dps_res)
    tg = np.asarray(tg, dtype=float)
    b = np.zeros(n)
    b[0] = math.sqrt(tg[0]) * phibar_inv_log(lF0) + b_off
    # cells j = 1..n-1 : (tg[j-1], tg[j]], midpoint s_j, b_s = (b[j-1] + b[j])/2 ; atom index 0 at s = tg[0], b[0]
    S_mid = np.zeros(n); S_mid[0] = tg[0]
    S_mid[1:] = (tg[:-1] + tg[1:]) / 2
    LW = np.zeros(n); LW[0] = lF0; LW[1:] = lw[1:]
    BS = np.zeros(n); BS[0] = b[0]
    for k in range(1, n):
        t = tg[k]
        sq = math.sqrt(t)
        mid = S_mid[:k + 1].copy(); lwk = LW[:k + 1]
        dt = np.sqrt(t - mid)
        bsk = BS[:k + 1].copy()
        comp = switch and (b[k - 1] < 0.0)

        def r(bk):
            bsk[k] = (b[k - 1] + bk) / 2
            if comp:
                rhs = np.logaddexp(lS[k], np.logaddexp.reduce(lwk + log_ndtr((bk - bsk) / dt)))
                return rhs - log_ndtr(bk / sq)
            rhs = np.logaddexp.reduce(lwk + log_ndtr(-(bk - bsk) / dt))
            return rhs - log_ndtr(-bk / sq)
        if k >= 2:
            step = (tg[k] - tg[k - 1]) / (tg[k - 1] - tg[k - 2]) * (b[k - 1] - b[k - 2])
            g = b[k - 1] + step
            d0 = max(3e-4, 0.3 * abs(step))
        else:
            g = b[0]; d0 = 1e-3
        rg = r(g)
        found = None
        d = d0
        for _ in range(16):
            for sgn in (-1, 1):
                x = g + sgn * d
                rx = r(x)
                if rx * rg <= 0:
                    found = (g, x); break
            if found:
                break
            d *= 2
        if not found:
            raise RuntimeError('no bracket at t=%g (r(g)=%g)' % (t, rg))
        bk = brentq(r, min(found), max(found), xtol=1e-14, rtol=1e-14)
        b[k] = bk
        BS[k] = (b[k - 1] + bk) / 2
    return b


def lin_test(h, t0, tmax):
    """linear boundary b = 1 - t/2: crossing law F(t) = Phibar((1 - t/2)/sqrt t) + e Phibar((1 + t/2)/sqrt t)."""
    dps = 80
    tg = [t0 + k * h for k in range(int(round((tmax - t0) / h)) + 1)]
    with mp.workdps(dps):
        def Fl(t):
            t = mp.mpf(t)
            Pb = lambda x: mp.erfc(x / mp.sqrt(2)) / 2
            return Pb((1 - t / 2) / mp.sqrt(t)) + mp.e * Pb((1 + t / 2) / mp.sqrt(t))
        Fv = [Fl(t) for t in tg]
        Sv = [1 - f for f in Fv]
    b = solve_boundary(tg, Fv, Sv, dps, switch=False)
    ex = np.array([1 - t / 2 for t in tg])
    return np.array(tg), b, ex


def legendre_b(t, dps=30):
    """L(t) = sup_s (psi(s) - s t), psi(s) = log xi(1/2 + sqrt s) - log xi(1/2), s > 0 real."""
    with mp.workdps(dps):
        t = mp.mpf(t)
        x0 = mp.re(logxi(mp.mpf(1) / 2))
        def g(x):
            s = mp.exp(mp.mpf(x))
            return float(mp.re(logxi(mp.mpf(1) / 2 + mp.sqrt(s)) - x0) - s * t)
        lo, hi = 0.0, 4.0 * math.log(10) * 40
        gr = (math.sqrt(5) - 1) / 2
        c = hi - gr * (hi - lo); d = lo + gr * (hi - lo); fc, fd = -g(c), -g(d)
        while hi - lo > 1e-7:
            if fc < fd:
                hi, d, fd = d, c, fc; c = hi - gr * (hi - lo); fc = -g(c)
            else:
                lo, c, fc = c, d, fd; d = lo + gr * (hi - lo); fd = -g(d)
        xm = (lo + hi) / 2
        s = mp.exp(mp.mpf(xm))
        # polish by Newton on psi'(s) = t with a numerical derivative
        ps = lambda ss: mp.re(logxi(mp.mpf(1) / 2 + mp.sqrt(ss)))
        sm = s
        L = ps(sm) - x0 - sm * t
        return L, sm


def task_smallt(R, quick):
    out = {}
    for e in (4, 8, 12, 16):
        t = mp.mpf(10) ** (-e)
        L, s = legendre_b(t, dps=30)
        L2, s2 = legendre_b(t, dps=45)
        b1 = mp.sqrt(2 * t * L)
        lp = -float(L)
        b2 = math.sqrt(float(t)) * phibar_inv_log(lp)
        ratio1 = float(b1) / (e * math.log(10)); ratio2 = b2 / (e * math.log(10))
        out['1e-%d' % e] = dict(L=float(L), s=float(s), relL=float(abs(L - L2) / L), ratio_sqrt2tL=ratio1, ratio_phibarinv=ratio2)
        log('smallt t=1e-%d: L=%.6g saddle s=%.4g  b/log(1/t): sqrt(2tL) -> %.4f, sqrt(t) Phibar^-1(e^-L) -> %.4f; relative change of L between 30 and 45 digits %.2e' %
            (e, float(L), float(s), ratio1, ratio2, float(abs(L - L2) / L)))
    out['limit_1_over_2sqrt2'] = 1 / (2 * math.sqrt(2))
    R['smallt'] = out


def task_formula(R, quick):
    g1 = float(mp.mpf(open(ZFILE).readline().strip()))
    c = math.sqrt(2) * g1
    out = dict(gamma1=g1, c=c, three_over_2c=3 / (2 * c))
    log('c = sqrt2 gamma_1 = %.6f, 3/(2c) = %.5f' % (c, 3 / (2 * c)))
    fm = {}
    for t in (0.8, 1.2, 1.6, 2.0):
        fm[t] = 3 / (2 * c) - 1 / (c * math.log(g1 * g1 * t))
    out['mu_formula_t(b%+c)'] = fm
    log('mu formula 3/(2c) - 1/(c log(gamma1^2 t)) at t=0.8,1.2,1.6,2: ' + ', '.join('%.4f' % v for v in fm.values()))
    fe = {t: 3 / (2 * math.sqrt(2)) - 1 / (math.sqrt(2) * math.log(t)) for t in (20, 50, 100, 150)}
    out['exp_formula'] = fe
    log('exponential law formula at t=20,50,100,150: ' + ', '.join('%.4f' % v for v in fe.values()))
    R['formula'] = out


def get_setup(quick):
    if quick:
        return dict(nz=60, split=0.012, t0=0.005, h=5e-4, tmax=0.06, dps_pole=80)
    return dict(nz=120, split=0.0064, t0=0.003, h=4e-5, tmax=0.12, dps_pole=None)


def build(quick, base_dps=None):
    """base poles at the working precision chosen from log10 F(split)."""
    S = get_setup(quick)
    with mp.workdps(40):
        law0 = Law([])
        v, _ = inv_vertical(LTmu, S['split'], P=0.04)
        lf = float(mp.log(v))
    dps = S['dps_pole'] or int(0.4343 * (-lf)) + 110
    log('F(%g) = exp(%.3f); residue precision %d digits' % (S['split'], lf, dps))
    zs, base = poles_base(S['nz'], dps)
    return S, dps, zs, base, lf


def task_cdf(R, quick, ctx):
    S, dps, zs, base, lf = ctx
    law = Law(base)
    out = dict(logF_split=lf, dps_res=dps)
    # Talbot against residues on the overlap, and the two P, two precisions
    rows = []
    split = S['split']
    for t in ([0.013, 0.015, 0.02] if quick else [0.0065, 0.007, 0.008, 0.0085, 0.01]):
        a, _ = inv_vertical(law.lt, t, P=0.04, dps=40)
        a2, sg = inv_vertical(law.lt, t, P=0.06, dps=55)
        with mp.workdps(dps):
            r = 1 - law.S(mp.mpf(t))
            rel = abs(a - r) / r; rel2 = abs(a - a2) / a
        rows.append(dict(t=t, F=mp.nstr(r, 20), rel_vertical_vs_residue=float(rel), rel_P_dps=float(rel2), sigma=sg, sigmaP=sg * 0.04 + float(mp.log(r))))
        log('overlap t=%g: F=%s, vertical line vs residues relative %.2e, vertical line P=0.04/40 digits vs P=0.06/55 digits %.2e (sigma P + log F = %.1f)' % (t, mp.nstr(r, 15), float(rel), float(rel2), rows[-1]['sigmaP']))
    out['overlap'] = rows
    # residue truncation: 120 vs 100 zeros at the split
    with mp.workdps(dps):
        l2 = Law(base[:-20])
        t = mp.mpf(split)
        out['trunc_20_zeros_rel'] = float(abs(law.S(t) - l2.S(t)) / (1 - law.S(t)))
    log('residue truncation: dropping the last 20 zeros changes F(split) by relative %.2e' % out['trunc_20_zeros_rel'])
    # mean
    with mp.workdps(40):
        g = lambda w: xi(mp.mpf(1) / 2 + w)
        mean = mp.diff(g, 0, 2) / (2 * xi(mp.mpf(1) / 2))
        zs_all = [mp.mpf(x) for x in open(ZFILE).read().split()]
        T = zs_all[-1]
        part = mp.fsum(1 / (z * z) for z in zs_all) + (mp.log(T / (2 * mp.pi)) + 1) / (2 * mp.pi * T)
    out['mean_xi_second_derivative'] = float(mean); out['mean_zero_sum'] = float(part)
    log('mean of mu: xi\'\'(1/2)/(2 xi(1/2)) = %.10f ; sum 1/gamma^2 over %d zeros + tail integral = %.10f' % (float(mean), len(zs_all), float(part)))
    # mass outside [0.0045, 0.12]
    ts = [0.0045, 0.12]
    Fs = []
    with mp.workdps(dps):
        for t in ts:
            if t <= split:
                v, _ = inv_vertical(law.lt, t, P=0.04, dps=40); Fs.append(mp.mpf(v))
            else:
                Fs.append(1 - law.S(mp.mpf(t)))
    out['F(0.0045)'] = mp.nstr(Fs[0], 6); out['S(0.12)'] = mp.nstr(1 - Fs[1], 6)
    log('mass outside [0.0045, 0.12]: F(0.0045) = %s, S(0.12) = %s, sum = %s' % (mp.nstr(Fs[0], 5), mp.nstr(1 - Fs[1], 5), mp.nstr(Fs[0] + 1 - Fs[1], 5)))
    R['cdf'] = out


def run_mu(ctx, h, t0, tmax, b_off=0.0, shift=0.0, law=None, P=0.04):
    S, dps, zs, base, lf = ctx
    law = law or Law(base)
    n = int(round((tmax - t0) / h))
    tg = [t0 + k * h for k in range(n + 1)]
    Fv, Sv = Fvals(law, tg, S['split'], dps, P=P, shift=shift)
    b = solve_boundary(tg, Fv, Sv, dps, b_off=b_off)
    return np.array(tg), b, (Fv, Sv)


def cmp_max(t1, b1, t2, b2, lo, hi):
    """max |b1 - b2| over common nodes in [lo, hi]."""
    d = {round(t, 9): v for t, v in zip(t2, b2)}
    m = 0.0
    for t, v in zip(t1, b1):
        if lo - 1e-12 <= t <= hi + 1e-12 and round(t, 9) in d:
            m = max(m, abs(v - d[round(t, 9)]))
    return m


def task_lin(R, quick):
    out = {}
    cfg = [(2e-3, 0.02, 1.0)] if quick else [(5e-4, 0.01, 1.0), (2.5e-4, 0.01, 1.0)]
    for h, t0, tm in cfg:
        t, b, ex = lin_test(h, t0, tm)
        e = float(np.max(np.abs(b - ex)))
        out['h=%g' % h] = dict(max_err=e, max_err_after_0p05=float(np.max(np.abs(b - ex)[t >= 0.05])))
        log('linear boundary b = 1 - t/2, t0=%g, h=%g on [%g, %g]: max |b - (1 - t/2)| = %.3e (t >= 0.05: %.3e)' % (t0, h, t0, tm, e, out['h=%g' % h]['max_err_after_0p05']))
    R['lin'] = out


def task_solve(R, quick, ctx):
    S, dps, zs, base, lf = ctx
    h, t0, tmax = S['h'], S['t0'], S['tmax']
    law = Law(base)
    t1, b1, (Fv, Sv) = run_mu(ctx, h, t0, tmax)
    t2, b2, _ = run_mu(ctx, h / 2, t0, tmax)
    lo = 0.0045 if not quick else 0.007
    out = dict(h=h, t0=t0)
    for name, a in (('0.0045', lo), ('0.006', 0.006 if not quick else 0.009), ('0.01', 0.01 if not quick else 0.012)):
        out['halving_max_from_' + name] = cmp_max(t1, b1, t2, b2, a, tmax)
    log('step halving h=%g -> %g: max change on [%g,%g] %.2e ; t>=%s: %.2e ; t>=%s: %.2e' % (h, h / 2, lo, tmax, out['halving_max_from_0.0045'], '0.006', out['halving_max_from_0.006'], '0.01', out['halving_max_from_0.01']))
    # sensitivity to t0 and to b(t0)
    t0b = t0 - (10 if not quick else 2) * h
    t3, b3, _ = run_mu(ctx, h, t0b, tmax)
    t4, b4, _ = run_mu(ctx, h, t0, tmax, b_off=0.01)
    out['t0_change_max'] = cmp_max(t1, b1, t3, b3, lo, tmax)
    out['boff_change_max'] = cmp_max(t1, b1, t4, b4, lo, tmax)
    log('t0 -> %g: max change %.2e on [%g, %g]; b(t0) + 0.01: max change %.2e' % (t0b, out['t0_change_max'], lo, tmax, out['boff_change_max']))
    # enclosure
    lo_b, hi_b, inside, wid = [], [], True, {}
    for t, b, F in zip(t1, b1, Fv):
        lp = float(mp.log(F))
        lb = math.sqrt(t) * phibar_inv_log(lp)
        ub = math.sqrt(t) * phibar_inv_log(lp + math.log(0.5) * 0 + math.log(0.5)) if False else math.sqrt(t) * phibar_inv_log(float(mp.log(F / 2)))
        if t >= lo - 1e-12:
            inside = inside and (lb - 1e-9 <= b <= ub + 1e-9)
            for tt in (0.0045, 0.01) if not quick else (0.007, 0.012):
                if abs(t - tt) < 1e-9:
                    wid[str(tt)] = ub - lb
    out['inside_enclosure_all_nodes'] = bool(inside); out['widths'] = wid
    out['decreasing'] = bool(np.all(np.diff(b1[t1 >= lo - 1e-12]) < 0))
    log('enclosure of Proposition smallb holds at every node: %s; widths %s; b decreasing: %s' % (inside, wid, out['decreasing']))
    # reference values
    out['b_at'] = {str(tt): float(b1[np.argmin(abs(t1 - tt))]) for tt in (0.0045, 0.01, 0.02, 0.05, 0.08, 0.12) if tt <= tmax}
    out['slope_end'] = float((b1[-1] - b1[-6]) / (t1[-1] - t1[-6]))
    log('b at nodes: %s ; slope at the end of the grid %.3f' % (out['b_at'], out['slope_end']))
    np.save(str(Path(R['_out']).with_suffix('')) + '_b_h%g.npy' % h, np.vstack([t1, b1]))
    R['solve'] = out
    return t1, b1


def task_resp(R, quick, ctx, tb):
    S, dps, zs, base, lf = ctx
    h, t0, tmax = S['h'], S['t0'], S['tmax']
    t1, b1 = tb if tb else (None, None)
    out = {}
    delta = 0.3085
    rows = [(1, 0.05), (1, 0.15), (1, 0.3085), (6, 0.3085), (20, 0.3085), (60, 0.3085)]
    if quick:
        rows = [(6, 0.3085), (20, 0.3085)]
    if ROWS is not None:
        rows = [(int(x.split(":")[0]), float(x.split(":")[1])) for x in ROWS.split(",") if x]
    ctl_cache = {}
    lo = 0.0045 if not quick else 0.007
    norm = {}
    for m, dl in rows:
        if m > len(zs):
            continue
        if m not in ctl_cache:
            with mp.workdps(dps):
                lawc = toy_law(base, zs, m, dl, "ctl")
            ctl_cache[m] = run_mu(ctx, h, t0, tmax, law=lawc)
        tc, bc, _ = ctl_cache[m]
        with mp.workdps(dps):
            lawt = toy_law(base, zs, m, dl, "toy")
        tt, bt, _ = run_mu(ctx, h, t0, tmax, law=lawt)
        hi = tmax if m != 60 else min(0.1, tmax)
        mx = 0.0
        for t, x, y in zip(tt, bt, bc):
            if lo - 1e-12 <= t <= hi + 1e-12:
                mx = max(mx, abs(x - y))
        # negativity onset of the density
        lawd = lawt
        grid = [S['split'] + 0.005 * i for i in range(1, 1200)] if not quick else [S['split'] + 0.01 * i for i in range(1, 600)]
        onset = None
        with mp.workdps(dps):
            for y in grid:
                if lawd.f(mp.mpf(y)) < 0:
                    onset = y; break
        neg_small = None
        if onset is None:
            for tq in ([0.004, 0.005, 0.006] if not quick else [0.006, 0.0095]):
                v, _ = inv_vertical(lawd.lt, tq, dens=True, P=0.04)
                if v < 0:
                    neg_small = tq
        key = 'gamma_%d_delta_%g' % (m, dl)
        out[key] = dict(max_diff=mx, onset_negative_from=onset, negative_below_split=neg_small, hi=hi)
        log('tab:ifpt row zero %d delta %g: max |b_toy - b_GGC| on [%g, %g] = %.3e ; density negative from y = %s (grid step 0.005 above %g; below it %s)' % (m, dl, lo, hi, mx, onset, S['split'], neg_small))
        # normalised differences for tab:resp
        if dl == delta and m in (6, 20, 60):
            vals = []
            for tq in (0.01, 0.02, 0.05, 0.08):
                if tq > tmax:
                    continue
                i = int(np.argmin(abs(tt - tq)))
                g = float(zs[m - 1])
                vals.append((tq, (bt[i] - bc[i]) * g ** 4 / (6 * delta ** 2)))
            norm[m] = vals
            log('tab:resp zero %d normalised (b_toy - b_ctl) g^4/(6 delta^2) at t=0.01,0.02,0.05,0.08: %s' % (m, ['%.2f' % v for _, v in vals]))
    out['normalised'] = {str(k): v for k, v in norm.items()}
    if not DOK:
        R["resp"] = out
        return
    # K: response to a shift of the law, Delta = 2e-6 and 4e-6
    Ks = {}
    for Dl in (2e-6, 4e-6):
        tK, bK, _ = run_mu(ctx, h, t0, tmax, shift=Dl)
        K = (bK - b1) / Dl
        Ks[Dl] = (tK, K)
    for tq in (0.01, 0.02, 0.05, 0.08):
        if tq > tmax:
            continue
        i = int(np.argmin(abs(t1 - tq)))
        log('K(%g) = %.2f (Delta 2e-6), %.2f (Delta 4e-6)' % (tq, Ks[2e-6][1][i], Ks[4e-6][1][i]))
    sel = (t1 >= lo - 1e-12)
    out['K'] = {str(tq): [float(Ks[2e-6][1][int(np.argmin(abs(t1 - tq)))]), float(Ks[4e-6][1][int(np.argmin(abs(t1 - tq)))])] for tq in (0.01, 0.02, 0.05, 0.08) if tq <= tmax}
    out['maxabsK_2e-6'] = float(np.max(np.abs(Ks[2e-6][1][sel]))); out['maxabsK_4e-6'] = float(np.max(np.abs(Ks[4e-6][1][sel])))
    log('max |K| over [%g, %g]: %.1f (Delta 2e-6), %.1f (Delta 4e-6); argmax at t = %.4f' % (lo, tmax, out['maxabsK_2e-6'], out['maxabsK_4e-6'], t1[sel][np.argmax(np.abs(Ks[2e-6][1][sel]))]))
    R['resp'] = out


def solve_nonuniform(tg, Fv, Sv, dps, switch=True):
    return solve_boundary(tg, Fv, Sv, dps, switch=switch)


def task_tail(R, quick, ctx):
    """mu to t = 2 by the complementary form, t(b' + c); exponential law to t = 150."""
    S, dps, zs, base, lf = ctx
    out = {}
    c = math.sqrt(2) * float(zs[0])
    h1 = S['h']
    hs = (2e-3, 1e-3) if not quick else (1e-2, 5e-3)
    t_end = 2.0 if not quick else 0.5
    res = {}
    for h2 in hs:
        tg = [S['t0'] + k * h1 for k in range(int(round((0.12 - S['t0']) / h1)) + 1)] if not quick else [S['t0'] + k * h1 for k in range(int(round((0.06 - S['t0']) / h1)) + 1)]
        t_last = tg[-1]
        tg += [t_last + k * h2 for k in range(1, int(round((t_end - t_last) / h2)) + 1)]
        Fv, Sv = Fvals(Law(base), tg, S['split'], dps)
        b = solve_boundary(np.array(tg), Fv, Sv, dps)
        res[h2] = (np.array(tg), b)
        log('mu tail solved, h=%g, nodes %d, b(%g) = %.6f' % (h2, len(tg), t_end, b[-1]))
    tg, b = res[hs[1]]
    ts = (0.8, 1.2, 1.6, 2.0) if not quick else (0.2, 0.3, 0.4, 0.5)
    vals = {}
    for tq in ts:
        for (name, (tgx, bx)) in (('fine', res[hs[1]]), ('coarse', res[hs[0]])):
            i = int(np.argmin(abs(tgx - tq)))
            k = max(1, int(round(0.02 / hs[1] if name == 'fine' else 0.02 / hs[0])))
            sl = (bx[i + k] - bx[i - k]) / (tgx[i + k] - tgx[i - k]) if i + k < len(bx) else (bx[i] - bx[i - k]) / (tgx[i] - tgx[i - k])
            vals[(tq, name)] = (tgx[i] * (sl + c), sl)
        log('mu: t(b\'+c) at t=%g: %.4f (fine h=%g), %.4f (coarse h=%g)' % (tq, vals[(tq, 'fine')][0], hs[1], vals[(tq, 'coarse')][0], hs[0]))
    out['mu'] = {str(tq): [vals[(tq, 'fine')][0], vals[(tq, 'coarse')][0]] for tq in ts}
    out['c'] = c
    # slope on [0.03, 0.12]
    tgx, bx = res[hs[1]]
    sel = (tgx >= 0.03) & (tgx <= 0.12) if not quick else (tgx >= 0.02) & (tgx <= 0.06)
    out['slope_0.03_0.12'] = float(np.polyfit(tgx[sel], bx[sel], 1)[0])
    log('slope of b on [0.03, 0.12] (least squares): %.3f' % out['slope_0.03_0.12'])
    # exponential law
    tex = 150.0 if not quick else 20.0
    ex = {}
    for hh in ((0.02, 0.01) if not quick else (0.05, 0.025)):
        t0e = hh
        tg = [t0e * k for k in range(1, int(round(tex / hh)) + 1)]
        with mp.workdps(60):
            Fv = [1 - mp.exp(-mp.mpf(t)) for t in tg]
            Sv = [mp.exp(-mp.mpf(t)) for t in tg]
        b = solve_boundary(np.array(tg), Fv, Sv, 60)
        ex[hh] = (np.array(tg), b)
        i = int(np.argmax(b))
        log('exponential law, h=%g: max of b = %.4f at t = %.3f' % (hh, b[i], tg[i]))
    qs = (20, 50, 100, 150) if not quick else (5, 10, 15, 20)
    exo = {}
    for tq in qs:
        row = []
        for hh, (tgx, bx) in ex.items():
            i = int(np.argmin(abs(tgx - tq))); k = max(1, int(round(0.1 * tq / hh / 10)))
            sl = (bx[i + k] - bx[i - k]) / (tgx[i + k] - tgx[i - k]) if i + k < len(bx) else (bx[i] - bx[i - k]) / (tgx[i] - tgx[i - k])
            row.append(tgx[i] * (sl + math.sqrt(2)))
        exo[str(tq)] = row
        log('exponential law t(b\'+sqrt2) at t=%g: %s (h = %s)' % (tq, ['%.4f' % v for v in row], list(ex.keys())))
    out['exp'] = exo
    R['tail'] = out


def main():
    global LOG, ZFILE, ROWS, DOK
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', default='all')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default=str(HERE.parent / 'data' / 'ht_curved'))
    ap.add_argument('--zeros', default=str(HERE.parent / 'data' / 'zeros_1700.txt'))
    ap.add_argument("--rows", default=None, help="resp rows m:delta,... (default all); an empty string runs no rows")
    ap.add_argument("--k", type=int, default=1, help="1: also compute K")
    a = ap.parse_args()
    ZFILE = a.zeros
    ROWS = a.rows
    DOK = bool(a.k)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    LOG = open(a.out + '.log', 'w')
    R = {'_out': a.out, 'quick': a.quick, 'task': a.task}
    tasks = ['smallt', 'formula', 'lin', 'cdf', 'solve', 'resp', 'tail'] if a.task == 'all' else a.task.split(',')
    mp.mp.dps = 30
    ctx = None
    tb = None
    for tk in tasks:
        log('=== task', tk)
        if tk == 'smallt': task_smallt(R, a.quick)
        elif tk == 'formula': task_formula(R, a.quick)
        elif tk == 'lin': task_lin(R, a.quick)
        else:
            if ctx is None:
                ctx = build(a.quick)
            if tk == 'cdf': task_cdf(R, a.quick, ctx)
            elif tk == 'solve': tb = task_solve(R, a.quick, ctx)
            elif tk == 'resp':
                if tb is None and DOK:
                    S = ctx[0]; t1, b1, _ = run_mu(ctx, S['h'], S['t0'], S['tmax']); tb = (t1, b1)
                task_resp(R, a.quick, ctx, tb)
            elif tk == 'tail': task_tail(R, a.quick, ctx)
    json.dump(R, open(a.out + '.json', 'w'), indent=1, default=str)
    log('done')


if __name__ == '__main__':
    main()
