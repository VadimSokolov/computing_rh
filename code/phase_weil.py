# code/phase_weil.py (October 2026): Section cp:sec:weil ("The Weil Gram matrix on compactly supported test functions") of Chapter 14
# ("The Thorin phase diagram"): Table cp:tab:weil, Table cp:tab:weilcheck, Figure cp:fig:weil (fig/cp_weil.pdf), and the two statements of
# that section that had no output (the Gaussian check of the normalisation of the explicit formula, and the error of the r space integral
# truncated at r = 400). Written from the mathematics, not from the authors' scripts part3.py, part3_check.py, ef_test.py, figs_fa.py.
# Weil's formula (Section ch:clock) for h with Fourier transform hhat(r) = int h(x) e^{irx} dx, rho = 1/2 + i gamma:
#   sum_rho hhat(gamma) = hhat(i/2) + hhat(-i/2) - sum_{n>=2} Lambda(n) n^{-1/2} (h(log n) + h(-log n)) + (1/2pi) int hhat(r) Re psi(1/4 + ir/2) dr - h(0) log pi.
# Test functions: hat functions g_j(x) = max(0, 1 - |x - c_j|/w) at c_j = -L/2 + j w, j = 1..N, N = L/w - 1, h_jk = g_j * g_k~. Then
#   h_jk(x) = Phi(x - d),  d = c_j - c_k = (j - k) w,  Phi(s) = w A(|s|/w),  A(a) = 2/3 - a^2 + a^3/2 (a <= 1), (2 - a)^3/6 (1 <= a <= 2), 0 beyond,
#   hhat_jk(r) = e^{ird} Phihat(r),  Phihat(r) = w^2 (sin(rw/2)/(rw/2))^4,
# so the Gram matrix G_jk = W(h_jk) is a Toeplitz matrix, G_jk = T_{|j-k|}, T_m = W(Phi(. - m w)) = pole_m - prime_m + arch_m, with
#   pole_m  = 2 cosh(mw/2) w^2 (sinh(w/4)/(w/4))^4                                  (hhat(i/2) + hhat(-i/2))
#   prime_m = sum_{n = p^k >= 2} Lambda(n) n^{-1/2} (Phi(log n - mw) + Phi(log n + mw))   (a finite sum: n <= e^L, no truncation)
#   arch_m  = int_0^inf [Phi(-d) e^{-2x}/x - 2 h_e(x) e^{-x/2}/(1 - e^{-2x})] dx - Phi(d) log pi,  h_e(x) = (Phi(x - d) + Phi(x + d))/2,
# the archimedean term in x space (Gauss's integral for psi), which is exact; it is integrated cell by cell between the knots of the cubic
# pieces (multiples of w), where the integrand is analytic (the singularity at x = 0 is removable), plus Phi(d) E_1(2 (m + 2) w) beyond the support.
# G(L) is the leading N x N section of one Toeplitz sequence, so lambda_min decreases with L (Cauchy interlacing) and all L at a resolution share T_m.
# The script computes (all in mpmath, working precision dps + 25 digits, dps = 40, or 25 with --quick):
#  gaussian_check: hhat(r) = e^{-r^2/50} against the first 60 zeros (and all 1700): both sides to about dps digits, and the x space against the r space
#              archimedean term (cp6.tex:11, "to twenty digits");
#  rspace_check:   the archimedean diagonal term T_0 by the r integral truncated at r = 400 against the x space value (cp6.tex:15, "of order 1e-3");
#  toeplitz:       pole, prime, archimedean parts of T_m for each resolution w, the quadrature error estimates, a second quadrature (tanh-sinh at
#              higher precision) for a few m, and T_m against the zero sum over 1700 zeros plus the smooth tail estimate (an entry level check);
#  gram:           for w = 0.25, 0.1 and L = 1..8: lambda_min, lambda_max, the pole, -prime, archimedean parts of v*Gv for the minimising unit vector v
#              (Table cp:tab:weil), the cancellation depth, the same parts for the vector v_j ~ c_j (the pole term is negative there), and
#              check: the sum of |ghat_min(gamma)|^2 over the first 1500 and 1700 zeros, the shares of the zeros, and the estimate of the omitted zeros
#              from the smooth zero density N(T) = (T/2pi) log(T/2pi e) + 7/8 (Table cp:tab:weilcheck, cp6.tex:62);
#  rates:          lambda_min(L)/lambda_min(L+1) at each resolution ("by a factor between about 2.5 and 20 per unit of L");
#  control:        the pair of zeros (1,2), (11,12), (51,52) moved to the quadruple 1/2 +- 1/4 +- ig, g = the mean of the pair (Figure cp:fig:weil, right):
#              lambda_min of G + Delta, Delta_m = sum over the quadruple of e^{i gamma d} Phihat(gamma) minus the same for the four zeros removed,
#              at w = 0.2, 0.1, 0.05, 0.025 for L = 1..6 and w = 0.0125 for L = 1..4 (the figure takes w = 0.05 and L <= 4, as the chapter does; the other
#              resolutions show how the detection depends on w, the law "w <~ 1/g" of cp6.tex:62; control.first_L_detected lists the smallest L detecting each pair).
# Error control. The prime sum is finite, so the Table entries have no truncation error: only the quadrature of arch_m (Gauss-Legendre cell by cell with
# the mp.quad error estimates, recorded, and a second rule, tanh-sinh at 15 more digits, on a few entries) and rounding (a guard of 25 digits; the
# eigenvalues come from mp.eigsy, with the residual |Gv - lambda v| recorded). The zero sum is a check, not an input: its truncation at 1700 zeros is
# estimated, not bounded (a bound for the minimising vector through |ghat| <= 4 sum|v_j|/(w g^2) is of order N/w^2 times 1e-9, far above lambda_min).
# The estimate integrates 2 Phihat(g) |sum_j v_j e^{i g c_j}|^2 against the smooth density log(g/2pi)/2pi of the zeros, from T_K with N(T_K) = K + 1/2
# (Simpson's rule up to g = Gamma = 40000, then the mean of sin^4(g w/2)|.|^2 over a period, which is |v * (1,-2,1)|^2/16, not |v|^2 3/8: the vector is
# smooth and the zeros of Phihat at g = 2 pi k/w kill the aliasing peaks of |.|^2). It is calibrated on the zeros 1501..1700, which are in the file
# (tail_est_1500_to_1700 against tail_actual_1500_to_1700), and what it leaves unexplained is reported (deficit_over_tail = 1 would be a perfect estimate).
# Run from the repository root (or from a folder holding zeros_1700.txt): python3 code/phase_weil.py [--quick] [--out FILE] [--fig FILE] [--workers K];
# writes data/phase_weil.json and fig/cp_weil.pdf; --quick (a reduced case, about 15 seconds on one process) writes /tmp/phase_weil_quick.json and .pdf;
# --figure redraws the figure from the JSON given by --out without computing. The full case (about 10 minutes of CPU, 3 to 4 minutes on 8 workers) is
# for the Hopper cluster (misc/reviews/nick-ch14/compute-C3-run.sh); --workers K runs the independent tasks (T_m, Gram matrices, control eigenvalues) in K processes.
import argparse
import json
import math
import os
import platform
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAIRS = [(1, 2), (11, 12), (51, 52)]   # the zeros (1 based) of the three displaced pairs of Figure cp:fig:weil
KS = (1500, 1700)                      # zeros summed in the check: 1500 as in the caption of Table cp:tab:weilcheck, 1700 = the whole file
CFG = None
PP = None
ZMP = None


def config(quick):
    if quick:
        return dict(quick=True, dps=25,
                    gram={'0.25': [1, 2, 3, 4, 5, 6], '0.1': [1, 2, 3, 4, 5, 6]},
                    check=[('0.25', 1), ('0.25', 3), ('0.25', 6), ('0.1', 2), ('0.1', 4), ('0.1', 6)],
                    ctrl={'0.2': [1, 2], '0.1': [1, 2, 3, 4], '0.05': [1, 2]}, fig_w='0.1',
                    rspace_w=['0.25', '0.1'], toep=[('0.25', [0, 1, 2, 3]), ('0.1', [0, 1, 5])], second_quad=[('0.25', [0, 1, 3])],
                    Gamma=4000.0, h=0.1)
    gram = {'0.25': list(range(1, 9)), '0.1': list(range(1, 9))}
    return dict(quick=False, dps=40, gram=gram,
                check=[(ws, L) for ws in gram for L in gram[ws] if L <= 6],
                ctrl={'0.2': [1, 2, 3, 4, 5, 6], '0.1': [1, 2, 3, 4, 5, 6], '0.05': [1, 2, 3, 4, 5, 6], '0.025': [1, 2, 3, 4, 5, 6],
                      '0.0125': [1, 2, 3, 4]}, fig_w='0.05',
                rspace_w=['0.25', '0.1', '0.05', '0.01'], toep=[('0.25', [0, 1, 2, 3, 5, 10]), ('0.1', [0, 1, 2, 5, 20, 50])],
                second_quad=[('0.25', [0, 1, 3, 10]), ('0.1', [0, 1, 2, 20])], Gamma=40000.0, h=0.05)


def fr2mp(fr):
    return mp.mpf(fr.numerator) / fr.denominator


def nodes_of(ws, Ls):
    w, L = Fraction(ws), Fraction(str(Ls))
    return int(L / w) - 1


def ser(x):
    return mp.nstr(x, mp.mp.dps + 3)


def js(x, n=25):
    return mp.nstr(x, n, min_fixed=-3)


def init(cfg):
    global CFG, PP
    CFG = cfg
    mp.mp.dps = cfg['dps'] + 25
    Lmax = max(max(v) for v in cfg['gram'].values())
    Lmax = max([Lmax] + [max(v) for v in cfg['ctrl'].values()])
    nmax = int(math.exp(Lmax)) + 2
    sieve = bytearray([1]) * (nmax + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(nmax ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, nmax + 1, i)))
    pp = []
    for p in range(2, nmax + 1):
        if sieve[p]:
            q = p
            while q <= nmax:
                pp.append((q, p))
                q *= p
    pp.sort()
    PP = dict(n=[a for a, _ in pp], lam=[mp.log(b) for _, b in pp], ln=[mp.log(a) for a, _ in pp],
              lnf=np.array([math.log(a) for a, _ in pp]))


def zeros_mp():
    global ZMP
    if ZMP is None:
        with open(CFG['zeros']) as fh:
            ZMP = [mp.mpf(s) for s in fh.read().split()]
    return ZMP


def Ahat(a):
    """(unit hat) * (unit hat)~ at lag a >= 0, in units of w: piecewise cubic."""
    a = mp.mpf(a)
    if a <= 1:
        return mp.mpf(2) / 3 - a * a + a ** 3 / 2
    if a <= 2:
        return (2 - a) ** 3 / 6
    return mp.mpf(0)


def Phihat(r, w):
    x = r * w / 2
    if x == 0:
        return w * w
    return w * w * (mp.sin(x) / x) ** 4


def pole_term(w, m):
    return 2 * mp.cosh(m * w / 2) * w * w * (mp.sinh(w / 4) / (w / 4)) ** 4


def prime_term(w, m):
    d = m * w
    s = mp.mpf(0)
    lnf = PP['lnf']
    lo, hi = float(d - 2 * w) - 1e-9, float(d + 2 * w) + 1e-9
    i0, i1 = int(np.searchsorted(lnf, lo, 'left')), int(np.searchsorted(lnf, hi, 'right'))
    i2 = int(np.searchsorted(lnf, float(2 * w - d) + 1e-9, 'right')) if d < 2 * w else 0
    for i in sorted(set(range(i0, i1)) | set(range(0, i2))):
        ln = PP['ln'][i]
        a = Ahat(abs(ln - d) / w) + Ahat((ln + d) / w)
        if a != 0:
            s += PP['lam'][i] / mp.sqrt(PP['n'][i]) * w * a
    return s


def arch_term(w, m, method='gauss-legendre'):
    """arch_m = int_0^inf [Phi(d) e^{-2x}/x - 2 h_e(x) e^{-x/2}/(1 - e^{-2x})] dx - Phi(d) log pi; returns (value, sum of the quad error estimates)."""
    d = m * w
    h0 = w * Ahat(m)

    def f(x):
        he = w * (Ahat(abs(x - d) / w) + Ahat((x + d) / w)) / 2
        return h0 * mp.exp(-2 * x) / x - 2 * he * mp.exp(-x / 2) / (-mp.expm1(-2 * x))
    knots = sorted(set(abs(m + k) for k in range(-2, 3)))
    tot, err = mp.mpf(0), mp.mpf(0)
    for a, b in zip(knots[:-1], knots[1:]):
        v, e = mp.quad(f, [a * w, b * w], method=method, error=True)
        tot += v
        err += e
    if h0 != 0:
        tot += h0 * (mp.e1(2 * (m + 2) * w) - mp.log(mp.pi))
    return tot, err


def task_T(a):
    t0 = time.time()
    ws, m = a
    w = fr2mp(Fraction(ws))
    ar, er = arch_term(w, m)
    return ws, m, ser(pole_term(w, m)), ser(prime_term(w, m)), ser(ar), js(er, 5), time.time() - t0


def task_T2(a):
    """T_m by tanh-sinh at 15 more digits: the second quadrature."""
    t0 = time.time()
    ws, m = a
    w = fr2mp(Fraction(ws))
    with mp.extradps(15):
        w = fr2mp(Fraction(ws))
        ar, er = arch_term(w, m, method='tanh-sinh')
    return ws, m, ser(ar), time.time() - t0


def toeplitz(T, N):
    G = mp.matrix(N, N)
    for j in range(N):
        for k in range(N):
            G[j, k] = T[abs(j - k)]
    return G


def form(X, v):
    """v^T X v for the Toeplitz matrix with first row X, through the autocorrelation of v."""
    N = len(v)
    s = X[0] * mp.fsum(x * x for x in v)
    for m in range(1, N):
        s += 2 * X[m] * mp.fsum(v[j] * v[j + m] for j in range(N - m))
    return s


def smooth_T(Nt):
    """T with (T/2pi) log(T/2pi e) + 7/8 = Nt."""
    lo, hi = 10.0, 1e6
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid / (2 * math.pi) * math.log(mid / (2 * math.pi * math.e)) + 0.875 < Nt:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def tail_integral(weight, M, wf, Ta, Tb, Gam, h):
    """int_{Ta}^{Tb} 2 Phihat(g) weight(g) rho(g) dg, rho = log(g/2pi)/2pi (the smooth density of the zeros), by Simpson's rule; Tb None: up to Gam, plus the
    remainder beyond Gam, in which sin^4(g w/2) weight(g) is replaced by its mean M over a period 2 pi/w, using 2 Phihat(g) = 32 sin^4(g w/2)/(w^2 g^4).
    For weight = |sum_j v_j e^{i g c_j}|^2, M = |v * (1,-2,1)|^2/16, small for a smooth v; for weight = cos(g m w), M = 3/8, -1/4, 1/16 for m = 0, 1, 2, else 0."""
    top = Gam if Tb is None else Tb
    n = int(math.ceil((top - Ta) / h))
    n += n % 2
    g = Ta + (top - Ta) * np.arange(n + 1) / n
    x = g * wf / 2
    f = 2 * wf * wf * (np.sin(x) / x) ** 4 * weight(g) * np.log(g / (2 * math.pi)) / (2 * math.pi)
    hh = (top - Ta) / n
    val = hh / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum())
    if Tb is None:
        val += M * 32 / (2 * math.pi * wf ** 2) * (math.log(Gam / (2 * math.pi)) / (3 * Gam ** 3) + 1 / (9 * Gam ** 3))
    return val


def hat_mean(vf):
    """mean over a period of sin^4(g w/2) |sum_j v_j e^{i g c_j}|^2 = |v * (1, -2, 1)|^2 / 16 (sin^4(x) = |1 - e^{2ix}|^4/16)."""
    return float((np.convolve(vf, [1.0, -2.0, 1.0]) ** 2).sum() / 16)


def vec_weight(vf, wf):
    def weight(g):
        z = np.exp(1j * g * wf)
        acc = np.full(g.shape, vf[-1], dtype=complex)
        for j in range(len(vf) - 2, -1, -1):
            acc = acc * z + vf[j]
        return acc.real ** 2 + acc.imag ** 2
    return weight


def zero_contrib(v, w, K):
    """contrib_k = |ghat(gamma_k)|^2 + |ghat(-gamma_k)|^2 = 2 Phihat(gamma_k) |sum_j v_j e^{i gamma_k c_j}|^2 for k = 1..K (v real)."""
    N = len(v)
    out = []
    for g in zeros_mp()[:K]:
        z = mp.expj(g * w)
        acc = mp.mpc(v[N - 1])
        for j in range(N - 2, -1, -1):
            acc = acc * z + v[j]
        out.append(2 * Phihat(g, w) * (acc.real ** 2 + acc.imag ** 2))
    return out


def task_gram(a):
    t0 = time.time()
    ws, L, comp, do_check = a
    w, Lf = Fraction(ws), Fraction(str(L))
    N = int(Lf / w) - 1
    wm = fr2mp(w)
    P = [mp.mpf(c[0]) for c in comp[:N]]
    Q = [mp.mpf(c[1]) for c in comp[:N]]
    Ar = [mp.mpf(c[2]) for c in comp[:N]]
    T = [P[m] - Q[m] + Ar[m] for m in range(N)]
    G = toeplitz(T, N)
    t1 = time.time()
    E, V = mp.eigsy(G)
    t_eig = time.time() - t1
    order = sorted(range(N), key=lambda i: E[i])
    lam = [E[i] for i in order]
    v = [V[j, order[0]] for j in range(N)]
    nrm = mp.sqrt(mp.fsum(x * x for x in v))
    v = [x / nrm for x in v]
    res = max(abs(mp.fsum(G[j, k] * v[k] for k in range(N)) - lam[0] * v[j]) for j in range(N))
    pole_p, prime_p, arch_p = form(P, v), -form(Q, v), form(Ar, v)
    c = [-fr2mp(Lf) / 2 + (j + 1) * wm for j in range(N)]
    cn = mp.sqrt(mp.fsum(x * x for x in c))
    u = [x / cn for x in c]
    rec = dict(w=ws, L=str(L), N=N, lambda_min=js(lam[0]), lambda_2=js(lam[1]) if N > 1 else None, lambda_max=js(lam[-1]),
               pole=js(pole_p), minus_prime=js(prime_p), arch=js(arch_p), parts_sum=js(pole_p + prime_p + arch_p),
               cancellation_digits=float(-mp.log10(abs(lam[0]) / max(abs(pole_p), abs(prime_p), abs(arch_p)))),
               eig_residual=js(res, 5),
               linear_vector=dict(pole=js(form(P, u)), minus_prime=js(-form(Q, u)), arch=js(form(Ar, u)), form=js(form(T, u))))
    if do_check:
        Kmax = max(KS)
        contrib = zero_contrib(v, wm, Kmax)
        s = {K: mp.fsum(contrib[:K]) for K in KS}
        tot = s[1700]
        share = [x / tot for x in contrib]
        top = sorted(range(Kmax), key=lambda i: -contrib[i])[:3]
        cum, half = mp.mpf(0), None
        for i, x in enumerate(contrib):
            cum += x
            if half is None and cum >= tot / 2:
                half = i + 1
        vf = np.array([float(x) for x in v])
        wf = float(wm)
        wt = vec_weight(vf, wf)
        Ta, Tb = smooth_T(1500.5), smooth_T(1700.5)
        est_mid = tail_integral(wt, 1.0, wf, Ta, Tb, CFG['Gamma'], CFG['h'])
        est_inf = tail_integral(wt, hat_mean(vf), wf, Tb, None, CFG['Gamma'], CFG['h'])
        act_mid = mp.fsum(contrib[1500:1700])
        zc = s[1700] + mp.mpf(est_inf)
        rec['check'] = dict(
            zero_sum_K1500=js(s[1500]), zero_sum_K1700=js(s[1700]), ratio_K1500=js(s[1500] / lam[0], 12), ratio_K1700=js(s[1700] / lam[0], 12),
            tail_est_1500_to_1700=js(est_mid, 6), tail_actual_1500_to_1700=js(act_mid, 6), tail_est_after_1700=js(est_inf, 6),
            zero_sum_K1700_plus_tail=js(zc), ratio_corrected=js(zc / lam[0], 12),
            deficit_over_tail=js((lam[0] - s[1700]) / est_inf, 6) if est_inf != 0 else None,
            half_mass_zero=half,
            top_zeros=[dict(index=i + 1, gamma=js(zeros_mp()[i], 12), share=js(share[i], 6)) for i in top])
    rec['seconds'] = dict(total=time.time() - t0, eigsy=t_eig)
    return rec


def task_ctrl(a):
    t0 = time.time()
    ws, L, k, Ts, zA, zB = a
    w, Lf = Fraction(ws), Fraction(str(L))
    N = int(Lf / w) - 1
    wm = fr2mp(w)
    T = [mp.mpf(s) for s in Ts[:N]]
    if k is not None:
        A, B = mp.mpf(zA), mp.mpf(zB)
        g = (A + B) / 2
        Pp = Phihat(mp.mpc(g, mp.mpf(1) / 4), wm)
        PA, PB = Phihat(A, wm), Phihat(B, wm)
        for m in range(N):
            d = m * wm
            B1 = mp.exp(-d / 4) * mp.expj(g * d) + mp.exp(d / 4) * mp.expj(-g * d)
            T[m] += 2 * (Pp * B1).real - 2 * PA * mp.cos(A * d) - 2 * PB * mp.cos(B * d)
    E = mp.eigsy(toeplitz(T, N), eigvals_only=True)
    return ws, str(L), k, js(min(E[i] for i in range(N))), time.time() - t0


def gaussian_check(zs):
    t0 = time.time()
    dps = CFG['dps']
    out = {}
    with mp.workdps(dps + 10):
        a = mp.mpf(1) / 50
        hhat = lambda r: mp.exp(-a * r * r)
        h = lambda x: mp.exp(-x * x / (4 * a)) / (2 * mp.sqrt(mp.pi * a))
        gam = [mp.mpf(s) for s in zs]
        lhs60, lhsall = 2 * mp.fsum(hhat(g) for g in gam[:60]), 2 * mp.fsum(hhat(g) for g in gam)
        pole = hhat(mp.mpc(0, mp.mpf(1) / 2)).real + hhat(mp.mpc(0, -mp.mpf(1) / 2)).real
        prime = mp.mpf(0)
        for i, n in enumerate(PP['n']):
            if n < 400:
                prime += PP['lam'][i] / mp.sqrt(n) * 2 * h(PP['ln'][i])
        h0 = h(0)
        f = lambda x: h0 * mp.exp(-2 * x) / x - 2 * h(x) * mp.exp(-x / 2) / (-mp.expm1(-2 * x))
        ax = mp.quad(f, [0, 1, 3, 6], method='gauss-legendre') + h0 * mp.e1(12) - h0 * mp.log(mp.pi)
        g2 = lambda r: hhat(r) * mp.re(mp.digamma(mp.mpf(1) / 4 + mp.mpc(0, r / 2)))
        ar = 2 * mp.quad(g2, [0, 10, 20, 30, 40, 60, 80], method='gauss-legendre') / (2 * mp.pi) - h0 * mp.log(mp.pi)
        rhs_x, rhs_r = pole - prime + ax, pole - prime + ar
        dig = lambda x, y: float(-mp.log10(abs(x - y) / abs(x))) if x != y else None
        out = dict(hhat='exp(-r^2/50)', zeros_60=js(lhs60, 30), zeros_all=js(lhsall, 30), pole=js(pole, 30), prime=js(prime, 30),
                   arch_xspace=js(ax, 30), arch_rspace=js(ar, 30), arithmetic_xspace=js(rhs_x, 30), arithmetic_rspace=js(rhs_r, 30),
                   digits_zeros60_vs_xspace=dig(lhs60, rhs_x), digits_zeros60_vs_rspace=dig(lhs60, rhs_r), digits_xspace_vs_rspace=dig(ax, ar))
    out['seconds'] = time.time() - t0
    return out


def rspace_check(ws, R=400):
    w = fr2mp(Fraction(ws))
    h0 = 2 * w / 3
    ax, _ = arch_term(w, 0)
    f = lambda r: w * w * (mp.sin(r * w / 2) / (r * w / 2)) ** 4 * mp.re(mp.digamma(mp.mpf(1) / 4 + mp.mpc(0, r / 2)))
    ar = mp.quad(f, [25 * i for i in range(R // 25 + 1)], method='gauss-legendre') / mp.pi - h0 * mp.log(mp.pi)
    return dict(w=ws, R=R, arch_xspace=js(ax), arch_rspace_truncated=js(ar), error=js(ar - ax, 6))


def toeplitz_zero_check(ws, ms, comp, zs):
    """T_m against 2 sum_k cos(gamma_k m w) Phihat(gamma_k) over the 1700 zeros plus the smooth tail estimate."""
    w = fr2mp(Fraction(ws))
    wf = float(w)
    gam = [mp.mpf(s) for s in zs]
    rows = []
    for m in ms:
        T = mp.mpf(comp[m][0]) - mp.mpf(comp[m][1]) + mp.mpf(comp[m][2])
        zsum = 2 * mp.fsum(mp.cos(g * m * w) * Phihat(g, w) for g in gam)
        weight = (lambda mm: (lambda g: np.cos(g * mm * wf)))(m)
        tail = tail_integral(weight, {0: 3 / 8, 1: -1 / 4, 2: 1 / 16}.get(m, 0.0), wf, smooth_T(1700.5), None, CFG['Gamma'], CFG['h'])
        rows.append(dict(m=m, arithmetic=js(T), zero_sum_1700=js(zsum), tail_est=js(tail, 6), difference=js(T - zsum, 6),
                         difference_minus_tail=js(T - zsum - tail, 6)))
    return dict(w=ws, rows=rows)


def run(cfg, workers):
    t00 = time.time()
    init(cfg)
    zs = open(cfg['zeros']).read().split()
    pool = Pool(workers, initializer=init, initargs=(cfg,)) if workers > 1 else None

    def pmap(fn, tasks):
        if pool is None:
            return [fn(t) for t in tasks]
        return list(pool.imap_unordered(fn, tasks, chunksize=1))
    sec, cpu = {}, {}
    res = dict(meta=dict(script='code/phase_weil.py', quick=cfg['quick'], dps=cfg['dps'], working_dps=cfg['dps'] + 25, workers=workers,
                         zeros_file=os.path.basename(cfg['zeros']), n_zeros=len(zs), KS=list(KS), python=platform.python_version(),
                         mpmath=mp.__version__, numpy=np.__version__, tail_Gamma=cfg['Gamma'], tail_h=cfg['h'], fig_w=cfg['fig_w'],
                         command=' '.join(sys.argv), pairs=PAIRS))
    # 1. the Toeplitz sequences
    t = time.time()
    need = {}
    for dct in (cfg['gram'], cfg['ctrl']):
        for ws, Ls in dct.items():
            need[ws] = max(need.get(ws, 0), nodes_of(ws, max(Ls)))
    comp = {ws: [None] * n for ws, n in need.items()}
    qerr = {ws: mp.mpf(0) for ws in need}
    cpu['T'] = 0.0
    for ws, m, p, q, ar, er, dt in pmap(task_T, [(ws, m) for ws, n in need.items() for m in range(n)]):
        comp[ws][m] = (p, q, ar)
        qerr[ws] = max(qerr[ws], mp.mpf(er))
        cpu['T'] += dt
    sec['toeplitz'] = time.time() - t
    # second quadrature on a few entries
    t = time.time()
    sq = []
    for ws, ms in cfg['second_quad']:
        for ws2, m, ar2, dt in pmap(task_T2, [(ws, m) for m in ms]):
            d = abs(mp.mpf(ar2) - mp.mpf(comp[ws2][m][2]))
            sq.append(dict(w=ws2, m=m, arch_glegendre=js(mp.mpf(comp[ws2][m][2])), arch_tanhsinh=js(mp.mpf(ar2)), difference=js(d, 5)))
    res['toeplitz'] = {ws: dict(nmax=need[ws], max_quad_error_estimate=js(qerr[ws], 5),
                                pole=[js(mp.mpf(c[0])) for c in comp[ws]], prime=[js(mp.mpf(c[1])) for c in comp[ws]],
                                arch=[js(mp.mpf(c[2])) for c in comp[ws]]) for ws in need}
    res['second_quadrature'] = sq
    sec['second_quadrature'] = time.time() - t
    # 2. the Gram matrices
    t = time.time()
    gt = [(ws, L, comp[ws], (ws, L) in set(cfg['check'])) for ws, Ls in cfg['gram'].items() for L in Ls]
    gt.sort(key=lambda a: -nodes_of(a[0], a[1]))
    gram, cpu['gram'], cpu['eigsy_vec'] = {}, 0.0, 0.0
    for rec in pmap(task_gram, gt):
        gram['%s,%s' % (rec['w'], rec['L'])] = rec
        cpu['gram'] += rec['seconds']['total']
        cpu['eigsy_vec'] += rec['seconds']['eigsy']
    gram = {k: gram[k] for k in sorted(gram, key=lambda s: (-float(s.split(',')[0]), float(s.split(',')[1])))}
    res['gram'] = gram
    sec['gram'] = time.time() - t
    # 3. rates, interlacing, positivity
    rates = {}
    for ws, Ls in cfg['gram'].items():
        lm = {L: mp.mpf(gram['%s,%s' % (ws, L)]['lambda_min']) for L in Ls}
        Lc = [L for L in Ls if L + 1 in lm]
        rates[ws] = dict(L=Lc, ratio=[js(lm[L] / lm[L + 1], 8) for L in Lc],
                         decreasing=all(lm[L] > lm[L + 1] for L in Lc), all_positive=all(x > 0 for x in lm.values()))
    r15 = [float(mp.mpf(r)) for ws, d in rates.items() for L, r in zip(d['L'], d['ratio']) if L <= 5]
    rates['ratio_min_max_L_le_5'] = [min(r15), max(r15)]
    res['rates'] = rates
    # 4. the control
    t = time.time()
    Tt = {ws: [js(mp.mpf(c[0]) - mp.mpf(c[1]) + mp.mpf(c[2]), mp.mp.dps) for c in comp[ws]] for ws in need}
    zg = [mp.mpf(s) for s in zs]
    pinfo = [dict(index=[a, b], A=zs[a - 1], B=zs[b - 1], g=js((zg[a - 1] + zg[b - 1]) / 2, 12)) for a, b in PAIRS]
    ct = []
    for ws, Ls in cfg['ctrl'].items():
        for L in Ls:
            for k in [None, 0, 1, 2]:
                if k is None and '%s,%s' % (ws, L) in gram:
                    continue
                a, b = PAIRS[k] if k is not None else (1, 2)
                ct.append((ws, L, k, Tt[ws], zs[a - 1], zs[b - 1]))
    ct.sort(key=lambda a: -nodes_of(a[0], a[1]))
    ctrl, cpu['ctrl'] = {}, 0.0
    for ws, L, k, lam, dt in pmap(task_ctrl, ct):
        ctrl[(ws, L, k)] = lam
        cpu['ctrl'] += dt
    rows = {}
    for ws, Ls in cfg['ctrl'].items():
        for L in Ls:
            base = gram['%s,%s' % (ws, L)]['lambda_min'] if '%s,%s' % (ws, L) in gram else ctrl[(ws, str(L), None)]
            row = dict(w=ws, L=str(L), N=nodes_of(ws, L), lambda_min=base)
            for k in range(3):
                lp = mp.mpf(ctrl[(ws, str(L), k)])
                row['pair%d' % (k + 1)] = dict(lambda_min=js(lp), detected=bool(lp < 0), ratio_to_base=js(lp / mp.mpf(base), 8),
                                               w_times_g=float(mp.mpf(ws) * mp.mpf(pinfo[k]['g'])))
            rows['%s,%s' % (ws, L)] = row
    first = {ws: {'pair%d' % (k + 1): next((L for L in Ls if rows['%s,%s' % (ws, L)]['pair%d' % (k + 1)]['detected']), None) for k in range(3)}
             for ws, Ls in cfg['ctrl'].items()}
    res['control'] = dict(pairs=pinfo, rows=rows, first_L_detected=first)
    sec['control'] = time.time() - t
    # 5. the checks that do not need a pool
    t = time.time()
    res['gaussian_check'] = gaussian_check(zs)
    res['rspace_check'] = [rspace_check(ws) for ws in cfg['rspace_w']]
    res['toeplitz_zero_check'] = [toeplitz_zero_check(ws, ms, comp[ws], zs) for ws, ms in cfg['toep']]
    sec['gaussian_rspace_toeplitz_checks'] = time.time() - t
    if pool is not None:
        pool.close()
        pool.join()
    sec['total'] = time.time() - t00
    res['meta']['seconds'] = {k: round(v, 2) for k, v in sec.items()}
    res['meta']['cpu_seconds'] = {k: round(v, 2) for k, v in cpu.items()}
    return res


def draw_figure(res, path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception as e:
        print('no matplotlib, figure not drawn:', e)
        return
    plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
    sty = {'0.25': dict(ls='-', marker='o', color='#0072B2'), '0.1': dict(ls='--', marker='s', color='#D55E00')}
    for ws, kw in sty.items():
        pts = sorted((float(r['L']), float(r['lambda_min'])) for r in res['gram'].values() if r['w'] == ws)
        if pts:
            ax[0].semilogy([p[0] for p in pts], [p[1] for p in pts], ms=4, lw=1.3, label=r'$w=%s$' % ws, **kw)
    ax[0].set_xlabel(r'support length $L$')
    ax[0].set_ylabel(r'$\lambda_{\min}$ of the Weil Gram matrix')
    ax[0].legend(frameon=False)
    fw = res['meta']['fig_w']
    pst = [dict(ls='-', marker='o', color='#0072B2'), dict(ls='--', marker='s', color='#D55E00'), dict(ls='-.', marker='^', color='#009E73')]
    for k, p in enumerate(res['control']['pairs']):
        pts = sorted((float(r['L']), float(r['pair%d' % (k + 1)]['lambda_min'])) for r in res['control']['rows'].values()
                     if r['w'] == fw and float(r['L']) <= 4)
        y = [math.copysign(math.log10(abs(v) / 1e-9), v) for _, v in pts]
        ax[1].plot([x for x, _ in pts], y, ms=4, lw=1.3, label=r'pair at $g=%.1f$' % float(p['g']), **pst[k])
    ax[1].axhline(0, color='k', ls=':', lw=0.8)
    for a in ax:
        lo_, hi_ = a.get_xlim()
        a.set_xticks(list(range(int(math.ceil(lo_ - 1e-9)), int(math.floor(hi_ + 1e-9)) + 1)))
    ax[1].set_xlabel(r'support length $L$')
    ax[1].set_ylabel(r'sign$(\lambda_{\min})\,\log_{10}(|\lambda_{\min}|/10^{-9})$')
    ax[1].set_title(r'pair moved to $\frac{1}{2}\pm\frac{1}{4}\pm ig$, $w=%s$' % fw, fontsize=9)
    ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path)
    print('wrote', path)


def find_zeros(arg):
    for p in ([arg] if arg else []) + ['zeros_1700.txt', 'data/zeros_1700.txt', '../data/zeros_1700.txt', os.path.join(HERE, '..', 'data', 'zeros_1700.txt')]:
        if p and os.path.exists(p):
            return os.path.abspath(p)
    raise SystemExit('zeros_1700.txt not found')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quick', action='store_true', help='reduced case, about 30 seconds on one process; writes under /tmp by default')
    ap.add_argument('--out', default=None, help='JSON path (default data/phase_weil.json, or /tmp/phase_weil_quick.json with --quick)')
    ap.add_argument('--fig', default=None, help='figure path (default fig/cp_weil.pdf, or the --out path with .pdf with --quick)')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--zeros', default=None)
    ap.add_argument('--figure', action='store_true', help='redraw the figure from the JSON given by --out, no computation')
    a = ap.parse_args()
    out = a.out or ('/tmp/phase_weil_quick.json' if a.quick else 'data/phase_weil.json')
    fig = a.fig or (os.path.splitext(out)[0] + '.pdf' if a.quick else 'fig/cp_weil.pdf')
    if a.figure:
        draw_figure(json.load(open(out)), fig)
        return
    cfg = config(a.quick)
    cfg['zeros'] = find_zeros(a.zeros)
    res = run(cfg, a.workers)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as fh:
        json.dump(res, fh, indent=1)
    print('wrote', out, res['meta']['seconds'], res['meta']['cpu_seconds'])
    draw_figure(res, fig)


if __name__ == '__main__':
    main()
