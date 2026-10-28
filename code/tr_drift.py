# tr_drift.py: the drift law of Chapter 28, Section tr: "The drift hierarchy" (thm:drift and its table), the values M(0.3), M_40^(1)(0.3),
# M_200^(0)(0.3), the local accuracy of the mean corrected first zero, and the test of the finding M2 of the review: the coefficient of
# n^-(2m+1) is delta_m only when kappa_{m+1}(R) = o(n^-(2m+1)); for m = 2 with R a gamma variable matching the mean and the variance of T_n
# the claim is (1 - kappa_3(R)/kappa_3(T_n)) delta_2 = (4/9) delta_2.
# Conventions. omega = sum_k g_k/lambda_k, lambda_k = 2 pi^2 (k - 1/2)^2, omega_n = first n atoms, T_n = tail, kappa_j(T_n) =
# (j-1)! (2 pi^2)^-j zeta(2j, n + 1/2) (Hurwitz). M_n^R(z) = E (omega_n + R)^(-z), M(z) = (2/sqrt(pi)) (2 pi)^z Lambda(2z). R = 0 (m = 0),
# R = c_n = kappa_1(T_n) (m = 1), R ~ Gamma(shape k, scale th) with k = kappa_1^2/kappa_2, th = kappa_2/kappa_1 (m = 2, "gamma2").
# Method (no A_k cancellation). E X^(-z) = (1/Gamma(z)) int_0^inf s^(z-1) E exp(-s X) ds, E exp(-s omega_n) = prod_{k<=n} (1 + s/lambda_k)^-1,
# E exp(-s R) = 1, exp(-s c_n), (1 + th s)^(-k). With s = e^u the integrand e^(zu) F(e^u) is analytic in |Im u| < pi (poles at s = -lambda_k,
# -1/th), so the trapezoid rule of step h has error about exp(-2 pi^2/h) (h = 1/4: 1e-34; h = 1/8: 1e-68). F is tabulated once per (m, n).
# Zero z_n near z_0 = rho_1/2 found by the secant method; delta^(n) = n^(2m+1) (z_n - z_0). Every number is computed twice: (dps 40, h = 1/4)
# and (dps 50, h = 1/8); the difference of the zeros is logged and written. Prediction: thm:drift at the same zero; for m = 2 also the
# corrected prediction with the exact kappa_3(R)/kappa_3(T_n) at each n.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
# Output: data/tr_drift.json, log data/tr_drift.log. Usage: python3 tr_drift.py [--quick] [--out FILE]; --quick: n <= 32, seconds. Full run
# on Hopper (about 5 minutes): bash misc/tools/hopper_run.sh -c 1 -t 60 -g tr_drift.json -g tr_drift.log code/tr_drift.py -- --out tr_drift.json
import argparse, json, time
import mpmath as mp

def L(s): return mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4) - mp.zeta(s, mp.mpf(3) / 4))
def Lam(s): return (4 / mp.pi) ** (s / 2) * mp.gamma((s + 1) / 2) * L(s)
def Mtrue(z): return 2 / mp.sqrt(mp.pi) * (2 * mp.pi) ** z * Lam(2 * z)
def lam(k): return 2 * mp.pi ** 2 * (k - mp.mpf(1) / 2) ** 2
def kap(j, n): return mp.factorial(j - 1) * (2 * mp.pi ** 2) ** (-j) * mp.zeta(2 * j, n + mp.mpf(1) / 2)

class Trunc:
    def __init__(self, kind, n, dps, h, zre):
        mp.mp.dps = dps
        self.h = mp.mpf(h); self.kind = kind; self.n = n
        lams = [lam(k) for k in range(1, n + 1)]
        if kind == 'mean': c = kap(1, n)
        if kind == 'gamma2':
            shape = kap(1, n) ** 2 / kap(2, n); th = kap(2, n) / kap(1, n)
        tol = mp.mpf(10) ** (-dps - 4)
        def F(u):
            s = mp.exp(u); p = mp.mpf(1)
            for l in lams: p *= l / (l + s)
            if kind == 'mean': p *= mp.exp(-c * s)
            if kind == 'gamma2': p *= (1 + th * s) ** (-shape)
            return p
        self.us = []; self.Fs = []
        j = 0
        while True:   # downwards: integrand ~ e^(zre u)
            u = -j * self.h; f = F(u)
            self.us.append(u); self.Fs.append(f)
            if mp.exp(zre * u) * f < tol: break
            j += 1
        j = 1
        while True:
            u = j * self.h; f = F(u)
            self.us.append(u); self.Fs.append(f)
            if mp.exp(zre * u) * f < tol: break
            j += 1
    def M(self, z):
        return self.h / mp.gamma(z) * mp.fsum(mp.exp(z * u) * f for u, f in zip(self.us, self.Fs))

def solve(kind, n, z0, dps, h):
    T = Trunc(kind, n, dps, h, mp.re(z0))
    z = mp.findroot(T.M, (z0, z0 * (1 + mp.mpf(10) ** -6)), solver='secant', tol=mp.mpf(10) ** (-2 * dps + 10), maxsteps=60)
    return z

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--quick', action='store_true'); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or ('/tmp/tr_drift.json' if a.quick else 'data/tr_drift.json')
    logf = open(out.replace('.json', '.log'), 'w')
    def log(s): print(s); logf.write(s + '\n'); logf.flush()
    t0 = time.time(); res = {}
    mp.mp.dps = 50
    g1 = mp.findroot(lambda x: mp.re(mp.exp(1j * (x / 2 * mp.log(4 / mp.pi) + mp.im(mp.loggamma((mp.mpf(3) / 2 + 1j * x) / 2)))) * L(mp.mpf(1) / 2 + 1j * x)), 6.02)
    rho = mp.mpf(1) / 2 + 1j * g1; z0 = rho / 2
    Lp = lambda s: -mp.log(4) * L(s) + mp.mpf(4) ** (-s) * (mp.zeta(s, mp.mpf(1) / 4, 1) - mp.zeta(s, mp.mpf(3) / 4, 1))
    lp = (4 / mp.pi) ** (rho / 2) * mp.gamma((rho + 1) / 2) * Lp(rho)
    def delta(m): return (-1) ** (m + 1) * mp.rf(z0, m + 1) * Lam(rho + 2 * m + 2) / (2 * (m + 1) * (2 * m + 1) * mp.pi ** (m + 1) * lp)
    pred = {m: delta(m) for m in (0, 1, 2)}
    log(f"gamma_1 = {mp.nstr(g1, 15)}; predicted delta_0 = {mp.nstr(pred[0], 8)}, delta_1 = {mp.nstr(pred[1], 8)}, delta_2 = {mp.nstr(pred[2], 8)}")
    # values at z = 0.3
    mp.mp.dps = 40
    Mt = Mtrue(mp.mpf('0.3')); log(f"M(0.3) = {mp.nstr(Mt, 10)}")
    for kind, n in ([('mean', 40), ('hard', 200)] if not a.quick else [('mean', 10), ('hard', 20)]):
        v1 = Trunc(kind, n, 40, 0.25, mp.mpf('0.3')).M(mp.mpf('0.3')); v2 = Trunc(kind, n, 50, 0.125, mp.mpf('0.3')).M(mp.mpf('0.3'))
        log(f"M_{n}^({kind})(0.3) = {mp.nstr(v1, 10)}; second run differs by {mp.nstr(abs(v1 - v2), 3)}")
        res[f'val_{kind}_{n}'] = str(v1)
    cases = {'hard': [8, 16, 32, 64, 128, 256], 'mean': [4, 8, 16, 32, 64, 128], 'gamma2': [16, 32, 64, 128]}
    if a.quick: cases = {'hard': [8, 16], 'mean': [4, 8], 'gamma2': [8, 16]}
    mm = {'hard': 0, 'mean': 1, 'gamma2': 2}
    for kind, ns in cases.items():
        m = mm[kind]
        for n in ns:
            zA = solve(kind, n, z0, 40, 0.25); zB = solve(kind, n, z0, 50, 0.125)
            mp.mp.dps = 40
            d = n ** (2 * m + 1) * (zA - z0)
            line = f"{kind} n={n}: n^{2*m+1}(z_n - z0) = {mp.nstr(d, 8)}; Re z_n - 1/4 = {mp.nstr(mp.re(zA) - mp.mpf(1) / 4, 6)}; |dz| between runs = {mp.nstr(abs(zA - zB), 2)}; delta_{m} = {mp.nstr(pred[m], 8)}; ratio d/delta = {mp.nstr(d / pred[m], 6)}"
            if kind == 'gamma2':
                shape = kap(1, n) ** 2 / kap(2, n); k3R = 2 * kap(2, n) ** 2 / kap(1, n)
                fac = 1 - k3R / kap(3, n)
                line += f"; 1 - k3(R)/k3(T_n) = {mp.nstr(fac, 6)}; d/(fac*delta_2) = {mp.nstr(d / (fac * pred[2]), 6)}"
            log(line)
            res[f'{kind}_{n}'] = {'d': [str(d.real), str(d.imag)], 'rez_minus_quarter': str(mp.re(zA) - mp.mpf(1) / 4), 'run_diff': str(abs(zA - zB))}
    res['pred'] = {str(m): [str(v.real), str(v.imag)] for m, v in pred.items()}; res['gamma1'] = str(g1)
    mp.mp.dps = 20
    log(f"hard cutoff: n for Re z_n - 1/4 = 8e-7 from Re delta_0/n: {mp.nstr(mp.re(pred[0]) / mp.mpf('8e-7'), 3)}")
    json.dump(res, open(out, 'w')); log(f"time {time.time() - t0:.1f} s")
if __name__ == '__main__':
    main()
