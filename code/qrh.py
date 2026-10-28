"""Computations for the chapter on the quasi Riemann hypothesis proof.
(a) Poisson summation on the Eisenstein lattice Z[omega] for a Gaussian (covolume sqrt3/2).
(b) Mellin transform of the Bessel K kernel: int_0^inf K_nu(x) x^{s-1} dx = 2^{s-2} Gamma((s-nu)/2) Gamma((s+nu)/2), nu = 1/3.
(c) the principal signal for zeta: f(Z) = (1/2 pi i) int_(2) Z^{x-2/3} e^{(x-5/6)^2} / zeta(x) dx
    = Z^{-2/3}/(2 sqrt pi) sum_n mu(n) (Z/n)^{5/6} exp(-(log(Z/n))^2/4)   (Gaussian smoothed Mobius function in log scale);
    both sides, the shifted line Re x = 0.6, and the size of the first zero's residue.
(d) the 2 sin^2 splitting for the Dedekind zeta of F = Q(sqrt -3): drift threshold alpha_0^F and merged threshold alpha_1^F.
Usage: python3 qrh.py PART"""
import mpmath as mp, numpy as np, json, sys
part = sys.argv[1]; res = {}
import os
if len(sys.argv) > 2: os.chdir(sys.argv[2])
if part == "a":
    mp.mp.dps = 20
    w = mp.e**(2j*mp.pi/3)
    def S(a):  # sum over Z[omega] of exp(-pi a |lambda|^2)
        return mp.fsum(mp.e**(-mp.pi*a*abs(m + n*w)**2) for m in range(-30, 31) for n in range(-30, 31))
    # dual lattice under <x,y> = Re(x conj y): basis dual to 1, w
    B = mp.matrix([[1, mp.re(w)], [0, mp.im(w)]]); Bd = (B**-1).T
    def Sd(a):
        return mp.fsum(mp.e**(-mp.pi*(( (Bd[0,0]*m + Bd[0,1]*n)**2 + (Bd[1,0]*m + Bd[1,1]*n)**2))/a) for m in range(-30, 31) for n in range(-30, 31))
    cov = mp.im(w)
    for a in [0.5, 1, 2]:
        lhs = S(a); rhs = Sd(a)/(cov*a)
        print("Poisson", a, lhs, rhs); res[str(a)] = [float(lhs), float(rhs)]
    res["covolume"] = float(cov)
if part == "b":
    mp.mp.dps = 25; nu = mp.mpf(1)/3
    for s in [mp.mpf(2), mp.mpf(3.5), mp.mpc(2, 1.5)]:
        lhs = mp.quad(lambda x: mp.besselk(nu, x)*x**(s-1), [0, 1, 10, mp.inf])
        rhs = 2**(s-2)*mp.gamma((s-nu)/2)*mp.gamma((s+nu)/2)
        print("Mellin K", s, lhs, rhs); res[str(s)] = [str(lhs), str(rhs)]
if part == "c":
    mp.mp.dps = 20
    N = 4*10**7
    isp = np.ones(N+1, bool); isp[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if isp[p]: isp[p*p::p] = False
    mu = np.ones(N+1, np.int8)
    for p in np.nonzero(isp)[0]:
        p = int(p); mu[p::p] *= -1
        if p*p <= N: mu[p*p::p*p] = 0
    mu[0] = 0
    nz = np.nonzero(mu)[0]; mv = mu[nz].astype(float); ln = np.log(nz.astype(float))
    def f_mob(Z):
        u = np.log(Z) - ln
        return Z**(-2/3)/(2*np.sqrt(np.pi))*np.sum(mv*np.exp(5/6*u - u*u/4))
    def f_line(Z, c):
        g = lambda t: mp.re(mp.mpf(Z)**(c + 1j*t - mp.mpf(2)/3)*mp.e**((c + 1j*t - mp.mpf(5)/6)**2)/mp.zeta(c + 1j*t))
        return mp.quad(g, [-12, -4, 0, 4, 12])/(2*mp.pi)
    rows = []
    for Z in [10.0, 100.0, 1e3, 1e4, 1e5]:
        a = f_mob(Z); b = f_line(Z, 2); c = f_line(Z, mp.mpf("0.6"))
        rows.append([Z, a, float(b), float(c)]); print("signal", Z, a, float(b), float(c), flush=True)
    rho = mp.mpc(0.5, 14.134725141734693)
    amp = mp.e**((rho - mp.mpf(5)/6)**2)/mp.zeta(rho, derivative=1)
    res = {"rows": rows, "first_zero_residue_abs": str(abs(amp)), "window_N": N}
    print("first zero residue size", abs(amp))
if part == "d":
    mp.mp.dps = 25
    L3 = lambda s: mp.dirichlet(s, [0, 1, -1])
    zl = lambda s: mp.zeta(s, derivative=1)/mp.zeta(s) + mp.diff(L3, s)/L3(s)
    # completed Dedekind zeta: xi_F(s) = s(s-1) (sqrt3/2pi)^s Gamma(s) zeta_F(s);  1/s + psi(s) = psi(1+s)
    DF = lambda a: mp.psi(0, 1 + a) + mp.log(mp.sqrt(3)/(2*mp.pi)) + zl(a)
    a0 = mp.findroot(DF, 3)
    P = lambda a, y: (a - 1)/((a - 1)**2 + y**2)
    G = lambda a, y: mp.re(mp.psi(0, 1 + a + 1j*y)) - mp.psi(0, 1 + a)
    def m(a):
        f = lambda y: P(a, y) + G(a, y)
        y0 = min((mp.mpf(i)/40 for i in range(1, 600)), key=f)
        y0 = mp.findroot(lambda y: mp.diff(f, y), y0)
        return f(y0) + DF(a)
    a1 = mp.findroot(m, 2.5)
    res = {"alpha0_F": mp.nstr(a0, 12), "alpha1_F": mp.nstr(a1, 12), "D_F(2)": mp.nstr(DF(2), 8), "D_F(1.5)": mp.nstr(DF(1.5), 8)}
    print(res)
json.dump(res, open("qrh_%s.json" % part, "w"), indent=1, default=str)
if part == "e":
    # trivial zero asymptotics of the zeta principal signal, crossover with the nontrivial zeros, and dVP in Thorin form
    mp.mp.dps = 25
    c1 = mp.e**((mp.mpf(17)/6)**2)/mp.zeta(-2, derivative=1)       # residue coefficient at x = -2
    res["trivial_coeff"] = float(c1)
    f_line = lambda Z: mp.quad(lambda t: mp.re(mp.mpf(Z)**(mp.mpf(0.6) + 1j*t - mp.mpf(2)/3)*mp.e**((mp.mpf(0.6) + 1j*t - mp.mpf(5)/6)**2)/mp.zeta(mp.mpf(0.6) + 1j*t)), [-12, -4, 0, 4, 12])/(2*mp.pi)
    rows = []
    for Z in [1e3, 1e4, 1e5, 1e6]:
        v = f_line(Z); tr = c1*mp.mpf(Z)**(-mp.mpf(8)/3)
        rows.append([Z, float(v), float(tr), float(v/tr)]); print("trivial", rows[-1], flush=True)
    res["trivial_rows"] = rows
    rho = mp.mpc(0.5, 14.134725141734693)
    r1 = abs(mp.e**((rho - mp.mpf(5)/6)**2)/mp.zeta(rho, derivative=1))
    Zc = (abs(c1)/(2*r1))**(mp.mpf(2)/5)
    res["first_zero_coeff"] = str(r1); res["crossover_log10Z"] = float(mp.log10(Zc))
    print("crossover log10 Z", float(mp.log10(Zc)))
    # de la Vallee Poussin in Thorin form at alpha = 1.5: 4 Z(y) + Z(2y) <= 8 A, A = -zeta'/zeta(alpha)
    a = mp.mpf(1.5); A = -mp.zeta(a, derivative=1)/mp.zeta(a)
    Zp = lambda y: mp.re(mp.zeta(a + 1j*y, derivative=1)/mp.zeta(a + 1j*y)) + A
    ys = [mp.mpf(i)/20 for i in range(1, 2000)]
    ratio = max((4*Zp(y) + Zp(2*y))/(8*A) for y in ys); triv = max(Zp(y)/(2*A) for y in ys)
    res["dVP_ratio_max"] = float(ratio); res["one_freq_ratio_max"] = float(triv); res["A"] = float(A)
    print("dVP ratio", float(ratio), "one frequency ratio", float(triv))
    json.dump(res, open("qrh_e.json", "w"), indent=1, default=str)
