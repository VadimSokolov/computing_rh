# bsy_tailbound.py: the bounds for the Balazard, Saias and Yor defect D_theta(0) above the verified range, Chapter 35 (Balazard, Saias and Yor
# meet Thorin), Section bsy:sec:bsy (remark local versus global) and the sentence after Table bsy:tab:yks: claims 7.4e-13, 1.85e-13 and 4.6e-14 for
# D_theta(0) at theta = 1/2, 3/4, 7/8 (task C11 of reconcile-M.md). Added for the book in October 2026, written from the mathematics (the
# collaborator's scripts were not read, run or copied).
# Conventions. D_theta(0) = sum over the zeros rho = beta + i gamma with beta > theta (both signs of gamma, multiplicity) of G_theta(1, rho), where
# G_theta(s, rho) = log|(s - 2 theta + conj(rho))/(s - rho)|. At s = 1, |1 - 2 theta + conj(rho)|^2 - |1 - rho|^2 = 4 (1 - theta)(beta - theta), so
# G = (1/2) log(1 + 4 (1 - theta)(beta - theta)/((1 - beta)^2 + gamma^2)) <= 2 (1 - theta)(beta - theta)/gamma^2 < 2 (1 - theta)^2/gamma^2 (beta < 1).
# Zeros with |gamma| <= T_PT = 3 000 175 332 800 (Platt and Trudgian 2021, Theorem 1) lie on the line, so they add nothing (beta = 1/2 <= theta). The zero
# set with gamma > 0 is invariant under rho -> 1 - conj(rho) (same gamma, same multiplicity), and beta > theta >= 1/2 excludes beta' = 1 - beta > theta, so
# at most half of the zeros of ordinate gamma have beta > theta; both signs of gamma double it back. Hence D_theta(0) <= 2 (1 - theta)^2 Sigma, with
# Sigma = sum_{gamma > T_PT} gamma^(-2) (all zeros, multiplicity).
# Sigma. With nu(x) = N(x) - N(T), N(x) = M(x) + R(x), M(x) = x/(2 pi) log(x/(2 pi e)) + 7/8, |R(x)| <= E(x) (Hasanalizade, Shen and Wong 2022 as stated in
# ch/basepoint.tex:275-280: E(x) = 0.1095 log x + 0.2042 log log x + 3.0305 + 1/(50 x) for x > X0 = 30 610 046 000, and E(x) < 0.3 log x there),
# partial summation gives Sigma = 2 int_T^inf nu(x) x^(-3) dx <= B = 2 int_T^inf [M(x) - M(T) + E(x) + E(T)] x^(-3) dx. The M part is exact:
# 2 int_T^inf M(x) x^-3 dx = log(T/2 pi)/(pi T) + 7/(8 T^2); the E part is integrated by mp.quad (tanh-sinh, x = T e^y, 40 digits; checked by a second
# step count) and, as a closed form cross check, with E(x) replaced by 0.3 log x. Also logged: Sigma at T0 = 3e12 as used in ch/almost.tex (1.48e-12).
# Output: data/bsy_tailbound.json and the log bsy_tailbound.log. Usage: python3 bsy_tailbound.py [--quick] [--out FILE]; takes under a second.
import argparse, json, sys, time
from pathlib import Path
import mpmath as mp

ap = argparse.ArgumentParser()
ap.add_argument("--quick", action="store_true")
ap.add_argument("--out", default=None)
args = ap.parse_args()
here = Path(__file__).resolve().parent
out = Path(args.out) if args.out else (Path("/tmp/bsy_tailbound_quick.json") if args.quick else here.parent / "data" / "bsy_tailbound.json")
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)

def sigma_bound(T, dps, use_E="hsw"):
    mp.mp.dps = dps
    T = mp.mpf(T); L = mp.log; pi = mp.pi
    M = lambda x: x/(2*pi)*L(x/(2*pi*mp.e)) + mp.mpf(7)/8
    if use_E == "hsw":
        E = lambda x: mp.mpf("0.1095")*L(x) + mp.mpf("0.2042")*L(L(x)) + mp.mpf("3.0305") + 1/(50*x)
        # int_T^inf E(x) x^-3 dx with x = T e^y
        IE = mp.quad(lambda y: E(T*mp.e**y)*(T*mp.e**y)**-2, [0, 1, 5, 20, mp.inf])
    else:
        E = lambda x: mp.mpf("0.3")*L(x)
        IE = mp.mpf("0.3")*(L(T)/(2*T**2) + 1/(4*T**2))
    IM = L(T/(2*pi))/(2*pi*T) + mp.mpf(7)/(16*T**2)         # int_T^inf M x^-3
    B = 2*(IM - M(T)/(2*T**2) + E(T)/(2*T**2) + IE)
    return B

res = {}
Tpt = 3000175332800
for dps in (30, 50):
    B = sigma_bound(Tpt, dps)
    P("dps", dps, "Sigma bound (HSW22 E) at T_PT:", mp.nstr(B, 12))
res["Sigma_HSW"] = mp.nstr(sigma_bound(Tpt, 40), 12)
res["Sigma_HSW_dps30_minus_dps50"] = mp.nstr(sigma_bound(Tpt, 30) - sigma_bound(Tpt, 50), 3)
res["Sigma_03logx"] = mp.nstr(sigma_bound(Tpt, 40, "03"), 12)
mp.mp.dps = 40
main = (mp.log(Tpt/(2*mp.pi)) + 1)/(2*mp.pi*Tpt)
P("main term (log(T/2pi)+1)/(2 pi T):", mp.nstr(main, 12), " error part of the bound:", mp.nstr(mp.mpf(res["Sigma_HSW"]) - main, 3))
res["main_term"] = mp.nstr(main, 12)
P("Sigma with E = 0.3 log x:", res["Sigma_03logx"])
T0 = 3*10**12
res["Sigma_T0_3e12"] = mp.nstr(sigma_bound(T0, 40), 12)
P("Sigma at T0 = 3e12 (ch/almost.tex uses < 1.48e-12):", res["Sigma_T0_3e12"])
S = mp.mpf(res["Sigma_HSW"])
res["D_bounds"] = {}
for th, chap in ((mp.mpf(1)/2, "7.4e-13"), (mp.mpf(3)/4, "1.85e-13"), (mp.mpf(7)/8, "4.6e-14")):
    d = 2*(1 - th)**2*S
    d148 = 2*(1 - th)**2*mp.mpf("1.48e-12")
    P("theta", mp.nstr(th, 3), " bound 2(1-theta)^2 Sigma =", mp.nstr(d, 6), " with Sigma = 1.48e-12:", mp.nstr(d148, 6), " chapter:", chap)
    res["D_bounds"][str(th)] = {"bound": mp.nstr(d, 6), "bound_with_1.48e-12": mp.nstr(d148, 6), "chapter": chap}
# the size per zero at the first zero above T_PT, as in the remark: G at theta, beta -> 1
for th in (mp.mpf(1)/2, mp.mpf(3)/4, mp.mpf(7)/8):
    P("max per zero (beta -> 1, gamma = T_PT), both signs: ", mp.nstr(2*2*(1 - th)**2/Tpt**2, 4))
Path(out).write_text(json.dumps(res, indent=1))
Path(str(out).replace(".json", ".log")).write_text("\n".join(log) + "\n")
