# Function field model. For E: y^2 = x^3 + A x + B over F_p, a_p = p + 1 - #E(F_p). The completed zeta in x = s - 1/2 is
# L_p(x) = 2 cosh(x log p) - a_p/sqrt(p) (up to a constant). Its zeros are on the imaginary axis iff |a_p| <= 2 sqrt(p)
# (Hasse), and its Thorin density u_y(v) = (1/pi) Re L'/L(y+iv) is then a single width Cauchy mixture.
import numpy as np, json, math
def count(A, B, p):
    sq = [0]*p
    for y in range(p): sq[y*y % p] += 1
    return 1 + sum(sq[(x*x*x + A*x + B) % p] for x in range(p))
primes = [p for p in range(5, 1000) if all(p % d for d in range(2, int(p**0.5)+1))]
rows = []
for p in primes:
    if (4*1**3 + 27*1**2) % p == 0: continue
    a = p + 1 - count(1, 1, p)
    rows.append((p, a, a/(2*math.sqrt(p))))
print("curve y^2 = x^3 + x + 1: primes", len(rows), " max |a_p|/(2 sqrt p) =", max(abs(r[2]) for r in rows))
json.dump(rows, open("curve_ap.json", "w"))
# Thorin densities for a genuine Frobenius factor and for a fake one (|a| > 2 sqrt q)
def u(y, v, c, lq):  # L = 2 cosh(x lq) - 2c
    x = y + 1j*v
    return (lq*2*np.sinh(x*lq)/(2*np.cosh(x*lq) - 2*c)).real/np.pi
p = 101; lq = math.log(p)
ap = [r for r in rows if r[0] == p][0][1]; c_true = ap/(2*math.sqrt(p)); c_fake = 1.3
v = np.linspace(-np.pi/lq, np.pi/lq, 2001)
out = dict(p=p, a=ap, c_true=c_true, c_fake=c_fake, v=v.tolist(), lq=lq,
           dens={f"{kind},{y}": u(y, v, c, lq).tolist() for kind, c in [("true", c_true), ("fake", c_fake)] for y in [0.02, 0.1, 0.25]})
out["delta_fake"] = math.acosh(c_fake)/lq
json.dump(out, open("curve_thorin.json", "w"))
print("p =", p, "a_p =", ap, "c =", round(c_true, 4), " fake strip half width delta =", round(out["delta_fake"], 4))
for kind, c in [("true", c_true), ("fake", c_fake)]:
    mins = {y: float(np.min(u(y, v, c, lq))) for y in [0.02, 0.05, 0.1, 0.25]}
    print(kind, "min Thorin density by level:", {k: round(w, 4) for k, w in mins.items()})
