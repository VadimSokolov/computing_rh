"""Small numerical checks for track C, using no zeta zeros.

Run: timeout 30 python3 -I tail_rates.py
The displayed numbers are ordinary high precision computations, not ball bounds.
The tail inequalities themselves are proved in report.md.
"""

from pathlib import Path
import mpmath as mp

mp.mp.dps = 60
xi0 = -mp.mpf(1) / 8 * mp.pi ** (-mp.mpf(1) / 4) * mp.gamma(mp.mpf(1) / 4) * mp.zeta(mp.mpf(1) / 2)
ratio = 16 * mp.exp(-3 * mp.pi)
c0 = 4 * mp.pi ** 2 / (xi0 * (1 - ratio))


def log_tail_bound(radius, cutoff):
    y = mp.pi * mp.exp(2 * cutoff)
    a = radius / 2 + mp.mpf(9) / 4
    if y <= a - 1:
        return mp.inf
    return mp.log(c0 / mp.pi) + (radius + mp.mpf(5) / 2) * cutoff - y - mp.log1p(-(a - 1) / y)


def cutoff_for_log_bound(radius, target):
    lo, hi = mp.mpf(0), mp.mpf(1)
    while log_tail_bound(radius, hi) > target:
        hi *= 2
    for _ in range(220):
        mid = (lo + hi) / 2
        if log_tail_bound(radius, mid) > target:
            lo = mid
        else:
            hi = mid
    return hi


lines = [
    "Numerical checks for track C. No zeta zeros used. Not ball arithmetic.",
    "precision_decimal_digits = 60",
    "xi(1/2) = " + mp.nstr(xi0, 35),
    "geometric_ratio = " + mp.nstr(ratio, 18),
    "density_bound_constant_C0 = " + mp.nstr(c0, 18),
    "Tail bound: log10 T_R(B)",
]
for cutoff in [mp.mpf("0.5"), mp.mpf(1), mp.mpf("1.5"), mp.mpf(2), mp.mpf(3)]:
    row = ["B=" + str(cutoff)]
    for radius in [mp.mpf(0), mp.mpf("0.5"), mp.mpf(1), mp.mpf(5)]:
        row.append("R=" + str(radius) + ": " + mp.nstr(log_tail_bound(radius, cutoff) / mp.log(10), 14))
    lines.append("; ".join(row))

lines.append("Cutoff B with T_(1/2)(B) <= 10^(-D)")
for digits in [6, 12, 30, 100]:
    cutoff = cutoff_for_log_bound(mp.mpf("0.5"), -digits * mp.log(10))
    lines.append("D=" + str(digits) + "; B=" + mp.nstr(cutoff, 16))

lines.append("Reference scale exp(-pi*T/4); this is not a contour lower bound.")
for height in [100, 1000, 10000]:
    target = -mp.pi * height / 4
    cutoff = cutoff_for_log_bound(mp.mpf("0.5"), target)
    digits = -target / mp.log(10)
    central_digits = -(target - cutoff / 2 - mp.log(2)) / mp.log(10)
    lines.append("T=" + str(height) + "; scale_digits=" + mp.nstr(digits, 14)
                 + "; cutoff=" + mp.nstr(cutoff, 14)
                 + "; sufficient_TV_digits_for_central_term=" + mp.nstr(central_digits, 14))

output = "\n".join(lines) + "\n"
Path(__file__).with_name("tail_numbers.txt").write_text(output)
print(output)
