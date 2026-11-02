#!/usr/bin/env python3
"""D2 proposal P3, NOT RUN (Hopper only; syntax checked, never executed).

Exact sign of Weil's even form at X with a planted quadruple 1/2 +- delta +- i gamma added to the true zeros, gamma the midpoint
of the ordinates k and k+1. The matrix comes from the prime side (code/wl_weil.py build), so no ordinate enters except gamma itself.
With A > 0 the form on the true zeros and Q = 4 (v_r v_r^T - v_i v_i^T), v = (F_n(z)), z = gamma - i delta, inertia additivity gives:
A + Q has a negative eigenvalue if and only if the 2 x 2 matrix S = [[1/4 + M11, M12], [M12, -1/4 + M22]] is positive definite,
where M = V^T A^{-1} V, V = [v_r, v_i]. Since S11 > 0 this is det S > 0 together with S22 > 0.
Ball mode (--ball) gives a certificate of either sign when the balls of S exclude 0; the default strips radii in the factorisation.

Place beside code/wl_weil.py and data/zeros_8000.txt:
  python3 D2_hopper_P3.py 164 1517 0.01 --digits 1300 --out p3_164.json
  python3 D2_hopper_P3.py 181 1517 0.01 --digits 1300 --out p3_181.json
  python3 D2_hopper_P3.py 200 1517 0.01 --digits 1300 --out p3_200.json
Size: N about 1700 modes (N = ceil(L (gamma+45)/pi) + 2), about 5 core hours each at 1300 digits and 3 GB (factor 3 uncertain);
ball mode at digits = 3.5 N + 1920 about 4 core days each.
"""
import argparse, json, math, time
from flint import arb, acb, ctx
import wl_weil as W


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('X', type=int)
    ap.add_argument('k', type=int)
    ap.add_argument('delta')
    ap.add_argument('--digits', type=int, default=1300)
    ap.add_argument('--ball', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    t0 = time.time()
    lines = [s.strip() for s in open('zeros_8000.txt') if s.strip()]
    L0 = 0.5 * math.log(a.X)
    gm0 = 0.5 * (float(lines[a.k - 1]) + float(lines[a.k]))
    N = math.ceil(L0 * (gm0 + 45) / math.pi) + 2
    D = int(3.5 * N + 1920) if a.ball else a.digits
    A, g, info = W.build(a.X, N, D)
    ctx.prec = W.bits_for(D)
    gm = (arb(lines[a.k - 1]) + arb(lines[a.k])) / 2
    z = acb(gm, -arb(a.delta))
    L = g.L
    v = [g.c[n] * (((z - g.a[n]) * L).sin() / (z - g.a[n]) + ((z + g.a[n]) * L).sin() / (z + g.a[n])) for n in range(N + 1)]
    vr = [x.real for x in v]
    vi = [x.imag for x in v]
    f = W.ldl(A, a.ball)
    if f is None:
        json.dump(dict(X=a.X, N=N, D=D, status='A not positive in this precision'), open(a.out, 'w'))
        return
    Dd, Lr = f
    yr = W.ldl_solve(Dd, Lr, vr)
    yi = W.ldl_solve(Dd, Lr, vi)
    M11 = sum((p * q for p, q in zip(vr, yr)), arb(0))
    M12 = sum((p * q for p, q in zip(vr, yi)), arb(0))
    M22 = sum((p * q for p, q in zip(vi, yi)), arb(0))
    q = arb(1) / 4
    S11, S22 = q + M11, M22 - q
    det = S11 * S22 - M12 * M12
    out = dict(X=a.X, N=N, D=D, k=a.k, delta=a.delta, ball=a.ball, S22=S22.str(20), det=det.str(20),
               negative_eigenvalue=bool(S22 > 0 and det > 0), positive_definite=bool(det < 0 or S22 < 0), seconds=round(time.time() - t0, 1))
    json.dump(out, open(a.out, 'w'), indent=0)
    print(out)


if __name__ == '__main__':
    main()
