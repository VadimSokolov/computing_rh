# tlz_k.py: TLZ.1 of reconcile-TLZ.md. Zeros of zeta(z) - zeta(z+2) (the sum of four copies, k = 4) and of zeta(z) - 5 zeta(z+2) + 4 zeta(z+4) (k = 6),
# fm:res:kzeros (Proposition fm:prop:kspins: D_{Y_4}(w) = (2pi)^(-2w) Gamma(2w+4)/3 [zeta(2w+1) - zeta(2w+3)], D_{Y_6} likewise, z = 2w+1).
# Method: Newton at 50 digits (mpmath) from the plan's starting points; then the argument principle in ball arithmetic (python-flint acb zeta, 128 bits,
# tlz_lib.Counter) on the square of half side 2^-27 about the root rounded to the grid 2^-60: the winding number is 1, so the square holds exactly one zero
# (computed in ball arithmetic: a proof, as the balls exclude 0 on the contour). Also a second run at 192 bits. 10 digits printed.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import sys, json, argparse, time
import mpmath as mp
from flint import acb, arb, ctx
import tlz_lib

FUN = {
 'k4': (lambda z: z.zeta() - (z + 2).zeta(), lambda z: mp.zeta(z) - mp.zeta(z + 2)),
 'k6': (lambda z: z.zeta() - 5 * (z + 2).zeta() + 4 * (z + 4).zeta(), lambda z: mp.zeta(z) - 5 * mp.zeta(z + 2) + 4 * mp.zeta(z + 4)),
}
START = {'k4': [(1.76891972665, 23.3078297413), (0.650238826497, 31.7124838244)],
         'k6': [(3.18634599215, 31.979112615), (1.51363208762, 41.215007487)]}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='tlz_k.json'); ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    mp.mp.dps = 50
    res = {}
    for key in ('k4', 'k6'):
        fa, fm = FUN[key]
        res[key] = []
        for (x, y) in START[key]:
            z = mp.findroot(fm, mp.mpc(x, y), tol=1e-45, maxsteps=100)
            row = {'re': mp.nstr(z.real, 20), 'im': mp.nstr(z.imag, 20), 'resid': mp.nstr(abs(fm(z)), 3)}
            cr, ci = tlz_lib.grid(mp.nstr(z.real, 25)), tlz_lib.grid(mp.nstr(z.imag, 25))
            r = arb(2) ** -27
            for prec in ((128,) if a.quick else (128, 192)):
                ctx.prec = prec
                C = tlz_lib.Counter(fa)
                n, w = C.winding(cr - r, cr + r, ci - r, ci + r)
                row['winding_%d' % prec] = n
                row['minabs_lower_%d' % prec] = str(C.minabs.str(5)) if C.minabs is not None else None
                row['evals_%d' % prec] = C.evals
                # the root lies in the square
                row['in_square'] = bool(abs(float(z.real) - float(cr.mid())) < float(r.mid()) and abs(float(z.imag) - float(ci.mid())) < float(r.mid()))
            ctx.prec = 128
            row['print10'] = '%s + %si' % (mp.nstr(z.real, 12), mp.nstr(z.imag, 12))
            row['plan_start'] = [x, y]
            row['diff_from_plan'] = [float(z.real) - x, float(z.imag) - y]
            print(key, json.dumps(row), flush=True)
            res[key].append(row)
    json.dump(res, open(a.out, 'w'), indent=1)

if __name__ == '__main__':
    main()
