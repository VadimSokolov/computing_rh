# tlz_onecos.py: TLZ.2 of reconcile-TLZ.md. The pair of two copies of the law with density (1+cos t)/(2 pi) on [-pi, pi] (fm:res:onecos):
# D(w) = D0(w) + 2 D0(w+1) + D0(w+2), D0(w) = -2 pi (2w+1) zeta(-2w) / sin(pi w) (fm:eq:Dres, continued meromorphically).
# Counts of zeros in the boxes -1 < Re w < 1 and -3 < Re w < -2, 1/2 < Im w < 30 by the argument principle in ball arithmetic (python-flint acb zeta and sin,
# tlz_lib.Counter, 128 bits and 192 bits): there are no poles off the real axis, so the winding number is the number of zeros (computed in ball arithmetic).
# Control: the same count for D0 in the first box (the 13 zeros -1/4 + i gamma_k/2, gamma_k < 60). The zeros in the second box are then located: Newton at 50 digits
# (mpmath) from the starting points of the plan, each isolated by a winding number 1 on a square of half side 2^-27. Checks of the closed form: D at -0.3+4i.
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import json, argparse
import mpmath as mp
from flint import acb, arb, ctx
import tlz_lib

def D0a(w):
    return -2 * arb.pi() * (2 * w + 1) * (-2 * w).zeta() / (arb.pi() * w).sin()
def Fa(w):
    return D0a(w) + 2 * D0a(w + 1) + D0a(w + 2)
def D0m(w):
    return -2 * mp.pi * (2 * w + 1) * mp.zeta(-2 * w) / mp.sin(mp.pi * w)
def Fm(w):
    return D0m(w) + 2 * D0m(w + 1) + D0m(w + 2)

START = [(-2.444, 15.838), (-2.977, 19.239), (-2.479, 21.119), (-2.273, 24.427), (-2.862, 26.207), (-2.757, 28.554)]

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', default='tlz_onecos.json'); ap.add_argument('--quick', action='store_true')
    a = ap.parse_args()
    out = {}
    mp.mp.dps = 50
    out['D0_at_-0.3+4i'] = str(mp.nstr(D0m(mp.mpc(-0.3, 4)), 12))
    boxes = {'control_D0_box1': (D0a, -1, 1), 'TLZ.2a_box1': (Fa, -1, 1), 'TLZ.2b_box2': (Fa, -3, -2)}
    for name, (F, x0, x1) in boxes.items():
        for prec in ((128,) if a.quick else (128, 192)):
            ctx.prec = prec
            C = tlz_lib.Counter(F)
            n, w = C.winding(arb(x0), arb(x1), arb('0.5'), arb(30))
            out['%s_%d' % (name, prec)] = {'count': n, 'winding_ball': w.str(8), 'minabs_lower': C.minabs.str(5), 'evals': C.evals, 'maxdepth': C.maxdepth_used}
            print(name, prec, json.dumps(out['%s_%d' % (name, prec)]), flush=True)
    ctx.prec = 128
    zs = []
    for (x, y) in START:
        z = mp.findroot(Fm, mp.mpc(x, y), tol=1e-45, maxsteps=100)
        cr, ci = tlz_lib.grid(mp.nstr(z.real, 25)), tlz_lib.grid(mp.nstr(z.imag, 25)); r = arb(2) ** -27
        C = tlz_lib.Counter(Fa); n, w = C.winding(cr - r, cr + r, ci - r, ci + r)
        row = {'re': mp.nstr(z.real, 20), 'im': mp.nstr(z.imag, 20), 'print6': '%s + %si' % (mp.nstr(z.real, 7), mp.nstr(z.imag, 7)), 'winding': n, 'resid': mp.nstr(abs(Fm(z)), 3)}
        print('zero', json.dumps(row), flush=True); zs.append(row)
    out['zeros_box2'] = zs
    # distinct and inside the box
    out['distinct'] = len(set((r_['re'], r_['im']) for r_ in zs))
    out['inside_box2'] = all(-3 < float(r_['re']) < -2 and 0.5 < float(r_['im']) < 30 for r_ in zs)
    json.dump(out, open(a.out, 'w'), indent=1)

if __name__ == '__main__':
    main()
