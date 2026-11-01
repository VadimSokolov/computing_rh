# wp3_ly.py: work package 3 (run book Section sec:spin): certify the two spin Lee and Yang property at couplings w < 1/2 on windows of height, and compare the cost with the one variable test.
# Theorem cm:thm:pair (a): the pair has the Lee and Yang property at w if no zero r of L has Re r >= w. Written for Xi(z) = xi(1/2+iz) (zero z, r = -iz...), this is: no zero of Xi with |Im z| >= w, so the window version is
#   LY(w, [X0,X1]):  Xi has no zero in the rectangle X0 <= Re z <= X1, w <= Im z <= 1.5  (Im z > 1.5 is sigma > 2, where zeta has no zero; the lower half plane is the mirror image by evenness).
# It is certified by the argument principle in ball arithmetic: the winding number of Xi around the boundary of the rectangle, each boundary segment [a,b] being accepted only when a ball enclosing Xi on the
# segment, divided by Xi(a), has real part > 0 (so the argument change is the principal argument of Xi(b)/Xi(a)); otherwise the segment is halved. The sum divided by 2 pi must be within 1e-6 of an integer: that integer is the number of zeros.
# One variable test at the same w (ferromagnetism of cosh(wx)p in the window, Theorem cm:thm:tilt): all zeros of A_w(z) = Xi(z+iw) + Xi(z-iw) in the box X0 <= Re z <= X1, |Im z| <= 1 are real: the number of
# zeros of A_w in the box (winding number, same routine) must equal the number of certified sign changes of A_w on the real grid of step 0.04. (A zero of A_w outside the box is not excluded by either count; the Lee and Yang rectangle needs no box.)
# Planted pair (nc_cauchy_*): F(z) = Xi(z)(z^2-c^2)(z^2-cbar^2)/((z^2-g_k^2)(z^2-g_{k+1}^2)), c = m + i eta, m = (g_k+g_{k+1})/2, in place of Xi. ordinates: certified acb.zeta_zeros.
# Cost is counted in evaluations of Xi (one acb zeta each) and seconds, 96 bit balls.
# Usage: python3 wp3_ly.py [--quick] [--workers 16] [--slice wp3_slice.json] [--out wp3_ly.json] [zeros_8000.txt]
import argparse, json, math, sys, time
from multiprocessing import Pool
from flint import acb, arb, ctx
ctx.prec = 96
CNT = [0]
def Xi(z):
    CNT[0] += 1
    sg = 1 if float(z.imag.mid()) >= 0 else -1
    s = acb(arb(0.5)) - acb(0, sg) * z          # Xi even: Re s = 1/2 + |Im z|
    return s * (s - 1) / 2 * (-(s / 2) * arb.pi().log()).exp() * (s / 2).gamma() * s.zeta()
def make_f(kind, planted, w):
    if planted:
        k, eta = planted
        gs = acb.zeta_zeros(k, 2); g1, g2 = gs[0].imag, gs[1].imag
        m = arb(float(((g1 + g2) / 2).mid())); c = acb(m, arb(eta)); cb = acb(m, -arb(eta))
        def base(z):
            zz = z * z
            return Xi(z) * (zz - c * c) * (zz - cb * cb) / ((zz - g1 * g1) * (zz - g2 * g2))
    else:
        base = Xi
    if kind == 'ly': return base
    iw = acb(0, arb(w))
    return lambda z: base(z + iw) + base(z - iw)
def winding(f, pts, maxlen=0.5):
    cache = {}
    def val(p):
        if p not in cache: cache[p] = f(acb(arb(p[0]), arb(p[1])))
        return cache[p]
    def seg(a, b, depth):
        fa, fb = val(a), val(b)
        cx, cy = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        box = acb(arb(cx, abs(b[0] - a[0]) / 2), arb(cy, abs(b[1] - a[1]) / 2))
        V = f(box)
        q = V / fa
        if q.real > 0:
            return float((fb / fa).arg().mid())
        if depth > 45: raise RuntimeError(f"segment {a}-{b} not resolved")
        c = (cx, cy)
        return seg(a, c, depth + 1) + seg(c, b, depth + 1)
    tot = 0.0
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        n = max(1, math.ceil(max(abs(b[0] - a[0]), abs(b[1] - a[1])) / maxlen))
        for j in range(n):
            p = (a[0] + (b[0] - a[0]) * j / n, a[1] + (b[1] - a[1]) * j / n)
            q = (a[0] + (b[0] - a[0]) * (j + 1) / n, a[1] + (b[1] - a[1]) * (j + 1) / n)
            tot += seg(p, q, 0)
    v = tot / (2 * math.pi)
    if abs(v - round(v)) > 1e-6: raise RuntimeError(f"winding {v} not an integer")
    return int(round(v))
def rect(x0, x1, y0, y1): return [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
def sign_changes(f, x0, x1, step=0.04, fine=None):
    # coarse grid of step 0.04; if fine = (center, halfwidth, step), also a fine grid there (the two real zeros born from a planted pair at w near eta are closer than 0.04)
    n = int((x1 - x0) / step); xs = [x0 + i * step for i in range(n + 1)]
    if fine:
        c, hw, st = fine
        xs += [c - hw + i * st for i in range(int(2 * hw / st) + 1)]
    xs = sorted(set(xs)); prev = None; cnt = 0; bad = 0
    for x in xs:
        v = f(acb(arb(x), arb(0))).real
        if v > 0: s = 1
        elif v < 0: s = -1
        else: bad += 1; continue
        if prev is not None and s != prev: cnt += 1
        prev = s
    return cnt, bad
def job(a):
    kind, planted, w, x0, x1 = a
    t0 = time.time(); CNT[0] = 0
    f = make_f(kind, planted, w)
    row = dict(kind=kind, planted=planted, w=w, x0=x0, x1=x1)
    try:
        if kind == 'ly':
            row['zeros_in_rectangle'] = winding(f, rect(x0, x1, w, 1.5))
        else:
            row['zeros_in_box'] = winding(f, rect(x0, x1, -1.0, 1.0))
            fine = None
            if planted:
                gs = acb.zeta_zeros(planted[0], 2); fine = (float(((gs[0].imag + gs[1].imag) / 2).mid()), 0.4, 0.0005)
            row['real_sign_changes'], row['unresolved_grid_points'] = sign_changes(f, x0, x1, fine=fine)
            row['all_real_in_box'] = row['zeros_in_box'] == row['real_sign_changes'] and row['unresolved_grid_points'] == 0
    except Exception as e:
        row['error'] = str(e)
    row['evals'] = CNT[0]; row['seconds'] = time.time() - t0
    print(json.dumps(row), flush=True)
    return row
def gapmid(G, t, side):
    import bisect
    j = bisect.bisect_right(G, t)
    return (G[j - 1] + G[j]) / 2 if side == 'lo' else (G[j] + G[j + 1]) / 2
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--slice', default=None); ap.add_argument('--planted-only', dest='po', action='store_true'); ap.add_argument('--out', default='wp3_ly.json'); ap.add_argument('zeros', nargs='?', default='zeros_8000.txt')
    a = ap.parse_args()
    G = [float(l.split()[-1]) for l in open(a.zeros) if l.strip()]
    wins = [(10.0, 110.0), (1000.0, 1100.0)] if not a.quick else [(10.0, 40.0)]
    # window ends at midpoints of gaps
    wins = [(gapmid(G, lo, 'lo') if lo > G[0] else lo, gapmid(G, hi, 'lo')) for lo, hi in wins]
    if not a.quick:   # the window of ordinates 9900 to 10000 (heights about 9778 to 9878), ordinates from the certified acb.zeta_zeros
        ctx.prec = 96; zz = [float(z.imag.mid()) for z in acb.zeta_zeros(9900, 102)]
        wins.append(((zz[0] + zz[1]) / 2, (zz[100] + zz[101]) / 2))
    ws = [0.5, 0.4, 0.3, 0.2, 0.1, 0.05] if not a.quick else [0.4, 0.1]
    tasks = []
    for x0, x1 in ([] if a.po else wins):
        for w in ws:
            tasks.append(('ly', None, w, x0, x1)); tasks.append(('onevar', None, w, x0, x1))
    planted = []
    wst = {}
    if a.slice:
        for r in json.load(open(a.slice)):
            if r['place'] == 'midpoint': wst[(r['k'], r['eta'])] = r['w_refined']
    if not a.quick:
        for k, eta in [(8, 0.3), (300, 0.3), (1000, 0.1), (9999, 0.1)]:
            if k < len(G) - 5: x0, x1 = gapmid(G, G[k - 1] - 8, 'lo'), gapmid(G, G[k] + 8, 'lo')
            else:
                ctx.prec = 96; zz = [float(z.imag.mid()) for z in acb.zeta_zeros(k - 8, 20)]; x0, x1 = (zz[0] + zz[1]) / 2, (zz[16] + zz[17]) / 2
            wmid = (wst[(k, eta)] + eta) / 2 if (k, eta) in wst else eta * 0.99
            for w in sorted({eta * 0.9, wst.get((k, eta), eta * .98) - 1e-4 * eta, wmid, eta * 1.05}):
                planted.append(((k, eta), w, x0, x1))
        for (kk, w, x0, x1) in planted:
            tasks.append(('ly', kk, w, x0, x1)); tasks.append(('onevar', kk, w, x0, x1))
    with Pool(a.workers) as p: rows = p.map(job, tasks, chunksize=1)
    json.dump(dict(wstar=[dict(k=k, eta=e, w_star=v) for (k, e), v in wst.items()], rows=rows), open(a.out, 'w'), indent=0)
