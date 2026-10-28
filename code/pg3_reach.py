# pg3_reach.py: the reach of the Cauchy circle for the heat trace W (Remark of Section ch:heat added in the merge of note X3, Chapter 12,
# end of the paragraph on the floor; placeholder N1). Added for the book in October 2026, written from the mathematics (the collaborator's
# scripts were not read, run or copied).
# Claim. The circle applied to W bounds the prime terms of order k by a quantity below the archimedean term only for t <~ c_P/k, with
# c_P = (log 2)^2/6, while the zeros above H stay below the first zero only for t >~ 2 k log(H/gamma_1)/H^2. The two ranges meet at
# k = K(H) = H (c_P/(2 log(H/gamma_1)))^(1/2); the script prints K(H) for H = 3e12 and the ratio K/H (heuristic, not certified).
# Method. gamma_1 = 14.1347251417... from mpmath.zetazero(1) at 30 digits and, as a check, from the first line of data/zeros_8000.txt;
# K is computed at 30 and at 50 digits and the relative difference is logged. Options: --quick (same, nothing to shorten), --out FILE.
import sys, json, argparse
import mpmath as mp

def K(H, dps):
    mp.mp.dps = dps
    g1 = mp.zetazero(1).imag
    cP = mp.log(2) ** 2 / 6
    H = mp.mpf(H)
    return g1, cP, H * mp.sqrt(cP / (2 * mp.log(H / g1)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--out', default='pg3_reach.json')
    a = ap.parse_args()
    L = []
    def log(s):
        print(s); L.append(s)
    H = 3 * 10 ** 12
    g30, c30, K30 = K(H, 30)
    g50, c50, K50 = K(H, 50)
    rel = abs(K30 - K50) / K50
    log('gamma_1 (30 digits) = %s' % mp.nstr(g30, 25))
    log('c_P = (log 2)^2/6 = %s' % mp.nstr(c30, 20))
    log('H = 3e12: K(H) = %s, K/H = %s (30 digits); relative difference of the 30 and 50 digit runs %s'
        % (mp.nstr(K30, 12), mp.nstr(K30 / H, 8), mp.nstr(rel, 3)))
    log('1.29 (H-23)^2 = %s, the order of Theorem thm:ch12:window exceeds K(H) by the factor %s'
        % (mp.nstr(mp.mpf('1.29') * (H - 23) ** 2, 6), mp.nstr(mp.mpf('1.29') * (H - 23) ** 2 / K30, 6)))
    res = dict(H=H, gamma_1=mp.nstr(g30, 30), cP=mp.nstr(c30, 30), K=mp.nstr(K30, 20), K_over_H=mp.nstr(K30 / H, 15),
               N1=mp.nstr(K30, 5), rel_diff_precision=mp.nstr(rel, 3), certified=False)
    json.dump(res, open(a.out, 'w'), indent=1)
    open(a.out.replace('.json', '.log'), 'w').write('\n'.join(L) + '\n')
