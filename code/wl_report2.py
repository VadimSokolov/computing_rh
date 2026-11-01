# wl_report2.py: writes misc/reviews/nick-runs/checks-WL2.md from data/wl_*.json (the second part of Work package 1: prime dependence, Davenport and Heilbronn, accuracy law, odd modes).
# Reads json only. Run: timeout 120 python3 -I code/wl_report2.py
import json, os, math
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/'
D = R + 'data/'


def J(f):
    p = D + f
    return json.load(open(p)) if os.path.exists(p) else None


def sg(c):
    return 'POS' if c.get('certified_positive') else 'NEG' if c.get('certified_negative') else '??'


def fe(x, n=3):
    try:
        return ('%.' + str(n) + 'e') % float(x)
    except Exception:
        return str(x)


out = []
w = out.append
w('# Checks WL2 (Work package 1, second part: prime dependence, Davenport and Heilbronn, accuracy law, odd modes)\n')
w('Written from the mathematics (code/wl_dh.py header: the explicit formula of the Davenport and Heilbronn function, the odd modes; code/wl_primes.py; code/wl_accuracy.py; this file by code/wl_report2.py). Nick\'s files in misc/reviews were not read, run or copied. Nothing was committed. Every number below is copied from the json file named in its section. Status words as in checks-WL.md: POS = certified positive (ball LDL^T, every pivot ball positive), NEG = certified negative (ball Rayleigh quotient of the float eigenvector, its upper bound below 0, so the true matrix has a negative eigenvalue), ?? = sign not certified.\n')

# ---------------------------------------------------------------- validation
V = J('wl_dh_validate.json')
Z = J('wl_dh_zeros.json')
w('## 1. Files and runs\n')
w('| file | what |\n|---|---|')
for f, t in (('code/wl_dh.py', 'general builder (zeta/DH, even/odd), zeros of f, validation, cell routine'), ('code/wl_primes.py', 'truncated prime sum grid'), ('code/wl_accuracy.py', 'accuracy grid and fit'),
             ('code/wl_report2.py', 'this note'), ('data/wl_dh_zeros.json', 'zeros of f (line zeros by sign changes, off line zeros by a complex secant search)'),
             ('data/wl_dh_validate.json', 'explicit formula against the zeros'), ('data/wl_dh_scan1.json', 'DH scan: signs and eigenvalues'), ('data/wl_dh_scan2.json', 'DH failing cells: zeros of xi_hat, decomposition over the zeros'),
             ('data/wl_primes_scan.json', 'truncated primes'), ('data/wl_acc.json, data/wl_accuracy_fit.json', 'accuracy grid and fits'), ('data/wl_odd_scan.json', 'odd modes, zeta')):
    w('| %s | %s |' % (f, t))
w('')
w('Runs on Hopper (AMD nodes, python-flint 0.9, mpmath 1.3), run folders wl_dh_z4, wl_dh_val, wl_dh_scan1, wl_dh_scan2, wl_primes3, wl_acc4, wl_odd3 under /scratch/vsokolov/rh_book_repro/agents (the job wl_grid3 was not touched).\n')
if Z:
    w('## 2. Zeros of the Davenport and Heilbronn function (wl_dh_zeros.json)\n')
    w('kappa = %s; functional equation Lambda_f(s) = Lambda_f(1 - s): max relative residual %s at three points; |Im Lambda_f|/|Re Lambda_f| on the line at three heights: %s.' % (Z['kappa'], Z['checks']['functional_equation_max_rel'], Z['checks']['Lambda_imag_on_line_max_rel']))
    w('Zeros up to height T = %d: %d on the line (sign changes of Lambda_f(1/2 + it) on a grid of step %s), %d off the line in the upper half plane (each with its three images). Computation %s s.\n' % (Z['T'], Z['n_online'], Z['step'], Z['n_offline'], Z['seconds']))
    w('| Re rho | Im rho |\n|---|---|')
    for a, b in Z['offline']:
        w('| %s | %s |' % (a[:22], b[:22]))
    if 'cf_check' in Z:
        w('\nAgainst the certified count of data/cf_dhzeros.json (T2: %s sign changes on the line below 200, four off line zeros): this run found %d line zeros and %d off line zeros; the first 20 digits of the off line zeros agree with the certified ones (cf: %s).\n' % (Z['cf_check']['cf_sign_changes'], Z['n_online'], Z['n_offline'], '; '.join('%s+%si' % tuple(t) for t in Z['cf_check']['cf_offline'])))
if V:
    w('## 3. Validation (wl_dh_validate.json)\n')
    w('Test functions in the span of the first modes with high order contact at +-L: f = cos^4(pi u/2L) (even) and f = sin x cos^3 x, x = pi u/2L (odd). W(f) from the ball matrix against sum over the zeros of ghat(gamma) = F(gamma)F(-gamma) (line zeros: 2 ghat(t); off line zeros: the four images +-t +-i delta), T = 200 for DH, the 1700 ordinates of data/zeros_1700.txt for zeta.\n')
    w('| test | kind | parity | X | N | W from the matrix | zero sum | rel diff | rel diff without the off line zeros |\n|---|---|---|---|---|---|---|---|---|')
    for r in V:
        if r['test'] == 'explicit_formula':
            w('| formula | %s | %s | %d | %d | %s | %s | %s | %s |' % (r['kind'], r['parity'], r['X'], r['N'], r['W_matrix'], r['zero_sum_all'], fe(r['rel_diff']), fe(r['rel_diff_without_offline']) if 'rel_diff_without_offline' in r else '--'))
        else:
            w('| builder vs wl_weil.build | zeta | even | %d | %d | entries overlap: %s | max mid difference 10^%s | | |' % (r['X'], r['N'], r['all_overlap'], r['max_mid_diff_log10']))
    w('')

# ---------------------------------------------------------------- DH scan
S1 = J('wl_dh_scan1.json')
S2 = J('wl_dh_scan2.json')
if S1:
    w('## 4. The Davenport and Heilbronn ladder (wl_dh_scan1.json%s)\n' % (', wl_dh_scan2.json' if S2 else ''))
    ev = {(c['X'], c['N']): c for c in S1 + (S2 or []) if c['parity'] == 'even'}
    Xs = sorted(set(k[0] for k in ev))
    Ns = sorted(set(k[1] for k in ev))
    w('Least eigenvalue of the DH form on the N+1 even modes (mpmath eigsy at 3.5N + 60 digits and above), with the certified sign.\n')
    w('| X \\ N | ' + ' | '.join(str(n) for n in Ns) + ' |\n|---|' + '---|' * len(Ns))
    for X in Xs:
        w('| %d | ' % X + ' | '.join(('%s %s' % (sg(ev[(X, n)]), fe(ev[(X, n)]['eps'], 2)) if (X, n) in ev else '') for n in Ns) + ' |')
    w('')
    w('First failing N at each X, and first failing X at each N (certified negative, N and X on this grid): ')
    fx = []
    for X in Xs:
        bad = [n for n in Ns if (X, n) in ev and ev[(X, n)].get('certified_negative')]
        fx.append('X=%d: %s' % (X, min(bad) if bad else 'none'))
    w('; '.join(fx) + '.')
    fn = []
    for n in Ns:
        bad = [X for X in Xs if (X, n) in ev and ev[(X, n)].get('certified_negative')]
        fn.append('N=%d: %s' % (n, min(bad) if bad else 'none'))
    w('; '.join(fn) + '.\n')
    odd = [c for c in S1 if c['parity'] == 'odd']
    if odd:
        w('DH odd modes: ' + '; '.join('X=%d N=%d %s %s' % (c['X'], c['N'], sg(c), fe(c['eps'], 3)) for c in sorted(odd, key=lambda c: (c['X'], c['N']))) + '.\n')
if S2:
    w('### 4b. The failing cells: zeros of xi_hat and the decomposition over the zeros of f (wl_dh_scan2.json)\n')
    for c in sorted(S2, key=lambda c: (c['X'], c['N'])):
        line = '- X=%d N=%d %s eps=%s, %s negative eigenvalues, %s zeros of xi_hat below pi N/L = %s' % (c['X'], c['N'], sg(c), fe(c['eps'], 3), c.get('n_negative'), c.get('nroots'), fe(c.get('top', float('nan')), 3))
        if c.get('err_log10'):
            line += '; distances to the first DH line zeros (log10): ' + ', '.join('%.1f' % e if e is not None else '--' for e in c['err_log10'][:6]) + '; matched to 1e-3: %s of %s' % (c.get('matched_1e_3'), c.get('n_ordinates'))
        w(line)
        if c.get('decomp'):
            d = c['decomp']
            w('  - eigenvector norm 1: W = eigenvalue %s; line zeros contribute %s; off line quadruples: ' % (fe(c['eps'], 4), d['line_sum']) + '; '.join('%s+%si: %s' % (o['re'], o['im'], o['quad']) for o in d['offline']))
    w('')

# ---------------------------------------------------------------- primes
P = J('wl_primes_scan.json')
if P:
    w('## 5. Truncated prime sum for zeta (wl_primes_scan.json)\n')
    w('Weil form of zeta on the N+1 even modes of [1/lambda, lambda], X = lambda^2, with the prime sum cut at Y <= X (all prime powers n <= Y; Y = 1: no primes). eps = least eigenvalue; sign certified as above.\n')
    from collections import defaultdict
    g = defaultdict(list)
    for c in P:
        g[(c['X'], c['N'])].append(c)
    for (X, N) in sorted(g):
        cs = sorted(g[(X, N)], key=lambda c: c['Y'])
        full = [c for c in cs if c['Y'] == X][0]
        pos = [c for c in cs if c['certified_positive']]
        neg = [c for c in cs if c['certified_negative']]
        unc = [c['Y'] for c in cs if not c['certified_positive'] and not c['certified_negative']]
        w('### X = %d, N = %d (%d cells)\n' % (X, N, len(cs)))
        w('| Y | ' + ' | '.join(str(c['Y']) for c in cs if c['Y'] in (1, 2, 3, 5, 7, 11, 13, 17, 23, 30, 50, 100, 300, 600, X) or c is cs[-1]) + ' |\n|---|' + '---|' * len([1 for c in cs if c['Y'] in (1, 2, 3, 5, 7, 11, 13, 17, 23, 30, 50, 100, 300, 600, X) or c is cs[-1]]))
        sel = [c for c in cs if c['Y'] in (1, 2, 3, 5, 7, 11, 13, 17, 23, 30, 50, 100, 300, 600, X) or c is cs[-1]]
        w('| eps | ' + ' | '.join('%s %s' % (sg(c), fe(c['eps'], 2)) for c in sel) + ' |')
        w('| matched 1e-3 | ' + ' | '.join(str(c.get('matched_1e_3')) for c in sel) + ' |')
        w('')
        pm = min(pos, key=lambda c: c['Y']) if pos else None
        nm = max(neg, key=lambda c: c['Y']) if neg else None
        s = 'largest Y with certified negative eigenvalue: %s (eps %s); smallest Y with certified positive form: %s' % (nm['Y'] if nm else None, fe(nm['eps'], 3) if nm else None, pm['Y'] if pm else None)
        if pm:
            s += ' (eps %s against %s for the full sum, ratio %s; matched %s against %s)' % (fe(pm['eps'], 3), fe(full['eps'], 3), fe(float(pm['eps']) / float(full['eps']), 3), pm.get('matched_1e_3'), full.get('matched_1e_3'))
        s += '; Y values positive: %d of %d; not certified: %s; non monotone (a negative cell above a positive one): %s.' % (len(pos), len(cs), unc, any(n['Y'] > p['Y'] for n in neg for p in pos))
        w(s + '\n')

# ---------------------------------------------------------------- accuracy
A = J('wl_acc.json')
Fit = J('wl_accuracy_fit.json')
if A:
    w('## 6. Accuracy of the zeros of xi_hat (wl_acc.json, wl_accuracy_fit.json)\n')
    w('e_k = log10 of the distance from the ordinate gamma_k to the nearest real zero of xi_hat (zeros searched in (0, pi N/L)); eigenvector from mpmath eigsy at 3.5N + 60 digits and above; ordinates by mp.zetazero at the working precision. eps = log10 of the least eigenvalue.\n')
    d = {(c['X'], c['N']): c for c in A if 'err_log10' in c}
    Xs = sorted(set(k[0] for k in d))
    Ns = sorted(set(k[1] for k in d))
    for X in Xs:
        w('X = %d (L = %.4f)\n' % (X, math.log(X) / 2))
        w('| N | log10 eps | pi N/L | ordinates | e_1 | e_2 | e_5 | e_10 | e_20 | e_last | e_1 - log10 eps |\n|---|---|---|---|---|---|---|---|---|---|---|')
        for N in Ns:
            if (X, N) not in d:
                continue
            c = d[(X, N)]
            es = c['err_log10']
            le = math.log10(abs(float(c['eps'])))
            g_ = lambda k: ('%.1f' % es[k - 1]) if len(es) >= k and es[k - 1] is not None else '--'
            w('| %d | %.1f | %.1f | %d | %s | %s | %s | %s | %s | %s | %s |' % (N, le, c['top'], len(es), g_(1), g_(2), g_(5), g_(10), g_(20), g_(len(es)) if es else '--', ('%.2f' % (es[0] - le)) if es and es[0] is not None else '--'))
        w('')
    if Fit:
        w('Fits (wl_accuracy_fit.json; least squares of e on the features, all entries with -D + 8 < e < -6; n = number of (X, N, k) entries):\n')
        w('| model | n | coefficients | rms of the residual | max abs residual | std of e |\n|---|---|---|---|---|---|')
        for f in [f for f in Fit.get('fits_all_k', []) if '(k = 1 only)' not in f['model']] + Fit.get('fits_k1', []):
            w('| %s | %d | %s | %.2f | %.2f | %.2f |' % (f['model'], f['n'], ', '.join('%.4g' % v for v in f['coef']), f['rms'], f['max_abs'], f['y_rms']))
        w('')
        for key, title in (('k1_slope_in_N_by_X', 'k = 1: digits gained per added mode (linear fit of e_1 against N at fixed X)'), ('k1_slope_in_L_by_N', 'k = 1: digits gained per doubling of lambda (fit of e_1 against log2 lambda at fixed N)')):
            if key in Fit:
                w(title + ': ' + '; '.join('%s: %.2f' % (k, v.get('digits_per_mode', v.get('digits_per_doubling_of_lambda'))) for k, v in Fit[key].items()) + '.\n')
        if Fit.get('by_k'):
            w('Digits lost per ordinate (slope of e_k against k at fixed X, N): ' + '; '.join('%s: %.2f (k up to %d, e_1 %.1f, e_last %.1f)' % (k, v['digits_per_ordinate'], v['klast'], v['e_k1'], v['e_klast']) for k, v in Fit['by_k'].items()) + '.\n')
        f = Fit.get('k1_NL_N_le_40')
        if f:
            w('k = 1, N <= 40 (before the saturation of eps): %s, n = %d, coefficients %s, rms %.2f, max abs %.2f.\n' % (f['model'], f['n'], ', '.join('%.4g' % v for v in f['coef']), f['rms'], f['max_abs']))
        for extra in ('summary',):
            if Fit.get(extra):
                w(Fit[extra] + '\n')

# ---------------------------------------------------------------- odd
O = J('wl_odd_scan.json')
if O:
    w('## 7. Odd modes of zeta (wl_odd_scan.json)\n')
    w('Odd modes psi_n = L^(-1/2) sin(pi n u/L), n = 1..N (they vanish at +-L), the pole term is -2 F_j(i/2) F_k(i/2) (negative semidefinite, see the header of wl_dh.py), the rest as for the even modes. eps = least eigenvalue; zeros of the odd transform G (F = iG) against the ordinates.\n')
    Xs = sorted(set(c['X'] for c in O))
    Ns = sorted(set(c['N'] for c in O))
    d = {(c['X'], c['N']): c for c in O}
    w('| X \\ N | ' + ' | '.join(str(n) for n in Ns) + ' |\n|---|' + '---|' * len(Ns))
    for X in Xs:
        w('| %d | ' % X + ' | '.join('%s %s' % (sg(d[(X, n)]), fe(d[(X, n)]['eps'], 2)) for n in Ns) + ' |')
    w('')
    w('| X | N | e_1 (log10 distance of the first zero of G to gamma_1) | matched / ordinates below pi N/L |\n|---|---|---|---|')
    for X in Xs:
        for n in Ns:
            c = d[(X, n)]
            es = c.get('err_log10') or []
            w('| %d | %d | %s | %s / %s |' % (X, n, ('%.1f' % es[0]) if es and es[0] is not None else '--', c.get('matched_1e_3'), c.get('n_ordinates')))
    w('\nAll %d cells certified positive: %s.\n' % (len(O), all(c.get('certified_positive') for c in O)))


# ---------------------------------------------------------------- brief (inserted after the title)
br = ['## 0. Results in brief (all numbers from the json files named below)\n']
if S1:
    evb = {(c['X'], c['N']): c for c in S1 + (S2 or []) if c['parity'] == 'even'}
    first = {}
    for X in sorted(set(k[0] for k in evb)):
        bad = sorted(n for (x, n), c in evb.items() if x == X and c.get('certified_negative'))
        okn = sorted(n for (x, n), c in evb.items() if x == X and c.get('certified_positive') and (not bad or n < bad[0]))
        first[X] = (bad[0] if bad else None, max(okn) if okn else None)
    br.append('- Davenport and Heilbronn (wl_dh_scan1.json, wl_dh_scan2.json): the form is certified positive on every even cell with N <= 52 (X up to 2000, N up to 80 for X = 30); it first fails (certified negative eigenvalue) at, per X: ' + '; '.join('X=%d: N=%s (last certified positive N=%s)' % (X, f[0], f[1]) for X, f in first.items() if f[0]) + '. X = 30 is positive up to N = 80 and X = 1000, 2000 up to N = 60 on this grid. DH odd modes: all %d cells positive.' % len([c for c in S1 if c['parity'] == 'odd']))
if S2:
    br.append('- How (section 4b, decomposition of W(f) over the zeros of f for the least eigenvector): in the last positive cells the zero at 0.8085 + 85.6993 i contributes a negative term of about the size of the line zero sum (it eats it); at the first negative cell its term exceeds the line sum and the eigenvector no longer reproduces the DH line zeros (matched 0 of the ordinates below pi N/L).')
if P:
    br.append('- Truncating the primes of zeta (wl_primes_scan.json): the form is certified negative for every Y below (nearly) X: the largest negative Y and smallest positive Y are given per (X, N) in section 5; eps is O(1) negative (not small) as soon as the prime sum is cut.')
if Fit and Fit.get('summary'):
    br.append('- Accuracy (wl_accuracy_fit.json): ' + Fit['summary'])
if O:
    br.append('- Odd modes (wl_odd_scan.json): %d of %d cells certified positive.' % (sum(1 for c in O if c.get('certified_positive')), len(O)))
out[2:2] = br + ['']

open(R + 'misc/reviews/nick-runs/checks-WL2.md', 'w').write('\n'.join(out) + '\n')
print('written', len(out), 'lines')
