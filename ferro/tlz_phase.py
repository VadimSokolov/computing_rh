# tlz_phase.py: TLZ.5 of reconcile-TLZ.md. The window of Proposition ev:prop:lbphase on the first 200 zeros: with theta the Riemann Siegel theta function, S_mid = n - 3/2 - theta(gamma_n)/pi
# (the midpoint of the jump of S at gamma_n: N = theta/pi + 1 + S), Delta = arctan tanh(pi gamma/2), Theta - theta = Delta + arg zeta(3/2 + i gamma) - arctan(3/(2 gamma)) = pi/4 + pi a(gamma),
# Re c_D(rho_n) < 0 exactly when sin(Theta - theta - pi S_mid) > 0 (sign of Z'(gamma_n) = (-1)^(n+1)); window: left when -3/4 + a < S_mid < 1/4 + a, right when 1/4 + a < S_mid < 5/4 + a, mod 2.
# Compared with the certified signs of Re c_D(rho_n) in data/lb_c.json (code/lb_c.py, ball arithmetic). mpmath at 30 digits, zeros from data/zeros_1700.txt (first 200).
# Written from the mathematics (the collaborator's scripts were not read, run or copied).
import json, sys, argparse
import mpmath as mp
mp.mp.dps = 30
ap = argparse.ArgumentParser(); ap.add_argument('zeros'); ap.add_argument('lbc'); ap.add_argument('--out', default='tlz_phase.json'); a = ap.parse_args()
gam = [mp.mpf(l.split()[0]) for l in open(a.zeros) if l.strip()][:200]
rows = json.load(open(a.lbc))['rows']
out, bad_sign, bad_window, amax = [], [], [], 0
for n, g in enumerate(gam, 1):
    th = mp.siegeltheta(g)
    S = n - mp.mpf(3) / 2 - th / mp.pi
    Delta = mp.atan(mp.tanh(mp.pi * g / 2))
    az = mp.arg(mp.zeta(mp.mpf(3) / 2 + 1j * g))
    pa = Delta + az - mp.atan(3 / (2 * g)) - mp.pi / 4
    A = pa / mp.pi
    amax = max(amax, abs(A))
    left_pred = mp.sin(mp.pi / 4 + pa - mp.pi * S) > 0
    Sm = (S + mp.mpf(3) / 4 - A) % 2 - mp.mpf(3) / 4 + A          # S_mid reduced into (-3/4 + a, 5/4 + a)
    left_window = Sm < mp.mpf(1) / 4 + A
    r = rows[n - 1]; assert r['n'] == n
    left_true = bool(r['neg'])
    if left_pred != left_true: bad_sign.append(n)
    if left_window != left_true: bad_window.append(n)
    margin = min(abs(Sm - (mp.mpf(1) / 4 + A)), abs(Sm - (-mp.mpf(3) / 4 + A)), abs(Sm - (mp.mpf(5) / 4 + A)))
    out.append({'n': n, 'gamma': mp.nstr(g, 12), 'S_mid': mp.nstr(S, 8), 'S_red': mp.nstr(Sm, 8), 'a': mp.nstr(A, 8), 'quarter_plus_a': mp.nstr(mp.mpf(1) / 4 + A, 8), 'margin_to_window_edge': mp.nstr(margin, 4), 'left_true_certified': left_true, 'left_pred_phase_law': bool(left_pred), 'left_window': bool(left_window)})
pos = [91, 63, 64, 65, 92, 93, 115, 126, 136, 160, 165, 183, 187]
res = {'n_checked': 200, 'sign_mismatches_phase_law': bad_sign, 'sign_mismatches_window': bad_window, 'max_abs_a': mp.nstr(amax, 6),
       'min_margin': min(float(o['margin_to_window_edge']) for o in out), 'rows': out}
json.dump(res, open(a.out, 'w'), indent=1)
print('mismatches phase law', bad_sign, 'window', bad_window, 'max|a|', res['max_abs_a'], 'min margin', res['min_margin'])
for o in out:
    if o['n'] in (63, 64, 65, 92, 93, 115, 126, 136, 160, 165, 183, 187, 1):
        print(o['n'], o['gamma'], 'S_mid', o['S_mid'], 'red', o['S_red'], '1/4+a', o['quarter_plus_a'], 'margin', o['margin_to_window_edge'], 'true_left', o['left_true_certified'], 'window_left', o['left_window'])
