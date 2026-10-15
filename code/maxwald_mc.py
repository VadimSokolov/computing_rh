# maxwald_mc.py: the simulation that the maximum of exponentials and the weighted sum of exponentials are the same law, for Section
# sec:maxwald (ch/maxwald.tex): the simulation paragraph of mw:sec:max and mw:fig:mc. Added for the book in October 2026 from the checks of
# the authors' chapter (misc/reviews/nick-maxwald/), recomputed from the mathematics with its own random streams (numpy, SeedSequence).
# M = max{E_n/(pi n), E'_m/(pi(2m-1)): n, m <= K} (K factors of each kind) is compared with M' = pi^{-1} sum_{n<=K} E''_n/n^2 + the mean
# sum_{n>K} 1/(pi n^2) of the remainder (Hurwitz zeta); both have the distribution function theta_4(t), mean pi/6 and variance pi^2/90.
#  1. N_REP independent replicates of the chapter's size (n = 2 x 10^5 draws, K = 400): the two sample Kolmogorov Smirnov statistic and
#     its p-value, the two means against pi/6, all nine deciles of each sample, the one sample statistic of each against theta_4.
#  2. For the first replicate (the chapter's configuration) the table of the nine deciles against the exact deciles, solutions of
#     theta_4(q) = p at 30 digits (mpmath), with the standard error sqrt(p(1-p)/n)/f(q) of a sample decile (f = theta_4') and the z score of
#     the difference of the two samples; the digits in which the deciles agree, -log10 of the largest relative difference.
#  3. One pooled run of N_BIG = 10^7 draws (the same K): the same table, to show how the agreement of the deciles grows with the sample
#     size (the replicates of 2 x 10^5 agree to about two digits in the deciles; the run of 14 October 2026 gave 2.9 digits at 10^7).
#  4. fig/maxwald_mc.pdf: histograms (density, 120 bins on [0, 2]) of M and M' for the first replicate, with the exact density theta_4'.
# Run from a folder holding this script: python3 maxwald_mc.py [--quick]; writes maxwald_mc.json and fig/maxwald_mc.pdf (created if absent).
# The tasks (one per replicate and one per 2 x 10^5 draws of the pooled run) go to a process pool with SLURM_CPUS_PER_TASK workers (one
# worker if unset). --quick is the laptop test (n = 2 x 10^4, K = 60, 2 replicates, pooled run 2 x 10^5). On Hopper (about 10 minutes, 8 cores):
# bash misc/tools/hopper_run.sh -c 8 -t 90 -g maxwald_mc.json -g fig/maxwald_mc.pdf code/maxwald_mc.py
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import mpmath as mp
from scipy.stats import ks_2samp, kstest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QUICK = '--quick' in sys.argv
K = 60 if QUICK else 400
N_CHAP = 20000 if QUICK else 200000
N_REP = 2 if QUICK else 20
N_BIG = 200000 if QUICK else 10000000
SEED = 20261014
DPS = 30
PS = [round(0.1 * j, 1) for j in range(1, 10)]
T0 = time.time()


def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)


def theta4_np(t):
    """theta_4 in double precision: series for t >= 0.4, Poisson summation below."""
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    big = t >= 0.4
    tb = t[big]
    s = np.ones_like(tb)
    for n in range(1, 14):
        s += 2 * (-1) ** n * np.exp(-np.pi * n * n * tb)
    out[big] = s
    ts = t[~big]
    s = np.zeros_like(ts)
    for n in range(0, 14):
        s += 2 * np.exp(-np.pi * (n + 0.5) ** 2 / ts)
    out[~big] = s / np.sqrt(ts)
    return out


def dtheta4_np(t):
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    big = t >= 0.4
    tb = t[big]
    s = np.zeros_like(tb)
    for n in range(1, 14):
        s += (-1) ** n * n * n * np.exp(-np.pi * n * n * tb)
    out[big] = -2 * np.pi * s
    ts = t[~big]
    s = np.zeros_like(ts)
    d = np.zeros_like(ts)
    for n in range(0, 14):
        e = 2 * np.exp(-np.pi * (n + 0.5) ** 2 / ts)
        s += e
        d += e * np.pi * (n + 0.5) ** 2 / ts ** 2
    out[~big] = d / np.sqrt(ts) - s / (2 * ts ** 1.5)
    return out


def simulate(task):
    """One block of n draws of M (maximum of 2K exponentials) and of M' (sum of K exponentials plus the mean of the remainder)."""
    seq, n, kk = task
    rng = np.random.default_rng(seq)
    k = np.arange(1, kk + 1, dtype=float)
    r1 = 1.0 / (np.pi * k)
    r2 = 1.0 / (np.pi * (2 * k - 1))
    r3 = 1.0 / (np.pi * k ** 2)
    rem = float(mp.zeta(2, kk + 1)) / np.pi
    mmax = np.empty(n)
    msum = np.empty(n)
    ch = 5000
    for b in range(0, n, ch):
        m = min(ch, n - b)
        e1 = rng.standard_exponential((m, kk))
        e1 *= r1
        e2 = rng.standard_exponential((m, kk))
        e2 *= r2
        mmax[b:b + m] = np.maximum(e1.max(axis=1), e2.max(axis=1))
        e3 = rng.standard_exponential((m, kk))
        e3 *= r3
        msum[b:b + m] = e3.sum(axis=1) + rem
    return mmax, msum


def exact_deciles():
    mp.mp.dps = DPS

    def th(t):
        t = mp.mpf(t)
        if t >= 1:
            return 1 + 2 * mp.nsum(lambda n: (-1) ** int(n) * mp.exp(-mp.pi * n * n * t), [1, mp.inf])
        return mp.nsum(lambda n: mp.exp(-mp.pi * (n + mp.mpf(1) / 2) ** 2 / t), [-mp.inf, mp.inf]) / mp.sqrt(t)

    res = []
    for p in PS:
        q = mp.findroot(lambda t: th(t) - mp.mpf(str(p)), (mp.mpf('0.2'), mp.mpf('0.35')) if p < 0.3 else (mp.mpf('0.4'), mp.mpf('0.6')),
                        solver='secant', tol=mp.mpf(10) ** (-DPS + 6))
        f = mp.diff(th, q)
        res.append((p, q, f))
    return res


def decile_table(mmax, msum, ex):
    n = len(mmax)
    qa = np.quantile(mmax, PS)
    qb = np.quantile(msum, PS)
    rows = []
    worst_rel = 0.0
    for i, (p, q, f) in enumerate(ex):
        qe, fe = float(q), float(f)
        se = np.sqrt(p * (1 - p) / n) / fe
        rel = abs(qa[i] - qb[i]) / qe
        worst_rel = max(worst_rel, rel)
        rows.append({'p': p, 'exact': mp.nstr(q, 20), 'max_sample': float(qa[i]), 'sum_sample': float(qb[i]),
                     'max_minus_sum': float(qa[i] - qb[i]), 'rel_diff': float(rel), 'se_of_a_sample_decile': float(se),
                     'z_two_samples': float((qa[i] - qb[i]) / (np.sqrt(2) * se)),
                     'z_max_vs_exact': float((qa[i] - qe) / se), 'z_sum_vs_exact': float((qb[i] - qe) / se)})
    return rows, worst_rel


def stats(mmax, msum):
    ks = ks_2samp(mmax, msum)
    k1 = kstest(mmax, theta4_np)
    k2 = kstest(msum, theta4_np)
    n = len(mmax)
    pp = np.array(PS)
    sd = np.sqrt(pp * (1 - pp) / n)
    zmean = np.pi / np.sqrt(90) / np.sqrt(n)
    return {'n': n, 'mean_max': float(mmax.mean()), 'mean_sum': float(msum.mean()), 'ks_two_sample': float(ks.statistic),
            'ks_two_sample_pvalue': float(ks.pvalue), 'ks_max_vs_theta4': float(k1.statistic), 'ks_sum_vs_theta4': float(k2.statistic),
            'std_max': float(mmax.std()), 'std_sum': float(msum.std()),
            'z_deciles_max_vs_theta4': [float(z) for z in (theta4_np(np.quantile(mmax, pp)) - pp) / sd],
            'z_deciles_sum_vs_theta4': [float(z) for z in (theta4_np(np.quantile(msum, pp)) - pp) / sd],
            'z_mean_max': float((mmax.mean() - np.pi / 6) / zmean), 'z_mean_sum': float((msum.mean() - np.pi / 6) / zmean)}


def main():
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    ntask_big = max(1, N_BIG // N_CHAP)
    seqs = np.random.SeedSequence(SEED).spawn(N_REP + ntask_big)
    tasks = [(seqs[i], N_CHAP, K) for i in range(N_REP + ntask_big)]
    log('K =', K, 'n =', N_CHAP, 'replicates', N_REP, 'pooled blocks', ntask_big, 'workers', nproc)
    reps = []
    big_a, big_b = [], []
    first = None
    if nproc > 1:
        pool = Pool(nproc)
        it = pool.imap(simulate, tasks)
    else:
        pool = None
        it = map(simulate, tasks)
    for i, (a, b) in enumerate(it):
        if i < N_REP:
            if i == 0:
                first = (a, b)
            reps.append(stats(a, b))
            reps[-1]['decile_max_abs_diff'] = float(np.max(np.abs(np.quantile(a, PS) - np.quantile(b, PS))))
            reps[-1]['decile_max_rel_diff'] = float(np.max(np.abs(np.quantile(a, PS) - np.quantile(b, PS)) / np.quantile(a, PS)))
            log('replicate', i, reps[-1]['ks_two_sample'], reps[-1]['mean_max'], reps[-1]['mean_sum'])
        else:
            big_a.append(a)
            big_b.append(b)
    if pool is not None:
        pool.close()
        pool.join()
    big_a, big_b = np.concatenate(big_a), np.concatenate(big_b)
    ex = exact_deciles()
    mp.mp.dps = DPS
    out = {'working_precision': 'exact deciles at 30 digits (mpmath); simulations in float64', 'K_factors_of_each_kind': K,
           'n_chapter_size': N_CHAP, 'n_replicates': N_REP, 'n_pooled': len(big_a), 'seed': SEED, 'quick': QUICK,
           'pi_over_6': float(np.pi / 6), 'std_exact_pi_over_sqrt90': float(np.pi / np.sqrt(90))}
    tab, w = decile_table(first[0], first[1], ex)
    out['first_replicate'] = {'stats': reps[0], 'deciles': tab, 'digits_of_agreement_of_deciles': float(-np.log10(w)),
                              'max_abs_decile_diff': max(abs(r['max_minus_sum']) for r in tab)}
    out['replicates'] = reps
    arr = lambda key: np.array([r[key] for r in reps])
    out['replicates_summary'] = {'ks_two_sample_mean': float(arr('ks_two_sample').mean()), 'ks_two_sample_min': float(arr('ks_two_sample').min()),
                                 'ks_two_sample_max': float(arr('ks_two_sample').max()),
                                 'ks_two_sample_pvalue_min': float(arr('ks_two_sample_pvalue').min()),
                                 'mean_max_range': [float(arr('mean_max').min()), float(arr('mean_max').max())],
                                 'mean_sum_range': [float(arr('mean_sum').min()), float(arr('mean_sum').max())],
                                 'decile_max_abs_diff_mean': float(arr('decile_max_abs_diff').mean()),
                                 'decile_max_abs_diff_max': float(arr('decile_max_abs_diff').max()),
                                 'expected_ks_scale_sqrt_2_over_n': float(np.sqrt(2.0 / N_CHAP)),
                                 'mean_se_of_a_sample_mean': float(np.pi / np.sqrt(90) / np.sqrt(N_CHAP))}
    tab, w = decile_table(big_a, big_b, ex)
    out['pooled_run'] = {'stats': stats(big_a, big_b), 'deciles': tab, 'digits_of_agreement_of_deciles': float(-np.log10(w)),
                         'max_abs_decile_diff': max(abs(r['max_minus_sum']) for r in tab)}
    zmax = np.array([r['z_deciles_max_vs_theta4'] for r in reps])
    zsum = np.array([r['z_deciles_sum_vs_theta4'] for r in reps])
    out['replicates_z_summary'] = {'note': 'z = (theta_4(sample decile) - p)/sqrt(p(1-p)/n); an unbiased simulation has mean 0 and standard deviation 1',
                                   'mean_z_per_decile_max': [float(x) for x in zmax.mean(0)], 'std_z_per_decile_max': [float(x) for x in zmax.std(0)],
                                   'mean_z_per_decile_sum': [float(x) for x in zsum.mean(0)], 'std_z_per_decile_sum': [float(x) for x in zsum.std(0)],
                                   'mean_z_of_the_sample_mean_max_sum': [float(arr('z_mean_max').mean()), float(arr('z_mean_sum').mean())],
                                   'std_z_of_the_sample_mean_max_sum': [float(arr('z_mean_max').std()), float(arr('z_mean_sum').std())]}
    f0 = reps[0]
    rng_ks = [float(arr('ks_two_sample').min()), float(arr('ks_two_sample').max())]
    se_mean = float(np.pi / np.sqrt(90) / np.sqrt(N_CHAP))
    digits_rep = -np.log10(np.maximum(arr('decile_max_rel_diff'), 1e-300))
    out['chapter_comparison'] = {
        'two sample KS statistic 0.0026 (L129)': {'chapter': 0.0026, 'ours_first_replicate': f0['ks_two_sample'], 'ours_replicates_min_max': rng_ks,
                                                  'typical_size_0.87_sqrt(2/n)': float(0.87 * np.sqrt(2.0 / N_CHAP)),
                                                  'chapter_value_inside_the_range_of_our_replicates': bool(rng_ks[0] <= 0.0026 <= rng_ks[1])},
        'means 0.52330 and 0.52382 against pi/6 = 0.52360 (L129)': {'chapter': [0.52330, 0.52382], 'pi_over_6': float(np.pi / 6),
                                                                    'ours_first_replicate': [f0['mean_max'], f0['mean_sum']], 'se_of_a_sample_mean': se_mean,
                                                                    'chapter_means_within_2_se_of_pi_over_6': bool(abs(0.52330 - np.pi / 6) < 2 * se_mean and abs(0.52382 - np.pi / 6) < 2 * se_mean)},
        'pi/6 = 0.52360 (L129)': {'chapter': 0.52360, 'ours': float(np.pi / 6), 'agree': bool(abs(np.pi / 6 - 0.52360) <= 0.5e-5)},
        'deciles agree to three digits (L129)': {'ours_digits_of_agreement_first_replicate': out['first_replicate']['digits_of_agreement_of_deciles'],
                                                 'ours_digits_mean_over_replicates': float(digits_rep.mean()),
                                                 'ours_replicates_with_at_least_3_digits': int((digits_rep >= 3).sum()), 'replicates': len(reps),
                                                 'ours_pooled_run_digits': out['pooled_run']['digits_of_agreement_of_deciles'], 'pooled_n': len(big_a),
                                                 'review_1_M7_quantiles_10_50_90_percent_of_the_chapters_run': {'max': [0.20748, 0.43568, 0.95292], 'sum': [0.20801, 0.43650, 0.95519]}}}
    # figure: the first replicate
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'axes.spines.top': False, 'axes.spines.right': False})
    os.makedirs('fig', exist_ok=True)
    bins = np.linspace(0, 2, 121)
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    ax.hist(first[0], bins=bins, density=True, histtype='step', color='black', lw=1.5, label='maximum of exponentials')
    ax.hist(first[1], bins=bins, density=True, histtype='step', color='#777777', lw=1.3, ls='--',
            label=r'$\pi^{-1}\sum_n E_n/n^2$ (sum)')
    xx = np.linspace(0.003, 2, 400)
    ax.plot(xx, dtheta4_np(xx), ':', color='#333333', lw=1.0, label=r'exact density $\theta_4^{\prime}$')
    ax.set_xlabel(r'$m$')
    ax.set_ylabel('density')
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig('fig/maxwald_mc.pdf')
    out['figure'] = {'file': 'fig/maxwald_mc.pdf', 'bins': 120, 'range': [0, 2], 'sample_size': N_CHAP}
    out['seconds'] = round(time.time() - T0, 1)
    with open('maxwald_mc.json', 'w') as f:
        json.dump(out, f, indent=1)
    log('done; wrote maxwald_mc.json and fig/maxwald_mc.pdf')


if __name__ == '__main__':
    main()
