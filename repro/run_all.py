# Run the book's scripts in dependency order, in parallel where possible, and log each run.
# Layout under BASE, the environment variable RH_REPRO_BASE or by default /scratch/vsokolov/rh_book_repro: the scripts of
# code/ and repro/zeros.py run in run/; the five programs of verification_reciprocal_zeta/ (enclosures.py,
# spot_checks.py, certify_k9_k2.py, harmonic_decay.py, sawtooth_truncation.py) run in verif/;
# code_bosch/verify_bosch_xi_v7.py runs in bosch/; and the scripts of ferro/ (Chapter 13) run in ferro/, with a copy
# of data/zeros_1700.txt. compare.py then compares the outputs with the archived ones, copied to
# book_data/ (data/), book_verif/ (verification_reciprocal_zeta/) and book_bosch/ (code_bosch/). Not run here: the
# programs of verification_reciprocal_zeta/torus/, which run as Slurm arrays by the steps of torus/README.md, and the
# packages of Part VI and of the Riesz program, which have their own drivers, run6.py and run_riesz.py.
import subprocess, time, json, os, sys
from concurrent.futures import ThreadPoolExecutor
BASE = os.environ.get('RH_REPRO_BASE', '/scratch/vsokolov/rh_book_repro')
RUN = f'{BASE}/run'; LOG = f'{BASE}/logs'
PY = '/scratch/vsokolov/rh_book_repro/venv/bin/python3'
for d in (LOG, f'{RUN}/fig', f'{RUN}/book/fig'):  # most plotting scripts do not create their folders
    os.makedirs(d, exist_ok=True)
TO = 24 * 3600
EPS = ['1.0', '0.5', '0.25', '0.1', '0.05', '0.02', '0.01', '0.005']
T = {}
def task(name, cmd, deps=(), cwd=RUN):
    T[name] = dict(cmd=cmd, deps=list(deps), cwd=cwd)
task('zeros', [PY, '-u', 'zeros.py'])
task('part1', [PY, '-u', 'part1.py'])
for e in EPS: task(f'part2_{e}', [PY, '-u', 'part2.py', e])
P2 = [f'part2_{e}' for e in EPS]
task('analyse', [PY, '-u', 'analyse.py'], ['zeros', 'part1'] + P2)
task('turing_ig', [PY, '-u', 'turing_ig.py'], P2)
task('deconv', [PY, '-u', 'deconv.py'], ['zeros', 'part2_0.1'])
task('clockgrid', [PY, '-u', 'clockgrid.py'], ['part1'])
task('hazard_spec', [PY, '-u', 'hazard_spec.py'], ['part1'])
task('rsboundary', [PY, '-u', 'rsboundary.py'], ['hazard_spec'])
task('hitting', [PY, '-u', 'hitting.py'], ['hazard_spec'])
task('hitting2', [PY, '-u', 'hitting2.py'], ['hazard_spec'])
task('volterra', [PY, '-u', 'volterra.py'], ['hazard_spec', 'clockgrid', 'rsboundary'])
task('mc', [PY, '-u', 'mc.py'], ['rsboundary', 'hazard_spec'])
task('mc2', [PY, '-u', 'mc2.py'], ['volterra', 'rsboundary', 'hazard_spec'])
task('dbn', [PY, '-u', 'dbn.py'], ['part1'])
task('dbnplot', [PY, '-u', 'dbnplot.py'], ['dbn'])
task('zeros_hp', [PY, '-u', 'zeros_hp.py'], ['zeros', 'part1'])
task('arith', [PY, '-u', 'arith.py'], ['zeros_hp'])
for s in ['bdlp_ladder', 'branch_coeffs', 'cm_orders', 'firstfail', 'gaussweil', 'heat_ladder', 'pakes_disc', 'pakes_sizebias', 'phasefree', 'vdpair', 'wigfree_upsilon']:
    task(s, [PY, '-u', f'{s}.py'], ['zeros_hp'])
task('firstfail_arb', [PY, '-u', 'firstfail_arb.py'], ['firstfail'])
# ball arithmetic certificates that read firstfail.json and compute their zeros with Arb
task('gsmooth_arb', [PY, '-u', 'gsmooth_arb.py'], ['firstfail'])
task('phasefree_arb', [PY, '-u', 'phasefree_arb.py'], ['firstfail'])
for s in ['hankel', 'kenttail', 'kenttail2', 'kmono', 'perturb', 'selberg', 'spacing']:
    task(s, [PY, '-u', f'{s}.py'], ['zeros'])
for s in ['biane_clock', 'biane_flow', 'biane_polya', 'biane_traces', 'dbeta_minimum', 'dbeta_thorin', 'dbeta_window', 'dbeta_zeros', 'dbn_collisions', 'dbn_delta', 'dbn_hankel_long', 'hcm', 'idsd_check', 'idsd_height', 'jacobi_dbn', 'jacobi_zeros', 'ladder_bound', 'lchi12', 'lehmer', 'lehmer_dbn', 'lfun', 'li', 'lsin2', 'maxexp_checks', 'maxexp_flow', 'maxexp_hb', 'maxexp_psizeros', 'maxlaw', 'pakes_lambert', 'primeside', 'saddle', 'sensitivity', 'sinh_bench', 'sqfree_atoms', 'sqfree_constants', 'sqfree_gammanu', 'sqfree_quadratic', 'sqfree_sieve', 'sqfree_smooth', 'turing1000', 'wigbiane_ineq', 'wigbiane_poly', 'wigfree_gamma2', 'wigfree_heat', 'wigfree_jacobi', 'williams_mc', 'tilted']:
    task(s, [PY, '-u', f'{s}.py'])
# Chapters 22 and 23 (scripts of code/mi_*.py and code/ki_*.py; no dependencies checked against the data they read)
for s in ['mi_friable', 'mi_humps', 'mi_survey', 'mi_depart', 'mi_miyatake', 'mi_gb', 'mi_ladder', 'mi_dh', 'mi_figs', 'ki_check']:
    task(s, [PY, '-u', f'{s}.py'])
# Chapter 24 (code/ht_*.py; the Hopper scripts are listed for completeness, see code/README.md; ht_lcert needs its nine --mj jobs first)
for s in ['ht_kent', 'ht_residues', 'ht_resid', 'ht_jensen', 'ht_jlaws', 'ht_dh', 'ht_ev', 'ht_euler32', 'ht_lcert', 'ht_dhsharp', 'ht_hcm2exp', 'ht_curved', 'ht_misc', 'ht_primebuilt', 'ht_figs']:
    task(s, [PY, '-u', f'{s}.py'])
# Chapter 25 (code/ev_*.py, default parts; ev_rouche5 and the counts part of ev_dbn are Hopper jobs, see code/README.md; ev_evt_lib.py is a library)
for s in ['ev_g1_fold', 'ev_g1_string', 'ev_g1_tail', 'ev_rouche5', 'ev_dbn', 'ev_prony2', 'ev_frac', 'ev_dev', 'ev_figs']:
    task(s, [PY, '-u', f'{s}.py'])
# Chapters 26 to 29 (code/wt_*.py, cg_checks.py, tr_*.py, lv_*.py; wt_*, tr_zeros, tr_counts, tr_hermite and lv_gram are Hopper jobs, see code/README.md; tr_counts takes 4 to 5 hours)
for s in ['wt_heat', 'wt_window', 'wt_rig', 'tr_zeros', 'tr_drift', 'tr_counts', 'tr_hermite', 'lv_sech_tables', 'lv_gram']:
    task(s, [PY, '-u', f'{s}.py'])
task('cg_checks', [PY, '-u', 'cg_checks.py', '--out', 'cg_checks.json'])
# Chapters 30 to 33 (code/lo_*.py, cx_*.py, es_*.py, ex_zeros3e12.py, cd_slopes.py; lo_*, cx_u_certify, cx_min, cx_em, es_pencil, es_zeros1e13 and ex_zeros3e12 are Hopper jobs, see code/README.md and the checks-*.md of misc/reviews/nick-b4; es_pair, es_bf, es_detect* and es_window are listed there and not run here; lo_common.py is a library)
task('lo_flow', [PY, '-u', 'lo_flow.py', '--procs', '8', '--out', 'lo_flow.json'])
task('lo_em_table', [PY, '-u', 'lo_em.py', 'table', '--procs', '8', '--out', 'lo_em_table.json'])
task('lo_em_inject', [PY, '-u', 'lo_em.py', 'inject', '--procs', '8', '--out', 'lo_em_inject.json'])
task('lo_gbc', [PY, '-u', 'lo_gbc.py', '--procs', '8', '--out', 'lo_gbc.json'])
for s in ['cx_u_certify', 'cx_slope_X', 'cx_min', 'cx_em', 'cd_slopes']:
    task(s, [PY, '-u', f'{s}.py'])
# Chapters 34 and 35 (code/cf_*.py and bsy_*.py; Hopper jobs, see code/README.md and the headers of the scripts; each has a --quick mode for a laptop, not used here; cf_eval.py and bsy_zdh.py are libraries, bsy_figs.py only reads data files and is not run here)
for s in ['cf_dhzeros', 'cf_lp', 'cf_arith', 'cf_moments', 'cf_xi', 'cf_info', 'cf_ns']:
    task(s, [PY, '-u', f'{s}.py'])
task('bsy_cm', [PY, '-u', 'bsy_cm.py', '--out', 'bsy_cm.json'])
task('bsy_levy', [PY, '-u', 'bsy_levy.py', '--out', 'bsy_levy.json'])
for s in ['bsy_count', 'bsy_defects', 'bsy_strip', 'bsy_defl', 'bsy_yks', 'bsy_mc', 'bsy_tailbound', 'bsy_riesz2d', 'bsy_defl_dh', 'bsy_vertical', 'bsy_negw']:
    task(s, [PY, '-u', f'{s}.py'])
task('bsy_halfdisk', [PY, '-u', 'bsy_halfdisk.py', '--procs', '16', '--out', 'bsy_halfdisk.json'])
task('bsy_wmc', [PY, '-u', 'bsy_wmc.py', '--procs', '12', '--out', 'bsy_wmc.json'])
task('bsy_cuts', [PY, '-u', 'bsy_cuts.py'])
task('bsy_isbound', [PY, '-u', 'bsy_isbound.py'])
task('bsy_cm_5000', [PY, '-u', 'bsy_cm.py', '--rows', 'zeta:0.5:5000', '--out', 'bsy_cm_5000.json'])
# bsy_trunc.py runs once per theta and the two outputs are merged (see its header)
task('bsy_trunc_t75', [PY, '-u', 'bsy_trunc.py', '--theta', '0.75', '--Ts', '10,40,160,640,1280', '--out', 'bsy_trunc_t75.json'])
task('bsy_trunc_t50', [PY, '-u', 'bsy_trunc.py', '--theta', '0.5', '--Ts', '10,40,160,640', '--out', 'bsy_trunc_t50.json'])
task('bsy_trunc', [PY, '-u', 'bsy_trunc.py', '--merge', 'bsy_trunc_t75.json', 'bsy_trunc_t50.json', '--out', 'bsy_trunc.json'], ['bsy_trunc_t75', 'bsy_trunc_t50'])
task('ki_tilt_book', [PY, '-u', 'ki_tilt.py', '--norm', 'book', '--out', 'ki_tilt_book.json'])
task('ki_tilt_chapter', [PY, '-u', 'ki_tilt.py', '--norm', 'chapter', '--out', 'ki_tilt_chapter.json'])
# sqfree_sieve.py needs about 8 GB of memory
# wigbiane_flow.py reads dbn_collisions.json
task('wigbiane_flow', [PY, '-u', 'wigbiane_flow.py'], ['dbn_collisions'])
# pakes_filter.py takes about 15 minutes on 32 cores; pakes_filter_cert.py (python-flint) and pakes_characters.py a few minutes
task('pakes_filter', [PY, '-u', 'pakes_filter.py', '32'])
task('pakes_filter_cert', [PY, '-u', 'pakes_filter_cert.py', '16'])
task('pakes_characters', [PY, '-u', 'pakes_characters.py', '8'])
# ttest.py and ttest2.py take the command lines given in their headers (repro/reconstruct/ttest2/)
task('ttest', [PY, '-u', 'ttest.py', '0.1', '200', '0.004', '4'])
task('ttest2_0.10', [PY, '-u', 'ttest2.py', '0.10', '200', '0.004', '4', '200', '500', '1000'])
task('ttest2_0.02', [PY, '-u', 'ttest2.py', '0.02', '220', '0.001', '4.2', '1000', '2000', '5000'])
task('sinh_bench_closed', [PY, '-u', 'sinh_bench.py', 'closed'])
task('hcm_arb', [PY, '-u', 'hcm_arb.py'], ['hcm'])
task('lplots', [PY, '-u', 'lplots.py'], ['lsin2', 'lfun'])
task('maxexp_figs', [PY, '-u', 'maxexp_figs.py'], ['maxexp_psizeros', 'maxexp_flow'])
task('wigfree_plots', [PY, '-u', 'wigfree_plots.py'], ['wigfree_jacobi', 'wigfree_upsilon', 'wigfree_heat'])
# biane_figs.py reads biane_polya.json, biane_flow.json, biane_clock.json and biane_traces.json; biane_flow.py takes about 22 minutes
task('biane_figs', [PY, '-u', 'biane_figs.py'], ['biane_polya', 'biane_flow', 'biane_clock', 'biane_traces'])
# dbeta_figs.py reads dbeta_zeros.json
task('dbeta_figs', [PY, '-u', 'dbeta_figs.py'], ['dbeta_zeros'])
task('plots', [PY, '-u', 'plots.py'], ['analyse', 'hcm', 'primeside'])
task('plots2', [PY, '-u', 'plots2.py'], ['hazard_spec', 'clockgrid', 'rsboundary', 'mc'])
task('threezones', [PY, '-u', 'threezones.py'], ['volterra', 'rsboundary', 'mc2'])
task('plots3', [PY, '-u', 'plots3.py'], ['zeros', 'arith', 'primeside', 'analyse', 'turing_ig', 'hankel', 'spacing', 'sinh_bench'])
task('plots4', [PY, '-u', 'plots4.py'], ['analyse'])
task('deconvplot', [PY, '-u', 'deconvplot.py'], ['deconv'])
task('lehmerflowplot', [PY, '-u', 'lehmerflowplot.py'])
for s in ['enclosures', 'spot_checks', 'certify_k9_k2', 'harmonic_decay', 'sawtooth_truncation']:
    task(s, [PY, '-u', f'{s}.py'], cwd=f'{BASE}/verif')
task('bosch_v7', [PY, '-u', 'verify_bosch_xi_v7.py'], ['zeros'], cwd=f'{BASE}/bosch')
FERRO = f'{BASE}/ferro'
for s in ['fm_polya', 'fm_interval', 'fm_hcm', 'fm_sizebias']:
    task(s, [PY, '-u', f'{s}.py'], cwd=FERRO)
task('fm_hermite', [PY, '-u', 'fm_hermite.py'] + [str(k) for k in (1, 2, 3, 5, 10, 16, 17, 20, 30, 50, 100, 150, 200, 300, 400, 600, 800)], cwd=FERRO)
task('fm_hermite_all', [PY, '-u', 'fm_hermite.py', 'all', '200'], cwd=FERRO)
task('fm_maxwell', [PY, '-u', 'fm_maxwell.py', '10000000'], cwd=FERRO)
# ball arithmetic certificates of Chapter 13; fm_qcert.py zeros writes fm_zeros_hp.txt (800 zeros and Var X from Arb)
task('fm_cutoff', [PY, '-u', 'fm_cutoff.py'], cwd=FERRO)
task('fm_qcert_zeros', [PY, '-u', 'fm_qcert.py', 'zeros', '800', '9600'], cwd=FERRO)
task('fm_qcert_table', [PY, '-u', 'fm_qcert.py', 'table'], ['fm_qcert_zeros', 'fm_hermite'], cwd=FERRO)
task('fm_qcert_all', [PY, '-u', 'fm_qcert.py', 'all', '200'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_qcert_minq', [PY, '-u', 'fm_qcert.py', 'minq', '1-200,300,400,600,800'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_tvcert', [PY, '-u', 'fm_tvcert.py'], ['fm_qcert_zeros'], cwd=FERRO)
task('fm_intcert', [PY, '-u', 'fm_intcert.py'], cwd=FERRO)
# note W (multipliers on Polya's kernel): wf_fit.py writes wf_fit.json, which the other scripts read; wf_tvcert.py also reads fm_tvcert.json
task('wf_fit', [PY, '-u', 'wf_fit.py'], cwd=FERRO)
task('wf_fit_printed', [PY, '-u', 'wf_fit_printed.py'], ['wf_fit'], cwd=FERRO)
task('wf_tvcert', [PY, '-u', 'wf_tvcert.py'], ['wf_fit', 'fm_tvcert'], cwd=FERRO)
task('wf_clock', [PY, '-u', 'wf_clock.py'], ['wf_fit'], cwd=FERRO)
task('wf_zeros', [PY, '-u', 'wf_zeros.py'], ['wf_fit', 'wf_tvcert'], cwd=FERRO)
# Note FC (ferro/fc_*.py: Stern counts, Perron table, Fisher zeros, Farey chain; Hopper jobs, fc_stern about 4 minutes on 8 cores, fc_fisher about 2 hours on 16 cores, all need python-flint and numpy; fc_perron, fc_fisher, fc_extra read fc_births.npz written by fc_stern and zeros_1700.txt)
task('fc_stern', [PY, '-u', 'fc_stern.py', '--nmax', '34', '--nproc', '8', '--out', 'fc_stern.json'], cwd=FERRO)
task('fc_yk', [PY, '-u', 'fc_stern.py', '--yk', '--ykmax', '8', '--out', 'fc_yk.json'], cwd=FERRO)
task('fc_perron', [PY, '-u', 'fc_perron.py', '--nmax', '34', '--nproc', '8', '--out', 'fc_perron.json'], ['fc_stern'], cwd=FERRO)
task('fc_fisher', [PY, '-u', 'fc_fisher.py', '--nmax', '34', '--nproc', '16', '--certmax', '26', '--argmax', '28', '--out', 'fc_fisher.json'], ['fc_stern'], cwd=FERRO)
task('fc_wind', [PY, '-u', 'fc_extra.py', '--mode', 'wind', '--nlist', '30,32,34', '--nmax', '34', '--npts', '160', '--nproc', '12', '--out', 'fc_wind.json'], ['fc_stern'], cwd=FERRO)
task('fc_means', [PY, '-u', 'fc_extra.py', '--mode', 'means', '--nmax', '34', '--nproc', '8', '--out', 'fc_means.json'], ['fc_stern'], cwd=FERRO)
task('fc_ly', [PY, '-u', 'fc_ly.py', '--nmax', '30', '--certmax', '26', '--nproc', '8', '--tables', 'fc_ly_tables.npz', '--out', 'fc_ly.json'], cwd=FERRO)
task('fc_ly_hi', [PY, '-u', 'fc_ly.py', '--nmin', '28', '--nmax', '36', '--certmax', '32', '--nproc', '8', '--tables', 'fc_ly_hi_tables.npz', '--out', 'fc_ly_hi.json'], cwd=FERRO)
task('fc_q1', [PY, '-u', 'fc_extra.py', '--mode', 'q1', '--tables', 'fc_ly_tables.npz', 'fc_ly_hi_tables.npz', '--nproc', '2', '--out', 'fc_q1.json'], ['fc_ly', 'fc_ly_hi'], cwd=FERRO)
task('rx_euler', [PY, '-u', 'rx_euler.py'])
task('rx_centre', [PY, '-u', 'rx_centre.py'])
task('rx_tail', [PY, '-u', 'rx_tail.py'])
task('rx_maslanka', [PY, '-u', 'rx_maslanka.py'])
task('rx_blocks', [PY, '-u', 'rx_blocks.py'])
# Note P (code/hx_*.py: the shifted Gauss Hermite evaluator of xi_t and the checks of the merged text of Section cp:sec:minflow; hx_shift reads xg_tilt.json and hx_half zeros_8000.txt from data/; hx_scan and hx_far use all cores, about 10 to 30 minutes on 40 cores; hx_cramer needs python-flint)
for s in ['hx_shift', 'hx_half', 'hx_edge', 'hx_cramer']:
    task(s, [PY, '-u', f'{s}.py'])
for s in ['hx_scan', 'hx_far']:
    task(s, [PY, '-u', f'{s}.py'], ['hx_shift'])
# Note X3 (code/pg3_*.py: Li cumulants, velocity of the zeros, triangle mixture of the Polya law; they read zeros_8000.txt from the run folder; pg3_velocity (--cert) and pg3_fejer need python-flint, pg3_li about 40 s)
for s in ['pg3_reach', 'pg3_li', 'pg3_fejer']:
    task(s, [PY, '-u', f'{s}.py'])
task('pg3_velocity', [PY, '-u', 'pg3_velocity.py', '--cert'])
task('pg3_velocity_flow', [PY, '-u', 'pg3_velocity.py', '--flow'])
# Notes on pole cancellation and the excursion maximum (code/px_*.py; px_lib.py is their library; px_cbox, px_xiw and px_dyck are Hopper jobs on 8 cores, about 1 minute each; px_const, px_cbox, px_xiw and px_dyck need python-flint, px_xiw and px_dyck read zeros_1700.txt from the run folder)
task('px_const', [PY, '-u', 'px_const.py'])
task('px_cbox', [PY, '-u', 'px_cbox.py'])
for s in ['px_xiw', 'px_dyck']:
    task(s, [PY, '-u', f'{s}.py'], ['zeros'])
# Lattice bias block of Chapter 25 (code/lb_*.py; lb_c is a Hopper job on 1 core, 1 second, lb_zero a Hopper job on 3 to 8 cores, 2 seconds to 2 minutes; both need python-flint, read zeros_1700.txt from the run folder and import px_lib.py; lb_zero needs lb_c.py beside it)
task('lb_c', [PY, '-u', 'lb_c.py', 'zeros_1700.txt', '--out', 'lb_c.json'], ['zeros'])
task('lb_zero3', [PY, '-u', 'lb_zero.py', 'zeros_1700.txt', '--ns', '6400,25600', '--ks', '1-13', '--workers', '8', '--out', 'lb_zero3.json'], ['zeros'])
task('lb_zero2', [PY, '-u', 'lb_zero.py', 'zeros_1700.txt', '--ns', '100000000,30000000,10000000,3000000,1000000,300000,100000', '--ks', '63,64,65', '--cont', '--workers', '3', '--out', 'lb_zero2.json'], ['zeros'])
# Note TLZ (ferro/tlz_*.py: Section ch:fmmellin of Chapter 13 and, for tlz_phase, Chapter 25; tlz_lib.py is their library; tlz_k and tlz_onecos are Hopper jobs on 2 cores, about 1 and 4 minutes, tlz_deform and tlz_bessel Hopper jobs on 1 core, about 1 minute; all need python-flint and mpmath; tlz_phase runs on a laptop in under 30 seconds and reads zeros_1700.txt and lb_c.json from the run folder)
for s in ['tlz_k', 'tlz_onecos', 'tlz_deform', 'tlz_bessel']:
    task(s, [PY, '-u', f'{s}.py', '--out', f'{s}.json'], cwd=FERRO)
task('tlz_phase', [PY, '-u', 'tlz_phase.py', f'{RUN}/zeros_1700.txt', f'{RUN}/lb_c.json', '--out', 'tlz_phase.json'], ['zeros', 'lb_c'], cwd=FERRO)
# Chapters nc_* and Chapter 22 (code/nc_*.py, mi_levy.py, mi_reweight.py, bc_*.py; the nc_cauchy and nc_route scripts write their json into the run folder; floating point laptop scripts except bc_cert (Hopper, 8 cores, 35 s) and bc_margin (Hopper, 1 core); not run here, listed in code/README.md: the work packages wl_*.py, wp2_*.py, wp3_*.py, wp4_hankel.py, wp5_*.py (Hopper array jobs and report scripts, see misc/reviews/nick-runs/checks-*.md) and nc_cauchy_two_pairs_fig.py, nc_route_toeplitz_fig.py (figures from json of other scripts))
for s in ['mi_levy', 'mi_reweight', 'nc_barrier_check', 'nc_qrh_numbers', 'nc_zdist_checks', 'nc_riemanncont_checks', 'nc_riemanncont_half', 'nc_riemanncont_scan', 'nc_route_genus2', 'nc_route_hankel_order', 'nc_route_hankel_order_hi', 'nc_cauchy_twospin', 'bc_margin']:
    task(s, [PY, '-u', f'{s}.py'])
task('bc_cert', [PY, '-u', 'bc_cert.py', '--workers', '8', '--out', 'bc_cert.json'])
for s in ['curves', 'cw_planted', 'exact_criterion', 'field_zero', 'gap_theory', 'multifield', 'multiplicity', 'phase', 'planted_pos', 'planted', 'strip_ly', 'strip_search', 'survey', 'thorin_compare', 'tilt_zeros', 'track', 'two_pairs_close', 'two_pairs', 'uniform']:
    task(f'nc_cauchy_{s}', [PY, '-u', f'nc_cauchy_{s}.py'])
status = {}
def run(name):
    t = T[name]
    for d in t['deps']:
        if status.get(d, {}).get('rc') != 0:
            status[name] = dict(rc=None, skipped=f'dependency {d} failed'); return
    t0 = time.time()
    with open(f'{LOG}/{name}.out', 'w') as out:
        try:
            rc = subprocess.run(t['cmd'], cwd=t['cwd'], stdout=out, stderr=subprocess.STDOUT, timeout=TO).returncode
        except subprocess.TimeoutExpired:
            rc = 'timeout'
    status[name] = dict(rc=rc, seconds=round(time.time() - t0, 1))
    print(f'{name}: rc={rc} {status[name]["seconds"]}s', flush=True)
    json.dump(status, open(f'{BASE}/status.json', 'w'), indent=1)
# simple scheduler: repeatedly launch tasks whose dependencies are finished
from concurrent.futures import wait, FIRST_COMPLETED
pending = set(T); running = {}
with ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))) as ex:
    while pending or running:
        ready = [n for n in pending if all(d in status for d in T[n]['deps'])]
        for n in ready:
            pending.discard(n); running[ex.submit(run, n)] = n
        if not running: break
        done, _ = wait(list(running), return_when=FIRST_COMPLETED)
        for f in done: running.pop(f)
json.dump(status, open(f'{BASE}/status.json', 'w'), indent=1)
print('ALL DONE', flush=True)
