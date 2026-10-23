# Scripts, data and run times of the book

This file lists the scripts, the data files and the run times behind the tables, figures and certificates of the book *Computing the Riemann Hypothesis*. It holds the inventory that the book does not print: Appendix A of the book says where the code and the data are and describes the two Pólya evaluators, and this file describes everything else, chapter by chapter.

The sections follow the chapters of the book. The section of a chapter holds a table with one row for each section of the chapter that has scripts: the scripts, with the data file that each writes in brackets, and the tables, figures and results of the book that they produce. Chapters, sections, theorems, tables and figures are numbered as in the PDF of the same release. A script is named relative to the folder given in the text, `code/` unless another folder is named. The sections after the chapters give the reruns of the archive, the gaps that remain and the references.

## Folders of the repository

A number in this book can be checked once three things are known: the evaluator that produced it, the precision at which it ran, and the script that turned its output into a table or a figure. The book's archive of scripts and data files is in its repository, <https://github.com/VadimSokolov/computing_rh>, beside the LaTeX source of the book. The scripts of Chapters 1 to 5 and 8 are in the folder `code/`, and the data files they write are in `data/`. The interval enclosures, spot checks and other scripts of Chapters 6 and 7 are in `verification_reciprocal_zeta/`; the verification script of Section 4.3 is in `code_bosch/`, with an identical copy beside the authors' notes in `bosch_hcm_xi/`; the scripts of Chapters 9 to 11 are in `basepoint_one/`, `sato_bondesson/` and `riesz_program/`; and the reruns and reconstructions made for the book are in `repro/`.

The scripts of Chapters 12 and 14 to 21 are in `code/` as well, and their data files in `data/`, except those of Chapter 13, which are in `ferro/`.

## Libraries, the cluster and how the scripts run

The scripts run with Python 3.11 and the libraries of `requirements.txt` in the repository root: python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.6, scipy 1.17.1, sympy 1.14.0 and matplotlib 3.11.2, the versions of the reruns of September and October 2026. Some later runs used mpmath 1.4.1, as stated under Chapter 21. The long computations ran on the Hopper cluster of George Mason University; the drivers in `repro/` and the SLURM files name the paths of that cluster and need adapting elsewhere.

Ball arithmetic comes from `flint`, the Python interface to FLINT and Arb (Johansson 2017), which also gives the seventy digit ordinates of the first $`1700`$ zeros used by `arith.py`. Reference values of $`\zeta`$, $`\Gamma`$ and $`\xi`$, the other ordinates of the zeros, the interval arithmetic of Chapters 6 and 7 and the Talbot inversion of Laplace transforms (Talbot 1979) come from `mpmath` (mpmath development team 2023); `numpy` and `scipy` supply arrays, quadrature, interpolation, root finding, differential equations and nonnegative least squares, `sympy` the primes and the Möbius function, and `matplotlib` the figures.

## How the scripts read and write their files

The scripts run from `code/` and read and write their data files in the working directory; the data files named in the tables of the chapters below are in `data/`, except the reports and JSON files of Chapters 6 and 7, which are in `verification_reciprocal_zeta/` and its subfolder `torus/`. The Pólya evaluators use no zeros. Ten scripts, `analyse.py`, `deconv.py`, `deconvplot.py`, `hankel.py`, `kenttail.py`, `kenttail2.py`, `kmono.py`, `perturb.py`, `plots3.py` and `selberg.py`, load the first $`1700`$ ordinates of the zeros from `g_1_400.npy`, `g_401_1000.npy` and `g_1001_1700.npy`. These lists are not in `code/`: `repro/zeros.py` computes them with `mpmath.zetazero` at thirty digits, and identical copies are in `riesz_program/zeros/` and in several of the folders of the notes (see Chapter 9 below). `spacing.py` reads the lists when they are present and otherwise computes and saves them, and `zeros_hp.py`, which writes the seventy digit ordinates of `zeros_1700.txt` that `arith.py`, `cm_orders.py`, `firstfail.py`, `firstfail_arb.py`, `gaussweil.py`, `phasefree.py` and `vdpair.py` read, uses the lists and `part1.json` only to check its output. Two scripts come first, because others read their output. `part1.py` finds the sign changes of $`\Xi`$ below height $`200`$ on a grid of step $`0.05`$, refines the zeros by Newton's method and computes $`N(200)`$ by the argument principle. `part2.py` evaluates $`\pi\rho_\varepsilon(\theta)=\operatorname{Re}(\xi'/\xi)(\frac12+\varepsilon+i\theta)`$ and the continuous argument of $`\xi(\frac12+\varepsilon+i\theta)`$ for $`0\le\theta\le200`$, at each resolution $`\varepsilon`$ given on its command line, and saves them in the files `eps_*.npz` with the largest argument increment, which it also prints. The book uses the eight resolutions $`1`$, $`0.5`$, $`0.25`$, $`0.1`$, $`0.05`$, $`0.02`$, $`0.01`$ and $`0.005`$, and the largest increment over all of them, $`0.29`$ at $`\varepsilon=1`$, is the one quoted in Section 2.1.2. The distances of Section 3.2 (`analyse.py`), the zero density identity of Section 3.4 (`turing_ig.py`) and the deconvolution of Section 8.6 (`deconv.py` and `deconvplot.py`) start from these files, and `analyse.py` writes `plotcurves.npz`, which `plots.py`, `plots3.py` and `plots4.py` read; the `.npz` files are not in the archive. `analyse.py`, `clockgrid.py`, `hazard_spec.py` and `dbn.py` read `part1.json`, and the boundary computations of Section 8.2 start from `hazard_spec.json`. `hcm_arb.py`, which encloses in ball arithmetic the $`216`$ values that `hcm.py` computes, reads `hcm.json` to compare them. Likewise `firstfail_arb.py`, which certifies in ball arithmetic the negative values that `firstfail.py` finds, reads `firstfail.json`, and so do `gsmooth_arb.py` and `phasefree_arb.py`, for the pairs and their heights; these two compute the zeros they need with Arb instead of reading them.

## Figures and reconstructed scripts

The figures are in colour in the digital edition, the PDF at the DOI on the copyright page, and in grey in the paperback; where a figure tells its curves apart by colour alone, the digital edition is the one to read. `plots.py` takes no arguments; it reads `plotcurves.npz`, `analysis.json`, `hcm.json` and `results_prime.json`. `dbnplot.py` reads `dbn.json`. `plots.py`, `plots2.py`, `perturb.py` and `dbnplot.py` save their figures in `fig/` below the working directory, and `plots3.py`, `plots4.py`, `kmono.py`, `deconvplot.py`, `lplots.py` and `lehmerflowplot.py` in `book/fig/`; of these only `perturb.py` and `dbnplot.py` create their folders. Four scripts of `code/` are reconstructions, written for the book where the archive held a data file or a module but no script; the header of each says so and gives its command line, its inputs and what it reproduces. `hankel.py` and `spacing.py` reproduce `hankel.json` and `spacing.json` byte for byte, and `spacing.py` also draws `spacing.pdf` in the working directory. `sinh_bench.py` reproduces `sinh_bench.json`, whose lattice sums are cut at $`|k|\le400000`$, to a relative $`1.1\times10^{-12}`$, and `python3 sinh_bench.py closed` writes the closed forms, which Sections 1.3 and 5.4 quote, to `sinh_bench_closed.json`. `polya.py` wraps `setup(420, 200, 3)` and `xi_pair` of `polya_flint.py` for `primeside.py`, which imports it.

## Constants used in several sections

The table below collects the constants used in several sections and names the file of the archive in which each can be checked; the printed values agree with these files to every digit shown.

| quantity | value | file |
|---|---|---|
| $`\xi(\frac12)`$ | 0.4971207781883141 | `williams_mc.json` |
| $`\sum_ka_k^{-1}=\mathbb{E}T`$ | 0.0231050 | `kenttail2.json` |
| $`\sum_ka_k^{-2}=\mathrm{Var}\,T`$ | $`3.71726\times10^{-5}`$ | `kenttail2.json` |
| $`\gamma_1^2`$ | 199.790454832386859463 | `hazard_spec.json` |
| $`\gamma_2^2-\gamma_1^2`$ | 242.1357 | `hazard_spec.json` |
| $`\lambda_1=b_1=\sum_\rho\rho^{-1}`$ | 0.023095708966121033814 | `li.json` |
| $`c=\zeta'(\frac12)/\zeta(\frac12)`$, the constant of (6.1.6) | 2.6860917096 | `spot_checks_report.md` |
| $`\Xi'(200)/\Xi(200)`$ | $`-1.2506879`$ | `extra.json` |
| $`\kappa_{\rm cell}`$, the cell constant of the $`79`$ cells below height $`200`$ | 5.43665 | `extra.json` |
| $`N(200)`$ | 79 | `part1.json` |
| $`N(1000)`$ | 649 | `turing1000.json` |
| de Bruijn Newman collision time of $`\gamma_{71},\gamma_{72}`$, the closest pair below height $`200`$ | $`-0.2880`$ | `dbn.json` |
| collision time of the Lehmer pair near height $`7005`$ | $`-7.11\times10^{-4}`$ | `lehmer_dbn.json` |
| $`N_{\chi_{12}}(100)`$, the number of zeros of $`L(s,\chi_{12})`$ with $`0<\operatorname{Im}\rho\le100`$ | 67 | `lchi12.json` |

## Chapter 1: Introduction and foundations

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 1.1 | `plots4.py` | Figure 1.1 |
| 1.2 | `primeside.py`, `polya.py` (`results_prime.json`), `arith.py` (`arith.json`), `vdpair.py` (`vdpair.json`), `plots3.py` | Table 1.5; Figure 1.2 |
| 1.3 | `sinh_bench.py` (`sinh_bench.json`), `plots3.py` | Figure 1.3 |

## Chapter 2: Evaluating the Xi function

### The Pólya evaluators

The first evaluator is `code/polya_flint.py`, the evaluator of Section 2.1 on the real $`u`$ axis. Its function `setup` computes the node values $`\Phi_j=\Phi(jh)`$, $`0\le j\le J`$, from the theta series (2.1.1) in ball arithmetic and stores the polynomials $`P`$ and $`Q`$ of Section 2.1.1 as `PA` and `PB`; `xi_pair` returns $`P(E)+P(E^{-1})`$ and $`Q(E)-Q(E^{-1})`$ with $`E=e^{sh}`$, the trapezoid approximations to $`\xi(\frac12+s)`$ and $`\xi'(\frac12+s)`$. The defaults are $`240`$ bits, $`h=1/200`$ and $`|u|\le3`$, so that $`J=600`$. Most runs call `setup(420, 200, 3)`, the parameters of Section 2.1.1; the Keiper Li coefficients of Section 8.3 use $`600`$ bits, and the Selberg computation of Section 3.5 uses $`300`$.

The second evaluator is `code/tilted.py`, which implements Proposition 2.1.1. With $`\eta=\pi/4-\delta`$ and nodes $`u_j=x_j+i\eta`$, $`x_j=jh`$, $`|j|\le J`$, on the line $`\operatorname{Im}u=\eta`$, where $`J=\lceil x_{\max}/h\rceil`$ for the cut $`|x|\le x_{\max}`$, so that there are $`2J+1`$ nodes, it returns $`e^{i\eta s}h\sum_j\Phi(u_j)e^{sx_j}`$ and $`e^{i\eta s}h\sum_ju_j\Phi(u_j)e^{sx_j}`$, the trapezoid approximations to $`\xi(\frac12+s)`$ and $`\xi'(\frac12+s)`$ on the tilted line; `Phi_c` sums the theta series at complex $`u`$, using the evenness of $`\Phi`$ when $`\operatorname{Re}u<0`$. Of the two evaluators in the file, `xi_pair` applies Horner's rule and is called only by `ttest.py`, which exhibits the growth of the ball radii described in Section 2.1.5; the production scripts evaluate each exponential directly, through `xi_pair_direct` or an equivalent loop. The Lehmer pair and its de Bruijn Newman flow use $`\delta=0.01`$, $`256`$ bits, $`h=6\times10^{-4}`$ and $`|x|\le3.9`$; the count at height $`1000`$ in Section 2.2.3.2 uses $`\delta=0.1`$, $`200`$ bits, $`h=0.004`$ and $`|x|\le4.2`$; the module `clock.py`, which evaluates the Laplace transform of the clock for Sections 8.1 and 8.2, switches to the tilted evaluator with $`\delta=0.1`$, $`320`$ bits, $`h=0.003`$ and $`|x|\le4.5`$ when $`|\operatorname{Im}s|\ge120`$; and the runs of Table 2.2 take $`\delta`$, the precision, $`h`$, the cut and the heights from the command line of `ttest2.py`. The command lines, recorded at the top of the two test scripts, are `python3 ttest2.py 0.10 200 0.004 4 200 500 1000`, with $`2001`$ nodes, for the rows with $`\delta=0.1`$, `python3 ttest2.py 0.02 220 0.001 4.2 1000 2000 5000`, with $`8401`$ nodes, for the rows with $`\delta=0.02`$, and `python3 ttest.py 0.1 200 0.004 4`, which prints `nan` at the heights $`200`$ to $`2000`$, where the Horner ball of $`\xi`$ contains zero.

In both evaluators the radii of the returned balls bound the rounding errors of the finite sums, and nothing else. Three errors lie outside them: the discretisation error of the trapezoid rule, bounded in Section 2.1.1; the cut of the axis at $`|u|\le3`$, where $`\Phi<10^{-500}`$, or at $`|x|\le x_{\max}`$ on the tilted line; and the truncation of the theta series, whose loop stops at the first term with midpoint below $`2^{-b-20}`$ at working precision $`b`$ and does not add the tail to the radius. The stopping test compares double precision numbers, so for $`b\ge1055`$ its threshold underflows to zero and the test is never met. The untilted route needs more than $`1054`$ bits for thirty digits above height about $`840`$ (Table 2.3), and there the test has to be made in ball arithmetic. The runs in this book use at most $`600`$ bits and are not affected.

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 2.1 | `polya_flint.py`, `part1.py` (`part1.json`), `part2.py`, `tilted.py`, `ttest.py`, `ttest2.py`, `plots.py`, `plots3.py` | the zeros and $`N(200)`$ below height $`200`$; Table 2.2; Figure 2.1 |
| 2.2 | `lehmer.py` (`lehmer.json`), `turing1000.py` (`turing1000.json`) | the Lehmer pair near height $`7005`$; the count of Section 2.2.3.2 |
| 2.3 | `polya_flint.py`, `plots3.py` | Figure 2.3 |

## Chapter 3: The Thorin measure at finite resolution

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 3.1 | `primeside.py` (`results_prime.json`), `plots3.py` | Table 3.1 |
| 3.2 | `part2.py`, `analyse.py` (`analysis.json`), `plots.py` (`extra.json`) | Table 3.2; Figure 3.1 |
| 3.3 | `sensitivity.py` (`results_sens.json`), `plots.py` (`sens2.json`) | Figure 3.2 |
| 3.4 | `turing_ig.py` (`turing_ig.json`), `plots3.py` | Table 3.3 |
| 3.5 | `selberg.py` (`selberg.json`) | figures and tables no longer in the book |

## Chapter 4: Weak limits, hyperbolic complete monotonicity and a Bosch type theorem

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 4.1 | `hcm.py` (`hcm.json`), `hcm_arb.py` (`hcm_arb.json`), `plots.py` | Figure 4.1; the $`216`$ values of Section 4.1.2 and their certified signs |
| 4.3 | `code_bosch/verify_bosch_xi_v7.py`, `pakes_disc.py` (`pakes_disc.json`) | Tables 4.2, 4.3; Figure 4.2; the checks of Theorem 4.3.9 and Remark 4.3.10; the table of blocks below |

### The verification script of the Bosch type theorem

The script behind the numerical checks of Section 4.3, `code_bosch/verify_bosch_xi_v7.py`, is in the book's archive with its output, `verify_bosch_xi_v7_output.txt`, and the line numbers below are those of these two files. The script consists of thirteen blocks that run one after another in one process, each with its own imports; blocks 8 and 9 use `numpy` as well as `mpmath`, block 8 passes quadrature data to block 9, and block 13 calls the functions that block 12 defines. Run in an empty folder, the blocks write their results to nine text files, one per block except that block 8 writes none, blocks 5 to 7 share a file and so do blocks 9 and 10, and the output file joins these files in block order. The table below gives, for each block, its lines in the script and in the output and the statements of Section 4.3 that it checks.

The blocks of the verification script: their lines in the script and in the output, what they check, and the statements of Section 4.3 they support. Output lines 37 to 39 and 109 come from code that is not in the script, and lines 93 to 102 from an earlier version of block 13.

| block | script | output | what is checked | statement |
|---|---|---|---|---|
| 1 | 1 to 50 | 1 to 10 | HCM test by the signs of the Taylor coefficients of pair functions, for Bosch's family and for a toy with zeros $`14e^{\pm0.1i}`$ and $`21`$; minimum of $`{\operatorname{Im}\varepsilon_\varphi(re^{i\theta})/r}`$ over $`885`$ points | Theorems 4.3.3, 4.3.5; Corollary 4.3.7 |
| 2 | 51 to 86 | 11 to 24 | minimum of $`\operatorname{Im}\varepsilon_h`$ over $`49`$ radii and $`199`$ angles, for Bosch's family and for the toy | Theorems 4.3.2, 4.3.3, 4.3.5 |
| 3 | 87 to 116 | 25 to 36 | $`u_p`$ from (4.3.2) against (4.3.3) with $`30`$ zeros; $`\log\varphi(s^p)`$ recovered from $`u_p`$ | Theorem 4.3.12 |
| none |  | 37 to 39 | (4.3.3) with $`100`$ and $`200`$ zeros; $`u_p>0`$ at $`148`$ points | Theorem 4.3.12; Corollary 4.3.14 |
| 4 | 117 to 148 | 40 to 49 | density of $`{X_C=\sum_ke_k/(2k-1)^2}`$ against $`(4/\pi)\eta(4it/\pi)^3`$, and its Laplace transform at $`s=2`$; a sufficient inequality in the prime region | Theorem 4.3.20; Proposition 4.3.18 |
| 5, 6 | 149 to 199 | 50 to 60 | residue coefficients $`c_n`$ from $`200`$ zeros, their signs, sizes and gaps; the residue series against Talbot inversion | Theorem 4.3.22 and the numerical paragraph after it |
| 7 | 200 to 219 | 61 to 66 | the bound $`m_p(t)`$ of (4.3.4); $`{-\zeta'/\zeta(\sigma)<1/(\sigma-1)}`$ | Proposition 4.3.18; Table 4.2 |
| 8 to 10 | 220 to 318 | 67 to 80 | $`g`$ by quadrature and by series; Pick scan of $`-g'/g`$; $`(\log g)''`$ on $`[0,1]`$ | Proposition 4.3.26; Theorem 4.3.27 |
| 11 | 319 to 343 | 81 to 86 | the pair kernel; constants for $`{T_0=3\times10^{12}}`$ | Theorem 4.3.25 |
| 12, 13 | 344 to 419 | 87 to 108 | Pick strip $`h_*`$ and threshold $`\beta_*=\pi/h_*`$ in logistic models with frequencies removed | Table 4.3 |
| none |  | 109 | gaps among the first $`1000`$ zeros | Question 4.3.28 |

The script runs as it stands. In an empty folder it ends after about $`400`$ seconds, and its nine files, joined in block order, give $`95`$ lines, each of which occurs in the output file in the same order and agrees with it byte for byte, apart from the timings in seconds that some lines report. The output file is nevertheless not the record of a single run, and its other $`14`$ lines come from code that is not in the script. Lines 37 to 39 (the zero sums with $`100`$ and $`200`$ zeros and the positivity of $`u_p`$ at $`148`$ points) and line 109 (the gaps among the first $`1000`$ zeros) are reproduced exactly by reconstructions kept in the archive (see Reruns and gaps below), although the grid of $`148`$ points is a guess. Lines 93 to 102 come from an earlier version of block 13 with ten cases. Lines 103 to 107 run five of these cases again with the present block 13: for the removed frequencies $`\{3\}`$ they repeat line 102, and for $`\{2,3\}`$, $`\{2,3,5,6\}`$, $`\{3,4\}`$ and $`\{4,5\}`$ they give $`\pi/h_*=2.99997`$, $`3.0`$, $`1.48313`$ and $`1.31671`$ where lines 94, 97, 99 and 100 give $`1.0`$, the value printed when the scan finds no failure below $`\pi`$. Lines 103 to 106 therefore supersede those four lines; the present block 13, run on all ten cases, reproduces lines 93, 95, 96, 98, 101 and 102, from which Table 4.3 takes its other rows. Line 108 adds the logistic model with no frequency removed.

For the toy with zeros $`14e^{\pm0.1i}`$ and $`21`$, neither test locates the threshold $`p_*=0.9363`$. The sign test finds no violation at $`p=0.95`$ (line 8), and the scan of block 2 misses the poles at $`p=0.94`$ and $`p=0.97`$, which lie between its grid points or beyond its last ray (lines 22 and 23, where the flag for poles is computed from their position, not detected); both tests see the failure only at $`p=1`$ (lines 9 and 24). At $`t=0.003`$ the terms of the residue series reach $`2.7\times10^{20}`$ while the density is about $`3.3\times10^{-72}`$, so about $`92`$ digits cancel, more than the $`25`$ carried in line 56 and the $`90`$ carried in line 60: the series values at $`t=0.003`$ in these lines are rounding error, and the Talbot values there are negative. At $`150`$ digits the residue series over the first $`250`$ zeros and Talbot inversion agree on the density $`3.32007507753\times10^{-72}`$ at $`t=0.003`$; the reconstructions described under Reruns and gaps contain this computation.

## Chapter 5: The heat trace, the arithmetic and the limits of finite information

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 5.1 | `idsd_check.py` (`idsd_check.json`), `idsd_height.py` (`idsd_height.json`), `ladder_bound.py` (`ladder_bound.json`), `jacobi_zeros.py` (`jacobi_zeros.json`), `cm_orders.py` (`cm_orders.json`), `cm_orders_k116e25.py` (`cm_orders_k116e25.json`), `prime_side_check.py` (`prime_side_check.json`), `window_orders.py` (`window_orders.json`), `hankel.py` (`hankel.json`), `saddle.py` (`saddle.json`), `kmono.py` (`kmono.json`), `gaussweil.py` (`gaussweil.json`), `plots3.py`, `xg_defect.py` (`xg_defect.json`) | the proof of Theorem 5.1.1, the remark after it and the one range proof; the check of (5.1.1), Lemma 5.1.6, Theorem 5.1.7 and the numbers after it; Theorem 5.1.10; Theorem 5.1.13 and the unnumbered table after it; Tables 5.1, 5.2, 5.3; Figures 5.1, 5.2; the rate bound of Theorem 5.1.11 |
| 5.2 | `perturb.py` (`perturb.json`), `firstfail.py` (`firstfail.json`), `firstfail_arb.py` (`firstfail_arb.json`), `gsmooth_arb.py` (`gsmooth_arb.json`), `xg_local.py` (`xg_local.json`) | Table 5.4; Figure 5.3; the failure times of $`\mathcal G`$; the table of Section 5.2.4 |
| 5.3 | `zeros_hp.py` (`zeros_1700.txt`), `arith.py` (`arith.json`), `plots3.py` | Table 5.6 |
| 5.4 | `spacing.py` (`spacing.json`), `plots.py` (`extra.json`), `plots3.py` | the cell constant $`\kappa_{\rm cell}`$; Figures 5.4, 5.5 |

## Chapter 6: The reciprocal zeta function

In `verification_reciprocal_zeta/`, `enclosures.py` computes the enclosures E1 to E6 of Chapters 6 and 7 in outward rounded interval arithmetic and writes `enclosures_report.md`. `certify_k9_k2.py` encloses in ball arithmetic the value of $`\operatorname{Re}\kappa(\rho_{213})`$ used in Corollary 7.2.4 (both from its definition and from the formula of Theorem 7.2.3), the sign of $`\operatorname{Re}\kappa`$ at the first $`460`$ zeros, and the values of $`S`$ quoted in and after the proof of Theorem 7.1.3, and writes `certify_k9_k2.json`; it runs in a few seconds. `spot_checks.py` recomputes the constants of these chapters that `mpmath` and the first few hundred zeros can settle and prints each beside the value in the text. It writes no file: its report is its output, `python3 spot_checks.py > spot_checks_report.md`, which takes about eight minutes, and it checks the nine rows of Table 7.7. Two scripts compute numbers quoted in Section 7.1, each in about a minute. `sawtooth_truncation.py` computes the root mean square of $`S`$ and of the errors of its truncated sawtooth and classical series on uniform grids of steps $`10^{-3}`$, $`5\times10^{-4}`$ and $`2.5\times10^{-4}`$ in $`[10,50]`$, and the sums $`\frac12\sum_{N\le N_0}|B_N|`$, quoted after Proposition 7.1.10, and writes `sawtooth_truncation.json`. `harmonic_decay.py` computes the maxima of $`\log|\hat f_{\kappa^*}(v)|/v`$ over the windows $`[0.95V,1.05V]`$ quoted after Proposition 7.1.17, by root finding in each interval between consecutive zeros of $`\hat f_{\kappa^*}`$, where the quotient has a single maximum, rather than on a grid, and writes `harmonic_decay.json`. The subfolder `torus/` holds the branch and bound of Remark 7.2.12. It ran as arrays of tasks on a cluster, and its proofs in ball arithmetic took about twenty CPU hours; `torus/README.md` says what each script does and how to rerun it. Tables that collect statements, or evaluate the formulas given in their captions, need no script; the section Reruns and gaps below names the other tables and figures of Chapters 1 to 8 that the tables of the chapters below do not list.

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 6.1 to 7.2 | `verification_reciprocal_zeta/enclosures.py`, `spot_checks.py` (`enclosures_report.md`, `spot_checks_report.md`), `sawtooth_truncation.py` (`sawtooth_truncation.json`), `harmonic_decay.py` (`harmonic_decay.json`), `certify_k9_k2.py` (`certify_k9_k2.json`), `torus/` (`certified_bounds.json`, `arb_thresholds.json` and six other JSON files), `code/xg_kernel.py` (`data/xg_kernel.json`) | the enclosures E1 to E6; Table 7.7; the numbers after Propositions 7.1.10 and 7.1.17; $`\operatorname{Re}\kappa(\rho_{213})`$ in the proof of Corollary 7.2.4 and the values of $`S`$ quoted with Theorem 7.1.3; Tables 7.3 and 7.4 and the second condition of Theorem 7.2.10; the numbers after Proposition 7.1.23 and in Examples 7.1.24 and 7.1.25 |

## Chapter 7: Infinite divisibility and truncations of the reciprocal zeta function

The scripts of this chapter are listed under Chapter 6.

## Chapter 8: Probability, hitting times and other criteria

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 8.1 | `clock.py`, `clockgrid.py` (`clockgrid.json`), `hazard_spec.py` (`hazard_spec.json`), `kenttail.py` (`kenttail.json`), `kenttail2.py` (`kenttail2.json`), `williams_mc.py` (`williams_mc.json`), `plots2.py` | Tables 8.3, 8.4, 8.2; Figure 8.1 |
| 8.2 | `turing_ig.py`, `rsboundary.py` (`rsboundary.json`), `mc.py` (`mc.json`), `volterra.py` (`volterra.json`), `mc2.py` (`mc2.json`), `hitting.py`, `hitting2.py`, `plots2.py`, `threezones.py` | the check of Proposition 8.2.1; Figure 8.2; the boundary of Section 8.2.2.2 |
| 8.3 | `li.py` (`li.json`), `plots3.py` | Table 8.5; Figure 8.3 |
| 8.4 | `dbn.py` (`dbn.json`), `lehmer_dbn.py` (`lehmer_dbn.json`), `jacobi_dbn.py` (`jacobi_dbn.json`), `dbn_collisions.py` (`dbn_collisions.json`), `dbn_hankel_long.py` (`dbn_hankel_long.json`), `wigbiane_ineq.py` (`wigbiane_ineq.json`), `wigbiane_poly.py` (`wigbiane_poly.json`), `wigbiane_flow.py` (`wigbiane_flow.json`), `dbnplot.py`, `lehmerflowplot.py`, `xg_tilt.py` (`xg_tilt.json`) | Tables 8.6, 8.7, 8.8, 8.9, 8.10, 8.11, 8.12; the Hankel test along the flow and its counts to $`t=-20`$; the collisions of the first $`200`$ gaps; the tests of the Biane and deadline bounds; the table and the numbers of Section 8.4.7 |
| 8.4.9 | `newton_central_ratio.py` (`newton_central_ratio.json`) | the table after Proposition 8.4.21: the central ratio $`c_n`$ at $`t=\tfrac12`$ and $`t=0`$ for $`n=10^2`$ to $`10^6`$ against $`\tfrac18\mathcal W(n/\pi)^2`$ |
| 8.5 | `lchi12.py` (`lchi12.json`), `lfun.py` (`lfun.json`), `lsin2.py` (`lsin2.json`), `lplots.py` (`lrho.json`), `maxlaw.py` (`maxlaw.json`), `dbn_delta.py` (`dbn_delta.json`), `maxexp_checks.py` (`maxexp_checks.json`), `maxexp_psizeros.py` (`maxexp_psizeros.json`), `maxexp_hb.py` (`maxexp_hb.json`), `maxexp_flow.py` (`maxexp_flow.json`), `maxexp_figs.py`, `sqfree_quadratic.py` (`sqfree_quadratic.json`), `dbeta_zeros.py` (`dbeta_zeros.json`), `dbeta_thorin.py` (`dbeta_thorin.json`), `dbeta_minimum.py` (`dbeta_minimum.json`), `dbeta_figs.py` (`dbeta_figs.json`) | Table 8.14; a table that the scripts write but the book no longer prints (it is reproduced below); Tables 8.15, 8.16, 8.17, 8.18, 8.19; Figures 8.4, 8.5, 8.6, 8.7, 8.8 |
| 8.6 | `deconv.py` (`deconv.json`), `deconvplot.py` | Tables 8.20, 8.21; Figure 8.9 |
| 8.7 | `biane_polya.py` (`biane_polya.json`), `biane_clock.py` (`biane_clock.json`), `biane_flow.py` (`biane_flow.json`), `biane_traces.py` (`biane_traces.json`), `biane_figs.py` (`biane_figs.json`) | Tables 8.22, 8.23, 8.24; Figures 8.10, 8.11, 8.12; in Section 10.1, Tables 10.3, 10.4 and Figure 10.2 |
| 8.1 | `bn_tail.py` (`bn_tail_a.json`, `bn_tail_b.json`, `bn_tail_c.json`, `bn_tail_d.json`; reads `zeros500.json`) | the saddlepoint tails of $`H`$ and the two expansions of its density, the Talbot comparison of the Möbius series at $`30`$ digits and the terms of the prime $`2`$ (the numbers of the subsection on the reciprocal gamma factor) |

### The cost of the three evaluators

The table below compares the evaluator for $`L(s,\Delta)`$ of Section 8.5.2 (step $`1/80`$, truncation $`|u|\le4`$) with those for $`\zeta`$ and $`\chi_{12}`$ at the same precision. In the variables of that section their steps are $`1/40`$ and $`1/160`$, and their truncations $`|u|\le6`$ and $`|u|\le2.5`$. The settings are not matched: the truncation is far longer than $`60`$ digits require for $`\zeta`$, and too short for them for $`\chi_{12}`$, whose kernel is still about $`10^{-16}`$ at $`|u|=2.5`$ (Section 8.5.3).

| | $`\zeta`$ | $`L(s,\chi_{12})`$ | $`L(s,\Delta)`$ |
|---|---|---|---|
| kernel | theta series with terms of both signs | $`\eta`$, a positive product | $`\eta^{24}`$, a positive product |
| decay of $`\xi`$ on the critical line | $`e^{-\pi t/4}`$ | $`e^{-\pi t/4}`$ | $`e^{-\pi t/2}`$ |
| pole of $`L`$ | at $`s=1`$ | none | none |
| degree, conductor | $`1`$, $`1`$ | $`1`$, $`12`$ | $`2`$, $`1`$ |
| zeros below height 50 | 10 | 28 | 20 |
| first zero | 14.134725 | 3.804628 | 9.222379 |
| digits lost per unit height | 0.31 | 0.36 | 0.59 |
| digits lost per zero near height 40 | about 1 | about 0.5 | about 1 |
| nodes | 481 | 801 | 641 |
| time per point | 9.5 ms | 12.8 ms | 10.7 ms |

The trapezoid evaluators for $`\zeta`$, $`L(s,\chi_{12})`$ and $`L(s,\Delta)`$ at $`60`$ digits, with the settings given above. The decay of $`\xi(\frac12+it)`$ is shown up to a power of $`t`$; the digits lost per unit height are the slope of $`\log_{10}|\xi(\frac12+it)|`$ between heights $`10.3`$ and $`40.3`$; the times are per point, in `mpmath`.

The kernels are evaluated once, before the timing, so the times measure an exponential sum over the nodes and follow the node counts of the table rather than any intrinsic difference: the theta series, the eta series and the $`q`$ expansion of $`\Delta`$ all converge geometrically, with ratio at most $`e^{-2\pi}`$, once the modular symmetry is used. Per unit height $`\Delta`$ is the most expensive, because its gamma factor $`\Gamma(s+\frac{11}2)`$ decays like $`e^{-\pi t/2}`$ on the critical line, twice as fast as the $`e^{-\pi t/4}`$ of $`\zeta`$ and $`\chi_{12}`$: the untilted evaluator loses asymptotically $`\pi/(2\log10)=0.68`$ digits per unit height against $`\pi/(4\log10)=0.34`$. Between heights $`10`$ and $`40`$ the powers of $`t`$ in the completed functions and the size of $`|L|`$ at the two ends give the measured rates of the table. Tilting the contour to within $`\delta`$ of the edge of the strip reduces the loss to about $`\delta\theta/\log10`$ digits at height $`\theta`$ (Proposition 2.1.1, whose proof applies to both eta kernels), so the wider strip of the $`\Delta`$ kernel offsets its faster decay up to a factor polynomial in $`1/\delta`$; we have not run the tilted evaluator for $`\Delta`$. Per zero $`L(s,\chi_{12})`$ is the cheapest: its conductor packs the zeros at low height.

## Chapter 9: The basepoint one

### The folders of the notes

The scripts of Chapters 9 to 11 stay in the folders of the authors' notes, and each reads and writes its data files in its own folder; the exceptions are Sections 9.1.6, 10.3, 11.2 and 11.5, whose scripts are in `code/` and whose data files are in `data/`, and Section 10.1.5, computed by the scripts of Section 8.7 in the table of Chapter 8 above. In `basepoint_one/`, the folder `xi` holds the scripts for $`\xi`$ of Chapter 9, and `eta` those for $`L(s,\chi_{12})`$ of Section 11.3; `scale_mixture` and `wiener` hold those of the two halves of Section 10.1, and `process_map` those of the map of processes of Section 11.1 and of the boundaries across the basepoints of Section 9.2.6. The folder `sato_bondesson/` holds the other scripts of Section 11.1, and `riesz_program/code/` those of Sections 10.2 and 11.4, whose data files are in `riesz_program/data/`; the script `setup.sh` links these files and the zero lists of `riesz_program/zeros/` into `code/`, from which the scripts run. Several names recur across the folders, so the tables below that have a folder column give the folder of every script: `comp.py` is in `scale_mixture`, `process_map` and `sato_bondesson`, `figs.py` in `xi`, `eta`, `process_map` and `sato_bondesson`, `improve.py` and `figs2.py` in `xi` and `scale_mixture`, and `salem_riesz.py` in `wiener` and `riesz_program/code`, while `mc.py` of `xi` and `arith.py` of `eta` are not the scripts of those names in `code/`. The scripts take no arguments, with three exceptions: `bnd_alpha.py` takes the basepoint, $`0.6`$ or $`0.75`$; `weil.py` takes the weights $`\log p`$ in double precision, or at thirty digits when called as `python3 weil.py mp`; and `riesz_large.py` takes $`x`$, the truncation $`K`$ and the twist $`h`$, as in `python3 riesz_large.py 1e10 200 50`.

In `xi`, `compute.py` writes the Thorin density at the basepoint one, the sine squared check with the primes up to $`10^7`$ and the values of Proposition 9.2.1 to `alpha1.json`, and each other script that computes writes the data file of its own name, except `rs1c.py`, which adds the tangent boundary to `rs1.json`. `mc.py` also writes `surv.json` and a sample `X_sample.npy`, and `rs1d.py` and `errformula.py` read `rs1.json` and `below.json`. The sample is not in the archive, so `figs.py`, which draws a figure no longer in the book from it, runs only after `mc.py`. The script `heat_primes.py` was added for the book: it reads no file, computes the thresholds of Table 9.5 in floating point and certifies Proposition 9.1.17 in ball arithmetic with `python-flint`, in about a minute and a half on $`32`$ cores. In `eta`, `zeros12.py` computes the $`297`$ zeros of $`L(s,\chi_{12})`$ below height $`340`$ and saves them in `g12.npy`, which `base.py`, `arith.py` and `hit.py` read; `hit2.py` adds to `hit.json`, and `figs.py` reads `alpha1.json` of `xi` as well as the files of its own folder.

In `scale_mixture`, `dens2.py` adds the densities of $`Y_\alpha`$ to `strip.json`, and `figs2.py` draws Figure 10.1 from `strip.json` and `comp.json`. In `wiener`, two scripts were added for the book: `figs_wiener.py`, which draws Figure 10.3, whose drawing script the authors' package did not contain, and `riesz_explicit.py`, a reconstruction that writes `riesz_explicit.json` from `salem_riesz.json`, byte for byte as archived. In `process_map`, `comp.py` and `comp2.py` read a copy of `alpha1.json`, and `figs.py` draws Figure 11.2 from `ou.json`, the files `bnd_0.6.json` and `bnd_0.75.json`, and copies of `rs1.json` from `xi` and of `volterra.json` from `data/`. In `sato_bondesson/`, `comp.py` reads `alpha1.json`, `rs1.json` and `improve.json` of `basepoint_one/xi` and `clockgrid.json` and `hazard_spec.json` of `data/`, and `weil.py` prints its results and writes no file; `figs.py` holds as constants the five values of the explicit formula that a figure no longer in the book showed as dots.

In `riesz_program/code/`, `killing2.py` adds to `killing.json`, which `rmt_ratios.py` reads, `more2.py` reads `almost.json`, and `figs_riesz.py` reads `riesz_thorin.json`, `salem_riesz.json` and `riesz_explicit.json`, the last two identical to the files of `basepoint_one/wiener`. Seven scripts of the folder were added for the book: `v1min.py`, which certifies Proposition 9.1.4; `explicit_bounds.py`, which certifies in ball arithmetic the numbers used in the proofs of Theorems 10.1.3(3) and 11.4.4; `almost_dip.py`, `check_large_twisted.py`, `ggc_offset.py` and `rmt_ratios.py`, which compute numbers of Sections 10.2 and 11.4 that had no archived output; and `figs_almost.py`, which redraws Figures 11.9 and 11.8 in the notation of the book. Two checks were added in `repro/`: `ks_ak.py` recomputes the arithmetic factors that `ks.py` holds as constants, and `zeros_mp.py` and `v1_grid_mp.py`, run as arrays of tasks on a cluster and collected by `v1_collect.py` in `v1_check.json`, check Proposition 9.1.4 independently with `mpmath`. The figures of Chapters 9 to 11 are copies, renamed in `fig/`, of the PDF files that these scripts write: `bp1_*.pdf` from `xi`, `scale_mixture`, `wiener` and `process_map`, `bp1eta_*.pdf` from `eta`, `sb_*.pdf` from `sato_bondesson`, and `ri_*.pdf` and `ag_*.pdf` from `riesz_program/code`, where `figs_strategy.py` writes the panels of Figure 10.4. The labels inside the other plots keep the notation of the notes.

Scripts of Chapter 9 and Section 10.1, the data files they write, and the tables and figures they produce. The folders `xi`, `process_map`, `scale_mixture` and `wiener` are in `basepoint_one/`; the scripts of `riesz_program` are in its folder `code/` and their data files in its folder `data/`, and those of Section 9.1.6 are in `code/` and `data/` of the archive.

Scripts of Sections 10.2 and 10.3 and Chapter 11, as in the note above; the folders `process_map` and `eta` are in `basepoint_one/`, and the scripts of Sections 10.3, 11.2 and 11.5 are in `code/` and their data files in `data/`.

Scripts of the notes with the folder that holds each, the data files they write in brackets, and the tables and figures they produce.

| section | folder | scripts (data files) | tables and figures |
|---|---|---|---|
| 9.1 | `xi` | `compute.py` (`alpha1.json`), `figs.py` | Table 9.2; Figure 9.1 |
| 9.1 | `xi` | `levy.py` (`levy.json`), `figs3.py` | Figure 9.2 |
| 9.1 | `xi` | `below.py` (`below.json`), `errformula.py` (`errformula.json`), `figs4.py` | Table 9.3; Figure 9.3 |
| 9.1 | `xi` | `heat_primes.py` (`heat_primes.json`) | Table 9.5; Proposition 9.1.17 |
| 9.1 | `xi` | `heat_thresholds.py` (`heat_thresholds.json`) | the thresholds of Table 9.5, in ball arithmetic |
| 9.1 | `process_map` | `comp.py` (`comp.json`) | Table 9.1 |
| 9.1 | `riesz_program` | `v1min.py` (`v1min.json`) | Proposition 9.1.4 and the values after its proof |
| 9.1 | `repro` | `zeros_mp.py`, `v1_grid_mp.py`, `v1_collect.py` (`v1_check.json`) | the check of Proposition 9.1.4 |
| 9.1 | `repro` | `valpha_min.py` (`valpha_min.json`) | Proposition 9.1.5 and the values after its proof |
| 9.1.6 | `code` | `dbeta_minimum.py` (`dbeta_minimum.json`), `dbeta_window.py` (`dbeta_window.json`) | Table 9.4; the window count in the proof of Theorem 9.1.11(v) |
| 9.2 | `xi` | `compute.py` (`alpha1.json`) | Proposition 9.2.1 |
| 9.2 | `xi` | `mc.py` (`mc.json`, `surv.json`), `figs.py` | figures and tables no longer in the book |
| 9.2 | `xi` | `improve.py` (`improve.json`), `figs2.py` | Table 9.6; Figure 9.4 |
| 9.2 | `xi` | `rs1.py`, `rs1c.py` (`rs1.json`), `rs1d.py` (`rs1d.json`), `figs5.py` | Figure 9.5 |
| 9.2 | `process_map` | `bnd_alpha.py` (`bnd_0.6.json`, `bnd_0.75.json`), with `rs1.py` of `xi` | Table 9.7 |

## Chapter 10: The critical strip and Riesz's criterion

Scripts of the notes with the folder that holds each, the data files they write in brackets, and the tables and figures they produce.

| section | folder | scripts (data files) | tables and figures |
|---|---|---|---|
| 10.1 | `scale_mixture` | `comp.py` (`comp.json`) | Table 10.1 |
| 10.1 | `scale_mixture` | `strip.py`, `dens2.py` (`strip.json`), `improve.py` (`improve.json`), `figs2.py` | Table 10.2; Figure 10.1; Section 10.1.4 |
| 10.1 | `wiener` | `wiener.py` (`wiener.json`), `figs_wiener.py` | Table 10.5 |
| 10.1 | `wiener` | `salem_riesz.py` (`salem_riesz.json`), `riesz_explicit.py` (`riesz_explicit.json`), `figs_wiener.py` | Figure 10.3 |
| 10.1 | `riesz_program` | `explicit_bounds.py` (`explicit_bounds.json`) | the numbers in the proof of Theorem 10.1.3(3) |
| 10.2 | `riesz_program` | `twisted.py` (`twisted.json`), `sens.py` (`sens.json`), `figs_strategy.py` | Figure 10.4 |
| 10.2 | `riesz_program` | `riesz_large.py`, `check_large.py` (`check_large.json`), `check_large_twisted.py` (`check_large_twisted.json`) | the Riesz function up to $`x=10^{12}`$ in Section 10.2.1 |
| 10.2 | `riesz_program` | `riesz_thorin.py` (`riesz_thorin.json`), `figs_riesz.py` | Table 10.6; Figure 10.5 |
| 10.2 | `riesz_program` | `ggc_coef.py` (`ggc_coef.json`), `ggc_offset.py` (`ggc_offset.json`) | figures and tables no longer in the book |
| 10.2 | `riesz_program` | `killing.py`, `killing2.py` (`killing.json`) | Table 10.7 |
| 10.3 | `code` | `sqfree_constants.py` (`sqfree_constants.json`) | the constants of the section |
| 10.3 | `code` | `sqfree_sieve.py` (`sqfree_sieve.json`) | Tables 10.8, 10.9 |
| 10.3 | `code` | `sqfree_atoms.py` (`sqfree_atoms.json`) | Figure 10.6 |
| 10.3 | `code` | `sqfree_smooth.py` (`sqfree_smooth.json`) | Table 10.10 |
| 10.3 | `code` | `sqfree_gammanu.py` (`sqfree_gammanu.json`) | Table 10.11 |

## Chapter 11: Infinite divisibility, the eta analogue and almost generalised gamma convolutions

Scripts of the notes with the folder that holds each, the data files they write in brackets, and the tables and figures they produce.

| section | folder | scripts (data files) | tables and figures |
|---|---|---|---|
| 11.1 | `sato_bondesson` | `weil.py` | Table 11.1 |
| 11.1 | `sato_bondesson` | `comp.py` (`comp.json`), `figs.py` | Table 11.2; Figure 11.1 |
| 11.1 | `process_map` | `comp2.py` (`comp2.json`) | figures and tables no longer in the book |
| 11.1 | `process_map` | `ou.py` (`ou.json`), `figs.py` | Figure 11.2 |
| 11.2 | `code` | `wigfree_upsilon.py` (`wigfree_upsilon.json`) | Tables 11.5, 11.6 and the numbers of Section 11.2.3 |
| 11.2 | `code` | `wigfree_jacobi.py` (`wigfree_jacobi.json`) | Table 11.7 |
| 11.2 | `code` | `wigfree_gamma2.py` (`wigfree_gamma2.json`) | Table 11.8 |
| 11.2 | `code` | `wigfree_heat.py` (`wigfree_heat.json`) | the numbers of Section 11.2.6 |
| 11.2 | `code` | `wigfree_plots.py` | Figures 11.3 and 11.4 |
| 11.3 | `eta` | `zeros12.py` (`g12.npy`), `base.py` (`base.json`) | Table 11.9 |
| 11.3 | `eta` | `arith.py` (`arith.json`) | Tables 11.10, 11.11 |
| 11.3 | `eta` | `hit.py`, `hit2.py` (`hit.json`), `figs.py` | Figures 11.5, 11.6, 11.7 |
| 11.4 | `riesz_program` | `almost.py` (`almost.json`), `almost_dip.py` (`almost_dip.json`) | Table 11.14 |
| 11.4 | `riesz_program` | `explicit_bounds.py` (`explicit_bounds.json`) | the numbers in Theorem 11.4.4 and its proof |
| 11.4 | `riesz_program` | `more2.py` (`more.json`), `high.py` (`high.json`), `figs_almost.py` (`figs_almost.json`) | Table 11.13; Figures 11.8, 11.9 |
| 11.4 | `riesz_program` | `sens.py` (`sens.json`), `figs_strategy.py` | Table 11.15 |
| 11.4 | `riesz_program` | `killing2.py` (`killing.json`), `rmt_ratios.py` (`rmt_ratios.json`) | the exact ratios after a table no longer in the book |
| 11.4 | `riesz_program` | `ks.py` (`ks.json`) | Table 11.16 |
| 11.4 | `repro` | `ks_ak.py` (`ks_ak.json`) | the factors $`a_1`$, $`a_2`$, $`a_3`$ of Table 11.16 |
| 11.5 | `code` | `pakes_sizebias.py` (`pakes_sizebias.json`) | the numbers after Proposition 11.5.1 |
| 11.5 | `code` | `pakes_lambert.py` (`pakes_lambert.json`) | the numbers of Section 11.5.3 |
| 11.5 | `code` | `bdlp_ladder.py` (`bdlp_ladder.json`) | Theorem 11.5.8 and the scan after it |
| 11.5 | `code` | `branch_coeffs.py` (`branch_coeffs.json`) | Theorem 11.5.11 and the numbers after it |
| 11.5 | `code` | `heat_ladder.py` (`heat_ladder.json`) | Table 11.17 and the numbers of Section 11.5.4.3 |
| 11.5 | `code` | `pakes_filter.py` (`pakes_filter.json`) | Tables 11.18, 11.19; the infinite divisibility counts, the survey of discriminants, the family modulo $`7`$ and the values $`e^{\pi/q}\vartheta_a(1)`$ of Section 11.5.5.4 |
| 11.5 | `code` | `pakes_filter_cert.py` (`pakes_filter_cert.json`) | the certified entries of Tables 11.18, 11.19; Propositions 11.5.15, 11.5.18; the sums $`S_a`$ after Lemma 11.5.17 |
| 11.5 | `code` | `pakes_characters.py` (`pakes_characters.json`) | the checks after Propositions 11.5.13 and 11.5.15; the constants, zeros and functional equation of $`f_{\rm DH}`$; Proposition 11.5.19 |

## Chapter 12: Conclusions and open problems

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 12.2 | `phasefree.py` (`phasefree.json`), `phasefree_arb.py` (`phasefree_arb.json`) | the thresholds $`t_0(u)`$ of Section 12.2.6, and $`\mathcal G(u,6)>0`$ for $`u\ge10^{213}`$ |

## Chapter 13: Ferromagnetic random variables

The scripts of Chapter 13 are in the folder `ferro/`, each with the data file of its own name and the log of its run on the Hopper cluster; the exception is `xg_cauchy.py` of Section 13.4.6, which is in `code/` with `xg_cauchy.json` in `data/`. They work in the units of (1.2.3), and those that need the zeros read the first $`1700`$ ordinates, to $`70`$ digits, from a copy of `data/zeros_1700.txt` in the working folder. `fm_hermite.py` runs in two modes: with a list of values of $`K`$ it brackets the thresholds $`v^*(K)`$ of Table 13.5 and writes `fm_hermite.json`, and with the arguments `all 200` it checks $`f_K`$ itself for every $`K\le200`$ and writes `fm_hermite_all.json`. `fm_maxwell.py` takes the number of draws, $`10^7`$ for the book; `fm_maxwell_1e6.json` is an earlier run with $`10^6`$ draws. `fm_hcm.py` computes in ball arithmetic with python-flint (Johansson 2017), and its output is the certificate of Proposition 13.5.8, in the variable $`t/2`$. `fm_cutoff.py` also works in ball arithmetic, with the zeros isolated by Arb instead of read from the file, and so do `fm_qcert.py` and `fm_tvcert.py`: in the mode `zeros 800 9600`, `fm_qcert.py` writes $`\operatorname{Var}X`$ and the first $`800`$ ordinates, computed by Arb at $`9600`$ bits, to `fm_zeros_hp.txt` (not in the archive), which both scripts then read; its modes `table`, `all 200` and `minq` give the brackets of $`v^*(K)`$, the positivity of $`f_K`$ for every $`K\le200`$ and the minima of $`Q_{K,2\beta_K}`$. The table below lists the scripts and what they produce.

Scripts of Chapter 13, in the folder `ferro/` except `xg_cauchy.py`, which is in `code/`, the data files they write, and the tables and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 13.2 | `fm_interval.py` (`fm_interval.json`) | the zeros after Theorem 13.2.1 and after Corollary 13.2.4 |
|  | `fm_intcert.py` (`fm_intcert.json`) | the same, in ball arithmetic |
| 13.3 | `fm_polya.py` (`fm_polya.json`) | Tables 13.3, 13.4; the values in Sections 13.3.1 and 13.3.2 |
| 13.4 | `fm_polya.py` (`fm_polya.json`) | the zeros in Section 13.4.1; floating point values of the distances |
|  | `fm_cutoff.py` (`fm_cutoff.json`) | the inequality $`3s_{16}<\operatorname{Var}X<3s_{17}`$ in Theorem 13.4.3(i), in ball arithmetic |
|  | `fm_hermite.py` (`fm_hermite.json`, `fm_hermite_all.json`) | the grid search that starts the brackets of $`v^*(K)`$ |
|  | `fm_qcert.py` (`fm_qcert_table.json`, `fm_qcert_all.json`, `fm_qcert_minq.json`) | Table 13.5: $`f_K>0`$, $`v^*(K)`$ and $`\min Q_{K,2\beta_K}`$, in ball arithmetic |
|  | `fm_tvcert.py` (`fm_tvcert.json`) | the distances, variances and values in Sections 13.4.1 and 13.4.2 and Table 13.5, in ball arithmetic |
|  | `fm_maxwell.py` (`fm_maxwell.json`) | the simulation after Theorem 13.4.5 |
|  | `xg_cauchy.py` (`xg_cauchy.json`), in `code/` and `data/` | the table and the numbers after Theorem 13.4.7, with a numerical check of its proof |
| 13.5 | `fm_sizebias.py` (`fm_sizebias.json`) | Table 13.6; the checks of Theorems 13.5.3, 13.5.6 and 13.5.7 |
|  | `fm_hcm.py` (`fm_hcm.json`) | Proposition 13.5.8 |

## Chapter 14: The Pólya density as a maximum order statistic and the Thorin measure of the clock

Scripts in `code/` unless another folder is named, with the data files they write in brackets, and the tables and figures they produce.

| section | scripts (data files) | tables and figures |
|---|---|---|
| 14.2 to 14.16 | `maxwald_maxcdf.py` (`maxwald_maxcdf.json`), `maxwald_mc.py` (`maxwald_mc.json`), `maxwald_offline.py` (`maxwald_offline.json`), `maxwald_weak.py` (`maxwald_weak.json`), `maxwald_basepoint.py` (`maxwald_basepoint.json`), `maxwald_basepoint_fig.py`, `maxwald_FT.py` (`maxwald_FT.json`), `maxwald_thorin.py` (`maxwald_thorin.json`), `maxwald_finite.py` (`maxwald_finite.json`), `maxwald_tstar.py` (`maxwald_tstar.json`), `maxwald_bondesson.py` (`maxwald_bondesson.json`), `maxwald_lagarias.py` (`maxwald_lagarias.json`) | Tables 14.2, 14.3, 14.4, 14.5, 14.6, 14.7; Figures 14.1, 14.2, 14.3, 14.4; the bound of Proposition 14.7.2 and the orders of Theorem 14.13.2 |

## Chapter 15: The Thorin phase diagram

The scripts of Chapter 15 are in the folder `code/`, each with the data file of its own name in `data/` and the log of its run on the Hopper cluster. `phase_rig.py` and `phase_collide.py` compute in ball arithmetic with python-flint (Johansson 2017) by the method of the evaluator `polya_flint.py` of Chapter 2, the trapezoid rule on the real axis, with the weight $`e^{tu^2}`$ and an enclosure of the error of the rule; the other scripts that evaluate $`\xi_t`$ or $`H_t`$ do so in mpmath, by the trapezoid rule of Section 2.1.1 applied to the formula for $`H_t`$ of Section 8.4.1. The computations were first made for the working report (Polson and Sokolov 2026) and recomputed for the book. The table below lists the scripts and what they produce.

Scripts of Chapter 15, all in the folder `code/`, the data files they write in `data/`, and the tables and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 15.4 | `phase_thorin.py` (`phase_thorin.json`) | Table 15.2; Proposition 15.4.1; Figures 15.1 and 15.2 |
| 15.5 | `phase_cert.py` (`phase_cert.json`) | Table 15.3; the floating point march and its control |
|  | `phase_rig.py` (`phase_rig.json`) | Theorem 15.5.8, in ball arithmetic |
| 15.6.1 | `phase_weil.py` (`phase_weil.json`) | Tables 15.4 and 15.5; Figure 15.3 |
| 15.6.2 | `phase_gas.py` (`phase_gas.json`) | parts (i) and (ii) of Proposition 15.6.1; the data of Figure 15.4 |
|  | `phase_collide.py` (`phase_collide.json`) | part (iii) of Proposition 15.6.1, in ball arithmetic; Figure 15.4 |
| 15.6.3 | `phase_scaling.py` (`phase_scaling.json`) | Proposition 15.6.2; Figure 15.5 |

## Chapter 16: Pólya and the zeros of entire functions

The scripts of Chapter 16 are in the folder `code/`, each with the data file of the same name in `data/`; the heavy ones were run on the Hopper cluster of George Mason University. They compute in mpmath at 20 to about 300 digits and in numpy; `polya_clock.py` evaluates $`\xi`$ at height with python-flint (Johansson 2017), in ball arithmetic at 100 bits with the midpoints kept, and falls back to mpmath when python-flint is missing or fails its self test. The statements of the chapter were first computed by its author and recomputed for the book from the mathematics, by scripts written without reading, running or copying the author's; where the chapter printed a value, the data file records the value computed beside it. Nothing here is interval arithmetic: a certificate of the chapter is a floating point computation at the precision stated, with explicit analytic bounds. The table below lists the scripts and what they produce.

Scripts of Chapter 16, all in the folder `code/`, the data files they write in `data/`, and the tables and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 16.4 | `polya_coeffs.py` (`polya_coeffs.json`) | Proposition 16.4.1: the coefficients of $`\xi`$ at $`\frac12`$, the Turán ratios and the Jensen polynomials; Problem 16.4: the shifted Jensen polynomials and the zeros of $`\xi'`$ |
| 16.6 | `polya_clock.py` (`polya_clock.json`) | Table 16.1; Theorem 16.11.1, Proposition 16.11.10 and the counts of Theorem 16.11.15; the left panel of Figure 16.1 |
| 16.7 | `polya_xistar.py` (`polya_xistar.json`) | Propositions 16.7.1 and 16.7.2; the right panel of Figure 16.1 |
| 16.8 | `polya_mellin.py` (`polya_mellin.json`), `polya_mellin_z1.py` (`polya_mellin_z1.json`) | the zeros of the Mellin transform of the Riemann kernel and their count; its zero $`z_1`$ of Theorem 5.2.2, with $`\|M'(z_1)\|`$, the winding number on a circle about it and a second quadrature |
| 16.9 | `polya_bessel.py` (`polya_bessel.json`) | Table 16.2; the compensated condition of Theorem 16.11.7 for the Bessel model |
| 16.11 | `polya_density.py` (`polya_density.json`) | Theorems 16.11.5 and 16.11.6, Lemma 16.11.3, Proposition 16.11.14; the left panel of Figure 16.2 |
|  | `polya_primes.py` (`polya_primes.json`) | Proposition 16.11.8, Lemma 16.11.12, the numerical statements of Theorem 16.11.13, and the direct computations of $`D_\alpha`$ in Section 16.12 |
| 16.12 | `polya_qrh.py` (`polya_qrh.json`) | Lemma 16.12.2, Theorem 16.12.3, Table 16.3 and the threshold $`\alpha_*`$ |
|  | `polya_qrhpair.py` (`polya_qrhpair.json`) | the constants of Theorem 16.12.4 |
|  | `qrh.py` (`qrh_c.json`, part c) | the Gaussian smoothed Möbius sum $`M_{\mathrm g}`$ of the subsection of that name: the value at $`10^3`$ from the sum, from the lines $`\operatorname{Re}v=2`$ and $`0.6`$, the ratios to the leading term and the size of the first zero's residue (parts a, b, d are not quoted in the book) |
|  | `polya_fig.py` (`fig/polya_fig1.pdf`, `fig/polya_fig2.pdf`) | Figures 16.1 and 16.2, drawn from the data files above and from the ordinates of `data/zeros_1700.txt` |
|  | `pol_qrh_edge.py` (`fig/pol_qrh_edge.pdf`, `data/pol_qrh_edge.json`) | Figure 16.3, $`\operatorname{Re}\xi'/\xi(\alpha+i\tau)`$ at $`\alpha=7/8,1,3/2`$ on $`0\le\tau\le60`$, from mpmath with no zeros |

## Chapter 17: Yoshida's Hermitian forms, Weil positivity and the Thorin programme

The scripts of Chapter 17 are in the folder `code/`, each with the data file and the log of the same name in `data/`, and the figures are in `fig/`. They were written for the book from the mathematics, by editors who had not read, run or copied the scripts of the chapter's author; where the draft of the chapter printed a value, the data file records the value computed beside it, with the verdict that it agrees to the printed digits, differs, or cannot be reproduced. The ordinates are the first $`1700`$ zeros of `data/zeros_1700.txt`, to seventy digits (height $`2197.26`$); the gap lemma and the symbol also read the summary of the $`1800`$ zeros certified to height $`2304.3`$ in `riesz_program/data/v1min.json` and the constants of `data/dbeta_minimum.json`. Arithmetic is mpmath at 20 to 60 digits (more on the line $`\alpha=\frac12`$, where the least eigenvalues are below $`10^{-15}`$ from $`a=0.8`$ on), numpy in double precision where the text says so, and python-flint (Johansson 2017) ball arithmetic in `yo_gaps.py`, whose comparisons hold for whole balls; everything else is a floating point computation at the precision stated, with the analytic bounds of the text. The lower bound $`L_N`$ of Section 17.9 is a branch and bound in double precision, checked at its minimiser in mpmath at $`30`$ digits.

Each script is run from the repository root, for instance `timeout 60 python3 -I code/yo_tent.py`; `--quick` (every script except `yo_matrix.py`, which has only `--selftest`) runs a reduced case in a few seconds and writes to the temporary folder unless `--out` is given, except `yo_damped.py`, whose `--quick` takes about $`30`$ s and updates `data/yo_damped.json` in place, and the header of each script lists its checks in the order of the chapter. The table below gives in the last column the time in seconds of the full case in one process on a laptop, except for `yo_damped.py`, whose complete run was made with $`8`$ processes on the cluster. The figures are drawn by the scripts from the data they compute, in black and white. The extended precision values of `yo_damped.py` on the line $`\alpha=\frac12`$ are finite sections at the sizes of Table 17.2, and the section for $`a=4`$ ($`N=96`$) was computed with $`640`$ digits.

<!-- CHECK: optional reruns of the other scripts on Hopper are in checks-Y1.md (yo_explicit, yo_tent) and checks-Y4.md (yo_symbol); nothing else is pending. -->

Scripts of Chapter 17, all in the folder `code/`, the data files they write in `data/` (and the figures in `fig/`), the statements they check, and the time of the full case in seconds; each script also writes a log of the same name.

| script (files) | statements of the chapter it checks | seconds |
|---|---|---:|
| `yo_matrix.py` (`yo_matrix.json`) | The prime side matrix of Weil's damped form on a window, in the Fourier, Legendre and cosine and sine bases, even and odd sectors, in double or extended precision; imported by `yo_damped.py` and by the scripts of Chapters 18 and 20. Twelve self tests (`--selftest`): Weil's functional on a Gaussian against the zero side, the closed form matrix against an independent quadrature, the three bases against each other, the multiplier identity of Proposition 17.8.1 on trigonometric vectors, double precision against $`30`$ digits, and the floor row $`2b_\alpha`$ of Table 17.2 | 10 |
| `yo_explicit.py` (`yo_explicit.json`) | $`p_1`$ and $`p_2`$ from the Taylor coefficients of $`\log\xi`$ at $`\frac12`$, against the zeros and the sums $`\sum a_k^{-1}`$, $`\sum a_k^{-2}`$ of the book; the heat trace $`W(0.003)`$ from the primes and from the zeros; $`V_1`$ of the Gaussian by both forms of (17.1.2); $`V_1(\Delta_1)`$ by the closed form of Proposition 17.6.1 and by the residue series; Table 17.1 (prime side, zero side with $`1500`$ and $`1700`$ zeros, the tails, Landau's term, $`V_1(\Delta_x)`$ two ways); the coefficient $`0.18338`$ of Yoshida's matrix; five entries of the prime side matrix at $`a=1`$ against Theorem 17.7.1; $`\xi'/\xi(1)`$; the bounds of Theorem 17.8.5 on $`2\xi'(\alpha)/\xi(\alpha)`$ | 16 |
| `yo_tent.py` (`yo_tent.json`, `fig/yo_D.pdf`, `fig/yo_off.pdf`) | The tent coefficient: its least value on $`[0.05,16]`$ and its mean on $`[2,16]`$ from the primes by an exact Riesz mean; the first negative values of $`D`$ with a hypothetical pair, on a grid and in the continuum, and Figure 17.2; the constants of Proposition 17.9.1, the lower bound $`L_N`$ of $`D_{\rm on}`$ and the levels $`46.69`$, $`47.10`$ and $`47.48`$ of Section 17.9 | 7 |
| `yo_damped.py` (`yo_damped.json`, `fig/yo_lam.pdf`) | Table 17.2: the least eigenvalue of the finite sections on $`\|n\|\le N`$ for $`N=12`$ to $`384`$, the rate of convergence and the limits by Richardson extrapolation; the line $`\alpha=\frac12`$ in extended precision and the decay in $`a`$; the floor row $`2b_\alpha`$; Figure 17.3; the integral of the symbol against the prime side entry at $`\alpha=1`$, $`a=0.8`$ (Section 17.8); the printed values at $`\alpha=\frac12`$ of the tables of Chapters 18 and 20; the complete run (`--part all --procs 8`) took $`607`$ s on the cluster, $`571`$ s of them for the line $`\alpha=\frac12`$ | 607 |
| `yo_gaps.py` (`yo_gaps.json`) | Lemma 17.8.7: the gaps between the ordinates as balls ($`\gamma_2-\gamma_1=6.8873`$, the largest later gap $`5.4141`$); the explicit bounds on $`S(T)`$ of Hasanalizade, Shen and Wong (2022) quoted in Proposition 9.1.4, the thresholds in $`T`$ for $`h=8.3`$, $`6.8874`$ and $`5.36`$, and the step above $`X_0=3.061\times10^{10}`$ by monotonicity, with a grid in ball arithmetic as far as $`10^{3000}`$; the bound of Trudgian (2014) for comparison; the algebra and the margins in the proof of Theorem 17.8.5; it needs python-flint | 1.4 |
| `yo_ladder.py` (`yo_ladder.json`, `fig/yo_W.pdf`) | The stress test of Theorem 17.3.1 (Section 17.3): $`\sec\theta_*-1=1.1065\times10^{-4}`$, $`y_*=3660.05`$, the first failing order $`K_0=6349`$, the $`3000`$ orders above it and the $`395`$ orders $`6349+997j`$ that fail, the $`3000`$ orders below it that do not; the sign of $`T_k`$ for $`k=K_0-1`$ (positive on fine grids of $`x`$, with the real zeros dominating outside them) and for $`k=K_0`$ (negative), confirmed at one point in mpmath at $`60`$ digits; the zeros outside the window; $`K_0`$ for thirteen further pairs; Figure 17.1 | 11 |
| `yo_symbol.py` (`yo_symbol.json`, `fig/yo_symbol.pdf`) | The symbol $`v_\alpha`$ on $`[0,2000]`$ for $`\alpha=1,\,0.75,\,0.6,\,0.55`$ by the Hadamard sum over the $`1700`$ zeros with the Riemann von Mangoldt tail, against the direct value from $`\Gamma`$ and $`\zeta'/\zeta`$: its minima above $`\tau=20`$, the local minimum at $`17.49`$ and the minimum at $`\tau=0`$; the Laguerre quantity $`-(\log\Xi)''`$ and its minima; $`p_1`$, $`p_2`$, $`b_\alpha`$; the margins of the proof of Theorem 17.8.5; Figure 17.4 | 16 |

## Chapter 18: Zhu's compact windows and the Landau Widom horizon

The scripts of Chapter 18 are in the folder `code/`, each with the data file of the same name in `data/` and the log of its run (`zh_zeros.py` writes `data/zeros_8000.txt`, as said below); the figure is in `fig/`. The statements of the chapter were first computed by its author and recomputed for the book from the mathematics, by scripts written without reading, running or copying the author's; where the chapter printed a value, the data file records the value computed beside it, with the verdict that it agrees to the printed digits, differs, or cannot be reproduced. The ordinates are the first $`1700`$ zeros of `data/zeros_1700.txt`, to height $`2197.26`$, and above them the tails of the zero sums are bounded by the explicit Riemann von Mangoldt bound of Hasanalizade, Shen and Wong (2022) as the book uses it (`zh_comb.py`, `zh_product_bound.py`) or estimated from the density of the ordinates (`zh_window_eigs.py`). Arithmetic is mpmath at $`40`$ digits in `zh_comb.py`; double precision logarithms for the zero sums of `zh_product_bound.py`, with the optimum of every window recomputed at $`50`$ digits and the geometric side of the explicit formula at $`100`$ and $`140`$ digits (the two runs agree to $`15`$ digits or more); and $`60`$ digits in `zh_window_eigs.py`, which works at $`75`$ digits inside the module `yo_matrix.py` shared with Chapter 17. Raising the precision to $`90`$ digits changes the least eigenvalues of that script by a relative $`2\times10^{-69}`$ at most at $`a=0.5`$ and by $`8\times10^{-59}`$ at most at $`a=0.8`$, and $`1.5`$ times the number of Gauss Legendre nodes by $`3\times10^{-69}`$ and $`2\times10^{-58}`$, so that the quadrature and rounding error is some $`55`$ digits below anything printed in Table 18.1. The detection scripts use double precision Gram matrices, with python-flint (Johansson 2017) ball arithmetic for the weights, for the exact sign decisions (inertia additivity, precision doubling from $`192`$ to $`4096`$ bits) and for the ordinates of `zh_zeros.py`. Nothing here is interval arithmetic except those decisions; every other sign or number is a floating point computation at the precision stated, and a finite section value is an upper bound for the least eigenvalue and not a certificate. Every record in `data/` was made on one core of a laptop, except the finite sections of $`192`$ to $`400`$ modes and the damping ladder at $`64`$ and $`96`$ modes of `data/zh_window_eigs.json`, which were made on four cores of the Hopper cluster of George Mason University on 17 October 2026, and the files of the detection experiment (`zeros_8000.txt`, `zh_detect.json`, `zh_detect.log`) with the figure `fig/zh_horizon.pdf`, which come from a run on the same cluster on the same day: one job of `zh_zeros.py` and an array of $`81`$ tasks of `zh_detect.py` on four cores each. Table 18.3, Figure 18.1 and the statements of Section 18.6 that rest on them stand on that record, which uses the first $`2100`$ ordinates of `zeros_8000.txt` except where $`4000`$ or $`8000`$ are stated. The extrapolated values of Table 18.6 are fits to the finite sections of $`256`$, $`320`$ and $`400`$ modes. The first table below lists the scripts and what they produce, and the second their usage and run times.

<!-- CHECK: the optional rerun of code/zh_product_bound.py at its default settings (100 and 140 digits, 40 and 50 nodes per piece, all geometric cases) is PENDING on the cluster and replaces data/zh_product_bound.json and data/zh_product_bound.log, so keep the local copies until it has succeeded. Commands, from the repository root, after `ssh -O check hopper` reports "Master running": cp data/zh_product_bound.json /tmp/zh_product_bound_local.json; cp data/zh_product_bound.log /tmp/zh_product_bound_local.log; bash misc/tools/hopper_run.sh -c 1 -t 60 -n zh_product_bound -g zh_product_bound.json code/zh_product_bound.py data/zeros_1700.txt -- --geometric --out zh_product_bound.json > data/zh_product_bound.log; mv zh_product_bound.json data/zh_product_bound.json -->

Scripts of Chapter 18, all in the folder `code/`, the data files they write in `data/`, and the tables and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 18.5 | `zh_comb.py` (`zh_comb.json`) | Table 18.5: the $`42`$ values of $`A_a(\alpha)`$ and $`\log_{10}T_1`$, and the values of its caption ($`A_{0.8}(\frac12)`$, $`T_1`$ at $`a=0.8`$, $`1.6`$, $`2`$, and $`-2\zeta'/\zeta(\alpha)`$ at $`\alpha=1.1`$, $`1.5`$, $`2`$); the ratios of (18.5.2) to the exact sums; the heights $`2\pi e^{2a}`$, the counts $`N(T_*)`$ and the half length $`13.446`$ (printed $`13.4`$ in Section 18.6 and $`13.45`$ in Section 18.4); $`p_1=0.0231049931\ldots`$ by two routes and the sandwich of Theorem 17.8.5(iv); the first two columns of Table 18.3 |
| 18.3 | `zh_product_bound.py` (`zh_product_bound.json`) | the product table: the best bound, the number of atoms killed and the height of the last one, for $`a=0.5`$, $`0.8`$, $`1.0`$ (and $`0.35`$, $`1.2`$); the evaluation of each optimum on the geometric side of the explicit formula, with no hypothesis on the zeros |
| 18.4 | `zh_detect.py` (`zh_detect.json`, `fig/zh_horizon.pdf`), `zh_zeros.py` (`data/zeros_8000.txt`) | Table 18.3 and Figure 18.1; the tent levels $`x_*/2`$ and their ratios to the thresholds, the count preserving variant, the runs with $`4000`$ and $`8000`$ ordinates, with the margin $`90`$ and with bands of $`75`$ to $`600`$ modes, the refined thresholds and $`\delta=0.001`$, $`0.003`$ of Section 18.4 |
| 18.2, 18.5.1 | `zh_window_eigs.py` (`zh_window_eigs.json`; uses `yo_matrix.py`) | Table 18.1 and its convergence in the number of modes up to $`400`$; the second even eigenvalue at $`a=0.8`$; Zhu's reference values of his Table 1; the check of the normalisation of Proposition 18.2.1 by the zero side; the numbers after Theorem 18.5.2 and Table 18.6 |

Usage and run times of the scripts of the first table. With `--quick` a script writes to `/tmp` unless `--out` is given.

| usage (from the repository root) | run time |
|---|---|
| `python3 -I code/zh_comb.py [--quick] [--out FILE] [--zeros FILE]` | $`0.3`$ s on one core; `--quick` $`1`$ s |
| `python3 -I code/zh_product_bound.py [--quick] [--out FILE] [--zeros FILE] [--geometric] [--geom-cases CASES] [--Ls L1,...] [--geom-Ls L1,...] [--kmax K] [--mspan M] [--dps D] [--ngl N] [--no-extend]` | $`16`$ s for the stored case (zero side for five windows, geometric side at the optimum of each at $`100`$ digits and $`12`$ nodes per piece); `--quick` a few seconds |
| `python3 -I code/zh_window_eigs.py [--part PART] [--L L1,...] [--K N1,...] [--dps D] [--procs P] [--out FILE] [--quick]` | the record was made in fifteen local chunks of at most $`25`$ s on one core (listed in the header of the script) and two cluster jobs on four cores (the table part for $`192`$ to $`400`$ modes, $`593`$ s, and the damping part for $`64`$ and $`96`$ modes, $`136`$ s), with PART one of table, prec, zeros, damping, claims; `--quick` $`2`$ s |
| `python3 -I code/zh_detect.py [--selftest] [--quick] [--task K] [--list-tasks] [--gather] [--figure]` | `--selftest` $`2.5`$ s; `--quick` (the height $`51.37`$, $`a=0.9`$ to $`1.4`$) $`25`$ s; the $`81`$ tasks of `--task K` ran on the cluster, in $`3`$ to $`107`$ s for the $`24`$ main tasks and up to $`1637`$ s for the band of $`600`$ modes; python-flint is needed |
| `python3 -I code/zh_zeros.py [--count 8000] [--prec 192] [--digits 50] [--out FILE] [--check FILE] [--quick]` | about one minute for $`8000`$ zeros; python-flint is needed |

Four points of the computations are not in the chapter's tables. First, the search of `zh_product_bound.py` runs over $`K\le15`$ atoms and even exponents $`m`$ from $`2K+2`$ to $`2K+40`$, and the box is enlarged to $`K\le25`$ and $`m\le2K+80`$ where the optimum sits on its edge, which happens at $`a=1.2`$: the optimum there is $`K=16`$, $`m=46`$, with bound $`8.14\times10^{-21}`$, and the chapter prints no value at that window. The geometric value agrees with the zero side value to all digits shown for $`a\ge0.8`$ and to the central estimate of the tail beyond the $`1700`$ ordinates ($`3.5\times10^{-7}`$ relative) at $`a=0.5`$, which tests the normalisation of the explicit formula in the pole, archimedean and comb terms together. Second, the normalisation of Proposition 18.2.1 is checked by the zero side for the eight bottom eigenvectors of `zh_window_eigs.py` ($`a=0.5`$ and $`0.8`$, $`26`$ and $`40`$ modes, even and odd): the partial sum $`2\sum_{k\le K}F(\gamma_k)^2`$ over the first $`K`$ ordinates, plus a tail estimate from the density $`\frac1{2\pi}\log\frac t{2\pi}`$ of the ordinates, reproduces the eigenvalue to within $`0.16`$ percent of its value for $`K=300`$, $`400`$ and $`1700`$ in all eight cases, and to within $`2.1`$ percent for $`K=100`$ in the six cases with a tail estimate (no tail is estimated at $`a=0.5`$ with $`40`$ modes, where $`\gamma_{100}`$ lies below $`1.2`$ times the largest mode frequency). Third, the convergence in the number of modes is estimated by a three point power law $`\lambda_N=\lambda_\infty+CN^{-p}`$ through consecutive sizes; it is a fit and not a theorem, and the run of $`p`$ and $`\lambda_\infty`$ over the triples is written to the data file so that its drift can be judged. Fourth, in `zh_detect.py` the half length is scanned on a grid of step $`0.05`$ and refined by bisection to $`0.005`$ for each of four criteria: the float64 eigenvalue below $`-10^{-11}`$, the float64 eigenvalue below minus the Weyl bound of its rounding error (which certifies $`\lambda_{\min}<0`$), and the exact decisions $`\lambda_{\min}<0`$ and $`\lambda_{\min}<-10^{-11}`$. The exact decisions use that the part of the matrix from the on line zeros is positive definite once the number of zeros is at least the size of the parity sector, and that the off line quadruple is a perturbation of rank two in each parity sector, so that by Haynsworth's inertia additivity the sign is read from a $`2\times2`$ matrix.

## Chapter 19: Three families of generalised gamma convolutions and the deconvolution obstruction

The scripts of Chapter 19 are in the folder `code/`, each with the data file of the same name in `data/` and the log of its run, and the figures are in `fig/`. They compute in mpmath at $`20`$ to $`40`$ digits, except the sieve and the prime sums of `tf_primes.py`, which are double precision sums in numpy added with `math.fsum`, and the double precision Gaussian average of `tf_heat.py`, which is compared with the $`40`$ digit one. The zeros are the first $`1700`$ ordinates of `data/zeros_1700.txt`, to seventy digits, and above them the ordinates computed by `mpmath.zetazero`. The two heavier scripts, `tf_deltac.py` and `tf_window.py`, use eight processes on the Hopper cluster of George Mason University; the command for the cluster, with `misc/tools/hopper_run.sh`, is in the header of each script, and the others run on one core of a laptop in seconds. The statements of the chapter were first computed by its author and recomputed for the book from the mathematics, by scripts written without reading, running or copying the author's; where the chapter printed a value, the log records the value computed beside it, and each script prints its error controls: a recomputation at higher precision, the closed forms against their definitions, the imaginary parts that must vanish. Nothing here is interval arithmetic. The minima of Table 19.3 are minima over a grid of $`41`$ points with a local refinement in the script, not certified, and every other sign or number is a floating point computation at the precision stated. The first table below lists the scripts and what they produce, and the second their usage and run times.

Scripts of Chapter 19, all in the folder `code/`, the data files they write in `data/`, and the tables, figures and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 19.2 | `tf_primes.py` (`tf_primes.json`) | Table 19.1 with its slopes; Proposition 19.2.3 against its right side, built from the zeros of $`L(s,\chi_4)`$ below height $`62`$ ($`26`$ with positive ordinate), in all $`56`$ cases of $`\alpha`$, $`N`$ and $`b`$, with the largest difference, $`5.4\times10^{-14}`$ |
| 19.3 | `tf_deltac.py` (`tf_deltac.json`) | Table 19.2 in the three configurations, and the statements after it |
|  | `tf_window.py` (`tf_window.json`, `fig/tf_window.pdf`) | Table 19.3, the critical levels $`0.1893`$ and $`0.1895`$, and Figure 19.1 |
|  | `tf_visibility.py` (`tf_visibility.json`, `fig/tf_visibility.pdf`) | the lattice constants and heights; Proposition 19.3.2 with its closed forms against direct summation and the limits; Figure 19.2; Corollary 19.3.3 and Table 19.4 |
|  | `tf_xi.py` | the formulas for $`\xi'/\xi`$, $`\mathcal L_\Xi`$ and $`\mathcal L_0`$ shared by `tf_deltac.py` and `tf_window.py`, imported and never run |
| 19.4 | `tf_heat.py` (`tf_heat.json`) | the normalisation $`H_0(z)=\Xi(z/2)/8`$, the Gaussian average of Theorem 19.4.1(b) at $`z=0,10,28.27`$, the bound on $`\Phi_{\rm dBN}`$ beyond $`u=1.5`$, and the arithmetic of the bounds on $`\Lambda`$ |

Usage and run times of the scripts of the first table. With `--quick` a script writes to `/tmp` unless `--out` is given.

| usage (from the folder `code/`) | run time |
|---|---|
| `python3 tf_primes.py [--quick] [--out FILE]` | $`6`$ s on one core; `--quick` $`2`$ s |
| `python3 tf_visibility.py [--quick] [--out FILE] [--fig FILE] [--zeros FILE] [--fig-only]` | $`7`$ s on one core; `--fig-only` redraws the figure from the JSON file; `--quick` a few seconds |
| `python3 tf_heat.py [--quick] [--out FILE] [--dps D]` | $`2`$ s on one core |
| `python3 tf_deltac.py [--quick] [--out FILE] [--workers W] [--zeros FILE] [--dps D] [--heights N0,...] [--count 20]` | six rows on $`8`$ processes, $`47`$ s; `--quick`, the four rows $`n_0\le1000`$ in one process, $`4`$ s; one row at $`n_0=10^4`$ or $`10^5`$, $`19`$ s or $`23`$ s |
| `python3 tf_window.py [--quick] [--out FILE] [--workers W] [--zeros FILE] [--dps D] [--gaps N,...] [--deltas D,...] [--ks K,...] [--jstep J] [--no-figure]`, and `python3 tf_window.py --plot [--data FILE] [--fig FILE] [--root DIR]` | two gaps on $`8`$ processes, $`252`$ s; `--plot` $`2`$ s; `--quick`, the gap at height $`100`$, $`3`$ s |

Three points of the computations are not in the chapter's tables. The nodes and weights of the Gauss Legendre rule with $`200`$ nodes on $`[0,1.5]`$ and of the Gauss Hermite rule with $`80`$ nodes in `tf_heat.py` are polished by Newton's method at the working precision: those of `numpy.polynomial` (`leggauss` and `hermgauss`) have weights that differ from the polished ones by up to $`2.2\times10^{-11}`$ and $`2.7\times10^{-14}`$ in relative value, which gives an absolute error of $`2.7\times10^{-15}`$, against $`2.9\times10^{-18}`$ with accurate nodes, in the Gaussian average at $`z=28.27`$ in double precision. The bound $`\Phi_{\rm dBN}(u)\le4\pi^2\exp(9u-\pi e^{4u})`$ for $`u\ge\tfrac32`$, which decreases because $`9-4\pi e^{4u}<0`$, is proved in the script; its value at $`u=\tfrac32`$ is $`1.07\times10^{-543}`$, which bounds the kernel beyond the cut of the quadrature, and the kernel itself is $`5.36\times10^{-544}`$ there. In `tf_xi.py` the contribution $`\mathcal L_0`$ of the other zeros at an on line zero $`\tfrac12+i\gamma`$ is the closed form $`-1/s^2-1/(s-1)^2+\tfrac14\psi'(s/2)+\zeta'''/(3\zeta')-(\zeta''/\zeta')^2/4`$ at $`s=\tfrac12+i\gamma`$, which avoids the cancellation of a numerical average of $`\mathcal L_\Xi(\gamma\pm\tau)-\tau^{-2}`$; it agrees with that average at $`\tau=10^{-4}`$ to within $`8.3\times10^{-7}`$, the error being $`O(\tau^2)`$.

## Chapter 20: Deconvolution as a damped Weil criterion

The scripts of Chapter 20 are in the folder `code/`, each with the data file of the same name in `data/` and the log of its run. `dc_weil.py` computes in mpmath at 30 to 70 digits (the working precision of a part is 10 or 20 digits above the `dps` recorded in its data file) and in numpy; its parts W3 and W4 evaluate an integral over $`\theta`$ of $`\operatorname{Re}\xi'/\xi`$ against an explicit weight, by Gauss Legendre panels with a smooth taper. `dc_ladder.py` and `dc_rate.py` compute the least eigenvalue of the damped form on $`[-a,a]`$ from the matrices of the shared module `yo_matrix.py` of Chapter 17, used unchanged, in numpy at double precision on a Legendre and a Fourier basis, with mpmath for the control at thirty digits, for $`\xi'/\xi`$ and its Taylor coefficients at forty digits, and for the rung at $`\varepsilon=0`$ at $`a=0.75`$, $`0.8`$ and $`1`$ at sixty and eighty digits; `dc_rate.py` takes its sizes and the extrapolation from `dc_ladder.py`. The statements of the chapter were first computed by its author and recomputed for the book from the mathematics, by scripts written without reading, running or copying the author's. Nothing here is interval arithmetic: an eigenvalue is a variational upper bound in floating point, and a converged value is an extrapolation with an error estimate. The table below lists the scripts and what they produce.

Scripts of Chapter 20, all in the folder `code/`, the data files they write in `data/`, and the results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 20.1.1 | `dc_weil.py`, part W1 (`dc_weil.json`, `dc_weil.log`) | the Gaussian check of the second expression of Theorem 20.1.1 at widths $`\frac12`$ and $`1`$: the prime side against $`2\sum_{\gamma>0}\hat g(\gamma)`$ over the zeros of `data/zeros_1700.txt` |
|  | part W4 | the first identity of Theorem 20.1.1 for $`\hat\varphi`$ a triangle on $`[-2,2]`$, at $`\varepsilon_0=\frac12`$ and $`1`$, against the right side and against the zero sum |
| 20.1.2 | part W3 | formula (20.1.1) for $`h(u)=\cos(\pi u/2a)`$ at $`a=1`$, $`\varepsilon=\frac12`$ and $`\frac15`$: the prime side against the integral over $`\theta`$ |
| 20.2.1 | part W2 | the constants $`2+\gamma_E-\log4\pi`$, $`2b_{1/2+\varepsilon}`$ at $`\varepsilon=0.5`$, $`0.375`$, $`0.2`$, $`0.1`$, $`0.05`$, $`V_\xi`$, $`p_1`$, and the column $`\pi^2v_{1/2+\varepsilon}''(0)/4`$ of Table 20.3 |
| 20.2.1 | `dc_ladder.py`, parts grid and cross (`dc_ladder.json`, `dc_ladder.log`; the spot checks at $`N=1536`$ in `dc_ladder_spot.json`, `dc_ladder_spot.log`) | Table 20.1: the converged values for $`\varepsilon=0.5`$ to $`0.05`$ with their uncertainties, the Legendre and Fourier sizes, and the comparison with the converged values of Table 17.2 |
|  | parts y0 and y0x | Table 20.2, the entries of Table 20.1 at $`\varepsilon=0`$, the orders of magnitude before and after the prime $`2`$ enters, and the extended precision upper bounds at $`a=0.75`$, $`0.8`$ and $`1`$ |
|  | parts check and fig (`fig/dc_ladder.pdf`) | the double precision control against thirty digits and against $`1.5`$ times the nodes; Figure 20.1 |
|  | `dc_rate.py`, parts const, table and fit (`dc_rate.json`, `dc_rate.log`) | Table 20.3 and the ratios to its constants $`\pi^2v_{1/2+\varepsilon}''(0)/4`$, and the fit of Problem 20.2 |
|  | parts trial and lower | the cosine trial function and the numerical bounds of the proof of Theorem 20.2.1(b): Table 20.4 |

At $`\varepsilon=\frac15`$ the identity presupposes that no zero has $`\beta\ge0.7`$, which is unknown, so that run checks the formula and its normalisation, not the hypothesis. Parts W3 and W4 treat the Dirichlet terms $`n\le N`$ of $`\zeta'/\zeta`$, with $`N=100`$ or $`1000`$, in closed form against the weight (hybrid mode) and, with $`N=0`$, use no prime at all (pure mode, which needs the long windows). The usage is

```
python3 dc_weil.py [--quick] [--out FILE] [--zeros FILE] [--workers W]
    [--parts W1,W2,W3,W4] [--nzeros N] [--w3y 0.5,0.2] [--w4y 0.5,1]
```

and `--quick` uses small sizes: about 3 s for parts W1 and W2, 16 s for W3 and 14 s for W4 in one process, so the parts are run one at a time to keep each local run under 30 s, and the output goes to `/tmp` unless `--out` is given. The full sizes, with windows to $`\theta=600`$ for W3 and $`\theta=1500`$ for W4, eighteen nodes per panel and about $`108000`$ evaluations of $`\zeta'`$, were run on the cluster:

```
bash misc/tools/hopper_run.sh -c 16 -t 120 -g dc_weil.json -g dc_weil.log \
    code/dc_weil.py data/zeros_1700.txt -- --workers 16 --out dc_weil.json
```

which writes `dc_weil.json` and `dc_weil.log`, moved to `data/`. The output holds numbers as strings, with the parameters, per part, and a part written again replaces only its own entry.

The scripts `dc_ladder.py` and `dc_rate.py` are run from the root of the repository, and their usage is

```
timeout 60 python3 -I code/dc_ladder.py \
    --part check|grid|y0|y0x|cross|fig|all [--quick] [--out FILE] \
    [--fig FILE] [--procs P] [--y 0.5,0.2] [--a 0.25,1] \
    [--legendre 24,32] [--fourier 24,48] [--xa 0.75,1] [--xleg 24] \
    [--xfou 12,24] [--xdps 60]
timeout 60 python3 -I code/dc_rate.py \
    --part const|table|fit|trial|lower|all [--quick] [--out FILE] \
    [--procs P] [--y 0.5,0.2] [--a 1,2] [--afit 2,3,4] \
    [--atrial 1,2] [--alower 1,2] [--fourier 96,192] [--legendre 24,32]
```

and `--quick` runs a small case of every part in under 30 s in one process and writes to `/tmp` unless `--out` is given. With the default sizes each part takes under 20 s on one core of a laptop (the longest are the grid of Table 20.1, 16.5 s, and the fit of Problem 20.2, 19.1 s), but the tables of the chapter hold the larger runs, with the Fourier size $`1536`$ for the whole ladder grid, the Legendre sizes $`24`$ to $`96`$ at extended precision for $`a=0.75`$, $`0.8`$ and $`1`$, and $`a=7`$ in the fits ($`a=8`$ for the trial function), which were run on the cluster:

```
bash misc/tools/hopper_run.sh -c 16 -t 240 -g dc_ladder.json \
    -g dc_ladder.log -g dc_ladder.pdf code/dc_ladder.py \
    code/yo_matrix.py data/yo_damped.json -- --part all --procs 16 \
    --fourier 24,48,96,192,384,768,1536 --xleg 24,32,48,64,96 \
    --xfou 12,24,36,48,60 --out dc_ladder.json --fig dc_ladder.pdf
bash misc/tools/hopper_run.sh -c 16 -t 240 -m 48G -g dc_rate.json \
    -g dc_rate.log code/dc_rate.py code/dc_ladder.py code/yo_matrix.py \
    -- --part all --procs 16 --fourier 96,192,384,768,1536 \
    --afit 1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7 \
    --atrial 1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8 --out dc_rate.json
bash misc/tools/hopper_run.sh -c 16 -t 300 -g dc_ladder.json \
    -g dc_ladder.log code/dc_ladder.py code/yo_matrix.py \
    data/yo_damped.json misc/cluster-out/dc_ladder/dc_ladder.json \
    -- --part cross --procs 16 \
    --fourier 24,48,96,192,384,768,1536 --out dc_ladder.json
```

The first command ran on 16 cores until its limit of 240 minutes and was stopped inside part cross: the parts check, grid, y0 and y0x had finished (in $`3.1`$, $`10480.5`$, $`1826.7`$ and $`1888.1`$ s) and the figure part had not started. The second command ended with return code $`0`$ after $`10635`$ s ($`3269`$ s for the table of Table 20.3 and $`7353`$ s for the fit of Problem 20.2). The third command ran part cross alone, on the JSON of the first, and ended with return code $`0`$ after $`5455.6`$ s; it added the entries `cross` and `cross_summary` to that JSON and left every other entry unchanged. The three commands wrote `dc_ladder.json`, `dc_ladder.log`, `dc_rate.json` and `dc_rate.log`, moved to `data/` (the JSON of `dc_ladder.py` being that of the first command with the entries of the third merged into it, and its log that of the first command followed by that of the third), and `fig/dc_ladder.pdf` was then drawn by part fig from `data/dc_ladder.json`. Each data file holds numbers as strings, per part, and a part written again replaces only its own entries.

## Chapter 21: Thorin representations, Dirichlet $`L`$ functions and the Linnik constant

The scripts of Chapter 21 are in the folder `code/`, each with the data file of the same name in `data/` and the log of its run. They compute in mpmath (mpmath development team 2023) at $`25`$ to $`200`$ digits, and `gl_dirichlet.py` and `gl_mixture_mod5.py` also in ball arithmetic with python-flint (Johansson 2017), at $`110`$ to $`128`$ bits (the runs used python-flint 0.9.0 and mpmath 1.4.1, except the cluster run of `gl_sinsq_chi4.py`, which needs no python-flint and used mpmath 1.3.0). The zeros of $`\zeta`$ are the first $`1700`$ ordinates of `data/zeros_1700.txt`, to seventy digits. The statements of the chapter were first computed by its author and recomputed for the book from the mathematics, by scripts written without reading, running or copying the author's; where the chapter printed a value, the log records the value computed beside it and says whether the two agree, and the text of the chapter gives the computed value. Nothing here is interval arithmetic except the certificates named in the chapter: the counts of the zeros of $`L(s,\chi)`$ below height $`300`$ for $`\chi_4`$ and $`\chi_5`$ and of those of $`L(s,\psi)`$ below height $`12`$, the enclosures of the ordinates by two certified signs, the positivity of $`L(\sigma,\chi)`$ on $`[\tfrac12,2]`$, the brackets and the number of zeros of the mixture modulo $`5`$ in its two intervals, and the enclosures of $`L(\tfrac12,\chi_5)`$ and $`L(\tfrac12,\psi)`$ away from $`0`$. Every other sign or number is a floating point computation at the precision stated. The first table below lists the scripts and what they produce, and the second their usage and run times.

The script `gl_pole_lock.py` (data file `gl_pole_lock.json`, run as `python3 gl_pole_lock.py [--N N] [--out FILE]` from `code/`) computes the ratios $`0.402`$, $`0.510`$, $`0.645`$ and $`0.735`$ of the prime sum to $`e^{t/4}`$ at $`t=0.5`$, $`1`$, $`2`$ and $`3`$ in the pole lock heuristic of Section 21.6.3, with $`\Phi(\sqrt{t/2})`$, the part below the first prime and the remainder; it takes $`1`$ s on one core for the default $`N=6\times10^{6}`$.

All the runs so far were local, in one process, except the full runs of `gl_zeta_moments.py` and `gl_sinsq_chi4.py`, which were made on the Hopper cluster of George Mason University. The values of Tables 21.1 and 21.2 and of the paragraph after the first, the tilt by one quarter and the values of the transfer function in Section 21.3.3, the negative moments and the tilt $`u_0=1`$ come from the first of them, at $`100`$ digits (`data/gl_zeta_moments.json` and `data/gl_zeta_moments.log`), and so do the statistics over the cuts $`N=50`$ to $`1699`$ quoted in Section 21.3.3; the same run also evaluates the continuation of the negative moments at heights up to $`80`$ on $`\operatorname{Re}s=-2`$ at $`200`$ digits. The numerical check for $`\chi_4`$ in Section 21.5.6 comes from the second of them (`data/gl_sinsq_chi4.json` and `data/gl_sinsq_chi4.log`), in which the script finds the zeros of $`L(s,\chi_4)`$ itself.

<!-- CHECK: code/gl_dirichlet.py and code/gl_mixture_mod5.py ran locally within 30 seconds each (full size, python-flint 0.9.0) and are complete; the cluster run would only repeat numbers that are already in data/gl_dirichlet.json and data/gl_mixture_mod5.json. The part central of gl_mixture_mod5.py (L(1/2, chi_5), L(1/2, psi), L(1/2, psibar) at 128 bits, finding M5 of misc/reviews/nick-linnik/final-X21.md) was added and run locally on 17 October 2026 (/tmp/zh_venv python-flint 0.9.0, under 0.1 s) and is stored in data/gl_mixture_mod5.json under the key central. Optional provenance rerun, from the repository root after ssh hopper: bash misc/tools/hopper_run.sh -c 1 -t 90 -g gl_dirichlet.json -g gl_dirichlet.log code/gl_dirichlet.py -- --out gl_dirichlet.json and bash misc/tools/hopper_run.sh -c 1 -t 60 -g gl_mixture_mod5.json -g gl_mixture_mod5.log code/gl_mixture_mod5.py -- --out gl_mixture_mod5.json (misc/reviews/nick-linnik/checks-L2.md Section 5). -->

Scripts of Chapter 21, all in the folder `code/`, the data files they write in `data/` (each with a log of the same name), and the tables and results they produce.

| section | scripts (data files) | tables and results |
|---|---|---|
| 21.3.3 | `gl_zeta_moments.py` (`gl_zeta_moments.json`) | Tables 21.1 and 21.2; the classical constant and the tilt by $`\tfrac14`$ of Proposition 21.3.5, with equation (21.3.11); the values of the transfer function in Proposition 21.3.6; the negative moments of Proposition 21.3.10 and the residue ratio of Proposition 21.3.11; the tilt $`u_0=1`$ of Section 21.6.1 |
| 21.3.3.5 | `gl_density_pf.py` (`gl_density_pf.json`) | the density of $`Y`$ and its values at $`0`$ and $`2`$, the log concavity grid and the minors of Table 21.3; the Stirling numbers after Lemma 21.4.1; the smooth model and the elasticities of Problem 21.4 |
| 21.6.1 | `gl_dirichlet.py`, parts `zeros`, `real`, `classical`, `table`, `transfer`, `negmom` (`gl_dirichlet.json`) | the certified counts $`203`$ and $`214`$ and the ordinates; the grid and the certificate for $`L(\sigma,\chi)>0`$; Table 21.5 with the fluctuation $`S(300,\chi)`$; the transfer function and the Hadamard product at $`s=3`$; the negative moments, the classical constants and the series of equation (21.3.11) for $`\chi`$ |
| 21.5 | `gl_dirichlet.py`, parts `complex`, `dedekind`, `levy` (`gl_dirichlet.json`) | Proposition 21.5.2 (the ratio $`\xi(\tfrac12+s,\psi)/\xi(\tfrac12-s,\psi)`$ at $`s=0.7`$ for $`\psi(2)=i`$); Proposition 21.5.3 (the ratio $`\xi_F/(\xi\,\xi(\cdot,\chi))`$ with $`\zeta_F`$ computed independently, and the first power sum for $`\mathbb Q(i)`$); Proposition 21.5.4 (equation (21.5.4) at three pairs $`(\alpha,s)`$) |
| 21.6.2 | `gl_mixture_mod5.py` (`gl_mixture_mod5.json`) | the ordinates of the three $`L`$ functions below height $`12`$ with certified counts, the counts of atoms below the ends of the two intervals, the central values $`L(\tfrac12,\chi_5)`$ and $`L(\tfrac12,\psi)`$ (part `central`), Table 21.6 with the certificates of one zero per interval, and the sizes near the poles |
| 21.5.6 | `gl_sinsq_chi4.py` (`gl_sinsq_chi4.json`) | the numerical check for $`\chi_4`$ (the kernel identity, the prime side against the direct value, the zero side, the minima over $`b`$, the peaks), and the differences of the plain prime sum quoted in Problem 21.7 |

Usage and run times of the scripts of the first table. With `--quick` a script writes to `/tmp` unless `--out` is given. `gl_dirichlet.py` and `gl_mixture_mod5.py` write their data file to the current folder by default, so from `code/` they are run with `--out ../data/gl_dirichlet.json` and `--out ../data/gl_mixture_mod5.json`. For `gl_sinsq_chi4.py` the zeros of $`L(s,\chi_4)`$ are computed by default, and `--zeros` reads them from a file such as `data/dbeta_zeros.json`.

| usage (from the folder `code/`) | digits and run time |
|---|---|
| `python3 gl_zeta_moments.py [--quick] [--out FILE] [--zeros FILE] [--dps D] [--line-dps D] [--taus LIST]` | `--quick`: $`50`$ digits, $`9`$ s on one core; the full run ($`100`$ digits, every cut $`N`$ from $`50`$ to $`1699`$, $`200`$ digits on $`\operatorname{Re}s=-2`$) took $`671`$ s on one core of the cluster |
| `python3 gl_density_pf.py [--quick] [--out FILE] [--zeros FILE] [--dps D]` | $`70`$ digits, trapezoid steps $`\tfrac14`$ and $`\tfrac18`$, $`8`$ s on one core; `--quick`: $`40`$ digits |
| `python3 gl_dirichlet.py [--quick] [--out FILE] [--parts P1,P2,...] [--T 300] [--h 0.03125] [--P 10000000] [--N 2000000]` | full size ($`T=300`$, prime cutoff $`2\times10^{6}`$, ideal cutoff $`10^{7}`$): $`15.8`$ s for the parts zeros, real and classical and $`11.5`$ s for the other six, in two runs; `--quick`: $`12.8`$ s |
| `python3 gl_mixture_mod5.py [--quick] [--out FILE] [--parts zeros,central,mixture] [--lambdas L1,L2,...] [--no-rect]` | $`8.9`$ s with the rectangle certificates; `--quick`: $`0.6`$ s |
| `python3 gl_sinsq_chi4.py [--quick] [--out FILE] [--parts 1,2,3,4,5] [--nmax N] [--tmax T] [--zeros compute\|book\|FILE] [--bmax B] [--nq N] [--workers W] [--merge]` | `--quick`: $`4`$ s; the full run ($`203`$ zeros below height $`300`$, $`b\le150`$, prime powers up to $`2\times10^{6}`$) took $`162`$ s on $`8`$ processes of the cluster |

Seven points of the computations are not in the chapter's tables. First, $`L(s,\chi)`$ is evaluated from the Hurwitz zeta function and exact character tables, checked against the character tables of python-flint and against its `dirichlet_l` at $`\tfrac12+17.25i`$; the real function $`\xi(\tfrac12+it,\chi)`$ is scanned for sign changes on a grid of step $`\tfrac1{32}`$, and each is refined by Illinois iteration at $`110`$ bits and enclosed by two certified signs. The count below height $`300`$ is certified by the argument principle on the half contour $`\tfrac12\to2\to2+iT\to\tfrac12+iT`$, in adaptive boxes on each of which the ball enclosure of $`L`$ excludes $`0`$ and has radius below half its modulus; its equality with the number of sign changes shows that the zeros are simple and on the line. Second, the tails are the Riemann von Mangoldt density as the chapter states it, the smooth density $`\vartheta'(t)/\pi`$ (in closed form from its asymptotic series, checked once against quadrature of the digamma function) and the anchored tail of equation (21.3.8); sums of other functions of $`\gamma^2`$, such as the tilts and the transfer function, use the exact series in $`\gamma^{-2}`$, so that no quadrature enters the zero side. Third, the Taylor coefficients of $`\xi(\tfrac12+s)`$ are trapezoid sums on circles of radius $`3`$ and $`2`$, whose agreement is the error control, and for the characters Cauchy integrals on circles of radius $`2`$ and $`3`$ with $`96`$ points at $`60`$ digits; half integer radii are avoided because they meet the poles of the gamma factor, which cancel against the trivial zeros of $`L`$. Fourth, the two forms of equation (21.3.14) are one integral after $`u=w^2`$, so their difference measures the quadrature at the singular endpoint, and the value used is that of the substitution $`u=v^4`$, which has a smooth integrand.

Fifth, in equation (21.5.4) the expression $`e^{-sx}-1+sx`$ is replaced by its Taylor series when $`sx<0.01`$, because the quadrature otherwise samples points near $`0`$ where the expression cancels; the Binet integrals agree with their closed forms to $`3\times10^{-31}`$, the sums over prime ideals run to norm $`30000`$ and $`2\times10^{6}`$, and $`\zeta_F`$ is computed independently of $`\zeta\,L(s,\chi)`$, by the Euler product over the prime ideals and by an approximate functional equation built from the ideal counts and the gamma factor of $`F`$ alone, so that its comparison with $`\xi\,\xi(\cdot,\chi)`$ tests the gamma factors and the constants. Sixth, the zeros of the mixture modulo $`5`$ in an interval free of poles are the zeros of the entire function $`\lambda c_AD_B+(1-\lambda)c_BD_A`$, where $`D_A=\xi(\tfrac12+iy,\chi_5)`$, $`D_B=\xi_{\psi\bar\psi}(\tfrac12+iy)`$ and $`c_A=D_A(0)`$, $`c_B=D_B(0)`$, which stays of size one at the poles; it is evaluated at $`128`$ bits on a grid with geometric offsets from $`10^{-9}`$ to one half of the interval length from both ends, each sign change is a certified bracket checked to contain no pole and refined by bisection to width $`10^{-22}`$, and the argument principle on the rectangle with the interval as base and half height $`\tfrac14`$ counts one zero. Seventh, the oscillatory integral of the kernel identity (21.5.7) is done over half periods plus the asymptotic expansion of the tail, because the quadrature of mpmath for oscillatory integrals is wrong by $`8\times10^{-5}`$ on this integrand, and the integral over $`b`$ of each prime term is the kernel identity itself, so that the integral of the truncated density is known term by term.

## Chapter 22: Friable approximants of Xi after Miyatake and the cosh ladder

The scripts of Chapter 22 are in the folder `code/`, each with the data files of the same name in `data/`. They were written for the book from the mathematics, without reading or running the author's scripts, and they compute in mpmath at the precision stated in the chapter. The first table lists the scripts of the chapter's first half, the second those of the second half.

| script | data | produces | time |
|---|---|---|---|
| `code/mi_friable.py` | `data/mi_friable_*.json` | approximants, kinks, Lorentz and curvature bounds, shifts (tasks `approx`, `kink`, `lorentz`, `curv`, `shifts`) | under 30 s per task |
| `code/mi_humps.py` | `data/mi_humps_*.json`, `data/mi_check215_0.json` to `data/mi_check215_7.json` | hump count, humps, scans, 215 digit check | parts of under 30 s on the laptop |
| `code/mi_survey.py` | `data/mi_survey.json` | survey table | 4 s |
| `code/mi_depart.py` | `data/mi_depart.json` | departed zeros, counts, chain (in 11 segments) | the chain runs on the laptop in segments |
| `code/mi_miyatake.py` | `data/mi_miyatake_{norm,kink,tail,xi2,xi2trap,cusp,cscan_a,cscan_b,cscan_big,interp}.json` | identity, R, Q, kinks, count for xi^{(2)}, tail constant, interpolation, C_b, cusp ratios, b_0 | each task under 30 s on the laptop (interp 35 s) |
| `code/mi_gb.py` | `data/mi_gb_*.json` | g_b, clocks, zero counts, thresholds, g^{(2)}, g_{9/2}, escape law | each task under 30 s |
| `code/mi_ladder.py` | `data/mi_ladder_{path,density,arith}.json`, `data/mi_ladder_c1_counts.json`, `data/mi_ladder_c2_counts.json` | Thorin path, interlacing, phase count, arithmetic of the zeta bound | under 30 s |
| `code/mi_dh.py` | `data/mi_dh_*.json` | Davenport Heilbronn definition, zeros, ladder, order, counts, Table thr, heat flow | each call under 30 s (the variables `MI_*` select subsets) |

## Chapter 23: Haseo Ki and the Thorin programme

The scripts of Chapter 23 are in the folder `code/`, with their data files in `data/` and their logs beside them.

| script | data | produces | time |
|---|---|---|---|
| `code/ki_tilt.py` | `data/ki_tilt_book.json`, `data/ki_tilt_chapter.json` (and `.log`) | moment tilt of Phi by u^(2n): cumulants, six zeros, Edgeworth shifts, Laplace error (run with `--norm book` or `--norm chapter`, `--out FILE`) | about 2 s per n (mpmath at 40 and 50 digits) |
| `code/ki_check.py` | `data/ki_check.log` | maximum of x^2/(2 cosh x), saddle formulas, percent claims, normalisation of Phi | seconds |

## Reruns and gaps

In September 2026 every script of the archive was run again from a fresh copy on the Hopper cluster of George Mason University, and every data file it wrote was compared with the archived one; the reports are `compare_report.md`, `compare_part6_report.md` and `compare_riesz_report.md` in `repro/`. Apart from the running times that four of the files record, the rerun reproduced every archived data file that a script writes to a relative $`10^{-9}`$, most of them bit for bit, and the report of the enclosures of Chapters 6 and 7 byte for byte. There are two exceptions: in `selberg.json`, $`31`$ of $`4020`$ values differ by up to $`2.1\times10^{-8}`$ relative, and in `improve.json` of `scale_mixture` the densities agree to $`2.9\times10^{-8}`$ in absolute terms but not in relative terms in the far tail, where they are below $`10^{-7}`$ and no printed number depends on them. The enclosure program was corrected afterwards to bound the remainders of the sums in E3 and E5 and to compare every bound, rounded outward, with the number used in the text; its report replaced the archived one. The scripts changed since then and those added for the book were run on the cluster too, and the data files in the archive are their output or agree with it. On 28 September, after the last of these changes, the driver behind `compare_report.md`, which by then ran every script of `code/`, the Bosch script and the five programs of `verification_reciprocal_zeta/` outside `torus/`, was run once more in an empty folder. All $`65`$ runs ended without error, and every data file that a script writes agreed with the archived one to a relative $`10^{-9}`$, apart from the running times that six of the files record, one host name and the same $`31`$ values of `selberg.json`; the report of the enclosures and the output of `spot_checks.py` were reproduced byte for byte (`repro/compare_report_0928.md`). The scripts `vdpair.py`, `idsd_height.py`, `threezones.py` and `cm_orders.py`, added on 29 September after that run, were run on the cluster on their own, `vdpair.json`, `idsd_height.json`, `fig/threezones.pdf` and `cm_orders.json` are their output, and the driver now includes them. The scripts of Chapters 14 and 15, added in October 2026 after that run, were run on the cluster on their own, and the data files of the archive are their output. So were the scripts of Chapter 16, added later in October, except `polya_qrhpair.py` and `polya_fig.py`, which take seconds and were run on a laptop.

The rerun also showed what the archive lacked, and each item was followed up on the cluster in a folder of `repro/reconstruct/` with its scripts, its logs and a report. `primeside.py` imports a module `polya` that was not in the archive; with the reconstructed `polya.py`, its `results_prime.json` agrees with the archived file in all $`162`$ values to $`7.4\times10^{-16}`$ relative, $`148`$ of them bit for bit, and `plots.py` run on it regenerates `extra.json` and `sens2.json` byte for byte. `ttest.py` and `ttest2.py` stopped because they were given no arguments; the driver now passes the command lines given under Chapter 2, and with them `ttest2.py` reproduces the node counts and the errors of Table 2.2, while its times per point vary with the machine. Two lines of the Bosch script executed a file `gaptest.py` that is not in the archive, which defined the functions of block 12 again and an empty list `out`; they are replaced by `out=[]`, and the script runs as it stands; the part of Chapter 4 on the Bosch script says what it reproduces. `spot_checks.py` writes no file, so the rerun had nothing to compare; run as described under Chapter 6, it reproduces the report in the archive byte for byte, and the $`112`$ numbers and statements of Chapters 6 and 7, Section 5.1 and this file that the report bears on agree with it at the printed precision. Of the five data files that no script wrote, `hankel.json`, `spacing.json` and `sinh_bench.json` are reproduced by the reconstructions described under Figures and reconstructed scripts above and `riesz_explicit.json` by `riesz_explicit.py` of `basepoint_one/wiener`; the fifth is `results_main.json`.

Four gaps remain in these items. The printed output of the Bosch script is not the record of a single run: output lines 37 to 39 and 109 come from code that is not in the archive, and lines 93 to 102 from an earlier version of block 13; reconstructions in `repro/reconstruct/bosch_v7` reproduce lines 37 to 39 and 109 and the density at $`t=0.003`$ at $`150`$ digits, but the output file is kept as the authors made it. `primeside.py` runs only with the reconstructed `polya.py`, and its Pólya column then differs from the archived one by up to four units in the last place, a difference that none of $`196`$ tested conventions removes and that no printed number shows. `spot_checks.py` still prints, beside the values it computes, three reference values from an earlier version of Chapters 6 and 7, and its docstring gives two minutes for a run that takes about eight; the values it computes agree with the text. Finally, `data/results_main.json` has no script that writes or reads it. A reconstruction in `repro/reconstruct/primeside` reproduces $`117`$ of its $`135`$ numbers, but not its distances `tv` and `sup`, which rest on a tail model that cannot be recovered, nor its timings; the book takes from it only $`\gamma_{1700}=2197.2587`$, which the zero lists also give.

Some computations of the book have no script in the archive. In Chapters 1 to 5 and 8 these are Table 2.1, whose values `plots3.py` holds as constants, Table 8.1, and Figures 2.2 and 3.3, most of whose data files exist although the code that drew them does not; the intermediate values of Section 2.3, which come from calls of `polya_flint.py` at its defaults; the check with the step $`h`$ halved and the time per evaluation in Section 2.1.1; the double precision run with $`\delta=0.15`$ of Section 2.1.4; and the numerical checks in the text of Section 4.2, made with `mpmath` as that section says. The Horner experiment of Section 2.1.5 is reproduced by `horner100.py` in `repro/reconstruct/ttest2`, outside `code/`; `lehmer.py` prints the timings of Section 2.2.3.1 but does not save them; and `kmono.json` holds only minima that underflow to zero, while `kmono.py` draws Figure 5.2 from the values it computes. In Chapters 6 and 7 no script reproduces Tables 7.1 and 7.6, the zeros of the truncations in Table 7.2, or the seven figures in `rzfig/`. Other computed statements of Section 7.1 have no script either: the grid evaluation of Proposition 7.1.5, with the mode of $`f`$ and the agreement of Talbot and de Hoog inversion quoted after it, the values of $`w_\sigma`$ in Section 7.1.2.3 and the data of Proposition 7.1.18. In Chapters 9 to 11 no script produces the column of Table 9.2 that sums over the first $`1700`$ zeros, `ou.json` does not store the batch standard errors of the autocorrelations quoted in Section 11.1, and `comp.py` of `process_map` prints the sample mean $`0.6664`$ of $`S_2`$ given with Table 9.1 but does not store it in `comp.json`.

## References

* Hasanalizade, E., Shen, Q. and Wong, P.-J. (2022). Counting zeros of the Riemann zeta function. Journal of Number Theory 235, pages 219 to 241.
* Johansson, F. (2017). Arb: efficient arbitrary-precision midpoint-radius interval arithmetic. IEEE Transactions on Computers 66, pages 1281 to 1292.
* The mpmath development team (2023). mpmath: a Python library for arbitrary-precision floating-point arithmetic, version 1.3.0.
* Polson, N. G. and Sokolov, V. (2026). Four attacks on the computational path: a working report. Manuscript, 7 October 2026.
* Talbot, A. (1979). The accurate numerical inversion of Laplace transforms. Journal of the Institute of Mathematics and its Applications 23, pages 97 to 120.
* Trudgian, T. S. (2014). An improved upper bound for the argument of the Riemann zeta-function on the critical line II. Journal of Number Theory 134, pages 280 to 292.
