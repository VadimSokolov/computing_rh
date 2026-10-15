# Appendix A. Code and reproducibility

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), appendix A. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23400373. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This appendix says how the numbers, tables and figures of the book were computed and how to check them. The archive of scripts and data files it describes is in the book's repository https://github.com/VadimSokolov/computing_rh. It prints the two Pólya evaluators of the Xi function, `polya_flint.py` on the real axis and `tilted.py` on a tilted line, which apply the trapezoid rule to the Pólya kernel in ball arithmetic (FLINT and Arb through `python-flint`), and states what their ball radii certify: the rounding errors of the finite sums, not the trapezoid error, the cut of the axis or the truncation of the theta series. Tables list, in book order, the scripts and data files behind the tables and figures of the chapters. The last section reports the September 2026 rerun of every script on the Hopper cluster of George Mason University, which reproduced every archived data file that a script writes to a relative $10^{-9}$ apart from running times and two stated exceptions, and lists the computations that have no script and the gaps that remain.

## From the literature

The appendix presents no result of the literature; the scripts use Arb (Johansson 2017), `mpmath` (mpmath development team 2023) and Talbot's inversion of Laplace transforms (Talbot 1979).

## New results

This appendix states no theorem or proposition.

<!-- results:begin -->
No new results of this appendix are catalogued.
<!-- results:end -->

## Key terms

reproducibility, code archive, book repository, `computing_rh`, `python-flint`, FLINT, Arb, ball arithmetic, certified enclosure, mpmath, Pólya evaluator, Polya evaluator, trapezoid rule, Xi function, Talbot inversion, Hopper cluster, data files
