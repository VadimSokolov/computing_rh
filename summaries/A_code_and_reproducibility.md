# Appendix A. Code and reproducibility

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), appendix A. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23595861. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This appendix says where the code and the data of the book are and how its two Pólya evaluators work. The scripts and data files are not printed: they are released with the PDF in the book's repository https://github.com/VadimSokolov/computing_rh, and the file `code/README.md` of the repository holds the inventory of the scripts, their usage and run time tables, the reruns made for the book and the computations that have no script. The appendix itself has one section, on the two Pólya evaluators of the Xi function, `polya_flint.py` on the real axis and `tilted.py` on a tilted line, which apply the trapezoid rule to the Pólya kernel in ball arithmetic (FLINT and Arb through `python-flint`, Python 3.11). It gives their default parameters and those of the main runs, and states what their ball radii certify: the rounding errors of the finite sums, and not the trapezoid error, the cut of the axis or the truncation of the theta series.

## From the literature

The appendix presents no result of the literature; the evaluators use Arb (Johansson 2017) and `mpmath` (mpmath development team 2023) for reference values.

## New results

This appendix states no theorem or proposition.

<!-- results:begin -->
No new results of this appendix are catalogued.
<!-- results:end -->

## Key terms

reproducibility, code archive, book repository, `computing_rh`, `python-flint`, FLINT, Arb, ball arithmetic, certified enclosure, mpmath, Pólya evaluator, Polya evaluator, trapezoid rule, Xi function, tilted evaluator, Hopper cluster, data files
