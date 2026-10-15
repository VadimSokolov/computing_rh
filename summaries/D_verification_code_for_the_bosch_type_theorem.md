# Appendix D. Verification code for the Bosch type theorem

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), appendix D. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23400373. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This appendix documents the numerical checks behind the Bosch type theorem of Chapter 4: an elasticity criterion that reproves Bosch's theorem and characterizes, in terms of the angles of the zeros, for which $p$ the reciprocal Xi function at $x^p$ is hyperbolically completely monotone (HCM). Its checks are one Python script, `code_bosch/verify_bosch_xi_v7.py`, with its output `verify_bosch_xi_v7_output.txt`, in the book's archive; it has thirteen blocks using mpmath and numpy and runs in about 400 seconds. A table gives, for each block, its lines in the script and the output, what it checks and the statements it supports: HCM tests for Bosch's family and a toy example, the Thorin density of the clock at $s^p$, the residue series against Talbot inversion, the prime bound, the condition on the first gap between zeros and the constants for $T_0=3\times10^{12}$. The checks confirm formulas and constants and prove nothing. The output file is not the record of one run: 14 of its 109 lines come from code that is not in the script.

## From the literature

The checks concern Bosch's family of densities, which contains the half Cauchy density (Bosch 2015), and use the verification of the Riemann hypothesis to height $3\times10^{12}$ (Platt and Trudgian 2021).

## New results

This appendix states no theorem or proposition.

<!-- results:begin -->
No new results of this appendix are catalogued.
<!-- results:end -->

## Key terms

Bosch type theorem, Bosch's theorem, half Cauchy distribution, hyperbolic complete monotonicity, HCM, elasticity criterion, Pick function, Thorin measure, reciprocal Xi function, Talbot inversion, mpmath
