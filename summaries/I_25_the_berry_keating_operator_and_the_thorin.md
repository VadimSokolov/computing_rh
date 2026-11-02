# Volume I, Chapter 25. The Berry Keating operator and the Thorin programme

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 25. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function*, edition 2, https://doi.org/10.5281/zenodo.23891944. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter asks what is true about the Berry Keating operator, the symmetrised quantisation $xp$ of the dilation generator, and how it meets the Thorin measure results of the book. It computes the spectral measure of Pólya's kernel for the dilation generator through a Mellin transform, proves that no finite rank change of the resolvent of a Gram regular operator, in particular no change of finitely many boundary conditions, can give the zeros as spectrum, and writes the Cauchy deconvolution of the Thorin density with its comb on the primes (conditional on RH). Under RH the squared ordinates of the zeros are also the spectral measure of a Krein string without initial gap, with explicit Weyl function, length and infinite mass at its end. The chapter gives no route to the Riemann Hypothesis, and nothing in it separates $\zeta$ from the Davenport and Heilbronn function: the results use the explicit formula and the functional equation, not the Euler product beyond them. The string statement is RH by inclusion. Bead data of truncated strings and a negative Langevin computation are reported as computed.

## From the literature

- the Berry Keating proposal and the Planck cell count of the zeros (Berry and Keating 1999)
- the operator has Lebesgue spectrum of multiplicity two, and a restriction with defect indices $(2,2)$ (Bagarello and Kuzel 2026, an unrefereed preprint; the chapter corrects the form of its extension statements)
- the zero counting formula with the Riemann Siegel phase, and $S(T)$ unbounded (Selberg 1946)
- the explicit formula (Weil 1952; Guinand 1948)
- Krein's correspondence between strings and spectral measures (Dym and McKean 1976; Kotani and Watanabe 1982), and hitting times of strings (Kent 1980)
- a positive proportion of the zeros are simple (Alpöge and Furman 2026)

## New results

<!-- results:begin -->
- [Spectral measure of Pólya's kernel for the Berry Keating dilation generator](../RESULTS.md#bk_spectral_measure_of_polya_kernel_for_dilation_generator) (Volume I, Proposition 25.3.1, page 713; proved): The function built from Pólya's kernel has spectral measure |xi(1/2+iE)|^2 dE/(4 pi) for the dilation generator, so the zeros appear as zeros of a density and not as eigenvalues. New in this book.
- [Finite rank changes cannot give the zeros as spectrum of a Gram regular operator](../RESULTS.md#bk_finite_rank_obstruction_for_berry_keating_boundary_condition) (Volume I, Proposition 25.4.1, page 713; proved): A finite rank change of the resolvent of a Gram regular operator, in particular a change of finitely many boundary conditions, cannot produce the squared ordinates of the zeros as spectrum, unconditionally. New in this book.
- [Cauchy form of the Thorin density and its comb on the prime side (conditional on RH)](../RESULTS.md#bk_cauchy_deconvolution_comb_conditional_on_rh) (Volume I, Proposition 25.5.1, page 714; conditional on the Riemann Hypothesis): Under RH the Thorin density at level sigma is the zero measure smoothed by a Cauchy kernel of width sigma minus one half, and the explicit formula gives the transform of the zero measure as a comb on the primes. New in this book.
- [Krein string with spectral measure the squared ordinates of the zeros (conditional on RH)](../RESULTS.md#bk_krein_string_of_the_zeros_conditional_on_rh) (Volume I, Theorem 25.6.1, page 715; conditional on the Riemann Hypothesis): Under RH the squared ordinates of the zeros are the spectral measure of a Krein string with no initial gap, whose Weyl function, length 0.02310499311541897 and infinite mass at its end are explicit; this restates a known equivalence and gives no route to RH. New in this book.
<!-- results:end -->

## Key terms

Berry Keating operator, xp Hamiltonian, Hilbert and Pólya conjecture, dilation generator, Mellin transform, spectral measure, Pólya kernel, Gram regular operator, finite rank perturbation, boundary condition, defect indices, Krein string, Weyl function, hitting time, explicit formula, Cauchy deconvolution, Thorin measure, Davenport and Heilbronn function, Riemann Hypothesis, RH
