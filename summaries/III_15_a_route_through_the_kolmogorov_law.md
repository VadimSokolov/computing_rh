# Chapter III.15. A route through the Kolmogorov law

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 15. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Can the Kolmogorov law of the supremum of a Brownian bridge be used as a route to the Riemann Hypothesis? This chapter of Volume III puts Ramanujan's pole removing combination, Dirichlet's eta function, the Kolmogorov Smirnov limit law, Riemann's kernel and the Thorin measures of the reciprocal Xi function on one line. It proves that the positive even kernel of the Kolmogorov law equals twice Riemann's kernel convolved with two tilted indicator kernels, derives the zeros of the resulting cosine transform and its log concavity under a hypothesis, and proves that the logarithmic derivative of xi is the tilted mean of log K plus an elementary term, where K is the supremum of the absolute Brownian bridge. RH is then equivalent to one inequality for that tilted mean. This is a reformulation, not a reduction, and the chapter says what blocks three attacks (probabilistic, analytic, descent on the ladder). It ends with three open problems: a condition on Thorin measures putting the Mellin zeros of an exponential series on one line, a monotone functional along the ladder, and a pole removing kernel with no extra zeros. Nothing here proves RH.

## From the literature

The Kolmogorov Smirnov limit law (Kolmogorov 1933) and its representation as a hitting time of a Bessel process, with the Mellin transform through xi (Biane, Pitman and Yor 2001); the criterion that RH holds if and only if the real part of the logarithmic derivative of xi is positive to the right of the critical line (Lagarias 1999; Garunkštis 2002); the monotonicity reformulation (Sondow and Dumitrescu 2010); the Lee and Yang property of Fourier transforms (Newman 1976); the Prékopa theorem on log concave functions (Prékopa 1973); the concavity of log Phi of the square root of t, used for the log concavity of Riemann's kernel and not checked against the paper (Csordas, Norfolk and Varga 1986); Bohr's almost periodicity; verification of RH to height 3 times 10 to the 12 (Platt and Trudgian 2021).

## New results

The results below are from an unrefereed note of Polson (November 2026).

<!-- results:begin -->
- [Riemann's kernel smoothed into the Kolmogorov kernel](../RESULTS.md#riemann_kernel_kolmogorov_smoothing) (Volume III, Theorem 15.3.1, page 287; proved): The positive even kernel built from the Kolmogorov law is twice Riemann's kernel convolved with a pair of exponentially tilted indicator kernels, so its law is that of Riemann's variable plus the difference of two independent tilted uniform variables. New in this book.
- [Zeros of the smoothed kernel transform and log concavity of Phi_R](../RESULTS.md#smoothing_zeros_and_log_concavity) (Volume III, Corollary 15.3.2, page 287; proved): The cosine transform of the Kolmogorov kernel has the zeros of Xi plus an extra row of zeros on the edge of the strip, and it is log concave whenever Riemann's kernel is, by Prekopa's theorem. New in this book.
- [Tilted mean of log K and the missing step, equivalent to RH](../RESULTS.md#kolmogorov_tilted_mean_missing_step) (Volume III, Theorem 15.5.1, page 288; proved): The logarithmic derivative of xi is the tilted mean of the log of the Kolmogorov supremum plus an elementary term, so RH is equivalent to an inequality for that tilted mean; this is a reformulation of RH, not a reduction. New in this book.
- [Open problem, Thorin condition for Mellin zeros of an exponential series on one line](../RESULTS.md#kolmogorov_mellin_zeros_exponential_series_problem) (Volume III, Problem 15.1, page 290; open): Find a condition on the Thorin measure of an exponential series under which the zeros of its Mellin transform in the strip lie on one line, satisfied by the Kolmogorov case. New in this book.
- [Open problem, a monotone functional along the ladder](../RESULTS.md#ladder_monotone_functional_problem) (Volume III, Problem 15.2, page 290; open): Find a functional of the Thorin density of the rung that is monotone in the rung height and has a fixed value at one half, and decide whether the ratio of the shift to the rung function is monotone. New in this book.
- [Open problem, a pole removing kernel with no extra zeros](../RESULTS.md#pole_removing_kernel_problem) (Volume III, Problem 15.3, page 290; open): Decide whether an infinite combination, a Dirichlet series, can remove the pole of zeta without adding zeros off the critical line, which no finite combination can do. New in this book.
<!-- results:end -->

## Key terms

Kolmogorov law, Kolmogorov Smirnov statistic, Brownian bridge supremum, Ramanujan combination, Dirichlet eta function, Riemann kernel, Xi function, tilted mean, logarithmic derivative of xi, Lagarias criterion, Prékopa, log concavity, Thorin measure, generalised gamma convolution, generalized gamma convolution, ladder, rung, Mellin zeros, Riemann hypothesis, RH
