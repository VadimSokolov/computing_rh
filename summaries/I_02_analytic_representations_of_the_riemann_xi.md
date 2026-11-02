# Volume I, Chapter 2. Analytic representations of the Riemann Xi function

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 2. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function*, edition 2, https://doi.org/10.5281/zenodo.23891944. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter collects the analytic representations of the completed Riemann zeta function $\xi$ and the Riemann Xi function $\Xi$ on which the rest of the book rests. It asks how the zeros of $\xi$, the Pólya kernel and density, the heat trace and the Thorin density at a shifted basepoint are related, and when the reciprocal of $\xi(\frac12+\sqrt u)/\xi(\frac12)$ is the Laplace transform of a generalised gamma convolution (GGC). From the Hadamard product and the partial fraction expansion of $\xi'/\xi$ it defines the clock variable $T$ with Laplace transform $1/G$, and proves a rigidity proposition: the clock is a GGC if and only if the Riemann Hypothesis (RH) holds, with Thorin measure the sum of unit atoms at the squared ordinates. It derives the Pólya kernel $\Phi$ from the theta function, shows that it is even, positive and analytic in a strip, and gives the Pólya density $p=\Phi/\Xi(0)$ with moment generating function $\xi(\frac12+s)/\xi(\frac12)$, so that the clock and the characteristic function of the Pólya variable form a van Dantzig pair exactly when RH holds. The heat trace $W(t)$ of the zeros is given in spectral and arithmetic form and checked numerically. The chapter proves a Poisson structure formula writing the Thorin density at basepoint $\alpha$ as Cauchy kernels centred at the zeros, which is a GGC for every $\alpha$ to the right of all zeros (in particular $\alpha\ge1$, unconditionally), and derives the same density from the primes through the sine squared lemma, compared in floating point with the Pólya route. Statements about $\alpha<1$ are conditional on RH or on a zero free half plane.

## From the literature

- Hadamard factorisation and the partial fraction expansion of $\xi'/\xi$ (classical; Titchmarsh)
- The Pólya kernel representation of $\Xi$ through the theta function (Pólya 1926; Titchmarsh)
- Generalised gamma convolutions and Thorin measures (Thorin 1977; Bondesson 1992)
- Van Dantzig pairs and the equivalence of RH with the Pólya law pairing (van Dantzig; Lukacs 1968; Konstantopoulos, Patie and Sarkar 2024)
- The hitting time construction of $2\xi(s)$ as a Mellin transform (Williams 1990; Biane, Pitman and Yor 2001)
- The sine squared lemma and the Thorin density at a shifted basepoint, the rigidity of the clock and the unconditional heat trace identity (Polson 2017; Polson 2026)

## New results

Results of this chapter, with their status.

<!-- results:begin -->
- [Rigidity of the clock, a GGC if and only if RH](../RESULTS.md#rigidity_clock_ggc_iff_rh) (Volume I, Proposition 2.1.1, page 18; proved): The clock, the reciprocal of the normalised completed zeta function in the variable $u=s^2$, is the Laplace transform of a generalised gamma convolution exactly when the Riemann Hypothesis holds, and then its Thorin measure is the sum of unit atoms at the squared ordinates of the zeros. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.
- [Shifted clock is a GGC to the right of every zero](../RESULTS.md#shifted_clock_ggc_right_of_every_zero) (Volume I, Theorem 2.1.2, page 22; proved): To the right of every zero of the completed zeta function, in particular at every basepoint $\alpha\ge 1$, the shifted clock is the Laplace transform of a generalised gamma convolution with an explicit positive Thorin density, and it is one exactly when $\alpha$ is at least the supremum of the real parts of the zeros. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.
- [Positivity of the Thorin density on one line, open question](../RESULTS.md#one_line_positivity_exclusion_open) (Volume I, Remark 2.1.4, page 23; open): It is not known, without assuming the Riemann Hypothesis, whether nonnegativity of the real part of the logarithmic derivative of the completed zeta function on a single vertical line excludes zeros to the right of that line. New in this book.
<!-- results:end -->

## Key terms

Riemann Xi function, completed zeta function, xi function, Hadamard product, partial fraction expansion, clock variable, reciprocal xi function, generalised gamma convolution, generalized gamma convolution, GGC, Thorin measure, Thorin density, rigidity, Pólya kernel, Polya kernel, Pólya density, van Dantzig pair, heat trace, Cauchy flow, Poisson kernel, shifted basepoint, sine squared lemma, Riemann Hypothesis, RH
