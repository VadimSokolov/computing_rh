# Chapter III.10. The classical sources: von Mangoldt, Hardy, Littlewood, Ramanujan and Pólya

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 10. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Does the oscillator of Chapter III.6 have classical ancestors, and what do they say about it? The chapter reads four sources for what they say about the oscillator f(L), the sum over the zeros of 2(1 - cos(gamma L))/gamma^2 whose values lie in the band [0, W] under RH: von Mangoldt's proof of Riemann's explicit formula for psi(x) (1895), Ramanujan's paper on highly composite numbers (1915) with its band for the sum of x^rho over rho squared, the Acta paper of Hardy and Littlewood (1916) with the Riesz and Hardy Littlewood criteria, and the Inequalities of Hardy, Littlewood and Polya (1934) with the Hilbert inequality of Montgomery and Vaughan. Each reading is checked by computation from the primes to 1e9. The chapter proves an exact formula for Ramanujan's zero sum in terms of the prime sum of Lambda(n) log(x/n), shows by Landau's theorem that Ramanujan's band is a criterion equivalent to RH, and applies the mean value theorem of Montgomery and Vaughan on finite windows of the oscillator. It explains why exponential kernels are numerically blind to zeros at height 1000 while the 2 sin squared form keeps every zero in view.

## From the literature

von Mangoldt's explicit formula for psi(x) (von Mangoldt 1895), Ramanujan's band under RH and his highly composite numbers (Ramanujan 1915, 1997), the Riesz criterion (Riesz 1916), the Hardy and Littlewood criterion (Hardy and Littlewood 1916), Littlewood's sign change theorem (Littlewood 1914), the Hilbert inequality and its finite window mean value form (Hardy, Littlewood and Polya 1934; Montgomery and Vaughan 1974), the closed form of the sum of rho to the minus 2 (Keiper 1992; Lehmer 1988), Robin's criterion (Robin 1984), Landau's theorem on integrals with nonnegative integrand (Landau; see Ingham 1932), and Polya's integral representation of xi (Polya 1926).

## New results



<!-- results:begin -->
- [Ramanujan's zero sum from the primes (exact formula)](../RESULTS.md#ramanujan_sum_from_primes) (Volume III, Proposition 10.2.1, page 195; proved): The sum over zeta zeros of x^rho over rho squared equals an explicit expression in the prime sum of Lambda(n) log(x/n), with an explicit constant S_2 given by Stieltjes constants, as in Keiper and Lehmer; the formula is von Mangoldt's explicit formula integrated term by term. New in this book.
- [Ramanujan's band as a criterion for RH](../RESULTS.md#ramanujan_band_criterion_equivalence) (Volume III, Theorem 10.2.2, page 196; proved): RH is equivalent to Ramanujan's symmetric bound on the zero sum of x^rho over rho squared, and also to a one sided lower bound of order square root of x, by Landau's theorem on integrals of nonnegative functions. New in this book.
<!-- results:end -->

## Key terms

oscillator, von Mangoldt formula, explicit formula, Ramanujan, highly composite numbers, Ramanujan band, Riesz criterion, Hardy Littlewood criterion, Montgomery Vaughan inequality, Hilbert inequality, Landau theorem, 2 sin squared kernel, Chebyshev psi function, Robin criterion, colossally abundant numbers
