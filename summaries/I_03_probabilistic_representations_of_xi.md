# Volume I, Chapter 3. Probabilistic representations of Xi

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 3. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function*, edition 2, https://doi.org/10.5281/zenodo.23861971. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter asks whether the clock variable $T$ of the reciprocal Xi function, the positive random variable with Laplace transform $\xi(\tfrac12)/\xi(\tfrac12+\sqrt u)$ for the completed Riemann zeta function $\xi$, can be represented by Brownian motion itself and not only through its spectrum, the set $\{\gamma^2\}$ of squared ordinates of the zeros. Under the Riemann Hypothesis (RH) with simple zeros, Kent's theorem makes $T$ the hitting time of a level by a diffusion, but that diffusion is described only through the zeros. The chapter proves that unconditionally $T$ is not the first passage time of Brownian motion with drift $c\ge0$ to an independent random level, by a growth argument along a parabola in the $u$ plane using Stirling's formula. It then treats the inverse first passage problem: by theorems of Ekström and Janson a curved Brownian boundary exists whose crossing time has the law of the clock, and the chapter computes it from the Volterra integral equation of the first kind, in three zones (small times, the bulk, and a right tail with a logarithmic correction), and compares it with the Roberts Shortland tangent boundary and with simulation. Wald's identities constrain the exact boundary. Finally Williams' representation of the Pólya kernel through hitting times of the three dimensional Bessel process (Biane, Pitman and Yor) explains the tilted contour, supplies a Gamma(2) benchmark with a closed Thorin measure, and gives Maślanka's series for $(\sigma-1)\zeta(\sigma)$ as Laguerre moments of Williams' variable; none of this transfers positivity to the Thorin measure of the clock. The boundary computations are floating point.

## From the literature

- Kent's theorem: under RH with simple zeros the clock is a hitting time of a diffusion (Kent 1978)
- The inverse first passage problem, and existence and uniqueness of the continuous boundary (Cheng, Chadam, Cheng and Sheu 2011; Ekström and Janson 2016), with the Volterra equation of the first kind (Peskir 2002)
- The tangent approximation of a Brownian first passage hazard (Roberts and Shortland 1995)
- De Bruijn's exponential Tauberian theorem (Bingham, Goldie and Teugels 1987)
- Williams' representation of $2\xi(s)$ through Bessel hitting times and the variable $S_2$ (Williams 1990; Biane, Pitman and Yor 2001)
- Maślanka's series for the zeta function (Maślanka 2001), proved rigorously by Báez Duarte (2010)
- Exact sampling of Gamma(2) type sums (Devroye 1986)

## New results

The status words below refer to the proved statements on $T$; the boundary numbers are computed in floating point.

<!-- results:begin -->
- [The clock is not a Brownian first passage time to an independent random level](../RESULTS.md#clock_not_first_passage_time_to_random_level) (Volume I, Proposition 3.1.1, page 25; proved): The clock variable of the Riemann xi function is not, for any constant drift and any independent random level, the first passage time of Brownian motion with that drift to that level, and the proof uses no hypothesis on the zeros. New in this book.
- [Maślanka coefficients as Laguerre moments of Williams' variable](../RESULTS.md#maslanka_coefficients_as_laguerre_moments_of_williams_variable) (Volume I, Proposition 3.1.2, page 29; proved): The coefficients of Maślanka's series for (sigma minus 1) zeta(sigma) are Laguerre moments of Williams' variable, with an explicit generating function and an integral formula, with no zeros and no primes, for (sigma minus 1) zeta(sigma) on the strip from zero to two. New in this book.
<!-- results:end -->

## Key terms

clock variable, reciprocal Xi function, Brownian motion, first passage time, hitting time, inverse first passage problem, Brownian boundary, curved boundary, Volterra equation, Ekström and Janson, Roberts Shortland tangent boundary, inverse Gaussian law, Kent's theorem, Wald's identities, Williams' representation, Bessel process, Biane Pitman Yor, Gamma(2) structure, Maslanka series, Maślanka series, Baez Duarte, Laguerre moments, Pólya kernel, Polya kernel, Thorin measure
