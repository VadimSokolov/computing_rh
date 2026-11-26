# Chapter III.31. Linear response of the rates and the pair formula

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 31. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

How does each rate of the birth and death chain of the zeros respond when a zero moves, and what does a pair of zeros off the line do to the rates? The chapter reports first steps along route B to RH, the rates of the Stieltjes continued fraction of the Riemann xi function. It proves an exact linear response formula for the derivative of the log of every orthogonal polynomial norm, hence of every rate, with respect to every atom position of a discrete measure (by the envelope theorem), and computes the norm of the response at the smooth zeros of zeta at eight levels between 100 and 1900 (between 0.99 and 1.12, floating point), with a conjecture that it stays bounded. It uses the rank update (Uvarov) for the Hankel determinants of a measure with added complex atoms to prove that, if each update determinant changes sign once near a level, exactly two rates in the window around that level are negative, and reports an experiment in which a pair of zeros of zeta is collided and pushed off the line. The path asks for the string of xi, the inverse half of the Gamma chain law, and a positive formula for the rates. The experiments are from an unrefereed note of Polson (November 2026); nothing here proves RH.

## From the literature

The Uvarov transformation of a measure by point masses and the determinant identity for a low rank update (Uvarov 1969); Krein's spectral theory of strings (Krein 1952; Kotani and Watanabe 1982); the Stieltjes continued fraction (Stieltjes 1894); the Thorin representation of generalised gamma convolutions (Thorin 1977); pair correlation of the zeros (Montgomery 1973; Rodgers and Tao 2020).

## New results

<!-- results:begin -->
- [Exact linear response of the Hankel norms (and rates) to the zeros](../RESULTS.md#linear_response_of_rates_to_zeros) (Volume III, Theorem 31.2.1, page 650; proved): An exact formula gives the derivative of the log of every orthogonal polynomial norm, and hence of every rate, with respect to every atom position. New in this book.
- [Bounded linear response norm at the smooth zeros](../RESULTS.md#linear_response_sup_norm_conjecture) (Volume III, Conjecture 31.2.7, page 655; open): The norm of the linear response of the log rates to displacements of the smooth zeros is conjectured to stay bounded in the level, as computed at eight levels between 100 and 1900. New in this book.
- [Added atoms with one sign change make exactly two rate factors negative](../RESULTS.md#two_negative_rates_per_complex_pair) (Volume III, Proposition 31.3.3, page 657; proved): If the two rank update determinants for added atoms each change sign once near level N, then among the rates at levels 2N-4 to 2N+1 exactly two of the factors relative to the base measure are negative, at positions fixed by where the two determinants change sign. New in this book.
<!-- results:end -->

## Key terms

linear response, Jacobian of the rates, orthogonal polynomials, Hankel determinants, rank update, Uvarov transformation, pair of zeros, complex conjugate pair, negative rates, Krein string, route B, envelope theorem, collision of zeros, de Bruijn Newman, backward heat flow, Stieltjes continued fraction
