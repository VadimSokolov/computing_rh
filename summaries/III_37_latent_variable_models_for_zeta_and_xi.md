# Chapter III.37. Latent variable models for zeta and Xi

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 37. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter asks whether the finite field proof of the Riemann hypothesis (Frobenius, the method of Stepanov) can be read as a latent variable model, with the zeta function as partition function and weights that make the functional equation visible. For curves the dictionary is exact (augmentation, the Schmidt weight, Toeplitz positivity). For the integers it gives Chebyshev's identity as a model with two marginals (prime and integer portfolios, with the Moebius hedge as the signed way back), the ramp as the measure of a real cofactor, and the statistic of Xi, symmetric under x to 1/x. The central limit theorem holds in scale only, since the limit law is bounded and not normal. The second part builds the latent model of Xi from the Biane, Pitman and Yor representation, studies its finite systems (truncation, folding at a complex angle, reflection, tilting, coupling) and the Lee and Yang property; the first part also shows that real temperature is easy (log concave tilted Polya laws) while complex temperature is a sign problem. It proves neither (M) nor the Riemann hypothesis.

## From the literature

Edwards and Sokal (1988); Jerrum and Sinclair (1993); the proof of F. K. Schmidt; Stepanov's method and Bombieri's proof (Bombieri 1974); the Weil theorem (1949); Kurlberg and Rudnick (2009); Katz and Sarnak (1999); the heuristic of Good and Churchhouse (1968); Wintner (1935) and Rubinstein and Sarnak (1994); Lagarias (1999); Applegate and Kannan (1991); Biane, Pitman and Yor (2001); Newman (1974, 1976); de Bruijn (1950); Rodgers and Tao (2020); the circle theorem (Lee and Yang 1952); Polya's kernel (Polya 1926); Spira (1967); Gonek and Montgomery (2013); Haglund (2009); Hejhal (1990); Ingham (1942); Deninger (1992); Platt and Trudgian (2021).

## New results

Finite system tables, zero counts and sampling costs are computations, not theorems.

<!-- results:begin -->
- [Chebyshev's identity as a latent variable model with two marginals](../RESULTS.md#lvz_chebyshev_identity_as_a_model_with_two_marginals) (Volume III, Proposition 37.3.1, page 751; proved): The prime portfolio and the integer portfolio are the two marginals of one joint weight on pairs (m,r), the conditional law is a probability, and the Moebius sequence is the only exact way back. New in this book.
- [The statistic of Xi is symmetric under x to 1/x with an explicit Laplace transform](../RESULTS.md#lvz_statistic_of_xi_even_under_inversion_and_its_laplace_transform) (Volume III, Proposition 37.6.1, page 754; proved): The sum over the zeros alone of a_rho x^rho satisfies G(1/x) = G(x)/x without hypothesis, and its Laplace transform in log x is given by the logarithmic derivative of xi. New in this book.
- [Real temperature is easy, log concavity, spectral gap and sampling for the tilted Polya laws](../RESULTS.md#lvz_log_concave_tilted_polya_laws_real_temperature_is_easy) (Volume III, Proposition 37.9.1, page 756; proved): If minus the logarithm of the Polya density has second derivative at least kappa, every tilt is log concave with a Poincare inequality, the second logarithmic derivative of xi is at most 1/kappa, and ratios of xi at real arguments are estimated by sampling. New in this book.
- [Open problem, a positive way back from the integers to the primes with power saving loss](../RESULTS.md#lvz_positive_way_back_from_integers_to_primes_open) (Volume III, Problem 37.1, page 758; open): Find a positive representation of the prime portfolio through the integer portfolio whose loss is a power saving and not the factor two of the sieve. New in this book.
- [Complex Gamma weights, folding Riemann's second proof at a point on the unit circle](../RESULTS.md#lxi_complex_gamma_weights_for_the_folded_representation_of_xi) (Volume III, Proposition 37.11.1, page 761; proved): Folding the theta integral at the point c = e^(i phi) with the angle phi below pi/2 in modulus gives xi as a sum of complex incomplete Gamma weights plus two explicit boundary terms, real on the critical line for every truncation. New in this book.
<!-- results:end -->

## Key terms

latent variable model, partition function, augmentation, Edwards and Sokal, Lee and Yang property, ferromagnet problem, Polya density, Xi function, xi function, Chebyshev identity, Moebius hedge, Mobius function, central limit theorem, log concave, sign problem, folding, de Bruijn Newman constant, generalized gamma convolution, generalised gamma convolution, inequality (M)
