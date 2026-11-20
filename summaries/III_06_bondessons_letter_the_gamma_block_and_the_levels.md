# Chapter III.6. Bondesson's letter, the gamma block and the levels three quarters and seven eighths

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 6. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

In a letter of 23 August 2016 Lennart Bondesson showed that e to the minus gamma root s over Gamma(1 plus root s), with gamma Euler's constant, is the Laplace transform of a generalised gamma convolution (GGC). This chapter proves the sign condition of the letter, which Bondesson checked with Maple plots, by a damped sine lemma, gives the Thorin density in closed form, and generalises the letter to every centred infinitely divisible law with only negative jumps, which is the sine squared lemma for a Lévy measure. It shows that Gamma(a) exp(minus beta root s) over Gamma(a plus root s) is a GGC transform exactly when beta is at least minus digamma of a. The level criterion says that the shifted clock xi(alpha) over xi(alpha plus root s) is a GGC transform exactly when xi has no zeros to the right of alpha; it is proved through the Hadamard product, gives positivity for alpha at least 1 and, by a theorem of OpenAI (2026), for every alpha at least 7/8, and a descent criterion turns Poisson deconvolution into an exact test of RH. Computations at levels 3/4 and 7/8 (EM and Tikhonov deconvolution, a planted off line pair) are reported. Below level 1 the letter helps computation, not proof.

## From the literature

- Bondesson's unpublished letter (2016) and the Pick characterisation of generalised gamma convolutions (Bondesson 1992).
- The Thorin representation of generalised gamma convolutions (Thorin 1977; Bondesson 1992); the sine squared lemma is from Polson (2017).
- The EM algorithm for Poisson deconvolution (Richardson 1972; Lucy 1974; Shepp and Vardi 1982; Vardi and Lee 1993).
- RH as positivity of the real part of xi prime over xi (Lagarias 1999; Garunkštis 2002; Sondow and Dumitrescu 2010).
- The level criterion of Polson 2026, here with a new proof.
- The zero free half plane Re s > 7/8, a theorem of OpenAI (2026).

## New results

<!-- results:begin -->
- [Generalisation of Bondesson's letter to infinitely divisible laws with negative jumps](../RESULTS.md#bondesson_letter_generalised_to_centred_infinitely_divisible_laws_with_negative_jumps) (Volume III, Theorem 6.2.1, page 96; proved): The reciprocal of the moment generating function of a centred infinitely divisible variable with only negative jumps, evaluated at the square root of s, is a generalised gamma convolution transform with an explicit Thorin measure. New in this book.
- [Sharp condition for the gamma block with drift to be a GGC transform](../RESULTS.md#sharp_condition_for_gamma_block_with_drift_to_be_ggc_transform) (Volume III, Proposition 6.2.3, page 97; proved): The function Gamma(a) exp(minus beta root s) over Gamma(a plus root s) is a GGC transform exactly when beta is at least minus digamma of a, for every positive a. New in this book.
- [Level criterion, a new proof by the Hadamard product (GGC at level alpha iff no zeros to the right)](../RESULTS.md#level_criterion_ggc_iff_zero_free_half_plane_new_proof_by_hadamard_product) (Volume III, Theorem 6.3.1, page 98; proved): For alpha at least one half, the shifted clock xi(alpha) over xi(alpha plus root s) is a GGC transform exactly when xi has no zeros with real part above alpha, and then the Thorin density is a sum of Cauchy kernels centred at the ordinates. New in this book.
- [Consequences of the level criterion at levels 1, 1/2 and 7/8](../RESULTS.md#consequences_of_level_criterion_levels_one_one_half_seven_eighths_and_poisson_convolution) (Volume III, Corollary 6.3.2, page 99; proved): The Thorin measure U_alpha is positive for every alpha at least 1 unconditionally and for every alpha at least 7/8 by a theorem of OpenAI (2026), GGC at level one half is equivalent to RH, and under RH U_alpha is the Poisson smoothing of the sum of atoms at the ordinates. New in this book.
- [Harmonic structure of the archimedean and arithmetic densities](../RESULTS.md#harmonic_structure_of_archimedean_and_arithmetic_densities_across_levels) (Volume III, Proposition 6.4.1, page 100; proved): The archimedean density is harmonic in the half plane to the right of minus 2 and so moves between levels by the Poisson extension, while the arithmetic density is harmonic away from the zeros of zeta. New in this book.
- [Descent criterion, RH if and only if U_alpha is a Poisson smoothing of a signed measure](../RESULTS.md#descent_criterion_rh_iff_poisson_deconvolution_of_thorin_measure_exists) (Volume III, Proposition 6.6.1, page 108; proved): If xi has no zeros to the right of theta, then RH holds exactly when the Thorin measure U_alpha at a level alpha at least theta is the Poisson kernel convolved with some locally finite signed measure, which is then the sum of atoms at the ordinates. New in this book.
<!-- results:end -->

## Key terms

Bondesson letter, generalised gamma convolution, generalized gamma convolution, GGC, Thorin measure, HCM, Pick function, Stieltjes inversion, sine squared lemma, infinitely divisible, Lévy measure, level criterion, shifted clock, Hadamard product, zero free half plane, Cauchy kernel, Poisson kernel, descent, deconvolution, Tikhonov, level 7/8, Riemann Hypothesis, RH
