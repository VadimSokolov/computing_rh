# Chapter III.12. Zeros from the primes: positive definiteness, the truncation bound and blind location

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 12. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Can the zeros of the Riemann zeta function be located from the primes alone? The chapter works with the Riesz mean of the von Mangoldt function, whose normalised remainder f(u), u = log x, is computed from the primes up to x. Under the Riemann Hypothesis f is a positive definite cosine sum over the ordinates of the zeros, with weights 2/(1/4 + gamma^2) and f(0) = 2 + gamma_E - log 4 pi. The chapter tests positive definiteness of f from the primes below 10^7 (smallest eigenvalue of the Toeplitz matrix, computed, with a proved lower bound that uses the verified height 3 x 10^12). It proves, under RH, a truncation bound for the level alpha density and a consistency theorem for the cosine transform at rate 1 over log X, and unconditionally a propagation bound for location error and a loss theorem showing that linear Fourier division amplifies noise by exp(a Omega) = X^a. The atomic prior (Gauss and Newton joint least squares from a peak search) beats this law in computation: all 138 zeros below height 300 are found from the primes below 10^7, and at level 0.51 the density is obtained to about four figures. This is a method of computation, not a route to a proof of RH.

## From the literature

Landau's oscillation theorem, by which one sided boundedness of the remainder already implies RH (Landau 1906); positivity of the level alpha density under RH (Lagarias 1999, correction 2005); RH verified up to height 3 x 10^12 (Platt and Trudgian 2021); EM facts (Csiszar and Tusnady 1984; Vardi, Shepp and Kaufman 1985); the deconvolution path for mixtures (Madrid Padilla, Polson and Scott 2018).

## New results

The results are a truncation bound and a consistency theorem under RH, a propagation bound and a loss theorem, all with computed tables.

<!-- results:begin -->
- [Truncation error of the level alpha density from the primes below X (under RH)](../RESULTS.md#xfp_truncation_error_of_level_alpha_density_from_primes_under_rh) (Volume III, Theorem 12.3.1, page 215; conditional on the Riemann Hypothesis): Under RH, truncating the damped error transform at the primes below X changes the level alpha data by at most a constant times X to the power minus (alpha minus one half), the law that makes descending towards the critical line costly. New in this book.
- [Consistency of the cosine transform of the prime remainder at rate 1 over log X (under RH)](../RESULTS.md#xfp_consistency_of_cosine_transform_of_prime_remainder_under_rh) (Volume III, Theorem 12.4.1, page 215; conditional on the Riemann Hypothesis): Under RH the finite cosine transform of the Riesz mean remainder over [0, log X] recovers each ordinate atom, with its multiplicity and weight, at rate 1 over log X, and tends to zero at every non ordinate (a mean value argument of Bohr type). New in this book.
- [Propagation of ordinate location error into the level alpha density](../RESULTS.md#xfp_propagation_of_ordinate_location_error_into_cauchy_kernel) (Volume III, Proposition 12.5.1, page 217; proved): Moving one ordinate by at most eta changes its Cauchy kernel by at most eta times an explicit bound, so near a zero the relative error is about 0.65 eta over a and away from the zeros it is of order a eta. New in this book.
- [Loss of linear reconstruction by Fourier division, factor exp(a Omega)](../RESULTS.md#xfp_loss_of_fourier_division_factor_exp_a_omega) (Volume III, Theorem 12.6.1, page 218; proved): Recovering a band limited convolution of a measure from Cauchy smoothed data with noise by Fourier division amplifies the noise in L2 by exactly e to the a Omega over pi at worst, and every linear estimator that is exact on noise free data has the same worst case error, so with Omega = log X the factor is X to the a. New in this book.
<!-- results:end -->

## Key terms

zeros from the primes, Riesz mean of the von Mangoldt function, Chebyshev function, explicit formula, positive definite function, Toeplitz matrix eigenvalue, truncation error, level alpha density, Cauchy deconvolution, Fourier division, band limited reconstruction, ill posed inverse problem, cosine transform, Bohr mean value, line spectrum estimation, joint least squares, Gauss Newton, zero location, critical line, Riemann Hypothesis, Riemann zeta zeros, ordinates, EM algorithm, deconvolution path
