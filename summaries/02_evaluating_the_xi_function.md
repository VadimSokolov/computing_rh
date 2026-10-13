# Chapter 2. Evaluating the Xi function

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 2. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23331091. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter asks how to evaluate the Riemann xi function $\xi$ and its derivative at complex arguments, the raw material for computing the Thorin density at every resolution. The Pólya route applies the trapezoid rule to Pólya's integral $\Xi(\theta)=\int\Phi(u)\cos(\theta u)\,du$. It needs neither primes nor zeros and is spectrally accurate because the Pólya kernel $\Phi$ is analytic in the strip $|\operatorname{Im}u|<\pi/4$. The evaluator is two fixed polynomials in Arb ball arithmetic; below height 200 it gives about 60 correct digits and recovers the 79 zeros, with the count read from the argument principle. Cancellation costs about $0.341\,\theta$ digits at height $\theta$; tilting the contour into the strip reduces the loss to about $\delta\theta/\log 10$ digits and gives 25 digits of the Thorin density at height 5000. The chapter compares the route with the Riemann Siegel formula, the Odlyzko Schönhage algorithm and Turing's method, resolves the Lehmer pair near height 7005, reads the count $N(1000)=649$ from the Thorin mass, and ends with a worked computation from the theta series to the first zero. Ball radii certify rounding only, so the results are computed, not certified.

## From the literature

Pólya's integral for $\xi$ (Pólya 1926); the exponentially convergent trapezoid rule (Trefethen and Weideman 2014); Arb ball arithmetic (Johansson 2017); the Riemann Siegel formula (Siegel 1932; Gabcke 1979); Gram points and Gram's law (Gram 1903), Rosser's rule (Rosser, Yohe and Schoenfeld 1969) and their failures (Brent 1979; Trudgian 2011); Backlund's count (Backlund 1918; Hutchinson 1925); Turing's method (Turing 1953; Lehman 1970; Trudgian 2011); the Odlyzko Schönhage algorithm (Odlyzko and Schönhage 1988); the Lehmer pair (Lehmer 1956); the verification of RH to height $3\times10^{12}$ (Platt and Trudgian 2021).

## New results

The tilted Pólya representation is new in the book. Digit counts, timings and zero counts are computed, with ball arithmetic certifying rounding only.

<!-- results:begin -->
- [Tilted Pólya representation of xi and its derivative](../RESULTS.md#tilted_polya_representation_of_xi) (Proposition 2.1.1, page 34; proved): Moving the line of integration of Pólya's integral into the strip of analyticity of the Pólya kernel gives exact integral formulas for $\xi$ and $\xi'$ on a line at distance $\delta$ from the edge of the strip, with kernel bounds of order $\delta^{-5/2}$. New in this book.
<!-- results:end -->

## Key terms

Riemann xi function, Pólya route, Polya route, Pólya kernel, trapezoid rule, spectral accuracy, tilted contour, ball arithmetic, Arb, Riemann Siegel formula, Odlyzko Schönhage algorithm, Gram points, Rosser's rule, Turing's method, Backlund's method, argument principle, zero counting, Lehmer pair, Thorin density, Thorin mass
