# Volume II, Appendix A. The Pólya evaluators

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation* by Nicholas G. Polson and Vadim Sokolov (2026), appendix A. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation*, edition 2, https://doi.org/10.5281/zenodo.23861976. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This appendix describes the two evaluators of $\xi$ and $\xi'$ built on the Pólya kernel $\Phi$ of the Riemann Xi function, and states exactly what their error balls do and do not bound. The question it answers is how the numbers of the book, the Thorin density $\rho_\varepsilon$, counts of zeros on the critical line, the Laplace transform of the clock and the Keiper Li coefficients, are computed, and which errors a ball arithmetic enclosure leaves outside the stated radius. The untilted evaluator computes node values $\Phi(jh)$ from the theta series and forms trapezoid approximations to $\xi(\frac12+s)$ and $\xi'(\frac12+s)$ through two polynomials in $e^{sh}$; its defaults are 240 bits, step $1/200$ and $|u|\le3$, and most runs use 420 bits. The tilted evaluator integrates along a line shifted into the strip of analyticity of $\Phi$, which keeps the cost manageable at large height, and the appendix lists the shift, precision, step and cut used for the Lehmer pair and its de Bruijn Newman flow, for the count of zeros at height 1000 and for the clock transform. The radii of the returned balls bound only the rounding errors of the finite sums. Three errors lie outside them: the discretisation error of the trapezoid rule (bounded in the main text), the cut of the axis, and the truncation of the theta series, whose stopping test does not add the tail to the radius and, because it compares double precision numbers, never triggers at working precision of 1055 bits or more. The untilted route needs more than 1054 bits for thirty digits above height about 840, so there the stopping test has to be made in ball arithmetic. The runs reported in the book use at most 600 bits and are not affected. Results computed with these evaluators are computed, and certified only where the three outside errors are bounded separately.

## From the literature

- Ball arithmetic as the basis of rigorous numerical enclosures, and arbitrary precision reference values in floating point, as cited by the appendix
- The theta series and the Pólya kernel representation of $\Xi$ (Pólya; Titchmarsh)

## New results

None are claimed beyond the error accounting.

<!-- results:begin -->
No new results of this appendix are catalogued.
<!-- results:end -->

## Key terms

Pólya evaluator, Polya evaluator, Pólya kernel, tilted contour, trapezoid rule, theta series truncation, ball arithmetic, interval arithmetic, rounding error, working precision, Riemann Xi function, xi function, Thorin density, Keiper Li coefficients, de Bruijn Newman flow, floating point reference values
