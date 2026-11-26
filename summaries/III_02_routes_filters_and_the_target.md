# Chapter III.2. Routes, filters and the target

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 2. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter assembles the results of Volumes I and II into one route to the Riemann Hypothesis (RH), says which steps are theorems, and states the one that is missing. It does not prove RH. The route is the descent in level of the Thorin measure, read in the language of Lee and Yang. Its target has six equivalent forms and five more. The top of the descent, level one, is proved with no hypothesis in each form, and the zero free half plane Re s above 7/8 is a theorem of OpenAI (2026), so positivity at every level above 7/8 is proved. Three filters rule out finite computations, arguments that ignore the Euler product (the Davenport and Heilbronn function has zeros off the line), and properties preserved by Gaussian damping. In the heat time there is no room. Four roads are closed, and routes A (descent in level) and C (heat flow) are weighed. Curves over finite fields are decided by Hasse and Weil positivity through Toeplitz sections, and certified ball arithmetic runs on zeta are reported. The missing step is a positive form built from the primes with the values of Weil's functional.

## From the literature

- Weil's positivity criterion (Weil 1952)
- Lee and Yang theorem and ferromagnetic couplings (Lee and Yang 1952; Lieb and Sokal 1981, as restated by Kozitsky)
- Hasse and Weil: the Riemann Hypothesis for curves over finite fields
- the closed half line of heat times with real zeros (de Bruijn 1950; Newman 1976) and the bound Lambda at least 0 (Rodgers and Tao 2020)
- the zero free half plane Re s above 7/8 (OpenAI 2026)
- Schoenberg's theorem on Pólya frequency functions (Schoenberg 1951)
- the Davenport and Heilbronn function with zeros off the critical line (Davenport and Heilbronn 1936)
- the Hilbert and Pólya programmes (Connes; Deninger)

## New results

<!-- results:begin -->
- [The Riemann Hypothesis for a curve over a finite field is equivalent to positive semidefinite Toeplitz matrices](../RESULTS.md#weil_riemann_hypothesis_for_curves_iff_toeplitz_positive) (Volume III, Proposition 2.4.1, page 16; proved): The Riemann Hypothesis for a curve over a finite field holds if and only if the Toeplitz matrices built from its normalised point counts are positive semidefinite for every size. New in this book.
- [A Curie Weiss coupling of strength a shrinks the zero strip of Xi type laws from delta squared to delta squared minus a](../RESULTS.md#curie_weiss_coupling_shrinks_strip) (Volume III, Proposition 2.5.1, page 18; proved): Coupling laws of Xi type with a Curie Weiss weight of strength a puts the zeros of the partition function in a strip whose squared half width is the maximum of delta squared minus a and zero. New in this book.
- [Problem: construct from the primes a positive semidefinite form for Weil's criterion](../RESULTS.md#question_positive_form_from_primes_for_weil_criterion) (Volume III, Problem 2.1, page 23; open): The book names as the missing piece a positive semidefinite Hermitian form built from the primes that equals the sum over zeros on a dense space of test functions, equivalently monotonicity of the modulus of xi in the real part or lowering the Pólya spin coupling to zero. New in this book.
- [Problem: lower bounds for the margins that the runs only observe](../RESULTS.md#question_margin_laws_at_fixed_windows) (Volume III, Problem 2.2, page 24; open): The book asks for proofs at fixed windows and bounded heights of the margin laws seen in the runs: the least eigenvalue of Weil's form, the Thorin margin and the curve margin. New in this book.
- [Problem: strip shrinking for pair couplings and the multiplicity rate](../RESULTS.md#question_quantitative_strip_shrinking_beyond_curie_weiss) (Volume III, Problem 2.3, page 24; open): The book asks for a shrinking bound for pair couplings that do not act through the magnetisation and for the rate delta squared minus N a for N identical copies. New in this book.
- [Problem: the gap law for clusters as wide as the spacing of their neighbours](../RESULTS.md#question_gap_law_for_many_zero_pairs) (Volume III, Problem 2.4, page 24; open): The book asks for the gap law of the one variable ferromagnetism test when clusters of zeros are as wide as the spacing of their neighbours, where a certificate for the tilt excludes at most the zeros with real part above 1/2 + w/cos(pi/4) unless simplicity of the zeros is also certified. New in this book.
- [Problem: a family along the Euler product with propagating Lee and Yang property](../RESULTS.md#question_family_along_the_euler_product) (Volume III, Problem 2.5, page 25; open): The book asks for an integer indexed family of laws or spins built from the Euler product of zeta whose Lee and Yang property or Toeplitz positivity propagates from one member to the next and whose limit is the Pólya spin. New in this book.
- [Problem: the two spin Lee and Yang property for every w above some w_0 below one half](../RESULTS.md#question_two_spin_lee_yang_below_one_half_coupling) (Volume III, Problem 2.6, page 25; open): The book asks for a proof of the two spin Lee and Yang property for every w above some w_0 below one half, which is a zero free half plane Re s above one half plus w_0, the quasi Riemann hypothesis, proved at w_0 equal to 3/8 by OpenAI (2026). New in this book.
<!-- results:end -->

## Key terms

Riemann Hypothesis, RH, route to RH, filters, zero free half plane, quasi Riemann hypothesis, heat flow, Pólya frequency function, Hankel determinants, Stieltjes continued fraction, closed roads, Thorin measure, descent in level, backward Cauchy flow, deconvolution, Weil positivity, Weil criterion, Lee and Yang theorem, ferromagnet, two spin system, Pólya law, Laguerre Pólya class, Hilbert Polya conjecture, Hasse Weil, curves over finite fields, Hodge index theorem, Toeplitz test, Hankel test, Davenport and Heilbronn function, ball arithmetic, de Bruijn Newman constant
