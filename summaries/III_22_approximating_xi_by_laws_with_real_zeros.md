# Chapter III.22. Approximating Xi by laws with real zeros

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 22. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

By Hurwitz's theorem RH would follow from functions with only real zeros converging locally uniformly to Xi(t) equal to xi(1/2+it). This chapter studies that idea in probabilistic form, with the laws of Volume II as approximants. The characteristic function Xi(t)/Xi(0) belongs to a law V, a tilted logarithm of the Biane Pitman Yor variable, with Wald partner H and van Dantzig partner. The chapter proves that under RH any independent factor of V has a canonical product over a split of the ordinates as characteristic function, and that one real zero of Xi can always be divided out and still leave a characteristic function, without hypothesis. It reports numerical evidence for splitting the zeros in halves, fits products of symmetrised generalised inverse Gaussian (GIG) laws with Bessel function characteristic functions to Xi on windows (each factor follows only two or three zeros), and proves the structure of the closure of such products: a Gaussian factor times a countable product of limits of single components, with the rigidity corollary that Xi in the closure forces RH and an exact split of the ordinates. This is a reduction and the expected failure of the route rests on an open problem, so it is a heuristic and not a theorem.

## From the literature

- Hurwitz's theorem and Hadamard factorisation of entire functions of order one
- the variable of Biane, Pitman and Yor with moments given by xi (Biane, Pitman and Yor 2001)
- van Dantzig pairs (Lukacs 1968)
- real zeros of Pólya's kernel and of Bessel type transforms (Pólya 1926)
- the Paley and Zygmund inequality and the Lyapunov central limit estimate
- monotonicity of Riemann's kernel (Wintner, source not cited in the book)
- Ingham's linear independence hypothesis (Ingham 1942) and Wald couples (Roynette and Yor 2005)

## New results

The removal of a zero and the closure results come from an unrefereed note of Polson (November 2026), and the fits are computed in double precision.

<!-- results:begin -->
- [Independent factors of the law V are canonical products over a split of the ordinates (under RH)](../RESULTS.md#factors_of_riemann_law_are_canonical_products_over_split_of_ordinates_under_rh) (Volume III, Proposition 22.2.1, page 413; conditional on The Riemann Hypothesis): Under RH, any decomposition of the law V with characteristic function Xi over Xi(0) into two independent nondegenerate summands splits the ordinates of the zeros into two infinite sets, each summand having the canonical product over its set as characteristic function. New in this book.
- [Removing a zero of Xi leaves a characteristic function](../RESULTS.md#one_real_zero_of_xi_can_be_divided_out_leaving_a_characteristic_function) (Volume III, Theorem 22.2.3, page 414; proved): Dividing the normalised Xi function by one factor 1 minus t squared over gamma squared, for any real zero gamma, still gives a characteristic function, with no hypothesis. New in this book.
- [Limits of single symmetrised GIG components](../RESULTS.md#limits_of_single_gig_component_normal_uniform_and_two_point_laws) (Volume III, Proposition 22.4.1, page 418; proved): The closure of the family of symmetrised generalised inverse Gaussian laws contains the normal law, the uniform law and the symmetric two point law. New in this book.
- [Structure of limits of sums of independent GIG components](../RESULTS.md#structure_of_limits_of_finite_sums_of_independent_gig_components) (Volume III, Theorem 22.4.3, page 419; proved): Every limit in law of finite sums of independent symmetrised GIG components has characteristic function a Gaussian factor times a countable product of limits of single components. New in this book.
- [Rigidity, if Xi over Xi(0) is in G then RH holds and the ordinates split into complete zero sets](../RESULTS.md#rigidity_if_xi_is_a_limit_of_products_of_gig_laws_then_rh_and_ordinates_split) (Volume III, Corollary 22.4.4, page 420; proved): If the normalised Xi function is a limit of products of GIG laws then RH holds and the ordinates of the zeros are a disjoint union of the complete zero multisets of nondegenerate members of the closure, each factor being the canonical product over its zeros. New in this book.
- [Open problem, is the list of limit laws of GIG components complete](../RESULTS.md#open_problem_completeness_of_list_of_limit_laws_of_gig_components) (Volume III, Problem 22.1, page 421; open): Prove that the listed limit laws (normal, products of cosine and Gaussian, truncated cosh laws) exhaust the closure of the family of symmetrised GIG laws and determine the zero set of each. New in this book.
- [Open problem, iterating the removal of zeros and the conjecture on splits](../RESULTS.md#open_problem_iterating_zero_removal_and_split_conjecture) (Volume III, Problem 22.2, page 421; open): Prove that the density left after removing one zero of Xi is again nonincreasing on the half line, so the removal can be iterated, and then prove the conjecture that sets of ordinates spread out in proportion at least theta_0 below one third give characteristic functions. New in this book.
- [Open problem, the ordinates of zeta are not a union of complete zero sets of GIG limit laws](../RESULTS.md#open_problem_ordinates_of_zeta_not_a_union_of_gig_zero_sets) (Volume III, Problem 22.3, page 421; open): Show unconditionally that the multiset of ordinates of the zeros of zeta is not a disjoint union of complete zero sets of nondegenerate members of the closure, which would prove that the normalised Xi function is not in G and close this route to RH. New in this book.
- [Open problem, GIG product fits on windows beyond 76](../RESULTS.md#open_problem_window_beyond_76_for_gig_product_fits) (Volume III, Problem 22.4, page 422; open): Settle the fit of products of GIG laws to Xi on the window up to 76 and beyond with a better optimiser, and measure how the number of components must grow with the window length. New in this book.
- [Open problem, the string of the Wald partner of a GIG component](../RESULTS.md#open_problem_string_of_wald_partner_of_gig_component) (Volume III, Problem 22.5, page 422; open): Identify the Krein string of the Wald partner of a symmetrised GIG component and compare its mass distribution with the string of H computed from the first 1545 zeros. New in this book.
<!-- results:end -->

## Key terms

Hurwitz theorem, Xi function, characteristic function with real zeros, ferromagnetic law, Lee Yang, generalised inverse Gaussian, GIG, symmetrised GIG, generalized inverse Gaussian, Bessel function K, Pólya approximant, Polya approximant, canonical product, Hadamard product, ordinates of zeros, van Dantzig pair, Wald partner, Thorin measure, splitting the zeros, removing a zero, limit in law, closure, rigidity, Riemann Hypothesis, RH
