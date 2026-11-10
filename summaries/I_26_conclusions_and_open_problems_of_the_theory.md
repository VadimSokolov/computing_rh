# Volume I, Chapter 25. Conclusions and open problems of the theory

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 26. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function*, edition 2, https://doi.org/10.5281/zenodo.24086101. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

How far does the Thorin approach to the Riemann Hypothesis (RH) reach, and what can computation add? The chapter summarises the theory and sets out a path along the basepoint $\alpha$ of the shifted clock of the reciprocal Xi function, with Laplace transform $\xi(\alpha)/\xi(\alpha+\sqrt u)$. For $\alpha\ge1$ it is a generalised gamma convolution (GGC) with positive Thorin measure, with no hypothesis. For $\frac12<\alpha<1$ the GGC property is equivalent to a zero free half plane $\operatorname{Re}s>\alpha$, and at the centre it is RH (Lagarias positivity). The gap is one positivity statement uniform in the height and the resolution. The chapter discusses the zero free half plane $\operatorname{Re}s>\frac78$ announced in October 2026 and not refereed, the sine squared lemma, five computational work packages (among them certified Thorin densities in ball arithmetic, twisted Riesz filters, margins against the circular unitary ensemble, the almost GGC defect), closed routes, the reach of finite information (a perturbation theorem with off line zeros) and what a proof must supply. It assembles the route from the whole theory and states the missing step as a single positivity problem. It poses open problems from free probability, on Jensen polynomials, Laguerre expansions, Jacobi parameters, the Hankel test and zero collisions along the de Bruijn Newman flow, and the heat trace ladder, and sorts the open problems and lessons of the Thorin programme. Nothing here proves RH.

## From the literature

- The Thorin representation of generalised gamma convolutions (Thorin 1977; Bondesson 1992)
- The criteria of Riesz (1916), Weil (1952), Li (1997; Keiper 1992) and Lagarias (1999)
- The de Bruijn Newman constant (de Bruijn 1950; Rodgers and Tao 2020; Platt and Trudgian 2021)
- Conrey and Li (2000) on a false positivity condition
- Zero density estimates (Ingham 1940; Huxley 1972; Guth and Maynard 2026)
- RH verified to height $3\times10^{12}$ (Platt and Trudgian 2021); Odlyzko's tables of zeros (1987; 2001)
- The finite free convolution (Marcus, Spielman and Srivastava 2022) and the Upsilon transform (Barndorff-Nielsen and Thorbjørnsen 2004)
- The half plane $\operatorname{Re}s>\frac78$ is announced by OpenAI (2026) and is not refereed, so results using it are conditional

## New results

This chapter summarises results proved earlier; its catalogued entries are open problems, among them the range of orders of the derivative signs of the heat trace.
<!-- results:begin -->
- [Jensen polynomials of the Pólya moments as finite free products (open problem)](../RESULTS.md#open_jensen_polynomials_as_finite_free_products) (Volume I, Problem 14.1, page 352; open): Find factorisations of the normalised coefficients of the Jensen polynomials of the Pólya moment sequence into two real rooted polynomials, which would give real zeros, or show that none exist with factors of a given kind. New in this book.
- [Laguerre and Meixner Pollaczek expansions of the Pólya density (open problem)](../RESULTS.md#open_laguerre_meixner_pollaczek_expansions_of_polya_density) (Volume I, Problem 14.2, page 353; open): Find a variable in which the Pólya density or the Xi function has an expansion in Laguerre or Meixner Pollaczek polynomials with positive coefficients and truncations with only real zeros. New in this book.
- [Jacobi parameters of the zero measure and the local semicircle model (open problem)](../RESULTS.md#open_jacobi_parameters_of_zero_measure_and_local_semicircle_model) (Volume I, Problem 14.3, page 353; open): Prove under the Riemann Hypothesis that the Jacobi parameters of the measure of the inverse squared ordinates satisfy alpha_n ~ 2 beta_n ~ L_n/2, a spectral form of the Riemann von Mangoldt formula, and explain the scatter of the ratios. New in this book.
- [Order at which a zero collision is registered by the Hankel test (open problem)](../RESULTS.md#open_order_at_which_a_zero_collision_is_registered_by_hankel_test) (Volume I, Problem 14.4, page 353; open): Find the rule giving the order of the Hankel matrix at which the collision of a pair of zeros along the de Bruijn Newman flow is registered, and follow the failure times of the Hankel test to orders of a few hundred. New in this book.
- [Unconditional brackets for the collision times of zeros along the de Bruijn Newman flow (open problem)](../RESULTS.md#open_unconditional_collision_brackets_for_de_bruijn_newman_flow) (Volume I, Problem 14.5, page 353; open): Find a priori control of the nonreal zeros of the heat evolved Xi function at negative times so that the book's brackets for collision times become statements computed from the Xi function at time zero alone. New in this book.
- [Extending the derivative sign ladder of the heat trace without new information on the zeros (open problem)](../RESULTS.md#open_moving_derivative_sign_ladder_of_heat_trace_without_new_zeros) (Volume I, Problem 14.6, page 353; open): Find a mechanism, for instance a functional equation for the law of the clock under the Upsilon transform, that extends the proved ranges of derivative signs of the heat trace without further information on the zeros. New in this book.
<!-- results:end -->

## Key terms

Riemann Hypothesis, RH, Thorin measure, generalised gamma convolution, generalized gamma convolution, GGC, shifted clock, reciprocal Xi function, basepoint, sine squared lemma, zero free half plane, de Bruijn Newman constant, Riesz criterion, Lagarias positivity, Jensen polynomials, ball arithmetic, open problems, path to RH
