# Volume I, Chapter 4. Infinite divisibility, GGC, EGGC and the Bondesson hierarchy

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 4. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume I: Representations of Xi, Thorin Measures and the Reciprocal Xi Function*, edition 2, https://doi.org/10.5281/zenodo.23861971. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter sets out the class of generalised gamma convolutions (GGC) in which the Thorin condition for the Riemann Hypothesis (RH) lives, and places it among the infinitely divisible laws. A random variable $H\ge0$ is a GGC when its Laplace transform is $\exp\{-du-\int\log(1+u/z)\,U(dz)\}$ with drift $d\ge0$ and a Thorin measure $U$. The chapter gives three readings of this definition: through the Lévy measure and a completely monotone canonical function $k$, through the Stieltjes transform and Stieltjes inversion, and as a product over the atoms of $U$. For the clock $1/G$, with $G(s^2)=\xi(\frac12+s)/\xi(\frac12)$ built from the completed zeta function, it identifies the Thorin condition with complete monotonicity of the heat trace and writes the Thorin density at resolution $\varepsilon$ as a regularised Stieltjes inversion along a parabola. It then lays out the Bondesson hierarchy: HCM, GGC, self decomposable (SD), the Bondesson class (BO) and infinitely divisible (ID) laws, with Bondesson's theorem that a transform is that of a GGC exactly when it is hyperbolically completely monotone, extended GGCs (EGGC) as Gaussian mixtures over a GGC time, and Urbanik's classes of multiply self decomposable laws. A table records which sign conditions on the heat trace are proved for the clock (with the zeros verified to a stated height) and which are equivalent to RH and open. A second table lists diagnostics for the Thorin condition and what a failure of each proves, and the chapter closes with explicit Thorin measures and the sinh benchmark, where distances to the lattice of atoms are computed in closed form.

## From the literature

- generalised gamma convolutions and the Thorin measure (Thorin 1977; Bondesson 1992)
- extended generalised gamma convolutions and the representation of a symmetric EGGC as $Z\sqrt{2H}$ (Thorin 1978; Bondesson 1992)
- Bondesson's theorem that a Laplace transform with value one at the origin belongs to a GGC exactly when it is hyperbolically completely monotone (Bondesson 1992)
- Frullani's integral, Bernstein's theorem and Stieltjes inversion, as they enter the Lévy and Stieltjes readings (James, Roynette and Yor 2008)
- the classes SD and BO and the examples separating them (Polson 2017)
- Urbanik's classes of multiply self decomposable laws
- the hyperbolic laws as Brownian exit times (Biane, Pitman and Yor 2001)
- the Thorin condition as a reformulation of RH, and the diagnostics for it (Polson 2017, 2026)

## New results

Statements are proved unless a table marks them as computed or equivalent to RH, in which case RH is not assumed; the tables of distances for the sinh clock are closed forms.

<!-- results:begin -->
- [The variable Z sqrt(2T) of the van Dantzig pair is self decomposable with no hypothesis, with an explicit canonical function](../RESULTS.md#self_decomposability_of_z_sqrt_2t_and_its_canonical_function_for_the_clock) (Volume I, Result 4.1.2, page 33; proved): With no hypothesis, $Z\sqrt{2T}$ (standard normal $Z$ independent of the clock $T$) is self decomposable, with canonical function $k(r)=\sum_je^{-r\sqrt{a_j}}$ obtained by subordinating the Lévy measure of $T$ to the heat kernel; it is an extended generalised gamma convolution exactly when RH holds. New in this book.
<!-- results:end -->

## Key terms

generalised gamma convolution, generalized gamma convolution, GGC, Thorin measure, Thorin condition, infinite divisibility, infinitely divisible, self decomposable, Bondesson class, hyperbolically completely monotone, HCM, extended generalised gamma convolution, EGGC, Urbanik class, Lévy measure, canonical function, completely monotone, Stieltjes transform, Stieltjes inversion, Poisson smoothing, heat trace, clock, sinh benchmark, Riemann Hypothesis, RH
