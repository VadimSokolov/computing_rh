# Chapter 23. Haseo Ki and the Thorin programme

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 23. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23595861. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

What does the work of Haseo Ki and his coauthors say about the Thorin condition on the reciprocal of the Riemann Xi function? Ki and his coauthors studied Fourier transforms of rapidly decaying positive kernels, their zeros under the heat flow (the de Bruijn Newman flow) and under repeated differentiation, and zeta functions with the functional equation of $\zeta$ but no Euler product. The chapter reads this through the Thorin dictionary: a zero of an even entire function is an atom of the Thorin measure of its reciprocal, and the reciprocal is a generalized gamma convolution (GGC) exactly when the zeros are real. It proves that after any positive time the reciprocal of the flowed function is a GGC transform up to a finite rational factor (the theorem of Ki, Kim and Lee), with the atoms of the factor confined to confocal parabolas that shrink with time. With the theorem of Rodgers and Tao this gives an obstruction: no property of the Pólya density that survives a Gaussian tilt can imply RH. For repeated differentiation, the rescaled reciprocals converge to the Laplace transform of the exit time of Brownian motion from an interval, a GGC, and a conjecture gives the rate. Epstein type zeta functions and the Davenport Heilbronn functions indicate that a Thorin argument must use the Euler product. Nothing here proves RH.

## From the literature

- the heat flow on Fourier transforms of positive kernels, the de Bruijn Newman constant, and the zeros of $\Xi_t$ (de Bruijn 1950; Newman 1976; Ki and Kim 2002, 2003; Ki, Kim and Lee 2009)
- the de Bruijn Newman constant is nonnegative (Rodgers and Tao 2020)
- Newman's characterisation of the Lee and Yang single spin laws by the Laguerre Pólya class (Newman 1974)
- Ki's theorem on repeated differentiation, the zeros of $\zeta'$, and approximation theorems for Epstein type zeta functions and the Davenport Heilbronn functions (Ki and coauthors)

## New results

Statements are proved unless marked: the Edgeworth rate is an open conjecture supported by computation, and the remark on what a Thorin proof must use is a heuristic and is not catalogued.

<!-- results:begin -->
- [Eventually Thorin: after any positive time the reciprocal of the flowed Xi function is a GGC transform times a finite rational factor](../RESULTS.md#after_any_positive_time_the_reciprocal_of_the_flowed_function_is_a_ggc_transform_up_to_a_rational_factor) (Proposition 23.2.1, page 794; proved): For every $t>0$ the ratio $\Xi_t(0)/\Xi_t(i\sqrt s)$ equals a generalised gamma convolution Laplace transform times a finite product $\prod_{\eta\in E_t}(1+s/\eta)^{-1}$ over a conjugation invariant multiset off the positive axis, which is empty exactly when $\Xi_t$ has only real zeros. New in this book.
- [Shrinking parabolas: Thorin atoms of the flowed Xi function lie to the right of a parabola that collapses at t = 1/2](../RESULTS.md#thorin_atoms_of_the_flowed_function_lie_right_of_a_shrinking_parabola) (Corollary 23.2.2, page 794; proved): For $0<t<\frac12$ every Thorin atom of $\Xi_t$ lies in the closed region to the right of the parabola $x=y^2/(4\Delta^2)-\Delta^2$ with $\Delta^2=\frac14-\frac t2$; the parabolas share the focus $0$ and collapse onto $[0,\infty)$ at $t=\frac12$. New in this book.
- [No slack in the Gaussian direction: no class of kernels closed under Gaussian damping and containing Phi has only real zeros](../RESULTS.md#no_property_closed_under_a_gaussian_tilt_of_the_polya_density_can_imply_rh) (Proposition 23.2.3, page 795; proved): If a class of positive even kernels is closed under multiplication by $e^{-bu^2}$ for every $b>0$ and every member has a Fourier transform with only real zeros, then Riemann's kernel $\Phi$ is not in the class; this follows from the theorem of Rodgers and Tao that the de Bruijn Newman constant is nonnegative. New in this book.
- [Brownian exit time limit: the reciprocals of the even derivatives of Xi, rescaled, converge to 1 over cosh of the square root of 2s](../RESULTS.md#reciprocals_of_derivatives_of_xi_converge_to_the_laplace_transform_of_the_brownian_exit_time) (Proposition 23.3.1, page 796; proved): For $\psi_n(s)=\Xi^{(2n)}(0)/\Xi^{(2n)}(i\sqrt{2s}/m_n)$, $\psi_n(s)\to1/\cosh\sqrt{2s}$ locally uniformly off the poles, the Laplace transform of the exit time of Brownian motion from $(-1,1)$, whose Thorin atoms $c_k^2/2$ are approached by the smallest atoms of $\psi_n$; under RH each $\psi_n$ is a GGC transform and the clocks converge in law. New in this book.
- [Conjecture: the Edgeworth rate of convergence of the zeros and clocks of the derivatives of Xi](../RESULTS.md#conjecture_on_the_edgeworth_rate_of_convergence_of_the_derivative_clocks) (Remark 23.3.3, page 797; open): The $k$th positive zero of $\Xi^{(2n)}$ is conjectured to satisfy $z_{k,n}m_n=c_k+\kappa_nc_k^3/(6m_n^3)(1+o(1))$ and the distance of $\psi_n$ from $1/\cosh\sqrt{2s}$ is $\sup_{x\le\sqrt{2S}}\frac{x^2}{2\cosh x}\cdot\sigma_n^2/m_n^2\,(1+o(1))$; the constant $0.53222$ is computed. New in this book.
<!-- results:end -->

## Key terms

Haseo Ki, de Bruijn Newman flow, heat flow, Thorin measure, generalised gamma convolution, generalized gamma convolution, GGC, eventually Thorin, confocal parabolas, Rodgers and Tao, Gaussian tilt, derivative of Xi, Brownian exit time, Edgeworth rate, Davenport Heilbronn function, Epstein zeta function, Euler product, Laguerre Pólya class
