# Chapter III.15. Maślanka, Haglund and the missing positivity

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 15. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Can the positive Maslanka sequence be written for the xi function itself, and does Haglund's incomplete gamma identity give a route to the Thorin form? The chapter defines g_k, the binomial differences of xi'/xi at even integers, shows that it is a pure Beta kernel sum over the zeros and equals minus the remainder R_k minus 1/(2k(k+1)), and proves a one sided Landau criterion: a zero free half plane Re s greater than theta is equivalent to either one sided bound on g_k of order k to the minus 1 plus theta over 2, with RH the case theta = 1/2 and the theorem of OpenAI (2026) giving theta = 7/8. The 2 sin squared lemma gives a margin of order k to the minus 1/2 where k to the minus 3/4 is needed. The second half writes the Polya density as a positive mixture of Haglund components, studies the kink at the origin, shows by a tolerance theorem that complete monotonicity of the xi clock survives a proportion of off line zeros, and gives the Hankel form of the total positivity ladder with certified signs for planted off line pairs in a sixty zero model. Open problems are listed.

## From the literature

the Maslanka expansion (Maslanka 2001), the Baez Duarte criterion (Baez Duarte 2003, 2010), Haglund's approximants of Xi by incomplete gamma functions (Haglund 2009), Polya's approximant and integral representation (Polya 1926), Ingham's one sided theorem (Ingham 1932), Conrey's two fifths theorem (Conrey 1989), verification of RH to height 3e12 (Platt and Trudgian 2021), total positivity (Karlin 1968), Schoenberg's theorem on Polya frequency functions (Schoenberg 1951), the divided difference form (Curry and Schoenberg 1966), hierarchical mixtures (Polson and Scott 2016) and Li's criterion (Li 1997).

## New results



<!-- results:begin -->
- [Maslanka sequence of xi as a sum over zeros](../RESULTS.md#maslanka_sequence_of_xi_zero_sum) (Volume III, Theorem 15.1.1, page 243; proved): The binomial differences g_k of the logarithmic derivative of xi at even integers equal a Beta kernel sum over the zeros and minus R_k minus 1 over 2k(k+1), with g_k = O(k to the minus 9/16) unconditionally. New in this book.
- [One sided Landau criterion for the Maslanka sequence of xi](../RESULTS.md#landau_criterion_maslanka_sequence_xi) (Volume III, Theorem 15.2.1, page 245; proved): A zero free half plane Re s greater than theta is equivalent to either one sided bound g_k at most or at least plus or minus k to the minus 1 plus theta over 2 plus epsilon, and RH is the case theta equal 1/2. New in this book.
- [The 2 sin squared lemma gives only a k to the minus 1/2 margin](../RESULTS.md#sin2_lemma_margin_maslanka_sequence) (Volume III, Proposition 15.3.1, page 246; proved): The nonnegativity d_k at least Re D_k(2+iq) holds with margin of order k to the minus 1/2, while RH needs k to the minus 3/4, so the 2 sin squared lemma is too weak by a factor k to the 1/4. New in this book.
- [Polya density as a positive mixture of Haglund components](../RESULTS.md#haglund_identity_positive_mixture) (Volume III, Proposition 15.4.1, page 247; proved): The Polya law with characteristic function Xi(t)/Xi(0) is a mixture of symmetrised laws whose densities are nonnegative components phi_n of the Polya kernel, with weights p_n given by the transforms of the components at zero. New in this book.
- [One real rate absorbs a complex pair in complete monotonicity](../RESULTS.md#one_real_rate_absorbs_complex_pair) (Volume III, Proposition 15.7.1, page 249; proved): The reciprocal of (1+u/r)(1+u/a)(1+u/conjugate a) is completely monotone exactly when the real rate r is at most the real part of a. New in this book.
- [Complete monotonicity of the xi clock tolerates a proportion of off line pairs](../RESULTS.md#complete_monotonicity_tolerates_off_line_pairs) (Volume III, Theorem 15.7.2, page 249; proved): If at every height T the number of off line pairs is at most the number of zeros on the line up to height T minus 1, then xi(1/2) over xi(1/2 plus sqrt u) is completely monotone, so existence of the Wald partner cannot separate zeta from a function that only behaves like it. New in this book.
- [Hankel determinants of a finite product of rates and their signs](../RESULTS.md#hankel_determinants_finite_product_asymptotics) (Volume III, Proposition 15.8.1, page 251; proved): For a finite product with distinct rates the Hankel determinant H_r is an explicit sum over r subsets of rates, and as x tends to infinity it oscillates in sign when a complex pair follows r minus 1 real rates. New in this book.
<!-- results:end -->

## Key terms

Maslanka sequence, Baez Duarte criterion, Landau criterion, one sided criterion, Haglund approximants, incomplete gamma function, Polya density, Wald couple, van Dantzig pair, complete monotonicity, Thorin form, kink, total positivity, Hankel determinant, planted off line zeros, Xi function
