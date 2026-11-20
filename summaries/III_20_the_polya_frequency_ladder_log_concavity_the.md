# Chapter III.20. The Pólya frequency ladder: log concavity, the detection rung and the heat flow

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 20. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Is the Riemann density a Pólya frequency function? The chapter takes the total positivity formulation of the Riemann Hypothesis (RH), in which the Riemann density F is a Pólya frequency function exactly when RH holds, and asks what can be proved and computed about each rung PF_k. It proves the far field of rung k, for each fixed configuration, from the first k simple zeros, certifies in Arb that F is strictly log concave (PF_2), proves a detection rung theorem (an off line quadruple above m simple real zeros breaks rung m+1) and checks it on the Davenport Heilbronn function, and embeds F in the de Bruijn Newman heat flow F_t, where for t below 0 some rung fails and an upper bound on the first failing rung is computed from the zero flow and Lehmer pairs of zeros. It also proves positivity of the Riemann density without computation, from a lower bound on the Lévy density, complete monotonicity of K_{1/2} to order 10^11 given the verification of RH to height 3 times 10^12, and an unconditional Thorin density of the half product.

## From the literature

- the Schoenberg characterisation of Pólya frequency functions (Schoenberg 1951) and total positivity (Karlin 1968)
- strong unimodality of log concave laws (Ibragimov 1956)
- the heat flow, the de Bruijn Newman constant and its characterisation of RH (de Bruijn 1950; Newman 1976)
- the bounds 0 at most Lambda at most 0.2 (Rodgers and Tao 2020; Polymath 2019; Platt and Trudgian 2021)
- verification of RH up to height 3 times 10^12 (Platt and Trudgian 2021)
- the Csordas Smith Varga bound from Lehmer pairs (Csordas, Smith and Varga 1994)
- the Davenport Heilbronn function as a control (Davenport and Heilbronn 1936)
- the Chebyshev coefficient bound for Gauss quadrature (Trefethen 2019)
- Arb ball arithmetic (Johansson 2017)

## New results

The theorems below are the authors' (the approach is from an unrefereed note of Polson, November 2026; the certificate was rebuilt and rerun).

<!-- results:begin -->
- [Residues of the Riemann density alternate in sign](../RESULTS.md#residues_of_riemann_density_alternate_in_sign_on_first_simple_zeros) (Volume III, Proposition 20.2.1, page 366; proved): If the first K positive zeros of Xi are real and simple, the residues c_n of the reciprocal Xi expansion alternate in sign, starting with c_1 positive. New in this book.
- [Far field of rung k of the Riemann density for each fixed configuration](../RESULTS.md#far_field_of_every_rung_for_each_fixed_configuration_of_riemann_density) (Volume III, Theorem 20.2.2, page 366; proved): If the first k zeros of xi are simple and on the line, every fixed k by k minor of the Riemann density is positive after a large enough shift, which is unconditional for k up to about 10^13 by the verification of RH. New in this book.
- [The Riemann density is strictly log concave (certified)](../RESULTS.md#riemann_density_strictly_log_concave_certified) (Volume III, Theorem 20.4.1, page 368; proved): The Riemann density F, with characteristic function xi(1/2)/xi(1/2+t), is strictly log concave on the whole real line, a result certified in ball arithmetic that settles the PF_2 part of the Pólya frequency question without assuming RH. New in this book.
- [Positivity of the Riemann density without computation](../RESULTS.md#riemann_density_positive_from_first_zero_bound_on_levy_density) (Volume III, Proposition 20.5.1, page 373; proved): The function K_{1/2}, for which K_{1/2}(|x|)/|x| is the Lévy density of the Riemann variable, is at least (1 minus 10^-7) times exp(minus gamma_1 x), so the Riemann variable contains a symmetrised gamma component and its density is positive. New in this book.
- [Complete monotonicity of the Lévy density to order 10^11](../RESULTS.md#levy_density_of_riemann_variable_completely_monotone_to_order_10_11_given_verified_zeros) (Volume III, Theorem 20.5.2, page 374; proved): Given the verification that every zero of height up to 3 times 10^12 lies on the critical line, the first 10^11 derivatives of K_{1/2} alternate in sign on the whole positive half line. New in this book.
- [Thorin density of the half product without hypothesis on the zeros](../RESULTS.md#unconditional_thorin_density_of_half_product_of_xi) (Volume III, Proposition 20.5.3, page 375; proved): The half product over zeros in the upper half plane has a Laplace transform with a positive Thorin density D(theta)/pi, with no hypothesis on the zeros, so positivity of this Thorin measure is free while positivity for the clock is RH. New in this book.
- [Detection rung theorem](../RESULTS.md#detection_rung_theorem_off_line_quadruple_breaks_next_rung) (Volume III, Theorem 20.6.1, page 376; proved): A single simple quadruple of zeros off the line, lying above m simple real zeros, makes the associated density fail Pólya frequency rung m+1 while rungs up to m hold in the far field. New in this book.
- [Detection rung for the Riemann density](../RESULTS.md#riemann_density_not_pf_m_plus_1_if_lowest_off_line_height_is_a_simple_quadruple) (Volume III, Corollary 20.6.2, page 377; proved): If RH fails with the lowest height off the line carrying one simple quadruple above m simple zeros on the line, the Riemann density is not PF_{m+1}, equivalently PF_k forces m at least k. New in this book.
- [The heat flow densities F_t are Pólya frequency functions exactly for t at least Lambda](../RESULTS.md#heat_flow_density_is_polya_frequency_function_iff_t_at_least_de_bruijn_newman_constant) (Volume III, Theorem 20.7.1, page 379; proved): The deformed density F_t is a Pólya frequency function if and only if t is at least the de Bruijn Newman constant, so every rung holds for t at least 0.2, RH is the case t equals 0, and every F_t with t below 0 fails some finite rung. New in this book.
- [Critical point test for nonreal zeros](../RESULTS.md#critical_point_test_for_zeros_off_the_real_axis_of_heat_flow_function) (Volume III, Proposition 20.7.2, page 381; proved): A real critical point of H at which the modulus has a local minimum and H is nonzero proves that H has zeros off the real axis. New in this book.
- [Detection scale, sign changes of the failing minor are spaced by pi over b](../RESULTS.md#detection_scale_sign_changes_of_minor_spaced_by_pi_over_b) (Volume III, Corollary 20.7.4, page 385; proved): In the setting of the detection rung theorem the failing minor of order m+1 equals a decaying exponential times a cosine of frequency b, so its sign changes are asymptotically spaced by pi over b. New in this book.
<!-- results:end -->

## Key terms

Pólya frequency function, Polya frequency function, PF_k, PF_infinity, totally positive, total positivity, log concavity, log concave, strongly unimodal, Riemann density, reciprocal xi, characteristic function, ladder of rungs, far field, near field, residue expansion, detection rung, off line zeros, Lehmer pair, de Bruijn Newman constant, heat flow, backward heat equation, Laguerre Pólya class, Davenport Heilbronn, ball arithmetic, Arb, certified, Lévy density, Levy density, complete monotonicity, Thorin density, generalised gamma convolution, generalized gamma convolution, Riemann Hypothesis, RH
