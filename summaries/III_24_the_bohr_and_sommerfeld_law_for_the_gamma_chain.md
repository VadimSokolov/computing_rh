# Chapter III.24. The Bohr and Sommerfeld law for the gamma chain

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 24. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24377757. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

Can the Bohr and Sommerfeld law for the rates of the Gamma chain be proved? The Gamma chain is the birth and death chain, the Jacobi matrix, of the smooth zeros of the Riemann xi function; its rates were found numerically to follow c_j = 1/(4 g_j^2) with (g_j/2) log(g_j/4 pi) + c_0 = j and c_0 = 7/4. The chapter proves the leading order of the direct half: for a zero diagonal Jacobi matrix with regularly decreasing weights the eigenvalue counting function equals the arccos sum (semiclassical count) up to O(sqrt(n log n)), by cutting into blocks, freezing the weight, Weyl inequalities and an overlap lemma, and the chain of the law has the main term of the Riemann and von Mangoldt zero count (Lambert W function, Abel transform) with error O(sqrt(Gamma) log Gamma). The constant 7/4 is not derived and the death level correction is open. A conditional inverse half at leading order assumes the weights are decreasing and regularly varying, hypothesis (R). Computations check the midpoint law and a spacing weighted Gamma chain. Nothing here involves the primes or bears on RH; it is a statement about the Gamma factor. The law is from an unrefereed note of Polson (November 2026).

## From the literature

Weyl's inequalities for eigenvalue counts (Weyl 1912); semiclassical counting of eigenvalues of Jacobi matrices with varying recurrence coefficients (Kuijlaars and Van Assche 1999); the Lambert W function (Corless, Gonnet, Hare, Jeffrey and Knuth 1996); the Riemann and von Mangoldt zero count (Riemann 1859; von Mangoldt 1905); regular variation (Bingham, Goldie and Teugels 1987).

## New results

<!-- results:begin -->
- [Semiclassical eigenvalue count for Jacobi matrices with slowly decreasing weights](../RESULTS.md#semiclassical_count_jacobi_arccos_sum) (Volume III, Theorem 24.3.1, page 447; proved): For a zero diagonal Jacobi matrix with regularly decreasing weights, the number of eigenvalues above mu equals the arccos sum of the weights up to an error of order the square root of n log n. New in this book.
- [Direct half of the Bohr and Sommerfeld law, leading order](../RESULTS.md#direct_half_gamma_chain_law_leading_order) (Volume III, Theorem 24.4.2, page 449; proved): The chain defined by the Gamma chain law has the main term of the Riemann and von Mangoldt zero count, with error of order the square root of Gamma times log Gamma; the constant term is not identified. New in this book.
- [Conditional inverse half of the Gamma chain law, leading order](../RESULTS.md#conditional_inverse_gamma_chain_leading_order) (Volume III, Proposition 24.7.1, page 455; conditional on The weights w_j=sqrt(c_j) of the Gamma chain are decreasing and satisfy hypothesis (R), k/j <= (w_j-w_{j+1})/w_j <= K/j for j >= j_0): If the weights of the Gamma chain are decreasing and regularly varying in the sense of (R), then g_j/2 log g_j is asymptotic to j and c_j is asymptotic to (log j)^2/(16 j^2); the book gives a proof sketch only. New in this book.
<!-- results:end -->

## Key terms

Gamma chain, Bohr and Sommerfeld quantisation, Bohr and Sommerfeld quantization, semiclassical eigenvalue count, Jacobi matrix, Lambert W function, Abel transform, Riemann and von Mangoldt formula, smooth zeros, turning point, Plancherel and Rotach, Weyl inequality, birth and death chain, rates, inverse spectral problem
