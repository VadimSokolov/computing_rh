# Chapter III.3. The toolkit and a solved model: the cosine

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 3. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume III: Is There a Path? Routes to RH from the Thorin Programme*, edition 2, https://doi.org/10.5281/zenodo.24565633. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter collects the short proofs on which the ladders of finite tests for the Riemann Hypothesis rest, and runs every tool on a solved model. The four languages for one statement about the zeros of the Riemann xi function are the density, the clock with its Thorin measure, the chain (Stieltjes continued fraction with birth and death rates) and the heat flow. The tools are the ladder of total positivity (Pólya frequency functions of order k, PF_k, with the Cauchy and Binet composition formula and the variation diminishing property), moments, Hankel determinants and continued fractions, and the heat flow of de Bruijn and Newman. The chapter proves that the density is log concave if and only if it is PF_2, that RH implies the Riemann density is PF_infinity (the easy half of Schoenberg's theorem), and states the detection rung theorem used to read a passed or failed rung. The worked example replaces xi by the cosine, where the density, Thorin measure, chain and heat flow have closed forms, and then plants a quadruple of zeros off the axis at height 3: the second log derivative of the density oscillates in sign and PF_2 fails, as the theorem predicts. The Davenport and Heilbronn control is introduced as the test that every argument must fail.

## From the literature

- Pólya frequency functions, totally positive kernels and the variation diminishing property (Schoenberg 1951; Karlin 1968)
- the Pólya frequency characterisation of RH and the Laguerre and Pólya class (Pólya 1926; Schoenberg 1951)
- the Cauchy and Binet composition formula for determinants
- generalised gamma convolutions and the Thorin measure (Thorin 1977; Bondesson 1992)
- the heat flow and the de Bruijn Newman constant (de Bruijn 1950; Newman 1976; Rodgers and Tao 2020)
- Lambert's continued fraction for the tangent

## New results

The catalogued result is the planted model, which is proved with computed constants.

<!-- results:begin -->
- [Cosine model with a planted quadruple of zeros fails log concavity](../RESULTS.md#planted_quadruple_cosine_model_fails_log_concavity) (Volume III, Proposition 3.7.1, page 43; proved): In a solved model where the cosine replaces the Xi function and four zeros are planted off the axis at height 3, the second log derivative of the density oscillates in sign with exact exponential rate and computed amplitude and phase, so the density is not PF2. New in this book.
<!-- results:end -->

## Key terms

total positivity, Pólya frequency function, PF_k, PF_infinity, log concavity, variation diminishing, Schoenberg theorem, Riemann density, clock, Thorin measure, Hankel determinants, Stieltjes continued fraction, birth and death chain, heat flow, de Bruijn Newman, detection rung theorem, cosine model, planted zeros, Lambert continued fraction, Davenport and Heilbronn control, generalized gamma convolution
