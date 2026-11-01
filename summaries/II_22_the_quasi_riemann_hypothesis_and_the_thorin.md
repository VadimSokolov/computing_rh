# Volume II, Chapter 22. The quasi Riemann hypothesis and the Thorin programme

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 22. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation*, edition 2, https://doi.org/10.5281/zenodo.23861976. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

What are the zero free half planes announced in October 2026, who claims them with what status, how far could their arguments be followed here, and what would such a half plane do for the Thorin programme? The chapter records the claims: $\operatorname{Re}s>\frac78$ and $\frac{11}{12}$ for every Dirichlet $L$ function, zeta included, and for Hecke $L$ functions over $\mathbb Q(\sqrt{-3})$ (OpenAI 2026, with a Lean formalisation of the $\frac78$ statement by its authors), and the constant $0.874957\ldots$ of Liu, with later lower constants. All are unrefereed preprints, so every statement using them is conditional on the hypothesis named $(\mathrm Q_\theta)$, that no Dirichlet $L$ function has a zero with real part above $\theta$. It sets out the structure of the argument as far as read: a continuation criterion turning a power saving for smoothed Möbius sums into a half plane (whose short proof is checked), Poisson summation in a sextic twist, Kubota's cubic theta function and a large sieve type recursion, with an exponent bookkeeping that bounds one step at $\frac67$. Under $(\mathrm Q_\theta)$ the Thorin measure of reciprocal xi is positive and has its minimum at the origin at every level above $\theta$, the shifted clock being a generalised gamma convolution. The sine squared lemma adds nothing to the half plane because the canonical measure is signed. The descent in level from $\theta$ to $\frac12$ is analysed in four readings, none a reduction, and five open problems are posed. The chapter concludes that the announcement gives no path to the Riemann Hypothesis.

## From the literature

- Zero free half planes $\operatorname{Re}s>\frac78$ and $\frac{11}{12}$ (OpenAI 2026, unrefereed) and an absolute lower bound for $1-\beta$ at real zeros of real Dirichlet $L$ functions of conductor $q$, in the form $(1-\beta)\log q\ge c$ (OpenAI 2026)
- Improved constants near $0.875$ and a further lowering (Liu 2026; Xu and Sang 2026), unrefereed
- Littlewood's equivalence between bounds for the Mertens function and zero free half planes, in a smoothed form
- Kubota's cubic theta function and Poisson summation in twists
- The Thorin programme for reciprocal xi (Polson 2026), the Biane, Pitman and Yor laws, and the verification of RH to height $3\times10^{12}$ (Platt and Trudgian 2021)
- A report of the reception in the press (Quanta Magazine, 19 October 2026)

## New results

Everything that uses a half plane is conditional on $(\mathrm Q_\theta)$; the arithmetic checks of the quoted constants and the computed condition numbers are floating point or exact rational computations, not ball arithmetic.

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

quasi Riemann hypothesis, zero free half plane, seven eighths, eleven twelfths, Dirichlet L function, Hecke L function, Landau Siegel zero, Möbius sums, Mertens function, sextic twist, cubic theta function, Kubota, Poisson summation, large sieve, Lean formalisation, Thorin programme, Thorin measure, generalised gamma convolution, generalized gamma convolution, GGC, reciprocal xi, shifted clock, level, sine squared lemma, canonical measure, backward Cauchy flow, deconvolution, descent in level, OpenAI, Liu, Riemann Hypothesis, RH
