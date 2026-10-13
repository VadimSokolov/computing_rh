# Chapter 3. The Thorin measure at finite resolution

<!-- cite:begin -->
From *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 3. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1, https://doi.org/10.5281/zenodo.23331091. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter asks what finite computations of the Thorin density $\rho_\varepsilon$, the Poisson smoothing of the zeros of $\xi$, can establish about the Riemann Hypothesis. On the window below height 200, where the 79 zeros are known, it measures the distance from the computed Thorin measure to the discrete measure demanded by RH. Raw total variation is maximal, matched total variation sits at a numerical floor of about $3.7\times10^{-9}\varepsilon$, and cell distance, mass defect and heat trace error vanish linearly in the resolution $\varepsilon$. The chapter proves the first order mass defect, closed form signatures of an off line zero, a Littlewood identity, and a zero density identity that counts the zeros to the right of any line from the cumulative Thorin mass (at $\varepsilon=\frac12$, Backlund's count with the vertical side on $\operatorname{Re}s=1$). It shows that the sine squared lemma and Selberg's smoothed explicit formula give the Thorin density from the primes only at a cost exponential in $1/\varepsilon$, and explains why a finite computation can refute RH by certifying a zero off the line but cannot measure a distance to RH.

## From the literature

Lagarias's criterion (Lagarias 1999); the argument principle count, the Riemann von Mangoldt formula and Littlewood's lemma (see Titchmarsh 1986); Backlund's method (Backlund 1918) and Turing's method (Turing 1953; Lehman 1970; Trudgian 2011); Selberg's smoothed explicit formula (Selberg 1946); Weil's positivity criterion (Weil 1952); zero density estimates (Ingham 1940); the verification of RH to height $3\times10^{12}$ (Platt and Trudgian 2021). The sine squared lemma is an earlier result of Polson (2017); the canonical measure of $\log\xi$ that the sine squared decomposition uses and the square root law for a displaced zero are earlier results of Polson (2026).

## New results

The propositions below are new in the book. The tables of distances, masses and counts at height 200 are computed, not certified.

<!-- results:begin -->
- [Mass defect of the Thorin measure on a window](../RESULTS.md#thorin_mass_first_order_defect) (Proposition 3.2.2, page 58; proved): When every zero up to height $\Theta$ lies on the critical line, the Thorin mass of the window $(0,\Theta^2]$ equals the zero count $N(\Theta)$ minus $(\varepsilon/\pi)\,\Xi'(\Theta)/\Xi(\Theta)$ up to a remainder of order $\varepsilon^2$. New in this book.
- [Signatures of an off line zero in the Thorin density](../RESULTS.md#off_line_zero_signatures) (Proposition 3.3.1, page 61; proved): The effect of an off line pair of zeros on the Thorin density is computed in closed form: a displacement that keeps the count has total variation $\frac{3\sqrt{3}}{4\pi}(\delta/\varepsilon)^2$ per displaced zero to leading order for $\delta<\varepsilon$, an omitted pair has total variation $\frac4\pi\arcsin(\varepsilon/\delta)$, and a displaced zero has negative mass $\frac1\pi\arcsin(\varepsilon/\delta)$, the last two for $\delta>\varepsilon$; the quadratic growth of the first is the square root law of Polson (2026). New in this book.
- [Littlewood identity for the cumulative Thorin mass](../RESULTS.md#littlewood_identity_thorin_mass) (Proposition 3.4.1, page 66; proved): The cumulative Thorin mass at resolution $\varepsilon$ is $\pi^{-1}$ times the derivative in $\varepsilon$ of the integral of $\log|\xi|$ along the line $\operatorname{Re}s=\frac12+\varepsilon$, which links the mass to the integral of $\log|\xi|$ on the critical line. New in this book.
- [Zero density identity from the Thorin mass](../RESULTS.md#zero_density_identity_thorin_mass) (Proposition 3.4.2, page 67; proved): The cumulative Thorin mass at resolution $\varepsilon$ plus a short horizontal correction equals the zero count $N(T)$ minus twice the number of zeros to the right of the line $\operatorname{Re}s=\frac12+\varepsilon$, and at $\varepsilon=\frac12$ it is Backlund's unconditional count of $N(T)$. New in this book.
<!-- results:end -->

## Key terms

Thorin density, Thorin mass, resolution, Poisson smoothing, Cauchy kernel, total variation, matched total variation, cell distance, mass defect, Richardson extrapolation, off line zero, square root law, Littlewood identity, zero density, argument principle, Backlund count, Turing's method, Selberg explicit formula, sine squared lemma, variogram, Weil criterion
