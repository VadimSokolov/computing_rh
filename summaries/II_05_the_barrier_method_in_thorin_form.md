# Volume II, Chapter 5. The barrier method in Thorin form

<!-- cite:begin -->
From *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation* by Nicholas G. Polson and Vadim Sokolov (2026), chapter 5. Cite as Polson and Sokolov (2026), *Computing the Riemann Hypothesis, Volume II: Cauchy Deconvolution and Certified Computation*, edition 2, https://doi.org/10.5281/zenodo.23891946. Book page: https://vsokolov.org/books/computing_rh/. Source: https://github.com/VadimSokolov/computing_rh.
<!-- cite:end -->

This chapter states the barrier method of Polymath15, as used by Platt and Trudgian, translates it into the language of Thorin measures, and runs a small barrier check. The method proves upper bounds for the de Bruijn Newman constant $\Lambda$ by combining a verification of the Riemann Hypothesis (RH) to a height, an analytic zero free region for $H_{t_0}$ at large heights, and a barrier, a window of heights through which no zero of $H_t$ can pass for any time $t\in[0,t_0]$. It gave $\Lambda\le0.22$ (Polymath15) and $\Lambda\le0.2$ (Platt and Trudgian). The criterion is stated as Polymath15 print it, with de Bruijn's strip theorem and the canopy. In Thorin language its conclusion is that the shifted clock at $\alpha_0=\frac12+\frac{y_0}2$ is a generalised gamma convolution (GGC) transform, so the reciprocal of $H_t$ is a GGC transform for every $t\ge t_0+\frac12y_0^2$ (proved); the barrier becomes the statement that no Thorin atom crosses a window bounded by four confocal parabolas with focus at the origin (proved), and a positive density $\operatorname{Re}(\xi_t'/\xi_t)$ on the barrier rectangle at every time is a sufficient certificate (proved). A check at height about $497.4$ for $t\in[0,0.2]$ and $y_0=0.2$ is run in floating point (computed) and certified in ball arithmetic for every $t\in[0,0.2]$. The Thorin form adds a geometry and a link to the phase diagram but no analytic input, and the method cannot reach $t=0$: every certified point lies on a level set $t+2(\alpha-\frac12)^2=c$ with $c>0$, and the corner of RH needs a zero free region for $\zeta$ itself.

## From the literature

- the barrier criterion for upper bounds on $\Lambda$ (Polymath15 2019)
- $\Lambda\le0.2$ by a barrier at larger height (Platt and Trudgian 2021)
- de Bruijn's strip theorem (de Bruijn 1950)
- $\Lambda\ge0$ (Rodgers and Tao 2020)

## New results

<!-- results:begin -->
No new results of this chapter are catalogued.
<!-- results:end -->

## Key terms

barrier method, Polymath15, de Bruijn Newman constant, Lambda, heat flow, H_t, canopy, strip theorem, Thorin measure, generalised gamma convolution, GGC, shifted clock, confocal parabolas, phase diagram, Platt and Trudgian, ball arithmetic, zero free region, Riemann Hypothesis, RH
