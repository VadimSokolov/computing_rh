# Results by chapter <small>Computing the Riemann Hypothesis</small>

The results of [Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times](https://vsokolov.org/books/computing_rh/) that are new in the book, that come from the authors' own papers, or that correct one of those papers, as the book states them: 278 entries (203 proved, 11 conditional, 61 open, 3 counterexample). Results of other authors that the book presents are credited in the [chapter summaries](https://github.com/VadimSokolov/computing_rh/tree/main/summaries).

Status: proved; conditional (proved assuming the named statements, which the book does not prove); open (stated as an open question); counterexample (a proved counterexample); withdrawn (found wrong after release; the entry says what was wrong and which version corrects it).

## Chapter 1. Introduction and foundations

<a id="rigidity_clock_ggc_iff_rh"></a>

### Rigidity of the clock, a GGC if and only if RH

Proposition 1.2.1, page 12 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s - 1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed zeta function, an entire function whose zeros are the nontrivial zeros of $\zeta$, and put $G(u)=\xi(\frac12+\sqrt{u})/\xi(\frac12)$, an entire function of $u$. The clock is the function $1/G$. A random variable $H\ge 0$ is a generalised gamma convolution (GGC) if $\mathbb{E}\,e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, with a drift $d\ge 0$ and a positive measure $U$ on $(0,\infty)$, the Thorin measure, such that $\int_{(0,1)}|\log z|\,U(dz)+\int_{(1,\infty)}z^{-1}U(dz)<\infty$; the drift and the Thorin measure are unique. Then $1/G$ is the Laplace transform of a GGC if and only if the Riemann Hypothesis (RH) holds, that is, if and only if every zero of $\xi$ has real part $\frac12$. In that case the drift is zero, the Thorin measure is $U_{\mathrm{RH}}=\sum_{\gamma>0}\delta_{\gamma^2}$, where $\frac12+i\gamma$ runs over the zeros of $\xi$ with $\gamma>0$, counted with multiplicity, and the clock variable is $T=\sum_{\gamma>0}e_\gamma/\gamma^2$ with independent standard exponential random variables $e_\gamma$.

Source: [`ch/ch02.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch02.tex), label `prop:rigid`.

<a id="shifted_clock_ggc_right_of_every_zero"></a>

### Shifted clock is a GGC to the right of every zero

Theorem 1.2.2, page 17 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi$ be the completed zeta function, and for real $\alpha$ and real $r$ put $v_\alpha(r)=\operatorname{Re}(\xi'/\xi)(\alpha+ir)$. If $\alpha>\operatorname{Re}\rho$ for every zero $\rho$ of $\xi$, in particular for every $\alpha\ge 1$, then $u\mapsto\xi(\alpha)/\xi(\alpha+\sqrt{u})$ is the Laplace transform of a generalised gamma convolution $T_\alpha$ (a random variable $H\ge 0$ with $\mathbb{E}\,e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, a drift $d\ge 0$ and a positive Thorin measure $U$) whose Thorin measure has density $v_\alpha(\sqrt{z})/(2\pi\sqrt{z})$ with respect to $dz$ on $(0,\infty)$, and $v_\alpha>0$. The function $u\mapsto\xi(\alpha)/\xi(\alpha+\sqrt{u})$ is the Laplace transform of a generalised gamma convolution if and only if $\alpha\ge\sup_\rho\operatorname{Re}\rho$.

Source: [`ch/ch02.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch02.tex), label `thm:poisson`.

<a id="one_line_positivity_exclusion_open"></a>

### Positivity of the Thorin density on one line, open question

Remark 1.2.4, page 18 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed zeta function and, for real $\alpha>\frac12$ and real $r$, put $v_\alpha(r)=\operatorname{Re}(\xi'/\xi)(\alpha+ir)$, so that $v_\alpha(r)/\pi$ is the Thorin density at resolution $\alpha - \frac12$. If some zero of $\xi$ has real part greater than $\alpha$, the shifted clock $\xi(\alpha)/\xi(\alpha+\sqrt{u})$ is not the Laplace transform of a generalised gamma convolution. The question is whether, without assuming the Riemann Hypothesis, $v_\alpha(r)\ge 0$ for all real $r$ implies that $\xi$ has no zero with real part greater than $\alpha$. The book states that this is not known for $\xi$. What is known is that if $v_\alpha\ge 0$ for every $\alpha>\sigma_0$ then no zero has real part greater than $\sigma_0$ (for $\sigma_0=\frac12$ this is Lagarias's criterion for RH), and that for an explicit entire function with the order, symmetry and reality of $\xi$ the answer is negative.

Source: [`ch/ch02.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch02.tex), label `rem:ch2:oneline`.

## Chapter 2. Evaluating the Xi function

<a id="tilted_polya_representation_of_xi"></a>

### Tilted Pólya representation of xi and its derivative

Proposition 2.1.1, page 34 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi(u)=2\sum_{n\ge 1}\bigl(2\pi^2n^4e^{9u/2} - 3\pi n^2e^{5u/2}\bigr)\exp(-\pi n^2e^{2u})$ be the Pólya kernel, an even function that is analytic in the strip $|\operatorname{Im}u|<\pi/4$, and let $\xi$ be the completed zeta function, so that $\xi(\frac12+s)=\int_{\mathbb{R}}\Phi(u)e^{su}\,du$. Let $0<\delta<\pi/4$ and $\eta=\pi/4 - \delta$. Then for every complex $s$, $\xi(\frac12+s)=e^{i\eta s}\int_{\mathbb{R}}\Phi(x+i\eta)\,e^{sx}\,dx$ and $\xi'(\frac12+s)=e^{i\eta s}\int_{\mathbb{R}}(x+i\eta)\,\Phi(x+i\eta)\,e^{sx}\,dx$. For real $x$, $|\Phi(x+i\eta)|\le C_\delta$ and $|(x+i\eta)\Phi(x+i\eta)|\le C_\delta$, with $C_\delta=O(\delta^{-5/2})$ as $\delta\to 0$.

Source: [`ch/ch04.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch04.tex), label `prop:tilt`.

## Chapter 3. The Thorin measure at finite resolution

<a id="thorin_mass_first_order_defect"></a>

### Mass defect of the Thorin measure on a window

Proposition 3.2.2, page 58 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed zeta function, $\Xi(\theta)=\xi(\frac12+i\theta)$, and for $\varepsilon>0$ let $\mu_\varepsilon$ be the signed measure on $(0,\infty)$ whose density in the variable $\theta=\sqrt{z}$ is the Thorin density $\rho_\varepsilon(\theta)=\pi^{-1}\operatorname{Re}(\xi'/\xi)(\frac12+\varepsilon+i\theta)$, so that $\mu_\varepsilon((0,\Theta^2])=\int_0^\Theta\rho_\varepsilon(\theta)\,d\theta$. Let $N(\Theta)$ be the number of zeros $\rho$ of $\xi$ with $0<\operatorname{Im}\rho\le\Theta$, counted with multiplicity. Suppose that $\Theta$ is not the ordinate of a zero and that every zero with $0<\operatorname{Im}\rho\le\Theta$ lies on the line $\operatorname{Re}s=\frac12$. Then $\mu_\varepsilon((0,\Theta^2])=N(\Theta) - \frac{\varepsilon}{\pi}\,\frac{\Xi'(\Theta)}{\Xi(\Theta)}+O(\varepsilon^2)$ as $\varepsilon\downarrow 0$. The book also shows that the remainder is in fact $O(\varepsilon^3)$.

Source: [`ch/ch07.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch07.tex), label `prop:mass`.

<a id="off_line_zero_signatures"></a>

### Signatures of an off line zero in the Thorin density

Proposition 3.3.1, page 61 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $P_a(x)=a/(\pi(a^2+x^2))$ for $a>0$ be the Poisson (Cauchy) kernel, $P_{-a}=-P_a$, and $P_0$ the unit point mass at the origin, and write $\|\cdot\|_{\mathrm{TV}}$ for the total variation norm on $\mathbb{R}$. For $\varepsilon>0$ the Thorin density $\rho_\varepsilon(\theta)=\pi^{-1}\operatorname{Re}(\xi'/\xi)(\frac12+\varepsilon+i\theta)$ of the completed zeta function $\xi$ is a sum over the zeros $\rho=\beta+i\gamma$ with $\gamma>0$ of kernels $P_{1/2+\varepsilon - \beta}(\theta - \gamma)$ together with their reflections at $-\gamma$: a zero on the line $\operatorname{Re}s=\frac12$ with ordinate $\gamma$ contributes $P_\varepsilon$ centred at $\gamma$, and a pair of off line zeros $\frac12\pm\delta+i\gamma$ with $0<\delta<\frac12$ contributes $P_{\varepsilon+\delta}+P_{\varepsilon - \delta}$ centred at $\gamma$. Let $\varepsilon>0$ and $0<\delta<\frac12$. (i) If two line zeros at height $\gamma$ are replaced by the quadruple $\frac12\pm\delta\pm i\gamma$, the density near $\gamma$ changes by $2g$, where $g=\frac12(P_{\varepsilon+\delta}+P_{\varepsilon - \delta}) - P_\varepsilon$ is the change per displaced zero, and for $\delta<\varepsilon$, $\|g\|_{\mathrm{TV}}=\frac{3\sqrt{3}}{4\pi}(\delta/\varepsilon)^2+O(\delta^4/\varepsilon^4)$. (ii) If the pair $\frac12\pm\delta+i\gamma$ exists but a list of line zeros omits it, the difference is $P_{\varepsilon+\delta}+P_{\varepsilon - \delta}$, and its total variation norm is $2$ for $\delta\le\varepsilon$ and $\frac8\pi\arctan\sqrt{(\delta+\varepsilon)/(\delta - \varepsilon)} - 2=\frac4\pi\arcsin(\varepsilon/\delta)$ for $\delta>\varepsilon$. (iii) For $\delta>\varepsilon$ the negative part of $\frac12(P_{\varepsilon+\delta}+P_{\varepsilon - \delta})$, the density of one displaced zero, has mass $\frac2\pi\arctan\sqrt{(\delta+\varepsilon)/(\delta - \varepsilon)} - \frac12=\frac1\pi\arcsin(\varepsilon/\delta)=\frac{\varepsilon}{\pi\delta}+O(\varepsilon^3/\delta^3)$. The expansion $g=\frac{\delta^2}{2}\,\partial_a^2P_a|_{a=\varepsilon}+O(\delta^4)$ behind (i), by which the signature grows like the square of the displacement, is Proposition 19.7 of Polson (2026).

Source: [`ch/ch08.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch08.tex), label `prop:signature`.

<a id="littlewood_identity_thorin_mass"></a>

### Littlewood identity for the cumulative Thorin mass

Proposition 3.4.1, page 66 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed zeta function and, for $\varepsilon>0$ and $T>0$, let $M_\varepsilon(T)=\pi^{-1}\int_0^T\operatorname{Re}(\xi'/\xi)(\frac12+\varepsilon+it)\,dt$ be the cumulative Thorin mass at resolution $\varepsilon$ up to height $T$; it equals $\pi^{-1}\arg\xi(\frac12+\varepsilon+iT)$, with the argument continued along the vertical line from the real axis, where $\xi>0$. Let $T>0$. For every $\varepsilon>0$ such that $\xi$ has no zero on the segment from $\frac12+\varepsilon$ to $\frac12+\varepsilon+iT$, $\pi M_\varepsilon(T)=\frac{\partial}{\partial\varepsilon}\int_0^T\log\bigl|\xi(\frac12+\varepsilon+it)\bigr|\,dt$, and for every $\varepsilon>0$, $\int_0^T\log\bigl|\xi(\frac12+\varepsilon+it)\bigr|\,dt=\int_0^T\log\bigl|\xi(\frac12+it)\bigr|\,dt+\pi\int_0^\varepsilon M_{\varepsilon'}(T)\,d\varepsilon'$, where $M_{\varepsilon'}(T)$ is defined except for the finitely many $\varepsilon'$ at which the segment meets a zero.

Source: [`ch/ch09.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch09.tex), label `prop:littlewood`.

<a id="zero_density_identity_thorin_mass"></a>

### Zero density identity from the Thorin mass

Proposition 3.4.2, page 67 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed zeta function. For $\varepsilon>0$ and $T>0$ let $M_\varepsilon(T)=\pi^{-1}\arg\xi(\frac12+\varepsilon+iT)$ be the cumulative Thorin mass, the argument being continued along the vertical line from the real axis, where $\xi>0$, and let $\tau_\varepsilon(T)$ be $\pi^{-1}$ times the change of $\arg\xi$ along the horizontal segment from $\frac12+\varepsilon+iT$ to $\frac12+iT$. Let $N(T)$ be the number of zeros $\rho$ with $0<\operatorname{Im}\rho\le T$, and $\nu(\sigma,T)$ the number of those with $\operatorname{Re}\rho>\sigma$, both counted with multiplicity. Let $\varepsilon>0$ and $T>0$, with $T$ not the ordinate of a zero and no zero of $\xi$ on the segment from $\frac12+\varepsilon$ to $\frac12+\varepsilon+iT$. Then $M_\varepsilon(T)+\tau_\varepsilon(T)=N(T) - 2\,\nu(\frac12+\varepsilon,T)$. The segment condition holds at $\varepsilon=\frac12$, and $\nu(1,T)=0$, which gives the unconditional count $N(T)=M_{1/2}(T)+\tau_{1/2}(T)$, and $\nu(\frac12+\varepsilon,T)=\frac12\,[N(T) - M_\varepsilon(T) - \tau_\varepsilon(T)]$ for every $\varepsilon$ satisfying the hypotheses.

Source: [`ch/ch09.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch09.tex), label `prop:zerodensity`.

## Chapter 4. Weak limits, hyperbolic complete monotonicity and a Bosch type theorem

<a id="rh_iff_clock_is_weak_limit_of_ggc_transforms"></a>

### RH as a weak limit criterion for generalized gamma convolutions

Corollary 4.1.2, page 77 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi(s)=\tfrac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function and $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, so that $1/G$ is the reciprocal Xi function (the clock). A generalized gamma convolution (GGC) is a law on $[0,\infty)$ whose Laplace transform is $\exp\bigl(-bx-\int_{(0,\infty)}\log(1+x/t)\,U(dt)\bigr)$ with drift $b\ge0$ and Thorin measure $U\ge0$. The Riemann Hypothesis holds if and only if there are GGC random variables $X_N$ with $\mathbb{E}e^{-uX_N}\to1/G(u)$ for all $u$ in a set with an accumulation point in $(0,\infty)$.

Source: [`ch/ch11.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch11.tex), label `cor:weak:closure`.

<a id="reciprocal_xi_hcm_iff_rh"></a>

### The reciprocal Xi function is HCM if and only if RH holds

Corollary 4.1.3, page 80 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\tfrac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function. A positive continuous function $f$ on $(0,\infty)$ is hyperbolically completely monotone (HCM) if for every $u>0$ the function $v\mapsto f(uv)f(u/v)$ is completely monotone in $w=v+v^{-1}$ on $[2,\infty)$. (i) For every $\alpha\ge1$ the function $\varphi_\alpha(u)=\xi(\alpha)/\xi(\alpha+\sqrt u)$ is HCM, unconditionally. (ii) The reciprocal Xi function $\varphi_{1/2}(u)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt u)$, the clock, is HCM if and only if the Riemann Hypothesis holds. Part (ii) is noted in versions 21 and 23 of Polson's arXiv paper 1708.02653; here (i) and (ii) follow from Theorem 3.7 and Proposition 3.4 of Polson (2026) with Bondesson's theorem.

Source: [`ch/ch11.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch11.tex), label `cor:hcm`.

<a id="deformed_reciprocal_xi_hcm_iff_time_at_least_lambda"></a>

### Deformed reciprocal Xi function is HCM exactly for t at least Lambda

Corollary 4.2.5, page 88 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, and $\Phi(u)=2\sum_{n\ge1}\bigl(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\bigr)\exp(-\pi n^2e^{2u})$ the Pólya kernel, so that $\Xi(\theta)=\int_{\mathbb{R}}\Phi(u)\cos(\theta u)\,du$. For real $t$ let $H_t(z)=\tfrac12\int_0^\infty e^{tu^2}\Phi(2u)\cos(zu)\,du$, so that $H_0(z)=\tfrac18\Xi(z/2)$, and let $\varphi_t(s)=H_t(0)/H_t(2i\sqrt s)$ for $s>0$, so that $\varphi_0(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$. Let $\Lambda$ be the de Bruijn Newman constant, defined by the property that $H_t$ has only real zeros exactly when $t\ge\Lambda$. Then $\varphi_t$ is HCM if and only if $t\ge\Lambda$. Hence $\varphi_t$ is HCM for $t\ge0.22$ (by the bound $\Lambda\le0.22$ of the Polymath project), $\varphi_t$ is not HCM for $t<0$ (by the theorem of Rodgers and Tao that $\Lambda\ge0$), and RH holds if and only if $\varphi_t$ is HCM for every $t>0$.

Source: [`ch/boschA.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschA.tex), label `bh:cor:dbn`.

<a id="question_hcm_powers_of_mellin_variable_w"></a>

### Question: for which beta is the density of W^beta HCM

Question 4.2.7, page 89 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function and $\Xi(z)=\xi(\tfrac12+iz)$. Let $W$ be the positive random variable, which exists unconditionally, with $\mathbb{E}W^{it}=\xi(\tfrac12)/\xi(\tfrac12+t)$ for real $t$ and $\mathbb{E}W^z=\Xi(0)/\Xi(z)$ for $|\operatorname{Re}z|<\gamma_1$, where $\gamma_1\approx14.1347$ is the smallest positive ordinate of a zero of $\xi$. For which real $\beta$ is the density of $W^\beta$ HCM (hyperbolically completely monotone)? On the Cauchy side, for $|C|$ with $C$ a standard Cauchy variable, so that $\mathbb{E}|C|^z=\sec(\pi z/2)$, the answer is $|\beta|\ge2$ (Bosch's theorem).

Source: [`ch/boschA.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschA.tex), label `bh:q:beta`.

<a id="question_single_power_hcm_implies_rh"></a>

### Question: does HCM of W^beta for a single beta imply RH

Question 4.2.8, page 89 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function and $\Xi(z)=\xi(\tfrac12+iz)$, and let $W$ be the positive random variable, which exists unconditionally, with $\mathbb{E}W^z=\Xi(0)/\Xi(z)$ for $|\operatorname{Re}z|<\gamma_1$, where $\gamma_1\approx14.1347$ is the smallest positive ordinate of a zero of $\xi$. Does HCM (hyperbolic complete monotonicity) of the density of $W^\beta$ for a single real $\beta$ imply the Riemann Hypothesis?

Source: [`ch/boschA.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschA.tex), label `bh:q:single`.

<a id="elasticity_criterion_for_hcm_functions"></a>

### Elasticity criterion for HCM functions (Bondesson's representation restated as a Pick condition)

Theorem 4.3.2, page 91 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Omega=\mathbb{C}\setminus(-\infty,0]$ and let $h$ be a positive function on $(0,\infty)$. Call $h$ HCM if it is continuous and for every $u>0$ the function $v\mapsto h(uv)h(u/v)$ is completely monotone in $w=v+v^{-1}$ on $[2,\infty)$. Then $h$ is HCM if and only if $h$ continues to a zero free analytic function on $\Omega$ and its elasticity $\varepsilon_h(x)=-xh'(x)/h(x)$ is a Pick function, that is, $\operatorname{Im}\varepsilon_h(x)\ge0$ for $\operatorname{Im}x>0$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:crit`.

<a id="reciprocal_xi_power_hcm_threshold"></a>

### Bosch type theorem for the reciprocal Xi function

Theorem 4.3.6, page 94 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, and let $\gamma$ run over the zeros of $\Xi$ with $\operatorname{Re}\gamma>0$, counted with multiplicity. Let $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)=\prod_\gamma(1+s/\gamma^2)^{-1}$ be the reciprocal Xi function, let $\delta_*=\max_\gamma|\arg\gamma|$ be the largest angle of these zeros (the maximum is attained if some zero is off the critical line, and $\delta_*=0$ if and only if the Riemann Hypothesis holds), and put $\alpha_\xi=2\delta_*/\pi$ and $p_*=1-\alpha_\xi$. Unconditionally, for real $p$ the function $x\mapsto\varphi(x^p)$ on $(0,\infty)$ is HCM (hyperbolically completely monotone) if and only if $|p|\le p_*=1-\alpha_\xi$. Equivalently, $\varphi(x^{1/\beta})$ is HCM if and only if $|\beta|\ge1/(1-\alpha_\xi)$. The Riemann Hypothesis is the statement $\alpha_\xi=0$, that is, $p_*=1$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:power`.

<a id="rh_iff_hcm_of_all_powers_below_one_and_pick_elasticity"></a>

### Equivalent forms of RH: HCM of all powers, Pick elasticity

Corollary 4.3.7, page 94 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ the reciprocal Xi function and $\varepsilon_\varphi(s)=-s\varphi'(s)/\varphi(s)=\tfrac{\sqrt s}{2}\,\frac{\xi'}{\xi}\bigl(\tfrac12+\sqrt s\bigr)$ its elasticity. The following are equivalent. (i) The Riemann Hypothesis. (ii) $\varphi$ is HCM. (iii) $\varphi(x^p)$ is HCM for every $|p|<1$. (iv) $\varepsilon_\varphi$ is a Pick function (analytic on the upper half plane with nonnegative imaginary part). (v) $\operatorname{Im}\bigl[z\,(\xi'/\xi)(\tfrac12+z)\bigr]\ge0$ for all $z$ with $\operatorname{Re}z>0$ and $\operatorname{Im}z>0$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:cor:rh`.

<a id="reciprocal_xi_power_hcm_unconditional_range"></a>

### Unconditional HCM range for the reciprocal Xi function composed with a power

Corollary 4.3.8, page 95 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma$ the zeros of $\Xi$ with $\operatorname{Re}\gamma>0$ counted with multiplicity, $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ the reciprocal Xi function with elasticity $\varepsilon_\varphi(s)=-s\varphi'(s)/\varphi(s)$, and $\delta_*=\max_\gamma|\arg\gamma|$. Let $T_0=3\times10^{12}$. RH has been verified for all zeros with $0<\operatorname{Re}\gamma\le T_0$ (Platt and Trudgian 2021). Hence $\delta_*<\arctan(1/(2T_0))$ and $\varphi(x^p)$ is HCM for all $|p|\le1-\frac{1}{\pi T_0}$, where $1-\frac{1}{\pi T_0}>1-1.1\times10^{-13}$. Unconditionally, $\operatorname{Im}\varepsilon_\varphi(s)\ge0$ on the sector $0<\arg s<\pi-1/T_0$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:cor:uncond`.

<a id="stieltjes_sign_of_clock_outside_parabolic_region"></a>

### Stieltjes sign of minus the logarithmic derivative of the clock outside a thin parabolic region

Theorem 4.3.10, page 96 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma$ the zeros of $\Xi$ with $\operatorname{Re}\gamma>0$ counted with multiplicity, $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ the reciprocal Xi function and $\psi(s)=-\varphi'(s)/\varphi(s)=\sum_\gamma(s+\gamma^2)^{-1}$. A Stieltjes function is a function $a/s+b+\int_{(0,\infty)}\sigma(dt)/(s+t)$ with $a,b\ge0$ and $\sigma\ge0$ a measure satisfying $\int(1+t)^{-1}\sigma(dt)<\infty$. Let $T_0=3\times10^{12}$, so that every zero $\gamma$ with $0<\operatorname{Re}\gamma\le T_0$ is real (Platt and Trudgian 2021), and let $\mathcal E_{T_0}=\bigl\{x+iy:\ x<-T_0^2+T_0+\tfrac14,\ y^2<\tfrac12-x\bigr\}$. Split $\psi=S_{\mathrm{on}}+R_{\mathrm{off}}$, where $S_{\mathrm{on}}(s)=\sum_{\gamma\ \mathrm{real}}(s+\gamma^2)^{-1}$ is the sum over the zeros on the critical line and $R_{\mathrm{off}}$ is the same sum over the nonreal $\gamma$. For a nonreal zero $\gamma$ let $D_\gamma$ be the closed disc with diameter the segment from $-\gamma^2$ to $-\bar\gamma^2$. Then (i) $S_{\mathrm{on}}$ is a Stieltjes function, with measure $\sum_{\gamma\ \mathrm{real}}\delta_{\gamma^2}$; (ii) $R_{\mathrm{off}}$ is meromorphic on $\mathbb{C}$ and nonnegative on $[0,\infty)$, for every nonreal $\gamma$ the disc $D_\gamma$, and with it the poles $-\gamma^2$ and $-\bar\gamma^2$ of $R_{\mathrm{off}}$ on its boundary, lies in $\mathcal E_{T_0}$, and $\operatorname{Im}R_{\mathrm{off}}\le0$ at every point of the upper half plane outside these discs; (iii) $\psi>0$ on $[0,\infty)$, and $\operatorname{Im}\psi(s)<0$ at every $s=x+iy$ with $y>0$ outside $\mathcal E_{T_0}$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:local`.

<a id="question_hcm_exponent_deficit_negative_heat_time"></a>

### Question: HCM exponent of the deformed reciprocal Xi function for negative time

Question 4.3.12, page 99 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\Phi$ the Pólya kernel (with $\Xi(\theta)=\int_{\mathbb{R}}\Phi(u)\cos(\theta u)\,du$) and, for real $t$, $H_t(z)=\tfrac12\int_0^\infty e^{tu^2}\Phi(2u)\cos(zu)\,du$ the de Bruijn Newman deformation, so that $H_0(z)=\tfrac18\Xi(z/2)$. Let $\delta_*(t)$ be the largest angle $|\arg z|$ over the zeros $z$ of $H_t$ with $\operatorname{Re}z>0$ and $p_*(t)=1-\tfrac2\pi\delta_*(t)$; whenever the zeros of $H_t$ satisfy $|\arg z|<\pi/4$, the function $x\mapsto H_t(0)/H_t(2i\sqrt{x^p})$ is HCM if and only if $|p|\le p_*(t)$. Describe $p_*(t)$ for $t<0$. By Rodgers and Tao, $p_*(t)<1$ for every $t<0$. Assume RH; how fast does $1-p_*(t)\to0$ as $t\uparrow0$?

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:four`.

<a id="thorin_density_of_reciprocal_xi_power"></a>

### Thorin measure of the reciprocal Xi function composed with a power

Theorem 4.3.13, page 99 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma$ the zeros of $\Xi$ with $\operatorname{Re}\gamma>0$ counted with multiplicity, $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ the reciprocal Xi function with elasticity $\varepsilon_\varphi(s)=-s\varphi'(s)/\varphi(s)=\tfrac{\sqrt s}{2}\frac{\xi'}{\xi}(\tfrac12+\sqrt s)=\sum_\gamma\frac{s}{s+\gamma^2}$, $\delta_*=\max_\gamma|\arg\gamma|$ and $p_*=1-2\delta_*/\pi$. A generalized gamma convolution (GGC) is a law on $[0,\infty)$ whose Laplace transform is $\exp\bigl(-bs-\int_{(0,\infty)}\log(1+s/t)\,U(dt)\bigr)$ with drift $b\ge0$ and Thorin measure $U\ge0$. Let $0<p\le p_*$. Then $\varphi(s^p)=\mathbb{E}e^{-sX_p}$ for a GGC variable $X_p$ with no drift. If $p<p_*$, its Thorin measure has the density $u_p(\tau)=\frac{p}{\pi\tau}\operatorname{Im}\varepsilon_\varphi\bigl(\tau^pe^{ip\pi}\bigr)=\frac{p}{2\pi\tau}\operatorname{Im}\Bigl[z\,\frac{\xi'}{\xi}\bigl(\tfrac12+z\bigr)\Bigr]_{z=\tau^{p/2}e^{ip\pi/2}}\ge0$ for $\tau>0$, which in terms of the zeros is $u_p(\tau)=\frac{p}{\pi\tau}\sum_\gamma\operatorname{Im}\frac{\tau^pe^{ip\pi}}{\tau^pe^{ip\pi}+\gamma^2}$. If $p=p_*$, its Thorin measure is $u_{p_*}(\tau)\,d\tau$ plus a unit atom at $|\gamma|^{2/p_*}$ for every zero $\gamma$ with $\arg\gamma=-\delta_*$, counted with multiplicity, that is, one atom for each conjugate pair of zeros at the extreme angle; when $p_*=1$ there is an atom at $\gamma^2$ for every zero, and $u_1=0$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:ggc`.

<a id="thorin_measure_of_reciprocal_xi_power_under_rh_lamperti_kernel"></a>

### Thorin measure of X_p under RH: scaled Lamperti laws (Bosch's variable)

Corollary 4.3.14, page 100 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on the Riemann Hypothesis.

Statement. Assume the Riemann Hypothesis. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma$ the (real, positive) zeros of $\Xi$ counted with multiplicity, $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$, and for $0<p<1$ let $X_p$ be the generalized gamma convolution (GGC) variable with $\mathbb{E}e^{-sX_p}=\varphi(s^p)$; a GGC is a law on $[0,\infty)$ whose Laplace transform is $\exp\bigl(-bs-\int_{(0,\infty)}\log(1+s/t)\,U(dt)\bigr)$ with drift $b\ge0$ and Thorin measure $U\ge0$. Let $L_p=Z_p/\tilde Z_p$ be the ratio of independent positive $p$ stable variables, whose density $\frac{\sin p\pi}{\pi}\,\frac{t^{p-1}}{t^{2p}+2t^p\cos p\pi+1}$ on $(0,\infty)$ is Lamperti's; its power $L_p^{\,p}$ is Bosch's variable $T_\alpha$ with $\alpha=p$. Writing $\mathcal L(Y)$ for the law of a random variable $Y$, the Thorin measure of $X_p$ is $U_p=p\sum_\gamma\mathcal L\bigl(\gamma^{2/p}L_p\bigr)$, and $X_p\stackrel{d}{=}S_pX^{1/p}$, where $X=\sum_\gamma e_\gamma/\gamma^2$ is the clock variable, with independent standard exponential variables $e_\gamma$, and $S_p$ is a positive $p$ stable variable independent of $X$ with $\mathbb{E}e^{-\lambda S_p}=e^{-\lambda^p}$. At $p=\tfrac12$, $L_{1/2}\stackrel{d}{=}C^2$ for a standard Cauchy variable $C$, and $U_{1/2}=\tfrac12\sum_\gamma\mathcal L(\gamma^4C^2)$. As $p\uparrow1$, $L_p\to1$ and $U_p\to\sum_\gamma\delta_{\gamma^2}$ vaguely, the Thorin measure of $\varphi$ itself.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:cor:lamperti`.

<a id="rh_iff_reciprocal_xi_power_is_ggc_transform_for_all_p"></a>

### RH is equivalent to the existence of the GGC variable X_p for every p in (0, 1)

Corollary 4.3.15, page 101 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $\varphi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ the reciprocal Xi function. A generalized gamma convolution (GGC) is a law on $[0,\infty)$ whose Laplace transform is $\exp\bigl(-bs-\int_{(0,\infty)}\log(1+s/t)\,U(dt)\bigr)$ with drift $b\ge0$ and Thorin measure $U\ge0$. For $0<p<1$ say that $X_p$ exists if $\varphi(s^p)=\mathbb{E}e^{-sX_p}$ for a GGC random variable $X_p$. Then the Riemann Hypothesis holds if and only if $X_p$ exists for every $p\in(0,1)$. Unconditionally, using the verification of the Riemann Hypothesis to height $T_0=3\times10^{12}$ (Platt and Trudgian 2021), $X_p$ exists for $0<p\le1-1/(\pi T_0)$, in particular for $0<p\le1-1.1\times10^{-13}$, and its Thorin measure then has the density $u_p(\tau)=\frac{p}{2\pi\tau}\operatorname{Im}\Bigl[z\,\frac{\xi'}{\xi}\bigl(\tfrac12+z\bigr)\Bigr]_{z=\tau^{p/2}e^{ip\pi/2}}$ for $\tau>0$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:cor:rhX`.

<a id="question_construct_x_half_without_zeros"></a>

### Question: construct X_{1/2} from Brownian motion without the zeros

Question 4.3.16, page 101 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function and, for $0<p<1$, let $X_p$ be the variable with $\mathbb{E}e^{-sX_p}=\xi(\tfrac12)/\xi(\tfrac12+s^{p/2})$, a generalized gamma convolution whenever it exists; unconditionally it exists for $0<p\le1-1.1\times10^{-13}$, and RH holds if and only if it exists for every $p\in(0,1)$. In particular $X_{1/2}$ has Laplace transform $\xi(\tfrac12)/\xi(\tfrac12+s^{1/4})$ and exists without any hypothesis. Construct $X_{1/2}$ from Brownian motion, Bessel processes or the Biane, Pitman and Yor variables without reference to the zeros. A construction that survived $p\uparrow1$ would prove RH.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:five`.

<a id="question_thorin_density_positive_for_all_powers"></a>

### Question: prove positivity of u_p for all p without locating zeros

Question 4.3.17, page 102 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function. For $0<p<1$ and $\tau>0$ put $u_p(\tau)=\frac{p}{2\pi\tau}\operatorname{Im}\Bigl[z\,\frac{\xi'}{\xi}\bigl(\tfrac12+z\bigr)\Bigr]_{z=\tau^{p/2}e^{ip\pi/2}}$, the Thorin density of the reciprocal Xi function composed with $s^p$; the formula makes sense whatever the zeros are. Prove $u_p\ge0$ for every $p\in(0,1)$ directly from this formula, without locating the zeros. (Positivity of $u_p$ for all $p\in(0,1)$ is equivalent to the Riemann Hypothesis.)

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:one`.

<a id="thorin_density_of_power_positive_beyond_strip_from_primes"></a>

### Thorin density u_p is positive beyond a finite point, from the primes

Proposition 4.3.19, page 103 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $0<p<1$, $\vartheta=p\pi/2$, and for $\tau>0$ put $t=\tau^{p/2}$, $z=te^{i\vartheta}$, $w=\tfrac12+z$ and $u_p(\tau)=\frac{p}{2\pi\tau}\operatorname{Im}\bigl[z\,\frac{\xi'}{\xi}(w)\bigr]$, the Thorin density of the reciprocal Xi function composed with $s^p$ (the formula makes sense whatever the zeros are). Let $\operatorname{Re}w>1$ and $\sigma=\operatorname{Re}w$. Then $\frac{2\pi\tau}{p}\,u_p(\tau)\ \ge\ m_p(t):=\operatorname{Im}\Bigl[z\Bigl(\frac{1}{2w}+\frac{1}{w-1}+\frac12\log\frac{w}{2\pi}\Bigr)\Bigr]-\frac{t}{6\sigma^2}-t\Bigl(-\frac{\zeta'}{\zeta}\Bigr)(\sigma)$. Every term of $m_p$ is explicit, no zero of $\zeta$ enters, and $m_p(t)=\frac t2\sin\vartheta\,\log t\,(1+o(1))$ as $t\to\infty$. Consequently for each $p\in(0,1)$ there is $\tau_1(p)<\infty$ with $u_p(\tau)>0$ for all $\tau>\tau_1(p)$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:prop:prime`.

<a id="cauchy_clock_density_eta_cubed_natural_boundary"></a>

### Cauchy analogue of the clock: eta cubed density with natural boundary, not HCM

Theorem 4.3.21, page 104 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X_C=\sum_{k\ge1}e_k/(2k-1)^2$ with independent standard exponential variables $e_k$, so that $\mathbb{E}e^{-sX_C}=\operatorname{sech}(\pi\sqrt s/2)$. The density of $X_C$ is $f_C(t)=\frac4\pi\sum_{k\ge1}(-1)^{k-1}(2k-1)e^{-(2k-1)^2t}=\frac4\pi\,\eta\Bigl(\frac{4it}{\pi}\Bigr)^3$ for $t>0$, with $\eta$ the Dedekind eta function. The line $\operatorname{Re}t=0$ is a natural boundary of $f_C$. Hence $X_C$ is a generalized gamma convolution whose Laplace transform is HCM, but the density of $X_C$ is not HCM.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:eta`.

<a id="question_natural_boundary_of_clock_density"></a>

### Question: natural boundary of the clock density

Question 4.3.22, page 105 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Assume the Riemann Hypothesis and that all zeros are simple. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma_1<\gamma_2<\cdots$ the positive zeros of $\Xi$, $X=\sum_ne_n/\gamma_n^2$ the clock variable with independent standard exponential variables $e_n$ (its Laplace transform is $\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$), and $A_n=\prod_{m\ne n}\bigl(1-\gamma_n^2/\gamma_m^2\bigr)^{-1}$. Does the density $f_X(t)=\sum_nA_n\gamma_n^2e^{-\gamma_n^2t}$ of $X$ have $\operatorname{Re}t=0$ as a natural boundary? If so, $X$ is not HCM.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:two`.

<a id="clock_density_natural_boundary_under_gap_hypotheses"></a>

### Natural boundary of the clock density under RH, simple zeros and gap hypotheses

Theorem 4.3.23, page 105 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on the Riemann Hypothesis; all zeros of the zeta function are simple; hypothesis (G), $\inf_n(\gamma_{n+1}^2-\gamma_n^2)>0$; hypothesis (D), $\log(1/|\zeta'(\tfrac12+i\gamma_n)|)=o(\gamma_n^2)$.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, and assume the Riemann Hypothesis and that all zeros are simple; let $\gamma_1<\gamma_2<\cdots$ be the positive zeros of $\Xi$ (so $\tfrac12+i\gamma_n$ are the zeros of $\xi$ in the upper half plane), and let $X=\sum_ne_n/\gamma_n^2$ be the clock variable, with independent standard exponential variables $e_n$ and Laplace transform $\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$. Assume (G) $\inf_n(\gamma_{n+1}^2-\gamma_n^2)>0$ and (D) $\log\bigl(1/|\zeta'(\tfrac12+i\gamma_n)|\bigr)=o(\gamma_n^2)$. Then for $t>0$ the density of $X$ is $f_X(t)=\sum_nc_ne^{-\gamma_n^2t}$ with $c_n=A_n\gamma_n^2=-\frac{2\gamma_n\Xi(0)}{\Xi'(\gamma_n)}$, where $A_n=\prod_{m\ne n}(1-\gamma_n^2/\gamma_m^2)^{-1}$, the series has abscissa of convergence $0$, and $\operatorname{Re}t=0$ is a natural boundary of $f_X$. In particular the density of $X$ is not HCM.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:nb`.

<a id="question_which_powers_of_clock_variable_are_hcm"></a>

### Question: for which p are X_p and its Cauchy analogue HCM

Question 4.3.25, page 107 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function and, for $0<p<1$, let $X_p$ be the variable with $\mathbb{E}e^{-sX_p}=\xi(\tfrac12)/\xi(\tfrac12+s^{p/2})$, a generalized gamma convolution whenever it exists (unconditionally for $0<p\le1-1.1\times10^{-13}$, and for every $p<1$ under RH), and let $X_{C,p}$ be the variable with Laplace transform $\operatorname{sech}(\pi s^{p/2}/2)$, the Cauchy family. Determine $\{p:\ X_p\ \text{is HCM}\}$, and the same set for the Cauchy family $X_{C,p}$. (A random variable is HCM if it has an HCM density.)

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:three`.

<a id="mellin_variable_w_density_lower_bound_unconditional"></a>

### Lower bound for the density of log W: existence of W (second proof)

Theorem 4.3.26, page 108 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\gamma$ the zeros of $\Xi$ with $\operatorname{Re}\gamma>0$ counted with multiplicity, $\chi(t)=\xi(\tfrac12)/\xi(\tfrac12+t)=\prod_\gamma(1+t^2/\gamma^2)^{-1}$ for real $t$, and $g(y)=\frac1\pi\int_0^\infty\chi(t)\cos(ty)\,dt$. Unconditionally, by the verification of RH to height $T_0=3\times10^{12}$ (Platt and Trudgian 2021), $g(y)\ge0.9998\,g_R(y)>0$ for all real $y$, where $g_R$ is the density of $\sum_{\gamma\ \mathrm{real}}L_\gamma/\gamma$, the sum over the real zeros $\gamma$ with independent standard Laplace variables $L_\gamma$ (density $\tfrac12e^{-|x|}$). Hence $\chi$ is a characteristic function and there is a positive random variable $W$ with $\mathbb{E}W^{it}=\chi(t)$ and $\mathbb{E}W^z=\Xi(0)/\Xi(z)$ for $|\operatorname{Re}z|<\gamma_1$, where $\gamma_1\approx14.1347$ is the smallest positive ordinate of a zero of $\xi$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:gpos`.

<a id="elasticity_criterion_for_powers_of_mellin_variable_w"></a>

### HCM criterion for powers of the Mellin variable W

Proposition 4.3.27, page 109 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\chi(t)=\xi(\tfrac12)/\xi(\tfrac12+t)$ and $g(y)=\frac1\pi\int_0^\infty\chi(t)\cos(ty)\,dt$, the density of $\log W$ for the positive random variable $W$ with $\mathbb{E}W^{it}=\chi(t)$; $g$ is an entire function because $\chi(t)=\exp\bigl(-\tfrac12|t|\log|t|\,(1+o(1))\bigr)$. Let $\beta\in\mathbb{R}\setminus\{0\}$. Then $W^\beta$ is HCM (has an HCM density) if and only if $g$ has no zeros in the strip $|\operatorname{Im}u|<\pi/|\beta|$ and $\operatorname{Im}\bigl(-g'(u)/g(u)\bigr)\ge0$ for $0<\operatorname{Im}u<\pi/|\beta|$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:prop:Wcrit`.

<a id="first_gap_necessary_condition_for_hcm_of_w_power"></a>

### First gap theorem: HCM of W^beta forces |beta| at least gamma_2 minus gamma_1

Theorem 4.3.28, page 109 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on the Riemann Hypothesis; all zeros of the zeta function are simple; $\log(1/|\zeta'(\rho_n)|)=O(\gamma_n)$, where $\rho_n=\tfrac12+i\gamma_n$ are the zeros of zeta in the upper half plane.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\chi(t)=\xi(\tfrac12)/\xi(\tfrac12+t)$, and $g(u)=\frac1\pi\int_0^\infty\chi(t)\cos(tu)\,dt$ (an entire function), the density of $\log W$ for the positive random variable $W$ with $\mathbb{E}W^{it}=\chi(t)$. Assume the Riemann Hypothesis, that all zeros are simple, and that $\log(1/|\zeta'(\rho_n)|)=O(\gamma_n)$, where $\rho_n=\tfrac12+i\gamma_n$ with $0<\gamma_1<\gamma_2<\cdots$ are the zeros of $\zeta$ in the upper half plane. Then $u_0=\frac\pi4+\limsup_n\gamma_n^{-1}\log(1/|\zeta'(\rho_n)|)$ is finite, and for $\operatorname{Re}u>u_0$ one has $g(u)=-\Xi(0)\sum_n\frac{e^{-\gamma_nu}}{\Xi'(\gamma_n)}=\sum_nb_ne^{-\gamma_nu}$ with $b_1>0>b_2$. If $W^\beta$ is HCM (has an HCM density), then $|\beta|\ge\gamma_2-\gamma_1=6.8873\ldots$

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:thm:gap`.

<a id="question_first_gap_conjecture"></a>

### Question: first gap conjecture for HCM powers of W

Question 4.3.29, page 110 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $W$ be the positive random variable, which exists unconditionally, with $\mathbb{E}W^z=\Xi(0)/\Xi(z)$ for $|\operatorname{Re}z|<\gamma_1$, where $\Xi(z)=\xi(\tfrac12+iz)$, $\xi$ is the completed Riemann zeta function and $0<\gamma_1<\gamma_2$ are the two smallest positive ordinates of zeros of $\xi$. Is $W^\beta$ HCM (does it have an HCM density) exactly when $|\beta|\ge\gamma_2-\gamma_1$? On the Cauchy side the analogous statement, with frequencies $2k-1$, is Bosch's theorem at $\alpha=1/2$: $3-1=2$.

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:gap`.

<a id="question_first_gap_conjecture_refined"></a>

### Refined first gap conjecture for sums of Laplace variables

Question 4.3.30, page 112 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $Y=\sum_nL_n/\lambda_n$ (a convergent series) with independent standard Laplace variables $L_n$ and $0<\lambda_1<\lambda_2<\cdots$, and let $W=e^Y$. If the characteristic function of $Y$ decays faster than every exponential, or if the first gap $\lambda_2-\lambda_1$ is at least every later gap, is the HCM threshold of $W^\beta$ exactly $\lambda_2-\lambda_1$? In particular, is it $\gamma_2-\gamma_1$ for the Riemann variable $W$ with $\mathbb{E}W^z=\Xi(0)/\Xi(z)$?

Source: [`ch/boschB.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/boschB.tex), label `bs:q:gaprefined`.

## Chapter 5. The heat trace, the arithmetic and the limits of finite information

<a id="clock_infinitely_divisible_and_self_decomposable"></a>

### The clock is infinitely divisible and self decomposable, unconditionally

Theorem 5.1.1, page 116 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, and let $W(t)=\sum_je^{-a_jt}$ be the heat trace, where $a_j=-(\rho_j-\tfrac12)^2$ and $\rho_j$ runs over one zero of each pair $\{\rho,1-\rho\}$ of zeros of $\xi$, counted with multiplicity. Unconditionally, that is, using only the verification of RH below height $1000$ (Platt and Trudgian 2021), $W(t)>0$ and $W'(t)<0$ for every $t>0$. Consequently the clock $\xi(\tfrac12)/\xi(\tfrac12+\sqrt u)=\exp\bigl(-\int_0^\infty(1-e^{-ut})\,W(t)\,dt/t\bigr)$ is the Laplace transform of an infinitely divisible, self decomposable law on $[0,\infty)$ with Lévy density $W(t)/t$, and $\bigl(\Xi(\theta)/\Xi(0),\ \xi(\tfrac12)/\xi(\tfrac12+\theta)\bigr)$ is a van Dantzig pair, that is, a pair $(\varphi_1,\varphi_2)$ of characteristic functions with $\varphi_2(\theta)=1/\varphi_1(i\theta)$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:idsd`.

<a id="hankel_positivity_at_one_time_iff_rh"></a>

### Hankel positivity of the heat trace at one time is equivalent to RH

Theorem 5.1.2, page 119 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $W(t)=\sum_je^{-a_jt}$ be the heat trace of the completed Riemann zeta function $\xi$ ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity) and $c_k(t)=(-1)^kW^{(k)}(t)=\sum_ja_j^ke^{-a_jt}$. For any fixed $t>0$ the following are equivalent: both infinite Hankel matrices $[c_{i+j}(t)]_{i,j\ge0}$ and $[c_{i+j+1}(t)]_{i,j\ge0}$ are positive semidefinite; Weil's quadratic form $f\mapsto\sum_\rho\hat f(\gamma_\rho)\overline{\hat f(\bar\gamma_\rho)}$, where $\rho=\tfrac12+i\gamma_\rho$ runs over the nontrivial zeros of $\zeta$ and $\hat f(r)=\int f(x)e^{irx}\,dx$, is positive semidefinite on the span of the test functions $f$ with $\hat f(r)=r^me^{-r^2t/2}$, $m=0,1,2,\dots$; the Riemann Hypothesis.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:hankel`.

<a id="heat_trace_hankel_blocks_positive_definite_to_size_3e7"></a>

### Positive definite Hankel blocks of the heat trace up to size 3 x 10^7

Theorem 5.1.3, page 120 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $W(t)=\sum_je^{-a_jt}$ its heat trace ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity) and $c_k(t)=(-1)^kW^{(k)}(t)$, so that the two Hankel matrices of the heat trace at time $t$ are $[c_{i+j}(t)]_{i,j\ge0}$ and $[c_{i+j+1}(t)]_{i,j\ge0}$. For $N\ge1$ let $T_N=20000N\sqrt{1+\log N}$. If all zeros of $\xi$ below height $T_N$, that is, with ordinate in $(0,T_N)$, lie on the critical line, then the leading $N\times N$ blocks $[c_{i+j}(t)]_{i,j<N}$ and $[c_{i+j+1}(t)]_{i,j<N}$ are positive definite at every $t>0$. By the verification of the Riemann Hypothesis to height $3\times10^{12}$ (Platt and Trudgian 2021) this holds unconditionally for every $N\le3\times10^7$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:blocks`.

<a id="heat_trace_derivative_signs_to_order_1e17"></a>

### Derivative signs of the heat trace up to order 10^17

Theorem 5.1.4, page 120 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $W(t)=\sum_je^{-a_jt}$ be the heat trace of the completed Riemann zeta function $\xi$ ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity). Unconditionally, that is, using only the verification of RH below height $10^{12}$ (Platt and Trudgian 2021), $(-1)^kW^{(k)}(t)>0$ for every $t>0$ and every integer $0\le k\le10^{17}$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:cm`.

<a id="hankel_inertia_counts_zeros_off_critical_line"></a>

### Hankel matrices of Stieltjes moments count the zeros off the critical line (eigenvalue form of the Grommer and Chebotarev criterion)

Theorem 5.1.5, page 122 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $W(t)=\sum_je^{-a_jt}$ be the heat trace of the completed Riemann zeta function $\xi$ ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity), let $m_n=\frac1{n!}\int_0^\infty t^nW(t)\,dt=\sum_ja_j^{-n-1}$ be its Stieltjes moments, and for $N\ge1$ let $\mathcal H_N=[m_{i+j}]_{i,j<N}$ and $\mathcal H'_N=[m_{i+j+1}]_{i,j<N}$ be the Hankel matrices, with $\iota_N$ and $\iota'_N$ the numbers of negative eigenvalues of $\mathcal H_N$ and $\mathcal H'_N$. Let $K\in\{0,1,2,\dots\}\cup\{\infty\}$ be the number of distinct zeros $\rho$ of $\xi$ with $\operatorname{Re}\rho>\tfrac12$ and $\operatorname{Im}\rho>0$. Then $\iota_N$ and $\iota'_N$ are nondecreasing in $N$, $\iota_N\le K$ and $\iota'_N\le K$ for every $N$, and $\iota_N\to K$ and $\iota'_N\to K$ as $N\to\infty$. In particular RH holds if and only if $\iota_N=0$ for every $N$, and if it fails, $\iota_N\ge1$ from some $N$ on; the same holds for $\iota'_N$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:inertia`.

<a id="heat_trace_completely_monotone_near_origin"></a>

### Heat trace completely monotone on (0, 1.29]; every order up to t = 0.019

Theorem 5.1.7, page 125 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $W(t)=\sum_je^{-a_jt}$ be the heat trace of the completed Riemann zeta function $\xi$ ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity). (a) For every integer $k\ge1$ and every $t\in(0,0.019]$, $(-1)^kW^{(k)}(t)>0$; the proof uses no information on the zeros. (b) Unconditionally, that is, using only the verification of RH below height $10^{12}$ (Platt and Trudgian 2021), $W$ is completely monotone on $(0,1.29]$: $(-1)^kW^{(k)}(t)>0$ for every $t\in(0,1.29]$ and every integer $k\ge0$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:cmint`.

<a id="heat_trace_derivative_along_ray_is_gaussian_smoothing_of_zeros"></a>

### Heat trace derivative along a ray is the Gaussian smoothing of the zeros

Proposition 5.1.8, page 130 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $W(t)=\sum_je^{-a_jt}$ its heat trace ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity). Take each $\rho_j$ with positive imaginary part, write $\rho_j=\tfrac12+\delta+i\gamma$ with $\gamma>0$ and $\delta$ real, and put $\tilde\gamma=\gamma-i\delta$, so that $a_j=\tilde\gamma^2$. Let $N(r)$ be the number of zeros of $\xi$ with ordinate in $(0,r]$, counted with multiplicity. Let $u\ge2$, $t>0$ and $k=u^2t$, not necessarily an integer, let $c_k(t)=\sum_ja_j^ke^{-a_jt}$ with $a_j^k$ the principal power (for an integer $k$ this is $(-1)^kW^{(k)}(t)$), and let $f_0=k(\log u^2-1)$. Then $e^{-f_0}c_k(t)=\sum_{|\tilde\gamma-u|\le u/2}e^{-2t(\tilde\gamma-u)^2}\,(1+\varepsilon_\rho)+R$, where the sum runs over the zeros with $|\tilde\gamma-u|\le u/2$, $\varepsilon_\rho$ is the relative error of the Gaussian bump of the zero $\rho$, $R$ is the contribution of the other zeros, and $|\varepsilon_\rho|\le\exp\bigl\{\tfrac{4t}{3u}|\tilde\gamma-u|^3\bigr\}-1$. With $\ell_\pm=u\pm\tfrac12\sqrt{u^2-1}$, $H(r)=u^2\log\frac{r^2+1/4}{u^2}-r^2+\tfrac14+u^2$ and $\eta=1-u^2/(\ell_+^2+\tfrac14)$, which is positive, $|R|\le N(\ell_-)\,e^{tH(\ell_-)}+\bigl(\ell_+^2+\frac1{t\eta}\bigr)e^{tH(\ell_+)}$. No hypothesis on the zeros is used.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `prop:ch12:ray`.

<a id="rh_iff_gaussian_smoothing_of_zeros_nonnegative"></a>

### RH is equivalent to nonnegativity of the Gaussian smoothing of the zeros

Proposition 5.1.9, page 131 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function. Write each zero of $\xi$ with positive imaginary part, counted with multiplicity, as $\rho=\tfrac12+\delta+i\gamma$ with $\gamma>0$ and $\delta$ real, and put $\tilde\gamma=\gamma-i\delta$. For real $u$ and $t>0$ let $\mathcal G(u,t)=\sum e^{-2t(\tilde\gamma-u)^2}$, the sum over these zeros, the Gaussian smoothing of the zeros at resolution $1/(2\sqrt t)$ (it is real, because zeros off the critical line occur in pairs with conjugate values of $\tilde\gamma$). Then the Riemann Hypothesis holds if and only if $\mathcal G(u,t)\ge0$ for every real $u$ and every $t>0$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `prop:ch12:gauss`.

<a id="heat_trace_derivative_signs_to_order_1e25"></a>

### Derivative signs of the heat trace up to order 1.16 x 10^25

Theorem 5.1.10, page 136 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $W(t)=\sum_je^{-a_jt}$ be the heat trace of the completed Riemann zeta function $\xi$ ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity). Let $H=3\times10^{12}$, the height to which RH is verified (Platt and Trudgian 2021). Then $(-1)^kW^{(k)}(t)>0$ for every $t\ge1.29$ and every integer $k$ with $0\le k\le t(H-23)^2$. Consequently $(-1)^kW^{(k)}(t)>0$ for every $t>0$ and every integer $k$ with $0\le k\le1.29(H-23)^2$, and $1.29(H-23)^2>1.16\times10^{25}$.

Source: [`ch/ch12.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch12.tex), label `thm:ch12:window`.

<a id="off_line_perturbation_passing_finite_heat_trace_tests"></a>

### Four off line zeros that pass the finite heat trace tests

Theorem 5.2.1, page 137 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $W(t)=\sum_je^{-a_jt}$ its heat trace ($a_j=-(\rho_j-\tfrac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, counted with multiplicity); the clock is $\xi(\tfrac12)/\xi(\tfrac12+\sqrt u)$. (a) There are ordinates $1000<A<B\le2000$ of zeros on the critical line, with $B-A\le2$, such that, with $g=(A+B)/2$ and $w=\tfrac14+ig$, $F(s)=\frac{\xi(\frac12+s)}{\xi(\frac12)}\,\frac{(1-s^2/w^2)(1-s^2/\bar w^2)}{(1+s^2/A^2)(1+s^2/B^2)}$ is an even entire function with the order, type and reflection symmetry of $\xi(\tfrac12+s)$ and with four zeros off the critical line, and the following hold. The reciprocal clock $1/F(\sqrt u)$ is the Laplace transform of an infinitely divisible, self decomposable law, and its heat trace $W_F(t)=W(t)+2\operatorname{Re}e^{-at}-e^{-A^2t}-e^{-B^2t}$, $a=-w^2=g^2-\tfrac1{16}-\tfrac{ig}2$, which is $t$ times its Lévy density, is strictly log convex. $W_F$ has the same small time expansion as $W$ through the constant term, and $1/F$ has the same nearest poles as the reciprocal Xi function. (b) Let $N\ge1$ and $R\ge T_N$, where $T_N=20000N\sqrt{1+\log N}$, and suppose that every zero of $\xi$ with ordinate in $(0,R+50]$ lies on the critical line. Then the pair can instead be chosen with $R<A\le B\le R+50$ and $B-A<2$ (if $A=B$, a multiple zero supplies both), so that $F$ keeps every property listed in (a) and the leading $N\times N$ Hankel blocks $[c_{i+j}]_{i,j<N}$ and $[c_{i+j+1}]_{i,j<N}$ of $W_F$, with $c_k=(-1)^kW_F^{(k)}$, remain positive definite at every $t>0$. If moreover $R\ge10^{12}$, the derivative signs are also preserved: $(-1)^kW_F^{(k)}(t)>0$ for every $t>0$ and every integer $0\le k\le10^{17}$. For every $N\le3\times10^7$ the verification of RH to height $3\times10^{12}$ (Platt and Trudgian 2021) establishes the hypothesis on the zeros for $R=\max(T_N,10^{12})$, so that such an $F$ with both properties exists unconditionally.

Source: [`ch/ch13.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch13.tex), label `thm:ch13:perturb`.

<a id="bernstein_pick_route_closed_by_certified_zero"></a>

### Certified zero of a Mellin transform closes the Bernstein Pick route to RH

Theorem 5.2.2, page 143 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Xi(z)=\xi(\tfrac12+iz)$, $\Phi(u)=2\sum_{n\ge1}\bigl(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\bigr)\exp(-\pi n^2e^{2u})$ the Pólya kernel, so that $\Xi(\theta)=\int_{\mathbb{R}}\Phi(u)\cos(\theta u)\,du$, and let $X$ be a random variable with the Pólya density $\Phi/\Xi(0)$. Put $M(s)=\mathbb{E}|X|^{2s}$, analytic for $\operatorname{Re}s>-\tfrac12$ and continued meromorphically to $\operatorname{Re}s<0$ by $M(s)=\frac{2\,\Gamma(2s+1)\sin(-\pi s)}{\pi\,\Xi(0)}\int_0^\infty\Xi(\theta)\,\theta^{-2s-1}\,d\theta$ for $-\tfrac12<\operatorname{Re}s<0$, and let $\varphi_I(s)=4(s-\tfrac12)M(s-1)/M(s)$, the canonical interpolant of the moment ratios. Then $M$ has exactly one zero in the disc $|s-s_1|<10^{-8}$ with centre $s_1=-8.494182256992559+2.937767068567584i$; the zero is simple, and $M(s-1)$ does not vanish in the disc. Hence $\varphi_I$ has a nonreal pole, and every Bernstein interpolant of the moment ratios (a Bernstein function that takes the values $\varphi_I(n)$ at the integers $n\ge1$) agrees with $\varphi_I$ on $(1,\infty)$ and cannot be continued through that pole, so none is a complete Bernstein function. For every $0<q<1$ the sequence $n\varphi_I(n)/(n-q)$ has no Bernstein interpolant either.

Source: [`ch/ch13.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch13.tex), label `thm:ch13:mellin`.

## Chapter 6. The reciprocal zeta function

<a id="euler_factor_thorin_type_density"></a>

### Thorin type density of an Euler factor

Proposition 6.1.2, page 159 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $0<q<1$ and $\lambda>0$, and let $E(s)=1-q\,e^{-\lambda\sqrt s}$, where $\sqrt s$ is the principal square root on the cut plane $\mathbb{C}\setminus(-\infty,0]$. Then $E$ is analytic and has no zeros on the cut plane, and there $\log\frac{E(s)}{1-q}=-\int_0^\infty\log\bigl(1+\frac st\bigr)\,u(t)\,dt$ with $u(t)=-\frac{q\lambda}{2\pi\sqrt t}\,\frac{\cos(\lambda\sqrt t)-q}{1-2q\cos(\lambda\sqrt t)+q^2}=-\frac{\lambda}{2\pi\sqrt t}\sum_{k\ge1}q^k\cos(k\lambda\sqrt t)$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:prop:euler`.

<a id="lai_vere_jones_identity"></a>

### Lai and Vere-Jones identity (the exponential as a Möbius sum of hyperbolic secants)

Proposition 6.1.3, page 160 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\mu$ be the Möbius function and $\chi_4$ the nontrivial Dirichlet character modulo $4$. For $\operatorname{Re}x>0$, $\frac12\operatorname{sech}x=\sum_{m\ \mathrm{odd}}\chi_4(m)\,e^{-mx}$ and $e^{-x}=\frac12\sum_{m\ge1}\mu(m)\chi_4(m)\operatorname{sech}(mx)$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:prop:lvj`.

<a id="vere_jones_representation_reciprocal_zeta"></a>

### Vere-Jones representation of the reciprocal zeta function by hyperbolic secants

Theorem 6.1.5, page 160 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\zeta$ be the Riemann zeta function, $\mu$ the Möbius function, $\chi_4$ the nontrivial Dirichlet character modulo $4$, $\sqrt s$ the principal square root, and $T_\lambda$ the exit time of $\sqrt2B$ from $(-\lambda,\lambda)$ for a standard Brownian motion $B$ started at $0$, so that $\mathbb{E}\,e^{-sT_\lambda}=\operatorname{sech}(\lambda\sqrt s)$. For $s$ in the region $H_{1/2}=\{s:\operatorname{Re}\sqrt s>\frac12\}=\{x+iy:\ x>\frac14-y^2\}$, $\frac{1}{\zeta(\frac12+\sqrt s)}=1+\sum_{N\ge2}A_N\operatorname{sech}(\sqrt s\log N)=1+\sum_{N\ge2}A_N\,\mathbb{E}\,e^{-sT_{\log N}}$, where $A_N=\frac12\,\mu(m)\chi_4(m)\mu(n)\,n^{-1/2}$ if $N=n^m$ with $n\ge2$ squarefree and $m$ odd, and $A_N=0$ otherwise. In particular $A_N=\mu(N)/(2\sqrt N)$ for squarefree $N\ge2$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:thm:sech`.

<a id="exponential_as_eta_weighted_sum_of_log_blocks"></a>

### Logarithmic twin of the Lai and Vere-Jones identity (the exponential as a sum of logarithms with Dirichlet eta weights)

Proposition 6.1.6, page 161 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\mu$ be the Möbius function and $\eta(z)=\sum_{k\ge1}(-1)^{k+1}k^{-z}=(1-2^{1-z})\zeta(z)$ Dirichlet's eta function. Define $d_n$ by $\sum_nd_nn^{-v}=1/\eta(1+v)$, that is $d_n=\sum_{2^jm=n}\mu(m)/m$. Then $|d_n|\le1$ for all $n\ge1$, and for $\operatorname{Re}x>0$, $e^{-x}=\sum_{m\ge1}d_m\log\bigl(1+e^{-mx}\bigr)$.

Supersedes: Proposition 27.6 of Polson (2026) states the same identity with the weaker bound $|d_n|\le2$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:prop:eta`.

<a id="sawtooth_thorin_measure_of_vere_jones_block"></a>

### Sawtooth Thorin measure of a Vere-Jones block

Lemma 6.1.7, page 161 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\lambda>0$, let $s$ lie in the cut plane $\mathbb{C}\setminus(-\infty,0]$ with $\sqrt s$ the principal square root, and put $\sigma_k(\lambda)=\bigl((2k+1)\pi/\lambda\bigr)^2$ for $k\ge0$. Then $\log\bigl(1+e^{-\lambda\sqrt s}\bigr)=\log2-\int\log\bigl(1+\frac st\bigr)\,U_\lambda(dt)$ with $U_\lambda(dt)=\frac{\lambda}{4\pi\sqrt t}\,dt-\sum_{k\ge0}\delta_{\sigma_k(\lambda)}$, where $\delta_a$ is the unit point mass at $a$. Equivalently $U_\lambda(dt)=-d_t\Sigma_\lambda(\sqrt t)$, where $\Sigma_\lambda(\tau)=\#\{k\ge0:(2k+1)\pi/\lambda\le\tau\}-\lambda\tau/(2\pi)$ is a sawtooth function with $|\Sigma_\lambda|\le\frac12$. So the reciprocal block $2/(1+e^{-\lambda\sqrt s})$ has the signed Thorin type measure $-U_\lambda=\sum_{k\ge0}\delta_{\sigma_k(\lambda)}-\frac{\lambda}{4\pi\sqrt t}\,dt$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:lem:saw`.

<a id="vere_jones_representation_log_zeta"></a>

### Vere-Jones representation of log zeta by sawtooth blocks

Theorem 6.1.8, page 162 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Lambda$ be the von Mangoldt function and $\mu$ the Möbius function, and put $a_n=\Lambda(n)/(\sqrt n\log n)$, $d_m=\sum_{2^jk=m}\mu(k)/k$ (so that $\sum_m d_mm^{-v}=1/\eta(1+v)$ for Dirichlet's eta function $\eta$, and $|d_m|\le1$), $B_N=\sum_{n^m=N}d_m\,a_n$ (the sum over $n\ge2$ and $m\ge1$ with $n^m=N$) and $\sigma_k(\lambda)=\bigl((2k+1)\pi/\lambda\bigr)^2$. For $s$ in the region $H_{1/2}=\{s:\operatorname{Re}\sqrt s>\frac12\}$, with $\sqrt s$ the principal square root, $\log\zeta\bigl(\tfrac12+\sqrt s\bigr)=\sum_{N\ge2}B_N\log\bigl(1+N^{-\sqrt s}\bigr)$ and $\frac{d}{ds}\log\zeta\bigl(\tfrac12+\sqrt s\bigr)=\sum_{N\ge2}B_N\Bigl[\sum_{k\ge0}\frac{1}{s+\sigma_k(\log N)}-\frac{\log N}{4\sqrt s}\Bigr]$. The first values are $B_p=p^{-1/2}$ for primes $p$, $B_4=\frac14+\frac1{2\sqrt2}$, $B_8=-\frac1{6\sqrt2}$ and $B_9=\frac16+\frac{1}{2\sqrt3}$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:thm:saw`.

<a id="thorin_measure_normalised_reciprocal_zeta_under_rh"></a>

### Thorin type measure of the normalised reciprocal zeta function under RH

Theorem 6.1.10, page 163 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: conditional on The Riemann Hypothesis: every nontrivial zero of the Riemann zeta function has real part $\frac12$.

Statement. Assume the Riemann Hypothesis. Let $\zeta$ be the Riemann zeta function, $\varphi(s)=\zeta(\frac12)/\{(1-4s)\,\zeta(\frac12+\sqrt s)\}$ (so $\varphi(0)=1$), $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ the Riemann Siegel theta function, and $N(T)$ the number of zeros $\frac12+i\gamma$ of $\zeta$ with $0<\gamma\le T$, counted with multiplicity, as are the sums over $\gamma$ below. Then for $s$ in the cut plane $\mathbb{C}\setminus(-\infty,0]$, $\log\varphi(s)=-\int_{(0,\infty)}\log\bigl(1+\tfrac st\bigr)\,dM(t)$ with $M(t)=N(\sqrt t)-\theta(\sqrt t)/\pi$, that is $dM(t)=\sum_{\gamma>0}\delta_{\gamma^2}-\frac{\theta'(\sqrt t)}{2\pi\sqrt t}\,dt$ and $-\frac{\varphi'}{\varphi}(s)=\int\frac{dM(t)}{s+t}$. By the Riemann von Mangoldt formula $M(t)=1+S(\sqrt t)$ for $\sqrt t$ not an ordinate, where $S(T)=\pi^{-1}\arg\zeta(\frac12+iT)$ is the argument function and $S(0^+)=-1$. In particular $dM$ is not a positive measure, so $\varphi$ is not the Laplace transform of a GGC even under RH.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:thm:thorin`.

<a id="canonical_function_normalised_reciprocal_zeta"></a>

### Canonical function of the normalised reciprocal zeta function, prime and zero forms

Theorem 6.1.11, page 164 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\zeta$ be the Riemann zeta function, $\Lambda$ the von Mangoldt function, $\varphi(s)=\zeta(\frac12)/\{(1-4s)\,\zeta(\frac12+\sqrt s)\}$ the normalised reciprocal zeta function, and $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ the Riemann Siegel theta function. Define, with no hypothesis, $w(x)=e^{x/4}-P(x)$ with $P(x)=\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,e^{-(\log n)^2/4x}$ for $x>0$. Then $-\varphi'/\varphi(s)=\int_0^\infty e^{-sx}w(x)\,dx$ for real $s>\frac14$, and for every $x>0$, $w(x)=W(x)-W_0(x)$ with $W(x)=\frac12\sum_\rho e^{x(\rho-\frac12)^2}$ and $W_0(x)=\frac1\pi\int_0^\infty e^{-x\tau^2}\theta'(\tau)\,d\tau$, the sum running over all nontrivial zeros $\rho$ of $\zeta$ with multiplicity. Under RH, $W(x)=\sum_{\gamma>0}e^{-x\gamma^2}$ and $w(x)=\int e^{-tx}\,dM(t)$, where $dM(t)=\sum_{\gamma>0}\delta_{\gamma^2}-\frac{\theta'(\sqrt t)}{2\pi\sqrt t}\,dt$ and $\gamma>0$ runs over the ordinates of the zeros $\frac12+i\gamma$ (counted with multiplicity). The identity $w=W-W_0$ is Weil's explicit formula for the Gaussian test function, in the arithmetic form of the heat trace stated as Theorem 4.4 of Polson (2026). Moreover $w(0^+)=1$ and $w(x)\sim c/(2\sqrt{\pi x})$ as $x\to\infty$, with $c=-\theta'(0)=\zeta'(\frac12)/\zeta(\frac12)\approx2.68609$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:thm:w`.

<a id="theta_series_for_canonical_function_reciprocal_zeta"></a>

### Theta function series for the canonical function w

Proposition 6.1.13, page 166 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Lambda$ be the von Mangoldt function, $\mu$ the Möbius function, $d_m=\sum_{2^jk=m}\mu(k)/k$, $B_N=\sum_{n^m=N}d_m\Lambda(n)/(\sqrt n\log n)$ (the sum over $n\ge2$ and $m\ge1$ with $n^m=N$), $\sigma_k(\lambda)=\bigl((2k+1)\pi/\lambda\bigr)^2$, and $h_\lambda(x)=\sum_{k\ge0}e^{-x\sigma_k(\lambda)}-\frac{\lambda}{4\sqrt{\pi x}}=\frac{\lambda}{2\sqrt{\pi x}}\sum_{j\ge1}(-1)^je^{-j^2\lambda^2/4x}$. Let $w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,e^{-(\log n)^2/4x}$. Then for every $x>0$, with absolutely convergent sums, $w(x)=e^{x/4}+\sum_{N\ge2}B_N\,h_{\log N}(x)$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:prop:theta`.

<a id="sawtooth_expansion_argument_function_gaussian_average"></a>

### Sawtooth expansion of the argument function after Gaussian averaging, under RH

Corollary 6.1.14, page 167 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: conditional on The Riemann Hypothesis: every nontrivial zero of the Riemann zeta function has real part $\frac12$.

Statement. Assume the Riemann Hypothesis. Let $S(\tau)=\pi^{-1}\arg\zeta(\frac12+i\tau)$ be the argument function of the Riemann zeta function, let $B_N$ be the coefficients of the Vere-Jones representation of $\log\zeta$ ($B_N=\sum_{n^m=N}d_m\Lambda(n)/(\sqrt n\log n)$ with $d_m=\sum_{2^jk=m}\mu(k)/k$, $\Lambda$ the von Mangoldt function and $\mu$ the Möbius function), let $\Sigma_\lambda(\tau)=\#\{k\ge0:(2k+1)\pi/\lambda\le\tau\}-\lambda\tau/(2\pi)$ be the sawtooth function, and let $h_\lambda(x)=\sum_{k\ge0}e^{-x((2k+1)\pi/\lambda)^2}-\frac{\lambda}{4\sqrt{\pi x}}$. Then for every $x>0$, $\int_0^\infty S(\tau)\,2\tau e^{-x\tau^2}\,d\tau=\frac{e^{x/4}-1}{x}+\sum_{N\ge2}B_N\int_0^\infty\Sigma_{\log N}(\tau)\,2\tau e^{-x\tau^2}\,d\tau$, and each integral on the right equals $h_{\log N}(x)/x$.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:cor:Ssaw`.

<a id="rh_criterion_from_canonical_function"></a>

### Criterion for RH from the canonical function w

Theorem 6.1.15, page 167 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Lambda$ be the von Mangoldt function and $w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,e^{-(\log n)^2/4x}$ for $x>0$. The following are equivalent. (i) The Riemann Hypothesis. (ii) There is a real signed Borel measure $m$ on $[0,\infty)$ with $\int|m|(dt)/(1+t)<\infty$ and $w(x)=\int e^{-tx}m(dt)$ for all $x>0$. In that case $m=dM$, where $dM(t)=\sum_{\gamma>0}\delta_{\gamma^2}-\frac{\theta'(\sqrt t)}{2\pi\sqrt t}\,dt$, the sum is over the ordinates $\gamma>0$ of the zeros $\frac12+i\gamma$ of the Riemann zeta function counted with multiplicity, and $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ is the Riemann Siegel theta function.

Source: [`ch/rz1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz1.tex), label `rz:thm:crit`.

<a id="thorin_measure_of_the_gamma_factor"></a>

### Thorin measure of the gamma factor as a renormalised product of Poisson kernels

Theorem 6.2.1, page 170 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $G(s)=\Gamma(\frac14+\frac{\sqrt s}2)\,\pi^{-\sqrt s/2}/\Gamma(\frac14)$ be the gamma factor, $c_n=2n+\frac12$ for $n\ge0$, $b=\frac12(\gamma_E+\log\pi)$ with $\gamma_E$ Euler's constant, and $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ the Riemann Siegel theta function. For $s$ in the cut plane $\mathbb{C}\setminus(-\infty,0]$, $G(s)=e^{-b\sqrt s}\prod_{n\ge0}\frac{e^{\sqrt s/(2n+2)}}{1+\sqrt s/c_n}$, and the Thorin type measure $U_G$ of $G$, defined by $G(s)=\exp\bigl(-\int_{(0,\infty)}\log(1+\frac st)\,U_G(dt)\bigr)$, is $U_G(dt)=\frac{1}{2\pi\sqrt t}\Bigl\{b+\sum_{n\ge0}\Bigl(\frac{c_n}{t+c_n^2}-\frac{1}{2n+2}\Bigr)\Bigr\}dt=-\frac{\theta'(\sqrt t)}{2\pi\sqrt t}\,dt$. Each factor $(1+\sqrt s/\kappa)^{-1}$ with $\kappa>0$ is the Laplace transform of a generalized gamma convolution (GGC) with positive Thorin density $q_\kappa(t)=\frac{\kappa}{2\pi\sqrt t\,(t+\kappa^2)}$, a Poisson kernel in $\sqrt t$; each factor $e^{a\sqrt s}$ with $a>0$ has the negative Thorin density $-a/(2\pi\sqrt t)$.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:thm:gamma`.

<a id="heat_trace_closed_form_and_formal_thorin_measure_of_reciprocal_xi"></a>

### Closed form heat trace over the zeros and formal Thorin measure of the reciprocal Xi function

Theorem 6.2.3, page 172 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\zeta$ be the Riemann zeta function, $\xi(z)=\frac12z(z-1)\pi^{-z/2}\Gamma(z/2)\zeta(z)$ Riemann's xi function, $W(x)=\frac12\sum_\rho e^{x(\rho-\frac12)^2}$ for $x>0$ (the sum over the nontrivial zeros $\rho$ of $\zeta$ with multiplicity), $\Lambda$ the von Mangoldt function and $\mu$ the Möbius function. Put $d_m=\sum_{2^jk=m}\mu(k)/k$, $B_N=\sum_{n^m=N}d_m\Lambda(n)/(\sqrt n\log n)$, $\sigma_k(\lambda)=\bigl((2k+1)\pi/\lambda\bigr)^2$, $h_\lambda(x)=\sum_{k\ge0}e^{-x\sigma_k(\lambda)}-\frac{\lambda}{4\sqrt{\pi x}}$, $b=\frac12(\gamma_E+\log\pi)$ with $\gamma_E$ Euler's constant, and $c_n=2n+\frac12$. Then for every $x>0$, with absolutely convergent sums, $W(x)=e^{x/4}+\sum_{N\ge2}B_N\,h_{\log N}(x)-\frac{b}{2\sqrt{\pi x}}+\sum_{n\ge0}\Bigl\{\frac{1}{(4n+4)\sqrt{\pi x}}-\frac12e^{c_n^2x}\operatorname{erfc}\bigl(c_n\sqrt x\bigr)\Bigr\}$. Equivalently, $W(x)$ is the Laplace transform, taken term by term, of the formal series $U_R=\delta_{-1/4}+\sum_{N\ge2}B_N\Bigl\{\sum_{k\ge0}\delta_{\sigma_k(\log N)}-\frac{\log N}{4\pi\sqrt t}\,dt\Bigr\}-\frac{1}{2\pi\sqrt t}\Bigl\{b+\sum_{n\ge0}\Bigl(\frac{c_n}{t+c_n^2}-\frac1{2n+2}\Bigr)\Bigr\}dt$, the formal Thorin measure of $R(s)=\xi(\frac12)/\xi(\frac12+\sqrt s)$. Under RH, $R$ is a GGC (generalized gamma convolution), and its Thorin measure $\sum_{\gamma>0}\delta_{\gamma^2}=dN(\sqrt t)$, where $N(T)$ counts the zeros $\frac12+i\gamma$ with $0<\gamma\le T$, is a positive measure whose Laplace transform is also $W$.

Supersedes: Theorem 27.18 of Polson (2026) states that the series $U_R$ is the Thorin measure of $R$ in the weak sense of Laplace transforms on $x>0$, and that $U_R=\sum_{\gamma>0}\delta_{\gamma^2}$ under RH.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:thm:RLVJ`.

<a id="sine_squared_lemma_polson_2017"></a>

### Sine squared lemma (Polson 2017, Lemma 1)

Lemma 6.2.4, page 174 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Suppose that $f$ is a function, $\alpha$ a real number and $\mu$ a measure on $(0,\infty)$ with $\int_0^\infty(1\wedge x)\,e^{-\alpha x}\,\mu(dx)<\infty$ such that $\log\frac{f(\alpha+s)}{f(\alpha)}-s\frac{f'(\alpha)}{f(\alpha)}=\int_0^\infty\bigl(e^{-sx}-1+sx\bigr)e^{-\alpha x}\,\frac{\mu(dx)}x$ for $s\ge0$. Then, for $s\ge0$, $\log\frac{f(\alpha)}{f(\alpha+s)}+s\frac{f'(\alpha)}{f(\alpha)}=-\int_0^\infty\bigl(1-e^{-\frac12s^2t}\bigr)\frac{\nu_\alpha(t)}t\,dt$ with $\nu_\alpha(t)=\frac{1}{\sqrt{2\pi}}\int_0^\infty e^{-tz}\Bigl\{\int_0^\infty2\sin^2\bigl(x\sqrt{z/2}\bigr)e^{-\alpha x}\,\mu(dx)\Bigr\}\frac{dz}{\sqrt{\pi z}}$. The map $\mu\mapsto\nu_\alpha$ is linear and positivity preserving.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:lem:P1`.

<a id="sine_squared_lemma_clock_variable"></a>

### Sine squared lemma in the variable s = u squared (Polson 2017, Lemma 1)

Lemma 6.2.5, page 174 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2017, version 24 of 2026), On Hilbert's 8th problem, arXiv:1708.02653.

Status: proved.

Statement. Let $\tilde\mu$ be a signed measure on $(0,\infty)$ with $\int(1\wedge x)\,|\tilde\mu|(dx)<\infty$, and let $f$ be a function with $\log\frac{f(u)}{f(0)}-u\frac{f'(0)}{f(0)}=\int_0^\infty\bigl(e^{-ux}-1+ux\bigr)\frac{\tilde\mu(dx)}x$ for real $u\ge0$. Then for $s$ in the cut plane $\mathbb{C}\setminus(-\infty,0]$, $\log\frac{f(0)}{f(\sqrt s)}+\sqrt s\,\frac{f'(0)}{f(0)}=-\int_0^\infty\log\Bigl(1+\frac st\Bigr)U(t)\,dt$ with $U(t)=\frac{1}{2\pi\sqrt t}\int_0^\infty2\sin^2\Bigl(\frac{x\sqrt t}2\Bigr)\tilde\mu(dx)$. If $\tilde\mu\ge0$ then $U\ge0$, and $s\mapsto f(0)\,e^{\sqrt s f'(0)/f(0)}/f(\sqrt s)$ is the Laplace transform of a GGC (generalized gamma convolution). Since $\sqrt s=\int\log(1+s/t)\,\frac{dt}{2\pi\sqrt t}$, the Thorin type measure of $f(0)/f(\sqrt s)$ itself is $U(t)\,dt+\frac{f'(0)}{f(0)}\frac{dt}{2\pi\sqrt t}$.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:lem:sin2`.

<a id="reciprocal_zeta_right_of_line_one_tilted_ggc"></a>

### Reciprocal zeta to the right of the line Re z = 1 is a GGC after a stable tilt

Proposition 6.2.6, page 175 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\zeta$ be the Riemann zeta function and $\Lambda$ the von Mangoldt function. For $\beta>1$ and $s$ in the cut plane $\mathbb{C}\setminus(-\infty,0]$, $\frac{\zeta(\beta)}{\zeta(\beta+\sqrt s)}\,e^{\sqrt s\,\zeta'(\beta)/\zeta(\beta)}=\exp\Bigl(-\int_0^\infty\log\Bigl(1+\frac st\Bigr)U_\beta(t)\,dt\Bigr)$ with $U_\beta(t)=\frac{1}{2\pi\sqrt t}\sum_{n\ge2}\frac{\Lambda(n)}{n^\beta}\,2\sin^2\Bigl(\frac{\sqrt t\log n}2\Bigr)\ge0$, so the left side is the Laplace transform of a generalized gamma convolution (GGC) with Thorin density $U_\beta$. Adding the tilt, the Thorin type measure of $\zeta(\beta)/\zeta(\beta+\sqrt s)$ is $-\frac{1}{2\pi\sqrt t}\sum_{n\ge2}\Lambda(n)n^{-\beta}\cos(\sqrt t\log n)\,dt$.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:zetaright`.

<a id="shifted_vere_jones_block_positive_thorin_density"></a>

### Positive Thorin density of a shifted Vere-Jones block

Proposition 6.2.7, page 176 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\lambda>0$, $\alpha>0$ and $q=e^{-\lambda\alpha}$, let $f(u)=1+e^{-\lambda(\alpha+u)}$, for which $\log\frac{f(u)}{f(0)}-u\frac{f'(0)}{f(0)}=\int_0^\infty\bigl(e^{-ux}-1+ux\bigr)\frac{\tilde\mu_{\lambda,\alpha}(dx)}x$ with the alternating lattice $\tilde\mu_{\lambda,\alpha}=\lambda\sum_{k\ge1}(-1)^{k+1}q^k\delta_{k\lambda}$, and let $\sigma_k(\lambda)=\bigl((2k+1)\pi/\lambda\bigr)^2$. The sine squared transform of $\tilde\mu_{\lambda,\alpha}$ is positive: $U_{\lambda,\alpha}(t)=\frac{1}{2\pi\sqrt t}\int_0^\infty2\sin^2\bigl(\frac{x\sqrt t}2\bigr)\tilde\mu_{\lambda,\alpha}(dx)=\frac{\lambda}{2\pi\sqrt t}\,\frac{q(1-q)\bigl(1-\cos\lambda\sqrt t\bigr)}{(1+q)\bigl(1+2q\cos\lambda\sqrt t+q^2\bigr)}\ge0$, so $\frac{1+e^{-\lambda\alpha}}{1+e^{-\lambda(\alpha+\sqrt s)}}\,e^{-\sqrt s\,\lambda q/(1+q)}$ is the Laplace transform of a generalized gamma convolution (GGC) with Thorin density $U_{\lambda,\alpha}$. As $\alpha\downarrow0$, $U_{\lambda,\alpha}(t)\,dt\to\sum_{k\ge0}\delta_{\sigma_k(\lambda)}$ vaguely, the tilt density tends to $-\frac{\lambda}{4\pi\sqrt t}$, and the Thorin type measure of $2/(1+e^{-\lambda\sqrt s})$ is recovered as $\sum_{k\ge0}\delta_{\sigma_k(\lambda)}-\frac{\lambda}{4\pi\sqrt t}\,dt$.

Supersedes: Proposition 28.4 of Polson (2026) states that $\frac{1+e^{-\lambda(\alpha+\sqrt s)}}{1+e^{-\lambda\alpha}}\,e^{-\sqrt s\,\lambda q/(1+q)}$, the block rather than its reciprocal, is the Laplace transform of a generalized gamma convolution (GGC), and that $U_{\lambda,\alpha}(t)\,dt\to\sum_{k\ge0}\delta_{\sigma_k(\lambda)}$ weakly as $\alpha\downarrow0$, where the book has vaguely because the limit has infinite mass.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:block`.

<a id="cosine_levy_series_of_vere_jones_blocks_diverges"></a>

### Divergence of the cosine Lévy series of the Vere-Jones blocks

Proposition 6.2.9, page 177 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For a measure $U$ on $(0,\infty)$ put $K_U(x)=2\int e^{-\sqrt t\,x}\,U(dt)$, let $U_\lambda(dt)=\frac{\lambda}{4\pi\sqrt t}\,dt-\sum_{k\ge0}\delta_{((2k+1)\pi/\lambda)^2}$ be the sawtooth measure, for which $K_{U_\lambda}(x)=\frac{\lambda}{\pi x}-\frac{1}{\sinh(\pi x/\lambda)}$, and let $B_N=\sum_{n^m=N}d_m\Lambda(n)/(\sqrt n\log n)$ with $d_m=\sum_{2^jk=m}\mu(k)/k$ be the coefficients of the Vere-Jones representation of $\log\zeta$ ($\Lambda$ the von Mangoldt function, $\mu$ the Möbius function). Then for every $x>0$ the series $\sum_{N\ge2}B_N\,K_{U_{\log N}}(x)$, summed in order of increasing $N$, diverges.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:diverge`.

<a id="meixner_density_not_polya_frequency_function"></a>

### The Meixner density of the gamma factor is not a Pólya frequency function of infinite order

Proposition 6.2.10, page 178 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $p(\tau)=|\Gamma(\frac14+\frac{i\tau}2)|^2/(2^{3/2}\pi^{3/2})$ be the Meixner density, a probability density on $\mathbb{R}$ with characteristic function $\int p(\tau)e^{iv\tau}\,d\tau=\operatorname{sech}(v)^{1/2}$. Then $p$ is not a Pólya frequency function of infinite order, that is, in the sense of Schoenberg, the kernel $p(x-y)$ is not totally positive of every order. The proof applies Schoenberg's theorem: the reciprocal of the characteristic function of a Pólya frequency density of infinite order is an entire function (of the Laguerre Pólya class), whereas here the reciprocal is $\cosh(v)^{1/2}$, which has branch points at the simple zeros $v=i(k+\frac12)\pi$ of $\cosh$ and is not entire.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:meixstruct`.

<a id="gamma_factor_bernstein_function_and_meixner_score"></a>

### Bernstein function of the gamma factor: theta derivative and Meixner score

Proposition 6.2.11, page 179 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\psi=\Gamma'/\Gamma$ be the digamma function, $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ the Riemann Siegel theta function, $p(\tau)=|\Gamma(\frac14+\frac{i\tau}2)|^2/(2^{3/2}\pi^{3/2})$ the Meixner density, $c_k=2k+\frac12$, and $\nu(dy)=\frac{e^{y/2}}{\sinh y}\,dy=2\sum_{k\ge0}e^{-c_ky}\,dy$ on $(0,\infty)$. Then $\Phi_\Gamma(\lambda)=\psi(\frac14+\frac\lambda2)-\psi(\frac14)=\int_0^\infty(1-e^{-\lambda y})\,\nu(dy)$ is a Bernstein function, and its boundary values on the imaginary axis are $\operatorname{Re}\Phi_\Gamma(i\tau)=2\bigl(\theta'(\tau)-\theta'(0)\bigr)=\int_0^\infty(1-\cos\tau y)\,\nu(dy)$ and $\operatorname{Im}\Phi_\Gamma(i\tau)=-\frac{d}{d\tau}\log p(\tau)=\int_0^\infty\sin(\tau y)\,\nu(dy)$. Consequently $2\theta(\tau)+\tau\log\pi$ and $\log p(\tau)$ are, up to constants, the imaginary and real parts of $2\log\Gamma(\frac14+\frac{i\tau}2)$, the boundary value of a function analytic in the right half plane, so $2\theta(\tau)+\tau\log\pi$ is the harmonic conjugate of the Meixner log density. The smooth heat trace is the heat trace of $\nu$: $W_0(x)=\frac1\pi\int_0^\infty e^{-x\tau^2}\theta'(\tau)\,d\tau=\frac{\theta'(0)}{2\sqrt{\pi x}}+\frac{1}{4\sqrt{\pi x}}\int_0^\infty\bigl(1-e^{-y^2/4x}\bigr)\,\nu(dy)$.

Supersedes: Proposition 28.16 of Polson (2026) says that the smooth part $\theta/\pi$ of the zero counting function, rather than $2\theta(\tau)+\tau\log\pi$, is the harmonic conjugate of the Meixner log density.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:bernstein`.

<a id="tilted_reciprocal_gamma_factor_is_ggc"></a>

### The reciprocal gamma factor is a GGC after a stable tilt

Theorem 6.2.12, page 180 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\theta(\tau)=\operatorname{Im}\log\Gamma(\frac14+\frac{i\tau}2)-\frac\tau2\log\pi$ be the Riemann Siegel theta function, $\psi=\Gamma'/\Gamma$, $c=-\theta'(0)=\frac12\bigl(\log\pi-\psi(\frac14)\bigr)\approx2.68609$ (which equals $\zeta'(\frac12)/\zeta(\frac12)$), $G(s)=\Gamma(\frac14+\frac{\sqrt s}2)\,\pi^{-\sqrt s/2}/\Gamma(\frac14)$ the gamma factor, and $c_k=2k+\frac12$. Then for $s$ in the cut plane $\mathbb{C}\setminus(-\infty,0]$, $\frac{e^{-c\sqrt s}}{G(s)}=\frac{\Gamma(\frac14)\,\pi^{\sqrt s/2}\,e^{-c\sqrt s}}{\Gamma(\frac14+\frac{\sqrt s}2)}=\exp\Bigl(-\int_0^\infty\log\Bigl(1+\frac st\Bigr)U_\Gamma(t)\,dt\Bigr)$ with $U_\Gamma(t)=\frac{\theta'(\sqrt t)-\theta'(0)}{2\pi\sqrt t}=\frac{1}{2\pi\sqrt t}\sum_{k\ge0}\frac{t}{c_k\,(c_k^2+t)}\ge0$. So $e^{-c\sqrt s}/G(s)$ is the Laplace transform of a generalized gamma convolution (GGC), and its Thorin distribution function is $\bigl(\theta(\sqrt t)-\theta'(0)\sqrt t\bigr)/\pi$, a convex increasing function of $\sqrt t$.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:thm:gammaggc`.

<a id="polya_autocorrelation_meixner_factorisation"></a>

### Factorisation of the Pólya autocorrelation through the Meixner law

Proposition 6.2.14, page 183 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi$ be Riemann's xi function, $\Xi(\tau)=\xi(\frac12+i\tau)$, and $\Phi$ the Pólya kernel, the even function with $\Xi(\tau)=\int_{\mathbb{R}}\Phi(u)e^{i\tau u}\,du$. Let $p(\tau)=|\Gamma(\frac14+\frac{i\tau}2)|^2/(2^{3/2}\pi^{3/2})$ be the Meixner density, with characteristic function $\int p(\tau)e^{iv\tau}d\tau=\operatorname{sech}(v)^{1/2}$, and put $\mathcal{C}(v)=\int p(\tau)\,|\zeta(\frac12+i\tau)|^2e^{iv\tau}\,d\tau$. Then the integral $\mathcal{C}(v)$ converges absolutely, and $(\Phi*\Phi)(v)=2^{-3/2}\Bigl(\frac{d^2}{dv^2}-\frac14\Bigr)^2\mathcal{C}(v)$. Formally $\mathcal{C}=\operatorname{sech}(\cdot)^{1/2}*\mathcal{A}$ with $\mathcal{A}=\sum_{m,n\ge1}\frac{\delta_{\log(m/n)}}{\sqrt{mn}}$, where $\operatorname{sech}(v)^{1/2}$ is the characteristic function of the Meixner law and $\mathcal{A}$, whose Fourier transform is formally $|\zeta(\frac12+i\tau)|^2$, is the autocorrelation of the measure $\sum_nn^{-1/2}\delta_{\log n}$. The measure $\mathcal{A}$ is positive and positive definite, and the series for $\mathcal{C}$ diverges: truncated at $m,n\le X$ it equals $X\int\operatorname{sech}(v-r)^{1/2}e^{-|r|/2}\,dr+O(X^{1/2})$. This growth comes from the pole of $\zeta$, and the factor $(\frac{d^2}{dv^2}-\frac14)^2$ does not remove it, since $(\frac{d^2}{dv^2}-\frac14)e^{-|v|/2}=-\delta_0$.

Supersedes: Proposition 28.22 of Polson (2026) states the factorisation formally with the constant $\frac14\pi^{-1/2}2^{3/2}\pi^{3/2}$, and says that the convolution $\operatorname{sech}(\cdot)^{1/2}*\mathcal{A}$ diverges logarithmically and that the factor $\bigl(\frac{d^2}{dv^2}-\frac14\bigr)^2$ removes the divergence.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:prop:autocorr`.

<a id="floor_reciprocal_xi_infinitely_divisible_ggc_iff_rh"></a>

### Splitting of the canonical function of the floor of the reciprocal Xi function, and GGC if and only if RH

Theorem 6.2.15, page 185 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi(z)=\frac12z(z-1)\pi^{-z/2}\Gamma(z/2)\zeta(z)$ be Riemann's xi function, $c=\zeta'(\frac12)/\zeta(\frac12)=-\theta'(0)\approx2.68609$ with $\theta$ the Riemann Siegel theta function, and $W(x)=\frac12\sum_\rho e^{x(\rho-\frac12)^2}$ for $x>0$, the sum over the nontrivial zeros $\rho$ of $\zeta$ with multiplicity. Unconditionally, $e^{-c\sqrt s}\,\xi(\frac12)/\xi(\frac12+\sqrt s)$ is the Laplace transform of an infinitely divisible law on $(0,\infty)$ with Lévy density $\frac{W(x)+c/(2\sqrt{\pi x})}{x}$, and $W(x)+c/(2\sqrt{\pi x})>0$ for all $x>0$. The following are equivalent: (i) the Riemann Hypothesis; (ii) this law is a generalized gamma convolution (GGC); (iii) $W(x)+c/(2\sqrt{\pi x})$ is completely monotone; (iv) $W$ is completely monotone. The proof uses the positivity of $w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}e^{-(\log n)^2/4x}$ for all $x>0$ ($\Lambda$ the von Mangoldt function), which rests on the verification of RH to height $50$, an explicit bound for the zero counting function (Hasanalizade, Shen and Wong 2022) and three interval enclosures. That the floor is infinitely divisible also follows from Theorem A of Polson (2026), which gives $W>0$ and $W'<0$.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:thm:floor`.

<a id="floor_canonical_function_derivative_signs_to_order_1e10"></a>

### Derivative signs of the canonical function of the floor to order ten billion, a separate proof (Williamson measures)

Theorem 6.2.17, page 187 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $k(y)=W(y)+c/(2\sqrt{\pi y})$ with $W(y)=\frac12\sum_\rho e^{y(\rho-\frac12)^2}$ (the sum over the nontrivial zeros $\rho$ of the Riemann zeta function $\zeta$ with multiplicity) and $c=\zeta'(\frac12)/\zeta(\frac12)$. For every integer $n$ with $0\le n\le10^{10}$ and every $y>0$, $(-1)^n\,k^{(n)}(y)>0$. Consequently, by Williamson's theorem, for each $1\le n\le10^{10}$ there is a positive measure $\Omega_n$ on $(0,\infty)$ with $k(y)=\int(1-ty)_+^{\,n-1}\,\Omega_n(dt)$. For $n\ge1$ the proof uses the verification of RH to height $T_0=3\times10^{12}$ (Platt and Trudgian 2021), an explicit bound for the zero counting function (Hasanalizade, Shen and Wong 2022) and one interval enclosure, and no other information about the zeros. This is a separate proof of a smaller range: the same signs hold for every $n\le10^{17}$ by Theorem 7.10 of Polson's arXiv paper 1708.02653 (version 24, 2026), which proves them for $W$, and the book proves them for every $n\le1.16\times10^{25}$ by sharper estimates of the heat trace.

Source: [`ch/rz2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz2.tex), label `rz:thm:nmono`.

## Chapter 7. Infinite divisibility and truncations of the reciprocal zeta function

<a id="reciprocal_zeta_infinitely_divisible_unconditional"></a>

### The normalised reciprocal zeta function is infinitely divisible

Theorem 7.1.1, page 192 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\zeta$ be the Riemann zeta function, $\sqrt s$ the principal square root on $\mathbb{C}\setminus(-\infty,0]$, $\Lambda$ the von Mangoldt function and \[ \varphi(s)=\frac{\zeta(\frac12)}{(1-4s)\,\zeta(\frac12+\sqrt s)},\qquad w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,e^{-(\log n)^2/4x}\quad(x>0), \] so that $\varphi$ is the normalised reciprocal zeta function (the factor $1-4s$ removes the pole of $\zeta(\frac12+\sqrt s)$ at $s=\frac14$) and $w$ is its canonical function: $-\varphi'(s)/\varphi(s)=\int_0^\infty e^{-sx}w(x)\,dx$ for $s>\frac14$. Then $w(x)>0$ for every $x>0$. Consequently, unconditionally (with no hypothesis on the zeros of $\zeta$), $\varphi$ is the Laplace transform of an infinitely divisible probability law on $(0,\infty)$ with zero drift and Lévy density $w(x)/x$: \[ \frac{\zeta(\frac12)}{(1-4s)\,\zeta(\frac12+\sqrt s)}=\exp\Bigl(-\int_0^\infty\bigl(1-e^{-sx}\bigr)\frac{w(x)}{x}\,dx\Bigr),\qquad \operatorname{Re}s>0 . \] Equivalently, the positivity of $w$ is the Gaussian smoothed Chebyshev bound $\sum_{n\ge2}\Lambda(n)\,n^{-1/2}e^{-(\log n)^2/4x}<\sqrt{4\pi x}\,e^{x/4}$ for all $x>0$. The proof uses the verification of the Riemann Hypothesis up to height $50$, an explicit Riemann von Mangoldt bound for the zero counting function, and three interval enclosures computed in outward rounded interval arithmetic.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:thm:id`.

<a id="gaussian_smoothed_prime_number_theorem_error_negative_bias"></a>

### Negative bias of the Gaussian smoothed prime number theorem error

Proposition 7.1.2, page 193 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Lambda$ be the von Mangoldt function, $\psi_C(y)=\sum_{n\le y}\Lambda(n)$ Chebyshev's function, $\Delta(y)=\psi_C(y)-y$ and $\operatorname{erfc}$ the complementary error function. For $x>0$ put \[ w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,e^{-(\log n)^2/4x},\qquad I(x)=\tfrac12e^{x/4}\operatorname{erfc}\bigl(\tfrac{\sqrt x}2\bigr), \] \[ Q_\Delta(x)=\frac{1}{\sqrt{4\pi x}}\int_{1^-}^\infty y^{-1/2}e^{-(\log y)^2/4x}\,d\Delta(y), \] where $w$ is the canonical function of the normalised reciprocal zeta function $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\,\zeta(\frac12+\sqrt s)\bigr)$ (so $-\varphi'/\varphi$ is the Laplace transform of $w$) and $Q_\Delta$, a Stieltjes integral, is the Gaussian smoothed error of the prime number theorem. Then $w(x)=I(x)-Q_\Delta(x)$ for $x>0$, and as $x\to\infty$, $I(x)\sim(\pi x)^{-1/2}$ and $w(x)\sim\frac c2(\pi x)^{-1/2}$, where $c=\frac12\bigl(\log\pi-\psi(\frac14)\bigr)\approx2.68609$ and $\psi$ is the digamma function. Hence $Q_\Delta(x)\sim-(\frac c2-1)(\pi x)^{-1/2}$ with $\frac c2-1\approx0.343$, a definite negative bias. On the other hand, any bound $|\Delta(y)|\le A\sqrt y\log^2y$ for $y\ge2$, including the one that the Riemann Hypothesis implies, yields only $|Q_\Delta(x)|\le A\,\bigl(x/2+O(\sqrt x)\bigr)$ as $x\to\infty$ after integration by parts.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:bias`.

<a id="reciprocal_zeta_bondesson_hierarchy_classification"></a>

### The reciprocal zeta function in the Bondesson hierarchy: not self decomposable, not a GGC

Theorem 7.1.3, page 195 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\,\zeta(\frac12+\sqrt s)\bigr)$ be the normalised reciprocal zeta function ($\zeta$ the Riemann zeta function, $\sqrt s$ the principal square root). For an infinitely divisible law (ID) on $[0,\infty)$ with Lévy density $k(x)/x$, the class SD of self decomposable laws is given by $k$ nonincreasing, the Bondesson class BO by the Lévy density $k(x)/x$ completely monotone, and the class GGC of generalised gamma convolutions by $k$ completely monotone; these classes are nested as $\mathrm{GGC}\subset\mathrm{SD}\subset\mathrm{ID}$ and $\mathrm{GGC}\subset\mathrm{BO}\subset\mathrm{ID}$. Then, with no hypothesis on the zeros of $\zeta$, the law with Laplace transform $\varphi$ is in ID but not in SD and not in BO, and it is not a mixture of exponentials (its density is not completely monotone, equivalently $\varphi$ is not a Stieltjes function). In particular it is not a generalised gamma convolution. The ID clause is the infinite divisibility theorem for $\varphi$; the SD and GGC clauses follow from the canonical function of $\varphi$ having positive derivative $\frac14$ at the origin; the BO clause combines the fact that membership in BO would force the Riemann Hypothesis with Selberg's central limit theorem for the argument function $S$ of $\zeta$, applied under that hypothesis.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:thm:class`.

<a id="reciprocal_zeta_density_rises_from_origin"></a>

### Density of the reciprocal zeta law at the origin

Theorem 7.1.4, page 197 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\,\zeta(\frac12+\sqrt s)\bigr)$ be the normalised reciprocal zeta function and $X_\varphi$ a random variable with $\mathbb{E}e^{-sX_\varphi}=\varphi(s)$ (it is infinitely divisible, by the positivity of the canonical function $w$ of $\varphi$). Here $w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\Lambda(n)\,n^{-1/2}e^{-(\log n)^2/4x}$ with $\Lambda$ the von Mangoldt function, and $w(0^+)=1$, $w'(0^+)=\frac14$. The law of $X_\varphi$ is absolutely continuous, with a density $f$ that can be taken continuous on $(0,\infty)$ and satisfies Steutel's equation $x f(x)=\int_0^x f(x-y)\,w(y)\,dy$. Then $f(0^+)=-\zeta(\frac12)/4\approx0.36509$ and $f'(0^+)=f(0^+)\,w'(0^+)=-\zeta(\frac12)/16\approx0.09127>0$. In particular $f$ is not monotone, so $X_\varphi$ is not a mixture of uniform laws on intervals $[0,a]$, and $f$ is not log convex and not completely monotone.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:thm:shape`.

<a id="reciprocal_zeta_law_not_built_from_ggc_operations"></a>

### The reciprocal zeta law is not a sum, product or limit of GGC variables

Corollary 7.1.6, page 199 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $X_\varphi$ be a random variable whose Laplace transform is the normalised reciprocal zeta function $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\,\zeta(\frac12+\sqrt s)\bigr)$. Then, with no hypothesis on the zeros of $\zeta$, $X_\varphi$ is not a sum or a product of independent generalised gamma convolution (GGC) variables, nor a weak limit of such. (The class of GGC laws is closed under independent sums and weak limits, and by a theorem of Bondesson under independent products, so any variable assembled from GGC variables by these operations is a GGC, while $X_\varphi$ is not a GGC.)

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:cor:hm`.

<a id="ggc_class_not_closed_under_harmonic_mean"></a>

### The GGC class is not closed under the harmonic mean

Proposition 7.1.7, page 199 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. The statement (H) says: if $V_1$ and $V_2$ are independent generalised gamma convolution (GGC) variables, then so is the harmonic mean $2V_1V_2/(V_1+V_2)$. (H) is false. If $E_1$ and $E_2$ are independent standard exponential variables, then $H=2E_1E_2/(E_1+E_2)$ has density \[ f_H(h)=e^{-h}\int_0^\infty\Bigl(1+\frac h{2u}\Bigr)e^{-u-h^2/4u}\,du,\qquad h>0, \] and $H$ is not a GGC.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:hmfalse`.

<a id="vere_jones_weights_are_necklace_polynomials"></a>

### Vere-Jones weights of the zeta function are necklace polynomials

Theorem 7.1.8, page 201 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For $k\ge1$ let $M_k(q)=\frac1k\sum_{d\mid k}\mu(k/d)\,q^d$ be the necklace polynomial ($\mu$ the Möbius function) and, for $n\ge1$, $b_n(q)=\sum_{i\ge0,\ 2^i\mid n}M_{n/2^i}(q)$. Then for $|y|<1$ and $0<q\le1$ \[ \frac1{1-qy}=\prod_{n\ge1}\bigl(1+y^n\bigr)^{b_n(q)}, \] so that $\zeta(\frac12+u)=\prod_p\prod_{n\ge1}\bigl(1+p^{-nu}\bigr)^{b_n(p^{-1/2})}$ for $\operatorname{Re}u>\frac12$, the product over the primes $p$. Equivalently, the weights $B_N$ of the logarithmic Vere-Jones representation $\log\zeta(\frac12+u)=\sum_{N\ge2}B_N\log\bigl(1+N^{-u}\bigr)$ for $\operatorname{Re}u>\frac12$, which are defined by $B_N=\sum_{n^m=N}d_m\,\Lambda(n)/(\sqrt n\,\log n)$ with $\Lambda$ the von Mangoldt function, $\sum_{m\ge1}d_mm^{-v}=1/\eta(1+v)$ and $\eta(s)=(1-2^{1-s})\zeta(s)$, satisfy $B_N=0$ unless $N=p^n$ is a prime power, and then $B_{p^n}=b_n(p^{-1/2})$. At $q=1$, $b_n(1)=1$ if $n$ is a power of $2$ and $b_n(1)=0$ otherwise, and the identity becomes Euler's binary identity $\frac1{1-y}=\prod_{j\ge0}\bigl(1+y^{2^j}\bigr)$.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:thm:necklace`.

<a id="necklace_weights_sign_rule"></a>

### Sign rule for the necklace weights of the zeta function

Proposition 7.1.9, page 201 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For a prime $p$ put $q=p^{-1/2}$, and let $b_n(q)=\sum_{i\ge0,\ 2^i\mid n}M_{n/2^i}(q)$, where $M_k(q)=\frac1k\sum_{d\mid k}\mu(k/d)\,q^d$ is the necklace polynomial and $\mu$ the Möbius function (these numbers $b_n(p^{-1/2})$ are the Vere-Jones weights $B_{p^n}$ of the zeta function at prime powers). Then $b_{2^a}(q)>0$ for all $a\ge1$, and for every odd $n>1$, $b_n(q)\ne0$ and the sign of $b_n(q)$ is $(-1)^{\omega(n)}$, where $\omega(n)$ is the number of distinct prime factors of $n$. The book proves that $b_n(q)\ne0$ for every odd $n>1$; the earlier statement of the sign rule, in Polson's 2026 book On Hilbert's 8th Problem, is for the odd $n>1$ with $b_n(q)\ne0$ only.

Supersedes: Proposition 29.2 of Polson (2026) states the sign rule only for odd $n>1$ with $b_n(q)\ne0$.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:signs`.

<a id="levy_representation_reciprocal_zeta_across_critical_strip"></a>

### Lévy representation of the reciprocal zeta function across the critical strip

Theorem 7.1.12, page 204 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Fix $\sigma\in(0,1)$, let $\zeta$ be the Riemann zeta function, $\sqrt s$ the principal square root, and put \[ \varphi_\sigma(s)=\frac{\zeta(\sigma)}{\bigl(1-s/(1-\sigma)^2\bigr)\,\zeta(\sigma+\sqrt s)}, \] so that $\varphi_{1/2}$ is the normalised reciprocal zeta function. For $s>(1-\sigma)^2$, $-\varphi_\sigma'(s)/\varphi_\sigma(s)=\int_0^\infty e^{-sx}w_\sigma(x)\,dx$, where $w_\sigma$ has the three equivalent forms \begin{align*} w_\sigma(x)&=e^{(1-\sigma)^2x}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\frac{\Lambda(n)}{n^{\sigma}}\,e^{-(\log n)^2/4x}\\ &=e^{(1-\sigma)^2x}+\sum_p\sum_{n\ge1}b_n(p^{-\sigma})\,h_{n\log p}(x)\\ &=\tfrac12e^{(1-\sigma)^2x}\operatorname{erfc}\bigl((1-\sigma)\sqrt x\bigr)+\sum_{\gamma>0}\operatorname{Re}\Bigl[e^{a_\gamma^2x}\operatorname{erfc}(a_\gamma\sqrt x)\Bigr]\\ &\qquad+\frac{\gamma_E+\log\pi}{4\sqrt{\pi x}}+\sum_{k\ge0}\Bigl[\tfrac12e^{a_k^2x}\operatorname{erfc}(a_k\sqrt x)-\frac1{4(k+1)\sqrt{\pi x}}\Bigr]. \end{align*} Here $\Lambda$ is the von Mangoldt function, $\gamma_E$ is Euler's constant, $\operatorname{erfc}$ is the complementary error function, $b_n(q)=\sum_{i\ge0,\ 2^i\mid n}M_{n/2^i}(q)$ with $M_k(q)=\frac1k\sum_{d\mid k}\mu(k/d)q^d$ the necklace polynomials, $h_\lambda(x)=\sum_{k\ge0}e^{-x((2k+1)\pi/\lambda)^2}-\frac{\lambda}{4\sqrt{\pi x}}=\frac{\lambda}{\sqrt{4\pi x}}\sum_{j\ge1}(-1)^je^{-j^2\lambda^2/4x}$ is the theta function of a block, $a_\gamma=\sigma-\frac12-i\gamma$ and $a_k=2k+2+\sigma$. The last form is written under the Riemann Hypothesis, with $\gamma>0$ the ordinates of the zeros $\frac12+i\gamma$ of $\zeta$; in general $\gamma$ runs over $-i(\rho-\frac12)$ for the zeros $\rho$. As $x\to\infty$, $w_\sigma(x)\sim\dfrac{\zeta'(\sigma)/\zeta(\sigma)}{2\sqrt{\pi x}}$, with $\zeta'(\sigma)/\zeta(\sigma)>0$ on $(0,1)$. If $w_\sigma\ge0$, then $\varphi_\sigma$ is the Laplace transform of an infinitely divisible law on $(0,\infty)$ with Lévy density $w_\sigma(x)/x$: \[ \log\varphi_\sigma(s)=-\int_0^\infty\bigl(1-e^{-sx}\bigr)\frac{w_\sigma(x)}x\,dx,\qquad\operatorname{Re}s>0 . \]

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:thm:levystrip`.

<a id="shifted_reciprocal_zeta_never_ggc_on_the_strip"></a>

### The shifted reciprocal zeta function is never a GGC on a line of the strip

Proposition 7.1.13, page 206 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For every $\sigma\in(0,1)$ let $\varphi_\sigma(s)=\zeta(\sigma)/\bigl((1-s/(1-\sigma)^2)\,\zeta(\sigma+\sqrt s)\bigr)$ ($\zeta$ the Riemann zeta function, $\sqrt s$ the principal square root). Then $\varphi_\sigma$ is not the Laplace transform of a generalised gamma convolution, and if $\varphi_\sigma$ is the Laplace transform of an infinitely divisible law, that law is not self decomposable. No hypothesis on the zeros of $\zeta$ and no computation is used: the canonical function $w_\sigma(x)=e^{(1-\sigma)^2x}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\Lambda(n)\,n^{-\sigma}e^{-(\log n)^2/4x}$ of $\varphi_\sigma$, with $-\varphi_\sigma'/\varphi_\sigma(s)=\int_0^\infty e^{-sx}w_\sigma(x)\,dx$ for $s>(1-\sigma)^2$ and $\Lambda$ the von Mangoldt function, satisfies $w_\sigma(x)=1+(1-\sigma)^2x+o(x)$ as $x\downarrow0$ and so is strictly increasing near the origin. The book proves this for every $\sigma\in(0,1)$ and includes the clause on self decomposability; the earlier statement, in Polson's 2026 book On Hilbert's 8th Problem, proves that $\varphi_\sigma$ is not a generalised gamma convolution for $\sigma\le\frac12$ only (for $\frac12<\sigma<1$ it is verified numerically) and has no such clause.

Supersedes: Proposition 29.5 of Polson (2026) proves that the shifted reciprocal zeta function is not a generalised gamma convolution for $\sigma\le\frac12$ only (for $\frac12<\sigma<1$ it is verified numerically and not proved), and has no clause on self decomposability.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:onlyhalf`.

<a id="laguerre_polya_class_integer_atom_ggc"></a>

### Laguerre Pólya class as integer atom GGCs

Proposition 7.1.15, page 209 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $G$ be a real even entire function of order less than $2$ with $G(0)\ne0$. Then $G$ is in the Laguerre Pólya class (the locally uniform limits of real polynomials with only real zeros, which for such $G$ means that all zeros of $G$ are real) if and only if $s\mapsto G(0)/G(i\sqrt s)$ is the Laplace transform of a generalised gamma convolution whose Thorin measure is $\sum_km_k\delta_{\alpha_k^2}$ with integers $m_k\ge1$, where $\pm\alpha_k$ are the zeros of $G$ with multiplicities $m_k$.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:lpthorin`.

<a id="lee_yang_property_integer_atom_ggc"></a>

### Lee and Yang single spin laws as integer atom GGCs

Proposition 7.1.16, page 211 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\rho$ be a symmetric probability law on $\mathbb{R}$ with $\int e^{bt^2}\rho(dt)<\infty$ for every $b$. Say that $\rho$ has the Lee and Yang property if ferromagnetic systems of $\rho$ spins have partition function zeros only on the imaginary field axis; by a theorem of Newman this holds exactly when $z\mapsto\int e^{izt}\rho(dt)$ is in the Laguerre Pólya class. Then $\rho$ has the Lee and Yang property if and only if $s\mapsto1/\int e^{\sqrt s\,t}\rho(dt)$ is the Laplace transform of a generalised gamma convolution with zero drift whose Thorin measure has integer atoms. The atoms are at the squares $\alpha_k^2$ of the zeros $\pm\alpha_k$ of $\int e^{izt}\rho(dt)$, with the multiplicities of the zeros.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:newmanthorin`.

<a id="random_harmonic_series_matches_polya_kernel_strip"></a>

### Random harmonic series and the strip of the Pólya kernel

Proposition 7.1.17, page 212 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Phi(t)=2\sum_{n\ge1}\bigl(2\pi^2n^4e^{9t/2}-3\pi n^2e^{5t/2}\bigr)\exp(-\pi n^2e^{2t})$ be the Pólya kernel, an even function whose Fourier transform $\int_{\mathbb{R}}\Phi(t)\cos(\tau t)\,dt$ is $\xi(\frac12+i\tau)$, with $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$; it is analytic in the strip $|\operatorname{Im}t|<\pi/4$ and $\log\Phi(t)=-\pi e^{2t}+\frac92t+O(1)$ as $t\to+\infty$. For $\kappa>-1$ let $\varepsilon_k$ be independent fair signs $\pm1$ and $H_\kappa=\sum_{k\ge1}\varepsilon_k/(2(k+\kappa))$, with characteristic function $\hat f_\kappa(v)=\prod_{k\ge1}\cos\bigl(v/(2(k+\kappa))\bigr)$ and density $f_\kappa$. (i) As $s\to\infty$, $\log\mathbb{E}e^{sH_\kappa}=\frac s2\log s+A(\kappa)s+o(s)$ with $A(\kappa)=A(0)-\frac12\bigl(\psi(1+\kappa)+\gamma_E\bigr)$, where $\psi$ is the digamma function, $\gamma_E$ is Euler's constant and $A(0)=-0.1485757\ldots$. Hence $\log\Pr(H_\kappa>u)=-\frac12e^{2u-2A(\kappa)-1}(1+o(1))$ as $u\to\infty$, and the double exponential terms of $H_\kappa$ and of $\Phi$ agree exactly when $A(\kappa)=-\frac12(1+\log2\pi)$, that is at $\kappa^*=6.61845\ldots$. (ii) $\limsup_{|v|\to\infty}|v|^{-1}\log|\hat f_\kappa(v)|=-\pi/4$ for every $\kappa>-1$. Hence $f_\kappa$ is analytic in the strip $|\operatorname{Im}t|<\pi/4$, which is exactly the strip of analyticity of $\Phi$.

Source: [`ch/rz3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz3.tex), label `rz:prop:harmonic`.

<a id="truncated_lattice_zeros_first_order_expansion"></a>

### First order expansion of the zeros of the truncated lattice transforms

Proposition 7.2.2, page 220 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $e_k,e_k'$ ($k\ge1$) be independent standard exponential random variables, $\mathcal S_n=\sum_{k\le n}k^{-2}(e_k+e_k')$ and $F_n(w)=\mathbb{E}\,\mathcal S_n^{\,w}$ for $\operatorname{Re}w>-2$. As $n\to\infty$, $F_n\to F$ locally uniformly on $\operatorname{Re}w>-2$, where $F(w)=\mathbb{E}\bigl[(\sum_{k\ge1}k^{-2}(e_k+e_k'))^w\bigr]=2\pi^w\xi(2w)$ and $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ is the completed zeta function, so the zeros of $F$ are the points $\rho/2$ for the nontrivial zeros $\rho$ of $\zeta$. Let $\rho$ be a simple zero of $\xi$ and put \[ \kappa(\rho)=\frac{\rho\,\xi(\rho-2)}{2\pi\,\xi'(\rho)}=\frac{\rho\,\xi(3-\rho)}{2\pi\,\xi'(\rho)} . \] There is a constant $C$ (depending on $\rho$) such that for all large $n$, $F_n$ has exactly one zero $w_n(\rho)$ in the disc $\bigl|w-\frac\rho2-\frac{\kappa(\rho)}n\bigr|<C/n^2$; thus $w_n(\rho)=\frac\rho2+\frac{\kappa(\rho)}n+O(n^{-2})$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prop:firstorder`.

<a id="truncated_lattice_displacement_formula_hardy_function"></a>

### Displacement constant of the truncated zeros in terms of Hardy's function

Theorem 7.2.3, page 221 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\rho=\frac12+i\gamma$ be a simple zero of $\zeta$ on the critical line, $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, and $\kappa(\rho)=\rho\,\xi(\rho-2)/(2\pi\,\xi'(\rho))$ the constant of the first order displacement of the zeros of the truncated lattice transforms $F_n(w)=\mathbb{E}\bigl[(\sum_{k\le n}k^{-2}(e_k+e_k'))^w\bigr]$, with $e_k,e_k'$ independent standard exponentials. Let $\theta(t)=\arg\Gamma(\frac14+\frac{it}2)-\frac t2\log\pi$ be the Riemann Siegel theta function, $Z(t)=e^{i\theta(t)}\zeta(\frac12+it)$ Hardy's function, and $a(\gamma)=e^{i\theta(\gamma)}\zeta(\frac52+i\gamma)$. Then \begin{align*} \operatorname{Re}\kappa(\rho)&=\frac{(\gamma^2-\frac{15}4)\operatorname{Im}a(\gamma)-4\gamma\operatorname{Re}a(\gamma)}{4\pi^2Z'(\gamma)},\\ \operatorname{Im}\kappa(\rho)&=\frac{4\gamma\operatorname{Im}a(\gamma)+(\gamma^2-\frac{15}4)\operatorname{Re}a(\gamma)}{4\pi^2Z'(\gamma)} . \end{align*}

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:kappaexact`.

<a id="conjecture_t_false_truncated_zeros_right_of_the_line"></a>

### Conjecture T is false for large n

Corollary 7.2.4, page 222 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: counterexample, proved.

Statement. Let $e_k,e_k'$ be independent standard exponential variables and $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ for $\operatorname{Re}w>-2$. Conjecture T says that for every $n\ge1$, $F_n$ has no zeros with $\operatorname{Re}w>\frac14$ (which, if it held for infinitely many $n$, would imply the Riemann Hypothesis by Hurwitz's theorem, since $F_n\to2\pi^w\xi(2w)$). Conjecture T is false: for all sufficiently large $n$, $F_n$ has a zero with $\operatorname{Re}w>\frac14$ near $\rho_{213}/2$, where $\rho_{213}=\frac12+i\gamma_{213}$ and $\gamma_{213}=415.4552\ldots$ is the ordinate of the 213th zero of the zeta function in the upper half plane. The proof uses that the first $1.5\times10^9$ zeros lie on the critical line and are simple (van de Lune, te Riele and Winter), the first order displacement $\kappa(\rho)=\rho\,\xi(\rho-2)/(2\pi\,\xi'(\rho))$ of the zero of $F_n$ near $\rho/2$, and a ball arithmetic enclosure $233.95<\operatorname{Re}\kappa(\rho_{213})<233.96$. The statement is for all sufficiently large $n$ and gives no explicit $n$; the book notes that the first order regime requires $n$ large compared with $|\kappa(\rho_{213})|\approx2311$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:cor:Tfalse`.

<a id="symmetrised_truncations_symmetry_and_limit"></a>

### Symmetrised truncations: symmetry, limit and persistence of simple zeros on the line

Proposition 7.2.6, page 223 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $e_k,e_k'$ be independent standard exponential variables, $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ for $\operatorname{Re}w>-2$, and define the symmetrised truncation and its restriction to the line by \[ K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\tfrac12-w),\qquad Z_n(y)=K_n(\tfrac14+iy). \] Then $K_n(\frac12-w)=K_n(w)$ and $K_n(\bar w)=\overline{K_n(w)}$, so $Z_n$ is real for real $y$ and the zeros of $K_n$ off the line $\operatorname{Re}w=\frac14$ come in quadruples. As $n\to\infty$, $K_n\to4\,\xi(2w)$ locally uniformly on $-\frac74<\operatorname{Re}w<\frac94$, where $\xi$ is the completed zeta function. Every simple zero $\frac14+i\frac\gamma2$ of $\xi(2w)$ on the line $\operatorname{Re}w=\frac14$ is, for all large $n$, the limit of a zero of $K_n$ lying exactly on that line.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prop:symtrunc`.

<a id="symmetrised_truncations_real_rooted_implies_rh"></a>

### Real rooted symmetrised truncations for infinitely many n imply the Riemann Hypothesis

Theorem 7.2.7, page 224 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $e_k,e_k'$ be independent standard exponential variables, $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ for $\operatorname{Re}w>-2$, and $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\tfrac12-w)$. If for infinitely many $n$ all zeros of $K_n$ in the strip $0<\operatorname{Re}w<\frac12$ lie on the line $\operatorname{Re}w=\frac14$, then the Riemann Hypothesis holds.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:symhurwitz`.

<a id="symmetrised_truncations_single_term_bound_explicit_heights"></a>

### Single term bound: no off line zeros of the symmetrised truncation above explicit heights

Theorem 7.2.8, page 225 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For $n\ge1$ let $\alpha_k,\beta_k$ ($k\le n$) be defined by the partial fraction expansion $\prod_{j\le n}(1+\lambda/j^2)^{-2}=\sum_{k\le n}k^{-2}\bigl\{\alpha_k(1+\lambda/k^2)^{-1}+\beta_k(1+\lambda/k^2)^{-2}\bigr\}$ (so $\beta_k=k^2\prod_{j\le n,\,j\ne k}(1-k^2/j^2)^{-2}$ and $\alpha_k=-2\beta_k\sum_{j\le n,\,j\ne k}k^2/(j^2-k^2)$), and put $Q_n(u)=\sum_{k\le n}k^{-2u}(\alpha_k+\beta_ku)$. Then $F_n(u-1)=\Gamma(u)Q_n(u)$ for $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ with $e_k,e_k'$ independent standard exponentials, and the symmetrised truncation is $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\frac12-w)$. For $y>|\alpha_1|/|\beta_1|$ let \[ \delta_n(y)=\sum_{k=2}^nk^{-2}\frac{|\beta_k|y+|\alpha_k|}{|\beta_1|y-|\alpha_1|},\qquad \delta_n'(y)=\sum_{k=2}^nk^{-2}\frac{2\log k\,(|\beta_k|y+|\alpha_k|)+|\beta_k|}{|\beta_1|y-|\alpha_1|}+\delta_n(y)\,\frac{|\beta_1|}{|\beta_1|y-|\alpha_1|}, \] where $\delta_n'$ is a bound for a derivative used in the proof and is not the derivative of $\delta_n$. If $\delta_n(y)<1$ and \[ \log\frac y\pi-\frac1y-\frac{|\beta_1|}{|\beta_1|y-|\alpha_1|}-\frac{\delta_n'(y)}{1-\delta_n(y)}>0, \] then $K_n$ has no zeros off the line $\operatorname{Re}w=\frac14$ in $0<\operatorname{Re}w<\frac12$ at any height (imaginary part) at least $y$. The condition holds for all $y\ge Y_n$, with $Y_1=4.81$, $Y_2=5.85$, $Y_3=7.11$, $Y_4=9.08$, $Y_5=12.8$, $Y_6=20.8$, $Y_8=104.1$ and $Y_{10}=21232$ (values obtained by evaluating the bounds with the exact coefficients).

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:Yn`.

<a id="conjecture_s_small_n_above_explicit_heights"></a>

### Conjecture S for n up to 6 and n = 8 above explicit heights

Corollary 7.2.9, page 226 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\frac12-w)$ be the symmetrised truncation, with $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ and $e_k,e_k'$ independent standard exponential variables. Conjecture S says that all zeros of $K_n$ in the strip $0<\operatorname{Re}w<\frac12$ lie on the line $\operatorname{Re}w=\frac14$. Proved: for $n=1,2,3,4,5,6,8$ the zeros of $K_n$ in the strip at height at least $Y_n$ lie on the line, with $Y_1=4.81$, $Y_2=5.85$, $Y_3=7.11$, $Y_4=9.08$, $Y_5=12.8$, $Y_6=20.8$, $Y_8=104.1$ (the heights of the single term bound). Below these heights the book verifies the same statement for these $n$ only numerically, by floating point argument principle counts and not by an enclosure (up to heights $6,7,8$ for $n=1,2,3$, up to $70$ for $4\le n\le6$ and up to $140$ for $n=8$); that part is not a certified result.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:cor:Ssmall`.

<a id="symmetrised_truncations_torus_bound"></a>

### Torus bound for the symmetrised truncations

Theorem 7.2.10, page 226 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\alpha_k,\beta_k$ ($k\le n$) be defined by $\prod_{j\le n}(1+\lambda/j^2)^{-2}=\sum_{k\le n}k^{-2}\bigl\{\alpha_k(1+\lambda/k^2)^{-1}+\beta_k(1+\lambda/k^2)^{-2}\bigr\}$, so that $F_n(u-1)=\Gamma(u)Q_n(u)$ with $Q_n(u)=\sum_{k\le n}k^{-2u}(\alpha_k+\beta_ku)=uP_n(u)+R_n(u)$, $P_n(u)=\sum_{k\le n}\beta_kk^{-2u}$, $R_n(u)=\sum_{k\le n}\alpha_kk^{-2u}$, for $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ with $e_k,e_k'$ independent standard exponentials, and let $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\frac12-w)$. For real $x$ and $\omega=(\omega_p)_{p\le n}$ in the torus $\mathbb{T}^{\pi(n)}$ of phases indexed by the primes $p\le n$, put $P_n(x,\omega)=\sum_{k\le n}\beta_kk^{-2x}\exp\bigl(-i\sum_{p\le n}v_p(k)\,\omega_p\bigr)$, where $v_p(k)$ is the exponent of $p$ in $k$, let $P_n'(x,\omega)=-\sum_{k\le n}2\log k\,\beta_kk^{-2x}\exp\bigl(-i\sum_{p\le n}v_p(k)\,\omega_p\bigr)$, $r(x)=\sum_{k\le n}|\alpha_k|k^{-2x}$ and $r'(x)=\sum_{k\le n}2\log k\,|\alpha_k|k^{-2x}$. Let $m_n=\min|P_n(x,\omega)|$, $M_n=\max\{-\operatorname{Re}(P_n'/P_n)(x,\omega)\}$ and $\ell_n=\max|P_n'/P_n|(x,\omega)$ over $\omega\in\mathbb{T}^{\pi(n)}$ and the band $1.2\le x\le1.3$, and \[ B_n(y)=\max_{0.3\le\sigma\le\frac12,\ \omega}\ \frac{\log\bigl(|P_n(\tfrac32-\sigma,\omega)|+r(\tfrac32-\sigma)/y\bigr)-\log\bigl(|P_n(1+\sigma,\omega)|-r(1+\sigma)/y\bigr)}{2\sigma-\frac12}. \] If $m_n>r(1.2)/y$, \[ \log\frac y\pi-\frac2y-M_n-\frac{\ell_n\,r(1.2)+r'(1.2)+r(1.2)/y}{y\,\bigl(m_n-r(1.2)/y\bigr)}>0\quad\text{and}\quad\log\frac y\pi-\frac1y>B_n(y), \] then $K_n$ has no zeros off the line $\operatorname{Re}w=\frac14$ in $0<\operatorname{Re}w<\frac12$ at height $y$, and the same holds at every larger height.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:torus`.

<a id="conjecture_s_above_explicit_thresholds_n_to_36"></a>

### Conjecture S for n up to 24, 32 and 36 above explicit thresholds

Corollary 7.2.11, page 227 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\frac12-w)$ be the symmetrised truncation, with $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ and $e_k,e_k'$ independent standard exponential variables. Conjecture S says that all zeros of $K_n$ in the strip $0<\operatorname{Re}w<\frac12$ lie on the line $\operatorname{Re}w=\frac14$. Proved, by torus bounds whose final boxes are enclosed in ball arithmetic: for $n=16,17,\dots,24,32,36$ the zeros of $K_n$ in the strip at height at least the threshold lie on the line, with thresholds $30,34,37,41,46,50,56,62,69$ for $n=16,\dots,24$ and $213$, $558$ for $n=32$, $36$; and the same holds for $n=7$ and $9\le n\le15$ above the thresholds $14,24,31,41,57,84,131,229$ for $n=7,9,10,11,12,13,14,15$. For all these $n$ the statement is also verified numerically below the thresholds, by floating point argument principle counts, which are not certified. Together with the case $n=1,\dots,6,8$ above explicit heights, this covers every $n\le24$ and $n=32,36$. The torus bounds behind the thresholds are proved in the book; in Polson's 2026 book On Hilbert's 8th Problem the corresponding corollary supposes that the extremes of the torus function, estimated by global optimisation and not enclosed, are bounds.

Supersedes: Corollary 30.11 of Polson (2026) is conditional: it supposes that the torus extremes, estimated by global optimisation and not enclosed, are bounds.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:cor:Storus`.

<a id="thorin_positivity_of_reciprocal_xi_is_weil_criterion"></a>

### Thorin positivity of R is Weil's criterion

Proposition 7.2.13, page 230 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $R(s)=\xi(\frac12)/\xi(\frac12+\sqrt s)$. Let $D(x)=\sum_\rho e^{(\rho-\frac12)x}$ be Weil's distribution, summed over the nontrivial zeros $\rho$ of $\zeta$, whose action on a test function $h$ is $\sum_\rho\hat h(\gamma_\rho)$ with $\gamma_\rho=-i(\rho-\frac12)$ and $\hat h$ the Fourier transform of $h$. Then $R$ is the Laplace transform of a generalised gamma convolution (whose Thorin measure is then $\sum_{\gamma>0}\delta_{\gamma^2}$ over the ordinates of the zeros) if and only if $D$ is positive definite, that is, $\sum_\rho\hat h(\gamma_\rho)\ge0$ for every $h=g*\tilde g$ with $g$ smooth and compactly supported and $\tilde g(x)=\overline{g(-x)}$. Both conditions are equivalent to the Riemann Hypothesis: the first by Proposition 3.4 of Polson's 2026 book On Hilbert's 8th Problem, the second by Weil's criterion (Weil 1952; Bombieri 2000). The proposition combines the two; the correspondence between Thorin measures and Weil's distribution behind it is formal.

Supersedes: Theorem 30.13 of Polson (2026) states the correspondence between Thorin positivity of $R$ and positive definiteness of Weil's distribution as a theorem.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:weil`.

<a id="weil_positivity_on_bernstein_cone"></a>

### Weil positivity on the Bernstein cone

Theorem 7.2.14, page 233 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. For $y>0$ let $W(y)=\frac12\sum_\rho e^{y(\rho-\frac12)^2}$, the sum over all nontrivial zeros $\rho$ of the Riemann zeta function (the heat trace of the zeros), and put $\gamma_\rho=-i(\rho-\frac12)$. Let $\hat h(r)=\int_0^\infty e^{-yr^2}g(y)\,dy$ with $g\ge0$, $g\not\equiv0$ and $\int_0^1g(y)\,y^{-1/2}\bigl(1+\log(1/y)\bigr)\,dy+\int_1^\infty g(y)\,dy<\infty$, so that $\hat h$ is a completely monotone function of $r^2$. Then \[ \sum_\rho\hat h(\gamma_\rho)=2\int_0^\infty W(y)\,g(y)\,dy>0, \] the sum being over all nontrivial zeros. In particular Weil's functional $\sum_\rho\hat h(\gamma_\rho)$ is positive for every test function of Thorin type $\hat h(r)=\prod_k(1+r^2/t_k)^{-1}$ (a nonempty product with $t_k>0$ and $\sum_k1/t_k<\infty$), for $\hat h(r)=\operatorname{sech}(cr)$ with $c>0$, and for $\hat h(r)=(1+r^2/a)^{-m}$ with $a>0$ and $m>\frac12$. No hypothesis on the zeros is made. The positivity of $W(y)$ for every $y>0$, which with the identity gives the result, is Theorem A of Polson's 2026 book On Hilbert's 8th Problem; the book gives a second proof of it, which uses the verification of the Riemann Hypothesis up to height $50$ and interval enclosures, and corrects the integrability condition and the examples of Theorem 30.14 of that book.

Supersedes: Theorem 30.14 of Polson (2026) has the integrability condition with $\log(1/y)$ in place of $1+\log(1/y)$, Thorin products without nonemptiness or $\sum_k1/t_k<\infty$, and $(1+r^2/a)^{-m}$ without $m>\frac12$, and its proof uses the verification of the Riemann Hypothesis to height $3\times10^{12}$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:bernsteincone`.

<a id="debruijn_newman_flow_ggc_path"></a>

### The de Bruijn Newman flow as a path of GGCs

Theorem 7.2.15, page 235 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}\bigl(2\pi^2n^4e^{9u}-3\pi n^2e^{5u}\bigr)\exp(-\pi n^2e^{4u})$, $H_t(z)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$ (so that $H_0(z)=\frac18\xi(\frac12+\frac{iz}2)$ for the completed zeta function $\xi$) and $R_t(s)=H_t(0)/H_t(2i\sqrt s)$, so that $R_0(s)=\xi(\frac12)/\xi(\frac12+\sqrt s)=R(s)$. Let $\Lambda$ be the de Bruijn Newman constant, so that $H_t$ has only real zeros exactly when $t\ge\Lambda$. Then for every $t\ge\Lambda$, $R_t$ is the Laplace transform of a generalised gamma convolution with Thorin measure $\sum_k\delta_{(z_k(t)/2)^2}$, where $z_k(t)>0$ are the positive zeros of $H_t$ listed with multiplicity. In particular $R_t$ is a GGC for every $t\ge0.2$, unconditionally (using the bound $\Lambda\le0.2$ of Platt and Trudgian, which improves the bound $0.22$ of the Polymath project). The Riemann Hypothesis holds if and only if $R_t$ is a GGC for arbitrarily small $t>0$ (and, since $\Lambda\ge0$ by Rodgers and Tao, if and only if $\Lambda=0$).

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:thm:dbn`.

<a id="no_ggc_multiplier_can_create_positivity"></a>

### No GGC multiplier can create the positivity of R

Proposition 7.2.16, page 236 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, $R(s)=\xi(\frac12)/\xi(\frac12+\sqrt s)$, $G(s)=\Gamma(\frac14+\frac{\sqrt s}2)\,\pi^{-\sqrt s/2}/\Gamma(\frac14)$ and $\varphi(s)=R(s)G(s)=\zeta(\frac12)/\bigl((1-4s)\zeta(\frac12+\sqrt s)\bigr)$, and let $c=\frac12(\log\pi-\psi(\frac14))\approx2.68609$ with $\psi$ the digamma function. Let $V$ be the Laplace transform of any generalised gamma convolution (GGC). Then $V\cdot R$ is the Laplace transform of a GGC if and only if the Riemann Hypothesis holds. For $\varphi$ only one direction holds for every $V$: if $V\cdot\varphi$ is the Laplace transform of a GGC, then the Riemann Hypothesis holds. The converse depends on $V$: it fails for $V\equiv1$ (since $\varphi$ is not a GGC, unconditionally), and holds for $V(s)=e^{-c\sqrt s}/G(s)$, which is the Laplace transform of a GGC and for which $V\varphi=e^{-c\sqrt s}R$ is the floor of the reciprocal Xi function. In Polson's 2026 book On Hilbert's 8th Problem the same is claimed with $\varphi$ in place of $R$, that is, both directions for every multiplier, which fails for $V\equiv1$; the book keeps one direction.

Supersedes: Proposition 30.16 of Polson (2026) states that the same is true with $\varphi$ in place of $R$, that is, both directions for every generalised gamma convolution multiplier, which fails for the multiplier $1$ by Theorem 28.9 of that book.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prop:nomult`.

<a id="open_conjecture_s_for_infinitely_many_n"></a>

### Conjecture S for infinitely many n

Problem 7.1, page 240 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $e_k,e_k'$ be independent standard exponential variables, $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$ for $\operatorname{Re}w>-2$, and $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-\frac12}F_n(\frac12-w)$ the symmetrised truncation. Conjecture S: for every $n\ge1$ all zeros of $K_n$ in the strip $0<\operatorname{Re}w<\frac12$ lie on the line $\operatorname{Re}w=\frac14$. Problem: prove that this holds for infinitely many $n$; by a Hurwitz argument this would imply the Riemann Hypothesis. The book proves it above explicit heights for $n\le24$ and $n=32,36$ and checks it numerically below those heights. Since the one sided Conjecture T (no zero of $F_n$ with $\operatorname{Re}w>\frac14$) and the Hermite Biehler inequality for $\pi^{-w}F_n(w)$ fail for large $n$, a proof must establish real rootedness without a one sided inequality, for instance by an interlacing argument in $n$, or by following the zeros of $K_n(\frac14+iy)$ continuously in a parameter interpolating between $n$ and $n+1$ and excluding collisions of zeros on the line.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:S`.

<a id="open_torus_method_beyond_n_36"></a>

### Torus method beyond n = 36

Problem 7.2, page 240 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $\alpha_k,\beta_k$ ($k\le n$) be the partial fraction coefficients of the Laplace transform $\prod_{j\le n}(1+\lambda/j^2)^{-2}$ of $\sum_{k\le n}k^{-2}(e_k+e_k')$, with $\beta_k=k^2\prod_{j\le n,\,j\ne k}(1-k^2/j^2)^{-2}$, and let $P_n(x,\omega)=\sum_{k\le n}\beta_kk^{-2x}\exp\bigl(-i\sum_{p\le n}v_p(k)\,\omega_p\bigr)$ be the torus function of the Dirichlet polynomial $P_n(u)=\sum_{k\le n}\beta_kk^{-2u}$ (phases $\omega_p$ indexed by the primes $p\le n$, $v_p(k)$ the exponent of $p$ in $k$). The torus method proves Conjecture S, that all zeros of the symmetrised truncation $K_n$ in $0<\operatorname{Re}w<\frac12$ lie on $\operatorname{Re}w=\frac14$, above explicit heights for $n=7$, $9\le n\le24$ and $n=32,36$ (the cases $n\le6$ and $n=8$ come from a single term bound); at $n=40$ the threshold, estimated by optimisation and not certified, is already near $4900$, and with the band $1.2\le x\le1.3$ the method ends between $n=40$ and $n=44$, where the torus minimum of $|P_n|$ in the band reaches zero. Problem: extend the torus bound with narrower bands around $x=\frac54$, and find the first $n$ at which $P_n$ has torus zeros with $x\ge\frac54$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:torus`.

<a id="open_finite_order_to_all_orders_n_monotonicity_of_the_floor"></a>

### From finite order to all orders for the floor

Problem 7.3, page 240 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $W(y)=\frac12\sum_\rho e^{y(\rho-\frac12)^2}$, the sum over the nontrivial zeros $\rho$ of the Riemann zeta function, let $c=\frac12(\log\pi-\psi(\frac14))\approx2.68609$, and let $k(y)=W(y)+c/(2\sqrt{\pi y})$ be the canonical function of the floor $e^{-c\sqrt s}\,\xi(\frac12)/\xi(\frac12+\sqrt s)$ of the reciprocal Xi function. The book proves, by a two range argument, that $(-1)^nk^{(n)}(y)>0$ for every $y>0$ and every $n\le10^{10}$, the bound on $n$ in that argument coming only from the height $3\times10^{12}$ to which the Riemann Hypothesis is verified, and the small $y$ branch of the proof holds for every $n$ on $(0,0.05/n]$. This is a separate proof of a smaller range: the same signs hold for every $n\le10^{17}$ by Theorem 7.10 of Polson's arXiv paper 1708.02653 (version 24, 2026), which proves them for $W$, and the book proves them for every $n\le1.16\times10^{25}$ by sharper estimates of the heat trace; the signs for every $n$ would make $k$ completely monotone, which is equivalent to the Riemann Hypothesis. Problem: can the two regimes be made to overlap for all $n$? This asks for a prime side estimate for $(-1)^nk^{(n)}(y)$ on a range $y\le y_n$ with $y_n$ growing at least linearly in $n$. In the proof the obstruction is the factor $2^n$ in the Cauchy estimate on the circle of radius $y/2$; a circle of radius close to $y$ loses only a factor of order $n$, but that argument also stops at a finite order.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:allorders`.

<a id="open_largest_cone_of_unconditional_weil_positivity"></a>

### Beyond the Bernstein cone

Problem 7.4, page 241 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Weil's functional $\sum_\rho\hat h(\gamma_\rho)$, with $\gamma_\rho=-i(\rho-\frac12)$ for the nontrivial zeros $\rho$ of the Riemann zeta function, is unconditionally positive for every test function $\hat h$ completely monotone in $r^2$ under an integrability condition (the Bernstein cone). Problem: find the largest cone of test functions, closed under the operations that preserve the Thorin structure, on which positivity of Weil's functional can be proved unconditionally, and determine how far it can localise in height.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:cone`.

<a id="open_deforming_the_binary_product"></a>

### Deforming the binary product

Problem 7.5, page 241 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. With the necklace weights $b_n(q)$ of the Vere-Jones factorisation, each Euler factor $(1-p^{-\sigma-u})^{-1}$ of $\zeta(\sigma+u)$ equals $\prod_{n\ge1}(1+p^{-nu})^{b_n(p^{-\sigma})}$, with integer multiplicities $0$ and $1$ at $\sigma=0$ (Euler's binary product) and necklace multiplicities $b_n(p^{-1/2})$ at $\sigma=\frac12$. The pseudo spectral trace formula writes the heat trace of the zeros of $\zeta$ (under the Riemann Hypothesis $\sum_{\gamma>0}e^{-x\gamma^2}$, and in general $\frac12\sum_\rho e^{x(\rho-\frac12)^2}$), net of its smooth part, as $e^{x/4}$ plus the sum over primes $p$ and $n\ge1$ of $b_n(p^{-1/2})$ times the heat trace $\sum_{k\ge0}e^{-x((2k+1)\pi/(n\log p))^2}$ of the arithmetic progression $(2k+1)\pi/(n\log p)$, net of its mean $n\log p/(4\sqrt{\pi x})$. Problem: follow this formula along the path $q=p^{-\sigma}$, $0\le\sigma\le\frac12$, and determine where and how the positivity of the interference pattern of the progressions can fail.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:binary`.

<a id="open_regularising_the_autocorrelation_of_the_polya_kernel"></a>

### Regularising the autocorrelation of the Pólya kernel

Problem 7.6, page 241 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $p(\tau)=|\Gamma(\frac14+\frac{i\tau}2)|^2/(2^{3/2}\pi^{3/2})$ be the Meixner density, $\Phi$ the Pólya kernel (an even function with Fourier transform $\Xi(\tau)=\xi(\frac12+i\tau)$) and $\mathcal C(v)=\int p(\tau)\,|\zeta(\frac12+i\tau)|^2e^{iv\tau}\,d\tau$. Then $(\Phi*\Phi)(v)=2^{-3/2}\bigl(\frac{d^2}{dv^2}-\frac14\bigr)^2\mathcal C(v)$. Formally $\mathcal C=\operatorname{sech}(\cdot)^{1/2}*\mathcal A$, where $\operatorname{sech}(v)^{1/2}$ is the characteristic function of the Meixner law and $\mathcal A=\sum_{m,n\ge1}(mn)^{-1/2}\delta_{\log(m/n)}$ is a positive, positive definite measure. The lattice series diverges: truncated at $m,n\le X$ it grows linearly in $X$, through the square of the pole term of the partial sums of $\zeta$, and the factor $(\frac{d^2}{dv^2}-\frac14)^2$ does not remove this growth. Problem: find the regularised form of $\operatorname{sech}(\cdot)^{1/2}*\mathcal A$ and decide whether it has a Thorin or Lee and Yang structure in $v$. In Polson's 2026 book On Hilbert's 8th Problem the problem is posed up to a logarithmic divergence removed by the pole; the book finds linear growth that the factor does not remove.

Supersedes: Open Problem 30.22 of Polson (2026) poses the problem up to a logarithmic divergence removed by the pole.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:meixner`.

<a id="open_unimodality_of_the_canonical_function_and_density"></a>

### Unimodality of w and of the density f

Problem 7.7, page 241 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. Corrects a statement of Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $w(x)=e^{x/4}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\Lambda(n)\,n^{-1/2}e^{-(\log n)^2/4x}$ be the canonical function of the normalised reciprocal zeta function $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\zeta(\frac12+\sqrt s)\bigr)$ and $f$ the density of the law with Laplace transform $\varphi$, which satisfies Steutel's equation $xf(x)=\int_0^xf(x-y)\,w(y)\,dy$. The book proves that $w'>0$ on $(0,0.004]$ and $w'<0$ on $[4.52,\infty)$; a grid evaluation (computed, not proved) finds exactly one sign change of $w'$ on $[0.004,4.6]$, at $x_m=0.015851\ldots$, and numerically $f$ rises from $0.365089$ to its maximum $0.366883$ at $x=0.024046$ and then decreases. Problem: prove that $w$ is unimodal, and prove that $f$ is unimodal with mode near $0.0240$, for instance from Steutel's equation and the unimodality of $w$. Yamazato's theorem, that every self decomposable law is unimodal, does not apply to the law of $\varphi$, which is infinitely divisible and not self decomposable. In Polson's 2026 book On Hilbert's 8th Problem the unimodality of $w$ is stated as a theorem, whose proof uses a bound on $w''$ that is not derived, and the problem asks only for the unimodality of $f$; the book proves only the two end ranges of $w'$ and asks for both.

Supersedes: Theorem 28.12 of Polson (2026) states that $w$ is unimodal, with a proof that uses a bound on $w''$ not derived there, and Open Problem 30.23 asks only for the unimodality of $f$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:unimodal`.

<a id="open_levy_density_positive_on_every_line_of_the_strip"></a>

### Infinite divisibility across the strip

Problem 7.8, page 241 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For $\sigma\in(0,1)$ let $\varphi_\sigma(s)=\zeta(\sigma)/\bigl((1-s/(1-\sigma)^2)\zeta(\sigma+\sqrt s)\bigr)$ and $w_\sigma(x)=e^{(1-\sigma)^2x}-\frac{1}{\sqrt{4\pi x}}\sum_{n\ge2}\Lambda(n)\,n^{-\sigma}e^{-(\log n)^2/4x}$, its canonical function, with $\Lambda$ the von Mangoldt function. Problem: is $w_\sigma(x)>0$ for every $\sigma\in(0,1)$ and every $x>0$? Positivity gives infinite divisibility of $\varphi_\sigma$ with Lévy density $w_\sigma(x)/x$. The book proves $w_{1/2}>0$ and checks $w_\sigma>0$ numerically only for $x\le10^3$ and at nine values of $\sigma$. For $\frac12<\sigma<1$ the method of the case $\sigma=\frac12$ may apply, since the terms of the zeros below the verified height are Poisson kernels of positive width, so the question reduces to the sign of the explicit pole and archimedean terms at large $x$, whose limit is governed by $\zeta'(\sigma)/\zeta(\sigma)>0$. For the completed functions $R_\sigma(s)=\xi(\sigma)/\xi(\sigma+\sqrt s)$ the answer is known: $R_\sigma$ is infinitely divisible for every $\sigma>\frac12$ with no hypothesis; since $\varphi_\sigma$ differs from $R_\sigma$ by the gamma factor and elementary factors in $\sqrt s$, this does not settle the question for $w_\sigma$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:strip`.

<a id="open_meixner_pollaczek_truncations_real_rooted"></a>

### Meixner Pollaczek truncations with real zeros

Problem 7.9, page 242 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $Z$ be Hardy's function and $p(\tau)=|\Gamma(\frac14+\frac{i\tau}2)|^2/(2^{3/2}\pi^{3/2})$ the Meixner density. The truncations of the expansion of the real even function $(\tau^2+\frac14)Z(\tau)$ in the even Meixner Pollaczek polynomials $P^{(1/4)}_{2m}(\tau/2;\frac\pi2)$, which are orthogonal for $p$, are real even polynomials. The symmetric expansion captures the first ordinate at degree $40$ but converges slowly, and the real zeros of its truncations do not track the ordinates of the zeros of the zeta function. Problem: is there a weighting of its coefficients, or a Jensen type construction from them, whose polynomials are real rooted and whose zeros converge to the ordinates?

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:mptrunc`.

<a id="open_nyman_beurling_and_the_sech_expansion"></a>

### Nyman Beurling criterion and the sech expansion

Problem 7.10, page 242 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. For $\operatorname{Re}\sqrt s>\frac12$ the book writes $1/\zeta(\frac12+\sqrt s)=1+\sum_{N\ge2}A_N\operatorname{sech}(\sqrt s\log N)$ with explicit Möbius type weights $A_N$ supported on the powers $N=n^m$ with $n$ squarefree and $m$ odd, and $A_N=\mu(N)/(2\sqrt N)$ for squarefree $N$ ($\mu$ the Möbius function). Problem: relate the completeness question of Lai and Vere-Jones for the dilates $\{\operatorname{sech}(\lambda x)\}$ to the Nyman Beurling criterion for the Riemann Hypothesis, and decide whether the sech expansion admits an $L^2$ version on the boundary $\operatorname{Re}\sqrt s=\frac12$.

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:nymanbeurling`.

<a id="open_product_representation_of_the_reciprocal_xi_law"></a>

### Product representations of the law of R

Problem 7.11, page 242 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. First stated in Polson (2026), On Hilbert's 8th Problem: the reciprocal Xi clock, Thorin measures, van Dantzig pairs, zero free approximations, hitting times and the reciprocal zeta function, first edition, Independently published.

Status: open.

Statement. Let $X_R$ be a random variable with Laplace transform $R(s)=\xi(\frac12)/\xi(\frac12+\sqrt s)$ and $X_\varphi$ one with Laplace transform $\varphi(s)=\zeta(\frac12)/\bigl((1-4s)\zeta(\frac12+\sqrt s)\bigr)$. The class of generalised gamma convolutions (GGC) is closed under independent sums, independent products and weak limits, so an unconditional representation of $X_R$ built from independent GGC variables by sums, products and weak limits would make $R$ a GGC and prove the Riemann Hypothesis; for $X_\varphi$ the book proves that no such representation exists, and the harmonic mean is not a closure operation of the class. The factorisation of $\varphi$ into $R$ and the gamma factor $\Gamma(\frac14+\frac{\sqrt s}2)\pi^{-\sqrt s/2}/\Gamma(\frac14)$ gives, formally, $X_R$ as $X_\varphi$ minus a positive stable law of index $\frac12$ and a compensated series of independent variables $Y_{c_n}$, where $\mathbb{E}e^{-sY_\kappa}=\kappa/(\kappa+\sqrt s)$ and $c_n=2n+\frac12$ for $n\ge0$. Question: can the compensation be absorbed into a product?

Source: [`ch/rz4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/rz4.tex), label `rz:prob:xrproduct`.

## Chapter 8. Probability, hitting times and other criteria

<a id="clock_not_first_passage_time_to_random_level"></a>

### The clock is not a Brownian first passage time to an independent random level

Proposition 8.2.1, page 252 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's completed zeta function, and let $T$ be the clock variable, the nonnegative random variable whose Laplace transform is $\mathbb{E}e^{-uT}=\xi(\tfrac12)/\xi(\tfrac12+\sqrt{u})$ for $u\ge0$ (it exists with no hypothesis on the zeros of $\xi$). For a drift $c\ge0$ and a level $a>0$, the first passage time of standard Brownian motion with drift $c$ (unit variance per unit time) to the level $a$ has the inverse Gaussian law, with Laplace transform $\exp\{a(c-\sqrt{c^2+2u})\}$. Then for no drift $c\ge0$ and no probability law $F$ on $(0,\infty)$ is $T$ equal in law to the first passage time of Brownian motion with drift $c$ to an independent random level $A$ with law $F$. For a drift $c<0$ the first passage time is infinite with probability $1-\mathbb{E}e^{2cA}>0$, so equality is impossible as well.

Source: [`ch/ch19.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch19.tex), label `prop:nogo`.

<a id="keiper_li_coefficients_as_thorin_functional"></a>

### Keiper Li coefficients as values of the Thorin functional

Proposition 8.3.1, page 261 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be Riemann's completed zeta function and, for $n\ge1$, let $\lambda_n=\sum_\rho\big[1-(1-1/\rho)^n\big]$ be the Keiper Li coefficient, the sum over the zeros $\rho$ of $\xi$ taken in conjugate pairs. Let $U$ be the Thorin functional of the zeros, which gives a unit point mass at $-(\rho-\tfrac12)^2$ to each pair $\{\rho,1-\rho\}$ of zeros of $\xi$, on the critical line or off it, so that $\langle U,\varphi\rangle=\sum_{\{\rho,1-\rho\}}\varphi\big(-(\rho-\tfrac12)^2\big)$. With no hypothesis on the zeros, $\lambda_n=\langle U,\varphi_n\rangle$, where $\varphi_n(z)=2-\Big(\frac{2\sqrt z+i}{2\sqrt z-i}\Big)^n-\Big(\frac{2\sqrt z-i}{2\sqrt z+i}\Big)^n$ is a rational function of $z$ with a pole of order $n$ at $z=-\tfrac14$ and $\varphi_n(z)\sim n^2/z$ as $z\to\infty$; for example $\varphi_1(z)=(z+\tfrac14)^{-1}$ and $\varphi_2(z)=4z(z+\tfrac14)^{-2}$. For $z>0$, $\varphi_n(z)=4\sin^2\big(n\arctan\tfrac1{2\sqrt z}\big)\ge0$. Consequently, by Li's criterion, the Riemann Hypothesis holds if and only if $\langle U,\varphi_n\rangle\ge0$ for every $n\ge1$.

Source: [`ch/ch20.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch20.tex), label `prop:ch20:thorin`.

<a id="backward_heat_flow_critical_point_collision_bound"></a>

### Critical point lower bound for the collision time of zeros under the backward heat flow

Theorem 8.4.6, page 278 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$; for a real polynomial $p$ and $D=d/dz$ one such family is $f_s=e^{sD^2}p=\sum_{k\ge0}s^kp^{(2k)}/k!$. For a real entire function $f$ of genus at most one with zeros $z_j$, Laguerre's quantity is $\mathcal{L}_f(x)=-\frac{d^2}{dx^2}\log|f(x)|=\big(f'(x)^2-f(x)f''(x)\big)/f(x)^2=\sum_j(x-z_j)^{-2}$, a real sum because nonreal zeros come in conjugate pairs. Fix a gap between two adjacent simple real zeros $z_k(s)<z_{k+1}(s)$ of $f_s$, which move continuously until they meet, if they do. When all zeros of $f_s$ are real, $f_s'/f_s$ falls from $+\infty$ to $-\infty$ across the gap and vanishes at exactly one point $c(s)$, the critical point of the gap; put $\ell(s)=\mathcal{L}_{f_s}(c(s))$. Let $S>0$. If all zeros of $f_\sigma$ are real and simple for $0\le\sigma<S$, then $\dot\ell\le\ell^2$ on $[0,S)$ and $\ell(\sigma)\le\ell(0)/(1-\sigma\ell(0))$ for $0\le\sigma<\min(S,1/\ell(0))$. In particular, if $p$ is a real polynomial with only real simple zeros, and the first collision of any two zeros of $e^{sD^2}p$ happens at $s=s_1$, between the two zeros that bound the gap $(z_k,z_{k+1})$ of $p$ at $s=0$, then $s_1\ge1/\mathcal{L}_p(c_k)$, where $c_k$ is the critical point of $p$ in that gap. The bound is attained by $p(z)=z^2-a^2$.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `thm:wigbiane:riccati`.

<a id="backward_heat_flow_biane_form_collision_bound"></a>

### Biane form lower bound for the collision time of zeros under the backward heat flow

Theorem 8.4.8, page 279 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$; for a real polynomial $p$ and $D=d/dz$ one such family is $f_s=e^{sD^2}p=\sum_{k\ge0}s^kp^{(2k)}/k!$. For a real entire function $f$ of genus at most one with zeros $z_j$, Laguerre's quantity is $\mathcal{L}_f(x)=-\frac{d^2}{dx^2}\log|f(x)|=\big(f'(x)^2-f(x)f''(x)\big)/f(x)^2=\sum_j(x-z_j)^{-2}$, a real sum because nonreal zeros come in conjugate pairs. Fix a gap between two adjacent simple real zeros $z_k(s)<z_{k+1}(s)$ of $f_s$, which move continuously until they meet, if they do. The quantity $\mathcal{L}_{f_s}$ tends to $+\infty$ at both ends of the gap; let $m(s)$ be its minimum over the gap and $D^+m$ the upper right Dini derivative of $m$. Let $S>0$. If all zeros of $f_\sigma$ are real and simple for $0\le\sigma<S$, then $D^+m\le m^2$ on $[0,S)$ and $m(\sigma)\le m(0)/(1-\sigma m(0))$ for $0\le\sigma<\min(S,1/m(0))$. In particular, if $p$ is a real polynomial with only real simple zeros, and the first collision of any two zeros of $e^{sD^2}p$ happens at $s=s_1$, between the two zeros that bound the gap $(z_k,z_{k+1})$ of $p$ at $s=0$, then $s_1\ge1/\min_{(z_k,z_{k+1})}\mathcal{L}_p$. This bound is at least as strong as the critical point bound $1/\mathcal{L}_p(c_k)$, where $c_k$ is the critical point of $p$ in the gap, because the minimum does not exceed the value at $c_k$, and it is attained by $p(z)=z^2-a^2$.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `thm:wigbiane:biane2`.

<a id="backward_heat_flow_local_bound_nonreal_zeros"></a>

### Local lower bound for the closing of a gap with a correction from the nonreal zeros

Theorem 8.4.11, page 281 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$, for example $f_s=e^{sD^2}p$ with $D=d/dz$ and $p$ a real polynomial. Let $\mathcal{L}(x)=-\frac{d^2}{dx^2}\log|f_s(x)|=\sum_j(x-z_j)^{-2}$, over the zeros $z_j$ of $f_s$, be Laguerre's quantity (a real sum, because nonreal zeros come in conjugate pairs). Fix a gap between two adjacent simple real zeros of $f_s$, which move continuously until they meet, if they do; $\mathcal{L}\to+\infty$ at both ends of the gap. Let $m(s)$ be the minimum of $\mathcal{L}$ over the gap, $D^+m$ its upper right Dini derivative, and at a time $s$ let $u$ be a point of the gap where the minimum is attained, with $m=m(s)=\mathcal{L}(u)>0$. For the conjugate pairs $a\pm ib$, $b>0$, of nonreal zeros of $f_s$, each pair counted once, put $B_p=\sum b^2|u-a|^{-p}$, and assume that $B_4$, $B_5$ and $B_6$ are finite. Then $\mathcal{L}''(u)-3m^2\le\eta m^2$, where $\eta=\big(36\,mB_4+108\,B_4^2+72\,(m+6B_4)^{1/2}B_5+120\,B_6\big)/m^2$, and $D^+m(s)\le(1+\eta)m^2$. Hence, if the two zeros bounding the gap meet at a time $s_c$, and $m>0$ and $\eta\le\bar\eta$ on $[0,s_c)$ for a constant $\bar\eta$ (with $\eta$ evaluated at a minimiser at each time), then $s_c\ge1/\big((1+\bar\eta)\,m(0)\big)$. All of this holds with the coefficient $72$ in $\eta$ replaced by $36$, by the sharp inequality $S_4\le\tfrac12S_2^2+\tfrac12S_2^{1/2}|S_3|$ for the power sums $S_k=\sum_jx_j^k$ of a real square summable sequence $(x_j)$.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `thm:wigbiane:local`.

<a id="backward_heat_flow_delta_collision_bound"></a>

### Lower bound for the closing of a gap from the quantity delta, with no knowledge of the nonreal zeros

Corollary 8.4.12, page 282 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$, for example $f_s=e^{sD^2}p$ with $D=d/dz$ and $p$ a real polynomial, and let $\mathcal{L}(x)=-\frac{d^2}{dx^2}\log|f_s(x)|$ be Laguerre's quantity. Fix a gap between two adjacent simple real zeros of $f_s$, let $m(s)$ be the minimum of $\mathcal{L}$ over the gap, $D^+m$ its upper right Dini derivative, and at a time $s$ let $u$ be a point of the gap where the minimum is attained, with $m=m(s)=\mathcal{L}(u)>0$. For the conjugate pairs $a\pm ib$, $b>0$, of nonreal zeros of $f_s$, each pair counted once, put $B_p=\sum b^2|u-a|^{-p}$, assume that $B_4$, $B_5$ and $B_6$ are finite, and put $\eta=\big(36\,mB_4+108\,B_4^2+72\,(m+6B_4)^{1/2}B_5+120\,B_6\big)/m^2$. Put $\delta=\mathcal{L}''(u)/m^2-3$, computed from $f_s$ and its first four derivatives at $u$. Then $D^+m(s)\le(1+\delta)m^2$ and $\delta\le\eta$, and $\delta\le0$ when all zeros of $f_s$ are real. Hence, if the two zeros bounding the gap meet at a time $s_c$, and $m>0$ and $\delta\le\bar\delta$ on $[0,s_c)$ for a constant $\bar\delta>-1$, then $s_c\ge1/\big((1+\bar\delta)\,m(0)\big)$; if $\delta\le0$ throughout, the collision cannot come before Biane's time $1/m(0)$, exactly as for functions with only real zeros.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `cor:wigbiane:delta`.

<a id="backward_heat_flow_collision_deadline"></a>

### Upper bound (deadline) for the closing time of a gap between real zeros

Theorem 8.4.13, page 282 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$, for example $f_s=e^{sD^2}p$ with $D=d/dz$ and $p$ a real polynomial, and let $\mathcal{L}(x)=-\frac{d^2}{dx^2}\log|f_s(x)|=\sum_j(x-z_j)^{-2}$ be Laguerre's quantity (a real sum, because nonreal zeros come in conjugate pairs). Fix a gap between two adjacent simple real zeros $z_k<z_{k+1}$ of $f_s$, let $m(s)$ be the minimum of $\mathcal{L}$ over the gap, and let $u$ be a point where it is attained, so $m=\mathcal{L}(u)$. Write $m=P+Q$, where $P=(u-z_k)^{-2}+(u-z_{k+1})^{-2}$ is the share of the two zeros bounding the gap and $Q=\sum(u-z)^{-2}$ over all the other zeros $z$ of $f_s$, and let $W=\sum|u-z|^{-4}$ over the nonreal zeros $z$ of $f_s$. While the two zeros bounding the gap are real, the right derivative $m'(s^+)$ exists and is the least value of $\mathcal{L}''(u)-2m^2$ over the minimisers $u$, and at every minimiser $\mathcal{L}''(u)-2m^2\ge m^2-6mQ+3Q^2-6W$. Suppose that, for constants $\rho$ and $\bar W$ with $9\rho^2+6\bar W>0$, $Q\le\rho$ and $W\le\bar W$ at every minimiser for as long as the two zeros bounding the gap are real, and that $m(0)>r_1$, where $r_{1,2}=3\rho\pm(9\rho^2+6\bar W)^{1/2}$. Then the gap closes, one of the two zeros bounding it meeting another zero, at a time $s_c\le s^*=\frac{1}{r_1-r_2}\log\frac{m(0)-r_2}{m(0)-r_1}$. If all zeros are real, then $W=0$, so $\bar W=0$ may be taken, and for $\rho>0$ one has $s^*=(6\rho)^{-1}\log\big(m(0)/(m(0)-6\rho)\big)=1/m(0)+3\rho/m(0)^2+12\rho^2/m(0)^3+\cdots$.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `thm:wigbiane:upper`.

<a id="backward_heat_flow_collision_deadline_kappa_form"></a>

### Deadline for the closing of a gap in terms of the quantity kappa

Corollary 8.4.14, page 283 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f_s$, $s\ge0$, be real entire functions of genus at most one, smooth in $(z,s)$, with $\partial_sf_s=\partial_z^2f_s$, and let $\mathcal{L}(x)=-\frac{d^2}{dx^2}\log|f_s(x)|=\sum_j(x-z_j)^{-2}$ be Laguerre's quantity. Fix a gap between two adjacent simple real zeros $z_k<z_{k+1}$ of $f_s$, let $m(s)$ be the minimum of $\mathcal{L}$ over the gap, let $u$ be a point where it is attained, and put $x_j=1/(u-z_j)$ for the zeros $z_j$ of $f_s$, $P=x_k^2+x_{k+1}^2$, $Q=\sum_{j\ne k,k+1}x_j^2$, so that $m=P+Q$, and $W=\sum|x_j|^4$ over the nonreal zeros. Put $\kappa=\mathcal{L}''(u)-3P^2$, computed from $f_s$ at $u$ and the two zeros bounding the gap. Then $\kappa=3\big(x_k^2-x_{k+1}^2\big)^2+6\sum_{j\ne k,k+1}x_j^4\ge-6W$, with $\kappa\ge0$ when all zeros are real. Put $\bar\kappa=\max\big(0,\sup(-\kappa)\big)$, the supremum over the minimisers at the times before the collision. If $Q\le\rho$ at every minimiser for as long as the two zeros bounding the gap are real, with $9\rho^2+\bar\kappa>0$, and $m(0)>r_1$, where $r_{1,2}=3\rho\pm(9\rho^2+\bar\kappa)^{1/2}$, then the gap closes, one of the two zeros bounding it meeting another zero, at a time $s_c\le\frac{1}{r_1-r_2}\log\frac{m(0)-r_2}{m(0)-r_1}$. This bound is never weaker than the same bound with $6\bar W$ in place of $\bar\kappa$, where $\bar W$ is a constant with $W\le\bar W$ at every minimiser for as long as the two zeros bounding the gap are real, since then $\bar\kappa\le6\bar W$.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `cor:wigbiane:kappa`.

<a id="nonreal_zeros_control_along_backward_heat_flow"></a>

### Control of the nonreal zeros of the de Bruijn Newman function at negative times

Problem 8.1, page 286 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For real $t$ let $H_t(z)=\int_0^\infty e^{tu^2}\Phi(u)\cos(zu)\,du$ be the de Bruijn Newman function, with $\Phi(u)=\sum_{n\ge1}\big(2\pi^2n^4e^{9u}-3\pi n^2e^{5u}\big)\exp(-\pi n^2e^{4u})$, so that $H_0(z)=\tfrac18\Xi(z/2)$ with $\Xi(\theta)=\xi(\tfrac12+i\theta)$ for Riemann's completed zeta function $\xi$. As $t$ decreases from $0$ the zeros of $H_t$ move by the backward heat flow, and the book bounds from below and from above, in terms of $H_0$ alone, the time at which a gap between two adjacent real zeros closes, under hypotheses along the flow that the computation checks on a grid of times but cannot yet guarantee in advance. Problem: find a bound for how far the nonreal zeros of $H_t$, $t<0$, travel, and for how close to a gap between real zeros they come. No such bound is known, and such control would make the lower and upper bounds for the closing time of a gap unconditional statements about $H_t$, computed from $H_0$ alone.

Source: [`ch/ch21.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch21.tex), label `prob:wigbiane:control`.

<a id="max_representation_of_ramanujan_delta_l_function"></a>

### Max representation of the L function of Ramanujan's Delta

Theorem 8.5.2, page 296 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $E_{n,j}$, $n\ge1$, $1\le j\le24$, be independent exponential random variables, $E_{n,j}$ of rate $2\pi n$, and let $M^*=\sup_{n,j}E_{n,j}$ (almost surely attained and finite, with $\Pr(M^*\le y)=\prod_{n\ge1}(1-e^{-2\pi ny})^{24}$ for $y>0$). Let $\Delta(z)=\eta(z)^{24}=\sum_{n\ge1}\tau(n)e^{2\pi inz}$ be Ramanujan's cusp form, $\eta$ Dedekind's eta function, $L(s,\Delta)=\sum_{n\ge1}\tau(n)n^{-s-11/2}$ its L function, normalised so that the critical line is $\operatorname{Re}s=\frac12$, and $\xi(s,\Delta)=(2\pi)^{-s-11/2}\Gamma(s+\frac{11}2)L(s,\Delta)$ its completed function. With $\Gamma(a,x)=\int_x^\infty e^{-v}v^{a-1}\,dv$ and $Q(a,x)=\Gamma(a,x)/\Gamma(a)$, for every complex $s$, $\xi(s,\Delta)=(2\pi)^{-s-11/2}\,\mathbb{E}\,\Gamma(s+\frac{11}2,2\pi M^*)$ and $L(s,\Delta)=\mathbb{E}\,Q(s+\frac{11}2,2\pi M^*)$, and both expectations converge absolutely. In particular $L(s,\Delta)=0$ at the points $s=-\frac{11}2-k$, $k=0,1,2,\dots$, where $Q(-k,x)=0$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `thm:maxexp:max`.

<a id="gamma_mixture_and_mixed_poisson_forms_of_delta_l_values"></a>

### Gamma mixture and mixed Poisson forms of L values of Delta

Corollary 8.5.3, page 297 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M^*$ be the maximum of independent exponential random variables, $24$ of each rate $2\pi n$, $n\ge1$, with distribution function $F(y)=\prod_{n\ge1}(1-e^{-2\pi ny})^{24}$, and let $L(s,\Delta)$ be the L function of Ramanujan's cusp form $\Delta$, normalised so that the critical line is $\operatorname{Re}s=\frac12$. Let $G_a$ be a gamma variable of shape $a>0$ and scale one and $N$ a Poisson process of unit rate, both independent of $M^*$. For real $s>-\frac{11}2$, with $a=s+\frac{11}2$, $L(s,\Delta)=\Pr(G_a>2\pi M^*)=\mathbb{E}\,F(G_a/2\pi)=\mathbb{E}\prod_{n\ge1}(1-e^{-nG_a})^{24}$, and for every integer $k\ge1$, $L(k-\frac{11}2,\Delta)=\Pr\big(N(2\pi M^*)\le k-1\big)=\mathbb{E}\Big[e^{-2\pi M^*}\sum_{j<k}\frac{(2\pi M^*)^j}{j!}\Big]$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `cor:maxexp:forms`.

<a id="delta_and_chi12_l_values_as_infinitely_divisible_probabilities"></a>

### L values of Delta and chi12 as infinitely divisible probabilities

Proposition 8.5.7, page 303 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M^*$ be the maximum of independent exponential random variables, $24$ of each rate $2\pi n$, $n\ge1$, and for real $a>0$ let $G_a$ be a gamma variable of shape $a$ and scale one independent of $M^*$ and $W_a=G_a/(2\pi M^*)$. Let $L(s,\Delta)$ be the L function of Ramanujan's cusp form $\Delta$, normalised so that the critical line is $\operatorname{Re}s=\frac12$. Then $W_a$ is a scale mixture of Gamma($a$) laws, the mixing law being that of $1/(2\pi M^*)$, and $L(a-\frac{11}2,\Delta)=\Pr(W_a>1)$. For $0<a\le2$, that is $-\frac{11}2<s\le-\frac72$, $W_a$ is infinitely divisible. In the same way, with $M_1$ the maximum of independent exponential random variables, one of each rate $2\pi n$, $n\ge1$, $G_{s/2}$ a gamma variable of shape $s/2$ and scale one independent of $M_1$ and $\chi_{12}$ the real primitive character modulo $12$ (equal to $1$ at $\pm1$ and to $-1$ at $\pm5$ modulo $12$), $L(s,\chi_{12})=\Pr\big(12G_{s/2}/(\pi M_1)>1\big)$ for real $s>0$, where $12G_{s/2}/(\pi M_1)$ is infinitely divisible for $s\le4$, and is a mixture of exponential laws for $s\le2$, in particular at the centre $s=\frac12$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:maxexp:krist`.

<a id="laplace_transform_of_max_of_exponentials_m1_closed_form"></a>

### Laplace transform of the maximum M1 of exponentials in closed form

Proposition 8.5.8, page 304 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M_1$ be the maximum of independent exponential random variables, one of each rate $2\pi n$, $n\ge1$, so that $\Pr(M_1\le y)=\prod_{n\ge1}(1-e^{-2\pi ny})$ for $y>0$, and put $b_n=\pi(n^2-1)/12$. For $\theta>0$ let $w=\sqrt{1-12\theta/\pi}$, positive or purely imaginary. Then $\mathbb{E}e^{-\theta M_1}=\frac{\pi(1-w^2)\sin(\pi w/3)}{2\sqrt3\,w\cos(\pi w/2)}=\prod_{k\ge1}\frac{1+\theta/b_{6k}}{(1+\theta/b_{6k-1})(1+\theta/b_{6k+1})}$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:ch16:M1`.

<a id="renyi_type_decomposition_of_max_of_exponentials_m1"></a>

### Rényi type decomposition of the maximum M1 of exponentials

Corollary 8.5.9, page 305 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M_1$ be the maximum of independent exponential random variables, one of each rate $2\pi n$, $n\ge1$, and $b_n=\pi(n^2-1)/12$. Then $M_1\overset{d}{=}\sum_{k\ge1}Y_k$ with independent $Y_k$, where $Y_k$ is exponential of rate $b_{6k-1}$ with probability $p_k$ and of rate $b_{6k+1}$ with probability $1-p_k$, and $p_k=\big(1-\frac{b_{6k-1}}{b_{6k}}\big)\frac{b_{6k+1}}{b_{6k+1}-b_{6k-1}}=\frac{(12k-1)(3k+1)}{2(36k^2-1)}$, so that $p_1=\frac{22}{35}$, $p_2=\frac{161}{286}$, and $p_k\downarrow\frac12$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `cor:maxexp:renyi`.

<a id="levy_density_and_bondesson_function_of_max_of_exponentials_m1"></a>

### Lévy density and Bondesson function of the maximum M1

Proposition 8.5.11, page 306 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M_1$ be the maximum of independent exponential random variables, one of each rate $2\pi n$, $n\ge1$, and $b_n=\pi(n^2-1)/12$. Put $a(n)=1$ for $n\equiv\pm1\pmod6$, $a(n)=-1$ for $6\mid n$ and $a(n)=0$ otherwise, $k_1(t)=\sum_{n\ge2}a(n)e^{-b_nt}$, so that $-\log\mathbb{E}e^{-\theta M_1}=\int_0^\infty(1-e^{-\theta t})\,\nu_1(t)\,dt$ with the Lévy density $\nu_1(t)=k_1(t)/t$ (the canonical function of $M_1$ is $k_1$), and let $U_1=\sum_{n\ge2}a(n)\delta_{b_n}$ be the signed Thorin measure of $M_1$ and $N_1(z)=U_1((0,z])$ its Bondesson function. Then $k_1(t)=\frac12e^{\pi t/12}\Big(1+\frac{\vartheta_0-3\vartheta_{1/2}-4\vartheta_{1/3}}{\sqrt{3t}}\Big)-1$ with $\vartheta_\alpha=\sum_{m\in\mathbb{Z}}e^{-12\pi(m+\alpha)^2/t}$, so that $\nu_1(t)=(2\sqrt3)^{-1}t^{-3/2}-(2t)^{-1}+O(t^{-1/2})$ as $t\downarrow0$ and $\nu_1(t)\sim e^{-2\pi t}/t$ as $t\to\infty$. Moreover $N_1(z)\ge0$ for all $z>0$, and $N_1(z)\sim\frac13\sqrt{3z/\pi}$ as $z\to\infty$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:maxexp:levy`.

<a id="thorin_measure_of_m1_difference_of_two_ggc_thorin_measures"></a>

### Thorin measure of M1 as a difference of two GGC Thorin measures

Theorem 8.5.12, page 307 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $b_n=\pi(n^2-1)/12$, let $E_n$ and $E'_m$ be independent standard exponential variables, $S=\sum_{n\text{ odd},\,n\ge3}E_n/b_n$ and $Z=\sum_{m\ge1}E'_m/b_{3m}$, with $Z$ independent of $M_1$, the maximum of independent exponential random variables, one of each rate $2\pi n$, $n\ge1$. Then $S$ and $Z$ are GGCs (generalised gamma convolutions: laws on $[0,\infty)$ with Laplace transform $\exp(-d\theta-\int\log(1+\theta/z)\,U(dz))$ for a drift $d\ge0$ and a positive measure $U$, the Thorin measure), with $\mathbb{E}e^{-\theta S}=\frac\pi4\,\frac{1-w^2}{\cos(\pi w/2)}$ and $\mathbb{E}e^{-\theta Z}=\frac{\sqrt3\,w/2}{\sin(\pi w/3)}$, where $w=\sqrt{1-12\theta/\pi}$, and $M_1+Z\overset{d}{=}S$. Equivalently, the signed Thorin measure of $M_1$, the signed measure $U_1=\sum_{n\ge2}a(n)\delta_{b_n}$ with $-\log\mathbb{E}e^{-\theta M_1}=\int\log(1+\theta/z)\,U_1(dz)$, where $a(n)=1$ for $n\equiv\pm1\pmod6$, $a(n)=-1$ for $6\mid n$ and $a(n)=0$ otherwise, equals $U_S-U_Z$ with $U_S=\sum_{n\text{ odd},\,n\ge3}\delta_{b_n}$ and $U_Z=\sum_{m\ge1}\delta_{b_{3m}}$, the atoms at $b_{3m}$ with $m$ odd cancelling.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `thm:maxexp:ggc`.

<a id="certified_zeros_of_laplace_transform_of_max_of_exponentials_m_star"></a>

### Certified zeros of the Laplace transform of the maximum M* of exponentials

Proposition 8.5.13, page 309 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M^*$ be the maximum of independent exponential random variables, $24$ of each rate $2\pi n$, $n\ge1$, let $\tau(n)$ be Ramanujan's tau function (the coefficients of $\Delta=\eta^{24}$), and let $\psi^*(\theta)=\mathbb{E}e^{-\theta M^*}$, which extends to a meromorphic function of $\theta\in\mathbb{C}$ with poles only at $\theta=-2\pi m$, $m\ge1$ (simple, with residue $-2\pi m\,\tau(m+1)$). (1) For each of the $25$ values $m\le60$ with $\tau(m+1)\tau(m+2)>0$, $\psi^*$ changes sign in an interval of length $2\times10^{-20}$ inside the gap $(-2\pi(m+1),-2\pi m)$; the zeros begin $\theta/2\pi=-5.48830041$, $-8.66373325$, $-11.4203774$, $-13.2533839$, $-14.7973058$, and the last is $-60.7670708$. (2) $\psi^*$ has exactly one zero in the disc of radius $10^{-7}$ about each of the ten points $-155.5940822775511402804\pm2.7631188861212609656\,i$, $-215.6211277570694988944\pm2.9647292148515237348\,i$, $-307.4008549057072052743\pm3.2283287024583906394\,i$, $-322.6135369899940182007\pm3.0449254921259967464\,i$, $-347.8816009468099010890\pm2.7171007701068929577\,i$, and it lies within $10^{-19}$ of that point.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:maxexp:psizeros`.

<a id="max_of_exponentials_m_star_not_in_bondesson_class"></a>

### The maximum M* of exponentials is not in the Bondesson class

Corollary 8.5.14, page 310 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M^*$ be the maximum of independent exponential random variables, $24$ of each rate $2\pi n$, $n\ge1$. Then $M^*$ is not in the Bondesson class BO (the infinitely divisible laws on $[0,\infty)$ with a completely monotone Lévy density, and a drift $d\ge0$); in particular $M^*$ is not a generalised gamma convolution.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `cor:maxexp:notBO`.

<a id="open_infinite_divisibility_of_maxima_of_exponentials"></a>

### Open question on infinite divisibility of maxima of exponentials

Question 8.5.15, page 311 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Is every maximum of independent exponential random variables infinitely divisible? In particular, is $M^*$ infinitely divisible, where $M^*$ is the maximum of independent exponential random variables, $24$ of each rate $2\pi n$, $n\ge1$?

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:maxexp:id`.

<a id="open_real_zeros_of_cosh_kernel_with_exponential_multiplier"></a>

### Open question on real zeros of cosh kernels with exponential multipliers

Question 8.5.16, page 315 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For which $a>0$, $c_m\ge0$ and $\beta_m>0$ does the Fourier transform of $e^{-a\cosh u}\exp\Big(-\sum_mc_m\big(e^{-\beta_me^u}+e^{-\beta_me^{-u}}\big)\Big)$ have only real zeros? (The kernels of the completed L functions of Ramanujan's cusp form, of the character $\chi_{12}$ and of the elliptic curve of conductor $11$ have this form.)

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:maxexp:mult`.

<a id="polya_kernel_of_imaginary_quadratic_field"></a>

### Pólya kernel of an imaginary quadratic field

Theorem 8.5.17, page 318 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $d<0$ be a fundamental discriminant, $K=\mathbb{Q}(\sqrt d)$, $\chi_d$ the Kronecker symbol of $d$ (a real primitive odd character modulo $|d|$), $h$ the class number of $K$ and $w$ the number of roots of unity in $K$. Let $\zeta_K(s)=\zeta(s)L(s,\chi_d)=\sum_{n\ge1}a_nn^{-s}$, where $a_n=\sum_{m\mid n}\chi_d(m)\ge0$ is the number of ideals of norm $n$, and put $\Lambda_K(s)=(\sqrt{|d|}/2\pi)^s\Gamma(s)\zeta_K(s)$, $\xi_K(s)=2s(s-1)\Lambda_K(s)$, $\vartheta_K(x)=\sum_{n\ge1}a_ne^{-2\pi nx/\sqrt{|d|}}$, $\Psi(y)=y(y-2)e^{-y}$, $H_K(x)=(x^2\vartheta_K'(x))'$ and $\Phi_K(u)=2e^{u/2}H_K(e^u)$. Then (i) $\vartheta_K(1/x)=x\,\vartheta_K(x)+\frac hw(x-1)$ for $x>0$. (ii) $H_K(x)=\sum_{n\ge1}a_n\Psi(2\pi nx/\sqrt{|d|})$ and $H_K(1/x)=x\,H_K(x)$; the kernel $\Phi_K$ is even, analytic in the strip $|\operatorname{Im}u|<\frac\pi2$ and decays doubly exponentially. (iii) $\xi_K(\frac12+s)=\int_{\mathbb{R}}\Phi_K(u)e^{su}\,du$ for every $s\in\mathbb{C}$; in particular $\xi_K(\frac12+i\theta)=\int_{\mathbb{R}}\Phi_K(u)\cos(\theta u)\,du$. (iv) $\Phi_K>0$ on $\mathbb{R}$ if and only if $H_K>0$ on $[1,\sqrt{|d|}/\pi]$; in particular $\Phi_K>0$ when $|d|<\pi^2$, that is for $d=-3,-4,-7,-8$. (v) If $\Phi_K>0$, then $\zeta_K(\frac12)<0$ and $L(\frac12,\chi_d)>0$; in particular $L(\frac12,\chi_d)>0$ for $d=-3,-4,-7,-8$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `thm:sqfree:kernelK`.

<a id="polya_kernel_takes_negative_values_when_two_and_three_inert"></a>

### Pólya kernel of an imaginary quadratic field takes negative values when 2 and 3 are inert

Theorem 8.5.18, page 320 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. For a fundamental discriminant $d<0$ let $K=\mathbb{Q}(\sqrt d)$, $\chi_d$ the Kronecker symbol of $d$, $a_n=\sum_{m\mid n}\chi_d(m)$ the number of ideals of norm $n$ in $K$, $\vartheta_K(x)=\sum_{n\ge1}a_ne^{-2\pi nx/\sqrt{|d|}}$, $H_K(x)=(x^2\vartheta_K'(x))'$ and $\Phi_K(u)=2e^{u/2}H_K(e^u)$ the Pólya kernel of $K$ (the kernel with $\xi_K(\frac12+s)=\int_{\mathbb{R}}\Phi_K(u)e^{su}\,du$ for the completed Dedekind zeta function $\xi_K$), and let $K_{\chi_d}(u)=2e^{3u/2}\sum_{n\ge1}n\chi_d(n)e^{-\pi n^2e^{2u}/|d|}$ be the kernel of the odd character $\chi_d$, with $\xi(\frac12+s,\chi_d)=\int_{\mathbb{R}}K_{\chi_d}(u)e^{su}\,du$ for $\xi(s,\chi_d)=(|d|/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_d)$. Let $d<0$ be a fundamental discriminant with $\chi_d(2)=\chi_d(3)=-1$, that is $d\equiv5\pmod{24}$. Then $H_K(x)<0$ at $x=1.42\sqrt{|d|}/2\pi$, so $\Phi_K$ takes negative values, and so does the kernel $K_{\chi_d}$. There are infinitely many such $d$, for instance $d=-p$ for every prime $p\equiv19\pmod{24}$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `thm:sqfree:inert`.

<a id="open_density_of_positive_polya_kernels_imaginary_quadratic"></a>

### Open question on the density of positive Pólya kernels

Question 8.5.19, page 322 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For a fundamental discriminant $d<0$ let $\Phi_K$ be the Pólya kernel of $K=\mathbb{Q}(\sqrt d)$, the even kernel with $\xi_K(\frac12+s)=\int_{\mathbb{R}}\Phi_K(u)e^{su}\,du$ for the completed Dedekind zeta function $\xi_K(s)=2s(s-1)(\sqrt{|d|}/2\pi)^s\Gamma(s)\zeta_K(s)$. Does the proportion of fundamental discriminants $-X\le d<0$ with $\Phi_K>0$ tend to a limit as $X\to\infty$, and is the limit positive? Is $\Phi_K$ positive for infinitely many $d$, for instance among the $d$ in which $2$ and $3$ split (that is, $\chi_d(2)=\chi_d(3)=1$ for the Kronecker symbol $\chi_d$), where the proportion up to $20000$ is $676/767$?

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:sqfree:density`.

<a id="open_positivity_criterion_polya_kernel_real_quadratic"></a>

### Open question on positivity criteria for real quadratic fields

Question 8.5.20, page 322 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For a number field $K$ with $r_1$ real and $r_2$ complex places the completed Dedekind zeta function $\xi_K$ has the gamma factor $\Gamma_{\mathbb{R}}(s)^{r_1}\Gamma_{\mathbb{C}}(s)^{r_2}$, with $\Gamma_{\mathbb{R}}(s)=\pi^{-s/2}\Gamma(\frac s2)$ and $\Gamma_{\mathbb{C}}(s)=2(2\pi)^{-s}\Gamma(s)$, and its Pólya kernel is the kernel $\Phi_K$ with $\xi_K(\frac12+s)=\int_{\mathbb{R}}\Phi_K(u)e^{su}\,du$. For real quadratic fields and for fields of higher degree the gamma factor replaces the exponentials in the theta series of an imaginary quadratic field by other kernels, for $\Gamma_{\mathbb{R}}(s)^2$ by a Bessel function $K_0$. For an imaginary quadratic field of discriminant $d$, $\Phi_K>0$ on $\mathbb{R}$ if and only if $H_K>0$ on $[1,\sqrt{|d|}/\pi]$ for an explicit function $H_K$ of one variable. Is there a criterion for the positivity of the Pólya kernel of $\xi_K$ for real quadratic fields and for fields of higher degree that is as simple as this one?

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:sqfree:realquad`.

<a id="levy_form_of_dirichlet_beta_function_from_basepoint_one"></a>

### Lévy form of the Dirichlet beta function at basepoints from one

Proposition 8.5.23, page 326 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $L(s,\chi_4)=\sum_{n\ge0}(-1)^n(2n+1)^{-s}$ be the Dirichlet beta function ($\chi_4$ the nonprincipal character modulo $4$), $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ its completed function (entire, with $\xi(s,\chi_4)=\xi(1-s,\chi_4)$) and $\Lambda$ the von Mangoldt function. For $\alpha\ge1$, $\log\frac{\xi(\alpha+s,\chi_4)}{\xi(\alpha,\chi_4)}=b_\alpha s+\int_0^\infty\big(e^{-sx}-1+sx\big)\,\Pi_\alpha(dx)$ with $b_\alpha=\frac{\xi'}{\xi}(\alpha,\chi_4)$ and the signed Lévy measure $\Pi_\alpha(dx)=\frac{e^{-\alpha x}}{2x\sinh x}\,dx+\sum_{n\ge2}\frac{\chi_4(n)\Lambda(n)}{n^\alpha\log n}\,\delta_{\log n}(dx)$, whose canonical measure $x\Pi_\alpha(dx)$ has the positive density $e^{-\alpha x}/(2\sinh x)$ and the atoms $\chi_4(n)\Lambda(n)n^{-\alpha}$ at $\log n$. For $\alpha>1$ the form holds for $\operatorname{Re}s>1-\alpha$, with absolute convergence. At $\alpha=1$ it holds for $\operatorname{Re}s>0$, the prime sum converging conditionally, as a limit over $n\le X$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:levy`.

<a id="thorin_measure_of_dirichlet_beta_function_from_basepoint_one"></a>

### Thorin measure of the Dirichlet beta function from basepoint one

Proposition 8.5.24, page 328 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $L(s,\chi_4)=\sum_{n\ge0}(-1)^n(2n+1)^{-s}$ be the Dirichlet beta function, $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ its completed function and $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_4)$. For every $\alpha\ge1$, $u\mapsto\xi(\alpha,\chi_4)/\xi(\alpha+\sqrt u,\chi_4)$ is the Laplace transform of a generalised gamma convolution whose Thorin measure has, in the chart $\theta=\sqrt z$, the positive density $v_\alpha(\theta)/\pi$ (that is, it is the image of $v_\alpha(\theta)\,d\theta/\pi$ under $\theta\mapsto\theta^2$), and $\log\frac{\xi(\alpha+\sqrt u,\chi_4)}{\xi(\alpha,\chi_4)}=\frac1\pi\int_0^\infty\log\Big(1+\frac u{\theta^2}\Big)v_\alpha(\theta)\,d\theta$ for $u\ge0$. For $\alpha>\frac12$ it is the Laplace transform of a generalised gamma convolution if and only if no zero of $L(s,\chi_4)$ has real part greater than $\alpha$. GRH for $\chi_4$ holds if and only if $v_\alpha\ge0$ for every $\alpha>\frac12$, and the Thorin measures at resolution $\varepsilon=\alpha-\frac12$ converge vaguely, as $\varepsilon\downarrow0$, to the sum of the point masses $\delta_{\gamma^2}$ over the zeros $\frac12+i\gamma$ of $L(s,\chi_4)$ on the critical line (under GRH this sum is the Thorin measure of $\xi(\frac12,\chi_4)/\xi(\frac12+\sqrt u,\chi_4)$).

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:thorin`.

<a id="termwise_positivity_thresholds_dirichlet_beta_and_chi12"></a>

### Termwise positivity thresholds for the Dirichlet beta function and chi12

Proposition 8.5.25, page 329 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ be the completed Dirichlet beta function, $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_4)$, $\psi=\Gamma'/\Gamma$ and $\Lambda$ the von Mangoldt function. For $\alpha>1$ and real $\theta$, $v_\alpha(\theta)\ge A^\infty_4(\alpha,\theta)-S_4(\alpha)$, where $A^\infty_4(\alpha,\theta)=\frac12\log\frac4\pi+\frac12\operatorname{Re}\psi\Big(\frac{\alpha+1+i\theta}2\Big)$ and $S_4(\alpha)=\sum_{n\text{ odd}}\frac{\Lambda(n)}{n^\alpha}=-\frac{\zeta'}\zeta(\alpha)-\frac{\log2}{2^\alpha-1}$, and $A^\infty_4(\alpha,\theta)$ increases in $|\theta|$. The bound is positive at every $\theta$ exactly when $\alpha>\alpha_c(\chi_4)=2.2823106368$, the root of $A^\infty_4(\alpha,0)=S_4(\alpha)$. For the real primitive character $\chi_{12}$ modulo $12$, with $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(\frac s2)L(s,\chi_{12})$ and $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_{12})$, the same bound with $A^\infty_{12}(\alpha,\theta)=\frac12\log\frac{12}\pi+\frac12\operatorname{Re}\psi\Big(\frac{\alpha+i\theta}2\Big)$ and $S_{12}(\alpha)=S_4(\alpha)-\frac{\log3}{3^\alpha-1}$ is positive at every $\theta$ exactly when $\alpha>\alpha_c(\chi_{12})=1.8348743295$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:termwise`.

<a id="no_bound_from_prime_weight_sizes_beats_absolute_values_kronecker"></a>

### No bound from the sizes of the prime weights beats the absolute values

Proposition 8.5.26, page 330 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $L(s,\chi_4)=\sum_{n\ge0}(-1)^n(2n+1)^{-s}$ be the Dirichlet beta function, $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ its completed function, $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_4)$, $A^\infty_4(\alpha,\theta)=\frac12\log\frac4\pi+\frac12\operatorname{Re}\psi\big(\frac{\alpha+1+i\theta}2\big)$ with $\psi=\Gamma'/\Gamma$, and $S_4(\alpha)=\sum_{n\text{ odd}}\Lambda(n)n^{-\alpha}$ with $\Lambda$ the von Mangoldt function. For $\alpha>1$ and every $T>0$, $\inf_{\theta\ge T}\big(v_\alpha(\theta)-A^\infty_4(\alpha,\theta)\big)=\inf_{\theta\ge T}\operatorname{Re}\frac{L'}L(\alpha+i\theta,\chi_4)=-S_4(\alpha)$, so no lower bound for $\operatorname{Re}\frac{L'}L(\alpha+i\theta,\chi_4)$ better than $-S_4(\alpha)$ holds at all large heights. The bound $v_\alpha\ge A^\infty_4-S_4$ yields $v_\alpha(\theta)>0$ exactly for $|\theta|>t_0(\alpha)$, where $A^\infty_4(\alpha,t_0)=S_4(\alpha)$, for $1<\alpha\le\alpha_c$ (so that $t_0(\alpha_c)=0$), and at every real $\theta$ for $\alpha>\alpha_c$, with $\alpha_c=2.2823106368$ the threshold at which $A^\infty_4(\alpha,0)=S_4(\alpha)$, and $\log t_0(\alpha)=\frac2{\alpha-1}-2\gamma_E-\log\frac8\pi+O(\alpha-1)$ as $\alpha\downarrow1$, where $\gamma_E$ is Euler's constant and $-2\gamma_E-\log\frac8\pi=-2.0891$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:kron`.

<a id="thorin_density_of_dirichlet_beta_as_difference_over_gaussian_field"></a>

### Thorin density of the Dirichlet beta function as a difference over the Gaussian field

Proposition 8.5.27, page 332 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$, let $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ be the completed Dirichlet beta function, $K=\mathbb{Q}(i)$ with $\zeta_K(s)=\zeta(s)L(s,\chi_4)$, and $\xi_K(s)=2s(s-1)\pi^{-s}\Gamma(s)\zeta_K(s)=\xi(s)\,\xi(s,\chi_4)$. Write $v^K_\alpha(\theta)=\operatorname{Re}\frac{\xi_K'}{\xi_K}(\alpha+i\theta)$, $v^\xi_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta)$ and $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_4)$. On every vertical line, $v_\alpha=v^K_\alpha-v^\xi_\alpha$. For $\alpha>1$ the canonical measure of the Lévy form of $\xi_K$ at the basepoint $\alpha$ is $x\Pi^K_\alpha(dx)=\Big(\frac{e^{-(\alpha+1)x}}{1-e^{-x}}-e^{-(\alpha-1)x}\Big)dx+\sum_{n\ge2}\frac{(1+\chi_4(n))\Lambda(n)}{n^\alpha}\,\delta_{\log n}(dx)$, and the canonical measure $x\Pi_\alpha$ of the Lévy form of $\xi(\cdot,\chi_4)$ at $\alpha$ (with density $e^{-\alpha x}/(2\sinh x)$ and atoms $\chi_4(n)\Lambda(n)n^{-\alpha}$ at $\log n$) equals $x\Pi^K_\alpha-x\Pi^\xi_\alpha$, where $x\Pi^\xi_\alpha(dx)=\Big(\frac{e^{-(\alpha+2)x}}{1-e^{-2x}}-e^{-(\alpha-1)x}\Big)dx+\sum_{n\ge2}\frac{\Lambda(n)}{n^\alpha}\,\delta_{\log n}(dx)$ is that of $\xi$: the pole densities $-e^{-(\alpha-1)x}$ cancel, the gamma densities differ by $e^{-\alpha x}/(2\sinh x)$, and the atoms by $\chi_4(n)\Lambda(n)n^{-\alpha}$. Here $\Lambda$ is the von Mangoldt function.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:factor`.

<a id="grh_chi4_iff_gaussian_field_thorin_density_dominates_xi"></a>

### GRH for chi4 and domination of Thorin densities over Q(i)

Corollary 8.5.28, page 332 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on Part (iii) assumes the Riemann hypothesis for the Riemann zeta function; parts (i) and (ii) hold without any hypothesis.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$, $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ the completed Dirichlet beta function, $\xi_K(s)=\xi(s)\,\xi(s,\chi_4)$ the completed Dedekind zeta function of $K=\mathbb{Q}(i)$, and $v^K_\alpha(\theta)=\operatorname{Re}\frac{\xi_K'}{\xi_K}(\alpha+i\theta)$, $v^\xi_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta)$. (i) For every $\alpha\ge1$ there are independent GGC variables (generalised gamma convolutions) $H^\xi_\alpha$ and $H^{(4)}_\alpha$ with Laplace transforms $\xi(\alpha)/\xi(\alpha+\sqrt u)$ and $\xi(\alpha,\chi_4)/\xi(\alpha+\sqrt u,\chi_4)$, whose sum $H^K_\alpha=H^\xi_\alpha+H^{(4)}_\alpha$ is the GGC with Laplace transform $\xi_K(\alpha)/\xi_K(\alpha+\sqrt u)$ and Thorin density $v^K_\alpha/\pi$ (in the chart $\theta=\sqrt z$). (ii) GRH for $\chi_4$ holds if and only if $v^K_\alpha\ge v^\xi_\alpha$ on every line $\operatorname{Re}s=\alpha>\frac12$: the Thorin density of $\xi_K$ dominates that of $\xi$ at every basepoint. (iii) Under the Riemann hypothesis, GRH for $\chi_4$ holds if and only if, for every $\alpha>\frac12$, $\xi_K(\alpha)/\xi_K(\alpha+\sqrt u)$ is the Laplace transform of $H^\xi_\alpha+H$ for some GGC $H$ independent of $H^\xi_\alpha$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `cor:dbeta:dom`.

<a id="open_coupling_construction_over_gaussian_primes"></a>

### Open question on a coupling construction over the Gaussian primes

Question 8.5.29, page 333 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function, $\xi(s,\chi_4)$ the completed Dirichlet beta function, $\xi_K=\xi\cdot\xi(\cdot,\chi_4)$ the completed Dedekind zeta function of $K=\mathbb{Q}(i)$, and for $\alpha\ge1$ let $H^\xi_\alpha$ and $H^{(4)}_\alpha$ be independent generalised gamma convolutions with Laplace transforms $\xi(\alpha)/\xi(\alpha+\sqrt u)$ and $\xi(\alpha,\chi_4)/\xi(\alpha+\sqrt u,\chi_4)$, so that $H^K_\alpha=H^\xi_\alpha+H^{(4)}_\alpha$ has Laplace transform $\xi_K(\alpha)/\xi_K(\alpha+\sqrt u)$. Is there a construction of $H^K_\alpha$ from the Gaussian primes themselves, two of norm $p$ for each prime $p\equiv1\pmod4$ and one of norm $p^2$ for each prime $p\equiv3\pmod4$, that persists below the basepoint one, that is for $\frac12<\alpha<1$? Such a construction must see the phases $p^{-i\theta}$ of the primes, not only the sizes of the weights, and below $\alpha=1$ its existence is equivalent, under the Riemann hypothesis, to GRH for $\chi_4$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:dbeta:coupling`.

<a id="thorin_mass_of_dirichlet_beta_as_argument_of_completed_function"></a>

### Thorin mass of the Dirichlet beta function as an argument

Proposition 8.5.30, page 333 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $L(s,\chi_4)=\sum_{n\ge0}(-1)^n(2n+1)^{-s}$ be the Dirichlet beta function, $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ its completed function and $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta,\chi_4)$. Let $\alpha>\frac12$ be such that no zero of $\xi(s,\chi_4)$ lies on the line $\operatorname{Re}s=\alpha$; this holds for every $\alpha\ge1$. Then for $\theta>0$, $\frac1\pi\int_0^\theta v_\alpha(t)\,dt=\frac1\pi\arg\xi(\alpha+i\theta,\chi_4)$, with the argument continuous along the segment from $\alpha$ to $\alpha+i\theta$ and zero at $\alpha$; for $\alpha\ge1$ the left side is the mass of $(0,\theta^2]$ under the Thorin measure of $u\mapsto\xi(\alpha,\chi_4)/\xi(\alpha+\sqrt u,\chi_4)$. As $\theta\to\infty$, with $a=\frac{\alpha+1}2$, $\frac1\pi\arg\xi(\alpha+i\theta,\chi_4)=\frac\theta{2\pi}\log\frac{2\theta}{\pi e}+\frac\alpha4+\frac1\pi\arg L(\alpha+i\theta,\chi_4)-\frac{a(a-1)+\frac16}{\pi\theta}+O(\theta^{-3})$, with $\arg L$ continuous along the same segment. As $\alpha\downarrow\frac12$, for $\theta$ not an ordinate, the left side tends to the number of zeros on the critical line with $0<\gamma<\theta$, which is $N(\theta)$ under GRH; at $\alpha=\frac12$ the expansion is the counting formula $N(T)=\frac T{2\pi}\log\frac{2T}{\pi e}+\frac18+\frac1\pi\arg L(\frac12+iT,\chi_4)+\frac1{48\pi T}+O(T^{-3})$, where $N(T)$ counts the zeros of $L(s,\chi_4)$ with $0<\operatorname{Im}\rho\le T$, for $T$ not an ordinate, with $\arg L(\frac12+iT,\chi_4)$ continued from $L(2,\chi_4)>0$ along $2\to2+iT\to\frac12+iT$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:arg`.

<a id="grh_chi4_iff_argument_of_completed_function_nondecreasing"></a>

### GRH for chi4 and monotonicity of the argument of the completed function

Corollary 8.5.31, page 335 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ be the completed Dirichlet beta function. GRH for $\chi_4$ holds if and only if, for every $\alpha>\frac12$, $\theta\mapsto\arg\xi(\alpha+i\theta,\chi_4)$ (continued along the vertical line from the real point $\alpha$) is nondecreasing. With no hypothesis, the heights $\theta>0$ at which it decreases lie in the union of the windows $|\theta-\gamma|<\sqrt{\delta^2-(\alpha-\frac12)^2}$, each of length less than $1$, over the zeros $\frac12+\delta+i\gamma$ of $L(s,\chi_4)$ with $\delta>\alpha-\frac12$, and there are none at heights below $2.5\times10^7-\frac12$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `cor:dbeta:argmono`.

<a id="half_kernel_of_dirichlet_beta_not_a_gamma_mixture"></a>

### The half kernel of the Dirichlet beta function is not a gamma mixture

Proposition 8.5.32, page 336 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K_4(u)=2e^{3u/2}\eta(ie^{2u})^3$ be the kernel of the completed Dirichlet beta function, $\xi(\frac12+s,\chi_4)=\int_{\mathbb{R}}K_4(u)e^{su}\,du$, with $\eta$ Dedekind's eta function. For no $a>0$ is $u^{1-a}K_4(u)$ completely monotone on $u>0$; in particular the half kernel $K_4(u)\mathbf 1\{u>0\}$ is not a positive mixture of Gamma($a$) densities $u^{a-1}e^{-\lambda u}$, for any shape $a$. The theta function $\eta(ix)^3$ is not completely monotone on $x>0$.

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `prop:dbeta:nocm`.

<a id="open_reciprocal_side_representation_from_eta_product"></a>

### Open question on a reciprocal side representation from the eta product

Question 8.5.33, page 337 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $q=e^{-2\pi x}$ and let $\eta$ be Dedekind's eta function. The function $e^{-\pi x/4}/\eta(ix)^3=\prod_{m\ge1}(1-q^m)^{-3}=\exp\Big(3\sum_{N\ge1}\sigma_{-1}(N)e^{-2\pi Nx}\Big)$ of $x>0$, with $\sigma_{-1}(N)=\sum_{d\mid N}1/d$, is completely monotone, since its logarithm is, and so is the analogue $\prod_{m\ge1}(1-q^m)^{-1}$ for the kernel of $L(s,\chi_{12})$. Can this complete monotonicity be turned into a positive representation on the reciprocal side, in the spirit of Thorin's construction for the reciprocal of a GGC transform?

Source: [`ch/ch16.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch16.tex), label `q:dbeta:recip`.

<a id="laplace_transform_of_bessel_hitting_times_with_curvature"></a>

### Laplace transform of Bessel hitting times with curvature

Proposition 8.7.6, page 348 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. For real $\kappa$ let $s_\kappa(z)=\sin(\sqrt\kappa z)/\sqrt\kappa$ (so that $s_{-1}=\sinh$, $s_0(z)=z$ and $s_1=\sin$), and let $J^{(\kappa)}_a$ be the hitting time of the radius $a$ from the centre by the radial part of Brownian motion on the three dimensional space form of curvature $\kappa$, the diffusion with generator $\frac12\frac{d^2}{dx^2}+\sqrt\kappa\cot(\sqrt\kappa x)\frac d{dx}$, with $a<\pi/\sqrt\kappa$ when $\kappa>0$. For $\kappa=0$ this is the hitting time $J_a$ of $a$ by the three dimensional Bessel process started at $0$. For $u\ge0$, $\mathbb{E}e^{-uJ^{(\kappa)}_a}=\frac{\sqrt{2u-\kappa}\;s_\kappa(a)}{\sinh\big(a\sqrt{2u-\kappa}\big)}$ (read through $\sinh(ix)/i=\sin x$ when $\kappa>2u$), the Laplace transform of a GGC with Thorin measure $\sum_{n\ge1}\delta_{\pi^2n^2/2a^2-\kappa/2}$. Equivalently, the law of $J^{(\kappa)}_a$ is the exponential tilt $e^{\kappa t/2}\Pr(J_a\in dt)/\mathbb{E}e^{\kappa J_a/2}$ of the law of $J_a$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `prop:biane:skappa`.

<a id="exact_quantisation_rule_for_zeros_of_polya_bessel_function"></a>

### Exact quantisation rule for the zeros of Pólya's Bessel function

Proposition 8.7.8, page 352 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K_\mu$ and $I_\mu$ be the modified Bessel functions of the second and first kind, $P(\theta)=K_{i\theta/2}(2\pi)$, $\nu_1<\nu_2<\cdots$ the positive zeros of $\nu\mapsto K_{i\nu}(2\pi)$ (so that the zeros of $P$ are $\pm2\nu_n$), $N_P(T)$ the number of zeros of $P$ in $(0,T]$ and $N_0(T)=\frac T{2\pi}\log\frac T{2\pi e}$. For real $\nu$ put $\vartheta_P(\nu)=\operatorname{Im}\log\Gamma(1+i\nu)-\nu\log\pi-\arg{}_0F_1(;1+i\nu;\pi^2)$, with ${}_0F_1(;b;z)=\sum_{k\ge0}z^k/((b)_kk!)$ and the branches continuous in $\nu$ and equal to zero at $\nu=0$, where ${}_0F_1(;1;\pi^2)=I_0(2\pi)>0$, and $x(T)=\frac1\pi\vartheta_P(T/2)$. Then (i) for $\nu>0$, $I_{i\nu}(2\pi)\ne0$ and $K_{i\nu}(2\pi)=\frac{\pi\,|I_{i\nu}(2\pi)|\,\sin\vartheta_P(\nu)}{\sinh\pi\nu}$. (ii) $\vartheta_P(\nu)>0$ for $\nu>0$, $\vartheta_P(\nu_n)=n\pi$ for every $n$, and the number of zeros of $\nu\mapsto K_{i\nu}(2\pi)$ in $(0,\nu)$ is $\lceil\vartheta_P(\nu)/\pi\rceil-1$; equivalently $N_P(T)=\lfloor x(T)\rfloor$ for every $T>0$ that is not a zero of $P$. (iii) As $\nu\to\infty$, $\vartheta_P(\nu)=\nu\log\frac\nu\pi-\nu+\frac\pi4+\frac{\pi^2-\frac1{12}}\nu+O(\nu^{-3})$, so that $x(T)=N_0(T)+\frac14+\big(2\pi-\frac1{6\pi}\big)\frac1T+O(T^{-3})$ as $T\to\infty$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `prop:biane:phase`.

<a id="index_side_thorin_measure_of_polya_bessel_function"></a>

### Index side Thorin measure of Pólya's Bessel function

Theorem 8.7.9, page 355 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K_\mu(x)=\int_0^\infty e^{-x\cosh t}\cosh(\mu t)\,dt$ be the Macdonald function and $\nu_1<\nu_2<\cdots$ the positive zeros of $\nu\mapsto K_{i\nu}(2\pi)$ (all zeros of $\mu\mapsto K_\mu(2\pi)$ are purely imaginary and simple). Let $H$ be the first passage time to the level $\sqrt{2\pi}$ of Brownian motion with drift $\sqrt{2\pi}$, with density $\tau^{-3/2}e^{2\pi}e^{-\pi(\tau+1/\tau)}$ on $(0,\infty)$ and $\mathbb{E}H^s=2e^{2\pi}K_{s-1/2}(2\pi)$. Then for $u\ge0$, $\frac{K_0(2\pi)}{K_{\sqrt u}(2\pi)}=\frac{\mathbb{E}H^{1/2}}{\mathbb{E}H^{1/2+\sqrt u}}=\prod_{n\ge1}\Big(1+\frac u{\nu_n^2}\Big)^{-1}$. This is the Laplace transform of the GGC $T_P=\sum_{n\ge1}e_n/\nu_n^2$, with independent standard exponential variables $e_n$, no drift and Thorin measure $\sum_n\delta_{\nu_n^2}$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `thm:biane:indexP`.

<a id="clock_of_polya_bessel_function_hitting_time_of_explicit_diffusion"></a>

### Clock of Pólya's Bessel function as a hitting time of an explicit diffusion

Theorem 8.7.10, page 357 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K_0$, $K_1$ and $K_\mu$ be Macdonald functions and $R$ the diffusion on $(0,\infty)$ with generator $\mathcal{G}=\frac12\frac{d^2}{dr^2}+\Big(\frac1{2r}-\frac{K_1(r)}{K_0(r)}\Big)\frac d{dr}$, the radial part of planar Brownian motion killed at rate $\frac12$ and Doob transformed by $K_0$. For $r>b>0$ let $T_b$ be the hitting time of $b$ and $\tau_r=\frac12\int_0^{T_b}dt/R_t^2$ under the law of $R$ started at $r$. Then $T_b<\infty$ almost surely and, for $u\ge0$, $\mathbb{E}e^{-u\tau_r}=\frac{K_{\sqrt u}(r)}{K_0(r)}\cdot\frac{K_0(b)}{K_{\sqrt u}(b)}$. As $r\to\infty$, $\tau_r$ converges in law to $\tau_\infty$ with $\mathbb{E}e^{-u\tau_\infty}=K_0(b)/K_{\sqrt u}(b)$. For $b=2\pi$, $\tau_\infty\overset{d}{=}T_P$, the GGC $\sum_{n\ge1}e_n/\nu_n^2$ with independent standard exponential variables $e_n$ and Thorin measure $\sum_n\delta_{\nu_n^2}$, where $\nu_n$ are the positive zeros of $\nu\mapsto K_{i\nu}(2\pi)$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `thm:biane:clock`.

<a id="sine_squared_forms_of_bessel_clock_drifted_clock_and_polya_clock"></a>

### Sine squared forms of the Bessel clock, the drifted clock and Pólya's clock

Proposition 8.7.13, page 359 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K_\mu$ be the Macdonald function, $\nu_n$ the positive zeros of $\nu\mapsto K_{i\nu}(2\pi)$ and $\Theta_P(v)=\sum_{n\ge1}e^{-\nu_nv}=\operatorname{tr}e^{-v\sqrt{A_P}}$, where $A_P=-\frac{d^2}{dx^2}+e^{2x}$ acts in $L^2([\log2\pi,\infty))$ with a Dirichlet condition at $\log2\pi$ and has the eigenvalues $\nu_n^2$. For real $\lambda$ and $a,c,\ell>0$, $2\log\frac{\sinh\lambda a}{\lambda a}=\int_0^\infty(1-\cos\lambda v)\,\frac{4\,dv}{v\,(e^{\pi v/a}-1)}$, $\ell\big(\sqrt{c^2+\lambda^2}-c\big)=\int_0^\infty(1-\cos\lambda v)\,\frac{2\ell c}\pi K_1(cv)\,\frac{dv}v$, and $\log\frac{K_\lambda(2\pi)}{K_0(2\pi)}=\int_0^\infty4\sin^2\Big(\frac{\lambda v}2\Big)\Theta_P(v)\,\frac{dv}v$. The left sides are $-\log\mathbb{E}e^{-\lambda^2W/2}$ for $W=J_a+J'_a$, the sum of two independent hitting times of the level $a$ by the three dimensional Bessel process started at $0$, and for $W=H^{(c)}_\ell$, the first passage time to the level $\ell$ of Brownian motion with drift $c$, and $-\log\mathbb{E}e^{-\lambda^2T_P}$ for the GGC $T_P=\sum_ne_n/\nu_n^2$ with independent standard exponential variables $e_n$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `prop:biane:sin2`.

<a id="klein_gordon_criterion_for_real_zeros"></a>

### Klein Gordon criterion for real zeros

Theorem 8.7.14, page 361 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $g$ be even, positive and continuous on $\mathbb{R}$ with $\int g(w)e^{c|w|}\,dw<\infty$ for every $c$, and $G(t)=\int g(w)e^{itw}\,dw$. Suppose that there are $x_0\in\mathbb{R}$, a real continuous function $V$ on $[x_0,\infty)$ bounded below, and $F\in C^2([x_0,\infty)\times\mathbb{R})$ such that (a) $(\partial_x^2-\partial_w^2-V(x))F=0$ and $F(x_0,\cdot)=g$; (b) $F$ and its partial derivatives of order at most two are $O(e^{-c|w|})$ for every $c>0$, uniformly for $x$ in compact subsets of $[x_0,\infty)$, and $\Psi_\mu(x)=\int F(x,w)e^{\mu w}\,dw$ lies in $L^2((x_0,\infty))$ for every $\mu\in\mathbb{C}$; (c) $\Psi_\mu$ is not identically zero for any $\mu\in\mathbb{C}$. Then every zero of $G$ is real, and its square is an eigenvalue of $-\frac{d^2}{dx^2}+V$ on $[x_0,\infty)$ with a Dirichlet condition at $x_0$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `thm:biane:kg`.

<a id="de_bruijn_newman_constant_of_polya_function_is_finite"></a>

### De Bruijn Newman constant of Pólya's function is finite

Proposition 8.7.15, page 363 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. For real $t$ let $P_t(\theta)=\int_{\mathbb{R}}e^{tu^2/4}e^{-2\pi\cosh2u}\,e^{i\theta u}\,du$, so that $P_0(\theta)=K_{i\theta/2}(2\pi)$ is Pólya's Macdonald function. There is a constant $\Lambda_P$ with $-\infty<\Lambda_P\le0$ such that $P_t$ has only real zeros if and only if $t\ge\Lambda_P$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `prop:biane:lambdaP`.

<a id="gaussian_damping_property_cannot_imply_real_zeros"></a>

### A property that survives Gaussian damping cannot imply real zeros

Proposition 8.7.17, page 368 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi(u)=2\sum_{n\ge1}\big(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\big)\exp(-\pi n^2e^{2u})$ be the Pólya kernel of Riemann's xi function, with $\xi(\frac12+i\theta)=\int_{\mathbb{R}}\Phi(u)\cos(\theta u)\,du$, and $\Phi^*(u)=8\pi^2\cosh\frac{9u}2\,e^{-2\pi\cosh2u}$ Pólya's asymptotic replacement. Let $\mathcal{P}$ be a property of even, positive, integrable kernels such that (i) if $\mathcal{P}(g)$ holds, then $\int g(u)e^{i\theta u}\,du$ has only real zeros, and (ii) there is $\varepsilon>0$ such that $\mathcal{P}(g)$ implies $\mathcal{P}(e^{-\varepsilon u^2/4}g)$. (a) If (ii) holds for $g=\Phi$, then $\mathcal{P}(\Phi)$ is false. (b) If (ii) holds for all $g$ with the same $\varepsilon$, then $\mathcal{P}(g)$ is false for every kernel $g$ that is $O(e^{-|u|^b})$ for some $b>2$ and satisfies $g(u)e^{Cu^4}\to0$ as $|u|\to\infty$ for every $C$; in particular $\mathcal{P}$ fails for $\Phi$, for $\Phi^*$ and for $e^{-2\pi\cosh2u}$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `prop:biane:damping`.

<a id="open_diffusion_for_clock_of_riemann_xi_function"></a>

### Open question on a diffusion for the clock of the Riemann xi function

Question 8.7.18, page 370 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$ and let $T$ be the clock of $\xi$, the random time with Laplace transform $\mathbb{E}e^{-uT}=\xi(\frac12)/\xi(\frac12+\sqrt u)$ (it exists without any hypothesis). Is there a diffusion on a half line, with an entrance boundary, whose hitting time of a level from the entrance boundary has the Laplace transform $\xi(\frac12)/\xi(\frac12+2\sqrt u)$, that of $4T$, as the clock $\tau_\infty$ of Pólya's Bessel function has the transform $K_0(2\pi)/K_{\sqrt u}(2\pi)$? By Kent's theorem such a hitting time is $\sum_ne_n/\lambda_n$ over the Dirichlet eigenvalues $\lambda_n>0$ of the generator, so a positive answer implies the Riemann hypothesis; under the Riemann hypothesis and the simplicity of the zeros, Krein's inverse spectral theory gives a generalised diffusion with the eigenvalues $\gamma^2/4$, $\frac12+i\gamma$ running over the zeros with $\gamma>0$. The question asks for one constructed without the Riemann hypothesis, for example as a perturbation of the diffusion for Pólya's function, with generator $\frac{d^2}{dx^2}+2\frac{\psi'}\psi\frac d{dx}$ and $\psi(x)=K_0(e^x)$, whose potential $e^{2x}$ has the right Weyl law.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `q:biane:radial`.

<a id="open_klein_gordon_extension_of_polya_kernel"></a>

### Open question on a Klein Gordon extension of Pólya's kernel

Question 8.7.19, page 370 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\Phi(u)=2\sum_{n\ge1}\big(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\big)\exp(-\pi n^2e^{2u})$ be the Pólya kernel of Riemann's xi function and $g(w)=\Phi(w/2)$. Does $g$ admit a Klein Gordon extension, that is: a real continuous potential $V$ on some $[x_0,\infty)$, bounded below and asymptotic to $e^{2x}$ at $+\infty$, and a function $F\in C^2([x_0,\infty)\times\mathbb{R})$ with $(\partial_x^2-\partial_w^2-V(x))F=0$ and $F(x_0,\cdot)=g$, such that $F$ and its partial derivatives of order at most two are $O(e^{-c|w|})$ for every $c>0$, uniformly for $x$ in compact sets, and $\Psi_\mu(x)=\int F(x,w)e^{\mu w}\,dw$ lies in $L^2((x_0,\infty))$ for every complex $\mu$ and is not identically zero for any $\mu$? By the Klein Gordon criterion a positive answer implies the Riemann hypothesis. Under the Riemann hypothesis the question asks whether $\mu\mapsto2\xi(\frac12+2\mu)$ is, up to normalisation, the Jost function of such a potential whose Klein Gordon solution has the boundary value $g$.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `q:biane:kg`.

<a id="open_de_bruijn_newman_constants_of_polya_functions"></a>

### Open question on the de Bruijn Newman constants of Pólya's functions

Question 8.7.20, page 370 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For real $t$ let $P_t(\theta)=\int_{\mathbb{R}}e^{tu^2/4}e^{-2\pi\cosh2u}\,e^{i\theta u}\,du$ (so that $P_0(\theta)=K_{i\theta/2}(2\pi)$ is Pólya's Macdonald function) and $\Xi^*_t(\theta)=\int_{\mathbb{R}}e^{tu^2/4}\Phi^*(u)e^{i\theta u}\,du$ with Pólya's kernel $\Phi^*(u)=8\pi^2\cosh\frac{9u}2\,e^{-2\pi\cosh2u}$. The de Bruijn Newman constant $\Lambda_P$ of $P$, the least $t$ from which all zeros of $P_t$ are real, satisfies $-\infty<\Lambda_P\le0$, and the same holds for the constant $\Lambda_{\Xi^*}$ of $\Xi^*$. What are $\Lambda_P$ and $\Lambda_{\Xi^*}$? At what height do the first nonreal zeros of $P_t$ appear as $t$ decreases through $\Lambda_P$, and do they appear, as for the de Bruijn Newman flow of Riemann's $\Xi$, through the collision of two real zeros?

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `q:biane:lambda`.

<a id="open_family_joining_the_two_clocks"></a>

### Open question on a family joining the two clocks

Question 8.7.21, page 370 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $H$ be the first passage time of Brownian motion with drift $\sqrt{2\pi}$ to the level $\sqrt{2\pi}$, whose Mellin transform is $\mathbb{E}H^s=2e^{2\pi}K_{s-1/2}(2\pi)$ and for which $1/H$ has the law of $H$ size biased by $H$, a reciprocal symmetry that encodes the functional equation $\mathbb{E}H^s=\mathbb{E}H^{1-s}$. Let $S_2=J_1+J'_1$ be the sum of two independent hitting times of the level $1$ by the three dimensional Bessel process started at $0$, for which $2\xi(s)=\mathbb{E}\big[(\pi S_2/2)^{s/2}\big]$ with Riemann's xi function $\xi$, the functional equation coming from Jacobi's theta transformation. Is there a one parameter family of laws joining the drifted first passage time and the sum of two Bessel hitting times, each with a reciprocal symmetry that encodes a functional equation, along which the index side Thorin measure can be followed? The curvature family of Brownian hitting times, an exponential tilt of the Bessel hitting time, breaks the symmetry. Positivity of the index side Thorin measure at the Bessel end is the Riemann hypothesis.

Source: [`ch/biane.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/biane.tex), label `q:biane:curv`.

## Chapter 9. The basepoint one

<a id="compensated_sine_squared_formula_basepoint_one"></a>

### Compensated sine squared formula at the basepoint one

Proposition 9.1.2, page 380 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, with zeros $\rho$, and let $v_1(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(1+i\theta)$. Let $\Lambda$ be the von Mangoldt function, $\psi=\Gamma'/\Gamma$ the digamma function, $\gamma_E$ Euler's constant, and $\nu_1=\sum_n\Lambda(n)n^{-1}\delta_{\log n}-\mathbf 1_{(0,\infty)}(x)\,dx$ the prime measure compensated by Lebesgue measure. For every real $\theta\ne0$, \[v_1(\theta)=\frac1{1+\theta^2}-\frac12\log\pi+\frac12\operatorname{Re}\psi\Big(\frac{1+i\theta}2\Big)+\gamma_E+\lim_{X\to\infty}\int_0^{\log X}2\sin^2\Big(\frac{\theta x}2\Big)\nu_1(dx),\] that is, the limit is that of $\sum_{n\le X}\Lambda(n)n^{-1}\,2\sin^2(\tfrac12\theta\log n)-\big(\log X-\sin(\theta\log X)/\theta\big)$. The limit exists by the prime number theorem, and at $\theta=0$ the formula reduces to $v_1(0)=\sum_\rho\rho^{-1}$ (the sum over the nontrivial zeros taken in conjugate pairs).

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:sin2`.

<a id="levy_form_of_xi_at_basepoint_one"></a>

### Lévy form of the completed zeta function at the basepoint one

Proposition 9.1.3, page 382 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\psi=\Gamma'/\Gamma$ the digamma function, $\gamma_E$ Euler's constant and $\Lambda$ the von Mangoldt function. For $\operatorname{Re}s>0$, \[\log\frac{\xi(1+s)}{\xi(1)}=b_1s+\int_0^\infty\big(e^{-sx}-1+sx\big)\,\Pi_1(dx),\] with drift $b_1=-\frac12\log\pi+\frac12\psi(\frac32)+\gamma_E=\sum_\rho\rho^{-1}=0.0230957089661$ (the sum over the nontrivial zeros $\rho$ of $\zeta$ in conjugate pairs) and the signed Lévy measure \[\Pi_1(dx)=\frac{e^{-3x}}{1-e^{-2x}}\,\frac{dx}x-\frac{dx}x+\sum_{n\ge2}\frac{\Lambda(n)}{n\log n}\,\delta_{\log n}(dx),\] the integral against the last two pieces (the pole term $-dx/x$ and the prime atoms) being taken as $\lim_{X\to\infty}\int_0^{\log X}$.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:levy`.

<a id="thorin_density_minimum_at_origin_basepoint_one"></a>

### Minimum of the Thorin density at the origin and stable summand at the basepoint one

Proposition 9.1.4, page 385 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function with zeros $\rho$, $v_1(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(1+i\theta)$ and $b_1=v_1(0)=\frac{\xi'}{\xi}(1)=\sum_\rho\rho^{-1}=0.0230957089661$. Let $T_1$ be the generalized gamma convolution (GGC) with $\mathbb{E}e^{-uT_1}=\xi(1)/\xi(1+\sqrt u)$ for $u\ge0$, whose Thorin density in the variable $\theta=\sqrt z$ is $v_1(\theta)/\pi$. With no hypothesis on the zeros, $v_1(\theta)>v_1(0)=b_1$ for every real $\theta\ne0$; the proof uses a computation of the zeros below height $2000$ in ball arithmetic and, above that height, explicit bounds for the argument of the zeta function on the critical line (Hasanalizade, Shen and Wong 2022, or the classical bounds of Backlund 1918). Consequently $T_1\stackrel{d}{=}S+T_1^\circ$, where $S$ and $T_1^\circ$ are independent, $S$ is positive stable of index one half with $\mathbb{E}e^{-uS}=e^{-b_1\sqrt u}$, and $T_1^\circ$ is a GGC with no drift whose Thorin density in the chart $\theta=\sqrt z$ is $(v_1-b_1)/\pi$.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:v1min`.

<a id="thorin_density_minimum_at_origin_every_basepoint_from_one"></a>

### Minimum at the origin and stable summand at every basepoint from one

Proposition 9.1.5, page 388 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, and for real $\alpha$ and $\theta$ let $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta)$ and $b_\alpha=v_\alpha(0)=\frac{\xi'}{\xi}(\alpha)$. With no hypothesis on the zeros, $v_\alpha(\theta)>v_\alpha(0)=b_\alpha$ for every $\alpha\ge1$ and every real $\theta\ne0$; the proof uses the computation of the zeros below height $2000$ in ball arithmetic of the case $\alpha=1$, above that height explicit bounds for the argument of the zeta function on the critical line (Hasanalizade, Shen and Wong 2022), and a further computation in ball arithmetic. Consequently, for every $\alpha\ge1$, the generalized gamma convolution (GGC) $T_\alpha$ with $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt u)$ satisfies $T_\alpha\stackrel{d}{=}S_\alpha+T_\alpha^\circ$, where $S_\alpha$ and $T_\alpha^\circ$ are independent, $S_\alpha$ is positive stable of index one half with $\mathbb{E}e^{-uS_\alpha}=e^{-b_\alpha\sqrt u}$, which is the law of the first passage time of Brownian motion to the level $b_\alpha/\sqrt2$, and $T_\alpha^\circ$ is a GGC with no drift, with Thorin density $(v_\alpha-b_\alpha)/\pi$ in the chart $\theta=\sqrt z$ and mean $\mathbb{E}T_\alpha^\circ=\frac12(\log\xi)''(\alpha)$.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:valpha`.

<a id="levy_form_below_basepoint_one_and_zero_free_half_plane"></a>

### Lévy form below the basepoint one and zero free half planes

Proposition 9.1.6, page 390 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $\Lambda$ the von Mangoldt function, $\psi=\Gamma'/\Gamma$ the digamma function and $\Theta=\sup\{\operatorname{Re}\rho\}$ over the nontrivial zeros $\rho$ of $\zeta$. Let $\frac12<\alpha<1$, $b_\alpha=\frac{\xi'}{\xi}(\alpha)$ and $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(\alpha+i\theta)$, and let $\Pi_\alpha$ be the signed measure on $(0,\infty)$ with \[x\,\Pi_\alpha(dx)=\frac{e^{-(\alpha+2)x}}{1-e^{-2x}}\,dx+\Big(\sum_n\frac{\Lambda(n)}{n^{\alpha}}\,\delta_{\log n}(dx)-e^{(1-\alpha)x}\,dx\Big),\] the integral against the part in parentheses being taken as $\lim_{X\to\infty}\int_0^{\log X}$. (1) If $\Theta<\alpha$, the Lévy form $\log\frac{\xi(\alpha+s)}{\xi(\alpha)}=b_\alpha s+\int_0^\infty(e^{-sx}-1+sx)\,\Pi_\alpha(dx)$ holds for $\operatorname{Re}s>0$, the compensated limit converging, and the identity $e^{-sx}+sx-1=\frac{2x}\pi\int_0^\infty\log(1+s^2/\theta^2)\sin^2(\frac{x\theta}2)\,d\theta$ converts it into \[v_\alpha(\theta)=b_\alpha+A_\alpha(\theta)+\lim_{X\to\infty}\Big[\sum_{n\le X}\frac{\Lambda(n)}{n^\alpha}\,2\sin^2\Big(\frac{\theta\log n}2\Big)-\int_0^{\log X}2\sin^2\Big(\frac{\theta x}2\Big)e^{(1-\alpha)x}\,dx\Big],\] where $A_\alpha(\theta)=\frac12\big[\operatorname{Re}\psi\big(\frac{\alpha+2+i\theta}2\big)-\psi\big(\frac{\alpha+2}2\big)\big]$; this value of $v_\alpha(\theta)$ is positive. The error of the compensated sum with prime cutoff $X$ (the bracket truncated at $X$ minus its limit) is $O(X^{\Theta-\alpha+\delta})$ for every $\delta>0$. (2) Conversely, if the Lévy form holds for $\operatorname{Re}s>0$, then $\xi(\alpha+s)$ has no zeros with $\operatorname{Re}s>0$, so $\Theta\le\alpha$.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:below`.

<a id="thorin_density_ratio_comparison_across_vertical_lines"></a>

### Comparison of the Thorin density across vertical lines

Theorem 9.1.8, page 395 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be one of $\xi$, $\xi(\cdot,\chi_4)$ and $\xi(\cdot,\chi_{12})$, where $\xi$ is the completed Riemann zeta function and $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ and $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(\frac s2)L(s,\chi_{12})$ are the completed Dirichlet $L$ functions of the nonprincipal character modulo $4$ and of the real primitive character modulo $12$; each $f$ is entire of order one, real on the real axis, with $f(s)=f(1-s)$ and no real zeros. Write $v_\alpha(\theta)=\operatorname{Re}\frac{f'}{f}(\alpha+i\theta)$. Let $\frac12<\alpha<\alpha'$ and let $\theta$ be real. Suppose that every zero $\frac12+\delta+i\gamma$ of $f$ with $\delta>(\alpha-\frac12)/\sqrt2$ satisfies $|\theta-\gamma|\ge(1+\sqrt2)\delta$. Then \[\frac{v_\alpha(\theta)}{\alpha-\frac12}\ \ge\ \frac{v_{\alpha'}(\theta)}{\alpha'-\frac12}.\] The hypothesis holds at every $\theta$ under the Riemann Hypothesis (RH) for $\xi$, or the Generalized Riemann Hypothesis (GRH) for $\chi_4$ or $\chi_{12}$; at every $\theta$ when $\alpha\ge\frac12+\frac1{\sqrt2}$, with no hypothesis; and at every $\theta$ with $|\theta|\le H_0-\frac{1+\sqrt2}2$, where $H_0$ is a height below which every zero of $f$ is known to lie on the critical line: $H_0=3\times10^{12}$ for $\xi$ (Platt and Trudgian 2021), and $H_0=10^8/q$ for $\xi(\cdot,\chi_q)$, that is $2.5\times10^7$ for $\chi_4$ and $10^8/12>8.33\times10^6$ for $\chi_{12}$ (Platt 2016).

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `thm:dbeta:comp`.

<a id="explicit_lower_bound_of_thorin_density_by_comparison_with_euler_product"></a>

### Explicit lower bound of the Thorin density by comparison with the Euler product

Corollary 9.1.9, page 396 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be one of $\xi$, $\xi(\cdot,\chi_4)$ and $\xi(\cdot,\chi_{12})$, where $\xi$ is the completed Riemann zeta function and $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ and $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(\frac s2)L(s,\chi_{12})$ are the completed Dirichlet $L$ functions of the nonprincipal character modulo $4$ and of the real primitive character modulo $12$; each $f$ is entire of order one, real on the real axis, with $f(s)=f(1-s)$ and no real zeros. Write $v_\alpha(\theta)=\operatorname{Re}\frac{f'}{f}(\alpha+i\theta)$, let $\psi=\Gamma'/\Gamma$ be the digamma function and $\Lambda$ the von Mangoldt function, and let $q=1$, $4$, $12$ for the three functions in the order above. For $\sigma>1$ and real $\theta$, $v_\sigma(\theta)\ge A^\infty_f(\sigma,\theta)-S_f(\sigma)$, where \[A^\infty_\xi(\sigma,\theta)=-\tfrac12\log\pi+\tfrac12\operatorname{Re}\psi\big(\tfrac{\sigma+i\theta}2\big),\quad A^\infty_{\chi_4}(\sigma,\theta)=\tfrac12\log\tfrac4\pi+\tfrac12\operatorname{Re}\psi\big(\tfrac{\sigma+1+i\theta}2\big),\quad A^\infty_{\chi_{12}}(\sigma,\theta)=\tfrac12\log\tfrac{12}\pi+\tfrac12\operatorname{Re}\psi\big(\tfrac{\sigma+i\theta}2\big),\] and $S_f(\sigma)=\sum_{(n,q)=1}\Lambda(n)n^{-\sigma}$, so that $S_\xi=-\zeta'/\zeta$, $S_{\chi_4}(\sigma)=S_\xi(\sigma)-\log2/(2^\sigma-1)$ and $S_{\chi_{12}}(\sigma)=S_{\chi_4}(\sigma)-\log3/(3^\sigma-1)$; moreover $A^\infty_f(\sigma,\theta)$ increases in $|\theta|$. Hence, if $\frac12<\alpha<\alpha'$ with $\alpha'>1$, $\theta$ is real, and every zero $\frac12+\delta+i\gamma$ of $f$ with $\delta>(\alpha-\frac12)/\sqrt2$ satisfies $|\theta-\gamma|\ge(1+\sqrt2)\delta$ (the comparison hypothesis, which holds at every $\theta$ when all zeros of $f$ lie on the line $\operatorname{Re}s=\frac12$, that is under the Riemann Hypothesis for $\xi$ or the Generalized Riemann Hypothesis for $\chi_4$ or $\chi_{12}$, and at every $\theta$ when $\alpha\ge\frac12+\frac1{\sqrt2}$ with no hypothesis), then \[v_\alpha(\theta)\ \ge\ \frac{\alpha-\frac12}{\alpha'-\frac12}\Big(A^\infty_f(\alpha',\theta)-S_f(\alpha')\Big)\ =\ \frac{\alpha-\frac12}{2\alpha'-1}\log|\theta|+O(1),\] the $O(1)$ being as $|\theta|\to\infty$ with $\alpha$, $\alpha'$ and $f$ fixed.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `cor:dbeta:lower`.

<a id="thorin_density_minimum_at_origin_below_one_equivalent_to_rh"></a>

### Minimum of the Thorin density at the origin and its equivalence with RH

Theorem 9.1.11, page 397 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $f$ be one of $\xi$, $\xi(\cdot,\chi_4)$ and $\xi(\cdot,\chi_{12})$, where $\xi$ is the completed Riemann zeta function and $\xi(s,\chi_4)=(4/\pi)^{(s+1)/2}\Gamma(\frac{s+1}2)L(s,\chi_4)$ and $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(\frac s2)L(s,\chi_{12})$ are the completed Dirichlet $L$ functions of the nonprincipal character modulo $4$ and of the real primitive character modulo $12$; each $f$ is entire of order one, real on the real axis, with $f(s)=f(1-s)$ and no real zeros. Write $v_\alpha(\theta)=\operatorname{Re}\frac{f'}{f}(\alpha+i\theta)$, $b_\alpha=v_\alpha(0)=\frac{f'}{f}(\alpha)$ and $V_f=(\log f)''(\frac12)$, which is $0.04620998623$, $0.1560296499$ and $0.3329236082$ for the three functions in the order above. Let $H_0$ be a height below which every zero of $f$ is known to lie on the critical line, $H_0=3\times10^{12}$ (Platt and Trudgian 2021), $2.5\times10^7$ and $10^8/12$ (Platt 2016) for the three functions, and let $\alpha_{\max}=7.35$, $4.70$ and $3.10$ for the three functions. Then: (i) for every $\alpha>\frac12$, $0<b_\alpha\le(\alpha-\frac12)V_f$, and the ratio $b_\alpha/(\alpha-\frac12)$ decreases in $\alpha$ and tends to $V_f$ as $\alpha\downarrow\frac12$; (ii) for $\frac12<\alpha\le\alpha_{\max}$, $v_\alpha(\theta)>b_\alpha$ for $0<|\theta|\le H_0-\frac{1+\sqrt2}2$; (iii) for $\frac12+\frac1{\sqrt2}\le\alpha\le\alpha_{\max}$, $v_\alpha(\theta)>b_\alpha$ for every real $\theta\ne0$; (iv) for $f=\xi(\cdot,\chi_4)$ and $f=\xi(\cdot,\chi_{12})$, $v_\alpha(\theta)>b_\alpha$ for every $\alpha\ge\frac12+\frac1{\sqrt2}$ and every real $\theta\ne0$; (v) for $f=\xi(\cdot,\chi_4)$ and $f=\xi(\cdot,\chi_{12})$, $v_\alpha(\theta)>b_\alpha$ for every $\alpha\ge1$ and every real $\theta\ne0$; (vi) the Riemann Hypothesis (RH) for $\xi$, respectively the Generalized Riemann Hypothesis (GRH) for $\chi_4$ or for $\chi_{12}$ (all zeros of $f$ on the line $\operatorname{Re}s=\frac12$), holds if and only if $v_\alpha(\theta)>b_\alpha$ for every $\alpha\in(\frac12,1)$ and every real $\theta\ne0$.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `thm:dbeta:minall`.

<a id="heat_trace_at_basepoint_one_positive_certified_from_primes"></a>

### Certified positivity of the heat trace at the basepoint one

Proposition 9.1.17, page 404 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $v_1(\theta)=\operatorname{Re}\frac{\xi'}{\xi}(1+i\theta)$, and $T_1$ the generalized gamma convolution with $\mathbb{E}e^{-uT_1}=\xi(1)/\xi(1+\sqrt u)$. Let $k_1(t)=\pi^{-1}\int_0^\infty e^{-t\theta^2}v_1(\theta)\,d\theta$ for $t>0$ be the heat trace of $T_1$, that is, the Laplace transform in the variable $z=\theta^2$ of its Thorin measure. Then $2\sqrt{\pi t}\,k_1(t)>0.021$ for every $t>0$. The proof uses the primes up to $10^7$, explicit bounds for $\psi(y)-y$, where $\psi(y)=\sum_{n\le y}\Lambda(n)$ is Chebyshev's function and $\Lambda$ the von Mangoldt function (Büthe 2018; Broadbent, Kadiri, Lumley, Ng and Wilk 2021; Johnston and Yang 2023), and a computation in ball arithmetic.

Source: [`ch/basepoint.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/basepoint.tex), label `bp:prop:heatcert`.

<a id="no_brownian_first_passage_to_random_level_at_basepoint_one"></a>

### No Brownian first passage to a random level at the basepoint one

Proposition 9.2.1, page 407 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $T_1$ the generalized gamma convolution with $\mathbb{E}e^{-uT_1}=\xi(1)/\xi(1+\sqrt u)$. Let $B$ be a standard Brownian motion and, for $c\in\mathbb{R}$ and $a\ge0$, let $H^{(c)}_a=\inf\{t:B_t+ct\ge a\}$ be the first passage time of $B_t+ct$ to the level $a$. For no drift $c\ge0$ and no random level $L\ge0$ independent of the Brownian motion is $T_1$ equal in law to the first passage time $H^{(c)}_L$. (For $c<0$ the first passage time is infinite with positive probability, so equality is impossible.)

Source: [`ch/bphit.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bphit.tex), label `bp:prop:nogo`.

<a id="basepoint_one_clock_subordination_of_tilted_centre_clock"></a>

### Subordination of the tilted centre clock at the basepoint one

Theorem 9.2.2, page 408 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $T$ the centre clock, the positive random variable with $\mathbb{E}e^{-uT}=\xi(\frac12)/\xi(\frac12+\sqrt u)$ for $u\ge0$, which exists with no hypothesis on the zeros (its existence uses the verification of the Riemann Hypothesis below height $1000$). Let $\widetilde T$ have the tilted law $P(\widetilde T\in dt)=e^{-t/4}P(T\in dt)/\mathbb{E}e^{-T/4}$, and let $(\sigma_\ell)_{\ell\ge0}$ be a subordinator independent of $\widetilde T$ with Laplace exponent $u+\sqrt u$, that is, $\mathbb{E}e^{-u\sigma_\ell}=e^{-\ell(u+\sqrt u)}$. Then the shifted clock $T_1$, with $\mathbb{E}e^{-uT_1}=\xi(1)/\xi(1+\sqrt u)$, satisfies \[T_1\stackrel{d}{=}\sigma_{\widetilde T}\stackrel{d}{=}\widetilde T+H_{\widetilde T/\sqrt2},\] where $H_a$ is the first passage time of a standard Brownian motion to the level $a$, independent of $\widetilde T$.

Source: [`ch/bphit.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bphit.tex), label `bp:thm:sub`.

<a id="basepoint_one_clock_kent_string_hitting_time_under_rh"></a>

### Kent string hitting time of the shifted clock at the basepoint one under RH

Corollary 9.2.3, page 409 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on the Riemann Hypothesis (every nontrivial zero of the Riemann zeta function has real part one half) and the simplicity of every nontrivial zero.

Statement. Let $T$ be the centre clock, with $\mathbb{E}e^{-uT}=\xi(\frac12)/\xi(\frac12+\sqrt u)$, where $\xi$ is the completed Riemann zeta function, let $\widetilde T$ be the clock tilted by $e^{-t/4}$, $P(\widetilde T\in dt)=e^{-t/4}P(T\in dt)/\mathbb{E}e^{-T/4}$, and let $T_1$ be the shifted clock with $\mathbb{E}e^{-uT_1}=\xi(1)/\xi(1+\sqrt u)$, so that $T_1\stackrel{d}{=}\sigma_{\widetilde T}$ for an independent subordinator $\sigma$ with Laplace exponent $u+\sqrt u$, the inverse local time at zero of a Brownian motion $B$ that is sticky at zero. Assume RH and that all nontrivial zeros are simple. Then $\widetilde T$ is Kent's hitting time for a Krein string whose Dirichlet spectrum is $\{\gamma^2+\frac14\}=\{|\rho|^2\}$ ($\gamma>0$ running over the ordinates of the zeros $\rho=\frac12+i\gamma$), and $T_1$ is the first hitting time of the level $y_*$ by the process $D(\ell_t)$, where $D$ is that string diffusion started at its left end, $y_*$ is the level whose hitting time by $D$ is $\widetilde T$, and $\ell_t$ is the local time at zero of the independent sticky Brownian motion $B$. Equivalently, $T_1$ is the first hitting time of the line $\{y=y_*\}$ by the two dimensional Markov process $(D(\ell_t),B_t)$.

Source: [`ch/bphit.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bphit.tex), label `bp:cor:kent`.

<a id="shifted_clock_subordination_every_basepoint_above_centre"></a>

### Subordination of the centre clock at every basepoint above the centre

Theorem 9.2.5, page 415 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function with zeros $\rho=\beta+i\gamma$, and $T$ the centre clock, with $\mathbb{E}e^{-uT}=\xi(\frac12)/\xi(\frac12+\sqrt u)$. For every $\alpha>\frac12$, with no hypothesis on the zeros (the existence of the centre clock uses the verification of the Riemann Hypothesis below height $1000$), $\varphi_\alpha(u)=\xi(\alpha)/\xi(\alpha+\sqrt u)$ is the Laplace transform of $T_\alpha\stackrel{d}{=}\sigma^{(\alpha)}_{\widetilde T_\alpha}$, where $\widetilde T_\alpha$ is the centre clock tilted by $e^{-(\alpha-\frac12)^2T}$, that is, $P(\widetilde T_\alpha\in dt)=e^{-(\alpha-\frac12)^2t}P(T\in dt)/\mathbb{E}e^{-(\alpha-\frac12)^2T}$, and $\sigma^{(\alpha)}$ is an independent subordinator with exponent $u+2(\alpha-\frac12)\sqrt u$, that is, $T_\alpha\stackrel{d}{=}\widetilde T_\alpha+H_{\sqrt2(\alpha-\frac12)\widetilde T_\alpha}$ with $H_a$ the first passage time of a standard Brownian motion to the level $a$, independent of $\widetilde T_\alpha$. The law has infinite mean and $P(T_\alpha>x)\sim\frac{\xi'}{\xi}(\alpha)\,(\pi x)^{-1/2}$ as $x\to\infty$. By the threshold theorem of Polson (2026), the variable $T_\alpha$ is a generalized gamma convolution (GGC) if and only if no zero has $\beta>\alpha$, in particular for every $\alpha\ge1$, and when every zero has $\beta<\alpha$ its Thorin density in the chart $\theta=\sqrt z$ is $\frac1\pi\sum_\rho\frac{\alpha-\beta}{(\alpha-\beta)^2+(\theta-\gamma)^2}$.

Source: [`ch/bphit.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bphit.tex), label `bp:thm:alpha`.

<a id="shifted_clock_levy_measure_and_self_decomposability"></a>

### Lévy measure and self decomposability of the shifted clock

Corollary 9.2.6, page 416 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function with zeros $\rho$, let $G(u)=\xi(\frac12+\sqrt u)/\xi(\frac12)=\prod_j(1+u/a_j)$ with $a_j=-(\rho_j-\frac12)^2$, one zero $\rho_j$ from each pair $\{\rho,1-\rho\}$, and let $W(x)=\sum_je^{-a_jx}$, so that the centre clock $T$, with $\mathbb{E}e^{-uT}=1/G(u)$, is infinitely divisible with Lévy measure $W(x)\,dx/x$; with no hypothesis on the zeros (using the verification of the Riemann Hypothesis below height $1000$), $W>0$ and $W'<0$ on $(0,\infty)$. Let $\alpha>\frac12$, $a=\alpha-\frac12$, and let $T_\alpha$ be the shifted clock with $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt u)$. Then $T_\alpha$ is infinitely divisible, with no drift and Lévy measure $W_\alpha(x)\,dx/x$, where \[W_\alpha(x)=\frac{ax}{\sqrt\pi}\int_0^x(x-s)^{-3/2}e^{-a^2xs/(x-s)}\,W(s)\,ds=\mathbb{E}\Big[(1+4xV)^{-1/2}\,W\Big(\frac x{1+4xV}\Big)\Big],\] with $V=a^2/(2Z^2)$ and $Z$ a standard normal variable. Hence $W_\alpha>0$ and $W_\alpha'<0$ on $(0,\infty)$, so that $T_\alpha$ is self decomposable; and on every interval $(0,x_0]$ on which $W$ is completely monotone, so is $W_\alpha$.

Source: [`ch/bphit.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bphit.tex), label `bp:cor:levyalpha`.

## Chapter 10. The critical strip and Riesz's criterion

<a id="reciprocal_xi_real_axis_normal_scale_mixture"></a>

### Reciprocal Xi function on the real axis as a normal scale mixture

Proposition 10.1.1, page 423 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function and $b_\alpha=\xi'(\alpha)/\xi(\alpha)$. For $\alpha\ge\frac12$ let $T_\alpha$ be the positive random variable with Laplace transform $\mathbb{E}\,e^{-\lambda T_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, $\lambda\ge0$ (it exists with no hypothesis on the zeros of $\zeta$; $T_{1/2}$ is the centre clock), let $Z$ be a standard normal variable independent of $T_\alpha$, and put $Y_\alpha=\sqrt{2T_\alpha}\,Z$. Then, for every $\alpha\ge\frac12$, the function $u\mapsto\xi(\alpha)/\xi(\alpha+|u|)$ on $\mathbb{R}$ is the characteristic function of $Y_\alpha$, and the law of $Y_\alpha$ is symmetric, unimodal and infinitely divisible, with no hypothesis on the zeros. At the centre its variance is $2\,\mathbb{E}\,T_{1/2}=2\sum_k a_k^{-1}=0.046210$, where $a_k=-(\rho_k-\frac12)^2$ with one zero $\rho_k$ taken from each pair $\{\rho,1-\rho\}$ of nontrivial zeros of $\zeta$. For $\alpha>\frac12$, $\xi(\alpha)/\xi(\alpha+|u|)=1-b_\alpha|u|+O(u^2)$ as $u\to0$ (a corner at the origin), so $Y_\alpha$ lies in the domain of attraction of the Cauchy law, and $P(|Y_\alpha|>y)\sim 2b_\alpha/(\pi y)$ as $y\to\infty$.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `sm:prop:cf`.

<a id="reciprocal_xi_thorin_form_and_eggc_criterion"></a>

### Thorin form of the reciprocal Xi function and the EGGC criterion

Theorem 10.1.2, page 425 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, let $\rho=\beta+i\gamma$ run over the nontrivial zeros of $\zeta$, and for $\alpha\ge\frac12$ let $T_\alpha$ be the positive random variable with Laplace transform $\mathbb{E}\,e^{-\lambda T_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, $\lambda\ge0$, let $Z$ be a standard normal variable independent of $T_\alpha$, $Y_\alpha=\sqrt{2T_\alpha}\,Z$, and $v_\alpha(\theta)=\operatorname{Re}\,(\xi'/\xi)(\alpha+i\theta)$ for real $\theta$. Recall that a symmetric law is an extended generalised gamma convolution (EGGC) if and only if it is the law of $Z\sqrt{2T'}$ with $T'$ a generalised gamma convolution (GGC) independent of $Z$, where a GGC is a law on $[0,\infty)$ with Laplace transform $\exp\{-a\lambda-\int_{(0,\infty)}\log(1+\lambda/z)\,U(dz)\}$ for a drift $a\ge0$ and a Thorin measure $U$. Let $\alpha\ge1$, or more generally let $\alpha>\frac12$ be such that $\operatorname{Re}\rho<\alpha$ for every zero $\rho$. Then $v_\alpha>0$ and $\xi(\alpha)/\xi(\alpha+|u|)=\exp\{-\frac1\pi\int_0^\infty\log(1+u^2/\theta^2)\,v_\alpha(\theta)\,d\theta\}$ for real $u$, so that $T_\alpha$ is a GGC without drift whose Thorin measure is $(v_\alpha(\theta)/\pi)\,d\theta$ under the substitution $z=\theta^2$, and $Y_\alpha$ is a symmetric EGGC. With no hypothesis on the zeros, for $\frac12<\alpha<1$ the law of $Y_\alpha$ is an EGGC if and only if $\zeta$ has no zeros with real part greater than $\alpha$, and $Y_{1/2}$ is an EGGC if and only if the Riemann Hypothesis holds (that centre statement is the earlier rigidity theorem of Polson, 2026, for the centre clock $T_{1/2}$, read through $Y_{1/2}=\sqrt{2T_{1/2}}\,Z$). Under the Riemann Hypothesis, $\xi(\frac12)/\xi(\frac12+|u|)=\prod_{\gamma>0}(1+u^2/\gamma^2)^{-1}$, the product over the ordinates of the zeros $\frac12+i\gamma$ with $\gamma>0$, counted with multiplicity, and $Y_{1/2}$ has the law of $\sum_kL_k/\gamma_k$ with independent standard Laplace variables $L_k$ (density $\frac12e^{-|y|}$).

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `sm:thm:thorin`.

<a id="eggc_property_across_the_critical_strip"></a>

### EGGC property of the reciprocal Xi law across the critical strip

Theorem 10.1.3, page 428 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, let $\rho=\beta+i\gamma$ run over the nontrivial zeros of $\zeta$ (with multiplicity), and for $\alpha\ge\frac12$ let $T_\alpha$ be the positive random variable with $\mathbb{E}\,e^{-\lambda T_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, $\lambda\ge0$, $Z$ a standard normal variable independent of $T_\alpha$, $Y_\alpha=\sqrt{2T_\alpha}\,Z$ and $v_\alpha(\theta)=\operatorname{Re}\,(\xi'/\xi)(\alpha+i\theta)$. Recall that a symmetric law is an extended generalised gamma convolution (EGGC) if and only if it is the law of $Z\sqrt{2T'}$ with $T'$ a generalised gamma convolution independent of $Z$. Let $\frac12<\alpha<1$ and $\varepsilon=\alpha-\frac12$. (1) With no hypothesis on the zeros, $\xi(\alpha)/\xi(\alpha+|u|)$ is the characteristic function of the infinitely divisible, symmetric, unimodal law $Y_\alpha$. If every zero has real part less than $\alpha$, then $\xi(\alpha)/\xi(\alpha+|u|)=\exp\{-\frac1\pi\int_0^\infty\log(1+u^2/\theta^2)\,v_\alpha(\theta)\,d\theta\}$ holds with the positive density $v_\alpha(\theta)/\pi=\frac1\pi\sum_\rho\frac{\alpha-\beta}{(\alpha-\beta)^2+(\theta-\gamma)^2}$. (2) $Y_\alpha$ is an EGGC if and only if no zero has $\beta>\alpha$, that is, if and only if $\zeta$ has no zeros in the half plane $\operatorname{Re}s>\alpha$. (3) Let $T_0\ge2$ be a height such that all zeros $\rho$ with $0<\gamma\le T_0$ lie on the critical line. Then $\xi(\alpha)/\xi(\alpha+|u|)=\prod_{0<\gamma\le T_0}\frac{\varepsilon^2+\gamma^2}{(\varepsilon+|u|)^2+\gamma^2}\times\mathcal{E}_{\alpha,T_0}(u)$, where the finite product, over the ordinates $\gamma$ of the zeros $\frac12+i\gamma$ with $0<\gamma\le T_0$ (with multiplicity), is the characteristic function of an EGGC with the positive Thorin density $\sum_{0<\gamma\le T_0}[P_\varepsilon(\theta-\gamma)+P_\varepsilon(\theta+\gamma)]$, $P_\varepsilon(x)=\varepsilon/(\pi(\varepsilon^2+x^2))$, and the remainder satisfies $|\log\mathcal{E}_{\alpha,T_0}(u)|\le C(|u|+u^2)\frac{\log T_0}{T_0}$ for an absolute constant $C$ (the proof gives $C=4$).

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `sm:thm:strip`.

<a id="poisson_trace_completely_monotone_iff_rh"></a>

### Complete monotonicity of the Poisson trace is equivalent to RH

Corollary 10.1.4, page 431 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, let $\rho_k=\beta_k+i\gamma_k$ ($\gamma_k>0$) be the nontrivial zeros of $\zeta$ in the upper half plane, counted with multiplicity, and put $a_k=-(\rho_k-\frac12)^2$ and $\sqrt{a_k}=\gamma_k+i(\frac12-\beta_k)$ (the root with positive real part). Define the heat trace $W(t)=\sum_ke^{-a_kt}$ and the Poisson trace $K_{1/2}(y)=\sum_ke^{-y\sqrt{a_k}}$ for $t,y>0$. For $\alpha\ge1$ put $v_\alpha(\theta)=\operatorname{Re}\,(\xi'/\xi)(\alpha+i\theta)$, $b_\alpha=\xi'(\alpha)/\xi(\alpha)$ and $K_\alpha(y)=\frac1\pi\int_0^\infty e^{-\theta y}v_\alpha(\theta)\,d\theta$. Then the Riemann Hypothesis holds if and only if $K_{1/2}$ is completely monotone on $(0,\infty)$, and if and only if $W$ is completely monotone on $(0,\infty)$ (the equivalence with $W$ is the earlier result of Polson, 2026, which the corollary recalls). For every $\alpha\ge1$, with no hypothesis, $K_\alpha$ is completely monotone and $K_\alpha(y)\sim b_\alpha/(\pi y)$ as $y\to\infty$; consequently the Lévy density $K_\alpha(|y|)/|y|$ of the law with characteristic function $\xi(\alpha)/\xi(\alpha+|u|)$ is asymptotic to $b_\alpha/(\pi y^2)$, the Lévy density of the Cauchy component.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `sm:cor:poisson`.

<a id="poisson_trace_derivative_conditions_hold_below_verified_height"></a>

### Unconditional range of the complete monotonicity conditions of the Poisson trace

Theorem 10.1.6, page 433 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\rho_k=\beta_k+i\gamma_k$ ($\gamma_k>0$) be the nontrivial zeros of $\zeta$ in the upper half plane, counted with multiplicity, $a_k=-(\rho_k-\frac12)^2$, $\sqrt{a_k}=\gamma_k+i\delta_k$ with $\delta_k=\frac12-\beta_k$, $\varphi_k=\arctan(\delta_k/\gamma_k)$, and let $K_{1/2}(y)=\sum_ke^{-y\sqrt{a_k}}$ be the Poisson trace of the zeros. For $n=0,1,2,\dots$ and $y>0$, $(-1)^nK_{1/2}^{(n)}(y)=\sum_k|\sqrt{a_k}|^n\,e^{-\gamma_ky}\cos(n\varphi_k-\delta_ky)$. Suppose that every zero with $0<\gamma_k\le H$ lies on the critical line. Then for $0<y<\pi$ and every integer $n$ with $0\le n<H(\pi-y)$, every term of this series is positive, and $(-1)^nK_{1/2}^{(n)}(y)>0$. With $H=T_0=3\times10^{12}$, the height up to which Platt and Trudgian verified that all zeros lie on the critical line, this covers all $n<9.12\times10^{12}$ at $y=0.1$, all $n<6.42\times10^{12}$ at $y=1$ and all $n<4.24\times10^{11}$ at $y=3$.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `thm:biane:cmrange`.

<a id="poisson_trace_positive_definite_and_bounded_forms_of_rh"></a>

### Positive definite and bounded forms of the Poisson trace criterion for RH

Proposition 10.1.7, page 435 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\rho_k=\beta_k+i\gamma_k$ ($\gamma_k>0$) be the nontrivial zeros of $\zeta$ in the upper half plane, counted with multiplicity, and put $a_k=-(\rho_k-\frac12)^2$ and $\sqrt{a_k}=\gamma_k+i(\frac12-\beta_k)$. For complex $z$ with $\operatorname{Re}z>0$ let $K_{1/2}(z)=\sum_ke^{-z\sqrt{a_k}}$ (the Poisson trace of the zeros, absolutely and locally uniformly convergent, hence analytic). The following are equivalent: (i) the Riemann Hypothesis; (ii) $K_{1/2}$ is completely monotone on $(0,\infty)$; (iii) for every $x>0$ the function $y\mapsto K_{1/2}(x+iy)$ is positive definite on $\mathbb{R}$; (iv) for some $x>0$ the function $y\mapsto K_{1/2}(x+iy)$ is bounded on $[0,\infty)$. Under these conditions $|K_{1/2}(x+iy)|\le K_{1/2}(x)$ for $x>0$ and real $y$.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `prop:biane:pd`.

<a id="poisson_trace_versus_polya_trace_small_y_asymptotics"></a>

### Small y asymptotics of the Poisson trace and Pólya's trace

Proposition 10.1.8, page 437 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\rho_k=\beta_k+i\gamma_k$ ($\gamma_k>0$) be the nontrivial zeros of $\zeta$ in the upper half plane, counted with multiplicity, $a_k=-(\rho_k-\frac12)^2$, $\sqrt{a_k}=\gamma_k+i(\frac12-\beta_k)$, and $K_{1/2}(y)=\sum_ke^{-y\sqrt{a_k}}$ (the Poisson trace of the zeros). Let $P(\theta)=K_{i\theta/2}(2\pi)$, with $K_\nu$ the modified Bessel function of the second kind (Pólya's function), let $\theta_1\le\theta_2\le\cdots$ be its positive zeros, counted with multiplicity, and $K^P(y)=\sum_ne^{-\theta_ny}$. Let $\gamma_E$ be Euler's constant, $M(y)=\frac1{2\pi y}\big(\log\frac1{2\pi y}-\gamma_E\big)$ and $c=2\pi-\frac1{6\pi}-\frac1{48\pi}=6.2235022035$. As $y\to0+$: (i) $K_{1/2}(y)=M(y)+\frac78+O(y\log\frac1y)$, and under the Riemann Hypothesis $K_{1/2}(y)=M(y)+\frac78+\frac1{48\pi}y\log\frac1y+o(y\log\frac1y)$; (ii) $K^P(y)=M(y)-\frac14+\big(2\pi-\frac1{6\pi}\big)y\log\frac1y+O(y)$; (iii) $D(y)=K_{1/2}(y)-K^P(y)=\frac98+O(y\log\frac1y)$, and under the Riemann Hypothesis $D(y)=\frac98-c\,y\log\frac1y+o(y\log\frac1y)$.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `prop:biane:traces`.

<a id="tilted_polya_kernel_translates_obstruction_and_approximation_cost"></a>

### Explicit obstruction and approximation cost for translates of the tilted Pólya kernel

Proposition 10.1.9, page 440 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and let $\Phi(u)=2\sum_{n\ge1}\big(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\big)\exp(-\pi n^2e^{2u})$ be Pólya's kernel, an even, positive function that decays doubly exponentially and satisfies $\int_{\mathbb{R}}\Phi(u)e^{su}\,du=\xi(\frac12+s)$ for every complex $s$. For $\alpha\ge\frac12$ put $\Phi_\alpha(u)=e^{(\alpha-\frac12)u}\Phi(u)$, so that $\int_{\mathbb{R}}\Phi_\alpha(u)e^{i\theta u}\,du=\xi(\alpha+i\theta)$, and for $g\in L^1(\mathbb{R})$ put $\hat g(\theta)=\int_{\mathbb{R}}g(u)e^{i\theta u}\,du$. (1) The translates $\Phi_\alpha(\cdot-x)$, $x\in\mathbb{R}$, span a dense subspace of $L^1(\mathbb{R})$ if and only if $\zeta$ has no zeros on the line $\operatorname{Re}s=\alpha$. (2) At the centre ($\alpha=\frac12$) they are not dense: for every finite combination $h=\sum_kc_k\Phi(\cdot-x_k)$ and every $g\in L^1(\mathbb{R})$, $\|g-h\|_1\ge\sup_\gamma|\hat g(\gamma)|$, the supremum over the ordinates $\gamma$ of the zeros of $\zeta$ on the critical line. (3) At the basepoint one they are dense, with no hypothesis on the zeros. (4) For $\alpha_0\in[\frac12,1)$, $\zeta$ has no zeros in the half plane $\operatorname{Re}s>\alpha_0$ if and only if the translates of $\Phi_\alpha$ are dense in $L^1(\mathbb{R})$ for every $\alpha>\alpha_0$; in particular the Riemann Hypothesis holds if and only if they are dense for every $\alpha>\frac12$. (5) If the translates of $\Phi_\alpha$ are dense, $h=\sum_kc_k\Phi_\alpha(\cdot-x_k)$, $g\in L^1(\mathbb{R})$ and $\|g-h\|_1\le\delta$, then $\sum_k|c_k|\ge(|\hat g(\gamma)|-\delta)/|\xi(\alpha+i\gamma)|$ for every ordinate $\gamma$ of a zero of $\zeta$. Parts (1), (3) and (4) restate Wiener's closure theorem (Wiener 1932) for $\Phi_\alpha$; parts (2) and (5) are what the completeness form adds.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `wr:prop:wiener`.

<a id="polya_kernel_analogue_of_salem_criterion"></a>

### Pólya kernel analogue of Salem's criterion

Corollary 10.1.10, page 442 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, let $\Phi$ be Pólya's kernel, the even positive function with $\int_{\mathbb{R}}\Phi(u)e^{su}\,du=\xi(\frac12+s)$ for every complex $s$, and put $\Phi_\alpha(u)=e^{(\alpha-\frac12)u}\Phi(u)$ and $\Xi(\gamma)=\xi(\frac12+i\gamma)$. The Riemann Hypothesis holds if and only if, for every $\frac12<\alpha<1$, the equation $\int_{\mathbb{R}}\Phi_\alpha(x-y)f(y)\,dy=0$ for all $x\in\mathbb{R}$ has no bounded solution $f\ne0$. At the centre ($\alpha=\frac12$) the equation has the explicit bounded solutions $f(y)=\cos(\gamma y)$, one for each ordinate $\gamma$ of a zero of $\zeta$ on the critical line, because $\Xi(\gamma)=0$ there; at the basepoint one ($\alpha=1$) it has none. Salem (1953) proved the corresponding statement with the Fermi Dirac kernel in place of $\Phi_\alpha$.

Source: [`ch/strip.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/strip.tex), label `wr:cor:salem`.

<a id="twisted_riesz_function_and_basepoint_riesz_criterion"></a>

### Twisted Riesz function and the basepoint form of Riesz's criterion

Proposition 10.2.1, page 447 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Möbius function and, for real $h$ and $x>0$, let $R_h(x)=x\sum_{n\ge1}\frac{\mu(n)}{n^2}\big(\frac x{n^2}\big)^{ih}e^{-x/n^2}$; for $h=0$ this is Riesz's function $R(x)=x\sum_{n\ge1}\mu(n)n^{-2}e^{-x/n^2}=-\sum_{k\ge1}\frac{(-x)^k}{(k-1)!\,\zeta(2k)}$. For every real $h$ and $\frac12<\operatorname{Re}s<1$, $\int_0^\infty R_h(x)\,x^{-s-1}\,dx=\Gamma(1-s+ih)/\zeta(2s)$. For $\frac12\le\alpha<1$, $\zeta$ has no zeros in the half plane $\operatorname{Re}s>\alpha$ if and only if $R_h(x)=O(x^{\alpha/2+\delta})$ as $x\to\infty$ for every $\delta>0$. In particular the Riemann Hypothesis holds if and only if $R_h(x)=O(x^{1/4+\delta})$ for every $\delta>0$. For $h=0$ and $\alpha=\frac12$ this is Riesz's criterion (Riesz 1916).

Source: [`ch/riesz.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/riesz.tex), label `ri:prop:twisted`.

<a id="riesz_function_as_centre_clock_average"></a>

### The Riesz function as an average over the centre clock

Proposition 10.2.4, page 450 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $R(x)=x\sum_{n\ge1}\mu(n)n^{-2}e^{-x/n^2}=-\sum_{k\ge1}\frac{(-x)^k}{(k-1)!\,\zeta(2k)}$ be Riesz's function ($\mu$ the Möbius function), $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ the completed Riemann zeta function, and $T$ the centre clock, the positive random variable with $\mathbb{E}\,e^{-\lambda T}=\xi(\frac12)/\xi(\frac12+\sqrt\lambda)$, $\lambda\ge0$. For every $x>0$, with no hypothesis on the zeros of $\zeta$, $R(x)=-\frac1{\xi(1/2)}\,\mathbb{E}\Big[\sum_{k\ge1}k(2k-1)\Big(-\frac x\pi\Big)^ke^{-(2k-\frac12)^2T}\Big]$, the expectation of the series being the series of the expectations.

Source: [`ch/riesz.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/riesz.tex), label `ri:prop:avg`.

<a id="riesz_coefficients_as_ggc_laplace_values"></a>

### Riesz coefficients as values of GGC Laplace transforms

Proposition 10.2.5, page 452 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, and for $\alpha\ge\frac12$ let $T_\alpha$ be the positive random variable with Laplace transform $\mathbb{E}\,e^{-\lambda T_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, $\lambda\ge0$ (it exists with no hypothesis on the zeros). Put $\kappa_k=2k(2k-1)(k-1)!/\pi^k$ and $\rho_{1/2}(\theta)=\frac1\pi\operatorname{Re}\,(\xi'/\xi)(1+i\theta)$. For every $\alpha\in[\frac12,1]$ and every integer $k\ge1$, $\frac1{\zeta(2k)}=\frac{\kappa_k}{2\xi(\alpha)}\,\mathbb{E}\,e^{-(2k-\alpha)^2T_\alpha}$. In particular, for the centre clock $T_{1/2}$ and the basepoint one clock $T_1$, where $\xi(1)=\frac12$, $\frac1{\zeta(2k)}=\frac{\kappa_k}{2\xi(\frac12)}\,\mathbb{E}\,e^{-(2k-\frac12)^2T_{1/2}}=\kappa_k\,\mathbb{E}\,e^{-(2k-1)^2T_1}$. At $\alpha=1$, with no hypothesis on the zeros, $T_1$ is a generalised gamma convolution without drift whose Thorin measure is $\rho_{1/2}(\theta)\,d\theta$ under the substitution $z=\theta^2$, and $\frac1{\zeta(2k)}=\kappa_k\exp\Big\{-\int_0^\infty\log\Big(1+\frac{(2k-1)^2}{\theta^2}\Big)\rho_{1/2}(\theta)\,d\theta\Big\}$.

Source: [`ch/riesz.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/riesz.tex), label `ri:prop:ggccoef`.

<a id="riesz_criterion_as_ggc_property_of_basepoint_clocks"></a>

### GGC form of Riesz's criterion

Theorem 10.2.7, page 454 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, $R(x)=x\sum_{n\ge1}\mu(n)n^{-2}e^{-x/n^2}=-\sum_{k\ge1}\frac{(-x)^k}{(k-1)!\,\zeta(2k)}$ Riesz's function ($\mu$ the Möbius function), $\Lambda$ the von Mangoldt function and $v_\alpha(\theta)=\operatorname{Re}\,(\xi'/\xi)(\alpha+i\theta)$. For $\alpha>\frac12$ let $T_\alpha$ be the positive random variable with Laplace transform $\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, $\lambda\ge0$ (it exists with no hypothesis on the zeros); a generalised gamma convolution (GGC) is a law on $[0,\infty)$ with Laplace transform $\exp\{-a\lambda-\int_{(0,\infty)}\log(1+\lambda/z)\,U(dz)\}$ for a drift $a\ge0$ and a Thorin measure $U$. For $\frac12\le\alpha<1$ the following are equivalent. (1) $R(x)=O(x^{\alpha/2+\delta})$ as $x\to\infty$ for every $\delta>0$. (2) $\zeta$ has no zeros with $\operatorname{Re}s>\alpha$. (3) $T_{\alpha'}$ is a GGC for every $\alpha'>\alpha$. (4) For every $\alpha'>\alpha$, $v_{\alpha'}(\theta)\ge0$ at every real $\theta$ for which $\alpha'+i\theta$ is not a zero of $\xi$. (5) For every $\alpha'>\alpha$ and every real $\theta$ the limit $\lim_{L\to\infty}\Big[\sum_{\log n\le L}\frac{\Lambda(n)}{n^{\alpha'}}\,2\sin^2\big(\tfrac12\theta\log n\big)-\int_0^L2\sin^2\big(\tfrac12\theta x\big)\,e^{(1-\alpha')x}\,dx\Big]$ exists. In particular ($\alpha=\frac12$) the Riemann Hypothesis holds if and only if $T_\alpha$ is a GGC for every $\alpha>\frac12$. For $\alpha=\frac12$ the equivalence of (1) and (2) is Riesz's criterion (Riesz 1916), and that of (2) and (4) is Lagarias's positivity criterion (Lagarias 1999), whose argument gives the equivalence of (2) and (4) for every $\alpha$.

Source: [`ch/riesz.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/riesz.tex), label `ri:thm:ggcriesz`.

<a id="squarefree_reciprocal_sum_closed_form_constant"></a>

### Closed form constant of the squarefree reciprocal sum over number fields

Proposition 10.3.1, page 464 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K$ be a number field of degree $n\ge2$ with Dedekind zeta function $\zeta_K$, let $c_K=\lim_{s\to1}(s-1)\zeta_K(s)$, let $\gamma_K$ be its Euler Kronecker constant, defined by $\zeta_K(s)=c_K/(s-1)+c_K\gamma_K+O(s-1)$, and put $C_K=\gamma_K-2\zeta_K'(2)/\zeta_K(2)$. Then, as $x\to\infty$, $\sum_{N(\mathfrak a)\le x}\frac{\mu(\mathfrak a)^2}{N(\mathfrak a)}=\frac{c_K}{\zeta_K(2)}\big(\log x+C_K\big)+O\big(x^{-1/n}\log x\big)$, the sum running over the integral ideals $\mathfrak a$ of $K$ ($\mu$ the Möbius function of ideals, $N$ the norm).

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `prop:sqfree:alpha1`.

<a id="squarefree_dirichlet_series_as_gamma_factor_times_polya_ratio"></a>

### Squarefree Dirichlet series as a gamma factor times a ratio of Pólya moment generating functions

Proposition 10.3.3, page 466 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, let $\Phi$ be Pólya's kernel, the even positive function with $\int_{\mathbb{R}}\Phi(u)e^{su}\,du=\xi(\frac12+s)$ for every complex $s$, and let $X$ be the symmetric random variable with density $\Phi/\xi(\frac12)$, so that $\mathbb{E}\,e^{sX}=\xi(\frac12+s)/\xi(\frac12)$ for all complex $s$. Put $\Gamma_{\mathrm{sq}}(s)=\frac{2s-1}{s-1}\,2^s\pi^{-(s+1)/2}\,\Gamma\big(\frac{s+1}2\big)$. For every $s$ with $\xi(2s)\ne0$, $\sum_{n\ge1}\frac{\mu(n)^2}{n^s}=\Gamma_{\mathrm{sq}}(s)\,\frac{\mathbb{E}\,e^{(s-\frac12)X}}{\mathbb{E}\,e^{(2s-\frac12)X}}$, where $\mu$ is the Möbius function and the left side is read as the continuation $\zeta(s)/\zeta(2s)$ where the series diverges.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `prop:sqfree:polya`.

<a id="squarefree_constant_as_difference_of_tilted_polya_means"></a>

### Squarefree constant as a combination of two tilted Pólya means

Theorem 10.3.5, page 467 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, $b_\alpha=\xi'(\alpha)/\xi(\alpha)$ (the mean of the Pólya variable $V_\alpha$ tilted to the basepoint $\alpha$, with $\mathbb{E}\,e^{sV_\alpha}=\xi(\alpha+s)/\xi(\alpha)$), $a_k=-(\rho_k-\frac12)^2$ with one zero $\rho_k$ taken from each pair $\{\rho,1-\rho\}$ of nontrivial zeros of $\zeta$, $G(u)=\xi(\frac12+\sqrt u)/\xi(\frac12)$, $W(t)=\sum_ke^{-a_kt}$, and let $\gamma_E$ be Euler's constant. The constant $C=\gamma_E-2\zeta'(2)/\zeta(2)$ of $\sum_{n\le x}\mu(n)^2/n=\frac6{\pi^2}(\log x+C)+O(x^{-1/2})$ equals $2+\log2-\frac12\log\pi-\frac12\gamma_E+b_1-2b_2$. With no hypothesis on the zeros, $b_1-2b_2=\frac{G'}G\big(\frac14\big)-6\,\frac{G'}G\big(\frac94\big)=\sum_k\Big(\frac1{\frac14+a_k}-\frac6{\frac94+a_k}\Big)=\int_0^\infty\big(e^{-t/4}-6e^{-9t/4}\big)\,W(t)\,dt$, the sum running over the pairs $\{\rho,1-\rho\}$ and converging absolutely. Under the Riemann Hypothesis, $a_k=\gamma_k^2$ and $b_1-2b_2=\int\big((z+\frac14)^{-1}-6(z+\frac94)^{-1}\big)\,U_{\mathrm{RH}}(dz)$, where $U_{\mathrm{RH}}=\sum_{\gamma>0}\delta_{\gamma^2}$ (over the ordinates $\gamma$ of the zeros $\frac12+i\gamma$, $\gamma>0$, with multiplicity) is the Thorin measure of the clock $1/G$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `thm:sqfree:const`.

<a id="dedekind_clock_is_ggc_iff_rh_for_dedekind_zeta"></a>

### Dedekind zeta clock is a GGC if and only if RH holds for the field

Proposition 10.3.6, page 469 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K$ be a number field with $r_1$ real and $r_2$ pairs of complex embeddings and discriminant $d_K$, let $\zeta_K$ be its Dedekind zeta function, $\Gamma_{\mathbb{R}}(s)=\pi^{-s/2}\Gamma(s/2)$, $\Gamma_{\mathbb{C}}(s)=2(2\pi)^{-s}\Gamma(s)$, and $\xi_K(s)=s(s-1)\,|d_K|^{s/2}\,\Gamma_{\mathbb{R}}(s)^{r_1}\Gamma_{\mathbb{C}}(s)^{r_2}\,\zeta_K(s)$, an entire function of order one with $\xi_K(s)=\xi_K(1-s)$ whose zeros are the nontrivial zeros of $\zeta_K$. Suppose $\zeta_K(\frac12)\ne0$ and let $G_K(u)=\xi_K(\frac12+\sqrt u)/\xi_K(\frac12)$. Then $1/G_K$ is the Laplace transform of a generalised gamma convolution if and only if every zero of $\zeta_K$ in the critical strip lies on the line $\operatorname{Re}s=\frac12$ (the Riemann Hypothesis for $\zeta_K$). Its drift is then $0$ and its Thorin measure is $\sum_{\gamma>0}\delta_{\gamma^2}$, over the ordinates $\gamma$ of the zeros counted with multiplicity. This is the number field form of the rigidity theorem of Polson, 2026, with the same proof. A real zero $\beta\in(\frac12,1)$ of $\zeta_K$ would be a pole of $1/G_K$ at the point $(\beta-\frac12)^2$ of the positive axis.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `prop:sqfree:thorinK`.

<a id="squarefree_constant_of_a_number_field_through_zeros"></a>

### Squarefree constant of a number field through the zeros of its Dedekind zeta function

Theorem 10.3.7, page 470 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $K$ be a number field with $r_1$ real and $r_2$ pairs of complex embeddings and discriminant $d_K$, with Dedekind zeta function $\zeta_K$, residue $c_K=\lim_{s\to1}(s-1)\zeta_K(s)$ and Euler Kronecker constant $\gamma_K$ (defined by $\zeta_K(s)=c_K/(s-1)+c_K\gamma_K+O(s-1)$), and put $C_K=\gamma_K-2\zeta_K'(2)/\zeta_K(2)$. Let $\Gamma_{\mathbb{R}}(s)=\pi^{-s/2}\Gamma(s/2)$, $\Gamma_{\mathbb{C}}(s)=2(2\pi)^{-s}\Gamma(s)$ and $\xi_K(s)=s(s-1)\,|d_K|^{s/2}\,\Gamma_{\mathbb{R}}(s)^{r_1}\Gamma_{\mathbb{C}}(s)^{r_2}\,\zeta_K(s)$, and let $\gamma_E$ be Euler's constant and $A_K=2+\frac12\log|d_K|+r_1\big(\log2-\frac12\log\pi-\frac12\gamma_E\big)+r_2\big(2-\gamma_E-\log2\pi\big)$. Then $C_K=A_K+\frac{\xi_K'}{\xi_K}(1)-2\,\frac{\xi_K'}{\xi_K}(2)=A_K+\sum_k\Big(\frac1{\frac14+a_{K,k}}-\frac6{\frac94+a_{K,k}}\Big)$, the sum running over the pairs $\{\rho,1-\rho\}$ of zeros of $\xi_K$, with multiplicity, where $a_{K,k}=-(\rho_k-\frac12)^2$ and a zero of order $2m$ at $s=\frac12$ counts as $m$ pairs with $a_{K,k}=0$. Under the Riemann Hypothesis for $\zeta_K$ and $\zeta_K(\frac12)\ne0$, the zero part is $\frac{G_K'}{G_K}\big(\frac14\big)-6\,\frac{G_K'}{G_K}\big(\frac94\big)$ with $G_K(u)=\xi_K(\frac12+\sqrt u)/\xi_K(\frac12)$, the same combination of values of the Stieltjes transform of the Thorin measure $\sum_{\gamma>0}\delta_{\gamma^2}$ of $1/G_K$ as over $\mathbb{Q}$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `thm:sqfree:constK`.

<a id="smoothed_squarefree_counts_equivalent_to_rh"></a>

### Exponentially and Gamma(2) smoothed squarefree counts are equivalent to RH

Theorem 10.3.9, page 473 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Möbius function, and for $x>0$ let $Q_1(x)=\sum_{n\ge1}\mu(n)^2e^{-n/x}$ and $Q_2(x)=\sum_{n\ge1}\mu(n)^2\big(1+\frac{2n}x\big)e^{-2n/x}$, the number of squarefree integers up to a random point of mean $x$ with the exponential law and with the Gamma(2) law respectively. (i) The Riemann Hypothesis holds if and only if $Q_1(x)=\frac6{\pi^2}x+O(x^{1/4+\varepsilon})$ as $x\to\infty$ for every $\varepsilon>0$, and if and only if $Q_2(x)=\frac6{\pi^2}x+O(x^{1/4+\varepsilon})$ for every $\varepsilon>0$. (ii) Under the Riemann Hypothesis, if the zeros are simple, $Q_1(x)=\frac6{\pi^2}x+1+\sum_\rho A_\rho x^{\rho/2}+\frac{c_1\log x+c_2}x+O(x^{-2})$ and $Q_2(x)=\frac6{\pi^2}x+1+\sum_\rho B_\rho x^{\rho/2}+\frac{2c_1}x+O(x^{-2})$, where $A_\rho=\frac{\Gamma(\rho/2)\,\zeta(\rho/2)}{2\zeta'(\rho)}$, $B_\rho=\big(1+\frac\rho2\big)\frac{2^{-\rho/2}\,\Gamma(\rho/2)\,\zeta(\rho/2)}{2\zeta'(\rho)}$, the sums over the nontrivial zeros $\rho$ are taken over $|\operatorname{Im}\rho|<2T_j$ along a sequence $T_j\to\infty$, $c_1=-\frac{\zeta(-1)}{2\zeta'(-2)}=-1.36843277762021$ and $c_2=c_1\big(1-\gamma_E+\frac{\zeta'(-1)}{\zeta(-1)}-\frac{\zeta''(-2)}{\zeta'(-2)}\big)=-0.339381226129599$, with $\gamma_E$ Euler's constant.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `thm:sqfree:smooth`.

<a id="squarefree_error_lower_bound_order_fourth_root_of_x"></a>

### Lower bound for the squarefree error of the order of the fourth root of x

Corollary 10.3.11, page 476 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $Q(x)=\sum_{n\le x}\mu(n)^2$ be the number of squarefree integers up to $x$ and $\Delta(x)=Q(x)-\frac6{\pi^2}x$. Let $\rho_1$ be the zero of $\zeta$ with least positive imaginary part, and $R_{\rho_1}=\frac{\zeta(\rho_1/2)}{\rho_1\zeta'(\rho_1)}$, so that $|R_{\rho_1}|=0.100404$. With no hypothesis, $\limsup_{x\to\infty}\Delta(x)/x^{1/4}\ge|R_{\rho_1}|$ and $\liminf_{x\to\infty}\Delta(x)/x^{1/4}\le-|R_{\rho_1}|$. The same holds for $Q_2(x)-\frac6{\pi^2}x-1$, where $Q_2(x)=\sum_{n\ge1}\mu(n)^2\big(1+\frac{2n}x\big)e^{-2n/x}$ is the squarefree count smoothed by a Gamma(2) law of mean $x$, with $|B_{\rho_1}|=9.944\times10^{-5}$ in place of $|R_{\rho_1}|$, where $B_\rho=\big(1+\frac\rho2\big)\frac{2^{-\rho/2}\,\Gamma(\rho/2)\,\zeta(\rho/2)}{2\zeta'(\rho)}$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `cor:sqfree:omega`.

<a id="squarefree_gamma_n_smoothing_bound_under_rh"></a>

### Squarefree count at a Gamma(n) random point under RH

Proposition 10.3.12, page 477 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on the Riemann Hypothesis.

Statement. Let $Q(x)=\sum_{m\le x}\mu(m)^2$ be the number of squarefree integers up to $x$, and for $n\ge1$ let $G_{n,x}$ have the Gamma($n$) law with mean $x$ (shape $n$, rate $n/x$). Assume the Riemann Hypothesis. Then, for every $\varepsilon>0$ and uniformly in $n\ge1$, $\mathbb{E}\,Q(G_{n,x})=\frac6{\pi^2}\,x+O_\varepsilon\big(x^{1/4+\varepsilon}\,n^{1/8+\varepsilon}\big)$ as $x\to\infty$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `prop:sqfree:gamman`.

<a id="gamma2_smoothed_squarefree_error_as_mertens_average"></a>

### Gamma(2) smoothed squarefree error as a positive average of the Mertens function

Theorem 10.3.15, page 479 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Möbius function, $M(t)=\sum_{d\le t}\mu(d)$ the Mertens function, and $Q_2(x)=\sum_{n\ge1}\mu(n)^2\big(1+\frac{2n}x\big)e^{-2n/x}$, the squarefree count smoothed by a Gamma(2) law of mean $x$. For $y>0$ let $g_2(y)=\frac2y-\frac1{e^y-1}-\frac y{4\sinh^2(y/2)}$ and $f_2(y)=-2g_2'(y)=\frac4{y^2}-\frac{y\cosh(y/2)}{2\sinh^3(y/2)}$. Then $f_2>0$ (this is Lazarević's inequality $(\sinh z/z)^3>\cosh z$ with $z=y/2$), $f_2(y)=y^2/60+O(y^4)$ as $y\to0$ and $f_2(y)\sim4/y^2$ as $y\to\infty$, and $f_2$ is the density of a probability law $\mathcal{V}_2$ on $(0,\infty)$ with $P(\mathcal{V}_2\ge y)=2g_2(y)$ and $\mathbb{E}\,\mathcal{V}_2^s=-2\,\Gamma(s+2)\,\zeta(s)$ for $-3<\operatorname{Re}s<1$. For every $x>0$, $Q_2(x)=\frac6{\pi^2}\,x-\frac12\,\mathbb{E}\,M\Big(\sqrt{x\mathcal{V}_2/2}\Big)$. Moreover $f_2$ is twice the difference between the continuous and the integer mixtures of the Gamma(2) densities $\gamma_{2,\lambda}(y)=\lambda^2ye^{-\lambda y}$ over their rate: $f_2(y)=2\Big(\int_0^\infty\gamma_{2,\lambda}(y)\,d\lambda-\sum_{m\ge1}\gamma_{2,m}(y)\Big)$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `thm:sqfree:mertens`.

<a id="gamma_smoothing_threshold_shape_two"></a>

### Gamma smoothing threshold, shape two

Theorem 10.3.18, page 482 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\mu$ be the Möbius function, $Q(x)=\sum_{n\le x}\mu(n)^2$ the number of squarefree integers up to $x$ and $M(t)=\sum_{d\le t}\mu(d)$ the Mertens function. For a shape $\nu>0$ let $R_\nu(y)=\frac{y^{\nu+1}}{\Gamma(\nu+1)}\sum_{m\ge1}m^\nu e^{-my}$, the right Riemann sum with step $y$ of the Gamma($\nu+1$) density, and $f_\nu(y)=\frac{2\nu\,(1-R_\nu(y))}{y^2}$ for $y>0$. The function $f_\nu$ is a signed density of total mass $1$, and for $Y$ with the Gamma($\nu$) law of rate $\lambda$, $\mathbb{E}\,Q(Y)=\frac6{\pi^2}\,\mathbb{E}\,Y-\frac12\int_0^\infty M\big(\sqrt{y/\lambda}\big)f_\nu(y)\,dy$. Then $\{\nu>0:\ f_\nu(y)\ge0\text{ for all }y>0\}=(0,2]$. Thus two is the largest shape for which the gamma smoothed squarefree error is a positive average of the Mertens function, and the obstruction is the trivial zero of $\zeta$ at $-2$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `thm:sqfree:threshold`.

<a id="squarefree_error_pointwise_bound_under_rh_open"></a>

### Pointwise bound for the squarefree error under RH (open)

Question 10.3.19, page 484 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $Q(x)=\sum_{n\le x}\mu(n)^2$ be the number of squarefree integers up to $x$ and $\Delta(x)=Q(x)-\frac6{\pi^2}x$. Open question: prove, under the Riemann Hypothesis, that $\Delta(x)=O(x^{1/4+\varepsilon})$ for every $\varepsilon>0$ (by the book, the order $x^{1/4}$ cannot be lowered, and under the Riemann Hypothesis the best bound proved there is $O(x^{2/5+\varepsilon})$, weaker than known conditional exponents such as $7/22$ (Baker and Pintz 1985)). The book gives these forms: remove the loss $n^{1/8}$ from the bound $\mathbb{E}\,Q(G_{n,x})=\frac6{\pi^2}x+O(x^{1/4+\varepsilon}n^{1/8+\varepsilon})$ valid under the Riemann Hypothesis for the count at a Gamma($n$) random point $G_{n,x}$ of mean $x$ (the conjecture implies the bound with no loss, and a loss $n^\delta$ in place of $n^{1/8}$ would give $\Delta(x)=O(x^{1-3/(8\delta+4)+\varepsilon})$, whose exponent tends to $\frac14$ as $\delta\to0$); or show that for large $X$ no $x\in[X,2X]$ has $|\Delta(x)|>CX^{1/4+\varepsilon}$, for some constant $C$ (under the Riemann Hypothesis the book shows that for $C=1$ the set of such $x$ has measure $O(X^{1-\varepsilon})$), which would follow from a bound $O(X^{1/4+\varepsilon})$ for the oscillation of $\Delta$ on the subintervals of $[X,2X]$ of length $X^{1-\varepsilon/2}$.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `q:sqfree:pointwise`.

<a id="squarefree_error_fejer_constant_and_unboundedness_open"></a>

### Sharpening the lower bound constant for the squarefree error (open)

Question 10.3.20, page 485 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $Q(x)=\sum_{n\le x}\mu(n)^2$ be the number of squarefree integers up to $x$, $\Delta(x)=Q(x)-\frac6{\pi^2}x$, $\rho_1$ the zero of $\zeta$ with least positive imaginary part and $R_{\rho_1}=\frac{\zeta(\rho_1/2)}{\rho_1\zeta'(\rho_1)}$, so that $|R_{\rho_1}|=0.1004$ (to four digits) and, by the book, $\limsup\Delta(x)/x^{1/4}\ge|R_{\rho_1}|$ and $\liminf\Delta(x)/x^{1/4}\le-|R_{\rho_1}|$ with no hypothesis. Open questions: sharpen the constant $|R_{\rho_1}|$ by using many zeros at once with a Fejér kernel, as Ingham did for the Mertens function, and compare it with the observed $\max|\Delta(x)|/x^{1/4}$, about $1$ for $x$ up to $2\times10^8$; and decide whether $\limsup|\Delta(x)|/x^{1/4}=\infty$, as the linear independence of the ordinates of the zeros would suggest.

Source: [`ch/sqfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sqfree.tex), label `q:sqfree:fejer`.

## Chapter 11. Infinite divisibility, the eta analogue and almost generalised gamma convolutions

<a id="centre_clock_background_driving_subordinator"></a>

### Background driving subordinator of the centre clock

Proposition 11.1.1, page 489 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function and $G(u)=\xi(\frac{1}{2}+\sqrt{u})/\xi(\frac{1}{2})$. The centre clock $T$ is the positive random variable with Laplace transform $\mathbb{E}e^{-uT}=1/G(u)$, $u\ge0$. Let $W(t)=\sum_k e^{-a_kt}$, $t>0$, be the heat trace of the zeros, where $a_k=-(\rho_k-\frac{1}{2})^2$ and $\rho_k$ runs over the zeros of $\xi$ with positive imaginary part, counted with multiplicity; under the Riemann Hypothesis (RH), $W(t)=\sum_{\gamma>0}e^{-\gamma^2t}$ with $\gamma$ the ordinates of the zeros. With no hypothesis (the proof rests on the inequalities $W>0$ and $W'<0$, whose proof uses only the verification of RH below height $1000$), $T\overset{d}{=}\int_0^\infty e^{-t}\,dZ_t$ for a subordinator $Z$ with Lévy density $h(x)=-W'(x)$, which is $\sum_{\gamma>0}\gamma^2e^{-\gamma^2x}$ under RH and is given for every $x>0$ by the derivative of Weil's explicit formula for $W$ with the test function $e^{-tr^2}$. The condition $\mathbb{E}\log(1+Z_1)<\infty$ holds because $h$ decays exponentially. Moreover RH holds if and only if $h$ is completely monotone.

Source: [`ch/sato.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sato.tex), label `sb:prop:bdlp`.

<a id="centre_clock_and_shifted_clocks_unimodal"></a>

### Unimodality of the centre clock and the shifted clocks

Corollary 11.1.2, page 491 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function. The centre clock $T$ is the positive random variable with Laplace transform $\mathbb{E}e^{-uT}=\xi(\frac{1}{2})/\xi(\frac{1}{2}+\sqrt{u})$, and for $\alpha>\frac{1}{2}$ the shifted clock $T_\alpha$ is the positive random variable with Laplace transform $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt{u})$. With no hypothesis, the law of $T$ is unimodal, and so is the law of $T_\alpha$ for every $\alpha>\frac{1}{2}$. Here a law is unimodal if its distribution function is convex to the left of some point, a mode, and concave to the right of it.

Source: [`ch/sato.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sato.tex), label `sb:cor:unimodal`.

<a id="polya_tilt_basepoint_one_not_quasi_infinitely_divisible"></a>

### The tilted Pólya variable at the basepoint one is not quasi infinitely divisible

Proposition 11.1.3, page 495 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function and let $V_1$ be the tilted Pólya variable at the basepoint one, the random variable with moment generating function $\mathbb{E}e^{sV_1}=\xi(1+s)/\xi(1)$ (the Esscher tilt by $e^{x/2}$ of the Pólya variable $X$, which has $\mathbb{E}e^{sX}=\xi(\frac{1}{2}+s)/\xi(\frac{1}{2})$). A law is quasi infinitely divisible if its characteristic function is the quotient of two infinitely divisible characteristic functions, equivalently if it has a Lévy Khintchine form with a signed Lévy measure $\nu^+-\nu^-$ in which $\nu^+$ and $\nu^-$ are both Lévy measures. Then $V_1$ is not quasi infinitely divisible, and so, an infinitely divisible law being quasi infinitely divisible, it is not infinitely divisible. The reflection $x\mapsto-x$ of the measure $\Pi_1(dx)=\big(e^{-3x}(1-e^{-2x})^{-1}-1\big)\,dx/x+\sum_{n\ge2}\Lambda(n)(n\log n)^{-1}\delta_{\log n}(dx)$, with $\Lambda$ the von Mangoldt function, is the candidate signed Lévy measure of $V_1$; its pole term $-dx/x$ and its prime atoms at $\log n$ each have infinite mass on $(1,\infty)$, so its positive part and its negative part both have infinite mass there, and its Lévy form converges only as a compensated limit, which exists by the prime number theorem.

Source: [`ch/sato.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sato.tex), label `sb:prop:qid`.

<a id="thorin_bernstein_criterion_for_shifted_clock_exponent"></a>

### Thorin Bernstein form of the shifted clock criterion

Proposition 11.1.4, page 498 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $G(u)=\xi(\frac{1}{2}+\sqrt{u})/\xi(\frac{1}{2})$ and, for $\alpha>\frac{1}{2}$, $f_\alpha(u)=u+2(\alpha-\frac{1}{2})\sqrt{u}$. A Thorin Bernstein function is the Laplace exponent of a generalised gamma convolution (GGC), equivalently a Bernstein function whose derivative is a Stieltjes function. For every $\alpha\ge1$ the composition $\log G\big((\alpha-\frac{1}{2})^2+f_\alpha(u)\big)-\log G\big((\alpha-\frac{1}{2})^2\big)$, which is the Laplace exponent of the shifted clock $T_\alpha$ (the positive random variable with $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt{u})$), is a Thorin Bernstein function of $u$ with no hypothesis, while $\log G$ itself, the Laplace exponent of the centre clock, is a Thorin Bernstein function if and only if the Riemann Hypothesis holds. For $\frac{1}{2}<\alpha<1$ the composition is a Thorin Bernstein function if and only if $\zeta$ has no zeros with real part greater than $\alpha$.

Source: [`ch/sato.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/sato.tex), label `sb:prop:tbf`.

<a id="upsilon_image_and_free_counterpart_of_centre_clock"></a>

### Upsilon image and free counterpart of the centre clock

Theorem 11.2.3, page 513 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $G(u)=\xi(\frac{1}{2}+\sqrt{u})/\xi(\frac{1}{2})=\prod_k(1+u/a_k)$ with $a_k=-(\rho_k-\frac{1}{2})^2$, where $\rho_k$ runs over the zeros of $\xi$ with positive imaginary part, counted with multiplicity, and let $T$ be the centre clock, the positive random variable with $\mathbb{E}e^{-uT}=1/G(u)$. Put $W(t)=\sum_ke^{-a_kt}$, $h=-W'$, $m_n=\frac{1}{n!}\int_0^\infty t^nW(t)\,dt=\sum_ka_k^{-n-1}$, let $\gamma_1=14.1347\ldots$ be the least ordinate of a zero and $N(\theta)$ the number of zeros with ordinate in $(0,\theta]$. With no hypothesis $T$ is infinitely divisible with no drift and Lévy measure $W(t)\,dt/t$, and self decomposable (because $W>0$ and $W'<0$ on $(0,\infty)$). For a law $\mu$ on $[0,\infty)$ write $\Psi_\mu(s)=-\log\int e^{-sx}\mu(dx)$ for its Laplace exponent. A generalised gamma convolution (GGC) is a law on $[0,\infty)$ with $\Psi(s)=ds+\int\log(1+s/t)\,U(dt)$ for a drift $d\ge0$ and a Thorin measure $U$, and BO is the class of infinitely divisible laws on $[0,\infty)$ with a completely monotone Lévy density. The Upsilon transform $\Upsilon$ maps an infinitely divisible law $\mu$, the law of $X_1$ for a Lévy process $X$, to the law of $\int_0^1\log(1/t)\,dX_t$. The Bercovici Pata bijection $\Lambda_{\mathrm{BP}}$ maps an infinitely divisible law to the freely infinitely divisible law with the same Lévy Khintchine triplet; $G_\mu(z)=\int\mu(dx)/(z-x)$ is the Cauchy transform of $\mu$ and $\mathcal{C}_\mu(z)=z\varphi_\mu(1/z)$ its free cumulant transform, with $\varphi_\mu$ the Voiculescu transform. A freely infinitely divisible law $\mu$ is free regular if its free convolution powers $\mu^{\boxplus t}$ are carried by $[0,\infty)$ for every $t>0$ (these are the images under $\Lambda_{\mathrm{BP}}$ of the infinitely divisible laws on $[0,\infty)$). (a) With no hypothesis, $\Upsilon(T)$ is a GGC with no drift and Thorin measure $\widetilde{U}(du)=u^{-2}h(1/u)\,du$, so that $\widetilde{U}((0,u])=W(1/u)>0$ for every $u>0$. Its Laplace exponent is $\Psi_{\Upsilon(T)}(s)=\int_0^\infty\log G(sx)\,e^{-x}\,dx=\sum_ke^{a_k/s}E_1(a_k/s)$, $s>0$, with the principal branch of the exponential integral $E_1$; its Lévy density is the completely monotone function $\ell(z)=z^{-1}\int_0^\infty W(z/w)\,e^{-w}\,dw$; its cumulants are $\kappa_n(\Upsilon(T))=n!\,(n-1)!\,m_{n-1}$; and $\int u^{-n}\,\widetilde{U}(du)=n!\,m_{n-1}$ for $n\ge1$. Moreover $\mathbb{E}e^{\lambda\Upsilon(T)}=\infty$ for every $\lambda>0$, while $\mathbb{E}e^{\lambda T}<\infty$ exactly when $\lambda<\gamma_1^2$. (b) With no hypothesis, $\Lambda_{\mathrm{BP}}(T)$ is free regular, freely self decomposable (for every $c\in(0,1)$ it is the free convolution of its dilation by $c$ with some law) and unimodal. Its free cumulant transform is $\mathcal{C}_{\Lambda_{\mathrm{BP}}(T)}(z)=\int_0^\infty zW(x)/(1-zx)\,dx$, $z\in\mathbb{C}^-$ (the lower half plane), where the right side continues it analytically to $\mathbb{C}\setminus(0,\infty)$, with $\mathcal{C}_{\Lambda_{\mathrm{BP}}(T)}(-s)=-\Psi_{\Upsilon(T)}(s)$ for $s>0$; its free cumulants are $(n-1)!\,m_{n-1}$, the cumulants of $T$. Its Cauchy transform $G_{\Lambda_{\mathrm{BP}}(T)}$ satisfies $\Psi_{\Upsilon(T)}\big(-G_{\Lambda_{\mathrm{BP}}(T)}(z)\big)=1-z\,G_{\Lambda_{\mathrm{BP}}(T)}(z)$ for $z\in\mathbb{C}^+$ (the upper half plane), where $\Psi_{\Upsilon(T)}$ is continued analytically to $\mathbb{C}\setminus(-\infty,0)$ by $\Psi_{\Upsilon(T)}(s)=\int_0^\infty sW(y)/(1+sy)\,dy$, and its support lies in $[c,\infty)$, where $c=\max_{s>0}(\Psi_{\Upsilon(T)}(s)-1)/s>0$. (c) RH holds if and only if $T$ is the Upsilon image of an infinitely divisible law on $[0,\infty)$. In that case $T=\Upsilon(\mu_0)$, where $\mu_0$ is the self decomposable law with no drift and Lévy measure $N(x^{-1/2})\,dx/x$ (canonical function $N(v^{-1/2})$); $\mu_0$ is neither a GGC nor in BO. With no hypothesis, $T$ is not the Upsilon image of a law in BO.

Source: [`ch/wigfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/wigfree.tex), label `thm:wigfree:clock`.

<a id="gamma_layer_peeling_keeps_real_zeros_under_rh"></a>

### Peeling a gamma layer from the Pólya generating function keeps zeros real under RH

Proposition 11.2.8, page 525 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on The Riemann Hypothesis.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function and $G(u)=\xi(\frac{1}{2}+\sqrt{u})/\xi(\frac{1}{2})=\sum_{k\ge0}b_ku^k$, so that the Pólya variable $X$, with $\mathbb{E}e^{sX}=G(s^2)$, has $b_k=\mathbb{E}X^{2k}/(2k)!$. For $a>0$ let $(a)_k=a(a+1)\cdots(a+k-1)$ and $G_a(u)=\sum_{k\ge0}b_ku^k/(a)_k$. If the Riemann Hypothesis holds, then for every $a>0$ the function $G_a$ is in the Laguerre Pólya class (the real entire functions that are locally uniform limits of real polynomials with only real zeros) with only real negative zeros, and all its Jensen polynomials $\sum_{i=0}^{d}\binom{d}{i}(n+i)!\,b_{n+i}x^i/(a)_{n+i}$, $d,n\ge0$, have only real negative zeros. Multiplying the coefficients by $(a)_k$ instead does not preserve real zeros: it turns $(1+u)^2=1+2u+u^2$ into $1+2au+a(a+1)u^2$, whose discriminant is $-4a<0$.

Source: [`ch/wigfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/wigfree.tex), label `prop:wigfree:peel`.

<a id="no_gamma_scale_mixture_for_polya_and_williams_variables"></a>

### No gamma scale mixture for the Pólya variable and Williams' variable

Theorem 11.2.9, page 526 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function. The Pólya variable $X$ is the real random variable with $\mathbb{E}e^{sX}=\xi(\frac{1}{2}+s)/\xi(\frac{1}{2})$, $s\in\mathbb{C}$, and Williams' variable $S_2$ is the positive random variable with $\mathbb{E}e^{-\lambda S_2}=\big(\sqrt{2\lambda}/\sinh\sqrt{2\lambda}\big)^2$, $\lambda\ge0$ (the sum of two independent hitting times of the level one by a three dimensional Bessel process), for which $\mathbb{E}(\frac{\pi}{2}S_2)^{s/2}=2\xi(s)$. A law is a scale mixture of Gamma($a$) laws if it is the law of $\Theta\,\Gamma$ with $\Gamma$ gamma distributed with shape $a$ and unit scale and independent of a random variable $\Theta\ge0$. Then for no $a>0$ is $|X|$, $X^2$ or $S_2$ a scale mixture of Gamma($a$) laws.

Source: [`ch/wigfree.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/wigfree.tex), label `thm:wigfree:nomix`.

<a id="eta_analogue_thorin_measure_positive_at_basepoint_one"></a>

### Thorin measure of the Dedekind eta analogue at the basepoint one

Proposition 11.3.1, page 531 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\chi_{12}$ be the real primitive even Dirichlet character modulo $12$ ($\chi_{12}(n)=1$ for $n\equiv\pm1$ and $-1$ for $n\equiv\pm5 \pmod{12}$, and $0$ if $\gcd(n,12)>1$), $L(s,\chi_{12})$ its Dirichlet $L$ function and $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(s/2)L(s,\chi_{12})$ the completed function, which is entire (there is no pole to remove) and satisfies $\xi(s,\chi_{12})=\xi(1-s,\chi_{12})$. Put $v_1(\theta)=\operatorname{Re}(\xi'/\xi)(1+i\theta,\chi_{12})$. Then $\varphi_1(u)=\xi(1,\chi_{12})/\xi(1+\sqrt{u},\chi_{12})$ is the Laplace transform of a generalised gamma convolution (GGC) with no drift, whose Thorin measure, in the variable $\theta=\sqrt{z}$, has the density $\frac{1}{\pi}v_1(\theta)=\frac{1}{\pi}\sum_\rho\frac{1-\beta}{(1-\beta)^2+(\theta-\gamma)^2}>0$, the sum over the nontrivial zeros $\rho=\beta+i\gamma$ of $L(s,\chi_{12})$, in conjugate pairs. Positivity is unconditional because $L(1+it,\chi_{12})\ne0$ for real $t$.

Source: [`ch/bpeta.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bpeta.tex), label `bp:eta:prop:thorin`.

<a id="eta_analogue_levy_form_at_basepoint_one"></a>

### Lévy form of the Dedekind eta analogue at the basepoint one

Proposition 11.3.2, page 532 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\chi_{12}$ be the real primitive even Dirichlet character modulo $12$ ($\chi_{12}(n)=1$ for $n\equiv\pm1$ and $-1$ for $n\equiv\pm5 \pmod{12}$, and $0$ if $\gcd(n,12)>1$), $\xi(s,\chi_{12})=(12/\pi)^{s/2}\Gamma(s/2)L(s,\chi_{12})$ the completed function of its Dirichlet $L$ function, $b_1=(\xi'/\xi)(1,\chi_{12})$ (which equals $\sum_\rho\rho^{-1}$ over the nontrivial zeros and is $0.1650833123$), $\Lambda$ the von Mangoldt function and $\Lambda_{\chi_{12}}(n)=\chi_{12}(n)\Lambda(n)$. For $\operatorname{Re}s>0$, $\log\frac{\xi(1+s,\chi_{12})}{\xi(1,\chi_{12})}=b_1s+\int_0^\infty(e^{-sx}-1+sx)\,\Pi(dx)$, where $\Pi(dx)=\frac{e^{-x}}{x(1-e^{-2x})}\,dx+\sum_{n\ge2}\frac{\Lambda_{\chi_{12}}(n)}{n\log n}\,\delta_{\log n}(dx)$, the prime sum converging conditionally, as a limit over $n\le X$ as $X\to\infty$.

Source: [`ch/bpeta.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/bpeta.tex), label `bp:eta:prop:levy`.

<a id="almost_ggc_theorem"></a>

### Almost GGC theorem for the shifted clocks

Theorem 11.4.2, page 541 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, let $\rho=\beta+i\gamma$ run over its zeros (the nontrivial zeros of $\zeta$), counted with multiplicity, let $v_\alpha(\theta)=\operatorname{Re}(\xi'/\xi)(\alpha+i\theta)$ for $\alpha>\frac{1}{2}$ and let $P_w(x)=w/(\pi(w^2+x^2))$ be the Cauchy kernel of width $w>0$, with $P_0$ the unit point mass at $0$. Write $v_\alpha=v_\alpha^+-v_\alpha^-$ with $v_\alpha^+(\theta)/\pi=\sum_{\beta\le\alpha,\,\gamma>0}[P_{\alpha-\beta}(\theta-\gamma)+P_{\alpha-\beta}(\theta+\gamma)]$ and $v_\alpha^-(\theta)/\pi=\sum_{\beta>\alpha,\,\gamma>0}[P_{\beta-\alpha}(\theta-\gamma)+P_{\beta-\alpha}(\theta+\gamma)]$ (a zero on the line $\operatorname{Re}s=\alpha$ contributes a unit atom at its ordinate to the positive part), and put $m_\alpha^\pm(T)=\pi^{-1}\int_0^Tv_\alpha^\pm(\theta)\,d\theta$. Recall that the shifted clock $T_\alpha$, with Laplace transform $\xi(\alpha)/\xi(\alpha+\sqrt{u})$, is a generalised gamma convolution (GGC) exactly when no zero has $\beta>\alpha$, that is when $v_\alpha^-\equiv0$. For every $\frac{1}{2}<\alpha<1$, with no hypothesis, $m_\alpha^-(T)/m_\alpha^+(T)\le T^{e(\alpha)-1+o(1)}$ as $T\to\infty$, where $e(\alpha)=3(1-\alpha)/(2-\alpha)$ for $\frac{1}{2}<\alpha<\frac{7}{10}$ and $e(\alpha)=15(1-\alpha)/(3+5\alpha)$ for $\frac{7}{10}\le\alpha<1$, and $e(\alpha)<1$ throughout; so the negative part $v_\alpha^-/\pi$, the Cauchy kernels of the zeros to the right of $\alpha$, carries on $[0,T]$ a vanishing fraction of the mass of the positive part, with a power saving. For $\alpha\ge1$ no zero lies to the right of $\alpha$ and $m_\alpha^-=0$. Under the density hypothesis ($N(\sigma,T)\le T^{2(1-\sigma)+o(1)}$, where $N(\sigma,T)$ counts the zeros with $\beta\ge\sigma$ and $0<\gamma\le T$) the same bound holds with $e(\alpha)=2(1-\alpha)$. RH holds if and only if $v_\alpha^-$ vanishes identically for every $\alpha>\frac{1}{2}$. The proof combines the zero density estimates of Ingham, of Guth and Maynard and of Huxley.

Source: [`ch/almost.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/almost.tex), label `ag:thm:almost`.

<a id="canonical_factorisation_of_shifted_clock_transform"></a>

### Canonical factorisation of the shifted clock transform into a GGC and a defect

Proposition 11.4.3, page 542 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, let $\rho=\beta+i\gamma$ run over its zeros (the nontrivial zeros of $\zeta$), counted with multiplicity, let $P_w(x)=w/(\pi(w^2+x^2))$ be the Cauchy kernel of width $w>0$, with $P_0$ the unit point mass at $0$, and for $\alpha>\frac{1}{2}$ let $\varphi_\alpha(u)=\xi(\alpha)/\xi(\alpha+\sqrt{u})$ be the Laplace transform of the shifted clock $T_\alpha$ and $v_\alpha^+(\theta)/\pi=\sum_{\beta\le\alpha,\,\gamma>0}[P_{\alpha-\beta}(\theta-\gamma)+P_{\alpha-\beta}(\theta+\gamma)]$ (a zero on the line $\operatorname{Re}s=\alpha$ contributes a unit atom at its ordinate). For every $\alpha>\frac{1}{2}$ and $u\ge0$, $\varphi_\alpha(u)=\varphi_\alpha^+(u)\,D_\alpha(u)$, where $\varphi_\alpha^+(u)=\exp\big\{-\int_0^\infty\log\big(1+u/\theta^2\big)\,\frac{v_\alpha^+(\theta)}{\pi}\,d\theta\big\}$ and $D_\alpha(u)=\prod_{\beta>\alpha,\,\gamma>0}\frac{(\beta-\alpha)^2+\gamma^2}{(\sqrt{u}-\beta+\alpha)^2+\gamma^2}$. Here $\varphi_\alpha^+$ is the Laplace transform of a generalised gamma convolution (GGC) $T_\alpha^+$ with Thorin density $v_\alpha^+/\pi$ in the chart $\theta=\sqrt{z}$, and the product defining the defect $D_\alpha$ runs over the zeros to the right of the line $\operatorname{Re}s=\alpha$ and converges absolutely. The defect is identically one if and only if $\zeta$ has no zeros with $\operatorname{Re}s>\alpha$, and then $T_\alpha\overset{d}{=}T_\alpha^+$. Otherwise $D_\alpha(u)\to0$ as $u\to\infty$, and $T_\alpha^+$ is not the sum of $T_\alpha$ and an independent nonnegative variable.

Source: [`ch/almost.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/almost.tex), label `ag:prop:factor`.

<a id="explicit_defect_bounds_from_rh_verification"></a>

### Explicit bounds on the defect from the verification of RH

Theorem 11.4.4, page 543 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, let $\rho=\beta+i\gamma$ run over its zeros (the nontrivial zeros of $\zeta$), counted with multiplicity, let $P_w(x)=w/(\pi(w^2+x^2))$ be the Cauchy kernel of width $w>0$, with $P_0$ the unit point mass at $0$, and for $\alpha>\frac{1}{2}$ let $v_\alpha(\theta)=\operatorname{Re}(\xi'/\xi)(\alpha+i\theta)=v_\alpha^+(\theta)-v_\alpha^-(\theta)$ with $v_\alpha^+(\theta)/\pi=\sum_{\beta\le\alpha,\,\gamma>0}[P_{\alpha-\beta}(\theta-\gamma)+P_{\alpha-\beta}(\theta+\gamma)]$ and $v_\alpha^-(\theta)/\pi=\sum_{\beta>\alpha,\,\gamma>0}[P_{\beta-\alpha}(\theta-\gamma)+P_{\beta-\alpha}(\theta+\gamma)]$. The shifted clock $T_\alpha$ has Laplace transform $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt{u})=\varphi_\alpha^+(u)D_\alpha(u)$, where $\varphi_\alpha^+(u)=\exp\{-\int_0^\infty\log(1+u/\theta^2)\,v_\alpha^+(\theta)\,\pi^{-1}d\theta\}$ is the Laplace transform of a generalised gamma convolution with Thorin density $v_\alpha^+/\pi$ and $D_\alpha(u)=\prod_{\beta>\alpha,\,\gamma>0}\frac{(\beta-\alpha)^2+\gamma^2}{(\sqrt{u}-\beta+\alpha)^2+\gamma^2}$ is the defect. Let $T_0=3\times10^{12}$, below which all zeros lie on the critical line (the verification of Platt and Trudgian), and let $\Sigma_0=\sum_{\gamma>T_0}\gamma^{-2}$, the sum over all zeros above $T_0$ counted with multiplicity. By the Riemann von Mangoldt formula with the explicit error term of Hasanalizade, Shen and Wong, $\Sigma_0$ differs from $\frac{\log(T_0/2\pi)+1}{2\pi T_0}=1.4797\times10^{-12}$ by less than $2\times10^{-24}$, so $\Sigma_0<1.48\times10^{-12}$. Then for every $\alpha>\frac{1}{2}$ and every $u\ge0$: (1) $|\log D_\alpha(u)|\le(u+\sqrt{u})\,\Sigma_0$, so that $e^{-1.5\times10^{-12}(u+\sqrt{u})}\varphi_\alpha^+(u)\le\mathbb{E}e^{-uT_\alpha}\le e^{1.5\times10^{-12}(u+\sqrt{u})}\varphi_\alpha^+(u)$; (2) the negative part of the density satisfies $v_\alpha^-(\theta)/\pi\le2\Sigma_0/\pi<9.43\times10^{-13}$ for $0\le\theta\le T_0/2$; (3) consequently $v_\alpha(\theta)>0$ for all $0\le\theta\le T_0/2$ whenever $\alpha-\frac{1}{2}>6\times10^{-10}$. Part of the proof is certified in ball arithmetic.

Source: [`ch/almost.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/almost.tex), label `ag:thm:explicit`.

<a id="thorin_density_positive_on_verified_window"></a>

### Positivity of the Thorin density on the verified window

Corollary 11.4.5, page 545 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac{1}{2}s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, let $v_\alpha(\theta)=\operatorname{Re}(\xi'/\xi)(\alpha+i\theta)$ for $\alpha>\frac{1}{2}$ (the function whose quotient by $\pi$ is the Thorin density of the shifted clock $T_\alpha$, the positive random variable with $\mathbb{E}e^{-uT_\alpha}=\xi(\alpha)/\xi(\alpha+\sqrt{u})$, in the chart $\theta=\sqrt{z}$, when every zero lies to the left of the line $\operatorname{Re}s=\alpha$), let $T_0=3\times10^{12}$ (every zero off the critical line has $|\gamma|>T_0$, by the verification of Platt and Trudgian) and let $V_\xi=(\log\xi)''(\frac{1}{2})=0.04620998623$. (1) For every $\alpha>\frac{1}{2}$, $v_\alpha(\theta)>0$ for $|\theta|\le T_0-\frac{1}{2}$. (2) For every $\alpha\in(\frac{1}{2},7.35]$, $v_\alpha(\theta)>v_\alpha(0)$ for $0<|\theta|\le T_0-\frac{1+\sqrt{2}}{2}$, and $0<v_\alpha(0)\le(\alpha-\frac{1}{2})V_\xi$.

Source: [`ch/almost.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/almost.tex), label `cor:dbeta:window`.

<a id="size_biased_clock_decomposition_exponential_mixture_iff_rh"></a>

### Size biased clock and the exponential mixture form of RH

Proposition 11.5.1, page 560 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function and $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, an entire function of $u$ whose zeros are the numbers $-a_k$, $a_k=-(\rho_k-\tfrac12)^2$, one for each pair $\{\rho,1-\rho\}$ of zeros of $\xi$ counted with multiplicity (so $a_k=\gamma_k^2>0$ for a zero $\rho_k=\tfrac12+i\gamma_k$ on the critical line). With no hypothesis $1/G$ is the Laplace transform of a probability law on $[0,\infty)$, the law of the clock $T$, so that $\mathbb{E}e^{-uT}=1/G(u)$ for $u\ge0$. Let $W(t)=\sum_ke^{-a_kt}$, $t>0$, be the heat trace, let $m_0=\mathbb{E}T=\sum_ka_k^{-1}=0.0231049931\ldots$, and let $\widehat T$ be the size biased clock, the variable with law $t\,\mathbb{P}(T\in dt)/m_0$.
With no hypothesis, $\widehat T\overset{d}{=}T+A$, where $A$ is independent of $T$ and has density $W(t)/m_0$ on $(0,\infty)$. The Riemann Hypothesis (RH) holds if and only if $A$ is a mixture of exponential laws, that is, $A\overset{d}{=}E/V$ with $E$ a standard exponential variable independent of a variable $V>0$; and then $V$ has the law $m_0^{-1}\sum_{\gamma>0}\gamma^{-2}\delta_{\gamma^2}$, the ordinates $\gamma$ of the zeros $\tfrac12+i\gamma$ counted with multiplicity.

Source: [`ch/pakes.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes.tex), label `pk:prop:sizebias`.

<a id="centre_clock_urbanik_class_l2_not_l3"></a>

### Centre clock lies in Urbanik class L2 and not in L3

Theorem 11.5.8, page 569 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function and $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, an entire function of $u$ whose zeros are the numbers $-a_k$, $a_k=-(\rho_k-\tfrac12)^2$, one for each pair $\{\rho,1-\rho\}$ of zeros of $\xi$ counted with multiplicity. With no hypothesis $1/G$ is the Laplace transform of a probability law on $[0,\infty)$, the law of the clock $T$ (the centre clock), and $W(t)=\sum_ke^{-a_kt}$, $t>0$, is the heat trace. The ladder functions are $W^{[0]}=W$ and $W^{[j+1]}(t)=-t\,\frac{d}{dt}W^{[j]}(t)$. With no hypothesis the clock is self decomposable: $T\overset{d}{=}\int_0^\infty e^{-t}\,dZ_t$ for a subordinator $Z=Z^{[1]}$, its driving process, with Laplace exponent $u(\log G)'(u)$ and Lévy density $-W'(x)$; and $Z^{[j+1]}$ denotes the driving process of the law of $Z^{[j]}_1$ whenever that law is self decomposable. A law $\mu$ on $\mathbb{R}$ is self decomposable if for every $c\in(0,1)$ there is a law $\mu_c$ with $\widehat\mu(z)=\widehat\mu(cz)\widehat{\mu_c}(z)$ for the characteristic functions; Urbanik's class $L_0$ consists of the self decomposable laws and, for $n\ge1$, $L_n$ consists of the laws in $L_{n-1}$ all of whose components $\mu_c$ lie in $L_{n-1}$.
Unconditionally, that is, using only the verification of the Riemann Hypothesis below height $2197.3$ (Platt and Trudgian): $W^{[3]}(t)>0$ for every $t>0$, and $W^{[4]}(e^{-4.15})<0$ (in ball arithmetic $W^{[4]}(e^{-4.15})=-0.11089729616255\pm10^{-23}$). Consequently $W$, $W^{[1]}$ and $W^{[2]}$ are positive on $(0,\infty)$; the centre clock lies in $L_2$ and not in $L_3$; its driving process $Z=Z^{[1]}$ and the driving process $Z^{[2]}$ of $Z$ are self decomposable at every time (the laws of $Z_s$ and of $Z^{[2]}_s$ are self decomposable for every $s>0$); and the driving process $Z^{[3]}$ of $Z^{[2]}$ is self decomposable at no time.

Source: [`ch/pakes2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes2.tex), label `pk:thm:ladder`.

<a id="clock_first_two_driving_processes_self_decomposable_unimodal"></a>

### First two driving processes of the clock are self decomposable and unimodal

Corollary 11.5.10, page 572 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\tfrac12 s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be the completed Riemann zeta function, $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, and $T$ the clock, the random variable with $\mathbb{E}e^{-uT}=1/G(u)$ for $u\ge0$ (it exists with no hypothesis). With no hypothesis $T$ is self decomposable: $T\overset{d}{=}\int_0^\infty e^{-t}\,dZ_t$ for a subordinator $Z$, its driving process, with Laplace exponent $u(\log G)'(u)$; and $Z^{[2]}$ denotes the driving process of the self decomposable law of $Z_1$. A law $\mu$ on $\mathbb{R}$ is self decomposable if for every $c\in(0,1)$ there is a law $\mu_c$ with $\widehat\mu(z)=\widehat\mu(cz)\widehat{\mu_c}(z)$ for the characteristic functions; Urbanik's class $L_1$ consists of the self decomposable laws all of whose components $\mu_c$ are self decomposable.
With no hypothesis, for every $s>0$ the laws of $Z_s$ and of $Z^{[2]}_s$ are self decomposable, hence unimodal, and the law of $Z_s$ lies in $L_1$.

Source: [`ch/pakes2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes2.tex), label `pk:cor:unimodal`.

<a id="clock_not_supercritical_csbp_martingale_limit"></a>

### The clock is not the martingale limit of a supercritical continuous state branching process

Theorem 11.5.13, page 577 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function, $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, and $T$ the clock, $\mathbb{E}e^{-uT}=1/G(u)$ for $u\ge0$. With no hypothesis $T\overset{d}{=}\int_0^\infty e^{-t}\,dZ_t$ for a subordinator $Z$ (its driving process) with Laplace exponent $\Psi_Z(u)=u(\log G)'(u)$. Let $H$ be the inverse function of $\log G$ near $0$ and $g=\Psi_Z\circ H=H/H'=\sum_{k\ge1}g_kw^k$, so that $g_1=1$. A supercritical continuous state branching process with finite Malthusian parameter is obtained, by Lamperti's time change, from a spectrally positive Lévy process $Y$ with Laplace exponent $\Psi_Y(u)=-\log\mathbb{E}e^{-uY_1}$ and $\Psi_Y'(0)\in(0,\infty)$; started at $z>0$ it is a process $(\mathsf{B}_t)_{t\ge0}$ with $\mathbb{E}\mathsf{B}_t=ze^{\Psi_Y'(0)t}$, and the nonnegative martingale $e^{-\Psi_Y'(0)t}\mathsf{B}_t$ converges almost surely to a limit $\mathsf{M}$.
With no hypothesis, $(-1)^{15}g_{16}<0$: in ball arithmetic $g_{16}=1.49524474945112\times10^{-15}\pm2\times10^{-30}$. Consequently, for no $\lambda>0$, no such $Y$ and no initial state $z>0$ does $\lambda T$ have the law of $\mathsf{M}$: no positive multiple of the centre clock is the martingale limit of a supercritical continuous state branching process with finite Malthusian parameter, from any initial state.

Source: [`ch/pakes2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes2.tex), label `pk:thm:nobranch`.

<a id="criticality_filter_damped_polya_kernel"></a>

### Criticality filter for properties of the Pólya kernel

Proposition 11.5.14, page 579 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $\Phi(v)=2\sum_{n\ge1}(2\pi^2n^4e^{9v/2}-3\pi n^2e^{5v/2})\exp(-\pi n^2e^{2v})$ the Pólya kernel, an even positive function with $\xi(\tfrac12+i\theta)=\int_{\mathbb{R}}\Phi(v)\cos(\theta v)\,dv$; for $t\le0$ put $\Phi_t(v)=e^{tv^2/4}\Phi(v)$. For an even positive integrable function $g$ on $\mathbb{R}$ put $\widehat g(\theta)=\int_{\mathbb{R}}g(v)\cos(\theta v)\,dv$ and $R_g(s)=\widehat g(0)/\widehat g(i\sqrt s)=\int_{\mathbb{R}}g(v)\,dv\big/\int_{\mathbb{R}}g(v)\cosh(\sqrt s\,v)\,dv$ for the $s\ge0$ at which the last integral converges; for $g=\Phi$ this is $R_\Phi(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$, the Laplace transform of the clock. A GGC (generalised, or generalized, gamma convolution) is a random variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, where the drift $d\ge0$ and the Thorin measure $U$, a positive measure on $(0,\infty)$, satisfy $\int_{(0,1)}|\log z|\,U(dz)+\int_{(1,\infty)}z^{-1}U(dz)<\infty$.
Let $\mathcal{P}$ be a property of even positive integrable functions on $\mathbb{R}$ that $\Phi_t$ has for some $t<0$. Then $\mathcal{P}$ implies neither that $\widehat g$ has only real zeros nor that $R_g$ is the Laplace transform of a GGC, for the functions $g$ with $\mathcal{P}$: the function $g=\Phi_t$ has $\mathcal{P}$, while $\widehat g$ has nonreal zeros and $R_g$ is not the Laplace transform of a GGC. At every $t\le0$ the kernel $\Phi_t$ (i) is positive and even, and strictly decreasing on $(0,\infty)$; (ii) is log concave; (iii) satisfies $\log\Phi_t(v)=-\pi e^{2|v|}+O(v^2)$ as $|v|\to\infty$; (iv) extends analytically to the strip $|\operatorname{Im}v|<\pi/4$. So no combination of (i) to (iv) implies that $\widehat g$ has only real zeros, or that $R_g$ is the Laplace transform of a GGC. The proof uses the theorem of Rodgers and Tao that the de Bruijn Newman constant is nonnegative, which gives a nonreal zero of $\widehat{\Phi_t}$ for every $t<0$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prop:filter`.

<a id="damped_polya_law_arcsine_mixture"></a>

### Damped Pólya law as a mixture of arcsine laws

Proposition 11.5.15, page 581 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $\Phi(v)=2\sum_{n\ge1}(2\pi^2n^4e^{9v/2}-3\pi n^2e^{5v/2})\exp(-\pi n^2e^{2v})$ the Pólya kernel, an even positive function with $\int_{\mathbb{R}}\Phi(v)\cos(\theta v)\,dv=\xi(\tfrac12+i\theta)$. For $t\le0$ let $\Phi_t(v)=e^{tv^2/4}\Phi(v)$, let $p_t=\Phi_t/\int_{\mathbb{R}}\Phi_t$ be the corresponding probability density, which is even and decreasing on $(0,\infty)$, and let $X_t$ be a variable with density $p_t$ (so $X_0$ is the Pólya variable).
For every $t\le0$ the function $h_t(\rho)=-2\rho\int_\rho^\infty p_t'(r)\,dr/\sqrt{r^2-\rho^2}$, $\rho>0$, is a positive probability density, and $p_t(v)=\frac1\pi\int_{|v|}^\infty h_t(\rho)\,d\rho/\sqrt{\rho^2-v^2}$ for $v\ne0$. Equivalently, $X_t\overset{d}{=}\varrho_t\cos\Theta$ with $\Theta$ uniform on $[0,2\pi)$ and $\varrho_t$ of density $h_t$ independent of $\Theta$: $X_t$ is a mixture of arcsine laws, the first coordinate of the rotation invariant vector $\varrho_t(\cos\Theta,\sin\Theta)$ in $\mathbb{R}^2$ and also of a rotation invariant vector in $\mathbb{R}^3$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prop:arcsine`.

<a id="biane_pitman_yor_lognormal_tilt_not_ggc_near_zero"></a>

### Tilted Biane Pitman Yor variable is not a GGC near t = 0 and at ten certified values of t

Proposition 11.5.17, page 583 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be the completed Riemann zeta function and $f(b)=\sum_{m\ge1}(4\pi^2m^4b-6\pi m^2)e^{-\pi m^2b}$, $b>0$, the density of the variable $\mathcal{B}$ of Biane, Pitman and Yor, with $\mathbb{E}e^{-\lambda\mathcal{B}}=(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^2$ and $\mathbb{E}\mathcal{B}^{s/2}=2\xi(s)$; $\mathcal{B}$ is a GGC with Thorin measure $2\sum_{m\ge1}\delta_{\pi m^2}$. For real $t$ let $\mathcal{B}_t$ be the variable with density proportional to $e^{t(\log b)^2/16}f(b)$ on $(0,\infty)$, so that $\mathcal{B}_0=\mathcal{B}$. A GGC (generalised, or generalized, gamma convolution) is a random variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, with drift $d\ge0$ and Thorin measure $U$, a positive measure on $(0,\infty)$ with $\int_{(0,1)}|\log z|\,U(dz)+\int_{(1,\infty)}z^{-1}U(dz)<\infty$.
(i) $\mathcal{B}_t$ is not a GGC for $t=-20$, $-5$, $-1$, $-\tfrac12$, $-\tfrac1{10}$, $-\tfrac1{100}$, $\tfrac1{100}$, $\tfrac15$, $1$ and $5$. (ii) There is $t_0>0$ such that $\mathcal{B}_t$ is not a GGC for $0<|t|<t_0$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prop:bpytilt`.

<a id="theta_variable_character_minus_8_not_ggc"></a>

### Theta variable of the odd real character modulo 8 and its reciprocal are not GGCs

Proposition 11.5.20, page 587 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: counterexample, proved.

Statement. Let $\chi_{-8}$ be the real primitive Dirichlet character modulo $8$ with $\chi_{-8}(-1)=-1$ (so $\chi_{-8}(n)=1$ for $n\equiv1,3$, $\chi_{-8}(n)=-1$ for $n\equiv5,7 \pmod{8}$, and $\chi_{-8}(n)=0$ for even $n$), let $\vartheta(x)=\sum_{n\ge1}n\chi_{-8}(n)e^{-\pi n^2x/8}$ for $x>0$, a positive function with $\vartheta(1/x)=x^{3/2}\vartheta(x)$, and $\xi(s,\chi_{-8})=\int_0^\infty\vartheta(x)x^{(s+1)/2}\,dx/x=(8/\pi)^{(s+1)/2}\Gamma(\tfrac{s+1}2)L(s,\chi_{-8})$. Let $\mathcal{B}_{\chi_{-8}}$ be the variable with density $b^{-1/2}\vartheta(b)/\xi(0,\chi_{-8})$ on $(0,\infty)$, so that $\mathbb{E}\mathcal{B}_{\chi_{-8}}^{s/2}=\xi(s,\chi_{-8})/\xi(0,\chi_{-8})$; its reciprocal $1/\mathcal{B}_{\chi_{-8}}$ has density $\vartheta(x)/\xi(0,\chi_{-8})$. A GGC (generalised, or generalized, gamma convolution) is a random variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, with drift $d\ge0$ and Thorin measure $U$, a positive measure on $(0,\infty)$ with $\int_{(0,1)}|\log z|\,U(dz)+\int_{(1,\infty)}z^{-1}U(dz)<\infty$.
Neither $\mathcal{B}_{\chi_{-8}}$ nor $1/\mathcal{B}_{\chi_{-8}}$ is a GGC. (The proof is a certificate in ball arithmetic: the cumulants $\kappa_2,\kappa_3,\kappa_4$ of a GGC with a finite fourth moment satisfy $2\kappa_2\kappa_4\ge3\kappa_3^2$, and $3\kappa_3^2-2\kappa_2\kappa_4$ equals $7.55465911093\pm4\times10^{-12}$ for $\mathcal{B}_{\chi_{-8}}$ and $78.8854168657\pm2\times10^{-11}$ for $1/\mathcal{B}_{\chi_{-8}}$, both positive.)

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prop:chi8`.

<a id="theta_odd_character_reciprocal_laplace_transform"></a>

### Reciprocal theta variable of odd characters is a GGC for discriminants -3 and -4, not for -8, for d congruent to 1 mod 8 and for the function of Davenport and Heilbronn

Proposition 11.5.21, page 588 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $q\ge3$ and let $a=(a_n)_{n\ge1}$ be a real sequence of period $q$ with $a_1=1$ and $a_{q - n}=-a_n$, whose theta series $\vartheta_a(x)=\sum_{n\ge1}na_ne^{-\pi n^2x/q}$, $x>0$, satisfies $\vartheta_a(1/x)=x^{3/2}\vartheta_a(x)$ and is positive. Put $\xi(s,a)=\int_0^\infty\vartheta_a(x)x^{(s+1)/2}\,dx/x$, let $\mathcal{B}_a$ be the variable with density $b^{-1/2}\vartheta_a(b)/\xi(0,a)$ on $(0,\infty)$, so that $1/\mathcal{B}_a$ has density $\vartheta_a(x)/\xi(0,a)$, and for $\lambda>0$ put $w=\sqrt{\pi q\lambda}$ and $\Sigma_a(w)=\sum_{1\le n<q/2}a_n\sinh\frac{(q - 2n)w}{q}$. A GGC (generalised, or generalized, gamma convolution) is a random variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, with drift $d\ge0$ and Thorin measure $U$ (a positive measure on $(0,\infty)$ with $\int_{(0,1)}|\log z|\,U(dz)+\int_{(1,\infty)}z^{-1}U(dz)<\infty$); the Laplace transform of a GGC extends holomorphically and without zeros to $\mathbb{C}\setminus(-\infty,0]$ and cannot vanish on $(-\infty,0)$ where it extends meromorphically.
Then $\int_0^\infty e^{-\lambda x}\vartheta_a(x)\,dx=\sqrt q\,\Sigma_a(w)/\sinh w$ and $\mathbb{E}e^{-\lambda/\mathcal{B}_a}=\sqrt q\,\Sigma_a(w)/(\xi(0,a)\sinh w)$, and $\xi(0,a)=q^{-1/2}\sum_{1\le n<q/2}(q - 2n)a_n$; the first formula is formula (96) of Biane, Pitman and Yor in this normalisation. In particular: (i) for the character of discriminant $-4$ ($q=4$, $a=(1,0,-1,0)$), $\mathbb{E}e^{-\lambda/\mathcal{B}_a}=\operatorname{sech}\sqrt{\pi\lambda}$, and $1/\mathcal{B}_a$, which is $\tfrac\pi2$ times the exit time of Brownian motion from $(-1,1)$, is a GGC; (ii) for the character of discriminant $-3$ ($q=3$, $a=(1,-1,0)$), $\mathbb{E}e^{-\lambda/\mathcal{B}_a}=3\sinh(w/3)/\sinh w=\prod_{3\nmid m}(1+3\lambda/(\pi m^2))^{-1}$, so that $1/\mathcal{B}_a$ is a GGC with Thorin measure $\sum_{3\nmid m}\delta_{\pi m^2/3}$, and it is $\tfrac{3\pi}2$ times the hitting time of $1$ by a three dimensional Bessel process started at $\tfrac13$; (iii) for the character of discriminant $-8$ ($q=8$, $a=(1,0,1,0,-1,0,-1,0)$), $\mathbb{E}e^{-\lambda/\mathcal{B}_a}=\cosh(w/4)/\cosh(w/2)=\prod_{k\ge1}(1+\lambda/(4\alpha_k))/(1+\lambda/\alpha_k)$ with $\alpha_k=\tfrac\pi8(2k-1)^2$, so that $1/\mathcal{B}_a$ is infinitely divisible with Lévy density $\sum_k(e^{-\alpha_kx}-e^{-4\alpha_kx})/x>0$, and it is $\pi$ times the hitting time of $1$ by a reflected Brownian motion started at $\tfrac12$ (the identifications of (ii) and (iii) with hitting times are those of Biane, Pitman and Yor); it is not a GGC, because its Laplace transform, continued meromorphically, vanishes at $\lambda=-\pi/2$; (iv) if $a=\chi_d$ for a real primitive character of fundamental discriminant $d<0$ with $\chi_d(2)=1$, that is $d\equiv1 \pmod{8}$, and $q=|d|$, the Laplace transform of $1/\mathcal{B}_a$, continued meromorphically by the closed form, vanishes at a point of $(-4\pi/q,-\pi/q)$, so $1/\mathcal{B}_a$ is not a GGC; for $\chi_{-7}$ the point is $\lambda=-1.033610573$. The same holds for the function of Davenport and Heilbronn, with $q=5$, $a=(1,\varkappa,-\varkappa,-1,0)$ and $\varkappa=(\sqrt{10-2\sqrt5}-2)/(\sqrt5-1)$, at $\lambda=-2.046575671$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prop:reflected`.

<a id="open_ggc_class_theta_functions_zeros_on_line"></a>

### Open problem on theta type functions whose GGC property forces zeros onto the line

Problem 11.1, page 591 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $q\ge3$, $p\in\{0,1\}$ and let $a=(a_n)_{n\ge1}$ be a real sequence of period $q$ with $a_1=1$ and $a_{q - n}=(-1)^pa_n$, whose theta series $\vartheta_a(x)=\sum_{n\ge1}n^pa_ne^{-\pi n^2x/q}$, $x>0$, satisfies $\vartheta_a(1/x)=x^{p+1/2}\vartheta_a(x)$; put $\xi(s,a)=\int_0^\infty\vartheta_a(x)x^{(s+p)/2}\,dx/x=(q/\pi)^{(s+p)/2}\Gamma(\tfrac{s+p}2)\sum_{n\ge1}a_nn^{-s}$, an entire function with $\xi(s,a)=\xi(1-s,a)$. When $\vartheta_a>0$ (a positive kernel) let $\mathcal{B}_a$ be the variable with density $b^{p/2-1}\vartheta_a(b)/\xi(0,a)$ on $(0,\infty)$, so that $\mathbb{E}\mathcal{B}_a^{s/2}=\xi(s,a)/\xi(0,a)$. For the completed Riemann zeta function $\xi$ the corresponding variable is $\mathcal{B}$ with $\mathbb{E}e^{-\lambda\mathcal{B}}=(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^2$ and $\mathbb{E}\mathcal{B}^{s/2}=2\xi(s)$, which is a GGC (generalised, or generalized, gamma convolution: a variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$, $u>0$, for a drift $d\ge0$ and a positive measure $U$ on $(0,\infty)$, the Thorin measure). The test by Hankel matrices of cumulants: for a positive variable $\mathcal{V}$ with $\ell(\lambda)=\mathbb{E}e^{-\lambda\mathcal{V}}$ put $\kappa_k(\lambda)=(-1)^k\frac{d^k}{d\lambda^k}\log\ell(\lambda)$ and $c_k(\lambda)=\kappa_{k+1}(\lambda)/k!$; if $\mathcal{V}$ is a GGC, the Hankel matrices $\mathcal{K}_n(\lambda)=[c_{i+j}(\lambda)]_{i,j=0}^{n-1}$ and $\mathcal{K}'_n(\lambda)=[c_{i+j+1}(\lambda)]_{i,j=0}^{n-1}$ are positive semidefinite for all $n$ and $\lambda\ge0$, so a variable for which one of them has a negative eigenvalue is not a GGC, and is said to fail the test (applied numerically for $n\le12$ at $\lambda=0$, $\tfrac12$, $2$, $10$, $20$, $50$ and $200$).
Question. Find a natural class of functions $\xi(s,a)$ with positive kernels, containing $\xi$, in which the GGC property of $\mathcal{B}_a$ implies that every zero of $\xi(s,a)$ lies on the critical line $\operatorname{Re}s=\tfrac12$, or find a function in such a class whose variable $\mathcal{B}_a$ is a GGC and which has a zero off the line. Since $\mathcal{B}$ is a GGC, the first answer for a class containing $\xi$ proves the Riemann Hypothesis. No candidate for the second is known: the function of Davenport and Heilbronn (coefficients of period $5$, $(1,\varkappa,-\varkappa,-1,0)$ with $\varkappa=(\sqrt{10-2\sqrt5}-2)/(\sqrt5-1)$), the members with a positive kernel, among the $17$ tested, of the real one parameter family $f_r$ of combinations of $L(s,\chi_{-7})$, $L(s,\chi)$ and $L(s,\bar\chi)$, with $\chi$ and $\bar\chi$ the complex odd characters modulo $7$, that have real coefficients, $a_1=1$ and the functional equation of the class, and the tilts $\mathcal{B}_t$ of $\mathcal{B}$ with density proportional to $e^{t(\log b)^2/16}f(b)$ at the eleven values of $t<0$ tested ($f$ the density of $\mathcal{B}$) all fail the test by Hankel matrices of cumulants.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prob:ggc`.

<a id="open_biane_pitman_yor_tilt_ggc_for_no_t"></a>

### Open problem, is the tilted Biane Pitman Yor variable a GGC for no t other than 0

Problem 11.2, page 592 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\mathcal{B}$ be the variable of Biane, Pitman and Yor, a GGC (generalised, or generalized, gamma convolution) with density $f(b)=\sum_{m\ge1}(4\pi^2m^4b-6\pi m^2)e^{-\pi m^2b}$, $b>0$, $\mathbb{E}e^{-\lambda\mathcal{B}}=(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^2$ and $\mathbb{E}\mathcal{B}^{s/2}=2\xi(s)$ for the completed Riemann zeta function $\xi$, and for real $t$ let $\mathcal{B}_t$ be the variable with density proportional to $e^{t(\log b)^2/16}f(b)$ on $(0,\infty)$, so that $\mathcal{B}_0=\mathcal{B}$. Let $\Phi(v)=2\sum_{n\ge1}(2\pi^2n^4e^{9v/2}-3\pi n^2e^{5v/2})\exp(-\pi n^2e^{2v})$ be the Pólya kernel, an even positive function with $\xi(\tfrac12+i\theta)=\int_{\mathbb{R}}\Phi(v)\cos(\theta v)\,dv$, so that $f(b)=b^{-5/4}\Phi(\tfrac12\log b)$, and for real $t$ let $\Phi_t(v)=e^{tv^2/4}\Phi(v)$ be the damped kernel; then $\mathcal{B}_t$ has density proportional to $b^{-5/4}\Phi_t(\tfrac12\log b)$.
Question. Is $\mathcal{B}_t$ a GGC for no $t\ne0$? For $t<0$ a positive answer would mean that the GGC property of the variable $\mathcal{B}_t$ built from the damped Pólya kernel $\Phi_t(v)=e^{tv^2/4}\Phi(v)$ is lost at every $t<0$, although the Pólya kernel $\Phi$ itself has it ($\mathcal{B}_0=\mathcal{B}$ is a GGC). Known: $\mathcal{B}_t$ is not a GGC for $0<|t|<t_0$ for some $t_0>0$, and for $t=-20$, $-5$, $-1$, $-\tfrac12$, $-\tfrac1{10}$, $-\tfrac1{100}$, $\tfrac1{100}$, $\tfrac15$, $1$ and $5$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prob:tilt`.

<a id="open_ggc_theta_variables_of_passing_characters"></a>

### Open problem on the GGC property of theta variables of characters

Problem 11.3, page 592 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For a real primitive Dirichlet character $\chi_d$ of fundamental discriminant $d$ ($p=1$ for $d<0$, $p=0$ for $d>0$) with positive theta series $\vartheta_{\chi_d}(x)=\sum_{n\ge1}n^p\chi_d(n)e^{-\pi n^2x/|d|}$, let $\mathcal{B}_{\chi_d}$ be the variable with density proportional to $b^{p/2-1}\vartheta_{\chi_d}(b)$ on $(0,\infty)$. A GGC (generalised, or generalized, gamma convolution) is a random variable $H\ge0$ with $\mathbb{E}e^{-uH}=\exp\{-du-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, with drift $d\ge0$ and a positive Thorin measure $U$ on $(0,\infty)$. For a positive variable $\mathcal{V}$ with $\ell(\lambda)=\mathbb{E}e^{-\lambda\mathcal{V}}$ put $\kappa_k(\lambda)=(-1)^k\frac{d^k}{d\lambda^k}\log\ell(\lambda)$ and $c_k(\lambda)=\kappa_{k+1}(\lambda)/k!$; if $\mathcal{V}$ is a GGC, the Hankel matrices $\mathcal{K}_n(\lambda)=[c_{i+j}(\lambda)]_{i,j=0}^{n-1}$ and $\mathcal{K}'_n(\lambda)=[c_{i+j+1}(\lambda)]_{i,j=0}^{n-1}$ are positive semidefinite for all $n$ and $\lambda\ge0$, so a negative eigenvalue shows that $\mathcal{V}$ is not a GGC. A variable passes the test to order $12$ if no $\mathcal{K}_n(\lambda)$ or $\mathcal{K}'_n(\lambda)$ with $n\le12$ has a negative eigenvalue at $\lambda=0$, $\tfrac12$, $2$, $10$, $20$, $50$ and $200$ (a numerical test).
Question. Which of the variables $\mathcal{B}_{\chi_d}$ for $d=-4$, $-3$, $5$, $8$ and $12$, which pass the test to order $12$, are GGCs? (For $d=24$ the test fails at order $11$.) On the reciprocal side, $1/\mathcal{B}_{\chi_d}$ is a GGC for $d=-3$ and $d=-4$; for $d=-4$ it is $\tfrac\pi2$ times the exit time of Brownian motion from $(-1,1)$, whose density is not hyperbolically completely monotone (HCM), so the question for $\mathcal{B}_{\chi_{-4}}$ cannot be settled through HCM. For even characters the Laplace transform of the reciprocal variable is the series $\frac2{\xi(0,a)}\sum_na_nK_0\big(2n\sqrt{\pi\lambda/q}\big)$ of Bessel functions $K_0$ (with $q=|d|$, $a_n=\chi_d(n)$ and $\xi(0,a)$ the normalising constant of the density of $\mathcal{B}_{\chi_d}$), and no product formula for it is known to the authors.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prob:side`.

<a id="open_li_coefficients_unconditional_positivity_size_biased_clock"></a>

### Open problem on unconditional positivity of the Li coefficients

Problem 11.4, page 592 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be the completed Riemann zeta function, $G(u)=\xi(\tfrac12+\sqrt u)/\xi(\tfrac12)$, with zeros $-a_k$, $a_k=-(\rho_k-\tfrac12)^2$, one for each pair $\{\rho,1-\rho\}$ of zeros of $\xi$ counted with multiplicity. Let $W(t)=\sum_ke^{-a_kt}$ be the heat trace, $m_0=\sum_ka_k^{-1}$, and $A$ a variable with density $W(t)/m_0$ on $(0,\infty)$ (the size biased clock is the clock plus an independent copy of $A$). The Li coefficients are $\lambda_n=\sum_\rho\big(1-(1-1/\rho)^n\big)$, $n\ge1$, the sum over the zeros $\rho$ of $\xi$ grouped in pairs $\{\rho,1-\rho\}$; they satisfy $\lambda_n=\sum_k\varphi_n(a_k)$ with $\varphi_n(z)=2-\big(\frac{2\sqrt z+i}{2\sqrt z-i}\big)^n-\big(\frac{2\sqrt z-i}{2\sqrt z+i}\big)^n$, a rational function of $z$ with a pole of order $n$ at $z=-\tfrac14$ and no constant term in $(z+\tfrac14)^{-1}$ (for example $\varphi_1(z)=(z+\tfrac14)^{-1}$ and $\varphi_2(z)=4z(z+\tfrac14)^{-2}$). Let $k_n$ be the function with $\int_0^\infty e^{-zt}k_n(t)\,dt=\varphi_n(z)$; it has the form $k_n(t)=e^{-t/4}p_n(t)$ with $p_n$ a polynomial of degree $n-1$ (for example $k_1(t)=e^{-t/4}$ and $k_2(t)=(4-t)e^{-t/4}$).
With no hypothesis, $\lambda_n=\int_0^\infty k_n(t)W(t)\,dt=m_0\,\mathbb{E}k_n(A)$. Question. Which unconditional properties of $W$, or of $A$, force $\mathbb{E}k_n(A)\ge0$ for every $n$, which is the Riemann Hypothesis by Li's criterion? Positivity of $W$ gives $n=1$; complete monotonicity of $W$, which is equivalent to RH, gives every $n$, since $\varphi_n\ge0$ on $(0,\infty)$.

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prob:li`.

<a id="open_ladder_index_along_de_bruijn_newman_flow"></a>

### Open problem on the ladder index along the de Bruijn Newman flow

Problem 11.5, page 592 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\Phi$ be the Pólya kernel, with $\xi(\tfrac12+i\theta)=\int_{\mathbb{R}}\Phi(v)\cos(\theta v)\,dv$ for the completed Riemann zeta function $\xi$, let $\Phi_t(v)=e^{tv^2/4}\Phi(v)$, $\Xi_t(\theta)=\int_{\mathbb{R}}\Phi_t(v)\cos(\theta v)\,dv$ and $R_t(s)=\Xi_t(0)/\Xi_t(i\sqrt s)$, so that $R_0(s)=\xi(\tfrac12)/\xi(\tfrac12+\sqrt s)$ is the Laplace transform of the clock. The zeros of $\Xi_t$ are all real exactly when $t\ge\Lambda$, where $\Lambda$ is the de Bruijn Newman constant, $0\le\Lambda\le0.2$ (Rodgers and Tao; Platt and Trudgian). For $t\ge\Lambda$, in particular unconditionally for every $t\ge0.2$, let $\pm\theta_k(t)$, $0<\theta_1(t)<\theta_2(t)<\cdots$, be the zeros of $\Xi_t$ and $a_k(t)=\theta_k(t)^2$; then $R_t(s)=\prod_k(1+s/a_k(t))^{-1}$ is the Laplace transform of a GGC (generalised, or generalized, gamma convolution) with Thorin measure $\sum_k\delta_{a_k(t)}$. At $t=0$ put $a_k(0)=a_k=-(\rho_k-\tfrac12)^2$, one for each pair $\{\rho,1-\rho\}$ of zeros of $\xi$ (real when RH holds), the numbers of the clock. Put $W(x;t)=\sum_ke^{-a_k(t)x}$ and $W^{[j]}(x;t)=(-x\,\partial_x)^jW(x;t)$. A law with this Laplace transform lies in Urbanik's class $L_n$ exactly when $W^{[n+1]}(\cdot\,;t)\ge0$ on $(0,\infty)$ (the classes $L_n$: $L_0$ is the class of self decomposable laws and $L_n$ the class of laws in $L_{n-1}$ all of whose components lie in $L_{n-1}$); the index on the window $0\le-\log x\le8$ is the largest $n$ with $W^{[1]},\dots,W^{[n+1]}\ge0$ there.
Question. Describe the ladder index of $R_t$ along the de Bruijn Newman flow. A numerical computation (not certified for $t>0$) finds on the window the index two of the clock at the tested values of $t$ up to $27.8$ and the index one at the tested values $t=27.9$, $30$, $40$ and $50$. Is the index nonincreasing in $t$, and is it one for all large $t$?

Source: [`ch/pakes3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/pakes3.tex), label `pk:prob:ladder`.

## Chapter 12. Conclusions and open problems

<a id="open_jensen_polynomials_as_finite_free_products"></a>

### Jensen polynomials of the Pólya moments as finite free products (open problem)

Problem 12.1, page 609 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$; its moment generating function is $\mathbb{E}e^{sX}=\xi(\frac12+s)/\xi(\frac12)$. Put $b_j=\mathbb{E}X^{2j}/(2j)!$, so that $\Xi(t)/\Xi(0)=\sum_{j\ge0}(-1)^jb_jt^{2j}$ with $\Xi(t)=\xi(\frac12+it)$, and for $d,n\ge0$ let $J_{d,n}(x)=\sum_{i=0}^d\binom di(n+i)!\,b_{n+i}x^i$ be the Jensen polynomial of degree $d$ and shift $n$ of the sequence $j!\,b_j$. Pólya proved that the Riemann Hypothesis holds if and only if every $J_{d,n}$ has only real zeros. The finite free multiplicative convolution $\boxtimes_d$ of polynomials of degree $d$ multiplies their normalised coefficients (the coefficient of $x^{d - i}$ divided by $\binom di$, up to the sign $(-1)^i$) and keeps roots real and of one sign when both factors have real roots of one sign. The normalised coefficients of the reversed Jensen polynomial $x^dJ_{d,n}(1/x)$ are, up to the signs $(-1)^i$, the numbers $(n+i)!\,b_{n+i}$ for $0\le i\le d$. A factorisation $(n+i)!\,b_{n+i}=\alpha_i\beta_i$ into the normalised coefficients of two polynomials with only real roots of one sign would therefore give real zeros of $J_{d,n}$. Open problem: find such factorisations for every shift $n$ at a fixed degree $d$, or show that none exists with factors of a given kind. The obvious split, a Gaussian multiplier sequence $q^{i^2}$ with $0<q\le1$ times a residual close to the coefficients of $(1+x)^d$, fails, because a root of multiplicity $d$ splits into nonreal roots under a generic perturbation. Dividing the normalised coefficients of a polynomial with real roots of one sign by $(a)_i=a(a+1)\cdots(a+i-1)$ keeps its roots real, but multiplying them by $(a)_i$ does not (it turns $(1+x)^2$ into $1+2ax+a(a+1)x^2$, whose discriminant is $-4a<0$). A factorisation has to put real structure into both factors.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:jensen`.

<a id="open_laguerre_meixner_pollaczek_expansions_of_polya_density"></a>

### Laguerre and Meixner Pollaczek expansions of the Pólya density (open problem)

Problem 12.2, page 609 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and let $\Xi(t)=\xi(\frac12+it)$. The book proves that for no $a>0$ is $|X|$ or $X^2$ a scale mixture of Gamma$(a)$ laws, and that, under the Riemann Hypothesis, removing a gamma layer keeps zeros real while adding one does not: the function $G_a(u)=\sum_kb_ku^k/(a)_k$, with $b_k=\mathbb{E}X^{2k}/(2k)!$ and $(a)_k=a(a+1)\cdots(a+k-1)$, lies in the Laguerre Pólya class (the entire functions that are locally uniform limits of real rooted polynomials) for every $a>0$, whereas multiplying the coefficients by $(a)_k$ instead does not preserve real zeros. Open problem: find a variable in which the Pólya density, or $\Xi$ itself, has an expansion in Laguerre or Meixner Pollaczek polynomials with positive coefficients and truncations with only real zeros. The symmetric Meixner Pollaczek expansions treated in the book converge slowly, and a real rooted construction from their coefficients is a separate open problem of the chapter on the reciprocal zeta function.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:laguerre`.

<a id="open_jacobi_parameters_of_zero_measure_and_local_semicircle_model"></a>

### Jacobi parameters of the zero measure and the local semicircle model (open problem)

Problem 12.3, page 609 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\gamma_1\le\gamma_2\le\cdots$ be the positive ordinates of the zeros $\frac12\pm i\gamma_j$ of $\zeta$, put $a_j=\gamma_j^2$ and let $\nu=\sum_ja_j^{-1}\delta_{1/a_j}$, a measure whose moments are computed from the Taylor coefficients at the centre of the Xi function $\Xi(t)=\xi(\frac12+it)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ is Riemann's xi function. Under the Riemann Hypothesis $\nu$ is a positive measure and all its Jacobi parameters exist: write $\alpha_n$ for the diagonal entries and $\beta_n>0$ for the off diagonal entries of its Jacobi matrix, the recurrence coefficients of its orthogonal polynomials. A semicircle law is a local model for them: constant parameters $\alpha_n=\alpha$, $\beta_n=\beta$ give the semicircle law on $[\alpha-2\beta,\alpha+2\beta]$, and a local model with right end $L_n$ of the interval for row $n$ has $\alpha_n=\frac12L_n$ and $\beta_n=\frac14L_n$. Let $N$ be the smooth Riemann von Mangoldt count of the ordinates, $n(L)=\int_{(0,L^{-1/2})}dN(T)/\sqrt{1-LT^2}$, and define $L_n$ by $n(L_n)=n$. Open problem: under the Riemann Hypothesis, prove that $\alpha_n\sim2\beta_n\sim\frac12L_n$ as $n\to\infty$, a spectral form of the Riemann von Mangoldt formula. A table in the book supports it to a few per cent for $n\le49$, but the local model is a heuristic, outside the scaling of Kuijlaars and Van Assche (1999). Also decide between the rules $n(L_n)=n$ and $n(L_n')=n-\frac12$, which the data do not separate, and explain the scatter of the ratios, about six per cent, through the individual gaps between the zeros.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:jacobi`.

<a id="open_order_at_which_a_zero_collision_is_registered_by_hankel_test"></a>

### Order at which a zero collision is registered by the Hankel test (open problem)

Problem 12.4, page 609 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $\Xi(t)=\xi(\frac12+it)$ the Xi function, and $H_t$ the de Bruijn Newman heat flow in the normalisation of Rodgers and Tao, which the book uses: $8H_t(2\theta)=\int_{\mathbb R}\Phi(v)e^{tv^2/4}\cos(\theta v)\,dv$, where $\Phi(v)=\sum_{n\ge1}(4\pi^2n^4e^{9v/2}-6\pi n^2e^{5v/2})\exp(-\pi n^2e^{2v})$, so that $8H_0(2\theta)=\Xi(\theta)$. Write $\pm\theta_j(t)$ for the zeros of $\theta\mapsto H_t(2\theta)$, so that $\theta_j(0)=\gamma_j$, where $\frac12\pm i\gamma_j$ are the nontrivial zeros of $\zeta$, numbered by height. For every $t<0$ the function $H_t$ has nonreal zeros, by the theorem of Rodgers and Tao, and as $t$ decreases from $0$ pairs of consecutive real zeros $\theta_k(t),\theta_{k+1}(t)$, which start at $(\gamma_k,\gamma_{k+1})$, collide and leave the real axis. Let $\mathcal H_N(t)=[m_{i+j}(t)]_{i,j<N}$ be the Hankel matrix of order $N$ of the moments $m_n(t)=\sum_j\theta_j(t)^{-2n-2}$ (at $t=0$ these are $m_n(0)=\sum_k\gamma_k^{-2n-2}$). Under the Riemann Hypothesis $\mathcal H_N(0)$ is positive definite. Let $n_{\rm fail}(t)$ be the least $N$ for which $\mathcal H_N(t)$ is not positive definite and $t^*_N=\sup\{t<0:n_{\rm fail}(t)\le N\}$, the time, going back from $t=0$, at which the Hankel test of order $N$ first fails. Under the Riemann Hypothesis every $t^*_N$ is negative, $t^*_N\le t^*_{N+1}$ and $t^*_N\to0$. The book computes $t^*_N$ for $N\le50$ and explains the steps of the staircase by collisions. The rule suggested in an unpublished note of N. G. Polson, that the collision of $(\gamma_k,\gamma_{k+1})$ is registered from order about $2k+2$ on, holds for $N=30$ and $40$ and fails at $N=50$: it gives the collision of $(\gamma_{24},\gamma_{25})$ at $t=-1.144$, while the matrix of order $50$ registers $(\gamma_{27},\gamma_{28})$, for which $2k+2=56$, at $t=-1.0325$, before $(\gamma_{24},\gamma_{25})$ has collided. Open problem: find the rule that gives the order at which a pair is registered, as a function of the index $k$ of the pair and of its distance from the axis, and follow the staircase of $t^*_N$ to orders of a few hundred, towards the limit $0$ that it has under the Riemann Hypothesis.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:window`.

<a id="open_unconditional_collision_brackets_for_de_bruijn_newman_flow"></a>

### Unconditional brackets for the collision times of zeros along the de Bruijn Newman flow (open problem)

Problem 12.5, page 610 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $\Xi(t)=\xi(\frac12+it)$ the Xi function, and $H_t$ the de Bruijn Newman heat flow in the normalisation of Rodgers and Tao, which the book uses: $8H_t(2\theta)=\int_{\mathbb R}\Phi(v)e^{tv^2/4}\cos(\theta v)\,dv$, where $\Phi(v)=\sum_{n\ge1}(4\pi^2n^4e^{9v/2}-6\pi n^2e^{5v/2})\exp(-\pi n^2e^{2v})$, so that $8H_0(2\theta)=\Xi(\theta)$; as $t$ decreases below $0$, pairs of adjacent real zeros of $H_t$ collide and leave the real axis, and $H_t$ has nonreal zeros for every $t<0$ by the theorem of Rodgers and Tao. The book brackets the collision times of such pairs by quantities computed along the flow and through the nonreal zeros of $H_t$, and checks the hypotheses of the brackets on a grid of times for the $58$ collisions of a census computed in the book. Open problem: find a priori control of the nonreal zeros of $H_t$ for $t<0$, of how far they travel and how close to a gap between real zeros they come, which would make the brackets statements about $H_t$ computed from $H_0$ alone.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:collide`.

<a id="open_moving_derivative_sign_ladder_of_heat_trace_without_new_zeros"></a>

### Extending the derivative sign ladder of the heat trace without new information on the zeros (open problem)

Problem 12.6, page 610 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function, $\gamma_k$ the ordinates of its zeros $\frac12\pm i\gamma_k$ in the upper half plane, $a_k=\gamma_k^2$ and $W(t)=\sum_ke^{-ta_k}$ the heat trace of the zeros. Let $T$ be the clock, the positive random variable with Laplace transform $1/G(u)$, $G(u)=\xi(\frac12+\sqrt{u})/\xi(\frac12)$, which has Lévy measure $W(x)\,dx/x$. The Upsilon transform sends an infinitely divisible law $\mu$ with Lévy process $(X_t)$, $X_1\sim\mu$, to the law of $\int_0^1\log(1/t)\,dX_t$, and the class BO is the class of infinitely divisible laws on $[0,\infty)$ whose Lévy measure has a completely monotone density. The book proves the derivative signs $(-1)^kW^{(k)}>0$ on ranges that grow with the height to which the Riemann Hypothesis (RH) is verified, and complete monotonicity of $W$ on $(0,\infty)$ is RH. Open problem: find a mechanism that extends these ranges with no new information on the zeros, for instance a functional equation for the law of $T$ under the Upsilon transform, as suggested in an unpublished note of N. G. Polson. The book restricts the form such an equation can take: with no hypothesis $T$ is not the Upsilon image of a law in BO, and under RH $T$ is the Upsilon image of the law $\mu_0$, the self decomposable law with no drift and Lévy measure $N(x^{-1/2})\,dx/x$, where $N(T)$ is the number of ordinates $\gamma_k\le T$, which is neither a generalised gamma convolution (GGC) nor in BO.

Source: [`ch/ch22.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/ch22.tex), label `prob:wigfree:ladder`.

## Chapter 13. Ferromagnetic random variables

<a id="ferromagnetic_law_clock_and_van_dantzig_partner"></a>

### Clock and van Dantzig partner of a ferromagnetic law

Theorem 13.1.4, page 618 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable $X$ is called ferromagnetic if it is symmetric, $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function $L_X(z)=\mathbb{E}e^{zX}$, which is then entire, vanishes only on the imaginary axis; equivalently every zero of the characteristic function $\varphi_X(\theta)=L_X(i\theta)$ is real. By Newman's theorem such an $X$ has the Hadamard form $L_X(s)=e^{bs^2}\prod_k(1+s^2/c_k^2)$ with $b\ge0$, $c_k>0$ and $\sum_kc_k^{-2}<\infty$, each $c_k$ repeated according to the multiplicity of the zero $ic_k$. Let $X$ be ferromagnetic with this form, let $e_k$ be independent standard exponential variables and $Z$ a standard normal variable independent of them, and put $T_X=b+\sum_ke_k/c_k^2$ and $Y_X=Z\sqrt{2T_X}$; the series converges almost surely and in mean. Then (i) $\mathbb{E}e^{-uT_X}=1/L_X(\sqrt{u})$ for $u\ge0$, and $T_X$ is a generalised gamma convolution (GGC, a law on $[0,\infty)$ with Laplace transform $\exp\{-b\lambda-\int\log(1+\lambda/t)\,\mu(dt)\}$ for a drift $b\ge0$ and a positive Thorin measure $\mu$) with drift $b$ and Thorin measure $\sum_k\delta_{c_k^2}$, one unit atom for each pair of zeros $\pm ic_k$. (ii) $\mathbb{E}e^{i\theta Y_X}=1/L_X(\theta)$ for real $\theta$, so that $\varphi_{Y_X}(\theta)\,\varphi_X(i\theta)=1$ and $(\varphi_X,\varphi_{Y_X})$ is a van Dantzig pair, and $Y_X\overset{d}{=}\sqrt{2b}\,Z+\sum_kL_k/c_k$ with independent standard Laplace variables $L_k$, whose characteristic function is $1/(1+\theta^2)$. (iii) The moment generating function of $Y_X$ is $1/\varphi_X(t)$ for $|t|<c_1=\min_kc_k$ and infinite for $|t|\ge c_1$, unless $X$ is Gaussian. (iv) Unless $X$ is degenerate, $Y_X$ has a Pólya frequency density of infinite order. (v) $\kappa_m(T_X)=(-1)^{m+1}m!\,\kappa_{2m}(X)/(2m)!$ for $m\ge1$, where $\kappa_j$ denotes the $j$th cumulant; in particular $\mathbb{E}T_X=\frac12\operatorname{Var}X$ and $\operatorname{Var}T_X=-\kappa_4(X)/12$.

Source: [`ch/fm1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm1.tex), label `fm:thm:dual`.

<a id="ferromagnetic_laws_closed_under_convergence_in_law"></a>

### Ferromagnetic laws are closed under convergence in law

Theorem 13.1.5, page 620 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable $X$ is called ferromagnetic if it is symmetric, $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function $L_X(z)=\mathbb{E}e^{zX}$, which is then entire, vanishes only on the imaginary axis; Gaussian laws are allowed. Let $X_n$ be ferromagnetic and suppose that $X_n\to X$ in law. Then $X$ is ferromagnetic, $\operatorname{Var}X_n\to\operatorname{Var}X$, and $L_{X_n}\to L_X$ locally uniformly on $\mathbb{C}$.

Source: [`ch/fm1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm1.tex), label `fm:thm:closed`.

<a id="rh_iff_polya_law_is_a_limit_in_law_of_ferromagnetic_laws"></a>

### RH holds exactly when the Pólya law is a limit in law of ferromagnetic laws

Corollary 13.1.7, page 621 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable $X$ is called ferromagnetic if it is symmetric, $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function $\mathbb{E}e^{zX}$, which is then entire, vanishes only on the imaginary axis. The Pólya law is the law of the symmetric random variable $X$ with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$; its moment generating function is $\mathbb{E}e^{sX}=\xi(\frac12+s)/\xi(\frac12)$. The Riemann Hypothesis holds if and only if the Pólya law is a limit in law of ferromagnetic laws.

Source: [`ch/fm1.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm1.tex), label `fm:cor:limit`.

<a id="mixture_of_uniform_law_with_atom_at_origin_is_not_ferromagnetic"></a>

### A uniform law mixed with an atom at the origin is not ferromagnetic

Theorem 13.2.1, page 622 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: counterexample, proved.

Statement. A real random variable is called ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function, which is then entire, vanishes only on the imaginary axis. Let $U$ be uniform on $[-1,1]$, so that $\mathbb{E}e^{sU}=\sinh s/s$. For every $0<q<1$ the mixture $(1-q)\,\mathrm{Unif}[-1,1]+q\,\delta_0$ of the uniform law on $[-1,1]$ with weight $1-q$ and the point mass at $0$ with weight $q$ is symmetric and unimodal but not ferromagnetic: its moment generating function $(1-q)\sinh s/s+q$ has no real zero, only finitely many zeros on the imaginary axis, and infinitely many zeros in all, so that all but finitely many of its zeros are off the imaginary axis. Consequently, for every $0<q<1$ and every small enough $a>0$, the mixture of the uniform laws on $[-1,1]$ and on $[-a,a]$ with weights $1-q$ and $q$ is not ferromagnetic. Both components of these mixtures are ferromagnetic, so the class of ferromagnetic laws is not closed under mixtures.

Source: [`ch/fm2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm2.tex), label `fm:thm:atom`.

<a id="shoulder_theorem_for_ferromagnetic_laws_on_an_interval"></a>

### Shoulder theorem for symmetric laws on an interval

Theorem 13.2.3, page 623 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable is called ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function, which is then entire, vanishes only on the imaginary axis; equivalently every zero of its characteristic function $\varphi_X$ is real. Let $p$ be a symmetric probability density on $[-1,1]$ that is absolutely continuous and nonincreasing on $[0,1]$, with $p(1)=0$, let $[-c,c]$ with $0<c\le1$ be its support, and let $X$ have density $p$. (i) If $p$ is concave on $[0,1]$, then $X$ is ferromagnetic. (ii) If $p$ is convex on $[0,1]$ and not linear on $[0,c]$, then $X$ is not ferromagnetic; moreover $\varphi_X\ge0$ on $\mathbb{R}$, and its real zeros, if any, are double and lie in a lattice $(2\pi/s_0)\mathbb{Z}$ with $0<s_0<c$.

Source: [`ch/fm2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm2.tex), label `fm:thm:shoulder`.

<a id="ferromagnetic_laws_among_power_kernels_and_products_of_uniform_laws"></a>

### Ferromagnetic laws among the kernels 1 minus |x|^b, k monotone kernels and products of uniform laws

Corollary 13.2.4, page 624 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable is called ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function, which is then entire, vanishes only on the imaginary axis. (i) The law with density proportional to $1-|x|^b$ on $[-1,1]$, $b>0$, is ferromagnetic if and only if $b\ge1$. In Khinchine's form, $UR$ with $U$ uniform on $[-1,1]$ and $R$ independent with the Beta$(a,1)$ law (density $ar^{a-1}$ on $[0,1]$) is ferromagnetic if and only if $a\ge2$. (ii) For an integer $k\ge1$ the symmetric $k$ monotone kernel $f_k(x)=\frac k2(1-|x|)_+^{k-1}$ (the uniform density for $k=1$ and the triangle density for $k=2$) is ferromagnetic if and only if $k\le2$. (iii) The laws with density proportional to $(1-x^2)^m$ on $[-1,1]$ with $0<m\le1$, which include the semicircle law, are ferromagnetic. (iv) The product $U_1U_2$ of two independent uniform variables on $[-1,1]$, with density $\frac12\log(1/|x|)$ on $[-1,1]$, is not ferromagnetic.

Source: [`ch/fm2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm2.tex), label `fm:cor:shoulder`.

<a id="shell_mixtures_of_uniform_laws_are_ferromagnetic"></a>

### Densities nondecreasing in |x| and mixtures of uniform laws on shells are ferromagnetic

Theorem 13.2.5, page 625 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. A real random variable is called ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function, which is then entire, vanishes only on the imaginary axis. Every symmetric probability density on $[-1,1]$ that is nondecreasing in $|x|$ is ferromagnetic. Equivalently, every mixture of the uniform laws on the shells $\{a\le|x|\le1\}$, $0\le a<1$, is ferromagnetic.

Source: [`ch/fm2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm2.tex), label `fm:thm:ushaped`.

<a id="polya_law_has_no_gaussian_component_and_no_finite_zero_summand"></a>

### The Pólya law has no Gaussian component and no ferromagnetic summand with finitely many zeros

Theorem 13.3.1, page 631 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be a symmetric random variable with $\mathbb{E}e^{sX}<\infty$ for every real $s$ and $\log\mathbb{E}e^{sX}=o(s^2)$ as $s\to\infty$, and suppose that $X\overset{d}{=}D+R$ with $D$ and $R$ independent. Then $\log\mathbb{E}e^{sD}=o(s^2)$ as $|s|\to\infty$. In particular, if $D$ is ferromagnetic (symmetric, with $\mathbb{E}e^{b_0D^2}<\infty$ for some $b_0>0$ and a moment generating function vanishing only on the imaginary axis) and not $0$, then in its Hadamard form $\mathbb{E}e^{sD}=e^{bs^2}\prod_k(1+s^2/c_k^2)$ the Gaussian factor has $b=0$, and the characteristic function of $D$ has infinitely many zeros. The Pólya variable, the symmetric random variable with density $\Phi/\xi(\frac12)$ and moment generating function $\xi(\frac12+s)/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, satisfies the hypothesis with no assumption on the zeros of $\zeta$. So it has no Gaussian component (it is not the sum of an independent normal variable of positive variance and another variable), and it is not a finite or convergent infinite sum of independent ferromagnetic variables of which one is not $0$ and has a characteristic function with finitely many zeros.

Source: [`ch/fm3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm3.tex), label `fm:thm:nosqrt`.

<a id="cosine_and_gaussian_tilts_of_polya_law_ferromagnetic_if_no_zero_beyond_seven_eighths"></a>

### Cosine and Gaussian tilts of the Pólya law with variance at least 9/64 if there is no zero beyond 7/8

Proposition 13.4.2, page 633 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on The Riemann zeta function has no zero with real part greater than 7/8 (announced by OpenAI in October 2026, not refereed, and not proved in the book).

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$. For $\lambda_1,\dots,\lambda_N>0$ and $\mu\ge0$ put $m(t)=e^{\mu t^2}\prod_{k\le N}\cosh(\lambda_kt)$ and $c(m)=2\mu+\sum_{k\le N}\lambda_k^2$, so that $m(t)=\mathbb{E}e^{tY}$ for $Y=\sqrt{2\mu}\,Z+\sum_k\lambda_k\varepsilon_k$ with $Z$ standard normal and $\varepsilon_k$ independent fair random signs, and $c(m)=\operatorname{Var}Y$; the same applies with $N=\infty$ when $\sum_k\lambda_k^2<\infty$. A law is ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function vanishes only on the imaginary axis. Assume that $\zeta$ has no zero with real part greater than $\frac78$. Then the zeros of $z\mapsto\int\Phi(t)m(t)e^{izt}\,dt$ lie in the strip $|\operatorname{Im}z|\le(\max(\frac9{64}-c(m),0))^{1/2}$. In particular the law with density proportional to $m\,p$ is ferromagnetic whenever $\operatorname{Var}Y\ge\frac9{64}$; for example $\cosh(at)\,p(t)$ is ferromagnetic for $a\ge\frac38$ and $\frac{\sinh(ct)}{ct}\,p(t)$ for $c\ge\frac{3\sqrt3}8$. With no hypothesis on the zeros the same statements hold with $\frac14$ in place of $\frac9{64}$, by theorems of de Bruijn (1950), so that $\cosh(at)\,p(t)$ is ferromagnetic for $a\ge\frac12$.

Source: [`ch/fm4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm4.tex), label `fm:prop:bruijn`.

<a id="zero_completed_laws_maxwell_construction_and_rh_criterion"></a>

### Zero completed laws of the Pólya variable, the Maxwell construction and a criterion for RH

Theorem 13.4.3, page 635 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, with characteristic function $\varphi_X$. Let $\gamma_1\le\gamma_2\le\cdots$ be the positive ordinates of the zeros $\frac12\pm i\gamma_k$ of $\zeta$ on the critical line, with multiplicity, and put $s_K=\sum_{k\le K}\gamma_k^{-2}$, $\beta_K=\frac12\operatorname{Var}X - s_K$ and $\varphi_K(\theta)=e^{-\beta_K\theta^2}\prod_{k\le K}(1-\theta^2/\gamma_k^2)$. Let $M_1,M_2,\dots$ be independent symmetric Maxwell variables (density $x^2e^{-x^2/2}/\sqrt{2\pi}$, characteristic function $(1-\theta^2)e^{-\theta^2/2}$, variance $3$), each ferromagnetic, and let $Z$ be a standard normal variable independent of them. A law is ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function vanishes only on the imaginary axis. (i) With no hypothesis, $\beta_K>0$ for every $K$. If $3s_K\le\operatorname{Var}X$, then $\varphi_K$ is the characteristic function of $X_K=(\operatorname{Var}X-3s_K)^{1/2}Z+\sum_{k\le K}M_k/\gamma_k$, whose law $f_K$ is ferromagnetic, with $\operatorname{Var}X_K=\operatorname{Var}X$ and $\kappa_4(X_K)=-12\sum_{k\le K}\gamma_k^{-4}$. The condition $3s_K\le\operatorname{Var}X$ holds for $K\le16$ and fails for $K=17$, a certified computation in ball arithmetic ($3s_{16}=0.0456803$, $\operatorname{Var}X=0.0462100$, $3s_{17}=0.0463006$). (ii) For every $K$, if the inverse Fourier transform $f_K$ of $\varphi_K$ is nonnegative, then $f_K$ is a ferromagnetic probability density. (iii) The Riemann Hypothesis holds if and only if $\varphi_K\to\varphi_X$ pointwise as $K\to\infty$. Under the Riemann Hypothesis, $f_K$ converges in law to the Pólya law along every sequence of $K$ for which $f_K\ge0$.

Source: [`ch/fm4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm4.tex), label `fm:thm:zc`.

<a id="rh_iff_polya_variable_smoothed_by_gaussian_equals_maxwell_series_in_law"></a>

### RH is equivalent to an identity in law between the smoothed Pólya variable and a series of Maxwell variables

Theorem 13.4.5, page 638 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and let $Z$ be a standard normal variable independent of $X$. Let $\gamma_1\le\gamma_2\le\cdots$ be the positive ordinates of the zeros $\frac12\pm i\gamma_k$ of $\zeta$ on the critical line, with multiplicity, let $M_1,M_2,\dots$ be independent symmetric Maxwell variables with density $x^2e^{-x^2/2}/\sqrt{2\pi}$, and put $\mathcal M=\sum_kM_k/\gamma_k$, a series that converges almost surely because its terms are independent and centred with summable variances $3\gamma_k^{-2}$. Then the Riemann Hypothesis holds if and only if $X+\big(\tfrac12\operatorname{Var}X\big)^{1/2}Z\overset{d}{=}\mathcal M$.

Source: [`ch/fm4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm4.tex), label `fm:thm:maxwell`.

<a id="polya_law_within_total_variation_two_millionths_of_a_ferromagnetic_law"></a>

### Total variation distance from the Pólya law to ferromagnetic laws is below two millionths, with no hypothesis

Proposition 13.4.6, page 639 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$. A law is ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0Y^2}<\infty$ for some $b_0>0$, and its moment generating function vanishes only on the imaginary axis. Let $\delta$ be the infimum of the total variation distances $\sup_A|P(X\in A)-P(Y\in A)|$ from $X$ to ferromagnetic variables $Y$; when $Y$ has a density $f$ the distance is $\operatorname{TV}(f,p)=\frac12\int|f - p|$. Let $\gamma_1\le\gamma_2\le\cdots$ be the positive ordinates of the zeros $\frac12\pm i\gamma_k$ of $\zeta$ on the critical line, with multiplicity, put $s_K=\sum_{k\le K}\gamma_k^{-2}$ and $\beta_K=\frac12\operatorname{Var}X - s_K$, and let $f_K$ be the function with Fourier transform $\varphi_K(\theta)=e^{-\beta_K\theta^2}\prod_{k\le K}(1-\theta^2/\gamma_k^2)$. With no hypothesis on the zeros: $f_{16}$ is a ferromagnetic probability density, an explicit sum of a normal variable and Maxwell variables, and $\operatorname{TV}(f_{16},p)<1.266\times10^{-4}$, so that $\delta<1.266\times10^{-4}$; and $f_{150}$ is a ferromagnetic probability density by a computer assisted proof in ball arithmetic (the polynomial in the Hermite basis that multiplies the Gaussian is positive on the real line) with $\operatorname{TV}(f_{150},p)<1.932\times10^{-6}$, so that $\delta<1.932\times10^{-6}$. The distances are enclosed in ball arithmetic, with explicit bounds for the discretisation and the tails. The Riemann Hypothesis holds if and only if $\delta=0$.

Source: [`ch/fm4.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm4.tex), label `fm:prop:tvcert`.

<a id="powers_of_exponential_of_polya_variable_not_infinitely_divisible_below_two"></a>

### Powers of the exponential of the Pólya variable are not infinitely divisible for exponent below 2

Proposition 13.5.1, page 640 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$. For every real $\beta$ with $0<|\beta|<2$ the random variable $e^{\beta X}$ is not infinitely divisible. No hypothesis on the zeros of $\zeta$ is used.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:prop:notid`.

<a id="exponential_of_polya_variable_is_not_a_generalised_gamma_convolution"></a>

### The exponential of the Pólya variable is not a generalised gamma convolution

Theorem 13.5.3, page 642 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, and let $B_0=\sum_{m\ge1}\Gamma_m/(\pi m^2)$ with independent Gamma$(2,1)$ variables $\Gamma_m$; this is the variable $\frac\pi2S_2$ of Biane, Pitman and Yor, a generalised gamma convolution (GGC) with Laplace transform $(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^2$ and $\mathbb{E}B_0^\sigma=2\xi(2\sigma)$ for every complex $\sigma$. For $0<q<1$ let $B_q$ be $B_0$ size biased by $B_0^q$, that is $\mathbb{E}g(B_q)=\mathbb{E}[B_0^qg(B_0)]/\mathbb{E}B_0^q$, so that $\mathbb{E}B_q^\sigma=\xi(2\sigma+2q)/\xi(2q)$. A GGC is a law on $[0,\infty)$ with Laplace transform $\exp\{-b\lambda-\int\log(1+\lambda/t)\,\mu(dt)\}$ for a drift $b\ge0$ and a positive Thorin measure $\mu$. Then for every $0<q<1$ the variable $B_q$ is not a GGC. In particular $e^{2X}$ is not a GGC, where $X$ is the Pólya variable, the symmetric random variable with density $\Phi/\xi(\frac12)$, $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and moment generating function $\xi(\frac12+s)/\xi(\frac12)$, because $e^{2X}\overset{d}{=}B_{1/4}$. No hypothesis on the zeros is used.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:thm:notggc`.

<a id="size_bias_principle_for_sums_of_gamma_variables_with_integer_shapes"></a>

### Size bias principle for sums of gamma variables with integer shapes

Theorem 13.5.5, page 643 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $Y=\sum_{n=1}^K\Gamma_n/\tau_n$ with $0<\tau_1<\dots<\tau_K$ and independent gamma variables $\Gamma_n$ with integer shapes $t_n\ge1$ and unit scale, let $0<q<1$, and let $Y_q$ be $Y$ size biased by $Y^q$, that is $\mathbb{E}g(Y_q)=\mathbb{E}[Y^qg(Y)]/\mathbb{E}Y^q$. Then $Y_q+G\overset{d}{=}\sum_{n=1}^K\Gamma_n'/\tau_n$, where the $\Gamma_n'$ are independent gamma variables with shapes $t_n+q$ and unit scale, and $G$, independent of $Y_q$, is a generalised gamma convolution without drift whose Thorin measure has a continuous distribution function $\theta_G$ that vanishes on $(0,\tau_1]$, increases by exactly $q$ on each interval $[\tau_k,\tau_{k+1}]$ and is constant on $[\tau_K,\infty)$.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:thm:sbp`.

<a id="factorisation_of_exponential_of_polya_variable_against_biane_pitman_yor_variable"></a>

### Factorisation of the exponential of the Pólya variable against a Biane Pitman Yor variable

Theorem 13.5.6, page 645 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$, let $B_0=\sum_{m\ge1}\Gamma_m/(\pi m^2)$ with independent Gamma$(2,1)$ variables $\Gamma_m$ (the variable $\frac\pi2S_2$ of Biane, Pitman and Yor, with $\mathbb{E}B_0^\sigma=2\xi(2\sigma)$), and for $0<q<1$ let $B_q$ be $B_0$ size biased by $B_0^q$, so that $\mathbb{E}B_q^\sigma=\xi(2\sigma+2q)/\xi(2q)$. Let $\frac\pi2S_{2+q}=\sum_{m\ge1}\Gamma_m'/(\pi m^2)$ with independent gamma variables $\Gamma_m'$ of shape $2+q$ and unit scale, with Laplace transform $(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^{2+q}$. Then $B_q+G_q\overset{d}{=}\frac\pi2S_{2+q}$, where $G_q$ is independent of $B_q$ and is a generalised gamma convolution (a law on $[0,\infty)$ with Laplace transform $\exp\{-b\lambda-\int\log(1+\lambda/t)\,\mu(dt)\}$) without drift whose Thorin measure has a continuous distribution function $\theta_q$ with $\theta_q=0$ on $(0,\pi]$ and $\theta_q(\pi m^2)=(m-1)q$ for every $m\ge1$. In particular $e^{2X}+G_{1/4}\overset{d}{=}\frac\pi2S_{9/4}$, where $X$ is the Pólya variable, the symmetric random variable with density $\Phi/\xi(\frac12)$, $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and moment generating function $\xi(\frac12+s)/\xi(\frac12)$, because $e^{2X}\overset{d}{=}B_{1/4}$. No hypothesis on the zeros is used.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:thm:factor`.

<a id="exponential_of_polya_variable_self_decomposable_with_theta_levy_density"></a>

### The exponential of the Pólya variable is self decomposable with a Lévy density built from the theta function

Theorem 13.5.7, page 646 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\omega(x)=\sum_{n\ge1}e^{-\pi n^2x}$, let $B_q$, $0<q<1$, be the variable $B_0=\sum_{m\ge1}\Gamma_m/(\pi m^2)$ (independent Gamma$(2,1)$ variables $\Gamma_m$, the variable $\frac\pi2S_2$ of Biane, Pitman and Yor) size biased by $B_0^q$, and let $\theta_q$ be the distribution function of the Thorin measure of the generalised gamma convolution $G_q$ in the identity $B_q+G_q\overset{d}{=}\frac\pi2S_{2+q}$, where $\frac\pi2S_{2+q}$ has Laplace transform $(\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda})^{2+q}$; $\theta_q$ is continuous, vanishes on $(0,\pi]$, and satisfies $\theta_q(\pi m^2)=(m-1)q$. Then $B_q$ is self decomposable, with no drift and Lévy density $k_q(x)/x$, where $k_q(x)=(2+q)\,\omega(x)-\int_0^\infty e^{-tx}\,d\theta_q(t)$ and $2\,\omega(x)\le k_q(x)\le(2+q)\,\omega(x)$. In particular $e^{2X}$, with $X$ the Pólya variable (the symmetric random variable with density $\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and with moment generating function $\xi(\frac12+s)/\xi(\frac12)$), which has the law of $B_{1/4}$, is self decomposable, hence infinitely divisible and unimodal; its Lévy density lies between $2\omega(x)/x$ and $\frac94\omega(x)/x$, and it is not a generalised gamma convolution. No hypothesis on the zeros is used.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:thm:sd`.

<a id="hyperbolically_completely_monotone_powers_of_exponential_of_polya_variable"></a>

### Hyperbolically completely monotone densities of the powers of the exponential of the Pólya variable

Proposition 13.5.8, page 648 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$. A density $f$ on $(0,\infty)$ is hyperbolically completely monotone (HCM) if for every $u>0$ the function $f(uv)f(u/v)$ of $w=v+1/v$ is completely monotone; an HCM law is a generalised gamma convolution (Bondesson). The variable $e^{\beta X}$ has an HCM density, and so is a generalised gamma convolution, for $|\beta|\ge5.5791422$, and it does not have an HCM density for $|\beta|\le5.5791395$. The proof of the first statement is computer assisted, in ball arithmetic. No hypothesis on the zeros is used.

Source: [`ch/fm5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm5.tex), label `fm:prop:hcm`.

<a id="open_mixed_shoulders_sine_transform_of_slope"></a>

### Mixed shoulders and the sine transform of the slope (open problem)

Problem 13.1, page 650 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. A real random variable is called ferromagnetic if it is symmetric, has $\mathbb{E}e^{b_0X^2}<\infty$ for some $b_0>0$, and its moment generating function vanishes only on the imaginary axis. Let $p$ be a symmetric probability density on $[-1,1]$, absolutely continuous and nonincreasing on $[0,1]$ with $p(1)=0$, and put $g=-p'\ge0$ on $(0,1)$. For a variable $X$ with density $p$ and characteristic function $\varphi_X$, integration by parts gives $z\,\varphi_X(z)=2\int_0^1g(y)\sin(zy)\,dy$, so $X$ is ferromagnetic if and only if the sine transform of the slope $g$ has only real zeros. The book decides the case in which $g$ is monotone: $X$ is ferromagnetic if $p$ is concave on $[0,1]$, and not ferromagnetic if $p$ is convex on $[0,1]$ and not linear on its support. Open problem: decide the slopes that first rise and then fall. One class is settled. If $g$ is symmetric about $\frac12$, then $p(0)=\int_0^1g=1$, so that $u\mapsto g(\frac12+u)$ is a probability density on $[-\frac12,\frac12]$, and $X\overset{d}{=}\frac12U+B$ with $U$ uniform on $[-1,1]$ and $B$ independent with that density; since $\sin(z/2)/(z/2)$ has only real zeros, $X$ is then ferromagnetic exactly when $B$ is. For the beta slopes $g(y)\propto y^a(1-y)^b$ with $a,b\ge0$ this decides the diagonal $a=b$, where $B$ is a symmetric beta law and $X$ is ferromagnetic, and the shoulder theorem decides the edges: $X$ is ferromagnetic for $b=0$, where $p$ is concave, and not for $a=0<b$, where $p$ is convex and not linear. The rest of the quadrant is not decided. On the half line the slope of the Pólya law, the symmetric law with density $\Phi/\xi(\frac12)$ where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, rises and then falls, and the Riemann Hypothesis is the statement that its sine transform has only real zeros.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:shoulder`.

<a id="open_polya_law_in_griffiths_simon_class"></a>

### Does the Pólya law belong to the Griffiths Simon class (open problem)

Problem 13.2, page 650 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. The Griffiths Simon class consists of the weak limits of suitably scaled sums of ferromagnetically coupled Ising spins; its laws are ferromagnetic, that is symmetric with a moment generating function that vanishes only on the imaginary axis. The class of Ellis, Monroe and Newman consists of the laws with density $e^{-V(x)}$, $V$ even and $V'$ convex on $[0,\infty)$, and their weak limits; by a theorem of Newman the Pólya law (the symmetric law with density $\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and with moment generating function $\xi(\frac12+s)/\xi(\frac12)$) belongs to it. Open problem: does the Pólya law belong to the Griffiths Simon class? A positive answer would imply the Riemann Hypothesis. Under Ingham's hypothesis, that the positive ordinates of the zeros of $\zeta$ are linearly independent over the rationals, the Pólya law is not a sum of independent scaled fair signs, so the spins of any approximating system would have to be coupled.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:gs`.

<a id="open_splitting_the_zeros_into_two_characteristic_functions"></a>

### Splitting the zeros of the Pólya characteristic function into two characteristic functions (open problem)

Problem 13.3, page 650 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Assume the Riemann Hypothesis, and let $X$ be the Pólya variable (the symmetric random variable with density $\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and with moment generating function $\xi(\frac12+s)/\xi(\frac12)$), so that its characteristic function is $\varphi_X(t)=\prod_\gamma(1-t^2/\gamma^2)$ over the positive ordinates $\gamma$ of the zeros of $\zeta$, with no Gaussian factor. Open problem: can the ordinates be split into two nonempty sets $A$ and $A'$ such that $\prod_{\gamma\in A}(1-t^2/\gamma^2)$ and $\prod_{\gamma\in A'}(1-t^2/\gamma^2)$ are both characteristic functions? This is the same as writing $X=D+R$ with $D$ and $R$ independent, ferromagnetic and not $0$, and both sets are then infinite. Under Ingham's hypothesis, that the positive ordinates are linearly independent over the rationals, no three ordinates are in arithmetic progression, so no piece is a scaled fair sign or a scaled uniform variable, whose zeros form arithmetic progressions.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:split`.

<a id="open_positivity_of_zero_completed_laws_and_the_threshold_v_star"></a>

### Positivity of the zero completed laws and the threshold near K = 5000 (open problem)

Problem 13.4, page 650 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, with characteristic function $\varphi_X$, let $\gamma_1\le\gamma_2\le\cdots$ be the positive ordinates of the zeros $\frac12\pm i\gamma_k$ of $\zeta$ on the critical line, with multiplicity, put $s_K=\sum_{k\le K}\gamma_k^{-2}$ and $\beta_K=\frac12\operatorname{Var}X - s_K>0$, and let $f_K$ be the function with Fourier transform $\varphi_K(\theta)=e^{-\beta_K\theta^2}\prod_{k\le K}(1-\theta^2/\gamma_k^2)$. For $v>0$ the function $e^{-v\theta^2/2}\prod_{k\le K}(1-\theta^2/\gamma_k^2)$ is the Fourier transform of $f_{K,v}(x)=g_v(x)\sum_{j=0}^Ke_jv^{-j}\mathrm{He}_{2j}(x/\sqrt{v})$, where $g_v$ is the $N(0,v)$ density, $\mathrm{He}_n(y)=\mathbb{E}(y+iZ)^n$ with $Z$ standard normal is the Hermite polynomial and $e_j$ is the $j$th elementary symmetric function of $\gamma_1^{-2},\dots,\gamma_K^{-2}$; there is a threshold $v^*(K)$ with $f_{K,v}\ge0$ if and only if $v\ge v^*(K)$, and $f_K=f_{K,2\beta_K}$. The book proves that $f_K$ is a ferromagnetic probability density for $K\le16$ by the Maxwell construction, and for every $K\le200$ and for $K=300$, $400$, $600$ and $800$ by a computer assisted proof in ball arithmetic; the ratio $v^*(K)/2\beta_K$ increases slowly, to $0.733$ at $K=800$. Open problem: compute the threshold $v^*(K)$ near $K=5000$, where a straight line in $\log K$ through the last computed points, which reaches $1$ near $K=1.8\times10^4$, and the power law $1-v^*(K)/2\beta_K\approx1.65K^{-0.27}$, which never does, part ways; the cancellation there calls for a working precision of tens of thousands of bits. Is $v^*(K)<2\beta_K$ for every $K$? If so, every $f_K$ is a ferromagnetic probability density, and with no hypothesis their limit in law, with characteristic function $e^{-\beta_\infty\theta^2}\prod_k(1-\theta^2/\gamma_k^2)$ over the zeros on the critical line, where $\beta_\infty=\lim_K\beta_K$ is $0$ under the Riemann Hypothesis, is a ferromagnetic law; it is the Pólya law exactly when the Riemann Hypothesis holds.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:zc`.

<a id="open_cosine_constant_omega_interval_and_improved_bounds"></a>

### The cosine constant Omega, whether the good tilts form an interval, and improved bounds (open problem)

Problem 13.5, page 651 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $p$ be the Pólya density $\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$, and let $\Omega$ be the infimum of the $a_0\ge0$ such that the law with density proportional to $\cosh(at)\,p(t)$ is ferromagnetic (symmetric with a moment generating function vanishing only on the imaginary axis) for every $a\ge a_0$. By de Bruijn's theorem, with no hypothesis, $\Omega\le\frac12$; $\Omega\le\frac38$ if the half plane $\operatorname{Re}s>\frac78$ is free of zeros of $\zeta$, as announced by OpenAI in October 2026 and not refereed; and the Riemann Hypothesis holds if and only if $\Omega=0$. Open problem: is the set of $a\ge0$ for which $\cosh(at)\,p$ is ferromagnetic an interval? Improve these bounds by computation, as computations with the zeros have improved de Bruijn's bound $\frac12$ for the de Bruijn Newman constant $\Lambda$ to $\Lambda\le0.2$ (Platt and Trudgian). A first case is $a=\sqrt{0.1}$, below both bounds: the law $\cosh(\sqrt{0.1}\,t)\,p$ is at total variation distance $1.094\times10^{-3}$ from $p$, as close as the de Bruijn Newman flow at time $\lambda=0.2$, the law with density proportional to $e^{\lambda t^2/4}p(t)$ (distance $1.097\times10^{-3}$), which is ferromagnetic, and deciding it is a computation with the zeros of $\frac12(\xi(s+a)+\xi(s - a))$, of the kind that gave the bound $0.2$ for the flow.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:omega`.

<a id="open_size_bias_principle_for_real_shapes"></a>

### The size bias principle for real shapes (open problem)

Problem 13.6, page 651 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $Y=\sum_{n=1}^K\Gamma_n/\tau_n$ with $0<\tau_1<\dots<\tau_K$ and independent gamma variables $\Gamma_n$ of shapes $t_n>0$ and unit scale, let $0<q<1$, and let $Y_q$ be $Y$ size biased by $Y^q$, that is $\mathbb{E}g(Y_q)=\mathbb{E}[Y^qg(Y)]/\mathbb{E}Y^q$. The size bias principle proved in the book states that for integer shapes $t_n\ge1$, $Y_q+G\overset{d}{=}\sum_n\Gamma_n'/\tau_n$, where the $\Gamma_n'$ are independent gamma variables of shapes $t_n+q$ and unit scale and $G$, independent of $Y_q$, is a generalised gamma convolution without drift whose Thorin measure has a continuous distribution function that vanishes on $(0,\tau_1]$, increases by exactly $q$ on each interval $[\tau_k,\tau_{k+1}]$ and is constant on $[\tau_K,\infty)$. Open problem: does it hold for all real shapes $t_n>0$? For $K=1$ it is the gamma identity, with $G=0$. Integer shapes enter the proof twice, through the divided difference formula for the density of $Y$ and through a sign lemma, and together they show that on each gap $(\tau_k,\tau_{k+1})$ the phase, a continuous determination of $-\frac1\pi$ times the argument of the boundary value at $-r+i0$ of the Laplace transform of $Y_q$, is nonincreasing and falls by exactly $q$. Is this still true for real shapes? With it, the rest of the proof goes through, provided the estimates of the transform near the points $-\tau_n$ and at infinity can be checked as well.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:sbp`.

<a id="open_powers_of_exponential_of_polya_variable_below_hcm_threshold"></a>

### Infinite divisibility of powers of the exponential of the Pólya variable below the HCM threshold (open problem)

Problem 13.7, page 651 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $X$ be the Pólya variable, the symmetric random variable with density $p=\Phi/\xi(\frac12)$, where $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ and $\Phi(t)=\sum_{n\ge1}(4\pi^2n^4e^{9t/2}-6\pi n^2e^{5t/2})\exp(-\pi n^2e^{2t})$. The variable $e^{\beta X}$ is not infinitely divisible for $0<|\beta|<2$, is self decomposable and not a generalised gamma convolution (GGC) at $|\beta|=2$, and has a hyperbolically completely monotone (HCM) density, so is a GGC, for $|\beta|\ge5.5791422$; it has no HCM density for $|\beta|\le5.5791395$. The HCM threshold, the least $|\beta|$ for which $e^{\beta X}$ has an HCM density, lies between $5.5791395$ and $\pi/0.563096=5.57914219\ldots$. Open problem: for $2<|\beta|$ below the HCM threshold, is $e^{\beta X}$ infinitely divisible, self decomposable or a GGC? Below the threshold it has no HCM density, which leaves the GGC property open, and self decomposability is established only at $|\beta|=2$. Is the HCM threshold exactly $\pi/v_m=5.5791396$, where $v_m=0.5630962646$ is the point at which $F(v)=\Phi(iv)$ has its local maximum on the imaginary axis? That requires $h\ge0$ on the open strip $0<v<v_m$, where $h(x,v)=\operatorname{Im}\big(-\Phi'/\Phi\big)(x+iv)$ and $h(0,v)\to0$ as $v\uparrow v_m$, so that a certificate must replace the balls with positive lower bounds near $iv_m$ by a local argument.

Source: [`ch/fm6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/fm6.tex), label `fm:prob:powers`.

## Chapter 14. The Pólya density as a maximum order statistic and the Thorin measure of the clock

<a id="polya_kernel_from_sum_of_two_maxima_of_exponentials"></a>

### Pólya kernel from a sum of two maxima of exponential variables

Proposition 14.2.1, page 657 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$ be Riemann's xi function and $Y$ the random variable with $\mathbb{E}e^{sY}=2\xi(s)$ for every complex $s$. Let $M_4$ be the maximum of independent exponential random variables with one variable of each rate $\pi n$, $n\ge1$, and one more of each rate $\pi(2m-1)$, $m\ge1$, so that $\Pr(M_4\le t)=\prod_{n\ge1}(1-e^{-\pi nt})\prod_{m\ge1}(1-e^{-\pi(2m-1)t})=\sum_{n\in\mathbb{Z}}(-1)^ne^{-\pi n^2t}$ for $t>0$, and let $M_4'$ be an independent copy of $M_4$. Then (i) $e^{2Y}$ has the same law as $M_4+M_4'$. (ii) The density of $e^Y=\sqrt{M_4+M_4'}$ is $f(y)=4\sum_{n\ge1}\big(2\pi^2n^4y^3-3\pi n^2y\big)e^{-\pi n^2y^2}$ for $y>0$, and the Pólya kernel $\Phi(u)=2\sum_{n\ge1}\big(2\pi^2n^4e^{9u/2}-3\pi n^2e^{5u/2}\big)\exp(-\pi n^2e^{2u})$, the even function with $\int_{\mathbb{R}}\Phi(u)e^{i\theta u}\,du=\xi(\frac12+i\theta)$, satisfies $\Phi(u)=\frac12e^{3u/2}f(e^u)$.

Source: [`ch/maxwald.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald.tex), label `mw:prop:sum`.

<a id="khinchin_factorisation_of_tilted_maximum_of_exponentials"></a>

### Khinchin factorisation of the tilted maximum of exponential variables

Theorem 14.3.1, page 659 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Fix $\sigma>1$. Let $M_4$ be the maximum of independent exponential random variables with one variable of each rate $\pi n$, $n\ge1$, and one more of each rate $\pi(2m-1)$, $m\ge1$, and let $\widetilde M_4$ have the tilted law $\Pr(\widetilde M_4\in dm)\propto m^{\sigma/2}\Pr(M_4\in dm)$. Let $K$, $N_\sigma$ and $\Gamma_{1+\sigma/2}$ be independent of each other and of $\widetilde M_4$, where $\Pr(K=k)=(1-a)a^k$ for $k\ge0$ with $a=2^{1-\sigma}$; $N_\sigma=\prod_pp^{G_p}$ is zeta distributed with parameter $\sigma$, that is $\Pr(N_\sigma=n)=n^{-\sigma}/\zeta(\sigma)$ for $n\ge1$, the product running over the primes $p$ with independent geometric variables $G_p$, $\Pr(G_p=k)=(1-p^{-\sigma})p^{-k\sigma}$; and $\Gamma_{1+\sigma/2}$ is gamma distributed with shape $1+\sigma/2$ and scale $1$. Then $\widetilde M_4\,4^{-K}$ has the same law as $\Gamma_{1+\sigma/2}/(\pi N_\sigma^2)=\frac{\Gamma_{1+\sigma/2}}\pi\prod_pp^{-2G_p}$.

Source: [`ch/maxwald.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald.tex), label `mw:thm:khinchin`.

<a id="wald_couples_of_maximum_of_exponentials_and_logistic_law"></a>

### Wald couples of the maximum of exponential variables and the logistic law

Theorem 14.4.2, page 660 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Call a pair $(\tau,\varsigma)$ of independent random variables with $\tau\ge0$ a Wald couple if $\mathbb{E}e^{s\varsigma}<\infty$ and $\mathbb{E}e^{s\varsigma}\,\mathbb{E}e^{-s^2\tau/2}=1$ for all real $s$; equivalently, with $Z$ standard normal and independent of $\tau$, the characteristic functions $\varphi_1$ of $\varsigma$ and $\varphi_2$ of $\sqrt\tau\,Z$ form a van Dantzig pair, $\varphi_2(t)=1/\varphi_1(it)$. Let $M_4$ be the maximum of independent exponential random variables with one variable of each rate $\pi n$, $n\ge1$, and one more of each rate $\pi(2m-1)$, $m\ge1$, so that $\mathbb{E}e^{-\lambda M_4}=\sqrt{\pi\lambda}/\sinh\sqrt{\pi\lambda}$, let $M_4'$ be an independent copy, let $\upsilon,\upsilon_1,\upsilon_2$ be uniform on $(-1,1)$, independent of each other and of $M_4,M_4'$, and let $c=\sqrt{\pi/2}$. Then $(M_4,c\upsilon)$ and $(M_4+M_4',c(\upsilon_1+\upsilon_2))$ are Wald couples. The associated van Dantzig pairs are $(\sin(ct)/(ct),\ ct/\sinh(ct))$, the characteristic functions of $c\upsilon$ and of $\sqrt{M_4}\,Z$ with $Z$ standard normal and independent of $M_4$, and their squares, and $\sqrt{M_4}\,Z$ has the logistic law, with characteristic function $ct/\sinh(ct)$.

Source: [`ch/maxwald.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald.tex), label `mw:thm:wald`.

<a id="stieltjes_criterion_for_rh_through_smoothed_thorin_density_of_clock"></a>

### Stieltjes criterion for RH through the smoothed Thorin density of the clock

Theorem 14.6.1, page 661 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$, let $\rho$ run over the nontrivial zeros of $\zeta$ with $\operatorname{Im}\rho>0$, counted with multiplicity, put $a_\rho=-(\rho-\frac12)^2$ and $G(\lambda)=\xi(\frac12+\sqrt\lambda)/\xi(\frac12)=\prod_\rho(1+\lambda/a_\rho)$. Let $T$ be the clock, the random variable with $\mathbb{E}e^{-\lambda T}=1/G(\lambda)$, which exists with no hypothesis on the zeros. Put $W(x)=\sum_\rho e^{-a_\rho x}$ for $x>0$, so that $(G'/G)(\lambda)=\sum_\rho(\lambda+a_\rho)^{-1}=\int_0^\infty e^{-\lambda x}W(x)\,dx$, and $r_\eta(x)=\frac1\pi\operatorname{Im}\frac{G'}G(-x-i\eta)$ for $x\in\mathbb{R}$ and $\eta>0$, wherever $-x-i\eta$ is not a pole of $G'/G$. The following are equivalent. (a) The Riemann Hypothesis (RH). (b) $T$ is a generalised gamma convolution (GGC), that is $\mathbb{E}e^{-\lambda T}=\exp\big(-a\lambda-\int_{(0,\infty)}\log(1+\lambda/x)\,U(dx)\big)$ for some $a\ge0$ and positive measure $U$, the Thorin measure; then $T$ has the same law as $\sum_{\gamma>0}e_\gamma/\gamma^2$ with independent standard exponential variables $e_\gamma$, and its Thorin measure is $U_{\mathrm{RH}}=\sum_{\gamma>0}\delta_{\gamma^2}$, the ordinates $\gamma=\operatorname{Im}\rho$ counted with multiplicity. (c) $G'/G$ is a Stieltjes function, that is of the form $b+\int_{[0,\infty)}\sigma(dx)/(\lambda+x)$ with $b\ge0$ and $\sigma$ a positive measure. (d) $r_\eta(x)\ge0$ for all $x\in\mathbb{R}$ and $\eta>0$ for which $-x-i\eta$ is not a pole of $G'/G$. (e) $W$ is completely monotone. In that case $r_\eta(x)=\sum_{\gamma>0}\eta/(\pi((x-\gamma^2)^2+\eta^2))$ is the Poisson smoothing of $U_{\mathrm{RH}}$ at scale $\eta$ and $U_{\mathrm{RH}}((0,x])=N(\sqrt x)$, where $N(t)$ is the number of zeros $\rho$ with $0<\operatorname{Im}\rho\le t$, counted with multiplicity.

Source: [`ch/maxwald.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald.tex), label `mw:thm:stieltjes`.

<a id="off_line_zero_makes_smoothed_thorin_density_of_clock_negative"></a>

### An off line zero makes the smoothed Thorin density of the clock negative

Lemma 14.6.2, page 662 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$, let $\rho'$ run over the nontrivial zeros of $\zeta$ with $\operatorname{Im}\rho'>0$, counted with multiplicity, put $a_{\rho'}=-(\rho'-\frac12)^2$, $G(\lambda)=\xi(\frac12+\sqrt\lambda)/\xi(\frac12)=\prod_{\rho'}(1+\lambda/a_{\rho'})$, so that $1/G$ is the Laplace transform of the clock, and $r_\eta(x)=\frac1\pi\operatorname{Im}\frac{G'}G(-x-i\eta)$ for $x\in\mathbb{R}$ and $\eta>0$. Let $\rho=\beta+i\gamma$ be a zero of multiplicity $m$ with $\beta>\frac12$ and $\gamma>0$, and put $\delta=\beta-\frac12$, $x_\rho=\gamma^2-\delta^2$ and $b=2\delta\gamma$, so that $a_\rho=x_\rho-ib$. Then $r_\eta(x_\rho)\to-\infty$ as $\eta\uparrow b$. More precisely, for $0<\eta<b$ such that $-x_\rho-i\eta$ is not a pole of $G'/G$, $\pi r_\eta(x_\rho)=\eta\Big(A_\eta-\frac{2m}{b^2-\eta^2}\Big)$, where $A_\eta=\frac1\eta\sum_{\rho'}\operatorname{Im}\,(a_{\rho'}-x_\rho-i\eta)^{-1}$, the sum running over the zeros $\rho'$ with $\operatorname{Im}\rho'>0$, counted with multiplicity, whose $a_{\rho'}$ is neither $a_\rho$ nor $\overline{a_\rho}$. Hence $r_\eta(x_\rho)<0$ at every $\eta\in(0,b)$ with $A_\eta<2m/b^2$. If all the other zeros lie on the critical line, $A_\eta\to A_0=\sum_{\gamma'}(\gamma'^2-x_\rho)^{-2}$ as $\eta\downarrow0$, the sum running over the ordinates $\gamma'>0$ of those zeros, counted with multiplicity.

Source: [`ch/maxwald.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald.tex), label `mw:lem:detect`.

<a id="explicit_ggc_within_laplace_distance_of_the_clock"></a>

### An explicit GGC within a stated Laplace distance of the clock

Proposition 14.7.2, page 664 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$, let $\rho_k=\beta_k+i\gamma_k$, $k\ge1$, be the nontrivial zeros of $\zeta$ with $\gamma_k>0$ in nondecreasing order of $\gamma_k$, counted with multiplicity, and put $a_k=-(\rho_k-\frac12)^2$, $G(\lambda)=\xi(\frac12+\sqrt\lambda)/\xi(\frac12)=\prod_k(1+\lambda/a_k)$ and $G_K(\lambda)=\prod_{k\le K}(1+\lambda/\gamma_k^2)$. Let $T$ be the clock, with $\mathbb{E}e^{-\lambda T}=1/G(\lambda)$. For laws $P,Q$ on $[0,\infty)$ with Laplace transforms $L_P,L_Q$ let $d_L(P,Q)=\sup_{\lambda\ge0}|L_P(\lambda)-L_Q(\lambda)|$. Let $H=3\times10^{12}$ and let $N$ be the number of zeros with $0<\gamma_k\le H$, all of which lie on the critical line by the verification of RH up to height $H$ of Platt and Trudgian. Let $T^{[N]}=\sum_{k\le N}e_k/\gamma_k^2$ with independent standard exponential variables $e_k$, the GGC with Laplace transform $1/G_N$, and put $\varepsilon_H=\sum_{k>N}1/|a_k|$ and $S_K=\sup_{\lambda\ge0}\lambda/G_K(\lambda)$ for $K\le N$. Then $d_L\big(T,T^{[N]}\big)\le S_K\,\varepsilon_H$ for every $K\le N$, with $\varepsilon_H<1.4798\times10^{-12}$ and $S_{1700}<16.85$, so that $d_L\big(T,T^{[N]}\big)<2.5\times10^{-11}$. The proof uses only the verification of RH up to height $H$ and makes no assumption on the zeros above $H$.

Source: [`ch/maxwald2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald2.tex), label `mw:prop:approx`.

<a id="term_of_one_zero_in_thorin_density_of_basepoint_family_lambda_chart"></a>

### Term of one zero in the Thorin density of the basepoint family in the lambda chart

Proposition 14.8.1, page 666 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be Riemann's xi function, $\alpha>\frac12$, $P_\delta(\theta)=\delta/(\pi(\delta^2+\theta^2))$ the Poisson kernel of width $\delta$ (so that $P_{-a}=-P_a$), and $u_\alpha(\lambda)=\operatorname{Re}\frac{\xi'}\xi(\alpha+i\sqrt\lambda)/(2\pi\sqrt\lambda)$ for $\lambda>0$, the density in the variable $\lambda$ of the Thorin measure of the shifted clock $T_\alpha$, the law with Laplace transform $\xi(\alpha)/\xi(\alpha+\sqrt\lambda)$, when every zero has real part below $\alpha$. Let $\rho=\beta+i\gamma$ be a nontrivial zero of $\zeta$ with $\gamma>0$ and put $\delta=\alpha-\beta\ne0$ and $c=\alpha-\rho=\delta-i\gamma$. The term of $\rho$ in $u_\alpha=\sum_{\gamma>0}u_\rho$, where the zeros $\rho$ and $\bar\rho$ are paired, is $u_\rho(\lambda)=\big(P_\delta(\sqrt\lambda-\gamma)+P_\delta(\sqrt\lambda+\gamma)\big)/(2\sqrt\lambda)$ for $\lambda>0$, and it satisfies the following. (a) $u_\rho(\lambda)=\kappa_\rho(\lambda)C_\rho(\lambda)=\delta(\lambda+|c|^2)/(\pi\sqrt\lambda\,|c^2+\lambda|^2)$, where $\kappa_\rho(\lambda)=(\lambda+\gamma^2+\delta^2)/(2\gamma\sqrt\lambda)$ and $C_\rho(\lambda)=\frac1\pi\cdot2\delta\gamma/\big((\lambda-\gamma^2+\delta^2)^2+4\delta^2\gamma^2\big)$; for $\delta>0$, $C_\rho$ is the Cauchy density centred at $\gamma^2-\delta^2$ with half width $2\delta\gamma$. (b) $\int_0^\infty u_\rho(\lambda)\,d\lambda=\operatorname{sign}\delta$. (c) $\kappa_\rho(\lambda)\ge\sqrt{1+\delta^2/\gamma^2}>1$ for every $\lambda>0$, $\kappa_\rho(\gamma^2-\delta^2)=(1-\delta^2/\gamma^2)^{-1/2}$ when $|\delta|<\gamma$, and for $K\ge1$ and $|\delta|\le\gamma/(2K)$, $|\kappa_\rho(\lambda)-1|\le(K^2/2+1)\,\delta^2/\gamma^2$ for $|\lambda-\gamma^2|\le K|\delta|\gamma$. Hence near its peak the term is a Lorentzian of full width $4|\delta|\gamma$ at half height, with the sign of $\delta$, up to a relative error $O(\delta^2/\gamma^2)$, while $u_\rho(\lambda)\sim\delta/(\pi(\gamma^2+\delta^2)\sqrt\lambda)$ as $\lambda\downarrow0$ and $u_\rho(\lambda)\sim\delta/(\pi\lambda^{3/2})$ as $\lambda\to\infty$.

Source: [`ch/maxwald2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald2.tex), label `mw:prop:lambda`.

<a id="clock_is_not_a_maximum_of_independent_exponential_variables"></a>

### The clock is not a maximum of independent exponential variables

Proposition 14.9.1, page 671 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be Riemann's xi function, $T$ the clock, the random variable with $\mathbb{E}e^{-\lambda T}=\xi(\frac12)/\xi(\frac12+\sqrt\lambda)$, and $F_T(t)=\Pr(T\le t)$. Then $-\log F_T$ is not completely monotone on $(0,\infty)$. Hence $T$ does not have the law of a maximum of independent exponential variables $E_i/r_i$, with $E_i$ standard exponential and deterministic rates $r_i>0$ with $\sum_ie^{-r_it}<\infty$ for $t>0$, for any choice of the rates, nor of the maximum of a Poisson cloud of exponential clocks, that is of the variables $E_i/x_i$ attached to the points $x_i$ of a Poisson process of intensity $\nu$ on $(0,\infty)$; for both of these laws $-\log F(t)=\int e^{-xt}\,\nu(dx)$ for a positive measure $\nu$, which is completely monotone in $t$. The proof uses the verification of RH up to height $3\times10^{12}$ and the simplicity of the zeros of height below $200$.

Source: [`ch/maxwald2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald2.tex), label `mw:prop:nomax`.

<a id="gumbel_limit_for_maxima_of_independent_copies_of_the_clock"></a>

### Gumbel limit for maxima of independent copies of the clock and of the maximum M_4

Proposition 14.10.1, page 672 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi$ be Riemann's xi function, $\Xi(\theta)=\xi(\frac12+i\theta)$, $\gamma_1$ the least positive ordinate of a zero of $\zeta$ (the zero $\frac12+i\gamma_1$ is simple) and $R_1=-2\gamma_1\Xi(0)/\Xi'(\gamma_1)$, the residue at $\lambda=-\gamma_1^2$ of $\xi(\frac12)/\xi(\frac12+\sqrt\lambda)$, the Laplace transform of the clock $T$. Let $M_4$ be the maximum of independent exponential random variables with one variable of each rate $\pi n$, $n\ge1$, and one more of each rate $\pi(2m-1)$, $m\ge1$, and let $M_4^{(j)}$ and $T^{(j)}$, $j\le n$, be independent copies of $M_4$ and of $T$. (i) $\Pr(M_4>t)\sim2e^{-\pi t}$ as $t\to\infty$, and $\Pr\big(\pi\max_{j\le n}M_4^{(j)}-\log(2n)\le x\big)\to\exp(-e^{-x})$ as $n\to\infty$, for every real $x$. (ii) $\Pr(T>t)\sim(R_1/\gamma_1^2)\,e^{-\gamma_1^2t}$ as $t\to\infty$, $\Pr\Big(\gamma_1^2\max_{j\le n}T^{(j)}-\log\frac{R_1n}{\gamma_1^2}\le x\Big)\to\exp(-e^{-x})$ for every real $x$, and $\max_{j\le n}T^{(j)}/\log n\to1/\gamma_1^2$ in probability.

Source: [`ch/maxwald2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald2.tex), label `mw:prop:gnedenko`.

<a id="finite_maxima_converge_to_maximum_and_finite_stages_to_xi"></a>

### Convergence of the finite maxima and of the finite stages of xi

Theorem 14.12.1, page 676 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $(E_n)_{n\ge1}$ and $(E'_m)_{m\ge1}$ be independent standard exponential random variables, for $N\ge1$ let $M_4^{(N)}=\max\big\{E_n/(\pi n)\ (n\le N),\ E'_m/(\pi(2m-1))\ (m\le N)\big\}$, and let $M_4=\lim_NM_4^{(N)}$ be the maximum of all the variables $E_n/(\pi n)$ and $E'_m/(\pi(2m-1))$, $n,m\ge1$. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$ be Riemann's xi function, and define the finite stage $\xi_N(s)=\frac{(s-1)\,\mathbb{E}(M_4^{(N)})^{s/2}}{2(1-2^{1-s})}$ and its symmetrisation $\xi_N^{\mathrm{sym}}(s)=\frac12\big(\xi_N(s)+\xi_N(1-s)\big)$, so that the same expression with $M_4$ in place of $M_4^{(N)}$ is $\xi(s)$ (here $\xi_N$ is the finite stage, not $\xi$ evaluated at $N$). Then (i) $M_4^{(N)}\uparrow M_4$ almost surely, and $\Pr(M_4^{(N)}\ne M_4)\le C\sqrt N\,e^{-\pi\sqrt N/2}$ for a constant $C$ independent of $N$. (ii) For every compact set $K$ in the open strip $0<\operatorname{Re}s<1$ there is a constant $C_K$ with $\sup_K|\xi_N-\xi|\le C_K\,N^{1/4}e^{-\pi\sqrt N/4}$ and $\sup_K|\xi_N^{\mathrm{sym}}-\xi|\le C_K\,N^{1/4}e^{-\pi\sqrt N/4}$ for all $N$.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:thm:finite`.

<a id="hurwitz_zeros_of_finite_stages_and_sufficient_condition_for_rh"></a>

### Zeros of the symmetrised finite stages and a sufficient condition for RH

Corollary 14.12.2, page 676 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $M_4^{(N)}$ be the maximum of the independent exponential variables $E_n/(\pi n)$ ($n\le N$) and $E'_m/(\pi(2m-1))$ ($m\le N$), with $E_n,E'_m$ independent standard exponential variables, let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$ be Riemann's xi function, $\xi_N(s)=\frac{(s-1)\,\mathbb{E}(M_4^{(N)})^{s/2}}{2(1-2^{1-s})}$ the finite stage and $\xi_N^{\mathrm{sym}}(s)=\frac12\big(\xi_N(s)+\xi_N(1-s)\big)$ its symmetrisation. (a) Every zero of $\xi$ in the open strip $0<\operatorname{Re}s<1$ is a limit of zeros of $\xi_N$, and of zeros of $\xi_N^{\mathrm{sym}}$, with multiplicity, and every limit of zeros of $\xi_N$ inside the open strip is a zero of $\xi$. (b) If $\rho=\frac12+i\gamma$ is a simple zero of $\xi$, then for all large $N$ the function $\xi_N^{\mathrm{sym}}$ has exactly one zero in a small disc about $\rho$, and that zero lies on the critical line. (c) Suppose that there are $N_j\to\infty$ and heights $h_j\to\infty$ such that $\xi^{\mathrm{sym}}_{N_j}$ has no zeros in $\{\frac12<\operatorname{Re}s<1,\ |\operatorname{Im}s|<h_j\}$ for every $j$. Then the Riemann Hypothesis holds. In particular this is so if there are $c>0$ and $N_0$ such that for every $N\ge N_0$ the function $\xi^{\mathrm{sym}}_N$ has no zeros in $\{\frac12<\operatorname{Re}s<1,\ |\operatorname{Im}s|<c\sqrt N\}$.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:cor:hurwitz`.

<a id="finite_order_ggc_bondesson_and_hcm_classes_of_the_clock"></a>

### Finite order GGC, Bondesson and HCM classes of the clock, and RH as the infinite order case

Theorem 14.13.2, page 683 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(\frac s2)\zeta(s)$ and let $T$ be the clock, with Laplace transform $\varphi(\lambda)=\xi(\frac12)/\xi(\frac12+\sqrt\lambda)$ and Lévy density $W(x)/x$, where $W(x)=\sum_\rho e^{(\rho-1/2)^2x}$ is summed over the nontrivial zeros $\rho$ of $\zeta$ with $\operatorname{Im}\rho>0$, counted with multiplicity. A function $f$ on an interval is $k$ monotone if $(-1)^jf^{(j)}\ge0$ for $0\le j\le k$, and completely monotone if it is $k$ monotone for every $k$. For $k\ge1$ let $\mathrm{GGC}_k$ be the class of laws on $[0,\infty)$ with a Lévy density $\ell$ such that $x\ell(x)$ is $k$ monotone, and $\mathrm{BO}_k$ the class with $\ell$ itself $k$ monotone, so that $\mathrm{GGC}_k\subset\mathrm{BO}_k$, $\mathrm{GGC}_1$ is the self decomposable class, $\bigcap_k\mathrm{GGC}_k$ is the class GGC of generalised gamma convolutions and $\bigcap_k\mathrm{BO}_k=\mathrm{BO}$ is the Bondesson class. A positive function $\varphi$ on $(0,\infty)$ is HCM (hyperbolically completely monotone) of order $k$ if for every $u>0$ the function $w\mapsto\varphi(uv)\varphi(u/v)$, $w=v+1/v\ge2$, is $k$ monotone on $[2,\infty)$, and HCM if it is HCM of every order. Let $H=3\times10^{12}$. Statements (i) and (ii) are unconditional: they use the verification of RH up to height $H$ of Platt and Trudgian. (i) $T\in\mathrm{GGC}_k$, and hence $T\in\mathrm{BO}_k$, for every $k\le1.16\times10^{25}$. (ii) $\varphi$ is HCM of order $k$ for every $k\le9\times10^{12}$. (iii) With no hypothesis, the following are equivalent: RH; $T$ is a GGC; $T\in\mathrm{BO}$; $\varphi$ is HCM; $\varphi$ is HCM of every order.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:thm:bondesson`.

<a id="finite_arithmetic_form_of_basepoint_density_below_basepoint_one"></a>

### Finite arithmetic form of the basepoint density below basepoint one

Proposition 14.14.1, page 687 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: conditional on Part (b) assumes the Riemann hypothesis. Part (a) holds without any hypothesis, and part (c) holds without any hypothesis when delta is at least one half and otherwise needs that every zero lies within delta of the critical line, which the Riemann hypothesis gives for every delta above zero.

Statement. Let $\xi$ be Riemann's xi function, $\Lambda$ the von Mangoldt function, $\psi$ the digamma function, $\alpha=\frac12+\delta$ with $\delta>0$, $t\in\mathbb{R}$ and $X>0$. Put $g_X(x)=e^{-\delta|x|}\cos(tx)\,(1-|x|/X)_+$ and $h_X(r)=\int_{\mathbb{R}}g_X(x)e^{irx}\,dx$, so that $h_X\ge0$ on $\mathbb{R}$: it is $(2\pi)^{-1}$ times the convolution of the two sided Poisson kernel $h_\infty(r)=\frac\delta{\delta^2+(t-r)^2}+\frac\delta{\delta^2+(t+r)^2}$ with the Fejér kernel $4\sin^2(Xr/2)/(Xr^2)$. Define $A_X(\delta,t)=\frac12\Big[h_X(\tfrac i2)+h_X(-\tfrac i2)-\log\pi+\frac1{2\pi}\int_{\mathbb{R}}h_X(r)\operatorname{Re}\psi\big(\tfrac14+\tfrac{ir}2\big)\,dr-2\sum_{n<e^X}\frac{\Lambda(n)}{n^{\alpha}}\cos(t\log n)\Big(1-\frac{\log n}X\Big)\Big]$, where $h_X(\tfrac i2)+h_X(-\tfrac i2)=2\operatorname{Re}\big[J(-\delta+\tfrac12+it)+J(-\delta-\tfrac12+it)\big]$ with $J(c)=-1/c+(e^{cX}-1)/(Xc^2)$. For a nontrivial zero $\rho=\beta+i\gamma$ of $\zeta$ put $\gamma_\rho=(\rho-\frac12)/i$. Then (a) unconditionally, $A_X(\delta,t)=\frac12\sum_\rho h_X(\gamma_\rho)$, the sum over all nontrivial zeros. (b) Under RH, $A_X(\delta,t)\ge0$ for every $X$, $\delta$ and $t$; in particular a single triple $(\delta,t,X)$ with $A_X(\delta,t)<0$ would disprove RH. (c) If every zero has $|\beta-\frac12|<\delta$ (under RH, or unconditionally when $\delta\ge\frac12$), then $A_X(\delta,t)\to\operatorname{Re}\frac{\xi'}\xi(\alpha+it)$ as $X\to\infty$.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:prop:AX`.

<a id="open_zero_freeness_of_symmetric_finite_stages_of_the_maximum"></a>

### Zero freeness of the symmetric finite stages of the maximum

Problem 14.1, page 690 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function. For $N\ge1$ let $M_4^{(N)}$ be the maximum of the independent exponential variables $E_n/(\pi n)$ ($n\le N$) and $E'_m/(\pi(2m-1))$ ($m\le N$), with $E_n,E'_m$ independent standard exponential variables, so that $\Pr(M_4^{(N)}\le t)=\sum_{m\ge0}p^{(N)}_me^{-\pi mt}$ where $\sum_{m\ge0}p^{(N)}_mq^m=\prod_{n\le N}(1-q^n)\prod_{m\le N}(1-q^{2m-1})$, and put $D_N(s)=\sum_{m\ge1}p^{(N)}_mm^{-s/2}$, a Dirichlet polynomial with integer coefficients. Let $F_N(s)=\frac{s-1}2\,\mathbb{E}\big(M_4^{(N)}\big)^{s/2}=-\frac{s-1}2\,\Gamma(1+\frac s2)\pi^{-s/2}D_N(s)$ and $F^{\mathrm{sym}}_N(s)=\frac12\big(F_N(s)+F_N(1-s)\big)$; as $N\to\infty$, $F^{\mathrm{sym}}_N\to\mathfrak q\,\xi$ uniformly on compact subsets of the open strip $0<\operatorname{Re}s<1$, with $\mathfrak q(s)=1-2^{-s}-2^{s-1}$. Find heights $h_N\to\infty$ such that for infinitely many $N$ the function $F^{\mathrm{sym}}_N$ has no zeros with $\frac12<\operatorname{Re}s<1$ and $|\operatorname{Im}s|<h_N$, for instance $h_N=c\sqrt N$ with some $c>0$. Either statement implies the Riemann Hypothesis, by Hurwitz's theorem, because a zero of $\xi$ in that region would be a zero of the limit $\mathfrak q\,\xi$. Numerical counts of zeros by the argument principle show, for $N=8$, $16$, $32$, $48$, $64$ and $96$, that all zeros of $F^{\mathrm{sym}}_N$ with $-1\le\operatorname{Re}s\le2$ and $0<\operatorname{Im}s\le T^*(N)$ lie on the critical line, for a height $T^*(N)$ of order $N$ (at $N=128$ a pair outside the strip $0<\operatorname{Re}s<1$ appears at height $60.7$, while the strip has no off line zero below $130$), and the sharper conjecture is $T^*(N)\ge\kappa N$ for some $\kappa>0$ and all $N$, which would also imply RH. A mechanism for $T^*(N)$, presumably through the ratio $|D_N(\sigma+it)/D_N(1-\sigma+it)|$, is the first thing to find.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:prob:hurwitz`.

<a id="open_termwise_positive_continuation_of_prime_form_below_basepoint_one"></a>

### A termwise positive continuation of the prime form below basepoint one

Problem 14.2, page 691 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function. For $\alpha>1$ the density $v_\alpha(\theta)=\operatorname{Re}\frac{\xi'}\xi(\alpha+i\theta)$ has the sine squared prime form, through a signed canonical measure (atoms $\Lambda(n)n^{-\alpha}$ at $\log n$, with $\Lambda$ the von Mangoldt function, an archimedean density and a pole density) whose prime and gamma terms are nonnegative and whose pole terms are negative; its positivity, which is Lagarias positivity for $\alpha>1$, comes from the absence of zeros with real part at least one, not from the terms one by one. The book continues this form below $\alpha=1$ by compensation, and shows that the compensated sine squared sum converges for every $\alpha'>\alpha$ exactly when $\zeta$ has no zero with $\operatorname{Re}s>\alpha$; a smoothed finite form, built from a Fejér tapered sum over the primes with explicit pole terms and a gamma integral, continues the identity below $\alpha=1$ but is not positive term by term. Is there a continuation of the prime form below $\alpha=1$ that is positive term by term? The prime sum diverges there, but smoothed versions of it, with $n^{-\alpha}$ replaced by a weight that cuts off at $n\approx e^{1/\varepsilon}$, remain positive termwise, and the question is whether the loss from the cutoff can be controlled by $(\xi'/\xi)(\alpha)$, the gamma sum and the pole terms. A positive answer on $\alpha\in(\alpha_0,1]$ would prove the zero free region $\beta\le\alpha_0$, that is, that $\zeta$ has no zero with real part above $\alpha_0$.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:prob:lagarias`.

<a id="open_existence_of_half_flow_clock_for_small_positive_times"></a>

### Existence of the half flow clock for small positive times

Problem 14.3, page 691 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function, $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}\big(2\pi^2n^4e^{9u}-3\pi n^2e^{5u}\big)\exp(-\pi n^2e^{4u})$ the Pólya kernel in the rescaled variable of Rodgers and Tao, so that $H_0(z)=\int_0^\infty\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du=\frac18\xi(\frac12+\frac{iz}2)$, and $H_t(z)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$ its de Bruijn Newman deformation, with de Bruijn Newman constant $\Lambda$, the least $t$ from which all zeros of $H_t$ are real, known to satisfy $0\le\Lambda\le0.2$ by Rodgers and Tao, Platt and Trudgian, and the Polymath project. Let $R_t(u)=H_t(0)/H_t(2i\sqrt u)$, the Laplace transform of the clock of the flow, with $R_0(u)=\xi(\frac12)/\xi(\frac12+\sqrt u)$. For $t\ge\Lambda$, $R_t$ is the transform of a GGC with Thorin measure $\sum_k\delta_{(z_k(t)/2)^2}$, where $z_k(t)>0$ are the positive zeros of $H_t$, each pair $\pm z$ giving one atom; so the law $T_t$ with Laplace transform $\lambda\mapsto R_t(\lambda/2)$ exists and is a GGC for $t\ge0.2$ unconditionally and, under RH, for every $t\ge0$, and $T_0$ is half the clock, which exists unconditionally. For $0<t<0.2$, does $T_t$ exist as a probability law, and is it infinitely divisible? This is the analogue along the de Bruijn Newman flow of the existence of the clock; it would need an explicit formula for the heat trace of the flow, which the Gaussian tilt destroys.

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:prob:flow`.

<a id="open_khinchin_analogue_of_tilted_maximum_at_one_half"></a>

### An analogue of the Khinchin factorisation at one half

Problem 14.4, page 691 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi$ be Riemann's xi function and $\Xi(\theta)=\xi(\frac12+i\theta)$. For $\sigma>1$ the maximum $M_4$ of independent exponential variables with one variable of each rate $\pi n$ and one more of each rate $\pi(2m-1)$, tilted to the law $\propto x^{\sigma/2}\Pr(M_4\in dx)$ and multiplied by an independent factor $4^{-K}$ with $\Pr(K=k)=(1-a)a^k$ for $k\ge0$ and $a=2^{1-\sigma}$, has the law of $\Gamma_{1+\sigma/2}/(\pi N_\sigma^2)$, where $\Gamma_{1+\sigma/2}$ is a gamma variable of shape $1+\sigma/2$ and scale $1$ and $N_\sigma$ is an independent zeta distributed integer, $\Pr(N_\sigma=n)=n^{-\sigma}/\zeta(\sigma)$, built from the primes by the Euler product. Is there an analogue of this factorisation at $\sigma=\frac12$, that is, a factorisation of the tilted law of the clock $T$ (the random variable with Laplace transform $\xi(\frac12)/\xi(\frac12+\sqrt\lambda)$) or of the Pólya variable (the variable with characteristic function $\Xi(\theta)/\Xi(0)$) through a law built from the primes, with the explicit formula replacing the Euler product?

Source: [`ch/maxwald3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/maxwald3.tex), label `mw:prob:khinchin`.

## Chapter 15. The Thorin phase diagram

<a id="threshold_theorem_shifted_clock_along_heat_flow"></a>

### Threshold theorem along the flow for the shifted clocks of the heat flow function

Theorem 15.2.3, page 700 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$, $H_t(z)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$ (so $H_0(z)=\frac18\xi(\frac12+\frac{iz}2)$) and $\xi_t(s)=8H_t(-2i(s-\frac12))=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, so that $\xi_0=\xi$. For $t\ge0$ and $\alpha\ge\frac12$ let $\varphi_{t,\alpha}(u)=\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ for $u\ge0$, the shifted time $t$ clock (for $\alpha=\frac12$ it is $R_t(u)=H_t(0)/H_t(2i\sqrt u)$, the transform of the centre clock of the flow). A generalized gamma convolution (GGC) is a random variable $Y\ge0$ with $\mathbb{E}e^{-uY}=\exp\{-au-\int_{(0,\infty)}\log(1+u/z)\,U(dz)\}$ for $u>0$, a drift $a\ge0$ and a Thorin measure $U$. Then $\varphi_{t,\alpha}$ is the Laplace transform of a GGC if and only if $\alpha\ge\sup_\rho\operatorname{Re}\rho$ over the zeros $\rho$ of $\xi_t$. In that case there is no drift, and the Thorin measure, written in the variable $\theta=\sqrt z$ on $\theta>0$, is $\sum_{\rho=\beta+i\gamma,\ \gamma>0,\ \beta<\alpha}\bigl[P_{\alpha-\beta}(\theta-\gamma)+P_{\alpha-\beta}(\theta+\gamma)\bigr]\,d\theta+\sum_{\rho=\alpha+i\gamma,\ \gamma>0}m_\rho\,\delta_\gamma$, where the sums run over the zeros $\rho=\beta+i\gamma$ of $\xi_t$ with $\gamma>0$ (in the first sum each zero is repeated according to its multiplicity $m_\rho$), $P_a(x)=a/(\pi(a^2+x^2))$ is the Cauchy density of scale $a>0$ and $\delta_\gamma$ is the point mass at $\gamma$. For $t=0$ this is the threshold theorem for the shifted clock of $\xi$.

Source: [`ch/cp2.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp2.tex), label `cp:thm:transfer`.

<a id="phase_diagram_theorem_de_bruijn_strip_theorem_in_thorin_form"></a>

### Phase diagram theorem, de Bruijn's strip theorem in Thorin form

Theorem 15.3.2, page 702 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ be Riemann's xi function, $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$, $H_t(z)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$ (so $H_0(z)=\frac18\xi(\frac12+\frac{iz}2)$) and $\xi_t(s)=8H_t(-2i(s-\frac12))$, so that $\xi_0=\xi$ and a zero $x+iy$ of $H_t$ is the zero $\frac12-\frac y2+\frac{ix}2$ of $\xi_t$. For $t\ge0$ and $\alpha\ge\frac12$ let $\varphi_{t,\alpha}(u)=\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ for $u\ge0$ (the shifted time $t$ clock), and let the GGC region $\mathcal{G}$ be the set of points $(t,\alpha)$ with $t\ge0$ and $\alpha\ge\frac12$ at which $\varphi_{t,\alpha}$ is the Laplace transform of a generalized gamma convolution; by the threshold theorem for the clocks of $\xi_t$ this holds exactly when $\xi_t$ has no zero with real part above $\alpha$, that is, when every zero $x+iy$ of $H_t$ has $|y|\le2(\alpha-\frac12)$. Let $\Lambda$ be the de Bruijn Newman constant, so that $H_t$ has only real zeros exactly when $t\ge\Lambda$. (i) If $(t_0,\alpha_0)\in\mathcal{G}$, then $(t,\alpha)\in\mathcal{G}$ for every $t\ge t_0$ and every $\alpha\ge\frac12$ with $(\alpha-\frac12)^2+\frac t2\ge(\alpha_0-\frac12)^2+\frac{t_0}2$. (ii) $\Lambda=\min\{t+2(\alpha-\frac12)^2:(t,\alpha)\in\mathcal{G}\}$. (iii) The Riemann Hypothesis (every zero of $\xi$ has real part $\frac12$) holds if and only if $(0,\frac12)\in\mathcal{G}$, and if and only if for every $\delta>0$ there is $t\in[0,\delta]$ with $(t,\frac12+\sqrt{(\delta-t)/2})\in\mathcal{G}$. Part (i) is de Bruijn's strip theorem (if the zeros of $H_{t_0}$ lie in $|\operatorname{Im}z|\le\Delta$, then the zeros of $H_t$ lie in $(\operatorname{Im}z)^2\le\max(\Delta^2-2(t-t_0),0)$ for $t\ge t_0$) read through the threshold theorem, and parts (ii) and (iii) also use the theorem of Rodgers and Tao that $\Lambda\ge0$.

Source: [`ch/cp3.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp3.tex), label `cp:thm:phase`.

<a id="winding_number_certificate_for_point_0_2_0_6_of_phase_diagram"></a>

### Winding number certificate for the point (0.2, 0.6) of the phase diagram

Theorem 15.5.6, page 708 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, the de Bruijn Newman deformation of Riemann's xi function $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ (so $\xi_0=\xi$), and $H_t(z)=\frac18\xi_t(\frac12+\frac{iz}2)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$. The kernel $\Phi_{\mathrm{dBN}}$ continues analytically to $|\operatorname{Im}u|<\pi/8$. For $0<\eta<\pi/8$ and real $x$ put $K_\eta(0.2,x)=4\int_{-\infty}^{\infty}\bigl|e^{0.2(u+i\eta)^2}\bigr|\,|\Phi_{\mathrm{dBN}}(u+i\eta)|\,e^{2xu}\,du$ and $\bar K_\eta(0.2,R)=\max_{|x|\le R}K_\eta(0.2,x)$; they give the bound $|\xi_{0.2}(\frac12+x+i\theta)|\le K_\eta(0.2,x)e^{-2\eta\theta}$ for $\theta\ge0$. Let $X>0$ be a barrier for the point $(0.2,0.6)$ in the variable of $H_{0.2}$, as Polymath15 provide one (no script of the book checks it), that is, a number such that $H_{0.2}$ has no zero $z$ with $\operatorname{Re}z\ge X$ and $\operatorname{Im}z>0.2$, so that $\xi_{0.2}$ has no zero with $\operatorname{Re}s>0.6$ and $\operatorname{Im}s\ge X/2$, and let $Q=\{s:0.6\le\operatorname{Re}s\le0.9,\ 0\le\operatorname{Im}s\le X/2\}$. Let $s_0,\dots,s_M$ be points of the left, top and right edges of $Q$, in order along the boundary of $Q$ from one bottom corner, $s_0$, to the other, $s_M$, and including the two top corners, and let $B_j$ be ball enclosures of $\xi_{0.2}(s_j)$ with radii $\varepsilon_j$. Suppose that for every $j$ there are $r>0$ and a tilt $\eta\in(0,\pi/8)$ (both may depend on $j$) such that $|s_{j+1}-s_j|\,\frac{\bar K_\eta(0.2,\,0.4+r)}{r}\,e^{-2\eta(\theta_j-r)}+\varepsilon_j+\varepsilon_{j+1}<\min|B_j|$, where $\theta_j=\min(\operatorname{Im}s_j,\operatorname{Im}s_{j+1})$, $\min|B_j|$ is the smallest modulus of a point of the ball $B_j$ and the factor $e^{-2\eta(\theta_j-r)}$ may be replaced by $1$ when $\theta_j\le r$, and that $\sum_j\operatorname{Arg}(B_{j+1}/B_j)$ (principal values), evaluated in ball arithmetic, lies within $(-\pi,\pi)$. (The bottom edge needs no points: $\xi_t(\sigma)>0$ for every real $\sigma$ and every $t\ge0$, so the argument of $\xi_{0.2}$ does not change along it.) Then $\xi_{0.2}$ has no zero in $Q$, hence no zero with real part above $0.6$: the point $(0.2,0.6)$ lies in the GGC region (the shifted clock $\xi_{0.2}(0.6)/\xi_{0.2}(0.6+\sqrt u)$ is the Laplace transform of a generalized gamma convolution), and the de Bruijn Newman constant satisfies $\Lambda\le0.22$.

Source: [`ch/cp5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp5.tex), label `cp:thm:cert`.

<a id="thorin_form_certificate_positive_log_derivative_on_rectangle"></a>

### Certificate in Thorin form, positive logarithmic derivative on a rectangle excludes zeros

Proposition 15.5.7, page 709 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, the de Bruijn Newman deformation of Riemann's xi function, and for real $\sigma,\theta$ with $\xi_t(\sigma+i\theta)\ne0$ put $v_t(\sigma,\theta)=\operatorname{Re}\bigl(\xi_t'(\sigma+i\theta)/\xi_t(\sigma+i\theta)\bigr)$ (when the shifted clock $\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ is a generalized gamma convolution, $v_t(\alpha,\theta)/\pi$ is its Thorin density in the variable $\theta=\sqrt z$, away from atoms). Let $t\ge0$, $\frac12<\alpha<\sigma_1$ and $T>0$, and suppose $v_t(\sigma,\theta)>0$ for every $(\sigma,\theta)\in[\alpha,\sigma_1]\times[0,T]$ at which $\xi_t(\sigma+i\theta)\ne0$. Then $\xi_t$ has no zero $\beta+i\gamma$ with $\alpha<\beta\le\sigma_1$ and $0<\gamma<T$.

Source: [`ch/cp5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp5.tex), label `cp:prop:thorincert`.

<a id="xi_time_0_2_no_zero_in_rectangle_up_to_height_160_certified"></a>

### No zero of the time 0.2 function in the rectangle up to height 160, certified

Theorem 15.5.8, page 711 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, the de Bruijn Newman deformation of Riemann's xi function $\xi(s)=\frac12s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$ (so $\xi_0=\xi$). Then $\xi_{0.2}$ has no zero $s$ with $0.6\le\operatorname{Re}s\le0.9$ and $0\le\operatorname{Im}s\le160$. This is certified by a computation in ball arithmetic (Arb, called from Python): $\xi_{0.2}$ is evaluated by the trapezoid rule for its integral with an enclosure of the error of the rule, and the boundary of the rectangle $[0.6,0.9]\times[0,160]$ is walked counterclockwise, up the right edge, along the top and down the left edge, in $23{,}219$ steps, with the step condition of the winding number certificate for the point $(0.2,0.6)$ (the step length times an explicit bound for $|\xi_{0.2}'|$ along the step, plus the radii of the enclosures at its two ends, is below the smallest modulus of the enclosure at its start) verified at every step. The bottom edge is not walked, since $\xi_{0.2}$ is real and positive on the real axis. The winding number of $\xi_{0.2}$ around the rectangle is enclosed in $[0\pm3.47\times10^{-7}]$, so it is $0$.

Source: [`ch/cp5.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp5.tex), label `cp:thm:rig`.

<a id="heat_flow_zero_pair_near_height_185_still_real_at_minus_0_287979_certified"></a>

### Pair of zeros near height 185 still real at time minus 0.287979, certified

Proposition 15.6.1, page 715 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: proved.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and, for real $t$ of either sign, $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, the de Bruijn Newman deformation of Riemann's xi function (so $\xi_0$ is the xi function), and $H_t(z)=\frac18\xi_t(\frac12+\frac{iz}2)$. A zero $\frac12+i\gamma$ of $\xi_t$ on the critical line is the real zero $z=2\gamma$ of $H_t$, and while the real zeros are simple they move by $dz_k/dt=2\sum_{j\ne k}(z_k-z_j)^{-1}$. Let $\gamma_{71}<\gamma_{72}$ be the ordinates of the 71st and 72nd zeros of $\xi$ in the upper half plane, the pair near height $185$. At $t=-0.287979$ the signs of $\xi_t(\frac12+i\theta)$ at the $131$ points $\theta=184.6+0.01k$, $k=0,\dots,130$, certified in ball arithmetic (Arb) at $640$ bits from the trapezoid sum for the integral with an enclosure of its discretisation error, show two sign changes; hence the real function $\theta\mapsto\xi_t(\frac12+i\theta)$ has at least two zeros in $[184.6,185.9]$, and the pair $(\gamma_{71},\gamma_{72})$ is still real at $t=-0.287979$. This is the certified upper end of the book's bracket $[-0.28803,-0.287979]$ for the time at which the pair collides when $t$ decreases from $0$; the lower end $-0.28803$ rests on the asymptotic pair law $D^2\approx8(t-t_*)$ for the distance $D$ of the pair in the variable $z$ near its collision time $t_*$, and is not certified. Newton's method on the same trapezoid sum puts the double zero of the pair at $z_*=370.78087$, $t_*=-0.287979837009$; a Krawczyk test certifies a real double zero of the trapezoid sum itself (a statement about the sum, not about $H_t$) in a box of half width $2^{-400}$ about this point, and at $t=t_*+1.2\times10^{-12}$ certified signs $+,-,+$ of $H_t$ at $z_*-4\times10^{-6}$, $z_*$ and $z_*+4\times10^{-6}$ show that the pair is still real there, so the collision lies below $-0.287979837008$. The floating point integration of the motion of the first $600$ zeros that the book reports with this result (with the monotonicity of two gap functionals at the sampled times) is computed, not certified, and is not claimed here.

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prop:gas`.

<a id="open_thorin_densities_at_time_0_and_0_2_near_polymath15_barrier"></a>

### Open problem on the candidate Thorin densities at times 0 and 0.2 near the Polymath15 barrier

Problem 15.1, page 718 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$, the de Bruijn Newman deformation of Riemann's xi function (so $\xi_0$ is the xi function), with $H_t(z)=\frac18\xi_t(\frac12+\frac{iz}2)$, and put $v_t(\alpha,\theta)=\operatorname{Re}\bigl(\xi_t'(\alpha+i\theta)/\xi_t(\alpha+i\theta)\bigr)$; when the shifted clock $\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ is the Laplace transform of a generalized gamma convolution, $v_t(\alpha,\theta)/\pi$ is its Thorin density in the variable $\theta=\sqrt z$, away from atoms. Problem: compute $v_t(\alpha,\theta)$ at $t=0$ and at $t=0.2$ near the barrier of Polymath15, at heights $\theta$ near $3\times10^{10}$, with a Riemann Siegel type evaluator of $H_{0.2}$ such as the effective approximation of Polymath15, and measure how the difference between the functions $v_{0.2}/\pi$ and $v_0/\pi$ grows with the height. The book computes this difference only for $\theta\le100$, from the Pólya integral.

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prob:high`.

<a id="open_boundary_of_ggc_region_for_times_up_to_0_2"></a>

### Open problem on the boundary of the GGC region for 0 at most t at most 0.2

Problem 15.2, page 718 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$ be the de Bruijn Newman deformation of Riemann's xi function, with $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$. The GGC region $\mathcal{G}$ is the set of points $(t,\alpha)$, $t\ge0$, $\alpha\ge\frac12$, at which the shifted clock $\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ is the Laplace transform of a generalized gamma convolution; this holds exactly when $\alpha\ge\alpha^*(t)$, where $\alpha^*(t)=\max(\frac12,\sup\operatorname{Re}\rho)$ over the zeros $\rho$ of $\xi_t$. Problem: describe the boundary $\alpha^*(t)$ of $\mathcal{G}$ for $0\le t\le0.2$ as far as computation and the known zero free regions allow; is $\alpha^*$ continuous, and is $(\alpha^*(t)-\frac12)^2+\frac t2$ strictly decreasing in $t$ where $\alpha^*(t)>\frac12$? (By the phase diagram theorem it is nonincreasing there.)

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prob:boundary`.

<a id="open_phase_diagram_coordinate_for_symmetrised_truncations"></a>

### Open problem on a coordinate on the phase diagram for the symmetrised truncations

Problem 15.3, page 718 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\xi_t(s)=8\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cosh((2s-1)u)\,du$ be the de Bruijn Newman deformation of Riemann's xi function, with $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$, and let the GGC region $\mathcal{G}$ be the set of points $(t,\alpha)$, $t\ge0$, $\alpha\ge\frac12$, at which $\xi_t$ has no zero with real part above $\alpha$, equivalently at which the shifted clock $\xi_t(\alpha)/\xi_t(\alpha+\sqrt u)$ is the Laplace transform of a generalized gamma convolution. The symmetrised truncations are $K_n(w)=\pi^{-w}F_n(w)+\pi^{w-1/2}F_n(\frac12-w)$ with $F_n(w)=\mathbb{E}\bigl[\bigl(\sum_{k\le n}k^{-2}(e_k+e_k')\bigr)^w\bigr]$, where $e_k,e_k'$ are independent standard exponential variables; Conjecture S says that, for every $n$, all zeros of $K_n$ in the strip $0<\operatorname{Re}w<\frac12$ lie on the line $\operatorname{Re}w=\frac14$, and at $t=0$, $\alpha=\frac12$ the Hurwitz argument shows that this for infinitely many $n$ would imply the Riemann Hypothesis. The known part of $\mathcal{G}$ consists of the points $(0,1)$ (the critical strip), $(0.2,0.6)$ (Polymath15) and a point on the level set $t+2(\alpha-\frac12)^2=0.2$ (Polymath15 with the verification of Platt and Trudgian), conditionally on an announced zero free half plane also $(0,\frac78)$, and every point $(t,\alpha)$ with $\alpha\ge\frac12$, $t\ge t_0$ and $(\alpha-\frac12)^2+\frac t2\ge(\alpha_0-\frac12)^2+\frac{t_0}2$ for one of these points $(t_0,\alpha_0)$ (de Bruijn's strip theorem in Thorin form). Problem: give the symmetrised truncations a coordinate on the phase diagram: for some $(t,\alpha)$ outside the known part, is there a family of explicit truncations of $\xi_t$, built from finitely many terms of the Pólya integral or of a Dirichlet series and making no use of the zeros, as the $K_n$ are, that converges to $\xi_t$ uniformly on compact sets, together with heights $h_n\to\infty$ such that for infinitely many $n$ the $n$th truncation has no zero with $\operatorname{Re}s>\alpha$ and $|\operatorname{Im}s|<h_n$? Hurwitz's theorem would then make $(t,\alpha)$ a new point of $\mathcal{G}$; at $t=0$ and $\alpha=\frac12$ the question for the $K_n$ is Conjecture S, and for the finite stages of the Pólya density as a maximum of exponential variables it is the open problem posed with those stages.

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prob:trunc`.

<a id="open_cholesky_factor_of_weil_gram_matrix_on_hat_functions"></a>

### Open problem on the Cholesky factor of the Weil Gram matrix

Problem 15.4, page 719 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. For $h$ smooth enough, with $\hat h(r)=\int h(x)e^{irx}\,dx$, let $W(h)=\hat h(\frac i2)+\hat h(-\frac i2)-\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\bigl(h(\log n)+h(-\log n)\bigr)+\frac1{2\pi}\int_{\mathbb{R}}\hat h(r)\operatorname{Re}\psi\bigl(\frac14+\frac{ir}2\bigr)\,dr-h(0)\log\pi$, where $\Lambda$ is von Mangoldt's function and $\psi=\Gamma'/\Gamma$; by Weil's explicit formula $W(h)=\sum_\rho\hat h(\gamma_\rho)$ over the zeros $\rho=\frac12+i\gamma_\rho$ of Riemann's xi function. For $g$ supported in $[-L/2,L/2]$ and $h=g*\tilde g$ with $\tilde g(x)=\overline{g(-x)}$, under the Riemann Hypothesis (every $\gamma_\rho$ real) $W(h)=\sum_\rho|\hat g(\gamma_\rho)|^2\ge0$ (Weil's criterion). Take the hat functions $g_j(x)=\max(0,1-|x-c_j|/w)$ at nodes $c_j$ spaced $w$ apart inside $[-L/2,L/2]$, $N=L/w-1$ functions in all, and the Weil Gram matrix $G_{jk}=W(g_j*\tilde g_k)$, which is Hermitian. Problem: compute the Cholesky factor of $G$ for decreasing $w$ and increasing $L$, and decide whether its entries converge to an object built from the arithmetic, as the first test of a factorisation $W(g*\tilde g)=\|Ag\|^2$ with $A$ built from the arithmetic.

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prob:chol`.

<a id="open_global_monotone_functional_of_zeros_along_heat_flow"></a>

### Open problem on a global functional of the zeros that is monotone along the flow

Problem 15.5, page 719 of Polson and Sokolov (2026), *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, edition 1. New in this book.

Status: open.

Statement. Let $\Phi_{\mathrm{dBN}}(u)=\sum_{n\ge1}(2\pi^2n^4e^{9u}-3\pi n^2e^{5u})\exp(-\pi n^2e^{4u})$ and $H_t(z)=\int_0^\infty e^{tu^2}\Phi_{\mathrm{dBN}}(u)\cos(zu)\,du$, the de Bruijn Newman deformation of Riemann's xi function (so $H_0(z)=\frac18\xi(\frac12+\frac{iz}2)$); while simple, its zeros $z_k(t)$ move by $dz_k/dt=2\sum_{j\ne k}(z_k-z_j)^{-1}$, like a log gas. Problem: find a global functional of the zeros of $H_t$, insensitive to individual close pairs, that is monotone along the flow and whose positivity at $t=0$ excludes collisions at positive time. Biane's gap functional, the minimum over one gap of Laguerre's quantity, is local; the functional asked for should instead be one quantity for all the zeros.

Source: [`ch/cp6.tex`](https://github.com/VadimSokolov/computing_rh/blob/v1.1.0/ch/cp6.tex), label `cp:prob:gas`.

***

<small>License: Text, figures and PDF: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: [MIT](https://opensource.org/license/mit).</small>

<small>This page is generated from the book's repository ([main](https://github.com/VadimSokolov/computing_rh/tree/main)). The numbers and pages are those of the PDF, [version 1.1.0](https://github.com/VadimSokolov/computing_rh/releases/tag/v1.1.0).</small>

<!-- Generated by bookrepo.py from book.yaml and results.yaml. Edit those files, not this one. -->
