"""Assemble the report and plots from saved results and checked proofs."""
from pathlib import Path
import json
import os
import hashlib
import numpy as np

HERE=Path(__file__).resolve().parent
base=json.loads((HERE/'results/fits.json').read_text())
ext=json.loads((HERE/'results_extended/fits.json').read_text())
diag=json.loads((HERE/'results/diagnostics.json').read_text())
zero=json.loads((HERE/'results/zero_validation.json').read_text())

def table(headers,rows):
    values=[[str(x) for x in headers]]+[[str(x) for x in row] for row in rows]
    widths=[max(len(r[j]) for r in values) for j in range(len(headers))]
    return '~~~text\n'+'\n'.join('  '.join(x.ljust(widths[j]) for j,x in enumerate(row)) for row in values)+'\n~~~\n'

def fmt(x):
    return f'{x:.7g}'

merged={(r['family'],r['target'],r['N']):r for r in base['fits']}
for r in ext['fits']:
    key=(r['family'],r['target'],r['N'])
    if key not in merged or r['D192']<merged[key]['D192']:
        merged[key]=r
unpaired=[r for (fam,name,n),r in merged.items() if not all((fam,law,n) in merged for law in ['polya','dh'])]
merged={key:r for key,r in merged.items() if all((key[0],law,key[2]) in merged for law in ['polya','dh'])}
exactw={(r['family'],r['target'],r['N']):r for r in diag['w1']}

def matched_rows(fam):
    ns=sorted({n for f,t,n in merged if f==fam})
    rows=[]
    for n in ns:
        p,d=[merged[(fam,t,n)] for t in ['polya','dh']]
        wp=exactw.get((fam,'polya',n),{}).get('W1',p['W1_grid'])
        wd=exactw.get((fam,'dh',n),{}).get('W1',d['W1_grid'])
        rows.append([n,fmt(p['D192']),fmt(d['D192']),fmt(wp),fmt(wd)])
    return rows

prefix=r'''1. No ferromagnetic weak limit equal to the Pólya law was constructed; the unrestricted representation remains open.
2. A precise obstruction excludes bounded critical-window mean-field rules for signs and fixed splines, and a fixed number of critical Ising blocks.
3. Nonnegative, zero-field Ising fits become very close in a smooth law distance for both Pólya and the Davenport and Heilbronn control.
4. The control decreases at least as fast in the pooled regressions; this finite-range difference is not a sustained advantage for Pólya or evidence of a limit.
5. The explicit quartic limit suggested by homogeneous couplings is wrong: its zeros disagree at the first zero and beyond the fitting heights.

**Scope and conventions.** This report follows BRIEF.md, rh_directions.txt, astra/E/report.md and the spline addition to track_E.txt. No opus/ or fable/ content was read, no book file was edited, and no prohibited review script was used. All writes are in this output folder. The numerical results are exploratory floating-point calculations, not ball-arithmetic certificates. One independent proof check was used; its complete argument is incorporated below. The bounded laptop runs were expressly authorised in the task. No Hopper job was submitted.

Let
\[
\Phi(x)=\sum_{k\ge1}(4\pi^2k^4e^{9x/2}-6\pi k^2e^{5x/2})e^{-\pi k^2e^{2x}},
\quad p(x)=\Phi(x)/\xi(1/2).
\]
The full sum is even and positive, and its moment generating function is
\(L(z)=\xi(1/2+z)/\xi(1/2)\). Numerics evaluate the positive half-line and reflect the full kernel. The control is exactly the law used in E:
\[
\kappa=\frac{\sqrt{10-2\sqrt5}-2}{\sqrt5-1},\qquad
(a_1,a_2,a_3,a_4,a_5)=(1,\kappa,-\kappa,-1,0),
\]
with periodic coefficients, kernel
\[
K_D(x)=2e^{3x/2}\sum_{k\ge1}k a_k e^{-\pi k^2e^{2x}/5},
\]
and normalised transform
\[
L_D(z)=\Lambda_D(1/2+z)/\Lambda_D(1/2),\quad
\Lambda_D(s)=(5/\pi)^{(s+1)/2}\Gamma((s+1)/2)\sum_{k\ge1}a_k k^{-s}.
\]
Evenness and positivity of this kernel were proved in E. Thus the control is a probability law, not a signed inverse transform. Its off-line zeros are a known obstruction to any ferromagnetic limit. The coefficient formula and the zero near \(0.808517+85.699348i\) appear on the indexed first page of [Balanzario and Sánchez-Ortiz, Mathematics of Computation (2007)](https://www.ams.org/mcom/2007-76-260/S0025-5718-07-01999-0/S0025-5718-07-01999-0.pdf). Direct full-PDF access failed; the Firecrawl fallback had no credits. A separate post-fit calculation below verifies this numerical zero from the stated Dirichlet function.

The operational class here is the weak closure of nonnegative linear combinations of finite, zero-field Ising spins with nonnegative pair couplings. This avoids assuming that every definition of the Griffiths and Simon class, including extra tail requirements, is literally closed in all probability laws. The Pólya law itself has all Gaussian exponential moments. The imported classical input is the finite ferromagnetic Lee and Yang theorem and its product form, as stated in [Williams, section 7](https://www.statslab.cam.ac.uk/~grg/books/hammfest/22-dw.pdf) and [Newman (1975)](https://fig.if.usp.br/~marchett/rg/inequalities-ising-models-field-theories-which-obey-lee-yang-theorem_newman.pdf). The limit obstructions and the tightness argument below are proved explicitly.

'''
proof=(HERE/'obstruction_check.md').read_text()
proof=proof.replace('**Independent mathematical check of the proposed obstruction.**','**Precise statements and full proofs.**',1)
proof=proof.replace('The \\(n=1\\) term in the theta kernel dominates as \\(x\\to+\\infty\\), giving',r'''For completeness, put \(y=e^{2x}\). The leading term is
\(4\pi^2e^{9x/2-\pi y}(1-3/(2\pi y))\). For \(x\ge0\), its parenthesis is bounded away from zero. The ratio of the remaining positive terms to it is bounded by
\(C\sum_{k\ge2}k^4e^{-\pi(k^2-1)y}\).
The ratio of consecutive summands is at most \((3/2)^4e^{-5\pi y}<1/2\), so the sum is at most \(32e^{-3\pi y}\) and tends to zero. Taking logarithms gives''')

method=r'''
**What the Euler form as splines does, and what coupling changes.** For a fixed spline \(V=\sum_{j=1}^m a_j U_j\), with independent \(U_j\) uniform on \([-1,1]\),
\[
\mathbb E e^{zV}=\prod_{j=1}^m\frac{\sinh(a_jz)}{a_jz}
=\prod_{j=1}^m\prod_{k\ge1}\left(1+\frac{a_j^2z^2}{\pi^2k^2}\right).
\]
These are integer frequencies, not primes. E's independent-sum obstruction applies to this product and its weak closure. The positive quadratic tilt in the proof genuinely escapes factorisation into independent transforms: it introduces the displayed nonnegative couplings between binary Ising approximations of the uniform coordinates. Thus the present obstruction concerns an interacting family that E's independent-sum determinant does not cover.

The conclusion is nevertheless negative for this specified mean-field rule. Its entropy expansion leaves a quartic limit, even when the single-site variable is a fixed spline. Maxima, a random mixture, logarithms, and size biases do not preserve the linear magnetisation construction merely because they produce the correct probability law. E's Brownian-range representation remains an exact nonlinear representation, with no new Griffiths and Simon preservation theorem supplied here. A growing number of spline constituents, an increasing number of critical blocks, or a different interaction geometry is outside the proved obstruction. Being outside it is not a construction.

**Numerical experiment, definition of the distance and the fits.** Every fitted system is
\[
P_J(\sigma)=Z_J^{-1}\exp\left(\sum_{i<j}J_{ij}\sigma_i\sigma_j\right),
\quad J_{ij}\ge0,\quad \sigma_i\in\{-1,1\}.
\]
All fields are exactly zero. Put \(S=\sum_i\sigma_i\) and choose the common spin weight
\[
a(J)=\sqrt{\operatorname{Var}(X_{\rm target})/\mathbb E_J S^2}.
\]
Thus the physical output \(a(J)S\) matches the target variance. Distances compare \(Y=S/\sqrt{\mathbb E_JS^2}\) with the variance-one target \(X/\sigma_X\), so the different physical scales of the two targets cannot create the comparison.

The fitting distance is the Gaussian Fourier distance
\[
D(\mu,\nu)^2=\sqrt{\frac2\pi}\int_0^\infty
|\widehat\mu(t)-\widehat\nu(t)|^2e^{-t^2/2}\,dt.
\]
This is a distance between laws: zero integral gives equal characteristic functions by continuity and positivity of the weight, then equality of probability laws. Numerically, the integral is cut at 12, with 64 Gauss-Legendre fitting nodes and 192 evaluation nodes. Since each characteristic function has modulus at most one, the norm of the omitted part is at most \(\sqrt{4\operatorname{erfc}(12/\sqrt2)}\). It is not a moment-only objective, although it strongly favours the central part of the law and low frequencies.

I also report the standardised Wasserstein distance \(W_1=\int|F_\mu-F_\nu|\,dx\), plus relative errors in moments 4 through 12 held out of the objective. To obtain physical \(W_1\), multiply by the relevant target standard deviation. Total variation is unsuitable here: every finite spin sum is discrete and each target has a density, so their total variation distance is exactly one, despite possible weak convergence.

There are three fitted families. For Curie-Weiss, \(J_{ij}=\beta/N\), so
\[
P(S=2k-N)\ \propto\ {N\choose k}\exp[\beta(2k-N)^2/(2N)],\quad\beta\ge0.
\]
For two blocks of size \(n=N/2\), the three parameters give
\[
J_{ij}=\theta_1/n\ \hbox{within block 1},\quad
J_{ij}=\theta_2/n\ \hbox{within block 2},\quad
J_{ij}=\theta_{12}/n\ \hbox{between blocks},
\]
all nonnegative. The exact probability of the two sums has binomial multiplicities and exponent
\(\theta_1S_1^2/(2n)+\theta_2S_2^2/(2n)+\theta_{12}S_1S_2/n\).
For the dense family, every one of the \(N(N-1)/2\) couplings is fitted separately. Enumeration fixes the first spin to +1 and uses global sign symmetry for every even observable. There is no Monte Carlo error.

The complete-graph search profiles a 61-point grid on \([0,3]\) and refines the best interval. This is a numerical minimum, not a global-minimum certificate. The other fits use bounded L-BFGS-B, four starts initially, twelve starts for the larger two-block comparison, analytic gradients including differentiation of the variance normalisation, at most 700 iterations and 1500 function evaluations per start. The fitted block and dense parameters lie in \([0,4]\); these boxes are search restrictions. Seeds and all termination messages are saved. No distance below is advertised as the global optimum over all ferromagnets.

'''

numeric='**Best distances found as N grows.** P means Pólya; D-H means the Davenport and Heilbronn control. The W1 entries for larger complete graphs and two-block systems use the refined CDF calculation. Smaller entries use the original fine-grid approximation.\n\n'
if unpaired:
    numeric+='The extended run reached the 600-second cap while fitting the control at N=512. Its completed Pólya N=512 record is saved in the raw output, but the paired tables and rate regressions stop at N=256 for this family. The larger paired comparison is a Hopper proposal.\n\n'
for fam,label in [('cw','Homogeneous complete graph'),('block','Two equal blocks'),('dense','Unrestricted pair couplings')]:
    numeric+=f'**{label}.**\n\n'+table(['N','D(P)','D(D-H)','W1(P)','W1(D-H)'],matched_rows(fam))+'\n'

numeric+='**The fitted couplings and the possible rules.** These parameters specify every edge and the output normalisation; all individual fits and all starts are in [fits.json](results/fits.json) and [the extended fits](results_extended/fits.json).\n\n'
rows=[]
for r in base['fits']:
    if r['family']=='cw' and r['N'] in [16,64,256,1024,4096]:
        rows.append([r['target'],r['N'],fmt(r['theta'][0]),fmt(r['b_window']),fmt(r['c_window'])])
numeric+=table(['law','N','beta','sqrt(N)*(beta-1)','N^(3/4)*spin_weight'],rows)
numeric+='\nThe homogeneous pattern is a bounded critical window. A fixed, explicit rule suggested by it is\n\n'
numeric+=r'''\[
J_{ij}^{(N)}=\frac{\max(0,1+b/\sqrt N)}{N},\qquad
X_N=cN^{-3/4}\sum_i\sigma_i.
\]
The limit of this rule is proved above:
\[
q_{b,c}(x)=\frac1{cC_b}\exp\left[-\frac{x^4}{12c^4}+\frac{b x^2}{2c^2}\right].
\]
The following values were fitted directly in the limiting quartic family, using the same law distance and no zeros. They are finite-precision fitted constants, not exact identities for the target. For the held-out zero calculation they are frozen as the exact displayed decimal parameters.

'''
numeric+=table(['law','b','c','D to target'],[[q['target'],repr(q['b']),repr(q['c']),fmt(q['D192'])] for q in base['quartic']])
numeric+='\nBoth limiting discrepancies are nonzero. The theorem proves that no other choice of these two constants can give Pólya either; the numerical minimum is not needed for that conclusion.\n\n'
rows=[]
for r in diag['block_starts']:
    lo,hi=sorted(r['theta'][:2])
    n=r['N']/2
    rows.append([r['target'],r['N'],fmt(lo),fmt(hi),fmt(r['theta'][2]),fmt((lo-1)*np.sqrt(n)),fmt((hi-1)*np.sqrt(n)),fmt(r['theta'][2]*np.sqrt(n))])
numeric+=table(['law','N','theta_low','theta_high','theta_cross','B_low','B_high','B_cross'],rows)
numeric+='\nBlock labels were exchanged where needed for comparison. The finite fitted B values motivate the fixed-block obstruction, but do not prove that an infinite fitted sequence has bounded B. If a proposed two-block rule does keep these three B parameters bounded, it is excluded by the theorem, whatever constants the fit suggests. No stable prime-based rule was identified.\n\n'
rows=[]
for name,r in diag['dense_couplings'].items():
    rows.append([name,r['N'],fmt(r['minimum']),fmt(r['maximum']),r['positive_count'],r['total_edges'],', '.join(fmt(x) for x in r['sorted_degrees'])])
numeric+=table(['law','N','min J','max J','J>1e-8','edges','sorted weighted degrees'],rows)
numeric+='\nThe complete fitted N=12 coupling matrices are in [dense_couplings.json](results/dense_couplings.json); the upper-triangle arrays for every smaller N are in fits.json. These heterogeneous finite solutions do not determine a unique graph: relabelling spins preserves the fitted law, and other nonidentifiabilities remain. They supply no extrapolation rule beyond the measured sizes.\n\n'
edge_rows=[]
for name,r in diag['dense_couplings'].items():
    mat=np.asarray(r['matrix'])
    for i in range(r['N']):
        for j in range(i+1,r['N']):
            if mat[i,j]>1e-8:
                edge_rows.append([name,i+1,j+1,f'{mat[i,j]:.10g}'])
numeric+='All nonzero N=12 edges, using the saved spin labels, are:\n\n'+table(['law','spin i','spin j','J_ij'],edge_rows)
numeric+=r'''
Here a concrete finite pattern is visible. Pólya's N=12 fit has a path of three interacting spins and nine isolated spins. The control's fit has a complete component of five spins and seven isolated spins. Their edge strengths are given above; every other edge is zero to the stated threshold. This is a sparse-component observation at small N, not evidence of an arithmetic graph rule.

There is a further elementary obstruction to its simplest extrapolation. Suppose only m_N=o(sqrt(N)) vertices participate in any nonzero coupling, while the other N-m_N spins are isolated, and use equal spin weights. Write S_N=R_N+B_N, where R_N is a sum of N-m_N independent fair signs, independent of B_N, and |B_N|<=m_N. Then, for every fixed real t,
\[
\mathbb E e^{itS_N/\sqrt N}
=\cos(t/\sqrt N)^{N-m_N}\mathbb E e^{itB_N/\sqrt N}
\longrightarrow e^{-t^2/2}.
\]
Indeed log cos(t/sqrt(N))=-t^2/(2N)+O(N^(-2)), while the absolute difference of the second factor from one is at most |t| m_N/sqrt(N). The continuity theorem gives S_N/sqrt(N) converging to a standard normal. The same tightness and subsequence scaling argument used above shows that every nondegenerate limit under an arbitrary common positive normalisation is a dilated Gaussian. Its quadratic log density cannot equal the Pólya log density in (5). Thus retaining a fixed correlated component and adding isolated spins cannot turn either observed small-N pattern into the desired limit. This result does not apply when the correlated part grows sufficiently rapidly or when spin weights become nonuniform.

'''

numeric+='**Rates and uncertainty.** The table gives the descriptive exponent alpha in D proportional to N^(-alpha), obtained by least-squares regression of log distance on log N over the stated sizes. It is not a convergence theorem or a statistical confidence interval.\n\n'
rows=[]
for fam in ['cw','block','dense']:
    p=diag['rates'][fam+'_polya']; d=diag['rates'][fam+'_dh']
    rows.append([fam,','.join(map(str,p['N'])),fmt(p['alpha_D']),fmt(d['alpha_D']),fmt(d['alpha_D']-p['alpha_D'])])
numeric+=table(['family','N in regression','alpha(P)','alpha(D-H)','alpha(D-H)-alpha(P)'],rows)
rows=[]
for fam in ['cw','block']:
    for name in ['polya','dh']:
        selected=sorted([r for r in diag['w1'] if r['family']==fam and r['target']==name],key=lambda r:r['N'])
        alpha=-np.polyfit(np.log([r['N'] for r in selected]),np.log([r['W1'] for r in selected]),1)[0]
        rows.append([fam,name,fmt(alpha)])
numeric+='\nFor the refined Wasserstein calculation, the corresponding descriptive exponents are:\n\n'+table(['family','law','alpha(W1)'],rows)
rows=[]
for i,first in enumerate([1729,1733,1737]):
    p=diag['rates']['block_polya']['optimizer_batches'][i]['alpha']
    d=diag['rates']['block_dh']['optimizer_batches'][i]['alpha']
    rows.append([f'{first} to {first+3}',fmt(p),fmt(d),fmt(d-p)])
numeric+='\nThree disjoint batches of four starts give the following two-block rate sensitivity. These are algorithmic replicates on deterministic laws, not independent observations of the target distribution.\n\n'+table(['seed batch','alpha(P)','alpha(D-H)','difference'],rows)
rows=[]
for p,d in zip(diag['rates']['block_polya']['adjacent_alpha'],diag['rates']['block_dh']['adjacent_alpha']):
    rows.append([str(p['N']),fmt(p['alpha']),fmt(d['alpha']),fmt(d['alpha']-p['alpha'])])
numeric+='\nThe sensitivity to the fitted size interval is much larger:\n\n'+table(['adjacent N','alpha(P)','alpha(D-H)','difference'],rows)
rows=[]
for p,d in zip(diag['rates']['dense_polya']['per_start_alpha'],diag['rates']['dense_dh']['per_start_alpha']):
    rows.append([p['seed'],fmt(p['alpha']),fmt(d['alpha']),fmt(d['alpha']-p['alpha'])])
numeric+='\nFor the unrestricted dense fits, using one initialisation at a time gives this optimisation sensitivity over N=4,6,8,10,12:\n\n'+table(['seed','alpha(P)','alpha(D-H)','difference'],rows)
rows=[]
for r in diag['block_starts']:
    rows.append([r['target'],r['N'],fmt(r['Dmin']),fmt(r['Dmedian']),fmt(r['Dmax']),f"{r['successes']}/{r['total']}"])
numeric+='\nThe full twelve-start spread is:\n\n'+table(['law','N','best D','median D','worst D','reported success'],rows)
numeric+=r'''
The measured pooled slopes favour the control in all three families. This is a finite-range observation, with the explicit qualifications in the tables. In the dense small-N family the control also fits better than Pólya at every displayed size, but the slope difference is sensitive to initialisation. In the complete-graph family Pólya has smaller errors over most finite sizes, while both explicit quartic candidates have almost the same positive discrepancy. In the two-block family, the best-of-four batches reproduce the small pooled exponent difference closely; it is therefore not explained by the displayed numerical quadrature error. But changing from one adjacent size interval to the next reverses the comparison, and the control error is nonmonotone. Local minima also create a wide best-to-worst spread, so the best values remain upper bounds on unknown optima. The Wasserstein decreases largely reflect a finer spin lattice and occur for both laws. There is no sustained Pólya advantage, and these family-dependent rates do not distinguish membership in the closed class or exploit the Euler product.

There is no sampling standard error to report: target quadrature and state enumeration are deterministic. Quadrature variation, truncation bounds, optimisation spread and the choice of the N window are separate uncertainties. A regression standard error would not quantify any of the unresolved infinite-size questions. Small D can coexist with much larger W1 because the Gaussian Fourier weight strongly smooths away lattice structure and high-frequency errors.

'''
numeric+='**Held-out moments and numerical checks.** Relative errors below are model moment divided by target moment, minus one, for the largest fitted system of each family. Variance is matched by construction.\n\n'
rows=[]
for fam in ['cw','block','dense']:
    for name in ['polya','dh']:
        rr=max((r for r in merged.values() if r['family']==fam and r['target']==name),key=lambda r:r['N'])
        rows.append([fam,name,rr['N']]+[fmt(x) for x in rr['moment_relative_errors'][1:]])
numeric+=table(['family','law','N','mu4','mu6','mu8','mu10','mu12'],rows)
numeric+='\nThe largest measured changes/check residuals are:\n\n'
checks={**base['validation'],**diag['validation']}
checks['W1_mesh_change_max']=max(r['mesh_change'] for r in diag['w1'])
numeric+=table(['check','value'],[[k,fmt(v) if isinstance(v,(float,int)) else v] for k,v in checks.items()])
numeric+=r'''
The target integrals use 512 Gauss-Legendre nodes on [0,4], compared against 1024, with 24 Pólya terms and 48 control terms. Applying the explicit geometric theta-tail and spatial-tail bounds of E with these cutoffs makes those two truncations smaller than 10^(-600) for the moments used here. The visible numerical uncertainty comes from floating-point quadrature and optimisation, not those cutoffs. Agreement under doubling quadrature is an error diagnostic, not a rigorous enclosure. The control positivity and kernel identity do not depend on the numerical quadrature.

Block multiplicities were independently compared with direct enumeration of all configurations modulo global sign at N=12, including variance normalisation. The analytic gradient was also checked by centred finite differences. All fitted probabilities are normalised from log weights, and all fields and coupling constraints are enforced in the code.

'''

zeros=r'''**The candidate fails the zero extrapolation.** No zero of zeta or of the control was available to the fitting script. Only after freezing its output did validate_results.py read zero data. It used 80-decimal-digit arithmetic and 256 versus 384 Gauss-Legendre nodes for the quartic transform, with a sign scan in physical height and local root refinement. An initial 192-node scan was insufficient for the highest control-candidate roots, so it was discarded for the final root validation. These are numerical zeros of the specified candidate, not certified root counts or a claim that its parameters are known to 80 digits.

'''
zeros+=f"The Pólya fit integrated standardised frequencies up to 12, corresponding to physical height {zero['polya']['fit_physical_height_cutoff']:.8f}. In particular, rows 12, 15 and 20 below test beyond that range. The same-index comparison is already wrong at index 1.\n\n"
rows=[]
for r in zero['polya']['comparisons']:
    rows.append([r['index'],f"{float(r['known_zeta_ordinate']):.9f}",f"{float(r['candidate_ordinate']):.9f}",f"{float(r['signed_difference']):+.9f}"])
zeros+=table(['index','zeta ordinate','quartic candidate','candidate minus zeta'],rows)
qchange=max(float(r['quadrature_root_change']) for name in ['polya','dh'] for r in zero[name]['roots'])
zeros+=f'\nAcross the scanned roots of both quartic candidates, the largest 256-to-384-node root change was {qchange:.3g}. This is accuracy for frozen decimal parameters, not uncertainty in fitting a proposed exact rule. The mismatches in the table are much larger.\n\n'
zeros+=r'''The control has a still sharper failure: every quartic candidate is a genuine Lee and Yang limit and therefore has no off-axis zero. The completed Davenport and Heilbronn function has an off-axis zero. The independent post-fit Hurwitz-zeta calculation gives
\[
s\simeq0.8085171824566374+85.69934848537759i.
\]
It evaluates \(5^{-s}\sum_{r=1}^4a_r\zeta(s,r/5)\), without truncating a nonconvergent Dirichlet series. The gamma factor is nonzero there, so this is also a zero of the completed function. The calculation corroborates the published zero and is not a new interval certificate. Its residual and all the scanned real candidate roots are saved in [zero_validation.json](results/zero_validation.json).

At a simple zero, a small transform perturbation displaces the zero approximately by the error divided by the derivative. The derivative can be extremely small at high heights. The completed zeta transform includes the factor exp(-pi |t|/4), up to powers and the zeta factor. Consequently an accurate central-law fit cannot control all zero heights. The explicit discrepancies above are a direct falsification of this candidate, independently of the tail proof.

**What a good fit does and does not show.** Every fitted finite system is a legitimate ferromagnet. Its D value shows only that one probability law approximates selected smooth observables of the target. Finitely many N, even with tiny errors and every low moment close, do not establish a sequence with an error tending to zero. Weak closure works in the important direction: an actual ferromagnetic sequence converging to Pólya would imply RH by the proof above. Hence any such apparent success must be presumed wrong until membership, infinite-size convergence and all analytic passages are verified.

There is also a useful compactness statement. Variance-one ferromagnetic laws have uniformly Gaussian tail bounds, as proved above. Any sequence has a weakly convergent subsequence, whose limit is again in the operational closed class and still has variance one. If their Gaussian Fourier distances to a fixed target went to zero, dominated convergence on the frequency integral would identify every subsequential limit with that target. Thus the known non-Lee-and-Yang control has a strictly positive infimum distance from this variance-one class. The experiment does not compute that infimum or a lower bound for it. Very small fitted control distances are consistent with that positive gap.

**The four ground rules and the stopping point.** The obstruction is infinite: it classifies limits of an explicitly stated sequence family, including arbitrary positive normalisations and fixed-block projections. The numerical tables alone are finite. No construction here reaches Pólya, and no claim about an unrefereed 7/8 half-plane result is used.

The spline Euler products are not the prime Euler product. The theta formula supplies zero-free input data to the fits, but the fitting rule does not extract multiplicativity or any other condition that excludes the control. Indeed a_6=1 while a_2 a_3=-kappa^2. A finite-prime calculation without a theorem controlling the remaining arithmetic cannot establish the required infinite limit. One can avoid enumerating all primes; one cannot replace all infinite arithmetic control by a finite fit.

Gaussian damping preserves the homogeneous quartic family: multiplying q_(b,c)(x) by exp(-epsilon x^2) replaces b by b-2 epsilon c^2. By contrast, the Pólya damping at any negative de Bruijn and Newman time has nonreal Fourier zeros, using [Rodgers and Tao](https://arxiv.org/abs/1801.05914). Thus a hypothesis stable under this damping cannot simply identify Pólya as one of these quartic Lee and Yang laws. The tail obstruction makes that failure explicit.

Finally, the positive quadratic interaction in the spline construction adds couplings. It proves admissibility of the approximating model, not removal of the coupling needed in the authors' two-spin route. Neither fitting nonnegative J nor letting a fixed mean-field parameter approach its critical value supplies that missing removal theorem.

What remains open is the unrestricted coupled-spin representation. The proof does not rule out increasing block complexity, changing base laws, different critical scaling, a lattice geometry, or a prime-driven rule with a new convergence theorem. No such rule emerged from these fits. Claiming a general obstruction from the quartic-family result would be incorrect.

**Reproducibility and the next computation.** [fit_spins.py](fit_spins.py) contains all target and spin-law code. [validate_results.py](validate_results.py) performs the independent law and zero checks. [empirical.md](empirical.md) records planned and completed runs; [numbers.txt](results/numbers.txt), [diagnostics.json](results/diagnostics.json), [the extended ledger](results_extended/numbers.txt), and the JSON parameter arrays contain the numerical provenance. [run_all.sh](run_all.sh) regenerates the final fits, checks, tables, plot and report, with each Python process bounded by timeout. The initial shorter-frequency trial is preserved separately and is not used in the final tables.

The larger dense search is [hopper_fit.py](hopper_fit.py) with [hopper_proposal.sbatch](hopper_proposal.sbatch), not submitted. It proposes 96 CPU array tasks: N=14,16,18, two laws, sixteen paired initialisations per size and law, at most eight concurrent tasks, one CPU and 4 GB per task, two hours per task. At N=18 there are 131072 configurations after global-sign reduction and 153 edge parameters; the feature matrix alone is 160432128 bytes. The script uses exact enumeration, analytic gradients, 96 fitting frequencies and 192 validation frequencies. [hopper_blocks.sbatch](hopper_blocks.sbatch) separately proposes the paired N=512,1024 block continuation, twelve starts, one CPU, 4 GB and a four-hour wall cap. Neither proposal has been submitted. The scratch directory and hopper_logs/ must be created and the code copied before a user elects to submit them.

This proposed search would test whether the small dense-graph fits develop a reproducible coupling pattern and whether their rate comparison survives a wider size range and more starts. It would still establish neither convergence nor RH. The next mathematical requirement is an explicit nonnegative rule outside the proved families, a uniform tail argument, and convergence to the exact theta law. Testing the same complete rule on the control, and testing its zeros beyond every fitting height, should precede any success claim.

'''

runtime=table(['run','seconds'],[['final-grid fits',fmt(base['runtime'])],['extended block attempt','600, timeout' if unpaired else fmt(ext['runtime'])],['law diagnostics',fmt(diag['runtime_seconds'])],['post-fit zeros',fmt(zero['runtime_seconds'])]])
report=prefix+proof+'\n'+method+numeric+zeros+'Recorded bounded-run durations:\n\n'+runtime
assert '\u2014' not in report and '---' not in report
(HERE/'report.md').write_text(report)
(HERE/'results/dense_couplings.json').write_text(json.dumps(diag['dense_couplings'],indent=2))

# A static research figure accompanies the numeric tables.
os.environ['MPLCONFIGDIR']=str(HERE/'fig/mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
(HERE/'fig').mkdir(exist_ok=True)
fig,axes=plt.subplots(1,3,figsize=(12,3.6))
for ax,fam,label in zip(axes,['cw','block','dense'],['Complete graph','Two blocks','Unrestricted edges']):
    for target,color,label2 in [('polya','#1b6ca8','Polya'),('dh','#c25932','Davenport and Heilbronn')]:
        rr=sorted([r for r in merged.values() if r['family']==fam and r['target']==target],key=lambda r:r['N'])
        ax.loglog([r['N'] for r in rr],[r['D192'] for r in rr],'o-',color=color,label=label2)
        if fam=='cw':
            qq=next(q for q in base['quartic'] if q['target']==target)
            ax.axhline(qq['D192'],color=color,linestyle=':',alpha=.6)
    ax.set_title(label)
    ax.set_xlabel('Number of spins N')
    ax.grid(True,which='both',alpha=.2)
axes[0].set_ylabel('Gaussian Fourier distance, variance one')
axes[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(HERE/'fig/distances.png',dpi=180)
plt.close(fig)
report=report.replace('**The fitted couplings and the possible rules.**','![Distance versus number of spins for both laws](fig/distances.png)\n\nThe dotted lines in the complete-graph panel are the fitted positive quartic-limit discrepancies. All plotted points are achieved distances, not certified optima.\n\n**The fitted couplings and the possible rules.**')
(HERE/'report.md').write_text(report)

# Preserve all number sources and a final manifest, without overwriting per-run ledgers.
ledger=[]
def walk(value,path):
    if isinstance(value,dict):
        for k,v in value.items(): walk(v,path+'.'+k)
    elif isinstance(value,list):
        for i,v in enumerate(value): walk(v,path+f'[{i}]')
    elif isinstance(value,(float,int)) and not isinstance(value,bool):
        ledger.append(f'{value}  {path}')
for name,data in [('results/fits.json',base),('results_extended/fits.json',ext),('results/diagnostics.json',diag),('results/zero_validation.json',zero)]:
    walk(data,name)
(HERE/'results/all_numbers.txt').write_text('\n'.join(ledger)+'\n')
import scipy
import mpmath
(HERE/'requirements.txt').write_text(f'numpy=={np.__version__}\nscipy=={scipy.__version__}\nmpmath=={mpmath.__version__}\nmatplotlib=={matplotlib.__version__}\n')
files=['fit_spins.py','validate_results.py','build_report.py','obstruction_check.md','report.md','results/fits.json','results_extended/fits.json','results/diagnostics.json','results/zero_validation.json']
manifest={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}
(HERE/'results/manifest.json').write_text(json.dumps(manifest,indent=2))
print('report.md written,',len(report.split()),'words')
print(runtime)
print('rate summary',json.dumps(diag['rates'],indent=2))
print('maximum refined W1 mesh change',checks['W1_mesh_change_max'])
print('maximum zero quadrature root change',qchange)
