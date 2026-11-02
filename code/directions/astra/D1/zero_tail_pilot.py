"""Zero-side and elementary-test pilot, with no RH assumption for the tail bound."""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 40
L = mp.log(3)/2
a = mp.log(2)
c = mp.sqrt(2)*a


def linear_test():
    norm = 2*L**3/3
    def g(x):
        return (norm-L**2*x+x**3/6)/norm
    arch = -(mp.euler+mp.log(mp.pi)+mp.log(1-mp.exp(-4*L)))
    arch += mp.quad(lambda x: 2*(mp.exp(-2*x)-mp.exp(-x/2)*g(x))/(-mp.expm1(-2*x)), [0, 2*L])
    s = mp.quad(lambda x: x*mp.sinh(x/2), [-L,L])/mp.sqrt(norm)
    pole = -2*s*s
    prime = c*g(a)
    return {'arch_gamma': str(arch), 'pole': str(pole), 'A': str(arch+pole),
            'c_g_log2': str(prime), 'Q': str(arch+pole-prime)}


def refined_zero_tail(H, nmax):
    H = mp.mpf(H)
    # Twice integration by parts. An elementary polynomial Markov bound is used.
    A = mp.sqrt((nmax*nmax*mp.cosh(L)+nmax)/L)
    # The second integration bounds endpoint derivatives and ||phi_n''||_1.
    B2 = mp.mpf(0)
    for n in range(nmax):
        k = mp.sqrt((2*n+1)/(2*L))
        boundary_derivative = k*n*(n+1)/L
        second_variation = 2*L*k*max(n,1)**4/L**2 if n else mp.mpf(0)
        B2 += (mp.exp(L/2)*(boundary_derivative+second_variation))**2
    B = mp.sqrt(B2)
    # With z=t+i*delta, ordinary integration by parts leaves no exponential derivatives.
    # The bound for F(z) is (A+B/H)/t, uniformly over a unit coefficient vector.
    C2 = (A+B/H)**2
    # HSW implies |N(t)-M(t)| <= log(t)+10 for t>=100, with 7/8 in M.
    # Include the boundary term in Stieltjes integration by parts.
    tail_sum = (mp.log(H/(2*mp.pi))+1)/(2*mp.pi*H) + (2*mp.log(H)+mp.mpf('20.5'))/H**2
    return {'H': str(H), 'nmax': nmax, 'A': str(A), 'B': str(B),
            'S2_upper': str(tail_sum), 'zero_matrix_tail_upper': str(2*C2*tail_sum),
            'zero_count_asymptotic': str(H/(2*mp.pi)*mp.log(H/(2*mp.pi))-H/(2*mp.pi)+mp.mpf(7)/8)}


if __name__ == '__main__':
    result = {'status': 'computed, not certified', 'odd_linear_test': linear_test(),
              'finite_zero_tail': [refined_zero_tail(H, n) for H,n in [('3e12',64),('1e14',80)]]}
    # Smooth cosine trial: exact transform, all sampled zeros are on the line.
    k = mp.pi/(2*L)
    def F(z):
        return 2*k*mp.cos(L*z)/(mp.sqrt(L)*(k*k-z*z))
    zeros = [mp.im(mp.zetazero(j)) for j in range(1,161)]
    partials = {}
    for J in (40,80,160):
        H = zeros[J-1]+mp.mpf('0.01')
        S = 2*mp.fsum(F(t)**2 for t in zeros[:J])
        # |F(t+i delta)| <= 2*k*cosh(L/2)/(sqrt(L)*(t^2-k^2-1/4)).
        C = 2*k*mp.cosh(L/2)/(mp.sqrt(L)*(1-(k*k+mp.mpf('0.25'))/H**2))
        # N(t)<=t*log t gives sum gamma^-4 <=4*(log H/3+1/9)/H^3.
        err = 2*C*C*4*(mp.log(H)/3+mp.mpf(1)/9)/H**3
        partials[str(J)] = {'H':str(H), 'sum':str(S), 'absolute_tail_bound':str(err),
                            'lower':str(S-err), 'upper':str(S+err)}
    result['cosine_trial_zero_sum'] = partials
    def g(x):
        return ((2*L-x)*mp.cos(k*x)+mp.sin(k*x)/k)/(2*L)
    gamma = -(mp.euler+mp.log(mp.pi)+mp.log(1-mp.exp(-4*L)))
    gamma += mp.quad(lambda x: 2*(mp.exp(-2*x)-mp.exp(-x/2)*g(x))/(-mp.expm1(-2*x)), [0, 2*L])
    result['cosine_trial_prime_side'] = str(mp.re(2*F(mp.j/2)**2)+gamma-c*g(a))
    Path('zero_tail_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
