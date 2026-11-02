# Q3: reach of the step grid. Necessary condition for one quadruple (gamma,beta) to change the sign on window a, from S6;
# row 2 of the note's table at delta=0.005; delays of Table 2 by delta in the note's two conventions.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
def A_d(d): return math.log(1/d)/12+0.9375482543/(2*math.pi**2)
def lam_model(a,d): return (d/2)*math.pi**2*d*d*A_d(d)/a**2
def norm(a,b,d,g): return 4*math.sinh(a*b)/math.tanh(b*d/2)/g**2
print('A. quadruple norm bound 4 sinh(a beta) coth(beta delta/2)/gamma^2 against lam_min(T), a=12 and a=6, delta=1/40: no effect above gamma_0')
for a,lam in [(12.0,3.19e-7),(6.0,5.64e-7),(3.0,3.12e-6)]:
    print('  a=%g:'%a,' '.join('beta=%g: %.3g'%(b,math.sqrt(4*math.sinh(a*b)/math.tanh(b*0.0125)/lam)) for b in [0.5,0.3,0.2,0.1,0.05,0.02,0.01]))
print('B. smallest a with norm >= floor model, gamma=3e12 (necessary condition for a zero above the verified height to be seen); true rung for comparison a=log(gamma/4 pi)')
G=3e12
for d in [0.025,0.003125]:
    for b in [0.5,0.3,0.1,0.05,0.01]:
        lo,hi=1.0,1e5
        for _ in range(200):
            m=math.sqrt(lo*hi)
            x=m*b
            ok=(math.log(4)+ (x-math.log(2)) -math.log(math.tanh(b*d/2)) -2*math.log(G) >= math.log(lam_model(m,d)))
            if ok: hi=m
            else: lo=m
        print('  delta=%g beta=%g: a>=%.1f'%(d,b,hi))
print('  true rung reaches height 3e12 at a=%.2f, X=%.3g, margin digits %.3g'%(math.log(G/(4*math.pi)),G/(4*math.pi),G/math.log(10)))
ns,Ls=prime_powers(4000); ln=np.log(ns); u=Ls/np.sqrt(ns)
print('C. the note\'s rows on delta=0.005: last PD, first non PD')
for N in [1,7,31,1000]:
    i0=int(np.searchsorted(ns,N,side='right')); d=0.005
    c=cov_arch_zeta(d,int((ln[i0]+0.3)/d)+2); add_primes(c,d,ln[:i0],u[:i0]); ks,vs,f=levinson(c)
    print('  N=%d: last PD %.4f first non PD %.4f delay(last) %.4f delay(first) %.4f'%(N,f*d,(f+1)*d,f*d-ln[i0],(f+1)*d-ln[i0]))
