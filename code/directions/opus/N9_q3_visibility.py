# Q3: planted quadruple 1/2 +- beta +- i gamma on the note's grid: which heights does a positive result see?
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
t0=time.time()
ns,Ls=prime_powers(250000); ln=np.log(ns); u=Ls/np.sqrt(ns)
rng=np.random.default_rng(9)
def plant(d,m,gam,beta):
    z=gam+1j*beta
    w=(1-np.cos(z*d))/z**2
    return 4*np.real(w*np.cos(np.arange(m)*z*d))
def scan(a,d,betas,cent,nj):
    m=int(round(a/d)); c=cov_zeta(d,m+1,ln,u)[:m]; T=toeplitz(c)
    lam0=eigvalsh(T,subset_by_index=[0,0])[0]
    print(' a=%g delta=%g m=%d Nyquist %.1f lam_min(T)=%.3e'%(a,d,m,math.pi/d,lam0))
    for b in betas:
        hit=[]
        for g0 in cent:
            cnt=0
            for j in range(nj):
                g=g0*(1+0.02*rng.random())
                D=plant(d,m,g,b)
                l=eigvalsh(T+toeplitz(D),subset_by_index=[0,0])[0]
                cnt+= l<0
            hit.append(cnt/nj)
        hit=np.array(hit)
        g50=cent[np.nonzero(hit>=0.5)[0][-1]] if (hit>=0.5).any() else float('nan')
        gany=cent[np.nonzero(hit>0)[0][-1]] if (hit>0).any() else float('nan')
        low=cent[np.nonzero(hit>=0.5)[0][0]] if (hit>=0.5).any() else float('nan')
        bound=math.sqrt(math.sinh(a*b)/math.tanh(b*d/2)*2/lam0)
        print('  beta=%.3f: detected in >=50%% of positions up to gamma=%.3g (lowest such %.3g), in any position up to %.3g; norm bound: no effect above %.3g'%(b,g50,low,gany,bound))
cent=np.exp(np.linspace(math.log(15),math.log(3e6),70))
scan(12,0.025,[0.5,0.3,0.2,0.1,0.05,0.02,0.01],cent,5)
print('time %.1f'%(time.time()-t0))
scan(6,0.025,[0.5,0.3,0.2,0.1,0.05,0.02,0.01],cent,5)
print('time %.1f'%(time.time()-t0))
scan(3,0.025,[0.5,0.2,0.1,0.02],cent,5)
cent2=np.exp(np.linspace(math.log(15),math.log(3e6),36))
scan(12,0.0125,[0.5,0.1,0.02],cent2,3)
print('total %.1f'%(time.time()-t0))
