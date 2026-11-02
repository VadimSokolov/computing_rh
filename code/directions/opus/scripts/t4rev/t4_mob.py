import sys; sys.path.insert(0,'/tmp/t4rev')
from lib import *
import mpmath as mp
mp.mp.dps=20
g=zeros8000(); T=100.; gT=g[g<=T]; g1=g[0]
S=float((gT**-2.).sum()); S4=float((gT**-4.).sum())
# var X = (log xi)''(1/2)
def xi(s): return s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
varX=float(mp.diff(lambda s:mp.log(xi(s)),mp.mpf(.5),2))
sig2=varX-2*S; print('S_T',S,'S4',S4,'varX',varX,'sigma2',sig2,'sum_gamma>0 g^-2',varX/2, 'n',len(gT))
# tails: sum_{gamma>T}gamma^-4, estimated with zeros+integral
gg=g[g>T]; tail_pos=(gg**-4.).sum()+ (lambda T2: (np.log(T2/2/np.pi)+1/3)/(2*np.pi*3*T2**3))(g[-1])
print('sum_{gamma>100}gamma^-4 (true)',tail_pos,'both signs',2*tail_pos,'  note bound 4.87e-7')
# phi_X by Hadamard product is what we want independent of spins: use Riemann kernel density + mp phi on grid
dt=0.01; tg=np.arange(0,120+dt/2,dt)
# phi via exact: exp(sum log(1-t^2/g^2)) sign-tracked, tail via series with mean density
def phi(t):
    # product over 8000 zeros and tail expo correction
    p=np.ones_like(t)
    for x in g: p*=1-(t/x)**2
    s2=varX/2-(g**-2.).sum(); s4=tail4=(np.log(g[-1]/2/np.pi)+1/3)/(2*np.pi*3*g[-1]**3)
    return p*np.exp(-s2*t*t-s4*t**4/2)
ph=phi(tg)
print('check phi(20)',ph[2000],float(mp.re(xi(mp.mpf(.5)+20j)/xi(mp.mpf(.5)))))
def Phi(u):
    u=abs(u); k=np.arange(1,14)[:,None]
    return (2*(2*np.pi**2*k**4*np.exp(4.5*u)-3*np.pi*k**2*np.exp(2.5*u))*np.exp(-np.pi*k**2*np.exp(2*u))).sum(axis=0)
xs=np.linspace(0,3,301); xih=float(xi(mp.mpf(.5))); fX=Phi(xs)/xih
print('fX(0)',fX[0],'inversion check',np.trapezoid(ph,tg)/np.pi)
def mob(N):
    m=np.ones(N+1,int); s=np.ones(N+1,bool)
    for p in range(2,N+1):
        if s[p]:
            s[2*p::p]=False; m[p::p]*=-1; m[p*p::p*p]=0
    return m
mu=mob(2000); mK=6/np.pi**2-np.cumsum(np.where(np.arange(2001)>0,mu/np.maximum(np.arange(2001),1)**2.,0))
z4=np.pi**4/90
Sig4_used=4.873e-7
rows=[]
Ks=list(range(1,61))+[70,78,80,90,100,110,120,150,200]
for K in Ks:
    U=min((K+1)*g1/2,T/2); t=tg[tg<=U]; ps=np.exp(-sig2*t*t/2)
    lp=np.zeros_like(t);sgn=np.ones_like(t)
    for k in range(1,K+1):
        if mu[k]==0: continue
        for x in gT:
            sc=np.sinc(t/(k*x))
            if mu[k]>0: ps=ps*sc
            else:
                ps=ps/np.where(np.abs(sc)<1e-300,1e-300,sc)
    w=np.full_like(t,dt); w[0]=w[-1]=dt/2
    fK=(np.cos(np.outer(xs,t))*(w*ps)).sum(axis=1)/np.pi
    act=np.abs(fK-fX).max()
    A=np.pi**2*t*t/6*abs(mK[K])*S+2*z4/(9*K**3)*t**4*S4; B=2/3*t**4*Sig4_used
    aph=np.abs(ph[tg<=U]); band=np.trapezoid(aph*(np.exp(A+B)-1),t)/np.pi
    tail=np.trapezoid(np.abs(ph[tg>=U]),tg[tg>=U])/np.pi
    rows.append((K,U,mK[K],act,band+tail,(band+tail)/act))
    if K in (1,2,3,5,10,30): print(f'K={K} U={U:.2f} m={mK[K]:+.6f} act={act:.5g} refined={band+tail:.5g} ratio={(band+tail)/act:.3f}')
r=np.array(rows); bad=r[r[:,5]>2]
print('ratio>2 at K:',[(int(a),round(float(c),6),round(float(b),2)) for a,c,b in zip(bad[:,0],bad[:,2],bad[:,5])][:40]); print('max ratio',r[:,5].max(),'at K',r[r[:,5].argmax(),0],'min ratio',r[:,5].min(),'at K',r[r[:,5].argmin(),0])
print('smallest act',r[:,3].min(),'last rows',[tuple(np.round(x,6)) for x in r[-5:]])
print('m30',mK[30])
