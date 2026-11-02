import sys; sys.path.insert(0,'/tmp/t4rev')
from lib import *
import mpmath as mp
from scipy.integrate import quad
g=zeros8000(); T=100.; gT=g[g<=T]; g1=g[0]; S=(gT**-2.).sum(); S4=(gT**-4.).sum(); varX=0.04620998623083794; sig2=varX-2*S
for nm,(a,b,c) in {'script .11,.29,2.29':(.11,.29,2.29),'Trudgian14 .112,.278,2.51 (memory)':(.112,.278,2.51),'HSW22 .1095,.2042,3.0305':(.1095,.2042,3.0305)}.items():
    Nup=lambda x:x/(2*np.pi)*np.log(x/(2*np.pi*np.e))+7/8+a*np.log(x)+b*np.log(np.log(x))+c+1/(50*x)
    print(nm,'4 int N^+ x^-5 =',quad(lambda x:4*Nup(x)*x**-5,100,np.inf)[0])
dt=.01; tg=np.arange(0,120+dt/2,dt)
p=np.ones_like(tg)
for x in g: p*=1-(tg/x)**2
s2=varX/2-(g**-2.).sum(); s4=(np.log(g[-1]/2/np.pi)+1/3)/(2*np.pi*3*g[-1]**3); ph=p*np.exp(-s2*tg*tg-s4*tg**4/2)
def Phi(u):
    k=np.arange(1,14)[:,None]; return (2*(2*np.pi**2*k**4*np.exp(4.5*u)-3*np.pi*k**2*np.exp(2.5*u))*np.exp(-np.pi*k**2*np.exp(2*u))).sum(axis=0)
xs=np.linspace(0,3,301); xih=0.4971207781883141; fX=Phi(xs)/xih
def mob(N):
    m=np.ones(N+1,int); s=np.ones(N+1,bool)
    for q in range(2,N+1):
        if s[q]:
            s[2*q::q]=False; m[q::q]*=-1; m[q*q::q*q]=0
    return m
mu=mob(300); mK=6/np.pi**2-np.cumsum(np.where(np.arange(301)>0,mu/np.maximum(np.arange(301),1)**2.,0)); z4=np.pi**4/90
for K in [1,2,3,5,10,30,11,17,29,42,150]:
    U=min((K+1)*g1/2,T/2); t=tg[tg<=U]; ps=np.exp(-sig2*t*t/2)
    for k in range(1,K+1):
        if mu[k]==0: continue
        for x in gT:
            sc=np.sinc(t/(k*x)); ps=ps*sc if mu[k]>0 else ps/np.where(np.abs(sc)<1e-300,1e-300,sc)
    w=np.full_like(t,dt); w[0]=w[-1]=dt/2
    act=np.abs((np.cos(np.outer(xs,t))*(w*ps)).sum(axis=1)/np.pi-fX).max()
    aph=np.abs(ph[tg<=U]); B=2/3*t**4*4.873e-7
    tail=np.trapezoid(np.abs(ph[tg>=U]),tg[tg>=U])/np.pi
    A1=np.pi**2*t*t/6*abs(mK[K])*S; A2=2*z4/(9*K**3)*t**4*S4
    ref=np.trapezoid(aph*(np.exp(A1+A2+B)-1),t)/np.pi+tail
    crude_code=np.trapezoid(aph*(np.exp(2*np.pi**2*t*t/(9*K)*S+B)-1),t)/np.pi+tail
    crude_text=np.trapezoid(aph*(np.exp(np.pi**2*t*t/6/K*S+A2+B)-1),t)/np.pi+tail
    print(f'K={K} m={mK[K]:+.3e} act={act:.4e} refined={ref:.4e} ratio={ref/act:.2f} parts: m-term only={np.trapezoid(aph*(np.exp(A1)-1),t)/np.pi:.3e} K^-3 term only={np.trapezoid(aph*(np.exp(A2)-1),t)/np.pi:.3e} B only={np.trapezoid(aph*(np.exp(B)-1),t)/np.pi:.3e} tail={tail:.3e} crude(code 4/(3K))={crude_code:.3f} crude(text 1/K)={crude_text:.3f}')
# R(eta) expansion check
for d,L in [(.5,14.5087),(.05,14.5087)]:
    x=np.linspace(0,L,200001);w=np.exp(-d*x)*(1-x/L);m2=np.trapezoid(w*x*x,x)/np.trapezoid(w,x)
    print('d',d,'L',L,'E[x^2]',m2,'1/d^2*2',2/d/d,'L^2/6',L*L/6,'R(.01)-1',np.trapezoid(w*np.cosh(.01*x),x)/np.trapezoid(w,x)-1,'eta^2/d^2',1e-4/d/d,'eta^2 L^2/12',1e-4*L*L/12)
