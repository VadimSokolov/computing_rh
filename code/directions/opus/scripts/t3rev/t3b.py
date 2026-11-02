# DH control. Run: timeout 600 python3 -I /tmp/t3rev/t3b.py
import numpy as np, mpmath as mp
mp.mp.dps=20
KAPPA=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
A=[mp.mpf(1),KAPPA,-KAPPA,mp.mpf(-1)]
def f_dh(s): return mp.mpf(5)**(-s)*sum(a*mp.zeta(s,mp.mpf(r+1)/5) for r,a in enumerate(A))
a=np.array([float(x) for x in A]+[0.0]); N=3000
n=np.arange(1,N+1); an=a[(n-1)%5]; ln=np.log(n)
# Cesaro smoothing for a bounded-partial-sum series
wt=an*(1-n/(N+1))
best=[]
for sig in (1.02,1.06,1.1,1.15,1.2,1.3):
    ts=np.arange(0,3000,0.05); vals=[]
    for i in range(0,len(ts),4000):
        t=ts[i:i+4000]; vals.append(np.abs((wt*n**(-sig)) @ np.exp(-1j*np.outer(ln,t)).T if False else np.exp(-1j*np.outer(t,ln)) @ (wt*n**(-sig))))
    v=np.concatenate(vals); j=np.argsort(v)[:5]
    print('sigma',sig,'min |f| ',[(round(float(ts[k]),2),round(float(v[k]),4)) for k in j[:3]])
    best+= [(float(v[k]),sig,float(ts[k])) for k in j[:3]]
best.sort(); roots=[]
for _,sig,t in best[:9]:
    try:
        r=mp.findroot(f_dh,mp.mpc(sig,t),tol=1e-12,verify=False)
        if abs(f_dh(r))<1e-8 and r.real>1 and r.imag>0 and all(abs(r-q)>1e-4 for q in roots): roots.append(r)
    except Exception as e: pass
print('roots beta>1',[mp.nstr(r,10) for r in roots])
def xif(s): return (mp.mpf(5)/mp.pi)**(s/2)*mp.gamma((s+1)/2)*f_dh(s)
dl=lambda s: mp.diff(lambda z: mp.log(xif(z)),s)
print('xi_DH real on real axis',mp.nstr(xif(mp.mpf(1)),10),'at 1/2',mp.nstr(xif(mp.mpf(.5)),10),' FE check',mp.nstr(abs(xif(mp.mpc(.3,7))/xif(mp.mpc(.7,7))),10))
if roots:
    r=max(roots,key=lambda z:z.real); be,ga=r.real,r.imag
    th=1-be
    print('use zero',mp.nstr(r,10),'theta at alpha=1',mp.nstr(th,6))
    print(' single factor at small s, s=1e-6:',mp.nstr((th**2+ga**2)/((th+mp.mpf('1e-3'))**2+ga**2),15),'(above 1 means not a Laplace transform)')
    v=mp.re(dl(mp.mpc(1,ga)))
    print(' u_1(gamma)=Re xi_f prime/xi_f(1+i gamma)/pi =',mp.nstr(v/mp.pi,6))
    print(' b_1 =',mp.nstr(mp.re(dl(mp.mpf(1))),8))
    al=be+mp.mpf('0.02'); print(' alpha above beta',mp.nstr(al,5),' u_alpha(gamma)',mp.nstr(mp.re(dl(mp.mpc(al,ga)))/mp.pi,6))
# first DH ordinates on line and doubling test
g=[]; t=mp.mpf(1); 
L=lambda t: mp.re((mp.mpf(5)/mp.pi)**(mp.mpc(.5,t)/2)*mp.gamma((mp.mpc(.5,t)+1)/2)*f_dh(mp.mpc(.5,t)))
prev=L(t)
while len(g)<6 and t<30:
    t2=t+0.1; cur=L(t2)
    if prev*cur<0: g.append(mp.findroot(L,(t,t2),solver='anderson'))
    t,prev=t2,cur
print('DH on-line ordinates',[mp.nstr(x,8) for x in g])
print('|f(1/2+2 i g1)|',mp.nstr(abs(f_dh(mp.mpc(.5,2*g[0]))),6),'|f(1/2+3 i g1)|',mp.nstr(abs(f_dh(mp.mpc(.5,3*g[0]))),6))
