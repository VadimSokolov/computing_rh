# Krein string with spectral measure sum_{n<=N} w delta_{lam_n}: h(-s)=w G'/G, G=prod(1+s/lam). Euclid on polynomials gives
# h=l0+1/(m1 s+1/(l1+1/(m2 s+...))). Calibration: lam_n=((n-1/2)pi)^2, w=2 is the truncated homogeneous string M(x)=x on [0,1].
import mpmath as mp, sys, bisect
Z=[float(l.split()[-1]) for l in open('./data/zeros_8000.txt') if l.strip() and l.strip()[0].isdigit()]
print('zeros',len(Z),Z[0])
def string(lam,w,dps):
    mp.mp.dps=dps
    lam=[mp.mpf(a) for a in lam]
    # G coefficients ascending
    g=[mp.mpf(1)]
    for a in lam:
        c=1/a;new=[mp.mpf(0)]*(len(g)+1)
        for i,x in enumerate(g):
            new[i]+=x;new[i+1]+=x*c
        g=new
    p=[w*(i+1)*g[i+1] for i in range(len(g)-1)]  # numerator
    A,B=g,p; m=[];l=[]
    # A/B = m s + rem/B ; then B/rem = l + rem2/rem
    while True:
        q=A[-1]/B[-1]; m.append(q)
        # A - q s B : degree drops by >=2 in top? A deg d, B deg d-1
        R=[A[i]-(q*B[i-1] if i>0 else 0) for i in range(len(A))]
        R=R[:-1]
        # R deg d-1 generically equal to deg B
        if len(R)==0 or len(B)==1 and all(abs(x)<mp.mpf(10)**(-dps//2)*abs(A[-1]) for x in R): 
            l.append(mp.mpf(0)); break
        q2=B[-1]/R[-1]; l.append(q2)
        S=[B[i]-q2*R[i] for i in range(len(B))]
        S=S[:-1]
        if len(S)==0: break
        A,B=R,S
    return m,l
def report(lam,w,dps,L,name):
    m,l=string(lam,w,dps)
    x=[];c=mp.mpf(0)
    for li in l: x.append(c); c+=li
    # beads at positions x_i (first at 0), lengths after it
    Lt=c; Ltrue=sum(w/mp.mpf(a) for a in lam)
    s=mp.mpf(1000)
    v=mp.mpf(0)
    for mi,li in zip(reversed(m),reversed(l)):
        v=li+v if False else v
    # backward evaluation h: start at the right end (Dirichlet): Y=0; Y->Y+l; Y->1/(1/Y+s m)
    Y=mp.mpf(0)
    for mi,li in zip(reversed(m),reversed(l)):
        Y=Y+li; Y=1/(1/Y+s*mi)
    ex=sum(w/(mp.mpf(a)+s) for a in lam)
    mean=sum(mi*(Lt-xi) for mi,xi in zip(m,x)); 
    print(name,'N',len(lam),'beads',len(m),'lengthsum',mp.nstr(Lt,10),'Sum w/lam',mp.nstr(Ltrue,10),'h(-1000) err',mp.nstr(abs(Y-ex)/ex,3),'int(L-y)dM',mp.nstr(mean,8),'Sum 1/lam',mp.nstr(sum(1/mp.mpf(a) for a in lam),8))
    xf=[float(a) for a in x];mf=[float(a) for a in m]
    out=[]
    Ms=0;cum=[]
    for a in mf: Ms+=a;cum.append(Ms)
    for fr in (.01,.1,.25,.5,.75,.9,.95,.99):
        X=fr*L;i=bisect.bisect_right(xf,X)
        big=max(mf[:i]) if i else 0
        out.append((fr,i,float('%.4g'%(cum[i-1] if i else 0)),float('%.3g'%big)))
    print('  (x/L, beads<=x, M(x), largest bead<=x)',out)
    print('  last 3 beads (x/L,mass)',[(round(xf[-k]/L,5),float('%.4g'%mf[-k])) for k in (3,2,1)], 'total mass',float('%.5g'%Ms))
L=0.023104993115419
for N in (50,100,200,400,800):
    report([Z[i]**2 for i in range(N)],mp.mpf(1),max(400,5*N),L,'zeta')
    report([((n-.5)*3.141592653589793238)**2 for n in range(1,N+1)],mp.mpf(2),max(400,5*N),1.0,'homog')
