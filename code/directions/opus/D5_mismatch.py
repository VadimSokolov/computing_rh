# D5 experiment 8: multiplicity mismatch. Laws psi^a (sum of a copies of spin one with strip s) and psi^b, on the integer lattice, coupling J_phys=J*s^2, common field. mpmath 60 digits, not certified.
# Local limit s->0: Q1=(w^2-1)^a, Q2=(w^2-1)^b, Z=sum_k J^k/k! Q1^(k) Q2^(k) in units delta=1.
import mpmath as mp
mp.mp.dps=60
def conv(a,b):
    c=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c
def power(p,m):
    out=[mp.mpf(1)]
    for _ in range(m): out=conv(out,p)
    return out
def Ylaw(a,b,s,Jp):
    q=1/(2*(1+mp.cosh(s))); p=[q,1-2*q,q]
    P1,P2=power(p,a),power(p,b); L1,L2=a,b
    c=[mp.mpf(0)]*(len(P1)+len(P2)-1)
    for i,x in enumerate(P1):
        for j,y in enumerate(P2): c[i+j]+=x*y*mp.exp(Jp*(i-L1)*(j-L2))
    m=max(abs(x) for x in c); c=[x/m for x in c]
    r=mp.polyroots(c[::-1],maxsteps=4000,extraprec=800)
    return max(abs(mp.log(abs(z))) for z in r)
def Ychain3(s,Jp):
    q=1/(2*(1+mp.cosh(s))); p={-1:q,0:1-2*q,1:q}
    c=[mp.mpf(0)]*7
    for a in p:
        for b in p:
            for d in p: c[a+b+d+3]+=p[a]*p[b]*p[d]*mp.exp(Jp*(a*b+b*d))
    m=max(abs(x) for x in c); c=[x/m for x in c]
    r=mp.polyroots(c[::-1],maxsteps=4000,extraprec=800)
    return max(abs(mp.log(abs(z))) for z in r)
print("(Y/s)^2 for laws, s=strip of spin one; J_phys=J s^2")
print("s       J    (1,1)    f(J)    (2,1)   chain3   (3,1)   (4,2)  (2,2)")
for s in (mp.mpf(1),mp.mpf('0.3'),mp.mpf('0.1'),mp.mpf('0.03')):
    for J in (mp.mpf('0.5'),mp.mpf(1),mp.mpf(4),mp.mpf(16)):
        Jp=J*s*s
        vals=[(Ylaw(a,b,s,Jp)/s)**2 for (a,b) in ((1,1),(2,1))]
        ch=(Ychain3(s,Jp)/s)**2
        v2=[(Ylaw(a,b,s,Jp)/s)**2 for (a,b) in ((3,1),(4,2),(2,2))]
        f=(1-2*J+mp.sqrt(1+2*J*J))/2 if J<2 else mp.mpf(0)
        print(f"{float(s):5.2f} {float(J):5.1f} {float(vals[0]):8.5f} {float(f):7.4f} {float(vals[1]):8.5f} {float(ch):8.5f} {float(v2[0]):8.5f} {float(v2[1]):8.5f} {float(v2[2]):8.5f}")
print("local root of (u-1)^3+8Ju(u-1)+J^2(12u-4) in (1/3,1):")
for J in (0.5,1,4,16,100,1e4):
    u=mp.findroot(lambda u:(u-1)**3+8*J*u*(u-1)+J**2*(12*u-4),(mp.mpf(1)/3,mp.mpf(1)),solver='anderson')
    print(f"  J={J}: u={mp.nstr(u,8)}")
