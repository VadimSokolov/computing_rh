# nearest approach to rational ratio p/q, q<=12, among first 8000 ordinates, float128-free via longdouble
import numpy as np
from decimal import Decimal
R='./data/zeros_8000.txt'
G=np.array([float(l) for l in open(R) if l.strip()])
from math import gcd
best={}
for q in range(1,13):
    bq=(1e9,None)
    for p in range(q+1,int(q*G[-1]/G[0])+2):
        if gcd(p,q)!=1: continue
        tg=G*p/q
        mask=tg<=G[-1]+1
        t=tg[mask]
        if t.size==0: break
        i=np.searchsorted(G,t); i=np.clip(i,1,len(G)-1)
        d=np.minimum(abs(G[i]-t),abs(G[i-1]-t))
        j=int(np.argmin(d))
        if d[j]<bq[0]: bq=(float(d[j]),(p,q,j+1))
    best[q]=bq
    print(q,bq)
