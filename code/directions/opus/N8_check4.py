#!/usr/bin/env python3
"""N8 check 4: (A) true level table at both cutoffs; (B) literal scale of the wall; (C) rank = span on finite sections, built explicitly.
Uses the functions of N8_check1.py by exec. Zeros are input in part C by construction. Run: timeout 300 python3 -I N8_check4.py"""
import io, contextlib, math
import numpy as np
H='./misc/research/directions/opus/'
g={}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(H+'N8_check1.py').read(),g)
Z,two,base,gram,lev_c,theta=g['Z'],g['two'],g['base'],g['gram'],g['lev_c'],g['theta']
print('gamma_400 = %.6f'%Z[399])
print('(A) levels up to cutoff: name, cutoff, intervals, without exactly one zero, first failure, counts 0/1/2, first level, last level')
gs=gram[gram>=17]
for nm,L0 in (('Gram from g_0',gs),('two term true from 18.717',two[two>=14]),('two term true with 7.604',two),('base true',base)):
    for top in (Z[399],680.0):
        L=L0[L0<=top]; c=np.searchsorted(Z,L[1:])-np.searchsorted(Z,L[:-1]); b=np.where(c!=1)[0]
        print('  %-26s %8.3f %4d %4d %7.2f %s %.3f %.3f'%(nm,top,len(L)-1,len(b),L[b[0]],np.bincount(c).tolist(),L[0],L[-1]))
print('(B) literal count in t=sqrt(mu) for the wall: nu<=T has about theta(2T)/pi+3/8 points, against theta(T)/pi')
for T in (100.,340.,680.):
    a2=np.ravel(theta(np.array([2*T])))[0]; a1=np.ravel(theta(np.array([T])))[0]; print('  T=%g  theta(2T)/pi=%.2f  theta(T)/pi=%.2f  ratio %.3f'%(T,a2/math.pi,a1/math.pi,a2/a1))
def chain(alpha,beta):
    """alpha, beta increasing, same length. Returns F with spec(diag(alpha)+F)=beta and rank p+q, from clipped counts."""
    n=len(alpha); U=np.unique(np.concatenate([alpha,beta])); Na=np.searchsorted(alpha,U,side='right'); Nb=np.searchsorted(beta,U,side='right'); M=np.minimum(Na,Nb)
    p=int((Na-Nb).max()); q=int((Nb-Na).max())
    def seq(C): return np.array([U[np.searchsorted(C,j)] for j in range(1,n+1)])
    path=[seq(np.minimum(Na,M+s)) for s in range(p,-1,-1)]+[seq(np.minimum(Nb,M+s)) for s in range(1,q+1)]
    Q=np.eye(n); cur=path[0].copy()
    assert np.array_equal(cur,alpha) and np.array_equal(path[-1],beta)
    A=np.diag(alpha).astype(float); F=np.zeros((n,n))
    for nxt in path[1:]:
        ci=np.setdiff1d(cur,nxt,assume_unique=True); ni=np.setdiff1d(nxt,cur,assume_unique=True)
        if len(ci)==0: cur=nxt; continue
        sigma=1.0 if ni[0]>ci[0] else -1.0
        w2=np.zeros(len(ci))
        for i,c in enumerate(ci):
            num=np.abs(ni-c); den=np.abs(np.delete(ci,i)-c); w2[i]=np.exp(np.log(num).sum()-np.log(den).sum())
        v=np.zeros(n); idx=np.searchsorted(cur,ci); v[idx]=np.sqrt(w2)
        B=np.diag(cur)+sigma*np.outer(v,v); ev,R=np.linalg.eigh(B)
        Qv=Q @ v; F+=sigma*np.outer(Qv,Qv); Q=Q @ R; cur=ev.copy(); assert np.abs(ev-nxt).max()<1e-6
        cur=nxt
    return F,p,q
print('(C) explicit F on [30,680]: name, n, D range, singular values of F (largest few and the first negligible ones), max |sqrt(eig)-gamma|')
tt=np.arange(30,680,0.0005); 
for nm,L0 in (('Gram',gs),('base true',base),('two term true',two)):
    L=L0[(L0>=30)&(L0<=680)]; zz=Z[Z>L[0]]; n=min(len(L),len(zz)); a=L[:n]**2; b=zz[:n]**2
    F,p,q=chain(a,b); sv=np.linalg.svd(F,compute_uv=False); ev=np.linalg.eigvalsh(np.diag(a)+F)
    print('  %-14s n=%d, N_nu-N_mu in [%d,%d], span %d; singular values'%(nm,n,-p,q,p+q),['%.3e'%x for x in sv[:p+q+2]],'; recovery error %.2e'%np.abs(np.sqrt(ev)-np.sqrt(b)).max())
