# Q5: is a prime power a surprise to the Levinson predictor? Grid chosen so that the whole spike sits on one lag.
import sys; sys.path.append('./misc/research/directions/opus')
from N9_common import *
ns,Ls=prime_powers(200000); ln=np.log(ns); u=Ls/np.sqrt(ns); isp=set(int(x) for x in ns)
print('A. f=(c_no-pred)/s: fraction of the spike forecast from the past; admissible interval is pred +- v; s/v is the naive partial correlation; k_no is the partial correlation if the spike is omitted')
for d0 in [0.025,0.00625]:
    print(' d0=%g\n  n k d s c_no pred f v s/v k_with k_no'%d0)
    for n0 in [2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49,53,59,61,64,67,71,73,79,81,83,89,97,101,127,128,251,256,509,1009,2003,4001,8009]:
        k=int(round(math.log(n0)/d0)); d=math.log(n0)/k; i=int(np.searchsorted(ns,n0))
        A=cov_arch_zeta(d,k+1); add_primes(A,d,ln[:i],u[:i]); ks,vs,f=levinson(A)
        if f is not None and f<k: print('  %d past break'%n0); continue
        v=vs[k-1]; pred=A[k]-ks[k]*v; s=0.5*u[i]*d
        print('  %5d %5d %.5f %.5f %+.5f %+.5f %6.3f %.5f %.3f %+.3f %+.3f'%(n0,k,d,s,A[k],pred,(A[k]-pred)/s,v,s/v,(A[k]-s-pred)/v,(A[k]-pred)/v))
print('B. control: n not a prime power. g=(c-pred)/s_n, s_n the spike a prime at n would carry. g near 0 means no false alarm')
for d0 in [0.025,0.00625]:
    out=[]
    for n0 in [6,10,12,14,15,18,20,21,22,24,26,28,30,33,34,35,36,38,39,40,42,44,45,46,48,50,51,52,54,55,56,57,58,60,62,63,65,66,68,69,70,72,74,75,76,77,78,80,82,84,85,86,87,88,90,91,92,93,94,95,96,98,99,100,250,500,1000,2000,4000,8000]:
        k=int(round(math.log(n0)/d0)); d=math.log(n0)/k; i=int(np.searchsorted(ns,n0))
        A=cov_arch_zeta(d,k+1); add_primes(A,d,ln[:i],u[:i]); ks,vs,f=levinson(A)
        if f is not None and f<k: out.append('%d:break'%n0); continue
        out.append('%d:%.2f'%(n0,ks[k]*vs[k-1]/(0.5*math.log(n0)/math.sqrt(n0)*d)))
    print(' d0=%g '%d0+' '.join(out))
print('C. share of the information -log(v_K/v_0) by lags: spike lags against all lags')
for d,amax in [(0.025,5.0),(0.025,3.0),(0.00625,5.0),(0.00625,3.0)]:
    K=int(amax/d); c=cov_arch_zeta(d,K+1); add_primes(c,d,ln,u); ks,vs,f=levinson(c)
    I=-np.log(1-ks[1:K+1]**2)
    sp=set()
    for l in ln[ln<amax]:
        m0=int(math.floor(l/d)); sp.add(m0); sp.add(m0+1)
    idx=np.array(sorted(x for x in sp if 1<=x<=K)); mask=np.zeros(K+1,bool); mask[idx]=True
    k2=int(math.log(2)/d)
    print(' d=%g a<=%g: total %.4f; before log 2 %.4f (v/v0 at log 2 = %.4f); on %d spike lags of %d: %.4f; share of what falls after log 2: %.3f; lag share after log 2 %.3f; v_K/v_0=%.4f'%(d,amax,I.sum(),I[:k2].sum(),vs[k2]/vs[0],len(idx),K,I[mask[1:]].sum(),I[mask[1:]].sum()/I[k2:].sum(),len(idx)/(K-k2),vs[K]/vs[0]))
print('D. DH peel: subtract the listed zeros and find the next break, delta on three grids, a<=14')
N=int(math.exp(14.2))+2; Lam=dh_coeffs(N)
def plant(d,m,gam,beta):
    z=gam+1j*beta; w=(1-np.cos(z*d))/z**2; return 4*np.real(w*np.cos(np.arange(m)*z*d))
def peaks(uu,d,k=3):
    n=1<<int(math.ceil(math.log2(16*len(uu)+1))); F=np.abs(np.fft.rfft(uu,n))
    idx=np.nonzero((F[1:-1]>F[:-2])&(F[1:-1]>F[2:]))[0]+1; idx=idx[np.argsort(F[idx])[::-1]][:k]
    return [round(float(2*math.pi*i/n/d),2) for i in idx]
L=[(85.69934849,.30851718),(240.404672,.369531),(114.16334273,.15083008),(176.70246124,.22425769),(166.47930591,.07435605),(983.752193,.262499),(1012.01998149,.19140481),(1018.79576231,.19809797),(1029.004692,.086137)]
for d in [0.025,0.02,0.0175]:
    M=int(14.0/d); c=cov_dh(d,M,Lam)
    for j in range(len(L)+1):
        ks,vs,f=levinson(c)
        if f is None: print(' d=%g removed %d: PD to %.2f'%(d,j,M*d))
        else:
            m=f+1; w,V=np.linalg.eigh(toeplitz(c[:m])); print(' d=%g removed %d (last %s): first non PD %.4f; peaks %s; aliases of next zero on the list: %s'%(d,j,('%.2f'%L[j-1][0]) if j else '-',m*d,peaks(V[:,0],d),('%.2f'%abs(((L[j][0]+math.pi/d)%(2*math.pi/d))-math.pi/d)) if j<len(L) else '-'))
        if j<len(L): c=c-plant(d,M,L[j][0],L[j][1])
