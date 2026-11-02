import sys,time; sys.path.insert(0,'/tmp/t4rev')
from lib import *
g=zeros8000()
# Resolution law: zero side (exact identity under verified RH for these zeros), kernel delta=1/L.
def analyse(T0,W,L,d,dw=None):
    if dw is None: dw=min(0.02,0.25/L*2)
    wg=np.arange(T0,T0+W,dw)
    sel=g[(g>T0-400)&(g<T0+W+400)]
    # far part via mean density; near via zeros
    O=np.empty_like(wg)
    for i in range(0,len(wg),200):
        ww=wg[i:i+200][:,None]
        O[i:i+200]=np.real(H(ww-sel,d,L)).sum(axis=1)
    # far zeros: H~2(d+1/L)/r^2, density log(x/2pi)/2pi, both sides; plus H(w+g) negligible
    dens=np.log(max(T0,10)/(2*np.pi))/(2*np.pi)
    O+=2*(d+1/L)*dens*(1/(wg-(T0-400))*(T0>400)+1/(T0+W+400-wg))
    H0=float(np.real(H(0.,d,L)))
    N=O/H0
    idx=np.nonzero((N[1:-1]>N[:-2])&(N[1:-1]>=N[2:]))[0]+1
    z=g[(g>T0+1)&(g<T0+W-1)]
    gaps=np.diff(g[(g>T0)&(g<T0+W)])
    pk=wg[idx]; hp=N[idx]
    # matching: every zero has a peak within a quarter of its nearest gap
    nearest=np.array([min(x-g[g<x].max(),g[g>x].min()-x) for x in z])
    jpk=np.array([np.argmin(np.abs(pk-x)) for x in z])
    ok=np.abs(pk[jpk]-z)<np.minimum(nearest/2,2*np.pi/L/2)
    distinct=len(set(jpk[ok]))==ok.sum()
    allz=ok.all() and distinct
    zero_h=N[np.searchsorted(wg,z)]
    # spurious: peak (inside) not matched to a zero, height
    ins=(pk>T0+1)&(pk<T0+W-1)
    matched=np.zeros(len(pk),bool); matched[jpk[ok]]=True
    spur=hp[ins&~matched]
    minzero=hp[jpk[ok]].min() if ok.any() else np.nan
    return dict(sup=N.max(),mean=N.mean(),minzp=minzero,spurmax=spur.max() if len(spur) else 0,allz=allz,fr=ok.mean(),errmed=np.median(np.abs(pk[jpk[ok]]-z[ok])) if ok.any() else np.nan,mingap=gaps.min(),sep=(minzero>(spur.max() if len(spur) else 0)))
Ts=[(14,136),(150,150),(500,200),(1000,300),(2000,300),(4000,300),(8000,300)]
Ls=[2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,24,28,32,40,50,60,80,100,150]
print('T0 W logT/2pi mingap | L_all (all zeros resolved, threshold separates) C=L/log | L_2: sup<2 | L_med: 90% resolved')
for T0,W in Ts:
    lg=np.log((T0+W/2)/(2*np.pi)); La=None;L2=None;L9=None;rows=[]
    for L in Ls:
        r=analyse(T0,W,L,1/L); rows.append((L,r))
        if L9 is None and r['fr']>=.9: L9=L
        if L2 is None and r['sup']<2: L2=L
        if La is None and r['allz'] and r['sep']: La=L
    print(T0,W,round(lg,3),round(rows[0][1]['mingap'],4),'L_all',La,round(La/lg,2) if La else None,'pi/mingap*2',round(2*np.pi/rows[0][1]['mingap'],1),'L_sup<2',L2,round(L2/lg,2) if L2 else None,'L_90',L9,round(L9/lg,2) if L9 else None)
    for L,r in rows:
        if L in (4,8,14,20,32,60):
            print('   L',L,'mean',round(r['mean'],3),'(log/L',round(lg/L,3),') sup',round(r['sup'],3),'minzpk',round(r['minzp'],3),'spurmax',round(r['spurmax'],3),'frac',round(r['fr'],3),'errmed',round(r['errmed'],5))
# planted pairs at T=150,1000,8000: at L, delta: normalised height of own pair plus background vs sup
print('planted pair: replace nothing, add pair at gamma0 = midpoint of the largest gap in window; report peak/H0 and excess over on-line double (same place)')
for T0,W in [(100,50),(1000,100),(8000,100)]:
    z=g[(g>T0)&(g<T0+W)]; k=np.argmax(np.diff(z)); g0=(z[k]+z[k+1])/2
    sel=g[(g>T0-400)&(g<T0+W+400)]
    for L in [14.5,25,50]:
        d=1/L; wg=np.arange(g0-1,g0+1,0.002)
        base=np.real(H(wg[:,None]-sel,d,L)).sum(axis=1)
        H0=float(np.real(H(0.,d,L)))
        out=[]
        for eta in [0,0.001,0.003,0.01,0.03,0.05,0.1,0.2]:
            own=np.real(H(wg-g0-1j*eta,d,L)+H(wg-g0+1j*eta,d,L))
            N=(base+own)/H0; out.append((eta,round(float(N.max()),4)))
        # noise scale: sd of the peak heights of normal single zeros
        ww=z[:,None]; hz=np.real(H(ww-sel,d,L)).sum(axis=1)/H0
        print(T0,'g0',round(g0,3),'L',L,'single-zero heights mean',round(hz.mean(),4),'sd',round(hz.std(),4),'max',round(hz.max(),4),'pair heights',out)
