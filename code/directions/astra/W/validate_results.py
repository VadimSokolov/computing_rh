"""Post-fit diagnostics. Zero data enter only this held-out validation script."""
from pathlib import Path
import importlib.util
import json
import time
import hashlib
import numpy as np
import mpmath as mp
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fit_spins',HERE/'fit_spins.py')
fs=importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)

def exact_cdf_w1(model,theta,target,step=0.0002):
    y,p,v=model.law(np.array(theta))
    y,p=np.r_[y,-y],np.r_[p,p]/2
    atoms,idx=np.unique(y,return_inverse=True)
    mass=np.bincount(idx,weights=p)
    cp=np.cumsum(mass)
    cp=np.clip(cp/cp[-1],0,1)
    hi=max(12,float(np.max(np.abs(atoms)))+1)
    grid=np.linspace(-hi,hi,int(np.ceil(2*hi/step))+1)
    dens=target.sd*fs.kernel(target.name,np.abs(grid)*target.sd)/target.norm
    f=cumulative_simpson(dens,x=grid,initial=0)
    den=f[-1]
    f=np.clip(f/den,0,1)
    fm=cumulative_simpson(grid*dens,x=grid,initial=0)/den
    g=grid*f-fm
    gx=PchipInterpolator(grid,g)
    # Do not differentiate inverse-CDF intervals at underflow-scale masses.
    # Their total contribution is negligible at the displayed W1 precision.
    mask=(f>1e-14)&(f<1-1e-14)
    unique,inds=np.unique(f[mask],return_index=True)
    inv=PchipInterpolator(np.r_[0,unique,1],np.r_[-hi,grid[mask][inds],hi])
    a,b=atoms[:-1],atoms[1:]
    c=cp[:-1]
    r=np.clip(inv(c),a,b)
    between=gx(a)+gx(b)-2*gx(r)+c*(2*r-a-b)
    return float(np.sum(between)+gx(atoms[0])+gx(atoms[-1])-atoms[-1])

def slope(ns,ds):
    return float(-np.polyfit(np.log(ns),np.log(ds),1)[0])

def diagnostics():
    base=json.loads((HERE/'results/fits.json').read_text())
    ext=json.loads((HERE/'results_extended/fits.json').read_text())
    out={'validation':{},'rates':{},'block_starts':[],'dense_couplings':{},'w1':[],'moments':[]}
    block=fs.SpinModel('block',12)
    dense=fs.SpinModel('dense',12)
    th=np.array([0.35,0.47,0.18])
    ii,jj=np.triu_indices(12,1)
    jjval=np.where(jj<6,th[0]/6,np.where(ii>=6,th[1]/6,th[2]/6))
    tt,_=fs.frequencies(192)
    bc=block.cf_and_jac(th,tt)[0]
    dc=dense.cf_and_jac(jjval,tt)[0]
    out['validation']['block_vs_direct_enumeration_cf_error']=float(np.max(np.abs(bc-dc)))
    assert out['validation']['block_vs_direct_enumeration_cf_error']<1e-12
    for name in ['polya','dh']:
        target=fs.Target(name)
        for fam,rows,ns in [('cw',base['fits'],[64,128,256,512,1024,4096]),('block',ext['fits'],[64,128,256,512]),('dense',base['fits'],[4,6,8,10,12])]:
            # Rate comparisons use only N with a completed result for both laws.
            ns=[n for n in ns if all(any(r['family']==fam and r['target']==law and r['N']==n for r in rows) for law in ['polya','dh'])]
            chosen=[r for r in rows if r['target']==name and r['family']==fam and r['N'] in ns]
            chosen.sort(key=lambda r:r['N'])
            out['rates'][fam+'_'+name]={'N':[r['N'] for r in chosen],'alpha_D':slope([r['N'] for r in chosen],[r['D192'] for r in chosen]),'alpha_W1grid':slope([r['N'] for r in chosen],[r['W1_grid'] for r in chosen])}
            out['rates'][fam+'_'+name]['adjacent_alpha']=[{'N':[a['N'],b['N']],'alpha':slope([a['N'],b['N']],[a['D192'],b['D192']])} for a,b in zip(chosen,chosen[1:])]
            if fam=='dense':
                out['rates'][fam+'_'+name]['per_start_alpha']=[{'seed':seed,'alpha':slope([r['N'] for r in chosen],[next(s['D'] for s in r['starts'] if s['seed']==seed) for r in chosen])} for seed in range(1729,1733)]
            if fam=='block':
                batch=[]
                for first in [1729,1733,1737]:
                    dd=[min(s['D'] for s in r['starts'] if first<=s['seed']<first+4) for r in chosen]
                    batch.append({'first_seed':first,'D':dd,'alpha':slope([r['N'] for r in chosen],dd)})
                out['rates'][fam+'_'+name]['optimizer_batches']=batch
                for r in chosen:
                    m=r['N']/2
                    aa,bb,cc=r['theta']
                    out['block_starts'].append({'target':name,'N':r['N'],'Dmin':min(s['D'] for s in r['starts']),'Dmax':max(s['D'] for s in r['starts']),'Dmedian':float(np.median([s['D'] for s in r['starts']])),'successes':sum(s['success'] for s in r['starts']),'total':len(r['starts']),'theta':r['theta'],'B_window':[(aa-1)*np.sqrt(m),(bb-1)*np.sqrt(m),cc*np.sqrt(m)],'c_window':r['spin_weight']*m**0.75})
            for r in chosen:
                out['moments'].append({'target':name,'family':fam,'N':r['N'],'relative_errors':r['moment_relative_errors']})
                if fam=='cw' or (fam=='block' and r['N'] in [64,128,256,512]):
                    model=fs.SpinModel(fam,r['N'])
                    w=exact_cdf_w1(model,r['theta'],target)
                    fine=exact_cdf_w1(model,r['theta'],target,0.0001)
                    out['w1'].append({'target':name,'family':fam,'N':r['N'],'W1':fine,'mesh_change':abs(w-fine)})
            if fam=='dense':
                r=chosen[-1]
                a=np.asarray(r['theta'])
                mat=np.zeros((r['N'],r['N']))
                ii,jj=np.triu_indices(r['N'],1)
                mat[ii,jj]=a
                mat+=mat.T
                out['dense_couplings'][name]={'N':r['N'],'matrix':mat.tolist(),'minimum':float(a.min()),'maximum':float(a.max()),'positive_count':int((a>1e-8).sum()),'total_edges':len(a),'sorted_degrees':np.sort(mat.sum(axis=0)).tolist()}
    p=out['rates']['block_polya']['optimizer_batches']
    d=out['rates']['block_dh']['optimizer_batches']
    out['rates']['block_paired_batch_alpha_dh_minus_polya']=[dd['alpha']-pp['alpha'] for pp,dd in zip(p,d)]
    out['validation']['frequency_quadrature_max_change']=max(abs(r['D']-r['D192']) for r in base['fits']+ext['fits'])
    out['validation']['tail_norm_bound']=str(mp.sqrt(4*mp.erfc(12/mp.sqrt(2))))
    return out

class QuarticTransform:
    def __init__(self,b,c,nodes):
        self.b,self.c=mp.mpf(str(b)),mp.mpf(str(c))
        x,w=mp.gauss_quadrature(nodes,'legendre')
        self.x=[4*(x[j]+1) for j in range(nodes)]
        self.w=[4*w[j]*mp.exp(-self.x[j]**4/12+self.b*self.x[j]**2/2) for j in range(nodes)]
        self.norm=mp.fsum(self.w)
    def __call__(self,t):
        t=mp.mpf(t)*self.c
        return mp.fsum(w*mp.cos(t*x) for x,w in zip(self.x,self.w))/self.norm

def zero_checks():
    data=json.loads((HERE/'results/fits.json').read_text())
    mp.mp.dps=80
    record={'precision_digits':80,'quadrature_nodes':[256,384],'integration_limit_u':8,'scan_step':'0.25','fit_sha256':hashlib.sha256((HERE/'results/fits.json').read_bytes()).hexdigest()}
    for name in ['polya','dh']:
        pars=next(q for q in data['quartic'] if q['target']==name)
        f=QuarticTransform(pars['b'],pars['c'],256)
        roots=[]
        oldt=mp.mpf(0)
        old=f(oldt)
        maxheight=90 if name=='polya' else 100
        for j in range(1,4*maxheight+1):
            t=mp.mpf(j)/4
            value=f(t)
            if old*value<0:
                root=mp.findroot(f,(oldt,t),tol=mp.mpf('1e-65'))
                assert oldt<root<t
                roots.append(root)
            oldt,old=t,value
        g=QuarticTransform(pars['b'],pars['c'],384)
        checked=[]
        for r in roots:
            rr=mp.findroot(g,(r-mp.mpf('.001'),r+mp.mpf('.001')),tol=mp.mpf('1e-65'))
            assert abs(r-rr)<mp.mpf('1e-10'), 'Quadrature does not resolve this zero.'
            checked.append({'root':str(rr),'quadrature_root_change':str(abs(r-rr)),'residual':str(abs(g(rr)))})
        entry={'b':pars['b'],'c':pars['c'],'roots':checked,'fit_physical_height_cutoff':float(12/np.sqrt(data['targets'][name]['variance']))}
        if name=='polya':
            comparisons=[]
            for k in [1,2,3,5,10,12,15,20]:
                z=mp.im(mp.zetazero(k))
                qr=mp.mpf(checked[k-1]['root'])
                comparisons.append({'index':k,'known_zeta_ordinate':str(z),'candidate_ordinate':str(qr),'signed_difference':str(qr-z),'candidate_cf_at_zeta_zero':str(g(z))})
            entry['comparisons']=comparisons
        record[name]=entry
        print('zeros complete',name,'count',len(roots),flush=True)
    kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
    co=[1,kap,-kap,-1]
    def dh(s):
        return mp.power(5,-s)*mp.fsum(co[j]*mp.zeta(s,mp.mpf(j+1)/5) for j in range(4))
    dz=mp.findroot(dh,(mp.mpc('.80','85.69'),mp.mpc('.82','85.71')),tol=mp.mpf('1e-65'))
    record['dh_off_line_validation']={'s':str(dz),'dirichlet_function_residual':str(abs(dh(dz))),'input_start_source':'Balanzario and Sanchez-Ortiz (2007), after all fits were frozen'}
    return record

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--zeros',action='store_true')
    args=parser.parse_args()
    tic=time.monotonic()
    if args.zeros:
        out=zero_checks()
        stem='zero_validation'
    else:
        out=diagnostics()
        stem='diagnostics'
    out['runtime_seconds']=time.monotonic()-tic
    out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'results'/f'{stem}.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':
    main()
