"""Zero-field ferromagnetic fits. All target inputs are theta kernels, never zeros."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS', '1')
from pathlib import Path
import argparse
import json
import platform
import time
import hashlib
import numpy as np
import scipy
from scipy.special import gammaln, roots_legendre, logsumexp
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'
RESULTS.mkdir(exist_ok=True)

def kernel(name, x):
    x = np.asarray(x)
    n = np.arange(1, 25 if name == 'polya' else 49, dtype=float)
    e = np.exp(2*x)[..., None]
    if name == 'polya':
        u = np.pi * n*n * e
        return np.sum((4*u*u - 6*u)*np.exp(x[..., None]/2 - u), axis=-1)
    kap = (np.sqrt(10-2*np.sqrt(5))-2)/(np.sqrt(5)-1)
    a = np.array([0, 1, kap, -kap, -1])[n.astype(int) % 5]
    return 2*np.exp(1.5*x)*np.sum(n*a*np.exp(-np.pi*n*n*e/5), axis=-1)

def quadrature(n, hi):
    q, w = roots_legendre(n)
    return (q+1)*hi/2, w*hi/2

class Target:
    def __init__(self, name, nodes=512):
        self.name = name
        x, w = quadrature(nodes, 4)
        w = 2*w*kernel(name, x)
        assert np.all(w >= 0)
        self.norm = float(w.sum())
        self.p = w/w.sum()
        self.var = float(self.p @ x**2)
        self.sd = np.sqrt(self.var)
        self.y = x/self.sd
        self.moments = np.array([self.p @ self.y**k for k in range(2, 13, 2)])
    def cf(self, t):
        return np.cos(np.outer(t, self.y)) @ self.p

def frequencies(n=64):
    t, w = quadrature(n, 12)
    return t, w*np.sqrt(2/np.pi)*np.exp(-t*t/2)

class SpinModel:
    def __init__(self, family, n):
        self.family, self.n = family, n
        if family == 'cw':
            k = np.arange(n+1)
            self.s = (2*k-n).astype(float)
            self.logmult = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
            self.features = (self.s*self.s/(2*n))[:, None]
        elif family == 'block':
            assert n % 2 == 0
            m = n//2
            k = np.arange(m+1)
            s1, s2 = np.meshgrid(2*k-m, 2*k-m, indexing='ij')
            s1, s2 = s1.ravel().astype(float), s2.ravel().astype(float)
            lb = gammaln(m+1)-gammaln(k+1)-gammaln(m-k+1)
            self.logmult = (lb[:,None]+lb[None,:]).ravel()
            self.s = s1+s2
            self.features = np.stack([s1*s1/(2*m), s2*s2/(2*m), s1*s2/m], axis=1)
        elif family == 'dense':
            # Fix the first spin using global sign symmetry of energy and cos.
            ints = np.arange(1 << (n-1), dtype=np.uint64)
            spins = np.ones((len(ints), n))
            spins[:,1:] = 2*((ints[:,None] >> np.arange(n-1, dtype=np.uint64)) & 1).astype(float)-1
            ii, jj = np.triu_indices(n, 1)
            self.features = spins[:,ii]*spins[:,jj]
            self.s = spins.sum(axis=1)
            self.logmult = np.zeros(len(ints))
        else:
            raise ValueError(family)
        self.s2 = self.s*self.s
        self.dim = self.features.shape[1]
    def law(self, theta):
        lp = self.logmult + self.features @ theta
        p = np.exp(lp-logsumexp(lp))
        v = p @ self.s2
        return self.s/np.sqrt(v), p, v
    def cf_and_jac(self, theta, t):
        y, p, v = self.law(theta)
        ty = np.outer(t,y)
        c = np.cos(ty)
        cf = c @ p
        ef = p @ self.features
        covcf = (c*p) @ self.features - cf[:,None]*ef
        covv = (p*self.s2) @ self.features - v*ef
        jac = covcf + ((ty*np.sin(ty)) @ p)[:,None]*covv[None,:]/(2*v)
        return cf, jac

def distance(model, theta, target, nq=64):
    t, w = frequencies(nq)
    y, p, v = model.law(theta)
    delta = np.cos(np.outer(t,y)) @ p-target.cf(t)
    return float(np.sqrt(w @ (delta*delta)))

def w1_grid(y, p, target, step=0.0002):
    # Symmetrise the sign-fixed enumeration before computing its CDF.
    y = np.concatenate([y,-y])
    p = np.concatenate([p,p])/2
    order = np.argsort(y)
    y, p = y[order], p[order]
    cp = np.cumsum(p)
    hi = max(12.0, float(np.max(np.abs(y)))+step)
    x = np.arange(-hi, hi+step, step)
    dens = target.sd*kernel(target.name, np.abs(x)*target.sd)/target.norm
    cdf = cumulative_trapezoid(dens, x, initial=0)
    idx = np.searchsorted(y,x,side='right')-1
    mcdf = np.where(idx>=0,cp[np.maximum(idx,0)],0)
    return float(np.trapezoid(np.abs(cdf-mcdf),x))

def fit(model, target, starts=4):
    t,w = frequencies()
    sqrtw = np.sqrt(w)
    targetcf = target.cf(t)
    def fun(theta):
        cf,jac = model.cf_and_jac(theta,t)
        r=(cf-targetcf)*sqrtw
        return 1e6*float(r @ r), 2e6*(jac*sqrtw[:,None]).T @ r
    out=[]
    if model.family=='cw':
        grid=np.linspace(0,3,61)
        costs=np.array([fun(np.array([b]))[0] for b in grid])
        j=int(np.argmin(costs))
        opt=minimize_scalar(lambda b: fun(np.array([b]))[0], bounds=(grid[max(0,j-1)],grid[min(len(grid)-1,j+1)]), method='bounded',options={'xatol':2e-13})
        cand=[(float(opt.fun),float(opt.x)),(float(costs[0]),0.0),(float(costs[-1]),3.0)]
        val,b=min(cand)
        out=[{'theta':[b], 'D':np.sqrt(val/1e6),'success':bool(opt.success),'message':str(opt.message),'nit':int(opt.nfev),'seed':None}]
    else:
        for seed in range(1729,1729+starts):
            rng=np.random.default_rng(seed)
            if model.family=='block':
                initial=np.array([0.3,0.3,0.3]) if seed==1729 else rng.uniform(0,1.1,3)
                upper=4.0
            else:
                initial=np.full(model.dim,0.3/model.n) if seed==1729 else rng.exponential(0.6/model.n,model.dim)
                upper=4.0
            opt=minimize(fun,initial,jac=True,bounds=[(0,upper)]*model.dim,method='L-BFGS-B',options={'maxiter':700,'maxfun':1500,'ftol':1e-15,'gtol':1e-8,'maxls':40,'maxcor':15})
            out.append({'theta':opt.x.tolist(),'D':np.sqrt(opt.fun/1e6),'success':bool(opt.success),'message':str(opt.message),'nit':int(opt.nit),'seed':seed})
    best=min(out,key=lambda r:r['D'])
    theta=np.asarray(best['theta'])
    y,p,v=model.law(theta)
    mm=np.array([p @ y**k for k in range(2,13,2)])
    ans={'target':target.name,'family':model.family,'N':model.n,'D':best['D'],'D192':distance(model,theta,target,192),'theta':theta.tolist(),'spin_weight':float(target.sd/np.sqrt(v)),'raw_variance':float(v),'standard_moments':mm.tolist(),'moment_relative_errors':(mm/target.moments-1).tolist(),'W1_grid':w1_grid(y,p,target),'starts':out}
    if model.family=='cw':
        ans['b_window']=float((theta[0]-1)*np.sqrt(model.n))
        ans['c_window']=float(ans['spin_weight']*model.n**0.75)
    return ans

def quartic_law(b, nodes=512):
    x,w=quadrature(nodes,12)
    lp=-x**4/12+b*x*x/2
    p=w*np.exp(lp-np.max(lp))
    p/=p.sum()
    v=p @ x**2
    return x/np.sqrt(v),p,v

def fit_quartic(target):
    t,w=frequencies()
    tc=target.cf(t)
    def obj(b):
        y,p,v=quartic_law(b)
        r=np.cos(np.outer(t,y)) @ p-tc
        return w @ (r*r)
    opt=minimize_scalar(obj,bounds=(-50,5),method='bounded',options={'xatol':1e-11})
    y,p,v=quartic_law(opt.x)
    tt,ww=frequencies(192)
    r=np.cos(np.outer(tt,y)) @ p-target.cf(tt)
    return {'target':target.name,'b':float(opt.x),'c':float(target.sd/np.sqrt(v)),'D':float(np.sqrt(opt.fun)),'D192':float(np.sqrt(ww @ (r*r))),'standard_moments':[float(p @ y**k) for k in range(2,13,2)]}

def emit_summary(data):
    lines=['Floating-point exploratory fits, not interval certificates.', 'family target N D64 D192 W1_grid parameters spin_weight b_window c_window']
    for r in data['fits']:
        pars=','.join(f'{a:.10g}' for a in r['theta']) if r['family']!='dense' else '[full array in fits.json]'
        lines.append(f"{r['family']:6s} {r['target']:5s} {r['N']:5d} {r['D']:.12g} {r['D192']:.12g} {r['W1_grid']:.12g} {pars} {r['spin_weight']:.12g} {r.get('b_window',float('nan')):.12g} {r.get('c_window',float('nan')):.12g}")
    for r in data['quartic']:
        lines.append('quartic '+json.dumps(r,sort_keys=True))
    for key in ['targets','validation','runtime','environment']:
        if key in data: lines.append(key+' '+json.dumps(data[key],sort_keys=True))
    (RESULTS/'numbers.txt').write_text('\n'.join(lines)+'\n')
    (RESULTS/'fits.json').write_text(json.dumps(data,indent=2))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--targets',action='store_true')
    parser.add_argument('--families',default='cw,block,dense')
    parser.add_argument('--starts',type=int,default=4)
    parser.add_argument('--dense-n',default='4,6,8,10,12')
    parser.add_argument('--block-n',default='8,16,32,64,128')
    parser.add_argument('--output-dir',default='results')
    args=parser.parse_args()
    global RESULTS
    RESULTS=HERE/args.output_dir
    RESULTS.mkdir(exist_ok=True)
    tic=time.monotonic()
    targets={name:Target(name) for name in ['polya','dh']}
    data={'targets':{},'fits':[],'quartic':[],'validation':{},'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    for name,targ in targets.items():
        other=Target(name,1024)
        data['targets'][name]={'norm':targ.norm,'variance':targ.var,'standard_moments':targ.moments.tolist()}
        data['validation'][name+'_moment_quadrature_change']=float(np.max(np.abs(other.moments/targ.moments-1)))
        data['validation'][name+'_cf_quadrature_change']=float(np.max(np.abs(other.cf(np.linspace(0,12,101))-targ.cf(np.linspace(0,12,101)))))
    # Direct finite-difference check of derivatives, including the weight change.
    mod=SpinModel('block',12)
    th=np.array([0.35,0.47,0.18])
    tt,ww=frequencies()
    cf,jac=mod.cf_and_jac(th,tt)
    numerical=np.stack([(mod.cf_and_jac(th+1e-5*np.eye(3)[j],tt)[0]-mod.cf_and_jac(th-1e-5*np.eye(3)[j],tt)[0])/2e-5 for j in range(3)],axis=1)
    data['validation']['jacobian_max_error']=float(np.max(np.abs(jac-numerical)))
    assert data['validation']['jacobian_max_error']<1e-7
    print(json.dumps({'targets':data['targets'],'validation':data['validation'],'environment':data['environment']},indent=2),flush=True)
    if args.targets:
        (RESULTS/'targets_probe.json').write_text(json.dumps(data,indent=2))
        return
    for name,targ in targets.items():
        qr=fit_quartic(targ)
        data['quartic'].append(qr)
        print('quartic',qr,flush=True)
    sizes={'cw':[4,8,16,32,64,128,256,512,1024,4096],'block':[int(n) for n in args.block_n.split(',')],'dense':[int(n) for n in args.dense_n.split(',')]}
    for fam in args.families.split(','):
        for n in sizes[fam]:
            mod=SpinModel(fam,n)
            for name,targ in targets.items():
                r=fit(mod,targ,args.starts)
                data['fits'].append(r)
                data['runtime']=time.monotonic()-tic
                emit_summary(data)
                print(f"{fam} {name} N={n} D={r['D']:.10g} D192={r['D192']:.10g} W1={r['W1_grid']:.10g} min/max_start={min(s['D'] for s in r['starts']):.8g}/{max(s['D'] for s in r['starts']):.8g} seconds={data['runtime']:.2f}",flush=True)
    print('complete',data['runtime'],flush=True)

if __name__=='__main__':
    main()
