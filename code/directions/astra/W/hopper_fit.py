"""Unsubmitted Hopper proposal: one dense-graph fit per array task."""
import argparse
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fit_spins',HERE/'fit_spins.py')
fs=importlib.util.module_from_spec(spec)
spec.loader.exec_module(fs)

parser=argparse.ArgumentParser()
parser.add_argument('--task',type=int,required=True)
args=parser.parse_args()
assert 0<=args.task<96
rep=args.task%16
target_name=['polya','dh'][(args.task//16)%2]
n=[14,16,18][args.task//32]
seed=9000+rep
model=fs.SpinModel('dense',n)
target=fs.Target(target_name)
t,w=fs.frequencies(96)
sqrtw=np.sqrt(w)
tc=target.cf(t)

def loss(theta):
    cf,jac=model.cf_and_jac(theta,t)
    r=(cf-tc)*sqrtw
    return 1e6*float(r @ r),2e6*(jac*sqrtw[:,None]).T @ r

initial=np.random.default_rng(seed).exponential(0.6/n,model.dim)
tic=time.monotonic()
opt=minimize(loss,initial,jac=True,bounds=[(0,4)]*model.dim,method='L-BFGS-B',options={'maxiter':700,'maxfun':1500,'ftol':1e-15,'gtol':1e-8,'maxls':40,'maxcor':15})
y,p,v=model.law(opt.x)
record={'target':target_name,'N':n,'seed':seed,'theta':opt.x.tolist(),'spin_weight':float(target.sd/np.sqrt(v)),'D_fit':float(np.sqrt(opt.fun/1e6)),'D192':fs.distance(model,opt.x,target,192),'success':bool(opt.success),'message':str(opt.message),'nit':int(opt.nit),'seconds':time.monotonic()-tic}
out=HERE/'hopper_results'
out.mkdir(exist_ok=True)
(out/f'task_{args.task:03d}.json').write_text(json.dumps(record,indent=2))
print(json.dumps(record),flush=True)
