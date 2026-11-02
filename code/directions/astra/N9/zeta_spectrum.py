import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import json
import math
import numpy as np
from scipy.optimize import minimize_scalar
from zeta_review import psi_grid, covariance

d=.025
m=480
c=covariance(psi_grid(d,m))
k=np.arange(1,m)
w=(1+np.cos(math.pi*k/m))/2
def spectrum(r):
    return c[0]+2*np.dot(w*c[1:],np.cos(r*k*d))
data=json.loads(Path('zero_sums.json').read_text())
rows=[]
for s in data['zeros'][:13]:
    g=float(s)
    opt=minimize_scalar(lambda r:-spectrum(r),bounds=(g-.3,g+.3),method='bounded',options={'xatol':1e-11})
    rows.append({'gamma':g,'peak':float(opt.x),'error':float(opt.x-g)})
Path('spectrum_results.json').write_text(json.dumps({'hann_lag_definition':'w(m)=0.5*(1+cos(pi*m/M)), M=480','rows':rows,'max_abs_error':max(abs(r['error']) for r in rows)},indent=2))
print(rows)
