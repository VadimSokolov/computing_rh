"""Direct zero sums, computed with mpmath, not a zero certificate."""
import json
import time
from pathlib import Path
import mpmath as mp
mp.mp.dps=20
ts=[mp.mpf(s) for s in ['.3','.7','1','2','3','5','8']]
sums=[mp.mpf(0) for t in ts]
zeros=[]
start=time.monotonic()
for n in range(1,401):
    gamma=mp.zetazero(n).imag
    zeros.append(str(gamma))
    for j,t in enumerate(ts):
        sums[j]+=2*(1-mp.cos(gamma*t))/(gamma*gamma)
    if n%50==0:
        print(n, 'height',gamma, 'seconds',time.monotonic()-start,flush=True)
        Path('zero_sums.json').write_text(json.dumps({'N':n,'gamma_N':str(gamma),'sum':{str(t):str(v) for t,v in zip(ts,sums)},'zeros':zeros},indent=2))
