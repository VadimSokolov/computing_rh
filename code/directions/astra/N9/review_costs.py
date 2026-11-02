"""Arithmetic cost estimates for the stated heuristic resolution scale."""
import json
import math
from pathlib import Path
rows=[]
for a in [1,2,math.log(9),3,4,6,8,10,12]:
    x=math.exp(a)
    m=math.ceil(2*a*x)
    digits=4*math.pi*x/math.log(10)
    bits=math.ceil(4*math.pi*x/math.log(2))
    rows.append({'a':a,'X':x,'T_star':2*math.pi*x,'max_delta':1/(2*x),'m_at_least':m,'decimal_digits_heuristic':digits,'bits_heuristic':bits,'dense_double_bytes':8*m*m,'dense_high_precision_bytes_heuristic':m*m*math.ceil(bits/8),'toeplitz_high_precision_bytes_heuristic':m*math.ceil(bits/8)})
Path('cost_results.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(rows[-1],indent=2))
