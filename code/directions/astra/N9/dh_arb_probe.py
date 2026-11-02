import sys
sys.path.insert(0, '/private/tmp/d1lib')
from flint import arb, arb_mat, ctx
ctx.prec = 128
print('constants', [s for s in dir(arb) if s in ['pi', 'const_euler', 'const_catalan', 'sqrt', 'log', 'exp', 'lower', 'upper', 'mid', 'rad']])
print('matrix', [s for s in dir(arb_mat) if any(k in s for k in ['chol', 'eig', 'solve', 'inv', 'positive'])])
print('numbers', arb.pi(), arb.const_euler(), arb.const_catalan(), arb('0.25').sqrt())
