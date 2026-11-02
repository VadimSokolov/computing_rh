#!/usr/bin/env python3
"""N8: first 400 zeta ordinates, cached by this script alone.  Independent of misc/reviews.
Run: timeout 600 python3 -I N8_zeros.py"""
import mpmath as mp, json, sys, os, time
HERE='./misc/research/directions/opus'
mp.mp.dps=20
t0=time.time()
Z=[]
for n in range(1,421):
    Z.append(str(mp.im(mp.zetazero(n))))
    if time.time()-t0>560: break
json.dump(Z,open(os.path.join(HERE,'N8_zeros.json'),'w'))
print(len(Z),Z[0],Z[-1],round(time.time()-t0,1))
