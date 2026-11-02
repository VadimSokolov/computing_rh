import re
R='./'
t=open(R+'misc/reviews/nick-ferro-thorin/ferromagnetic_thorin_for_vadim/three_computations.tex').read()
r=open(R+'misc/research/directions/opus/T4.md').read()
for s in re.findall(r'^TEXT: "(.*)"$',r,flags=re.M): print(t.count(s),s[:70])
for s in ["\\frac{\\alpha+3}2\\log x_\\varepsilon+O(1)","$\\Sigma_T^{(4)}\\le4.87\\times10^{-7}$","The refined bound is within a factor of two of the truth throughout."]: print(t.count(s),s[:60])
print('emdash',r.count('—'),r.count('---'))
