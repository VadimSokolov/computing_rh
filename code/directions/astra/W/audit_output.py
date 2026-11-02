"""Final artifact and model sanity checks, no new optimisation."""
from pathlib import Path
import json
import re
import numpy as np

p=Path(__file__).resolve().parent
report=(p/'report.md').read_text()
assert '\u2014' not in report and '---' not in report
assert len(report.splitlines()[:5])==5
assert all(report.splitlines()[i].startswith(str(i+1)+'. ') for i in range(5))
links=re.findall(r'\]\(([^)]+)\)',report)
for link in links:
    if not link.startswith(('https://','http://')):
        assert (p/link).exists(),link
for path in ['results/fits.json','results_extended/fits.json']:
    data=json.loads((p/path).read_text())
    for fit in data['fits']:
        assert np.all(np.asarray(fit['theta'])>=0)
        assert fit['spin_weight']>0
        assert abs(fit['standard_moments'][0]-1)<1e-11
        assert fit['D192']>0
        assert abs(fit['D']-fit['D192'])<1e-10
z=json.loads((p/'results/zero_validation.json').read_text())
assert z['quadrature_nodes']==[256,384]
for name in ['polya','dh']:
    roots=[float(r['root']) for r in z[name]['roots']]
    assert all(a<b for a,b in zip(roots,roots[1:]))
    assert max(float(r['quadrature_root_change']) for r in z[name]['roots'])<1e-10
assert float(z['polya']['comparisons'][-1]['known_zeta_ordinate'])>z['polya']['fit_physical_height_cutoff']
d=json.loads((p/'results/diagnostics.json').read_text())
assert d['rates']['block_polya']['N']==d['rates']['block_dh']['N']
assert d['validation']['block_vs_direct_enumeration_cf_error']<1e-12
for name,description in [('polya',(3,2)),('dh',(5,10))]:
    mat=np.asarray(d['dense_couplings'][name]['matrix'])
    active=np.sum(np.any(mat>1e-8,axis=0))
    edges=np.sum(np.triu(mat,1)>1e-8)
    assert (int(active),int(edges))==description
print('PASS: five-line verdict, no em dashes, local links, coupling signs, variance matching, distance quadrature, paired rates, root checks, sparse components.')
print('Final root quadrature maximum:',max(float(r['quadrature_root_change']) for name in ['polya','dh'] for r in z[name]['roots']))
print('Report words:',len(report.split()))
