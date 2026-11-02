"""Local formatting, provenance, and evidence checks for this review."""
from pathlib import Path
import json

root=Path(__file__).resolve().parent
p=root/'report.md'
s=p.read_text()
s=s.replace('\tfrac',chr(92)+'tfrac')
lines=s.splitlines()
lines=[line.rstrip(chr(92))+chr(92)*2 if line.startswith('&') and line.endswith(chr(92)) else line for line in lines]
s='\n'.join(lines)+'\n'
p.write_text(s)
for name in ['gram_cusp_audit.md','reachability_audit.md']:
    q=root/name
    v=q.read_text().replace(':---',':--')
    q.write_text(v)
for name in ['report.md','gram_cusp_audit.md','reachability_audit.md','empirical.md']:
    v=(root/name).read_text()
    assert chr(0x2014) not in v, name
    assert '-'*3 not in v, name
    assert '\t' not in v, name
    assert v.count(chr(92)+'[')==v.count(chr(92)+']'), name
    assert v.count(chr(92)+'(')==v.count(chr(92)+')'), name
print('Formatting: no em dashes, triple hyphens, tabs, or unbalanced math delimiters.')
for name in ['moments','zeros','dh','folded','extra','contour','winding']:
    d=json.loads((root/'results'/(name+'.json')).read_text())
    print(name, 'seconds', d.get('seconds'))
cont=json.loads((root/'results/contour.json').read_text())
assert cont['1']['count']==336 and cont['2']['count']==335
low=json.loads((root/'results/contour_low.json').read_text())
assert low['1']['count']==415 and low['1']['real_bounds']==[1,700]
wind=json.loads((root/'results/winding.json').read_text())
for key,item in wind.items():
    if key=='seconds': continue
    assert abs(item['sampled_winding']-(336 if key.startswith('(1,') else 335))<1e-8
extra=json.loads((root/'results/extra.json').read_text())
assert len(extra['gram_interval_assignments_ball_checked'])==414
z=json.loads((root/'results/zeros.json').read_text())
assert len(z['N1']['brackets'])==415
assert z['N1']['sign_changes_to_100']==29 and z['N1']['sign_changes_to_260']==114
print('Key result files and counts agree.')
print('Report words:',len(s.split()))
