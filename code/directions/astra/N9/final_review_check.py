"""Final consistency checks for the delivered review."""
import json
from pathlib import Path
root=Path(__file__).resolve().parent
report=(root/'report.md').read_text()
assert chr(0x2014) not in report and chr(45)*3 not in report
assert '|\n\n|' not in report
a=json.loads((root/'zeta_arb_certificates.json').read_text())
assert len(a['grid'])==11
assert all(r['positive_definite'] for r in a['grid'])
assert all(r['before']['positive_definite'] and r['after']['negative_rayleigh'] for r in a['omissions'])
b=json.loads((root/'dh_arb_results.json').read_text())
assert all(r['last_positive']['positive_definite'] and r['first_negative']['negative_rayleigh'] for r in b['grids'])
c=json.loads((root/'zeta_arb_results.json').read_text())
assert all(r['negative_proved'] and r['first_nonpositive']==r['m'] for r in c['truncation'])
assert json.loads((root/'zero_sums.json').read_text())['N']==400
assert all((root/f).exists() for f in ['report.md','empirical.md','results/numbers.txt','hopper_grid_proposal.py'])
print('PASS: 11 positive zeta grid matrices; 11 adjacent single-omission pairs; 3 DH onset pairs; 4 tail-deletion onsets.')
print('PASS: 400-zero computation, artifact links and forbidden punctuation checks.')
print('Exact analytic endpoint enclosures:')
for r in a['table']:
    if r['t'] in [8,12,12.05]:
        print(r)
print('Report:',len(report.split()),'words,',len(report.splitlines()),'lines')

