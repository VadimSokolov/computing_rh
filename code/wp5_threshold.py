# wp5_threshold.py: where does the decay of eps_N in N start?  tau = 2L/log q around 2g, genus 1, 2, 3 at q = 13, N = 8, 16, 32 (laptop probe, mpmath, computed only)
import sys, json, math
sys.path.insert(0, '../code')
import wp5_curves as w
import mpmath as mp
counts = json.load(open('wp5_counts.json'))
cur = {1: [r for r in counts if r['g'] == 1 and r['q'] == 13 and not r.get('rejected')][0]['a'],
       2: [r for r in counts if r.get('named') and r['p'] == 13][0]['a'],
       3: [r for r in counts if r['g'] == 3 and r['q'] == 13 and not r.get('rejected') and 'named' not in r][0]['a']}
out = []
for g in (1, 2, 3):
    for dt in (-0.5, 0.0, 0.25, 0.5, 1.0, 1.5, 2.5):
        tau = 2 * g + dt
        row = []
        for N in (8, 16, 32):
            o = w.rung(13, g, cur[g], tau, N, 60 + int(0.16 * N * tau * math.log(13)) + 40, certify=False)
            row.append(o['log10eps'] if not o['neg'] else None)
        c = -(row[2] - row[1]) / 16 if None not in row else None
        out.append(dict(g=g, tau=tau, log10eps=row, c=c))
        print(g, tau, ['%.2f' % x if x is not None else 'NEG' for x in row], 'c=%s' % ('%.3f' % c if c is not None else None), flush=True)
json.dump(out, open('wp5_threshold.json', 'w'))
