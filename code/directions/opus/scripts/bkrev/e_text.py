import re,sys
c=open('./misc/reviews/nick-berry-keating/chapter_berry_keating.tex').read()
bad=0
for l in open('/tmp/bkrev/BK.md'):
    if l.startswith('TEXT: "'):
        t=l[7:l.rstrip().rfind('"')]
        n=c.count(t)
        if n!=1: print(n,t[:90]);bad+=1
print('bad',bad, 'emdash',open('/tmp/bkrev/BK.md').read().count('—'))
