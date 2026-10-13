#!/usr/bin/env python3
"""Write version.tex, the git commit the book is typeset from.

latexmkrc runs this before every build. main.tex prints the commit on the copyright page of the KDP interior and on the title page of the draft, and stores the full id in the PDF information under the key Commit.

The id gets -dirty appended when a tracked input of the build differs from the commit: the .tex files at the top (main.tex, main-kdp.tex, bookmeta.tex), latexmkrc, the folders ch/, fig/ and rzfig/ (the only folders the source reads figures from), and code/polya_flint.py and code/tilted.py, which the book typesets as listings. The publishing files (book.yaml, results.yaml, book.lock.json, summaries/, site/, the generated README block and the like), kdp.md, misc/, cover/ and the code are not inputs and are not checked, so refreshing them after a build does not mark it; untracked files are not checked either.
"""
import datetime
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS = [':(glob)*.tex', 'latexmkrc', 'ch', 'fig', 'rzfig', 'code/polya_flint.py', 'code/tilted.py']


def git(*args):
    return subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True, check=True).stdout.strip()


def day(d):
    return f'{d.day} {d:%B} {d.year}'


full = git('rev-parse', 'HEAD')
short = git('rev-parse', '--short=7', 'HEAD')
dirty = '-dirty' if git('status', '--porcelain', '--untracked-files=no', '--', *INPUTS) else ''
committed = datetime.date.fromisoformat(git('log', '-1', '--format=%cs', 'HEAD'))
lines = [
    '% Written by tools/stamp_version.py before every build (see latexmkrc); not under version control.',
    f'\\def\\bookcommit{{{short}{dirty}}}',
    f'\\def\\bookcommitfull{{{full}{dirty}}}',
    f'\\def\\bookcommitdate{{{day(committed)}}}',
    f'\\def\\bookbuilddate{{{day(datetime.date.today())}}}',
    f'\\def\\bookdirty{{{1 if dirty else 0}}}',
]
with open(os.path.join(ROOT, 'version.tex'), 'w') as f:
    f.write('\n'.join(lines) + '\n')
print(f'version.tex: commit {short}{dirty} of {day(committed)}')
