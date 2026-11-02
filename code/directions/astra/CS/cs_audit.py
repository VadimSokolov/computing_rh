"""Mechanical checks on the delivered Track CS research artifacts."""
from pathlib import Path

root = Path(__file__).resolve().parent
for name in ('report.md', 'empirical.md', 'cs_certify.py', 'cs_probe.py'):
    content = (root/name).read_text()
    bad = [(i, ord(c)) for i, c in enumerate(content) if ord(c) < 32 and c not in '\n\t']
    print(name, 'control_characters', bad)
    print(name, 'forbidden_dash_count', content.count(chr(0x2014)) + content.count('-'*3))
    assert not bad
    assert chr(0x2014) not in content and '-'*3 not in content
print('Report lines', len((root/'report.md').read_text().splitlines()))
