"""Normalize review math delimiters and check local report hygiene."""
from pathlib import Path
import re

root = Path(__file__).resolve().parent
p = root / "report.md"
s = p.read_text()
s = s.replace("\\\\(", "$").replace("\\\\)", "$")
bt = chr(96)
protected = re.compile(
    "(" + bt*3 + r"[\s\S]*?" + bt*3
    + r"|\\\[[\s\S]*?\\\]|" + bt + "[^" + bt + "]*" + bt
    + r"|\$[^\n$]*\$|\[[^\]]*\]\([^)]*\))"
)

def normalize(segment):
    stack = []
    pairs = []
    for j, ch in enumerate(segment):
        if ch == "(":
            stack.append(j)
        elif ch == ")" and stack:
            start = stack.pop()
            if not stack:
                body = segment[start+1:j]
                if ("\\" in body or "_" in body or "^" in body or "=" in body
                    or body in {"V", "X", "Y", "Z", "N", "K", "P", "R", "M", "t", "u", "z", "s", "x", "y", "g", "m", "k", "c", "v"}
                    or re.fullmatch(r"-?\d+\.\d+", body)):
                    pairs.append((start,j,body))
    for a,b,body in reversed(pairs):
        segment = segment[:a]+"$"+body+"$"+segment[b+1:]
    return segment

parts = protected.split(s)
for j in range(0,len(parts),2):
    parts[j] = normalize(parts[j])
s = "".join(parts)
s = s.replace("(u\\downarrow0$", "$u\\downarrow0$")
s = s.replace("$U_*'''$0.01$=412.6470<438$", "$U_*'''(0.01)=412.6470<438$")
s = s.replace("(F(z))", "$F(z)$").replace("(F)", "$F$")
s = s.replace("(K(z))", "$K(z)$").replace("(K_X)", "$K_X$")
s = s.replace("(K(3))", "$K(3)$").replace("(U'')", "$U''$")
s = s.replace("(O(Nm))", "$O(Nm)$").replace("(m>0)", "$m>0$")
s = s.replace("The notes do not prove the converse, either generally or for $V$.",
              "The general converse is false, as the positive Gaussian-polynomial counterexample in N1.8 shows. The converse for the particular Pólya law remains open here.")
s = s.replace("Statement (e) is only a sufficient condition. The title should not suggest",
              "Statement (e) is only a sufficient condition. The general question in line 142 has a negative answer, since N1.8 exhibits a member of the real-zero class violating GHS; the Pólya-specific converse remains open. The title should not suggest")
s = s.replace("The natural open domain of absolute convergence is $\\Re z>1$; a claim about every boundary point requires separate care.",
              "The natural open domain of absolute convergence is $\\Re z>1$; a claim about every boundary point requires separate care. Line 243 must also avoid asserting that all zeros of the completed transform are on $\\Re z=0$, since that would assume RH.")
s = s.replace("The Ding–Song–Sun extension invoked in line 88 is a separate recent result, with the required even-potential and convex-derivative hypotheses.",
"""The Ding–Song–Sun extension invoked in line 88 is a separate recent result. Its theorem also assumes a positive-definite quadratic interaction matrix. For the superquadratic Pólya potential this is repaired by choosing $a$ so that $aI-J$ is positive definite and replacing the site potential by $U(u)-au^2/2$; the total Hamiltonian is unchanged. The modified potential remains even, unbounded above and has the same third derivative.""")
p.write_text(s)

issues = []
if "\u2014" in s or "\u002d\u002d\u002d" in s:
    issues.append("forbidden dash sequence")
in_code = False
in_display = False
for line_no,line in enumerate(s.splitlines(),1):
    if line.startswith(bt*3):
        in_code = not in_code
    if line.strip() == r"\[":
        in_display = True
    if not in_code and not in_display and line.count("$") % 2:
        issues.append(f"odd inline math delimiters at line {line_no}: {line}")
    if line.strip() == r"\]":
        in_display = False
for target in re.findall(r"\]\(([^)]+)\)",s):
    if "://" not in target and not (root/target).exists():
        issues.append("missing local link: "+target)
if not s.endswith("\n"):
    issues.append("missing final newline")
print("report lines:", len(s.splitlines()), "words:", len(s.split()))
print("hygiene issues:", len(issues))
for issue in issues:
    print(issue)
