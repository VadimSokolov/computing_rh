#!/usr/bin/env python3
"""Split Computing the Riemann Hypothesis (v1.4.1) into two volumes.

Run from the root of the computing_rh repository:

    python3 split_volumes.py              # write main_vol1.tex, main_vol2.tex, vol/ and split_report.md
    python3 split_volumes.py --build      # also typeset both volumes (pdflatex, three passes each)
    python3 split_volumes.py --check      # only check the plan against main.tex and the sources

Nothing in ch/ is modified. The script
  1. reads the chapter plan below (the plan of the note rh_two_volumes.pdf),
  2. checks that every file input by main.tex is assigned exactly once and that labels are unique,
  3. writes one wrapper file per new chapter into vol/ : a new \\chapter heading, then the old
     chapter openings (their \\chapter commands removed or turned into unnumbered sections),
     then \\input of the old section files in the new order,
  4. writes stub files vol/new/*.tex for the new writing (only if they do not exist yet, so
     drafts are never overwritten),
  5. writes main_vol1.tex and main_vol2.tex from the preamble and title pages of main.tex,
     loading xr-hyper so each volume resolves references into the other volume,
  6. writes split_report.md: chapter list, page estimates, every reference from Volume I into
     Volume II with file and line, and the labels to lift into Volume I.

Edit PLAN to change the order; rerun. The old main.tex still builds the single volume edition.
"""

import argparse
import collections
import json
import os
import re
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# The plan. Each chapter: (title, short running head or None, [source files], new stub or None)
# A source file name is a file of ch/ without .tex. Header files cNN carry the old chapter
# opening. ('PART', title) starts a part. ('APPENDIX',) starts the appendices.
# ---------------------------------------------------------------------------

VOLUMES = {
    "1": {
        "main": "main_vol1",
        "other": "main_vol2",
        "title": "Computing the Riemann Hypothesis I",
        "volume": "I",
        "cover": r"Cover: the Thorin density of $\xi(\frac12)/\xi(\frac12+\sqrt u)$ at eight resolutions, sharpening onto the zeros.",
        "subtitle": r"Representations of $\Xi$, Thorin Measures\\[0.2em] and the Reciprocal $\Xi$ Function",
        "subtitle_flat": r"Representations of $\Xi$, Thorin Measures and the Reciprocal $\Xi$ Function",
        "plan": [
            ("PART", r"Foundations: $\Xi$, its representations and the reciprocal $\Xi$ function"),
            (r"Introduction: the Riemann Hypothesis as the positivity of a measure", "Introduction", ["c01", "ch01", "@stub", "c01ex"], "I_01_overture"),
            (r"Analytic representations of the Riemann $\Xi$ function", r"Analytic representations of $\Xi$", ["ch02"], None),
            (r"Probabilistic representations of $\Xi$", None, ["ch19"], None),
            (r"Infinite divisibility, GGC, EGGC and the Bondesson hierarchy", r"ID, GGC, EGGC and the Bondesson hierarchy", ["ch03"], None),
            (r"The reciprocal $\Xi$ function", None, ["c24", "ch10", "ch18", "ht1", "ht4", "ex_sin2more"], None),
            (r"Riemann's continuation and the generalised gamma convolution question", r"Riemann's continuation and GGC", ["nc_riemanncont"], None),
            (r"The reciprocal $\Xi$ function as an extended generalised gamma convolution", r"The reciprocal $\Xi$ function as an EGGC", ["nc_zdist"], None),
            (r"The heat trace and the hierarchy of positivity", r"The heat trace", ["c05", "ch12", "ch14", "ch15", "c05ex"], None),
            (r"Weak limits, HCM and a Bosch type threshold", None, ["c04", "ch11", "boschA", "boschB", "c04ex"], None),
            ("PART", r"The reciprocal zeta function"),
            (r"Thorin type forms of the reciprocal zeta function", r"The reciprocal zeta function", ["c06", "rz1", "rz2"], None),
            (r"Infinite divisibility, necklaces and solved models", None, ["c07", "rz3"], None),
            ("PART", r"Basepoints and the critical strip"),
            (r"The basepoint one and hitting times at the basepoints", r"The basepoint one", ["c09", "basepoint", "bphit"], None),
            (r"The critical strip, Riesz's criterion and the squarefree series", r"The critical strip", ["c10", "strip", "riesz", "sqfree"], None),
            (r"Below the basepoint one: a map of processes, free probability, almost GGCs and the Lambert law", r"Below the basepoint one", ["c11", "sato", "wigfree", "almost", "pakes", "pakes2", "pakes3"], None),
            ("PART", r"Beyond zeta"),
            (r"The Dedekind eta function, $L(s,\chi_{12})$ and Ramanujan's $\Delta$", r"The Dedekind eta function", ["ch16", "bpeta"], None),
            (r"Dirichlet $L$ functions and the Linnik constant", None, ["c21", "gl1", "gl2", "gl3", "gl4", "gl5", "gl6", "gl7", "gl8"], None),
            (r"Completeness of the secant dilates and the M\"obius problem", r"Completeness of the secant dilates", ["c29", "lv1", "lv2"], None),
            ("PART", r"Probability, hitting times and statistical mechanics"),
            (r"Hitting times, Krein strings and the P\'olya Bessel clock", r"Hitting times and Krein strings", ["c08", "biane", "c08ex"], None),
            (r"The P\'olya density as a maximum and as an extreme value law", r"Maxima and extreme values", ["c14", "maxwald", "maxwald2", "maxwald3", "ev1"], None),
            (r"Ferromagnetic random variables", None, ["c13", "fm1", "fm2", "fm2m", "fm3", "fm5", "fm5k", "fm6"], None),
            (r"Cauchy mixtures, the reverse sine squared lemma and a two spin Lee and Yang form of RH", r"Cauchy mixtures and Lee and Yang", ["nc_cauchy"], None),
            (r"P\'olya and the zeros of entire functions", None, ["c16", "polya1", "polya2", "polya3", "polya4", "polya5", "polya6", "polya7", "polya8"], None),
            (r"Haseo Ki, derivatives and the Keiper Li coefficients", r"Ki and the Keiper Li coefficients", ["c23", "ki1", "ki2", "ki3", "ch20"], None),
            ("PART", r"Synthesis"),
            (r"A route to the Riemann Hypothesis", None, ["nc_route"], None),
            (r"Conclusions and open problems of the theory", r"Conclusions and open problems", ["c12", "path", "ch22"], None),
        ],
    },
    "2": {
        "main": "main_vol2",
        "other": "main_vol1",
        "title": "Computing the Riemann Hypothesis II",
        "volume": "II",
        "cover": r"Cover: the backward Cauchy flow, from the level one density down to the zeros.",
        "subtitle": r"Cauchy Deconvolution\\[0.2em] and Certified Computation",
        "subtitle_flat": r"Cauchy Deconvolution and Certified Computation",
        "plan": [
            ("PART", r"What computation can decide"),
            (r"What computation can and cannot decide", None, ["ch06", "ch13", "appB"], "II_01_overture"),
            (r"Evaluating the $\Xi$ function in ball arithmetic", r"Evaluating $\Xi$", ["c02", "ch04", "ch05", "ch05b", "c02ex"], None),
            (r"The Thorin measure at finite resolution", None, ["c03", "ch07", "ch08", "ch09", "c03ex"], None),
            (r"The de Bruijn Newman flow computed from the P\'olya density", r"The de Bruijn Newman flow", ["ch21"], None),
            (r"The barrier method in Thorin form", None, ["nc_barrier"], None),
            (r"The Thorin phase diagram", None, ["c15", "cp1", "cp2", "cp3", "cp4", "cp5", "cp6", "c15ex"], None),
            ("PART", r"Cauchy deconvolution"),
            (r"The Cauchy deconvolution operator", None, ["c33", "cd1"], None),
            (r"Three families of generalised gamma convolutions and the deconvolution obstruction", r"Three families and deconvolution", ["c19", "tf1", "tf2", "tf3", "tf4"], None),
            (r"Deconvolution as a damped Weil criterion", None, ["c20", "dc1", "dc2", "dc3"], None),
            (r"The level one measure and the backward Cauchy flow", r"Level one and the backward flow", ["c30", "lo1", "lo2", "ch17"], None),
            (r"Complexity and acceleration of the backward flow", r"Complexity and acceleration", ["c31", "cx1", "cx2", "cx3"], None),
            (r"The explicit formula as a deconvolution", None, ["c32", "es1", "es2"], None),
            (r"From level one to the critical line: particle filters and nested sampling", r"Cauchy deconvolution to the line", ["c34", "cf1", "cf2", "cf3", "cf4"], None),
            (r"Balazard, Saias and Yor: Monte Carlo along vertical lines", r"Balazard, Saias and Yor", ["c35", "bsy1", "bsy2", "bsy3"], None),
            ("PART", r"Weil positivity in windows"),
            (r"Weil positivity and Yoshida's Hermitian forms: certificates and barriers", r"Weil positivity and Yoshida's forms", ["c26", "wt1", "wt2", "wt3", "c17", "yo1", "yo2", "yo3", "yo4"], None),
            (r"Zhu's compact windows and the Landau Widom horizon", r"Zhu's compact windows", ["c18", "zhu1", "zhu2", "zhu3", "zhu4"], None),
            ("PART", r"Approximants and finite tests"),
            (r"Truncations: the Biane Pitman Yor lattice, the beta lattice and Hermite coefficients", r"Truncations and Hermite coefficients", ["rz4", "c28", "tr1", "tr2"], None),
            (r"Friable approximants of $\Xi$ after Miyatake and the cosh ladder", r"Friable approximants after Miyatake", ["c22", "mi1", "mi2", "mi3", "mi4", "mi5"], None),
            (r"Ferromagnetic and extreme value approximants with certified real zeros", r"Ferromagnetic and extreme value approximants", ["fm4", "c25", "ev2", "ev3", "ev4"], None),
            (r"The Comtet Majumdar family: damped kernels and a solved model", r"The Comtet Majumdar family", ["c27", "cg1", "cg2", "cg3"], None),
            (r"Approximating clocks: Jensen, Euler product and boundary hitting times", r"Approximating clocks", ["ht2", "ht3", "ht5"], None),
            ("PART", r"Conditional results and synthesis"),
            (r"The quasi Riemann hypothesis and the Thorin programme", r"The quasi Riemann hypothesis", ["nc_qrh"], None),
            (r"Conditional certificates under a quasi Riemann hypothesis at $7/8$", r"Conditional certificates", ["polya9"], None),
            (r"What would close the gap: open problems of the computational programme", r"Open problems of the computation", [], "II_22_open_problems"),
            ("APPENDIX",),
            (r"The P\'olya evaluators", None, ["appA"], None),
        ],
    },
}

# Statements to lift into Volume I (label -> destination chapter of Volume I). The script cannot
# move prose; it lists these in the report with the file and line of each definition.
LIFT = {
    "ch:dbn": "I.2", "sec:ch21:flow": "I.2", "rz:thm:dbn": "I.2", "prop:tilt": "I.2",
    "ch:polya": "I.2", "sec:ch21:hankel": "I.2 (statement)", "sec:ch5:turing": "I.2 (statement)",
    "eq:split": "I.5", "prop:zerodensity": "I.5", "prop:signature": "I.5", "tf:thm:tilt": "I.4",
    "thm:ch13:perturb": "I.6 (statement)", "rz:sec:finite": "I.9 (statement)",
    "rz:thm:kappaexact": "I.9 (statement)", "ht:thm:32": "I.16", "ht:sec:jensen": "I.16",
    "fm:thm:bessel": "I.18", "fm:sec:simple": "I.18",
}

STUB_OUTLINES = {
    "I_01_overture": ["A map of the two volumes", "Status conventions", "Dependencies between the chapters"],
    "I_02_analytic": ["The completed zeta function and Hadamard's product", "Riemann's theta integral and P\\'olya's kernel",
                       "The de Bruijn Newman family", "A table of normalisations", "Hermite, friable and Miyatake representations",
                       "The Poisson smoothing of the zeros", "The L\\'evy form, the Euler product and the explicit formula"],
    "I_03_probabilistic": ["The P\\'olya variable and the clock", "Maxima and extreme values", "Spins, strings and hitting times",
                           "Normal scale mixtures and the Meixner law", "A table of representations"],
    "I_04_toolkit": ["Infinite divisibility and self decomposability", "Generalised gamma convolutions",
                     "The Bondesson class and hyperbolic complete monotonicity", "Extended generalised gamma convolutions",
                     "Van Dantzig pairs, Wald couples, z and Meixner laws", "P\\'olya frequency functions and the Laguerre P\\'olya class",
                     "Kent's theorem and Krein strings"],
    "I_05_reciprocal_xi": ["The clock, its existence and its infinite divisibility", "The $2\\sin^2$ lemma, stated once",
                           "Kent's formula and the spectral expansion"],
    "I_21_open_problems": ["Open problems of Volume I, merged and numbered"],
    "II_01_overture": ["What a finite computation can and cannot decide"],
    "II_06_operator_opening": ["The operator and its place in the volume"],
    "II_21_qrh_quarantine": ["The announced half plane and what rests on it"],
    "II_22_open_problems": ["Open problems of Volume II, merged and numbered"],
}

# ---------------------------------------------------------------------------


def strip_comments(s):
    return re.sub(r"(?m)(?<!\\)%.*$", "", s)


def match_brace(s, i):
    """s[i] == '{' ; return index just after the matching '}'."""
    depth = 0
    for j in range(i, len(s)):
        c = s[j]
        if c == "\\":
            continue
        if c == "{" and (j == 0 or s[j - 1] != "\\"):
            depth += 1
        elif c == "}" and s[j - 1] != "\\":
            depth -= 1
            if depth == 0:
                return j + 1
    raise ValueError("unbalanced braces")


def take_command(s, name):
    """Find \\name[opt]{arg} (or \\name*{...}); return (start, end, arg) or None."""
    m = re.search(r"\\" + name + r"\*?\s*(\[[^\]]*\])?\s*\{", s)
    if not m:
        return None
    end = match_brace(s, m.end() - 1)
    return m.start(), end, s[m.end():end - 1]


def pdfsafe(title):
    """Wrap inline math in \\texorpdfstring so PDF bookmarks keep a readable title."""
    def plain(m):
        x = m.group(1)
        t = re.sub(r"\\([A-Za-z]+)", r"\1", x).replace("{", "").replace("}", "").replace("_", "")
        return "\\texorpdfstring{$%s$}{%s}" % (x, t)
    return re.sub(r"\$([^$]+)\$", plain, title)


def demote_header(text, first):
    """Old chapter opening: drop \\chapter and \\chaptermark; if not the first opening of the new
    chapter, keep the old title as an unnumbered section."""
    out = text
    t = take_command(out, "chapter")
    title = None
    if t:
        a, b, title = t
        repl = "" if first else "\\section*{%s}" % title
        out = out[:a] + repl + out[b:]
    t = take_command(out, "chaptermark")
    if t:
        out = out[:t[0]] + out[t[1]:]
    return out, title


def promote_appendix(text):
    """A file with its own \\chapter (appA, appB) placed inside a chapter becomes a section."""
    t = take_command(text, "chapter")
    if not t:
        return text
    a, b, title = t
    return text[:a] + "\\section{%s}" % title + text[b:]


SRC_FOR_HOIST = {}


def split_extra_chapters(text):
    """Return (first part, [later \\chapter blocks]) for a file holding several chapters, such as
    appB.tex, which carries both the short history and the Notation appendix."""
    starts = [m.start() for m in re.finditer(r"(?m)^\\chapter\*?[\[{]", text)]
    if len(starts) < 2:
        return text, []
    return text[:starts[1]], [text[a:b] for a, b in zip(starts[1:], starts[2:] + [len(text)])]


EXTRAS = []  # (source file, index, text) of chapters split off from section files; set appendices of both volumes


def load_sources(root):
    main = open(os.path.join(root, "main.tex"), encoding="utf8").read()
    files = re.findall(r"\\input\{ch/(\w+)\}", main)
    src = {}
    for f in set(files):
        src[f] = open(os.path.join(root, "ch", f + ".tex"), encoding="utf8").read()
    # new chapters ch/nc_<key>.tex are written for the two volumes only; load those listed as checked
    # in vol/nc_ready.txt (one file name per line), so that a chapter still being written is never built
    ready_path = os.path.join(root, "vol", "nc_ready.txt")
    ready = set(open(ready_path).read().split()) if os.path.exists(ready_path) else set()
    for V in VOLUMES.values():
        for item in V["plan"]:
            if item[0] in ("PART", "APPENDIX"):
                continue
            for f in item[2]:
                if f.startswith("nc_") and f in ready and os.path.exists(os.path.join(root, "ch", f + ".tex")):
                    src[f] = open(os.path.join(root, "ch", f + ".tex"), encoding="utf8").read()
                    files.append(f)
                    nb = os.path.join(root, "ch", f + "_bib.tex")
                    if os.path.exists(nb) and "bib" in src:
                        src["bib"] = src["bib"].replace("\\end{thebibliography}",
                                                        open(nb, encoding="utf8").read().rstrip() + "\n\\end{thebibliography}")
    for V in VOLUMES.values():
        V["plan"] = [item for item in V["plan"] if item[0] in ("PART", "APPENDIX")
                     or not all(f.startswith("nc_") and f not in src for f in item[2] if f != "@stub") or not item[2]]
    SRC_FOR_HOIST.update(src)
    return main, files, src


def check_plan(files, src):
    assigned = collections.Counter()
    for V in VOLUMES.values():
        for item in V["plan"]:
            if item[0] in ("PART", "APPENDIX"):
                continue
            for f in item[2]:
                if f != "@stub":
                    assigned[f] += 1
    problems = []
    special = {"preface", "bib"}
    for f in files:
        if f in special:
            continue
        if assigned[f] == 0:
            problems.append("not assigned: %s" % f)
        if assigned[f] > 1:
            problems.append("assigned %d times: %s" % (assigned[f], f))
    for f in assigned:
        if f not in src:
            problems.append("in the plan but not input by main.tex: %s" % f)
    labels = collections.Counter()
    where = {}
    for f, s in src.items():
        for l in re.findall(r"\\label\{([^}]+)\}", strip_comments(s)):
            labels[l] += 1
            where[l] = f
    dup = [l for l, n in labels.items() if n > 1]
    if dup:
        problems.append("labels defined twice (xr needs unique labels): " + ", ".join(sorted(dup)[:20]))
    return problems, where


def write_volume(root, key, V, src, vol_of):
    os.makedirs(os.path.join(root, "vol", "new"), exist_ok=True)
    body = []
    n = 0
    for item in V["plan"]:
        if item[0] == "PART":
            body.append("\\part{%s}" % pdfsafe(item[1]))
            continue
        if item[0] == "APPENDIX":
            body.append("\\appendix")
            continue
        title, short, fs, stub = item
        n += 1
        name = "v%s_ch%02d" % (key, n)
        lines = ["%% Volume %s, Chapter %d, generated by split_volumes.py; edit PLAN, not this file" % (key, n)]
        lines.append("\\chapter{%s}\\label{vol%s:ch%d}" % (pdfsafe(title), key, n))
        if short:
            lines.append("\\chaptermark{%s}" % short)
        first = True
        if stub and "@stub" not in fs:
            lines.append("\\input{vol/new/%s}" % stub)
            stubpath = os.path.join(root, "vol", "new", stub + ".tex")
            if not os.path.exists(stubpath):
                secs = STUB_OUTLINES.get(stub, ["Draft"])
                with open(stubpath, "w", encoding="utf8") as fh:
                    fh.write("%% New writing for %s. This file is never overwritten by split_volumes.py.\n" % title)
                    fh.write("% Outline from the note rh_two_volumes.pdf; replace the placeholders.\n")
                    for sct in secs:
                        fh.write("\\section{%s}\n%% TODO\n\n" % sct)
        for f in fs:
            if f == "@stub":
                lines.append("\\input{vol/new/%s}" % stub)
                continue
            if re.fullmatch(r"c\d\d", f):
                text, _ = demote_header(src[f], first)
                hdr = os.path.join("vol", "hdr_%s.tex" % f)
                with open(os.path.join(root, hdr), "w", encoding="utf8") as fh:
                    fh.write("%% Opening of the old chapter %s, its \\chapter command removed; generated\n" % f)
                    fh.write(text)
                lines.append("\\input{vol/hdr_%s}" % f)
                first = False
            elif take_command(strip_comments(src[f]), "chapter"):
                first, extra = split_extra_chapters(src[f])
                for k, blk in enumerate(extra):
                    if (f, k) not in [(a, b) for a, b, _ in EXTRAS]:
                        EXTRAS.append((f, k, blk))
                text = promote_appendix(first)
                with open(os.path.join(root, "vol", "sec_%s.tex" % f), "w", encoding="utf8") as fh:
                    fh.write("%% %s with its \\chapter turned into a section; generated\n" % f)
                    fh.write(text)
                lines.append("\\input{vol/sec_%s}" % f)
            else:
                lines.append("\\input{ch/%s}" % f)
            first = False
        with open(os.path.join(root, "vol", name + ".tex"), "w", encoding="utf8") as fh:
            fh.write("\n".join(lines) + "\n")
        body.append("\\input{vol/%s}" % name)
    return body


def write_bibliography(root, key, V, src):
    """Write vol/bib_v<key>.tex holding only the entries of ch/bib.tex that the volume cites."""
    files = [f for item in V["plan"] if item[0] not in ("PART", "APPENDIX") for f in item[2] if f != "@stub"]
    texts = [src[f] for f in files] + [open(os.path.join(root, "vol", "new", item[3] + ".tex"), encoding="utf8").read()
                                       for item in V["plan"] if item[0] not in ("PART", "APPENDIX") and item[3]
                                       and os.path.exists(os.path.join(root, "vol", "new", item[3] + ".tex"))]
    texts += [blk for _, _, blk in EXTRAS] + [src["preface"]]
    cited = set()
    for t in texts:
        for g in re.findall(r"\\(?:no)?cite[pt]?\*?(?:\[[^\]]*\])?\{([^}]+)\}", strip_comments(t)):
            cited.update(x.strip() for x in g.split(","))
    lines = src["bib"].split("\n")
    head, entries, tail, buf, cur = [], [], [], [], None
    for line in lines:
        if line.startswith("\\bibitem"):
            if cur is not None:
                entries.append(cur)
            k = re.match(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", line).group(1)
            cur = [k, buf + [line]]
            buf = []
        elif line.startswith("\\end{thebibliography}"):
            if cur is not None:
                entries.append(cur)
                cur = None
            tail.append(line)
        elif cur is None and not entries and not line.startswith("%"):
            head.append(line)
        elif line.startswith("%") and (cur is not None or entries):
            buf.append(line)
        elif cur is not None:
            cur[1].append(line)
        else:
            head.append(line)
    kept = [e for e in entries if e[0] in cited]
    out = head + [l for e in kept for l in e[1]] + tail
    with open(os.path.join(root, "vol", "bib_v%s.tex" % key), "w", encoding="utf8") as fh:
        fh.write("%% bibliography of Volume %s: the %d of %d entries of ch/bib.tex it cites; generated\n" % (key, len(kept), len(entries)))
        fh.write("\n".join(out) + "\n")
    return len(kept), len(entries), sorted(cited - set(e[0] for e in entries))


def write_main(root, key, V, main, body):
    pre, rest = main.split("\\begin{document}", 1)
    # xr-hyper must come before hyperref; labels are unique across the book, so no prefix is needed
    pre = pre.replace("\\usepackage[hidelinks]{hyperref}",
                      "\\usepackage{xr-hyper}\n\\usepackage[hidelinks]{hyperref}\n\\externaldocument{%s_xr}" % V["other"])
    # theorem environments and macros that some chapter files define on first use (for example the
    # conjecture environment of c22): the new order may use them earlier, so define them in the preamble
    hoist = []
    for f in sorted(SRC_FOR_HOIST):
        for line in SRC_FOR_HOIST[f].split("\n"):
            t = line.strip()
            if (t.startswith("\\ifdefined") and "\\newtheorem" in t) or t.startswith("\\providecommand"):
                if t not in hoist:
                    hoist.append(t)
    if hoist:
        pre = pre + "%% hoisted from chapter files by split_volumes.py\n" + "\n".join(hoist) + "\n"
    front = rest.split("\\input{ch/preface}", 1)[0]
    # the title pages name the volume as the covers do ("Volume I"); the copyright page names it in the full title
    vol = V["volume"]
    front = front.replace(r"{\Huge\scshape Computing the\\[0.35em] Riemann Hypothesis\par}",
         r"{\Huge\scshape Computing the\\[0.35em] Riemann Hypothesis\par}" "\n" r"\vspace{0.5cm}" "\n" r"{\LARGE\scshape Volume %s\par}" % vol)
    front = front.replace(r"{\LARGE\scshape Computing the Riemann Hypothesis\par}",
         r"{\LARGE\scshape Computing the Riemann Hypothesis\par}\vspace{0.4cm}{\Large\scshape Volume %s\par}" % vol)
    front = front.replace(r"{\Huge\scshape Computing the Riemann Hypothesis\par}",
         r"{\Huge\scshape Computing the Riemann Hypothesis\par}\vspace{0.4cm}{\LARGE\scshape Volume %s\par}" % vol)
    front = front.replace(r"\textit{Computing the Riemann Hypothesis: ", r"\textit{Computing the Riemann Hypothesis, Volume %s: " % vol)
    front = front.replace("First edition, October 2026", "Second edition, in two volumes, November 2026")
    front = front.replace(r"Cover: the Thorin density of $\xi(\frac12)/\xi(\frac12+\sqrt u)$ at eight resolutions, computed by the authors.", V["cover"])
    pre = pre.replace(r"pdftitle={Computing the Riemann Hypothesis}", r"pdftitle={Computing the Riemann Hypothesis, Volume %s}" % vol)
    # the copyright page prints this volume's DOI: bookmeta.tex (bookrepo.py generate) picks it by \BookVolume
    pre = pre.replace(r"\input{bookmeta}", r"\def\BookVolume{%s}\input{bookmeta}" % vol)
    front = front.replace("Thorin Measures, P\\'olya Densities\\\\[0.2em] and Hitting Times", V["subtitle"])
    front = front.replace("Thorin Measures, P\\'olya Densities and Hitting Times", V["subtitle_flat"])
    out = [pre, "\\begin{document}", front,
           "% the preface is shared until each volume has its own (see the note, Section 5.6)",
           "\\input{ch/preface}", "\\tableofcontents", "\\mainmatter"]
    out += body
    # chapters split off section files (the Notation appendix of appB) become appendices of both volumes;
    # the second volume's copy drops its \\label commands so that xr sees each label once
    if EXTRAS:
        if not any(item[0] == "APPENDIX" for item in V["plan"]):
            out.append("\\appendix")
        for f, k, blk in EXTRAS:
            name = "extra_%s_%d_v%s" % (f, k, key)
            text = blk if key == "1" else re.sub(r"\\label\{[^}]*\}", "", blk)
            with open(os.path.join(root, "vol", name + ".tex"), "w", encoding="utf8") as fh:
                fh.write("%% %s, chapter %d of the file, placed as an appendix by split_volumes.py\n" % (f, k + 2))
                fh.write(text)
            out.append("\\input{vol/%s}" % name)
    # the lists of figures and tables go at the back of each volume, before the bibliography
    out += ["\\backmatter", "\\cleardoublepage\\phantomsection\\addcontentsline{toc}{chapter}{List of Figures}\\listoffigures", "\\cleardoublepage\\phantomsection\\addcontentsline{toc}{chapter}{List of Tables}\\listoftables", "\\input{vol/bib_v%s}" % key, "\\end{document}", ""]
    path = os.path.join(root, V["main"] + ".tex")
    with open(path, "w", encoding="utf8") as fh:
        fh.write("%% %s, generated by split_volumes.py from main.tex and PLAN; do not edit by hand\n" % V["main"])
        fh.write("\n".join(out))
    with open(os.path.join(root, V["main"] + "-kdp.tex"), "w", encoding="utf8") as fh:
        fh.write("%% KDP interior of %s, as main-kdp.tex does for main.tex\n\\def\\kdpbuild{}\n\\def\\kdpimprint{Independently published}\n\\input{%s}\n" % (V["main"], V["main"]))
    return path


def page_estimates(root):
    """Pages per source file from the PDF bookmarks, as in the audit; None if pypdf or the PDF is missing."""
    try:
        import pypdf
    except ImportError:
        return None
    pdf = os.path.join(root, "pdf", "computing_rh.pdf")
    if not os.path.exists(pdf):
        return None
    r = pypdf.PdfReader(pdf)
    out = []

    def walk(o, d):
        for x in o:
            if isinstance(x, list):
                walk(x, d + 1)
            else:
                out.append((d, r.get_destination_page_number(x) + 1, x.title))
    walk(r.outline, 0)
    secs = [p for d, p, t in out if d == 1]
    chs = [p for d, p, t in out if d == 0 and t not in ("Preface", "Notation")]
    main = open(os.path.join(root, "main.tex"), encoding="utf8").read()
    files = [f for f in re.findall(r"\\input\{ch/(\w+)\}", main) if f not in ("preface", "bib")]
    start = {}
    i = 0
    ci = 0
    for f in files:
        s = strip_comments(open(os.path.join(root, "ch", f + ".tex"), encoding="utf8").read())
        nsec = len(re.findall(r"\\section(?!\*)", s))
        if re.fullmatch(r"c\d\d|app[A-Z]", f) and ci < len(chs):
            start[f] = chs[ci]
            ci += 1
        elif nsec:
            start[f] = secs[i]
        i += nsec
    pages = {}
    known = [f for f in files if f in start]
    for k, f in enumerate(known):
        nxt = start[known[k + 1]] if k + 1 < len(known) else len(r.pages) - 34
        pages[f] = max(0, nxt - start[f])
    for f in files:
        pages.setdefault(f, 0)
    return pages


def report(root, src, where, vol_of, chap_of):
    lines = ["# Split report", ""]
    pages = page_estimates(root)
    for key, V in VOLUMES.items():
        lines.append("## Volume %s: %s" % ("I" if key == "1" else "II", V["subtitle_flat"]))
        n = 0
        total = 0
        for item in V["plan"]:
            if item[0] == "PART":
                lines.append("")
                lines.append("Part: " + item[1])
                continue
            if item[0] == "APPENDIX":
                lines.append("")
                lines.append("Appendices")
                continue
            n += 1
            p = sum(pages.get(f, 0) for f in item[2] if f != "@stub") if pages else None
            total += p or 0
            lines.append("%d. %s  [%s]%s" % (n, item[0], ", ".join(x for x in item[2] if x != "@stub") or "new",
                                            "  about %d pages now" % p if pages else ""))
        if pages:
            lines.append("")
            lines.append("Existing material: about %d pages, before new writing, lifts and trims." % total)
        lines.append("")
    refs_I_to_II = []
    for f, s in src.items():
        if vol_of.get(f) != "1":
            continue
        for ln, line in enumerate(s.split("\n"), 1):
            if line.lstrip().startswith("%"):
                continue
            for grp in re.findall(r"\\(?:eq|auto|c|C|page)?ref\{([^}]+)\}", line):
                for l in grp.split(","):
                    l = l.strip()
                    g = where.get(l)
                    if g and vol_of.get(g) == "2":
                        refs_I_to_II.append((f, ln, l, g))
    lines.append("## References from Volume I into Volume II (%d)" % len(refs_I_to_II))
    lines.append("")
    lines.append("xr resolves all of them, so both volumes build; each one should still become either a pointer"
                 " to a chapter of Volume II or a lifted statement. Lifted labels are marked LIFT.")
    lines.append("")
    by_label = collections.defaultdict(list)
    for f, ln, l, g in refs_I_to_II:
        by_label[(l, g)].append("%s.tex:%d" % (f, ln))
    for (l, g), locs in sorted(by_label.items(), key=lambda x: -len(x[1])):
        tag = "LIFT to %s" % LIFT[l] if l in LIFT else ""
        lines.append("- `%s` (defined in %s, %s) %s: %s" % (l, g, chap_of.get(g, "?"), tag, ", ".join(locs)))
    lines.append("")
    lines.append("## Labels to lift into Volume I")
    lines.append("")
    for l, d in LIFT.items():
        lines.append("- `%s` from ch/%s.tex to %s" % (l, where.get(l, "?"), d))
    lines.append("")
    path = os.path.join(root, "split_report.md")
    with open(path, "w", encoding="utf8") as fh:
        fh.write("\n".join(lines))
    return path, len(refs_I_to_II), len(by_label)


def write_xr(root, V):
    """Copy the book's own labels (all contain a colon) from a volume's .aux into <main>_xr.aux, which
    the other volume reads through xr. Citation keys are left out, so the shared bibliography does not
    produce multiply defined labels."""
    aux = os.path.join(root, V["main"] + ".aux")
    if not os.path.exists(aux):
        return
    keep = [line for line in open(aux, encoding="latin1")
            if line.startswith("\\newlabel{") and ":" in line[10:line.find("}")]]
    with open(os.path.join(root, V["main"] + "_xr.aux"), "w", encoding="latin1") as fh:
        fh.writelines(keep)


def build(root):
    """Three alternating passes so that each volume reads the other's .aux through xr."""
    if os.path.exists(os.path.join(root, "tools", "stamp_version.py")):
        subprocess.run([sys.executable, "tools/stamp_version.py"], cwd=root)
    cmd = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error"]
    for rnd in range(3):
        for V in VOLUMES.values():
            r = subprocess.run(cmd + [V["main"] + ".tex"], cwd=root, capture_output=True, text=True, errors="replace")
            if r.returncode != 0:
                log = r.stdout[-3000:]
                sys.exit("pdflatex failed on %s (pass %d):\n%s" % (V["main"], rnd + 1, log))
            write_xr(root, V)
    for V in VOLUMES.values():
        logp = os.path.join(root, V["main"] + ".log")
        log = open(logp, encoding="latin1").read()
        und = len(re.findall(r"Reference `[^']+' on page \d+ undefined", log))
        mult = len(re.findall(r"Label `[^']+' multiply defined", log))
        pages = re.search(r"Output written on .*?\((\d+) pages", log)
        print("%s.pdf: %s pages, %d undefined references, %d multiply defined labels" % (V["main"], pages.group(1) if pages else "?", und, mult))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="root of the computing_rh repository")
    ap.add_argument("--check", action="store_true", help="only check the plan")
    ap.add_argument("--build", action="store_true", help="typeset both volumes")
    ap.add_argument("--clean", action="store_true", help="remove generated files (keeps vol/new)")
    a = ap.parse_args()
    root = a.root
    if a.clean:
        for V in VOLUMES.values():
            for ext in (".tex", "-kdp.tex"):
                p = os.path.join(root, V["main"] + ext)
                if os.path.exists(p):
                    os.remove(p)
        for f in os.listdir(os.path.join(root, "vol")) if os.path.isdir(os.path.join(root, "vol")) else []:
            if f.endswith(".tex"):
                os.remove(os.path.join(root, "vol", f))
        print("removed generated files; vol/new kept")
        return
    main_tex, files, src = load_sources(root)
    problems, where = check_plan(files, src)
    if problems:
        print("plan problems:")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    print("plan ok: %d source files assigned once each, labels unique" % len(set(files) - {"preface", "bib"}))
    if a.check:
        return
    vol_of, chap_of = {}, {}
    for key, V in VOLUMES.items():
        n = 0
        for item in V["plan"]:
            if item[0] in ("PART", "APPENDIX"):
                continue
            n += 1
            for f in [x for x in item[2] if x != "@stub"]:
                vol_of[f] = key
                chap_of[f] = "%s.%d" % ("I" if key == "1" else "II", n)
    os.makedirs(os.path.join(root, "vol"), exist_ok=True)
    bodies = {key: write_volume(root, key, V, src, vol_of) for key, V in VOLUMES.items()}
    for key, V in VOLUMES.items():
        print("wrote", write_main(root, key, V, main_tex, bodies[key]))
        kept, total, missing = write_bibliography(root, key, V, src)
        print("  bibliography of volume %s: %d of %d entries%s" % (key, kept, total,
              "; cited but not in bib.tex: " + ", ".join(missing) if missing else ""))
    path, nref, nlab = report(root, src, where, vol_of, chap_of)
    print("wrote %s: %d references from Volume I into Volume II, to %d labels" % (path, nref, nlab))
    with open(os.path.join(root, "vol", "plan.json"), "w") as fh:
        json.dump({k: V["plan"] for k, V in VOLUMES.items()}, fh, indent=1)
    if a.build:
        build(root)


if __name__ == "__main__":
    main()
