P='./misc/research/directions/opus/T4.md'
r=open(P).read()
rep=[
("(mpmath at 110 digits for Section 8 of the small ball)","(mpmath at 110 digits for the small ball)"),
("4e20 primes at T = 3e12","3e20 primes at T = 3e12"),
("at 3e12 it needs 4e20 primes","at 3e12 it needs 3e20 primes"),
("| 150 | 3.2 | 6 | 400 | 80 |","| 150 | 3.2 | 6 | 400 | 70 |"),
("| 1e6 | 12.0 | 23 | 8e9 | 3.5e8 |","| 1e6 | 12.0 | 23 | 8e9 | 3.4e8 |"),
("| 3e12 | 27.0 | 51 | 2e22 | 4e20 |","| 3e12 | 26.9 | 51 | 1.6e22 | 3e20 |"),
("<= C(delta)/r^2 for L >= 1, uniformly in a in (0, 2 delta); and |phi^(r)| <= 2/a is not needed for large r. Each term tends to the Cauchy transform. Dominated convergence over the zeros (sum gamma^-2 < inf) gives","<= C(delta)/r^2 for L >= 1, uniformly in a in (0, 2 delta). Each term tends to the Cauchy transform, and only finitely many zeros have |w - gamma| < 1. Dominated convergence over the others (sum gamma^-2 < inf) gives"),
("- three_computations.tex:132  \"first 58 zeros\" is right; \"the first 58\" at [5,160] only.\n",""),
("It is below H(0) only when L is below about (2/(1-2 delta)) log(...), whereas resolving zeros near T needs L about 2 log(T/2 pi). The two are compatible only for 1/2 - delta of order log log T/log T, where the kernel is too wide to resolve.","It is below H(0) only when e^{(1/2-delta)L} is below about L^2 d^2. Resolving zeros near T needs L about 2 log(T/2 pi), so L d is of order one for neighbours and the two are compatible only for 1/2 - delta of order 1/log T, where the kernel, of width about 1/2, is too wide to resolve."),
("In units of Re xi'/xi(1+iw), with exp kernel, the infimum is at most 1/(1+u^2) per pair against 1/(1/4+u^2) on the line. The integrated density gives background (1/4) log(T/2 pi) after subtraction instead of (1/2) log(T/2 pi): the slope halves, the ceiling rises from about 2 pi e^8 to about 2 pi e^16, S(T) adds error,","In units of Re xi'/xi(1+iw), with exp kernel, the infimum is at most 1/(1+u^2) per pair against 1/(1/4+u^2) on the line. The integrated density leaves a background of at least (1/4) log(T/2 pi) after subtraction instead of (1/2) log(T/2 pi): the slope at best halves, the ceiling rises at best from about 2 pi e^8 = 1.9e4 to about 2 pi e^16 = 5.6e7, S(T) adds error,"),
("side lobes of height 2.9 and 3.0 at 239.84 and 240.98","side lobes of height 3.0 and 2.9 at 239.84 and 240.98"),
("Both report: a failure is a real certificate of off line zeros, and a pass says little.","Both reports agree: a failure is a proof of off line zeros, and a pass says little."),
("FIX for :214: \"I_\\alpha(\\varepsilon)=x_\\varepsilon^2\\varepsilon-\\frac{x_\\varepsilon}2+\\frac{\\alpha+3}2\\log x_\\varepsilon+c_\\alpha+o(1),\" with c_alpha as above. This is a mechanical change; the TEXT is \"\\frac{\\alpha+3}2\\log x_\\varepsilon+O(1)\".","TEXT (:214): \"I_\\alpha(\\varepsilon)=x_\\varepsilon^2\\varepsilon-\\frac{x_\\varepsilon}2+\\frac{\\alpha+3}2\\log x_\\varepsilon+O(1).\"\nFIX (:214): \"I_\\alpha(\\varepsilon)=x_\\varepsilon^2\\varepsilon-\\frac{x_\\varepsilon}2+\\frac{\\alpha+3}2\\log x_\\varepsilon+c_\\alpha+o(1),\" with c_alpha defined as above. ALSO AT :223 (same term, proof), :242, :259."),
("Delta = 0.5 (t4_dh.py): first at n about 1e4 the first peaks above 2 are 85.70 (2.76) and 114.16 (2.46), but on line clusters pass 2 at 145.44, 149.38 and above:","Delta = 0.5 (t4_dh.py): at n = 1e4 the lowest maxima above 2 are 85.70 (2.76) and 114.16 (2.46), but on line clusters pass 2 at 145.44, 149.38 and above:"),
("and past L about 30 it is dominated (Proposition C)","and by L about 30 it is dominated (Proposition C, kernel estimate, not run on the full sum)"),
("So at one fixed L = 14.51, usable to about T = 1e4.","So at one fixed L = 14.51, usable to about T = 1e4, close pairs aside (sup 2.03 at L = 14 on [4000,4300])."),
]
for a,b in rep:
    assert r.count(a)>=1,a[:50]
    r=r.replace(a,b)
open(P,'w').write(r)
t=open('./misc/reviews/nick-ferro-thorin/ferromagnetic_thorin_for_vadim/three_computations.tex').read()
print(t.count("I_\\alpha(\\varepsilon)=x_\\varepsilon^2\\varepsilon-\\frac{x_\\varepsilon}2+\\frac{\\alpha+3}2\\log x_\\varepsilon+O(1)."))
print(r.count('### [M'), r.count('—'), len(r.splitlines()))
