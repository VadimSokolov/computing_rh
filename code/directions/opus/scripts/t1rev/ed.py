import sys
p=sys.argv[1]; t=open(p).read()
R=[("every checked number except three semigroup entries","every checked number except two semigroup entries"),
('- F:136  sentence on pairing is awkward.\n',''),
('- F:335  "certificate" family: see M6; positivity only "on the grid" is correctly said.\n',''),
("- Table: all entries rerun by own code except at finer script step.","- Table: all entries rerun with own code; the semigroup rows at step 0.05, not the note's 0.025."),
("Proof: log(1+z^2/w^2) = log(w+iz) + log(w-iz) - 2 log w. Here w+iz is zero free in the upper half plane (its zero -iz has negative imaginary part), and w-iz in the lower one. The kernel P_a(w-gamma) evaluates the first at gamma+ia and the second at gamma-ia, since the combination decays like 1/w^2 and the singularity at 0 is log integrable. Split 2 log w = log_+ w + log_- w the same way.","Proof: on the real line log(1+z^2/w^2) = g_+(w) + g_-(w), with g_+ = log(w+iz) - log_+ w analytic in the upper half plane and g_- = log(w-iz) - log_- w analytic in the lower one, where log_+ and log_- are the boundary values of log from each side (their sum is 2 log|w|). Both decay like 1/|w| and have log integrable singularities at 0, so the Poisson kernel evaluates g_+ at gamma+ia and g_- at gamma-ia, which gives the right side."),
("The same three point law appears in ch/nc_cauchy.tex:200.","The square of this transform appears in ch/nc_cauchy.tex:200."),
("Every case tried stays negative to 1e6. The trend, about 0.37 per decade at most, extrapolates to positive near 3e7.","Every case tried stays negative to 1e6. The mean density gains 0.37 per decade (observed 0.31 and 0.10 here), which puts a sign change somewhere beyond 1e7, a heuristic extrapolation not computed."),
]
for a,b in R:
    assert t.count(a)==1,a
    t=t.replace(a,b)
open(p,'w').write(t)
