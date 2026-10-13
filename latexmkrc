# Before every build, write version.tex with the git commit of the source; main.tex prints it and stores it in the PDF information.
# The old stamp is removed first, so that a failed stamp shows the commit as unknown rather than a stale one.
unlink 'version.tex';
system('python3', 'tools/stamp_version.py') == 0 or warn "tools/stamp_version.py failed: the PDF will show the commit as unknown\n";
