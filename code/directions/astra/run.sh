#!/bin/zsh
cd "$(dirname "$0")"
for id in A B C D2 D5 D10 D1; do
  ( cd $id && codex exec -m gpt-6-astra -c model_reasoning_effort=max -s workspace-write --skip-git-repo-check --add-dir ../hilbert6 "$(cat prompt.txt)" > run.log 2>&1; echo "done $id rc=$?" ) &
done
wait
