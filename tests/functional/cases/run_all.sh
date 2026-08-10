#!/bin/bash
# Submit every functional case (f1..fN), wait for them, and print a PASS/FAIL summary.
# Runs on the login node — it only sbatch-submits, polls, and greps the RESULT lines.
#
#   bash run_all.sh
#
# Exit status is 0 only if every case passed (handy for scripting/CI).
# The per-case SLURM logs (<jobname>_<jobid>.out) land in this dir and are git-ignored.
set -uo pipefail   # NOT -e: we want to report a failing case, not abort on it

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"   # so %x_%j.out logs land here and can be found by job id

# Auto-discover cases (f1_*.sh, f2_*.sh, ...), version-sorted — picks up future fN without editing.
mapfile -t CASES < <(ls f[0-9]*.sh 2>/dev/null | sort -V | sed 's/\.sh$//')
[ "${#CASES[@]}" -gt 0 ] || { echo "no case scripts (f[0-9]*.sh) found in $HERE"; exit 1; }

declare -A JID
echo "submitting ${#CASES[@]} functional cases..."
for c in "${CASES[@]}"; do
  jid=$(sbatch --parsable "$c.sh") || { echo "  FAILED to submit $c"; continue; }
  JID[$c]=$jid
  printf "  %-26s job %s\n" "$c" "$jid"
done

echo "waiting for completion (Ctrl-C just stops watching; jobs keep running)..."
terminal() { case "$1" in COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL) return 0;; *) return 1;; esac; }
for _ in $(seq 1 240); do   # cap ~60 min; the suite normally finishes in ~2 min
  left=0
  for c in "${CASES[@]}"; do
    st=$(sacct -j "${JID[$c]:-x}" --format=State -n 2>/dev/null | head -1 | tr -d ' ')
    terminal "$st" || left=$((left+1))
  done
  [ "$left" -eq 0 ] && break
  printf "\r  still running: %d/%d   " "$left" "${#CASES[@]}"
  sleep 15
done
printf "\r%40s\r" " "

echo
echo "==================== RESULTS ===================="
pass=0; fail=0
for c in "${CASES[@]}"; do
  jid=${JID[$c]:-}
  [ -n "$jid" ] || { printf "  %-26s %-6s [not submitted]\n" "$c" "FAIL"; fail=$((fail+1)); continue; }
  st=$(sacct -j "$jid" --format=State -n 2>/dev/null | head -1 | tr -d ' ')
  out=$(ls -t ./*_"$jid".out 2>/dev/null | head -1)
  res=$(grep -h 'RESULT:' "$out" 2>/dev/null | tail -1)
  if   printf '%s' "$res" | grep -q 'PASS'; then mark="PASS"; pass=$((pass+1))
  elif printf '%s' "$res" | grep -q 'FAIL'; then mark="FAIL"; fail=$((fail+1))
  elif [ "$st" = "COMPLETED" ];              then mark="PASS"; pass=$((pass+1))
  else                                            mark="FAIL"; fail=$((fail+1)); fi
  printf "  %-26s %-6s [%s]  %s\n" "$c" "$mark" "${st:-?}" "${res#RESULT: }"
done
echo "================================================="
echo "  $pass passed, $fail failed  (of ${#CASES[@]})"
[ "$fail" -eq 0 ]
