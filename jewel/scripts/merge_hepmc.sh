#!/bin/bash
# Concatenates standalone HepMC2 files (e.g. a production's raw/ or filtered/
# chunks) into one valid HepMC2 file, optionally stopping after N events.
# Inputs are taken in version-sorted order and each must have its
# END_EVENT_LISTING footer.
#
#   scripts/merge_hepmc.sh <out.hepmc> [--max N] <in1.hepmc> [in2.hepmc ...]
#
# Example:
#   scripts/merge_hepmc.sh vacuum_10k.hepmc --max 10000 \
#     /sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc/test_v1/vacuum/filtered/*.hepmc

set -euo pipefail

[[ $# -ge 2 ]] || { echo "usage: $0 <out.hepmc> [--max N] <in.hepmc>..." >&2; exit 1; }
OUT=$1; shift
MAX=0
if [[ "$1" == "--max" ]]; then MAX=$2; shift 2; fi
[[ $# -ge 1 ]] || { echo "ERROR: no input files" >&2; exit 1; }

tmp="$OUT.part"; cnt="$OUT.count"
trap 'rm -f "$tmp" "$cnt"' EXIT
printf '\nHepMC::Version 2.06.05\nHepMC::IO_GenEvent-START_EVENT_LISTING\n' > "$tmp"

total=0; nfiles=0
for f in $(printf '%s\n' "$@" | sort -V); do
  tail -c 4096 "$f" | grep -q "END_EVENT_LISTING" || { echo "ERROR: $f has no footer (truncated?)" >&2; exit 1; }
  # Copy event lines only (drop HepMC:: header/footer and blank lines),
  # stopping before the event that would exceed the remaining budget.
  awk -v budget=$(( MAX > 0 ? MAX - total : -1 )) -v cnt="$cnt" '
        /^HepMC::/ || NF == 0 { next }
        /^E / { if (budget >= 0 && n >= budget) exit; n++ }
        { print }
        END { print n + 0 > cnt }' "$f" >> "$tmp"
  n=$(cat "$cnt")
  total=$((total + n)); nfiles=$((nfiles + 1))
  (( MAX > 0 && total >= MAX )) && break
done

printf 'HepMC::IO_GenEvent-END_EVENT_LISTING\n\n' >> "$tmp"
mv "$tmp" "$OUT"; rm -f "$cnt"
trap - EXIT
echo ">> wrote $OUT: $total events from $nfiles files"
(( MAX == 0 || total == MAX )) || echo "WARNING: only $total events available (asked for $MAX)"
