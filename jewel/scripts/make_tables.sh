#!/bin/bash
# Builds the cached integration tables (xsecs.dat, pdfs.dat, splitint.dat) for
# one config, into configs/<config>/tables/. Run once per config, and again
# whenever its params change (run_jewel.sh refuses stale tables).
#
#   scripts/make_tables.sh <config>
#
# The tables depend on the medium: xsecs.dat is tabulated over the Debye-mass
# range [m_D(T_c), m_D(T_max)], so vacuum (m_D = 0), TI = 0.26 and TI = 0.375
# each need their own. Never copy tables between configs.

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

[[ $# -eq 1 ]] || die "usage: $0 <config>"
load_config "$1"

tables="$CONFIG_DIR/tables"
work=$(mktemp -d "$CONFIG_DIR/.tables.XXXXXX")
trap 'rm -rf "$work"' EXIT
cd "$work"

[[ -f "$CONFIG_DIR/medium.params.dat" ]] && cp "$CONFIG_DIR/medium.params.dat" .
# NJOB 1, a single event: JEWEL integrates the tables first, then writes them out.
write_params params.dat 1 1 tables_test.hepmc tables.log

echo ">> integrating tables for '$CONFIG_NAME' with $(basename "$JEWEL_EXE") (can take several minutes)"
start=$SECONDS
"$JEWEL_EXE" params.dat < /dev/null > jewel.stdout 2>&1 || { tail -20 jewel.stdout; die "JEWEL failed"; }
echo ">> done in $((SECONDS - start)) s"

for f in xsecs.dat pdfs.dat splitint.dat; do
  [[ -s $f ]] || die "JEWEL did not write $f (see $work/tables.log)"
done
grep -q "HepMC::IO_GenEvent-END_EVENT_LISTING" tables_test.hepmc || die "test event file incomplete"

rm -rf "$tables"
mkdir -p "$tables"
mv xsecs.dat pdfs.dat splitint.dat tables.log "$tables/"
config_fingerprint > "$tables/fingerprint"

echo ">> wrote $tables"
echo ">> parameters JEWEL actually used (check these -- a value JEWEL can't parse silently keeps its default):"
sed -n '/NEVENT/,/^ *$/p' "$tables/tables.log" | head -60
if grep -q "TAUI" "$tables/tables.log"; then
  echo ">> medium parameters:"
  grep -E "^ *(TAUI|TI|TC|WOODSSAXON|CENTRMIN|CENTRMAX|NF|SIGMANN|MDFACTOR|MDSCALEFAC|A) *=" "$tables/tables.log"
fi
