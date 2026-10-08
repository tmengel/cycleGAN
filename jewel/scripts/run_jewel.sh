#!/bin/bash
# Runs one JEWEL job and copies its output to <outdir>. This is what each
# condor job executes, and it also works interactively for a quick test.
#
#   scripts/run_jewel.sh <config> <njob> <nevent> <outdir> [R ptmin etamax]
#
#   config   a directory name under configs/ (vacuum, medium, ...)
#   njob     JEWEL NJOB = random seed (seed = NJOB*1000); unique per job!
#   nevent   events to generate
#   outdir   gets raw/<config>_<njob>.hepmc and logs/<config>_<njob>.log;
#            with the optional filter args also filtered/<config>_<njob>.hepmc
#            (anti-kT R, keep events with a jet pT > ptmin, |eta| < etamax)
#
# Example:
#   scripts/run_jewel.sh vacuum 1 100 /tmp/jewel_test

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

[[ $# -eq 4 || $# -eq 7 ]] || die "usage: $0 <config> <njob> <nevent> <outdir> [R ptmin etamax]"
load_config "$1"
NJOB=$2
NEVENT=$3
OUTDIR=$4
check_njob "$NJOB"
check_tables

name="${CONFIG_NAME}_$(printf '%06d' "$NJOB")"
echo ">> $name on $(hostname) at $(date): $(basename "$JEWEL_EXE"), NEVENT=$NEVENT"

# Work in condor's scratch dir (local disk), or a temp dir interactively.
if [[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "${_CONDOR_SCRATCH_DIR}" ]]; then
  work="$_CONDOR_SCRATCH_DIR"
else
  work=$(mktemp -d)
  trap 'rm -rf "$work"' EXIT
fi
cd "$work"

cp "$CONFIG_DIR"/tables/{xsecs,pdfs,splitint}.dat .
[[ -f "$CONFIG_DIR/medium.params.dat" ]] && cp "$CONFIG_DIR/medium.params.dat" .
write_params params.dat "$NJOB" "$NEVENT" "$name.hepmc" "$name.log"

start=$SECONDS
set +e
"$JEWEL_EXE" params.dat < /dev/null > jewel.stdout 2>&1
ret=$?
set -e
echo ">> JEWEL exit code $ret after $((SECONDS - start)) s"

mkdir -p "$OUTDIR/logs"
cat params.dat "$name.log" > "$OUTDIR/logs/$name.log" 2>/dev/null || true

# A run killed mid-write (eviction, walltime) leaves no footer; never ship it.
if [[ $ret -ne 0 ]] || ! tail -c 4096 "$name.hepmc" 2>/dev/null | grep -q "END_EVENT_LISTING"; then
  tail -30 jewel.stdout
  die "JEWEL failed or $name.hepmc is truncated"
fi
nev=$(grep -c '^E ' "$name.hepmc")
# With a nuclear PDF JEWEL draws the pp/pn/np/nn event counts from Poisson
# distributions, so NEVENT +- sqrt(NEVENT) is normal.
echo ">> generated $nev events (asked for $NEVENT)"
(( nev * 100 >= NEVENT * 90 )) || echo "WARNING: far fewer events than asked for"

if [[ $# -eq 7 ]]; then
  "$JEWEL_ROOT/bin/jetFilter" "$name.hepmc" "$name.filtered.hepmc" "$5" "$6" "$7"
  mkdir -p "$OUTDIR/filtered"
  cp "$name.filtered.hepmc" "$OUTDIR/filtered/$name.hepmc"
  cp "$name.filtered.hepmc.stats" "$OUTDIR/filtered/$name.hepmc.stats"
fi

# Copy to a temp name then rename, so a half-copied file is never mistaken
# for a finished one.
mkdir -p "$OUTDIR/raw"
cp "$name.hepmc" "$OUTDIR/raw/.$name.hepmc.part"
mv "$OUTDIR/raw/.$name.hepmc.part" "$OUTDIR/raw/$name.hepmc"
echo ">> wrote $OUTDIR/raw/$name.hepmc"
