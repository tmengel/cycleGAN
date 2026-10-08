#!/bin/bash
# Runs one pass 1 job: GEANT4 on <nevents> events of a HepMC file, starting after <skip>.
# This is what each pass 1 condor job executes; it also works interactively.
#
#   scripts/run_pass1.sh <sample> <runnumber> <segment> <hepmc> <skip> <nevents> <outdir> [beam]
#
# Writes <outdir>/G4Hits/G4Hits_<sample>-<run>-<segment>.root. beam = fixed (default: vertex at
# 0,0,0), AuAu, pp or pp_zeroangle (see macros/Fun4All_G4_pass1.C). The random seed is fixed
# by (run, segment).

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

[[ $# -eq 7 || $# -eq 8 ]] || die "usage: $0 <sample> <runnumber> <segment> <hepmc> <skip> <nevents> <outdir> [beam]"
SAMPLE=$1; RUN=$2; SEG=$3; HEPMC=$4; SKIP=$5; NEV=$6; OUTDIR=$7; BEAM=${8:-fixed}
check_run_segment "$RUN" "$SEG"
[[ -f "$HEPMC" ]] || die "no such HepMC file: $HEPMC"
SEED=$(pass1_seed "$RUN" "$SEG")
DST=$(f4a_name G4Hits "$SAMPLE" "$RUN" "$SEG")
CDBTAG=${F4A_CDBTAG:-MDC2_ana.435}

echo ">> pass1 $SAMPLE run $RUN segment $SEG on $(hostname) at $(date)"
echo ">> $HEPMC events [$SKIP, $((SKIP + NEV))), beam $BEAM, seed $SEED, $F4A_RELEASE, CDB $CDBTAG"

enter_workdir
start=$SECONDS
set +e
root.exe -q -b "$F4A_ROOT/macros/Fun4All_G4_pass1.C($NEV,$SKIP,\"$HEPMC\",\"$DST\",$RUN,$SEG,$SEED,\"$BEAM\",\"$CDBTAG\")"
rc=$?
set -e
echo ">> root exit code $rc after $((SECONDS - start)) s"
[[ $rc -eq 0 ]] || die "pass1 root exited $rc"
[[ -s "$DST" ]] || die "pass1 did not write $DST"

# the DST must hold exactly the requested events (fewer = the HepMC file ran out)
nev=$(tree_entries "$DST" T)
echo ">> $DST: $nev events"
[[ "$nev" -eq "$NEV" ]] || die "$DST has $nev events, expected $NEV"

ship "$DST" "$OUTDIR/G4Hits"
echo ">> wrote $OUTDIR/G4Hits/$DST"
