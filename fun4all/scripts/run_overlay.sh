#!/bin/bash
# Runs one overlay job: adds HIJING 0-10% background images to a range of JEWEL image files.
# This is what each overlay condor job executes; it also works interactively.
#
#   scripts/run_overlay.sh <chunk list> <sample> <run> <jewel images dir> <outdir>
#
# <chunk list>: lines "<segment> <nevents> <first global event>" (made by submit_overlay.sh).
# Writes <outdir>/images_<sample>_hijing-<run>-<segment>.root; see macros/overlay_hijing.C.

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

[[ $# -eq 5 ]] || die "usage: $0 <chunk list> <sample> <run> <jewel images dir> <outdir>"
CHUNK=$1; SAMPLE=$2; RUN=$3; JEWEL=$4; OUTDIR=$5
HIJING_DIR=${F4A_HIJING_DIR:-/sphenix/user/yeonjugo/jetML/create_calo_images/macro/type4_run19_hijingAll_noNoise_0to10/histograms/cent0}
COUNTS=$F4A_ROOT/hijing/type4_run19_hijingAll_noNoise_0to10_cent0.counts

echo ">> overlay $SAMPLE run $RUN, $(wc -l < "$CHUNK") segments from $(basename "$CHUNK") on $(hostname) at $(date)"
enter_workdir
mkdir -p out
# interpreted on purpose (ACLiC would make every job write the same .so)
root.exe -q -b "$F4A_ROOT/macros/overlay_hijing.C(\"$CHUNK\",\"$JEWEL\",\"out\",\"$SAMPLE\",$RUN,\"$HIJING_DIR\",\"$COUNTS\")" || die "overlay_hijing failed"

n=0
for f in out/*.root; do
  [[ -e "$f" ]] || continue
  ship "$f" "$OUTDIR"
  n=$((n + 1))
done
echo ">> wrote $n files to $OUTDIR"
