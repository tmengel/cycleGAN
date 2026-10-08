#!/bin/bash
# Runs one pass 2 job on one pass 1 DST: trees (Fun4All_pass2.C), then images (make_images.C).
# This is what each pass 2 condor job executes; it also works interactively.
#
#   scripts/run_pass2.sh <G4Hits DST> <outdir>
#
# Sample, run number and segment are taken from the DST name
# (G4Hits_<sample>-<run>-<segment>.root). Writes
#   <outdir>/trees/trees_<sample>-<run>-<segment>.root    towertree, jettree, globaltree
#   <outdir>/images/images_<sample>-<run>-<segment>.root  one TH2D per event
# The random seed is fixed by (run, segment).

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

[[ $# -eq 2 ]] || die "usage: $0 <G4Hits DST> <outdir>"
DST=$1; OUTDIR=$2
[[ -f "$DST" ]] || die "no such DST: $DST"
base=$(basename "$DST" .root)
[[ "$base" =~ ^G4Hits_(.+)-([0-9]{10})-([0-9]{6})$ ]] || die "unexpected DST name $base"
SAMPLE=${BASH_REMATCH[1]}; RUN=$((10#${BASH_REMATCH[2]})); SEG=$((10#${BASH_REMATCH[3]}))
check_run_segment "$RUN" "$SEG"
SEED=$(pass2_seed "$RUN" "$SEG")
TREES=$(f4a_name trees "$SAMPLE" "$RUN" "$SEG")
IMAGES=$(f4a_name images "$SAMPLE" "$RUN" "$SEG")

echo ">> pass2 $SAMPLE run $RUN segment $SEG on $(hostname) at $(date): seed $SEED, $F4A_RELEASE"
nin=$(tree_entries "$DST" T)
echo ">> input $DST: $nin events"

enter_workdir
cp "$DST" input.root   # local copy: the tree stage is I/O bound on the shared filesystem

start=$SECONDS
set +e
root.exe -q -b "$F4A_ROOT/macros/Fun4All_pass2.C(\"input.root\",\"$TREES\",$SEED)"
rc=$?
set -e
echo ">> pass2 root exit code $rc after $((SECONDS - start)) s"
[[ $rc -eq 0 ]] || die "pass2 root exited $rc"

# every input event must come out (the writer fills all three trees once per event)
for t in towertree jettree globaltree; do
  n=$(tree_entries "$TREES" "$t")
  [[ "$n" -eq "$nin" ]] || die "$TREES: $t has $n entries, input had $nin"
done
echo ">> $TREES: $nin entries in each tree"

# interpreted on purpose: ACLiC ("+") would have every job write the same .so into macros/
root.exe -q -b "$F4A_ROOT/macros/make_images.C(\"$TREES\",\"$IMAGES\",$SEG)" || die "make_images failed"
nimg=$(root.exe -q -b -l -e 'TFile f("'"$IMAGES"'"); printf("%d\n", f.GetNkeys()); exit(0);' 2>/dev/null | tail -1)
[[ "$nimg" -eq "$nin" ]] || die "$IMAGES has $nimg images, expected $nin"

ship "$TREES" "$OUTDIR/trees"
ship "$IMAGES" "$OUTDIR/images"
echo ">> wrote $OUTDIR/trees/$TREES and $OUTDIR/images/$IMAGES"
