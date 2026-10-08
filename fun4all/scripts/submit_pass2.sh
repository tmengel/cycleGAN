#!/bin/bash
# Prepares (and by default submits) pass 2 on condor: one job per pass 1 DST of a sample.
#
#   scripts/submit_pass2.sh [--only-missing] [--skip-submitted] [--max-jobs N] [--dry-run] <sample> <tag>
#
#   --only-missing    queue only DSTs whose trees or images file is missing
#   --skip-submitted  also skip DSTs already submitted (condor/pass2_submitted.list), e.g. still queued
#   --max-jobs N      queue at most N jobs (the scheduler allows 100,000 queued jobs per user)
#   --dry-run       write the job file but don't submit
#
# Reads $F4A_OUTPUT_ROOT/<tag>/<sample>/G4Hits/*.root, writes trees/ and images/ next to it.

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

ONLYMISSING=0; SKIPSUB=0; MAXJOBS=0; DRYRUN=0
while [[ $# -gt 0 && "$1" == --* ]]; do
  case $1 in
    --only-missing) ONLYMISSING=1; shift ;;
    --skip-submitted) SKIPSUB=1; shift ;;
    --max-jobs) MAXJOBS=$2; shift 2 ;;
    --dry-run) DRYRUN=1; shift ;;
    *) die "unknown option $1" ;;
  esac
done
[[ $# -eq 2 ]] || die "usage: $0 [--only-missing] [--skip-submitted] [--max-jobs N] [--dry-run] <sample> <tag>"
SAMPLE=$1; TAG=$2
OUT="$F4A_OUTPUT_ROOT/$TAG/$SAMPLE"
[[ -d "$OUT/G4Hits" ]] || die "no pass 1 output in $OUT/G4Hits"
[[ -x "$F4A_INSTALL/lib/libjeweltreewriter.so" || -f "$F4A_INSTALL/lib/libjeweltreewriter.so" ]] || die "JewelTreeWriter not built -- run ./install.sh"

queue="$OUT/condor/pass2_queue.list"
mkdir -p "$OUT/condor" "$OUT/logs"
: > "$queue"
ndst=0
havetrees=$(mktemp); haveimages=$(mktemp); dsts=$(mktemp)
trap 'rm -f "$havetrees" "$haveimages" "$dsts"' EXIT
dir_names "$OUT/trees" > "$havetrees"
dir_names "$OUT/images" > "$haveimages"
dir_names "$OUT/G4Hits" | grep '^G4Hits_.*\.root$' | sed 's/^G4Hits_//; s/\.root$//' | sort -V > "$dsts"
ndst=$(wc -l < "$dsts")
awk -v only="$ONLYMISSING" -v dir="$OUT/G4Hits" '
    FILENAME == ARGV[1] {t[$1] = 1; next}
    FILENAME == ARGV[2] {i[$1] = 1; next}
    !(only == 1 && (("trees_" $1 ".root") in t) && (("images_" $1 ".root") in i)) {print dir "/G4Hits_" $1 ".root", $1}' \
    "$havetrees" "$haveimages" "$dsts" > "$queue"
submitted="$OUT/condor/pass2_submitted.list"
if [[ $SKIPSUB -eq 1 && -f "$submitted" ]]; then
  awk 'NR == FNR {s[$1] = 1; next} !($2 in s)' "$submitted" "$queue" > "$queue.tmp" && mv "$queue.tmp" "$queue"
fi
if (( MAXJOBS > 0 )); then
  head -n "$MAXJOBS" "$queue" > "$queue.tmp" && mv "$queue.tmp" "$queue"
fi
nq=$(wc -l < "$queue")
(( ndst > 0 )) || die "no G4Hits DSTs in $OUT/G4Hits"

jobfile="$OUT/condor/pass2.job"
cat > "$jobfile" <<EOF
universe            = vanilla
executable          = $F4A_ROOT/scripts/run_pass2.sh
transfer_executable = False
arguments           = \$(dst) $OUT
output              = $OUT/logs/pass2_\$(name).out
error               = $OUT/logs/pass2_\$(name).err
log                 = /tmp/${USER}_f4a_pass2_${TAG}_${SAMPLE}.log
request_memory      = 3000MB
PeriodicHold        = (NumJobStarts >= 1 && JobStatus == 1)
notification        = Never
queue dst, name from $queue
EOF

cat >> "$OUT/production.txt" <<EOF
# ---- pass2 $(date '+%F %T') on $(hostname)$([[ $DRYRUN -eq 1 ]] && echo " (prepared, not submitted)")
$ndst DSTs, $nq queued; rerun one: scripts/run_pass2.sh <G4Hits DST> <outdir>
seeds      pass2 = 1000000000 + run*100000 + segment
release    $F4A_RELEASE ($OFFLINE_MAIN), JewelTreeWriter from $F4A_INSTALL
EOF

echo ">> $SAMPLE / $TAG: $ndst DSTs, $nq to queue"
echo ">> job file: $jobfile"
if [[ $DRYRUN -eq 1 ]]; then
  echo ">> dry run: not submitted"
  exit 0
fi
(( nq > 0 )) || { echo ">> nothing to submit"; exit 0; }
condor_submit "$jobfile" && awk '{print $2}' "$queue" >> "$submitted"
