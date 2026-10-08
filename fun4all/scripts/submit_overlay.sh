#!/bin/bash
# Prepares (and by default submits) the HIJING-overlay images for one sample: every JEWEL image of
# <sample> plus one HIJING 0-10% background image (see macros/overlay_hijing.C, hijing/README.md).
#
#   scripts/submit_overlay.sh [--per-job N] [--max-events N] [--priority P] [--only-missing] [--dry-run]
#                             <sample> <tag> <runnumber>
#
#   --per-job N     JEWEL image files (segments) per condor job (default 100, i.e. ~1000 images)
#   --max-events N  only overlay JEWEL events 0..N-1 (e.g. the first 100k), default: all
#   --priority P    condor job priority among your own jobs (default 0; higher runs first)
#   --only-missing  queue only jobs with a missing output file
#   --dry-run       write the job files but don't submit
#
# JEWEL event e of the sample (events numbered in pass 1 segment order, from
# condor/pass1_jobs.list) gets HIJING event e, so the same index means the same background in
# every sample. No background is reused: events beyond the number of HIJING images (985,791) get
# no overlay image. Only segments whose JEWEL image file exists are queued.
#
# Output: $F4A_OUTPUT_ROOT/<tag>/<sample>_hijing/images/images_<sample>_hijing-<run>-<seg>.root

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

PERJOB=100; MAXEV=0; PRIO=0; ONLYMISSING=0; DRYRUN=0
while [[ $# -gt 0 && "$1" == --* ]]; do
  case $1 in
    --per-job) PERJOB=$2; shift 2 ;;
    --max-events) MAXEV=$2; shift 2 ;;
    --priority) PRIO=$2; shift 2 ;;
    --only-missing) ONLYMISSING=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    *) die "unknown option $1" ;;
  esac
done
[[ $# -eq 3 ]] || die "usage: $0 [--per-job N] [--max-events N] [--priority P] [--only-missing] [--dry-run] <sample> <tag> <runnumber>"
SAMPLE=$1; TAG=$2; RUN=$3
check_run_segment "$RUN" 0

IN="$F4A_OUTPUT_ROOT/$TAG/$SAMPLE"
OUT="$F4A_OUTPUT_ROOT/$TAG/${SAMPLE}_hijing"
P1LIST="$IN/condor/pass1_jobs.list"
[[ -f "$P1LIST" ]] || die "no pass 1 job list $P1LIST"
NHIJING=$(awk '$2 > 0 {s += $2} END {print s + 0}' "$F4A_ROOT/hijing/type4_run19_hijingAll_noNoise_0to10_cent0.counts")
mkdir -p "$OUT"/{condor/chunks,logs,images}

# segment -> (nevents, first global event); keep segments with an existing JEWEL image file and at
# least one event below NHIJING
all="$OUT/condor/overlay_segments.list"
LIMIT=$NHIJING
(( MAXEV > 0 && MAXEV < LIMIT )) && LIMIT=$MAXEV
# a segment straddling LIMIT is written with only its events below LIMIT when LIMIT = NHIJING
# (no reuse); with --max-events, segments are kept whole if they start below LIMIT
awk -v nh="$LIMIT" '{ if (e + 0 < nh) print $1, $4, e + 0; e += $4 }' "$P1LIST" > "$all"
: > "$OUT/condor/overlay_queue.list"
rm -f "$OUT"/condor/chunks/chunk_*.list
k=0; nseg=0; nmissing_in=0
chunk=""
# existing input / output segments, from one directory read each
havein=$(mktemp); haveout=$(mktemp)
trap 'rm -f "$havein" "$haveout"' EXIT
dir_names "$IN/images" | awk -v pre="images_${SAMPLE}-$(printf '%010d' "$RUN")-" 'index($0, pre) == 1 {print substr($0, length(pre) + 1, 6) + 0}' | sort -n > "$havein"
dir_names "$OUT/images" | awk -v pre="images_${SAMPLE}_hijing-$(printf '%010d' "$RUN")-" 'index($0, pre) == 1 {print substr($0, length(pre) + 1, 6) + 0}' | sort -n > "$haveout"
[[ $ONLYMISSING -eq 1 ]] || : > "$haveout"
awk 'FILENAME == ARGV[1] {inp[$1] = 1; next} FILENAME == ARGV[2] {out[$1] = 1; next}
     !($1 in inp) {print "MISSING_IN"; next} !($1 in out)' "$havein" "$haveout" "$all" > "$all.todo"
nmissing_in=$(grep -c '^MISSING_IN' "$all.todo" || true)
grep -v '^MISSING_IN' "$all.todo" > "$all.queue" || true
while read -r seg nev first; do
  if (( nseg % PERJOB == 0 )); then
    chunk="$OUT/condor/chunks/chunk_$(printf '%05d' $k).list"
    echo "$chunk" >> "$OUT/condor/overlay_queue.list"
    k=$((k + 1))
  fi
  echo "$seg $nev $first" >> "$chunk"
  nseg=$((nseg + 1))
done < "$all.queue"

jobfile="$OUT/condor/overlay.job"
cat > "$jobfile" <<EOF
universe            = vanilla
executable          = $F4A_ROOT/scripts/run_overlay.sh
transfer_executable = False
arguments           = \$(chunk) $SAMPLE $RUN $IN/images $OUT/images
output              = $OUT/logs/overlay_\$Fn(chunk).out
error               = $OUT/logs/overlay_\$Fn(chunk).err
log                 = /tmp/${USER}_f4a_overlay_${TAG}_${SAMPLE}.log
request_memory      = 2000MB
priority            = $PRIO
PeriodicHold        = (NumJobStarts >= 1 && JobStatus == 1)
notification        = Never
queue chunk from $OUT/condor/overlay_queue.list
EOF

cat >> "$OUT/production.txt" <<EOF
# ---- overlay $(date '+%F %T') on $(hostname)$([[ $DRYRUN -eq 1 ]] && echo " (prepared, not submitted)")
$SAMPLE + HIJING 0-10% ($NHIJING background images, no reuse): $nseg segments in $k jobs;
$nmissing_in segments skipped because the JEWEL image file does not exist (yet)
rerun one  scripts/run_overlay.sh <chunk list> $SAMPLE $RUN $IN/images <outdir>
EOF

echo ">> ${SAMPLE}_hijing / $TAG: $nseg segments in $k jobs ($nmissing_in segments without JEWEL images yet)"
echo ">> output: $OUT/images"
if [[ $DRYRUN -eq 1 ]]; then echo ">> dry run: not submitted"; exit 0; fi
(( k > 0 )) || { echo ">> nothing to submit"; exit 0; }
condor_submit "$jobfile"
