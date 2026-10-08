#!/bin/bash
# Prepares (and by default submits) pass 1 on condor for one sample.
#
#   scripts/submit_pass1.sh [options] <sample> <tag> <runnumber> <events_per_job> <max_events> <hepmc>...
#
#   sample          name used in file names (e.g. vacuum, medium)
#   tag             production name; output goes to $F4A_OUTPUT_ROOT/<tag>/<sample>/
#                   (default root /sphenix/tg/tg01/jets/tmengel/JEWEL_sim)
#   runnumber       100-9999, one per sample (convention: vacuum 1000, medium 2000); with the
#                   segment it fixes every job's random seeds
#   events_per_job  10 recommended (~27 s/event GEANT4 + ~3 min startup; ~22 MB DST per job)
#   max_events      stop after this many events in total (0 = all events in the files)
#   hepmc...        HepMC files, used in version-sorted order (e.g. a JEWEL production's raw/*.hepmc)
#
# options:
#   --beam B        fixed (default: vertex at 0,0,0), AuAu, pp or pp_zeroangle -- see the pass 1 macro
#   --only-missing  queue only the segments whose G4Hits DST does not exist yet (resubmit failures)
#   --extend        grow an existing production to more events: the existing job list is kept
#                   unchanged (its segments keep their events, DSTs and seeds) and new segments are
#                   appended, starting right after its last event. Combine with --only-missing to
#                   queue just the new segments.
#   --dry-run       write the job list and job file, but don't submit
#
# Job k (= segment k) simulates a fixed block of events [skip, skip + n) of one HepMC file; the
# list is written to condor/pass1_jobs.list and never changes for a production, so segment k
# always means the same events.
#
# Example (first 100k events of each JEWEL sample, 10 events/job -> ~10k jobs per sample):
#   H=/sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc/pthat16
#   scripts/submit_pass1.sh vacuum pthat16 1000 10 100000 $H/vacuum/raw/*.hepmc
#   scripts/submit_pass1.sh medium pthat16 2000 10 100000 $H/medium/raw/*.hepmc

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

BEAM=fixed; ONLYMISSING=0; DRYRUN=0; EXTEND=0
while [[ $# -gt 0 && "$1" == --* ]]; do
  case $1 in
    --beam) BEAM=$2; shift 2 ;;
    --only-missing) ONLYMISSING=1; shift ;;
    --extend) EXTEND=1; shift ;;
    --dry-run) DRYRUN=1; shift ;;
    *) die "unknown option $1" ;;
  esac
done
[[ $# -ge 6 ]] || die "usage: $0 [--beam B] [--only-missing] [--extend] [--dry-run] <sample> <tag> <runnumber> <events_per_job> <max_events> <hepmc>..."
SAMPLE=$1; TAG=$2; RUN=$3; NPER=$4; MAX=$5; shift 5
check_run_segment "$RUN" 0
[[ "$NPER" =~ ^[0-9]+$ && $NPER -gt 0 && "$MAX" =~ ^[0-9]+$ ]] || die "bad events_per_job/max_events"
case $BEAM in fixed|AuAu|pp|pp_zeroangle) ;; *) die "--beam must be fixed, AuAu, pp or pp_zeroangle" ;; esac

OUT="$F4A_OUTPUT_ROOT/$TAG/$SAMPLE"
mkdir -p "$OUT"/{condor,logs}

# --- job list: segment, hepmc, skip, nevents ----------------------------------------------
# Blocks are aligned to multiples of events_per_job within each file; a block is shorter only at
# the end of a file or where max_events cuts it.
list="$OUT/condor/pass1_jobs.list"
new=$(mktemp)
trap 'rm -f "$new"' EXIT
mapfile -t files < <(for f in "$@"; do readlink -f "$f"; done | sort -V)
seg=0; total=0; startfile=""; startpos=0
if [[ $EXTEND -eq 1 ]]; then
  [[ -f "$list" ]] || die "--extend needs an existing $list"
  cp "$list" "$new"
  read -r lseg startfile lskip ln < <(tail -1 "$list")
  seg=$((lseg + 1)); startpos=$((lskip + ln))
  total=$(awk '{s += $4} END {print s + 0}' "$list")
  # the existing list must cover a leading run of the given files, in the same order
  mapfile -t oldfiles < <(awk '!seen[$2]++ {print $2}' "$list")
  for i in "${!oldfiles[@]}"; do
    [[ "${files[$i]:-}" == "${oldfiles[$i]}" ]] || die "--extend: the given files do not start with the existing list's files"
  done
fi
started=$([[ -z "$startfile" ]] && echo 1 || echo 0)
for f in "${files[@]}"; do
  pos=0
  if [[ $started -eq 0 ]]; then
    [[ "$f" == "$startfile" ]] || continue
    started=1; pos=$startpos
  fi
  [[ -f "$f" ]] || die "no such file: $f"
  tail -c 4096 "$f" | grep -q END_EVENT_LISTING || die "$f has no END_EVENT_LISTING footer (truncated?)"
  n=$(grep -c '^E ' "$f")
  for (( skip = pos; skip < n; skip += take )); do
    next=$(( (skip / NPER + 1) * NPER ))
    take=$(( (next < n ? next : n) - skip ))
    (( MAX > 0 && total + take > MAX )) && take=$(( MAX - total ))
    (( take > 0 )) || break
    echo "$seg $f $skip $take" >> "$new"
    seg=$((seg + 1)); total=$((total + take))
  done
  (( MAX > 0 && total >= MAX )) && break
done
(( MAX == 0 || total == MAX )) || die "only $total events in the given files, asked for $MAX"
(( seg <= 1000000 )) || die "$seg jobs: segment numbers would exceed 999999; use more events per job"

if [[ -f "$list" ]] && ! cmp -s "$list" "$new"; then
  if [[ $EXTEND -eq 1 ]]; then
    cp "$list" "$list.before_extend.$(date +%Y%m%d%H%M%S)"
    cp "$new" "$list"
  else
    die "$list exists and differs from this request: a production's job list must never change
       (segment k = fixed events). Use --extend to add events, a new tag, or delete the old outputs and list."
  fi
elif [[ ! -f "$list" ]]; then
  cp "$new" "$list"
fi
chmod 644 "$list" "$OUT"/condor/pass1_jobs.list.before_extend.* 2>/dev/null || true

queue="$OUT/condor/pass1_queue.list"
if [[ $ONLYMISSING -eq 1 ]]; then
  dir_names "$OUT/G4Hits" | awk -v pre="G4Hits_${SAMPLE}-$(printf '%010d' "$RUN")-" '
      NR == FNR { if (index($0, pre) == 1) have[substr($0, length(pre) + 1, 6) + 0] = 1; next }
      !($1 in have)' - "$list" > "$queue"
else
  cp "$list" "$queue"
fi
nq=$(wc -l < "$queue")

jobfile="$OUT/condor/pass1.job"
cat > "$jobfile" <<EOF
universe            = vanilla
executable          = $F4A_ROOT/scripts/run_pass1.sh
transfer_executable = False
arguments           = $SAMPLE $RUN \$(segment) \$(hepmc) \$(skip) \$(nev) $OUT $BEAM
output              = $OUT/logs/pass1_\$(segment).out
error               = $OUT/logs/pass1_\$(segment).err
log                 = /tmp/${USER}_f4a_pass1_${TAG}_${SAMPLE}.log
request_memory      = 6000MB
PeriodicHold        = (NumJobStarts >= 1 && JobStatus == 1)
notification        = Never
queue segment, hepmc, skip, nev from $queue
EOF

cat >> "$OUT/production.txt" <<EOF
# ---- pass1 $(date '+%F %T') on $(hostname)$([[ $DRYRUN -eq 1 ]] && echo " (prepared, not submitted)")
sample $SAMPLE, run $RUN, beam $BEAM, $seg jobs x <= $NPER events = $total events, $nq queued
job list   condor/pass1_jobs.list (segment hepmc skip nevents)
rerun one  scripts/run_pass1.sh $SAMPLE $RUN <segment> <hepmc> <skip> <nevents> <outdir> $BEAM
seeds      pass1 = run*100000 + segment
release    $F4A_RELEASE ($OFFLINE_MAIN), CDB ${F4A_CDBTAG:-MDC2_ana.435}
EOF

echo ">> $SAMPLE / $TAG: $seg jobs, $total events (run $RUN, beam $BEAM), $nq to queue"
echo ">> output: $OUT"
echo ">> job file: $jobfile"
if [[ $DRYRUN -eq 1 ]]; then
  echo ">> dry run: not submitted"
  exit 0
fi
(( nq > 0 )) || { echo ">> nothing to submit"; exit 0; }
condor_submit "$jobfile"
