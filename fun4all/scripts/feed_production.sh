#!/bin/bash
# Keeps a production moving without babysitting: every 10 minutes it
#   1. tops up pass 2 for each sample to the scheduler's per-user cap (100,000 queued jobs),
#      never submitting a DST twice (submit_pass2.sh --skip-submitted),
#   2. once a sample's pass 2 outputs are complete, submits its HIJING overlay (all events),
#   3. resubmits failed jobs (outputs missing, nothing left in the queue) up to 2 times per stage.
# It stops when every overlay is complete, or after 7 days.
#
#   nohup scripts/feed_production.sh <tag> vacuum:1000 medium:2000 > feed.log 2>&1 &

set -uo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

TAG=$1; shift
SAMPLES=("$@")
CAP=98000          # stay a little under the 100,000-job limit
SLEEP=600
NHIJING=$(awk '$2 > 0 {s += $2} END {print s + 0}' "$F4A_ROOT/hijing/type4_run19_hijingAll_noNoise_0to10_cent0.counts")
log() { echo "[$(date '+%F %T')] $*"; }

queued_total() { timeout 300 condor_q "$USER" -totals 2>/dev/null | awk '/Total for query/ {print $4; exit}'; }
queued_for() { timeout 300 condor_q "$USER" -af Args 2>/dev/null | grep -c -- "$1" || true; }
count_root() { dir_names "$1" | grep -c '\.root$' || true; }

declare -A tries
START=$SECONDS
while (( SECONDS - START < 7 * 86400 )); do
  alldone=1
  Q=$(queued_total); Q=${Q:-$CAP}
  for s in "${SAMPLES[@]}"; do
    name=${s%%:*}; run=${s##*:}
    OUT="$F4A_OUTPUT_ROOT/$TAG/$name"; OVL="$F4A_OUTPUT_ROOT/$TAG/${name}_hijing"
    ndst=$(count_root "$OUT/G4Hits"); ntrees=$(count_root "$OUT/trees"); nimg=$(count_root "$OUT/images")
    log "$name: DSTs $ndst, trees $ntrees, images $nimg; queue total $Q"

    # ---- pass 2
    if (( ntrees < ndst || nimg < ndst )); then
      alldone=0
      room=$(( CAP - Q ))
      inq=$(queued_for "$OUT/G4Hits/")
      nsub=$(wc -l < "$OUT/condor/pass2_submitted.list" 2>/dev/null || echo 0)
      if (( nsub < ndst && room > 1000 )); then
        log "$name: submitting up to $room pass2 jobs"
        "$F4A_SCRIPTS/submit_pass2.sh" --only-missing --skip-submitted --max-jobs "$room" "$name" "$TAG" 2>&1 | grep -E "^>>|submitted|ERROR"
        Q=$(queued_total); Q=${Q:-$CAP}
      elif (( nsub >= ndst && inq == 0 )); then
        k="pass2_$name"; tries[$k]=$(( ${tries[$k]:-0} + 1 ))
        if (( tries[$k] > 2 )); then log "$name: pass2 still incomplete after 2 retries -- giving up on it"; continue; fi
        log "$name: pass2 queue empty but $(( ndst - ntrees )) outputs missing: retry ${tries[$k]}"
        "$F4A_SCRIPTS/submit_pass2.sh" --only-missing --max-jobs $(( room > 0 ? room : 1 )) "$name" "$TAG" 2>&1 | grep -E "^>>|submitted|ERROR"
        Q=$(queued_total); Q=${Q:-$CAP}
      fi
      continue
    fi

    # ---- overlay (pass 2 complete)
    expect=$(awk -v nh="$NHIJING" '{ if (e + 0 < nh) n++; e += $4 } END {print n + 0}' "$OUT/condor/pass1_jobs.list")
    novl=$(count_root "$OVL/images")
    log "${name}_hijing: images $novl / $expect"
    if (( novl < expect )); then
      alldone=0
      if [[ ! -f "$OVL/condor/overlay_full.launched" ]]; then
        log "${name}_hijing: submitting overlay for all events"
        "$F4A_SCRIPTS/submit_overlay.sh" --priority 10 --only-missing "$name" "$TAG" "$run" 2>&1 | grep -E "^>>|submitted|ERROR"
        touch "$OVL/condor/overlay_full.launched"
      elif (( $(queued_for "$OVL/images") == 0 )); then
        k="ovl_$name"; tries[$k]=$(( ${tries[$k]:-0} + 1 ))
        if (( tries[$k] > 2 )); then log "${name}_hijing: still incomplete after 2 retries -- giving up on it"; continue; fi
        log "${name}_hijing: queue empty but $(( expect - novl )) missing: retry ${tries[$k]}"
        "$F4A_SCRIPTS/submit_overlay.sh" --priority 10 --only-missing "$name" "$TAG" "$run" 2>&1 | grep -E "^>>|submitted|ERROR"
      fi
    fi
  done
  if (( alldone )); then log "all samples and overlays complete"; break; fi
  sleep $SLEEP
done
