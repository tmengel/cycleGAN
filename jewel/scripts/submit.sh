#!/bin/bash
# Prepares (and by default submits) a JEWEL production on condor: <njobs> jobs
# of <nevent> events each. The JEWEL seed is the job number: job k runs with
# NJOB = k (seed k*1000), k = first_job ... first_job+njobs-1 (first_job = 1
# unless --first-job is given). Rerunning job k with the same config and
# <nevent> reproduces its output byte for byte.
#
#   scripts/submit.sh <config> <tag> <njobs> <nevent> [--first-job N]
#                     [--filter R ptmin etamax] [--dry-run] [--reproduce]
#
#   --first-job N  extend a production: jobs N ... N+njobs-1
#   --dry-run      write the job files and manifest but don't submit or touch
#                  the ledger (rerun without it to submit)
#   --reproduce    deliberately regenerate jobs that seeds.ledger says were
#                  already produced with identical physics (exact rerun)
#
# Output goes to $JEWEL_OUTPUT_ROOT/<tag>/<config>/ (default root
# /sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc):
#   raw/             one standalone .hepmc per job: <config>_<job, 6 digits>.hepmc
#   filtered/        same, after the jet filter (only with --filter)
#   logs/            JEWEL params + log per job
#   config/          snapshot of the config (params, medium params, tables md5)
#   production.txt   everything needed to reproduce the production
#   condor/          job file, job list, stdout/stderr
#
# seeds.ledger records every submitted job range with the config's physics
# fingerprint. A range that overlaps one already produced by the same config
# with the same fingerprint would duplicate events, so it is refused (unless
# --reproduce). Different configs, or the same config with different physics,
# may reuse job numbers: vacuum and medium jobs with the same NJOB produce
# different events.
#
# Example (100 jobs x 10k events per config, seeds 1-100):
#   scripts/submit.sh vacuum pthat5-10 100 10000
#   scripts/submit.sh medium pthat5-10 100 10000

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/common.sh"

usage="usage: $0 <config> <tag> <njobs> <nevent> [--first-job N] [--filter R ptmin etamax] [--dry-run] [--reproduce]"
[[ $# -ge 4 ]] || die "$usage"
load_config "$1"
TAG=$2; NJOBS=$3; NEVENT=$4
shift 4
FIRST=1; FILTER=""; DRYRUN=0; REPRODUCE=0
while [[ $# -gt 0 ]]; do
  case $1 in
    --first-job) [[ $# -ge 2 ]] || die "$usage"; FIRST=$2; shift 2 ;;
    --filter)    [[ $# -ge 4 ]] || die "$usage"; FILTER="$2 $3 $4"; shift 4 ;;
    --dry-run)   DRYRUN=1; shift ;;
    --reproduce) REPRODUCE=1; shift ;;
    *) die "$usage" ;;
  esac
done
[[ "$NJOBS" =~ ^[0-9]+$ && "$NEVENT" =~ ^[0-9]+$ ]] || die "$usage"
LAST=$((FIRST + NJOBS - 1))
check_njob "$FIRST"; check_njob "$LAST"
check_tables
FP=$(config_fingerprint)

OUTDIR="$JEWEL_OUTPUT_ROOT/$TAG/$CONFIG_NAME"
LEDGER="$JEWEL_ROOT/seeds.ledger"

# --- seed bookkeeping --------------------------------------------------------
touch "$LEDGER"
if [[ $REPRODUCE -eq 0 ]]; then
  while read -r lo hi cfg tag _ rest; do
    [[ -z "$lo" || "$lo" == \#* ]] && continue
    [[ "$cfg" == "$CONFIG_NAME" && " $rest " == *" fp=$FP "* ]] || continue
    if (( FIRST <= hi && LAST >= lo )); then
      die "jobs $FIRST-$LAST overlap jobs $lo-$hi of $tag/$cfg, produced with identical physics (see $LEDGER).
       Use --first-job $((hi + 1)) to extend that production, or --reproduce to regenerate it exactly."
    fi
  done < "$LEDGER"
fi

# --- job files ---------------------------------------------------------------
mkdir -p "$OUTDIR"/{config,condor/out}
cp "$CONFIG_DIR"/config.sh "$CONFIG_DIR"/params.tmpl.dat "$CONFIG_DIR"/tables/fingerprint "$OUTDIR/config/"
[[ -f "$CONFIG_DIR/medium.params.dat" ]] && cp "$CONFIG_DIR/medium.params.dat" "$OUTDIR/config/"
(cd "$CONFIG_DIR/tables" && md5sum xsecs.dat pdfs.dat splitint.dat) > "$OUTDIR/config/tables.md5"
seq "$FIRST" "$LAST" > "$OUTDIR/condor/jobs_${FIRST}_${LAST}.list"

jobfile="$OUTDIR/condor/jewel_${FIRST}_${LAST}.job"
cat > "$jobfile" <<EOF
universe            = vanilla
executable          = $JEWEL_ROOT/scripts/run_jewel.sh
transfer_executable = False
arguments           = $CONFIG_NAME \$(job) $NEVENT $OUTDIR $FILTER
output              = $OUTDIR/condor/out/\$(job).out
error               = $OUTDIR/condor/out/\$(job).err
log                 = /tmp/${USER}_jewel_${TAG}_${CONFIG_NAME}.log
request_memory      = 2GB
notification        = Never
queue job from $OUTDIR/condor/jobs_${FIRST}_${LAST}.list
EOF

# Manifest: what to rerun, and what has to match for the rerun to be exact.
cat >> "$OUTDIR/production.txt" <<EOF
# ---- $(date '+%F %T') on $(hostname)$([[ $DRYRUN -eq 1 ]] && echo " (prepared, not yet submitted)")
command        scripts/submit.sh $CONFIG_NAME $TAG $NJOBS $NEVENT --first-job $FIRST${FILTER:+ --filter $FILTER}
jobs (= NJOB)  $FIRST-$LAST, $NEVENT events each; job k -> raw/${CONFIG_NAME}_<k, 6 digits>.hepmc
rerun one job  scripts/run_jewel.sh $CONFIG_NAME <k> $NEVENT <outdir>${FILTER:+ $FILTER}
physics        fingerprint $FP (config/params.tmpl.dat$([[ -f "$CONFIG_DIR/medium.params.dat" ]] && echo ", config/medium.params.dat"))
tables         config/tables.md5 (rebuild with scripts/make_tables.sh $CONFIG_NAME)
binary         $(md5sum < "$JEWEL_EXE" | cut -d' ' -f1)  $(readlink -f "$JEWEL_EXE")
environment    sPHENIX $JEWEL_SPHENIX_RELEASE, LHAPDF $LHAPDF6_PREFIX
pdf data path  $LHAPDF_DATA_PATH
EOF

echo ">> $CONFIG_NAME / $TAG: $NJOBS jobs x $NEVENT events, seed = job number $FIRST-$LAST${FILTER:+, filter (R ptmin etamax) = $FILTER}"
echo ">> output: $OUTDIR"
echo ">> job file: $jobfile"
if [[ $DRYRUN -eq 1 ]]; then
  echo ">> dry run: not submitted, ledger not updated. To submit, rerun without --dry-run."
  exit 0
fi
condor_submit "$jobfile"
echo "$FIRST $LAST $CONFIG_NAME $TAG $(date +%F) nevent=$NEVENT fp=$FP" >> "$LEDGER"
