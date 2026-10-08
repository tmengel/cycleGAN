#!/bin/bash
# Re-runs JetCenteredMaps.C over the v2_2026-09-14_genesis100k merged trees, with the same
# sample specs used by condor_qa.job, to refresh the maps under qa/jet_maps after macro edits
# (plotting/normalization only — inputs unchanged).
set -euo pipefail

# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.572

TAG=/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v2_2026-09-14_genesis100k
OUTDIR=/sphenix/user/tmengel/cycleGAN/qa/genesis
MACRO=/sphenix/user/tmengel/cycleGAN/qa/genesis/JetCenteredMaps.C

# USE_Y=1 builds the truth maps in Delta-y, which needs the tjet_constit_e branch that only the
# regenerated trees carry, so it also switches TREES to trees_pid/merged.
USE_Y=0
NORM_TO_LEAD=0
DIJET_VETO=0
INCLUSIVE=1
PT1MIN=20
PT1MAX=30
NEVENTS=-1

USE_Y=${USE_Y:-0}
if [ "$USE_Y" = "1" ]; then
  TREES=/sphenix/tg/tg01/jets/tmengel/production/JEWEL_pp200_signal_genesis/trees_pid/merged
else
  TREES=$TAG/trees/merged
fi

VACUUM="vacuum|Vacuum|$TREES/vacuum_trees_merged.root|$TAG/truthqa/truthqa_vacuum.root"
NORECOIL="norecoil|Medium w/o recoil|$TREES/norecoil_trees_merged.root|$TAG/truthqa/truthqa_norecoil.root"

# Event selection on the leading truth jet, applied to every map. PT1MAX <= 0 drops the upper
# edge; DIJET_VETO=1 also requires a back-to-back truth dijet (|dphi_12| > 3pi/4). Different
# selections overwrite each other's output, so change OUTDIR when scanning pT bins.
# NORM_TO_LEAD=1 normalises all four maps of a sample by its accepted leading truth jets
# (one per selected event) instead of by the jets in each map.
# INCLUSIVE=1 drops the leading/subleading split: every fiducial jet whose OWN pT is inside
# [PT1MIN, PT1MAX) goes into one map per level, and the pT window stops being an event veto.
# That makes it 2 figures instead of 4, named *_inclusive, and NORM_TO_LEAD then means "per
# truth jet in the window".


root -l -b -q "$MACRO(\"$NORECOIL\", \"$VACUUM\", \"$OUTDIR\", \"JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV\", $NEVENTS, $PT1MIN, $PT1MAX, $DIJET_VETO, $NORM_TO_LEAD, $USE_Y, $INCLUSIVE)"
