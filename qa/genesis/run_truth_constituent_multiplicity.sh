#!/bin/bash
# Truth-jet constituent multiplicity, vacuum vs no-recoil, in jet pT bins, over the
# v2_2026-09-14_genesis100k merged trees. Same sample specs as condor_qa.job.
set -euo pipefail

# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.572

TAG=/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v2_2026-09-14_genesis100k
OUTDIR=/sphenix/user/tmengel/cycleGAN/qa/genesis
MACRO=/sphenix/user/tmengel/cycleGAN/qa/genesis/TruthConstituentMultiplicity.C

VACUUM="vacuum|vacuum|$TAG/trees/merged/vacuum_trees_merged.root|$TAG/truthqa/truthqa_vacuum.root"
NORECOIL="norecoil|no-recoil (T_{i}=375 MeV)|$TAG/trees/merged/norecoil_trees_merged.root|$TAG/truthqa/truthqa_norecoil.root"

# Bins below 20 GeV are biased by the generator-level jet filter (see the macro header).
NEVENTS=${NEVENTS:--1}
PTBINS=${PTBINS:-10,15,20,30,40,60}

root -l -b -q "$MACRO(\"$VACUUM\", \"$NORECOIL\", \"$OUTDIR\", \"JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV\", $NEVENTS, \"$PTBINS\")"
