#!/bin/bash
# Runs JetCenteredMapsCalo.C over the genesis100k merged trees: inclusive reco jet maps split
# into EMCal and HCal towers, with the two samples and their ratio, the same sample specs
# condor_qa.job uses for JetCenteredMaps.C.
#
# The calo split reads jettree's rjet_constit_src branch, the Jet::SRC each tower constituent
# came from. Only trees written by a JewelTreeWriter carrying that branch can be split; against
# older trees the macro says so and stops. TREES therefore has to point at a production made
# after that branch was added.
set -euo pipefail

# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.572

TREES=/sphenix/tg/tg01/jets/tmengel/production/JEWEL_pp200_signal_genesis/trees_src/merged
OUTDIR=/sphenix/user/tmengel/cycleGAN/qa/genesis/jet_maps_calo
MACRO=/sphenix/user/tmengel/cycleGAN/qa/genesis/JetCenteredMapsCalo.C

# DIJET_VETO=1 also requires a back-to-back truth dijet (|dphi_12| > 3pi/4). The pT window is a
# per-jet cut on RECO pT, not an event veto -- this macro is inclusive-only. Different selections
# overwrite each other's output, so change OUTDIR when scanning pT bins.
DIJET_VETO=0
PT1MIN=20
PT1MAX=30
NEVENTS=-1

VACUUM="vacuum|Vacuum|$TREES/vacuum_trees_merged.root"
NORECOIL="norecoil|TESST w/o recoil|$TREES/norecoil_trees_merged.root"

root -l -b -q "$MACRO+(\"$NORECOIL\", \"$VACUUM\", \"$OUTDIR\", \"JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV\", $NEVENTS, $PT1MIN, $PT1MAX, $DIJET_VETO)"
