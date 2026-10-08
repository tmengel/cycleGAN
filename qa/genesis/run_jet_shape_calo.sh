#!/bin/bash
# Runs JetShapeCalo.C: radial profiles <dpT/dDeltaR> and rho(DeltaR) of inclusive reco jets,
# total and split into EMCal and HCal towers, plus the EMCal/HCal make-up of rho, for the two
# samples and their ratio. Needs trees carrying rjet_constit_src, i.e. trees_src/merged.
set -euo pipefail

# source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.572

TREES=/sphenix/tg/tg01/jets/tmengel/production/JEWEL_pp200_signal_genesis/trees_src/merged
OUTDIR=/sphenix/user/tmengel/cycleGAN/qa/genesis/jet_shape_calo
MACRO=/sphenix/user/tmengel/cycleGAN/qa/genesis/JetShapeCalo.C

# the pT window is a per-jet cut on RECO pT; DIJET_VETO=1 also requires a back-to-back truth dijet
DIJET_VETO=0
PTMIN=20
PTMAX=30
NEVENTS=-1

VACUUM="vacuum|Vacuum|$TREES/vacuum_trees_merged.root"
NORECOIL="norecoil|Medium w/o recoil|$TREES/norecoil_trees_merged.root"

root -l -b -q "$MACRO+(\"$VACUUM\", \"$NORECOIL\", \"$OUTDIR\", \"JEWEL 2.6.0, #sqrt{s_{NN}} = 200 GeV\", $NEVENTS, $PTMIN, $PTMAX, $DIJET_VETO)"
